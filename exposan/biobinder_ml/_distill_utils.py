# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""
Distillation runtime chooser utilities with VLE + design/cost sanity guardrails.
"""

import numpy as np

__all__ = (
    "safe_LHK_candidates",
    "find_best_Lr_Hr_runtime",
    "add_runtime_LHK_LrHr_spec",
)

def _vle_check(col, x_min=1e-12, bp_T_max=950.0, alpha_max=1e6, K_max=1e8):
    IDs = getattr(col, "_IDs_vle", None)
    LHK_vle = getattr(col, "_LHK_vle_index", None)
    dew = getattr(col, "_dew_point", None)
    bub = getattr(col, "_bubble_point", None)
    if (not IDs) or (LHK_vle is None) or (dew is None) or (bub is None):
        return False, {"why": "missing_vle_state"}

    Dstream, Bstream = col.outs
    zD = Dstream.get_normalized_mol(IDs)
    zB = Bstream.get_normalized_mol(IDs)

    try:
        if getattr(col, "_partial_condenser", True):
            dp = dew(zD, P=col.P)
        else:
            dp = bub(zD, P=col.P)
        bp = bub(zB, P=col.P)
    except Exception as e:
        return False, {"why": f"dp_bp_fail:{e}"}

    dp_min_x = float(np.min(dp.x))
    bp_min_x = float(np.min(bp.x))
    bp_T = float(bp.T)
    dp_T = float(dp.T)

    try:
        alpha_mean = col._estimate_mean_volatilities_relative_to_heavy_key()
        LK_vle = int(LHK_vle[0])
        alpha_LK = float(alpha_mean[LK_vle])
    except Exception as e:
        return False, {"why": f"alpha_fail:{e}"}

    # >>> ADD THIS BLOCK <<<
    eps = 1e-16
    Kdist = dp.y / np.maximum(dp.x, eps)
    Kbot  = bp.y / np.maximum(bp.x, eps)
    Kmax_obs = float(np.nanmax([np.nanmax(Kdist), np.nanmax(Kbot)]))

    ok = True
    why = []
    if (not np.isfinite(bp_T)) or (not np.isfinite(dp_T)):
        ok = False; why.append("nonfinite_T")
    if bp_T > bp_T_max:
        ok = False; why.append(f"bp_T>{bp_T_max}")
    if dp_min_x < x_min:
        ok = False; why.append(f"dp_min_x<{x_min}")
    if bp_min_x < x_min:
        ok = False; why.append(f"bp_min_x<{x_min}")
    if (not np.isfinite(alpha_LK)) or (alpha_LK <= 1.0) or (alpha_LK > alpha_max):
        ok = False; why.append(f"alpha_LK_bad:{alpha_LK}")
    if (not np.isfinite(Kmax_obs)) or (Kmax_obs > K_max):
        ok = False; why.append(f"Kmax>{K_max}:{Kmax_obs}")

    return ok, {
        "dp_T": dp_T,
        "bp_T": bp_T,
        "dp_min_x": dp_min_x,
        "bp_min_x": bp_min_x,
        "alpha_LK": alpha_LK,
        "Kmax": Kmax_obs,
        "why": "|".join(why) if why else "",
    }


def _design_cost_sanity(
    col,
    max_purchase_cost=1e9,
    allow_missing_purchase_costs=True,
):
    D = getattr(col, "design_results", None) or {}
    C = getattr(col, "baseline_purchase_costs", None) or {}

    # Names vary by BioSTEAM version/unit; check both common variants
    N = float(D.get("Theoretical stages", np.nan))
    H = float(D.get("Rectifier height", D.get("Height", np.nan)))
    Diam = float(D.get("Rectifier diameter", D.get("Diameter", np.nan)))
    W = float(D.get("Rectifier weight", D.get("Weight", np.nan)))
    tw = float(D.get("Rectifier wall thickness", D.get("Wall thickness", np.nan)))

    bad = []
    for k, v in (("N", N), ("H", H), ("Diam", Diam), ("W", W), ("tw", tw)):
        if (not np.isfinite(v)) or (v <= 0):
            bad.append(f"{k}_bad:{v}")

    # Purchase costs sanity (optional but very helpful)
    if C:
        try:
            tot = float(np.nansum(list(C.values())))
        except Exception:
            tot = np.nan
        if (not np.isfinite(tot)) or (tot <= 0) or (tot > max_purchase_cost):
            bad.append(f"purchase_cost_bad:{tot}")
    else:
        if not allow_missing_purchase_costs:
            bad.append("missing_purchase_costs")

    return (len(bad) == 0), "|".join(bad)

def safe_LHK_candidates(
    splitter,
    col_feed,
    idx,
    fallback=("4M-PHYNO", "INDOLE"),
    min_massfrac=1e-6,  #1e-6 for food/sludge, -3 for green/manure
):
    def massfrac(stream, cid):
        try:
            tot = float(stream.F_mass)
            if tot <= 0:
                return 0.0
            return float(stream.imass[cid]) / tot
        except Exception:
            return 0.0

    cands = []
    try:
        ks = splitter.keys
        if ks and len(ks) > idx and ks[idx]:
            lk, hk = ks[idx]
            if massfrac(col_feed, lk) > min_massfrac and massfrac(col_feed, hk) > min_massfrac:
                cands.append(ks[idx])
    except Exception:
        pass

    lk, hk = fallback
    if massfrac(col_feed, lk) > min_massfrac and massfrac(col_feed, hk) > min_massfrac:
        cands.append(fallback)
    else:
        cands.append(fallback)

    out, seen = [], set()
    for c in cands:
        if c not in seen:
            out.append(c)
            seen.add(c)
    return out

def find_best_Lr_Hr_runtime(
    col,
    target_ratio=None,
    tol=0.20,
    Lr_grid=None,
    Hr_grid=None,
    prefer="max_Lr_then_max_Hr",
    require_design=True,
    require_cost=False,
    vle_screen=True,
    x_min=1e-12,
    bp_T_max=950.0,
    alpha_max=1e6,
    design_sanity=True,
    max_purchase_cost=1e9,
):
    if Lr_grid is None:
        Lr_grid = [0.995, 0.99, 0.985, 0.98, 0.975, 0.97, 0.965, 0.96, 0.95, 0.93, 0.90]
    if Hr_grid is None:
        Hr_grid = [0.99, 0.98, 0.97, 0.95, 0.93, 0.90, 0.89, 0.87, 0.85]

    def top_ratio():
        tot = float(col.F_mass_out)
        if tot > 0:
            return float(col.outs[0].F_mass) / tot
        return np.nan

    feasible = []
    for Hr in Hr_grid:
        for Lr in Lr_grid:
            col.reset_cache()
            col.Hr = float(Hr)
            col.Lr = float(Lr)
            try:
                col._run()

                # Design first (your failures often happen here)
                if require_design:
                    col._design()
                    ok_d, why_d = _design_cost_sanity(col, max_purchase_cost=1e9)
                    if not ok_d:
                        continue

                # VLE screen (safe to apply after _design; dp/bp objects should exist)
                if vle_screen:
                    ok_vle, diag = _vle_check(col, x_min=x_min, bp_T_max=bp_T_max, alpha_max=alpha_max)
                    if not ok_vle:
                        continue

                # Design/cost sanity (prevents negative H/W/N, etc.)
                if design_sanity and require_design:
                    ok_dc, why_dc = _design_cost_sanity(col, max_purchase_cost=max_purchase_cost)
                    if not ok_dc:
                        continue

                if require_cost:
                    col._cost()
                    if design_sanity:
                        ok_dc, why_dc = _design_cost_sanity(col, max_purchase_cost=max_purchase_cost)
                        if not ok_dc:
                            continue

                r = top_ratio()
                if not np.isfinite(r):
                    continue

                if target_ratio is not None:
                    err = abs(r - float(target_ratio))
                    if err > float(tol):
                        continue
                else:
                    err = np.nan

                feasible.append((float(Hr), float(Lr), float(r), float(err) if np.isfinite(err) else np.nan))

            except Exception:
                try:
                    col.reset_cache()
                except Exception:
                    pass
                continue

    if not feasible:
        return None

    if prefer == "max_Hr_then_max_Lr":
        feasible.sort(key=lambda t: (t[0], t[1]), reverse=True)
    elif prefer == "max_Lr_then_max_Hr":
        feasible.sort(key=lambda t: (t[1], t[0]), reverse=True)
    elif prefer == "max_top_ratio":
        def q(v): return round(float(v), 6)
        feasible.sort(key=lambda t: (q(t[2]), t[0], t[1]), reverse=True)
    elif prefer == "min_err":
        feasible.sort(key=lambda t: (t[3], -t[0], -t[1]))
    else:
        raise ValueError(
            f"Unknown prefer='{prefer}'. Use one of: "
            "max_Hr_then_max_Lr, max_Lr_then_max_Hr, max_top_ratio, min_err"
        )

    best = feasible[0]
    print(f"[choose] prefer={prefer} -> Hr={best[0]}, Lr={best[1]}, top_ratio={best[2]}")
    return best

def add_runtime_LHK_LrHr_spec(
    splitter,
    col,
    idx,
    fallback_LHK=("4M-PHYNO", "INDOLE"),
    target_ratio=None,
    tol=0.20,
    Lr_grid=None,
    Hr_grid=None,
    prefer="max_top_ratio",
    require_design=True,
    require_cost=False,
    # Guardrails:
    vle_screen=True,
    x_min=1e-12,
    bp_T_max=950.0,
    alpha_max=1e6,
    design_sanity=True,
    max_purchase_cost=1e9,
):
    def spec():
        last_err = None

        feed = col.ins[0]
        if float(feed.F_mass) <= 0:
            if hasattr(col, "_runtime_best"):
                delattr(col, "_runtime_best")
            return

        for LHK in safe_LHK_candidates(splitter, feed, idx, fallback=fallback_LHK):
            try:
                col.LHK = LHK
                col.reset_cache()

                best = find_best_Lr_Hr_runtime(
                    col,
                    target_ratio=target_ratio,
                    tol=tol,
                    Lr_grid=Lr_grid,
                    Hr_grid=Hr_grid,
                    prefer=prefer,
                    require_design=require_design,
                    require_cost=require_cost,
                    vle_screen=vle_screen,
                    x_min=x_min,
                    bp_T_max=bp_T_max,
                    alpha_max=alpha_max,
                    design_sanity=design_sanity,
                    max_purchase_cost=max_purchase_cost,
                )
                if best is None:
                    raise RuntimeError("no feasible (Hr,Lr)")

                Hr, Lr, ratio, err = best
                col.Hr = float(Hr)
                col.Lr = float(Lr)

                # Final verification MUST re-apply guardrails (prevents rare cache/path differences)
                col.reset_cache()
                col._run()
                if require_design:
                    col._design()

                if vle_screen:
                    ok_vle, diag = _vle_check(col, x_min=x_min, bp_T_max=bp_T_max, alpha_max=alpha_max)
                    if not ok_vle:
                        raise RuntimeError(f"final_vle_fail:{diag.get('why','')}")

                if design_sanity and require_design:
                    ok_dc, why_dc = _design_cost_sanity(col, max_purchase_cost=max_purchase_cost)
                    if not ok_dc:
                        raise RuntimeError(f"final_design_fail:{why_dc}")

                if require_cost:
                    col._cost()
                    if design_sanity:
                        ok_dc, why_dc = _design_cost_sanity(col, max_purchase_cost=max_purchase_cost)
                        if not ok_dc:
                            raise RuntimeError(f"final_cost_fail:{why_dc}")

                col._runtime_best = {
                    "prefer": prefer,
                    "LHK": col.LHK,
                    "Hr": float(Hr),
                    "Lr": float(Lr),
                    "top_ratio": float(ratio),
                    "err": (float(err) if (err is not None and np.isfinite(err)) else np.nan),
                }
                return

            except Exception as e:
                last_err = e
                try:
                    col.reset_cache()
                except Exception:
                    pass
                continue

        if hasattr(col, "_runtime_best"):
            delattr(col, "_runtime_best")

        raise RuntimeError(f"No feasible LHK/Lr/Hr found. Last error: {last_err}")

    col.add_specification(spec)
    col.run_after_specifications = True

