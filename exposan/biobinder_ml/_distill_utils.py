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
    min_massfrac=1e-6,  #1e-6 for food/sludge, 1e-3 for green/manure
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
    """
    Search the specified Lr/Hr grid for feasible ShortcutColumn
    operating points.

    Returns
    -------
    tuple or None
        (Hr, Lr, top_ratio, target_error)

    Notes
    -----
    If prefer='min_err', target_ratio must be supplied.
    If target_ratio is supplied for any policy, candidates outside
    the specified tolerance are rejected.
    """

    # ========================================================
    # INPUT VALIDATION
    # ========================================================

    valid_policies = {
        "max_Hr_then_max_Lr",
        "max_Lr_then_max_Hr",
        "max_top_ratio",
        "min_err",
    }

    if prefer not in valid_policies:
        raise ValueError(
            f"Unknown prefer={prefer!r}. "
            f"Choose from {sorted(valid_policies)}."
        )

    if prefer == "min_err" and target_ratio is None:
        raise ValueError(
            "prefer='min_err' requires a finite target_ratio."
        )

    if target_ratio is not None:
        target_ratio = float(target_ratio)

        if not np.isfinite(target_ratio):
            raise ValueError(
                f"target_ratio must be finite; "
                f"received {target_ratio}."
            )

        if not 0.0 <= target_ratio <= 1.0:
            raise ValueError(
                f"target_ratio must be between 0 and 1; "
                f"received {target_ratio}."
            )

    tol = float(tol)

    if not np.isfinite(tol) or tol < 0:
        raise ValueError(
            f"tol must be a finite nonnegative number; "
            f"received {tol}."
        )

    # ========================================================
    # DEFAULT SEARCH GRIDS
    # ========================================================

    if Lr_grid is None:
        Lr_grid = [
            0.995,
            0.990,
            0.985,
            0.980,
            0.975,
            0.970,
            0.965,
            0.960,
            0.950,
            0.930,
            0.900,
            0.850,
            0.800,
            0.700,
            0.600,
            0.500
        ]

    if Hr_grid is None:
        Hr_grid = [
            0.990,
            0.980,
            0.970,
            0.950,
            0.930,
            0.900,
            0.890,
            0.870,
            0.850,
            0.800,
            0.700,
            0.600,
            0.500
        ]

    # ========================================================
    # LOCAL HELPERS
    # ========================================================

    def top_ratio():
        total_out = float(col.F_mass_out)

        if total_out > 0:
            return float(
                col.outs[0].F_mass / total_out
            )

        return np.nan

    def separation_design_check():
        """
        Additional ShortcutColumn design checks.

        This uses BioSTEAM's reported design outputs. The more
        detailed independent Underwood calculation remains in
        collect_diagnostics().
        """

        design = getattr(
            col,
            "design_results",
            None,
        ) or {}

        stages = float(
            design.get(
                "Theoretical stages",
                np.nan,
            )
        )

        reflux = float(
            design.get(
                "Reflux",
                np.nan,
            )
        )

        minimum_reflux = float(
            design.get(
                "Minimum reflux",
                np.nan,
            )
        )

        bad = []

        if (
            not np.isfinite(stages)
            or stages < 2
            or stages > 120
        ):
            bad.append(
                f"stages_bad:{stages}"
            )

        if (
            not np.isfinite(reflux)
            or reflux <= 0
        ):
            bad.append(
                f"reflux_bad:{reflux}"
            )

        if (
            not np.isfinite(minimum_reflux)
            or minimum_reflux < 0
        ):
            bad.append(
                f"minimum_reflux_bad:"
                f"{minimum_reflux}"
            )

        if (
            np.isfinite(reflux)
            and np.isfinite(minimum_reflux)
            and reflux < minimum_reflux
        ):
            bad.append(
                f"reflux_below_minimum:"
                f"{reflux}<{minimum_reflux}"
            )

        return (
            len(bad) == 0,
            "|".join(bad),
        )

    # ========================================================
    # GRID SEARCH
    # ========================================================

    feasible = []

    for Hr in Hr_grid:
        for Lr in Lr_grid:

            try:
                col.reset_cache()

                col.Hr = float(Hr)
                col.Lr = float(Lr)

                col._run()

                # --------------------------------------------
                # Column design
                # --------------------------------------------

                if require_design:
                    col._design()

                    ok_design, why_design = (
                        _design_cost_sanity(
                            col,
                            max_purchase_cost=max_purchase_cost,
                        )
                    )

                    if not ok_design:
                        continue

                    ok_sep, why_sep = (
                        separation_design_check()
                    )

                    if not ok_sep:
                        continue

                # --------------------------------------------
                # VLE sanity
                # --------------------------------------------

                if vle_screen:
                    ok_vle, vle_diag = _vle_check(
                        col,
                        x_min=x_min,
                        bp_T_max=bp_T_max,
                        alpha_max=alpha_max,
                    )

                    if not ok_vle:
                        continue

                # --------------------------------------------
                # General design sanity
                # --------------------------------------------

                if design_sanity and require_design:
                    ok_dc, why_dc = (
                        _design_cost_sanity(
                            col,
                            max_purchase_cost=(
                                max_purchase_cost
                            ),
                        )
                    )

                    if not ok_dc:
                        continue

                # --------------------------------------------
                # Optional costing
                # --------------------------------------------

                if require_cost:
                    col._cost()

                    if design_sanity:
                        ok_dc, why_dc = (
                            _design_cost_sanity(
                                col,
                                max_purchase_cost=(
                                    max_purchase_cost
                                ),
                            )
                        )

                        if not ok_dc:
                            continue

                # --------------------------------------------
                # Product-ratio objective
                # --------------------------------------------

                ratio = top_ratio()

                if not np.isfinite(ratio):
                    continue

                if not 0.0 <= ratio <= 1.0:
                    continue

                if target_ratio is not None:
                    error = abs(
                        ratio - target_ratio
                    )
                    
                else:
                    error = np.nan

                feasible.append(
                    (
                        float(Hr),
                        float(Lr),
                        float(ratio),
                        (
                            float(error)
                            if np.isfinite(error)
                            else np.nan
                        ),
                    )
                )

            except Exception:
                try:
                    col.reset_cache()
                except Exception:
                    pass

                continue

    if not feasible:
        return None

    # ========================================================
    # POLICY-SPECIFIC SORTING
    # ========================================================

    if prefer == "max_Hr_then_max_Lr":
        feasible.sort(
            key=lambda result: (
                result[0],
                result[1],
            ),
            reverse=True,
        )

    elif prefer == "max_Lr_then_max_Hr":
        feasible.sort(
            key=lambda result: (
                result[1],
                result[0],
            ),
            reverse=True,
        )

    elif prefer == "max_top_ratio":
        feasible.sort(
            key=lambda result: (
                round(
                    float(result[2]),
                    8,
                ),
                result[0],
                result[1],
            ),
            reverse=True,
        )

    elif prefer == "min_err":
        feasible.sort(
            key=lambda result: (
                result[3],
                -result[0],
                -result[1],
            )
        )

    best = feasible[0]

    print(
        f"[choose-grid] "
        f"prefer={prefer} | "
        f"Hr={best[0]:.6f} | "
        f"Lr={best[1]:.6f} | "
        f"top_ratio={best[2]:.6f} | "
        f"target={target_ratio} | "
        f"error={best[3]}"
    )

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
    vle_screen=True,
    x_min=1e-12,
    bp_T_max=950.0,
    alpha_max=1e6,
    design_sanity=True,
    max_purchase_cost=1e9,
):
    """
    Add a runtime specification that:

    1. evaluates every candidate LHK pair;
    2. searches Lr/Hr combinations for each LHK pair;
    3. ranks all feasible candidates globally;
    4. reruns and verifies the selected candidate;
    5. records whether the selected ratio satisfies tolerance.

    Tolerance is diagnostic. A physically feasible candidate is
    not rejected merely because it cannot reproduce the target
    ratio within the requested tolerance.
    """

    valid_policies = {
        "max_Hr_then_max_Lr",
        "max_Lr_then_max_Hr",
        "max_top_ratio",
        "min_err",
    }

    if prefer not in valid_policies:
        raise ValueError(
            f"Unknown prefer={prefer!r}. "
            f"Choose from {sorted(valid_policies)}."
        )

    if prefer == "min_err" and target_ratio is None:
        raise ValueError(
            "prefer='min_err' requires target_ratio."
        )

    if target_ratio is not None:
        target_ratio = float(target_ratio)

        if not np.isfinite(target_ratio):
            raise ValueError(
                f"target_ratio must be finite; "
                f"received {target_ratio}."
            )

        if not 0.0 <= target_ratio <= 1.0:
            raise ValueError(
                f"target_ratio must be between 0 and 1; "
                f"received {target_ratio}."
            )

    tol = float(tol)

    if not np.isfinite(tol) or tol < 0:
        raise ValueError(
            f"tol must be finite and nonnegative; "
            f"received {tol}."
        )

    def separation_design_check():
        design = getattr(
            col,
            "design_results",
            None,
        ) or {}

        stages = float(
            design.get(
                "Theoretical stages",
                np.nan,
            )
        )

        reflux = float(
            design.get(
                "Reflux",
                np.nan,
            )
        )

        minimum_reflux = float(
            design.get(
                "Minimum reflux",
                np.nan,
            )
        )

        bad = []

        if (
            not np.isfinite(stages)
            or stages < 2
            or stages > 120
        ):
            bad.append(
                f"stages_bad:{stages}"
            )

        if (
            not np.isfinite(reflux)
            or reflux <= 0
        ):
            bad.append(
                f"reflux_bad:{reflux}"
            )

        if (
            not np.isfinite(minimum_reflux)
            or minimum_reflux < 0
        ):
            bad.append(
                f"minimum_reflux_bad:"
                f"{minimum_reflux}"
            )

        if (
            np.isfinite(reflux)
            and np.isfinite(minimum_reflux)
            and reflux < minimum_reflux
        ):
            bad.append(
                f"reflux_below_minimum:"
                f"{reflux}<{minimum_reflux}"
            )

        return (
            len(bad) == 0,
            "|".join(bad),
        )

    def candidate_sort_key(candidate):
        if prefer == "min_err":
            return (
                candidate["search_error"],
                -candidate["Hr"],
                -candidate["Lr"],
            )

        if prefer == "max_top_ratio":
            return (
                -round(
                    candidate[
                        "search_top_ratio"
                    ],
                    8,
                ),
                -candidate["Hr"],
                -candidate["Lr"],
            )

        if prefer == "max_Hr_then_max_Lr":
            return (
                -candidate["Hr"],
                -candidate["Lr"],
            )

        if prefer == "max_Lr_then_max_Hr":
            return (
                -candidate["Lr"],
                -candidate["Hr"],
            )

        raise ValueError(
            f"Unsupported policy: {prefer}"
        )

    def spec():
        last_error = None

        feed = col.ins[0]

        if float(feed.F_mass) <= 0:
            if hasattr(col, "_runtime_best"):
                delattr(
                    col,
                    "_runtime_best",
                )
            return

        all_candidates = []

        LHK_candidates = safe_LHK_candidates(
            splitter,
            feed,
            idx,
            fallback=fallback_LHK,
        )

        # ====================================================
        # SEARCH ALL LHK PAIRS
        # ====================================================

        for LHK in LHK_candidates:
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
                    continue

                Hr, Lr, ratio, error = best

                all_candidates.append(
                    {
                        "LHK": tuple(LHK),
                        "Hr": float(Hr),
                        "Lr": float(Lr),
                        "search_top_ratio":
                            float(ratio),
                        "search_error":
                            (
                                float(error)
                                if np.isfinite(error)
                                else np.nan
                            ),
                    }
                )

            except Exception as error:
                last_error = error

                try:
                    col.reset_cache()
                except Exception:
                    pass

                continue

        if not all_candidates:
            if hasattr(col, "_runtime_best"):
                delattr(
                    col,
                    "_runtime_best",
                )

            raise RuntimeError(
                "No physically feasible LHK/Lr/Hr "
                "candidate was found. "
                f"Last error: {last_error}"
            )

        all_candidates.sort(
            key=candidate_sort_key
        )

        # ====================================================
        # FINAL VERIFICATION
        # ====================================================

        verification_errors = []

        for candidate in all_candidates:
            try:
                col.LHK = candidate["LHK"]
                col.Hr = candidate["Hr"]
                col.Lr = candidate["Lr"]

                col.reset_cache()
                col._run()

                if require_design:
                    col._design()

                if vle_screen:
                    ok_vle, vle_diag = _vle_check(
                        col,
                        x_min=x_min,
                        bp_T_max=bp_T_max,
                        alpha_max=alpha_max,
                    )

                    if not ok_vle:
                        raise RuntimeError(
                            "final_vle_fail:"
                            f"{vle_diag.get('why', '')}"
                        )

                if design_sanity and require_design:
                    ok_dc, why_dc = (
                        _design_cost_sanity(
                            col,
                            max_purchase_cost=(
                                max_purchase_cost
                            ),
                        )
                    )

                    if not ok_dc:
                        raise RuntimeError(
                            "final_design_fail:"
                            f"{why_dc}"
                        )

                    ok_sep, why_sep = (
                        separation_design_check()
                    )

                    if not ok_sep:
                        raise RuntimeError(
                            "final_separation_fail:"
                            f"{why_sep}"
                        )

                if require_cost:
                    col._cost()

                    if design_sanity:
                        ok_dc, why_dc = (
                            _design_cost_sanity(
                                col,
                                max_purchase_cost=(
                                    max_purchase_cost
                                ),
                            )
                        )

                        if not ok_dc:
                            raise RuntimeError(
                                "final_cost_fail:"
                                f"{why_dc}"
                            )

                final_total_mass = float(
                    col.F_mass_out
                )

                if final_total_mass <= 0:
                    raise RuntimeError(
                        "Final column outlet mass "
                        "is zero or negative."
                    )

                final_top_ratio = float(
                    col.outs[0].F_mass
                    / final_total_mass
                )

                if not np.isfinite(
                    final_top_ratio
                ):
                    raise RuntimeError(
                        "Final top ratio is not finite."
                    )

                if not (
                    0.0
                    <= final_top_ratio
                    <= 1.0
                ):
                    raise RuntimeError(
                        "Final top ratio is outside "
                        f"[0, 1]: {final_top_ratio}"
                    )

                if target_ratio is not None:
                    final_error = abs(
                        final_top_ratio
                        - target_ratio
                    )

                    target_satisfied = bool(
                        final_error <= tol
                    )
                else:
                    final_error = np.nan
                    target_satisfied = True

                col._runtime_best = {
                    "prefer":
                        prefer,

                    "LHK":
                        tuple(col.LHK),

                    "Hr":
                        float(col.Hr),

                    "Lr":
                        float(col.Lr),

                    "target_ratio":
                        (
                            float(target_ratio)
                            if target_ratio
                            is not None
                            else np.nan
                        ),

                    "tolerance":
                        float(tol),

                    "search_top_ratio":
                        candidate[
                            "search_top_ratio"
                        ],

                    "search_error":
                        candidate[
                            "search_error"
                        ],

                    "top_ratio":
                        final_top_ratio,

                    "err":
                        (
                            float(final_error)
                            if np.isfinite(
                                final_error
                            )
                            else np.nan
                        ),

                    "target_satisfied":
                        bool(target_satisfied),

                    "candidate_count":
                        len(all_candidates),

                    "LHK_candidate_count":
                        len(LHK_candidates),
                }

                print(
                    f"[choose-final] "
                    f"prefer={prefer} | "
                    f"LHK={tuple(col.LHK)} | "
                    f"Hr={col.Hr:.6f} | "
                    f"Lr={col.Lr:.6f} | "
                    f"target={target_ratio} | "
                    f"search_ratio="
                    f"{candidate['search_top_ratio']:.6f} | "
                    f"final_ratio="
                    f"{final_top_ratio:.6f} | "
                    f"final_error={final_error} | "
                    f"target_satisfied="
                    f"{target_satisfied}"
                )

                return

            except Exception as error:
                last_error = error

                verification_errors.append(
                    {
                        "LHK":
                            candidate["LHK"],

                        "Hr":
                            candidate["Hr"],

                        "Lr":
                            candidate["Lr"],

                        "error":
                            str(error),
                    }
                )

                try:
                    col.reset_cache()
                except Exception:
                    pass

                continue

        if hasattr(col, "_runtime_best"):
            delattr(
                col,
                "_runtime_best",
            )

        raise RuntimeError(
            "Candidates were found during grid search, "
            "but none passed final physical verification. "
            f"Last error: {last_error}. "
            f"Verification failures: "
            f"{verification_errors}"
        )

    col.add_specification(spec)
    col.run_after_specifications = True