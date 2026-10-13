# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import numpy as np

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return np.nan

def hx_duty_kJ_per_hr(hx):
    try:
        return float(hx.outs[0].H - hx.ins[0].H)  # kJ/hr
    except Exception:
        return np.nan

def safe_get_design(unit, key):
    try:
        v = unit.design_results.get(key, np.nan)
        return float(v) if np.isfinite(v) else np.nan
    except Exception:
        try:
            return float(unit.design_results[key])
        except Exception:
            return np.nan

def add_crudeheavy_dis_diagnostics(out, col, prefix="CHD_"):
    import numpy as _np

    feed = col.ins[0]
    Dstream, Bstream = col.outs

    # Basic feed state
    out[prefix+"feed_F_mass_kgph"] = safe_float(feed.F_mass)
    out[prefix+"feed_F_mol_kmolph"] = safe_float(getattr(feed, "F_mol", np.nan))
    out[prefix+"feed_T_K"] = safe_float(feed.T)
    out[prefix+"feed_P_Pa"] = safe_float(getattr(feed, "P", np.nan))
    out[prefix+"feed_H_kJph"] = safe_float(getattr(feed, "H", np.nan))

    # Nonzero components in feed
    try:
        nz = [ID for ID, m in zip(feed.chemicals.IDs, feed.mol) if m > 0]
        out[prefix+"feed_nonzero_n"] = len(nz)
        out[prefix+"feed_nonzero_IDs_head"] = ",".join(nz[:30])  # keep cell manageable
    except Exception:
        out[prefix+"feed_nonzero_n"] = np.nan
        out[prefix+"feed_nonzero_IDs_head"] = ""

    # Column settings / results
    out[prefix+"P_Pa"] = safe_float(getattr(col, "P", np.nan))
    out[prefix+"P_psi"] = safe_float(getattr(col, "P", np.nan) * 0.000145078) if np.isfinite(getattr(col, "P", np.nan)) else np.nan
    out[prefix+"is_vacuum"] = bool(out[prefix+"P_psi"] < 14.68) if np.isfinite(out[prefix+"P_psi"]) else np.nan

    # policy / runtime best
    best = getattr(col, "_runtime_best", None)
    out[prefix+"policy_found"] = bool(best is not None)
    if best:
        out[prefix+"policy_prefer"] = str(best.get("prefer", ""))
        out[prefix+"policy_LHK"] = str(best.get("LHK", ""))
        out[prefix+"policy_Lr"] = safe_float(best.get("Lr", np.nan))
        out[prefix+"policy_Hr"] = safe_float(best.get("Hr", np.nan))
        out[prefix+"policy_top_ratio"] = safe_float(best.get("top_ratio", np.nan))
        out[prefix+"policy_err"] = safe_float(best.get("err", np.nan))

    # final LHK/Lr/Hr
    try:
        out[prefix+"LHK_light"] = col.LHK[0]
        out[prefix+"LHK_heavy"] = col.LHK[1]
    except Exception:
        out[prefix+"LHK_light"] = str(getattr(col, "LHK", ""))
        out[prefix+"LHK_heavy"] = np.nan
    out[prefix+"Lr_final"] = safe_float(getattr(col, "Lr", np.nan))
    out[prefix+"Hr_final"] = safe_float(getattr(col, "Hr", np.nan))

    # top ratio actually achieved
    try:
        denom = float(getattr(col, "F_mass_out", np.nan))
        out[prefix+"column_top_ratio"] = float(col.outs[0].F_mass / denom) if np.isfinite(denom) and denom > 0 else np.nan
    except Exception:
        out[prefix+"column_top_ratio"] = np.nan

    # Design results (grab a lot)
    D = getattr(col, "design_results", {}) or {}
    out[prefix+"N_theoretical_raw"] = D.get("Theoretical stages", None)
    out[prefix+"N_feed_raw"] = D.get("Theoretical feed stage", None)
    out[prefix+"Reflux"] = safe_get_design(col, "Reflux")
    out[prefix+"Rmin"] = safe_get_design(col, "Minimum reflux")

    out[prefix+"Diameter_ft"] = safe_float(D.get("Rectifier diameter", D.get("Diameter", np.nan)))
    out[prefix+"Height_ft"] = safe_float(D.get("Rectifier height", D.get("Height", np.nan)))
    out[prefix+"twall_in"] = safe_float(D.get("Rectifier wall thickness", D.get("Wall thickness", np.nan)))
    out[prefix+"weight_lb"] = safe_float(D.get("Rectifier weight", D.get("Weight", np.nan)))

    # condenser/reboiler duties + costs
    try:
        out[prefix+"Q_cond_kJhr"] = hx_duty_kJ_per_hr(col.condenser)
        out[prefix+"cond_utility_cost_USDhr"] = safe_float(col.condenser.utility_cost)
    except Exception:
        out[prefix+"Q_cond_kJhr"] = np.nan
        out[prefix+"cond_utility_cost_USDhr"] = np.nan

    try:
        out[prefix+"Q_reb_kJhr"] = hx_duty_kJ_per_hr(col.reboiler)
        out[prefix+"reb_utility_cost_USDhr"] = safe_float(col.reboiler.utility_cost)
    except Exception:
        out[prefix+"Q_reb_kJhr"] = np.nan
        out[prefix+"reb_utility_cost_USDhr"] = np.nan

    # VLE internals (IDs, LK/HK indices)
    IDs = getattr(col, "_IDs_vle", None)
    LHK_vle = getattr(col, "_LHK_vle_index", None)
    out[prefix+"partial_condenser"] = bool(getattr(col, "_partial_condenser", True))
    out[prefix+"IDs_vle_n"] = len(IDs) if IDs else 0

    if IDs and (LHK_vle is not None):
        LK_vle, HK_vle = int(LHK_vle[0]), int(LHK_vle[1])

        # distillate/bottom compositions for IDs
        try:
            zD = Dstream.get_normalized_mol(IDs)
            zB = Bstream.get_normalized_mol(IDs)
            out[prefix+"zD_LK"] = safe_float(zD[LK_vle])
            out[prefix+"zD_HK"] = safe_float(zD[HK_vle])
            out[prefix+"zB_LK"] = safe_float(zB[LK_vle])
            out[prefix+"zB_HK"] = safe_float(zB[HK_vle])
            out[prefix+"zD_min"] = safe_float(np.min(zD))
            out[prefix+"zB_min"] = safe_float(np.min(zB))
        except Exception as e:
            out[prefix+"z_comp_error"] = str(e)

        dew = getattr(col, "_dew_point", None)
        bub = getattr(col, "_bubble_point", None)

        try:
            if getattr(col, "_partial_condenser", True):
                dp = dew(zD, P=col.P)
            else:
                dp = bub(zD, P=col.P)
            bp = bub(zB, P=col.P)

            out[prefix+"dp_T_K"] = safe_float(dp.T)
            out[prefix+"bp_T_K"] = safe_float(bp.T)

            out[prefix+"dp_min_x_raw"] = safe_float(_np.min(dp.x))
            out[prefix+"bp_min_x_raw"] = safe_float(_np.min(bp.x))

            out[prefix+"dp_x_LK"] = safe_float(dp.x[LK_vle])
            out[prefix+"dp_y_LK"] = safe_float(dp.y[LK_vle])
            out[prefix+"dp_x_HK"] = safe_float(dp.x[HK_vle])
            out[prefix+"dp_y_HK"] = safe_float(dp.y[HK_vle])

            out[prefix+"bp_x_LK"] = safe_float(bp.x[LK_vle])
            out[prefix+"bp_y_LK"] = safe_float(bp.y[LK_vle])
            out[prefix+"bp_x_HK"] = safe_float(bp.x[HK_vle])
            out[prefix+"bp_y_HK"] = safe_float(bp.y[HK_vle])

            eps = 1e-16
            Kdist = dp.y / _np.maximum(dp.x, eps)
            Kbot  = bp.y / _np.maximum(bp.x, eps)

            out[prefix+"Kdist_LK"] = safe_float(Kdist[LK_vle])
            out[prefix+"Kdist_HK"] = safe_float(Kdist[HK_vle])
            out[prefix+"Kbot_LK"]  = safe_float(Kbot[LK_vle])
            out[prefix+"Kbot_HK"]  = safe_float(Kbot[HK_vle])

            out[prefix+"alpha_dist_LK_over_HK"] = safe_float(Kdist[LK_vle] / max(Kdist[HK_vle], eps))
            out[prefix+"alpha_bot_LK_over_HK"]  = safe_float(Kbot[LK_vle] / max(Kbot[HK_vle], eps))

            alpha_mean = col._estimate_mean_volatilities_relative_to_heavy_key()
            out[prefix+"alpha_LK_used"] = safe_float(alpha_mean[LK_vle])
            out[prefix+"alpha_LK_le_1"] = bool(alpha_mean[LK_vle] <= 1.0)

            # Fenske ratio using actual LK/HK mole flows (shortcut internal indices)
            try:
                LK_i, HK_i = col._LHK_index
                LK_D = float(Dstream.mol[LK_i])
                HK_D = float(Dstream.mol[HK_i])
                LK_B = float(Bstream.mol[LK_i])
                HK_B = float(Bstream.mol[HK_i])
                ratio = (LK_D / max(HK_D, 1e-30)) * (HK_B / max(LK_B, 1e-30))
                out[prefix+"Fenske_ratio"] = safe_float(ratio)
                out[prefix+"log10_Fenske_ratio"] = safe_float(_np.log10(max(ratio, 1e-30)))
                out[prefix+"log10_alpha_LK"] = safe_float(_np.log10(max(alpha_mean[LK_vle], 1e-30)))
            except Exception as e:
                out[prefix+"Fenske_error"] = str(e)

            # Underwood diagnostics
            try:
                q = col.get_feed_quality()
                theta = col._solve_Underwood_constant(alpha_mean, alpha_mean[LK_vle])
                z_d = Dstream.get_normalized_mol(IDs)
                Rm = float((alpha_mean * z_d / (alpha_mean - theta)).sum() - 1.0)
                out[prefix+"q"] = safe_float(q)
                out[prefix+"theta"] = safe_float(theta)
                out[prefix+"Rm_calc"] = safe_float(Rm)
                out[prefix+"Rm_bad"] = bool((not _np.isfinite(Rm)) or (Rm < 0))
            except Exception as e:
                out[prefix+"Underwood_error"] = str(e)

        except Exception as e:
            out[prefix+"dp_bp_error"] = str(e)
    else:
        out[prefix+"dp_bp_error"] = "missing _IDs_vle or _LHK_vle_index"

    # Mass flows of top/bottom
    try:
        out[prefix+"F_top_kgph"] = safe_float(Dstream.F_mass)
        out[prefix+"F_bot_kgph"] = safe_float(Bstream.F_mass)
    except Exception:
        out[prefix+"F_top_kgph"] = np.nan
        out[prefix+"F_bot_kgph"] = np.nan

import flexsolve as flx

def _NPV_at_r(r, CF, dur):
    r = float(r)
    if (not np.isfinite(r)) or (r <= -0.999999):
        return np.nan
    return (CF / (1.0 + r)**dur).sum()

def solve_IRR_aitken_limit(tea, maxiter=200, xtol=1e-6, ytol=10.0, npv_tol_abs=1e3):
    CF = tea.cashflow_array
    dur = tea._get_duration_array()

    if (CF is None) or (dur is None) or (not np.isfinite(CF).all()) or (not np.isfinite(dur).all()):
        return np.nan, "CF_or_dur_nonfinite"

    args = (CF, dur)
    IRR_guess = float(getattr(tea, "IRR", 0.10) or 0.10)
    if (not np.isfinite(IRR_guess)) or (IRR_guess <= -0.999999):
        IRR_guess = 0.10

    x0 = IRR_guess
    x1 = 1.0001 * IRR_guess + 1e-3

    y0 = _NPV_at_r(x0, *args)
    y1 = _NPV_at_r(x1, *args)
    if (not np.isfinite(y0)) or (not np.isfinite(y1)):
        return np.nan, "endpoint_nonfinite"

    try:
        irr = flx.aitken_secant(
            _NPV_at_r,
            x0, x1,
            xtol=xtol, ytol=ytol, maxiter=maxiter,
            args=args, checkiter=False
        )
    except Exception as e:
        return np.nan, f"solver_exception:{type(e).__name__}"

    if (not np.isfinite(irr)) or (irr <= -0.999999):
        return np.nan, "solver_returned_invalid"

    y = _NPV_at_r(irr, *args)
    if not np.isfinite(y):
        return np.nan, "npv_at_solution_nonfinite"
    if abs(float(y)) > float(npv_tol_abs):
        return np.nan, f"npv_not_close:{float(y):.3g}"

    return float(irr), "ok"
