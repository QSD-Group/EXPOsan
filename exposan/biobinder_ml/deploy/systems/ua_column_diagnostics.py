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
    """Try Aitken/secant first, then scan brackets and solve with Brent.

    Every accepted rate must be real, finite, greater than -0.999999 and
    have absolute NPV <= npv_tol_abs (USD; default $1,000). The fallback
    searches -0.999999 < r <= 1e6 on a log(1+r) grid. Multiple distinct
    fallback roots are reported as unavailable rather than selecting one.
    A finite scan cannot prove uniqueness or detect every tangent root;
    an accepted primary-solver root is also not a uniqueness guarantee.
    """
    from scipy.optimize import brentq

    CF = tea.cashflow_array
    dur = tea._get_duration_array()
    if CF is None or dur is None:
        return np.nan, "CF_or_dur_nonfinite"
    CF = np.asarray(CF, dtype=float)
    dur = np.asarray(dur, dtype=float)
    if (CF.ndim != 1 or dur.shape != CF.shape or not CF.size
            or not np.isfinite(CF).all() or not np.isfinite(dur).all()):
        return np.nan, "CF_or_dur_nonfinite"
    if not np.any(CF):
        return np.nan, "IRR_indeterminate_zero_cashflows"
    if not (np.any(CF < 0) and np.any(CF > 0)):
        return np.nan, "no_cashflow_sign_change"
    if (not np.isfinite(npv_tol_abs) or npv_tol_abs < 0
            or not np.isfinite(xtol) or xtol <= 0):
        raise ValueError('NPV tolerance must be nonnegative; xtol must be positive')

    args = (CF, dur)
    IRR_guess = float(getattr(tea, "IRR", 0.10) or 0.10)
    if not np.isfinite(IRR_guess) or IRR_guess <= -0.999999:
        IRR_guess = 0.10

    def residual(r):
        with np.errstate(over='ignore', under='ignore', invalid='ignore', divide='ignore'):
            return float(_NPV_at_r(r, *args))

    def verified(r):
        if not np.isfinite(r) or r <= -0.999999:
            return False
        value = residual(r)
        return np.isfinite(value) and abs(value) <= npv_tol_abs

    try:
        x0, x1 = IRR_guess, 1.0001 * IRR_guess + 1e-3
        if np.isfinite(residual(x0)) and np.isfinite(residual(x1)):
            with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
                irr = flx.aitken_secant(
                    _NPV_at_r, x0, x1, xtol=xtol, ytol=ytol,
                    maxiter=maxiter, args=args, checkiter=False)
            if verified(irr):
                return float(irr), "ok"
    except Exception:
        # A failed primary iteration does not establish absence of an IRR.
        pass

    rates = np.expm1(np.linspace(np.log(1.000001e-6), np.log1p(1e6), 2049))
    rates = np.unique(np.append(rates, [0., IRR_guess]))
    values = [residual(r) for r in rates]
    roots = [float(r) for r, value in zip(rates, values)
             if value == 0. and verified(r)]
    rejected_residual = False
    failed_bracket = False
    for low, high, f_low, f_high in zip(rates[:-1], rates[1:], values[:-1], values[1:]):
        if not (np.isfinite(f_low) and np.isfinite(f_high)):
            continue
        if f_low == 0. or f_high == 0. or np.signbit(f_low) == np.signbit(f_high):
            continue
        try:
            # Tight rate tolerance is needed for steep NPV curves at negative
            # rates. Acceptance is determined by the USD residual, not xtol.
            root = brentq(residual, low, high, xtol=min(xtol, 1e-12),
                          rtol=4*np.finfo(float).eps, maxiter=maxiter)
            if not verified(root):
                root = brentq(residual, low, high, xtol=np.nextafter(0., 1.),
                              rtol=4*np.finfo(float).eps, maxiter=maxiter)
            if verified(root):
                roots.append(float(root))
            else:
                rejected_residual = True
        except (ValueError, RuntimeError, OverflowError):
            failed_bracket = True

    distinct = []
    for root in sorted(roots):
        if not distinct or not np.isclose(root, distinct[-1], rtol=1e-9, atol=1e-10):
            distinct.append(root)
    if len(distinct) > 1:
        return np.nan, f"multiple_bracketed_IRRs:{len(distinct)}"
    if rejected_residual or failed_bracket:
        return np.nan, "bracketed_search_incomplete_or_residual_not_close"
    if distinct:
        return distinct[0], "ok_bracketed"
    return np.nan, "no_verified_IRR_in_search_range"
