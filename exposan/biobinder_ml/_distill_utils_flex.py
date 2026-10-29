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
from datetime import datetime

_TRIALS = {}

def get_distillation_trials(col):
    return list(_TRIALS.get(col, []))


__all__ = (
    "safe_LHK_candidates",
    "find_best_Lr_Hr_runtime",
    "add_runtime_LHK_LrHr_spec",
)

def _vle_check(col, x_min=1e-12, bp_T_max=950.0, alpha_max=1e6, K_max=1e14):
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
            tot = float(np.sum(list(C.values())))
        except Exception:
            tot = np.nan
        if (not np.isfinite(tot)) or (tot <= 0) or (tot > max_purchase_cost):
            bad.append(f"purchase_cost_bad:{tot}")
    else:
        if not allow_missing_purchase_costs:
            bad.append("missing_purchase_costs")

    return (len(bad) == 0), "|".join(bad)

def safe_LHK_candidates(
    splitter, col_feed, idx, fallback=("4M-PHYNO", "INDOLE"),
    min_massfrac=1e-6, target_ratio=None, max_adjacent_pairs=6,
):
    """Keep legacy pairs and add adjacent boiling-point pairs near the mass cut.

    Uses normalized inlet mass fractions, so uniform flow scaling does not
    change pair selection. Normal boiling points generate candidates only;
    column VLE/design/cost checks still decide feasibility. The same rule
    applies to every feedstock. No thermodynamic checks are disabled.
    """
    total = float(col_feed.F_mass)
    if total <= 0:
        return []

    def massfrac(cid):
        try:
            return float(col_feed.imass[cid]) / total
        except (KeyError, IndexError, TypeError, ValueError):
            return 0.0

    pairs = []
    try:
        primary = tuple(splitter.keys[idx])
    except (AttributeError, IndexError, TypeError):
        primary = ()
    for pair in (primary, tuple(fallback) if fallback else ()):
        if len(pair) == 2 and all(massfrac(cid) > min_massfrac for cid in pair):
            pairs.append(pair)

    # Include every present component in ordering so that a small component
    # is not silently skipped between supposedly adjacent keys.
    present = []
    for chemical in col_feed.chemicals:
        tb = getattr(chemical, 'Tb', None)
        fraction = massfrac(chemical.ID)
        if (fraction > 0 and tb is not None and np.isfinite(tb)
                and getattr(chemical, 'locked_state', None) not in ('g', 'l', 's')):
            present.append((float(tb), chemical.ID, fraction))
    present.sort(key=lambda item: (item[0], item[1]))

    # Without a target, anchor around the original pair's boiling boundary.
    primary_tbs = [item[0] for item in present if item[1] in primary]
    anchor_tb = np.mean(primary_tbs) if primary_tbs else None
    adjacent = []
    cumulative = 0.0
    for light, heavy in zip(present, present[1:]):
        cumulative += light[2]
        if (light[2] <= min_massfrac or heavy[2] <= min_massfrac
                or heavy[0] <= light[0]):
            continue
        if target_ratio is not None:
            distance = abs(cumulative - float(target_ratio))
        elif anchor_tb is not None:
            distance = abs((light[0] + heavy[0]) / 2 - anchor_tb)
        else:
            continue
        adjacent.append((distance, light[1], heavy[1]))
    adjacent.sort()
    pairs.extend((lk, hk) for _, lk, hk in adjacent[:max_adjacent_pairs])
    return list(dict.fromkeys(pairs))

# Cost ranking uses the column's own baseline purchase-cost dictionary.
# Auxiliary equipment coverage is whatever this BioSTEAM version places there.
_POLICIES = {'max_Hr_then_max_Lr', 'max_Lr_then_max_Hr',
             'max_top_ratio', 'min_err', 'min_purchase_cost'}


def _fresh(col):
    col.reset_cache()
    # Clear previous candidate results on the column and auxiliary units.
    units = [col] + list(getattr(col, 'auxiliary_units', ()) or ())
    for unit in units:
        for name in ('design_results', 'baseline_purchase_costs',
                     'purchase_costs', 'installed_costs'):
            obj = getattr(unit, name, None)
            if hasattr(obj, 'clear'):
                obj.clear()


def _evaluate(col, LHK, Hr, Lr, options, phase):
    trial = dict(Column=col.ID, Timestamp=datetime.now().astimezone().isoformat(),
                 LHK_light=LHK[0], LHK_heavy=LHK[1], Hr=float(Hr), Lr=float(Lr),
                 Pressure_Pa=float(col.P), Target_ratio=options['target_ratio'],
                 Tolerance=options['tol'], Status='REJECT', Stage=phase,
                 Reason='', Top_ratio=np.nan, Purchase_cost=np.nan)
    try:
        col.LHK = tuple(LHK)
        col.Hr, col.Lr = float(Hr), float(Lr)
        _fresh(col)
        trial['Stage'] = phase + ':run'
        col._run()
        # VLE screening can update product temperatures: do it BEFORE design.
        if options['vle_screen']:
            trial['Stage'] = phase + ':VLE'
            ok, diag = _vle_check(col, x_min=options['x_min'],
                bp_T_max=options['bp_T_max'], alpha_max=options['alpha_max'],
                K_max=options['K_max'])
            trial.update({'VLE_' + k: v for k, v in diag.items()})
            if not ok:
                raise RuntimeError('VLE:' + str(diag.get('why', diag)))
        trial['Stage'] = phase + ':design'
        if options['require_design']:
            col._design()
            if options['design_sanity']:
                ok, reason = _design_cost_sanity(col)
                if not ok:
                    raise RuntimeError(reason)
        if options['require_cost']:
            trial['Stage'] = phase + ':cost'
            col._cost()
            costs = dict(col.baseline_purchase_costs)
            values = np.asarray(list(costs.values()), dtype=float)
            if not len(values) or not np.all(np.isfinite(values)) or np.any(values < 0):
                raise RuntimeError('missing/nonfinite/negative purchase-cost entry')
            cost = float(values.sum())
            trial['Purchase_cost'] = cost
            trial['Purchase_cost_items'] = repr(costs)
            if cost <= 0 or cost > options['max_purchase_cost']:
                raise RuntimeError(f'purchase_cost_bad:{cost}')
        total = float(col.F_mass_out)
        ratio = float(col.outs[0].F_mass) / total if total > 0 else np.nan
        trial['Top_ratio'] = ratio
        if not np.isfinite(ratio):
            raise RuntimeError('nonfinite_top_ratio')
        target = options['target_ratio']
        error = abs(ratio - float(target)) if target is not None else np.nan
        trial['Target_error'] = error
        if target is not None and error > options['tol']:
            raise RuntimeError('target_outside_tolerance')
        design = col.design_results
        trial.update(Stages=design.get('Theoretical stages', np.nan),
                     Reflux=design.get('Reflux', np.nan),
                     Status='FEASIBLE', Reason='')
        return dict(LHK=tuple(LHK), Hr=float(Hr), Lr=float(Lr),
                    top_ratio=ratio, err=error, purchase_cost=trial['Purchase_cost'])
    except Exception as error:
        trial['Reason'] = f'{type(error).__name__}: {error}'
        return None
    finally:
        _TRIALS.setdefault(col, []).append(trial)


def _rank(candidates, prefer, cost_band, cost_tiebreak):
    if prefer == 'min_purchase_cost':
        # Repeatedly form a band around the cheapest remaining candidate.
        pending = list(candidates)
        ranked = []
        while pending:
            minimum = min(c['purchase_cost'] for c in pending)
            band = [c for c in pending if c['purchase_cost'] <= minimum*(1+cost_band)]
            if cost_tiebreak == 'min_err':
                key = lambda c: (c['err'] if np.isfinite(c['err']) else np.inf,
                                  c['purchase_cost'], -c['top_ratio'])
            else:
                key = lambda c: (-c['top_ratio'], c['purchase_cost'])
            ranked.extend(sorted(band, key=key))
            ids = {id(c) for c in band}
            pending = [c for c in pending if id(c) not in ids]
        return ranked
    if prefer == 'max_top_ratio':
        key = lambda c: (-round(c['top_ratio'], 6), -c['Hr'], -c['Lr'])
    elif prefer == 'max_Hr_then_max_Lr':
        key = lambda c: (-c['Hr'], -c['Lr'])
    elif prefer == 'max_Lr_then_max_Hr':
        key = lambda c: (-c['Lr'], -c['Hr'])
    else:
        key = lambda c: (c['err'], -c['Hr'], -c['Lr'])
    return sorted(candidates, key=key)


def _options(target_ratio, tol, prefer, require_design, require_cost,
             vle_screen, x_min, bp_T_max, alpha_max, design_sanity,
             max_purchase_cost, K_max, cost_band, cost_tiebreak):
    if prefer not in _POLICIES:
        raise ValueError(f'Unknown prefer={prefer!r}; choose {sorted(_POLICIES)}')
    if prefer == 'min_err' and target_ratio is None:
        raise ValueError('min_err requires target_ratio')
    if cost_band < 0 or cost_tiebreak not in ('min_err', 'max_top_ratio'):
        raise ValueError('Invalid cost_band or cost_tiebreak')
    if tol < 0 or (target_ratio is not None and not 0 <= target_ratio <= 1):
        raise ValueError('Invalid target_ratio or tol')
    if prefer == 'min_purchase_cost':
        require_design = require_cost = True
    if require_cost:
        require_design = True
    return dict(target_ratio=target_ratio, tol=float(tol),
        require_design=require_design, require_cost=require_cost,
        vle_screen=vle_screen, x_min=x_min, bp_T_max=bp_T_max,
        alpha_max=alpha_max, design_sanity=design_sanity,
        max_purchase_cost=max_purchase_cost, K_max=K_max)


def _scan(col, LHK, options, Lr_grid, Hr_grid):
    if Lr_grid is None:
        Lr_grid = [0.995, .99, .985, .98, .975, .97, .965, .96, .95, .93, .90]
    if Hr_grid is None:
        Hr_grid = [.99, .98, .97, .95, .93, .90, .89, .87, .85]
    candidates = []
    for Hr in Hr_grid:
        for Lr in Lr_grid:
            if not (0 < Hr < 1 and 0 < Lr < 1):
                raise ValueError('Recoveries must be strictly between zero and one')
            candidate = _evaluate(col, LHK, Hr, Lr, options, 'grid')
            if candidate is not None:
                candidates.append(candidate)
    return candidates


def _verify_ranked(col, candidates, options, prefer, cost_band, cost_tiebreak,
                   cost_verify_rtol):
    for candidate in _rank(candidates, prefer, cost_band, cost_tiebreak):
        verified = _evaluate(col, candidate['LHK'], candidate['Hr'],
                             candidate['Lr'], options, 'final')
        if verified is None:
            continue
        if options['require_cost'] and not np.isclose(
                verified['purchase_cost'], candidate['purchase_cost'],
                rtol=cost_verify_rtol, atol=1.0):
            _TRIALS[col][-1].update(Status='REJECT', Reason='final_cost_changed',
                                   Grid_purchase_cost=candidate['purchase_cost'])
            continue
        _TRIALS[col][-1]['Status'] = 'FINAL_ACCEPTED'
        col._runtime_best = dict(verified, prefer=prefer)
        print(f"[choose] {col.ID}: prefer={prefer}, LHK={verified['LHK']}, "
              f"Hr={verified['Hr']}, Lr={verified['Lr']}, "
              f"top_ratio={verified['top_ratio']:.6f}, "
              f"purchase_cost={verified['purchase_cost']:,.2f}")
        return verified
    return None


def find_best_Lr_Hr_runtime(col, target_ratio=None, tol=.20,
        Lr_grid=None, Hr_grid=None, prefer='max_Lr_then_max_Hr',
        require_design=True, require_cost=False, vle_screen=True,
        x_min=1e-12, bp_T_max=950., alpha_max=1e6, design_sanity=True,
        max_purchase_cost=1e9, K_max=1e14, cost_band=0.0,
        cost_tiebreak='min_err', cost_verify_rtol=.01):
    """Search current LHK; return legacy (Hr, Lr, ratio, error) tuple."""
    options = _options(target_ratio, tol, prefer, require_design, require_cost,
        vle_screen, x_min, bp_T_max, alpha_max, design_sanity,
        max_purchase_cost, K_max, cost_band, cost_tiebreak)
    candidates = _scan(col, tuple(col.LHK), options, Lr_grid, Hr_grid)
    best = _verify_ranked(col, candidates, options, prefer, cost_band,
                          cost_tiebreak, cost_verify_rtol)
    return None if best is None else (best['Hr'], best['Lr'], best['top_ratio'], best['err'])


def add_runtime_LHK_LrHr_spec(splitter, col, idx,
        fallback_LHK=('4M-PHYNO', 'INDOLE'), target_ratio=None, tol=.20,
        Lr_grid=None, Hr_grid=None, prefer='max_top_ratio',
        require_design=True, require_cost=False, vle_screen=True,
        x_min=1e-12, bp_T_max=950., alpha_max=1e6, design_sanity=True,
        max_purchase_cost=1e9, K_max=1e14, cost_band=0.0,
        cost_tiebreak='min_err', cost_verify_rtol=.01):
    """Cost policy compares ALL legacy/adjacent pairs; old policies keep pair priority.

    cost_band=0 is strict minimum purchase cost. Set .01 to permit a 1% band.
    All ratios are mass-based top/feed fractions, not purity specifications.
    No TEA, prices, revenue or system simulation is evaluated here.
    """
    options = _options(target_ratio, tol, prefer, require_design, require_cost,
        vle_screen, x_min, bp_T_max, alpha_max, design_sanity,
        max_purchase_cost, K_max, cost_band, cost_tiebreak)
    def spec():
        _TRIALS[col] = []
        if hasattr(col, '_runtime_best'):
            delattr(col, '_runtime_best')
        feed = col.ins[0]
        if float(feed.F_mass) <= 0:
            for stream in col.outs:
                stream.empty()
            return
        all_candidates = []
        for LHK in safe_LHK_candidates(
                splitter, feed, idx, fallback=fallback_LHK,
                target_ratio=target_ratio):
            candidates = _scan(col, LHK, options, Lr_grid, Hr_grid)
            if prefer == 'min_purchase_cost':
                all_candidates.extend(candidates)
            else:
                best = _verify_ranked(col, candidates, options, prefer,
                    cost_band, cost_tiebreak, cost_verify_rtol)
                if best is not None:
                    return
        if prefer == 'min_purchase_cost':
            best = _verify_ranked(col, all_candidates, options, prefer,
                cost_band, cost_tiebreak, cost_verify_rtol)
            if best is not None:
                return
        reasons = [t['Reason'] for t in _TRIALS[col] if t.get('Reason')]
        raise RuntimeError('No verified feasible LHK/Lr/Hr found. '
                           'Recent rejections: ' + ' | '.join(reasons[-5:]))
    col.add_specification(spec)
    col.run_after_specifications = False

__all__ = ('safe_LHK_candidates', 'find_best_Lr_Hr_runtime',
           'add_runtime_LHK_LrHr_spec', 'get_distillation_trials')
