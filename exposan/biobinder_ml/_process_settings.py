#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Ali Ahmad <aa3056@scarletmail.rutgers.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import math
import qsdsan as qs
import numpy as np
from exposan.saf.utils import find_Lr_Hr
from exposan.htl import HTL_TEA
from exposan.saf import _process_settings
from exposan.saf._process_settings import *


__all__ = [i for i in _process_settings.__all__ if i != 'dry_flowrate']
__all__.extend([
    'BiobinderTEA',
    'central_dry_flowrate',
    'pilot_dry_flowrate',
    'cutoff_fracs_from_feed',
    'find_Lr_Hr',
    'TRACI_CATEGORIES',
    'feedprice_dct'
    ])
_ton_to_kg = 907.185
moisture = 0.7566 #0.403  # 0.7566
ash = (1-moisture)*0.0571 #0.05 #0.0571
feedstock_composition = {
    'Water': moisture,
    'Lipids': (1-moisture)*0.6245, #0.07, #0.6245
    'Proteins': (1-moisture)*0.0238, #0.07, #0.0238
    'Carbohydrates': (1-moisture)*0.2946, #0.81,#0.2946
    'Ash': ash,
    }
tea_indices = qs.utils.indices.tea_indices
central_dry_flowrate = dry_flowrate # 110 tpd converted to kg/hr
pilot_dry_flowrate = 11.46 # kg/hr

cost_year = 2020
PCE_indices = tea_indices['PCEPI']
FOG_year=2021
FOG_factor= PCE_indices[cost_year]/PCE_indices[FOG_year]
SnowdenSwan_year = 2016
SnowdenSwan_factor = PCE_indices[cost_year]/PCE_indices[SnowdenSwan_year]
Green_year= 2000
Green_factor= PCE_indices[cost_year]/PCE_indices[Green_year]
Badgett_year= 2019
Badgett_factor= PCE_indices[cost_year]/PCE_indices[Badgett_year]

feedprice_dct = {
    'fog': +390/1e3*FOG_factor,      # $/kg, $600/dry tonne*0.65 yellow grease USDA, FOG moisture 0.35, Ou et al. 2021
    'green': -16/_ton_to_kg *Green_factor, #Komilis et al 2000, 2000$ assumed, actual 1998 (PCEPI not available for 1998)
    'food': -39.7/1e3*SnowdenSwan_factor, # Ahmad et al. 2025
    'sludge': -40/1e3*Badgett_factor, #median, Badgett et al. 2019
    'manure': 0.0, #median, Badgett et al. 2019
}

# Salad dressing waste
HTL_yields = {
    'gas': 0.1756,   #8.50, 17.56
    'aqueous': 0.2925,  #36.70, 29.25
    'biocrude': 0.5219,  #24.27, 52.19
    'char': 1-0.1756-0.2925-0.5219, 
    }

_COEFS = {
    "const": 0.0013816662190209508,
    "Y":    -0.0852272885558131,   # yield (wt% of dry feed)
    "L":     0.07555880522176552,  # lipids (wt% dry)
    "P":     0.0905428196778164,   # proteins (wt% dry)
    "C":    -0.021574160701441924, # carbohydrates (wt% dry)
    "Ash":  -0.006360842296037167, # ash (wt% dry)
}

TRACI_CATEGORIES = {
'GlobalWarming': dict(
    alias='GlobalWarmingPotential',
    method='TRACI 2.1',
    category='environmental impact',
    unit='kg CO2-eq'
),
'Acidification': dict(
    method='TRACI 2.1',
    category='environmental impact',
    unit='kg SO2-eq'
),
'Ecotoxicity': dict(
    method='TRACI 2.1',
    category='environmental impact',
    unit='CTUe'
),
'Eutrophication': dict(
    method='TRACI 2.1',
    category='environmental impact',
    unit='kg N-eq'
),
'OzoneDepletion': dict(
    method='TRACI 2.1',
    category='environmental impact',
    unit='kg CFC-11-eq'
),
'PhotochemicalOxidation': dict(
    method='TRACI 2.1',
    category='environmental impact',
    unit='kg O3-eq'
),
'Carcinogenics': dict(
    method='TRACI 2.1',
    category='environmental impact',
    unit='CTUh'
),
'NonCarcinogenics': dict(
    method='TRACI 2.1',
    category='environmental impact',
    unit='CTUh'
),
'ParticulateMatterFormation': dict(
    method='TRACI 2.1',
    category='environmental impact',
    unit='kg PM2.5-eq'
)
}
def _frac_lt350(Y_wt, L_dry, P_dry, C_dry, Ash_dry):
    """Logistic predictor for fraction of biocrude boiling <350°C (0–1)."""
    z = (_COEFS["const"]
         + _COEFS["Y"]   * Y_wt
         + _COEFS["L"]   * L_dry
         + _COEFS["P"]   * P_dry
         + _COEFS["C"]   * C_dry
         + _COEFS["Ash"] * Ash_dry)
    return 1.0 / (1.0 + math.exp(-z))

def cutoff_fracs_from_feed(feedstock_composition, HTL_yields, dewatered=True):
    """
    Returns [light, mid, heavy]:
      light = volatiles/water fraction (small for dewatered crude)
      mid   = <350°C fuel-range fraction (predicted)
      heavy = >350°C fraction
    Expects feedstock_composition keys: 'Water','Lipids','Proteins','Carbohydrates','Ash'
    with WET-basis values. HTL_yields['biocrude'] is a fraction (0–1) of dry feed.
    """
    W = feedstock_composition['Water']  # wet-basis moisture (0–1)
    # convert to DRY-basis wt% for the predictor
    denom = max(1e-9, 1.0 - W)
    L = feedstock_composition['Lipids']        / denom * 100.0
    P = feedstock_composition['Proteins']      / denom * 100.0
    C = feedstock_composition['Carbohydrates'] / denom * 100.0
    Ash = feedstock_composition['Ash']         / denom * 100.0
    Y_wt = HTL_yields['biocrude'] * 100.0

    # predict <350°C fraction
    mid = _frac_lt350(Y_wt, L, P, C, Ash)

    # light = tiny volatile/water cut; fixed if dewatered, else tie to moisture
    if dewatered:
        light = 0.05  # 5%
    else:
        light = max(0.0, min(0.10, 0.01 + 0.2 * W))

    # heavy is whatever remains
    heavy = 1.0 - light - mid

    # clamp & renormalize for safety
    light = max(0.0, min(0.20, light))
    mid   = max(0.0, min(0.99, mid))
    heavy = max(0.0, heavy)
    s = light + mid + heavy if (light + mid + heavy) > 0 else 1.0
    return [light/s, mid/s, heavy/s]

# def find_Lr_Hr(unit,
#                Lr_trial_range=np.arange(0.5, 1, 0.01),
#                Hr_trial_range=np.arange(0.5, 1, 0.01)):

#     import numpy as np, pandas as pd

#     outs0, _ = unit.outs
#     results = {}
#     F_mass_in = max(unit.F_mass_in, 1e-9)

#     # Save original values
#     orig_Lr, orig_Hr = unit.Lr, unit.Hr

#     for Lr in Lr_trial_range:
#         unit.Lr = float(np.clip(Lr, 1e-6, 0.9999))
#         Hr_results = {}

#         for Hr in Hr_trial_range:
#             unit.Hr = float(np.clip(Hr, 1e-6, 0.9999))

#             try:
#                 # RUN COLUMN
#                 unit._run()

#                 # FEASIBILITY CHECK: design AND cost
#                 unit._design()
#                 unit._cost()

#                 # compute valid yield
#                 yield_val = outs0.F_mass / F_mass_in
#                 yield_val = float(yield_val) if yield_val > 0 else np.nan

#                 Hr_results[unit.Hr] = yield_val
#                 print(f"Lr={unit.Lr:.3f}, Hr={unit.Hr:.3f}, Yield={yield_val:.3f}")

#             except Exception as e:
#                 Hr_results[unit.Hr] = np.nan
#                 print(f"Lr={Lr:.3f}, Hr={Hr:.3f}, failed: {e}")

#         results[unit.Lr] = Hr_results

#     # restore original
#     unit.Lr, unit.Hr = orig_Lr, orig_Hr
#     try: unit.simulate()
#     except: pass

#     # ANALYZE RESULTS
#     results_df = pd.DataFrame(results, dtype=float)
#     flat = results_df.stack().dropna()

#     if flat.empty:
#         print("[Optimization] No feasible Lr/Hr.")
#         return results_df, None, None, None

#     best_idx = flat.idxmax()
#     best_Hr, best_Lr = best_idx
#     best_yield = flat.max()

#     print(f"[Optimization] Best feasible Lr={best_Lr:.3f}, Hr={best_Hr:.3f}, Yield={best_yield:.3f}")
#     return results_df, float(best_Lr), float(best_Hr), float(best_yield)
def find_Lr_Hr(unit,
                Lr_trial_range=np.linspace(0.6, 1.0, 9),
                Hr_trial_range=np.linspace(0.8, 1.0, 9)):
    """
    Find (Lr, Hr) that maximize the top-product (biofuel) mass fraction
    of a ShortcutColumn. Safe against simulation and dtype issues.
    Returns (results_df, best_Lr, best_Hr, max_yield).
    """
    import numpy as np, pandas as pd

    outs0, _ = unit.outs
    results = {}
    F_mass_in = max(unit.F_mass_in, 1e-9)
    _Lr, _Hr = unit.Lr, unit.Hr

    for Lr in Lr_trial_range:
        unit.Lr = float(np.clip(Lr, 1e-6, 0.9999))
        Hr_results = {}
        for Hr in Hr_trial_range:
            unit.Hr = float(np.clip(Hr, 1e-6, 0.9999))
            try:
                unit._run_flag = False
                unit._design_results = {}
                unit._baseline_purchase_costs = {}
                unit.simulate()
                val = outs0.F_mass / F_mass_in
                val = float(val) if np.isfinite(val) and val > 0 else np.nan
                Hr_results[unit.Hr] = val
                print(f"Lr={unit.Lr:.3f}, Hr={unit.Hr:.3f}, Biofuel yield={val:.3f}" if np.isfinite(val) else f"Lr={unit.Lr:.3f}, Hr={unit.Hr:.3f}, failed")
            except Exception as e:
                Hr_results[unit.Hr] = np.nan
                print(f"Lr={Lr:.3f}, Hr={Hr:.3f}, failed: {e}")
        results[unit.Lr] = Hr_results

    # Assemble clean numeric DataFrame
    results_df = pd.DataFrame(results, dtype=float)

    # Restore column to its original state
    unit.Lr, unit.Hr = _Lr, _Hr
    try:
        unit.simulate()
    except Exception:
        pass

    # Flatten numeric results
    flat = results_df.stack()
    flat = flat[flat.notna() & np.isfinite(flat)]
    if flat.empty:
        print("[Optimization] No successful Lr/Hr combination found.")
        return results_df, None, None, None

    # Explicit numeric conversion for safety
    flat = pd.Series(flat.values.astype(float),
                      index=pd.MultiIndex.from_tuples(flat.index))

    best_idx = flat.idxmax()
    best_Hr, best_Lr = best_idx
    max_yield = flat.max()

    print(f"[Optimization] Best Lr={best_Lr:.3f}, Hr={best_Hr:.3f}, "
          f"Biofuel yield={max_yield:.3f}")
    return results_df, float(best_Lr), float(best_Hr), float(max_yield)

class BiobinderTEA(HTL_TEA):
    
    @property
    def land(self):
        if callable(self._land): return self._land()
        return self._land
    @land.setter
    def land(self, i):
        self._land = i