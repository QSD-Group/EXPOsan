#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

#%% initialization

import numpy as np
from scipy.stats import qmc

#%% WWTP random sampling

# TODO: adjust typology if including AD and AeD (any any units before them)
# all use uniform distribution (not Monte Carlo)
# for now, assume WRRFs are only responsible for the transportation of e.g., solids (and ash) to landfills, biocrude to oil refineries
# for other products, assume they can be sold at the gate of WRRFs
typology_bounds = {
    # tonne/day
    'solids_mass_flow': (0, 250),
    # dry weight ratio
    'solids_dw_ash': (0.15, 0.45),
    # ash free dry weight ratio
    'solids_afdw_lipid': (0.05, 0.35),
    # ash free dry weight ratio
    'solids_afdw_protein': (0.3, 0.55),
    # $/kWh
    'electricity_cost': (0.05, 0.2),
    # kg CO2 eq/kWh
    'electricity_CI': (0, 1.5),
    # km
    'landfill_distance': (0, 1500),
    # km
    'biofuel_distance': (0, 1500)
    }

typology_keys = list(typology_bounds.keys())

# Latin Hypercube Sampling
sampler = qmc.LatinHypercube(d=len(typology_keys), seed=3221)
# TODO: decide the sample size based on the computation requirement
LHS_unit = sampler.random(n=10000)

typology_mins = np.array([typology_bounds[key][0] for key in typology_keys])
typology_maxs = np.array([typology_bounds[key][1] for key in typology_keys])
LHS_samples = qmc.scale(LHS_unit, typology_mins, typology_maxs)

typology = {key: LHS_samples[:, i] for i, key in enumerate(typology_keys)}

lipid = typology['solids_afdw_lipid']
protein = typology['solids_afdw_protein']
carbohydrate = 1 - lipid - protein
# no negative carb
mask = carbohydrate >= 0
typology = {k: v[mask] for k, v in typology.items()}
typology['solids_afdw_carbohydrate'] = carbohydrate[mask]

for i in range(len(LHS_samples)):
    solids_mass_flow = typology['solids_mass_flow'][i]
    solids_dw_ash = typology['solids_dw_ash'][i]
    solids_afdw_lipid = typology['solids_afdw_lipid'][i]
    solids_afdw_protein = typology['solids_afdw_protein'][i]
    solids_afdw_carbohydrate = typology['solids_afdw_carbohydrate'][i]
    electricity_cost = typology['electricity_cost'][i]
    electricity_CI = typology['electricity_CI'][i]
    landfill_distance = typology['landfill_distance'][i]
    biofuel_distance = typology['biofuel_distance'][i]
    
    # TODO: compare the baseline cost and GHG results from different systems
    # TODO: for each typology, identify a winner
    # TODO: plot 5th and 95th percentiles to represent the opportunity space

#%% preliminary utility check

from exposan.htl import (create_C1_system, create_C2_system, create_C3_system,
                         create_C4_system, create_C5_system, create_C6_system,
                         create_C7_system, create_C8_system, create_C9_system,
                         create_C10_system, create_C11_system, create_C12_system,
                         create_C13_system, create_C14_system, create_C15_system,
                         create_C16_system, create_C17_system, create_C18_system,
                         create_C19_system, create_C20_system, create_C21_system,
                         create_C22_system, create_C23_system, create_C24_system,
                         create_C25_system, create_T1_system, create_T2_system,
                         create_T3_system, create_T4_system, create_T5_system,
                         create_T6_system, create_T7_system, create_T8_system,
                         create_T9_system, create_T10_system, create_T11_system,
                         create_T12_system, create_T13_system, create_T14_system,
                         create_T15_system)

results = []
for function in (create_C1_system, create_C2_system, create_C3_system,
                 create_C4_system, create_C5_system, create_C6_system,
                 create_C7_system, create_C8_system, create_C9_system,
                 create_C10_system, create_C11_system, create_C12_system,
                 create_C13_system, create_C14_system, create_C15_system,
                 create_C16_system, create_C17_system, create_C18_system,
                 create_C19_system, create_C20_system, create_C21_system,
                 create_C22_system, create_C23_system, create_C24_system,
                 create_C25_system, create_T1_system, create_T2_system,
                 create_T3_system, create_T4_system, create_T5_system,
                 create_T6_system, create_T7_system, create_T8_system,
                 create_T9_system, create_T10_system, create_T11_system,
                 create_T12_system, create_T13_system, create_T14_system,
                 create_T15_system):
    sys = function(size=10)
    print('\n' + sys.ID)
    print(round(sys.get_cooling_duty()))
    print(round(sys.get_heating_duty()))
    print(sys.heat_utilities)

#%% preliminary TEA check

from exposan.htl import (create_C1_system, create_C2_system, create_C3_system,
                         create_C4_system, create_C5_system, create_C6_system,
                         create_C7_system, create_C8_system, create_C9_system,
                         create_C10_system, create_C11_system, create_C12_system,
                         create_C13_system, create_C14_system, create_C15_system,
                         create_C16_system, create_C17_system, create_C18_system,
                         create_C19_system, create_C20_system, create_C21_system,
                         create_C22_system, create_C23_system, create_C24_system,
                         create_C25_system, create_T1_system, create_T2_system,
                         create_T3_system, create_T4_system, create_T5_system,
                         create_T6_system, create_T7_system, create_T8_system,
                         create_T9_system, create_T10_system, create_T11_system,
                         create_T12_system, create_T13_system, create_T14_system,
                         create_T15_system)

TEA_results = []
for function in (create_C1_system, create_C2_system, create_C3_system,
                 create_C4_system, create_C5_system, create_C6_system,
                 create_C7_system, create_C8_system, create_C9_system,
                 create_C10_system, create_C11_system, create_C12_system,
                 create_C13_system, create_C14_system, create_C15_system,
                 create_C16_system, create_C17_system, create_C18_system,
                 create_C19_system, create_C20_system, create_C21_system,
                 create_C22_system, create_C23_system, create_C24_system,
                 create_C25_system, create_T1_system, create_T2_system,
                 create_T3_system, create_T4_system, create_T5_system,
                 create_T6_system, create_T7_system, create_T8_system,
                 create_T9_system, create_T10_system, create_T11_system,
                 create_T12_system, create_T13_system, create_T14_system,
                 create_T15_system):
    sys = function(size=10)
    
    print('\n' + sys.ID)
    
    print(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
    
    TEA_results.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
    
print(np.quantile(TEA_results[0:25], 0.05))
print(np.quantile(TEA_results[0:25], 0.5))
print(np.quantile(TEA_results[0:25], 0.95))

print(np.quantile(TEA_results[25:40], 0.05))
print(np.quantile(TEA_results[25:40], 0.5))
print(np.quantile(TEA_results[25:40], 0.95))

#%% preliminary LCA check

from exposan.htl import (create_C1_system, create_C2_system, create_C3_system,
                         create_C4_system, create_C5_system, create_C6_system,
                         create_C7_system, create_C8_system, create_C9_system,
                         create_C10_system, create_C11_system, create_C12_system,
                         create_C13_system, create_C14_system, create_C15_system,
                         create_C16_system, create_C17_system, create_C18_system,
                         create_C19_system, create_C20_system, create_C21_system,
                         create_C22_system, create_C23_system, create_C24_system,
                         create_C25_system, create_T1_system, create_T2_system,
                         create_T3_system, create_T4_system, create_T5_system,
                         create_T6_system, create_T7_system, create_T8_system,
                         create_T9_system, create_T10_system, create_T11_system,
                         create_T12_system, create_T13_system, create_T14_system,
                         create_T15_system)

LCA_results = []
for function in (create_C1_system, create_C2_system, create_C3_system,
                 create_C4_system, create_C5_system, create_C6_system,
                 create_C7_system, create_C8_system, create_C9_system,
                 create_C10_system, create_C11_system, create_C12_system,
                 create_C13_system, create_C14_system, create_C15_system,
                 create_C16_system, create_C17_system, create_C18_system,
                 create_C19_system, create_C20_system, create_C21_system,
                 create_C22_system, create_C23_system, create_C24_system,
                 create_C25_system, create_T1_system, create_T2_system,
                 create_T3_system, create_T4_system, create_T5_system,
                 create_T6_system, create_T7_system, create_T8_system,
                 create_T9_system, create_T10_system, create_T11_system,
                 create_T12_system, create_T13_system, create_T14_system,
                 create_T15_system):
    sys = function(size=10)
    
    print('\n' + sys.ID)
    
    print(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
    LCA_results.append(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
print(np.quantile(LCA_results[0:25], 0.05))
print(np.quantile(LCA_results[0:25], 0.5))
print(np.quantile(LCA_results[0:25], 0.95))

print(np.quantile(LCA_results[25:40], 0.05))
print(np.quantile(LCA_results[25:40], 0.5))
print(np.quantile(LCA_results[25:40], 0.95))