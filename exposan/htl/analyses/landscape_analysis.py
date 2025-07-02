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