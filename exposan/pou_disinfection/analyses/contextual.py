#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, pandas as pd
from qsdsan import PowerUtility, ImpactItem
from qsdsan.utils import load_data
from exposan.pou_disinfection import (
    create_system,
    data_path,
    ppl,
    results_path,
    update_number_of_householdsize,
    )

# Data for
water_path = os.path.join(data_path, '_raw_water.xlsx')
contextual_data = load_data(water_path, sheet='contextual')
results = pd.DataFrame(index=contextual_data.index)

# Create systems
sysA = create_system('A')
sysB = create_system('B')
sysC = create_system('C')
sysD = create_system('D')

all_syses = (sysA, sysB, sysC, sysD)
for sys in all_syses:
    sys.units[0].Mg = 0 # assume all hardness is Ca
    
e_item = ImpactItem.get_item('Electricity')

all_results = {}
for ID in (sys.ID for sys in all_syses):
    all_results[f'{ID}_cost'] = {}
    all_results[f'{ID}_gwp'] = {}

# Function to update data and record results
for community, data in contextual_data.iterrows():
    PowerUtility.price = data.e_cost
    e_item.CFs['GWP'] = data.e_gwp
    for sys in all_syses:
        update_number_of_householdsize(sys, household_size=data.household_size, ppl=ppl)
        RawWater = sys.units[0]
        RawWater.Ecoli = data.e_coli
        RawWater.turbidity = data.turbidity
        RawWater.Ca = data.hardness # assume all hardness is Ca
        sys.simulate()
        all_results[f'{sys.ID}_cost'][community] = sys.TEA.EAC/ppl
        all_results[f'{sys.ID}_gwp'][community] = sys.LCA.total_impacts['GWP']/(sys.LCA.lifetime*ppl)

results = pd.DataFrame.from_dict(all_results)
results.to_csv(os.path.join(results_path, 'contextual_analysis.csv'))
