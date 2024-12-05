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

# Temporarily ignoring warnings
import warnings
warnings.filterwarnings('ignore')

import os, numpy as np, pandas as pd, qsdsan as qs
from qsdsan.utils import time_printer
from exposan.saf import (
    config_baseline,
    config_EC,
    config_EC_future,
    create_system,
    dry_flowrate as default_dry_flowrate,
    get_GWP,
    get_MFSP,
    results_path,
    )


# %%

# 110 tpd sludge (default) is sludge from a WWTP of about 100 MGD in size

@time_printer
def evaluation_across_sizes(ratios, **config_kwargs):
    fuel_yields = []
    MFSPs = []
    GWPs = []
    for ratio in ratios:
        dry_flowrate = ratio * default_dry_flowrate
        sys = create_system(dry_flowrate=dry_flowrate, **config_kwargs)
        mixed_fuel = flowsheet.stream.mixed_fuel
        print(f'ratio: {ratio}; dry flowrate: {dry_flowrate:.0f} kg/hr.')
        try:
            sys.simulate()
            fuel_yield = mixed_fuel.F_mass/dry_flowrate
            MFSP = get_MFSP(sys, print_msg=False)
            GWP = get_GWP(sys, print_msg=False)
            print(f'Fuel yield: {fuel_yield:.2%}; MFSP: ${MFSP:.2f}/GGE; GWP: {GWP:.2f} kg CO2e/GGE.\n')
        except: 
            print('Simulation failed.\n')
            fuel_yield = MFSP = GWP = None
        fuel_yields.append(fuel_yield)
        MFSPs.append(MFSP)
        GWPs.append(GWP)

    return fuel_yields, MFSPs, GWPs
        
if __name__ == '__main__':
    # config_kwargs = config_baseline
    # config_kwargs = config_EC
    config_kwargs = config_EC_future
    
    flowsheet = qs.main_flowsheet
    dct = globals()
    dct.update(flowsheet.to_dict())
    
    # ratios = [1]
    # ratios = np.arange(1, 11, 1).tolist()
    ratios = np.arange(0.1, 1, 0.1).tolist() + np.arange(1, 11, 1).tolist()
    sizes_results = evaluation_across_sizes(ratios=ratios, **config_kwargs)
    sizes_df = pd.DataFrame()
    sizes_df['Ratio'] = ratios
    sizes_df['Fuel yields'] = sizes_results[0]
    sizes_df['MFSP'] = sizes_results[1]
    sizes_df['GWP'] = sizes_results[2]
    outputs_path = os.path.join(results_path, f'sizes_{flowsheet.ID}.csv')
    sizes_df.to_csv(outputs_path)
