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

# !!! Temporarily ignoring warnings
import warnings
warnings.filterwarnings('ignore')

import os, numpy as np, pandas as pd, qsdsan as qs
from qsdsan.utils import time_printer
from exposan.saf import (
    results_path,
    create_system,
    get_MFSP,
    )


# %%

# 110 tpd sludge (default) is about 100 MGD

@time_printer
def MFSP_across_sizes(ratios, **config_kwargs):
    # sys = create_system(**config_kwargs)
    # feedstock = sys.flowsheet.stream.feedstock
    # _default_mass = feedstock.F_mass
    _default_mass = 18980.787227243676

    MFSPs = []
    fuel_yields = []
    for ratio in ratios:
        print(f'ratio: {ratio}')
        # sys.reset_cache() # too many fails
        flowsheet_ID = 'saf' # new flowsheet doesn't help... same as old flowsheet but new system
        if config_kwargs['include_PSA']: flowsheet_ID += '_PSA'
        if config_kwargs['include_EC']: flowsheet_ID += '_EC'
        flowsheet_ID += f'_{str(ratio).replace(".", "_")}'
        flowsheet = qs.Flowsheet(flowsheet_ID)
        qs.main_flowsheet.set_flowsheet(flowsheet)
        sys = create_system(flowsheet=flowsheet, **config_kwargs)
        feedstock = flowsheet.stream.feedstock
        mixed_fuel = flowsheet.stream.mixed_fuel
        FeedstockCond = flowsheet.unit.FeedstockCond

        new_mass = ratio * _default_mass
        feedstock.F_mass = new_mass
        dry_feedstock = feedstock.F_mass-feedstock.imass['H2O']
        FeedstockCond.feedstock_dry_flowrate = dry_feedstock
        try:
            sys.simulate()
            MFSP = get_MFSP(sys, print_msg=False)
            fuel_yield = mixed_fuel.F_mass/dry_feedstock
            print(f'MFSP: ${MFSP:.2f}/GGE; fuel yields {fuel_yield:.2%}.\n')
        except: 
            print('Simulation failed.\n')
            MFSP = fuel_yield = None
        MFSPs.append(MFSP)
        fuel_yields.append(fuel_yield)

    return MFSPs, fuel_yields
        
if __name__ == '__main__':
    config = {'include_PSA': False, 'include_EC': False,}
    # config = {'include_PSA': True, 'include_EC': False,}
    # config = {'include_PSA': True, 'include_EC': True,}
    flowsheet = qs.main_flowsheet
    dct = globals()
    dct.update(flowsheet.to_dict())
    
    # ratios = [1]
    # ratios = np.arange(1, 11, 1).tolist()
    ratios = np.arange(0.1, 1, 0.1).tolist() + np.arange(1, 11, 1).tolist()
    sizes_results = MFSP_across_sizes(ratios=ratios, **config)
    sizes_df = pd.DataFrame()
    sizes_df['Ratio'] = ratios
    sizes_df['MFSP'] = sizes_results[0]
    sizes_df['Fuel yields'] = sizes_results[1]
    outputs_path = os.path.join(results_path, f'sizes_{flowsheet.ID}.csv')
    sizes_df.to_csv(outputs_path)