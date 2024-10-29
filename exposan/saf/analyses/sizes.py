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
from exposan.saf import (
    results_path,
    create_system,
    get_MFSP,
    )


# %%

# 110 tpd sludge (default) is about 100 MGD
# _default_size = 100
# MGDs = np.arange(10, 100, 10).tolist() + np.arange(100, 1300, 100).tolist()

def MFSP_across_sizes(ratios, **config_kwargs):
    sys = create_system(**config_kwargs)
    sys.simulate()
    
    feedstock = sys.flowsheet.stream.feedstock
    FeedstockCond = sys.flowsheet.unit.FeedstockCond
    _default_dry_mass = feedstock.F_mass
    
    MFSPs = []
    # for MGD in [10, 100, 1000]:
    for ratio in ratios:
        new_dry_mass = ratio * _default_dry_mass
        feedstock.F_mass = new_dry_mass
        FeedstockCond.feedstock_dry_flowrate = feedstock.F_mass-feedstock.imass['H2O']
        print(ratio, new_dry_mass)
        try:
            sys.simulate()
            MFSPs.append(get_MFSP(sys, True))
        except: 
            print('Simulation failed.')
            MFSPs.append(None)
    return MFSPs
        
if __name__ == '__main__':
    # config = {'include_PSA': False, 'include_EC': False,}
    config = {'include_PSA': True, 'include_EC': False,}
    # config = {'include_PSA': True, 'include_EC': True,}
    flowsheet = qs.main_flowsheet
    dct = globals()
    dct.update(flowsheet.to_dict())
    
    
    ratios = np.arange(0.1, 1, 0.1).tolist() + np.arange(1, 10, 1).tolist()
    sizes_results = MFSP_across_sizes(sizes=ratios, **config)
    sizes_df = pd.DataFrame()
    sizes_df['Ratio'] = ratios
    sizes_df['MFSP'] = sizes_results
    outputs_path = os.path.join(results_path, f'sizes_{flowsheet.ID}.csv')
    sizes_df.to_csv(outputs_path)