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

# Test if the module works
# from exposan import pm2_ecorecover
# pm2_ecorecover.load()
# sys = pm2_ecorecover.sys
# sys.simulate(t_span=(0,3), method='RK23')
# PBR = pm2_ecorecover.PBR20
# fig, axis = PBR.scope.plot_time_series(('S_P'))
# fig

# %%

import os, numpy as np, qsdsan as qs
from exposan import pm2_ecorecover as pm2
from exposan.phosphorus import (
    data_path,
    results_path,
    )

Temp = 286.08          # temperature [K]

def create_flowsheets():
    f_dct = dict.fromkeys(f'f{n}' for n in range(3, 16)) # inf=3-15 mg/L
    for k, v in f_dct.items():
        mgL = k[1:]
        path = os.path.join(data_path, f'biowin_effluent_ph{mgL}.tsv')
        # Create the EcoRecover system
        f = qs.Flowsheet(f'f{mgL}')
        pm2.create_system(flowsheet=f)
        SE = f.unit.SE
        SE._init_from_file(path)
        f_dct[k] = f
    return f_dct

f_dct = create_flowsheets()

# %%

def simulate_system(f_dct=None, t_step=1, t=25):
    f_dct = f_dct or create_flowsheets()
    for k, v in f_dct.items():
        mgL = k[1:]
        f = f_dct[f'f{mgL}']
        sys = f.system.sys
        sys.simulate(t_span=(0, t), method='RK23')
        path = os.path.join(results_path, f'ecorecover_results{mgL}.xlsx')
        sys.scope.export(path, t_eval=np.arange(0, t+t_step, t_step))
        print(f'finished {k}')
    return f_dct

f_dct = simulate_system(f_dct)

# f3 = f_dct['f3']
# PBR20 = f3.unit.PBR20
# fig, axis = PBR20.scope.plot_time_series(('S_P'))