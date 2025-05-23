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
    figures_path,
    results_path,
    )

Temp = 286.08          # temperature [K]
f_dct = None

def create_flowsheets(f_dct=f_dct):
    f_dct = f_dct or dict.fromkeys(f'f{n}' for n in range(3, 16)) # inf=3-15 mg/L
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

# Test if the inputs work
f_dct = create_flowsheets()
f3 = f_dct['f3']
sys = f3.system.sys
sys.simulate(t_span=(0,10), method='RK23')
PBR20 = f3.unit.PBR20
fig, axis = PBR20.scope.plot_time_series(('S_P'))

# TODO: find appropriate initial states for different inf P
#!!! Should look for effluent rather than PBR20, tho results seem to be the same.
S_P_dct = { # mg/L
    '3': 0.103786314929568,
    '4': 0.757374963832081,
    '5': 1.38383912772804,
    '6': 2.03682213914134,
    '7': 2.68256462847412,
    '8': 3.32964385968618,
    '9': 3.9775093368395,
    '10': 4.68,
    '11': 5.33,
    '12': 5.98,
    '13': 6.63,
    '14': 7.28,
    '15': 7.93,
    }
#%%
# Individual trial
mgL = '3'
t = 25 # 3 for trials, 25 for full-length
t_step = 1
f = f_dct[f'f{mgL}']
sys = f.system.sys
sys.reset_cache()
sys.simulate(t_span=(0,t), method='RK23')
PBR20 = f.unit.PBR20
fig, axis = PBR20.scope.plot_time_series(('S_P'))
fig.savefig(os.path.join(figures_path, f'ecorecover_results_{mgL}mgL_{t}d.jpg'))
sys.scope.export(os.path.join(results_path, f'ecorecover_results_{mgL}mgL_{t}d.xlsx'), t_eval=np.arange(0, t+t_step, t_step))


# %%

def simulate_system(f_dct=None, t_step=1, t=25):
    f_dct = create_flowsheets(f_dct)
    for k, v in f_dct.items():
        mgL = k[1:]
        f = f_dct[f'f{mgL}']
        sys = f.system.sys
        sys.simulate(t_span=(0, t), method='RK23')
        path = os.path.join(results_path, f'ecorecover_results_{mgL}mgL_{t}d.xlsx')
        sys.scope.export(path, t_eval=np.arange(0, t+t_step, t_step))
        print(f'finished {k}')
        
    return f_dct

# f_dct = simulate_system(f_dct=None, t_step=1, t=3)

