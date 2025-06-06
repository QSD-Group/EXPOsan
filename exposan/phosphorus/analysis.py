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

# Test if the ecorecovery module works
# from exposan import pm2_ecorecover
# pm2_ecorecover.load()
# sys = pm2_ecorecover.sys
# sys.simulate(t_span=(0,3), method='RK23')
# PBR = default_init_conds.PBR20
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

S_P_dct = { # steady state effluent P in mg/L, when Q=1000 m3/d
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


def create_flowsheets(f_dct=f_dct, inf_range=range(3, 16), Q=50): # inf=3-15 mg/L
    f_dct = f_dct or dict.fromkeys(f'f{n}' for n in inf_range) 
    for k, v in f_dct.items():
        mgL = k[1:]
        path = os.path.join(data_path, f'biowin_effluent_ph{mgL}.tsv')        
        df = qs.utils.load_data(path)
        df.Q = Q # update flow
        new_path = os.path.join(data_path, f'biowin_effluent_ph{mgL}_{Q}m3d.tsv')
        df.to_csv(new_path, sep='\t')
        # Create the EcoRecover system
        f = qs.Flowsheet(f'f{mgL}')
        pm2.create_system(flowsheet=f)
        SE = f.unit.SE
        SE._init_from_file(new_path)
        os.remove(new_path)
        f_dct[k] = f
    return f_dct

# f_dct = create_flowsheets() # this create all the flowsheets so takes some time


# %%

# Function to run all concentrations
def simulate_system(f_dct=None, inf_range=range(3, 16), Q=50, t_step=1, t=25, export=True):
    f_dct = create_flowsheets(f_dct, inf_range=inf_range, Q=Q)
    for k, v in f_dct.items():
        mgL = k[1:]
        f = f_dct[f'f{mgL}']
        sys = f.system.sys
        sys.reset_cache()
        sys.simulate(t_span=(0, t), method='RK23')
        inf = f.stream.Dynamic_influent
        Q = round(inf.get_total_flow('m3/d'))
        eff = f.stream.Effluent
        fig, axis = eff.scope.plot_time_series(('S_P'))
        if export:
            fig.savefig(os.path.join(figures_path, f'ecorecover_results_{mgL}mgL_{t}d_{Q}m3d.jpg'))
            sys.scope.export(os.path.join(results_path, f'ecorecover_results_{mgL}mgL_{t}d_{Q}m3d.xlsx'), t_eval=np.arange(0, t+t_step, t_step))
        print(f'\n\nfinished {k}\n\n')
    return f_dct


# Individual trial
f_dct = simulate_system(f_dct=None, inf_range=range(15, 16), Q=50, t_step=1, t=25, export=True)

# Run all results
# f_dct = simulate_system(f_dct=None, inf_range=range(3, 16), t_step=1, t=25, export=True)
