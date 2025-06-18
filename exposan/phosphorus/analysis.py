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

'''
Simulation Findings:
    
    - Changing the input step did not change decrease simulation time,
    either of the results can be used as the baseline for the given concentration.
    
    original input vs. t_step=1 input, inf=15 mg/L, Q=50 m3/d, t=25 d, 5-6 min
    (ecorecover_results_15mgL_25d_50m3d_2025-6-12-10-3-19
     ecorecover_results_15mgL_25d_50m3d_2025-6-12-10-9-22)
    
    original input vs. t=0/1000 input, inf=13 mg/L, Q=50 m3/d, t=5 d, ~1 min
    (ecorecover_results_13mgL_5d_50m3d_2025-6-12-10-24-34
    ecorecover_results_13mgL_5d_50m3d_2025-6-12-10-30-52)
    
    Now only has t=0 and t=1000 d influent data after 2025-6-12-10-56-00
    
    - Cycling occured for 100-d simulations, though biomass appears to be recovering,
    the general trend is still to decrease.
    (ecorecover_results_15mgL_100d_50m3d_2025-6-12-10-41-47, took 21 min).
    
    - Simulating 13 mg/L for 5 days took around 1 min, 100 days took 20 min.
    (ecorecover_results_13mgL_5d_50m3d_2025-6-12-10-30-52
     ecorecover_results_13mgL_100d_50m3d_2025-6-12-11-17-43)
    
    Also observed cycling. 5-day results are identical to the first 5 days of 100-day results.
    
    - Doubling V_mix leads to earlier effluent P increase and higher effluent P
    compared to the original settings,
    but the simulation time almost halved <3 min vs. 5-6 min).
    (ecorecover_results_15mgL_25d_50m3d_2025-6-12-11-27-11
     ecorecover_results_15mgL_25d_50m3d_2025-6-12-10-9-22)
    
    - Doubling V_pbr greatly improves treatment performance,
    but leads to longer simulation time (8 min vs. 5-6 min).
    (ecorecover_results_15mgL_25d_50m3d_2025-6-12-11-54-4
     ecorecover_results_15mgL_25d_50m3d_2025-6-12-10-9-22)
    
    Running for 100 days with doubled V_pbr more than doubled the simulation time,
    still observed biomass cycling, but the trend is going up instead of down,
    significant increase in effluent P still observed around day 73-76 and day 95 towards the end.
    (ecorecover_results_15mgL_100d_50m3d_2025-6-12-12-42-52
     ecorecover_results_15mgL_25d_50m3d_2025-6-12-10-9-22)
    
    - Setting V_pbr based on influent P led to much better treatment performance.
    (ecorecover_results_15mgL_25d_50m3d_2025-6-12-14-12-9
     ecorecover_results_15mgL_100d_50m3d_2025-6-12-15-7-38)
    
    25-d simulation took 9-10 min, 100-d simulation took 59 min.
    
    - When V_pbr was adjusted based on the influent P,
    changing Q did not affect results in a significant way

    (ecorecover_results_15mgL_5d_1000m3d_2025-6-12-14-51-6
     ecorecover_results_15mgL_5d_50m3d_2025-6-12-14-57-32)
     
     The simulation time nearly doubled for 1000 m3/d (90 min vs. 50 min for 50 m3/d)
     ecorecover_results_15mgL_100d_50m3d_2025-6-12-15-7-38
     ecorecover_results_15mgL_100d_1000m3d_2025-6-12-16-56-23)
    
    - Tried to adjust initial condition setting, but does not appear to be working, gave up.
    
    - When using the calibrated PM2 parameters, simulation time shortened by >2/3.
    
    From 1-2 min to 23 sec for 15 mg/L, 5 day, 1000 m3/d.
    (ecorecover_results_15mgL_5d_1000m3d_2025-6-16-11-19-54
     ecorecover_results_15mgL_5d_1000m3d_2025-6-12-14-51-6)
    
    From 90 min to 26-27 min for 15 mg/L, 100 day, 1000 m3/d.
    (ecorecover_results_15mgL_100d_1000m3d_2025-6-16-11-46-34
     ecorecover_results_15mgL_100d_1000m3d_2025-6-12-16-56-23)
    
    Biomass seems to increase monotonically, getting to a plateau toward the end,
    effluent P decreases to 0 at Day 12.
    
    - Even with calibrated parameter, still observed the same cycling (and ~non treatment)
    without adjusting for the volume.
    
    (ecorecover_results_15mgL_5d_1000m3d_2025-6-16-11-56-17
     ecorecover_results_15mgL_100d_1000m3d_2025-6-16-11-46-34)
    
    - Simulate time doesn't make much difference in SRT (5 vs. 100 d for 15 mg/L, 1000 m3/d),
    but influent P significantly affect SRT,
    probably because the amount of biomass harvested was not adjusted.
    
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

import os, datetime, numpy as np, qsdsan as qs
from qsdsan.utils import get_SRT, time_printer
from exposan import pm2_ecorecover as pm2
from exposan.pm2_ecorecover.validation_result import pm2 as pm2_cal
from exposan.phosphorus import (
    data_path,
    figures_path,
    results_path,
    )

V_pbr_default = pm2.V_pbr
V_pbr_indv_default = V_pbr_default / 20
Q_default = 449.06
all_infs = range(3, 16) # inf=3-15 mg/L

# Initial conditions in EcoRecover
default_init_conds = {
        'X_CHL':2.53, # Chlorophyll content of cells, g Chl/m3
        'X_ALG':505.97, # Carbon-accumulating mixotrophic organisms, g COD/m3
        'X_PG':22.99, # Stored carbohydrates, g COD/m3
        'X_TAG':93.78, # Stored lipids, g COD/m3
        'S_CO2':30.0, # Soluble carbon dioxide, g CO2/m3
        'S_A':5.0, # Extracellular dissolved organic carbon (acetate), g COD/m3
        'S_G':5.0, # Extracellular dissolved organic carbon (glucose), g COD/m3
        'S_O2':5.0, # Dissolved oxygen, g O2/m3
        'S_NH':35.80, # Dissolved ammonium, g N/m3
        'S_NO':0.7, # Dissolved nitrate/nitrite, g N/m3
        'S_P':0.36, # Dissolved phosphorus, g P/m3
        'X_N_ALG':3.23, # Stored nitrogen in microalgal cell, g N/m3
        'X_P_ALG':0.19, # Stored phosphorus in microalgal cell, g P/m3
    }

# Set initial conditions based on BioWin effluent
init_cond_dct = {}
for mgL in all_infs:
    path = os.path.join(data_path, f'biowin_effluent_ph{mgL}.tsv')
    df = qs.utils.load_data(path)
    init_cond_dct[str(mgL)] = init_conds = default_init_conds.copy()
    init_conds['S_NH'] = df.S_NH[0]
    init_conds['S_NO'] = df.S_NO[0]
    init_conds['S_P'] = df.S_P[0]

# BioWin effluent P
# S_P_dct = {
#     '3': 0.104,
#     '4': 0.78,
#     '5': 1.42,
#     '6': 2.08,
#     '7': 2.73,
#     '8': 3.38,
#     '9': 4.03,
#     '10': 4.68,
#     '11': 5.33,
#     '12': 5.98,
#     '13': 6.63,
#     '14': 7.28,
#     '15': 7.93,
#     }

# Steady state effluent P in mg/L, when Q=1000 m3/d
# S_P_dct = {
#     '3': 0.103786314929568,
#     '4': 0.757374963832081,
#     '5': 1.38383912772804,
#     '6': 2.03682213914134,
#     '7': 2.68256462847412,
#     '8': 3.32964385968618,
#     '9': 3.9775093368395,
#     '10': 4.68,
#     '11': 5.33,
#     '12': 5.98,
#     '13': 6.63,
#     '14': 7.28,
#     '15': 7.93,
#     }


def create_flowsheets(inf_range=all_infs, Q=1000, init_conds=default_init_conds):
    f_dct = dict.fromkeys(f'f{n}' for n in inf_range) 
    for k, v in f_dct.items():
        mgL = k[1:]
        path = os.path.join(data_path, f'biowin_effluent_ph{mgL}.tsv')        
        df = qs.utils.load_data(path)
        df.Q = Q # update flow
        new_path = os.path.join(data_path, f'biowin_effluent_ph{mgL}_{Q}m3d.tsv')
        df.to_csv(new_path, sep='\t')
        # Create the EcoRecover system
        f = qs.Flowsheet(f'f{mgL}')
        pm2.create_system(flowsheet=f, init_conds=init_cond_dct[mgL])
        # Update system
        # To account for the increased P conc. and flowrate,
        # 449 m3/d was EcoRecover baseline
        X = max(1, Q/449*init_cond_dct[mgL]['S_P']/default_init_conds['S_P'])
        for u in f.unit:
            if hasattr(u, 'suspended_growth_model'):
                if u.suspended_growth_model is not None:
                    u.suspended_growth_model = pm2_cal
                if 'PBR' in u.ID:
                    if round(X, 1) > 1.1: u.V_max = V_pbr_indv_default * X
                u._compile_ODE()
        
        SE = f.unit.SE
        SE._init_from_file(new_path)
        os.remove(new_path)
        f_dct[k] = f
    return f_dct


# %%

# Function to run all concentrations
def simulate_system(f_dct, t_step=1, t=25, export=True):
    simulation_kwargs = dict(t_step=t_step, t=t, export=export)
    @time_printer
    def single_simulation(flowsheet, **simulation_kwargs):
        sys = f.system.sys
        sys.reset_cache()
        sys.simulate(t_span=(0, t), method='RK23')
        inf = f.stream.Dynamic_influent
        Q = round(inf.get_total_flow('m3/d'))
        eff = f.stream.Effluent
        fig, axis = eff.scope.plot_time_series(('S_P'))
        if export:
            now = datetime.datetime.now()
            time = f'{now.year}-{now.month}-{now.day}-{now.hour}-{now.minute}-{now.second}'
            fig.savefig(os.path.join(figures_path, f'ecorecover_results_{mgL}mgL_{t}d_{Q}m3d_{time}.jpg'))
            sys.scope.export(os.path.join(results_path, f'ecorecover_results_{mgL}mgL_{t}d_{Q}m3d_{time}.xlsx'), t_eval=np.arange(0, t+t_step, t_step))
        print(f'\n\nfinished {k}\n\n')
    for k, v in f_dct.items():
        mgL = k[1:]
        f = f_dct[f'f{mgL}']
        single_simulation(flowsheet=f, **simulation_kwargs)
    return f_dct


# Individual trial
mgL = 3
f_dct = create_flowsheets(inf_range=range(mgL, mgL+1), Q=1000)
f = f_dct[f'f{mgL}']
sys = f.system.sys
f_dct = simulate_system(f_dct=f_dct, t_step=1, t=5, export=False)

biomass_IDs = ('X_ALG', 'X_PG', 'X_TAG', 'X_N_ALG', 'X_P_ALG')
SRT = get_SRT(sys, biomass_IDs=biomass_IDs,
        wastage=f.stream.Harvested_biomass,
        # active_unit_IDs=(u.ID for u in sys.units if u.get_retained_mass(biomass_IDs) is not None)
        active_unit_IDs=('MIX', *(f'PBR{i}' for i in range(1, 21)))
        )
print(SRT) # 2-3 d for 3 mg/L, 75 d for 15 mg/L

# Adjust SRT by changing the amount of harvested biomass
POST_MEM = f.unit.POST_MEM
POST_MEM.split = 0.95 # split to return, default 0.97, smaller number means harvesting, thus shorter SRT
sys.reset_cache()
sys.simulate(t_span=(0, 5), method='RK23')
fig, axis = f.stream.Effluent.scope.plot_time_series(('S_P'))

# Run all results
# f_dct = create_flowsheets(inf_range=range(3, 16), Q=1000)
# f_dct = simulate_system(f_dct=f_dct, t_step=1, t=25, export=True)
