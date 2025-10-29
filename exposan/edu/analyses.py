# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from qsdsan.utils import get_SRT
from exposan.edu import (
    create_system,
    biomass_IDs,
    active_unit_IDs,
    )

import numpy as np, pandas as pd

t = 20    
t_step = 1
method = "RK23"
Q_was_min=50
Q_was_max=11400
 
#%% At a single SRT 

def single_run(t, t_step, method=None, **kwargs):
    sys = create_system(**kwargs)
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=os.path.join(results_path, f'sol_{t}d_{method}_analyses_test.xlsx'),
        )
    
    # Simulation results
    srt = get_SRT(sys, biomass_IDs, sys.streams[-1], active_unit_IDs)
    # print(f'Estimated SRT assuming at steady state is {round(srt, 3)} days')

    ts = sys.scope.time_series
    ys = sys.scope._get_records()
    data = sys.scope._interpolate_eval(ts, ys, np.arange(0, t+t_step, t_step))
    df = pd.DataFrame(data.T, columns=sys.scope._get_headers())
    
    # Multiindex pre-processing
    lvl0 = df.columns.get_level_values(0)                   # WasteStream/SanUnit
    lvl1 = df.columns.get_level_values(1).to_series()       # Component + unit
    lvl1 = df.columns.get_level_values(1).to_series().str.replace(r"\s*\[.*?\]\s*", "", regex=True) # remove units
    df.columns = pd.MultiIndex.from_arrays([lvl0, lvl1])        # reassign
    
    # TSS/COD/BOD/TN/TP calculation
    idx_tracker = ['influent', 'effluent', 'recycle', 'wastage', 'O1', 'C1', 'S1']  #lvl0
    
    cmps = sys.units[0].components
    components = list(cmps.IDs)                                                     #lvl1

    X_comps = ['X_I','X_S','X_BH','X_BA','X_P','X_ND']         # TSS (only X components)
    w_mass  = np.array([0.75, 0.75, 0.75, 0.75, 0.75, 0.0])    # i_mass of X components
    i_COD = cmps.i_COD                      # COD
    f_BOD5_COD = cmps.f_BOD5_COD            # BOD
    i_N = cmps.i_N                          # TN
    i_P = cmps.i_P                          # TP

    # Calculated properties
    for stream in idx_tracker:
        # TSS (particulate components only)
        X = df.loc[:, pd.MultiIndex.from_product([[stream], X_comps])].to_numpy(dtype=float)
        df[(stream, "TSS")] = (X * w_mass).sum(axis=1)

        # COD, BO, TN, TP: particulate + soluble
        XS = df.loc[:, pd.MultiIndex.from_product([[stream], components])].to_numpy(dtype=float)
        df[(stream, "COD")] = (XS * i_COD).sum(axis=1)
        df[(stream, "BOD")] = (XS * (i_COD * f_BOD5_COD)).sum(axis=1)
        df[(stream, "TN")]  = (XS * i_N).sum(axis=1)
        df[(stream, "TP")]  = (XS * i_P).sum(axis=1)

    # Insert calculated properties under same upper level index
    cal_properties = ["TSS","COD","BOD","TN","TP"]
    new_cols = []
    for stream in idx_tracker:
        sub = list(df[stream].columns)
        base = [prop for prop in sub if prop not in cal_properties]
        ordered = base + cal_properties                    
        new_cols.extend([(stream, prop) for prop in ordered])
    df = df.reindex(columns=pd.MultiIndex.from_tuples(new_cols))

    return srt, df

#%% Multiple SRTs

def multiple_run(t, t_step, Q_was_min, Q_was_max, method=None, **kwargs):
    Q_WAS = []
    SRT = []
    effluent_substrate = []
    sludge_production = []
    mlss = []
    cal_yield = []

    for q_was in np.round(np.logspace(np.log10(Q_was_min), np.log10(Q_was_max), 8)).astype(int):
        srt, df = single_run(t, t_step, method=method, Q_was=q_was, **kwargs)
        last = lambda stream, comp: df[stream][comp].iloc[-1]
            
        Q_WAS.append(q_was)
        SRT.append(srt)
        effluent_substrate.append(last('effluent', 'S_S') + last('effluent', 'X_S'))         # mg COD/L
        sludge_production.append(q_was*0.75/1000*
                                 (last('wastage', 'X_BH') + last('wastage', 'X_BA') + last('wastage', 'X_I')))      # kg/d
        mlss.append(last('O1','X_BH') + last('O1','X_BA') + last('O1','X_I'))                                       # mg COD/L
        cal_yield.append((q_was*last('wastage', 'X_BH'))/
                         (last('influent','Q') * (last('influent','S_S') + last('influent','X_S')) -
                          last('effluent','Q') * (last('effluent','S_S') + last('effluent','X_S'))))   # unitless
        
    result = pd.DataFrame({'Q_WAS': Q_WAS,
                           'SRT': SRT,
                           'Effluent Substrate': effluent_substrate,
                           'Sludge Production': sludge_production,
                           'MLSS': mlss,
                           'Yield': cal_yield
                           })
    return result

#%% Run
single_result = single_run(t, t_step, method=method)                                        # Run simulation at a specific SRT
multiple_result = multiple_run(t, t_step, Q_was_min, Q_was_max, method=method)              # Run simulations at multiple SRTs

#%% Tested (what users can change)

# single_run(t, t_step, method=method)                                                      # pass (default)
# single_run(t, t_step, method=method, Q_was=500)                                           # pass (can change Q_was)
# single_run(t, t_step, method=method, Q_was=20000)                                         # pass (but too high Q_was will return error)
# single_run(t, t_step, method=method, Q_was=500, Temp=10)                                  # pass (can change temperature)
# single_run(t, t_step, method=method, Q_was=500, Temp=10, KLa=50)                          # pass (can change KLa)
# single_run(t, t_step, method=method, Q_was=500, Temp=50, KLa=50)                          # pass (temperature should be either 10 or 20)
# single_run(t, t_step, method=method, Q_was=500, Temp=10, KLa=50, asm_kwargs = dict(
#     Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
#     mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3, 
#     eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5, 
#     K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75, 
#     ))                                                                                    # pass (can change ASM1 parameters)
# single_run(t, t_step, method=method, Q_was=500, Temp=10, KLa=50, asm_kwargs = dict(
#     Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
#     mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3, 
#     eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5, 
#     K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75, 
#     ), inf_kwargs = {
#         'concentrations': {                                            
#             'S_S':0,
#             'X_BH':28.17,
#             'X_S':202.32,
#             'X_I':51.2,
#             'S_NH':31.56,
#             'S_I':30,
#             'S_ND':6.95,
#             'X_ND':10.59,
#             'S_ALK':7*12,
#           },
#         'units': ('m3/d', 'mg/L'),                
#         }
# )                                                                                         # pass (can change influent composition)