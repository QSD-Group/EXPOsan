# -*- coding: utf-8 -*-

# Import packages
# import qsdsan as qs
import numpy as np, pandas as pd
from qsdsan import sanunits as su, processes as pc, WasteStream, System
from qsdsan.utils import time_printer, load_data, get_SRT
from exposan.edu import data_path
import os

import warnings
warnings.filterwarnings('ignore')  
 
#%% Components
cmps = pc.create_asm1_cmps()          

#%% Parameters (flowrates, temperature)
Q_inf = 18446                               # influent flowrate [m3/d]
Q_was = 400                                 # sludge wastage flowrate [m3/d]
Q_ext = 18446                               # external recycle flowrate [m3/d]

Temp = 273.15+20                            # temperature [K]

#%% Create influent, effluent, recycle stream
influent = WasteStream('influent', T=Temp)
effluent = WasteStream('effluent', T=Temp)
solid = WasteStream('solid', T=Temp)

ext_recycle = WasteStream('recycle', T=Temp)
wastage = WasteStream('wastage', T=Temp)                     

#%% Set the influent composition
default_inf_kwargs = {                        # default influent composition
    'concentrations': {                                            
        'S_S':69.5,
        'X_BH':28.17,
        'X_S':202.32,
        'X_I':51.2,
        'S_NH':31.56,
        'S_I':30,
        'S_ND':6.95,
        'X_ND':10.59,
        'S_ALK':7*12,
      },
    'units': ('m3/d', 'mg/L'),                # ('input total flowrate', 'input concentrations')
    }                                                            

influent.set_flow_by_concentration(Q_inf, **default_inf_kwargs)   

#%% Parameters (volumes)
V_ae = 1333                                 # aerated tank volume [m3/d]
aer = pc.DiffusedAeration('aer', DO_ID='S_O', KLa=240, DOsat=8.0, V=V_ae)             # aeration model 

#%% ASM1
asm1 = pc.ASM1(Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
               mu_H=6.0, K_S=20.0, K_O_H=0.2, K_NO=0.5, b_H=0.62, 
               eta_g=0.8, eta_h=0.4, k_h=3.0, K_X=0.03, mu_A=0.8, 
               K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.08, fr_SS_COD=0.75)   # 20 degree C

# asm1 = pc.ASM1(Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
#                mu_H=3.0, K_S=20.0, K_O_H=0.2, K_NO=0.5, b_H=0.2, 
#                eta_g=0.8, eta_h=0.4, k_h=1.0, K_X=0.01, mu_A=0.3, 
#                K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.04, fr_SS_COD=0.75, )   # 10 degree C

# asm1 = pc.ASM1(Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
#                mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3, 
#                eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5, 
#                K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,  )   # bsm1 default

#%% System configuration
O1 = su.CSTR('O1', ins=[influent, ext_recycle], V_max=V_ae, aeration=aer,
             DO_ID='S_O', suspended_growth_model=asm1)

C1 = su.IdealClarifier('C1', ins=O1-0, outs=[effluent, solid],
                       sludge_flow_rate=Q_ext+Q_was, solids_removal_efficiency=1)

S1 = su.Splitter('S1', ins=C1-1, outs=[ext_recycle, wastage], split=1-Q_was/Q_ext)

# Create system
sys = System('example_system', path=(O1, C1, S1), recycle=(ext_recycle))    
sys.diagram(display=False)

#%% Import initial condition excel file
df = load_data(os.path.join(data_path, 'initial_conditions_asm1.xlsx'), sheet='default')

#%% Create a function to set initial conditions of the reactors
def batch_init(sys, df):
    dct = df.to_dict('index')                                        
    u = sys.flowsheet.unit                                            
    u.O1.set_init_conc(**dct['O1'])    

    # c1s = {k:v for k,v in dct['C1_s'].items() if v>0}                
    # c1x = {k:v for k,v in dct['C1_x'].items() if v>0}
    # tss = [v for v in dct['C1_tss'].values() if v>0]
    # u.C1.set_init_solubles(**c1s)                             # set solubles
    # u.C1.set_init_sludge_solids(**c1x)                        # set sludge solids
    # u.C1.set_init_TSS(tss)   

batch_init(sys,df)

#%% Simulation settings
sys.set_dynamic_tracker(influent, effluent, ext_recycle, wastage, O1, C1, S1)          
sys.set_tolerance(rmol=1e-6)

biomass_IDs = ('X_BH', 'X_BA')
active_unit_IDs = ("O1", )
# Simulation settings
t = 50                         # total time for simulation
t_step = 1                      # times at which to store the computed solution

# method = 'BDF'                  # integration method to use
# method = 'RK45'
method = 'RK23'
# method = 'DOP853'
# method = 'Radau'
# method = 'LSODA'

# Run simulation, this could take several minuates
sys.simulate(state_reset_hook='reset_cache',
             t_span=(0,t),
             t_eval=np.arange(0, t+t_step, t_step),
             method=method,
             export_state_to=f'results/sol_{t}d_{method}_ideal_{Q_was}_cal_properties.xlsx',               
            )

srt = get_SRT(sys, biomass_IDs, wastage, active_unit_IDs)
print(f'Estimated SRT assuming at steady state is {round(srt, 8)} days')
#%%
ts = sys.scope.time_series
ys = sys.scope._get_records()
data = sys.scope._interpolate_eval(ts, ys, np.arange(0, t+t_step, t_step))
df = pd.DataFrame(data.T, columns=sys.scope._get_headers())

lvl0 = df.columns.get_level_values(0)                   # WasteStream/SanUnit
lvl1 = df.columns.get_level_values(1).to_series()       # Component + unit
lvl1 = df.columns.get_level_values(1).to_series().str.replace(r"\s*\[.*?\]\s*", "", regex=True) # remove units
df.columns = pd.MultiIndex.from_arrays([lvl0, lvl1])        # reassign

#%%
idx_tracker = ['influent', 'effluent', 'recycle', 'wastage', 'O1', 'C1', 'S1']  #lvl0
components = list(cmps.IDs)                                                         #lvl1

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

#%%
# Q_was = 50, SRT=13.33 d
# Q_was = 100, SRT=6.665 d
# Q_was = 200, SRT=3.333 d
# Q_was = 400, SRT=1.667 d
# Q_was = 1000, SRT=0.667 d
# Q_was = 2000, SRT=0.335 d
# Q_was = 2200, SRT=0.305 d
# Q_was = 2500, SRT=0.269 d
# Q_was = 5000, SRT=0.138 d
# Q_was = 10000, SRT= 0.078 d
# Q_was = 11000, SRT= 0.074 d
# Q_was = 11300, SRT= 0.073 d

# Q_was = 11500, SRT= d (error: effluent empty)
# Q_was = 12000, SRT= d (error: effluent empty)
# Q_was = 15000, SRT= d (error: effluent empty)

#%%
# O1.scope.plot_time_series(('X_BH', 'X_BA')) 
# C1._ODE = None
# effluent.COD
# b=O1.scope.record
# c=O1.scope.time_series
# d=O1.scope.header
# C1.wastage
# C1._ODE

# X_BA (Autotrophic nitrifiers): oxidizing ammonia directly to nitrate under aerobic conditions.
# X_BH (Heterotrophic denitrifiers): reducing nitrate to nitrogen gas in anoxic (oxygen-deficient) conditions. Heterotrophs use the nitrate as an electron acceptor to break down organic carbon.