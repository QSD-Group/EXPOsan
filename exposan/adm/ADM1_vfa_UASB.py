# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 19:57:56 2023

@author: Junhyung Park
"""

import numpy as np
import matplotlib.pyplot as plt
from chemicals.elements import molecular_weight as get_mw
from qsdsan import processes as pc, WasteStream, System
# from qsdsan.utils import time_printer

from exposan.metab import UASB

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)        # to ignore Pandas future warning

#%%
# Components
cmps = pc.create_adm1_vfa_cmps()      # create state variables for ADM1_vfa
# cmps.show()                         # 30 components in ADM1_vfa + water

#%%
# Processes
adm1 = pc.ADM1_vfa()                     # create ADM1 processes
# adm1.show()                            # 22 processes in ADM1

# adm1.__dict__                          # adm1 is composed of...

# Petersen stoichiometric matrix
# adm1.stoichiometry

#%%
# Kinetics
# rhos = pc.rhos_adm1_vfa()
# Flow rate, temperature, HRT (R3G20)
# Q = 0.00007                                       # influent flowrate [m3/d]
Q = 7                                               #!!! increasing Q shouldn't affect process simulation, but it'd increase numerical stability
Temp = 273.15 + 40                                  # temperature [K]
HRT = 20                                            # HRT [d]

#%%
# WasteStream
inf = WasteStream('Influent', T=Temp)               # influent
eff = WasteStream('Effluent', T=Temp)               # effluent
gas = WasteStream('Biogas')                         # gas

#%%
# Set influent concentration (growth medium (glucose))
C_mw = get_mw({'C':1})        # molecular weight of carbon
N_mw = get_mw({'N':1})        # molecular weight of nitrogen

# Default values
'''
default_inf_kwargs = {
    'concentrations': {
        'S_su':0.01,                                              
        'S_aa':1e-3,
        'S_fa':1e-3,
        'S_la':1e-3,
        'S_et':1e-3,
        'S_va':1e-3,
        'S_bu':1e-3,
        'S_pro':1e-3,
        'S_ac':1e-3,
        'S_h2':1e-8,
        'S_ch4':1e-5,
        'S_IC':0.04*C_mw,                                             
        'S_IN':0.01*N_mw,
        'S_I':0.02,
        'X_c':2.0,
        'X_ch':5.0,
        'X_pr':20.0,
        'X_li':5.0,
        'X_aa':1e-2,
        'X_fa':1e-2,
        'X_la':1e-2,
        'X_et':1e-2,
        'X_c4':1e-2,
        'X_pro':1e-2,
        'X_ac':1e-2,
        'X_h2':1e-2,
        'X_I':25,
        'S_cat':0.04,
        'S_an':0.02,
        },
    'units': ('m3/d', 'kg/m3'),                                 
    }                                                           # concentration of each state variable in influent
'''
# Medium (Glucose)
default_inf_kwargs = {
    'concentrations': {
        'S_su':20.14,                                               # glucose 10 g/L = 10.7 kg/m3, 20 g/L = 20.14 kg/m3
        'S_aa':0,
        'S_fa':0,
        'S_la':0,
        'S_et':0,
        'S_va':0,
        'S_bu':0,
        'S_pro':0,
        'S_ac':0,
        'S_h2':0,
        'S_ch4':0,
        'S_IC':0*C_mw,                                             
        'S_IN':5.553*1e-5*N_mw,                                     #!!! S_IN: 4.283*1e-5 kmole N / m3? * N_mw, why different? Fixed value
        'S_I':0,
        'X_c':0,
        'X_ch':0,
        'X_pr':0,
        'X_li':0,
        'X_aa':0,
        'X_fa':0,
        'X_la':0,
        'X_et':0,
        'X_c4':0,
        'X_pro':0,
        'X_ac':0,
        'X_h2':0,
        'X_I':0,
        'S_cat':1e-5,
        'S_an':1e-5,
        },
    'units': ('m3/d', 'kg/m3'),                                 #!!! kg/m3 = g/L, Is it kg COD / m3?
    }                                                           # concentration of each state variable in influent
inf.set_flow_by_concentration(Q, **default_inf_kwargs)          # set influent concentration
# inf

#%%
# SanUnit
U1 = UASB('UASB', ins=inf, outs=(gas, eff), model=adm1,        # !!!Even though my model does not contain recirculation ratio, is it defined as CSTR? or PFR regarding HRT?
          V_liq=Q*HRT, V_gas=Q*HRT*0.1,                        # !!! Considering real experiments including either high recirculation rate or not
          T=Temp, pH_ctrl=False,                               # pH adjustment X
          fraction_retain=0.95,                                 # needs to set this value properly
          )                                                    

                                                               # fraction_retain : float, optional
                                                               # The assumed fraction of ideal retention of select components. The default is 0.95.
                                                               # To make all solids sent to effluent

# U1                                                             # anaerobic CSTR with influent, effluent, and biogas
                                                               # before running the simulation, 'outs' have nothing
# print(f"The liquid volume of the reactor is: {U1.V_liq} m^3")

# Set initial condition of the reactor (Cow manure (Inoculum) in bioreactor, 10% of total volume)
# Default values
'''
default_init_conds = {
    'S_su': 0.0124*1e3,                                        
    'S_aa': 0.0055*1e3,
    'S_fa': 0.1074*1e3,
    'S_la': 0.0,
    'S_et': 0.0,
    'S_va': 0.0123*1e3,
    'S_bu': 0.0140*1e3,                                  
    'S_pro': 0.0176*1e3,                                  
    'S_ac': 0.0893*1e3,                                      
    'S_h2': 2.5055e-7*1e3,
    'S_ch4': 0.0555*1e3,
    'S_IC': 0.0951*C_mw*1e3,
    'S_IN': 0.0945*N_mw*1e3,
    'S_I': 0.1309*1e3,
    'X_ch': 0.0205*1e3,
    'X_pr': 0.0842*1e3,
    'X_li': 0.0436*1e3,
    'X_su': 0.3122*1e3,
    'X_aa': 0.9317*1e3,
    'X_fa': 0.3384*1e3,
    'X_la': 0.0,
    'X_et': 0.0,
    'X_c4': 0.3258*1e3,
    'X_pro': 0.1011*1e3,
    'X_ac': 0.6772*1e3,
    'X_h2': 0.2848*1e3,
    'X_I': 17.2162*1e3
    }                   # in mg/L                         
'''
# Inoculum (Cow manure)
default_init_conds = {
    'S_su': 0.0124*1e3,                                        # fixed according to R3G20 (No glucose)
    'S_aa': 0.0055*1e3,
    'S_fa': 0.1074*1e3,
    'S_la': 0.0,
    'S_et': 0.0,
    'S_va': 0.0,
    'S_bu': 0.0,                                  # fixed according to R3G20
    'S_pro': 0,                                  # fixed according to R3G20
    'S_ac': 0,                                      # fixed according to R3G20
    'S_h2': 0,
    'S_ch4': 0,
    'S_IC': 0,
    'S_IN': 0,
    'S_I': 0,
    'X_ch': 0.5*1e2,
    'X_pr': 0.5*1e2,
    'X_li': 0.5*1e2,
    'X_su': 1.5*1e1,
    'X_aa': 1.5*1e1,
    'X_fa': 1.5*1e1,
    'X_la': 0.01*1e1,
    'X_et': 0.1*1e1,
    'X_c4': 2.0*1e2,
    'X_pro': 2.0*1e2,
    'X_ac': 4.0*1e3,
    'X_h2': 1.0*1e3,
    'X_I': 1.0*1e3
    }                   # in mg/L                         #!!! Is it also mg COD/L?

U1.set_init_conc(**default_init_conds)                          # set initial condition of AD

#%%
# System
sys = System('Anaerobic_Digestion', path=(U1,))                 # aggregation of sanunits
sys.set_dynamic_tracker(eff, gas, U1)                           # what you want to track changes in concentration
# sys                                                           # before running the simulation, 'outs' have nothing

#%%
# Simulation settings
t = 40                        # total time for simulation
t_step = 1                    # times at which to store the computed solution

method = 'BDF'                  # integration method to use
# method = 'RK45'
# method = 'RK23'
# method = 'DOP853'
# method = 'Radau'
# method = 'LSODA'

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html

# Run simulation
sys.simulate(state_reset_hook='reset_cache',
             t_span=(0,t),
             t_eval=np.arange(0, t+t_step, t_step),
             method=method,
             export_state_to=f'sol_{t}d_{method}_AD.xlsx',               # export simulation result as excel file
            )
#
# sys                                                                      # now you have 'outs' info.

#%%
eff.scope.plot_time_series(('S_aa', 'S_fa', 'S_la', 'S_et', 'S_va', 'S_bu', 'S_pro', 'S_ac'))  # you can plot how each state variable changes over time

eff.scope.plot_time_series(('S_su', 'S_et'))

eff.scope.plot_time_series(('S_IC'))

eff.scope.plot_time_series(('X_aa', 'X_fa', 'X_la', 'X_et', 'X_c4', 'X_pro', 'X_ac', 'X_h2'))

gas.scope.plot_time_series(('S_h2'))

gas.scope.plot_time_series(('S_h2','S_ch4', 'S_IC'))
#!!! Soluble biogas could be defined as biogas? Gas measurement in real experiments by Gas Chromatography, it could be partial pressure.

#%%
# VFA change over days
# Total VFAs = 'S_va' + 'S_bu' + 'S_pro' + 'S_ac'    (you can change the equations based on your assumption)
idx_vfa = cmps.indices(['S_va', 'S_la', 'S_bu', 'S_pro', 'S_ac']) #S_la as vfa
idx_h2 = cmps.indices(['S_h2'])
idx_ch4 = cmps.indices(['S_ch4'])
idx_co2 = cmps.indices(['S_IC'])
#idx_pH = cmps.indices(root.data['pH'])
t_stamp = eff.scope.time_series

vfa = eff.scope.record[:,idx_vfa]
total_vfa = np.sum(vfa, axis=1)

plt.figure(figsize=(10, 6))

plt.plot(t_stamp, total_vfa)
plt.xlabel("Time [day]")
plt.ylabel("Total VFA [mg/l]")

# max_vfa = np.max(total_vfa)
# print(f"Maximum total VFA during the simulation is: {max_vfa:.2f} mg/L")

#%%
#!!! pH change over days (right codes below?, if pH_ctrl has fixed value, below not work)
t = 40  # total simulation time in days
t_step = 1  # simulation time step in days
time_stamps = np.arange(0, t + t_step, t_step)

# Prepare a list to hold pH values over time
pH_values = []

# Loop over the simulation period, simulate, and collect pH values
for current_time in time_stamps:
    # Simulate for the current time step
    sys.simulate(state_reset_hook='reset_cache', t_span=(current_time, current_time + t_step), t_eval=[current_time])
    
    # Access and store the current pH value from the root data structure
    current_pH = sys.path[0].model.rate_function._params['root'].data['pH']  # sys.path[0] should be your UASB reactor
    pH_values.append(current_pH)

# Plotting the pH values over time
plt.figure(figsize=(10, 6))
plt.plot(time_stamps, pH_values, label='pH Level', marker='o', linestyle='-')
plt.xlabel("Time [days]")
plt.ylabel("pH Level")
plt.title("pH Levels Over Time")
plt.legend()
plt.grid(True)
plt.show()
#%%
#!!! Partial Pressure of gas over days (S_ch4 vs Partial Pressure of CH4)
# Simulation settings
t = 40  # total simulation time in days
t_step = 1  # simulation time step in days
time_stamps = np.arange(0, t + t_step, t_step)

# Prepare lists to hold partial pressures of H2, CH4, and IC over time
pH2_values = []
pCH4_values = []
pIC_values = []

# Loop over the simulation period, simulate, and collect partial pressures
for current_time in time_stamps:
    # Simulate for the current time step
    sys.simulate(state_reset_hook='reset_cache', t_span=(current_time, current_time + t_step), t_eval=[current_time])
    
    # Access and store the current partial pressures from the root data structure
    current_pH2 = sys.path[0].model.rate_function._params['root'].data['biogas_p_h2']
    current_pCH4 = sys.path[0].model.rate_function._params['root'].data['biogas_p_ch4']
    current_pIC = sys.path[0].model.rate_function._params['root'].data['biogas_p_IC']
    
    pH2_values.append(current_pH2)
    pCH4_values.append(current_pCH4)
    pIC_values.append(current_pIC)

# Plotting the partial pressures over time
plt.figure(figsize=(12, 8))
plt.plot(time_stamps, pH2_values, label='H2 Partial Pressure', marker='o', linestyle='-', color='blue')
plt.plot(time_stamps, pCH4_values, label='CH4 Partial Pressure', marker='x', linestyle='-', color='green')
plt.plot(time_stamps, pIC_values, label='IC Partial Pressure', marker='^', linestyle='-', color='red')

plt.xlabel("Time [days]")
plt.ylabel("Partial Pressure[Pa]")
plt.title("Partial Pressures Over Time")
plt.legend()
plt.grid(True)
plt.show()
#%%
# if pH levels varies
pH_levels = [4, 5, 6, 7, 8]

# Dictionaries to store VFA values and time stamps for each pH
VFA_values_for_each_pH = {}
t_stamps_for_each_pH = {}

h2_values_for_each_pH = {}
ch4_values_for_each_pH = {}
co2_values_for_each_pH = {}

# Simulate the system at each pH
for pH in pH_levels:
    # Reset the reactor to its initial conditions
    U1.set_init_conc(**default_init_conds)

    U1.pH_ctrl = pH
    sys.simulate(state_reset_hook='reset_cache',
                  t_span=(0, t),
                  t_eval=np.arange(0, t + t_step, t_step),
                  method=method,
                )

    # Calculate VFA values after this simulation
    vfa = eff.scope.record[:, idx_vfa]
    total_vfa = np.sum(vfa, axis=1)

    # Store the VFA values and time stamps in the dictionaries
    VFA_values_for_each_pH[pH] = total_vfa
    t_stamps_for_each_pH[pH] = eff.scope.time_series

    h2_values_for_each_pH[pH] = gas.scope.record[:, idx_h2]
    ch4_values_for_each_pH[pH] = gas.scope.record[:, idx_ch4]
    co2_values_for_each_pH[pH] = gas.scope.record[:, idx_co2]


# Plot the results for each pH on the same graph
# vfa
plt.figure(figsize=(10, 6))

for pH, vfa in VFA_values_for_each_pH.items():
    plt.plot(t_stamps_for_each_pH[pH], vfa, label=f"pH = {pH}")

plt.xlabel("Time [day]")
plt.ylabel("Total VFA [mg/l]")
plt.title("Total VFA production at varying pH levels")
plt.legend()
plt.grid(True)
plt.show()

# h2
plt.figure(figsize=(10, 6))

for pH, h2 in h2_values_for_each_pH.items():
    plt.plot(t_stamps_for_each_pH[pH], h2, label=f"pH = {pH}")

plt.xlabel("Time [day]")
plt.ylabel("H2 [mg/L]")
plt.legend()
plt.grid(True)
plt.show()

# ch4
plt.figure(figsize=(10, 6))

for pH, ch4 in ch4_values_for_each_pH.items():
    plt.plot(t_stamps_for_each_pH[pH], ch4, label=f"pH = {pH}")

plt.xlabel("Time [day]")
plt.ylabel("CH4 [mg/L]")
plt.legend()
plt.grid(True)
plt.show()

# co2
plt.figure(figsize=(10, 6))

for pH, co2 in co2_values_for_each_pH.items():
    plt.plot(t_stamps_for_each_pH[pH], co2, label=f"pH = {pH}")

plt.xlabel("Time [day]")
plt.ylabel("CO2 [mg/L]")
plt.legend()
plt.grid(True)
plt.show()
#%%

# Extracting indices for each component
idx_s_aa = cmps.indices(['S_aa'])[0]
idx_s_fa = cmps.indices(['S_fa'])[0]
idx_s_la = cmps.indices(['S_la'])[0]
idx_s_et = cmps.indices(['S_et'])[0]
idx_s_va = cmps.indices(['S_va'])[0]
idx_s_bu = cmps.indices(['S_bu'])[0]
idx_s_pro = cmps.indices(['S_pro'])[0]
idx_s_ac = cmps.indices(['S_ac'])[0]

# Extracting concentration data for each component
data_s_aa = eff.scope.record[:, idx_s_aa]
data_s_fa = eff.scope.record[:, idx_s_fa]
data_s_la = eff.scope.record[:, idx_s_la]
data_s_et = eff.scope.record[:, idx_s_et]
data_s_va = eff.scope.record[:, idx_s_va]
data_s_bu = eff.scope.record[:, idx_s_bu]
data_s_pro = eff.scope.record[:, idx_s_pro]
data_s_ac = eff.scope.record[:, idx_s_ac]

# Time stamps
time_stamps = eff.scope.time_series

# Stacking data for cumulative plotting
cumulative_data = np.vstack((data_s_aa, data_s_fa, data_s_la, data_s_et, data_s_va, data_s_bu, data_s_pro, data_s_ac)).T
cumulative_sums = np.cumsum(cumulative_data, axis=1)

# Plotting
plt.figure(figsize=(12, 8))
components = ['S_aa', 'S_fa', 'S_la', 'S_et', 'S_va', 'S_bu', 'S_pro', 'S_ac']
colors = plt.cm.tab10(np.linspace(0, 1, len(components)))  # Optional: specify a colormap for the bars

for i, (component, color) in enumerate(zip(components, colors)):
    plt.bar(time_stamps, cumulative_sums[:, i] if i == 0 else cumulative_sums[:, i] - cumulative_sums[:, i-1],
            bottom=None if i == 0 else cumulative_sums[:, i-1],
            label=component, color=color)

plt.xlabel('Time [day]')
plt.ylabel('Concentration [mg/L]')
plt.title('Cumulative Concentration of Components Over Time')
plt.legend()
plt.show()
