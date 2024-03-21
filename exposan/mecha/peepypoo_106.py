# vim: set tabstop=4
# s_competition.py
#!/usr/bin/env python3

""" This is a implementation of the Tilburg water line using the python package
peepypoo and can be used for the competition on wastewater modelling with
mechanistic, data-driven and hybrid models:
https://www.kaggle.com/competitions/dynamic-modeling-of-wastewater-treatment-process/overview
This script builds on the usage example of BSM1 in peepypoo:
https://datinfo.gitlab.io/PeePyPoo/api.html#bsm1 """

# Copyright (C) 2023 Mariane Yvonne Schneider
# Copyright (C) 2023 Saba Daneshgar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# Mariane Yvonne Schneider <mariane.schneider@ugent.be>
# Saba Daneshgar <saba.daneshgar@ugent.be>
# 19.09.2023

############
## Imports
# built-ins

# third party
import numpy as np
import pandas as pd
import os
from exposan.mecha import data_path

# user
import peepypoo # how to install: https://datinfo.gitlab.io/PeePyPoo/installation.html


# %%
# Stoichiometric parameters at 20 degrees celsius
f_P = 0.08  # fraction of biomass leading to particulate material (default 0.08 -)
i_XB = 0.08 # nitrogen fraction in biomass (default 0.086 gN/gCOD)
i_XP = 0.06 # nitrogen fraction in endogenous mass (default 0.01 gN/gCOD)

# Kinetic parameters
mu_H = 6.0  # maximum specific growth rate (default 6.0 1/day)
K_S = 20.0  # substrate saturation constant (default 20.0 gCOD/m3)
K_OH = 0.2  # Oxygen saturation constant (default 0.2 gO2/m3)
b_H = 0.62  # specific decay rate (default 0.62 1/day)
eta_g = 0.8 # anoxic growth correction factor (default 0.8 -)
eta_h = 0.4 # anoxic hydrolysis correction factor (default 0.4 -)
k_h = 3.0   # maximum specific hydrolysis rate (default 3.0 1/day)
K_X = 0.03  # half-saturation coefficient for hydrolysis of XS (default 0.03 -)

b_A = 0.1   # specific decay rate (default 0.1 1/day)
K_OA = 0.4  # oxygen saturation constant (default 0.4 gO2/m3)
k_a = 0.08  # ammonification rate constant (default 0.08 m3/(gCOD*day))

# increase
Y_A = 0.6  # autotrophic yield (default 0.24 gCOD/gCOD)
# Y_A = 0.24  # autotrophic yield (default 0.24 gCOD/gCOD)


# decrease
Y_H = 0.67  # heterotrophic yield (default 0.67 gCOD/gCOD)
mu_A = 0.8  # maximum specific growth rate (default 0.8 1/day)
# mu_A = 0.8  # maximum specific growth rate (default 0.8 1/day)
K_NO = 0.5  # nitrate saturation constant (default 0.5 gNO3-N/m3)

K_NH = 0.8  # ammonium saturation constant (default 1.0 gO2/m3)
# K_NH = 1.0  # ammonium saturation constant (default 1.0 gO2/m3)
#X_BA

# Volume all in m3
V_anox = 3*3640            # volume of the three paralell anoxic streets
V_facultative = 2820+2*3640 # volume of the three paralell facultative streets
V_aer1 = 11179              # sum of aerobic tanks 1/3
V_aer2 = 11179              # sum of aerobic tanks 2/3
V_aer3 = 11179              # sum of aerobic tanks 3/3

#%%
# the average of 6 months traingin data is 2281 m3/h all flows need to be of the same type.
# Therefore they need to be floats in this skrip.

# # static inflow
# inf_stab_data = pd.DataFrame([[0, 30, 69.50, 51.2, 202.32, 28.17, 0, 0, 0, 0, 31.56, 6.95, 10.59, 7.0, 24*2281.0],
#                               [100, 30, 69.50, 51.2, 202.32, 28.17, 0, 0, 0, 0, 31.56, 6.95, 10.59, 7.0, 18446.0]],
#                              columns=['t', 'Si', 'Ss', 'Xi', 'Xs', 'Xbh', 'Xba', 'Xp',
#                                       'So', 'Sno', 'Snh', 'Snd', 'Xnd', 'Salk', 'Q'])
#
# # here rain data could be added
# # inf_rain_data = pd.read_csv()
# # inf_rain_data.t += 100
# # inf_data = pd.concat([inf_stab_data, inf_rain_data])
# inf_data = inf_stab_data
#
# # Inflow, split into flow rate and concentrations
# inflow = peepypoo.Systems.Static.PiecewiseLinearInterpolation(
#     [inf_data.t.to_numpy()],  # Time vector
#     [inf_data.Q.to_numpy(), inf_data.drop(['t', 'Q'], axis=1).to_numpy().T], # Input vector as [Q, C]
#     name="Inflow")

# adapt this path depending if you are in the main folder or the
print("warning: please check if the input file has the same variables and the same order as in the ASM1 module of PeePyPoo defined")
print("it is: | S_I | S_S  | X_I | X_S | X_BH| X_BA| X_P | S_O  |S_NO  |S_NH  |S_ND  | X_ND| S_ALK|")
inf_dyn = pd.read_csv(os.path.join(data_path, 'peepypoo_training_dynamic_influent_kaggle.csv'), sep='\t', skiprows=[1])
# inf_dyn = pd.read_csv('data/peepypoo_training_dynamic_influent.csv', sep='\t', skiprows=[1])

#%%
# Select the inflow

inflow = peepypoo.Systems.Static.PiecewiseLinearInterpolation(t=[inf_dyn['Time'].values,
                                                              inf_dyn['Time'].values],
                                                              y=[inf_dyn['Q_in'].values,
                                                              inf_dyn[[c for c in inf_dyn.columns if c not in ('Time', 'Q_in', 'Unnamed: 0')]].values.T])

#%%
# Constant parameters for the recirculations and the clarifier sludge output
# Clarifier underflow rate (Qu = Qr + Qw)
Qu_rate = peepypoo.Systems.Static.Constant([24*1630.5], name="Qu_rate")
# Wastage Flow rate Qw
Qw_rate = peepypoo.Systems.Static.Constant([24*30.5], name="Qw_rate")
# Internal recycle flow rate Qa
Qa_rate = peepypoo.Systems.Static.Constant([24*9800.0], name="Qa_rate")

# Flow elements
# Flow unifier, unites inflow, internal and external recycle
unifier = peepypoo.Systems.Static.FlowElements.FlowUnifier(
    num_materials=13,
    num_flows=3,
    name='Unifier'
)
# Flow separator for the wastage pump
separator = peepypoo.Systems.Static.FlowElements.FlowSeparator(
    num_materials=13,
    num_flows=2,
    use_percentages=False,
    name='Wastage Pump'
)
# Flow separator for the internal recycle
internal_rec = peepypoo.Systems.Static.FlowElements.FlowSeparator(
    num_materials=13,
    num_flows=2,
    use_percentages=False,
    name='Internal Recycle'
)

# Clarifier
clarifier = peepypoo.Systems.Dynamic.Continuous.Settlers.SettlerTakacs(
    num_materials=13,
    #               | S_I | S_S  | X_I | X_S | X_BH| X_BA| X_P | S_O  |S_NO  |S_NH  |S_ND  | X_ND| S_ALK|
    solid_materials=[False, False, True, True, True, True, True, False, False, False, False, True, False],
    provided_outflow='sludge',
    num_layers=10,
    Xf=(lambda x: 0.75*sum(x[2:7])),  # Formula to calculate TSS from the states according to ASM1
    feedlayer=5,  # Stating to count at 0, thus this is the 6th layer
    A=1700, z=2.3, v0_max=250.0, v0=474.0, rh=0.00019, rp=0.0028, fns=0.001, Xt=3000.0,
    initial_state=np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]*10),
    name='Settler'
)

# System Parameters
params = peepypoo.Systems.Static.Constant(
    [np.array([Y_H, Y_A, f_P, i_XB, i_XP, mu_H, K_S, K_OH, K_NO, b_H, mu_A, K_NH, K_OA, b_A, eta_g, k_a, k_h, K_X, eta_h])]
)

# Reactors
# This is the anoxic tank that receives inflow, return sludge and internal recirculation. The flows and
#
r_anox1 = peepypoo.Systems.Dynamic.Continuous.Reactors.CSTR(
        v=V_anox,
        initial_state=np.array([30, 69.5, 1000, 202.32, 1000, 500, 100, 0.0, 1.0, 0.23, 6.95, 10.59, 7.0]),
        # initial_state=np.array([30, 69.5, 51.2, 202.32, 28.17, 1.0, 1.0, 0.0, 1.0, 31.56, 6.95, 10.59, 7.0]),
        stoichiometry_file='ASM1',
        name='Anoxic tank'
)

r_facultative = peepypoo.Systems.Dynamic.Continuous.Reactors.CSTR(
        v=V_facultative,
        initial_state=np.array([30, 69.5, 1000, 202.32, 1000, 500, 100, 0.2, 1.0, 0.23, 6.95, 10.59, 7.0]),
        # initial_state=np.array([30, 69.5, 51.2, 202.32, 28.17, 1.0, 1.0, 0.2, 1.0, 31.56, 6.95, 10.59, 7.0]),
        stoichiometry_file='ASM1',
        name='Facultative tank'
)

r_aer1 = peepypoo.Systems.Dynamic.Continuous.Reactors.CSTR(
        v=V_aer1,
        initial_state=np.array([30, 69.5, 1000, 202.32, 1000, 500, 100, 0.46, 1.0, 0.23, 6.95, 10.59, 7.0]),
        # initial_state=np.array([30, 69.5, 51.2, 202.32, 28.17, 1, 1, 0.46, 1, 31.56, 6.95, 10.59, 7.0]),
        stoichiometry_file='ASM1_constO',
        name='Aerobic tank 1'
)

r_aer2 = peepypoo.Systems.Dynamic.Continuous.Reactors.CSTR(
        v=V_aer2,
        initial_state=np.array([30, 69.5, 1000, 202.32, 1000, 500, 100, 0.29, 1.0, 0.23, 6.95, 10.59, 7.0]),
        # initial_state=np.array([30, 69.5, 51.2, 202.32, 28.17, 1.0, 1.0, 0.29, 1.0, 31.56, 6.95, 10.59, 7.0]),
        stoichiometry_file='ASM1_constO',
        name='Aerobic tank 2'
)

r_aer3 = peepypoo.Systems.Dynamic.Continuous.Reactors.CSTR(
        v=V_aer3,
        initial_state=np.array([30, 69.5, 1000, 202.32, 1000, 500, 100, 0.67, 1.0, 0.23, 6.95, 10.59, 7.0]),
        # initial_state=np.array([30, 69.5, 51.2, 202.32, 28.17, 1.0, 1.0, 0.67, 1.0, 31.56, 6.95, 10.59, 7.0]),
        stoichiometry_file='ASM1_constO',
        name='Aerobic tank 3'
)

# Connect the params to the reactors
r_anox1.add_input_connection(params, [2], [0])
r_facultative.add_input_connection(params, [2], [0])
r_aer1.add_input_connection(params, [2], [0])
r_aer2.add_input_connection(params, [2], [0])
r_aer3.add_input_connection(params, [2], [0])

# Sludge feedback rate
Qu_rate.add_output_connection(clarifier, [2], [0])
# Excess sludge rate
Qw_rate.add_output_connection(separator, [2], [0])
# internal recirculation
Qa_rate.add_output_connection(internal_rec, [2], [0])


# Connect the forward flows
inflow.add_output_connection(unifier, [0, 1], [0, 1])
unifier.add_output_connection(r_anox1, [0, 1], [0, 1])   # (system, ranox(flow, concentration-port), unifier)
r_anox1.add_output_connection(r_facultative, [0, 1], [0, 1])
r_facultative.add_output_connection(r_aer1, [0, 1], [0, 1])
r_aer1.add_output_connection(r_aer2, [0, 1], [0, 1])
r_aer2.add_output_connection(r_aer3, [0, 1], [0, 1])
r_aer3.add_output_connection(internal_rec, [0, 1], [0, 1])
internal_rec.add_output_connection(clarifier, [0, 1], [0, 1])

# connect the recycling flows
clarifier.add_output_connection(separator, [0, 1], [2, 3])
separator.add_output_connection(unifier, [2, 3], [0, 1])
internal_rec.add_output_connection(unifier, [4, 5], [2, 3])

# evaluation
wrapper = peepypoo.BDWrapper([inflow]) # adds the whole connected tree gets added automatically to the wrapper
# Simulate the water line of the plant in julia and returns a pandas dataframe
out = wrapper.evaluate_blockdiagram_julia(0, 270, np.arange(0, 270, 1/96))
out.to_excel('results/peepypoo_result_kaggle.xlsx')

# %%
# import matplotlib.pyplot as plt
# q_in = []
# for row in out[inflow]:
#     q_in.append(row[1][7])

# # (output from 'clarifier.outputs'):
# # e is effluent, s is sludge (underflow clarifier)
# # (q_out_e, (C_out_e0, C_out_e1, C_out_e2, C_out_e3, C_out_e4, C_out_e5,
# # C_out_e6, C_out_e7, C_out_e8, C_out_e9, C_out_e10, C_out_e11, C_out_e12),
# # q_out_s, (C_out_s0, C_out_s1, C_out_s2, C_out_s3, C_out_s4, C_out_s5,
# # C_out_s6, C_out_s7, C_out_s8, C_out_s9, C_out_s10, C_out_s11, C_out_s12))

# q_out = []
# for row in out[clarifier]:
#     q_out.append(row[1][7])


# plt.plot(q_in)
# plt.plot(q_out)
# plt.show()
#%%

result=pd.DataFrame()
result['t_stamp']=out.index

for i in range(25920):
    index = out[r_aer3].index[i]
    result.loc[i,'S_NO'] = out[r_aer3][index][1][8]   # S_NO
    result.loc[i,'S_NH'] = out[r_aer3][index][1][9]   # S_NH

result.to_excel('results/peepypoo_result_Aer3_NO_NH_only.xlsx')

answer = pd.read_excel('data/train_val_test_online_dataset_.xlsx', sheet_name=0, index_col='DATETIME')

from datetime import datetime
day0 = datetime.strptime('5-5-2021 11:57 AM', '%m-%d-%Y %I:%M %p') # set the first time stamp as 'day0'
calc_day = lambda t_stamp: (t_stamp-day0).total_seconds()/60/60/24

t_stamp = answer.index
t_intp = calc_day(t_stamp).to_numpy()

import matplotlib.pyplot as plt, matplotlib.ticker as ticker

def plot_raw_vs_smooth(x1, y1, x2, y2):
    fig, ax = plt.subplots(figsize=(15, 3)) # figure size
    l1, = ax.plot(x1, y1, label='Model') # first plot = raw data
    l2, = ax.plot(x2, y2, label='Answer') # second plot = smoothed result
    ax.xaxis.set_major_locator(ticker.MultipleLocator(30)) # x-axis tick interval
    ax.set_xlabel('Time [day]') # x-axis title
    ax.legend(handles=[l1,l2])
    return fig, ax

plot_raw_vs_smooth(t_intp, answer['NH4'], result['t_stamp'], result['S_NH'])
plot_raw_vs_smooth(t_intp, answer['NO3'], result['t_stamp'], result['S_NO'])

#%%

answer.drop(answer.tail(1).index, inplace=True)
answer.index=result.index

eq_predicted = 30 * result['S_NH'] + 10 * result['S_NO']
eq_measured = 30 * answer['NH4'] + 10 * answer['NO3']

rmse = (sum((eq_predicted - eq_measured)**2)/len(answer.index))**0.5

print(rmse)
