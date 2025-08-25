#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

#%% initialization

import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.mathtext import _mathtext as mathtext
from colorpalette import Color
from scipy.stats import qmc
from warnings import filterwarnings
from exposan.htl import (create_C1_system, create_C2_system, create_C3_system,
                         create_C4_system, create_C5_system, create_C6_system,
                         create_C7_system, create_C8_system, create_C9_system,
                         create_C10_system, create_C11_system, create_C12_system,
                         create_C13_system, create_C14_system, create_C15_system,
                         create_C16_system, create_C17_system, create_C18_system,
                         create_C19_system, create_C20_system, create_C21_system,
                         create_C22_system, create_C23_system, create_C24_system,
                         create_C25_system, create_T1_system, create_T2_system,
                         create_T3_system, create_T4_system, create_T5_system,
                         create_T6_system, create_T7_system, create_T8_system,
                         create_T9_system, create_T10_system, create_T11_system,
                         create_T12_system, create_T13_system, create_T14_system,
                         create_T15_system)

# color palette
b = Color('blue', (96, 193, 207)).HEX
g = Color('green', (121, 191, 130)).HEX
r = Color('red', (237, 88, 111)).HEX
o = Color('orange', (249, 143, 96)).HEX
y = Color('yellow', (243, 195, 84)).HEX
a = Color('gray', (144, 145, 142)).HEX
p = Color('purple', (162, 128, 185)).HEX

db = Color('dark_blue', (53, 118, 127)).HEX
dg = Color('dark_green', (77, 126, 83)).HEX
dr = Color('dark_red', (156, 75, 80)).HEX
do = Color('dark_orange', (167, 95, 62)).HEX
dy = Color('dark_yellow', (171, 137, 55)).HEX
da = Color('dark_gray', (78, 78, 78)).HEX
dp = Color('dark_purple', (76, 56, 90)).HEX

#%% WWTP random sampling

# TODO: adjust typology if including AD and AeD (any any units before them)
# all use uniform distribution (not Monte Carlo)
# for now, assume WRRFs are only responsible for the transportation of e.g., solids (and ash) to landfills, biocrude to oil refineries
# for other products, assume they can be sold at the gate of WRRFs
typology_bounds = {
    # tonne/day
    'solids_mass_flow': (0, 250),
    # dry weight ratio
    'solids_dw_ash': (0.15, 0.45),
    # ash free dry weight ratio
    'solids_afdw_lipid': (0.05, 0.35),
    # ash free dry weight ratio
    'solids_afdw_protein': (0.3, 0.55),
    # $/kWh
    'electricity_cost': (0.05, 0.2),
    # kg CO2 eq/kWh
    'electricity_CI': (0, 1.5),
    # km
    'landfill_distance': (0, 1500),
    # km
    'biofuel_distance': (0, 1500)
    }

typology_keys = list(typology_bounds.keys())

# Latin Hypercube Sampling
sampler = qmc.LatinHypercube(d=len(typology_keys), seed=3221)
# TODO: decide the sample size based on the computation requirement
LHS_unit = sampler.random(n=10000)

typology_mins = np.array([typology_bounds[key][0] for key in typology_keys])
typology_maxs = np.array([typology_bounds[key][1] for key in typology_keys])
LHS_samples = qmc.scale(LHS_unit, typology_mins, typology_maxs)

typology = {key: LHS_samples[:, i] for i, key in enumerate(typology_keys)}

lipid = typology['solids_afdw_lipid']
protein = typology['solids_afdw_protein']
carbohydrate = 1 - lipid - protein
# no negative carb
mask = carbohydrate >= 0
typology = {k: v[mask] for k, v in typology.items()}
typology['solids_afdw_carbohydrate'] = carbohydrate[mask]

for i in range(len(LHS_samples)):
    solids_mass_flow = typology['solids_mass_flow'][i]
    solids_dw_ash = typology['solids_dw_ash'][i]
    solids_afdw_lipid = typology['solids_afdw_lipid'][i]
    solids_afdw_protein = typology['solids_afdw_protein'][i]
    solids_afdw_carbohydrate = typology['solids_afdw_carbohydrate'][i]
    electricity_cost = typology['electricity_cost'][i]
    electricity_CI = typology['electricity_CI'][i]
    landfill_distance = typology['landfill_distance'][i]
    biofuel_distance = typology['biofuel_distance'][i]
    
    # TODO: compare the baseline cost and GHG results from different systems
    # TODO: for each typology, identify a winner
    # TODO: plot 5th and 95th percentiles to represent the opportunity space

#%% preliminary utility check

results = []
for function in (create_C1_system, create_C2_system, create_C3_system,
                 create_C4_system, create_C5_system, create_C6_system,
                 create_C7_system, create_C8_system, create_C9_system,
                 create_C10_system, create_C11_system, create_C12_system,
                 create_C13_system, create_C14_system, create_C15_system,
                 create_C16_system, create_C17_system, create_C18_system,
                 create_C19_system, create_C20_system, create_C21_system,
                 create_C22_system, create_C23_system, create_C24_system,
                 create_C25_system, create_T1_system, create_T2_system,
                 create_T3_system, create_T4_system, create_T5_system,
                 create_T6_system, create_T7_system, create_T8_system,
                 create_T9_system, create_T10_system, create_T11_system,
                 create_T12_system, create_T13_system, create_T14_system,
                 create_T15_system):
    sys = function(size=10)
    print('\n' + sys.ID)
    print(round(sys.get_cooling_duty()))
    print(round(sys.get_heating_duty()))
    print(sys.heat_utilities)

#%% preliminary TEA check

TEA_results = []
for function in (create_C1_system, create_C2_system, create_C3_system,
                 create_C4_system, create_C5_system, create_C6_system,
                 create_C7_system, create_C8_system, create_C9_system,
                 create_C10_system, create_C11_system, create_C12_system,
                 create_C13_system, create_C14_system, create_C15_system,
                 create_C16_system, create_C17_system, create_C18_system,
                 create_C19_system, create_C20_system, create_C21_system,
                 create_C22_system, create_C23_system, create_C24_system,
                 create_C25_system, create_T1_system, create_T2_system,
                 create_T3_system, create_T4_system, create_T5_system,
                 create_T6_system, create_T7_system, create_T8_system,
                 create_T9_system, create_T10_system, create_T11_system,
                 create_T12_system, create_T13_system, create_T14_system,
                 create_T15_system):
    sys = function(size=10)
    
    print('\n' + sys.ID)
    
    print(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
    
    TEA_results.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
    
print(np.quantile(TEA_results[0:25], 0.05))
print(np.quantile(TEA_results[0:25], 0.5))
print(np.quantile(TEA_results[0:25], 0.95))

print(np.quantile(TEA_results[25:40], 0.05))
print(np.quantile(TEA_results[25:40], 0.5))
print(np.quantile(TEA_results[25:40], 0.95))

plt.boxplot((TEA_results[0:25], TEA_results[25:40]))

#%% preliminary LCA check

LCA_results = []
for function in (create_C1_system, create_C2_system, create_C3_system,
                 create_C4_system, create_C5_system, create_C6_system,
                 create_C7_system, create_C8_system, create_C9_system,
                 create_C10_system, create_C11_system, create_C12_system,
                 create_C13_system, create_C14_system, create_C15_system,
                 create_C16_system, create_C17_system, create_C18_system,
                 create_C19_system, create_C20_system, create_C21_system,
                 create_C22_system, create_C23_system, create_C24_system,
                 create_C25_system, create_T1_system, create_T2_system,
                 create_T3_system, create_T4_system, create_T5_system,
                 create_T6_system, create_T7_system, create_T8_system,
                 create_T9_system, create_T10_system, create_T11_system,
                 create_T12_system, create_T13_system, create_T14_system,
                 create_T15_system):
    sys = function(size=10)
    
    print('\n' + sys.ID)
    
    print(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
    LCA_results.append(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)
    
print(np.quantile(LCA_results[0:25], 0.05))
print(np.quantile(LCA_results[0:25], 0.5))
print(np.quantile(LCA_results[0:25], 0.95))

print(np.quantile(LCA_results[25:40], 0.05))
print(np.quantile(LCA_results[25:40], 0.5))
print(np.quantile(LCA_results[25:40], 0.95))

plt.boxplot((LCA_results[0:25], LCA_results[25:40]))

#%% preliminary TEA and LCA animation

filterwarnings('ignore')

system_ID_list = []
system_type_list = []
size_list = []
TEA_result_list = []
LCA_result_list = []

for function in (create_C1_system, create_C2_system, create_C3_system,
                 create_C4_system, create_C5_system, create_C6_system,
                 create_C7_system, create_C8_system, create_C9_system,
                 create_C10_system, create_C11_system, create_C12_system,
                 create_C13_system, create_C14_system, create_C15_system,
                 create_C16_system, create_C17_system, create_C18_system,
                 create_C19_system, create_C20_system, create_C21_system,
                 create_C22_system, create_C23_system, create_C24_system,
                 create_C25_system, create_T1_system, create_T2_system,
                 create_T3_system, create_T4_system, create_T5_system,
                 create_T6_system, create_T7_system, create_T8_system,
                 create_T9_system, create_T10_system, create_T11_system,
                 create_T12_system, create_T13_system, create_T14_system,
                 create_T15_system):
    for size in np.linspace(1, 50, 50):
        sys = function(size=size)
        
        if size == 1:
            print(sys.ID)
        
        system_ID_list.append(sys.ID)
        system_type_list.append(sys.ID[7])
        size_list.append(size)
        TEA_result_list.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
        LCA_result_list.append(sys.LCA.get_total_impacts(operation_only=True,
                                                     exclude=(sys.flowsheet.raw_wastewater,),
                                                     annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)

animation_data = pd.DataFrame({'system_ID': system_ID_list,
                               'system_type': system_type_list,
                               'size': size_list,
                               'cost': TEA_result_list,
                               'CI': LCA_result_list})

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 38
plt.rcParams['ytick.labelsize'] = 38

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

fig, ax = plt.subplots(figsize=(14, 12))

ax = plt.gca()
ax.set_xlim(0, 4000)
ax.set_ylim(-1000, 2500)

ax.tick_params(direction='inout', length=20, width=3,
                bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 4500, 500), fontname='Arial')
plt.yticks(np.arange(-1000, 3000, 500), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
plt.xticks(np.arange(0, 4500, 500), fontname='Arial')

ax_top.tick_params(direction='in', length=10, width=3,
                   bottom=False, top=True, left=False, right=False,
                   labelcolor='none')

ax_bottom = ax.twinx()
ax_bottom.set_ylim(ax.get_ylim())
plt.xticks(np.arange(0, 4500, 500), fontname='Arial')

ax_bottom.tick_params(direction='in', length=10, width=3,
                      bottom=False, top=False, left=False, right=True,
                      labelcolor='none')

mathtext.FontConstantsBase.sup1 = 0.35

ax.set_xlabel(r'$\mathbf{Cost}$ ' + '[\$·day${^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{CI}$ ' + '[tonne CO${_2}$ eq·day${^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

scatter_C = ax.scatter(x=[], y=[], s=500, color=a, linewidths=2, edgecolors='k')
scatter_T = ax.scatter(x=[], y=[], s=500, color=r, linewidths=2, edgecolors='k')

flow_scenarios = animation_data['size'].unique()

def update(frame):
    fr = flow_scenarios[frame]
    animation_data_fr = animation_data[animation_data['size'] == fr]
    danimation_data_C = animation_data_fr[animation_data_fr['system_ID'].str.contains('C')]
    animation_data_T = animation_data_fr[animation_data_fr['system_ID'].str.contains('T')]

    scatter_C.set_offsets(list(zip(danimation_data_C['cost'], danimation_data_C['CI'])))
    scatter_T.set_offsets(list(zip(animation_data_T['cost'], animation_data_T['CI'])))
    
    ax.set_title(f'{int(fr)} ' + 'dry tonne·day${^{-1}}$', fontname='Arial', fontsize=45)
    return scatter_C, scatter_T

ani = FuncAnimation(fig, update, frames=len(flow_scenarios), interval=800, blit=True)

fig.tight_layout()

ani.save('preliminary_TEA_LCA_results.mp4', writer="ffmpeg", fps=2, savefig_kwargs={'bbox_inches':'tight'})

plt.show()