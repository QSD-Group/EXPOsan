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

import numpy as np, matplotlib.pyplot as plt
from colorpalette import Color
from matplotlib.mathtext import _mathtext as mathtext
from matplotlib.colors import to_hex
from matplotlib.gridspec import GridSpec
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

mathtext.FontConstantsBase.sup1 = 0.35

# color palette
b = Color('blue', (96, 193, 207)).HEX
g = Color('green', (121, 191, 130)).HEX
r = Color('red', (237, 88, 111)).HEX
o = Color('orange', (249, 143, 96)).HEX
y = Color('yellow', (243, 195, 84)).HEX
a = Color('gray', (144, 145, 142)).HEX
p = Color('purple', (162, 128, 185)).HEX

lb = to_hex((96/256*0.5 + 1*0.5, 193/256*0.5 + 1*0.5, 207/256*0.5 + 1*0.5))
lg = to_hex((121/256*0.5 + 1*0.5, 191/256*0.5 + 1*0.5, 130/256*0.5 + 1*0.5))
lr = to_hex((237/256*0.5 + 1*0.5, 88/256*0.5 + 1*0.5, 111/256*0.5 + 1*0.5))
lo = to_hex((249/256*0.5 + 1*0.5, 143/256*0.5 + 1*0.5, 96/256*0.5 + 1*0.5))
ly = to_hex((243/256*0.5 + 1*0.5, 195/256*0.5 + 1*0.5, 84/256*0.5 + 1*0.5))
la = to_hex((144/256*0.5 + 1*0.5, 145/256*0.5 + 1*0.5, 142/256*0.5 + 1*0.5))
lp = to_hex((162/256*0.5 + 1*0.5, 128/256*0.5 + 1*0.5, 185/256*0.5 + 1*0.5))

db = Color('dark_blue', (53, 118, 127)).HEX
dg = Color('dark_green', (77, 126, 83)).HEX
dr = Color('dark_red', (156, 75, 80)).HEX
do = Color('dark_orange', (167, 95, 62)).HEX
dy = Color('dark_yellow', (171, 137, 55)).HEX
da = Color('dark_gray', (78, 78, 78)).HEX
dp = Color('dark_purple', (76, 56, 90)).HEX

#%% maturity opportunity

filterwarnings('ignore')

system_ID_list = []
system_type_list = []
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
    sys = function(size=10)
    
    system_ID_list.append(sys.ID)
    system_type_list.append(sys.ID[7])
    TEA_result_list.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
    LCA_result_list.append(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)

for function in (create_T1_system, create_T2_system, create_T3_system,
                 create_T4_system, create_T5_system, create_T6_system,
                 create_T7_system, create_T8_system, create_T9_system,
                 create_T10_system, create_T11_system, create_T12_system,
                 create_T13_system, create_T14_system, create_T15_system):
    sys = function(size=10, FOAK=False)
    
    system_ID_list.append(sys.ID)
    system_type_list.append(sys.ID[7])
    TEA_result_list.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
    LCA_result_list.append(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)

fig = plt.figure(figsize=(13, 12))
gs = GridSpec(2, 2, width_ratios=[4, 1.25], height_ratios=[1.25, 4], wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 1800)
ax.set_ylim(-1000, 2500)

ax.tick_params(direction='inout', length=20, width=3,
                bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 2100, 300), fontname='Arial')
plt.yticks(np.arange(-1000, 3000, 500), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
plt.xticks(np.arange(0, 2100, 300), fontname='Arial')

ax_top.tick_params(direction='in', length=10, width=3,
                   bottom=False, top=True, left=False, right=False,
                   labelcolor='none')

ax_bottom = ax.twinx()
ax_bottom.set_ylim(ax.get_ylim())
plt.xticks(np.arange(0, 2100, 300), fontname='Arial')

ax_bottom.tick_params(direction='in', length=10, width=3,
                      bottom=False, top=False, left=False, right=True,
                      labelcolor='none')

ax.set_xlabel(r'$\mathbf{Cost}$ ' + '[\$·tonne${^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{CI}$ ' + '[kg CO${_2}$ eq·tonne${^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

scatter_C = ax.scatter(x=TEA_result_list[0:25], y=LCA_result_list[0:25], s=1200, color=a, linewidths=3, edgecolors='k')
scatter_T = ax.scatter(x=TEA_result_list[25:40], y=LCA_result_list[25:40], s=1200, color=r, linewidths=3, edgecolors='k')
scatter_TT = ax.scatter(x=TEA_result_list[40:55], y=LCA_result_list[40:55], s=1200, color=y, linewidths=3, edgecolors='k')

for i in range(25, 40):
    ax.annotate('', xy=(TEA_result_list[i+15], LCA_result_list[i+15]), xytext=(TEA_result_list[i], LCA_result_list[i]),
            arrowprops=dict(arrowstyle='-|>', lw=1, color='k'), size=15)

ax_right_boxplot = fig.add_subplot(gs[1, 1], sharey=ax)
ax_right_data = [LCA_result_list[0:25], LCA_result_list[25:40], LCA_result_list[40:55]]
bp = ax_right_boxplot.boxplot(ax_right_data, vert=True, whis=[5, 95], showfliers=False, widths=0.8, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=a, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=r, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=y, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

ax_right_boxplot.axis('off')

ax_top_boxplot = fig.add_subplot(gs[0, 0], sharex=ax)
ax_top_data = [TEA_result_list[40:55], TEA_result_list[25:40], TEA_result_list[0:25]]
bp_top = ax_top_boxplot.boxplot(ax_top_data, vert=False, whis=[5, 95], showfliers=False, widths=0.8, patch_artist=True)

bp_top['boxes'][0].set(color='k', facecolor=y, linewidth=3)
bp_top['boxes'][1].set(color='k', facecolor=r, linewidth=3)
bp_top['boxes'][2].set(color='k', facecolor=a, linewidth=3)

for whisker in bp_top['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp_top['medians']:
    median.set(color='k', linewidth=3)

for cap in bp_top['caps']:
    cap.set(color='k', linewidth=3)

ax_top_boxplot.axis('off')

# plt.savefig('/Users/jiananfeng/Desktop/maturity_opportunity.png', transparent=True, bbox_inches='tight')

#%% deployment opportunity

filterwarnings('ignore')

system_ID_list = []
system_type_list = []
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
    sys = function(size=10)
    
    system_ID_list.append(sys.ID)
    system_type_list.append(sys.ID[7])
    TEA_result_list.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
    LCA_result_list.append(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)

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
    sys = function(size=100)
    
    system_ID_list.append(sys.ID)
    system_type_list.append(sys.ID[7])
    TEA_result_list.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
    LCA_result_list.append(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)

fig = plt.figure(figsize=(13, 12))
gs = GridSpec(2, 2, width_ratios=[4, 1.25], height_ratios=[1.25, 4], wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 1800)
ax.set_ylim(-1000, 2500)

ax.tick_params(direction='inout', length=20, width=3,
                bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 2100, 300), fontname='Arial')
plt.yticks(np.arange(-1000, 3000, 500), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
plt.xticks(np.arange(0, 2100, 300), fontname='Arial')

ax_top.tick_params(direction='in', length=10, width=3,
                   bottom=False, top=True, left=False, right=False,
                   labelcolor='none')

ax_bottom = ax.twinx()
ax_bottom.set_ylim(ax.get_ylim())
plt.xticks(np.arange(0, 2100, 300), fontname='Arial')

ax_bottom.tick_params(direction='in', length=10, width=3,
                      bottom=False, top=False, left=False, right=True,
                      labelcolor='none')

ax.set_xlabel(r'$\mathbf{Cost}$ ' + '[\$·tonne${^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{CI}$ ' + '[kg CO${_2}$ eq·tonne${^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

scatter_C = ax.scatter(x=TEA_result_list[0:25], y=LCA_result_list[0:25], s=1200, color=a, linewidths=3, edgecolors='k')
scatter_T = ax.scatter(x=TEA_result_list[25:40], y=LCA_result_list[25:40], s=1200, color=r, linewidths=3, edgecolors='k')
scatter_CC = ax.scatter(x=TEA_result_list[40:65], y=LCA_result_list[40:65], s=1200, color=la, linewidths=3, edgecolors='k')
scatter_TT = ax.scatter(x=TEA_result_list[65:80], y=LCA_result_list[65:80], s=1200, color=y, linewidths=3, edgecolors='k')

for i in range(0, 40):
    ax.annotate('', xy=(TEA_result_list[i+40], LCA_result_list[i+40]), xytext=(TEA_result_list[i], LCA_result_list[i]),
            arrowprops=dict(arrowstyle='-|>', lw=1, color='k'), size=15)

ax_right_boxplot = fig.add_subplot(gs[1, 1], sharey=ax)
ax_right_data = [LCA_result_list[0:25], LCA_result_list[40:65], LCA_result_list[25:40], LCA_result_list[65:80]]
bp = ax_right_boxplot.boxplot(ax_right_data, vert=True, whis=[5, 95], showfliers=False, widths=0.8, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=a, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=la, linewidth=3)
bp['boxes'][2].set(color='k', facecolor=r, linewidth=3)
bp['boxes'][3].set(color='k', facecolor=y, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

ax_right_boxplot.axis('off')

ax_top_boxplot = fig.add_subplot(gs[0, 0], sharex=ax)
ax_top_data = [TEA_result_list[65:80], TEA_result_list[25:40], TEA_result_list[40:65], TEA_result_list[0:25]]
bp_top = ax_top_boxplot.boxplot(ax_top_data, vert=False, whis=[5, 95], showfliers=False, widths=0.8, patch_artist=True)

bp_top['boxes'][0].set(color='k', facecolor=y, linewidth=3)
bp_top['boxes'][1].set(color='k', facecolor=r, linewidth=3)
bp_top['boxes'][2].set(color='k', facecolor=la, linewidth=3)
bp_top['boxes'][3].set(color='k', facecolor=a, linewidth=3)

for whisker in bp_top['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp_top['medians']:
    median.set(color='k', linewidth=3)

for cap in bp_top['caps']:
    cap.set(color='k', linewidth=3)

ax_top_boxplot.axis('off')

# plt.savefig('/Users/jiananfeng/Desktop/deployment_opportunity.png', transparent=True, bbox_inches='tight')

#%% credit opportunity

filterwarnings('ignore')

system_ID_list = []
system_type_list = []
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
    sys = function(size=10)
    
    system_ID_list.append(sys.ID)
    system_type_list.append(sys.ID[7])
    TEA_result_list.append(-sys.TEA.solve_price(sys.flowsheet.raw_wastewater)*3785411.78)
    LCA_result_list.append(sys.LCA.get_total_impacts(operation_only=True,
                                                 exclude=(sys.flowsheet.raw_wastewater,),
                                                 annual=True)['GlobalWarming']/sys.flowsheet.raw_wastewater.F_vol/2.3141471786573806)

CI_baseline = 2500

for i in range(25, 40):
    # TODO: assume 180 $/tonne CO2 removal
    TEA_result_list.append(TEA_result_list[i] - (CI_baseline - LCA_result_list[i])/1000*180)
    LCA_result_list.append(CI_baseline)

fig = plt.figure(figsize=(13, 12))
gs = GridSpec(2, 2, width_ratios=[4, 1.25], height_ratios=[1.25, 4], wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])

plt.rcParams['axes.linewidth'] = 3
plt.rcParams['hatch.linewidth'] = 3
plt.rcParams['xtick.labelsize'] = 36
plt.rcParams['ytick.labelsize'] = 36
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 1800)
ax.set_ylim(-1000, 2500)

ax.tick_params(direction='inout', length=20, width=3,
                bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 2100, 300), fontname='Arial')
plt.yticks(np.arange(-1000, 3000, 500), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
plt.xticks(np.arange(0, 2100, 300), fontname='Arial')

ax_top.tick_params(direction='in', length=10, width=3,
                   bottom=False, top=True, left=False, right=False,
                   labelcolor='none')

ax_bottom = ax.twinx()
ax_bottom.set_ylim(ax.get_ylim())
plt.xticks(np.arange(0, 2100, 300), fontname='Arial')

ax_bottom.tick_params(direction='in', length=10, width=3,
                      bottom=False, top=False, left=False, right=True,
                      labelcolor='none')

ax.set_xlabel(r'$\mathbf{Cost}$ ' + '[\$·tonne${^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{CI}$ ' + '[kg CO${_2}$ eq·tonne${^{-1}}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

scatter_C = ax.scatter(x=TEA_result_list[0:25], y=LCA_result_list[0:25], s=1200, color=a, linewidths=3, edgecolors='k')
scatter_T = ax.scatter(x=TEA_result_list[25:40], y=LCA_result_list[25:40], s=1200, color=r, linewidths=3, edgecolors='k')

for i in range(25, 40):
    ax.annotate('', xy=(TEA_result_list[i+15], LCA_result_list[i+15]), xytext=(TEA_result_list[i], LCA_result_list[i]),
            arrowprops=dict(arrowstyle='-|>', lw=3, color=y), size=15)

ax_right_boxplot = fig.add_subplot(gs[1, 1], sharey=ax)
ax_right_data = [LCA_result_list[0:25], LCA_result_list[25:40]]
bp = ax_right_boxplot.boxplot(ax_right_data, vert=True, whis=[5, 95], showfliers=False, widths=0.8, patch_artist=True)
    
bp['boxes'][0].set(color='k', facecolor=a, linewidth=3)
bp['boxes'][1].set(color='k', facecolor=r, linewidth=3)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp['medians']:
    median.set(color='k', linewidth=3)

for cap in bp['caps']:
    cap.set(color='k', linewidth=3)

ax_right_boxplot.axis('off')

ax_top_boxplot = fig.add_subplot(gs[0, 0], sharex=ax)
ax_top_data = [TEA_result_list[25:40], TEA_result_list[0:25]]
bp_top = ax_top_boxplot.boxplot(ax_top_data, vert=False, whis=[5, 95], showfliers=False, widths=0.8, patch_artist=True)

bp_top['boxes'][0].set(color='k', facecolor=r, linewidth=3)
bp_top['boxes'][1].set(color='k', facecolor=a, linewidth=3)

for whisker in bp_top['whiskers']:
    whisker.set(color='k', linewidth=3)

for median in bp_top['medians']:
    median.set(color='k', linewidth=3)

for cap in bp_top['caps']:
    cap.set(color='k', linewidth=3)

ax_top_boxplot.axis('off')

# plt.savefig('/Users/jiananfeng/Desktop/credit_opportunity.png', transparent=True, bbox_inches='tight')