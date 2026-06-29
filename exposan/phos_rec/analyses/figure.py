#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Xuan Wang <easonbiubiu99@gmail.com>
    
    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

#%% initialization

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib.colors as colors
from qsdsan.utils import auom
from matplotlib.mathtext import _mathtext as mathtext
from colorpalette import Color
from matplotlib.colors import to_hex
from matplotlib.patches import Rectangle

_ton_to_tonne = auom('ton').conversion_factor('tonne')

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

folder = os.path.dirname(os.path.dirname(__file__))

Fig_1_path = os.path.join(folder, 'data/Fig_1_data.xlsx')
Fig_1a_data = pd.read_excel(Fig_1_path, 'VFA')
Fig_1b_data = pd.read_excel(Fig_1_path, 'pH')
Fig_1c_data = pd.read_excel(Fig_1_path, 'PO43-')
Fig_1d_data = pd.read_excel(Fig_1_path, 'Fe2+')
Fig_1e_data = pd.read_excel(Fig_1_path, 'FePmolar')

Fig_2c_3ab_data = pd.read_excel(os.path.join(folder, 'results/sludge_management_cost_CI_baseline_2026-06-29.xlsx'), header=[1])
Fig_2c_3ab_data = Fig_2c_3ab_data.drop(0)
Fig_3c_data = pd.read_excel(os.path.join(folder, 'results/FePO4_result_cost_credit_2026-06-29.xlsx'))
Fig_3d_data = pd.read_excel(os.path.join(folder, 'results/FePO4_result_CI_credit_2026-06-29.xlsx'))

Fig_4ab_data = pd.read_excel(os.path.join(folder, 'results/decision_heatmap_FePO4_2026-06-29.xlsx'))
Fig_4ab_data = Fig_4ab_data[Fig_4ab_data['ratio']<=1]

Fig_5a_data = pd.read_excel(os.path.join(folder, 'results/context_heatmap_2026-06-29.xlsx'))

Fig_S2_data = pd.read_excel(os.path.join(folder, 'data/SI_NH4_data.xlsx'))

Fig_S3_data = pd.read_excel(os.path.join(folder, 'results/decision_heatmap_sludge_2026-06-29.xlsx'))
Fig_S3_data = Fig_S3_data[Fig_S3_data['ratio']<=1]

Fig_S4a_data = pd.read_excel(os.path.join(folder, 'results/sludge_cost_IRR_2026-06-29.xlsx'))
Fig_S4b_data = pd.read_excel(os.path.join(folder, 'results/FePO4_cost_IRR_2026-06-29.xlsx'))

Fig_S5_data = pd.read_excel(os.path.join(folder, 'results/sludge_result_size_2026-06-29.xlsx'))

Fig_S6_data = pd.read_excel(os.path.join(folder, 'results/FePO4_result_size_2026-06-29.xlsx'))

Fig_S7_data = pd.read_excel(os.path.join(folder, 'results/FePO4_VFA_CI_2026-06-29.xlsx'))

#%% Fig. 1a

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 144)
ax.set_ylim(0, 5000)

ax.set_xlabel(r'$\mathbf{Fermentation\ time}$ [h]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{VFA}$ [mg COD·L$^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')
plt.yticks(np.arange(0, 6000, 1000), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 6000, 1000), fontname='Arial')

plt.errorbar(Fig_1a_data['Time/h'],
             Fig_1a_data['FW:sludge 0:3'],
             yerr=Fig_1a_data['std'],
             lw=5, color=da, linestyle='dotted',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1a_data['Time/h'],
             Fig_1a_data['FW:sludge 1:3'],
             yerr=Fig_1a_data['std'],
             lw=5, color=da, linestyle='dashdot',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1a_data['Time/h'],
             Fig_1a_data['FW:sludge 2:3'],
             yerr=Fig_1a_data['std'],
             lw=5, color=da, linestyle='dashed',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1a_data['Time/h'],
             Fig_1a_data['FW:sludge 3:3'],
             yerr=Fig_1a_data['std'],
             lw=5, color=da, linestyle='solid',
             capsize=12, capthick=5, zorder=0)

plt.scatter(Fig_1a_data['Time/h'],
            Fig_1a_data['FW:sludge 0:3'],
            s=500, color=a, lw=5, edgecolor=da, marker='o', zorder=1)

plt.scatter(Fig_1a_data['Time/h'],
            Fig_1a_data['FW:sludge 1:3'],
            s=500, color=a, lw=5, edgecolor=da, marker='s', zorder=1)

plt.scatter(Fig_1a_data['Time/h'],
            Fig_1a_data['FW:sludge 2:3'],
            s=500, color=a, lw=5, edgecolor=da, marker='^', zorder=1)

plt.scatter(Fig_1a_data['Time/h'],
            Fig_1a_data['FW:sludge 3:3'],
            s=500, color=a, lw=5, edgecolor=da, marker='D', zorder=1)

plt.savefig(os.path.join(folder, 'figures/Fig_1a.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 1b

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 144)
ax.set_ylim(4, 7)

ax.set_xlabel(r'$\mathbf{Fermentation\ time}$ [h]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{pH}$',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')
plt.yticks(np.arange(4, 7.5, 0.5), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(4, 7.5, 0.5), fontname='Arial')

plt.errorbar(Fig_1b_data['Time/h'],
             Fig_1b_data['FW:sludge 0:3'],
             yerr=Fig_1b_data['std'],
             lw=5, color=dg, linestyle='dotted',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1b_data['Time/h'],
             Fig_1b_data['FW:sludge 1:3'],
             yerr=Fig_1b_data['std'],
             lw=5, color=dg, linestyle='dashdot',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1b_data['Time/h'],
             Fig_1b_data['FW:sludge 2:3'],
             yerr=Fig_1b_data['std'],
             lw=5, color=dg, linestyle='dashed',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1b_data['Time/h'],
             Fig_1b_data['FW:sludge 3:3'],
             yerr=Fig_1b_data['std'],
             lw=5, color=dg, linestyle='solid',
             capsize=12, capthick=5, zorder=0)

plt.scatter(Fig_1b_data['Time/h'],
            Fig_1b_data['FW:sludge 0:3'],
            s=500, color=g, lw=5, edgecolor=dg, marker='o', zorder=1)

plt.scatter(Fig_1b_data['Time/h'],
            Fig_1b_data['FW:sludge 1:3'],
            s=500, color=g, lw=5, edgecolor=dg, marker='s', zorder=1)

plt.scatter(Fig_1b_data['Time/h'],
            Fig_1b_data['FW:sludge 2:3'],
            s=500, color=g, lw=5, edgecolor=dg, marker='^', zorder=1)

plt.scatter(Fig_1b_data['Time/h'],
            Fig_1b_data['FW:sludge 3:3'],
            s=500, color=g, lw=5, edgecolor=dg, marker='D', zorder=1)

plt.savefig(os.path.join(folder, 'figures/Fig_1b.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 1c

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 144)
ax.set_ylim(0, 160)

ax.set_xlabel(r'$\mathbf{Fermentation\ time}$ [h]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{PO_4^{3-}\!-P}$ [mg·L$^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')
plt.yticks(np.arange(0, 200, 40), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 200, 40), fontname='Arial')

plt.errorbar(Fig_1c_data['Time/h'],
             Fig_1c_data['FW:sludge 0:3'],
             yerr=Fig_1c_data['std'],
             lw=5, color=db, linestyle='dotted',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1c_data['Time/h'],
             Fig_1c_data['FW:sludge 1:3'],
             yerr=Fig_1c_data['std'],
             lw=5, color=db, linestyle='dashdot',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1c_data['Time/h'],
             Fig_1c_data['FW:sludge 2:3'],
             yerr=Fig_1c_data['std'],
             lw=5, color=db, linestyle='dashed',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1c_data['Time/h'],
             Fig_1c_data['FW:sludge 3:3'],
             yerr=Fig_1c_data['std'],
             lw=5, color=db, linestyle='solid',
             capsize=12, capthick=5, zorder=0)

plt.scatter(Fig_1c_data['Time/h'],
            Fig_1c_data['FW:sludge 0:3'],
            s=500, color=b, lw=5, edgecolor=db, marker='o', zorder=1)

plt.scatter(Fig_1c_data['Time/h'],
            Fig_1c_data['FW:sludge 1:3'],
            s=500, color=b, lw=5, edgecolor=db, marker='s', zorder=1)

plt.scatter(Fig_1c_data['Time/h'],
            Fig_1c_data['FW:sludge 2:3'],
            s=500, color=b, lw=5, edgecolor=db, marker='^', zorder=1)

plt.scatter(Fig_1c_data['Time/h'],
            Fig_1c_data['FW:sludge 3:3'],
            s=500, color=b, lw=5, edgecolor=db, marker='D', zorder=1)

plt.savefig(os.path.join(folder, 'figures/Fig_1c.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 1d

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 144)
ax.set_ylim(0, 300)

ax.set_xlabel(r'$\mathbf{Fermentation\ time}$ [h]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{Fe^{2+}}$ [mg·L$^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')
plt.yticks(np.arange(0, 350, 50), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 350, 50), fontname='Arial')

plt.errorbar(Fig_1d_data['Time/h'],
             Fig_1d_data['FW:sludge 0:3'],
             yerr=Fig_1d_data['std'],
             lw=5, color=dr, linestyle='dotted',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1d_data['Time/h'],
             Fig_1d_data['FW:sludge 1:3'],
             yerr=Fig_1d_data['std'],
             lw=5, color=dr, linestyle='dashdot',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1d_data['Time/h'],
             Fig_1d_data['FW:sludge 2:3'],
             yerr=Fig_1d_data['std'],
             lw=5, color=dr, linestyle='dashed',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1d_data['Time/h'],
             Fig_1d_data['FW:sludge 3:3'],
             yerr=Fig_1d_data['std'],
             lw=5, color=dr, linestyle='solid',
             capsize=12, capthick=5, zorder=0)

plt.scatter(Fig_1d_data['Time/h'],
            Fig_1d_data['FW:sludge 0:3'],
            s=500, color=r, lw=5, edgecolor=dr, marker='o', zorder=1)

plt.scatter(Fig_1d_data['Time/h'],
            Fig_1d_data['FW:sludge 1:3'],
            s=500, color=r, lw=5, edgecolor=dr, marker='s', zorder=1)

plt.scatter(Fig_1d_data['Time/h'],
            Fig_1d_data['FW:sludge 2:3'],
            s=500, color=r, lw=5, edgecolor=dr, marker='^', zorder=1)

plt.scatter(Fig_1d_data['Time/h'],
            Fig_1d_data['FW:sludge 3:3'],
            s=500, color=r, lw=5, edgecolor=dr, marker='D', zorder=1)

plt.savefig(os.path.join(folder, 'figures/Fig_1d.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 1e

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 144)
ax.set_ylim(0, 1.2)

ax.set_xlabel(r'$\mathbf{Fermentation\ time}$ [h]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{Fe:P\ molar\ ratio}$',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')
plt.yticks(np.arange(0, 1.4, 0.2), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 1.4, 0.2), fontname='Arial')

plt.errorbar(Fig_1e_data['Time/h'],
             Fig_1e_data['FW:sludge 0:3'],
             yerr=Fig_1e_data['std'],
             lw=5, color=dp, linestyle='dotted',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1e_data['Time/h'],
             Fig_1e_data['FW:sludge 1:3'],
             yerr=Fig_1e_data['std'],
             lw=5, color=dp, linestyle='dashdot',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1e_data['Time/h'],
             Fig_1e_data['FW:sludge 2:3'],
             yerr=Fig_1e_data['std'],
             lw=5, color=dp, linestyle='dashed',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_1e_data['Time/h'],
             Fig_1e_data['FW:sludge 3:3'],
             yerr=Fig_1e_data['std'],
             lw=5, color=dp, linestyle='solid',
             capsize=12, capthick=5, zorder=0)

plt.scatter(Fig_1e_data['Time/h'],
            Fig_1e_data['FW:sludge 0:3'],
            s=500, color=p, lw=5, edgecolor=dp, marker='o', zorder=1)

plt.scatter(Fig_1e_data['Time/h'],
            Fig_1e_data['FW:sludge 1:3'],
            s=500, color=p, lw=5, edgecolor=dp, marker='s', zorder=1)

plt.scatter(Fig_1e_data['Time/h'],
            Fig_1e_data['FW:sludge 2:3'],
            s=500, color=p, lw=5, edgecolor=dp, marker='^', zorder=1)

plt.scatter(Fig_1e_data['Time/h'],
            Fig_1e_data['FW:sludge 3:3'],
            s=500, color=p, lw=5, edgecolor=dp, marker='D', zorder=1)

plt.savefig(os.path.join(folder, 'figures/Fig_1e.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 2c Fe and P recovery

fig, ax = plt.subplots(figsize=(4, 8))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_ylim(0, 100)

ax.set_ylabel(r'$\mathbf{Recovery}$ [%]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=False, top=False, left=True, right=False)

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

bp = plt.boxplot([Fig_2c_3ab_data['P recovery [-]'].dropna()*100,
                  Fig_2c_3ab_data['Fe recovery [-]'].dropna()*100],
                 whis=[5, 95], showfliers=False, widths=0.6, patch_artist=True)

plt.xticks([1, 2], ('P','Fe'), fontname='Arial')

bp['boxes'][0].set(color='k', facecolor=b, linewidth=5)
bp['boxes'][1].set(color='k', facecolor=r, linewidth=5)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=5)

for median in bp['medians']:
    median.set(color='k', linewidth=5)
    
for cap in bp['caps']:
    cap.set(color='k', linewidth=5)

ax.plot((1.5, 1.5), (0, 100), color='k', linewidth=5, linestyle='--')

rectangle_fill = Rectangle((0.5, 73.1859), 1, 99.3-73.1859,
                           fc=a, alpha=0.5, zorder=0)
ax.add_patch(rectangle_fill)

plt.savefig(os.path.join(folder, 'figures/Fig_2c.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 3a sludge management cost

fig, ax = plt.subplots(figsize=(2.5, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0.5, 1.5)
ax.set_ylim(0, 1000)

ax.set_ylabel(r'$\mathbf{Sludge\ management\ cost}$'+'\n'+r'[$\$·tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5,
              linespacing=0.8)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=False, top=False, left=True, right=False,
               labelbottom=False)

plt.yticks(np.arange(0, 1200, 200), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 1200, 200), fontname='Arial')

bp = plt.boxplot(Fig_2c_3ab_data['Sludge management cost [$/tonne]'].dropna(),
                 whis=[5, 95], showfliers=False, widths=0.6, patch_artist=True)

for box in bp['boxes']:
    box.set(color='k', facecolor=o, linewidth=5)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=5)

for median in bp['medians']:
    median.set(color='k', linewidth=5)
    
for cap in bp['caps']:
    cap.set(color='k', linewidth=5)

# conventional sludge management: 100 to 800 $·ton-1
# Peccia, J. & Westerhoff, P. We Should Expect More out of Our Sewage Sludge. Environ. Sci. Technol. 49, 8271–8276 (2015).
rectangle_fill = Rectangle((0.5, 100/_ton_to_tonne), 1.5, (800-100)/_ton_to_tonne,
                           fc=a, alpha=0.5, zorder=0)
ax.add_patch(rectangle_fill)

plt.savefig(os.path.join(folder, 'figures/Fig_3a.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 3b sludge management CI

fig, ax = plt.subplots(figsize=(2.5, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0.5, 1.5)
ax.set_ylim(0, 4000)

ax.set_ylabel(r'$\mathbf{Sludge\ management\ CI}$'+'\n'+r'[$kg\ CO_2e·tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5,
              linespacing=0.8)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=False, top=False, left=True, right=False,
               labelbottom=False)

plt.yticks(np.arange(0, 5000, 1000), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 5000, 1000), fontname='Arial')

bp = plt.boxplot(Fig_2c_3ab_data['Sludge management GWP [kg_CO2_eq/tonne]'].dropna(),
                 whis=[5, 95], showfliers=False, widths=0.6, patch_artist=True)

for box in bp['boxes']:
    box.set(color='k', facecolor=o, linewidth=5)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=5)

for median in bp['medians']:
    median.set(color='k', linewidth=5)
    
for cap in bp['caps']:
    cap.set(color='k', linewidth=5)

rectangle_fill = Rectangle((0.5, 414), 1.5, 3618-414,
                           fc=a, alpha=0.5, zorder=0)
ax.add_patch(rectangle_fill)

plt.savefig(os.path.join(folder, 'figures/Fig_3b.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 3c FePO4 MSP

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 1000)
ax.set_ylim(-15, 15)

ax.set_xlabel(r'$\mathbf{Avoided\ management\ cost}$'+'\n'r'[$\$·tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0,
              linespacing=0.8)

ax.set_ylabel(r'$\mathbf{FePO_4\ MSP}$ [$\$·kg^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 1200, 200), fontname='Arial')
plt.yticks(np.arange(-15, 20, 5), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0, 1200, 200), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(-15, 20, 5), fontname='Arial')

x = Fig_3c_data['credit']
low = Fig_3c_data['FePO4_MSP_5th']
mid = Fig_3c_data['FePO4_MSP_50th']
high = Fig_3c_data['FePO4_MSP_95th']

ax.fill_between(x, low, high, color=p, alpha=0.5)

ax.plot(x, low, color=dp, linestyle='--', linewidth=5)
ax.plot(x, high, color=dp, linestyle='--', linewidth=5)
ax.plot(x, mid, color=dp, linewidth=5)

# conventional sludge management: 100 to 800 $·ton-1
# Peccia, J. & Westerhoff, P. We Should Expect More out of Our Sewage Sludge. Environ. Sci. Technol. 49, 8271–8276 (2015).
rectangle_fill = Rectangle((100/_ton_to_tonne, -15), (800-100)/_ton_to_tonne, 35, 
                           fc=a, alpha=0.5, zorder=0)
ax.add_patch(rectangle_fill)

rectangle_fill = Rectangle((0, 1.327), 1000, 3.89-1.327, 
                           fc=g, alpha=0.5, zorder=0)
ax.add_patch(rectangle_fill)

rectangle_fill = Rectangle((0, 4.55), 1000, 6.97-4.55, 
                           fc=dg, alpha=0.8, zorder=0)
ax.add_patch(rectangle_fill)

plt.savefig(os.path.join(folder, 'figures/Fig_3c.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 3d FePO4 CI

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 4000)
ax.set_ylim(-50, 50)

ax.set_xlabel(r'$\mathbf{Avoided\ management\ CI}$'+'\n'+r' [$kg\ CO_2e·tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0,
              linespacing=0.8)

ax.set_ylabel(r'$\mathbf{FePO_4\ CI}$ [$kg·kg^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 5000, 1000), fontname='Arial')
plt.yticks(np.arange(-50, 75, 25), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0, 5000, 1000), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(-50, 75, 25), fontname='Arial')

x = Fig_3d_data['credit']
low = Fig_3d_data['FePO4_GWP_5th']
mid = Fig_3d_data['FePO4_GWP_50th']
high = Fig_3d_data['FePO4_GWP_95th']

ax.fill_between(x, low, high, color=p, alpha=0.5)

ax.plot(x, low, color=dp, linestyle='--', linewidth=5)
ax.plot(x, high, color=dp, linestyle='--', linewidth=5)
ax.plot(x, mid, color=dp, linewidth=5)

rectangle_fill = Rectangle((414, -50), 3618-414, 100,
                           fc=a, alpha=0.5, zorder=0)
ax.add_patch(rectangle_fill)

ax.plot((0, 4000), (16.6, 16.6), color=g, alpha=0.5, linewidth=5)

ax.plot((0, 4000), (22.7, 22.7), color=dg, alpha=0.8, linewidth=5)

plt.savefig(os.path.join(folder, 'figures/Fig_3d.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 4a FePO4 MSP technological heatmap

fig, ax = plt.subplots(figsize=(12.5, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(1, 3)
ax.set_ylim(12, 132)

ax.set_xlabel(r'$\mathbf{Food\ waste:sludge\ ratio}$',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{Fermentation\ time}$ [h]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(1, 4, 1), ('1:3','2:3','3:3'), fontname='Arial')
plt.yticks(np.arange(12, 144, 24), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(1, 4, 1), ('1:3','2:3','3:3'), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(12, 144, 24), fontname='Arial')

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [dp, p, 'w'][::-1])

X = np.array(Fig_4ab_data['ratio'].map({1/3: 1,
                                        2/3: 2,
                                        1: 3}))
Y = np.array(Fig_4ab_data['fermentation_time'])
Z = np.array(Fig_4ab_data['FePO4_MSP_50th'])

fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)

fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(X, Y, Z, levels=7, linewidths=5, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=45)

plt.rcParams['hatch.color'] = g

ax.tricontourf(X, Y, Z, levels=[1.327, 3.89], colors=g, alpha=0.5)

plt.rcParams['hatch.color'] = dg

ax.tricontourf(X, Y, Z, levels=[4.55, 6.97], colors=dg, alpha=0.8)

plt.savefig(os.path.join(folder, 'figures/Fig_4a.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 4b FePO4 CI technological heatmap

fig, ax = plt.subplots(figsize=(12.5, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(1, 3)
ax.set_ylim(12, 132)

ax.set_xlabel(r'$\mathbf{Food\ waste:sludge\ ratio}$',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{Fermentation\ time}$ [h]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(1, 4, 1), ('1:3','2:3','3:3'), fontname='Arial')
plt.yticks(np.arange(12, 144, 24), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(1, 4, 1), ('1:3','2:3','3:3'), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(12, 144, 24), fontname='Arial')

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [dp, p, 'w'][::-1])


X = np.array(Fig_4ab_data['ratio'].map({1/3: 1,
                                        2/3: 2,
                                        1: 3}))
Y = np.array(Fig_4ab_data['fermentation_time'])
Z = np.array(Fig_4ab_data['FePO4_GWP_50th'])

fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)

fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(X, Y, Z, levels=7, linewidths=5, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=45)

ax.tricontour(X, Y, Z, levels=[16.6, 22.7], colors=[g, dg], linewidths=5)

plt.savefig(os.path.join(folder, 'figures/Fig_4b.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. 5a FePO4 MSP contextual heatmap

fig, ax = plt.subplots(figsize=(12.5, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 120)
ax.set_ylim(0, 3)

ax.set_xlabel(r'$\mathbf{Disposal\ cost}$ [$\$·tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{VFA-rich\ stream\ price}$'+'\n'+r'[$\$·m^{-3}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5,
              linespacing=0.8)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 140, 20),  fontname='Arial')
plt.yticks(np.arange(0, 3.5, 0.5), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0, 140, 20),  fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 3.5, 0.5), fontname='Arial')

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [dp, p, 'w'][::-1])

X = np.array(Fig_5a_data['residue_cost']*1000)
Y = np.array(Fig_5a_data['VFA_price']*1000)
Z = np.array(Fig_5a_data['FePO4_MSP_50th'])

fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)

fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(X, Y, Z, levels=6, linewidths=5, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=45)

plt.rcParams['hatch.color'] = g

ax.tricontourf(X, Y, Z, levels=[1.327, 3.89], colors=g, alpha=0.5)

plt.rcParams['hatch.color'] = dg

ax.tricontourf(X, Y, Z, levels=[4.55, 6.97], colors=dg, alpha=0.8)

plt.savefig(os.path.join(folder, 'figures/Fig_5a.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. S2 NH4+

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 144)
ax.set_ylim(0, 100)

ax.set_xlabel(r'$\mathbf{Fermentation\ time}$ [h]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{NH_4^{+}-N}$ [mg·L$^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')
plt.yticks(np.arange(0, 120, 20), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0, 156, 24), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

plt.errorbar(Fig_S2_data['Time/h'],
             Fig_S2_data['FW:sludge 0:3'],
             yerr=Fig_S2_data['std'],
             lw=5, color=dy, linestyle='dotted',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_S2_data['Time/h'],
             Fig_S2_data['FW:sludge 1:3'],
             yerr=Fig_S2_data['std'],
             lw=5, color=dy, linestyle='dashdot',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_S2_data['Time/h'],
             Fig_S2_data['FW:sludge 2:3'],
             yerr=Fig_S2_data['std'],
             lw=5, color=dy, linestyle='dashed',
             capsize=12, capthick=5, zorder=0)

plt.errorbar(Fig_S2_data['Time/h'],
             Fig_S2_data['FW:sludge 3:3'],
             yerr=Fig_S2_data['std'],
             lw=5, color=dy, linestyle='solid',
             capsize=12, capthick=5, zorder=0)

plt.scatter(Fig_S2_data['Time/h'],
            Fig_S2_data['FW:sludge 0:3'],
            s=500, color=y, lw=5, edgecolor=dy, marker='o', zorder=1)

plt.scatter(Fig_S2_data['Time/h'],
            Fig_S2_data['FW:sludge 1:3'],
            s=500, color=y, lw=5, edgecolor=dy, marker='s', zorder=1)

plt.scatter(Fig_S2_data['Time/h'],
            Fig_S2_data['FW:sludge 2:3'],
            s=500, color=y, lw=5, edgecolor=dy, marker='^', zorder=1)

plt.scatter(Fig_S2_data['Time/h'],
            Fig_S2_data['FW:sludge 3:3'],
            s=500, color=y, lw=5, edgecolor=dy, marker='D', zorder=1)

plt.savefig(os.path.join(folder, 'figures/Fig_S2.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. S3a sludge cost technological heatmap

fig, ax = plt.subplots(figsize=(12.5, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(1, 3)
ax.set_ylim(12, 132)

ax.set_xlabel(r'$\mathbf{Food\ waste:sludge\ ratio}$',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{Fermentation\ time}$ [h]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(1, 4, 1), ('1:3','2:3','3:3'), fontname='Arial')
plt.yticks(np.arange(12, 144, 24), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(1, 4, 1), ('1:3','2:3','3:3'), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(12, 144, 24), fontname='Arial')

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [do, o, 'w'][::-1])

X = np.array(Fig_S3_data['ratio'].map({1/3: 1,
                                       2/3: 2,
                                       1: 3}))
Y = np.array(Fig_S3_data['fermentation_time'])
Z = np.array(Fig_S3_data['sludge_cost_50th'])

fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)

fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(X, Y, Z, levels=7, linewidths=5, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=45)

plt.savefig(os.path.join(folder, 'figures/Fig_S3a.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. S3b sludge CI technological heatmap

fig, ax = plt.subplots(figsize=(12.5, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(1, 3)
ax.set_ylim(12, 132)

ax.set_xlabel(r'$\mathbf{Food\ waste:sludge\ ratio}$',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{Fermentation\ time}$ [h]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(1, 4, 1), ('1:3','2:3','3:3'), fontname='Arial')
plt.yticks(np.arange(12, 144, 24), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(1, 4, 1), ('1:3','2:3','3:3'), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(12, 144, 24), fontname='Arial')

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [do, o, 'w'][::-1])

X = np.array(Fig_S3_data['ratio'].map({1/3: 1,
                                       2/3: 2,
                                       1: 3}))
Y = np.array(Fig_S3_data['fermentation_time'])
Z = np.array(Fig_S3_data['sludge_GWP_50th'])

fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)

fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(X, Y, Z, levels=7, linewidths=5, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=45)

plt.savefig(os.path.join(folder, 'figures/Fig_S3b.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. S4a sludge cost IRR

fig, ax = plt.subplots(figsize=(2.5, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 0.05)
ax.set_ylim(0, 1000)

ax.set_xlabel(r'$\mathbf{IRR}$',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{Sludge\ management\ cost}$'+'\n'+r'[$\$·tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5,
              linespacing=0.8)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks([0, 0.05], fontname='Arial')
plt.yticks(np.arange(0, 1200, 200), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks([0, 0.05], fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 1200, 200), fontname='Arial')

x = Fig_S4a_data['IRR']
low = Fig_S4a_data['sludge_cost_5th']
mid = Fig_S4a_data['sludge_cost_50th']
high = Fig_S4a_data['sludge_cost_95th']

ax.fill_between(x, low, high, color=o, alpha=0.5)

ax.plot(x, low, color=do, linestyle='--', linewidth=5)
ax.plot(x, high, color=do, linestyle='--', linewidth=5)
ax.plot(x, mid, color=do, linewidth=5)

# conventional sludge management: 100 to 800 $·ton-1
# Peccia, J. & Westerhoff, P. We Should Expect More out of Our Sewage Sludge. Environ. Sci. Technol. 49, 8271–8276 (2015).
rectangle_fill = Rectangle((0, 100/_ton_to_tonne), 0.05, (800-100)/_ton_to_tonne,
                           fc=a, alpha=0.5, zorder=0)
ax.add_patch(rectangle_fill)

plt.savefig(os.path.join(folder, 'figures/Fig_S4a.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. S4a FePO4 MSP IRR

fig, ax = plt.subplots(figsize=(7.5, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0.05, 0.2)
ax.set_ylim(1, 7)

ax.set_xlabel(r'$\mathbf{IRR}$',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{FePO_4\ MSP}$ [$\$·kg^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0.05, 0.25, 0.05), fontname='Arial')
plt.yticks(np.arange(1, 8, 1), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0.05, 0.25, 0.05), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(1, 8, 1), fontname='Arial')

x = Fig_S4b_data['IRR']
low = Fig_S4b_data['FePO4_MSP_5th']
mid = Fig_S4b_data['FePO4_MSP_50th']
high = Fig_S4b_data['FePO4_MSP_95th']

ax.fill_between(x, low, high, color=p, alpha=0.5)

ax.plot(x, low, color=dp, linestyle='--', linewidth=5)
ax.plot(x, high, color=dp, linestyle='--', linewidth=5)
ax.plot(x, mid, color=dp, linewidth=5)

rectangle_fill = Rectangle((0.05, 1.327), 0.15, 3.89-1.327, 
                           fc=g, alpha=0.5, zorder=0)
ax.add_patch(rectangle_fill)

rectangle_fill = Rectangle((0.05, 4.55), 0.15, 6.97-4.55, 
                           fc=dg, alpha=0.8, zorder=0)
ax.add_patch(rectangle_fill)

plt.savefig(os.path.join(folder, 'figures/Fig_S4b.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. S5a sludge cost size

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(1, 256)
ax.set_ylim(0, 1000)

ax.set_xlabel(r'$\mathbf{Sludge\ throughput}$'+'\n'+r'[$dry\ tonne·day^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0,
              linespacing=0.8)

ax.set_ylabel(r'$\mathbf{Sludge\ management\ cost}$'+'\n'+r'[$\$·tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5,
              linespacing=0.8)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(fontname='Arial')
plt.yticks(np.arange(0, 1200, 200), fontname='Arial')

plt.xscale('log')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 1200, 200), fontname='Arial')

x = Fig_S5_data['size']
low = Fig_S5_data['sludge_cost_5th']
mid = Fig_S5_data['sludge_cost_50th']
high = Fig_S5_data['sludge_cost_95th']

ax.fill_between(x, low, high, color=o, alpha=0.5, zorder=1)

ax.plot(x, low, color=do, linestyle='--', linewidth=5)
ax.plot(x, high, color=do, linestyle='--', linewidth=5)
ax.plot(x, mid, color=do, linewidth=5)

# conventional sludge management: 100 to 800 $·ton-1
# Peccia, J. & Westerhoff, P. We Should Expect More out of Our Sewage Sludge. Environ. Sci. Technol. 49, 8271–8276 (2015).
rectangle_fill = Rectangle((1, 100/_ton_to_tonne), 255, (800-100)/_ton_to_tonne,
                           fc=a, alpha=0.5, zorder=0)
ax.add_patch(rectangle_fill)

plt.savefig(os.path.join(folder, 'figures/Fig_S5a.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. S5b sludge CI size

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(1, 256)
ax.set_ylim(0, 4000)

ax.set_xlabel(r'$\mathbf{Sludge\ throughput}$'+'\n'+r'[$dry\ tonne·day^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0,
              linespacing=0.8)

ax.set_ylabel(r'$\mathbf{Sludge\ management\ CI}$'+'\n'+r'[$kg\ CO_2e·tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5,
              linespacing=0.8)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(fontname='Arial')
plt.yticks(np.arange(0, 5000, 1000), fontname='Arial')

plt.xscale('log')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 5000, 1000), fontname='Arial')

x = Fig_S5_data['size']
low = Fig_S5_data['sludge_GWP_5th']
mid = Fig_S5_data['sludge_GWP_50th']
high = Fig_S5_data['sludge_GWP_95th']

ax.fill_between(x, low, high, color=o, alpha=0.5)

ax.plot(x, low, color=do, linestyle='--', linewidth=5)
ax.plot(x, high, color=do, linestyle='--', linewidth=5)
ax.plot(x, mid, color=do, linewidth=5)

rectangle_fill = Rectangle((1, 414), 255, 3618-414,
                           fc=a, alpha=0.5, zorder=0)
ax.add_patch(rectangle_fill)

plt.savefig(os.path.join(folder, 'figures/Fig_S5b.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. S6a FePO4 MSP size

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(1, 256)
ax.set_ylim(0, 18)

ax.set_xlabel(r'$\mathbf{Sludge\ throughput}$'+'\n'+r'[$dry\ tonne·day^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0,
              linespacing=0.8)

ax.set_ylabel(r'$\mathbf{FePO_4\ MSP}$ [$\$·kg^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(fontname='Arial')
plt.yticks(np.arange(0, 21, 3), fontname='Arial')

plt.xscale('log')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 21, 3), fontname='Arial')

x = Fig_S6_data['size']
low = Fig_S6_data['FePO4_MSP_5th']
mid = Fig_S6_data['FePO4_MSP_50th']
high = Fig_S6_data['FePO4_MSP_95th']

ax.fill_between(x, low, high, color=p, alpha=0.5)

ax.plot(x, low, color=dp, linestyle='--', linewidth=5)
ax.plot(x, high, color=dp, linestyle='--', linewidth=5)
ax.plot(x, mid, color=dp, linewidth=5)

rectangle_fill = Rectangle((1, 1.327), 255, 3.89-1.327,
                           fc=g, alpha=0.5, zorder=0)
ax.add_patch(rectangle_fill)

rectangle_fill = Rectangle((1, 4.55), 255, 6.97-4.55, 
                           fc=dg, alpha=0.8, zorder=0)
ax.add_patch(rectangle_fill)

plt.savefig(os.path.join(folder, 'figures/Fig_S6a.pdf'), transparent=True, bbox_inches='tight')

#%% Fig. S6b FePO4 CI size

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(1, 256)
ax.set_ylim(13, 25)

ax.set_xlabel(r'$\mathbf{Sludge\ throughput}$'+'\n'+r'[$dry\ tonne·day^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0,
              linespacing=0.8)

ax.set_ylabel(r'$\mathbf{FePO_4\ CI}$ [$kg\ CO_2e·kg^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(fontname='Arial')
plt.yticks(np.arange(13, 28, 3), fontname='Arial')

plt.xscale('log')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(13, 28, 3), fontname='Arial')

x = Fig_S6_data['size']
low = Fig_S6_data['FePO4_GWP_5th']
mid = Fig_S6_data['FePO4_GWP_50th']
high = Fig_S6_data['FePO4_GWP_95th']

ax.fill_between(x, low, high, color=p, alpha=0.5)

ax.plot(x, low, color=dp, linestyle='--', linewidth=5)
ax.plot(x, high, color=dp, linestyle='--', linewidth=5)
ax.plot(x, mid, color=dp, linewidth=5)

ax.plot((1, 256), (16.6, 16.6), color=g, alpha=0.5, linewidth=5)

ax.plot((1, 256), (22.7, 22.7), color=dg, alpha=0.8, linewidth=5)

plt.savefig(os.path.join(folder, 'figures/Fig_S6b.pdf'), transparent=True, bbox_inches='tight')

#%% optional FePO4 CI VFA

fig, ax = plt.subplots(figsize=(10, 10))

plt.rcParams['axes.linewidth'] = 5
plt.rcParams['hatch.linewidth'] = 5
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] = 45
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams.update({'mathtext.fontset':'custom'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.bf':'Arial: bold'})

ax.set_xlim(0, 10)
ax.set_ylim(-20, 30)

ax.set_xlabel(r'$\mathbf{VFA\ CI}$ [$kg\ CO_2e·m^{-3}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0,
              linespacing=0.8)

ax.set_ylabel(r'$\mathbf{FePO_4\ CI}$ [$kg\ CO_2e·kg^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(fontname='Arial')
plt.yticks(np.arange(-20, 40, 10), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(-20, 40, 10), fontname='Arial')

x = Fig_S7_data['VFA_GWP']*1000
low = Fig_S7_data['FePO4_GWP_5th']
mid = Fig_S7_data['FePO4_GWP_50th']
high = Fig_S7_data['FePO4_GWP_95th']

ax.fill_between(x, low, high, color=p, alpha=0.5)

ax.plot(x, low, color=dp, linestyle='--', linewidth=5)
ax.plot(x, high, color=dp, linestyle='--', linewidth=5)
ax.plot(x, mid, color=dp, linewidth=5)

ax.plot((0, 256), (16.6, 16.6), color=g, alpha=0.5, linewidth=5)

ax.plot((0, 256), (22.7, 22.7), color=dg, alpha=0.8, linewidth=5)