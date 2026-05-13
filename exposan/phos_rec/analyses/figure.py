#!/usr/bin/env python3 for the CI and Cost heatmap of food_sludge_ratio and HRT
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

import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib.colors as colors
from matplotlib.mathtext import _mathtext as mathtext
from colorpalette import Color
from matplotlib.colors import to_hex

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

Figure_2_file = r'D:\1_20250901-20260730\UIUC\writing_paper\EST-phosphorus_recovery_to_FePO4\Figure_2_data.xlsx'
Figure_2A_data = pd.read_excel(Figure_2_file, 'VFA')
Figure_2B_data = pd.read_excel(Figure_2_file, 'pH')
Figure_2C_data = pd.read_excel(Figure_2_file, 'PO43-')
Figure_2D_data = pd.read_excel(Figure_2_file, 'Fe2+')
Figure_2E_data = pd.read_excel(Figure_2_file, 'FePmolar')

Figure_3B_P_recovery_data = pd.read_excel(Figure_2_file, 'Precovery')
Figure_3B_Fe_recovery_data = pd.read_excel(Figure_2_file, 'Ferecovery')
Figure_3B_sludge_management_cost_data = pd.read_excel(Figure_2_file, 'SludgeManagementCost')
Figure_3B_sludge_management_CI_data = pd.read_excel(Figure_2_file, 'SludgeManagementCI')
Figure_3C_avoid_waste_sludge_management_cost_data = pd.read_excel(Figure_2_file, 'AvoidWMC_MSPFePO4')
Figure_3D_avoid_waste_sludge_management_CI_data = pd.read_excel(Figure_2_file, 'AvoidWMCI_CIFePO4')

Figure_4AB_data = pd.read_excel(Figure_2_file, 'FePO4CostCI')
Figure_4CD_data = pd.read_excel(Figure_2_file, 'sludgeCostCI')
Figure_4EF_data = pd.read_excel(Figure_2_file, 'VFA_residue_FePO4MSPCI')

# SI
Figure_S1_NH4_release_data = pd.read_excel(Figure_2_file, 'SI_NH4+')

#%% Figure 2A

fig, ax = plt.subplots(figsize=(14, 12))

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

ax.set_xlim(0, 132)
ax.set_ylim(0, 5000)
plt.subplots_adjust(left=0.20, bottom=0.13)

ax.set_xlabel(r'$\mathbf{Time}$ [hr]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{VFAs}$ [mg COD·L$^{-1}$]',
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

plt.errorbar(Figure_2A_data['Time/h'],
            Figure_2A_data['FW:sludge 0:3'],
            yerr=Figure_2A_data['std'],
            lw=5, color=da, linestyle='dotted',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2A_data['Time/h'],
            Figure_2A_data['FW:sludge 1:3'],
            yerr=Figure_2A_data['std'],
            lw=5, color=da, linestyle='dashdot',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2A_data['Time/h'],
            Figure_2A_data['FW:sludge 2:3'],
            yerr=Figure_2A_data['std'],
            lw=5, color=da, linestyle='dashed',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2A_data['Time/h'],
            Figure_2A_data['FW:sludge 3:3'],
            yerr=Figure_2A_data['std'],
            lw=5, color=da, linestyle='solid',
            capsize=12, capthick=5, zorder=0)

plt.scatter(Figure_2A_data['Time/h'],
            Figure_2A_data['FW:sludge 0:3'],
            s=500, color=a, lw=5, edgecolor=da, marker='o', zorder=1)

plt.scatter(Figure_2A_data['Time/h'],
            Figure_2A_data['FW:sludge 1:3'],
            s=500, color=a, lw=5, edgecolor=da,marker='s', zorder=1)

plt.scatter(Figure_2A_data['Time/h'],
            Figure_2A_data['FW:sludge 2:3'],
            s=500, color=a, lw=5, edgecolor=da,marker='^', zorder=1)

plt.scatter(Figure_2A_data['Time/h'],
            Figure_2A_data['FW:sludge 3:3'],
            s=500, color=a, lw=5, edgecolor=da,marker='D', zorder=1)

#%% Figure 2B

fig, ax = plt.subplots(figsize=(13, 12))

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

ax.set_xlim(0, 132)
ax.set_ylim(4, 7)
plt.subplots_adjust(left=0.15, bottom=0.13)

ax.set_xlabel(r'$\mathbf{Time}$ [hr]',
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

plt.errorbar(Figure_2B_data['Time/h'],
            Figure_2B_data['FW:sludge 0:3'],
            yerr=Figure_2B_data['std'],
            lw=5, color=dg, linestyle='dotted',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2B_data['Time/h'],
            Figure_2B_data['FW:sludge 1:3'],
            yerr=Figure_2B_data['std'],
            lw=5, color=dg, linestyle='dashdot',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2B_data['Time/h'],
            Figure_2B_data['FW:sludge 2:3'],
            yerr=Figure_2B_data['std'],
            lw=5, color=dg, linestyle='dashed',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2B_data['Time/h'],
            Figure_2B_data['FW:sludge 3:3'],
            yerr=Figure_2B_data['std'],
            lw=5, color=dg, linestyle='solid',
            capsize=12, capthick=5, zorder=0)

plt.scatter(Figure_2B_data['Time/h'],
            Figure_2B_data['FW:sludge 0:3'],
            s=500, color=g, lw=5, edgecolor=dg, marker='o', zorder=1)

plt.scatter(Figure_2B_data['Time/h'],
            Figure_2B_data['FW:sludge 1:3'],
            s=500, color=g, lw=5, edgecolor=dg,marker='s', zorder=1)

plt.scatter(Figure_2B_data['Time/h'],
            Figure_2B_data['FW:sludge 2:3'],
            s=500, color=g, lw=5, edgecolor=dg,marker='^', zorder=1)

plt.scatter(Figure_2B_data['Time/h'],
            Figure_2B_data['FW:sludge 3:3'],
            s=500, color=g, lw=5, edgecolor=dg,marker='D', zorder=1)

#%% Figure 2C

fig, ax = plt.subplots(figsize=(14, 12))

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

ax.set_xlim(0, 132)
ax.set_ylim(0, 160)
plt.subplots_adjust(left=0.17, bottom=0.13)

ax.set_xlabel(r'$\mathbf{Time}$ [hr]',
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

plt.errorbar(Figure_2C_data['Time/h'],
            Figure_2C_data['FW:sludge 0:3'],
            yerr=Figure_2C_data['std'],
            lw=5, color=db, linestyle='dotted',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2C_data['Time/h'],
            Figure_2C_data['FW:sludge 1:3'],
            yerr=Figure_2C_data['std'],
            lw=5, color=db, linestyle='dashdot',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2C_data['Time/h'],
            Figure_2C_data['FW:sludge 2:3'],
            yerr=Figure_2C_data['std'],
            lw=5, color=db, linestyle='dashed',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2C_data['Time/h'],
            Figure_2C_data['FW:sludge 3:3'],
            yerr=Figure_2C_data['std'],
            lw=5, color=db, linestyle='solid',
            capsize=12, capthick=5, zorder=0)

plt.scatter(Figure_2C_data['Time/h'],
            Figure_2C_data['FW:sludge 0:3'],
            s=500, color=b, lw=5, edgecolor=db, marker='o', zorder=1)

plt.scatter(Figure_2C_data['Time/h'],
            Figure_2C_data['FW:sludge 1:3'],
            s=500, color=b, lw=5, edgecolor=db,marker='s', zorder=1)

plt.scatter(Figure_2C_data['Time/h'],
            Figure_2C_data['FW:sludge 2:3'],
            s=500, color=b, lw=5, edgecolor=db,marker='^', zorder=1)

plt.scatter(Figure_2C_data['Time/h'],
            Figure_2C_data['FW:sludge 3:3'],
            s=500, color=b, lw=5, edgecolor=db,marker='D', zorder=1)

#%% Figure 2D

fig, ax = plt.subplots(figsize=(14, 12))

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

ax.set_xlim(0, 132)
ax.set_ylim(0, 300)
# TODO: add bottom to others, if necessary
plt.subplots_adjust(bottom=0.122, left=0.15)

ax.set_xlabel(r'$\mathbf{Time}$ [hr]',
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

plt.errorbar(Figure_2D_data['Time/h'],
            Figure_2D_data['FW:sludge 0:3'],
            yerr=Figure_2D_data['std'],
            lw=5, color=dr, linestyle='dotted',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2D_data['Time/h'],
            Figure_2D_data['FW:sludge 1:3'],
            yerr=Figure_2D_data['std'],
            lw=5, color=dr, linestyle='dashdot',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2D_data['Time/h'],
            Figure_2D_data['FW:sludge 2:3'],
            yerr=Figure_2D_data['std'],
            lw=5, color=dr, linestyle='dashed',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2D_data['Time/h'],
            Figure_2D_data['FW:sludge 3:3'],
            yerr=Figure_2D_data['std'],
            lw=5, color=dr, linestyle='solid',
            capsize=12, capthick=5, zorder=0)

plt.scatter(Figure_2D_data['Time/h'],
            Figure_2D_data['FW:sludge 0:3'],
            s=500, color=r, lw=5, edgecolor=dr, marker='o', zorder=1)

plt.scatter(Figure_2D_data['Time/h'],
            Figure_2D_data['FW:sludge 1:3'],
            s=500, color=r, lw=5, edgecolor=dr,marker='s', zorder=1)

plt.scatter(Figure_2D_data['Time/h'],
            Figure_2D_data['FW:sludge 2:3'],
            s=500, color=r, lw=5, edgecolor=dr,marker='^', zorder=1)

plt.scatter(Figure_2D_data['Time/h'],
            Figure_2D_data['FW:sludge 3:3'],
            s=500, color=r, lw=5, edgecolor=dr,marker='D', zorder=1)

#%% Figure 2E

fig, ax = plt.subplots(figsize=(13, 12))

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

ax.set_xlim(0, 132)
ax.set_ylim(0, 1.2)
plt.subplots_adjust(left=0.15, bottom=0.13)

ax.set_xlabel(r'$\mathbf{Time}$ [hr]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{Fe:P\ Molar\ Ratio}$',
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

plt.errorbar(Figure_2E_data['Time/h'],
            Figure_2E_data['FW:sludge 0:3'],
            yerr=Figure_2E_data['std'],
            lw=5, color=dp, linestyle='dotted',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2E_data['Time/h'],
            Figure_2E_data['FW:sludge 1:3'],
            yerr=Figure_2E_data['std'],
            lw=5, color=dp, linestyle='dashdot',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2E_data['Time/h'],
            Figure_2E_data['FW:sludge 2:3'],
            yerr=Figure_2E_data['std'],
            lw=5, color=dp, linestyle='dashed',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_2E_data['Time/h'],
            Figure_2E_data['FW:sludge 3:3'],
            yerr=Figure_2E_data['std'],
            lw=5, color=dp, linestyle='solid',
            capsize=12, capthick=5, zorder=0)

plt.scatter(Figure_2E_data['Time/h'],
            Figure_2E_data['FW:sludge 0:3'],
            s=500, color=p, lw=5, edgecolor=dp, marker='o', zorder=1)

plt.scatter(Figure_2E_data['Time/h'],
            Figure_2E_data['FW:sludge 1:3'],
            s=500, color=p, lw=5, edgecolor=dp,marker='s', zorder=1)

plt.scatter(Figure_2E_data['Time/h'],
            Figure_2E_data['FW:sludge 2:3'],
            s=500, color=p, lw=5, edgecolor=dp,marker='^', zorder=1)

plt.scatter(Figure_2E_data['Time/h'],
            Figure_2E_data['FW:sludge 3:3'],
            s=500, color=p, lw=5, edgecolor=dp,marker='D', zorder=1)

#%% Figure 3B P recovery

fig, ax = plt.subplots(figsize=(5, 12))

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
plt.subplots_adjust(left=0.45)

ax.set_ylabel(r'$\mathbf{P\ Recovery}$ [%]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=False, top=False, left=True, right=False,
               labelbottom=False)

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

bp = plt.boxplot(Figure_3B_P_recovery_data['P recovery'].dropna()*100,
                 whis=[5, 95], showfliers=False, widths=0.6, patch_artist=True)

for box in bp['boxes']:
    box.set(color='k', facecolor=b, linewidth=5)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=5)

for median in bp['medians']:
    median.set(color='k', linewidth=5)
    
for cap in bp['caps']:
    cap.set(color='k', linewidth=5)

#%% Figure 3B Fe recovery

fig, ax = plt.subplots(figsize=(5, 12))

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
plt.subplots_adjust(left=0.45)

ax.set_ylabel(r'$\mathbf{Fe\ Recovery}$ [%]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=False, top=False, left=True, right=False,
               labelbottom=False)

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0, 120, 20), fontname='Arial')

bp = plt.boxplot(Figure_3B_Fe_recovery_data['Fe recovery'].dropna()*100,
                 whis=[5, 95], showfliers=False, widths=0.6, patch_artist=True)

for box in bp['boxes']:
    box.set(color='k', facecolor=r, linewidth=5)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=5)

for median in bp['medians']:
    median.set(color='k', linewidth=5)
    
for cap in bp['caps']:
    cap.set(color='k', linewidth=5)

#%% Figure 3B sludge management cost

fig, ax = plt.subplots(figsize=(6, 12))

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

ax.set_ylim(200, 400)
plt.subplots_adjust(left=0.5)

ax.set_ylabel(r'$\mathbf{Sludge\ management\ cost}$'+'\n'+r'[$\$\cdot tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=False, top=False, left=True, right=False,
               labelbottom=False)

plt.yticks(np.arange(200, 440, 40), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(200, 440, 40), fontname='Arial')

bp = plt.boxplot(Figure_3B_sludge_management_cost_data['Sludge management cost'].dropna(),
                 whis=[5, 95], showfliers=False, widths=0.6, patch_artist=True)

for box in bp['boxes']:
    box.set(color='k', facecolor=o, linewidth=5)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=5)

for median in bp['medians']:
    median.set(color='k', linewidth=5)
    
for cap in bp['caps']:
    cap.set(color='k', linewidth=5)

#%% Figure 3B CI of sludge management

fig, ax = plt.subplots(figsize=(6, 12))

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

ax.set_ylim(1000, 1400)
plt.subplots_adjust(left=0.55)

ax.set_ylabel(r'$\mathbf{CI\ of\ Sludge\ management}$'+'\n'+r'[$kg\ CO_2\ eq\cdot tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=False, top=False, left=True, right=False,
               labelbottom=False)

plt.yticks(np.arange(1000, 1500, 100), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(1000, 1500, 100), fontname='Arial')

bp = plt.boxplot(Figure_3B_sludge_management_CI_data['Sludge management GWP'].dropna(),
                 whis=[5, 95], showfliers=False, widths=0.6, patch_artist=True)

for box in bp['boxes']:
    box.set(color='k', facecolor=y, linewidth=5)

for whisker in bp['whiskers']:
    whisker.set(color='k', linewidth=5)

for median in bp['medians']:
    median.set(color='k', linewidth=5)
    
for cap in bp['caps']:
    cap.set(color='k', linewidth=5)

#%% Figure 3C

fig, ax = plt.subplots(figsize=(16, 12))

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
ax.set_ylim(-10, 30)
plt.subplots_adjust(left=0.25, bottom=0.20)

ax.set_xlabel(r'$\mathbf{Avoid\ waste\ management\ cost}$'+'\n'r'[$\$\cdot tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{MSP\ of\ recovered\ FePO_4}$'+'\n'r'[$\$\cdot kg^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0, 1200, 200), fontname='Arial')
plt.yticks(np.arange(-10, 35, 5), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

x = Figure_3C_avoid_waste_sludge_management_cost_data['credit']
low = Figure_3C_avoid_waste_sludge_management_cost_data['FePO4_MSP_5th']
mid = Figure_3C_avoid_waste_sludge_management_cost_data['FePO4_MSP_50th']
high = Figure_3C_avoid_waste_sludge_management_cost_data['FePO4_MSP_95th']

ax.fill_between(x, low, high, color=p, alpha=0.85)

ax.plot(x, low, color=r,linestyle='--',linewidth=3)
ax.plot(x, high, color=r,linestyle='--',linewidth=3)

ax.plot(x, mid, color=r, linewidth=4)

#%% Figure 3D

fig, ax = plt.subplots(figsize=(16, 12))

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

ax.set_xlim(-300, 1900)
ax.set_ylim(-1200, 600)
plt.subplots_adjust(left=0.25, bottom=0.20)

ax.set_xlabel(r'$\mathbf{Avoid\ waste\ management\ CI}$'+'\n' r' [$kg\ CO_2\ eq\cdot tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{CI\ of\ recovered\ FePO_4}$'+'\n'+r'[$kg\ CO_2\ eq\cdot tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(-300, 2100, 300), fontname='Arial')
plt.yticks(np.arange(-1200, 800, 200), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

x = Figure_3D_avoid_waste_sludge_management_CI_data['credit']
low = Figure_3D_avoid_waste_sludge_management_CI_data['FePO4_GWP_5th']
mid = Figure_3D_avoid_waste_sludge_management_CI_data['FePO4_GWP_50th']
high = Figure_3D_avoid_waste_sludge_management_CI_data['FePO4_GWP_95th']

ax.fill_between(x, low, high, color=p, alpha=0.85)

ax.plot(x, low, color=r,linestyle='--',linewidth=3)
ax.plot(x, high, color=r,linestyle='--',linewidth=3)

ax.plot(x, mid, color=r, linewidth=4)

#%% Figure 4A

fig, ax = plt.subplots(figsize=(16, 12))

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
plt.subplots_adjust(bottom=0.15, left=0.15)

ax.set_xlabel(r'$\mathbf{Food\ waste:Sludge\ ratio}$',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{Fermentation\ time}$ [hr]',
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

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [dp, p, lp, 'w'][::-1])


X = np.array(Figure_4AB_data['FWsludgeratio'].map({'FW:sludge 1:3': 1,
                                                   'FW:sludge 2:3': 2,
                                                   'FW:sludge 3:3': 3}))
Y = np.array(Figure_4AB_data['HRT/h'])
Z = np.array(Figure_4AB_data['FePO4_MSP_50th'])

fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)

fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(X, Y, Z, levels=7, linewidths=5, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=45)

#%% Figure 4B

fig, ax = plt.subplots(figsize=(16, 12))

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
plt.subplots_adjust(bottom=0.15, left=0.15)

ax.set_xlabel(r'$\mathbf{Food\ waste:Sludge\ ratio}$',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{Fermentation\ time}$ [hr]',
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

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [dp, p, lp, 'w'][::-1])


X = np.array(Figure_4AB_data['FWsludgeratio'].map({'FW:sludge 1:3': 1,
                                                   'FW:sludge 2:3': 2,
                                                   'FW:sludge 3:3': 3}))
Y = np.array(Figure_4AB_data['HRT/h'])
Z = np.array(Figure_4AB_data['FePO4_GWP_50th'])

fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)

fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(X, Y, Z, levels=7, linewidths=5, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=45)

#%% Figure 4C

fig, ax = plt.subplots(figsize=(16, 12))

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
plt.subplots_adjust(bottom=0.15, left=0.15)

ax.set_xlabel(r'$\mathbf{Food\ waste:Sludge\ ratio}$',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{Fermentation\ time}$ [hr]',
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

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [do, o, lo, 'w'][::-1])


X = np.array(Figure_4CD_data['FWsludgeratio'].map({'FW:sludge 1:3': 1,
                                                   'FW:sludge 2:3': 2,
                                                   'FW:sludge 3:3': 3}))
Y = np.array(Figure_4CD_data['HRT/h'])
Z = np.array(Figure_4CD_data['sludge_cost_50th'])

fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)

fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(X, Y, Z, levels=7, linewidths=5, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=45)

#%% Figure 4D

fig, ax = plt.subplots(figsize=(16, 12))

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
plt.subplots_adjust(bottom=0.15, left=0.15)

ax.set_xlabel(r'$\mathbf{Food\ waste:Sludge\ ratio}$',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{Fermentation\ time}$ [hr]',
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

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [do, o, lo, 'w'][::-1])


X = np.array(Figure_4CD_data['FWsludgeratio'].map({'FW:sludge 1:3': 1,
                                                   'FW:sludge 2:3': 2,
                                                   'FW:sludge 3:3': 3}))
Y = np.array(Figure_4CD_data['HRT/h'])
Z = np.array(Figure_4CD_data['sludge_GWP_50th'])

fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)

fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(X, Y, Z, levels=7, linewidths=5, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=45)

#%% Figure 4E

fig, ax = plt.subplots(figsize=(18, 12))

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

ax.set_xlim(0.03, 0.09)
ax.set_ylim(0.001, 0.003)
plt.subplots_adjust(bottom=0.14, left=0.25)

ax.set_xlabel(r'$\mathbf{Landfill\ Credit}$ [$\$$ ·kg$^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{VFA\!-\!rich\ fermentation\ supernatant}$'
              '\n[$\\$$·L$^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0.03, 0.11, 0.02),  fontname='Arial')
plt.yticks(np.arange(0.001, 0.0035, 0.0005), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0.03, 0.09, 0.02), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0.001, 0.003, 0.0005), fontname='Arial')

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [do, o, lo, 'w'][::-1])


X = np.array(Figure_4EF_data['residue_price'])
Y = np.array(Figure_4EF_data['VFA_price'])
Z = np.array(Figure_4EF_data['FePO4_MSP_50th'])

fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)

fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(X, Y, Z, levels=7, linewidths=5, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=45)

#%% Figure 4F

fig, ax = plt.subplots(figsize=(18, 12))

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

ax.set_xlim(0.03, 0.09)
ax.set_ylim(0.001, 0.003)
plt.subplots_adjust(bottom=0.14, left=0.25)

ax.set_xlabel(r'$\mathbf{Landfill\ Credit}$ [$\$$ ·kg$^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.set_ylabel(r'$\mathbf{VFA\!-\!rich\ fermentation\ supernatant}$'
              '\n[$\\$$·L$^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(0.03, 0.11, 0.02),  fontname='Arial')
plt.yticks(np.arange(0.001, 0.0035, 0.0005), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

plt.xticks(np.arange(0.03, 0.09, 0.02), fontname='Arial')

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

plt.yticks(np.arange(0.001, 0.003, 0.0005), fontname='Arial')

color_map_Guest = colors.LinearSegmentedColormap.from_list('color_map_Guest', [do, o, lo, 'w'][::-1])


X = np.array(Figure_4EF_data['residue_price'])
Y = np.array(Figure_4EF_data['VFA_price'])
Z = np.array(Figure_4EF_data['FePO4_GWP_50th'])

fills = ax.tricontourf(X, Y, Z, levels=10000, cmap=color_map_Guest)

fig.colorbar(fills, ax=ax)

fig.delaxes(fig.axes[3])

lines = ax.tricontour(X, Y, Z, levels=2, linewidths=5, linestyles='solid', colors='k')

ax.clabel(lines, lines.levels, inline=True, fontsize=45)
#%% SI- NH4+

fig, ax = plt.subplots(figsize=(14, 12))

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

ax.set_xlim(0, 132)
ax.set_ylim(0, 100)
plt.subplots_adjust(left=0.20, bottom=0.13)

ax.set_xlabel(r'$\mathbf{Time}$ [hr]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{NH_4^{+}\!-N}$ [mg ·L$^{-1}$]',
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

plt.errorbar(Figure_S1_NH4_release_data['Time/h'],
            Figure_S1_NH4_release_data['FW:sludge 0:3'],
            yerr=Figure_S1_NH4_release_data['std'],
            lw=5, color=dy, linestyle='dotted',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_S1_NH4_release_data['Time/h'],
            Figure_S1_NH4_release_data['FW:sludge 1:3'],
            yerr=Figure_S1_NH4_release_data['std'],
            lw=5, color=dy, linestyle='dashdot',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_S1_NH4_release_data['Time/h'],
            Figure_S1_NH4_release_data['FW:sludge 2:3'],
            yerr=Figure_S1_NH4_release_data['std'],
            lw=5, color=dy, linestyle='dashed',
            capsize=12, capthick=5, zorder=0)

plt.errorbar(Figure_S1_NH4_release_data['Time/h'],
            Figure_S1_NH4_release_data['FW:sludge 3:3'],
            yerr=Figure_S1_NH4_release_data['std'],
            lw=5, color=dy, linestyle='solid',
            capsize=12, capthick=5, zorder=0)

plt.scatter(Figure_S1_NH4_release_data['Time/h'],
            Figure_S1_NH4_release_data['FW:sludge 0:3'],
            s=500, color=y, lw=5, edgecolor=dy, marker='o', zorder=1)

plt.scatter(Figure_S1_NH4_release_data['Time/h'],
            Figure_S1_NH4_release_data['FW:sludge 1:3'],
            s=500, color=y, lw=5, edgecolor=dy,marker='s', zorder=1)

plt.scatter(Figure_S1_NH4_release_data['Time/h'],
            Figure_S1_NH4_release_data['FW:sludge 2:3'],
            s=500, color=y, lw=5, edgecolor=dy,marker='^', zorder=1)

plt.scatter(Figure_S1_NH4_release_data['Time/h'],
            Figure_S1_NH4_release_data['FW:sludge 3:3'],
            s=500, color=y, lw=5, edgecolor=dy,marker='D', zorder=1)

leg = ax.legend(loc='upper left',
                bbox_to_anchor=(0.02, 0.98),
                frameon=False,
                fontsize=28,
                handlelength=2.5,
                handletextpad=0.5,
                labelspacing=0.4)

for text in leg.get_texts():
    text.set_fontname('Arial')

#%% SI-IRR-A

fig, ax = plt.subplots(figsize=(16, 12))

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

ax.set_xlim(-300, 1900)
ax.set_ylim(-1200, 600)
plt.subplots_adjust(left=0.25, bottom=0.20)

ax.set_xlabel(r'$\mathbf{Avoid\ waste\ management\ CI}$'+'\n' r' [$kg\ CO_2\ eq\cdot tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{CI\ of\ recovered\ FePO_4}$'+'\n'+r'[$kg\ CO_2\ eq\cdot tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(-300, 2100, 300), fontname='Arial')
plt.yticks(np.arange(-1200, 800, 200), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

x = Figure_3D_avoid_waste_sludge_management_CI_data['credit']
low = Figure_3D_avoid_waste_sludge_management_CI_data['FePO4_GWP_5th']
mid = Figure_3D_avoid_waste_sludge_management_CI_data['FePO4_GWP_50th']
high = Figure_3D_avoid_waste_sludge_management_CI_data['FePO4_GWP_95th']

ax.fill_between(x, low, high, color=p, alpha=0.85)

ax.plot(x, low, color=r,linestyle='--',linewidth=3)
ax.plot(x, high, color=r,linestyle='--',linewidth=3)

ax.plot(x, mid, color=r, linewidth=4)

#%% SI-IRR-B

fig, ax = plt.subplots(figsize=(16, 12))

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

ax.set_xlim(-300, 1900)
ax.set_ylim(-1200, 600)
plt.subplots_adjust(left=0.25, bottom=0.20)

ax.set_xlabel(r'$\mathbf{Avoid\ waste\ management\ CI}$'+'\n' r' [$kg\ CO_2\ eq\cdot tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=0)

ax.set_ylabel(r'$\mathbf{CI\ of\ recovered\ FePO_4}$'+'\n'+r'[$kg\ CO_2\ eq\cdot tonne^{-1}$]',
              fontname='Arial',
              fontsize=45,
              labelpad=5)

ax.tick_params(direction='inout', length=30, width=5,
               bottom=True, top=False, left=True, right=False)

plt.xticks(np.arange(-300, 2100, 300), fontname='Arial')
plt.yticks(np.arange(-1200, 800, 200), fontname='Arial')

ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())

ax_top.tick_params(direction='in', length=15, width=5,
                   bottom=False, top=True, left=False, right=False,
                   labeltop=False)

ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())

ax_right.tick_params(direction='in', length=15, width=5,
                     bottom=False, top=False, left=False, right=True,
                     labelright=False)

x = Figure_3D_avoid_waste_sludge_management_CI_data['credit']
low = Figure_3D_avoid_waste_sludge_management_CI_data['FePO4_GWP_5th']
mid = Figure_3D_avoid_waste_sludge_management_CI_data['FePO4_GWP_50th']
high = Figure_3D_avoid_waste_sludge_management_CI_data['FePO4_GWP_95th']

ax.fill_between(x, low, high, color=p, alpha=0.85)

ax.plot(x, low, color=r,linestyle='--',linewidth=3)
ax.plot(x, high, color=r,linestyle='--',linewidth=3)

ax.plot(x, mid, color=r, linewidth=4)

#%% legend

fig, ax = plt.subplots(figsize=(8,4))
ax.axis('off')

ax.text(0,1.05, 'Food Waste (FW): Sludge ratio',
        transform=ax.transAxes,
        fontsize=30,
        fontweight='bold',
        va='bottom',
        ha='left')
handles = [

    plt.Line2D([0],[0],
               color='k',
               lw=6,
               linestyle='dotted',
               marker='o',
               markersize=24,
               label='0:3'),

    plt.Line2D([0],[0],
               color='k',
               lw=6,
               linestyle='dashdot',
               marker='s',
               markersize=24,
               label='1:3'),

    plt.Line2D([0],[0],
               color='k',
               lw=6,
               linestyle='dashed',
               dashes=(8,6),
               marker='^',
               markersize=24,
               label='2:3'),

    plt.Line2D([0],[0],
               color='k',
               lw=6,
               linestyle='solid',
               marker='D',
               markersize=24,
               label='3:3'),
]

leg = ax.legend(handles=handles,
                loc='center',
                bbox_to_anchor=(0.52, 0.40),
                frameon=False,
                fontsize=28,
                handlelength=10,
                labelspacing=1.2)

for t in leg.get_texts():
    t.set_fontname('Arial')
    t.set_fontweight('bold')
    t.set_fontsize(36)

plt.tight_layout()
plt.show()