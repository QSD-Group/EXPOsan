# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import pandas as pd, numpy as np, matplotlib.pyplot as plt
import matplotlib 
from qsdsan.utils import load_data, ospath, palettes
from exposan.werf import data_path, results_path, figures_path

# import warnings
# warnings.filterwarnings("ignore")

Guest = palettes['Guest']

b = Guest.blue.HEX
g = Guest.green.HEX
r = Guest.red.HEX
o = Guest.orange.HEX
y = Guest.yellow.HEX
p = Guest.purple.HEX
a = Guest.gray.HEX
da = Guest.dgray.HEX
db = Guest.dblue.HEX
# %%
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams["figure.autolayout"] = False
plt.rcParams['xtick.minor.visible'] = False
plt.rcParams['ytick.minor.visible'] = True

thresholds = {
    'BOD': [1,18,10], # xmin, xmax, y
    'TSS': [1,18,15],
    'NH4 N': [7,18,2],
    'TN': [10,18,10],
    'TP': [10,18,2],
    }

ylabels = {
    'COD': '$\mathbf{COD}$ [mg·L$^{-1}$]', 
    'BOD': '$\mathbf{BOD}$ [mg·L$^{-1}$]', 
    'TSS': '$\mathbf{TSS}$ [mg·L$^{-1}$]', 
    'TN': '$\mathbf{TN}$ [mg·L$^{-1}$]', 
    'NH4 N': '$\mathbf{NH_4-N}$ [mg·L$^{-1}$]', 
    'TP': '$\mathbf{TP}$ [mg·L$^{-1}$]', 
    'Ortho P': '$\mathbf{Ortho-P}$ [mg·L$^{-1}$]', 
    'SRT': '$\mathbf{SRT}$ [d]', 
    'CH4 production': '$\mathbf{CH_4\ production\ rate}$ [tonne·d$^{-1}$]', 
    'CH4 content': '$\mathbf{CH_4\ content}$ [%]', 
    'Liquid aeration flowrate': '$\mathbf{Aeration\ flowrate}$ [10$^5$ m$^3$·d$^{-1}$]', 
    'Liquid aeration energy':'$\mathbf{Aeration\ energy}$ [MWh·d$^{-1}$]'
    }
# %%
data = load_data(
    ospath.join(results_path, 'baseline_unopt_performance.xlsx'), sheet=0,
    header=[0,1], index_col=0,
    )
data.columns = [b.split(' [')[0] for a,b in data.columns]
data['Liquid aeration flowrate'] /= 1e5
data['Liquid aeration energy'] *= 24/1e3
data['CH4 production'] *=24/1e3

gpsx = load_data(
    ospath.join(data_path, 'GPS-X_results.xlsx'), sheet='All',
    header=[0], index_col=0,
    )
gpsx['Liquid aeration flowrate'] /= 1e5
gpsx['Liquid aeration energy'] *= 24/1e3
gpsx['CH4 production'] *=24/1e3

# %%

def plot():
    for var in gpsx.columns:
        if var.startswith('CH4'): continue
        fig, ax = plt.subplots(figsize=(6.5,3.5))
        if var in thresholds:
            xmin, xmax, y = thresholds[var]
            ax.hlines(y, xmin-0.5, xmax+0.5, colors=r, linestyles='--', zorder=0)
            if var == 'TP':
                ax.text(xmax+0.5, y, 'treatment\n', ha='right', va='bottom', 
                        fontsize=11, fontstyle='italic', color='black', 
                        linespacing=0.2, zorder=0)
                ax.text(xmax+0.5, y, '\ntarget', ha='right', va='top', 
                        fontsize=11, fontstyle='italic', color='black', 
                        linespacing=0.2, zorder=0)
            else:
                ax.text(xmax+0.5, y, 'treatment target\n', ha='right', va='bottom', 
                        fontsize=11, fontstyle='italic', color='black', 
                        linespacing=0.2, zorder=0)
        x = np.arange(data.shape[0]) + 1
        if var in gpsx:
            ax.plot(x, gpsx[var], marker='D', ls='', 
                    mec='black', mfc='white', ms=3.5, zorder=2)
        ax.bar(x, data[var], width=0.4, color=db, zorder=1)
        ax.set_xticks(x, data.index.values)     
        ax.set_xlim((0, len(x)+1))
        ax.set_xlabel('WRRF configuration', fontname='Arial', weight='bold', fontsize=12)
        ax.set_ylabel(ylabels[var], fontname='Arial', fontsize=12)
        ax.tick_params(axis='y', which='major', direction='inout', length=8, labelsize=11)
        ax.tick_params(axis='y', which='minor', direction='inout', length=5)
        ax.tick_params(axis='x', which='major', direction='out', length=4, labelsize=10.5)
        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(axis='y', which='major', direction='in', length=4)
        ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
        ax2y.yaxis.set_major_formatter(plt.NullFormatter())
        fig.subplots_adjust(bottom=0.15)
        fig.savefig(ospath.join(figures_path, f'baseline_vs_gpsx/{var.replace(" ", "_")}.tif'),
                    dpi=300, transparent=True)

# %%
    fig, ax = plt.subplots(figsize=(6.5,3.5))


    yy = data[['CH4 production', 'CH4 content']].dropna()
    xx = np.arange(len(yy)) + 1
    ax.plot(xx-0.15, gpsx['CH4 production'][yy.index], marker='D', ls='', 
            mec=db, mfc='white', ms=4, zorder=2)
    ax.bar(xx-0.15, yy['CH4 production'], width=0.3, color=db, zorder=1)
    ax.set_xticks(xx, yy.index.values)
    ax.set_xlim((0.3, len(xx)+0.7))
    ax.set_xlabel('WRRF configuration', fontname='Arial', weight='bold', fontsize=12)
    ax.set_ylabel(ylabels['CH4 production'], fontname='Arial', fontsize=12, color=db)
    ax.tick_params(axis='y', which='major', direction='inout', length=8, 
                   labelsize=11, labelcolor=db)
    ax.tick_params(axis='y', which='minor', direction='inout', length=5)
    ax.tick_params(axis='x', which='major', direction='out', length=4, labelsize=10.5)

    ax2 = ax.twinx()
    ax2.plot(xx+0.15, gpsx['CH4 content'][yy.index], marker='D', ls='', 
            mec=o, mfc='white', ms=4, zorder=2)
    ax2.bar(xx+0.15, yy['CH4 content'], width=0.3, color=o, zorder=1)
    ax2.set_ylabel(ylabels['CH4 content'], fontname='Arial', fontsize=12, rotation=-90,
                   rotation_mode='anchor', va='bottom', color=o)
    ax2.tick_params(axis='y', which='major', direction='inout', length=8, 
                    labelsize=11, labelcolor=o)
    ax2.tick_params(axis='y', which='minor', direction='inout', length=5)
    
    fig.subplots_adjust(bottom=0.15)
    fig.savefig(ospath.join(figures_path, 'baseline_vs_gpsx/Biogas.tif'),
                dpi=300, transparent=True)

# %%
if __name__ == '__main__':
    plot()

