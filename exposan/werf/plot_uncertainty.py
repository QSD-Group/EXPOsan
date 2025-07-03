# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

import pandas as pd, numpy as np, matplotlib.pyplot as plt
from qsdsan.utils import load_data, ospath, palettes
from exposan.werf import results_path, figures_path

import warnings
warnings.filterwarnings("ignore")

Guest = palettes['Guest']

b = Guest.blue.HEX
g = Guest.green.HEX
r = Guest.red.HEX
o = Guest.orange.HEX
y = Guest.yellow.HEX
p = Guest.purple.HEX
a = Guest.gray.HEX
db = Guest.dblue.HEX
dg = Guest.dgreen.HEX
dy = Guest.dyellow.HEX

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams["figure.autolayout"] = False
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True

# %%
path = ospath.join(results_path, 'HA_F1_UA_opex_stats.xlsx')
data = load_data(path, sheet=None)
bl = []
for strength in ('low', 'mid', 'high'):
    bl.append(
        load_data(
            ospath.join(results_path, f'HA_F1_UA-{strength}.xlsx'), 
            header=[0,1], skiprows=[2,], sheet='no_HA'
            )[('OPEX', 'Total OPEX [USD/d]')]
        )
bl = pd.DataFrame(bl, index=['low', 'mid', 'high']).T

colors = {
    'low': (b, db),
    'mid': (g, dg),
    'high': (y, dy)
    }

# %%
def plot_linebox():
    fig, (ax, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(8, 6), width_ratios=[3,1])
    pos = 1
    for strength, quantiles in data.items():    
        x = [float(i)*100 for i in quantiles.index[1:]]
        color, dcolor = colors[strength]
        ax.fill_between(
            x, y1=quantiles.loc['0.50':, 0.05], y2=quantiles.loc['0.50':, 0.95], 
            alpha=0.3, color=color, lw=0, 
            )
        ax.fill_between(
            x, y1=quantiles.loc['0.50':, 0.25], y2=quantiles.loc['0.50':, 0.75], 
            alpha=0.5, color=color, lw=0,
            )
        ax.plot(x, quantiles.loc['0.50':, 0.5], color=dcolor)
        ax.plot(x, quantiles.loc['0.50':, [0.25, 0.75]], color=color, lw=1.5, ls='--')

        lineprops = dict(color=dcolor)
        bplot = ax1.boxplot(
            bl[strength], whis=(5, 95), 
            positions=[pos,], tick_labels=[strength,], widths=0.5,
            showfliers=False, patch_artist=True,
            whiskerprops=lineprops, capprops=lineprops, medianprops=lineprops
            )
        bplot['boxes'][0].set(facecolor=color, edgecolor=dcolor)
        pos += 1
        
    ax.set_xlim((x[0], x[-1]))
    xmin, xmax = ax1.get_xlim()
    ax1.set_xlim((xmin-0.25, xmax))
    ax.set_ylim((0, 1e4))
    
    major = dict(which='major', length=8, labelsize=12)
    minor = dict(which='minor', length=5)
    ax.tick_params(axis='both', direction='inout', **major)
    ax.tick_params(axis='both', direction='inout', **minor)
    ax1.tick_params(axis='both', direction='inout' **major)
    ax1.tick_params(axis='y', direction='inout', **minor)
    ax1.tick_params(axis='x', which='minor', length=0)
    
    ax2y = ax1.secondary_yaxis('right')
    ax2y.tick_params(axis='y', which='major', direction='in', length=4)
    ax2y.tick_params(axis='y', which='minor', direction='in', length=2.5)
    ax2y.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax.set_title('Configuration F1 + HA', fontsize=14, fontweight='bold')
    ax1.set_title('F1', fontsize=14, fontweight='bold')
    ax.set_xlabel('$\mathbf{HA\ NH_4\ recovery\ efficiency}$ [%]', fontname='Arial', fontsize=13)
    ax.set_ylabel('$\mathbf{OPEX}$ [USD·day${^{-1}}$]', fontname='Arial', fontsize=13)
    ax1.set_xlabel('Wastewater \n strength', fontname='Arial', fontweight='bold', fontsize=13)
    
    fig.subplots_adjust(wspace=0, bottom=0.25)
    
    fig.savefig(ospath.join(figures_path, 'HA_F1_UA.png'),
                dpi=300, transparent=True)
    
# %%
def compile_dopex(ID='F1', save_as=''):
    zz = {}
    xx = None
    for strength in ('low', 'mid', 'high'):
        dfs = load_data(ospath.join(results_path, f'HA_{ID}_UA-{strength}.xlsx'), 
                        header=[0,1], skiprows=[2,], sheet=None)
        dopex = []
        keys = []
        for k, df in dfs.items():
            if k != 'no_HA':
                dopex.append(bl[strength] - df[('OPEX', 'Total OPEX [USD/d]')])
                keys.append(k)
            else:
                if xx is None: 
                    xx = df[('System', 'Electricity price [cents/kWh]')]
                    yy = df[('System', 'Sludge disposal price [USD/wet tonne]')]
        dopex = pd.DataFrame(dopex, index=keys).T
        zz[strength] = dopex
    if save_as:
        with pd.ExcelWriter(save_as) as writer:
            for k, v in zz.items():
                v.to_excel(writer, sheet_name=k)
    return xx, yy, zz

# %%
# strength = 'low'
# fig, axes = plt.subplots(nrows=3, sharex=True, figsize=(4, 12))
# for strg, ax in zip(('low', 'mid', 'high'), axes):
#     zs = zz[strg]
    

# %%
    
if __name__ == '__main__':
    # plot_linebox()
    dopex = compile_dopex(save_as='HA_F1_dopex.xlsx')
