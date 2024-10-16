# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 13:47:07 2024

@author: joy_c
"""

import pandas as pd, numpy as np, matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as itp
from qsdsan.utils import load_data, ospath

mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams["figure.autolayout"] = False
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

folder = ospath.dirname(__file__)

#%%

# var = 'GHG'
def plot_lines(save=True, percentiles=[5, 25, 50, 75, 95], 
               var='Operational'):
    data = {}
    for conf in ('B1', 'metro',):
        path = ospath.join(folder, f'results/{conf}.xlsx')
        y = load_data(path, sheet='Pairwise '+var)
        x = np.float64(y.columns)
        xhat = np.arange(x.min(), x.max(), 0.1)
        f = itp.interp1d(x, y)
        yhat = f(xhat)
        data[conf] = {'x':xhat, 'y':np.percentile(yhat, percentiles, axis=0)}
        
    fig, ax = plt.subplots(figsize=(5, 5))
    # colors = ('#35767f', '#a75f3e')
    # shades = ('#60c1cf', '#f98f60')
    colors = ('#4c385a', '#4d7e53')
    shades = ('#a280b9', '#79bf82')
    confs = ''
    for (conf, xy), c, s in zip(data.items(), colors, shades):
        x, y = xy.values()
        confs += '_'+conf
        ax.axhline(y=0, color='black', lw=1)

        ax.fill_between(x, y[0], y[-1], fc=s, alpha=0.25)
        ax.fill_between(x, y[1], y[-2], fc=s, alpha=0.25)
        ax.plot(x, y[2], color=c, ls='-', lw=2)     
        ax.plot(x, y[[1, 3]].T, color=c, ls='--', lw=1.5)

    ax.set_xlim(x.min(), x.max())
    ax.tick_params(axis='y', which='major', direction='inout', length=10, labelsize=14)
    ax.tick_params(axis='y', which='minor', direction='inout', length=6)
    
    ax.tick_params(axis='x', which='major', direction='inout', length=10, labelsize=14)
    ax.tick_params(axis='x', which='minor', direction='inout', length=6)
    
    ax2y = ax.secondary_yaxis('right')
    ax2y.tick_params(axis='y', which='major', direction='in', length=5)
    ax2y.tick_params(axis='y', which='minor', direction='in', length=3)
    ax2y.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax2x = ax.secondary_xaxis('top')
    ax2x.tick_params(axis='x', which='major', direction='in', length=5)
    ax2x.tick_params(axis='x', which='minor', direction='in', length=3)
    ax2x.xaxis.set_major_formatter(plt.NullFormatter())

    if save:
        file = ospath.join(folder, f'figures/{var}_delta{confs}.png')
        fig.savefig(file, dpi=300, transparent=True)
    else:
        return fig, ax

#%%
colors = {
    'Secondary treatment': '#60c1cf',
    'Aeration': '#60c1cf',
    'Effluent discharge': '#35767f',
    'Pumping': '#f98f60',
    'Sludge disposal': '#a75f3e',
    'Electricity': '#90918e',
    }

hatch = {
    'CH4': '////',
    'N2O': '---',
    'CO2': ''
    }


# var = 'Cost'
def plot_breakdown(save=True, var='GHG_eb'):
    path = ospath.join(folder, 'results/Baseline_scenario_analysis.xlsx')
    fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(4,5))
    plt.subplots_adjust(left=0.2, wspace=0)
    if var == 'GHG_eb': header = [0,1]
    else: header=0
    for conf, ax in zip(('B1', 'metro'), axes):
        df = load_data(path, sheet=f'{var}_{conf}', index_col=0, header=header)
        y = df.iloc[-2:,]
        x = [0.28, 0.72]
        y0 = np.zeros(y.shape[0])
        for k, v in y.items():
            if isinstance(k, tuple):
                c = colors[k[0]]
                hat = hatch[k[1]]
            else: 
                c = colors[k]
                hat = ''
            ax.bar(x, v, width=0.25, bottom=y0, color=c, hatch=hat, 
                   ec='black', lw=0.4)
            y0 += v
        ax.set_xlim(0, 1)
        ax.tick_params(axis='y', which='major', direction='inout', length=10, labelsize=14)
        ax.tick_params(axis='y', which='minor', direction='inout', length=6)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(axis='y', which='major', direction='in', length=5)
        ax2y.tick_params(axis='y', which='minor', direction='in', length=3)
        ax2y.yaxis.set_major_formatter(plt.NullFormatter())
        ax.tick_params(bottom=False, which='both', labelbottom=False)
    if save:
        fig.savefig(ospath.join(folder, f'figures/{var}_breakdown'), 
                    dpi=300, transparent=True)
    else:
        return fig, axes

#%%
if __name__ == '__main__':
    # plot_lines()
    # plot_lines(var='GHG')
    plot_breakdown()
    plot_breakdown(var='Cost')
