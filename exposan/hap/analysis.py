# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import qsdsan as qs, numpy as np, pandas as pd
from exposan.hap import (
    results_path,
    figures_path
    )
from qsdsan.utils import load_data, ospath
import matplotlib as mpl, matplotlib.pyplot as plt, seaborn as sns

mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams["figure.autolayout"] = False
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True

#%%
boxprops = dict(alpha=1, facecolor='white', edgecolor='black')
meanprops = dict(marker='x', markersize=7, markeredgecolor='black')
flierprops = dict(marker='.', markersize=2, markerfacecolor='#4e4e4e')
medianprops = dict(color='black')

# _market_prices={            # USD/kg
#     'low range': [136,      # puriss https://www.sigmaaldrich.com/US/en/product/sial/04238
#                   194],     # purum  https://www.sigmaaldrich.com/US/en/product/sial/21223
#     'high range': [4000,    # https://www.sigmaaldrich.com/US/en/product/aldrich/900203
#                    9600]    # https://www.sigmaaldrich.com/US/en/product/aldrich/677418
#     }

def plot_mpsp(data=None, seed=None, save=True, 
              market_ranges=[[136, 194], [4000, 9600]],
              single_values=[50, 170]):
    if data is None:
        data = load_data(ospath.join(results_path, f'table_{seed}.xlsx'),
                          header=[0,1], skiprows=[2,])
    mpsp = data.loc[:,  ('TEA', 'MPSP [USD/kg]')]
    # mpsp = np.random.rand(100) * 100
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(2,6), height_ratios=[1,2])
    fig.subplots_adjust(left=0.3, hspace=0.05)
    for ll, ul in market_ranges:
        ax1.axhspan(ll, ul, alpha=0.3, fill=True, facecolor='#4e4e4e')
    for v in single_values:
        ax1.axhline(v, color='black', linewidth=0.5)
        ax2.axhline(v, color='black', linewidth=0.5)
    ax1.set_ylim(91, 11000)
    ax1.set_yscale('log')
    ax2.boxplot(mpsp,
                whis=(5,95), 
                widths=0.6,
                patch_artist=True,
                showmeans=True,
                boxprops=boxprops,
                meanprops=meanprops,
                flierprops=flierprops,
                medianprops=medianprops
                )
    ax2.set_ylim(0, 55)
    for ax in (ax1, ax2):
        ax.tick_params(axis='y', which='major', direction='inout', length=10, 
                       labelsize=13)
        ax.tick_params(axis='y', which='minor', direction='inout', length=6)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax2y = ax.secondary_yaxis('right')
        ax2y.tick_params(axis='y', which='major', direction='in', length=5)
        ax2y.tick_params(axis='y', which='minor', direction='in', length=3)
        ax2y.yaxis.set_major_formatter(plt.NullFormatter())
        ax.tick_params(bottom=False, which='both', labelbottom=False)
    ax1.spines.bottom.set_visible(False)
    ax2.spines.top.set_visible(False)
    ax1.tick_params(labeltop=False)
    
    d = 0.35  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=10,
                  linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
    if save:
        fig.savefig(ospath.join(figures_path, 'msps.png'), 
                    dpi=300, transparent=True)
    else:
        return fig, ax

#%%

patch_dct = {
    'CAPEX': ('white', r'\\\\'),
    'OPEX': ('white', ''),
    'HAp fermenter': ('#f98f60', ''),
    'Yeast production': ('#60c1cf', ''),
    'Transportation': ('#f3c354', ''),
    'Post processing': ('#a280b9', ''),
    'Other': ('#90918e', ''),
    }

def plot_area(df, figsize=(1.5, 6)):
    # df.sort_values(by=df.columns[0], inplace=True)
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=0.3)
    x = range(df.shape[0])
    yp = np.zeros(df.shape[0])
    yn = yp.copy()
    for k, v in patch_dct.items():
        if k in df.columns:
            c, hat = v
            y = df.loc[:,k]
            y_offset = (y>=0)*yp + (y<0)*yn
            ax.fill_between(x, y+y_offset, y_offset, facecolor=c, hatch=hat, 
                            zorder=0) 
                            # linewidth=0.25, ec='black')
            yp += (y>=0) * y
            yn += (y<0) * y
    ax.set_xlim(0, max(x))
    ax.set_xticks([])
    ax.set_ylim(0,100)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(labelsize=13)
    ax.tick_params(axis='y', which='major', direction='inout', length=8)
    ax.tick_params(axis='y', which='minor', direction='inout', length=4)
    ax2y = ax.secondary_yaxis('right')
    ax2y.set_yticks(ax.get_yticks())
    ax2y.tick_params(axis='y', which='major', direction='in', length=4)
    ax2y.tick_params(axis='y', which='minor', direction='in', length=2)
    ax2y.yaxis.set_ticklabels([])
    return fig, ax

def plot_npv_breakdown(data=None, seed=None):
    if data is None:
        data = load_data(ospath.join(results_path, f'table_{seed}.xlsx'),
                          header=[0,1], skiprows=[2,])
    # data.sort_values(by=('TEA', 'MPSP [USD/kg]'), inplace=True)
    npv = data.TEA[['CAPEX [% ANPV]', 'OPEX [% ANPV]']]
    for j, df in enumerate([npv, data.CAPEX, data.OPEX]):
        df.columns = [i.split(' [')[0] for i in df.columns]
        figsize= (1.5, 6)
        if j == 1: 
            df.columns = [i.rstrip(' CAPEX') for i in df.columns]
            figsize= (2, 1.8)
        if j == 2: 
            df.columns = [i.rstrip(' OPEX') for i in df.columns]
            figsize= (2, 3.9)
        by = df.columns[0]
        fig, ax = plot_area(df.sort_values(by), figsize)
        fig.savefig(ospath.join(figures_path, f'breakdown_{j}.png'), 
                    dpi=300, transparent=True)

color_dct = {
    'HAp fermenter': '#f98f60',
    'Yeast production': '#60c1cf',
    'Transportation': '#f3c354',
    'Post processing': '#a280b9',
    'Other': '#90918e',
    }

hatch_dct = {
    'CAPEX': r'\\\\',
    'OPEX': '',
    }

def plot_unified_breakdown(data=None, seed=None):
    if data is None:
        data = load_data(ospath.join(results_path, f'table_{seed}.xlsx'),
                          header=[0,1], skiprows=[2,])
    data.sort_values(by=('TEA', 'CAPEX [% ANPV]'), inplace=True)
    npv = data.TEA[['CAPEX [% ANPV]', 'OPEX [% ANPV]']]
    npv.columns = [i.split(' [')[0] for i in npv.columns]
    df = pd.concat([data.CAPEX.mul(npv.CAPEX/100, axis=0), 
                    data.OPEX.mul(npv.OPEX/100, axis=0)],
                   axis=1)
    df.columns = [i.split(' [')[0] for i in df.columns]
    cols = []
    for i in df.columns:
        col = i.split()
        cat = col.pop()
        cols.append((cat, ' '.join(col)))
    df.columns = cols
    fig, ax = plt.subplots(figsize=(6,6))
    x = range(df.shape[0])
    yp = np.zeros(df.shape[0])
    yn = yp.copy()
    for k, y in df.items():
        hat, c = k
        c = color_dct[c]
        hat = hatch_dct[hat]
        y_offset = (y>=0)*yp + (y<0)*yn
        ax.fill_between(x, y+y_offset, y_offset, facecolor=c, hatch=hat,
                        zorder=0)
                        # linewidth=0.2, 
                        # ec='black')
        yp += (y>=0) * y
        yn += (y<0) * y
    ax.set_xlim(min(x), max(x))
    ax.set_xticks([])
    ax.set_ylim(0,100)
    ax.tick_params(labelsize=13)
    ax.tick_params(axis='y', which='major', direction='inout', length=8)
    ax.tick_params(axis='y', which='minor', direction='inout', length=4)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax2y = ax.secondary_yaxis('right')
    ax2y.tick_params(axis='y', which='major', direction='in', length=4)
    ax2y.tick_params(axis='y', which='minor', direction='in', length=2)
    ax2y.yaxis.set_ticklabels([])
    fig.savefig(ospath.join(figures_path, 'breakdown_all.png'), 
                dpi=300, transparent=True)

 
#%%
def plot_sensitivity(seed=None, save=True):
    data = load_data(ospath.join(results_path, f'Spearman_{seed}.xlsx'),
                     header=[0,1], skiprows=[2,], index_col=[0,1], sheet=None)
    rho, p = data.values()
    col = ('TEA', 'MPSP [USD/kg]')
    r = rho.loc[:,col]
    sig = ['***' if i < 0.001 else
           '**' if i < 0.01 else
           '*' if i < 0.05 else 
           ''  for i in p.loc[:,col]]
    # colors = ['#4d7e53' if i<0 else '#9c4b50' for i in r]
    colors = ['#79bf82' if i<0 else '#ed586f' for i in r]    
    fig, ax = plt.subplots(figsize=(12, 6))
    y_pos = list(range(len(r)))
    ax.barh(y_pos, r, height=0.7, linewidth=0.7,
            color=colors, edgecolor='black')
    ax.axvline(x=0, color='black', linewidth=0.8)
    for s, x, y in zip(sig, r, y_pos):
        off = -0.5*len(s)-0.1 if x < 0 else 0.5*len(s)-1
        ax.annotate(s, (x,y), (off, -0.55), textcoords='offset fontsize',
                    fontsize='large')
    ax.set_xlim(-1, 1)
    ax.invert_yaxis()
    ax.tick_params(axis='x', which='major', direction='inout', length=12, labelsize=14)
    ax.tick_params(axis='x', which='minor', direction='inout', length=7)
    ax.set_xlabel('')
    ax.set_yticks([])
    ax2x = ax.secondary_xaxis('top')
    ax2x.tick_params(axis='x', which='major', direction='in', length=6)
    ax2x.tick_params(axis='x', which='minor', direction='in', length=3.5)
    ax2x.xaxis.set_major_formatter(plt.NullFormatter())
    ax.tick_params(right=False, which='both', labelright=False)
    if save:
        fig.savefig(ospath.join(figures_path, 'spearman.png'), 
                    dpi=300, transparent=True)
    else: return fig, ax

#%%
if __name__ == '__main__':
    seed = 292
    plot_mpsp(seed=seed)
    plot_npv_breakdown(seed=seed)
    plot_unified_breakdown(seed=seed)
    plot_sensitivity(seed)
