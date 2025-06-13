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
from qsdsan.utils import load_data, ospath
from exposan.werf import results_path, figures_path

plt.rcParams['font.sans-serif'] = 'Arial'

#%%

# # diff = load_data(ospath.join(results_path, 'performance.xlsx'), sheet='diff').T

thresholds = {
    'BOD': [10]*18,
    'NH4 N': [40]*6 + [2]*12,
    'TN': [40]*9 + [10]*9,
    'TP': [7]*9 + [2]*9,
    }
thresholds = pd.DataFrame.from_dict(
    thresholds, orient='index', 
    columns=('B1', 'B2', 'B3', 'C1', 'C2', 'C3', 
            'E2', 'E2P', 'F1', 
            'G1', 'G2', 'G3', 'H1', 'I1', 'I2', 'I3', 'N1', 'N2'))

eff_vars = ['COD', 'BOD', 'TN', 'NH4 N', 'TP']
op_vars = ['CH4 production', 'CH4 content', 'Sludge production', 
           'Liquid aeration flowrate', 'Sludge aeration flowrate']

#%%
def heatmap(data, colorbar=True, cmap='viridis', 
            show_ticklabels=True, row_labels=[], col_labels=[],
            annotate=None, valfmt="{x:.1f}",
            txtcolors=None, annotate_kw={}, 
            save_as=None, **kwargs):
    fig, ax = plt.subplots(figsize=(12,6))
    im = ax.imshow(data, cmap=cmap, **kwargs)
    if colorbar:
        cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, fraction=0.05)
        if cmap == 'bwr':
            cbar.ax.set_xticks(ticks=[-1,0,1], labels=['-100%', '0%', '+100%'], fontsize=12)
        else:
            cbar.ax.set_xticks(ticks=[0,1], labels=['min', 'max'], fontsize=12)

    # Show all ticks and label them with the respective list entries.
    if show_ticklabels:
        row_labels = row_labels or data.index
        col_labels = col_labels or data.columns
    else:
        row_labels = col_labels = []
        
    ax.set_xticks(
        range(data.shape[1]), 
        labels=col_labels, fontsize=12,
        # rotation=0, ha="right", rotation_mode="anchor"
        )
    ax.set_yticks(
        range(data.shape[0]), 
        labels=row_labels, fontsize=12,
        )

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)   

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    if annotate is not None:
        kw = dict(horizontalalignment="center",
                  verticalalignment="center",
                  size=11)
        kw.update(annotate_kw)
        
        if isinstance(valfmt, str):
            valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)
        
        for i in range(annotate.shape[0]):
            row = annotate.index[i]
            for j in range(annotate.shape[1]):
                col = annotate.columns[j]
                val = annotate.at[row, col]
                if txtcolors is None: color = 'black'
                else: 
                    fill = data.at[row, col]
                    color = txtcolors(row, col, val, fill)
                im.axes.text(j, i, valfmt(val), color=color, **kw)
    
    if save_as:
        fig.savefig(ospath.join(figures_path, save_as), 
                    dpi=300, transparent=True)
    else:
        return fig, ax, im

#%%
def plot_absolute(data=None, suffix='baseline_unopt'):
    if data is None:
        data = load_data(
            ospath.join(results_path, 'baseline_unopt_performance.xlsx'), sheet=0,
            header=[0,1], skiprows=[2,]
            )
        data.columns = [i[1].split(' [')[0] for i in data.columns]
    
    fill = pd.DataFrame()
    for var, col in data.items():
        fill[var] = (col-col.min())/(col.max()-col.min())
    fill.index = data.index
    fill = fill.T
    
    vals = pd.DataFrame()
    for var, col in data.items():
        if col.min() > 1e5:
            vals[var] = (col*1e-5).round(1)
        else:
            vals[var] = col.round(1)
    vals.index = data.index
    vals = vals.T

    def valfmtwnan(val):
        if str(val) == 'nan': return ''
        return f"{val:.1f}"

    def txtcolors(var, config, val, fill):
        if var in thresholds.index:
            if thresholds.at[var, config] < val:
                return 'red'
        if fill <= 0.6: return 'white'
        return 'black'

    heatmap(fill.loc[eff_vars,:], annotate=vals.loc[eff_vars,:], valfmt=valfmtwnan, 
            txtcolors=txtcolors, save_as=f'eff_{suffix}.png',
            row_labels=['COD', 'BOD', 'TN', 'NH$_4$-N', 'TP'])
    heatmap(fill.loc[op_vars,:], 
            show_ticklabels=False, colorbar=False,
            annotate=vals.loc[op_vars,:], valfmt=valfmtwnan, 
            txtcolors=txtcolors, save_as=f'op_{suffix}.png')


#%%
from math import nan
def cal_perc_diff(baseline, alternative):
    baseline.replace(0, nan, inplace=True)
    out = alternative.sub(baseline).div(baseline)
    return out

#%%
def plot_diff(data, suffix=''):
    fill = data.T.copy()
    fill[fill > 1] = 1
    fill[fill < -1] = -1
    
    vals = data.T    
    
    def valfmtpc(val):
        if str(val) == 'nan': return ''
        if abs(val) >= 0.095: return f"{val:.0%}"
        return f"{val:.1%}"

    def txtcolors(var, config, val, fill):
        if abs(fill) > 0.6: return 'white'
        return 'black'
    
    heatmap(fill.loc[eff_vars,:], cmap='bwr', vmin=-1, vmax=1,
            annotate=vals.loc[eff_vars,:], valfmt=valfmtpc, 
            txtcolors=txtcolors, annotate_kw=dict(size=10), 
            row_labels=['COD', 'BOD', 'TN', 'NH$_4$-N', 'TP'],
            save_as=f'deff_{suffix}.png')
    heatmap(fill.loc[op_vars,:], cmap='bwr', vmin=-1, vmax=1,
            show_ticklabels=False, colorbar=False,
            annotate=vals.loc[op_vars,:], valfmt=valfmtpc, 
            txtcolors=txtcolors, annotate_kw=dict(size=10), 
            save_as=f'dop_{suffix}.png')

#%%
data_handles = {
    'Baseline': ('baseline_unopt_performance', 0),   # file name, sheet
    'Adjusted': ('baseline_opt_performance', 'combined'),
    'UD10': ('UD_opt_performance', '10'),
    # 'UD30': ('UD_opt_performance', '30'),
    'UD100': ('UD_opt_performance', '100'),
    'HA': ('HA_opt_performance', 0),
    'ECS': ('ECS_opt_performance', 0)
    }

def diff_opex():
    dfs = {}
    for k,v in data_handles.items():
        file, sheet = v
        df = load_data(
            ospath.join(results_path, f'{file}.xlsx'), sheet=sheet,
            header=[0,1], skiprows=[2,]
            )
        df = df.loc[:,[('OPEX', 'Total OPEX [USD/d]')]]
        dfs[k] = df
    bl = dfs['Adjusted']
    data = pd.DataFrame()
    for k in ['UD10', 'UD100', 'HA', 'ECS']:
        alt = dfs[k]
        data[k] = cal_perc_diff(bl, alt).iloc[:,[0]]
    return data

def plot_opex_diff(data=None):
    if data is None: data = diff_opex()
    fill = pd.DataFrame()
    for var, col in data.items():
        fill[var] = (col-col.min())/(col.max()-col.min())
    fill.index = data.index
    fill = fill.T
    
    vals = data.T    
    def valfmtpc(val):
        if str(val) == 'nan': return ''
        if abs(val) >= 0.095: return f"{val:.0%}"
        return f"{val:.1%}"
    
    def txtcolors(var, config, val, fill):
        if val > 0: return 'red'
        if fill <= 0.25: return 'white'
        return 'black'
    
    fig, ax, im = heatmap(
        fill, colorbar=False, cmap='GnBu_r',
        annotate=vals, valfmt=valfmtpc, 
        txtcolors=txtcolors, annotate_kw=dict(size=10), 
        )
    
    cbar = ax.figure.colorbar(
        im, ax=ax, orientation='horizontal', pad=0.05, fraction=0.035,
        aspect=25, anchor=(0.6, 1.0)
        )
    cbar.ax.set_xticks(ticks=[0,1], labels=['max', 'min'], fontsize=12)
    cbar.ax.invert_xaxis()
    ax.text(x=2.4, y=4.25, s='OPEX Reduction', fontweight='bold', fontsize=13, )
    
    fig.savefig(ospath.join(figures_path, 'dopex.png'), 
                dpi=300, transparent=True)

#%%
if __name__ == '__main__':
    # df_baseline_opt = load_data(
    #     ospath.join(results_path, 'baseline_opt_performance.xlsx'), sheet='combined',
    #     header=[0,1], skiprows=[2,]
    #     )
    # df_baseline_opt.columns = [i[1].split(' [')[0] for i in df_baseline_opt.columns]
    # plot_absolute(df_baseline_opt, 'baseline_opt')
    # plot_absolute()
    
    # df_unopt = load_data(
    #     ospath.join(results_path, 'baseline_unopt_performance.xlsx'), sheet=0,
    #     header=[0,1], skiprows=[2,]
    #     )
    # df_unopt.columns = [i[1].split(' [')[0] for i in df_unopt.columns]

    # for tech in ('UD', 'HA', 'ECS'):
    #     alt = load_data(
    #         ospath.join(results_path, tech+'_performance.xlsx'), sheet=None,
    #         header=[0,1], skiprows=[2,]
    #         )
    #     for name, df in alt.items():
    #         if name == 'Sheet1': name = ''
    #         suffix = f'{name}{tech}'
    #         df.columns = [i[1].split(' [')[0] for i in df.columns]
    #         if 'NH4 recovery' in df.columns: df.drop(columns='NH4 recovery', inplace=True)
    #         diff = cal_perc_diff(df_unopt, df)
    #         plot_diff(diff, suffix)
        
    # df = load_data(
    #     ospath.join(results_path, 'UD_opt_performance.xlsx'), sheet='100',
    #     header=[0,1], skiprows=[2,]
    #     )
    # df.columns = [i[1].split(' [')[0] for i in df.columns]
    # diff = cal_perc_diff(df_baseline_opt, df)
    # plot_diff(diff, '100UD_opt')
    
    diff = diff_opex()
    plot_opex_diff(diff)
