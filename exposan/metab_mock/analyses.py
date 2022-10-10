# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from time import time
from copy import deepcopy
from qsdsan.utils import load_data, load_pickle, save_pickle, ospath
from qsdsan.stats import plot_uncertainties
from exposan.metab_mock import (
    create_modelA,
    create_modelB,
    create_modelC,
    run_model,
    run_modelB, 
    results_path, 
    figures_path
     )
import os
import numpy as np
import pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt, matplotlib.ticker as tk

mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams["figure.autolayout"] = True

#%% Global variables
mdlA = create_modelA()
mdlB = create_modelB()
mdlC = create_modelC()
cmps = mdlB._system.flowsheet.stream.Effluent_B.components
AnR1 = mdlB._system.flowsheet.unit.AnR1

N = 400
T = 120
t_step = 3

keys = ['X_ch', 'S_su', 'S_va', 'S_ac',
        'X_pr', 'S_aa', 'S_bu', 'S_h2',
        'X_li', 'S_fa', 'S_pro', 'S_ch4']
irows = [0,0,0,0,1,1,1,1,2,2,2,2]
icols = [0,1,2,3,0,1,2,3,0,1,2,3]
plots_wide = zip(keys, irows, icols)

# indices = [(4, 5), (0, 1)]
# names = [('H2', 'CH4'), ('rCOD_1', 'rCOD_2')]
indices = [(0,1), (2,3)]
names = [('eff_COD', 'rCOD'), ('H2', 'CH4')]
pairs = zip(indices, names)


def seed_RGT():
    files = os.listdir(results_path)
    seeds = [int(file_name[-3:]) for file_name in files if file_name.startswith('time_series_data')]
    seed = int(str(time())[-3:])
    if len(set(seeds)) >= 1000:
        raise RuntimeError('The program has run out of 3-digit seeds to use. Consider'
                           'clean up the results folder.')
    while seed in seeds:
        seed = (seed+1) % 1000
    return seed

#%% UA with time-series data
def single_state_var_getter(arr, unit, variable):
    j = 1
    if unit in ('AnR2', 'CH4E'): j += len(AnR1._state)
    i = cmps.index(variable)
    return arr[:, i+j]

def analyze_timeseries(variable_getter, N, folder='', todf=True, **kwargs):
    outputs = {}
    for sample_id in range(N):
        arr = np.load(ospath.join(folder, f'state_{sample_id}.npy'))
        outputs[sample_id] = variable_getter(arr, **kwargs)
    if todf: outputs = pd.DataFrame.from_dict(outputs)
    return outputs

def analyze_vars(seed, N, sys='B', pickle=True):
    folder = os.path.join(results_path, f'sys{sys}_time_series_data_{seed}')
    state_vars = keys
    dct = {}
    t_arr = np.arange(0, T+t_step, t_step)
    pcs = np.array([0, 5, 25, 50, 75, 95, 100])
    col_names = [f'{p}th percentile' for p in pcs]
    units = ('AnR1', 'AnR2') if sys == 'B' else ('H2E', 'CH4E')
    for u in units:
        dct[u] = {}
        for var in state_vars:
            df = analyze_timeseries(single_state_var_getter, N=N, folder=folder,
                                    variable=var, unit=u)
            quants = np.percentile(df, q = pcs, axis=1)
            df[col_names] = quants.T
            df['t'] = t_arr
            dct[u][var] = df
    if pickle:
        path = os.path.join(results_path, f'sys{sys}_state_vars_{seed}.pckl')
        save_pickle(dct, path)
    return dct

def plot_timeseries(seed, N, unit, sys='B', data=None):
    if data is None:
        try:
            path = ospath.join(results_path, f'sys{sys}_state_vars_{seed}.pckl')
            data = load_pickle(path)
        except:
            data = analyze_vars(seed, N, sys=sys)
    plts = deepcopy(plots_wide)    
    x_ticks = np.linspace(0, T, 5)
    fig, axes = plt.subplots(3, 4, sharex=True, figsize=(14,8))
    for var, ir, ic in plts:
        df = data[unit][var]
        ax = axes[ir, ic]
        # ax = axes[ic]
        ys = df.iloc[:,:N].values
        ax.plot(df.t, ys, color='grey', linestyle='-', alpha=0.05)
        l5, l95, = ax.plot(df.t, df.loc[:,['5th percentile', '95th percentile']].values,
                      color='blue', linestyle='--',
                      label='5th, 95th')
        l50, = ax.plot(df.t, df.loc[:,'50th percentile'],
                      color='black', linestyle='-.',
                      label='50th')
        lct = tk.MaxNLocator(nbins=3, min_n_ticks=1)
        y_ticks = lct.tick_values(np.min(ys), np.max(ys))
        ax.xaxis.set_ticks(x_ticks)
        ax.yaxis.set_ticks(y_ticks)
        ax.tick_params(axis='both', direction='inout', length=6, labelsize=14)
        ax2x = ax.secondary_xaxis('top')
        ax2x.xaxis.set_ticks(x_ticks)
        ax2x.tick_params(axis='x', direction='in')
        ax2x.xaxis.set_major_formatter(plt.NullFormatter())
        ax2y = ax.secondary_yaxis('right')
        ax2y.yaxis.set_ticks(y_ticks)
        ax2y.tick_params(axis='y', direction='in')
        ax2y.yaxis.set_major_formatter(plt.NullFormatter())
        # if ir == 0 and ic == 0: ax.legend(handles=[l5, l50, lbl])
        plt.ylim(y_ticks[0], y_ticks[-1])
        plt.subplots_adjust(wspace=0.2, hspace=0.1)
    fig.savefig(ospath.join(figures_path, f'{unit}_states.png'), dpi=300)
    return fig, axes

def plot_kde2d_metrics(seed, model, sys_ID):
    metrics = model.metrics
    if model.table is None:
        path = ospath.join(results_path, f'sys{sys_ID}_table_{seed}.xlsx')
        model.table = load_data(path=path, header=[0,1])
    for ids, names in pairs:
        i, j = ids
        x, y = names
        fig, ax = plot_uncertainties(model, x_axis=metrics[i], y_axis=metrics[j],
                                     kind='kde-box', center_kws={'fill':True},
                                     margin_kws={'width':0.5})
        ax0, ax1, ax2 = fig.axes # KDE, top box, right box
        ax0x = ax0.secondary_xaxis('top')
        ax0y = ax0.secondary_yaxis('right')
        for ax in (ax0, ax0x, ax0y):
            ax.tick_params(axis='both', which='both', direction='inout', width=1)

        for txt in ax0.get_xticklabels()+ax0.get_yticklabels():
            txt.set_fontsize(16)
        ax0x.set_xticklabels('')
        ax0y.set_yticklabels('')

        ax1.xaxis.set_visible(False)
        ax1.spines.clear()
        ax2.spines.clear()
        ax2.yaxis.set_visible(False)
        xl, xu = ax0.get_xlim()
        yl, yu = ax0.get_ylim()
        if yl < 0:
            ax0.set_ylim((0, yu))
            yl, yu = ax0.get_ylim()
        ax0.set_xlim((xl, xu))
        ax0.set_ylim((yl, yu))
        fig.savefig(ospath.join(figures_path, f'sys{sys_ID}_{x}_vs_{y}.png'), dpi=300)
        del fig, ax

def plot_scatter(seed, modelA, modelC):
    nx = len(modelA.parameters)
    ny = len(modelA.metrics)
    if modelA.table is None:
        path = ospath.join(results_path, f'sysA_table_{seed}.xlsx')
        modelA.table = load_data(path=path, header=[0,1])
    if modelC.table is None:
        path = ospath.join(results_path, f'sysC_table_{seed}.xlsx')
        modelC.table = load_data(path=path, header=[0,1])
    dfa_x = modelA.table.iloc[:, :nx]
    dfa_y = modelA.table.iloc[:, nx:]
    dfc_x = modelC.table.iloc[:, :nx]
    dfc_y = modelC.table.iloc[:, nx:]
    fig, axes = plt.subplots(ny, nx, sharex=False, sharey=False, 
                             figsize=(nx*2, ny*2))
    x_ticks = np.linspace(0.1, 0.9, 5)
    for j in range(ny):
        ya = dfa_y.iloc[:,j].values
        yc = dfc_y.iloc[:,j].values
        ylct = tk.MaxNLocator(nbins=3, min_n_ticks=1)
        y_ticks = ylct.tick_values(min(np.min(ya), np.min(yc)), 
                                   max(np.max(ya), np.max(yc)))
        for i in range(nx):
            xa = dfa_x.iloc[:,i].values
            xc = dfc_x.iloc[:,i].values
            xlct = tk.MaxNLocator(nbins=3, min_n_ticks=1)
            x_ticks = xlct.tick_values(np.min(xa), np.max(xa))
            ax = axes[j,i]            
            ax.scatter(xa, ya, marker='o', s=1, c='orange', alpha=0.7)
            ax.scatter(xc, yc, marker='^', s=1, c='blue', alpha=0.7)
            ax.tick_params(axis='both', direction='inout', length=4, labelsize=11)
            ax.xaxis.set_ticks(x_ticks)
            ax.yaxis.set_ticks(y_ticks)
            if i > 0: ax.yaxis.set_ticklabels([])
            if j < ny-1: ax.xaxis.set_ticklabels([])
            ax2x = ax.secondary_xaxis('top')
            ax2x.xaxis.set_ticks(x_ticks)
            ax2x.tick_params(axis='x', direction='in', length=2)
            ax2x.xaxis.set_ticklabels([])
            ax2y = ax.secondary_yaxis('right')
            ax2y.yaxis.set_ticks(y_ticks)
            ax2y.tick_params(axis='y', direction='in', length=2)
            ax2y.yaxis.set_ticklabels([])
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.savefig(ospath.join(figures_path, 'AvC_table.png'), dpi=300)
    return fig, axes

#%%
def run_UA_AvC(seed=None, N=N, T=T, t_step=t_step, plot=True):
    seed = seed or seed_RGT()
    # run_model(mdlA, N, T, t_step, method='BDF', sys_ID='A', seed=seed)
    run_model(mdlC, N, T, t_step, method='BDF', sys_ID='C', seed=seed)
    if plot:
        plot_scatter(seed, mdlA, mdlC)
    print(f'Seed used for uncertainty analysis of system A is {seed}.')
    for mdl in (mdlA, mdlC):
        for p in mdl.parameters:
            p.setter(p.baseline)
    return seed

def run_UA_sysB(seed=None, N=N, T=T, t_step=t_step, plot=True):
    seed = seed or seed_RGT()
    run_modelB(mdlB, N, T, t_step, method='BDF', seed=seed)
    out = analyze_vars(seed, N, sys='B')
    if plot: 
        for u in ('AnR1', 'AnR2'):
            plot_timeseries(seed, N, u, data=out, sys='B')
        plot_kde2d_metrics(seed, mdlB, sys_ID='B')
    print(f'Seed used for uncertainty analysis of system B is {seed}.')
    for p in mdlB.parameters:
        p.setter(p.baseline)
    return seed

#%%
if __name__ == '__main__':
    run_UA_AvC(seed=123, N=10)
    # seed = 952
    # run_modelB(mdl, N, T, t_step, method='BDF', seed=seed)
    # out = analyze_vars(seed, N, 'B')
    # plot_timeseries(seed, N, 'AnR1', 'B')
    # plot_timeseries(seed, N, 'AnR2', 'B')
    # plot_kde2d_metrics(seed, 'B')
