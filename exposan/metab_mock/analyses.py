# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 17:05:15 2022

@author: joy_c
"""
from time import time
from copy import deepcopy
from qsdsan.utils import load_data, load_pickle, save_pickle, ospath
from qsdsan.stats import plot_uncertainties
from exposan.metab_mock import model_ua as mdl, system as s,\
    run_uncertainty, results_path, figures_path
import os
import numpy as np
import pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt, matplotlib.ticker as tk

mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams["figure.autolayout"] = True

#%% Global variables

N = 400
T = 120
t_step = 3

keys = ['X_ch', 'S_su', 'S_va', 'S_ac',
        'X_pr', 'S_aa', 'S_bu', 'S_h2',
        'X_li', 'S_fa', 'S_pro', 'S_ch4']
irows = [0,0,0,0,1,1,1,1,2,2,2,2]
icols = [0,1,2,3,0,1,2,3,0,1,2,3]
plots_wide = zip(keys, irows, icols)

indices = [(4, 5), (0, 1)]
names = [('H2', 'CH4'), ('rCOD_1', 'rCOD_2')]
pairs = zip(indices, names)

cmps = s.cmps
H2E = s.H2E

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
    if unit == 'CH4E': j += len(H2E._state)
    i = cmps.index(variable)
    return arr[:, i+j]

def analyze_timeseries(variable_getter, N, folder='', todf=True, **kwargs):
    outputs = {}
    for sample_id in range(N):
        arr = np.load(ospath.join(folder, f'state_{sample_id}.npy'))
        outputs[sample_id] = variable_getter(arr, **kwargs)
    if todf: outputs = pd.DataFrame.from_dict(outputs)
    return outputs

def analyze_vars(seed, N, pickle=True):
    folder = os.path.join(results_path, f'time_series_data_{seed}')
    state_vars = keys
    dct = {}
    t_arr = np.arange(0, T+t_step, t_step)
    pcs = np.array([0, 5, 25, 50, 75, 95, 100])
    col_names = [f'{p}th percentile' for p in pcs]
    for u in ('H2E', 'CH4E'):
        dct[u] = {}
        for var in state_vars:
            df = analyze_timeseries(single_state_var_getter, N=N, folder=folder,
                                    variable=var, unit=u)
            quants = np.percentile(df, q = pcs, axis=1)
            df[col_names] = quants.T
            df['t'] = t_arr
            dct[u][var] = df
    if pickle:
        path = os.path.join(results_path, f'state_vars_{seed}.pckl')
        save_pickle(dct, path)
    return dct

def plot_timeseries(seed, N, unit, data=None):
    if data is None:
        try:
            path = ospath.join(results_path, f'state_vars_{seed}.pckl')
            data = load_pickle(path)
        except:
            data = analyze_vars(seed, N)
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

def plot_kde2d_metrics(seed):
    metrics = mdl.metrics
    if mdl.table is None:
        path = ospath.join(results_path, f'table_{seed}.xlsx')
        mdl.table = load_data(path=path, header=[0,1])
    for ids, names in pairs:
        i, j = ids
        x, y = names
        fig, ax = plot_uncertainties(mdl, x_axis=metrics[i], y_axis=metrics[j],
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
        fig.savefig(ospath.join(figures_path, f'{x}_vs_{y}.png'), dpi=300)
        del fig, ax

def UA_w_all_params(seed=None, N=N, T=T, t_step=t_step, plot=True):
    seed = seed or seed_RGT()
    run_uncertainty(mdl, N, T, t_step, method='BDF', seed=seed)
    out = analyze_vars(seed, N)
    if plot: 
        for u in ('H2E', 'CH4E'):
            plot_timeseries(seed, N, u, data=out)
        plot_kde2d_metrics(mdl)
    print(f'Seed used for uncertainty analysis with all parameters is {seed}.')
    for p in mdl.parameters:
        p.setter(p.baseline)
    return seed

#%%
if __name__ == '__main__':
    seed = 952
    # run_uncertainty(mdl, N, T, t_step, method='BDF', seed=seed)
    # out = analyze_vars(seed, N)
    plot_timeseries(seed, N, 'H2E')
    plot_timeseries(seed, N, 'CH4E')
    # plot_kde2d_metrics(seed)
