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
    create_systems,
    create_modelA,
    create_modelB,
    create_modelC,
    create_ss_model,
    run_model,
    run_modelB,
    run_ss_model,
    results_path, 
    figures_path
     )
import os
import numpy as np
import pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt, matplotlib.ticker as tk


mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams["figure.autolayout"] = True
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True
# mpl.rcParams['figure.facecolor'] = 'white'

#colors
#D81B60
#1E88E5
#FFC107

#%% Global variables
mdlA = create_modelA()
mdlB = create_modelB()
mdlC = create_modelC()
cmps = mdlB._system.flowsheet.stream.Effluent_B.components
AnR1 = mdlB._system.flowsheet.unit.AnR1

sysA, sysB, sysC = create_systems()
mssA = create_ss_model(sysA, R1_ID='H2E', R2_ID='CH4E')
mssC = create_ss_model(sysC, R1_ID='R1', R2_ID='R2')

N = 400
T = 200
t_step = 5

keys = ['X_ch', 'S_su', 'S_va', 'S_ac',
        'X_pr', 'S_aa', 'S_bu', 'S_h2',
        'X_li', 'S_fa', 'S_pro', 'S_ch4']
irows = [0,0,0,0,1,1,1,1,2,2,2,2]
icols = [0,1,2,3,0,1,2,3,0,1,2,3]
plots_wide = zip(keys, irows, icols)

# indices = [(4, 5), (0, 1)]
# names = [('H2', 'CH4'), ('rCOD_1', 'rCOD_2')]
indices = [(3,4), (5,6)]
names = [('eff_COD', 'rCOD'), ('H2', 'CH4')]
pairs = zip(indices, names)

S_keys = [ID for ID in cmps.IDs if ID.startswith('S_') and ID not in ('S_an', 'S_cat')]
X_keys = [ID for ID in cmps.IDs if ID.startswith('X_')]
gas_keys = ['S_h2_gas', 'S_ch4_gas', 'S_IC_gas']
irows = [0,0,0,1,1,1,2,2,2,3,3,3]
icols = [0,1,2,0,1,2,0,1,2,0,1,2]
s_plots = zip(S_keys, irows, icols)
x_plots = zip(X_keys, irows, icols)
g_plots = zip(gas_keys, [0,0,0], [0,1,2])


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
    if unit in ('CH4E', 'AnR2', 'R2'): j += len(AnR1._state_keys)
    if variable == 'S_IC_gas': i = 29
    elif variable == 'S_ch4_gas': i = 28
    elif variable == 'S_h2_gas': i = 27
    else: i = cmps.index(variable)
    return arr[:, i+j]

def analyze_timeseries(variable_getter, N, folder='', todf=True, **kwargs):
    outputs = {}
    for sample_id in range(N):
        arr = np.load(ospath.join(folder, f'state_{sample_id}.npy'))
        outputs[sample_id] = variable_getter(arr, **kwargs)
    if todf: outputs = pd.DataFrame.from_dict(outputs)
    return outputs

def analyze_vars(seed, N, prefix='sys', sys_ID='B', pickle=True):
    folder = os.path.join(results_path, f'{prefix}{sys_ID}_time_series_data_{seed}')
    state_vars = S_keys+X_keys+gas_keys
    dct = {}
    t_arr = np.arange(0, T+t_step, t_step)
    pcs = np.array([5, 25, 50, 75, 95])
    col_names = [f'{p}th percentile' for p in pcs]
    if sys_ID == 'A': units = ('H2E', 'CH4E')
    elif sys_ID == 'B': units = ('AnR1', 'AnR2')
    else: units = ('R1', 'R2')
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
        path = os.path.join(results_path, f'{prefix}{sys_ID}_state_vars_{seed}.pckl')
        save_pickle(dct, path)
    return dct

def plot_timeseries(seed, N, unit, sys_ID='B', data=None):
    if data is None:
        try:
            path = ospath.join(results_path, f'sys{sys_ID}_state_vars_{seed}.pckl')
            data = load_pickle(path)
        except:
            data = analyze_vars(seed, N, sys_ID=sys_ID)
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

def plot_scatter(seed, modelA, modelC, bl_metrics=None):
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
                             figsize=(nx*2.5, ny*2.5)) 
    xlct = tk.MaxNLocator(nbins=3, min_n_ticks=1)
    ylct = tk.MaxNLocator(nbins=3, min_n_ticks=1)
    for j in range(ny):
        ya = dfa_y.iloc[:,j].values
        yc = dfc_y.iloc[:,j].values
        try: yb = bl_metrics[j]
        except: yb = -1
        ymin = max(0, min(np.min(ya), np.min(yc), yb))
        ymax = max(np.max(ya), np.max(yc), yb)
        y_ticks = ylct.tick_values(ymin, ymax)
        for i in range(nx):
            xa = dfa_x.iloc[:,i].values
            xc = dfc_x.iloc[:,i].values
            x_ticks = xlct.tick_values(np.min(xa), np.max(xa))
            ax = axes[j,i]
            ax.axhline(y=yb, ls='-', lw=0.3, c='black')
            ax.scatter(xa, ya, marker='o', s=1, c='#D81B60')
            ax.scatter(xc, yc, marker='^', s=1, c='#1E88E5')
            ax.tick_params(axis='both', which='both', 
                           direction='inout', labelsize=10.5)
            ax.xaxis.set_ticks(x_ticks)
            ax.yaxis.set_ticks(y_ticks)
            if i > 0: ax.yaxis.set_ticklabels([])
            if j < ny-1: ax.xaxis.set_ticklabels([])
            ax2x = ax.secondary_xaxis('top')
            ax2x.xaxis.set_ticks(x_ticks)
            ax2x.tick_params(axis='x', which='both', direction='in')
            ax2x.xaxis.set_ticklabels([])
            ax2y = ax.secondary_yaxis('right')
            ax2y.yaxis.set_ticks(y_ticks)
            ax2y.tick_params(axis='y', which='both', direction='in')
            ax2y.yaxis.set_ticklabels([])
    fig.subplots_adjust(wspace=0., hspace=0.)
    fig.savefig(ospath.join(figures_path, f'AvC_table_{seed}.png'), dpi=300)
    return fig, axes

def plot_2d_scatter(seed, model, bl_metrics=None):
    npr = len(model.parameters)
    nmt = len(model.metrics)
    if model.table is None:
        path = ospath.join(results_path, f'sysC_table_{seed}.xlsx')
        model.table = load_data(path=path, header=[0,1])
    df_p = model.table.iloc[:, :npr]
    df_m = model.table.iloc[:, npr:]    
    nplot = ceil(nmt/2)
    nrow = ceil(nplot/2)
    xlct = tk.MaxNLocator(nbins=3, min_n_ticks=1)
    ylct = tk.MaxNLocator(nbins=3, min_n_ticks=1)
    for i in range(npr):
        fig, axes = plt.subplots(nrow, 2, sharex=False, sharey=False, 
                                 figsize=(6, nrow*2.5))
        c = df_p.iloc[:,i].values
        for j, ax in zip(range(nplot), axes.ravel()):
            x = df_m.iloc[:,j*2].values
            y = df_m.iloc[:,j*2+1].values
            try: xmin = np.min(x)
            except: breakpoint()
            xmax = np.max(x)
            ymin = np.min(y)
            ymax = np.max(y)
            im = ax.scatter(x, y, marker='o', s=1, c=c, cmap='viridis')
            if bl_metrics is not None:
                x_bl, y_bl = bl_metrics[j*2: (j*2+2)]
                ax.scatter(x_bl, y_bl, s=1.5, marker='^', c='black')
                if x_bl < xmin: xmin = x_bl
                elif x_bl > xmax: xmax = x_bl
                if y_bl < ymin: ymin = y_bl
                elif y_bl > ymax: ymax = y_bl
            ax.tick_params(axis='both', which='both', 
                           direction='inout', labelsize=10.5)
            x_ticks = xlct.tick_values(xmin, xmax)
            y_ticks = ylct.tick_values(ymin, ymax)
            ax.xaxis.set_ticks(x_ticks)
            ax.yaxis.set_ticks(y_ticks)
            ax2x = ax.secondary_xaxis('top')
            ax2x.xaxis.set_ticks(x_ticks)
            ax2x.tick_params(axis='x', which='both', direction='in')
            ax2x.xaxis.set_ticklabels([])
            ax2y = ax.secondary_yaxis('right')
            ax2y.yaxis.set_ticks(y_ticks)
            ax2y.tick_params(axis='y', which='both', direction='in')
            ax2y.yaxis.set_ticklabels([])
            cbar = fig.colorbar(im, ax=ax)
            cbar.ax.tick_params(labelsize=10.5)
        fig.subplots_adjust(wspace=0., hspace=0.)
        fig.savefig(ospath.join(figures_path, f'2dscatter_param{i}_{seed}.png'), dpi=300)
        del fig, axes

def plot_ss_convergence(seed, N, unit='CH4E', prefix='ss', sys_ID='A', data=None,
                        baseline_sys=None, baseline_unit_ID='AnR2'):
    if data is None:
        try:
            path = os.path.join(results_path, f'{prefix}{sys_ID}_state_vars_{seed}.pckl')
            data = load_pickle(path)
        except:
            data = analyze_vars(seed, N, sys_ID=sys_ID)
    if baseline_sys is not None:
         blu = getattr(baseline_sys.flowsheet.unit, baseline_unit_ID)
         bl_df = pd.DataFrame(blu.scope.record, columns=blu._state_keys)
         bl_ts = blu.scope.time_series
    fig, axes = plt.subplots(4, 3, sharex=True, figsize=(12,12))
    x_ticks = np.linspace(0, 200, 6)
    for group, par_size in [(s_plots, 's'), (x_plots, 'x'), (g_plots, 'g')]:
        if par_size == 'g':
            fig, axes = plt.subplots(1, 3, sharex=True, figsize=(12,3.2))
        else:
            fig, axes = plt.subplots(4, 3, sharex=True, figsize=(12,12))
        plts = deepcopy(group)        
        for var, ir, ic in plts:
            df = data[unit][var]
            if par_size == 'g': ax = axes[ic]
            else: ax = axes[ir, ic]
            ys = df.iloc[:,:N].values
            ax.plot(df.t, ys, color='grey', linestyle='-', alpha=0.25)
            if baseline_sys is not None:
                bl_y = bl_df.loc[:, var].values
                ax.plot(bl_ts, bl_y, color='#FFC107')
                ymin = np.min(bl_y)
                ymax = np.max(bl_y)
            else: 
                ymin = np.inf
                ymax = -1
            lct = tk.MaxNLocator(nbins=3, min_n_ticks=1)
            y_ticks = lct.tick_values(min(np.min(ys), ymin), max(np.max(ys), ymax))
            ax.xaxis.set_ticks(x_ticks)
            ax.yaxis.set_ticks(y_ticks)
            ax.tick_params(axis='both', which='both', 
                           direction='inout', labelsize=14)
            ax2x = ax.secondary_xaxis('top')
            ax2x.xaxis.set_ticks(x_ticks)
            ax2x.tick_params(axis='x', which='both', direction='in')
            ax2x.xaxis.set_ticklabels([])
            ax2y = ax.secondary_yaxis('right')
            ax2y.yaxis.set_ticks(y_ticks)
            ax2y.tick_params(axis='y', which='both', direction='in')
            ax2y.yaxis.set_ticklabels([])
        del plts
        fig.subplots_adjust(hspace=0., wspace=0.)
        fig.savefig(ospath.join(figures_path, f'{prefix}{sys_ID}_{par_size}t_{unit}.png'), dpi=300)
        # fig.savefig(ospath.join(figures_path, f'{prefix}{sys_ID}_{par_size}t_tau01.png'), dpi=300)
        del fig, axes

#%%
import seaborn as sns
from math import ceil

def f_grouping(cod_removal):
    if cod_removal > 50: return 'high' 
    else: return 'low'

def plot_dist_bygroup(seed, model, sys_ID='A'):
    group = model.table.loc[:,('System','Total COD removal [%]')].apply(f_grouping)
    n_param = len(model.parameters)
    ncol = 6
    nrow = ceil(n_param/ncol)
    x_df = model.table.iloc[:, :n_param]
    fig, axes = plt.subplots(nrow, ncol, sharey=True, figsize=(ncol*3, nrow*2.5))
    for col, ax in zip(x_df, axes.ravel()):
        unit = col[0].split('-')[-1]
        var = col[1].split(' 0')[0]
        sns.histplot(data=x_df, x=col, hue=group, bins=30, stat='probability',
                     cumulative=True,
                     common_norm=False, kde=True, legend=False, ax=ax)
        ax.set_xlabel(f'{unit} {var}_0')
    fig.subplots_adjust(hspace=0., wspace=0.)
    fig.savefig(ospath.join(figures_path, f'hist_ss{sys_ID}_{seed}'), dpi=300)


#%%
def run_UA_AvC(seed=None, N=N, T=T, t_step=t_step, plot=True):
    seed = seed or seed_RGT()
    run_model(mdlA, N, T, t_step, method='BDF', sys_ID='A', seed=seed)
    run_model(mdlC, N, T, t_step, method='BDF', sys_ID='C', seed=seed)
    outC = analyze_vars(seed, N, prefix='sys', sys_ID='C')
    for p in mdlC.parameters:
        if 'removal' in p.name: p.setter(0)
        else: p.setter(1e-4)
    temp_sysC = mdlC._system
    temp_sysC.simulate(t_span=(0, T), method='BDF')
    tmC = [m.getter() for m in mdlC.metrics]
    if plot:
        plot_scatter(seed, mdlA, mdlC, tmC)
        plot_ss_convergence(seed, N, unit='R1', prefix='sys', sys_ID='C', data=outC,
                            baseline_sys=temp_sysC, baseline_unit_ID='R1')
        plot_ss_convergence(seed, N, unit='R2', prefix='sys', sys_ID='C', data=outC,
                            baseline_sys=temp_sysC, baseline_unit_ID='R2')
    print(f'Seed used for uncertainty analysis of systems A & C is {seed}.')
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

def run_ss_AvC(seed=None, N=N, T=T, t_step=t_step, plot=True):
    seed = seed or seed_RGT()
    run_ss_model(mssA, N, T, t_step, sys_ID='A', R1_ID='H2E', R2_ID='CH4E', seed=seed)
    run_ss_model(mssC, N, T, t_step, sys_ID='C', R1_ID='R1', R2_ID='R2', seed=seed)
    outA = analyze_vars(seed, N, prefix='ss', sys_ID='A')
    outC = analyze_vars(seed, N, prefix='ss', sys_ID='C')
    if plot:
        plot_ss_convergence(seed, N, unit='CH4E', prefix='ss', sys_ID='A', data=outA)
        plot_ss_convergence(seed, N, unit='R2', prefix='ss', sys_ID='C', data=outC)
        plot_dist_bygroup(seed, mssA)
        plot_dist_bygroup(seed, mssC, sys_ID='C')
    print(f'Seed used for convergence test is {seed}.')
    return seed

#%%
# if __name__ == '__main__':
#     run_UA_AvC(seed=123, N=10)
    # run_ss_AvC(seed=223, N=100)
    # seed = 952
    # run_modelB(mdl, N, T, t_step, method='BDF', seed=seed)
    # out = analyze_vars(seed, N, 'B')
    # plot_timeseries(seed, N, 'AnR1', 'B')
    # plot_timeseries(seed, N, 'AnR2', 'B')
    # plot_kde2d_metrics(seed, 'B')
