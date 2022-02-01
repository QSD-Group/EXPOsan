#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

'''
#TODO
Other things that I think we might want to include:
    - A figure in the SI (or just on GitHub) showing that we can converge to
    similar steady-states conditions
'''

# from qsdsan.utils import ords
from time import time
from qsdsan.utils import load_data, load_pickle, save_pickle
from qsdsan.stats import plot_uncertainties, get_correlations
from exposan.bsm1 import model_bsm1 as mdl, model_2dv as mdl2, model_ss as mdls, \
    cmps, run_uncertainty, analyze_timeseries, run_wdiff_init, \
    results_path, figures_path
import os
import numpy as np
import matplotlib.pyplot as plt


N = 1000
T = 50
t_step = 1

#%% UA with time-series data

idx = mdl._system._load_state()[1]
def single_state_var_getter(arr, unit, variable):
    head, tail = idx[unit]
    if unit == 'C1':
        if variable.startswith('TSS'): 
            i_layer = int(variable.lstrip('TSS'))
            i = 1 + head + len(cmps) + i_layer
        elif variable == 'Q': i = 1 + head + len(cmps)
        else: 
            i_cmp = cmps.index(variable)
            i = 1 + head + i_cmp
    else:
        if variable == 'Q': i = tail
        else: 
            i_cmp = cmps.index(variable)
            i = 1 + head + i_cmp
    if variable == 'S_ALK': return arr[:, i]/12
    else: return arr[:, i]

def analyze_SE_vars(seed, N, pickle=True):
    folder = os.path.join(results_path, f'time_series_data_{seed}')
    state_vars = [ID for ID in cmps.IDs if ID not in ('S_N2', 'H2O', 'S_I')]
    dct = {}
    t_arr = np.arange(51)
    pcs = np.array([0, 5, 25, 50, 75, 95, 100])
    col_names = [f'{p}th percentile' for p in pcs]
    for var in state_vars:
        df = analyze_timeseries(single_state_var_getter, N=N, folder=folder, 
                                unit='S1', variable=var)
        quants = np.percentile(df, q = pcs, axis=1)
        df[col_names] = quants.T
        df['t'] = t_arr
        dct[var] = df
    if pickle:
        path = os.path.join(results_path, f'SE_vars_{seed}.pckl')
        save_pickle(dct, path)
    return dct

keys = ['S_S', 'S_O', 'S_ALK', 
        'S_NO', 'S_NH', 'S_ND', 
        'X_S', 'X_I', 'X_ND',
        'X_BH', 'X_BA', 'X_P']
irows = [0,0,0,1,1,1,2,2,2,3,3,3]
icols = [0,1,2,0,1,2,0,1,2,0,1,2]
plots = zip(keys, irows, icols)
def plot_SE_timeseries(seed, N, data=None):
    if data is None:
        try:
            path = os.path.join(results_path, f'SE_vars_{seed}.pckl')
            data = load_pickle(path)
        except:
            data = analyze_SE_vars(seed)
    # bm_data = load_data(os.path.join(results_path, 'matlab_exported_data.xlsx'),
    #                     sheet='Data_ASin_changed', header=1, skiprows=0, usecols='A,CU:DF')
    # bm_data.columns = [col.rstrip('.6') for col in bm_data.columns]
    bl_data = load_data(os.path.join(results_path, 'sol_50d_BDF.xlsx'),
                        skiprows=[0,2], header=0, usecols='B,DF:DQ')
    bl_data.columns = [col.split(' ')[0] for col in bl_data.columns]
    bl_data.S_ALK = bl_data.S_ALK/12
    fig, axes = plt.subplots(4, 3, sharex=True, figsize=(12,12))
    for var, ir, ic in plots:
        df = data[var]
        ax = axes[ir, ic]
        ax.plot(df.t, df.iloc[:,:N], 
                color='grey', linestyle='-', alpha=0.05)
        lbl, = ax.plot(df.t, bl_data.loc[:, var],
                     color='orange', linestyle='-', linewidth=2.5,
                      label='Base')
        l5, l95, = ax.plot(df.t, df.loc[:,['5th percentile', '95th percentile']],
                      color='blue', linestyle='--', 
                      label='5th, 95th')
        l50, = ax.plot(df.t, df.loc[:,'50th percentile'],
                      color='black', linestyle='-.',
                      label='50th')
        ax.tick_params(axis='both', direction='inout')
        ax2 = ax.secondary_yaxis('right')
        ax2.tick_params(axis='y', direction='in')
        ax2.yaxis.set_major_formatter(plt.NullFormatter())
        # lbm,  = ax.plot(df.t, bm_data.loc[:, var],
        #               color='black', linestyle='-.', label='Matlab/Simulink')
        # ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        if ir == 0 and ic == 0: ax.legend(handles=[l5, l50, lbl])
        # if ir == 0 and ic == 0: ax.legend(handles=[lbl, lbm])
    plt.subplots_adjust(hspace=0.1)
    fig.savefig(os.path.join(figures_path, 'effluent_yt.png'), dpi=300)
    return fig, ax

def UA_w_all_params():
    seed = int(str(time())[-3:])
    run_uncertainty(mdl, N, T, t_step, method='BDF', seed=seed)
    out = analyze_SE_vars(seed)
    fig, ax = plot_SE_timeseries(seed, N, data=out)
    
#%% 
#!!! DO NOT UNCOMMENT !!!
# f_b2cod = cmps.f_BOD5_COD
# icod = cmps.i_COD
# def SE_BOD5_getter(arr):
#     head, tail = idx['S1']
#     state = arr[-1, head:(tail-1)]
#     return sum(state * f_b2cod * icod)

# seed = 157
# folder = os.path.join(results_path, f'time_series_data_{seed}')
# dct_bod = analyze_timeseries(SE_BOD5_getter, N, folder=folder, todf=False)
# path = os.path.join(results_path, f'table_{seed}.xlsx')
# df = load_data(path, header=[0,1])
# df['Effluent', 'Effluent BOD5 [mg/L]'] = list(dct_bod.values())
# df.to_excel(path, header=[0,1])

#%% SA to narrow down the number of DVs for further analyses

indices = [(0, 2), (1, 3), (4, 5)]
names = [('COD', 'TN'), ('BOD5', 'TKN'), ('TSS', 'DSP')]
pairs = zip(indices, names)
def plot_kde2d_metrics(model):
    metrics = model.metrics
    for ids, names in pairs:
        i, j = ids
        x, y = names
        fig, ax = plot_uncertainties(model, x_axis=metrics[i], y_axis=metrics[j],
                                      kind='kde-hist', center_kws={'fill':True},
                                      margin_kws={'kde':True, 'fill':True})
        fig.savefig(os.path.join(figures_path, f'{x}_vs_{y}.png'))
    fig, ax = plot_uncertainties(model, x_axis=metrics[6], kind='hist', 
                                 center_kws={'kde':True})
    fig.savefig(os.path.join(figures_path, 'SRT_hist.png'))

thresholds = [100, 30, 19, 3, 30] # Effluent COD, BOD5, TN, TKN, TSS
def run_ks_test(model):
    metrics_df = model.table['Effluent']
    to_analyze =  [i for i in range(5) if sum(metrics_df.iloc[:,i] > thresholds[i]) > 10]
    key_metrics = [model.metrics[i] for i in to_analyze]
    thrsh = [thresholds[i] for i in to_analyze]
    path = os.path.join(results_path, 'KS_test.xlsx')
    Ds, ps = get_correlations(model, input_y=key_metrics, kind='KS', 
                              file=path, thresholds=thrsh)

def run_sensitivity(seed):
    mdl.table = load_data(os.path.join(results_path, f'table_{seed}.xlsx'), header=[0,1])
    plot_kde2d_metrics(mdl)
    run_ks_test(mdl)


#%% 2-DV heatmap plotting
def meshgrid_sample(p1, p2, n):
    xl, xu = p1.distribution.lower[0], p1.distribution.upper[0]
    yl, yu = p2.distribution.lower[0], p2.distribution.upper[0]
    x = np.linspace(xl, xu, n)
    y = np.linspace(yl, yu, n)
    xx, yy = np.meshgrid(x, y)
    samples = np.array([xx.reshape(n**2), yy.reshape(n**2)])
    return samples.T, xx, yy

irs = [0, 0, 1, 1, 2, 2] 
ics = [0, 1, 0, 1, 0, 1]
def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s}" if plt.rcParams["text.usetex"] else f"{s}"

from exposan.bsm1 import system as sys 
xbl = sys.aer1.Q_air
ybl = sys.Q_was
def plot_heatmaps(xx, yy, n, model=None, path=''):
    if model: data = model.table
    else: 
        path = path or os.path.join(results_path, 'table_2dv.xlsx')
        data = load_data(path, header=[0,1])
    zs = data.loc[:,['Effluent', 'WAS']].to_numpy(copy=True)
    zs = zs.T
    zz = zs.reshape((zs.shape[0], n, n))
    fig, axes = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(9,12))
    for z, ir, ic in zip(zz, irs, ics):
        ax = axes[ir, ic]
        pos = ax.imshow(z, aspect='auto', extent=(xx[0,0], xx[0,-1], yy[0,0], yy[-1,0]),
                        interpolation='spline16', origin='lower')
        fig.colorbar(pos, ax=ax)
        cs = ax.contour(xx, yy, z, colors='white', origin='lower', linestyles='dashed',
                        linewidths=1, extent=(xx[0,0], xx[0,-1], yy[0,0], yy[-1,0]))
        ax.clabel(cs, cs.levels, inline=True, fmt=fmt)
        ax.plot(xbl, ybl, marker='D', mec='black', mew=1.5, mfc='white')
    fig.savefig(os.path.join(figures_path, 'heatmaps.png'), dpi=300)
    return fig, ax
    
def run_mapping(model, n, T, t_step, method='BDF', mpath='', tpath=''):
    x = model.parameters[0]
    y = model.parameters[1]
    samples, xx, yy = meshgrid_sample(x, y, n)
    model.load_samples(samples)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    model.evaluate(
        state_reset_hook='reset_cache',
        t_span=t_span,
        t_eval=t_eval,
        method=method,
        export_state_to=tpath
        )
    if mpath: model.table.to_excel(mpath)
    return xx, yy

def dv_analysis(n=20, T=50, t_step=1, save_to='table_2dv.xlsx'):
    xx, yy = run_mapping(mdl2, n, T, t_step, 
                         mpath=os.path.join(results_path, save_to))
    plot_heatmaps(xx, yy, n, mdl2)
    # plot_heatmaps(xx, yy, n, path=os.path.join(results_path, save_to))
    
#%% 50-d simulations with different initial conditions
keys_wide = ['S_S', 'S_NO', 'X_S', 'X_BH', 
             'S_O', 'S_NH', 'X_I', 'X_BA',
             'S_ALK', 'S_ND', 'X_ND', 'X_P']
irows_wide = [0,0,0,0,1,1,1,1,2,2,2,2]
icols_wide = [0,1,2,3,0,1,2,3,0,1,2,3]
plots_wide = zip(keys_wide, irows_wide, icols_wide)

def plot_SE_yt_w_diff_init(seed, N, data=None):
    if data is None:
        try:
            path = os.path.join(results_path, f'SE_vars_{seed}.pckl')
            data = load_pickle(path)
        except:
            data = analyze_SE_vars(seed)
    bm_data = load_data(os.path.join(results_path, 'matlab_exported_data.xlsx'),
                        sheet='Data_ASin_changed', header=1, skiprows=0, usecols='A,CU:DF')
    bm_data.columns = [col.rstrip('.6') for col in bm_data.columns]
    fig, axes = plt.subplots(3, 4, sharex=True, figsize=(15, 8))
    for var, ir, ic in plots_wide:
        df = data[var]
        ax = axes[ir, ic]
        ax.plot(df.t, df.iloc[:,:N], 
                color='grey', linestyle='-', alpha=0.25)
        lbm, = ax.plot(df.t, bm_data.loc[:, var],
                       color='red', linestyle='-.', label='MATLAB/Simulink')
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        ax.tick_params(axis='both', direction='inout')
        ax2 = ax.secondary_yaxis('right')
        ax2.tick_params(axis='y', direction='in')
        ax2.yaxis.set_major_formatter(plt.NullFormatter())
        if var == 'X_ND': ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.3f'))
        if ir == 0 and ic == 0: ax.legend(handles=[lbm])
    plt.subplots_adjust(wspace=0.21, hspace=0.1)
    fig.savefig(os.path.join(figures_path, 'effluent_yt_wdiff_init_wide.png'), dpi=300)
    return fig, ax

def UA_w_diff_inits(N=100):
    seed = int(str(time())[-3:])    
    # run_uncertainty(mdls, N, T, t_step, method='BDF', seed=seed)
    run_wdiff_init(mdls, N, T, t_step, method='BDF', seed=seed)
    out = analyze_SE_vars(seed, N)
    fig, ax = plot_SE_yt_w_diff_init(seed, N, data=out)

#%%
if __name__ == '__main__':
    # UA_w_all_params()
    # run_sensitivity(157)
    # dv_analysis()
    # plot_SE_timeseries(seed=157, N=1000)
    plot_SE_yt_w_diff_init(seed=708, N=100)
    # UA_w_diff_inits()