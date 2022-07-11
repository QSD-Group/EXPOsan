# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# from time import time
from copy import deepcopy
from qsdsan.utils import ospath, load_pickle, save_pickle
# from qsdsan.stats import plot_uncertainties, get_correlations
from exposan.adm import cmps, results_path, figures_path
import numpy as np, pandas as pd, os, matplotlib as mpl, \
    matplotlib.pyplot as plt, matplotlib.ticker as tk

mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams["figure.autolayout"] = True

N = 100
T = 200
t_step = 5

#%%
S_keys = [ID for ID in cmps.IDs if ID.startswith('S_') and ID not in ('S_an', 'S_cat')]
X_keys = [ID for ID in cmps.IDs if ID.startswith('X_')]
gas_keys = ['S_h2_gas', 'S_ch4_gas', 'S_IC_gas']
irows = [0,0,0,1,1,1,2,2,2,3,3,3]
icols = [0,1,2,0,1,2,0,1,2,0,1,2]
s_plots = zip(S_keys, irows, icols)
x_plots = zip(X_keys, irows, icols)
g_plots = zip(gas_keys, [0,0,0], [0,1,2])

bm_ss = [0.01195483	, 0.00531474, 0.098621401, 0.011625006, 0.01325073, 0.015783666,
         0.197629717, 2.35945e-7, 0.055088776, 1.832134448, 1.823217421, 0.328697664, 
         0.308697664, 0.02794724, 0.102574106, 0.02948305, 0.420165982	, 1.179171799, 
         0.243035345, 0.431921106, 0.137305909, 0.760562658, 0.317022953, 25.61739533, 
         6.40065e-7	, 0.025400113, 0.014150535]
bm_ss = dict(zip(S_keys+X_keys+gas_keys, bm_ss))

def single_state_var_getter(arr, variable):
    if variable == 'S_IC_gas': i = 29
    elif variable == 'S_ch4_gas': i = 28
    elif variable == 'S_h2_gas': i = 27
    else:
        i = cmps.index(variable)
    return arr[:, i+1]
    
def analyze_timeseries(variable_getter, N, folder='', todf=True, **kwargs):
    outputs = {}
    for sample_id in range(N):
        arr = np.load(os.path.join(folder, f'state_{sample_id}.npy'))
        outputs[sample_id] = variable_getter(arr, **kwargs)
    if todf: outputs = pd.DataFrame.from_dict(outputs)
    return outputs

def analyze_vars(seed, N, pickle=True):
    folder = os.path.join(results_path, f'time_series_data_{seed}')
    state_vars = [ID for ID in cmps.IDs if ID not in ('H2O', 'S_cat', 'S_an')]
    state_vars += gas_keys
    dct = {}
    t_arr = np.arange(0, T+t_step, t_step)
    pcs = np.array([0, 5, 25, 50, 75, 95, 100])
    col_names = [f'{p}th percentile' for p in pcs]
    for var in state_vars:
        df = analyze_timeseries(single_state_var_getter, N=N, folder=folder,
                                variable=var)
        quants = np.percentile(df, q = pcs, axis=1)
        df[col_names] = quants.T
        df['t'] = t_arr
        dct[var] = df
    if pickle:
        path = os.path.join(results_path, f'state_vars_{seed}.pckl')
        save_pickle(dct, path)
    return dct

def plot_yt_w_diff_init(seed, N, data=None):
    if data is None:
        path = os.path.join(results_path, f'state_vars_{seed}.pckl')
        data = load_pickle(path)
    
    x_ticks = np.linspace(0, 200, 6)
    
    for group, prefix in [(s_plots, 's'), (x_plots, 'x'), (g_plots, 'g')]:
        if prefix == 'g':
            fig, axes = plt.subplots(1, 3, sharex=True, figsize=(12,3.2))
        else:
            fig, axes = plt.subplots(4, 3, sharex=True, figsize=(12,12))
        plts = deepcopy(group)        
        for var, ir, ic in plts:
            df = data[var]
            if prefix == 'g': ax = axes[ic]
            else: ax = axes[ir, ic]
            ys = df.iloc[:,:N].values
            ax.plot(df.t, ys, color='grey', linestyle='-', alpha=0.25)
            ax.hlines(y=bm_ss[var], xmin=0, xmax=T, linestyle='-.', color='r')
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
        del plts
        plt.subplots_adjust(hspace=0.1)
        fig.savefig(ospath.join(figures_path, f'{prefix}_yt_wdiff_init.png'), dpi=300)
        del fig, axes

#%%
if __name__ == '__main__':
    seed = 119
    ##### Uncertainty analysis results with different initial conditions #####
    out = analyze_vars(seed, N)
    plot_yt_w_diff_init(seed, N, data=out)
    # plot_yt_w_diff_init(seed, N)    