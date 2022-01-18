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
from exposan.bsm1 import model_bsm1 as mdl, cmps,\
    run_uncertainty, analyze_timeseries, results_path, figures_path
import os
import numpy as np
import matplotlib.pyplot as plt


N = 1000
T = 50
t_step = 1

#%%
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

def analyze_SE_vars(seed, pickle=True):
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
    
def plot_SE_timeseries(seed, data=None):
    if data is None:
        try:
            path = os.path.join(results_path, f'SE_vars_{seed}.pckl')
            data = load_pickle(path)
        except:
            data = analyze_SE_vars(seed)
    bm_data = load_data(os.path.join(results_path, 'matlab_exported_data.xlsx'),
                        header=1, skiprows=0, usecols='A,CU:DF')
    bm_data.columns = [col.rstrip('.6') for col in bm_data.columns]
    fig, axes = plt.subplots(4, 3, sharex=True, figsize=(12, 12))
    keys = ['S_S', 'S_O', 'S_ALK', 
            'S_NO', 'S_NH', 'S_ND', 
            'X_S', 'X_I', 'X_ND',
            'X_BH', 'X_BA', 'X_P']
    # units = ['mg COD/L', 'mg O/L', 'mmol/L'] + \
    #         ['mg N/L'] * 3 + \
    #         ['mg COD/L'] * 2 + ['mg N/L'] + \
    #         ['mg COD/L'] * 3
    irows = [0,0,0,1,1,1,2,2,2,3,3,3]
    icols = [0,1,2,0,1,2,0,1,2,0,1,2]
    for var, ir, ic in zip(keys, irows, icols):
        df = data[var]
        ax = axes[ir, ic]
        ax.plot(df.t, df.iloc[:,:N], 
                color='grey', linestyle='-', alpha=0.05)
        ax.plot(df.t, df.loc[:,['5th percentile', '95th percentile']],
                color='blue', linestyle='--')
        ax.plot(df.t, df.loc[:,'50th percentile'],
                color='orange', linestyle='-', linewidth=2)
        ax.plot(df.t, bm_data.loc[:, var],
                color='black', linestyle='-.')
        # ax.set(title=f'{var} [{unit}]')
    # axes[-1, 1].set_xlabel('Time [d]')
    fig.savefig(os.path.join(figures_path, 'effluent_yt.png'), dpi=300)
    return fig, ax

    
#%%
if __name__ == '__main__':
    # seed = int(str(time())[-3:])
    # run_uncertainty(mdl, N, T, t_step, method='BDF', seed=seed)
    # out = analyze_SE_vars(157)
    # fig, ax = plot_SE_timeseries(157, data=out)
    fig, ax = plot_SE_timeseries(157)
