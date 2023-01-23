# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import os
from time import time
# from copy import deepcopy
from qsdsan.utils import ospath, load_data
from qsdsan.stats import get_correlations
from exposan.pm2_batch import (
    results_path, 
    # figures_path, 
    create_model,
    run_uncertainty
    )
# import numpy as np, pandas as pd, os, matplotlib as mpl, \
#     matplotlib.pyplot as plt, matplotlib.ticker as tk

# mpl.rcParams['font.sans-serif'] = 'arial'
# mpl.rcParams["figure.autolayout"] = True

# N = 10000
# T = 20
# t_step = 1/24

N = 10000
T = 7
t_step = 0.25/24

#%%
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

#%%
def run_UA_SA(seed=None, N=N, T=T, t_step=t_step, rmse_thresholds=[]):
    seed = seed or seed_RGT()
    mdl = create_model()
    #!!! Need to define algorithm that determines an acceptable rmse for each metric
    run_uncertainty(mdl, N, T, t_step, seed=seed)
    D, p = get_correlations(mdl, kind='KS', thresholds=rmse_thresholds,
                            file=ospath.join(results_path, f'KS_test_{seed}.xlsx'))
    #!!! Can add functions to plot SA results
    return mdl

#%%
if __name__ == '__main__':
    seed = 119
    mdl = run_UA_SA(seed=seed)
    
    
#%%

# import experiemental data
folder = ospath.dirname(__file__)
file = ospath.join(folder, 'data/batch_exp_result.xlsx')
exp_data = load_data(file, sheet=None)


# import uncertainty & sensitivity analysis results as time-series data

simul_data = 
# one-to-one correspondence

corr_VSS = exp_data.merge(simul_data[['t_stamp', 'VSS']], on='t_stamp', how='left')
corr_NH = exp_data.merge(simul_data[['t_stamp', 'S_NH']], on='t_stamp', how='left')
corr_P = exp_data.merge(simul_data[['t_stamp', 'S_P']], on='t_stamp', how='left')

# calculation of RMSE
