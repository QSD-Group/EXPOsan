# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 16:44:59 2024

@author: Ga-Yeong Kim

to convert 'online' into 'dynamic influent'
"""

# Import packages
import os, pandas as pd, numpy as np, matplotlib.pyplot as plt
# from matplotlib import dates
from qsdsan.utils import load_data
from exposan.mecha import data_path

#%% Load dynamic influent data

# dynamic_inf_file = os.path.join(data_path, 'train_val_test_dynamic_influent_unit_deleted.xlsx')       # Time ['days'], irregular interval
# dynamic_inf = load_data(dynamic_inf_file, sheet='train_val_test_dynamic_influent')

dynamic_inf_file = os.path.join(data_path, 'evaluation_dynamic_influent_unit_deleted.xlsx')         # evaluation dataset
dynamic_inf = load_data(dynamic_inf_file, sheet='evaluation_dynamic_influent')

#%% Load online data

# online_file = os.path.join(data_path,'train_val_test_online_dataset_.xlsx')               # Time ['Datetime'], regular interval of 15 mins
# online = load_data(online_file, sheet='train_val_test_online_dataset')

online_file = os.path.join(data_path,'evaluation_online_dynamic_influent_integrated_precipitation_reconciled_induced_do.xlsx')                   # evaluation dataset
online = load_data(online_file, sheet='Sheet1')

#%%

from datetime import datetime
# day0 = datetime.strptime('5-5-2021 11:57 AM', '%m-%d-%Y %I:%M %p')                        # set the first time stamp as 'day0' - train
day0 = datetime.strptime('3-13-2022 10:57 PM', '%m-%d-%Y %I:%M %p')                        # set the first time stamp as 'day0' - eval

calc_day = lambda t_stamp: (t_stamp-day0).total_seconds()/60/60/24                        # calculate 'DATETIME' as day

t_stamp = online.index
t_online = calc_day(t_stamp).to_numpy()

t_intp = dynamic_inf.Time

#%%

def plot_raw_vs_intp(x1, y1, x2, y2):
    fig, ax = plt.subplots(figsize=(15, 3))                                          # figure size
    l1, = ax.plot(x1, y1, label='Raw')                                               # first plot = raw data
    l2, = ax.plot(x2, y2, label='Interpolated')                                      # second plot = smoothed result
    ax.legend(handles=[l1,l2])
    return fig, ax

#%%

do = pd.DataFrame()
do['t_stamp'] = t_intp
do['DO_1'] = np.interp(t_intp, t_online, online['DO_1'])

do['Q_in'] = np.interp(t_intp, t_online, online['Q_in'])  # evaluation only
do['X_BH'] = np.interp(t_intp, t_online, online['X_BH'])
do['X_I'] = np.interp(t_intp, t_online, online['X_I'])
do['X_ND'] = np.interp(t_intp, t_online, online['X_ND'])
do['X_P'] = np.interp(t_intp, t_online, online['X_P'])
do['X_S'] = np.interp(t_intp, t_online, online['X_S'])
do['S_ALK'] = np.interp(t_intp, t_online, online['S_ALK'])
do['S_I'] = np.interp(t_intp, t_online, online['S_I'])
do['S_ND'] = np.interp(t_intp, t_online, online['S_ND'])
do['S_NH'] = np.interp(t_intp, t_online, online['S_NH'])
do['S_NO'] = np.interp(t_intp, t_online, online['S_NO'])
do['S_S'] = np.interp(t_intp, t_online, online['S_S'])
do['X_BA'] = np.interp(t_intp, t_online, online['X_BA'])


plot_raw_vs_intp(t_online, online['DO_1'], do['t_stamp'], do['DO_1'])

# do.to_excel('results/train_val_test_DO_interpolated.xlsx')
do.to_excel('results/evaluation_DO_interpolated.xlsx')
