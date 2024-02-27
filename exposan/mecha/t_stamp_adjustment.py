# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 16:44:59 2024

@author: Ga-Yeong Kim

to convert 'dynamic influent' into 'online'
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

online_file = os.path.join(data_path,'evaluation_online_dataset_.xlsx')                   # evaluation dataset
online = load_data(online_file, sheet='evaluation_online_dataset')

#%%

from datetime import datetime
# day0 = datetime.strptime('5-5-2021 11:57 AM', '%m-%d-%Y %I:%M %p')                        # set the first time stamp as 'day0' - train
day0 = datetime.strptime('3-13-2022 10:57 PM', '%m-%d-%Y %I:%M %p')                        # set the first time stamp as 'day0' - eval

calc_day = lambda t_stamp: (t_stamp-day0).total_seconds()/60/60/24                        # calculate 'DATETIME' as day

t_stamp = online.index
t_intp = calc_day(t_stamp).to_numpy()

#%%

def plot_raw_vs_intp(x1, y1, x2, y2):
    fig, ax = plt.subplots(figsize=(15, 3))                                          # figure size
    l1, = ax.plot(x1, y1, label='Raw')                                               # first plot = raw data
    l2, = ax.plot(x2, y2, label='Interpolated')                                      # second plot = smoothed result
    ax.legend(handles=[l1,l2])
    return fig, ax

#%%
columns=['Q_in','X_BH','X_I','X_ND','X_P','X_S',
         'S_ALK','S_I','S_ND','S_NH','S_NO','S_O','S_S','X_BA']

for i in range(14):
    online[columns[i]] = np.interp(t_intp, dynamic_inf.Time, dynamic_inf[columns[i]])
    plot_raw_vs_intp(dynamic_inf.Time, dynamic_inf[columns[i]], t_intp, online[columns[i]])

# online.to_excel('results/train_val_test_online_dynamic_influent_integrated.xlsx')
online.to_excel('results/evaluation_online_dynamic_influent_integrated.xlsx')

#%%

# columns=['Q_in','X_BH','X_I','X_ND','X_P','X_S',
#          'S_ALK','S_I','S_ND','S_NH','S_NO','S_O','S_S','X_BA']

# dynamic_inf_intp=pd.DataFrame()
# dynamic_inf_intp['t_stamp']=t_intp

# for i in range(14):
#     dynamic_inf_intp[columns[i]] = np.interp(t_intp, dynamic_inf.Time, dynamic_inf[columns[i]])
#     plot_raw_vs_intp(dynamic_inf.Time, dynamic_inf[columns[i]], t_intp, dynamic_inf_intp[columns[i]])

# dynamic_inf_intp.to_excel('results/train_val_test_dynamic_influent_interpolated.xlsx')
# # dynamic_inf_intp.to_excel('results/evaluation_dynamic_influent_interpolated.xlsx')