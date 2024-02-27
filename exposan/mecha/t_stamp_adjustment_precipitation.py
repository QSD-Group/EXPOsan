# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 16:44:59 2024

@author: Ga-Yeong Kim
"""

# Import packages
import os, pandas as pd, numpy as np, matplotlib.pyplot as plt
# from matplotlib import dates
from qsdsan.utils import load_data
from exposan.mecha import data_path

#%% Load precipitation data

# precipitation_file = os.path.join(data_path, 'train_val_test_rijen_precipitation.xlsx')       # Time ['days'], irregular interval
# precipitation = load_data(precipitation_file, sheet='Rijen_Precipitation')

precipitation_file = os.path.join(data_path, 'evaluation_rijen_precipitation.xlsx')         # evaluation dataset
precipitation = load_data(precipitation_file, sheet='Rijen_Precipitation')

#%% Load online data (time stamp I'd like to have)

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

t_stamp_precipitation = precipitation.index
t_precipitation = calc_day(t_stamp_precipitation).to_numpy()
#%%

def plot_raw_vs_intp(x1, y1, x2, y2):
    fig, ax = plt.subplots(figsize=(15, 3))                                          # figure size
    l1, = ax.plot(x1, y1, label='Raw')                                               # first plot = raw data
    l2, = ax.plot(x2, y2, label='Interpolated')                                      # second plot = interpolated data
    ax.legend(handles=[l1,l2])
    return fig, ax

#%%

online['precipitation'] = np.interp(t_intp, t_precipitation, precipitation['precipitation'] )
plot_raw_vs_intp(t_precipitation, precipitation['precipitation'], t_intp, online['precipitation'])

# online.to_excel('results/train_val_test_online_dynamic_influent_integrated_precipitation.xlsx')
online.to_excel('results/evaluation_online_dynamic_influent_integrated_precipitation.xlsx')