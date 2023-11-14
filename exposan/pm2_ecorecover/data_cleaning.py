# -*- coding: utf-8 -*-
'''
Created on Fri Nov  4 19:45:06 2022

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''
import matplotlib.pyplot as plt, matplotlib.ticker as ticker, \
    pandas as pd, numpy as np
from matplotlib import dates
from qsdsan.utils import load_data, ospath

from skfda.misc.hat_matrix import (KNeighborsHatMatrix,
                                   LocalLinearRegressionHatMatrix,
                                   NadarayaWatsonHatMatrix,)

from skfda.preprocessing.smoothing import KernelSmoother
from skfda.representation.grid import FDataGrid

#%% Load raw data
folder = ospath.dirname(__file__)
file = ospath.join(folder, 'data/data_cleaning_rawdata.xlsx')
raw_data = load_data(file, sheet=None)
scada = raw_data['SCADA']

# Join AIMS data columns
aims = raw_data['AIMS_Inf_NH4'].join(raw_data['AIMS_Eff_NH4'], how='outer', sort=True)  # initial
for k, df in raw_data.items():
    if k.startswith('AIMS_') and not k.endswith('NH4'):
        aims = aims.join(df, how='outer', sort=True)    # append
del raw_data

# Delete '_aims' in column names
aims.columns = [c.rstrip('_aims') for c in aims.columns]

# Set starting point
from datetime import datetime
day0 = datetime.strptime('03-12-2022 12:00 AM', '%m-%d-%Y %I:%M %p')

#%% Difference between SCADA and AIMS
all_data = scada.join(aims, how='right', rsuffix='_aims', sort=True)
sheets = ['nh4_inf', 'nh4_eff', 'no3_eff', 'po4_inf', 'po4_eff', 'tss']
dct_compare = {}
for sheet in sheets:
    cols = [col for col in all_data.columns if col.startswith(sheet)]
    df = all_data.loc[:,cols]
    if sheet == 'tss':
        for i in ('_harv', '_postmix'):
            df[sheet+i+'_diff'] = df[sheet+i] - df[sheet+'_ret']
            dct_compare[sheet+i] = df.loc[:, [sheet+'_ret', sheet+i, sheet+i+'_diff']].dropna(axis=0)
    else:
        df[sheet+'_diff'] =  df[sheet+'_aims'] - df[sheet]
        dct_compare[sheet] = df.dropna(axis=0)

with pd.ExcelWriter(ospath.join(folder, 'results/data_cleaning/excel/scada_vs_aims.xlsx')) as writer:
    for sheet, df in dct_compare.items():
        df.to_excel(writer, sheet_name=sheet)

#%% Plotting functions

plot_path = ospath.join(folder, 'results/data_cleaning/plot/')
date0 = dates.date2num(day0)

@ticker.FuncFormatter
def fmt_day(t_stamp, pos):
    day = t_stamp - date0
    return f'{day:.0f}'

def plot_single(x, y, x_by_date=True, file_name=None):
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.plot(x, y, label=y.name)
    ax.xaxis.set_major_locator(dates.DayLocator(interval=5))
    fig.legend(loc='upper right', bbox_to_anchor=(1,0.9))
    if x_by_date:
        ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%m-%d'))
        ax.set_xlabel('Date')
        subfolder = 'raw/date'
    else:
        ax.xaxis.set_major_formatter(fmt_day)
        ax.set_xlabel('Time [day]')
        subfolder = 'raw/day'
    if file_name:
        fig.savefig(ospath.join(plot_path, f'{subfolder}/{file_name}.png'),
                    dpi=300, bbox_inches='tight', facecolor='white')
    else: return fig, ax

def plot_SCADA(x_by_date=True):
    for col, y in scada.items():
        plot_single(scada.index, y, x_by_date=x_by_date, file_name=col)

# plot_SCADA(x_by_date=True)    # date
# plot_SCADA(x_by_date=False)   # day

def plot_multiple(df1, df2, labels=None, smooth=False, x_by_date=True, file_name=None):
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.plot(df1)
    ax.plot(df2, 'o', ms=4, linestyle='-')
    labels = labels or [*df1.columns, *[col+'_aims' for col in df2.columns]]
    ax.xaxis.set_major_locator(dates.DayLocator(interval=5))
    if x_by_date:
        ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%m-%d'))
        ax.set_xlabel('Date')
        file = f'raw/date/{file_name}_vs.png'
    else:
        ax.xaxis.set_major_formatter(fmt_day)
        ax.set_xlabel('Time [day]')
        file = f'raw/day/{file_name}_vs.png'
    fig.legend(labels=labels, loc='upper right', bbox_to_anchor=(1,0.9))
    if file_name:
        fig.savefig(ospath.join(plot_path, file),
                    dpi=300, bbox_inches='tight', facecolor='white')
    else: return fig, ax

def plot_SCADA_vs_AIMS(x_by_date=True):
    for var in set(scada.columns) & set(aims.columns):
        df1 = scada.loc[:, [var]]
        df2 = aims.loc[:, [var]].dropna(axis=0)
        plot_multiple(df1, df2, x_by_date=x_by_date, file_name=var)
    df1 = scada.loc[:,['tss_ret']]
    df2 = aims.loc[:,['tss_harv', 'tss_postmix']]                   # this should be separated? or masked?
    labels = [*df1.columns, *df2.columns]
    plot_multiple(df1, df2, labels, x_by_date, file_name='tss')

#%% To export plots
for by in (True, False):
    plot_SCADA(by)
    plot_SCADA_vs_AIMS(by)

#%% Kernel smoothing function

# Calculate t_stamp as day
calc_day = lambda t_stamp: (t_stamp-day0).total_seconds()/60/60/24

def plot_raw_vs_smooth(x1, y1, x2, y2, save_as=''):
    fig, ax = plt.subplots(figsize=(15, 5))
    l1, = ax.plot(x1, y1, label='Raw data')
    l2, = ax.plot(x2, y2, label='Smoothed')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.set_xlabel('Time [day]')
    ax.legend(handles=[l1,l2], title=save_as)
    path = ospath.join(plot_path, f'smoothed/{save_as}.png')
    fig.savefig(path, dpi=300, facecolor='white')

def kernel_smoothing(data, method, kernel, bw, file_name=None, **kwargs):
    '''
    Smooths input data with kernel smoother.

    Parameters
    ----------
    data : pandas.DataFrame or pandas.Series
        Raw data.
    method : str
        One of ('knn', 'llr', 'nw').
    kernel : Callable
        Kernel function.
    bw : float or int
        Bandwidth or the number of neighbor to consider
    **kwargs : optional
        Other parameters passed to `KernelSmoother`.
    '''
    # To use skfda package: 1) pandas series -> array,
    #                       2) float64 -> float32,
    #                       3) use 'day' data, not 'date' data

    x = pd.Series(data.index).apply(calc_day).to_numpy()

    if method == 'knn': kernel_estimator = KNeighborsHatMatrix(n_neighbors=bw, kernel=kernel)
    elif method == 'llr': kernel_estimator = LocalLinearRegressionHatMatrix(bandwidth=bw, kernel=kernel)
    else: kernel_estimator = NadarayaWatsonHatMatrix(bandwidth=bw, kernel=kernel)

    smoother = KernelSmoother(kernel_estimator=kernel_estimator, **kwargs)
    xhat = kwargs.pop('output_points', None) or x               # in case times points other than the raw data are specified as output_points (output time point of the smoothed data doesn't have to be exactly the same as the raw data; If not specified, then it'll be the same.)
    out = pd.DataFrame({'day':xhat})
    if len(data.shape) == 1:
        y = data.to_numpy()
        fd = FDataGrid(y, x)                                    # creation of FDataGrid
        fd_smoothed = smoother.fit_transform(fd)                # smooth
        yhat = pd.Series(fd_smoothed.data_matrix.flatten())     # smoothed output
        if xhat is x:
            out['raw'] = y
        out['smoothed'] = yhat
        plot_raw_vs_smooth(x, y, xhat, yhat, f'{data.name}_{method}')
        out.to_excel(ospath.join(folder, f'results/data_cleaning/excel/smoothed/{data.name}_{method}.xlsx'))
    else:                                                       # multiple y
        for col, y in data.items():
            fd = FDataGrid(y.to_numpy(), x)
            fd_smoothed = smoother.fit_transform(fd)
            yhat = pd.Series(fd_smoothed.data_matrix.flatten())
            out[col+'_smoothed'] = yhat
            plot_raw_vs_smooth(x, y, xhat, yhat, f'{col}_{method}')
        out.to_excel(ospath.join(folder, f'results/data_cleaning/excel/smoothed/{file_name}_{method}.xlsx'))

#%% Kernel smoothing commands
from skfda.misc.kernels import (
    normal,
    # cosine,
    # epanechnikov,
    # tri_weight,
    # quartic,
    # uniform,
    )

#%% SCADA data smoothing (need to check one by one manually)

# Flow & volume
kernel_smoothing(data=scada.flow_mix, method='llr', kernel=normal, bw=0.1)
kernel_smoothing(data=scada.flow_pbr, method='llr', kernel=normal, bw=0.5)
kernel_smoothing(data=scada.flow_mem, method='llr', kernel=normal, bw=0.1)
kernel_smoothing(data=scada.flow_ret, method='llr', kernel=normal, bw=0.1)
kernel_smoothing(data=scada.flow_eff, method='llr', kernel=normal, bw=0.1)
# kernel_smoothing(data=scada.flow_cen, method='llr', kernel=normal, bw=0.1)    # use as it is

kernel_smoothing(data=scada.vol_mix, method='llr', kernel=normal, bw=0.2)
kernel_smoothing(data=scada.vol_mem, method='llr', kernel=normal, bw=0.5)
kernel_smoothing(data=scada.vol_ret, method='llr', kernel=normal, bw=0.2)

# Temperature & light
kernel_smoothing(data=scada.temp, method='llr', kernel=normal, bw=0.1)
kernel_smoothing(data=scada.light, method='llr', kernel=normal, bw=0.001)

# Aqueous species
kernel_smoothing(data=scada.do, method='llr', kernel=normal, bw=0.001)
kernel_smoothing(data=scada.nh4_inf, method='llr', kernel=normal, bw=0.5)
kernel_smoothing(data=scada.nh4_eff, method='llr', kernel=normal, bw=0.5)
kernel_smoothing(data=scada.no3_eff, method='llr', kernel=normal, bw=0.8)
kernel_smoothing(data=scada.po4_inf, method='llr', kernel=normal, bw=0.5)
kernel_smoothing(data=scada.po4_eff, method='llr', kernel=normal, bw=0.5)
kernel_smoothing(data=scada.tss_ret, method='llr', kernel=normal, bw=0.7)

#%% SCADA vs. AIMS differences smoothing (need to check one by one manually)

var = ['nh4_inf', 'nh4_eff', 'no3_eff', 'po4_inf', 'po4_eff',
       'tss_harv', 'tss_postmix']

kernel_smoothing(data=dct_compare[var[0]][var[0]+'_diff'], method='llr', kernel=normal, bw=0.8)
kernel_smoothing(data=dct_compare[var[1]][var[1]+'_diff'], method='llr', kernel=normal, bw=0.6)
kernel_smoothing(data=dct_compare[var[2]][var[2]+'_diff'], method='llr', kernel=normal, bw=0.8)
kernel_smoothing(data=dct_compare[var[3]][var[3]+'_diff'], method='llr', kernel=normal, bw=1)
kernel_smoothing(data=dct_compare[var[4]][var[4]+'_diff'], method='llr', kernel=normal, bw=1)
kernel_smoothing(data=dct_compare[var[5]][var[5]+'_diff'], method='llr', kernel=normal, bw=1)
kernel_smoothing(data=dct_compare[var[6]][var[6]+'_diff'], method='llr', kernel=normal, bw=2)

#%% Interpolate diff & correct smoothed SCADA

t_stamp = scada.index
var = ['nh4_inf', 'nh4_eff', 'no3_eff', 'po4_inf', 'po4_eff', 'tss_harv', 'tss_postmix']

def plot_raw_vs_corrected(df1, df2, df3, save_as=''):
    fig, ax = plt.subplots(figsize=(15, 5))
    l1, = ax.plot(df1, label='Raw SCADA')
    l2, = ax.plot(df2, label='Raw AIMS')
    l3, = ax.plot(df3, label='Corrected')
    ax.xaxis.set_major_locator(dates.DayLocator(interval=5))
    ax.xaxis.set_major_formatter(fmt_day)
    ax.set_xlabel('Time [day]')
    ax.legend(handles=[l1,l2,l3], title=save_as)
    path = ospath.join(plot_path, f'corrected/{save_as}.png')
    fig.savefig(path, dpi=300, facecolor='white')

def correct_scada(t_stamp, var, smoothing_method='llr'):
    t_intp = calc_day(t_stamp).to_numpy()
    diff_intp = pd.DataFrame()
    diff_intp.index = t_stamp
    diff_intp['t_interp'] = t_intp
    with pd.ExcelWriter(ospath.join(folder, 'results/data_cleaning/excel/scada_corrected.xlsx')) as writer:
        for i in var:
            if i == 'tss_harv':
                pass
            else:
                diff_smoothed = load_data(ospath.join(folder, f'results/data_cleaning/excel/smoothed/{i}_diff_{smoothing_method}.xlsx'))
                diff_smoothed_t = diff_smoothed.day
                diff_smoothed_df = diff_smoothed.smoothed
                diff_intp[i+'diff_smo_interp'] = diff = np.interp(t_intp, diff_smoothed_t, diff_smoothed_df)
                if i == 'tss_postmix':
                    scada_smoothed = load_data(ospath.join(folder, f'results/data_cleaning/excel/smoothed/tss_ret_{smoothing_method}.xlsx'))
                    scada_smoothed.index = t_stamp
                    corrected = scada_smoothed[i+'_corr'] = scada_smoothed.smoothed + diff
                    corrected[corrected < 0] = 0
                    scada_smoothed[i+'_corr'] = corrected
                    scada_smoothed.to_excel(writer, sheet_name='tss_ret')
                    raw_scada = scada.loc[:,'tss_ret']
                else:
                    scada_smoothed = load_data(ospath.join(folder, f'results/data_cleaning/excel/smoothed/{i}_{smoothing_method}.xlsx'))
                    scada_smoothed.index = t_stamp
                    corrected = scada_smoothed[i+'_corr'] = scada_smoothed.smoothed + diff
                    corrected[corrected < 0] = 0
                    scada_smoothed[i+'_corr'] = corrected
                    scada_smoothed.to_excel(writer, sheet_name=i)
                    raw_scada = scada.loc[:,i]
                raw_aims = aims.loc[:,i].dropna(axis=0)
                plot_raw_vs_corrected(raw_scada, raw_aims, corrected, i)
    diff_intp.to_excel(ospath.join(folder, 'results/data_cleaning/excel/interpolated_smoothed_diff.xlsx'))

correct_scada(t_stamp, var, smoothing_method='llr')

#%% HRT of smoothed SCADA data

flow = ['flow_mix', 'flow_pbr', 'flow_mem', 'flow_ret', 'flow_cen']
vol = ['vol_mix', 'vol_pbr', 'vol_mem', 'vol_ret']

def hrt(smoothing_method='llr'):
    with pd.ExcelWriter(ospath.join(folder, 'results/data_cleaning/excel/hrt.xlsx')) as writer:
        for i in range(4):
            flow_smoothed = load_data(ospath.join(folder, f'results/data_cleaning/excel/smoothed/{flow[i]}_{smoothing_method}.xlsx'))
            if i == 1:
                v = 77.49   # fixed
            else:
                v = load_data(ospath.join(folder, f'results/data_cleaning/excel/smoothed/{vol[i]}_{smoothing_method}.xlsx')).smoothed
            flow_smoothed[vol[i]] = v
            flow_smoothed['hrt'] = v/flow_smoothed.smoothed
            flow_smoothed.to_excel(writer, sheet_name=flow[i])

            plt.hist(flow_smoothed['hrt'], bins=15)
            plt.show()

hrt(smoothing_method='llr')

#%% Split ratio

def ratio(smoothing_method='llr'):
    flow_pbr_smoothed = load_data(ospath.join(folder, f'results/data_cleaning/excel/smoothed/flow_pbr_{smoothing_method}.xlsx')).smoothed
    flow_mem_smoothed = load_data(ospath.join(folder, f'results/data_cleaning/excel/smoothed/flow_mem_{smoothing_method}.xlsx')).smoothed
    flow_eff_smoothed = load_data(ospath.join(folder, f'results/data_cleaning/excel/smoothed/flow_eff_{smoothing_method}.xlsx')).smoothed
    flow_cen = scada['flow_cen']
    flow_cen.index = flow_mem_smoothed.index

    ratio_pbr_to_mem = flow_mem_smoothed/flow_pbr_smoothed                    # ME = 0.51
    ratio_pbr_to_ret = 1 - ratio_pbr_to_mem                                   # INT = 0.49

    ratio_mem_to_eff = flow_eff_smoothed/flow_mem_smoothed                    # TE = 0.39
    ratio_mem_to_postmem = 1 - ratio_mem_to_eff                               # RETEN = 0.61

    ratio_postmem_to_cen = flow_cen/(flow_mem_smoothed - flow_eff_smoothed)   # CE = 0.03
    ratio_postmem_to_ret = 1 - ratio_postmem_to_cen                           # RE = 0.97

    ratios = pd.DataFrame({'t_stamp_day': calc_day(t_stamp).to_numpy(),
                          'ratio_pbr_to_mem (ME)':ratio_pbr_to_mem,
                          'ratio_pbr_to_ret (INT)':ratio_pbr_to_ret,
                          'ratio_mem_to_eff (TE)':ratio_mem_to_eff,
                          'ratio_mem_to_postmem (RETEN)':ratio_mem_to_postmem,
                          'ratio_postmem_to_cen (CE)':ratio_postmem_to_cen,
                          'ratio_postmem_to_ret (RE)':ratio_postmem_to_ret})

    ratios.to_excel(ospath.join(folder, 'results/data_cleaning/excel/ratio.xlsx'))

    for i in range(6):
        plt.hist(ratios.iloc[:,i+1], bins=15)
        plt.show()

ratio(smoothing_method='llr')