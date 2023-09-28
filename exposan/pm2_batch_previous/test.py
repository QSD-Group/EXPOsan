#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 22:00:26 2023

@author: macbook
"""


import os, numpy as np, qsdsan as qs
from qsdsan import processes as pc
from scipy.interpolate import interp1d
# from warnings import warn
from chaospy import distributions as shape
from qsdsan.utils import DictAttrSetter, ospath, time_printer, load_data

from exposan.pm2_batch import (
    create_system,
    data_path,
    results_path,
    )



#%%
file = ospath.join(data_path, 'batch_exp_result.xlsx')

result_exp = load_data(file, sheet=0, index_col=None)


t_exp = result_exp.t_stamp.to_numpy()

vss_exp = result_exp.VSS.to_numpy()
snh_exp = result_exp.S_NH.to_numpy()
sp_exp = np.maximum(result_exp.S_P.to_numpy(), 0.015)

vss_range = max(vss_exp) - min(vss_exp)
snh_range = max(snh_exp) - min(snh_exp)
sp_range = max(sp_exp) - min(sp_exp)


#%%



file_simul = ospath.join(data_path, 'test.xlsx')

result_simul = load_data(file_simul, sheet=0, index_col=None)

t_simul = result_simul.t_stamp.to_numpy()

vss = result_simul.VSS.to_numpy()
snh = result_simul.S_NH.to_numpy()
sp = np.maximum(result_simul.S_P.to_numpy(), 0.015)

f = interp1d(t_simul, vss)
vss_simul = f(t_exp)

rmse_vss = (sum((vss_exp - vss_simul)**2)/len(t_exp))**0.5
nrmse_vss = rmse_vss/vss_range


f2 = interp1d(t_simul, snh)
snh_simul = f2(t_exp)

rmse_snh = (sum((snh_exp - snh_simul)**2)/len(t_exp))**0.5
nrmse_snh = rmse_snh/snh_range

f3= interp1d(t_simul, sp)
sp_simul = f3(t_exp)

rmse_sp = (sum((sp_exp - sp_simul)**2)/len(t_exp))**0.5
nrmse_sp = rmse_sp/sp_range

avg = (nrmse_vss+nrmse_snh+nrmse_sp)/3

print(avg)

#%%

file = ospath.join(data_path, 'batch_exp_result_chli.xlsx')

result_exp = load_data(file, sheet=0, index_col=None)

t_exp_x = result_exp.t_stamp_x.to_numpy()
t_exp_nh = result_exp.t_stamp_nh.to_numpy()
t_exp_p = result_exp.t_stamp_p.to_numpy()

t_exp_x=t_exp_x[~np.isnan(t_exp_x)]
t_exp_p=t_exp_p[~np.isnan(t_exp_p)]

vss_exp = result_exp.VSS.to_numpy()
ch_exp = result_exp.CH.to_numpy()
li_exp = result_exp.LI.to_numpy()
snh_exp = result_exp.S_NH.to_numpy()
sp_exp = np.maximum(result_exp.S_P.to_numpy(), 0.015)

vss_exp=vss_exp[~np.isnan(vss_exp)]
ch_exp=ch_exp[~np.isnan(ch_exp)]
li_exp=li_exp[~np.isnan(li_exp)]
sp_exp=sp_exp[~np.isnan(sp_exp)]

vss_range = max(vss_exp) - min(vss_exp)
snh_range = max(snh_exp) - min(snh_exp)
sp_range = max(sp_exp) - min(sp_exp)
ch_range = max(ch_exp) - min(ch_exp)
li_range = max(li_exp) - min(li_exp)



file_simul = ospath.join(data_path, 'test_include.xlsx')

result_simul = load_data(file_simul, sheet=0, index_col=None)

t_simul = result_simul.t_stamp.to_numpy()

vss = result_simul.VSS.to_numpy()
snh = result_simul.S_NH.to_numpy()
sp = np.maximum(result_simul.S_P.to_numpy(), 0.015)
ch = result_simul.X_CH.to_numpy()
li = result_simul.X_LI.to_numpy()

f = interp1d(t_simul, vss)
vss_simul = f(t_exp_x)

rmse_vss = (sum((vss_exp - vss_simul)**2)/len(t_exp_x))**0.5
nrmse_vss = rmse_vss/vss_range





f2 = interp1d(t_simul, snh)
snh_simul = f2(t_exp_nh)

rmse_snh = (sum((snh_exp - snh_simul)**2)/len(t_exp_nh))**0.5
nrmse_snh = rmse_snh/snh_range




f3= interp1d(t_simul, sp)
sp_simul = f3(t_exp_p)

rmse_sp = (sum((sp_exp - sp_simul)**2)/len(t_exp_p))**0.5
nrmse_sp = rmse_sp/sp_range



f4 = interp1d(t_simul, ch)
ch_simul = f4(t_exp_x)

rmse_ch = (sum((ch_exp - ch_simul)**2)/len(t_exp_x))**0.5
nrmse_ch = rmse_ch/ch_range



f5 = interp1d(t_simul, li)
li_simul = f5(t_exp_x)

rmse_li = (sum((li_exp - li_simul)**2)/len(t_exp_x))**0.5
nrmse_li = rmse_li/li_range

avg = (nrmse_vss+nrmse_snh+nrmse_sp+nrmse_ch+nrmse_li)/5

print(avg)