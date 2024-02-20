# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 20:17:13 2023

@author: Ga-Yeong Kim
"""

from qsdsan.utils import ospath, load_data
from scipy.interpolate import interp1d

path = ospath.dirname(__file__)

data_path_model_prediction = ospath.join(path, 'data/model_prediction.xlsx')
data_path_observed = ospath.join(path, 'data/observed.xlsx')

model_prediction = load_data(data_path_model_prediction, sheet=0, index_col=None)
observed = load_data(data_path_observed, sheet=0, index_col=None)

t_model = model_prediction.t_stamp.to_numpy()
t_observed = observed.t_stamp.to_numpy()
t_observed_nh = observed.t_stamp_NH.to_numpy()
t_observed_p = observed.t_stamp_P.to_numpy()

vss_model = model_prediction.VSS.to_numpy()
vss_observed = observed.VSS.to_numpy()

ch_model = model_prediction.CH.to_numpy()
ch_observed = observed.CH.to_numpy()

li_model = model_prediction.LI.to_numpy()
li_observed = observed.LI.to_numpy()

nh_model = model_prediction.NH.to_numpy()
nh_observed = observed.NH.to_numpy()

p_model = model_prediction.P.to_numpy()
p_observed = observed.P.to_numpy()

#%%

f_vss = interp1d(t_observed, vss_observed)
vss_observed_interp = f_vss(t_model)
rmse_vss = (sum((vss_model - vss_observed_interp)**2)/len(t_model))**0.5

f_ch = interp1d(t_observed, ch_observed)
ch_observed_interp = f_ch(t_model)
rmse_ch = (sum((ch_model - ch_observed_interp)**2)/len(t_model))**0.5

f_li = interp1d(t_observed, li_observed)
li_observed_interp = f_li(t_model)
rmse_li = (sum((li_model - li_observed_interp)**2)/len(t_model))**0.5

f_nh = interp1d(t_observed_nh, nh_observed)
nh_observed_interp = f_nh(t_model)
rmse_nh = (sum((nh_model - nh_observed_interp)**2)/len(t_model))**0.5

f_p = interp1d(t_observed_p, p_observed)
p_observed_interp = f_p(t_model)
rmse_p = (sum((p_model - p_observed_interp)**2)/len(t_model))**0.5
#%%

print(f'RMSE_VSS: {round(rmse_vss, 3)}')
print(f'RMSE_CH: {round(rmse_ch, 3)}')
print(f'RMSE_LI: {round(rmse_li, 3)}')
print(f'RMSE_SNH: {round(rmse_nh, 3)}')
print(f'RMSE_SP: {round(rmse_p, 3)}')

# RMSE_VSS: 27.822
# RMSE_CH: 6.024
# RMSE_LI: 4.825
# RMSE_SNH: 1.863
# RMSE_SP: 0.011
