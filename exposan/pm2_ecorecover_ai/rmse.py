# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 20:17:13 2023

@author: Ga-Yeong Kim
"""

from qsdsan.utils import ospath, load_data
from scipy.interpolate import interp1d

import numpy as np

path = ospath.dirname(__file__)

data_path_dynamic_scada = ospath.join(path, 'data/dynamic_scada_result_calibration.xlsx')
data_path_dynamic_simul = ospath.join(path, 'data/dynamic_simul_result_calibration.xlsx')

# data_path_dynamic_scada = ospath.join(path, 'data/dynamic_scada_result_validation.xlsx')
# data_path_dynamic_simul = ospath.join(path, 'data/dynamic_simul_result_validation.xlsx')

result_scada = load_data(data_path_dynamic_scada, sheet=0, index_col=None)
result_simul = load_data(data_path_dynamic_simul, sheet=0, index_col=None)

t_scada = result_scada.t_stamp.to_numpy()
t_simul = result_simul.t_stamp.to_numpy()

# vss_scada = result_scada.VSS.to_numpy()
# vss_simul = result_simul.VSS.to_numpy()

sp_scada = result_scada.SP.to_numpy()
sp_simul = result_simul.SP.to_numpy()

# f_vss = interp1d(t_simul, vss_simul)
# vss_simul_interp = f_vss(t_scada)

# rmse_vss = (sum((vss_scada - vss_simul_interp)**2)/len(t_scada))**0.5

f_sp = interp1d(t_simul, sp_simul)
sp_simul_interp = f_sp(t_scada)

rmse_sp = (sum((sp_scada - sp_simul_interp)**2)/len(t_scada))**0.5

# print(f'RMSE_VSS: {round(rmse_vss, 3)}')
print(f'RMSE_SP: {round(rmse_sp, 3)}')

#%%

# t_scada_yield = result_scada.t_stamp_yield.to_numpy()
# t_simul_yield = result_simul.t_stamp_yield.to_numpy()

# t_scada_yield = t_scada_yield[~np.isnan(t_scada_yield)]
# t_simul_yield = t_simul_yield[~np.isnan(t_simul_yield)]

# yield_scada = result_scada.Yield_1d.to_numpy()
# yield_simul = result_simul.Yield_1d.to_numpy()

# yield_scada = yield_scada[~np.isnan(yield_scada)]
# yield_simul = yield_simul[~np.isnan(yield_simul)]

# f_yield = interp1d(t_simul_yield, yield_simul)
# yield_simul_interp = f_yield(t_scada_yield)

# rmse_yield = (sum((yield_scada - yield_simul_interp)**2)/len(t_scada_yield))**0.5

# print(f'RMSE_Yield_1d: {round(rmse_yield, 3)}')

# calidation results
# RMSE_VSS: 96.204
# RMSE_SP: 0.043
# RMSE_Yield_1d: 43.544

# validation results
# RMSE_VSS: 187.149
# RMSE_SP: 0.019
# RMSE_Yield_1d: 308.746

# Janus quotient
# VSS: 1.9453349133092181
# SP: 0.4418604651162791
# RMSE_Yield_1d: 7.090437