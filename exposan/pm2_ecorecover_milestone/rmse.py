# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 20:17:13 2023

@author: Ga-Yeong Kim
"""

from qsdsan.utils import ospath, load_data
from scipy.interpolate import interp1d

path = ospath.dirname(__file__)

# data_path_dynamic_scada = ospath.join(path, 'data/dynamic_scada_result_calibration_2.xlsx')
# data_path_dynamic_simul = ospath.join(path, 'data/dynamic_simul_result_calibration_2.xlsx')

data_path_dynamic_scada = ospath.join(path, 'data/dynamic_scada_result_validation_2.xlsx')
data_path_dynamic_simul = ospath.join(path, 'data/dynamic_simul_result_validation_2.xlsx')

result_scada = load_data(data_path_dynamic_scada, sheet=0, index_col=None)
result_simul = load_data(data_path_dynamic_simul, sheet=0, index_col=None)

t_scada = result_scada.t_stamp.to_numpy()
t_simul = result_simul.t_stamp.to_numpy()

vss_scada = result_scada.VSS.to_numpy()
vss_simul = result_simul.VSS.to_numpy()

sp_scada = result_scada.SP.to_numpy()
sp_simul = result_simul.SP.to_numpy()

sn_scada = result_scada.SN.to_numpy()
sn_simul = result_simul.SN.to_numpy()

f_vss = interp1d(t_simul, vss_simul)
vss_simul_interp = f_vss(t_scada)

rmse_vss = (sum((vss_scada - vss_simul_interp)**2)/len(t_scada))**0.5

f_sp = interp1d(t_simul, sp_simul)
sp_simul_interp = f_sp(t_scada)

rmse_sp = (sum((sp_scada - sp_simul_interp)**2)/len(t_scada))**0.5

f_sn = interp1d(t_simul, sn_simul)
sn_simul_interp = f_sn(t_scada)

rmse_sn = (sum((sn_scada - sn_simul_interp)**2)/len(t_scada))**0.5

print(f'RMSE_VSS: {round(rmse_vss, 3)}')
print(f'RMSE_SP: {round(rmse_sp, 3)}')
print(f'RMSE_SN: {round(rmse_sn, 3)}')


# calidation results
# RMSE_VSS: 103.139
# RMSE_SP: 0.006
# RMSE_SN: 4.83

# validation results
# RMSE_VSS: 74.905
# RMSE_SP: 0.011
# RMSE_SN: 6.113

# Janus quotient
# VSS: 0.726
# SP: 1.833
# SN: 1.266