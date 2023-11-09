# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 20:17:13 2023

@author: Ga-Yeong Kim
"""
import numpy as np

from qsdsan.utils import ospath, load_data
from scipy.interpolate import interp1d

path = ospath.dirname(__file__)

# data_path_scada = ospath.join(path, 'data/yield_scada_calibration.xlsx')
# data_path_simul = ospath.join(path, 'data/yield_simul_calibration.xlsx')

data_path_scada = ospath.join(path, 'data/yield_scada_validation.xlsx')
data_path_simul = ospath.join(path, 'data/yield_simul_validation.xlsx')

yield_scada = load_data(data_path_scada, sheet=0, index_col=None)
yield_simul = load_data(data_path_simul, sheet=0, index_col=None)

yield_scada[np.isnan(yield_scada)]=0
yield_simul[np.isnan(yield_simul)]=0

t_scada = yield_scada.t_stamp.to_numpy()
t_simul = yield_simul.t_stamp.to_numpy()

raw_scada = np.maximum(yield_scada.raw.to_numpy(), 0.01)
raw_simul = np.maximum(yield_simul.raw.to_numpy(), 0.01)

five_scada = np.maximum(yield_scada.five.to_numpy(), 0.01)
five_simul = np.maximum(yield_simul.five.to_numpy(), 0.01)

three_scada = np.maximum(yield_scada.three.to_numpy(), 0.01)
three_simul = np.maximum(yield_simul.three.to_numpy(), 0.01)
# three_scada[np.isnan(three_scada)]=0
# three_simul[np.isnan(three_simul)]=0

one_scada = np.maximum(yield_scada.one.to_numpy(), 0.01)
one_simul = np.maximum(yield_simul.one.to_numpy(), 0.01)
# one_scada[np.isnan(one_scada)]=0
# one_simul[np.isnan(one_simul)]=0
#%%

f_raw = interp1d(t_simul, raw_simul)
raw_simul_interp = f_raw(t_scada)

mape_raw = sum(np.abs(raw_scada - raw_simul_interp)/np.abs(raw_scada))/len(t_scada) * 100

f_five = interp1d(t_simul, five_simul)
five_simul_interp = f_five(t_scada)

mape_five = sum(np.abs(five_scada - five_simul_interp)/np.abs(five_scada))/len(t_scada) * 100

f_three = interp1d(t_simul, three_simul)
three_simul_interp = f_three(t_scada)

mape_three = sum(np.abs(three_scada - three_simul_interp)/np.abs(three_scada))/len(t_scada) * 100

f_one = interp1d(t_simul, one_simul)
one_simul_interp = f_one(t_scada)

mape_one = sum(np.abs(one_scada - one_simul_interp)/np.abs(one_scada))/len(t_scada) * 100

print(f'MAPE_RAW: {round(mape_raw, 3)}')
print(f'MAPE_5d: {round(mape_five, 3)}')
print(f'MAPE_3d: {round(mape_three, 3)}')
print(f'MAPE_1d: {round(mape_one, 3)}')

# calibration results
# MAPE_RAW: 58.274
# MAPE_5d: 18.501
# MAPE_3d: 21.186
# MAPE_1d: 32.653

# validation results
# MAPE_RAW: 48.792
# MAPE_5d: 31.255
# MAPE_3d: 35.821
# MAPE_1d: 39.208