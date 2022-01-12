# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


from exposan.bsm1 import model_bsm1 as mdl
import numpy as np

np.random.seed(112)

samples = mdl.sample(N=10, rule='L')
mdl.load_samples(samples)
t = 10
t_step = 1
mdl.evaluate(t_span=(0,t),
             t_eval=np.arange(0, t+t_step, t_step),
             method='LSODA',
             export_state_to='results/test_time_series.xlsx')
mdl.table.to_excel('results/table.xlsx')