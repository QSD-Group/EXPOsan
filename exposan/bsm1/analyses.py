#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

'''
#TODO
Other things that I think we might want to include:
    - A figure in the SI (or just on GitHub) showing that we can converge to
    similar steady-states conditions
'''

# from qsdsan.utils import ords
from time import time
from exposan.bsm1 import model_bsm1 as mdl, \
    run_uncertainty#, analyze_timeseries
import os
# import numpy as np

bsm1_path = os.path.dirname(__file__)
N = 100
T = 50
t_step = 1

#%%
if __name__ == '__main__':
    for i in range(10):
        seed = int(str(time())[-3:])
        mt_path = os.path.join(bsm1_path, f'results/table_{seed}.xlsx')
        ts_path = os.path.join(bsm1_path, f'results/states_{seed}.xlsx')
        print(f'\nBatch ID: {seed}')
        print('-'*15)
        run_uncertainty(mdl, N, T, t_step, method='BDF', 
                        metrics_path=mt_path, timeseries_path=ts_path, 
                        seed=seed)
# out = analyze_timeseries(SE_var_getter)

