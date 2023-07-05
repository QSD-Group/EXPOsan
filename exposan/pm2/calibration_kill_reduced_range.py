# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan.utils import ospath, time_printer
from exposan.pm2 import (
    results_path,
    create_model,
    sensitive_params,
    )

import numpy as np, pandas as pd
from scipy.optimize import shgo

# from datetime import datetime
import time

__all__ = ('cali_setup', 'optimizer', 'objective_function')

#%%

mdl = create_model()

def cali_setup():
    params = []
    baseline = np.array([])
    boundary = []

    for k, v in sensitive_params.items():
        b, units, bounds = v
        params.append(k)
        baseline = np.append(baseline, b)
        boundary.append(bounds)
    boundary = tuple(boundary)

    return params, baseline, boundary

sense_params, opt_params, bnds = cali_setup()         # start from init_guess

'''
sense_params: list of sensitive parameters
opt_params: initial guess of sensitive parameters
bnds: min & max of sensitive parameters
'''

#%%
def optimizer():

    opt = shgo(objective_function, bounds=bnds, iters=1, minimizer_kwargs={'method':'SLSQP', 'ftol':1e-3})
    # opt = shgo(objective_function, bounds=bnds, iters=1, minimizer_kwargs={'method':'SLSQP', 'ftol':1e-3}, workers=-1)

    opt_as_series = pd.Series(opt)
    opt_as_series.to_excel(excel_writer=(ospath.join(results_path, 'calibration_result_kill_reduced_range.xlsx')))

#%%

start_time = time.time()
print('Start time:', start_time)

def time_track(t, y):
    global start_time

    elapsed_time = time.time() - start_time

    return elapsed_time-720

#%%
time_track.terminal = True

@time_printer
def objective_function(opt_params, *args):
    global start_time
    start_time = time.time()    # reset start_time to be here

    print('Renewed start time:', start_time, '--------------------------------')
    print('Parameters:', opt_params)

    try:
        # mdl._update_state(opt_params, t_span=(0, 50), t_eval = np.arange(0, 51, 1), method='RK23', state_reset_hook='reset_cache', print_t=True)
        mdl._update_state(opt_params, t_span=(0, 50), t_eval = np.arange(0, 51, 1), method='RK23', state_reset_hook='reset_cache', print_t=False, events=time_track)

        # when terminated, what to return?

    except:
        print('Fail & return 5.0')
        return 5

    out = [metric() for metric in mdl.metrics]
    obj = np.average(out)

    print('Objective function:', obj)
    return obj

optimizer()