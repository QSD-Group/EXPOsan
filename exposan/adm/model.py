# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs
from chaospy import distributions as shape
from qsdsan.utils import ospath, time_printer
from exposan.adm import sys, _init_conds as _ic, AD, results_path
# import pandas as pd
import numpy as np
import os

__all__ = ('model_ss',)

#%%

# =============================================================================
# UA with random initial conditions to test steady state
# =============================================================================

model_ss = qs.Model(system=sys, exception_hook='raise')
param_ss = model_ss.parameter
# metric_ss = model_ss.metric

get_uniform_w_frac = lambda b, frac: shape.Uniform(lower=b*(1-frac), upper=b*(1+frac))


# Set initial conditions of all bioreactors
for k, v in _ic.items():
    b = v
    D = get_uniform_w_frac(b, 0.5)
    @param_ss(name='initial '+k, element=AD, kind='coupled', units='mg/L',
              baseline=b, distribution=D)
    def ic_setter(conc): pass


@time_printer
def run_wdiff_init(model, N, T, t_step, method='BDF', 
                   metrics_path='', timeseries_path='', 
                   rule='L', seed=None, pickle=False):
    if seed: np.random.seed(seed)
    samples = model.sample(N=N, rule=rule)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    if timeseries_path: tpath = timeseries_path
    else:
        folder = ospath.join(results_path, f'time_series_data_{seed}')
        os.mkdir(folder)
        tpath = os.path.join(folder, 'state.npy')
    for i, smp in enumerate(samples):
        concs = dict(zip(_ic.keys(), smp))
        AD.set_init_conc(**concs)          
        model._system.simulate(
            state_reset_hook='reset_cache',
            t_span=t_span,
            t_eval=t_eval,
            method=method,
            export_state_to=tpath,
            sample_id=i,
            )
    AD.set_init_conc(**_ic)

#%%
if __name__ == '__main__':
    seed = 119
    t = 200
    t_step = 5
    n = 100
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    run_wdiff_init(model_ss, n, t, t_step, method=method,
                   seed=seed)
