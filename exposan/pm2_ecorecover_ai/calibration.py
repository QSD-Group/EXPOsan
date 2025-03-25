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
from exposan.pm2_ecorecover import (
    results_path,
    create_model,
    )

import numpy as np, pandas as pd
import optuna

import time

__all__ = ('objective')

#%%

start_time = time.time()
print('Start time:', start_time)

def time_track(t, y):
    global start_time

    elapsed_time = time.time() - start_time

    return elapsed_time-720

#%%
mdl = create_model()

time_track.terminal = True

@time_printer
def objective(trial):

    params = {
        'arr_e': trial.suggest_uniform('arr_e', 1000, 10000),
        'K_P': trial.suggest_uniform('K_P', 0.01, 100),
        'f_CH_max': trial.suggest_uniform('f_CH_max', 0.1, 10),
        'exponent': trial.suggest_uniform('exponent', 1, 10)
        }

    current_params = np.array([])

    for k, v in params.items():
        current_params = np.append(current_params, v)

    global start_time
    start_time = time.time()    # reset start_time to be here

    print('Renewed start time:', start_time)

    try:
        mdl._update_state(current_params, t_span=(0, 25), t_eval = np.arange(0, 26, 1), method='RK23', state_reset_hook='reset_cache', print_t=False, events=time_track)

    except:
        print('Fail & return 15')
        return 15

    out = [metric() for metric in mdl.metrics]
    obj = np.average(out)

    print('Objective function:', obj)
    return min(obj, 15.1)

if __name__ == '__main__':

    # sampler = optuna.samplers.TPESampler(seed=555)
    sampler = optuna.samplers.TPESampler(seed=777)

    study = optuna.create_study(sampler=sampler, direction='minimize')
    # study = optuna.create_study(sampler=sampler, direction='minimize', pruner=optuna.pruners.MedianPruner())
    # w/ or w/o lead to the same results & the same time required, but hyperparameter importances and optimization history changes

    study.optimize(objective, n_trials=500)

    print(study.best_params)

    df = study.trials_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] == 500

    df.to_excel(ospath.join(results_path, 'optuna_cali_seed777.xlsx'))
    # df.to_excel(ospath.join(results_path, 'conti_optuna_kill.xlsx'))

    # optuna.visualization.matplotlib.plot_param_importances(study)
    # optuna.visualization.matplotlib.plot_optimization_history(study)

    # optuna.visualization.matplotlib.plot_intermediate_values(study)

    # optuna.visualization.plot_param_importances(study).show(renderer='browser')
    # optuna.visualization.plot_optimization_history(study).show(renderer='browser')