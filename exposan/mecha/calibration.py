# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan.utils import ospath, time_printer
from exposan.mecha import (
    results_path,
    create_model,
    )

import numpy as np, pandas as pd
import optuna

# import time

__all__ = ('objective')

#%%

# start_time = time.time()
# print('Start time:', start_time)

# def time_track(t, y):
#     global start_time

#     elapsed_time = time.time() - start_time

#     return elapsed_time-720

#%%
mdl = create_model()

# time_track.terminal = True

@time_printer
def objective(trial):

    params = {
        'Y_A': trial.suggest_uniform('Y_A', 0.1, 1),
        'Y_H': trial.suggest_uniform('Y_H', 0.1, 1),
        'mu_H': trial.suggest_uniform('mu_H', 1, 10),
        'K_S': trial.suggest_uniform('K_S', 10, 50),
        'K_O_H': trial.suggest_uniform('K_O_H', 0.02, 2),
        'K_NO': trial.suggest_uniform('K_NO', 0.05, 5),
        'b_H': trial.suggest_uniform('b_H', 0.062, 6.2),
        'eta_g': trial.suggest_uniform('eta_g', 0.1, 1),
        'eta_h': trial.suggest_uniform('eta_h', 0.1, 1),
        'k_h': trial.suggest_uniform('k_h', 0.3, 30),
        'K_X': trial.suggest_uniform('K_X', 0.003, 0.3),
        'mu_A': trial.suggest_uniform('mu_A', 0.08, 8),
        'K_NH': trial.suggest_uniform('K_NH', 0.1, 10),
        'b_A': trial.suggest_uniform('b_A', 0.01, 1),
        'K_O_A': trial.suggest_uniform('K_O_A', 0.04, 4),
        'k_a': trial.suggest_uniform('k_a', 0.008, 0.8),
        }

    current_params = np.array([])

    for k, v in params.items():
        current_params = np.append(current_params, v)

    # global start_time
    # start_time = time.time()    # reset start_time to be here

    # print('Renewed start time:', start_time)

    try:
        mdl._update_state(current_params, t_span=(0, 50), t_eval = np.arange(0, 50, 51), method='BDF', state_reset_hook='reset_cache', print_t=False)
        # mdl._update_state(current_params, t_span=(0, 25), t_eval = np.arange(0, 26, 1), method='RK23', state_reset_hook='reset_cache', print_t=False, events=time_track)

    except:
        print('Fail & return 15')
        return 10000

    out = [metric() for metric in mdl.metrics]
    obj = np.average(out)

    print('Objective function:', obj)
    return min(obj, 10000.1)

if __name__ == '__main__':

    # sampler = optuna.samplers.TPESampler(seed=555)
    sampler = optuna.samplers.TPESampler(seed=777)

    study = optuna.create_study(sampler=sampler, direction='minimize')
    # study = optuna.create_study(sampler=sampler, direction='minimize', pruner=optuna.pruners.MedianPruner())
    # w/ or w/o lead to the same results & the same time required, but hyperparameter importances and optimization history changes

    study.optimize(objective, n_trials=300)

    print(study.best_params)

    df = study.trials_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] == 300

    df.to_excel(ospath.join(results_path, 'optuna_cali_seed777.xlsx'))
    # df.to_excel(ospath.join(results_path, 'conti_optuna_kill.xlsx'))

    # optuna.visualization.matplotlib.plot_param_importances(study)
    # optuna.visualization.matplotlib.plot_optimization_history(study)

    # optuna.visualization.matplotlib.plot_intermediate_values(study)

    # optuna.visualization.plot_param_importances(study).show(renderer='browser')
    # optuna.visualization.plot_optimization_history(study).show(renderer='browser')