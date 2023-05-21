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
#%% Import everything for now

from qsdsan.utils import ospath, time_printer
from exposan.pm2 import (
    results_path,
    create_model,
    sensitive_params,
    )

import numpy as np, pandas as pd
import optuna

# import winsound as sd
# from datetime import datetime

__all__ = ('objective')

# start_time = datetime.now()
# print('Started (hh:mm:ss.ms) {}'.format(start_time))

#%%

mdl = create_model()

#%%
@time_printer
def objective(trial):

    params = {
        'arr_e': trial.suggest_uniform('arr_e', 1000, 10000),
        'K_P': trial.suggest_uniform('K_P', 0.01, 100),
        'f_CH_max': trial.suggest_uniform('f_CH_max', 0.1, 10),
        'exponent': trial.suggest_uniform('exponent', 1, 10),
        'q_CH': trial.suggest_uniform('q_CH', 0.1, 10),
        'q_LI': trial.suggest_uniform('q_LI', 1.5, 50),
        'V_NH': trial.suggest_uniform('V_NH', 0.01, 1),
        'V_P': trial.suggest_uniform('V_P', 0.01, 1),
        }

    current_params = np.array([])

    for k, v in params.items():
        current_params = np.append(current_params, v)

    # mdl._update_state(current_params, t_span=(0, 0.25), t_eval = np.arange(0, 0.26, 0.01), method='BDF', state_reset_hook='reset_cache')

    try:
        mdl._update_state(current_params, t_span=(0, 50), t_eval = np.arange(0, 51, 1), method='RK23', state_reset_hook='reset_cache')
        # if time

    except:
        return 0.5

    out = [metric() for metric in mdl.metrics]
    obj = np.average(out)

    # for step in range(100):
    #     trial.report(obj, step)

    #     if trial.should_prune():
    #         raise optuna.TrialPruned()

    return obj

if __name__ == '__main__':

    sampler = optuna.samplers.TPESampler(seed=333)

    study = optuna.create_study(sampler=sampler, direction='minimize')
    # study = optuna.create_study(sampler=sampler, direction='minimize', pruner=optuna.pruners.MedianPruner())
    # w/ or w/o lead to the same results & the same time required, but hyperparameter importances and optimization history changes

    study.optimize(objective, n_trials=2)

    # study.optimize(objective, n_trials=100, n_jobs=-1)  # errors due to python's GIL?
    # parallelization -sql

    print(study.best_params)

    df = study.trials_dataframe()
    assert isinstance(df, pd.DataFrame)
    # assert df.shape[0] == 10000
    assert df.shape[0] == 2

    df.to_excel(ospath.join(results_path, 'conti_optuna.xlsx'))

    optuna.visualization.matplotlib.plot_param_importances(study)
    optuna.visualization.matplotlib.plot_optimization_history(study)
    # optuna.visualization.matplotlib.plot_intermediate_values(study)

    # optuna.visualization.plot_param_importances(study).show(renderer='browser')
    # optuna.visualization.plot_optimization_history(study).show(renderer='browser')

#%% print time & beep

# time_elapsed = datetime.now()-start_time
# print('Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

# sd.Beep(2000, 1000)