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
from exposan.pm2_batch import (
    results_path,
    create_model,
    sensitive_params,
    )

import numpy as np, pandas as pd
import optuna


from datetime import datetime
start_time = datetime.now()

__all__ = ('objective')

# start_time = datetime.now()
# print('Started (hh:mm:ss.ms) {}'.format(start_time))

#%%

mdl = create_model(kind='include', analysis='cali')   # with uniform distribution of sensitive params (in model.py)
# mdl = create_model(kind='exclude', analysis='cali')   # with uniform distribution of sensitive params (in model.py)

#%%
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
        mdl._update_state(current_params, t_span=(0, 7), t_eval = np.arange(0, 7.01, 0.01), method='BDF', state_reset_hook='reset_cache')

    except:
        return 0.5

    out = [metric() for metric in mdl.metrics]
    obj = np.average(out)

    for step in range(100):
        trial.report(obj, step)

        if trial.should_prune():
            raise optuna.TrialPruned()

    return obj

#%%
# if __name__ == '__main__':

#     sampler = optuna.samplers.TPESampler(seed=777) #unit: 1000:0.195 (3 min), 5000:0.1363 (30-35min) (best~2670), 10000:0.1239 (2hr) (best~8593)
#     sampler = optuna.samplers.TPESampler(seed=333) #unit: 1000:0.195 (3 min), 5000:0.1363 (30-35min) (best~2670), 10000:0.1239 (2hr) (best~8593)

#     # study = optuna.create_study(sampler=sampler, direction='minimize')
#     study = optuna.create_study(sampler=sampler, direction='minimize', pruner=optuna.pruners.MedianPruner())
#     # w/ or w/o lead to the same results & the same time required, but hyperparameter importances and optimization history changes

#     # study.optimize(objective, n_trials=10000)  # obj자체는 10000이 더 낫지만, 결과plot해보면 5000이 더 나음
#     study.optimize(objective, n_trials=5000)

#     # study.optimize(objective, n_trials=100, n_jobs=-1)  # errors due to python's GIL?
#     # parallelization -sql

#     print(study.best_params)

#     df = study.trials_dataframe()
#     assert isinstance(df, pd.DataFrame)
#     # assert df.shape[0] == 10000
#     assert df.shape[0] == 5000

#     df.to_excel(ospath.join(results_path, 'batch_exclude_kinetic_optuna.xlsx'))

#     optuna.visualization.matplotlib.plot_param_importances(study)
#     optuna.visualization.matplotlib.plot_optimization_history(study)
    # optuna.visualization.matplotlib.plot_intermediate_values(study)

    # optuna.visualization.plot_param_importances(study).show(renderer='browser')
    # optuna.visualization.plot_optimization_history(study).show(renderer='browser')

#%% parallelization

if __name__ == '__main__':

    sampler = optuna.samplers.TPESampler(seed=777) #unit: 1000:0.195 (3 min), 5000:0.1363 (30-35min) (best~2670), 10000:0.1239 (2hr) (best~8593)

    study = optuna.create_study(sampler=sampler, direction='minimize', pruner=optuna.pruners.MedianPruner(),
                                storage="mysql://root@localhost/example")

    study.optimize(objective, n_trials=1000)

    print(study.best_params)

    df = study.trials_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] == 1000

    df.to_excel(ospath.join(results_path, 'batch_include_kinetic_optuna_parallel.xlsx'))

    optuna.visualization.matplotlib.plot_param_importances(study)
    optuna.visualization.matplotlib.plot_optimization_history(study)



#%% print time & beep

# time_elapsed = datetime.now()-start_time
# print('Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))

# sd.Beep(2000, 1000)