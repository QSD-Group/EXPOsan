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

from qsdsan.utils import ospath
from exposan.pm2_batch import (
    results_path,
    create_model,
    )

import numpy as np, pandas as pd
import optuna

__all__ = ('objective')

#%%

mdl = create_model(kind='include', analysis='cali')
# mdl = create_model(kind='exclude', analysis='cali')

#%%
def objective(trial):

    params = {
        'q_CH': trial.suggest_uniform('q_CH', 0.1, 10),
        'q_LI': trial.suggest_uniform('q_LI', 1.5, 50),
        'V_NH': trial.suggest_uniform('V_NH', 0.01, 1),
        'V_P': trial.suggest_uniform('V_P', 0.01, 1),
        }  # sequential calibration

    current_params = np.array([])

    for k, v in params.items():
        current_params = np.append(current_params, v)

    try:
        mdl._update_state(current_params, t_span=(0, 0.25), t_eval = np.arange(0, 0.26, 0.01), method='BDF', state_reset_hook='reset_cache')
        # mdl._update_state(current_params, t_span=(0, 7), t_eval = np.arange(0, 7.01, 0.01), method='BDF', state_reset_hook='reset_cache')

    except:
        return 5

    out = [metric() for metric in mdl.metrics]
    obj = np.average(out)

    return obj

#%%
if __name__ == '__main__':

    sampler = optuna.samplers.TPESampler(seed=333)

    study = optuna.create_study(sampler=sampler, direction='minimize')
    # study = optuna.create_study(sampler=sampler, direction='minimize', pruner=optuna.pruners.MedianPruner())

    study.optimize(objective, n_trials=5000)    # takes about 30 min

    print(study.best_params)

    df = study.trials_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] == 5000

    df.to_excel(ospath.join(results_path, 'calibration_result_include_optuna_sequential_cali.xlsx'))
    # df.to_excel(ospath.join(results_path, 'calibration_result_exclude_optuna_sequential_cali.xlsx'))

    # optuna.visualization.matplotlib.plot_param_importances(study)
    # optuna.visualization.matplotlib.plot_optimization_history(study)