# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan.utils import ospath
from exposan.pm2asm2d import (
    results_path,
    create_model,
    )

import numpy as np, pandas as pd
import optuna

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

__all__ = ('objective')

#%%
mdl = create_model(kind='', analysis='cali')
# mdl = create_model(kind='include', analysis='cali')
# mdl = create_model(kind='exclude', analysis='cali')

#%%
def objective(trial):

    params = {
        'a_c': trial.suggest_uniform('a_c', 0.005, 0.5),
        'arr_e': trial.suggest_uniform('arr_e', 4000, 8000),
        'K_P': trial.suggest_uniform('K_P', 0.1, 10),
        'K_A': trial.suggest_uniform('K_A', 1, 10),        
        'f_CH_max': trial.suggest_uniform('f_CH_max', 0.1, 10),
        'f_LI_max': trial.suggest_uniform('f_LI_max', 1, 10),
        'V_NH': trial.suggest_uniform('V_NH', 0.01, 1),
        'V_NO': trial.suggest_uniform('V_NO', 0.01, 1),        
        'V_P': trial.suggest_uniform('V_P', 0.001, 0.1),
        'eta_fe': trial.suggest_uniform('eta_fe', 0.1, 1),
        'K_O2': trial.suggest_uniform('K_O2', 0.1, 1),
        'K_O2_H': trial.suggest_uniform('K_O2_H', 0.1, 1),
        'K_P_H': trial.suggest_uniform('K_P_H', 0.001, 0.1),
        'mu_AUT': trial.suggest_uniform('mu_AUT', 0.1, 5),
        }  

    current_params = np.array([])

    for k, v in params.items():
        current_params = np.append(current_params, v)

    try:
        mdl._update_state(current_params, t_span=(0, 9), t_eval = np.arange(0, 9.1, 0.1), method='BDF', state_reset_hook='reset_cache')
        # mdl._update_state(current_params, t_span=(0, 7), t_eval = np.arange(0, 7.01, 0.01), method='BDF', state_reset_hook='reset_cache')
        
        out = [metric() for metric in mdl.metrics]
        obj = np.average(out)
        return obj
    
    except:
        return 15
    
#%% some trial can't be fully integrated until 9, instead it stops before 9
# in that case, system returns interpolation range error, so we have to return 15
    # if mdl._system.units[0].scope.time_series[-1] == 9:
    #     out = [metric() for metric in mdl.metrics]
    #     obj = np.average(out)
    #     return obj
    
    # else:
    #     return 15

#%%
if __name__ == '__main__':

    sampler = optuna.samplers.TPESampler(seed=333)

    study = optuna.create_study(sampler=sampler, direction='minimize')
    # study = optuna.create_study(sampler=sampler, direction='minimize', pruner=optuna.pruners.MedianPruner())

    study.optimize(objective, n_trials=10000)    # takes about 30 min

    print(study.best_params)

    df = study.trials_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] == 10000

    df.to_excel(ospath.join(results_path, 'calibration_result_optuna.xlsx'))
    # df.to_excel(ospath.join(results_path, 'calibration_result_include_optuna_sequential_cali.xlsx'))
    # df.to_excel(ospath.join(results_path, 'calibration_result_exclude_optuna_sequential_cali.xlsx'))

    # optuna.visualization.matplotlib.plot_param_importances(study)
    # optuna.visualization.matplotlib.plot_optimization_history(study)