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
from exposan.pm2abaco2 import (
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
        'arr_e': trial.suggest_uniform('arr_e', 4000, 8000),
        'K_P': trial.suggest_uniform('K_P', 0.1, 10),
        'K_A': trial.suggest_uniform('K_A', 1, 10),      
        'rho': trial.suggest_uniform('rho', 1, 5),      
        'K_STO': trial.suggest_uniform('K_STO', 0.1, 10),      
        'm_ATP': trial.suggest_uniform('m_ATP', 1, 20),      
        'mu_max': trial.suggest_uniform('mu_max', 0.1, 5),      
        'q_LI': trial.suggest_uniform('q_LI', 0.1, 20),      
        'V_NH': trial.suggest_uniform('V_NH', 0.01, 1),
        'V_NO': trial.suggest_uniform('V_NO', 0.0007, 0.01),        
        'V_P': trial.suggest_uniform('V_P', 0.001, 0.1),
        'Y_G': trial.suggest_uniform('Y_G', 0.1, 1),
        'f_BAC': trial.suggest_uniform('f_BAC', 0.01, 0.5),
        'mu_max_HET': trial.suggest_uniform('mu_max_HET', 1, 5),
        'temp_min_NIT': trial.suggest_uniform('temp_min_NIT', 260, 273),
        'temp_max_NIT': trial.suggest_uniform('temp_max_NIT', 300, 320),
        'temp_opt_NIT': trial.suggest_uniform('temp_opt_NIT', 273, 310),
        'temp_max_HET': trial.suggest_uniform('temp_max_HET', 300, 320),
        'temp_opt_HET': trial.suggest_uniform('temp_opt_HET', 273, 310),
        'ph_max_NIT': trial.suggest_uniform('ph_max_NIT', 10, 14),
        'ph_min_HET': trial.suggest_uniform('ph_min_HET', 4, 7),
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

    study.optimize(objective, n_trials=3000)    
    # study.optimize(objective, n_trials=5000)    

    print(study.best_params)

    df = study.trials_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] == 3000

    df.to_excel(ospath.join(results_path, 'calibration_result_optuna_pm2abaco2.xlsx'))
    # df.to_excel(ospath.join(results_path, 'calibration_result_include_optuna_sequential_cali.xlsx'))
    # df.to_excel(ospath.join(results_path, 'calibration_result_exclude_optuna_sequential_cali.xlsx'))

    # optuna.visualization.matplotlib.plot_param_importances(study)
    # optuna.visualization.matplotlib.plot_optimization_history(study)