# -*- coding: utf-8 -*-
'''
This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Jooman Kim <jooman@illinois.edu>

'''
#%% Import packages
from qsdsan.utils import ospath
from exposan.pm2_ecorecover_ai import data_path, results_path

import pandas as pd, math
import matplotlib.pyplot as plt

from lightgbm import LGBMRegressor
from sklearn.metrics import mean_squared_error

import optuna
from optuna.samplers import TPESampler

import shap

#%% Data load
df = pd.read_csv(ospath.join(data_path, 'cleaned_data.csv'), sep='\t')

# Select y
y_col = 'nh4_eff'
# y_col = 'no3_eff'
# y_col = 'po4_eff'
# y_col = 'tss_mix'

X = df.drop(columns=[y_col]).set_index('t_stamp')
y = df.set_index('t_stamp')[[y_col]]

#%% Data split (6:2:2)
train_split = math.ceil(len(df) * .6)
val_split = math.ceil(len(df) * .8)

X_train = X.iloc[:train_split, :]
y_train = y.iloc[:train_split, :]
X_val = X.iloc[train_split:val_split, :]
y_val = y.iloc[train_split:val_split, :]
X_test = X.iloc[val_split:, :]
y_test = y.iloc[val_split:, :]

#%% LGBM Hyperparameter tuning
seed=10                 
sampler = TPESampler(seed=seed)

# Objective function
def objective(trial):

    # lgbm_param = {
    #     'objective': 'regression',
    #     'verbose': -1,
    #     'metric': 'mse',
    #     'num_leaves': trial.suggest_int('num_leaves', 2, 1024, step=1, log=True),
    #     'colsample_bytree': trial.suggest_uniform('colsample_bytree', 0.7, 1.0),
    #     'reg_alpha': trial.suggest_uniform('reg_alpha', 0.0, 1.0),
    #     'reg_lambda': trial.suggest_uniform('reg_lambda', 0.0, 10.0),
    #     'max_depth': trial.suggest_int('max_depth', 3, 15),
    #     'learning_rate': trial.suggest_loguniform("learning_rate", 1e-8, 1e-2),
    #     'n_estimators': trial.suggest_int('n_estimators', 100, 3000),
    #     'min_child_samples': trial.suggest_int('min_child_samples', 5, 100),
    #     'subsample': trial.suggest_loguniform('subsample', 0.4, 1),
    #     }     # no3_eff, po4_eff, tss_mix

    lgbm_param = {
        'objective': 'regression',
        'verbose': -1,
        'metric': 'mse',
        'num_leaves': trial.suggest_int('num_leaves', 16, 64, step=1, log=True),
        'colsample_bytree': trial.suggest_uniform('colsample_bytree', 0.7, 1.0),
        'reg_alpha': trial.suggest_uniform('reg_alpha', 0.0, 0.5),
        'reg_lambda': trial.suggest_uniform('reg_lambda', 0.0, 0.5),
        'max_depth': -1,
        'learning_rate': trial.suggest_loguniform("learning_rate", 1e-3, 1e-1),
        'n_estimators': trial.suggest_int('n_estimators', 100, 3000),
        'min_child_samples': trial.suggest_int('min_child_samples', 20, 100),
        'subsample': trial.suggest_loguniform('subsample', 0.4, 1),
        }     # nh4_eff                                                                        

    model_lgbm = LGBMRegressor(**lgbm_param, early_stopping_rounds=25)               
    model_lgbm = model_lgbm.fit(X_train, y_train, eval_set=[(X_val, y_val)])

    MSE = mean_squared_error(y_val, model_lgbm.predict(X_val))                      
    return MSE                                                                       

# Optimization
optuna_lgbm = optuna.create_study(direction='minimize', sampler=sampler)  
optuna_lgbm.optimize(objective, n_trials=50)                 

# Optimized parameters
lgbm_trial = optuna_lgbm.best_trial                        
lgbm_trial_params = lgbm_trial.params                         

#%% Model train & test
lgbm = LGBMRegressor(**lgbm_trial_params) 
         
lgbm_study = lgbm.fit(X_train, y_train)            
pred = lgbm_study.predict(X_test)           

#%% SHAP
explainer = shap.TreeExplainer(lgbm_study)
shap_values = explainer.shap_values(X_test)

fig = shap.summary_plot(shap_values, X_test, max_display=10, show=False)
_, h = plt.gcf().get_size_inches()
plt.gcf().set_size_inches(h*5/3, h)    
plt.savefig(ospath.join(results_path, 'lgbm_shap.jpg'), bbox_inches='tight', dpi=300)