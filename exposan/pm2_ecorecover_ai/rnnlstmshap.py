# -*- coding: utf-8 -*-
'''
This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Jooman Kim <jooman@illinois.edu>

'''
#%% Import packages
from qsdsan.utils import ospath
from exposan.pm2_ecorecover_ai import data_path, results_path

import pandas as pd, numpy as np, math, random, time

from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf, keras_tuner as kt

import shap

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

#%% Data load
df = pd.read_csv(ospath.join(data_path, 'cleaned_data.csv'), sep='\t')

# Select y
y_col = 'nh4_eff'
# y_col = 'no3_eff'
# y_col = 'po4_eff'
# y_col = 'tss_mix'

X = df.set_index('t_stamp')
y = df.set_index('t_stamp')[[y_col]]

#%% Data split (6:2:2)
train_split = math.ceil(len(df) * .6)
val_split = math.ceil(len(df) * .8)

#%% Data scaling
x_scaler = MinMaxScaler()
y_scaler = MinMaxScaler()

X_mat = X.values
x_scaler.fit(X_mat[:train_split])
X_scaled = x_scaler.transform(X_mat)

y_mat = y.values
y_scaler.fit(y_mat[:train_split])
y_scaled = y_scaler.transform(y_mat)

#%% Util function to make window
def multivariate_data(dataset, target, start_index, end_index, history_size,
                      target_size, step, single_step=False):      
    data = []  # x
    labels = []  # y

    start_index = start_index + history_size
    if end_index is None:
        end_index = len(dataset) - target_size    

    for i in range(start_index, end_index):
        indices = range(i-history_size, i, step)
        data.append(dataset[indices])

        if single_step:
            labels.append(target[i+target_size])
        else:
            labels.append(target[i:i+target_size])

    return np.array(data), np.array(labels)

#%% Make window
future_target = 1       # future 5 minute
past_history = 1440     # past 5 days
step = 12               # data interval = 1 hour

x_train_multi, y_train_multi = multivariate_data(X_scaled, y_scaled, 0,                 
                                                 train_split, past_history,
                                                 future_target, step)
x_val_multi, y_val_multi = multivariate_data(X_scaled, y_scaled,                        
                                             train_split, val_split, past_history,      
                                             future_target, step)
x_test_multi, y_test_multi = multivariate_data(X_scaled, y_scaled, 
                                               val_split-past_history-1, None, past_history,
                                               future_target, step)

#%% LSTM & RNN Hyperparameter tuning
seed = 33  
tf.random.set_seed(seed)
np.random.seed(seed)
random.seed(seed)

# Choose model
#use_lstm = True        # True: lstm
use_lstm = False       # False: rnn

def lstm_model_builder(hp, x_train_multi, y_train_multi):   # LSTM
    n_outputs = y_train_multi.shape[1]

    model = tf.keras.models.Sequential()
    model.add(tf.keras.layers.LSTM(hp.Int('input_unit', min_value=32, max_value=512, step=32), 
                                   return_sequences=True, input_shape= x_train_multi.shape[-2:]))
    for i in range(hp.Int('n_layers', 1, 4)):
        model.add(tf.keras.layers.LSTM(hp.Int(f'lstm_{i}_units', min_value=32, max_value=512, step=32),
                                       return_sequences=True))
    model.add(tf.keras.layers.LSTM(hp.Int('layer_2_neurons', min_value=32, max_value=512, step=32)))
    model.add(tf.keras.layers.Dropout(hp.Float('Dropout_rate', min_value=0, max_value=0.5, step=0.1)))
    model.add(tf.keras.layers.Dense(30, activation=hp.Choice('dense_activation', 
                                                             values=['relu', 'sigmoid'], default='relu')))
    model.add(tf.keras.layers.Dense(n_outputs, activation=hp.Choice('dense_activation', 
                                                                    values=['relu', 'sigmoid'], default='relu')))
    model.compile(loss='mean_squared_error', optimizer='adam', metrics = ['mse'])

    return model

def rnn_model_builder(hp, x_train_multi, y_train_multi):     # RNN
    n_timesteps, n_features = x_train_multi.shape[-2:]
    n_outputs = y_train_multi.shape[1]

    model = tf.keras.models.Sequential()
    model.add(tf.keras.layers.SimpleRNN(hp.Int('input_unit', min_value=50, max_value=150, step=20), 
                                        dropout=hp.Float('in_dropout', min_value=0, max_value=0.5, step=0.1), 
                                        input_shape=(n_timesteps, n_features),                                         
                                        return_sequences=True, seed=seed)) 
    model.add(tf.keras.layers.SimpleRNN(hp.Int('layer 1', min_value=50, max_value=150, step=20), 
                                        activation=hp.Choice("l1_activation", values=["tanh", "relu", "sigmoid"]),
                                        dropout=hp.Float('l1_dropout', min_value=0, max_value=0.5, step=0.1), 
                                        return_sequences=True, seed=seed))
    model.add(tf.keras.layers.SimpleRNN(hp.Int('layer 2', min_value=50, max_value=150, step=20), 
                                        activation=hp.Choice("l2_activation", values=["tanh", "relu", "sigmoid"]),
                                        dropout=hp.Float('l2_dropout', min_value=0, max_value=0.5, step=0.1), 
                                        return_sequences=True, seed=seed))
    model.add(tf.keras.layers.SimpleRNN(hp.Int('layer 3', min_value=20, max_value=150, step=20), 
                                        activation=hp.Choice("l3_activation", values=["tanh", "relu", "sigmoid"]),
                                        dropout=hp.Float('l3_dropout', min_value=0, max_value=0.5, step=0.1), 
                                        return_sequences=True, seed=seed))
    model.add(tf.keras.layers.SimpleRNN(hp.Int('layer 4',min_value=20,max_value=150,step=20), 
                                        activation=hp.Choice("l4_activation", values=["tanh", "relu", "sigmoid"]),
                                        dropout=hp.Float('l4_dropout', min_value=0, max_value=0.5, step=0.1),
                                        seed=seed))
    model.add(tf.keras.layers.Dense(n_outputs))
    model.compile(loss='mean_squared_error', optimizer='adam', metrics = ['mse'])

    return model

#%% SHAP
def get_single_timestep(data, offset):
    return data[:, -offset, :]

def slice_data(train, val, test, offset):        
    sliced_train = get_single_timestep(train, offset)
    sliced_val = get_single_timestep(val, offset)
    sliced_test = get_single_timestep(test, offset)

    sliced_train = np.expand_dims(sliced_train, axis=1)
    sliced_val = np.expand_dims(sliced_val, axis=1)
    sliced_test = np.expand_dims(sliced_test, axis=1)

    return sliced_train, sliced_val, sliced_test

def finetune(x_train_multi, y_train_multi, x_val_multi, y_val_multi, use_lstm=True):    
    if use_lstm:
        model_builder = lambda hp: lstm_model_builder(hp, x_train_multi, y_train_multi)   
    else:
        model_builder = lambda hp: rnn_model_builder(hp, x_train_multi, y_train_multi)

    experiment_id = 'ex_' + time.strftime('%Y%m%d_%H%M%S', time.localtime(time.time())) 
    
    tuner = kt.RandomSearch(model_builder, objective="mse", max_trials=20, executions_per_trial=1, 
                            directory = "./{}".format(experiment_id), seed=seed)
    tuner.search(x=x_train_multi, y=y_train_multi, epochs=1, batch_size=256, 
                 validation_data=(x_val_multi, y_val_multi))
    
    multi_step_model = tuner.get_best_models(num_models=1)[0]
    multi_step_model.summary()

    return multi_step_model

def run_pipeline(X_train, y_train, X_val, y_val, X_test, y_test):
    num_time_steps = 120                                   
    num_samples = 50                                        
    num_features = X_train.shape[-1]                             

    shap_values_3d = np.zeros((num_time_steps, num_samples, num_features))      
    test_samples_3d = np.zeros((num_time_steps, num_samples, num_features))

    for offset in range(1, num_time_steps+1, 10):                 
        # Preprocess Data
        X_train_sl, X_val_sl, X_test_sl = slice_data(X_train, X_val, X_test, offset)         
        
        # Finetune
        model = finetune(X_train_sl, y_train, X_val_sl, y_val)                

        random_indices = np.random.choice(X_train_sl.shape[0], size=200, replace=False)            
        background = X_train_sl[random_indices]

        test_indices = np.random.choice(X_test_sl.shape[0], size=num_samples, replace=False)  
        test_samples = X_test_sl[test_indices]
       
        explainer = shap.GradientExplainer(model, background)  
        shap_values = explainer.shap_values(test_samples)       
        
        test_samples_2d = np.squeeze(test_samples, axis=1)  
        shap_values_2d = np.squeeze(shap_values, axis=(1, 3))  
       
        shap_values_3d[offset-1] = shap_values_2d   
        test_samples_3d[offset-1] = test_samples_2d       

    return shap_values_3d, test_samples_3d                    

shap_values_3d, test_samples_3d = run_pipeline(x_train_multi, y_train_multi, 
                                                x_val_multi, y_val_multi, 
                                                x_test_multi, y_test_multi)

#%% SHAP rank
shap_ranks_dict = {}

for t in range(0, shap_values_3d.shape[0], 10):
    mean_abs_shap = np.mean(np.abs(shap_values_3d[t]), axis=0)  
    sorted_indices = np.argsort(mean_abs_shap)[::-1]  
    sorted_features = np.array(df.columns.drop("t_stamp").tolist())[sorted_indices]  
    sorted_shap_values = mean_abs_shap[sorted_indices]  
    ranks = np.arange(1, len(sorted_features) + 1)  

    shap_ranks_dict[f"t-{t+1}"] = {"Feature": sorted_features,
                                   "Mean Absolute SHAP Value": sorted_shap_values,
                                   "Rank": ranks
                                   }

df_shap_ranks = pd.concat(
    [pd.DataFrame(shap_ranks_dict[t]).assign(Time_Step=t) for t in shap_ranks_dict],
    ignore_index=True
    )

#%% SHAP plot
heatmap_data_full = df_shap_ranks.pivot(index='Feature', columns='Time_Step', values='Rank')

cols = [c for c in heatmap_data_full.columns if c not in ['t-101', 't-111']]
cols += ['t-101', 't-111']
heatmap_data_full = heatmap_data_full[cols]

sorted_features = heatmap_data_full.sort_values(by='t-1')[['t-1']].index
heatmap_data_sorted = heatmap_data_full.loc[sorted_features]

custom_colors = ["#ff0051","#fefefe", "#0086f9"]
custom_cmap = LinearSegmentedColormap.from_list("pastel_custom", custom_colors)

plt.figure(figsize=(8, len(sorted_features) * 0.3))
ax = sns.heatmap(
        heatmap_data_sorted,
        annot=False,
        fmt=".0f",
        cmap=custom_cmap,
        cbar_kws={'label': 'Rank', 'shrink': 0.3},
        linewidths=0.5,
        linecolor='white',
        )
colorbar = ax.collections[0].colorbar
colorbar.ax.invert_yaxis()

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(0.5)

ax.tick_params(axis='both', which='both', direction='inout', length=6, width=1)
ax.set_xlabel("")
ax.set_ylabel("")
plt.tight_layout()

plt.savefig(ospath.join(results_path, 'lstm_shap.jpg') if use_lstm else 'rnn_shap.jpg', dpi=300, bbox_inches='tight')