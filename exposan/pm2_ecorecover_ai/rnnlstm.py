# -*- coding: utf-8 -*-
'''
This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Jooman Kim <jooman@illinois.edu>

'''
#%% Import packages
from qsdsan.utils import ospath
from exposan.pm2_ecorecover_ai import data_path

import pandas as pd, numpy as np, math, random, time

from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf, keras_tuner as kt

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

#%% Convert data format
batch_size = 256           
buffer_size = 10000

seed = 42  
tf.random.set_seed(seed)
np.random.seed(seed)
random.seed(seed)

train_data_multi = tf.data.Dataset.from_tensor_slices((x_train_multi, y_train_multi))      
train_data_multi = train_data_multi.cache().shuffle(buffer_size, seed=seed,
                                                    reshuffle_each_iteration=False).batch(batch_size).repeat()             

val_data_multi = tf.data.Dataset.from_tensor_slices((x_val_multi, y_val_multi))
val_data_multi = val_data_multi.batch(batch_size).repeat()              

test_data_multi = tf.data.Dataset.from_tensor_slices((x_test_multi, y_test_multi))
test_data_multi = test_data_multi.batch(batch_size)       

#%% LSTM & RNN Hyperparameter tuning

# Choose model
#use_lstm = True        # True: lstm
use_lstm = False       # False: rnn

def lstm_model_builder(hp):   # LSTM
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

def rnn_model_builder(hp):     # RNN
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

model_builder = lstm_model_builder if use_lstm else rnn_model_builder   

experiment_id = 'ex_' + time.strftime('%Y%m%d_%H%M%S', time.localtime(time.time())) 

tuner = kt.RandomSearch(model_builder, objective="mse", max_trials=20, executions_per_trial=1, 
                        directory = "./{}".format(experiment_id), seed=seed)
tuner.search(x=x_train_multi, y=y_train_multi, epochs=1, batch_size=batch_size, 
             validation_data=(x_val_multi, y_val_multi))  

multi_step_model = tuner.get_best_models(num_models=1)[0]
multi_step_model.summary()

#%% Model test 
pred = multi_step_model.predict(test_data_multi)
pred_org = y_scaler.inverse_transform(pred)