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

from qsdsan.utils import ospath
from exposan.pm2_batch import (
    results_path,
    create_model,
    sensitive_params,
    )

import numpy as np, pandas as pd
from scipy.optimize import minimize

__all__ = ('cali_setup', 'optimizer', 'objective_function')

#%% Questions
#1 When we create_model(kind='exclude', analysis='cali'), it will have sensitive parameters values as distribution we set,
#  and the other parameters will be based on 'basically 'default_pm2_kwargs'

                                                        # in system.py
                                                        # default_pm2_kwargs = dict(
                                                        #     a_c=0.049, I_n=1500, arr_a=1.8e10, arr_e=6842, beta_1=2.90,
                                                        #     beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
                                                        #     K_N=0.1, K_P=1.0, K_A=6.3, K_F=6.3, rho=1.186, K_STO=1.566,
                                                        #     f_CH_max=0.819, f_LI_max=3.249, m_ATP=15.835,
                                                        #     mu_max=1.969, q_CH=0.594, q_LI=0.910,
                                                        #     Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163,
                                                        #     V_NH=0.254, V_NO=0.254, V_P=0.016, exponent=4,
                                                        #     Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
                                                        #     Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
                                                        #     Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
                                                        #     Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
                                                        #     Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317, n_dark=0.7,
                                                        #     path=None,
                                                        #     )

# opt_params are automatically updated? starting from init_guess?

#%%

mdl = create_model(kind='exclude', analysis='cali')   # with uniform distribution of sensitive params (in model.py)

def cali_setup():
    params = []
    baseline = np.array([])
    boundary = []

    for k, v in sensitive_params.items():
        b, units, bounds = v
        params.append(k)
        baseline = np.append(baseline, b)
        boundary.append(bounds)
    boundary = tuple(boundary)

    return params, baseline, boundary

sense_params, opt_params, bnds = cali_setup()         # start from init_guess

'''
sense_params: list of sensitive parameters
opt_params: initial guess of sensitive parameters
bnds: min & max of sensitive parameters
'''

# need to let it know which parameters to change
# in objective-function?

#%%
def optimizer():

    opt = minimize(objective_function, opt_params, method='SLSQP', bounds=bnds, tol=1e-2)

    # scipy.optimize.minimize(fun, x0, args=(), method=None, jac=None, hess=None, hessp=None, bounds=None,\
    #                         constraints=(), tol=None, callback=None, options=None)

    # how to retrieve output of each run?

    opt_as_series = pd.Series(opt)
    opt_as_series.to_excel(excel_writer=(ospath.join(results_path, 'calibration_result.xlsx')))

def objective_function(opt_params, *args):

    mdl._update_state(opt_params, t_span=(0, 7), t_eval = np.arange(0, 7.01, 0.01), method='BDF', state_reset_hook='reset_cache')

    # how to know which parameters to be changed? need to point

    # mdl.evaluate()???
    # where to define t_step?
    # need to add into args=()?

    out = [metric() for metric in mdl.metrics]
    obj = np.average(out)

    return obj

optimizer()










#%% Brian's codes

# def fit_setup(i, eqn_set):
#     fit_params = np.array([])  # , dtype=np.float64)
#     bnds = []
#     sens_params = params.loc[eqn_list[i]].params
#     # Takes values and bounds of sensitive parameters from list of all uncertain parameters
#     for parameter in sens_params:
#         fit_params = np.append(fit_params, uparams_values.loc[parameter])
#         bnds.append(uparams_bnds.loc[parameter])
#     bnds = tuple(bnds)

#     # return fitfunc(fit_params, eqn_set, sens_params, )
#     return fitting(fit_params, bnds, eqn_set, sens_params)

# def fitting(fit_params, bnds, eqn_set, sens_params):
#     """
#     :param fit_params: initial guesses of parameters that will be fitted
#     :param bnds: bounds for fit parameters
#     :param eqn_set: tuple of equation decisions (e.g., Droop, Monod, or Andrews)
#     :param sens_params: parameters that were determined to be sensitive for a given model
#     :return: calculated model fit formatted as a pandas Series (includes objective function value, jacobian,
#               number of evaluations, number of iterations, optimized parameter values, and time it took to fit the model
#               uses "minimize" function from scipy with the SLSQP method
#     """
#     uparams_val = copy(uparams_values)
#     start1 = codetime()
#     fit = minimize(fitfunc, fit_params, args=(eqn_set, sens_params, uparams_val), method='SLSQP', bounds=bnds, tol=1e-2)
#     fit_as_series = pd.Series(fit)
#     fit_as_series['parameters'] = sens_params
#     fit_as_series['init_values'] = fit_params
#     fit_as_series.name = str(eqn_set)
#     # print((codetime()-start1)/60, eqn_set)
#     return fit_as_series


# def fitfunc(fit_params, *args):
#     """
#     :param fit_params: initial guesses of parameters that will be fitted
#     :param args: needed to include this to pass the eqn_set because of how scipy minimize function works
#     :return: value of objective function (RMSE)

#     -This function runs the model from 0 to the specified time (see ppm_params)
#     -After the model has been run, model outputs are compared to experimental data and the objective function
#     (i.e., root mean square error (RMSE)) is calculated
#     """
#     # eqn_set, sens_params = args
#     eqn_set, sens_params, uparams_val = args
#     # uparams_val = copy(uparams_values)

#     pm2 = model_setup(fit_params, eqn_set, sens_params, uparams_val)
#     # writer = pd.ExcelWriter('calibration_dmebait_model_output.xlsx', engine='xlsxwriter')
#     # pm2 = pd.DataFrame(pm2)
#     # pm2.to_excel(writer)
#     # writer.save()
#     pm2_frame = np.array(pm2.t)
#     pm2_frame = pm2_frame.reshape((len(pm2_frame), 1))
#     sv_output = np.array(pm2.y[[(num_pbrs + 1) * 20 + 6, (num_pbrs + 1) * 20 + 7, (num_pbrs + 1) * 20 + 9], :])
#     sv_output = sv_output.transpose()
#     pm2_frame = np.append(pm2_frame, sv_output, axis=1)
#     pm2_frame = pd.DataFrame(pm2_frame, columns=['time', 's_o2', 's_nh', 's_p'])

#     calibration_fitting_data = pd.read_csv('calibration_fitting_data_clean.csv')
#     calibration_fitting_data.time = round(calibration_fitting_data.time, 9)
#     # validation_fitting_data = pd.read_csv('validation_fitting_data_clean.csv')

#     # Match model values with calibration data based on time
#     calibration_fitting_data = calibration_fitting_data.merge(pm2_frame, on='time')
#     # print(calibration_fitting_data.shape[0])

#     # Calculate RMSE
#     calibration_fitting_data['s_nh_diff'] = (calibration_fitting_data['s_nh'] -
#                                               calibration_fitting_data['eff_amm']) ** 2
#     calibration_fitting_data['s_p_diff'] = (calibration_fitting_data['s_p'] -
#                                             calibration_fitting_data['eff_po4']) ** 2
#     calibration_fitting_data['s_o2_diff'] = (calibration_fitting_data['s_o2'] -
#                                               calibration_fitting_data['do_pbr']) ** 2
#     obj = np.sqrt(calibration_fitting_data['s_nh_diff'].sum() / calibration_fitting_data.shape[0]) + \
#           np.sqrt(calibration_fitting_data['s_p_diff'].sum() / calibration_fitting_data.shape[0]) + \
#           np.sqrt(calibration_fitting_data['s_o2_diff'].sum() / calibration_fitting_data.shape[0])
#     print(obj)
#     return obj


#%%





# # model = qs.Model(system=sys)

# mdl = create_model()

# def
#     new_params = np.array([])  # output from optimization function

#     mdl._update_state(new_params, t_span=(0, 7), method='BDF',...)

#     outputs = mdl.metrics[0]()