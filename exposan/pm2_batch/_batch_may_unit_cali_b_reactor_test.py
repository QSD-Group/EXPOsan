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

import numpy as np

import qsdsan.processes as pc, qsdsan.sanunits as su
from qsdsan import System
from qsdsan.utils import ospath, time_printer, \
    ExogenousDynamicVariable as EDV

cmps = pc.create_pm2_cmps()

pm2 = pc.PM2(arr_e=6842, K_P=1.0, f_CH_max=0.819, exponent=4, q_CH=1.92792246509906, q_LI=26.1535941900048, V_NH=0.150722549179019, V_P=0.540050768528713,
              a_c=0.049, I_n=1500, arr_a=1.8e10, beta_1=2.90,
              beta_2=3.50, b_reactor=0.003, I_opt=2000, k_gamma=1e-5,
              K_N=0.1, K_A=6.3, K_G=6.3, rho=1.186, K_STO=1.566,
              f_LI_max=3.249, m_ATP=10,
              mu_max=1.969, Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163,
              V_NO=0.003, n_dark=0.7,
              Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
              Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
              Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
              Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
              Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317)   # sequential calibration, seed333, b_reactor=0.003

# pm2 = pc.PM2(arr_e=6842, K_P=1.0, f_CH_max=0.819, exponent=4, q_CH=1.92792246509906, q_LI=26.1535941900048, V_NH=0.150722549179019, V_P=0.540050768528713,
#               a_c=0.049, I_n=1500, arr_a=1.8e10, beta_1=2.90,
#               beta_2=3.50, b_reactor=0.3, I_opt=2000, k_gamma=1e-5,
#               K_N=0.1, K_A=6.3, K_G=6.3, rho=1.186, K_STO=1.566,
#               f_LI_max=3.249, m_ATP=10,
#               mu_max=1.969, Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163,
#               V_NO=0.003, n_dark=0.7,
#               Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
#               Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
#               Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
#               Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
#               Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317)   # sequential calibration, seed333, b_reactor=0.3

pm2_path = ospath.dirname(__file__)
data_path = ospath.join(pm2_path, 'data/exo_vars_batch_may_unit.xlsx')

T, I = EDV.batch_init(data_path, 'linear')

PBR = su.BatchExperiment('PBR', model=pm2, exogenous_vars=(T, I))

init_concs = {
    'X_CHL':2.81,
    'X_ALG':561.57,
    'X_PG':13.74,
    'X_TAG':62.22,
    'S_CO2':30.0,
    'S_A':5.0,
    'S_G':5.0,
    'S_O2':20.36,
    'S_NH':25,
    'S_NO':9.30,
    'S_P':0.383,
    'X_N_ALG':3.62,
    'X_P_ALG':12.60,
    }

PBR.set_init_conc(**init_concs)

sys = System('sys', path=(PBR,))
sys.set_dynamic_tracker(PBR)

@time_printer
def run(t, t_step, method=None, print_t=False, **kwargs):
    if method:
        sys.simulate(state_reset_hook='reset_cache',
                      t_span=(0,t),
                      t_eval=np.arange(0, t+t_step, t_step),
                      method=method,
                      # rtol=1e-2,
                      # atol=1e-3,
                      export_state_to=f'results/sol_{t}d_{method}_batch_may_unit_cali_optuna_sequential_cali_breactor_0.003.xlsx',
                      # export_state_to=f'results/sol_{t}d_{method}_batch_may_unit_cali_optuna_sequential_cali_breactor_0.3.xlsx',
                      print_t=print_t,
                      **kwargs)
    else:
        sys.simulate(state_reset_hook='reset_cache',
                      solver='odeint',
                      t=np.arange(0, t+t_step/30, t_step/30),
                      # export_state_to=f'results/sol_{t}d_odeint.xlsx',
                      print_msg=True,
                      print_t=print_t,
                      **kwargs)

if __name__ == '__main__':
    t = 0.25
    t_step = 0.01
    # method = 'RK45'
    # method = 'RK23'  # original
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    # method = None
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    run(t, t_step, method=method,
        print_t = True,
        )