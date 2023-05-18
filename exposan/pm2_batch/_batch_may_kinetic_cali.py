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

# pm2 = pc.PM2(arr_e=3249, K_P=25.03, f_CH_max=7.527, exponent=7.752, q_CH=2.575, q_LI=13.64, V_NH=0.2577, V_P=0.7527,
#              a_c=0.049, I_n=1500, arr_a=1.8e10, beta_1=2.90,
#              beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
#              K_N=0.1, K_A=6.3, K_F=6.3, rho=1.186, K_STO=1.566,
#              f_LI_max=3.249, m_ATP=10,
#              mu_max=1.969, Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163,
#              V_NO=0.003, n_dark=0.7,
#              Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
#              Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
#              Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
#              Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
#              Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317) # with shgo, iter 3, tot -3, ret 0.5, took 4 hrs

pm2 = pc.PM2(arr_e=5500, K_P=0.01297, f_CH_max=0.1, exponent=9.952, q_CH=4.927, q_LI=1.617, V_NH=0.09811, V_P=0.4950,
             a_c=0.049, I_n=1500, arr_a=1.8e10, beta_1=2.90,
             beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
             K_N=0.1, K_A=6.3, K_F=6.3, rho=1.186, K_STO=1.566,
             f_LI_max=3.249, m_ATP=10,
             mu_max=1.969, Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163,
             V_NO=0.003, n_dark=0.7,
             Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
             Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
             Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
             Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
             Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317)   #optimized from unit exp. using shgo

# modified_pm2_kwargs = dict(
#     a_c=0.049, I_n=1500, arr_a=1.8e10, arr_e=6842, beta_1=2.90,
#     beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
#     K_N=0.1, K_P=1.0, K_A=6.3, K_F=6.3, rho=1.186, K_STO=1.566,
#     f_CH_max=0.819, f_LI_max=3.249, m_ATP=10,
#     mu_max=1.969, q_CH=1, q_LI=15,
#     Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163,
#     V_NH=0.1, V_NO=0.003, V_P=0.2, exponent=4,
#     Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
#     Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
#     Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
#     Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
#     Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317, n_dark=0.7,
#     path=None,
#     )

# sensitive_params = {
#     'arr_e': (6842, 'K', (1000, 10000)),
#     'K_P': (1.0, 'g P/m^3', (0.01, 100)),
#     'f_CH_max': (0.819, 'g COD/g COD', (0.1, 10)),
#     'exponent': (4, '', (1, 10)),
#     'q_CH': (1, 'g COD/g COD/d', (0.1, 10)),
#     'q_LI': (15, 'g COD/g COD/d', (1.5, 50)),
#     'V_NH': (0.1, 'g N/g COD/d', (0.01, 1)),
#     'V_P': (0.2, 'g P/g COD/d', (0.01, 1))
#     }

######################## Previous optimization results #########################################
# pm2 = pc.PM2(arr_e=6832, K_P=1.363, f_CH_max=0.8239, exponent=3.980, q_CH=1.032, q_LI=15.02, V_NH=0.1022, V_P=0.2005)  # new base, new sens_params
# pm2 = pc.PM2(arr_e=4973, K_P=1.071, K_STO=1.297, exponent=4.002, m_ATP=1.489, q_CH=1.545, V_NH=0.2065, V_P=0.3299)  # 4rd, initial baseline as optimized value
# pm2 = pc.PM2(arr_e=4996, K_P=1.008, K_STO=1.018, exponent=4.003, m_ATP=1.048, q_CH=1.008, V_NH=0.09981, V_P=0.1047)  # 3rd_e-4
# pm2 = pc.PM2(arr_e=4996, K_P=1.009, K_STO=1.019, exponent=4.006, m_ATP=1, q_CH=1.009, V_NH=0.0999, V_P=0.1048)    # 2nd
# pm2 = pc.PM2(arr_e=6845, K_P=1.009, K_STO=1.584, exponent=4.006, m_ATP=15.87, q_CH=0.6032, V_NH=0.2586, V_P=0.001598)    _init
# pm2 = pc.PM2(mu_max=2, I_n=1000, I_opt=2000, V_NH=0.12, V_NO=0.0035, V_P=0.24, q_CH=1.5, q_LI=24, arr_a=4*10**10, f_CH_max=1.2, f_LI_max=3.9)

pm2_path = ospath.dirname(__file__)
data_path = ospath.join(pm2_path, 'data/exo_vars_batch_may_kinetic.xlsx')

T, I = EDV.batch_init(data_path, 'linear')

PBR = su.BatchExperiment('PBR', model=pm2, exogenous_vars=(T, I))

init_concs = {
    'X_CHL':3.91,
    'X_ALG':782.30,
    'X_CH':19.24,
    'X_LI':101.06,
    'S_CO2':30.0,
    'S_A':5.0,
    'S_F':5.0,
    'S_O2':20.36,
    'S_NH':27.84,
    'S_NO':10.72,
    'S_P':0.617,
    'X_N_ALG':0.92,
    'X_P_ALG':8.06,
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
                      export_state_to=f'results/sol_{t}d_{method}_batch_may_kinetic_cali_shgo_iter3_ftot-3_ret0.5_new_usingunit.xlsx',
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
    t = 7
    t_step = 0.1
    # method = 'RK45'
    # method = 'RK23'   # original
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