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

pm2 = pc.PM2(arr_e=3249, K_P=25.03, f_CH_max=7.527, exponent=7.752, q_CH=2.575, q_LI=13.64, V_NH=0.2577, V_P=0.7527,
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
              Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317) # from kinetic assay optimzed results



# pm2 = pc.PM2(arr_e=5500, K_P=0.01297, f_CH_max=0.1, exponent=9.952, q_CH=4.927, q_LI=1.617, V_NH=0.09811, V_P=0.4950,
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
#              Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317) # with shgo, iter 3, tot -3, ret 0.5

# pm2 = pc.PM2(arr_e=6842, K_P=0.8017, f_CH_max=0.4432, exponent=4.009, q_CH=2.021, q_LI=15.08, V_NH=0.1625, V_P=0.4) # with new baseline, new sens_params

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

# pm2 = pc.PM2(arr_e=6843, K_P=0.8056, K_STO=3.797, exponent=3.806, m_ATP=15.81, q_CH=2.594, V_NH=0.3537, V_P=0.3631) #minimize with original init
# pm2 = pc.PM2(mu_max=2, I_n=1000, I_opt=2000, V_NH=0.12, V_NO=0.0035, V_P=0.24, q_CH=1.5, q_LI=24, arr_a=4*10**10, f_CH_max=1.2, f_LI_max=3.9) original

pm2_path = ospath.dirname(__file__)
data_path = ospath.join(pm2_path, 'data/exo_vars_batch_may_unit.xlsx')

T, I = EDV.batch_init(data_path, 'linear')

PBR = su.BatchExperiment('PBR', model=pm2, exogenous_vars=(T, I))

init_concs = {
    'X_CHL':2.81,
    'X_ALG':561.57,
    'X_CH':13.74,
    'X_LI':62.22,
    'S_CO2':30.0,
    'S_A':5.0,
    'S_F':5.0,
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
                      export_state_to=f'results/sol_{t}d_{method}_batch_may_unit_cali_shgo_iter3_tol-3_ret0.5_new_kinetic_used.xlsx',
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


