# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np

import qsdsan.processes as pc, qsdsan.sanunits as su
from qsdsan import System
from qsdsan.utils import ospath, time_printer, \
    ExogenousDynamicVariable as EDV

cmps = pc.create_pm2asm2d_cmps()

pm2asm2d = pc.PM2ASM2d(a_c=0.0133531740980039, arr_e=4527.1446341968, K_P=7.24385749531272, K_A=1.19865676173447,
                       f_CH_max=0.440984103016091, f_LI_max=8.27782534104099,  
                       # V_NH=0.0269485896197009, V_NO=0.007, V_P=0.00529581539056508,                                   # low V_NO
                       V_NH=0.0269485896197009, V_NO=0.0100476771433462, V_P=0.00529581539056508, 

                       eta_fe=0.509973141656043, K_O2=0.160459502203415,  K_O2_H=0.135500798373871, 
                       K_P_H=0.00789252301210439, mu_AUT=3.32504961805362,

                       I_n=1500, arr_a=1.8e10, beta_1=2.90,
                       beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
                       K_N=0.1, K_F=6.3, rho=1.186, K_STO=1.566,
                       m_ATP=15.835, mu_max=1.969, q_CH=0.594, q_LI=0.910,
                       Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163, exponent=4,
                       Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
                       Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
                       Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
                       Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
                       Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317, n_dark=0.7,  
                       
                       f_SI=0.0, Y_H=0.625, f_XI_H=0.1, Y_A=0.24, f_XI_AUT=0.1,
                       K_h=3.0, eta_NO3=0.6, K_NO3=0.5, K_X=0.1,
                       mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4,K_F_H=4.0,
                       K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_ALK_H=0.1,
                       b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,  
                       )   # seed333, 

pm2asm2d_path = ospath.dirname(__file__)
data_path = ospath.join(pm2asm2d_path, 'data/exo_vars_batch_november.xlsx')

T, I = EDV.batch_init(data_path, 'linear')

PBR = su.BatchExperiment('PBR', model=pm2asm2d, exogenous_vars=(T, I))

init_conds = {
    'X_CHL':2.61,
    'X_ALG':521.53,
    'X_CH':12.83,
    'X_LI':67.37,
    'S_CO2':30.0,
    'S_A':5.0,
    'S_F':5.0,
    'S_O2':20.36,
    'S_NH':32.10,
    'S_NO':91.17,
    'S_P':2.86,
    'X_N_ALG':0,
    'X_P_ALG':13.23,
    'S_N2': 0,           # Concentration of dinitrogen
    'S_ALK': 100,        # Concentration of alkalinity
    'S_I': 5,            # Concentration of inert soluble organic material
    'X_I': 5,            # Concentration of inert particulate organic material
    'X_S': 5,            # Concentration of slowly biodegradable substrates
    'X_H': 20,           # Concentration of heterotrophic organisms (including denitrifer)
    'X_AUT': 10,  
    }

PBR.set_init_conc(**init_conds)

sys = System('sys', path=(PBR,))
sys.set_dynamic_tracker(PBR)
sys.set_tolerance(rmol=1e-6)

@time_printer
def run(t, t_step, method=None, print_t=True, **kwargs):
    if method:
        sys.simulate(state_reset_hook='reset_cache',
                      t_span=(0,t),
                      t_eval=np.arange(0, t+t_step, t_step),
                      method=method,
                      # rtol=1e-2,
                      # atol=1e-3,
                      export_state_to=f'results/sol_{t}d_{method}_batch_nov_cali_optuna_seed333.xlsx',
                      # export_state_to=f'results/sol_{t}d_{method}_batch_nov_cali_optuna_vno_0.007.xlsx',
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
    t = 9
    t_step = 0.1
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