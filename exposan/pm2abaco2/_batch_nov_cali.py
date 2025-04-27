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

cmps = pc.create_pm2abaco2_cmps()

pm2abaco2 = pc.PM2ABACO2(arr_e=6351.18337472904, K_P=8.02801655725152, K_A=9.01340931339147, 
                         rho=3.6371936397996, K_STO=7.54872301376646, m_ATP=8.45523356363203, 
                         mu_max=0.865772845547299, q_LI=11.9873181834793, 
                         V_NH=0.0100103044597615, V_NO=0.00449164584420843, V_P=0.00956331001640809, 
                         
                         Y_G=0.999904421654707, f_BAC=0.0523677736811379, mu_max_HET=1.7913342116121, 
                         temp_min_NIT=264.759553255619, temp_max_NIT=300.607160923099, temp_opt_NIT=286.477196748904, 
                         temp_max_HET=317.371181590779, temp_opt_HET=291.691302314683, 
                         ph_max_NIT=10.5491787462377, ph_min_HET=5.85138668399908, 
        
                         a_c=0.049, I_n=1500, arr_a=1.8e10, beta_1=2.90,
                         beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
                         K_N=0.1, K_G=6.3, f_CH_max=0.819, f_LI_max=3.249, 
                         q_CH=0.594, Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163, exponent=4, 
                         Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
                         Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
                         Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
                         Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
                         Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317, n_dark=0.7,  # modified baseline for pm2 (I_n, I_opt)
                         
                         Y_NH_NIT=0.18, Y_NO_NIT=0.19, Y_NH_HET=9.09,  Y_O2_NIT=0.09, Y_O2_HET=2.78, 
                         mu_max_NIT=0.75, 
                         temp_min_HET=270,  ph_min_NIT=2, ph_opt_NIT=9,
                         ph_max_HET=12, ph_opt_HET=9, K_S_O2_NIT=1.08, K_I_O2_NIT=104.9, K_S_O2_HET=1.98,
                         K_S_NH_NIT=1.0, K_S_NH_HET=0.5, K_S_G_HET=0.32, theta_NIT=274.1, theta_HET=274.07,
                         )

# pm2abaco2 = pc.PM2ABACO2(a_c=0.0140264685430579, arr_e=6101.68141831652, K_P=4.13107619555845, K_A=2.77768066194642,
#                        f_CH_max=2.91791702097432, f_LI_max=6.17346108004295,  
#                        V_NH=0.0274321398165965, V_NO=0.00399, V_P=0.003977073923784, 
#                        # V_NH=0.0274321398165965, V_NO=0.00399447778517413, V_P=0.003977073923784, 

#                        eta_fe=0.497814493730393, K_O2=0.95589810728547,  K_O2_H=0.129725944552192, 
#                        K_P_H=0.0999763470472844, mu_AUT=0.833426065767295,

#                        I_n=1500, arr_a=1.8e10, beta_1=2.90,
#                        beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
#                        K_N=0.1, K_G=6.3, rho=1.186, K_STO=1.566,
#                        m_ATP=15.835, mu_max=1.969, q_CH=0.594, q_LI=0.910,
#                        Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163, exponent=4,
#                        Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
#                        Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
#                        Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
#                        Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
#                        Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317, n_dark=0.7,  
                       
#                        f_SI=0.0, Y_H=0.625, f_XI_H=0.1, Y_A=0.24, f_XI_AUT=0.1,
#                        K_h=3.0, eta_NO3=0.6, K_NO3=0.5, K_X=0.1,
#                        mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4,K_F_H=4.0,
#                        K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_ALK_H=0.1,
#                        b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,  
#                        )   # seed333, n=5000, 'V_NO': trial.suggest_uniform('V_NO', 0.0007, 0.01)

# pm2asm2d = pc.PM2ASM2d(a_c=0.0133531740980039, arr_e=4527.1446341968, K_P=7.24385749531272, K_A=1.19865676173447,
#                        f_CH_max=0.440984103016091, f_LI_max=8.27782534104099,  
#                        V_NH=0.0269485896197009, V_NO=0.0100476771433462, V_P=0.00529581539056508, 
#                        eta_fe=0.509973141656043, K_O2=0.160459502203415,  K_O2_H=0.135500798373871, 
#                        K_P_H=0.00789252301210439, mu_AUT=3.32504961805362,

#                        I_n=1500, arr_a=1.8e10, beta_1=2.90,
#                        beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
#                        K_N=0.1, K_G=6.3, rho=1.186, K_STO=1.566,
#                        m_ATP=15.835, mu_max=1.969, q_CH=0.594, q_LI=0.910,
#                        Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163, exponent=4,
#                        Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
#                        Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
#                        Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
#                        Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
#                        Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317, n_dark=0.7,  
                       
#                        f_SI=0.0, Y_H=0.625, f_XI_H=0.1, Y_A=0.24, f_XI_AUT=0.1,
#                        K_h=3.0, eta_NO3=0.6, K_NO3=0.5, K_X=0.1,
#                        mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4,K_F_H=4.0,
#                        K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_ALK_H=0.1,
#                        b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,  
#                        )   # seed333, n=10000, 'V_NO': trial.suggest_uniform('V_NO', 0.01, 1)

pm2abaco2_path = ospath.dirname(__file__)
data_path = ospath.join(pm2abaco2_path, 'data/exo_vars_batch_november_ph.xlsx')

T, I, pH = EDV.batch_init(data_path, 'linear')

PBR = su.BatchExperiment('PBR', model=pm2abaco2, exogenous_vars=(T, I, pH))

init_conds = {
    'X_CHL':2.61,
    'X_ALG':521.53,
    'X_PG':12.83,
    'X_TAG':67.37,
    'S_CO2':30.0,
    'S_A':5.0,
    'S_G':5.0,
    'S_O2':20.36,
    'S_NH':32.10,
    'S_NO':91.17,
    'S_P':2.86,
    'X_N_ALG':0,
    'X_P_ALG':13.23,
    'X_NIT':10,             # copied from ASM2d X_AUT
    'X_HET':20,   
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
                      export_state_to=f'results/sol_{t}d_{method}_batch_nov_cali_optuna_calibrated.xlsx',
                      # export_state_to=f'results/sol_{t}d_{method}_batch_nov_cali_optuna_seed333.xlsx',
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