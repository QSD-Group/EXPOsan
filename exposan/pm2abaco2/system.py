# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np, qsdsan as qs
from qsdsan import processes as pc, sanunits as su, System
from qsdsan.utils import time_printer, ExogenousDynamicVariable as EDV
from exposan.pm2abaco2 import data_path

__all__ = (
    'create_system',
    'default_pm2abaco2_kwargs',
    'default_init_conds_nov',
    'T_nov', 'I_nov', 'pH_nov'
    # 'default_init_conds_inc',
    # 'default_init_conds_exc',
    # 'T_exc', 'I_exc',
    # 'T_inc', 'I_inc',
    )

#%%

# =============================================================================
# Parameters and initial conditions
# =============================================================================

T_nov, I_nov, pH_nov = EDV.batch_init(os.path.join(data_path, 'exo_vars_batch_november_ph.xlsx'), 'linear')

default_init_conds_nov = {
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
    'X_HET':20,             # copied from ASM2d X_H    
    # 'S_N2': 0,           # Concentration of dinitrogen
    # 'S_ALK': 100,        # Concentration of alkalinity
    # 'S_I': 5,            # Concentration of inert soluble organic material
    # 'X_I': 5,            # Concentration of inert particulate organic material
    # 'X_S': 5,            # Concentration of slowly biodegradable substrates
    # 'X_H': 20,           # Concentration of heterotrophic organisms (including denitrifer)
    # 'X_AUT': 10,         # Concentration of nitrifying organisms
    } 

# T_exc, I_exc = EDV.batch_init(os.path.join(data_path, 'exo_vars_batch_may_kinetic.xlsx'), 'linear')
# T_inc, I_inc = EDV.batch_init(os.path.join(data_path, 'exo_vars_batch_may_unit.xlsx'), 'linear')

# # default for 'include' - unit process
# default_init_conds_inc = {
#     'X_CHL':2.81,
#     'X_ALG':561.57,
#     'X_PG':13.74,
#     'X_TAG':62.22,
#     'S_CO2':30.0,
#     'S_A':5.0,
#     'S_G':5.0,
#     'S_O2':20.36,
#     'S_NH':25,
#     'S_NO':9.30,
#     'S_P':0.383,
#     'X_N_ALG':3.62,
#     'X_P_ALG':12.60,
#     }

# # default for 'exclude' - kinetic assay
# default_init_conds_exc = {
#     'X_CHL':3.91,
#     'X_ALG':782.30,
#     'X_PG':19.24,
#     'X_TAG':101.06,
#     'S_CO2':30.0,
#     'S_A':5.0,
#     'S_G':5.0,
#     'S_O2':20.36,
#     'S_NH':27.84,
#     'S_NO':10.72,
#     'S_P':0.617,
#     'X_N_ALG':0.92,
#     'X_P_ALG':8.06,
#     }

default_pm2abaco2_kwargs = dict(
    a_c=0.049, I_n=1500, arr_a=1.8e10, arr_e=6842, beta_1=2.90,
    beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
    K_N=0.1, K_P=1.0, K_A=6.3, K_G=6.3, rho=1.186, K_STO=1.566,
    f_CH_max=0.819, f_LI_max=3.249, m_ATP=15.835,
    mu_max=1.969, q_CH=0.594, q_LI=0.910,
    Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163,
    V_NH=0.254, V_NO=0.254, V_P=0.016, exponent=4,
    Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
    Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
    Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
    Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
    Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317, n_dark=0.7,  # modified baseline for pm2 (I_n, I_opt)
    
    Y_NH_NIT=0.18, Y_NO_NIT=0.19, Y_NH_HET=9.09, Y_G=0.45, Y_O2_NIT=0.09, Y_O2_HET=2.78, f_BAC=0.1,
    mu_max_NIT=0.75, mu_max_HET=3.4, temp_min_NIT=265, temp_max_NIT=311, temp_opt_NIT=293,
    temp_min_HET=270, temp_max_HET=315, temp_opt_HET=298, ph_min_NIT=2, ph_max_NIT=13.4, ph_opt_NIT=9,
    ph_min_HET=6, ph_max_HET=12, ph_opt_HET=9, K_S_O2_NIT=1.08, K_I_O2_NIT=104.9, K_S_O2_HET=1.98,
    K_S_NH_NIT=1.0, K_S_NH_HET=0.5, K_S_G_HET=0.32, theta_NIT=274.1, theta_HET=274.07,  # original baseline for abaco2
    
    # f_SI=0.0, Y_H=0.625, f_XI_H=0.1, Y_A=0.24, f_XI_AUT=0.1,
    # K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
    # mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F_H=4.0,
    # K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
    # mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,  # original baseline for asm2d
    path=None,
    )                                                  

# %%

# =============================================================================
# Batch system
# =============================================================================

def create_system(flowsheet=None, pm2abaco2_kwargs={}, init_conds={}, kind=''):
    flowsheet = flowsheet or qs.Flowsheet('pm2abaco2')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components
    pc.create_pm2abaco2_cmps()

    # Process models
    pm2abaco2_kwargs = pm2abaco2_kwargs or default_pm2abaco2_kwargs
    pm2abaco2 = pc.PM2ABACO2(**pm2abaco2_kwargs)

    # Create unit operations

    # if kind == 'include':
    #     PBR = su.BatchExperiment('PBR', model=pm2asm2d, exogenous_vars=(T_inc, I_inc))
    #     init_conds = init_conds or default_init_conds_inc
    
    # else:
    #     PBR = su.BatchExperiment('PBR', model=pm2asm2d, exogenous_vars=(T_exc, I_exc))
    #     init_conds = init_conds or default_init_conds_exc

    PBR = su.BatchExperiment('PBR', model=pm2abaco2, exogenous_vars=(T_nov, I_nov, pH_nov))
    init_conds = init_conds or default_init_conds_nov

    PBR.set_init_conc(**init_conds)

    # System setup
    sys = System('sys', path=(PBR,))

    sys.set_dynamic_tracker(PBR)
    sys.set_tolerance(rmol=1e-6)

    return sys

# %%

@time_printer
def run(t, t_step, method=None, **kwargs):
    sys = create_system()
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        # rtol=1e-2,
        # atol=1e-3,
        export_state_to=f'results/pm2abaco2_sol_{t}d_{method}.xlsx',
        **kwargs)

if __name__ == '__main__':
    t = 9
    t_step = 1
    # method = 'RK45'
    method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    # method = 'LSODA'
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    run(t, t_step, method=method)