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
from exposan.pm2asm2d import data_path

__all__ = (
    'create_system',
    'default_pm2asm2d_kwargs',
    'default_init_conds_inc',
    'default_init_conds_exc',
    'T_exc', 'I_exc',
    'T_inc', 'I_inc',
    )

#%%

# =============================================================================
# Parameters and initial conditions
# =============================================================================

T_exc, I_exc = EDV.batch_init(os.path.join(data_path, 'exo_vars_batch_may_kinetic.xlsx'), 'linear')
T_inc, I_inc = EDV.batch_init(os.path.join(data_path, 'exo_vars_batch_may_unit.xlsx'), 'linear')

# default for 'include' - unit process exp.
default_init_conds_inc = {
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
    'S_N2': 0,           # Concentration of dinitrogen
    'S_ALK':,          # Concentration of alkalinity
    'S_I':,            # Concentration of inert soluble organic material
    'X_I':,            # Concentration of inert particulate organic material
    'X_S':,            # Concentration of slowly biodegradable substrates
    'X_H':,            # Concentration of heterotrophic organisms (including denitrifer)
    'X_AUT':,          # Concentration of nitrifying organisms
    } 


# default_init_conds['asm2d'] = {
#     'S_F':5,
#     'S_A':2,
#     'X_I':1000,
#     'X_S':100,
#     'X_H':500,
#     'X_AUT':100,
#     #'X_P':100,
#     'S_O2':2,
#     'S_NO3':20,
#     'S_NH4':2,
#     'S_ALK':7*12,
#     }

# default_inf_kwargs ['asm2d'] = {
#     'concentrations': {
#       'S_I': 14,
#       'X_I': 26.5,
#       'S_F': 20.1,
#       'S_A': 94.3,
#       'X_S': 409.75,
#       'S_NH4': 31,
#       'S_N2': 0,
#       'S_NO3': 0.266, 
#       'S_PO4': 2.8,
#       'X_PP': 0.05,
#       'X_PHA': 0.5,
#       'X_H': 0.15,
#       'X_AUT': 0, 
#       'X_PAO': 0, 
#       'S_ALK':7*12,
#       },
#     'units': ('m3/d', 'mg/L'),
#     }

# default for 'exclude' - kinetic assay
default_init_conds_exc = {
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
    'S_N2':,
    'S_ALK':,
    'S_I':,
    'X_I':,
    'X_S':,
    'X_H':,
    'X_AUT':,    
    }

default_pm2asm2d_kwargs = dict(
    a_c=0.049, I_n=1500, arr_a=1.8e10, arr_e=6842, beta_1=2.90,
    beta_2=3.50, b_reactor=0.03, I_opt=2000, k_gamma=1e-5,
    K_N=0.1, K_P=1.0, K_A=6.3, K_F=6.3, rho=1.186, K_STO=1.566,
    f_CH_max=0.819, f_LI_max=3.249, m_ATP=15.835,
    mu_max=1.969, q_CH=0.594, q_LI=0.910,
    Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163,
    V_NH=0.254, V_NO=0.254, V_P=0.016, exponent=4,
    Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
    Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
    Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
    Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
    Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317, n_dark=0.7,
    
    f_SI=0.0, Y_H=0.625, f_XI_H=0.1, Y_A=0.24, f_XI_AUT=0.1,
    K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
    mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F_H=4.0,
    K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
    mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,
    path=None,
    ) # original baseline

# %%

# =============================================================================
# Batch system
# =============================================================================

def create_system(flowsheet=None, pm2asm2d_kwargs={}, init_conds={}, kind=''):
    flowsheet = flowsheet or qs.Flowsheet('pm2asm2d')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components
    pc.create_pm2asm2d_cmps()

    # Process models
    pm2asm2d_kwargs = pm2asm2d_kwargs or default_pm2asm2d_kwargs
    pm2asm2d = pc.PM2ASM2d(**pm2asm2d_kwargs)

    # Create unit operations

    if kind == 'include':
        PBR = su.BatchExperiment('PBR', model=pm2asm2d, exogenous_vars=(T_inc, I_inc))
        init_conds = init_conds or default_init_conds_inc

    else:
        PBR = su.BatchExperiment('PBR', model=pm2asm2d, exogenous_vars=(T_exc, I_exc))
        init_conds = init_conds or default_init_conds_exc

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
        export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)

if __name__ == '__main__':
    t = 7
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