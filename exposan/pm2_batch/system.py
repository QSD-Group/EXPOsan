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

import os, numpy as np, qsdsan as qs
from qsdsan import processes as pc, sanunits as su, System
from qsdsan.utils import time_printer, ExogenousDynamicVariable as EDV
from exposan.pm2 import data_path

__all__ = (
    'create_system',
    'default_pm2_kwargs',
    'default_init_conds',
    'T', 'I',
    )

#%%

# =============================================================================
# Parameters and util functions
# =============================================================================

T, I = EDV.batch_init(os.path.join(data_path, 'exo_vars.xlsx'), 'linear')

default_init_conds = {
    'X_CHL':2.31,
    'X_ALG':461.18,
    'X_CH':11.34,
    'X_LI':59.58,
    'S_CO2':30.0,
    'S_A':5.0,
    'S_F':5.0,
    'S_O2':5.0,
    'S_NH':35.80,
    'S_NO':0.7,
    'S_P':0.36,
    'X_N_ALG':2.94,
    'X_P_ALG':10.55,
    }

default_pm2_kwargs = dict(
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
    path=None,
    )

# %%

# =============================================================================
# Validation & Verification of PM2
# =============================================================================

def create_system(flowsheet=None, pm2_kwargs={}, init_conds={}):
    flowsheet = flowsheet or qs.Flowsheet('pm2_batch')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components
    pc.create_pm2_cmps()
    
    # Process models
    pm2_kwargs = pm2_kwargs or default_pm2_kwargs
    pm2 = pc.PM2(**pm2_kwargs)
    
    # Create unit operations
    PBR = su.BatchExperiment('PBR', model=pm2, exogenous_vars=(T, I))
    
    init_conds = init_conds or default_init_conds
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
    t = 20
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
    