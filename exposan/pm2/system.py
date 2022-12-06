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
from qsdsan import processes as pc, sanunits as su, WasteStream, System
from qsdsan.utils import time_printer, load_data, get_SRT, ExogenousDynamicVariable as EDV
from exposan.pm2 import data_path

__all__ = (
    'biomass_IDs',
    'create_system',
    'default_pm2_kwargs',
    'default_inf_kwargs',
    'default_init_conds',
    'Q', 'Temp', 'V_mix', 'V_pbr', 'V_mem', 'V_ret',
    )

#%%

# =============================================================================
# Parameters and util functions
# =============================================================================

Q = 445.47           # influent flowrate [m3/d]
Temp = 286.13        # temperature [K]

V_mix = 66  
V_pbr = 77.49
V_mem = 7.03
V_ret = 6.14
biomass_IDs = ('X_ALG',)

# T_pbr, I_pbr = EDV.batch_init(os.path.join(data_path, 'exo_vars.xlsx'))

# T_mix = T
# I_mix = EDV('light_I_mix', function=lambda t: 0)

default_inf_kwargs = {
    'concentrations': {
        'X_CHL':0.0,
        'X_ALG':0.0,
        'X_CH':0.0,
        'X_LI':0.0,
        'S_CO2':30.0,
        'S_A':5.0,
        'S_F':5.0,
        'S_O2':15.0,
        'S_NH':10.0,
        'S_NO':1.0,
        'S_P':0.32,
        'X_N_ALG':0.0,
        'X_P_ALG':0.0,
        },
    'units': ('m3/d', 'mg/L'),
    }

default_pm2_kwargs = dict(
    a_c=0.049, I_n=1000, arr_a=1.8e10, arr_e=6842, beta_1=2.90, 
    beta_2=3.50, b_reactor=0.03, I_opt=1500, k_gamma=1e-5, 
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
    path=os.path.join(data_path, '_pm2.tsv'),
    )

default_init_conds = {
        'X_CHL':2.1,
        'X_ALG':420.0,
        'X_CH':4.2,
        'X_LI':1.0,
        'S_CO2':30.0,
        'S_A':5.0,
        'S_F':5.0,
        'S_O2':15.0,
        'S_NH':10.0,
        'S_NO':1.0,
        'S_P':0.32,
        'X_N_ALG':50,
        'X_P_ALG':10,
    }

def batch_init(sys, path, sheet):                        
    df = load_data(path, sheet)                   
    dct = df.to_dict('index')                       
    u = sys.flowsheet.unit #unit registry
    for k in [u.MIX, 
              u.PBR1, u.PBR2, u.PBR3, u.PBR4, u.PBR5, 
              u.PBR6, u.PBR7, u.PBR8, u.PBR9, u.PBR10, 
              u.PBR11, u.PBR12, u.PBR13, u.PBR14, u.PBR15, 
              u.PBR16, u.PBR17, u.PBR18, u.PBR19, u.PBR20,
              u.MEV, u.RET]:
        k.set_init_conc(**dct[k._ID[:3]])  

# %%

# =============================================================================
# Validation & Verification of PM2
# =============================================================================

def create_system(flowsheet=None, inf_kwargs={}, pm2_kwargs={}, init_conds={}):
    flowsheet = flowsheet or qs.Flowsheet('pm2')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components
    cmps = pc.create_pm2_cmps()
    
    # Exogeneous dynamic variables
    T_pbr, I_pbr = EDV.batch_init(os.path.join(data_path, 'exo_vars.xlsx'))

    T_mix = T_pbr
    I_mix = EDV('light_I_mix', function=lambda t: 0)
    
    # Streams
    SE = WasteStream('Secondary_effluent', T=Temp)
    inf_kwargs = inf_kwargs or default_inf_kwargs
    SE.set_flow_by_concentration(Q, **inf_kwargs)

    PHO = WasteStream('to_PBR', T=Temp)

    ME = WasteStream('to_membrane', T=Temp)
    INT = WasteStream('internal_recycle', T=Temp)

    RETEN = WasteStream('retentate', T=Temp)

    RE = WasteStream('to_return_tank', T=Temp)
    CE = WasteStream('to_centrifuge', T=Temp)

    CEN = WasteStream('centrate', T=Temp)
    ALG = WasteStream('harvested_biomass', T=Temp)

    RAA = WasteStream('return_activated_algae', T=Temp)
    
    TE = WasteStream('Tertiary_effluent', T=Temp)
    
    # Process models
    pm2_kwargs = pm2_kwargs or default_pm2_kwargs
    pm2 = pc.PM2(**pm2_kwargs)
    
    # Create unit operations
    MIX = su.CSTR('MIX', ins=[SE, RAA], outs=[PHO], V_max=V_mix,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_mix, I_mix))

    PBR1 = su.CSTR('PBR1', ins=MIX-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR2 = su.CSTR('PBR2', ins=PBR1-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR3 = su.CSTR('PBR3', ins=PBR2-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR4 = su.CSTR('PBR4', ins=PBR3-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR5 = su.CSTR('PBR5', ins=PBR4-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR6 = su.CSTR('PBR6', ins=PBR5-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR7 = su.CSTR('PBR7', ins=PBR6-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR8 = su.CSTR('PBR8', ins=PBR7-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR9 = su.CSTR('PBR9', ins=PBR8-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR10 = su.CSTR('PBR10', ins=PBR9-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR11 = su.CSTR('PBR11', ins=PBR10-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR12 = su.CSTR('PBR12', ins=PBR11-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR13 = su.CSTR('PBR13', ins=PBR12-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR14 = su.CSTR('PBR14', ins=PBR13-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR15 = su.CSTR('PBR15', ins=PBR14-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR16 = su.CSTR('PBR16', ins=PBR15-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR17 = su.CSTR('PBR17', ins=PBR16-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR18 = su.CSTR('PBR18', ins=PBR17-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR19 = su.CSTR('PBR19', ins=PBR18-0, V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))
    PBR20 = su.CSTR('PBR20', ins=PBR19-0, outs=[ME, INT], split=[0.51, 0.49], V_max=V_pbr/20,
                  aeration=None, suspended_growth_model=pm2, exogenous_vars=(T_pbr, I_pbr))

    MEM = su.Splitter('MEM', PBR20-0, outs=[TE, RETEN], split=0.39*(1-cmps.x))

    MEV = su.CSTR('MEV', ins=MEM-1, V_max=V_mem, aeration=None, suspended_growth_model=None)

    POST_MEM = su.Splitter('POST_MEM', MEV-0, outs=[RE, CE], split=0.97)

    CENT = su.Splitter('CENT', POST_MEM-1, outs=[CEN, ALG], split=0.975*(1-cmps.x))

    RET = su.CSTR('RET', ins=[PBR20-1, POST_MEM-0, CENT-0], outs=[RAA], V_max=V_ret,
                  aeration=None, suspended_growth_model=None)
    
    # System setup
    sys = System('sys',
                 path=(MIX, PBR1, PBR2, PBR3, PBR4, PBR5, PBR6, PBR7, PBR8, PBR9, PBR10, 
                 PBR11, PBR12, PBR13, PBR14, PBR15, PBR16, PBR17, PBR18, PBR19, PBR20, 
                 MEM, MEV, POST_MEM, CENT, RET), 
                 recycle=(RAA,))

    init_conds = init_conds or default_init_conds
    batch_init(sys, os.path.join(data_path, 'initial_conditions.xlsx'), 'default')
    
    sys.set_dynamic_tracker(MIX, PBR1, PBR20, RET, TE, ALG)
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
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)
    srt = get_SRT(sys, biomass_IDs)
    print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')

if __name__ == '__main__':
    t = 50
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