# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from qsdsan import processes as pc, sanunits as su, System, WasteStream
from qsdsan.utils import time_printer

__all__ = (
    'biomass_IDs',
    'active_unit_IDs',
    'create_system',
    'default_asm1_kwargs',
    'default_inf_kwargs',
    'default_init_conds',
    'default_Q_was',
    'Q_inf', 'Q_ext', 'V_ae',
    )

#%%

# =============================================================================
# Parameters and util functions
# =============================================================================

Q_inf = 18446                               # influent flowrate [m3/d]
default_Q_was = 400                         # sludge wastage flowrate [m3/d]
Q_ext = 18446                               # external recycle flowrate [m3/d]

V_ae = 1333                                 # aerated tank volume [m3/d]

biomass_IDs = ('X_BH', 'X_BA')
active_unit_IDs = ("O1", )

valid_temps = (10, 20)
default_asm1_kwargs = dict.fromkeys(valid_temps)
default_asm1_kwargs[10] = dict(
    Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
    mu_H=3.0, K_S=20.0, K_O_H=0.2, K_NO=0.5, b_H=0.2, 
    eta_g=0.8, eta_h=0.4, k_h=1.0, K_X=0.01, mu_A=0.3, 
    K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.04, fr_SS_COD=0.75, 
    )    # at 10 degree C
default_asm1_kwargs[20] = dict(
    Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
    mu_H=6.0, K_S=20.0, K_O_H=0.2, K_NO=0.5, b_H=0.62, 
    eta_g=0.8, eta_h=0.4, k_h=3.0, K_X=0.03, mu_A=0.8, 
    K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.08, fr_SS_COD=0.75, 
    )    # at 20 degree C

# default_asm1_kwargs = dict(
#     Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
#     mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3, 
#     eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5, 
#     K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75, 
#     )    # bsm default

default_inf_kwargs = {
    'concentrations': {                                            
        'S_S':69.5,
        'X_BH':28.17,
        'X_S':202.32,
        'X_I':51.2,
        'S_NH':31.56,
        'S_I':30,
        'S_ND':6.95,
        'X_ND':10.59,
        'S_ALK':7*12,
      },
    'units': ('m3/d', 'mg/L'),                
    }

default_init_conds = {
        'S_I':30,
        'S_S':5,
        'X_I':1000,
        'X_S':100,
        'X_BH':500,
        'X_BA':100,
        'X_P':100,
        'S_O':2,
        'S_NO':20,
        'S_NH':2,
        'S_ND':1,
        'X_ND':1,
        'S_ALK':7*12,
    } 
    
# %%

# =============================================================================
# CSTR with ideal clarifier
# =============================================================================

def create_system(                  
        flowsheet=None, 
        inf_kwargs={},
        asm_kwargs={},
        init_conds={},
        Q_was=None,
        Temp=20,                     # default = 20 degree C
        KLa=240,
        ):    
    
    # Users can change influent composition, asm1 parameters, 
    # initial condition, Q_was, Temp, KLa

    flowsheet = flowsheet or qs.Flowsheet('edu')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and stream
    pc.create_asm1_cmps()
    
    if Temp == 10:
        asm_kwargs = asm_kwargs or default_asm1_kwargs[Temp]
        Temp = 273.15+10
    elif Temp == 20:
        asm_kwargs = asm_kwargs or default_asm1_kwargs[Temp]
        Temp = 273.15+20
    else:
        raise ValueError('`Temp` can only be either 10 or 20.')
   
    asm = pc.ASM1(**asm_kwargs)
       
    influent = WasteStream('influent', T=Temp)
    inf_kwargs = inf_kwargs or default_inf_kwargs
    influent.set_flow_by_concentration(Q_inf, **inf_kwargs)   

    effluent = WasteStream('effluent', T=Temp)
    solid = WasteStream('solid', T=Temp)
    ext_recycle = WasteStream('recycle', T=Temp)
    wastage = WasteStream('wastage', T=Temp)    
    
    Q_was = Q_was or default_Q_was
    if Q_was > 11400:
        raise ValueError('`Flow rate of wastage stream` can only be <11400.')
   
    # Process model
    aer = pc.DiffusedAeration('aer', DO_ID='S_O', KLa=KLa, DOsat=8.0, V=V_ae)
    
    # Create unit operations
    O1 = su.CSTR('O1', ins=[influent, ext_recycle], V_max=V_ae, aeration=aer,
                 DO_ID='S_O', suspended_growth_model=asm)
    
    C1 = su.IdealClarifier('C1', ins=O1-0, outs=[effluent, solid],
                           sludge_flow_rate=Q_ext+Q_was, solids_removal_efficiency=1)
    
    S1 = su.Splitter('S1', ins=C1-1, outs=[ext_recycle, wastage], split=1-Q_was/Q_ext)

    # System setup
    sys = System('example_system', path=(O1, C1, S1), recycle=(ext_recycle))    
    
    init_conds = init_conds or default_init_conds
    O1.set_init_conc(**init_conds)
    
    sys.set_dynamic_tracker(influent, effluent, ext_recycle, wastage, O1, C1, S1)          
    sys.set_tolerance(rmol=1e-6)
    
    return sys

#%%
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
    # srt = get_SRT(sys, biomass_IDs, wastage, active_unit_IDs)
    # print(f'Estimated SRT assuming at steady state is {round(srt, 8)} days')

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