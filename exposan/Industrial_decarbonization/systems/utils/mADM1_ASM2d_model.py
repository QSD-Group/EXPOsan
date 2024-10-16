#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 13:37:59 2024

@author: saumitrarai
"""

import os, qsdsan as qs
import numpy as np
from qsdsan import processes as pc, sanunits as su
from exposan.metab import create_system
from qsdsan.utils import (
    ospath, 
    time_printer,
    )
folder = ospath.dirname(__file__)

__all__ = (
    'create_WRRF_system',
    # 'default_asm1_kwargs', 'default_init_conds',
    'Temp', 
    )

# %%

# =============================================================================
# Parameters and util functions
# =============================================================================

Q_brewery_WERF = 785         # industrial wastewater [m3/d] flowrate is determined based on a constant

Temp = 273.15+20 # temperature [K]

bp3_metab = {'S_su': 9.472992246201887, 'S_aa': 4.024200833257111, 'S_fa': 63.96534516028177, 'S_va': 6.820554767531797, 'S_bu': 8.54620520480524, 'S_pro': 10.99219313097201, 'S_ac': 49.53408595467583, 'S_h2': 6.351890451616634e-05, 'S_ch4': 68.61711746845934, 'S_IC': 484.23496056838974, 'S_IN': -20.03079722956212, 'S_IP': 0, 'S_I': 97.38826145568594, 'X_ch': 74.4190824995273, 'X_pr': 52.19686027730525, 'X_li': 88.01751263684295, 'X_su': 504.58936901277445, 'X_aa': 95.26785584876515, 'X_fa': 101.69509191151278, 'X_c4': 306.07315898850646, 'X_pro': 90.49162763945084, 'X_ac': 253.14558016689116, 'X_h2': 191.54026396859808, 'X_I': 145.31476848616092, 'X_PHA': 0, 'X_PP': 0, 'X_PAO': 0, 'S_K': 0, 'S_Mg': 0, 'X_MeOH': 0, 'X_MeP': 0, 'S_cat': 40.0, 'S_an': 20.0}

#%%

def create_components():
     mADM1_cmps = pc.create_adm1_p_extension_cmps(False)
     cmps = qs.Components([*mADM1_cmps])
     cmps.compile()
     return cmps

def create_WRRF_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, sys_in = None, init_conds={}, lifetime=100, discount_rate=0.1,
                  aeration_processes=()):
   
    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    thermo_adm1 = qs.get_thermo()
    adm1 = qs.processes.ADM1_p_extension()
    
    ind_ww = qs.WasteStream('industrial_wastewater', T=Temp)
    
    ind_ww.set_flow_by_concentration(Q_brewery_WERF, 
                                      concentrations= bp3_metab,  
                                      units=('m3/d', 'mg/L'))
    
    cmps_asm2d = qs.processes.create_asm2d_cmps()
    qs.set_thermo(cmps_asm2d)
    thermo_asm2d = qs.get_thermo()
    asm2d = pc.ASM2d()
    
    J1 = su.mADM1toASM2d('J1', upstream = ind_ww,
                     # sys_in.outs[2], 
                     thermo=thermo_asm2d, isdynamic=True, adm1_model=adm1, asm2d_model=asm2d)
    sys = qs.System('mADM1_to_ASM2d', 
                    path=(J1,) )
    
    return sys
#%%

@time_printer
def run(t, method=None, sys_in=None, **kwargs):
    
    sys = create_WRRF_system(sys_in=sys_in)    
   
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        method=method,
        # print_t=True,
        **kwargs)
    
    S_O2, S_N2, S_NH4, S_NO3, S_PO4, S_F, S_A, S_I, S_ALK, X_I, X_S, X_H, X_PAO, X_PP, X_PHA, X_AUT, X_MeOH, X_MeP, H2O = sys.outs[0].conc 
    
    brewery_ww_asm2d = {'S_O2': S_O2,  
        'S_N2': S_N2,
        'S_NH4': S_NH4,     
        'S_NO3': S_NO3,
        'S_PO4': S_PO4,
        'S_F': S_F, 
        'S_A': S_A,
        'S_I': S_I, 
        'S_ALK': S_ALK, 
        'X_I': X_I, 
        'X_S': X_S, 
        'X_H': X_H, 
        'X_PAO': X_PAO, 
        'X_PP': X_PP,
        'X_PHA': X_PHA,
        'X_AUT': X_AUT,
        'X_MeOH': X_MeOH,
        'X_MeP': X_MeP
        }
    
    print(f'brewery_ww = {brewery_ww_asm2d}')
    
    return sys_in

    
if __name__ == '__main__':
    t = 1
    # method = 'RK45'
    method = 'RK23' 
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    # method = 'LSODA'
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    
    sys1 = create_WRRF_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, sys_in = None, init_conds={}, lifetime=100, discount_rate=0.1,
                      aeration_processes=())
    
    sys = run(t, method=method, sys_in = sys1)
    
    sys.diagram()