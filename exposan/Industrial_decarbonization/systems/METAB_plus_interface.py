#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 11:10:50 2024

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

Q_brewery = 785         # industrial wastewater [m3/d] flowrate is determined based on a constant

Temp = 273.15+20 # temperature [K]

baseline_minus_two_influent = {
    'S_su':0.18,
    'S_aa':0.18,
    'S_fa':0.18,
    'S_va':0.15,
    'S_bu':0.15,
    'S_pro':0.15,
    'S_ac':0.15,
    'S_IC':0.001,
    # 'S_IC':0.04*12,
    # 'S_IN':0.01*14,
    'S_I':0.024,
    'X_c':0.12,
    'X_ch':0.36,
    'X_pr':0.6,
    'X_li':0.3,
    'X_I':0.03, 
    'S_cat':0.04, 
    'S_an':0.02
    }

baseline_minus_two_metab =  {'S_su': 4.282270896827433, 'S_aa': 2.075051137773224, 'S_fa': 26.717857012142748, 'S_va': 3.470865426952685, 'S_bu': 4.171021239982134, 'S_pro': 4.561508532171465, 'S_ac': 18.822138719107855, 'S_h2': 2.968335491317047e-05, 'S_ch4': 42.87455256316345, 'S_IC': 185.45489347815246, 'S_IN': 55.86899262735419, 'S_I': 48.624383607039434, 'X_c': 84.2670061436853, 'X_ch': 10.847816256785055, 'X_pr': 17.5144829246916, 'X_li': 9.605057719225021, 'X_su': 49.35837838369387, 'X_aa': 57.94957534126758, 'X_fa': 4.773428276205071, 'X_c4': 31.2709036785768, 'X_pro': 11.237613626397083, 'X_ac': 20.191859875671863, 'X_h2': 20.99400813763462, 'X_I': 78.15562262059682, 'S_cat': 40.00000000000001, 'S_an': 20.000000000000004}
#%%

def create_components():
     adm1_cmps = pc.create_adm1_cmps(False)
     cmps = qs.Components([*adm1_cmps])
     cmps.compile()
     return cmps

def create_WRRF_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, sys_in = None, init_conds={}, lifetime=100, discount_rate=0.1,
                  aeration_processes=()):
   
    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    thermo_adm1 = qs.get_thermo()
    adm1 = qs.processes.ADM1()
    
    ind_ww = qs.WasteStream('industrial_wastewater', T=Temp)
    
    ind_ww.set_flow_by_concentration(Q_brewery, 
                                      concentrations= baseline_minus_two_metab,  
                                      units=('m3/d', 'mg/L'))

    
    cmps_asm1 = qs.processes.create_asm1_cmps()
    qs.set_thermo(cmps_asm1)
    thermo_asm1 = qs.get_thermo()
    
    J1 = su.ADMtoASM('J2', upstream = ind_ww,
                     # sys_in.outs[2], 
                     thermo=thermo_asm1, isdynamic=True, adm1_model=adm1)
    sys = qs.System('ADM_to_ASM', 
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
    
    S_I, S_S, X_I, X_S, X_BH, X_BA, X_P, S_O, S_NO, S_NH, S_ND, X_ND, S_ALK, S_N2, H2O = sys.outs[0].conc
    
    brewery_ww_after_METAB = { 'S_I' : S_I,  
        'S_S':S_S,
        'X_I':X_I,     
        'X_S': X_S,
        'X_BH': X_BH,
        'X_BA': X_BA, 
        'X_P': X_P,
        'S_O': S_O, 
        'S_NO': S_NO, 
        'S_NH': S_NH, 
        'S_ND': S_ND, 
        'X_ND': X_ND, 
        'S_ALK': S_ALK, 
        'S_N2': S_N2
        }
    
    print(f'brewery_ww_after_METAB = {brewery_ww_after_METAB}')
    
    return sys

def run_METAB():
    
    sys = create_system(tot_HRT=1, reactor_type='PB', inf_concs=baseline_minus_two_influent, Q=Q_brewery)
    sys.simulate(state_reset_hook='reset_cache', t_span=(0,400), method='BDF')
    
    return sys

    
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
    
    sys1 = run_METAB()
    sys1.diagram()
    
    S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I, X_c, X_ch, X_pr, \
    X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, S_cat, S_an, H2O = sys1.outs[2].conc
    
    ww_after_METAB = {'S_su': S_su, 
                      'S_aa': S_aa, 
                      'S_fa': S_fa, 
                      'S_va': S_va, 
                      'S_bu': S_bu, 
                      'S_pro': S_pro, 
                      'S_ac': S_ac, 
                      'S_h2': S_h2, 
                      'S_ch4': S_ch4, 
                      'S_IC': S_IC, 
                      'S_IN': S_IN - (sys1.outs[2].TN - sys1.ins[0].TN), 
                      'S_I': S_I, 
                      'X_c': X_c, 
                      'X_ch': X_ch, 
                      'X_pr': X_pr, 
                      'X_li': X_li, 
                      'X_su': X_su, 
                      'X_aa': X_aa, 
                      'X_fa': X_fa, 
                      'X_c4': X_c4, 
                      'X_pro': X_pro, 
                      'X_ac': X_ac, 
                      'X_h2': X_h2, 
                      'X_I': X_I, 
                      'S_cat': S_cat, 
                      'S_an': S_an}
    
    print(f'ww_after_METAB = {ww_after_METAB}')
    
    sys = run(t, method=method, sys_in = sys1)
    
    sys.diagram()

