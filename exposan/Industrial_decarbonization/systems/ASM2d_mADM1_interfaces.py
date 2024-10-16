#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 14:44:41 2024

@author: saumitrarai
"""

import qsdsan as qs
from qsdsan import processes as pc, sanunits as su
# , Model as mod
# import numpy as np
from qsdsan.utils import (
    ospath, 
    time_printer, 
    load_data, 
    )
# from chaospy import distributions as shape
from exposan.metab import UASB, flex_rhos_adm1

# from . import data_path
folder = ospath.dirname(__file__)


__all__ = (
    'create_system',
    'default_asm2d_kwargs', 'influent_ww', 'default_init_conds',
    'Q_influent', 'Temp', 
    )

# %%

# =============================================================================
# Parameters and util functions
# =============================================================================

Q_influent = 139.70 # influent flowrate [m3/d]
Temp = 273.15+35 # temperature [K]

# Using WAS from Metropolitan configuration
influent_ww = {
            'S_NH4': 2.359e+01, 
            'S_NO3': 6.170e+00, 
            'S_PO4': 2.558e+00, 
            'S_F': 1.437e+01, 
            'S_A': 6.485e+01, 
            'S_I': 1.496e+01, 
            'S_ALK': 8.304e+01, 
            'X_I': 9.225e+03, 
            'X_S': 5.276e+04, 
            'X_H': 9.410e+03, 
            'X_PAO': 2.621e+01, 
            'X_PP': 1.471e+01, 
            'X_PHA': 6.170e+01, 
            'X_AUT': 5.080e+02, 
                }

default_asm2d_kwargs = dict(iN_SI=0.01, iN_SF=0.03, iN_XI=0.02, iN_XS=0.04, iN_BM=0.07,
            iP_SI=0.0, iP_SF=0.01, iP_XI=0.01, iP_XS=0.01, iP_BM=0.02,
            iTSS_XI=0.75, iTSS_XS=0.75, iTSS_BM=0.9,
            f_SI=0.0, Y_H=0.625, f_XI_H=0.1,
            Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, f_XI_PAO=0.1,
            Y_A=0.24, f_XI_AUT=0.1,
            K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
            mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F=4.0,
            K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
            q_PHA=3.0, q_PP=1.5, mu_PAO=1.0, eta_NO3_PAO=0.6, b_PAO=0.2, b_PP=0.2,
            b_PHA=0.2, K_O2_PAO=0.2, K_NO3_PAO=0.5, K_A_PAO=4.0, K_NH4_PAO=0.05,
            K_PS=0.2, K_P_PAO=0.01, K_ALK_PAO=0.1,
            K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
            mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,
            k_PRE=1.0, k_RED=0.6, K_ALK_PRE=0.5, 
            # path=os.path.join(data_path, '_asm2d.tsv'),
            )

default_init_conds = {
        'S_F':5,
        'S_A':2,
        'X_I':1000,
        'X_S':100,
        'X_H':500,
        'X_AUT':100,
        #'X_P':100,
        'S_O2':2,
        'S_NO3':20,
        'S_NH4':2,
        'S_ALK':7*12,
    }

steady_state_ad_init_conds = {
      'S_su': 1.094e+01,
       # 'S_aa': 4.713e+00, # original ss
      # 'S_aa': 4.713e+01,
       'S_aa': 4.713e+03,
      'S_fa': 9.150e+01,
      'S_va': 9.512e+00,
      'S_bu': 1.255e+01,
      'S_pro': 1.486e+01,
      'S_ac': 3.749e+01,
       'S_h2': 2.193e-04,
      'S_ch4': 7.692e+01,
      'S_IC': 6.708e+02,
      'S_IN': 1.035e+03,
      # --
      'S_IP': 1.035e+00, #random
      # --
      'S_I': 2.325e+02,
      # --
      # 'X_c': 2.171e+02,
      # --
       # 'X_ch': 6.432e+01, # original ss
      # 'X_ch': 6.432e+00, 
       'X_ch': 6.432e-01,
      'X_pr': 6.719e+01,
       # 'X_li': 1.436e+02, # original ss
      # 'X_li': 1.436e+01,
       'X_li': 1.436e+00,
      'X_su': 1.116e+03,
      'X_aa': 8.755e+02,
      'X_fa': 1.274e+03,
      'X_c4': 3.730e+02,
      'X_pro': 1.751e+02,
      'X_ac': 1.405e+03,
      'X_h2': 6.820e+02,
      'X_I': 2.413e+04,
      # --
      'X_PHA': 10.000e-01, 
      'X_PP':  10.000e-01,  
      
      'X_PAO': 10.000e+00, 
      
      'S_K': 10.000e-001, 
      'S_Mg': 10.000e-001, 
      'X_MeOH': 10.000e-001, 
      'X_MeP': 10.000e-001, 
      # --
       'S_cat': 0.000e+00, 
       'S_an': 4.772e+00
      }

def batch_init(sys, path, sheet):
    # isa = isinstance
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    # u.DG.set_init_conc(**default_ad_init_conds)
    u.DG.set_init_conc(**steady_state_ad_init_conds)
    # for k in sys.units:
    #     if sheet.startswith('B1'):
    #         if k.ID.startswith('O'): k.set_init_conc(**dct[k.ID])
    #     else:
    #         if k.ID.startswith('O'):
    #             k.set_init_conc(**dct['O'])
    #         elif k.ID.startswith('A'):
    #             k.set_init_conc(**dct['A'])           
    # c1s = {k:v for k,v in dct['C1_s'].items() if v>0}
    # c1x = {k:v for k,v in dct['C1_x'].items() if v>0}
    # tss = [v for v in dct['C1_tss'].values() if v>0]
    # u.C1.set_init_solubles(**c1s)
    # u.C1.set_init_sludge_solids(**c1x)
    # u.C1.set_init_TSS(tss)

#%%

def create_components():
     cmps = pc.create_asm2d_cmps(False)
     cmps.X_I.i_N = 0.0600327162
     cmps.X_I.i_P = 0.01
     cmps.S_I.i_P = 0.01
     cmps.S_I.i_N = 0.0600327162
     cmps.refresh_constants()
     # asm2d = pc.ASM2d(iN_XI = 0.06)
     cmps.compile()
     return cmps

def create_system(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={}, lifetime=100, discount_rate=0.1,
                  aeration_processes=()):
    # flowsheet = flowsheet or qs.Flowsheet('bsm1')
    # qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    thermo_asm2d = qs.get_thermo()
    
    qs.set_thermo(cmps)
    sludge = qs.WasteStream('sludge', T=Temp)
    
    #inf_kwargs = inf_kwargs or default_inf_kwargs
    sludge.set_flow_by_concentration(Q_influent, 
                                     concentrations=influent_ww, 
                                     units=('m3/d', 'mg/L'))

    asm_kwargs = asm_kwargs or default_asm2d_kwargs
    
    cmps_adm1 = qs.processes.create_adm1_p_extension_cmps()
    cmps_adm1.X_PAO.i_N = 0.07
    cmps_adm1.X_PAO.i_P = 0.02
    cmps_adm1.refresh_constants()
    thermo_adm1 = qs.get_thermo()
    adm1 = qs.processes.ADM1_p_extension()
    
    # adm1 = qs.processes.ADM1_p_extension(flex_rate_function=flex_rhos_adm1)

    J1 = su.ASM2dtomADM1('J1', upstream= [sludge,], thermo=thermo_adm1, isdynamic=True, 
                      adm1_model=adm1)
    
    DG = su.AnaerobicCSTR(ID='DG', ins = J1.outs[0], outs= ['gas', 'sludge_DG'], 
                          model=adm1, thermo = thermo_adm1, V_liq = 3217, V_gas = 321.7) 
    
    # Volume of AD based on WERF report
    # DG = UASB(ID='DG', ins = J1.outs[0], outs= ['gas', 'sludge_DG'], 
    #                       model=adm1, thermo = thermo_adm1, V_liq = 3217, V_gas = 321.7, 
    #                       pH_ctrl = 6.8, fraction_retain=0) 
    
    cmps = create_components()
    qs.set_thermo(cmps)
    thermo_asm2d = qs.get_thermo()
    
    J2 = su.mADM1toASM2d('J2', upstream = DG-1, thermo=thermo_asm2d, isdynamic=True, adm1_model=adm1)

    sys = qs.System('metro_ASM2d', path=(J1, 
                                          DG, 
                                         J2),)
    
    return sys
#%%

@time_printer
def run(t, method=None, **kwargs):
    sys = create_system()    
    
    path = ospath.join(folder, "/Users/saumitrarai/Desktop/Research/MS_research/metro/data/initial_conditions_ASM2d.xlsx")    
    
    batch_init(sys, path, 
                sheet='ss_alt')
                # sheet='t=10')
        
    # RAS = sys.flowsheet.stream.RAS
    # C1 = sys.flowsheet.unit.C1
    # sys.set_dynamic_tracker(RAS, C1)
    
    sys.set_dynamic_tracker(*sys.products)
    
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        method=method,
        # print_t=True,
        **kwargs)
    
    # sys._setup()
    # sys.converge() 
    
    # fu = sys.flowsheet.unit
    
    # fu.DG.set_init_conc(**steady_state_ad_init_conds)
    
    # fu.DG.simulate(
    #     state_reset_hook='reset_cache',
    #     t_span=(0,t),
    #     method=method,
    #     # print_t=True,
    #     **kwargs)

    return sys
    
if __name__ == '__main__':
    t = 100
    # method = 'RK45'
    # method = 'RK23' 
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    sys = run(t, method=method)
    
    sys.diagram()
    fs = sys.flowsheet.stream
    fu = sys.flowsheet.unit