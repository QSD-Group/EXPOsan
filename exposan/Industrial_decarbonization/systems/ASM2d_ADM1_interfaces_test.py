#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:44:41 2024

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
    # get_SRT,
    # get_P_blower, get_power_utility, 
    # get_cost_sludge_disposal,
    # get_normalized_energy, 
    # get_daily_operational_cost, 
    # get_total_operational_cost, 
    # get_GHG_emissions_sec_treatment,
    # get_GHG_emissions_discharge,
    # get_GHG_emissions_electricity,
    # get_GHG_emissions_sludge_disposal,
    # get_CO2_eq_WRRF,
    # get_total_CO2_eq
    )
# from chaospy import distributions as shape

# from . import data_path
folder = ospath.dirname(__file__)


__all__ = (
    'biomass_IDs',
    'create_system',
    'default_asm2d_kwargs', 'domestic_ww', 'default_init_conds',
    'Q_domestic', 'Q_brewery', 'Q_ras', 'Q_was', 'Temp', 'V_an', 'V_ae', 
    )

# %%

# =============================================================================
# Parameters and util functions
# =============================================================================

Q_domestic = 605000 # influent flowrate [m3/d]
Q_brewery = 12500
Temp = 273.15+20 # temperature [K]
V_an = 25210     # anoxic zone tank volume [m3] 
V_ae = 25210     # aerated zone tank volume [m3] 
Q_was = 11356    # sludge wastage flowrate (0.04% of 69 MGD) [m3/day] 
Q_ras = 249837   # recycle sludge flowrate (0.96% of 69 MGD) [m3/day] 
biomass_IDs = ('X_H', 'X_AUT')
ammonia_ID = 'S_NH4'

domestic_ww = {
   'S_I': 14,
   'X_I': 26.5,
   'S_F': 20.1,
   'S_A': 94.3,
   'X_S': 395,
   'S_NH4': 31,
   'S_N2': 0,
   'S_NO3': 0.266, 
   'S_PO4': 2.8,
   'X_PP': 0.05,
   'X_PHA': 0.5,
   'X_H': 0.15,
   'X_AUT': 0, 
   'X_PAO': 0, 
   'S_ALK':7*12,
    }

bp2_metab = {'S_NH4': 170.73,
 'S_F': 41.47,
 'S_A': 39.406,
 'S_I': 61.3,
 'S_ALK': 452.7,
 'X_I': 250.965,
 'X_S': 1513.66}

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
      'S_I': 2.325e+02,
      'X_c': 2.171e+02,
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
    for k in sys.units:
        if sheet.startswith('B1'):
            if k.ID.startswith('O'): k.set_init_conc(**dct[k.ID])
        else:
            if k.ID.startswith('O'):
                k.set_init_conc(**dct['O'])
            elif k.ID.startswith('A'):
                k.set_init_conc(**dct['A'])           
    # c1s = {k:v for k,v in dct['C1_s'].items() if v>0}
    # c1x = {k:v for k,v in dct['C1_x'].items() if v>0}
    # tss = [v for v in dct['C1_tss'].values() if v>0]
    # u.C1.set_init_solubles(**c1s)
    # u.C1.set_init_sludge_solids(**c1x)
    # u.C1.set_init_TSS(tss)

#%%

def create_components():
     cmps = pc.create_asm2d_cmps(False)
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
    dom_ww = qs.WasteStream('domestic_wastewater', T=Temp)
    ind_ww = qs.WasteStream('industrial_wastewater', T=Temp)
    #inf_kwargs = inf_kwargs or default_inf_kwargs
    dom_ww.set_flow_by_concentration(Q_domestic, 
                                     concentrations=domestic_ww, 
                                     units=('m3/d', 'mg/L'))
    
    ind_ww.set_flow_by_concentration(Q_brewery, 
                                     concentrations= bp2_metab, 
                                     units=('m3/d', 'mg/L'))

    asm_kwargs = asm_kwargs or default_asm2d_kwargs
    
    Mixer = su.Mixer('Mixer', [dom_ww, ind_ww])
    
    cmps_adm1 = qs.processes.create_adm1_cmps()
    thermo_adm1 = qs.get_thermo()
    adm1 = qs.processes.ADM1()

    J1 = su.ASM2dtoADM1('J1', upstream= [Mixer-0], thermo=thermo_adm1, isdynamic=True, 
                      adm1_model=adm1)
    
    DG = su.AnaerobicCSTR(ID='DG', ins = J1.outs[0], outs= ['gas', 'sludge_DG'], 
                          model=adm1, thermo = thermo_adm1, V_liq = 3217, V_gas = 321.7) 
    
    J2 = su.ADM1toASM2d('J2', upstream = DG-1, thermo=thermo_asm2d, isdynamic=True, adm1_model=adm1)

    # thickener_perc and TSS_removal_perc based on WERF report
    DU = su.Centrifuge('DU', ins = J2.outs[0], outs = ['sludge_DU', 'eff_DU'], thermo = thermo_asm2d,
                        thickener_perc=23, TSS_removal_perc=95)
    
    sys = qs.System('metro_ASM2d', path=(Mixer,  J1, DG, J2, DU),)
    
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
    sys = run(t, method=method)
    
    sys.diagram()
    fs = sys.flowsheet.stream
    fu = sys.flowsheet.unit