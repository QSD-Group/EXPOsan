#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
This module is developed by:
    Siqi Tang <siqit@outlook.com>
    Yuyao Huang <yuyaoh2@illinois.edu>
    Aaron Marszewski <aaronpm3@illinois.edu>
    Rishabh Puri <rp34@illinois.edu>

This python file is used to perform uncertainty and sensitivity analysis for Enviroloo Clear Reinvented Toilet system.
'''

# %% 
import qsdsan as qs 
from qsdsan import (
#     Flowsheet, main_flowsheet,
#     WasteStream,
#     sanunits as su,
      ImpactItem, 
#     LCA, TEA, System,
     )
# from qsdsan.sanunits import Trucking
# from chaospy import distributions as shape
# from qsdsan.utils import clear_lca_registries
# from qsdsan.sanunits._excretion import ExcretionmASM2d

# from qsdsan.utils import load_components, set_thermo

from qsdsan.utils import (
    ospath, 
    time_printer, 
    load_data, 
    get_SRT,
    )

# from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
# from exposan.utils import add_V_from_rho
# from exposan.bwaise import create_components as create_bw_components
# import numpy as np
from qsdsan import Components, processes as pc, sanunits as su#, Model as mod
# from chaospy import distributions as shape


# from exposan.utils import add_fugitive_items
from exposan.enviroloo import _units as elu
from exposan.enviroloo import (
    data_path,
    # _load_components,
    _load_lca_data,
    # discount_rate,
    # get_decay_k,
    # get_toilet_users,
    # max_CH4_emission,
    # operator_daily_wage,
    ppl, baseline_ppl, scale_factor
    # get_tanker_truck_fee
    # update_resource_recovery_settings,
    )
# from exposan.enviroloo._EL_pumps import (
#     LiftPump, AgitationPump, DosingPump, ReturnPump, SelfPrimingPump, 
#     AirDissolvedPump, MicroBubblePump, ClearWaterPump,)
from exposan import enviroloo as el

folder = ospath.dirname(__file__)

__all__ = ('create_systemEL',)
#%%

'''
Name notes:

To make programming more convenient, we use the following names for the units in the system:

Control Room Housing: ELH 
Collection tank: CT
Primary clarifier: PC
Anoxic tank: AnoT
Aerobic tank: AerT
Membrane tank: MemT
Clear water tank: CWT
Pressure tank: PT
Primary clarifier return pump: P_PC_return
Glucose agitation pump: P_AnoT_agitation
Glucose dosing pump: P_AnoT_dosing
Anoxic mixing pump: P_AnoT_mixing.
PAC agitation pump: P_AerT_agitation
PAC dosing pump: P_AerT_dosing
Self-priming pump: P_MT_selfpriming
Lift pump: P_CT_lift
Aerobic blower: B_AerT
Membrane blower: B_MemT
Clear water pump: P_CWT
Ozone generator: O3_gen
Micro bubble pump (Ozone dosing pump): P_O3_dosing
Air dissolving pump: P_AirDissolved

'''
Temp = 273.15+20 # temperature [K]
# scale_factor = 10 #scale factor for flow, baseline is 100 ppl, imported from _init
Q_w = 3.63 * scale_factor # m3/day  #Assuming Toilet flows to be 60lpcd and Initial population to be 100 6m^3/day for 100 ppl, 60 for 1000 ppl
Q_ras = 1.815 * scale_factor # m3/day # 200 l/min for nitrate pump 3m^3/day for 100 ppl, 30 for 1000 ppl
Q_was = 0.05 * 24 # 0.05*12
# Q_was = 280 # m3/day # 200 l/min for sludge return pump
biomass_IDs = ('X_H', 'X_AUT', 'X_PAO')
toilet_waste={
    'S_NH4':  19, # 13.7
    # 'S_NO3': 0.2,     # Nitrate
    'S_NO3': 1.6,
    # 'S_PO4': 2.3,
    'S_PO4': 16.7, # 3.5
    'S_O2': 3.0,      # Dissolved oxygen
    'S_F': 65,  # 19.1    # Fermentable COD
    'S_I': 25,  # 12.8    # Inert soluble COD
    'X_S': 40,  #29.3      # Slowly biodegradable particulates = TSS
    'S_IC': 0.148,
    'S_K': 0.0694,
    'S_Mg': 0.00833,
    'S_Ca': 0.0117,
    'S_Na': 0.775,
    'S_Cl': 0.687,
    # 'H2O':   55.4
    }

#%%

'''
Excretion mASM2d Values
[0] urine
flow (g/hr): S_NH4  0.241
                S_PO4  0.0243
                S_F    0.302
                S_IC   0.148
                S_K    0.0694
                S_Mg   0.00833
                S_Ca   0.0117
                S_Na   0.775
                S_Cl   0.687
                H2O    55.4
    WasteStream-specific properties:
     pH         : 7.0
     COD        : 5198.4 mg/L
     BOD        : 3727.2 mg/L
     TC         : 4203.8 mg/L
     TOC        : 1655.3 mg/L
     TN         : 4317.1 mg/L
     TP         : 446.9 mg/L
     TK         : 1192.3 mg/L
[1] feces
phase: 'l', T: 298.15 K, P: 101325 Pa
flow (g/hr): S_NH4  0.00685
                S_PO4  0.0121
                S_A    0.472
                S_IC   0.0237
                S_K    0.0244
                S_Mg   0.0104
                X_S    0.817
                S_Ca   0.0792
                S_Na   0.124
                S_Cl   0.11
                H2O    8.85
    WasteStream-specific properties:
     pH         : 7.0
     COD        : 123414.2 mg/L
     BOD        : 77768.8 mg/L
     TC         : 44139.3 mg/L
     TOC        : 41869.9 mg/L
     TN         : 3278.4 mg/L
     TP         : 1591.1 mg/L
     TK         : 2332.9 mg/L
     TSS        : 58681.8 mg/L

'''
# default_masm2d_kwargs = dict(pH_ctrl=7.0, 
#                 f_SI=0.0, Y_H=0.625, Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, Y_A=0.24, 
#                 f_XI_H=0.1, f_XI_PAO=0.1, f_XI_AUT=0.1,
#                 k_h=3.0, mu_H=6.0, mu_PAO=1.0, mu_AUT=1.0, 
#                 q_fe=3.0, q_PHA=3.0, q_PP=1.5, 
#                 b_H=0.4, b_PAO=0.2, b_PP=0.2, b_PHA=0.2, b_AUT=0.15, 
#                 eta_NO3=0.6, eta_fe=0.4, eta_NO3_H=0.8, eta_NO3_PAO=0.6, 
#                 eta_NO3_Hl=0.5, eta_NO3_PAOl=0.33, eta_NO3_PPl=0.33, eta_NO3_PHAl=0.33, eta_NO3_AUTl=0.33,
#                 K_O2=0.2, K_O2_H=0.2, K_O2_PAO=0.2, K_O2_AUT=0.5, 
#                 K_NO3=0.5, K_NO3_H=0.5, K_NO3_PAO=0.5, K_NO3_AUT=0.5, 
#                 K_X=0.1, K_F=4.0, K_fe=4.0, K_A_H=4.0, K_A_PAO=4.0, 
#                 K_NH4_H=0.05, K_NH4_PAO=0.05, K_NH4_AUT=1.0, 
#                 K_P_H=0.01, K_P_PAO=0.01, K_P_AUT=0.01, K_P_S=0.2, 
#                 K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
#                 # k_mmp=(5.0, 300, 0.05, 150, 50, 1.0, 1.0),
#                 # pKsp=(6.45, 13.16, 5.8, 23, 7, 21, 26),
#                 # k_mmp=(0.024, 120, 0.024, 72, 0.024, 0.024, 0.024),  # Flores-Alsina 2016
#                 # pKsp=(8.3, 13.6, 18.175, 28.92, 7.46, 18.2, 37.76),  # Flores-Alsina 2016
#                 k_mmp=(8.4, 240, 1.0, 72, 1.0, 1.0e-5, 1.0e-5),              # MATLAB
#                 pKsp=(8.45, 13.5, 5.7, 29.1, 7.4, 18.2, 26.4),               # MINTEQ (except newberyite), 20 C    
#                 K_dis=(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
#                 K_AlOH=0.001, K_FeOH=0.001, 
#                 pKa=(14, 9.25, 6.37, 10.32, 2.12, 7.21, 12.32, 4.76),
#                 )
# toilet_waste={'S_NH4':  0.241,
#                 'S_PO4': 0.0243,
#                 'S_F':   0.302,
#                 'S_IC':   0.148,
#                 'S_K':   0.0694,
#                 'S_Mg':   0.00833,
#                 'S_Ca':   0.0117,
#                 'S_Na':   0.775,
#                 'S_Cl':   0.687,
#                 'H2O':   55.4
    
#     }
# feces_waste= {'S_NH4':  0.00685,
#                 'S_PO4':  0.0121,
#                 'S_A':    0.472,
#                 'S_IC':   0.0237,
#                 'S_K':   0.0244,
#                 'S_Mg':   0.0104,
#                 'X_S':   0.817,
#                 'S_Ca':   0.0792,
#                 'S_Na':   0.124,
#                 'S_Cl':   0.11,
#                 'H2O':    8.85
#                 }
#%% Create Universal Units and Functions
def batch_init(sys, path, sheet):
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    for k in sys.units:
        print(f'k={k}')
        if k.ID.startswith(('O', 'A', 'B')):
            k.set_init_conc(**dct[k.ID])

# %% Create EnviroLoo Clear system

def create_components(set_thermo = True
                      #adjust_MW_to_measured_as=False
                      ):
    # bw_cmps = create_bw_components(set_thermo=False)
    masm2d_cmps = pc.create_masm2d_cmps(set_thermo=False)
    cmps = Components([*masm2d_cmps, 
                       # Tissue, WoodAsh, H2O,
                       ])
    #cmps.compile()
    cmps.compile(ignore_inaccurate_molar_weight=True)
    if set_thermo: qs.set_thermo(cmps)
    return cmps
#%%
    
def create_systemEL(flowsheet=None, inf_kwargs={}, masm_kwargs={}, init_conds={},
                  aeration_processes=()):
    # Components and stream
    _load_lca_data()
    cmps = create_components()
    toilet_ins = qs.WasteStream('toilet_waste', T=Temp)
    toilet_ins.set_flow_by_concentration(Q_w*0.9, 
                                         concentrations=toilet_waste, 
                                         units=('m3/d', 'mg/L'))
    
   
    
    # PAC = qs.WasteStream('PAC_Dose', X_AlOH= 0.24, units='kg/hr', T=Temp) # 0.1207
    
    # Create effluent_CT WasteStream to match toilet_ins concentration
    effluent_CT = qs.WasteStream('effluent_CT', T=Temp)
    effluent_CT.set_flow_by_concentration(Q_w*0.9,
                                          concentrations=toilet_waste, 
                                          units=('m3/d', 'mg/L'))
    
   # WasteStream('CH4', phase='g', stream_impact_item=item)

    
    masm2d = pc.mASM2d(**masm_kwargs)
    # masm2d_A1 = pc.mASM2d(mu_H=6.0, eta_NO3_H=0.8, K_F=0.05, K_NO3_H=0.05, K_O2_H=0.01, K_O2=0.01, K_O2_PAO=0.01, K_O2_AUT=0.01)
    
    kwargs_O = dict(V_max=6.35, aeration=2, DO_ID='S_O2', suspended_growth_model=masm2d)
    kwargs_1 = dict(V_max=6.35, aeration=None, DO_ID=None, suspended_growth_model=masm2d)
    kwargs_2 = dict(V_max=2.89, aeration=2, DO_ID='S_O2', suspended_growth_model=masm2d)
    
    
    # CT = su.Mixer('CT', 
    #               # ins=(toilet_ins, 'sludge_PC', 'flushing_water_CT', 'spill_PC'), 
    #               ins=(toilet_ins, 'sludge_PC',), 
    #               outs=('effluent_CT'),
    #               )
    CT = elu.EL_CT(
        'CT', ins=(toilet_ins), outs=(effluent_CT),
        isdynamic=True, V_max=10, aeration=None, suspended_growth_model=None, ppl=ppl, baseline_ppl=baseline_ppl
        )
    
    # PC = su.PrimaryClarifier('PC', ins=[CT-0, 'sludge_AeroT'],
    #                          outs=('effluent_PC_total', sludge_PC), isdynamic=True,
    #                          init_with='WasteStream', thermo=thermo_masm2d,
    #                          sludge_flow_rate=280, 
    #                          solids_removal_efficiency=0.6)
    
    PC = elu.EL_PC('PC', ins=(CT-0, 'RAS_PC'), outs=('effluent_PC_total', 'sludge_PC'),
                   ppl=ppl, baseline_ppl=baseline_ppl,
                   solids_removal_efficiency=0.85,
                   isdynamic=True,
                   # thermo=thermo_masm2d,
                   sludge_flow_rate=Q_w*0.01)  # 0.00116 m3/hr
    

    # S3 = su.Splitter('S3', ins = PC-0, outs = ['effluent_PC', 3-CT], split= 0.8)
    
    mixing_ratio = 0.35 #kg/L for glucose and PAC
    
    #Glucose = qs.WasteStream('Glucose_Dose', S_F= 0.05, units='kg/hr', T=Temp) # 0.0805 0.05 kg COD/hr as S_F fermentable solute 
    Glucose = qs.WasteStream('Glucose_Dose', T=Temp) # 0.0805 0.05 kg COD/hr as S_F fermentable solute 
    
    
    A1 = elu.EL_Anoxic('A1', ins=(PC-0, 'RAS_A1', Glucose), outs=('effluent_AnoxT',),
                       isdynamic=True,  ppl = ppl, baseline_ppl = baseline_ppl, **kwargs_1 
                       # aeration=None, 
                       # DO_ID='S_O2', suspended_growth_model=masm2d,  **kwargs_1
                       # W_tank= 2.09,
                                       
                        # V_max= 7.33, 
                        )
    item = ImpactItem.get_item('Glucose_item').copy(f'A1_glucose_item', set_as_source=True)
    A1.ins[2].stream_impact_item = item      #effluent was
    
    PAC = qs.WasteStream('PAC_Dose', T=Temp) #  as Al_OH
    
    O1 = elu.EL_Aerobic('O1', ins=(A1-0, PAC), outs=('effluent_AeroT',),
                        isdynamic=True,  ppl = ppl, baseline_ppl = baseline_ppl, **kwargs_O
                        # aeration = 2, suspended_growth_model=asm2d,
                        #  # ppl = ppl, baseline_ppl = 100,
                        #  W_tank= 2.09,
                        #  V_max=7.33,**kwargs_0, 
                         )
    item = ImpactItem.get_item('PAC_item').copy(f'O1_PAC_item', set_as_source=True)
    O1.ins[1].stream_impact_item = item      #effluent was

    B1 = elu.EL_CMMBR('B1', ins=O1-0, outs=('effluent_MembT', 'sludge_MembT'),
                      isdynamic=True, ppl = ppl, baseline_ppl = baseline_ppl,
                      #V_max=2.9, 
                      #DO_ID='S_O2', aeration=2, suspended_growth_model=masm2d,
                      # pumped_flow=5, # after calculation # initial = 0.0001 # m3/hr
                      pumped_flow=(Q_ras*2+Q_was), # 0.026
                      **kwargs_2
                      # pumped_flow=Q_ras+Q_was,
                      # solids_capture_rate=0.999, 
                      )
    # breakpoint()
    # S2 = su.Splitter('S2', ins=B1-1, outs=['RAS', 'WAS'], split=0.95)
    S2 = su.Splitter('S2', ins=B1-1, outs=['RAS', 'WAS'], split=Q_ras*2/(Q_ras*2+Q_was))
    
    # # S2.run()
    
    S1 = su.Splitter('S1', ins=S2-0, outs=[1-PC, 1-A1], split=0.01) #0.05
    # S1 = su.Splitter('S1', ins=B1-1, outs=[1-PC, 1-A1], split=Q_was/(Q_ras+Q_was))
    
    CWT = elu.EL_CWT(
        'CWT', ins=(B1-0), outs=('effluent_CWT'),
        isdynamic=True, V_max=12, aeration=None, suspended_growth_model=None,
        ppl = ppl, baseline_ppl = baseline_ppl,
        )
    
    PV = elu.EL_WindSolar('PV', ins=CWT-0, outs='effluent',  isdynamic=True,
                          ppl = ppl, baseline_ppl = baseline_ppl,)
    
    ELH = elu.EL_Housing('ELH', ins=PV-0, outs='effluent',  isdynamic=True,
                         ppl = ppl, baseline_ppl = baseline_ppl,)
    
    # # S1.run()
    
    # CWT = elu.EL_CWT('CWT', ins=(B1-0), 
    #                 outs= ('effluent_CWT'),
    #                 isdynamic=False,
    #                 # thermo=thermo_masm2d
    #                 )
    
    # S4 = su.Splitter('S4', ins = ELH-0, outs = [1-CT, 'Reflushing'], split= 0.5)

    sys = qs.System('EL', path=(CT, PC, A1, O1, B1, S2, S1, CWT, PV, ELH))
    sys.set_dynamic_tracker(A1, O1, B1, B1-0, B1-1)
    
    
    batch_init(sys, 
               ospath.join(data_path, "units_data/bsm2p_init.xlsx"), 
               sheet='el')
    
    # sys = qs.System('EL', path=(CT, PC, S3, A1, O1, B1, S2, S1, CWT, S4),
    #                 recycle = [sludge_PC, sludge_MT_PC, sludge_MT_A1, flushing_water_CT],
    #                 ) # add flushing water
    
    GWP = qs.ImpactIndicator('GlobalWarming', alias='GWP', unit='kg CO2-eq')
    Ecosystems = qs.ImpactIndicator('H_Ecosystems', alias='Ecosystems', unit='points')
    Health = qs.ImpactIndicator('H_Health', alias='Health', unit='points')
    Resources = qs.ImpactIndicator('H_Resources', alias='Resources', unit='points')
    
    tea = qs.TEA(system=sys, discount_rate=0.05, lifetime=10, simulate_system=False)
    lca = qs.LCA(system=sys, lifetime=10, lifetime_unit='yr',
                 indicators=(GWP, Ecosystems, Health, Resources), simulate_system=False)  
    
    return sys


# %%

@time_printer
def run(t, method=None, **kwargs):
    sys = create_systemEL()    
    
    # batch_init(sys, "/Users/rishabhpuri/Desktop/bsm2p_init.xlsx", sheet='el')
    batch_init(sys, 
               ospath.join(data_path, "units_data/bsm2p_init.xlsx"), 
               sheet='el')
    
    
    # path = ospath.join(folder, "data/initial_conditions_ASM2d.xlsx")    
    # batch_init(sys, path, 
    #            sheet='el')
    # sys.set_dynamic_tracker(*sys.products)

    return sys

    
if __name__ == '__main__':
    t = 40
    # method = 'RK45'
    method = 'RK23' 
    # method = 'DOP853'
    # method = 'Radau'
    # method = 'BDF'
    # method = 'LSODA'
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    sysEL = run(t, method=method)
    
    fs = sysEL.flowsheet.stream
    fu = sysEL.flowsheet.unit

    
    sysEL.simulate(
        # state_reset_hook='reset_cache',
        t_span=(0,t),
        method=method,
        # print_t=True,
        )
    sysEL.diagram()
    act_units = [u.ID for u in sysEL.units if isinstance(u, elu.EL_Aerobic) or isinstance(u, elu.EL_Anoxic) or isinstance(u, elu.EL_CMMBR)]
    print(f'act_units = {act_units}')
    # act_units = [u.ID for u in sysEL.units if (u.ID.startswith('A1') or u.ID.startswith('O1'))]
    
    srt = get_SRT(sysEL, biomass_IDs, wastage= [fs.WAS, fs.effluent_MembT], active_unit_IDs=act_units)
    print(f'Estimated SRT assuming at steady state is {srt} days\n')
    
    # fig, axis = fs.effluent_AnoxT.scope.plot_time_series(('S_F', 'S_NO3', 'S_O2'))
    
    fig, axis = fs.effluent_MembT.scope.plot_time_series(('S_F', 'S_A', 'X_H','S_NH4', 'S_NO3', 'S_PO4', 'X_I', 'S_I', 'S_N2')) 

    
    # act_units = [u.ID for u in sys.units if isinstance(u, su.FlatBottomCircularClarifier) or u.ID.startswith('O')]
    # act_units = [u.ID for u in sys.units if u.ID.startswith('O')]
    
    # srt = get_SRT(sys, biomass_IDs, wastage= [fs.WAS, fs.effluent], active_unit_IDs=act_units)
    # print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days\n')

    qs.PowerUtility.price = 0     #    

    tea1 = qs.TEA(system=sysEL, discount_rate=0.05, lifetime=10)
    #tea1.show()
    el.get_TEA_metrics_breakdown(sysEL, include_breakdown=True)
  
    GWP = qs.ImpactIndicator('GlobalWarming', alias='GWP', unit='kg CO2-eq')
    Ecosystems = qs.ImpactIndicator('H_Ecosystems', alias='Ecosystems', unit='points/cap/yr')
    Health = qs.ImpactIndicator('H_Health', alias='Health', unit='points/cap/yr')
    Resources = qs.ImpactIndicator('H_Resources', alias='Resources', unit='points/cap/yr')
    
    lca1 = qs.LCA(system=sysEL, lifetime=10, lifetime_unit='yr',indicators=(GWP, Ecosystems, Health, Resources))
    #lca1 = qs.LCA(system=sysEL, lifetime=10, lifetime_unit='yr',)
    #lca1.get_total_impacts()
    #lca1.get_impact_table('Construction')
    lca1.show()
    #el.get_LCA_metrics_breakdown(sysEL, include_breakdown=True)





