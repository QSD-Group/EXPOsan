#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
This module is developed by:
    Siqi Tang <siqit@outlook.com>
    Yuyao Huang <yuyaoh2@illinois.edu>
    Aaron Marszewski <aaronpm3@illinois.edu>

This python file is used to perform uncertainty and sensitivity analysis for Enviroloo Clear Reinvented Toilet system.
'''

# %% 
import qsdsan as qs
from qsdsan import (
    Flowsheet, main_flowsheet,
    WasteStream,
    sanunits as su,
    ImpactItem, 
    LCA, TEA, System,
    )
from qsdsan.sanunits import Trucking
from chaospy import distributions as shape
from qsdsan.utils import clear_lca_registries
# from qsdsan.sanunits._excretion import ExcretionmASM2d

# from qsdsan.utils import load_components, set_thermo


import qsdsan as qs
import os

from qsdsan import (
    Flowsheet, main_flowsheet,
    WasteStream,
    sanunits as su,
    ImpactItem, 
    LCA, TEA, System,
    )

from qsdsan.utils import (
    ospath, 
    time_printer, 
    load_data, 
    get_SRT,
    )

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho
from exposan.bwaise import create_components as create_bw_components
import numpy as np
from qsdsan import processes as pc, sanunits as su, Model as mod
from chaospy import distributions as shape


from exposan.utils import add_fugitive_items
from exposan.enviroloo import _units as elu
from exposan.enviroloo import (
    data_path,
    _load_components,
    _load_lca_data,
    discount_rate,
    get_decay_k,
    get_toilet_users,
    max_CH4_emission,
    operator_daily_wage,
    ppl,
    #get_tanker_truck_fee
    update_resource_recovery_settings,
    )
from exposan.enviroloo._EL_pumps import (
    LiftPump, AgitationPump, DosingPump, ReturnPump, SelfPrimingPump, 
    AirDissolvedPump, MicroBubblePump, ClearWaterPump,)


folder = ospath.dirname(__file__)

__all__ = ('create_systemEL', 'create_system',)
#%%

'''
Name notes:

To make programming more convenient, we use the following names for the units in the system:

Toilet: Toilet
Collection tank: CT
Primary clarifier: PC
Anoxic tank: AnoT
Aerobic tank: AerT
Membrane tank: MemT
Clear water tank: CWT
Pressure tank: PT
Photovoltaic and Wind system: photovoltaic_wind
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
Q_w = 15 # m3/hr  # 250 l/min flow for lift pump
Q_ras = 280 # m3/day # 200 l/min for nitrate pump
Q_was = 0.15*24
# Q_was = 280 # m3/day # 200 l/min for sludge return pump

toilet_waste={
    'S_NH4':  13.7,
    # 'S_NO3': 0.2,     # Nitrate
    'S_NO3': 1.6,
    # 'S_PO4': 2.3,
    'S_PO4': 3.5,
    'S_O2': 3.0,      # Dissolved oxygen
    'S_F': 19.1,      # Fermentable COD
    'S_I': 12.8,      # Inert soluble COD
    'X_S': 29.3,        # Slowly biodegradable particulates = TSS
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
    # Tissue = Component('Tissue', MW=1, phase='s', particle_size='Particulate',
    #                     degradability='Undegradable', organic=False,
    #                     description='Tissue for toilet paper')
    # # 375 kg/m3 is the average of 250-500 for tissue from
    # # https://paperonweb.com/density.htm (accessed 2020-11-12)
    # add_V_from_rho(Tissue, 375)

    # WoodAsh = Component('WoodAsh', MW=1, phase='s', i_Mg=0.0224, i_Ca=0.3034,
    #                     particle_size='Particulate', degradability='Undegradable',
    #                     organic=False, description='Wood ash for desiccant')
    # add_V_from_rho(WoodAsh, 760)

    # for i in (Tissue, WoodAsh):
    #     i.copy_models_from(Chemical('Glucose'), ('Cn', 'mu'))
    
    # H2O = Component('H2O', phase='l', particle_size='Soluble',
    #                 degradability='Undegradable', organic=False)

    # PAC = Component('PAC', search_ID='10124-27-3', phase='s', particle_size='Particulate', degradability='Slowly', organic=False)
    # add_V_from_rho(PAC, rho=2800)
                    
    # # Glucose = Component('Glucose', search_ID='50-99-7', phase='s', particle_size='Particulate', degradability='Readily', organic=False)
    # # add_V_from_rho(Glucose, rho=1560)


    cmps = Components([*masm2d_cmps, 
                       # Tissue, WoodAsh, H2O,
                       ])
    #cmps.compile()
    cmps.compile(ignore_inaccurate_molar_weight=True)
    if set_thermo: qs.set_thermo(cmps)
    return cmps
# %%
def batch_create_streams(prefix, phases=('liq', 'sol')):
    item = ImpactItem.get_item('CH4_item').copy(f'{prefix}_CH4_item', set_as_source=True)
    WasteStream('CH4', phase='g', stream_impact_item=item)

    item = ImpactItem.get_item('N2O_item').copy(f'{prefix}_N2O_item', set_as_source=True)
    WasteStream('N2O', phase='g', stream_impact_item=item)

    price_dct = update_resource_recovery_settings()[0]
    for nutrient in ('N', 'P', 'K'):
        for phase in phases:
            original = ImpactItem.get_item(f'{nutrient}_item')
            new = original.copy(f'{phase}_{nutrient}_item', set_as_source=True)
            WasteStream(f'{phase}_{nutrient}', phase='l',
                        price=price_dct[nutrient], stream_impact_item=new)

    def create_stream_with_impact_item(stream_ID='', item_ID='', dct_key=''):
        item_ID = item_ID or stream_ID+'_item'
        dct_key = dct_key or stream_ID
        item = ImpactItem.get_item(item_ID).copy(f'{prefix}_{item_ID}', set_as_source=True)
        WasteStream(f'{stream_ID}', phase='s',
                    price=price_dct.get(dct_key) or 0., stream_impact_item=item)

    # create_stream_with_impact_item(stream_ID='ammonium')
    # create_stream_with_impact_item(stream_ID='struvite')
    # create_stream_with_impact_item(stream_ID='NaOH')
    # create_stream_with_impact_item(stream_ID='NaClO')
    # create_stream_with_impact_item(stream_ID='O3')
    # create_stream_with_impact_item(stream_ID='PAC')
    # create_stream_with_impact_item(stream_ID='Glucose')
    # create_stream_with_impact_item(stream_ID='air')
    WasteStream('Glucose_Dose', S_F= 0.9, units='kg/hr', T=Temp) # 0.0805
    WasteStream('PAC_Dose', X_AlOH= 0.12, units='kg/hr', T=Temp) # 0.1207

#%%
    
def create_systemEL(flowsheet=None, inf_kwargs={}, masm_kwargs={}, init_conds={},
                  aeration_processes=()):
    # Components and stream
    flowsheet = flowsheet or main_flowsheet
    streamEL = flowsheet.stream
    batch_create_streams('EL')
    
    cmps = create_components()
    toilet_ins = qs.WasteStream('toilet_waste', T=Temp)
    toilet_ins.set_flow_by_concentration(Q_w*0.9, 
                                         concentrations=toilet_waste, 
                                         units=('m3/hr', 'mg/L'))
    
    # Glucose = qs.WasteStream('Glucose_Dose', S_F= 0.9, units='kg/hr', T=Temp) # 0.0805
    # PAC = qs.WasteStream('PAC_Dose', X_AlOH= 0.12, units='kg/hr', T=Temp) # 0.1207
    
    masm2d = pc.mASM2d(**masm_kwargs)
    
    kwargs_O = dict(V_max=6.35, aeration=2, DO_ID='S_O2', suspended_growth_model=masm2d)
    kwargs_1 = dict(V_max=6.35, aeration=None, DO_ID=None, suspended_growth_model=masm2d)
    # kwargs_2 = dict(V_max=2.89, aeration=2, DO_ID='S_O2', suspended_growth_model=masm2d)
    

    # CT = su.Mixer('CT', 
    #               # ins=(toilet_ins, 'sludge_PC', 'flushing_water_CT', 'spill_PC'), 
    #               ins=(toilet_ins, 'sludge_PC',), 
    #               outs=('effluent_CT'),
    #               )
    CT = elu.EL_CT(
        'CT', ins=(toilet_ins, 'sludge_PC','flushing_water_CT'), outs=('effluent_CT'),
        isdynamic=True, V_max=10, aeration=None, suspended_growth_model=None
        )
    
    # PC = su.PrimaryClarifier('PC', ins=[CT-0, 'sludge_AeroT'],
    #                          outs=('effluent_PC_total', sludge_PC), isdynamic=True,
    #                          init_with='WasteStream', thermo=thermo_masm2d,
    #                          sludge_flow_rate=280, 
    #                          solids_removal_efficiency=0.6)
    
    PC = elu.EL_PC('PC', ins=(CT-0, 'RAS_PC'), outs=('effluent_PC_total', 1-CT),
                   ppl=ppl, baseline_ppl=100,
                   solids_removal_efficiency=0.85,
                   isdynamic=True,
                   # thermo=thermo_masm2d,
                   sludge_flow_rate=Q_w*0.1*24)  # 0.00116 m3/hr
    

    # S3 = su.Splitter('S3', ins = PC-0, outs = ['effluent_PC', 3-CT], split= 0.8)
    


    
    # A1 = su.CSTR('A1', ins=(S3-0), 
    #                         outs = ('effluent_AnoxT'), isdynamic=True, **kwargs_1
    #                         # W_tank= 2.09,
    #                         # # ppl = ppl, baseline_ppl = 100,
    #                         # aeration=None, DO_ID='S_O2', suspended_growth_model=asm2d, 
    #                         # V_max= 7.33, 
    #                         )
    
    A1 = elu.EL_Anoxic('A1', ins=(PC-0, 'RAS_A1', streamEL['Glucose_Dose']), outs=('effluent_AnoxT',),
                       isdynamic=True, **kwargs_1
                       # W_tank= 2.09,
                        # # ppl = ppl, baseline_ppl = 100,
                        # aeration=None, DO_ID='S_O2', suspended_growth_model=asm2d, 
                        # V_max= 7.33, 
                        )

    
    # O1 = su.CSTR('O1', ins=(A1-0), 
    #                        outs = ('effluent_AeroT', 1-PC),isdynamic=True, **kwargs_O
    #                        # aeration = 2, suspended_growth_model=asm2d,
    #                        #  # ppl = ppl, baseline_ppl = 100,
    #                        #  W_tank= 2.09,
    #                        #  V_max=7.33,**kwargs_0, 
    #                         )
    
    O1 = elu.EL_Aerobic('O1', ins=(A1-0, streamEL['PAC_Dose']), outs=('effluent_AeroT',),
                        isdynamic=True, **kwargs_O
                        # aeration = 2, suspended_growth_model=asm2d,
                        #  # ppl = ppl, baseline_ppl = 100,
                        #  W_tank= 2.09,
                        #  V_max=7.33,**kwargs_0, 
                         )


    B1 = elu.EL_CMMBR('B1', ins=O1-0, outs=('effluent_MembT', 'sludge_MembT'),
                      isdynamic=True, V_max=2.9, 
                      DO_ID='S_O2', aeration=2, suspended_growth_model=masm2d,
                      # pumped_flow=5, # after calculation # initial = 0.0001 # m3/hr
                      pumped_flow=(Q_ras*2+Q_was),
                      # pumped_flow=Q_ras+Q_was,
                      # solids_capture_rate=0.999, 
                      )
    # breakpoint()
    S2 = su.Splitter('S2', ins=B1-1, outs=['RAS', 'WAS'], split=Q_ras*2/(Q_ras*2+Q_was))
    
    # # S2.run()
    
    S1 = su.Splitter('S1', ins=S2-0, outs=[1-PC, 1-A1], split=0.5)
    # S1 = su.Splitter('S1', ins=B1-1, outs=[1-PC, 1-A1], split=Q_was/(Q_ras+Q_was))
    
    CWT = elu.EL_CWT(
        'CWT', ins=(B1-0), outs=('effluent_CWT'),
        isdynamic=True, V_max=12, aeration=None, suspended_growth_model=None
        )
    
    # # S1.run()
    
    # CWT = elu.EL_CWT('CWT', ins=(B1-0), 
    #                 outs= ('effluent_CWT'),
    #                 isdynamic=False,
    #                 # thermo=thermo_masm2d
    #                 )
    PV = elu.EL_WindSolar()
    
    S4 = su.Splitter('S4', ins = CWT-0, outs = [2-CT, 'Reflushing'], split= 0.5)
    
    
    

    sysEL = qs.System('EL', path=(CT, PC, A1, O1, B1, S2, S1, CWT, S4,))
    sysEL.set_dynamic_tracker(A1, O1, B1, B1-0, B1-1)
    # sysEL.simulate()
    sysEL.simulate(
        # state_reset_hook='reset_cache',
        t_span=(0,2),
        method='RK23',
        print_t=True,
        )
    
    teaEL = TEA(system=sysEL, discount_rate=discount_rate,
       start_year=2020, lifetime=20, uptime_ratio=1,
       # CEPCI = 567.5,
       # CAPEX = 2.00,  
       #lang_factor=None,
       lang_factor=None,
       annual_maintenance=0,
       # annual_labor=(operator_daily_wage*3*365),
       annual_labor=0
       )
    get_powerEL = lambda: sum([u.power_utility.rate for u in sysEL.units]) * (24 * 365 * teaEL.lifetime)
    LCA(system=sysEL, lifetime=20, lifetime_unit='yr', uptime_ratio=1.0, e_item=get_powerEL)
    
    # sys = qs.System('EL', path=(CT, PC, S3, A1, O1, B1, S2, S1, CWT, S4),
    #                 recycle = [sludge_PC, sludge_MT_PC, sludge_MT_A1, flushing_water_CT],
    #                 ) # add flushing water
    
    

    return sysEL

def create_system(system_ID='EL', flowsheet=None, 
                  #adjust_MW_to_measured_as=False
                  ):
    ID = system_ID.lower().lstrip('sys').upper()
    reload_lca = False

    #set flowsheet to avoid stream replacement warnings
    if flowsheet is None:
        flowsheet_ID = f'el{ID}'
        if hasattr(main_flowsheet.flowsheet, flowsheet_ID): # clear flowsheet
            getattr(main_flowsheet.flowsheet, flowsheet_ID).clear()
            clear_lca_registries()
            reload_lca = True
        flowsheet = Flowsheet(flowsheet_ID)
        main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    _load_lca_data(reload_lca)

    if system_ID == 'EL': f = create_systemEL
    elif system_ID == 'E': f = create_systemEL
    elif system_ID == 'L': f = create_systemEL
    else: raise ValueError(f'`system_ID` can only be "EL", "E", or "L", not "{ID}".')
    
    try: system = f(flowsheet)
    except:
        _load_components(reload=True)
        system = f(flowsheet)
    
    return system

# %%

# @time_printer
# def run(t, method=None, **kwargs):
#     sys = create_systemEL()    
    
#     # batch_init(sys, "/Users/rishabhpuri/Desktop/bsm2p_init.xlsx", sheet='el')
#     batch_init(sys, 
#                ospath.join(data_path, "units_data/bsm2p_init.xlsx"), 
#                sheet='el')
    
    
#     # path = ospath.join(folder, "data/initial_conditions_ASM2d.xlsx")    
#     # batch_init(sys, path, 
#     #            sheet='el')
#     # sys.set_dynamic_tracker(*sys.products)
    

#     return sys

    
# if __name__ == '__main__':
#     t = 2
#     # method = 'RK45'
#     method = 'RK23' 
#     # method = 'DOP853'
#     # method = 'Radau'
#     # method = 'BDF'
#     # method = 'LSODA'
#     msg = f'Method {method}'
#     print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
#     print(f'Time span 0-{t}d \n')
#     sys = run(t, method=method)
    
#     sys.diagram()
#     fs = sys.flowsheet.stream
#     fu = sys.flowsheet.unit
    
#     sys.simulate(
#         # state_reset_hook='reset_cache',
#         t_span=(0,t),
#         method=method,
#         print_t=True,
#         )
    
#     sys.diagram()
    
    '''teaEL = TEA(system=sys, discount_rate=discount_rate,
           start_year=2020, lifetime=20, uptime_ratio=1,
           # CEPCI = 567.5,
           # CAPEX = 2.00,  
           #lang_factor=None,
           lang_factor=None,
           annual_maintenance=0,
           # annual_labor=(operator_daily_wage*3*365),
           annual_labor=0
           )
    get_powerEL = lambda: sum([u.power_utility.rate for u in fu]) * (24 * 365 * teaEL.lifetime)
    LCA(system=sys, lifetime=20, lifetime_unit='yr', uptime_ratio=1.0, e_item=get_powerEL)
    
    
    # act_units = [u.ID for u in sys.units if isinstance(u, su.FlatBottomCircularClarifier) or u.ID.startswith('O')]
    # act_units = [u.ID for u in sys.units if u.ID.startswith('O')]
    
    # srt = get_SRT(sys, biomass_IDs, wastage= [fs.WAS, fs.effluent], active_unit_IDs=act_units)
    # print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days\n')
'''




