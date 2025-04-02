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
from qsdsan.sanunits._excretion import ExcretionmASM2d

# from qsdsan.utils import load_components, set_thermo


import qsdsan as qs
import os

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

__all__ = ('create_systemEL',)

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
# Q_urine = 3800
# Q_feces = 3800
# urine_waste={'S_NH4':  0.241,
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
                
                
    # c1s = {k:v for k,v in dct['C1_s'].items() if v>0}
    # c1x = {k:v for k,v in dct['C1_x'].items() if v>0}
    # tss = [v for v in dct['C1_tss'].values() if v>0]
    # u.C1.set_init_solubles(**c1s)
    # u.C1.set_init_sludge_solids(**c1x)
    # u.C1.set_init_TSS(tss)

# %% Create EnviroLoo Clear system

def create_components(set_thermo = True
                      #adjust_MW_to_measured_as=False
                      ):
    # bw_cmps = create_bw_components(set_thermo=False)
    masm2d_cmps = pc.create_masm2d_cmps(set_thermo=True)
    Tissue = Component('Tissue', MW=1, phase='s', particle_size='Particulate',
                        degradability='Undegradable', organic=False,
                        description='Tissue for toilet paper')
    # 375 kg/m3 is the average of 250-500 for tissue from
    # https://paperonweb.com/density.htm (accessed 2020-11-12)
    add_V_from_rho(Tissue, 375)

    WoodAsh = Component('WoodAsh', MW=1, phase='s', i_Mg=0.0224, i_Ca=0.3034,
                        particle_size='Particulate', degradability='Undegradable',
                        organic=False, description='Wood ash for desiccant')
    add_V_from_rho(WoodAsh, 760)

    for i in (Tissue, WoodAsh):
        i.copy_models_from(Chemical('Glucose'), ('Cn', 'mu'))
    
    H2O = Component('H2O', phase='l', particle_size='Soluble',
                    degradability='Undegradable', organic=False)

    # PAC = Component('PAC', search_ID='10124-27-3', phase='s', particle_size='Particulate', degradability='Slowly', organic=False)
    # add_V_from_rho(PAC, rho=2800)
                    
    # # Glucose = Component('Glucose', search_ID='50-99-7', phase='s', particle_size='Particulate', degradability='Readily', organic=False)
    # # add_V_from_rho(Glucose, rho=1560)


    cmps = Components((*masm2d_cmps, Tissue, WoodAsh, H2O,
                       ))
    cmps.compile(ignore_inaccurate_molar_weight=True)
    return cmps
#%%
    
def create_systemEL(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={},
                  aeration_processes=()):
    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    thermo_masm2d = qs.get_thermo()
    masm2d = pc.mASM2d()
    
    # urine = qs.WasteStream('urine', T=Temp)
    # feces = qs.WasteStream('feces', T=Temp)
    
    # urine.set_flow_by_concentration(Q_urine, 
    #                                   concentrations=urine_waste, 
    #                                   units=('m3/d', 'mg/L'))
    # feces.set_flow_by_concentration(Q_feces, 
    #                                   concentrations=feces_waste, 
    #                                   units=('m3/d', 'mg/L'))
    
    # effluent = qs.WasteStream('effluent', T=Temp)
    Recycle = qs.WasteStream('Recycle', T=Temp)
    sludge_PC = qs.WasteStream('sludge_PC', T=Temp)
    influent_PC = qs.WasteStream('influent_PC', T=Temp)
    sludge_MT = qs.WasteStream('sludge_MT', T=Temp)
    sludge_MT_PC = qs.WasteStream('sludge_MT_PC', T=Temp)
    sludge_MT_A1 = qs.WasteStream('sludge_MT_A1', T=Temp)
    flushing_water = qs.WasteStream('flushing_water', T=Temp)
    # RAS = qs.WasteStream('RAS', T=Temp)
    
    # eff_GT = qs.WasteStream(ID = 'effluent_GT')
    # eff_MT = qs.WasteStream(ID = 'effluent_MT')
    # eff_DU = qs.WasteStream(ID = 'effluent_DU')
    
    
    WasteWater = elu.EL_Excretion('WasteWater', outs=('urine','feces'))
    
    # WasteWater._run()
    
    
    '''
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

    Toilet = elu.EL_Toilet('Toilet',
                    ins=(WasteWater-0, WasteWater-1, 'toilet_paper', 'flushing_water','cleansing_water', 'desiccant'), # add flushing water
                    outs=('mixed_waste'),
                    N_user=100, N_tot_user=1,
                    # F_BM_default=1
                    # lifetime=10, if_include_front_end=True,
                    # if_toilet_paper=True, 
                    # if_flushing=True, 
                    # if_cleansing=False,
                    # if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,                 
                    # CAPEX=500*max(1, ppl/100), OPEX_over_CAPEX=0.06
                    )
    # M1 = su.Mixer('M1', ins=[WasteWater-0, WasteWater-1], outs='mixed_waste') # mixer
    # Toilet._run()
    
    
    # breakpoint()

    CT = elu.EL_CT('CT', ins=(Toilet-0, sludge_PC, Recycle), 
                    outs = ('effluent_CT'),
                    )
    
    # CT.run()
    
    PC = elu.EL_PC('PC', ins=(CT-0, sludge_MT_PC), outs=('effluent_PC', sludge_PC),
               ppl=ppl, baseline_ppl=100,
               solids_removal_efficiency=0.85,
               sludge_flow_rate=0.00116)
    
    # PC.run()
    # breakpoint()

    
    A1 = elu.EL_Anoxic('A1', ins=(PC-0, sludge_MT_A1), 
                            outs = ('effluent_AnoxT'),
                            # ppl = ppl, baseline_ppl = 100,
                            aeration=None, DO_ID='S_O2', suspended_growth_model=masm2d, 
                            V_max= 7.33
                            )

    # A1.run()
    # breakpoint()
    
    O1 = elu.EL_Aerobic('O1', ins=(A1-0), 
                           outs = ('effluent_AeroT'), 
                           aeration = 2, suspended_growth_model=masm2d,
                            # ppl = ppl, baseline_ppl = 100,
                            V_max=7.33
                            )

    # O1.run()
    # breakpoint()
    B1 = elu.EL_CMMBR('B1', ins=(O1-0), 
                        outs = (sludge_MT, 'effluent_MembT'),aeration=2, suspended_growth_model=masm2d,
                        pumped_flow=0.0001, 
                        solids_capture_rate=0.999, 
                        V_max=3.34, 
                        )
                        
    # B1 = su.CompletelyMixedMBR('B1', ins=(O1-0), 
    #                     outs = (sludge_MT, 'effluent_MembT'),DO_ID='S_O2', aeration=2, suspended_growth_model=masm2d,
    #                     pumped_flow=0.0001, 
    #                     solids_capture_rate=0.999, 
    #                     V_max=3.34,
                        
    #                     # ppl = ppl,
    #                     # baseline_ppl = 100,
    #                     )
    # B1.run()
    # breakpoint()
    
    S1 = su.Splitter('S1', ins = B1-0, outs = [sludge_MT_PC, sludge_MT_A1], split= 0.5)
    
    # S1.run()


    CWT = elu.EL_CWT('CWT', ins=(B1-1), 
                    outs= ('ClearWater'), 
                    # V_wf = 0.9, max_oveflow=15, 
                    # ppl = ppl, baseline_ppl = 100,
                    )
    
    # CWT.run()
    
    S2 = su.Splitter('S2', ins = CWT-0, outs = ['effluent', Recycle], split= 0.5)
   
    # S2.run()
    
    PT = su.Mixer('PT', ins=S2-0,)
   
    # PT = elu.EL_PT('PT', ins=S2-0, outs=(flushing_water), vessel_material = None, V_wf = None, # add flushing water
    #                     include_construction = True, length_to_diameter = None, 
    #                     F_BM_default = 1, kW_per_m3 = 0.1, vessel_type = None, tau = None, 
    #                     ppl = ppl, baseline_ppl = 100,
    #               )
    # PT.run(
    
    # breakpoint()
    
    sys = qs.System('EL', path=(PC, A1, O1, B1, S1, CWT, S2, PT, WasteWater, Toilet, CT))
                #recycle = [Recycle, sludge_PC, sludge_MT_PC, sludge_MT_A1, flushing_water]) # add flushing water

    # sys = qs.System('G1_WERF', path=(PC, S1, A1, A2, A3, A4, O1, O2, C1, GT, MT, M1, J1, DG, J2, DU, M2), 
    #                 recycle = [M2-0, 1-A3, RAS])
                    # recycle = [eff_GT, eff_MT, eff_DU, RAS])
    # sys.set_tolerance(rmol=1e-6)
    # sys.maxiter = 500
    return sys


    # sysEL_PCrecycle = System('sysEL_PCspill',
    #                  path = (WasteWaterGenerator, Toilet, CT, P_CT_lift, PC, P_PC_return),
    #                  recycle = PC-2
    #                  )
    
    # sysEL_CWTrecycle = System('sysEL_CWTrecycle',
    #                    path = (sysEL_PCrecycle, P_AnoxT_agitation, AnoxT, AeroT,MembT, P_NitrateReturn_PC, 
    #                            P_NitrateReturn_AnoxT, P_MT_selfpriming, CWT), 
    #                    recycle = CWT-1
    #                    )
    
    # sysEL = System('sysEL',
    #                          path = (sysEL_CWTrecycle, P_CWT, PT, 
    #                                  # Total_N2O, 
    #                                  # Total_CH4, 
    #                                  Pipeline_system),
    #                          recycle = PT-0
    #                          )
    
    
# def create_system(system_ID='EL', flowsheet=None, 
#                   #adjust_MW_to_measured_as=False
#                   ):
#     ID = system_ID.lower().lstrip('sys').upper()
#     reload_lca = False

#     #set flowsheet to avoid stream replacement warnings
#     if flowsheet is None:
#         flowsheet_ID = f'el{ID}'
#         if hasattr(main_flowsheet.flowsheet, flowsheet_ID): # clear flowsheet
#             getattr(main_flowsheet.flowsheet, flowsheet_ID).clear()
#             clear_lca_registries()
#             reload_lca = True
#         flowsheet = Flowsheet(flowsheet_ID)
#         main_flowsheet.set_flowsheet(flowsheet)
    
#     _load_components()
#     _load_lca_data(reload_lca)

#     if system_ID == 'EL': f = create_systemEL
#     elif system_ID == 'E': f = create_systemEL
#     elif system_ID == 'L': f = create_systemEL
#     else: raise ValueError(f'`system_ID` can only be "EL", "E", or "L", not "{ID}".')
    
#     try: system = f(flowsheet)
#     except:
#         _load_components(reload=True)
#         system = f(flowsheet)
    
#     return system

# %%

@time_printer
def run(t, method=None, **kwargs):
    sys = create_systemEL()    
    
    batch_init(sys, "/Users/rishabhpuri/Desktop/bsm2p_init.xlsx", sheet='el')
    
    
    # path = ospath.join(folder, "data/initial_conditions_ASM2d.xlsx")    
    # batch_init(sys, path, 
    #            sheet='el')
    # sys.set_dynamic_tracker(*sys.products)
    
    
    
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
    
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        method=method,
        # print_t=True,
        )

    
    # act_units = [u.ID for u in sys.units if isinstance(u, su.FlatBottomCircularClarifier) or u.ID.startswith('O')]
    # act_units = [u.ID for u in sys.units if u.ID.startswith('O')]
    
    # srt = get_SRT(sys, biomass_IDs, wastage= [fs.WAS, fs.effluent], active_unit_IDs=act_units)
    # print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days\n')





