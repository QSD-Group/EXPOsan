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
#%% Create Universal Units and Functions
def batch_init(sys, path, sheet):
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    for k in sys.units:
        print(f'k={k}')
        if k.ID.startswith(('O', 'A', 'M')):
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
    cmps.compile()    
    return cmps
#%%
    
def create_systemEL(flowsheet=None, inf_kwargs={}, asm_kwargs={}, init_conds={},
                  aeration_processes=()):
    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    thermo_masm2d = qs.get_thermo()
    masm2d = pc.mASM2d()
    
    effluent = qs.WasteStream('effluent', T=Temp)
    Recycle = qs.WasteStream('Recycle', T=Temp)
    sludge_PC = qs.WasteStream('sludge_PC', T=Temp)
    return_PC = qs.WasteStream('return_PC', T=Temp)
    sludge_MT = qs.WasteStream('sludge_MT', T=Temp)
    # RAS = qs.WasteStream('RAS', T=Temp)
    
    # eff_GT = qs.WasteStream(ID = 'effluent_GT')
    # eff_MT = qs.WasteStream(ID = 'effluent_MT')
    # eff_DU = qs.WasteStream(ID = 'effluent_DU')
    
    WasteWaterGenerator = elu.EL_Excretion('WasteWaterGenerator', outs=('urine','feces'))
    WasteWaterGenerator._run()
    # breakpoint()

    Toilet = elu.EL_Toilet('Toilet',
                    ins=(WasteWaterGenerator-0, WasteWaterGenerator-1, 'toilet_paper', 'flushing_water', 'cleansing_water', 'desiccant'),
                    outs=('mixed_waste'),
                    N_user=get_toilet_users(), N_tot_user=ppl,
                    # lifetime=10, if_include_front_end=True,
                    if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                    if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,                 
                    # CAPEX=500*max(1, ppl/100), OPEX_over_CAPEX=0.06
                    )
    Toilet._run()
    

    CT = elu.EL_CT('CT', ins=(Toilet-0, sludge_PC,return_PC, Recycle), 
                    outs = ('effluent_CT'),
                    )
    
    # CT.run()
    
    
    PC = elu.EL_PC('PC', ins=(CT-0, sludge_MT), outs=('effluent_PC', 1-CT , 2-CT),
                    ppl = ppl,  # The number of people served
                    baseline_ppl = 100,
                    solids_removal_efficiency = 0.85,  # The solids removal efficiency
                    sludge_flow_rate = .00116,  # Sludge flow rate
                    # max_oveflow = 15,
                    )
    
    # PC.run()

    
    A1 = elu.EL_Anoxic('A1', ins=(PC-0, 'sludge_MT'), 
                            outs = ('effluent_AnoxT'),
                            # ppl = ppl, baseline_ppl = 100,
                            aeration=None, DO_ID='S_O2', suspended_growth_model=masm2d, 
                            )

    # A1.run()
    
    O1 = elu.EL_Aerobic('O1', ins=(A1-0), 
                           outs = ('effluent_AeroT'), 
                            aeration = 2, suspended_growth_model=masm2d
                            # ppl = ppl, baseline_ppl = 100,
                            )

    # AeroT.run()
    # breakpoint()
    M1 = elu.EL_CMMBR('M1', ins=(O1-0), 
                        outs = ('effluent_MembT', 1-PC, 1-A1),
                        ppl = ppl,
                        baseline_ppl = 100,
                        )
    # M1.run()

    # CWT = elu.EL_CWT('CWT', ins=(M1-0), 
    #                 outs= ('ClearWater'), 
    #                 # V_wf = 0.9, max_oveflow=15, 
    #                 # ppl = ppl, baseline_ppl = 100,
    #                 )
    
    S1 = su.Splitter('S1', ins = M1-0, outs = ['effluent', 3-CT], split= 0.5)
   
   
    PT = elu.EL_PT('PT', ins=S1-0, outs=(3-Toilet), vessel_material = None, V_wf = None, 
                        include_construction = True, length_to_diameter = None, 
                        F_BM_default = 1, kW_per_m3 = 0.1, vessel_type = None, tau = None, 
                        ppl = ppl, baseline_ppl = 100,
                        )
    
    sys = qs.System('EL', path=(WasteWaterGenerator, Toilet, CT, PC, A1, O1, M1, S1, PT), 
                recycle = [Recycle, sludge_PC, return_PC, sludge_MT])

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
    sys.set_dynamic_tracker(*sys.products)
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        method=method,
        # print_t=True,
        **kwargs)
    
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
    
    # act_units = [u.ID for u in sys.units if isinstance(u, su.FlatBottomCircularClarifier) or u.ID.startswith('O')]
    act_units = [u.ID for u in sys.units if u.ID.startswith('O')]
    
    # srt = get_SRT(sys, biomass_IDs, wastage= [fs.WAS, fs.effluent], active_unit_IDs=act_units)
    # print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days\n')






