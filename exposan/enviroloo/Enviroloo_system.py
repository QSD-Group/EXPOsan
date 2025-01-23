#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Created by Yuyao Huang and Siqi Tang

This python file is used to perform uncertainty and sensitivity analysis for Enviroloo Clear Reinvented Toilet system.
'''

# %% 
import qsdsan as qs
from qsdsan import (
    Component, Components, get_thermo,
    Flowsheet, main_flowsheet,
    WasteStream, SanStream,
    sanunits as su,
    set_thermo as qs_set_thermo,
    ImpactItem, LCA, TEA,
    System,
    Model,
    )
from qsdsan.sanunits import Trucking

from chaospy import distributions as shape

from qsdsan.utils import clear_lca_registries
from exposan.utils import add_fugitive_items, get_generic_tanker_truck_fee as get_tanker_truck_fee
from exposan.enviroloo import (
    _load_components,
    _load_lca_data,
    discount_rate,
    get_decay_k,
    get_toilet_users,
    max_CH4_emission,
    operator_daily_wage,
    ppl,
    price_dct,
    update_resource_recovery_settings,
    )
from exposan.enviroloo._EL_pumps import (
    LiftPump, AgitationPump, DosingPump, 
    ReturnPump, SelfPrimingPump, 
    AirDissolvedPump, 
    MicroBubblePump, ClearWaterPump,)

__all__ = ('create_system',)

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
Anoxic mixing pump: P_AnoT_mixing
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

# %% Create Universal Units and Functions
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

    create_stream_with_impact_item(stream_ID='Glucose')
    create_stream_with_impact_item(stream_ID='PAC')
    create_stream_with_impact_item(stream_ID='O3')



# %% Create EnviroLoo Clear system
def create_system(flowsheet = None): # figure out the "components" and "flowsheet" arguments in Biogenic_refinery system (in Exposan repository)
    
    reload_lca = False;
    
    _load_components()
    _load_lca_data(reload_lca)

    flowsheet = flowsheet or main_flowsheet
    stream = flowsheet.stream
    batch_create_streams('EL')
    
    WasteWaterGenerator = su.Excretion('WasteWaterGenerator', outs=('urine', 'feces'))
    
    Toilet = su.MURT('Toilet', ins=(WasteWaterGenerator-0, WasteWaterGenerator-1, 'FlushingWater', 'ToiletPaper'), 
                    outs =('MixedWasteWater', 'Toilet_CH4', 'Toilet_N2O'),
                    decay_k_COD = get_decay_k(),
                    decay_k_N = get_decay_k(),
                    max_CH4_emission = max_CH4_emission,
                    N_user = 100, # two scenarios, i.e., 100 users served per day at household scale and 1000 users served per day at institutional scale
                    N_toilet = 3, # will be included for uncertainty analysis
                    if_flushing = True, if_desiccant = False, if_toilet_paper = True,
                    CAPEX = 0, # capital cost of a single toilet
                    OPEX_over_CAPEX = 0.07, # fraction of annual operating cost over total capital cost
                    )
    
    CT = su.EL_CT('CT', ins=(Toilet-0, 'ClearWaterTank_spill','PrimaryClar_spill', 'PrimaryClarP_return'), 
                  outs = ('TreatedWater', 'CT_CH4', 'CT_N20'),
                  V_wf=0.9, ppl=1000, baseline_ppl=30,
                  kW_per_m3=0.1,  # The power consumption per unit volume of the tank
                  )
    
    P_CT_lift = LiftPump('P_CT_lift', ins=CT-0, outs = 'TreatedWater',
                          working_factor=0.9,  # The ratio of the actual output and the design output
                          operation_time=12,  # Total run time of system or plant [h/d]
                          pump_lifetime=5,  # Lifetime of the pump [years]
                          pump_cost=200,
                          dP_design=0,
                          )
    
    PC = su.EL_PC('PC', ins=(LiftPump-0, 'NitrateReturn_MT'), outs=('TreatedWater', 2-CT , 3-CT, 'PC_CH4', 'PC_N2O'),
                 ppl = 1000,  # The number of people served
                 baseline_ppl = 30,
                 solids_removal_efficiency = 0.85,  # The solids removal efficiency
                 sludge_flow_rate=0.5,  # Sludge flow rate
                 max_oveflow=0.3,
                )

    P_PC_return = ReturnPump('P_PC_return', ins=PC-0, outs = 3-CT, 
                             #dP_design=405300,
                             working_factor=0.9,  # The ratio of the actual output and the design output
                             operation_time=12,  # Total run time of system or plant [h/d]
                             pump_lifetime=5,  # Lifetime of the pump [years]
                             pump_cost=200,
                             dP_design=0,
                             )

    P_Glu_agitation = AgitationPump('P_Glu_agitation', ins='Glucose', outs='GlucoseAgitation', 
                                    #dP_design=405300,
                                    working_factor=0.9,  # The ratio of the actual output and the design output
                                    operation_time=12,  # Total run time of system or plant [h/d]
                                    pump_lifetime=5,  # Lifetime of the pump [years]
                                    pump_cost=200,
                                    dP_design=0,
                                    )
    
    P_Glu_dosing = DosingPump('P_Glu_dosing', ins=P_Glu_agitation-0, outs='GlucoseDosing',
                              working_factor=0.9,  # The ratio of the actual output and the design output
                              operation_time=12,  # Total run time of system or plant [h/d]
                              pump_lifetime=5,  # Lifetime of the pump [years]
                              pump_cost=200,
                              dP_design=0,
                              ) 
    
    P_AnoxT_agitation = AgitationPump('P_AnoxT_agitation', ins=(PC-0, P_Glu_dosing-0), outs='AgitationWater', 
                                      #dP_design=405300,
                                      working_factor=0.9,  # The ratio of the actual output and the design output
                                      operation_time=12,  # Total run time of system or plant [h/d]
                                      pump_lifetime=5,  # Lifetime of the pump [years]
                                      pump_cost=200,
                                      dP_design=0,                                    
                                      )
    
    AnoxT = su.EL_Anoxic('AnoxT', ins=(P_AnoxT_agitation-0, 'NitrateReturn_MT'), 
                         outs = ('TreatedWater', 'AnoxT_CH4', 'AnoxT_N2O'),
                         degraded_components=('OtherSS',),  ppl = 1000, baseline_ppl = 30,
                         )

    P_PAC_agitation = AgitationPump('P_PAC_agitation', ins='PAC', outs='PACAgitation', 
                                    #dP_design=405300,
                                    working_factor=0.9,  # The ratio of the actual output and the design output
                                    operation_time=12,  # Total run time of system or plant [h/d]
                                    pump_lifetime=5,  # Lifetime of the pump [years]
                                    pump_cost=200,
                                    dP_design=0,                                      
                                    )
    
    P_PAC_dosing = DosingPump('P_PAC_dosing', ins=P_PAC_agitation-0, outs ='PACDosing')
    B_AeroT = su.EL_blower('B_AeroT', ins='Air', outs ='Air',
                            F_BM={
                                  'Blowers': 2.22,
                                  'Blower piping': 1,
                                  'Blower building': 1.11,
                                 },
                            lifetime=15, lifetime_unit='yr',
                            units={
                                  'Total gas flow': 'CFM',
                                  'Blower capacity': 'CFM',
                                  'Number of blowers': '',
                                  'Total blower power': 'kW',
                                 },
                            N_reactor=2, # the number of the reactors where the gas sparging modules will be installed
                            gas_demand_per_reactor=1, # gas demand per reactor
                            TDH=6, # total dynamic head for rhe blower, in psi
                            eff_blower=0.7, # efficiency of the blower in fraction
                            eff_motor=0.7, # efficiency of the motor in fraction
                            AFF=3.33, # air flow fraction
                            building_unit_cost=9, # unit cost of the building, in USD/ft2
                            ppl = 1000, baseline_ppl = 30,
                           )
    AeroT = su.EL_Aerobic('AeroT', ins=(AnoxT-0, P_PAC_dosing-0, B_AeroT-0), 
                          outs = ('TreatedWater', 'AeroT_CH4', 'AeroT_N2O'), 
                          ppl = 1000, baseline_ppl = 30,
                          )

    B_MembT = su.EL_blower('B_MembT', ins = 'Air', outs = 'Air', 
                            F_BM={
                                  'Blowers': 2.22,
                                  'Blower piping': 1,
                                  'Blower building': 1.11,
                                 },
                            lifetime=15, lifetime_unit='yr',
                            units={
                                  'Total gas flow': 'CFM',
                                  'Blower capacity': 'CFM',
                                  'Number of blowers': '',
                                  'Total blower power': 'kW',
                                 },
                            N_reactor=2, # the number of the reactors where the gas sparging modules will be installed
                            gas_demand_per_reactor=1, # gas demand per reactor
                            TDH=6, # total dynamic head for rhe blower, in psi
                            eff_blower=0.7, # efficiency of the blower in fraction
                            eff_motor=0.7, # efficiency of the motor in fraction
                            AFF=3.33, # air flow fraction
                            building_unit_cost=9, # unit cost of the building, in USD/ft2
                            ppl = 1000, baseline_ppl = 30,)
    P_NitrateReturn_PC = ReturnPump('P_NitrateReturn_PC', ins='MembT_return', outs=1-PC, dP_design=0,)

    P_NitrateReturn_AnoxT = ReturnPump('P_NitrateReturn_AnoxT', ins='MembT_return', outs=1-AnoxT, dP_design=0,) 
    MembT = su.EL_MBR('MembT', ins=(AeroT-0, B_MembT-0), 
                      outs = ('TreatedWater', 1-P_NitrateReturn_PC, 1-P_NitrateReturn_AnoxT, 'MemT_CH4', 'MemT_N2O'),
                      ppl = 1000,
                      baseline_ppl = 30,
                      )

    P_MT_selfpriming = SelfPrimingPump('P_MT_selfpriming', ins=MembT-0, outs='SelfPrimingWater', dP_design=405300,)
    
    # P_O3_gen = O3GenPump('P_O3_gen', ins=None, outs='O3', dP_design=405300); 
    P_O3_dosing = MicroBubblePump('P_O3_dosing', ins='P_O3_gen', outs='DosingO3', dP_design=0,)

    P_AirDissolved = AirDissolvedPump('P_AirDissolved', ins='CWTWater', outs='Water_With_Oxyen', dP_design=50000,)

    CWT = su.EL_CWT('CWT', ins=(P_MT_selfpriming-0, P_O3_dosing-0, P_AirDissolved-0), 
                    outs= ('ClearWater', 1-CT, 0-P_AirDissolved, 'CWT_CH4', 'CWT_N2O'), 
                    V_wf = 0.9, 
                    ppl = 1000, baseline_ppl = 30,)
    
    P_CWT = ClearWaterPump('P_CWT', ins=CWT-0, outs='ReuseWater', dP_design=405300,)
    
    PT = su.StorageTank('PT', ins=P_CWT-0, outs=2-Toilet, vessel_material = None, V_wf = None, 
                        include_construction = True, length_to_diameter = None, 
                        F_BM_default = 1, kw_per_m3 = None, vessel_type = None, tau = None, 
                        ppl = 1000, baseline_ppl = 30,)

    Total_CH4 = su.Mixer('Total_CH4', ins=(Toilet-1, CT-1, PC-3, AnoxT-1, AeroT-1, MembT-3, ), outs=stream['CH4'])
    Total_CH4.add_specification(lambda: add_fugitive_items(Total_CH4, 'CH4_item'))
    Total_CH4.line = 'fugitive CH4 mixer'

    Total_N2O = su.Mixer('Total_N2O', ins=(Toilet-2, CT-2, PC-4, AnoxT-2, AeroT-2, MembT-4,), outs=stream['N2O'])
    Total_N2O.add_specification(lambda: add_fugitive_items(Total_N2O, 'N2O_item'))
    Total_N2O.line = 'fugitive N2O mixer'
    
    # Other impacts and costs
    Other_system = su.EL_System('Other_system', ins=PT-0, outs=('housing',), 
                                ppl = 1000, baseline_ppl = 30, if_gridtied=True)
    Other_housing = su.EL_Housing('Other_housing', ins=(Other_system-0), outs='Transport', ppl=1000, baseline_ppl=30)
    Other_WasteTransport = Trucking('Other_WasteTransport', ins=Other_system-0, outs=('WasteTransport', 'ConveyanceLoss'), 
                                       load = 20, # transportation load per trip
                                       load_unit='kg',
                                       load_type='mass', # mass or volume
                                       distance=5.0, # transportation distance per trip 
                                       distance_unit='km', loss_ratio=0.02,
                                       interval= 365, # timeinterval between trips
                                       interval_unit='d',
                                       ) # here transport fee for waste clearance is considered
    Other_SetupShipping = Trucking('Other_SetupShipping', ins=(Other_system-0), outs=None,
                                       load = 20, # transportation load per trip
                                       load_unit='kg',
                                       load_type='mass', # mass or volume
                                       distance=5.0, # transportation distance per trip 
                                       distance_unit='km', loss_ratio=0.0,
                                       interval= 365, # timeinterval between trips
                                       interval_unit='d',
                                       )
                                                                           
    sysEL = System('sysEL', path=(WasteWaterGenerator,
                                  Toilet,
                                  CT,
                                  P_CT_lift,
                                  PC, P_PC_return,
                                  P_Glu_agitation,
                                  P_Glu_dosing,
                                  P_AnoxT_agitation,
                                  AnoxT,
                                  P_PAC_agitation,
                                  P_PAC_dosing,
                                  B_AeroT,
                                  AeroT,
                                  B_MembT,
                                  P_NitrateReturn_PC,
                                  P_NitrateReturn_AnoxT,
                                  MembT,
                                  P_MT_selfpriming,
                                  P_O3_dosing,
                                  P_AirDissolved,
                                  CWT,
                                  P_CWT,
                                  PT,
                                  Total_CH4,
                                  Total_N2O,
                                  Other_system,
                                  Other_housing,
                                  Other_WasteTransport,
                                  ))
    sysEL.simulate()
    
    teaEL = TEA(system = sysEL,
                   discount_rate=0.05,  
                   income_tax=0.05,  
                   CEPCI=567.5,  
                   start_year=2024,  
                   lifetime=10,  
                   uptime_ratio=1.0,  
                   CAPEX=0.0,  
                   lang_factor=None,  
                   annual_maintenance=0.0,  
                   annual_labor=0.0,  
                   system_add_OPEX={},
                   depreciation='SL',  
                   construction_schedule=(0, 1),  
                   accumulate_interest_during_construction=False,  
                   simulate_system=True,  
                   simulate_kwargs={},  
                   #**tea_kwargs, how to define **tea_kwargs here(???),
                 )
    
    get_power = lambda: sum([(u.power_utility.rate * u.uptime_ratio) for u in sysEL.units]) * (365 * teaEL.lifetime) * 12
    LCA(system = sysEL, lifetime = 10, lifetime_unit = 'yr', uptime_ratio = 1.0, e_item = get_power)

    return sysEL