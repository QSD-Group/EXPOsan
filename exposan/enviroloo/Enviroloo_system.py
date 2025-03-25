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
# from qsdsan.utils import load_components, set_thermo


import qsdsan as qs
import os
from qsdsan import Component, Components, set_thermo as qs_set_thermo
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

    #create_stream_with_impact_item(stream_ID='ammonium')
    #create_stream_with_impact_item(stream_ID='struvite')
    #create_stream_with_impact_item(stream_ID='NaOH')
    #create_stream_with_impact_item(stream_ID='NaClO')
    create_stream_with_impact_item(stream_ID='O3')
    create_stream_with_impact_item(stream_ID='PAC')
    create_stream_with_impact_item(stream_ID='Glucose')
    create_stream_with_impact_item(stream_ID='air')

def update_toilet_param(unit):
    # Use the private attribute so that the number of users/toilets will be exactly as assigned
    # (i.e., can be fractions)
    unit._N_user = get_toilet_users()
    unit._N_toilet = ppl / get_toilet_users()
    unit._run()

def update_carbon_COD_ratio(sys):
    for first_u in sys.units:
        if hasattr(first_u, 'carbon_COD_ratio'):
            carbon_COD_ratio = first_u.carbon_COD_ratio
            break
    for u in sys.units:
        if u is first_u: continue
        if hasattr(u, 'carbon_COD_ratio'): u.carbon_COD_ratio = carbon_COD_ratio

# %% Create EnviroLoo Clear system

def create_components(set_thermo = True
                      #adjust_MW_to_measured_as=False
                      ):
    # bw_cmps = create_bw_components(set_thermo=False)
    masm2d_cmps = pc.create_masm2d_cmps(set_thermo=True)

    # C = Component('C', phase='l', particle_size='Soluble', degradability='Undegradable', organic=False)
    
    # SolubleCH4 = Component('SolubleCH4', search_ID='CH4', phase='l', particle_size='Soluble', degradability='Slowly', organic=True)

    # PAC = Component('PAC', search_ID='10124-27-3', phase='s', particle_size='Particulate', degradability='Slowly', organic=False)
    # add_V_from_rho(PAC, rho=2800)
                    
    # # Glucose = Component('Glucose', search_ID='50-99-7', phase='s', particle_size='Particulate', degradability='Readily', organic=False)
    # # add_V_from_rho(Glucose, rho=1560)

    # Glucose = Component('Glucose', search_ID='50-99-7', phase='l', particle_size='Soluble', degradability='Readily', organic=True)
    # add_V_from_rho(Glucose, rho=1560)
          
    # O3 = Component('O3', search_ID='10028-15-6', phase='g', particle_size='Dissolved gas', degradability='Readily', organic=False)
          
          
    # NaOH = Component('NaOH', search_ID='1310-73-2', phase='s', particle_size='Particulate', degradability='Readily', organic=False)
    # add_V_from_rho(NaOH, rho=2130)
          
    # NaClO = Component('NaClO', search_ID='7681-52-9', phase='s', particle_size='Particulate', degradability='Readily', organic=False)
    # add_V_from_rho(NaClO, rho=1250)
    
    # # NO3 = Component('NO3', measured_as = 'N', phase='l', particle_size='Soluble', degradability='Undegradable', organic=False)
    # # add_V_from_rho(NO3, rho=1.15) # need check
          
    # #NH3_l = Component('NH3_l', measured_as = 'N', phase='l', particle_size='Soluble', degradability='Undegradable', organic=False)
          
    # #NonNH3 = Component('NonNH3', formula = 'N', measured_as = 'N', phase='l', particle_size='Soluble', degradability='Undegradable', organic=False, description='Non-NH3');
          
    # air = Component('air', search_ID='17778-88-0', phase='g', particle_size='Dissolved gas',
    #                 degradability='Readily', organic=False)
    #                 # 1.204 kg/m3, cited from https://en.wikipedia.org/wiki/Density_of_air#:~:text=Air%20density%2C%20like%20air%20pressure,International%20Standard%20Atmosphere%20(ISA).
    # add_V_from_rho(air, rho=1.204)

    #allowed_values = {
    #'particle_size': ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'),
    #'degradability': ('Readily', 'Slowly', 'Undegradable'),
    #'organic': (True, False)}
          
    cmps = Components((masm2d_cmps 
                       # C, SolubleCH4, 
                       # #H2O, CO2, CH4, N2O, NH3
                       # Glucose, O3, air, PAC, NaOH, NaClO
                       ))
    
    # for i in cmps:
    #     for attr in ('HHV', 'LHV', 'Hf'):
    #         if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()

    # cmps.set_alias('H2O', 'Water')
    # #cmps.set_alias('CO2', 'Carbon Dioxide')
    # cmps.set_alias('CH4', 'Methane')
    # if set_thermo: qs_set_thermo(cmps)

    return cmps



def create_systemEL(flowsheet = None):
    
    # Components and stream
    cmps = create_components()
    qs.set_thermo(cmps)
    thermo_masm2d = qs.get_thermo()
    
    flowsheet = flowsheet or main_flowsheet
    streamEL = flowsheet.stream
    batch_create_streams('EL')
    
    WasteWaterGenerator = elu.EL_Excretion('WasteWaterGenerator', outs=('urine', 'feces'))

    # Toilet = EL_MURT('Toilet', ins=(WasteWaterGenerator-0, WasteWaterGenerator-1, 'FlushingWater', 'ToiletPaper'), 
    #                 outs =('MixedWasteWater', 'Toilet_CH4', 'Toilet_N2O'),
    #                 decay_k_COD = get_decay_k(),
    #                 decay_k_N = get_decay_k(),
    #                 max_CH4_emission = max_CH4_emission,
    #                 N_user = get_toilet_users(), 
    #                 N_toilet = ppl / get_toilet_users(),
    #                 if_flushing = True, if_desiccant = False, if_toilet_paper = True,
    #                 # CAPEX = 0, # capital cost of a single toilet
    #                 CAPEX=500*max(1, ppl/100),
    #                 OPEX_over_CAPEX = 0.06, # fraction of annual operating cost over total capital cost
    #                 )
    # Toilet.add_specification(lambda: update_toilet_param(Toilet))
    
    Toilet = elu.EL_MURT('Toilet',
                    ins=(WasteWaterGenerator-0, WasteWaterGenerator-1, 'toilet_paper', 'flushing_water', 'cleansing_water', 'desiccant'),
                    outs=('mixed_waste', 'Toilet_CH4', 'Toilet_N2O'),
                    # N_user=get_toilet_users(), N_tot_user=ppl,
                    N_user=get_toilet_users(), N_tot_user=ppl,
                    lifetime=10, if_include_front_end=True,
                    if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                    if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                    CAPEX=500*max(1, ppl/100), OPEX_over_CAPEX=0.06,
                    decay_k_COD=get_decay_k(),
                    decay_k_N=get_decay_k(),
                    max_CH4_emission=max_CH4_emission
                    )
    Toilet.add_specification(lambda: update_toilet_param(Toilet))
    # Toilet = EL_MURT('Toilet', ins=(WasteWaterGenerator-0, WasteWaterGenerator-1, 'toilet_paper', 'flushing_water', 'cleansing_water', 'desiccant'),
    #             outs=('mixed_waste', 'Toilet_CH4', 'Toilet_N2O'),
    #             decay_k_COD=get_decay_k(),
    #             decay_k_N=get_decay_k(),
    #             max_CH4_emission=max_CH4_emission,
    #             N_user=25, N_tot_user=ppl, lifetime=10,
    #             if_flushing=True, if_desiccant=False, if_toilet_paper=True,
    #             CAPEX=500*max(1, ppl/100), OPEX_over_CAPEX=0.06)    

    CT = elu.EL_CT('CT', ins=(Toilet-0, 'PrimaryClarP_return','PrimaryClar_spill', 'ClearWaterTank_spill'), 
                    outs = ('TreatedWater'),
                    # V_wf = 0.9, ppl = ppl, baseline_ppl = 100,
                    # kW_per_m3=0.1,  # The power consumption per unit volume of the tank
                    )
    
    P_CT_lift = LiftPump('P_CT_lift', ins = CT-0, outs = 'TreatedWater',
                            working_factor = 0.9,  # The ratio of the actual output and the design output
                            operation_time = 12,  # Total run time of system or plant [h/d]
                            life_time = 5,  # Lifetime of the pump [years]
                            pump_cost = 104.0, # USD from Alibaba https://www.alibaba.com/product-detail/Jinmu-brand-professional-1-1kw-submersible_1600929097854.html?spm=a2700.details.you_may_like.3.44e35ef2MlxhPU
                            dP_design = 0,
                            )
    
    PC = elu.EL_PC('PC', ins=(P_CT_lift-0, 'NitrateReturn_MT'), outs=('TreatedWater', 'PC_return' , 2-CT),
                    ppl = ppl,  # The number of people served
                    baseline_ppl = 100,
                    solids_removal_efficiency = 0.85,  # The solids removal efficiency
                    sludge_flow_rate = 0.5,  # Sludge flow rate
                    max_oveflow = 15,
                    )
    
    P_PC_return = ReturnPump('P_PC_return', ins=PC-1, outs = 1-CT, 
                                working_factor = 0.9,  # The ratio of the actual output and the design output
                                operation_time = 12,  # Total run time of system or plant [h/d]
                                life_time = 5,  # Lifetime of the pump [years]
                                pump_cost = 123.76, # USD from Alibaba https://www.alibaba.com/product-detail/0-4kw-cast-iron-motor-housing_1600934836942.html?spm=a2700.galleryofferlist.normal_offer.d_title.3be013a0L7StzT
                                dP_design = 0,
                                )
    
    P_Glu_agitation = AgitationPump('P_Glu_agitation', 
                                    ins = streamEL['Glucose'], 
                                    #ins=Glucose_in-0,
                                    outs = 'GlucoseAgitation', 
                                    working_factor = 0.9,  # The ratio of the actual output and the design output
                                    operation_time = 12,  # Total run time of system or plant [h/d]
                                    life_time = 5,  # Lifetime of the pump [years]
                                    pump_cost = 696.30, # USD from https://www.grainger.com/product/DAYTON-Open-Drum-Mixer-115-230V-AC-32V133?opr=PLADS&analytics=FM%3APLA&a2c_sku_original=32V138&position=2
                                    dP_design = 0,
                                    )
    
    P_Glu_dosing = DosingPump('P_Glu_dosing', ins = P_Glu_agitation-0, outs = 'GlucoseDosing',
                                working_factor = 0.9,  # The ratio of the actual output and the design output
                                operation_time = 12,  # Total run time of system or plant [h/d]
                                life_time = 5,  # Lifetime of the pump [years]
                                pump_cost = 59, # USD from https://www.aliexpress.us/item/3256804645639765.html?src=google&gatewayAdapt=glo2usa
                                dP_design = 0,
                                ) 
    
    P_AnoxT_agitation = AgitationPump('P_AnoxT_agitation', ins= None, outs='AgitationWater', 
                                        working_factor = 0.9,  # The ratio of the actual output and the design output
                                        operation_time = 12,  # Total run time of system or plant [h/d]
                                        life_time = 5,  # Lifetime of the pump [years]
                                        pump_cost = 696.30, # USD from https://www.grainger.com/product/DAYTON-Open-Drum-Mixer-115-230V-AC-32V133?opr=PLADS&analytics=FM%3APLA&a2c_sku_original=32V138&position=2
                                        dP_design = 0,
                                        )
    
    AnoxT = elu.EL_Anoxic('AnoxT', ins=(PC-0, 'NitrateReturn_MT', P_Glu_dosing-0, P_AnoxT_agitation-0), 
                            outs = ('TreatedWater', 'AnoxT_CH4', 'AnoxT_N2O'),
                            degraded_components=('OtherSS',),  
                            ppl = ppl, baseline_ppl = 100,
                            )
    

    
    
    P_PAC_agitation = AgitationPump('P_PAC_agitation', ins=streamEL['PAC'], outs='PACAgitation', 
                                    working_factor = 0.9,  # The ratio of the actual output and the design output
                                    operation_time = 12,  # Total run time of system or plant [h/d]
                                    life_time = 5,  # Lifetime of the pump [years]
                                    pump_cost = 696.30, # USD from https://www.grainger.com/product/DAYTON-Open-Drum-Mixer-115-230V-AC-32V133?opr=PLADS&analytics=FM%3APLA&a2c_sku_original=32V138&position=2
                                    dP_design=0,                                      
                                    )
    
    P_PAC_dosing = DosingPump('P_PAC_dosing', ins=P_PAC_agitation-0, outs ='PACDosing', 
                                working_factor = 0.9,
                                operation_time = 12,
                                life_time = 5,
                                pump_cost = 59, # USD from https://www.aliexpress.us/item/3256804645639765.html?src=google&gatewayAdapt=glo2usa
                                dP_design = 0,
                                )
    # B_AeroT = elu.EL_blower('B_AeroT', ins=streamEL['air'], outs ='Air_to_aerobic',
    #                         # F_BM={
    #                         #       'Blowers': 2.22,
    #                         #       'Blower piping': 1,
    #                         #       'Blower building': 1.11,
    #                         #      },
    #                         F_BM = 2.2,
    #                         lifetime = 10, lifetime_unit='yr',
    #                         # units={
    #                         #       'Total gas flow': 'CFM',
    #                         #       'Blower capacity': 'CFM',
    #                         #       'Number of blowers': '',
    #                         #       'Total blower power': 'kW',
    #                         #      },
    #                         # N_reactor=2, # the number of the reactors where the gas sparging modules will be installed
    #                         # gas_demand_per_reactor=1, # gas demand per reactor
    #                         # TDH=6, # total dynamic head for rhe blower, in psi
    #                         # eff_blower=0.85, # efficiency of the blower in fraction
    #                         # eff_motor=0.95, # efficiency of the motor in fraction
    #                         # AFF=3.33, # air flow fraction
    #                         # building_unit_cost=9, # unit cost of the building, in USD/ft2
    #                         ppl = ppl, baseline_ppl = 100,
    #                        )
    # B_AeroT.line = 'Air to aerobic tank'
    
    AeroT = elu.EL_Aerobic('AeroT', ins=(AnoxT-0, P_PAC_dosing-0), 
                                         # B_AeroT-0), 
                            outs = ('TreatedWater', 'AeroT_CH4', 'AeroT_N2O'), 
                            ppl = ppl, baseline_ppl = 100,
                            )
    
    # B_MembT = elu.EL_blower('B_MembT', ins = streamEL['air'], outs = 'Air_to_membrane', 
    #                         # F_BM={
    #                         #       'Blowers': 2.22,
    #                         #       'Blower piping': 1,
    #                         #       'Blower building': 1.11,
    #                         #      },
    #                         F_BM = 2.2,
    #                         lifetime=10, lifetime_unit='yr',
    #                         # units={
    #                         #       'Total gas flow': 'CFM',
    #                         #       'Blower capacity': 'CFM',
    #                         #       'Number of blowers': '',
    #                         #       'Total blower power': 'kW',
    #                         #      },
    #                         # N_reactor=2, # the number of the reactors where the gas sparging modules will be installed
    #                         # gas_demand_per_reactor=1, # gas demand per reactor
    #                         # TDH=6, # total dynamic head for rhe blower, in psi
    #                         # eff_blower=0.85, # efficiency of the blower in fraction
    #                         # eff_motor=0.95, # efficiency of the motor in fraction
    #                         # AFF=3.33, # air flow fraction
    #                         # building_unit_cost=9, # unit cost of the building, in USD/ft2
    #                         ppl = ppl, baseline_ppl = 100,
    #                         )
    # B_MembT.line = 'Air to membrane tank'
    
    P_NitrateReturn_PC = ReturnPump('P_NitrateReturn_PC', ins='MembT_return', outs=1-PC,
                                    working_factor = 0.9,  # The ratio of the actual output and the design output
                                    operation_time = 12,  # Total run time of system or plant [h/d]
                                    life_time = 5,  # Lifetime of the pump [years]
                                    pump_cost = 123.76, # USD from Alibaba https://www.alibaba.com/product-detail/0-4kw-cast-iron-motor-housing_1600934836942.html?spm=a2700.galleryofferlist.normal_offer.d_title.3be013a0L7StzT
                                    dP_design = 0, 
                                    )

    P_NitrateReturn_AnoxT = ReturnPump('P_NitrateReturn_AnoxT', ins='MembT_return', outs=1-AnoxT, 
                                        working_factor = 0.9,  # The ratio of the actual output and the design output
                                        operation_time = 12,  # Total run time of system or plant [h/d]
                                        life_time = 5,  # Lifetime of the pump [years]
                                        pump_cost = 123.76, # USD from Alibaba https://www.alibaba.com/product-detail/0-4kw-cast-iron-motor-housing_1600934836942.html?spm=a2700.galleryofferlist.normal_offer.d_title.3be013a0L7StzT
                                        dP_design = 0,
                                        ) 
    MembT = elu.EL_CMMBR('MembT', ins=(AeroT-0), 
                                     # B_MembT-0), 
                        outs = ('TreatedWater', 0-P_NitrateReturn_PC, 0-P_NitrateReturn_AnoxT, 'MemT_CH4', 'MemT_N2O', 'Sludge'),
                        ppl = ppl,
                        baseline_ppl = 100,
                        )
    
    Solids_separation = su.ComponentSplitter('Solids_separation', ins=MembT-5,
                                             outs=(streamEL['sol_N'], streamEL['sol_P'], streamEL['sol_K'],
                                             'A_sol_non_fertilizers'),
                                             split_keys=(('NH3','NonNH3'), 'P', 'K')
                                             )
    
    P_MT_selfpriming = SelfPrimingPump('P_MT_selfpriming', ins=MembT-0, outs='SelfPrimingWater', 
                                        working_factor = 0.9,  # The ratio of the actual output and the design output
                                        operation_time = 12,  # Total run time of system or plant [h/d]
                                        life_time = 5,  # Lifetime of the pump [years]
                                        pump_cost = 155.24, # USD from Alibaba https://www.alibaba.com/product-detail/LEO-QDX-Series-Cast-Iron-Submersible_60671071414.html?spm=a2700.galleryofferlist.normal_offer.d_title.3be013a0BMhyxn
                                        dP_design = 0,
                                        )
    
    # P_O3_gen = O3GenPump('P_O3_gen', ins=None, outs='O3', dP_design=405300)
    P_O3_dosing = MicroBubblePump('P_O3_dosing', ins=streamEL['O3'], outs='DosingO3', 
                                    working_factor = 0.9,  # The ratio of the actual output and the design output
                                    operation_time = 12,  # Total run time of system or plant [h/d]
                                    life_time = 5,  # Lifetime of the pump [years]
                                    pump_cost = 160, # USD, assumption
                                    dP_design = 25331, # in Pa
                                    )

    P_AirDissolved = AirDissolvedPump('P_AirDissolved', ins=('CWTWater', streamEL['air']), outs='Water_With_Oxyen', 
                                        working_factor = 0.9,  # The ratio of the actual output and the design output
                                        operation_time = 12,  # Total run time of system or plant [h/d]
                                        life_time = 5,  # Lifetime of the pump [years]
                                        pump_cost = 155.24, # USD from https://www.alibaba.com/product-detail/LEO-QDX-Series-Cast-Iron-Submersible_60671071414.html?spm=a2700.galleryofferlist.normal_offer.d_title.3be013a0BMhyxn
                                        dP_design = 50000, # in Pa
                                        )

    CWT = elu.EL_CWT('CWT', ins=(P_MT_selfpriming-0, P_O3_dosing-0, P_AirDissolved-0), 
                    outs= ('ClearWater', 3-CT), 
                    # V_wf = 0.9, max_oveflow=15, 
                    # ppl = ppl, baseline_ppl = 100,
                    )
   
    P_CWT = ClearWaterPump('P_CWT', ins=CWT-0, outs='ReuseWater', 
                            working_factor = 0.9,  # The ratio of the actual output and the design output
                            operation_time = 12,  # Total run time of system or plant [h/d]
                            life_time = 5,  # Lifetime of the pump [years]
                            pump_cost = 118.46, # USD from https://www.alibaba.com/product-detail/SQD-0-75kw-1hp-high-pressure_1601032336743.html?spm=a2700.galleryofferlist.normal_offer.d_title.2b1413a0sqHM0o
                            dP_design = 202650, # in Pa
                            )
   
    PT = elu.EL_PT('PT', ins=P_CWT-0, outs=(3-Toilet, 'pipeline_system'), vessel_material = None, V_wf = None, 
                        include_construction = True, length_to_diameter = None, 
                        F_BM_default = 1, kW_per_m3 = 0.1, vessel_type = None, tau = None, 
                        ppl = ppl, baseline_ppl = 100,
                        )
    
    # CWT.add_specification(lambda: update_carbon_COD_ratio(sysEL))
    # CWT.run_after_specification = True
    
    Total_CH4 = su.Mixer('Total_CH4', ins=(Toilet-1, AnoxT-1, AeroT-1, MembT-3), outs=streamEL['CH4'])
    Total_CH4.add_specification(lambda: add_fugitive_items(Total_CH4, 'CH4_item'))
    Total_CH4.line = 'fugitive CH4 mixer'

    Total_N2O = su.Mixer('Total_N2O', ins=(Toilet-2, AnoxT-2, AeroT-2, MembT-4,), outs=streamEL['N2O'])
    Total_N2O.add_specification(lambda: add_fugitive_items(Total_N2O, 'N2O_item'))
    Total_N2O.line = 'fugitive N2O mixer'
    
    # Other impacts and costs
    Pipeline_system = elu.EL_System('Pipeline_system', ins=PT-1, 
                                # outs='PipelineConnection',
                                ppl = ppl, baseline_ppl = 100, if_gridtied=False)
    # #Other_housing = EL_Housing('Other_housing', ins=Other_system-0, outs='Transport', ppl = ppl, baseline_ppl = 30)
    # Other_WasteTransport = Trucking('Other_WasteTransport', ins = Other_system-0, outs = ('WasteTransport', 'ConveyanceLoss'), 
    #                                    load = 20, # transportation load per trip
    #                                    load_unit = 'kg',
    #                                    load_type = 'mass', # mass or volume
    #                                    distance = 5.0, # transportation distance per trip 
    #                                    distance_unit = 'km', loss_ratio=0.02,
    #                                    interval = 30, # timeinterval between trips
    #                                    interval_unit = 'd',
    #                                    ) # here transport fee for waste clearance is considered
    # Other_SetupShipping = Trucking('Other_SetupShipping', ins = (Other_system-0), outs = None,
    #                                    load = 53.472, # transportation load per trip, from BOM of Enviroloo
    #                                    load_unit = 'kg',
    #                                    load_type = 'mass', # mass or volume
    #                                    distance = 5000.0, # transportation distance per trip, assuming the delivery of setup installation
    #                                    distance_unit = 'km', loss_ratio=0.0,
    #                                    interval = 365, # timeinterval between trips
    #                                    interval_unit = 'd',
    #                                    )
    
    # ensure the same carbon_COD_ratio throughout the EL system
    # Other_SetupShipping.add_specification(lambda: update_carbon_COD_ratio(sysEL))
    # Other_SetupShipping.run_after_specification = True
                                                                           
    # sysEL_PCspill = System('sysEL_PCspill',
    #                  path = (WasteWaterGenerator, Toilet, CT, P_CT_lift, PC),
    #                  recycle = PC-2
    #                  )
    
    # sysEL_PCreturn = System('sysEL_PCreturn',
    #                   path = (sysEL_PCspill, P_PC_return),
    #                   recycle = P_PC_return-0
    #                   )
    
    # sysEL_PCreturn.simulate()
    # teaEL = TEA(system=sysEL_PCreturn, discount_rate=discount_rate,
    #        start_year=2024, lifetime=20, uptime_ratio=1,
    #        CEPCI = 567.5,
    #        CAPEX = 2.00,  
    #        #lang_factor=None,
    #        lang_factor=None,
    #        annual_maintenance=0,
    #        # annual_labor=(operator_daily_wage*3*365),
    #        annual_labor=0
    #        )
    # get_powerEL = lambda: sum([u.power_utility.rate for u in sysEL_PCreturn.units]) * (24 * 365 * teaEL.lifetime)
    # LCA(system=sysEL_PCreturn, lifetime=20, lifetime_unit='yr', uptime_ratio=1.0, e_item=get_powerEL)
    # return sysEL_PCreturn

    sysEL_PCrecycle = System('sysEL_PCspill',
                     path = (WasteWaterGenerator, Toilet, CT, P_CT_lift, PC, P_PC_return),
                     recycle = PC-2
                     )
    
    sysEL_CWTrecycle = System('sysEL_CWTrecycle',
                       path = (sysEL_PCrecycle, P_AnoxT_agitation, P_Glu_agitation, P_Glu_dosing, AnoxT, 
                               P_PAC_agitation, P_PAC_dosing, 
                               # B_AeroT, 
                               AeroT, 
                               # B_MembT, 
                               MembT, Solids_separation,
                               P_NitrateReturn_PC, P_NitrateReturn_AnoxT, P_MT_selfpriming, 
                               P_O3_dosing, P_AirDissolved, CWT), 
                       recycle = CWT-1
                       )
    
    sysEL = System('sysEL',
                             path = (sysEL_CWTrecycle, P_CWT, PT, Total_N2O, Total_CH4, Pipeline_system),
                             recycle = PT-0
                             )
    
    sysEL.simulate()
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
    return sysEL
    
    # sysEL_PCspill = System('sysEL_PCspill',
    #                  path = (sys1, P_CT_lift, PC),
    #                  recycle = PC-2
    #                  )
    # sysEL_PCspill.simulate()
    # teaEL = TEA(system=sysEL_PCspill, discount_rate=discount_rate,
    #        start_year=2024, lifetime=20, uptime_ratio=1,
    #        CEPCI = 567.5,
    #        CAPEX = 2.00,  
    #        #lang_factor=None,
    #        lang_factor=None,
    #        annual_maintenance=0,
    #        # annual_labor=(operator_daily_wage*3*365),
    #        annual_labor=0
    #        )
    # sysEL_PCreturn = System('sysEL_PCreturn',
    #                   path = (sysEL_PCspill, P_PC_return),
    #                   recycle = P_PC_return-0
    #                   )
    
    # sysEL_MBR_PC = System('sysEL_MBR_PC',
    #                 path = (sysEL_PCreturn, P_AnoxT_agitation, P_Glu_agitation, P_Glu_dosing, AnoxT, 
    #                         P_PAC_agitation,P_PAC_dosing, AeroT, MembT, P_NitrateReturn_PC),
    #                 recycle = P_NitrateReturn_PC-0
    #                 )
    
    # sysEL_MBR_AnoxT = System('sysEL_MBR_AnoxT',
    #                    path = (sysEL_MBR_PC, P_NitrateReturn_AnoxT),
    #                    recycle = P_NitrateReturn_AnoxT-0
    #                    )
    
    # sysEL_CWT = System('sysEL_PT',
    #              path = (sysEL_MBR_AnoxT, P_MT_selfpriming, P_O3_dosing, CWT),
    #              recycle = CWT-1
    #              )
    
    # sysEL = System('sysEL',
    #             path = (sysEL_CWT, P_CWT, PT),
    #             recycle = PT-0
    #             )
    
    # sysEL.simulate()
    
    # sysEL = System('sysEL', path=(WasteWaterGenerator,
    #                               Toilet,
    #                               CT,
    #                               P_CT_lift,
    #                               PC, P_PC_return,
    #                               P_Glu_agitation,
    #                               P_Glu_dosing,
    #                               P_AnoxT_agitation,
    #                               AnoxT,
    #                               P_PAC_agitation,
    #                               P_PAC_dosing,
    #                               B_AeroT,
    #                               AeroT,
    #                               B_MembT,
    #                               P_NitrateReturn_PC,
    #                               P_NitrateReturn_AnoxT,
    #                               MembT,
    #                               P_MT_selfpriming,
    #                               P_O3_dosing,
    #                               P_AirDissolved,
    #                               CWT,
    #                               P_CWT,
    #                               PT,
    #                               Total_CH4,
    #                               Total_N2O,
    #                               Pipeline_system,
    #                               # Other_housing,
    #                               # Other_WasteTransport,
    #                               ))
    # sysEL.simulate()
    
    # teaEL = TEA(system = sysEL,
    #                discount_rate = discount_rate,  
    #                income_tax = 0.05,  
    #                CEPCI = 567.5,  
    #                start_year = 2024,  
    #                lifetime = 10,  
    #                uptime_ratio = 1.0,  
    #                #CAPEX = 0.0,  
    #                lang_factor = None,  
    #                annual_maintenance = 0.0,  
    #                annual_labor = (operator_daily_wage * 3 * 365),  
    #                system_add_OPEX = {},
    #                depreciation = 'SL',  
    #                construction_schedule = (0, 1),  
    #                accumulate_interest_during_construction = False,  
    #                simulate_system = True,  
    #                simulate_kwargs = {},
    #                )
    
    # teaEL = TEA(system=sysEL, discount_rate=discount_rate,
    #        start_year=2024, lifetime=20, uptime_ratio=1,
    #        CEPCI = 567.5,
    #        CAPEX = 2.00,  
    #        #lang_factor=None,
    #        lang_factor=None,
    #        annual_maintenance=0,
    #        # annual_labor=(operator_daily_wage*3*365),
    #        annual_labor=0
    #        )
   
    # teaEL = TEA(system=sysEL_PT, discount_rate=discount_rate,
    #        start_year=2024, lifetime=20, uptime_ratio=1,
    #        CEPCI = 567.5,
    #        CAPEX = 2.00,  
    #        #lang_factor=None,
    #        lang_factor=None,
    #        annual_maintenance=0,
    #        # annual_labor=(operator_daily_wage*3*365),
    #        annual_labor=0
    #        )
   
    # get_powerEL = lambda: sum([u.power_utility.rate for u in sysEL.units]) * (24 * 365 * teaEL.lifetime)
    # #get_powerEL = lambda: sum([(getattr(u.power_utility, 'rate', 0) * u.uptime_ratio) for u in sysEL.units]) * (365 * teaEL.lifetime) * 12
    
    # LCA(system=sysEL, lifetime=20, lifetime_unit='yr', uptime_ratio=1.0, e_item=get_powerEL)

    # return sysEL

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