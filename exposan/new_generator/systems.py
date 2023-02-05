#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Shion Watabe <shionwatabe@gmail.com>
    
    Hannah Lohman <hlohman94@gmail.com>
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import (
    Flowsheet, main_flowsheet,
    WasteStream,
    sanunits as su,
    ImpactItem,
    System, SimpleTEA, LCA,
    )
from qsdsan.utils import clear_lca_registries
from exposan.utils import add_fugitive_items
try: from exposan import new_generator
except:
    raise RuntimeError(
        'No NEWgenerator unit operations available from QSDsan, '
        'most likely due to the use of the public QSDsan. '
        'Detailed information on the unit operations in the NEWgenerator system '
        'is protected by non-disclosure agreement, please reach out to the '
        'corresponding authors of the papers listed in README for full access.'
        )
from exposan.new_generator import (
    _load_components,
    _load_lca_data,
    default_ppl,
    discount_rate,
    get_decay_k,
    max_CH4_emission,
    update_resource_recovery_settings,
    )

__all__ = ('create_system',)


# =============================================================================
# Universal units and functions
# =============================================================================

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

    create_stream_with_impact_item(stream_ID='conc_NH3')
    create_stream_with_impact_item(stream_ID='biogas')
    create_stream_with_impact_item(stream_ID='GAC')
    create_stream_with_impact_item(stream_ID='zeolite')
    create_stream_with_impact_item(stream_ID='NaCl')
    create_stream_with_impact_item(stream_ID='NaCl1', dct_key='NaCl')
    create_stream_with_impact_item(stream_ID='NaOH')
    create_stream_with_impact_item(stream_ID='LPG')


# %%

# =============================================================================
# Scenario A (sysA):
# Original NEWgenerator, with front-end and photovoltaic power
# biogas is captured and flared
# no sludge trucking
# =============================================================================

def create_systemA(flowsheet=None, ppl=default_ppl):
    # Set flowsheet to avoid stream replacement warnings
    flowsheet = flowsheet or main_flowsheet
    streamA = flowsheet.stream
    batch_create_streams('A')

    ##### Human Inputs #####
    A1 = su.Excretion('A1', outs=('urine','feces'))

    ##### User Interface #####
    A2 = su.MURT('A2',
                 ins=(A1-0, A1-1, 'toilet_paper', 'flushing_water', 'cleansing_water', 'desiccant'),
                 outs=('mixed_waste', 'A2_CH4', 'A2_N2O'),
                 N_user=25, N_tot_user=ppl,
                 lifetime=10, if_include_front_end=True,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=500*max(1, ppl/100), OPEX_over_CAPEX=0.06,
                 decay_k_COD=get_decay_k(),
                 decay_k_N=get_decay_k(),
                 max_CH4_emission=max_CH4_emission)

    ##### Treatment #####
    A3 = su.NEWgeneratorAnMBR('A3', ins=A2-0, outs=('A3_treated', 'sludge', streamA.biogas,'A3_CH4', 'A3_N20'),
                              if_capture_biogas=True, if_N2O_emission=True,
                              decay_k_COD=get_decay_k(),
                              decay_k_N=get_decay_k(),
                              MCF_decay=0.1, N_max_decay=0.8, N2O_EF_decay=0.0005,
                              if_gridtied=False, if_HFMC_membrane=True, user_scale_up=ppl/100,
                              sludge_moisture_content=0.93) # communication with the design team

    A4 = su.SludgePasteurization('A4', ins=(A3-2, 'air', A3-1, streamA['LPG']),
                                 outs=('used', 'lost', 'wasted', 'treated_sludge'),
                                 if_biogas=True, heat_loss=0.1, target_MC=0.1,
                                 sludge_temp=283.15, temp_pasteurization=343.15,
                                 if_combustion=True, biogas_loss=0.1, biogas_eff=0.55,
                                 hhv_lpg=49.3, hhv_methane=55.5,
                                 user_scale_up=None, exponent_scale=1, if_sludge_service=True)
    # Clear power usage for the solar scenario
    old_A4_cost = A4._cost
    def A4_cost_no_power():
        old_A4_cost()
        A4.power_utility.empty()
    A4._cost = A4_cost_no_power

    A5 = su.NEWgeneratorIonExchange('A5', ins=(A3-0, streamA['zeolite'], streamA['GAC'], streamA['NaCl1'], streamA['NaOH']),
                                    outs=('A5_treated', 'spent_zeolite', 'spent_GAC',streamA['conc_NH3']),
                                    decay_k_COD=get_decay_k(),
                                    decay_k_N=get_decay_k(),
                                    if_gridtied=False, if_larger_NCS=False,
                                    if_zeolite_replacement_only=False, user_scale_up=ppl/100)

    A6 = su.NEWgeneratorChlorination('A6', ins=(A5-0, streamA['NaCl']), outs=('A6_treated', 'A6_solubleCH4'), if_gridtied=False, user_scale_up=ppl/100)

    # This unit has a capacity of 1160 W and a 14.7 kWh battery
    # enough to support the electricity consumption of the whole system (about 0.46 kW)
    # for ~32 hr (14.7/0.46), so over the night is OK
    A7 = su.NEWgeneratorPhotovoltaic('A7', ins=A6-0, outs='A7_in', if_lithium_battery=False, user_scale_up=ppl/100)
    A8 = su.NEWgeneratorControls('A8', ins=A7-0, outs='A8_in', if_gridtied=False)
    A9 = su.NEWgeneratorHousing('A9', ins=A8-0, outs='A9_in', if_cheaper_housing=False, if_larger_NCS_or_greater_150users=False, user_scale_up=ppl/100)
    A10 = su.NEWgeneratorPretreatment('A10', ins=A9-0, outs='A10_in', if_gridtied=False, user_scale_up=ppl/100)
    A11 = su.NEWgeneratorFoundation('A11', ins=A10-0, outs='A11_in', user_scale_up=ppl/100, if_larger_NCS=False)

    # CH4 emissions from MURT, AnMBR, lost CH4 during sludge pasteurization, CH4 remaining in effluent
    A12 = su.Mixer('A12', ins=(A2-1, A3-3, A4-1, A6-1), outs=streamA['CH4'])
    A12.add_specification(lambda: add_fugitive_items(A12, 'CH4_item'))
    A12.line = 'fugitive CH4 mixer'

    # N2O emissions from MURT and AnMBR
    A13 = su.Mixer('A13', ins=(A2-2, A3-4), outs=streamA['N2O'])
    A13.add_specification(lambda: add_fugitive_items(A13, 'N2O_item'))
    A13.line = 'fugitive N2O mixer'

    A14 = su.ComponentSplitter('A14', ins=A4-3,
                               outs=(streamA['sol_N'], streamA['sol_P'], streamA['sol_K'],
                                     'A_sol_non_fertilizers'),
                               split_keys=('N', 'P', 'K'))

    ##### Simulation, TEA, and LCA #####
    sysA = System('sysA', path=(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14))

    teaA = SimpleTEA(system=sysA, discount_rate=discount_rate,
                     start_year=2020, lifetime=25, uptime_ratio=1,
                     lang_factor=None, annual_maintenance=0,
                     annual_labor=0)

    get_powerA = sum([u.power_utility.rate for u in sysA.units]) * (24 * 365 * teaA.lifetime)

    LCA(system=sysA, lifetime=25, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerA)

    return sysA


# =============================================================================
# Scenario B (sysB):
# original NEWgenerator with front-end and grid-tied system
# biogas is captured and flared
# no sludge trucking
# =============================================================================

def create_systemB(flowsheet=None, ppl=default_ppl):
    flowsheet = flowsheet or main_flowsheet
    streamB = flowsheet.stream
    batch_create_streams('B')

    ##### Human Inputs #####
    B1 = su.Excretion('B1', outs=('urine','feces'))

    ##### User Interface #####
    B2 = su.MURT('B2',
                 ins=(B1-0, B1-1, 'toilet_paper', 'flushing_water', 'cleansing_water', 'desiccant'),
                 outs=('mixed_waste', 'B2_CH4', 'B2_N2O'),
                 N_user=25, N_tot_user=ppl,
                 lifetime=10, if_include_front_end=True,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=500*max(1, ppl/100), OPEX_over_CAPEX=0.06,
                 decay_k_COD=get_decay_k(),
                 decay_k_N=get_decay_k(),
                 max_CH4_emission=max_CH4_emission)

    ##### Treatment #####
    B3 = su.NEWgeneratorAnMBR('B3', ins=B2-0, outs=('B3_treated', 'sludge', streamB.biogas,'B3_CH4', 'B3_N20'),
                              if_capture_biogas=True, if_N2O_emission=True,
                              decay_k_COD=get_decay_k(),
                              decay_k_N=get_decay_k(),
                              MCF_decay=0.1, N_max_decay=0.8, N2O_EF_decay=0.0005,
                              if_gridtied=True, if_HFMC_membrane=True, user_scale_up=ppl/100,
                              sludge_moisture_content=0.93) # communication with the design team

    B4 = su.SludgePasteurization('B4', ins=(B3-2, 'air', B3-1, streamB['LPG']),
                                 outs=('used', 'lost', 'wasted', 'treated_sludge'),
                                 if_biogas=True, heat_loss=0.1, target_MC=0.1,
                                 sludge_temp=283.15, temp_pasteurization=343.15,
                                 if_combustion=True, biogas_loss=0.1, biogas_eff=0.55,
                                 hhv_lpg=49.3, hhv_methane=55.5,
                                 user_scale_up=None, exponent_scale=1, if_sludge_service=True)

    B5 = su.NEWgeneratorIonExchange('B5', ins=(B3-0, streamB['zeolite'], streamB['GAC'], streamB['NaCl1'], streamB['NaOH']),
                                    outs=('B5_treated', 'spent_zeolite', 'spent_GAC', streamB['conc_NH3']), user_scale_up=ppl/100,
                                    decay_k_COD=get_decay_k(),
                                    decay_k_N=get_decay_k(),
                                    if_gridtied=True, if_larger_NCS=False,
                                    if_zeolite_replacement_only=False)

    B6 = su.NEWgeneratorChlorination('B6', ins=(B5-0, streamB['NaCl']), outs=('B6_treated', 'B6_solubleCH4'), if_gridtied=True, user_scale_up=ppl/100)

    B7 = su.NEWgeneratorGridtied('B7', ins=B6-0, outs='B7_in', user_scale_up=ppl/100)
    B8 = su.NEWgeneratorControls('B8', ins=B7-0, outs='B8_in', if_gridtied=True)
    B9 = su.NEWgeneratorHousing('B9', ins=B8-0, outs='B9_in', if_cheaper_housing=False, if_larger_NCS_or_greater_150users=False, user_scale_up=ppl/100)
    B10 = su.NEWgeneratorPretreatment('B10', ins=B9-0, outs='B10_in', if_gridtied=True, user_scale_up=ppl/100)
    B11 = su.NEWgeneratorFoundation('B11', ins=B10-0, outs='B11_in', if_larger_NCS=False, user_scale_up=ppl/100)

    # CH4 emissions from MURT, AnMBR, lost CH4 during sludge pasteurization, CH4 remaining in effluent
    B12 = su.Mixer('B12', ins=(B2-1, B3-3, B4-1, B6-1), outs=streamB['CH4'])
    B12.add_specification(lambda: add_fugitive_items(B12, 'CH4_item'))
    B12.line = 'fugitive CH4 mixer'

    # N2O emissions from MURT and AnMBR
    B13 = su.Mixer('B13', ins=(B2-2, B3-4), outs=streamB['N2O'])
    B13.add_specification(lambda: add_fugitive_items(B13, 'N2O_item'))
    B13.line = 'fugitive N2O mixer'

    B14 = su.ComponentSplitter('B14', ins=B4-3,
                               outs=(streamB['sol_N'], streamB['sol_P'], streamB['sol_K'],
                                     'B_sol_non_fertilizers'),
                               split_keys=('N', 'P', 'K'))

    ##### Simulation, TEA, and LCA #####
    sysB = System('sysB', path=(B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14))

    teaB = SimpleTEA(system=sysB, discount_rate=discount_rate,
                     start_year=2020, lifetime=25, uptime_ratio=1,
                     lang_factor=None, annual_maintenance=0,
                     annual_labor=0)

    get_powerB = sum([u.power_utility.rate for u in sysB.units]) * (24 * 365 * teaB.lifetime)

    LCA(system=sysB, lifetime=25, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerB)

    return sysB


# %%

# =============================================================================
# Wrapper function
# =============================================================================

def create_system(system_ID='A', flowsheet=None, ppl=default_ppl):
    ID = system_ID.lower().lstrip('sys').upper()  # so that it'll work for "sysA"/"A"
    reload_lca = False

    # Set flowsheet to avoid stream replacement warnings
    if flowsheet is None:
        flowsheet_ID = f'ng{ID}'
        if hasattr(main_flowsheet.flowsheet, flowsheet_ID): # clear flowsheet
            getattr(main_flowsheet.flowsheet, flowsheet_ID).clear()
            clear_lca_registries()
            reload_lca = True
        flowsheet = Flowsheet(flowsheet_ID)
        main_flowsheet.set_flowsheet(flowsheet)

    _load_components()
    _load_lca_data(reload_lca)

    if system_ID == 'A': f = create_systemA
    elif system_ID == 'B': f = create_systemB
    else: raise ValueError(f'`system_ID` can only be "A" or "B", not "{ID}".')

    try: system = f(flowsheet, ppl=ppl)
    except:
        _load_components(reload=True)
        system = f(flowsheet, ppl=ppl)

    return system