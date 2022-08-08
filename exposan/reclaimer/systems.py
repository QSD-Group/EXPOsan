#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
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
from exposan.reclaimer import (
    _load_components,
    _load_lca_data,
    discount_rate,
    get_decay_k,
    max_CH4_emission,
    ppl,
    price_dct,
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

    create_stream_with_impact_item(stream_ID='MgOH2')
    create_stream_with_impact_item(stream_ID='KCl')
    create_stream_with_impact_item(stream_ID='GAC')
    create_stream_with_impact_item(stream_ID='zeolite')
    create_stream_with_impact_item(stream_ID='conc_NH3', dct_key='conc_NH3')
    create_stream_with_impact_item(stream_ID='LPG')


# %%

# =============================================================================
# Scenario A (sysA): Duke Reclaimer Design 2.0 Solids Treatment Only
# =============================================================================

def create_systemA(flowsheet=None):
    # Set flowsheet to avoid stream replacement warnings
    flowsheet = flowsheet or main_flowsheet
    streamA = flowsheet.stream
    batch_create_streams('A')

    ##### Human Inputs #####
    A1 = su.Excretion('A1', outs=('urine', 'feces'))

    ##### Treatment #####

    A3 = su.SepticTank('A3', ins=(A1-0, streamA['MgOH2']),
                       outs=('A3_treated', 'A3_CH4', 'A3_N2O', 'A3_sludge'),
                       decay_k_COD=get_decay_k(),
                       decay_k_N=get_decay_k(),
                       max_CH4_emission=max_CH4_emission,
                       ppl=ppl,
                       if_include_front_end=True,
                       if_generate_struvite=True,
                       if_struvite_in_sludge=True,
                       )

    # Biogas-related streams are empty as no biogas is generated in Reclaimer
    A4 = su.SludgePasteurization('A4', ins=('biogas', 'air', A3-3, streamA['LPG']),
                                 outs=('biogas_used', 'biogas_lost', 'biogas_wasted', 'treated_sludge'),
                                 heat_loss=0.1, target_MC=0.1, sludge_temp=283.15,
                                 temp_pasteurization=343.15, lhv_lpg=48.5,
                                 ppl=ppl, baseline_ppl=100,
                                 user_scale_up=None, exponent_scale=0.6,
                                 if_sludge_service=True)

    A8 = su.Mixer('A8', ins=(A3-1), outs=streamA['CH4'])
    A8.specification = lambda: add_fugitive_items(A8, 'CH4_item')
    A8.line = 'fugitive CH4 mixer'

    A9 = su.Mixer('A9', ins=(A3-2), outs=streamA['N2O'])
    A9.specification = lambda: add_fugitive_items(A9, 'N2O_item')
    A9.line = 'fugitive N2O mixer'

    ##### Simulation, TEA, and LCA #####
    sysA = System('sysA', path=(A1, A3, A4, A8, A9))

    teaA = SimpleTEA(system=sysA, discount_rate=discount_rate,
                     start_year=2020, lifetime=20, uptime_ratio=1,
                     lang_factor=None, annual_maintenance=0,
                     annual_labor=0)

    get_powerA = lambda: sum([u.power_utility.rate for u in sysA.units]) * (24 * 365 * teaA.lifetime)

    LCA(system=sysA, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerA)

    return sysA


# %%

# =============================================================================
# System B Duke Reclaimer 2.0 Full System on grid-tied energy
# =============================================================================

def create_systemB(flowsheet=None):
    flowsheet = flowsheet or main_flowsheet
    streamB = flowsheet.stream
    batch_create_streams('B')

    ##### Human Inputs #####
    B1 = su.Excretion('B1', outs=('urine','feces'))

    ##### User Interface #####
    # Reclaimer 2.0 can process ~30L/hr(net), 720L/24 hours of constant operation
    # flush volume of 6L per flush determines 120 flushes in one day (~20 users if they use toilet 6 times)
    # CAPEX of MURT is considered in the housing unit
    B2 = su.MURT('B2', ins=(B1-0, B1-1, 'toilet_paper', 'flushing_water', 'cleansing_water', 'desiccant'),
                 outs=('mixed_waste', 'B2_CH4', 'B2_N2O'),
                 decay_k_COD=get_decay_k(),
                 decay_k_N=get_decay_k(),
                 max_CH4_emission=max_CH4_emission,
                 N_user=ppl, N_tot_user=ppl, lifetime=8,
                 if_flushing=True, if_desiccant=False, if_toilet_paper=True,
                 CAPEX=500*max(1, ppl/100), OPEX_over_CAPEX=0.06)

    ##### Treatment #####
    B3 = su.SepticTank('B3', ins=(B2-0, streamB['MgOH2']),
                       outs=('B3_treated', 'B3_CH4', 'B3_N2O', 'B3_sludge'),
                       decay_k_COD=get_decay_k(),
                       decay_k_N=get_decay_k(),
                       max_CH4_emission=max_CH4_emission,
                       ppl=ppl,
                       if_include_front_end=True,
                       if_generate_struvite=True,
                       if_struvite_in_sludge=True,
                       )

    # Biogas-related streams are empty as no biogas is generated in Reclaimer
    B4 = su.SludgePasteurization('B4', ins=('biogas', 'air', B3-3, streamB['LPG']),
                                 outs=('biogas_used', 'biogas_lost', 'biogas_wasted', 'treated_sludge'),
                                 heat_loss=0.1, target_MC=0.1, sludge_temp=283.15,
                                 temp_pasteurization=343.15, lhv_lpg=48.5,
                                 ppl=ppl, baseline_ppl=100,
                                 user_scale_up=None, exponent_scale=0.6,
                                 if_sludge_service=True)

    B5 = su.ReclaimerUltrafiltration('B5', ins=(B3-0),
                                     outs=('B5_treated', 'retentate'),
                                     ppl=ppl, if_gridtied=True)

    B6 = su.ReclaimerIonExchange('B6', ins=(B5-0, streamB['zeolite'], streamB['GAC'], streamB['KCl']),
                                 outs=('B6_treated', 'spent_zeolite', 'spent_GAC', streamB['conc_NH3']),
                                 ppl=ppl)

    B7 = su.ReclaimerECR('B7', ins=(B6-0), outs='B7_treated',
                        ppl=ppl, if_gridtied=True)

    B8 = su.Mixer('B8', ins=(B2-1, B3-1), outs=streamB['CH4'])
    B8.specification = lambda: add_fugitive_items(B7, 'CH4_item')
    B8.line = 'fugitive CH4 mixer'

    B9 = su.Mixer('B9', ins=(B2-2, B3-2), outs=streamB['N2O'])
    B9.specification = lambda: add_fugitive_items(B8, 'N2O_item')
    B9.line = 'fugitive N2O mixer'

    ##### Non-Treatment #####
    B10 = su.ReclaimerHousing('B10', ins=(B7-0), outs='B10_out', ppl=ppl)
    B11 = su.ReclaimerSystem('B11', ins=(B10-0), outs='B11_out', ppl=ppl, if_gridtied=True)

    B12 = su.ComponentSplitter('B12', ins=B4-3,
                               outs=(streamB['sol_N'], streamB['sol_P'], streamB['sol_K'],
                                     'B_sol_non_fertilizers'),
                               split_keys=(('NH3', 'NonNH3'), ('P', 'Struvite'), 'K'))

    ##### Simulation, TEA, and LCA #####
    sysB = System('sysB', path=(B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12))

    teaB = SimpleTEA(system=sysB, discount_rate=discount_rate,
                     start_year=2020, lifetime=20, uptime_ratio=1,
                     lang_factor=None, annual_maintenance=0,
                     annual_labor=0)

    get_powerB = lambda: sum([u.power_utility.rate for u in sysB.units]) * (24 * 365 * teaB.lifetime)

    LCA(system=sysB, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerB)

    return sysB


# %%

# =============================================================================
# System C Duke Reclaimer 2.0 coupled with photovoltaic system
# =============================================================================

def create_systemC(flowsheet=None):
    flowsheet = flowsheet or main_flowsheet
    streamC = flowsheet.stream
    batch_create_streams('C')

    ##### Human Inputs #####
    C1 = su.Excretion('C1', outs=('urine','feces'))

    ##### User Interface #####
    # Reclaimer 2.0 can process ~30L/hr(net), 720L/24 hours of constant operation
    # flush volume of 6L per flush determines 120 flushes in one day (~20 users if they use toilet 6 times)
    C2 = su.MURT('C2',
                 ins=(C1-0, C1-1, 'toilet_paper', 'flushing_water', 'cleansing_water', 'desiccant'),
                 outs=('mixed_waste', 'C2_CH4', 'C2_N2O'),
                 decay_k_COD=get_decay_k(),
                 decay_k_N=get_decay_k(),
                 max_CH4_emission=max_CH4_emission,
                 N_user=ppl, N_tot_user=ppl, lifetime=8,
                 if_flushing=True, if_desiccant=False, if_toilet_paper=True,
                 CAPEX=500*max(1, ppl/100), OPEX_over_CAPEX=0.06)

    ##### Treatment #####
    C3 = su.SepticTank('C3', ins=(C2-0, streamC['MgOH2']),
                       outs=('C3_treated', 'C3_CH4', 'C3_N2O', 'C3_sludge'),
                       decay_k_COD=get_decay_k(),
                       decay_k_N=get_decay_k(),
                       max_CH4_emission=max_CH4_emission,
                       ppl=ppl,
                       if_include_front_end=True,
                       if_generate_struvite=True,
                       if_struvite_in_sludge=True,
                       )

    # Biogas-related streams are empty as no biogas is generated in Reclaimer
    C4 = su.SludgePasteurization('C4', ins=('biogas', 'air', C3-3, streamC['LPG']),
                                 outs=('biogas_used', 'biogas_lost', 'biogas_wasted', 'treated_sludge'),
                                 heat_loss=0.1, target_MC=0.1, sludge_temp=283.15,
                                 temp_pasteurization=343.15, lhv_lpg=48.5,
                                 ppl=ppl, baseline_ppl=100,
                                 user_scale_up=None, exponent_scale=0.6,
                                 if_sludge_service=True)

    C5 = su.ReclaimerUltrafiltration('C5', ins=(C3-0),
                                     outs=('C5_treated', 'retentate'),
                                     ppl=ppl,
                                     if_gridtied=False)

    C6 = su.ReclaimerIonExchange('C6', ins=(C5-0, streamC['zeolite'], streamC['GAC'], streamC['KCl']),
                                 outs=('C6_treated', 'spent_zeolite', 'spent_GAC',streamC['conc_NH3']),
                                 ppl=ppl)

    C7 = su.ReclaimerECR('C7', ins=(C6-0), outs='C7_treated',
                        ppl=ppl, if_gridtied=False)

    C8 = su.Mixer('C8', ins=(C3-1, C2-1), outs=streamC['CH4'])
    C8.specification = lambda: add_fugitive_items(C7, 'CH4_item')
    C8.line = 'fugitive CH4 mixer'

    C9 = su.Mixer('B9', ins=(C3-2, C2-2), outs=streamC['N2O'])
    C9.specification = lambda: add_fugitive_items(C8, 'N2O_item')
    C9.line = 'fugitive N2O mixer'

    ##### Non-Treatment #####
    C10 = su.ReclaimerHousing('C10', ins=(C7-0), outs='C10_out', ppl=ppl)
    C11 = su.ReclaimerSystem('C11', ins=(C10-0), outs='C11_out', ppl=ppl, if_gridtied=False)
    C12 = su.ReclaimerSolar('C12', ins=(C11-0), outs='C12_out')

    C13 = su.ComponentSplitter('C13', ins=C4-3,
                               outs=(streamC['sol_N'], streamC['sol_P'], streamC['sol_K'],
                                     'C_sol_non_fertilizers'),
                               split_keys=(('NH3', 'NonNH3'), ('P', 'Struvite'), 'K'))

    ##### Simulation, TEA, and LCA #####
    sysC = System('sysC', path=(C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13))

    teaC = SimpleTEA(system=sysC, discount_rate=discount_rate,
                     start_year=2020, lifetime=20, uptime_ratio=1,
                     lang_factor=None, annual_maintenance=0,
                     annual_labor=0)

    get_powerC = lambda: sum([u.power_utility.rate for u in sysC.units]) * (24 * 365 * teaC.lifetime)

    LCA(system=sysC, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerC)

    return sysC


# %%

# =============================================================================
# Scenario D (sysD): Targeted Nitrogen Removal Only
# =============================================================================

def create_systemD(flowsheet=None):
    flowsheet = flowsheet or main_flowsheet
    streamD = flowsheet.stream
    batch_create_streams('D')

    ##### Human Inputs #####
    D1 = su.Excretion('D1', outs=('urine','feces'))

    ##### Treatment #####
    D3 = su.SepticTank('D3', ins=(D1-0, streamD['MgOH2']),
                       outs=('D3_treated', 'D3_CH4', 'D3_N2O', 'D3_sludge'),
                       decay_k_COD=get_decay_k(),
                       decay_k_N=get_decay_k(),
                       max_CH4_emission=max_CH4_emission,
                       ppl=ppl,
                       if_include_front_end=True,
                       if_generate_struvite=True,
                       if_struvite_in_sludge=True,
                       )

    D4 = su.ReclaimerUltrafiltration('D4', ins=(D3-0),
                                     outs=('D4_treated', 'retentate'),
                                     ppl=ppl,
                                     if_gridtied=True)

    D5 = su.ReclaimerIonExchange('D5', ins=(D4-0, streamD['zeolite'], streamD['GAC'], streamD['KCl']),
                                 outs=('D5_treated', 'spent_zeolite', 'spent_GAC',streamD['conc_NH3']),
                                 ppl=ppl)

    D6 = su.Mixer('D6', ins=(D3-1), outs=streamD['CH4'])
    D6.specification = lambda: add_fugitive_items(D6, 'CH4_item')
    D6.line = 'fugitive CH4 mixer'

    D7 = su.Mixer('D7', ins=(D3-2), outs=streamD['N2O'])
    D7.specification = lambda: add_fugitive_items(D7, 'N2O_item')
    D7.line = 'fugitive N2O mixer'

    ##### Non-Treatment #####
    D8 = su.ReclaimerHousing('D8', ins=(D5-0), outs='D8_out', ppl=ppl)
    D9 = su.ReclaimerSystem('D9', ins=(D8-0), outs='D9_out', ppl=ppl, if_gridtied=True)

    ##### Simulation, TEA, and LCA #####
    sysD = System('sysD', path=(D1, D3, D4, D5, D6, D7, D8, D9))

    teaD = SimpleTEA(system=sysD, discount_rate=discount_rate,
                     start_year=2020, lifetime=20, uptime_ratio=1,
                     lang_factor=None, annual_maintenance=0,
                     annual_labor=0)

    get_powerD = lambda: sum([u.power_utility.rate for u in sysD.units]) * (24 * 365 * teaD.lifetime)

    LCA(system=sysD, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerD)

    return sysD


# %%

# =============================================================================
# Wrapper function
# =============================================================================

def create_system(system_ID='A', flowsheet=None):
    ID = system_ID.lower().lstrip('sys').upper() # so that it'll work for "sysA"/"A"
    reload_lca = False

    # Set flowsheet to avoid stream replacement warnings
    if flowsheet is None:
        flowsheet_ID = f're{ID}'
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
    elif system_ID == 'C': f = create_systemC
    elif system_ID == 'D': f = create_systemD
    else: raise ValueError(f'`system_ID` can only be "A", "B", "C", or "D", not "{ID}".')

    try: system = f(flowsheet)
    except:
        _load_components(reload=True)
        system = f(flowsheet)

    return system