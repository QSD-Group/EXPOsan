#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
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
from exposan.utils import add_fugitive_items, get_generic_tanker_truck_fee as get_tanker_truck_fee
from exposan.eco_san import (
    _load_components,
    _load_lca_data,
    discount_rate,
    get_decay_k,
    get_toilet_user,
    max_CH4_emission,
    operator_daily_wage,
    ppl,
    price_dct,
    )

__all__ = ('create_system',)


# %%

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
    create_stream_with_impact_item(stream_ID='struvite')
    create_stream_with_impact_item(stream_ID='salt')
    create_stream_with_impact_item(stream_ID='HCl_acid')


# %%

# =============================================================================
# Scenario A (sysA)
# original design tested in Durban
# anaerobic-aerobic-anoxic-aerobic-electrochemical reactor (A-O-A-O-ECR)
# =============================================================================

def create_systemA(flowsheet=None):
    flowsheet = flowsheet or main_flowsheet
    streamA = flowsheet.stream
    batch_create_streams('A')

    ##### Human Inputs #####
    A1 = su.Excretion('A1', outs=('urine','feces'))

    ##### User Interface #####
    A2 = su.MURT('A2',
                 ins=(A1-0, A1-1, 'toilet_paper', 'flushing_water', 'cleansing_water', 'desiccant'),
                 outs=('mixed_waste', 'A2_CH4', 'A2_N2O'),
                 decay_k_COD=get_decay_k(),
                 decay_k_N=get_decay_k(),
                 max_CH4_emission=max_CH4_emission,
                 N_user=300/7, N_toilet=7,
                 if_flushing=True, if_desiccant=False, if_toilet_paper=True,
                 CAPEX=0, OPEX_over_CAPEX=0.07)

    ##### Treatment #####
    A3 = su.EcoSanPrimary('A3', ins=(A2-0, streamA['MgOH2']), # MgOH2 is empty
                          outs=('A3_treated', 'A3_CH4', 'A3_N2O', 'sludge'),
                          decay_k_COD=get_decay_k(),
                          decay_k_N=get_decay_k(),
                          max_CH4_emission=max_CH4_emission,
                          if_with_MBR=True)

    A4 = su.EcoSanAnaerobic('A4', ins=A3-0, outs=('A4_treated', 'A4_CH4', 'A4_N2O'),
                            decay_k_COD=get_decay_k(),
                            decay_k_N=get_decay_k(),
                            max_CH4_emission=max_CH4_emission)
    A5 = su.EcoSanAerobic('A5', ins=A4-0, outs = ('A5_treated','A5_CH4', 'A5_N2O'),
                          decay_k_COD=get_decay_k(),
                          decay_k_N=get_decay_k(),
                          max_CH4_emission=max_CH4_emission)
    A6 = su.EcoSanAnoxic('A6', ins=A5-0, outs = ('A6_treated', 'A6_CH4', 'A6_N2O'),
                         decay_k_COD=get_decay_k(),
                         decay_k_N=get_decay_k(),
                         max_CH4_emission=max_CH4_emission) #double check
    A7 = su.EcoSanAerobic('A7', ins=A6-0, outs = ('A7_treated','A7_CH4', 'A7_N2O'),
                          decay_k_COD=get_decay_k(),
                          decay_k_N=get_decay_k(),
                          max_CH4_emission=max_CH4_emission)
    A8 = su.EcoSanECR('A8', ins=(A7-0, streamA['salt'], streamA['HCl_acid']),
                      outs = ('A8_treated'), if_after_MBR=False,
                      decay_k_COD=get_decay_k(),
                      decay_k_N=get_decay_k(),
                      max_CH4_emission=max_CH4_emission)

    A10 = su.Mixer('A10', ins=(A3-1, A2-1, A4-1, A5-1, A6-1, A7-1), outs=streamA['CH4'])
    A10.specification = lambda: add_fugitive_items(A10, 'CH4_item')
    A10.line = 'fugitive CH4 mixer'

    A11 = su.Mixer('A11', ins=(A3-2, A2-2, A4-2, A5-2, A6-2, A7-2), outs=streamA['N2O'])
    A11.specification = lambda: add_fugitive_items(A11, 'N2O_item')
    A11.line = 'fugitive N2O mixer'

    ##### Other impacts and costs #####
    A13 = su.EcoSanBioCost('A13', ins=A8-0, outs=('A13_out',))
    A14 = su.EcoSanSystem('A14', ins=A13-0, outs=('A14_out',))
    A15 = su.EcoSanSolar('A15', ins=A14-0, outs=('A15_out',))

    A16 = su.Trucking('A16', ins=A3-3, outs=('transported', 'conveyance_loss'),
                      load_type='mass', distance=5, distance_unit='km',
                      interval=365, interval_unit='d',
                      loss_ratio=0.02)
    def update_A16_param():
        A16._run()
        truck = A16.single_truck
        truck.interval = 365*24
        truck.load = A16.F_mass_in*truck.interval
        rho = A16.F_mass_in/A16.F_vol_in
        vol = truck.load/rho
        A16.fee = get_tanker_truck_fee(vol)
        A16._design()
    A16.specification = update_A16_param

    ##### Simulation, TEA, and LCA #####
    sysA = System('sysA', path=(A1, A2, A3, A4, A5, A6, A7, A8, A10, A11, A13, A14, A15, A16))

    teaA = SimpleTEA(system=sysA, discount_rate=discount_rate,
                      start_year=2020, lifetime=10, uptime_ratio=1,
                      lang_factor=None, annual_maintenance=0,
                      #!!! If multiplying by 52, shouldn't it be the weekily wage?
                      annual_labor=operator_daily_wage*52)

    get_powerA = lambda: sum([u.power_utility.rate for u in sysA.units]) * (24 * 365 * teaA.lifetime)
    LCA(system=sysA, lifetime=10, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerA)

    return sysA


# %%

# =============================================================================
# Scenario B (sysB)
# anaerobic-aerobic-membrane bioreactor (MBR)-ECR
# (anoxic and facultative aerobic with MBR)
# =============================================================================

def create_systemB(flowsheet=None):
    flowsheet = flowsheet or main_flowsheet
    streamB = flowsheet.stream
    batch_create_streams('B')

    ##### Human Inputs #####
    B1 = su.Excretion('B1', outs=('urine','feces'))

    ##### User Interface #####
    B2 = su.MURT('B2',
                 ins=(B1-0, B1-1,'toilet_paper', 'flushing_water', 'cleansing_water', 'desiccant'),
                 outs=('mixed_waste', 'B2_CH4', 'B2_N2O'),
                 decay_k_COD=get_decay_k(),
                 decay_k_N=get_decay_k(),
                 max_CH4_emission=max_CH4_emission,
                 N_user=300/7, N_toilet=7,
                 if_flushing=True, if_desiccant=False, if_toilet_paper=True,
                 CAPEX=0, OPEX_over_CAPEX=0.07)

    ##### Treatment #####
    B3 = su.EcoSanPrimary('B3', ins=(B2-0, streamB['MgOH2']), # MgOH2 is empty
                          outs=('B3_treated', 'B3_CH4', 'B3_N2O', 'sludge'),
                          decay_k_COD=get_decay_k(),
                          decay_k_N=get_decay_k(),
                          max_CH4_emission=max_CH4_emission,
                          if_with_MBR=False)

    B4 = su.EcoSanAnaerobic('B4', ins=B3-0, outs=('B4_treated', 'B4_CH4', 'B4_N2O'),
                        decay_k_COD=get_decay_k(),
                        decay_k_N=get_decay_k(),
                        max_CH4_emission=max_CH4_emission)
    B5 = su.EcoSanAerobic('B5', ins=B4-0, outs=('B5_treated','B5_CH4', 'B5_N2O'),
                        decay_k_COD=get_decay_k(),
                        decay_k_N=get_decay_k(),
                        max_CH4_emission=max_CH4_emission)

    B6 = su.EcoSanMBR('B6', ins=B5-0, outs=('B6_treated', 'B6_CH4', 'B6_N2O'),
                      decay_k_COD=get_decay_k(),
                      decay_k_N=get_decay_k(),
                      max_CH4_emission=max_CH4_emission)

    B8 = su.EcoSanECR('B8', ins=(B6-0, streamB['salt'], streamB['HCl_acid']),
                      outs=('B8_treated',), if_after_MBR=True,
                      decay_k_COD=get_decay_k(),
                      decay_k_N=get_decay_k(),
                      max_CH4_emission=max_CH4_emission)


    B10 = su.Mixer('B10', ins=(B2-1, B3-1, B4-1, B5-1, B6-1), outs=streamB['CH4'])
    B10.specification = lambda: add_fugitive_items(B10, 'CH4_item')
    B10.line = 'fugitive CH4 mixer'

    B11 = su.Mixer('B11', ins=(B2-2, B3-2, B4-2, B5-2, B6-2), outs=streamB['N2O'])
    B11.specification = lambda: add_fugitive_items(B11, 'N2O_item')
    B11.line = 'fugitive N2O mixer'

    ##### Other impacts and costs #####
    B14 = su.EcoSanSystem('B14', ins=(B8-0), outs = ('B14_out'))
    B15 = su.EcoSanSolar('B15', ins=(B14-0), outs = ('B15_out'))


    B16 = su.Trucking('B16', ins=B3-3, outs=('transported', 'conveyance_loss'),
                      load_type='mass', distance=5, distance_unit='km',
                      interval=365, interval_unit='d',
                      loss_ratio=0.02)
    def update_B16_param():
        B16._run()
        truck = B16.single_truck
        truck.interval = 365*24
        truck.load = B16.F_mass_in*truck.interval/7
        rho = B16.F_mass_in/B16.F_vol_in
        vol = truck.load/rho
        B16.fee = get_tanker_truck_fee(vol)
        B16._design()
    B16.specification = update_B16_param

    ##### Simulation, TEA, and LCA #####
    sysB = System('sysB', path=(B1, B2, B3, B4, B5, B6, B8, B10, B11, B14, B15, B16))

    teaB = SimpleTEA(system=sysB, discount_rate=discount_rate,
                      start_year=2020, lifetime=10, uptime_ratio=1,
                      lang_factor=None, annual_maintenance=0,
                      #!!! 52?
                      annual_labor=operator_daily_wage*12)

    get_powerB = lambda: sum([u.power_utility.rate for u in sysB.units]) * (24 * 365 * teaB.lifetime)
    LCA(system=sysB, lifetime=10, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerB)

    return sysB


# %%

# =============================================================================
# Scenario C (sysC)
# sysB with struvite precipitation
# MBR + struvite precipitation (Anoxic and Facultative Aerobic with MBR)
# =============================================================================

def create_systemC(flowsheet=None):
    flowsheet = flowsheet or main_flowsheet
    streamC = flowsheet.stream
    batch_create_streams('C')

    ##### Human Inputs #####
    C1 = su.Excretion('C1', outs=('urine','feces'))

    ##### User Interface #####
    C2 = su.MURT('C2',
                 ins=(C1-0, C1-1, 'toilet_paper', 'flushing_water', 'cleansing_water', 'desiccant'),
                 outs=('mixed_waste', 'C2_CH4', 'C2_N2O'),
                 decay_k_COD=get_decay_k(),
                 decay_k_N=get_decay_k(),
                 max_CH4_emission=max_CH4_emission,
                 N_user=get_toilet_user(), N_toilet=ppl/get_toilet_user(),
                 if_flushing=True, if_desiccant=False, if_toilet_paper=True,
                 CAPEX=0, OPEX_over_CAPEX=0.07)

    ##### Treatment #####
    C3 = su.EcoSanPrimary('C3', ins=(C2-0, streamC['MgOH2']),
                          outs=('C3_treated', 'C3_CH4', 'C3_N2O', 'sludge', streamC['struvite']),
                          decay_k_COD=get_decay_k(),
                          decay_k_N=get_decay_k(),
                          max_CH4_emission=max_CH4_emission,
                          if_with_MBR=True)
    C4 = su.EcoSanAnaerobic('C4', ins=C3-0, outs=('C4_treated', 'C4_CH4', 'C4_N2O'),
                            decay_k_COD=get_decay_k(),
                            decay_k_N=get_decay_k(),
                            max_CH4_emission=max_CH4_emission)
    C5 = su.EcoSanAerobic('C5', ins=C4-0, outs=('C5_treated','C5_CH4', 'C5_N2O'),
                          decay_k_COD=get_decay_k(),
                          decay_k_N=get_decay_k(),
                          max_CH4_emission=max_CH4_emission)
    C6 = su.EcoSanMBR('C6', ins=C5-0, outs=('C6_treated', 'C6_CH4', 'C6_N2O'),
                      decay_k_COD=get_decay_k(),
                      decay_k_N=get_decay_k(),
                      max_CH4_emission=max_CH4_emission)
    C8 = su.EcoSanECR('C8', ins=(C6-0, streamC['salt'], streamC['HCl_acid']),
                      outs = ('C8_treated'), if_after_MBR=True,
                      decay_k_COD=get_decay_k(),
                      decay_k_N=get_decay_k(),
                      max_CH4_emission=max_CH4_emission)

    C10 = su.Mixer('C10', ins=(C2-1, C3-1, C4-1, C5-1, C6-1), outs=streamC['CH4'])
    C10.specification = lambda: add_fugitive_items(C10, 'CH4_item')
    C10.line = 'fugitive CH4 mixer'

    C11 = su.Mixer('C11', ins=(C2-2, C3-2, C4-2, C5-2, C6-2), outs=streamC['N2O'])
    C11.specification = lambda: add_fugitive_items(C1, 'N2O_item')
    C11.line = 'fugitive N2O mixer'

    ##### Other impacts and costs #####
    C14 = su.EcoSanSystem('C14', ins=(C8-0), outs = ('C14_out'))

    C16 = su.Trucking('C16', ins=C3-3, outs=('transported', 'conveyance_loss'),
                      load_type='mass', distance=5, distance_unit='km',
                      interval=365, interval_unit='d',
                      loss_ratio=0.02)
    def update_C16_param():
        C16._run()
        truck = C16.single_truck
        truck.interval = 365*24
        truck.load = C16.F_mass_in*truck.interval/C2.N_toilet
        rho = C16.F_mass_in/C16.F_vol_in
        vol = truck.load/rho
        C16.fee = get_tanker_truck_fee(vol)
        C16._design()
    C16.specification = update_C16_param

    ##### Simulation, TEA, and LCA #####
    sysC = System('sysC', path=(C1, C2, C3, C4, C5, C6, C8, C10, C11, C14, C16))

    teaC = SimpleTEA(system=sysC, discount_rate=discount_rate,
                      start_year=2020, lifetime=10, uptime_ratio=1,
                      lang_factor=None, annual_maintenance=0,
                      #!!! 52?
                      annual_labor=operator_daily_wage*12)

    get_powerC = lambda: sum([u.power_utility.rate for u in sysC.units]) * (24 * 365 * teaC.lifetime)

    LCA(system=sysC, lifetime=10, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerC)

    return sysC


# %%

# =============================================================================
# Wrapper function
# =============================================================================

def create_system(system_ID='A', flowsheet=None):
    ID = system_ID.lower().lstrip('sys').upper()  # so that it'll work for "sysA"/"A"
    reload_lca = False

    # Set flowsheet to avoid stream replacement warnings
    if flowsheet is None:
        flowsheet_ID = f'es{ID}'
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
    else: raise ValueError(f'`system_ID` can only be "A" , "B", or "C", not "{ID}".')

    try: system = f(flowsheet)
    except:
        _load_components(reload=True)
        system = f(flowsheet)

    return system