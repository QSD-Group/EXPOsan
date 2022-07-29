#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Lane To <lane20@illinois.edu>

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
from exposan.biogenic_refinery import (
    _load_components,
    _load_lca_data,
    const_daily_wage,
    const_person_days,
    discount_rate,
    get_decay_k,
    get_handcart_and_truck_fee,
    get_ppl,
    get_tanker_truck_fee,
    get_toilet_user,
    max_CH4_emission,
    operator_daily_wage,
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

    create_stream_with_impact_item(stream_ID='polymer', dct_key='Polymer')
    create_stream_with_impact_item(stream_ID='resin', dct_key='Resin')
    create_stream_with_impact_item(stream_ID='filter_bag', dct_key='FilterBag')
    create_stream_with_impact_item(stream_ID='MgOH2')
    create_stream_with_impact_item(stream_ID='MgCO3')
    create_stream_with_impact_item(stream_ID='H2SO4')
    create_stream_with_impact_item(stream_ID='biochar')
    create_stream_with_impact_item(stream_ID='struvite')
    create_stream_with_impact_item(stream_ID='conc_NH3')


def update_toilet_param(unit):
    # Use the private attribute so that the number of users/toilets will be exactly as assigned
    # (i.e., can be fractions)
    unit._N_user = get_toilet_user()
    unit._N_toilet = get_ppl(unit.ID[0])/get_toilet_user()
    unit._run()


def update_carbon_COD_ratio(sys):
    for first_u in sys.units:
        if hasattr(first_u, 'carbon_COD_ratio'):
            carbon_COD_ratio = first_u.carbon_COD_ratio
            break
    for u in sys.units:
        if u is first_u: continue
        if hasattr(u, 'carbon_COD_ratio'): u.carbon_COD_ratio = carbon_COD_ratio


# %%

# =============================================================================
# Scenario A (sysA): pit latrine with 12,000 users
# =============================================================================

def create_systemA(flowsheet=None):
    # Set flowsheet to avoid stream replacement warnings
    flowsheet = flowsheet or main_flowsheet
    streamA = flowsheet.stream
    batch_create_streams('A')

    ##### Human Inputs #####
    A1 = su.Excretion('A1', outs=('urine','feces'))

    ##### User Interface #####
    A2 = su.PitLatrine('A2', ins=(A1-0, A1-1,
                                  'toilet_paper', 'flushing_water',
                                  'cleansing_water', 'desiccant'),
                       outs=('mixed_waste', 'leachate', 'A2_CH4', 'A2_N2O'),
                       N_user=get_toilet_user(), N_toilet=get_ppl('A')/get_toilet_user(),
                       if_flushing=False, if_desiccant=False, if_toilet_paper=False,
                       OPEX_over_CAPEX=0.05, lifetime=10,
                       decay_k_COD=get_decay_k(),
                       decay_k_N=get_decay_k(),
                       max_CH4_emission=max_CH4_emission
                       )
    A2.specification = lambda: update_toilet_param(A2)

    ##### Conveyance of Waste #####
    A3 = su.Trucking('A3', ins=A2-0, outs=('transported', 'conveyance_loss'),
                      load_type='mass', distance=5, distance_unit='km',
                      interval=A2.emptying_period, interval_unit='yr',
                      loss_ratio=0.02)
    def update_A3_param():
        A3._run()
        truck = A3.single_truck
        truck.interval = A2.emptying_period*365*24
        truck.load = A3.F_mass_in*truck.interval/A2.N_toilet
        rho = A3.F_mass_in/A3.F_vol_in
        vol = truck.load/rho
        A3.fee = get_tanker_truck_fee(vol)
        A3._design()
    A3.specification = update_A3_param

    ##### Treatment #####
    A4 = su.BiogenicRefineryControls('A4', ins=A3-0, outs='A4_out')

    A5 = su.BiogenicRefineryHousing('A5', ins=A4-0, outs='A5_out',
                                    const_wage=const_daily_wage,
                                    const_person_days=const_person_days)

    A6 = su.BiogenicRefineryScrewPress('A6', ins=(A5-0, streamA.polymer), outs=('liq', 'cake_sol'))

    A7 = su.LiquidTreatmentBed('A7', ins=A6-0,
                               outs=('liquid_bed_treated', 'A7_CH4', 'A7_N2O'),
                               decay_k_COD=get_decay_k(),
                               decay_k_N=get_decay_k(),
                               max_CH4_emission=max_CH4_emission)

    A8 = su.BiogenicRefineryCarbonizerBase('A8', outs=(streamA.biochar, 'A8_hot_gas', 'A8_N2O'))

    A9 = su.BiogenicRefineryPollutionControl('A9', ins=(A8-1, A8-2), outs=('A9_hot_gas_pcd', 'A9_N2O'))

    # Update uptime_ratio in all units to follow carbonizer base
    A9_old_cost = A9._cost
    def update_A9_uptime_ratio():
        A12.uptime_ratio = A11.uptime_ratio = A10.uptime_ratio = A9.uptime_ratio = A8.uptime_ratio
        A9_old_cost()
    A9._cost = update_A9_uptime_ratio

    A10 = su.BiogenicRefineryOHX('A10', ins=A9-0, outs='A10_hot_gas')
    A11 = su.BiogenicRefineryHHX('A11', ins=A10-0, outs='A11_hot_gas')
    A12 = su.BiogenicRefineryHHXdryer('A12', ins=(A6-1, A11-0), outs=('waste_out', 'A12_N2O', 'A12_CH4'))
    A12-0-A8

    A13 = su.Mixer('A13', ins=(A2-2, A7-1, A12-2), outs=streamA.CH4)
    A13.specification = lambda: add_fugitive_items(A13, 'CH4_item')
    A13.line = 'fugitive CH4 mixer'

    A14 = su.Mixer('A14', ins=(A2-3, A7-2, A9-1, A12-1), outs=streamA.N2O)
    A14.specification = lambda: add_fugitive_items(A14, 'N2O_item')
    A14.line = 'fugitive N2O mixer'

    A15 = su.ComponentSplitter('A15', ins=A7-0,
                               outs=(streamA.liq_N, streamA.liq_P, streamA.liq_K,
                                     'liq_non_fertilizers'),
                               split_keys=(('NH3', 'NonNH3'), 'P', 'K'))
    # Make sure the same carbon_COD_ratio is used throughout
    A15.specification = lambda: update_carbon_COD_ratio(sysA)
    A15.run_after_specification = True

    ##### Simulation, TEA, and LCA #####
    sysA = System('sysA', path=(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15))
    teaA = SimpleTEA(system=sysA, discount_rate=discount_rate,
                      start_year=2020, lifetime=20, uptime_ratio=1,
                      lang_factor=None, annual_maintenance=0,
                      annual_labor=(operator_daily_wage*3*365))

    # 12 is assuming the device is running 12 hr per day (50% of the time)
    # this isn't adjusted through `uptime_ratio` because other OPEX calculation
    # in this unit needs `uptime_ratio` to be 1
    get_powerA = lambda: sum([(u.power_utility.rate*u.uptime_ratio)
                              for u in sysA.units])*(365*teaA.lifetime)*12
    LCA(system=sysA, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerA)

    return sysA


# %%

# =============================================================================
# Scenario B (sysB): urine-diverting dry toilet (UDDT) with 12,000 users
# =============================================================================

def create_systemB(flowsheet=None):
    flowsheet = flowsheet or main_flowsheet
    streamB = flowsheet.stream
    batch_create_streams('B')

    ##### Human Inputs #####
    B1 = su.Excretion('B1', outs=('urine','feces'))

    ##### User Interface #####
    B2 = su.UDDT('B2', ins=(B1-0, B1-1,
                            'toilet_paper', 'flushing_water',
                            'cleaning_water', 'desiccant'),
                 outs=('liq_waste', 'sol_waste',
                       'struvite_scaled', 'HAP', 'B2_CH4', 'B2_N2O'),
                 N_user=get_toilet_user(), N_toilet=get_ppl('B')/get_toilet_user(),
                 if_flushing=False, if_desiccant=False, if_toilet_paper=False,
                 OPEX_over_CAPEX=0.1, lifetime=10,
                 decay_k_COD=get_decay_k(),
                 decay_k_N=get_decay_k(),
                 max_CH4_emission=max_CH4_emission)
    B2.specification = lambda: update_toilet_param(B2)

    ##### Conveyance of Waste #####
    # Liquid waste
    B3 = su.Trucking('B3', ins=B2-0, outs=('transported_l', 'loss_l'),
                     load_type='mass', distance=5, distance_unit='km',
                     loss_ratio=0.02)

    # Solid waste
    B4 = su.Trucking('B4', ins=B2-1, outs=('transported_s', 'loss_s'),
                     load_type='mass', load=1, load_unit='tonne',
                     distance=5, distance_unit='km',
                     loss_ratio=0.02)

    def update_B3_B4_param():
        B4._run()
        truck3, truck4 = B3.single_truck, B4.single_truck
        hr = truck3.interval = truck4.interval = B2.collection_period*24
        ppl_t = get_ppl('B') / B2.N_toilet
        truck3.load = B3.F_mass_in * hr / B2.N_toilet
        truck4.load = B4.F_mass_in * hr / B2.N_toilet
        rho3 = B3.F_mass_in/B3.F_vol_in
        rho4 = B4.F_mass_in/B4.F_vol_in
        B3.fee = get_handcart_and_truck_fee(truck3.load/rho3, ppl_t, True, B2)
        B4.fee = get_handcart_and_truck_fee(truck4.load/rho4, ppl_t, False, B2)
        B3._design()
        B4._design()
    B4.specification = update_B3_B4_param

    ##### Treatment #####
    B5 = su.BiogenicRefineryStruvitePrecipitation(
        'B5', ins=(B3-0, streamB.MgOH2, streamB.MgCO3, streamB.filter_bag),
        outs=('B5_treated', streamB.struvite))

    B6 = su.BiogenicRefineryIonExchange('B6', ins=(B5-0, streamB.resin, streamB.H2SO4),
                                        outs=('B6_treated', 'resin_out', streamB.conc_NH3))

    B7 = su.LiquidTreatmentBed('B7', ins=B6-0, outs=('liquid_bed_treated', 'B7_CH4', 'B7_N2O'),
                               decay_k_COD=get_decay_k(),
                               decay_k_N=get_decay_k(),
                               max_CH4_emission=max_CH4_emission)

    B8 = su.BiogenicRefineryControls('B8', ins=B7-0, outs='B8_out')

    B9 = su.BiogenicRefineryHousing('B9', ins=B8-0, outs='B9_out',
                                    const_wage=const_daily_wage,
                                    const_person_days=const_person_days)

    B10 = su.BiogenicRefineryGrinder('B10', ins=B4-0, outs='waste')

    B11 = su.BiogenicRefineryCarbonizerBase('B11', outs=(streamB.biochar, 'B11_hot_gas', 'B11_N2O'))

    B12 = su.BiogenicRefineryPollutionControl('B12', ins=(B11-1, B11-2), outs=('B12_hot_gas_pcd', 'B12_N2O'))

    # updating uptime_ratio in all units to follow carbonizer base
    B12_old_cost = B12._cost
    def update_B12_uptime_ratio():
        B15.uptime_ratio = B14.uptime_ratio = B13.uptime_ratio = B12.uptime_ratio = B11.uptime_ratio
        B12_old_cost()
    B12._cost = update_B12_uptime_ratio

    B13 = su.BiogenicRefineryOHX('B13', ins=B12-0, outs='B13_hot_gas')

    B14 = su.BiogenicRefineryHHX('B14', ins=B13-0, outs='B14_hot_gas')

    B15 = su.BiogenicRefineryHHXdryer('B15', ins=(B10-0, B14-0), outs=('waste_out', 'B15_N2O', 'B15_CH4'))
    B15-0-B11

    B16 = su.Mixer('B16', ins=(B2-4, B7-1, B15-2), outs=streamB.CH4)
    B16.specification = lambda: add_fugitive_items(B16, 'CH4_item')
    B16.line = 'fugitive CH4 mixer'

    B17 = su.Mixer('B17', ins=(B2-5, B7-2, B12-1, B15-1), outs=streamB.N2O)
    B17.specification = lambda: add_fugitive_items(B17, 'N2O_item')
    B17.line = 'fugitive N2O mixer'

    B18 = su.ComponentSplitter('B18', ins=B9-0,
                               outs=(streamB.liq_N, streamB.liq_P, streamB.liq_K,
                                     'liq_non_fertilizers'),
                               split_keys=(('NH3', 'NonNH3'), 'P', 'K'))
    # Make sure the same carbon_COD_ratio is used throughout
    B18.specification = lambda: update_carbon_COD_ratio(sysB)
    B18.run_after_specification = True

    ##### Simulation, TEA, and LCA #####
    sysB = System('sysB', path=(B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17, B18))
    teaB = SimpleTEA(system=sysB, discount_rate=discount_rate,
                      start_year=2020, lifetime=20, uptime_ratio=1,
                      lang_factor=None, annual_maintenance=0,
                      annual_labor=(operator_daily_wage*3*365))

    # 12 is assuming the device is running 12 hr per day (50% of the time)
    # this isn't adjusted through `uptime_ratio` because other OPEX calculation
    # in this unit needs `uptime_ratio` to be 1
    get_powerB = lambda: sum([(u.power_utility.rate*u.uptime_ratio)
                              for u in sysB.units])*(365*teaB.lifetime)*12
    LCA(system=sysB, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerB)

    return sysB


# %%

# =============================================================================
# Scenario C (sysC): pit latrine with 10,000 users
# =============================================================================

def create_systemC(flowsheet=None):
    flowsheet = flowsheet or main_flowsheet
    streamC = flowsheet.stream
    batch_create_streams('C')

    ##### Human Inputs #####
    C1 = su.Excretion('C1', outs=('urine','feces'))

    ##### User Interface #####
    C2 = su.PitLatrine('C2', ins=(C1-0, C1-1,
                                  'toilet_paper', 'flushing_water',
                                  'cleansing_water', 'desiccant'),
                       outs=('mixed_waste', 'leachate', 'C2_CH4', 'C2_N2O'),
                       N_user=get_toilet_user(), N_toilet=get_ppl('C')/get_toilet_user(),
                       if_flushing=False, if_desiccant=False, if_toilet_paper=False,
                       OPEX_over_CAPEX=0.05, lifetime=10,
                       decay_k_COD=get_decay_k(),
                       decay_k_N=get_decay_k(),
                       max_CH4_emission=max_CH4_emission)
    C2.specification = lambda: update_toilet_param(C2)

    ##### Conveyance of Waste #####
    C3 = su.Trucking('C3', ins=C2-0, outs=('transported', 'conveyance_loss'),
                      load_type='mass', distance=5, distance_unit='km',
                      interval=C2.emptying_period, interval_unit='yr',
                      loss_ratio=0.02)
    def update_C3_param():
        C3._run()
        truck = C3.single_truck
        truck.interval = C2.emptying_period*365*24
        truck.load = C3.F_mass_in*truck.interval/C2.N_toilet
        rho = C3.F_mass_in/C3.F_vol_in
        vol = truck.load/rho
        C3.fee = get_tanker_truck_fee(vol)
        C3._design()
    C3.specification = update_C3_param

    ##### Treatment #####

    C4 = su.BiogenicRefineryControls('C4', ins=C3-0, outs='C4_out')
    C5 = su.BiogenicRefineryHousing('C5', ins=C4-0, outs='C5_out',
                                    const_wage=const_daily_wage,
                                    const_person_days=const_person_days)

    C6 = su.BiogenicRefineryScrewPress('C6', ins=(C5-0, streamC.polymer), outs=('liq', 'cake_sol'))
    C7 = su.LiquidTreatmentBed('C7', ins=C6-0, outs=('liquid_bed_treated', 'C7_CH4', 'C7_N2O'),
                               decay_k_COD=get_decay_k(),
                               decay_k_N=get_decay_k(),
                               max_CH4_emission=max_CH4_emission)

    # Agricultural residue flowrate is proportional to sysA effluent to adjust for
    # the population difference between sysA and sysC
    #!!! Not sure what the 0.95 is, maybe adjusting for moisture content?
    #!!! Need to figure out a way to get the loading without needing A6
    ag_res_ratio = (get_ppl('A')-get_ppl('C'))/get_ppl('A') * (1/0.95)
    ag_res_kg_hr = 13.28 * ag_res_ratio # 13.28 is (A6.outs[1].F_mass-A6.outs[1].imass['H2O'])
    # ag_res_kg_hr = (A6.outs[1].F_mass-A6.outs[1].imass['H2O']) * ag_res_ratio
    rice_husk = WasteStream('rice_husk', RiceHusk=ag_res_kg_hr, units='kg/hr')

    C_ag_mixer = su.Mixer('C_ag_mixer', ins=(C6-1, rice_husk), outs=('C_ag_mixer_out'))

    C8 = su.BiogenicRefineryCarbonizerBase('C8', outs=(streamC.biochar, 'C8_hot_gas', 'C8_N2O'))

    C9 = su.BiogenicRefineryPollutionControl('C9', ins=(C8-1, C8-2), outs=('C9_hot_gas_pcd', 'C9_N2O'))

    # updating uptime_ratio in all units to follow carbonizer base
    C9_old_cost = C9._cost
    def update_C9_uptime_ratio():
        C12.uptime_ratio = C11.uptime_ratio = C10.uptime_ratio = C9.uptime_ratio = C8.uptime_ratio
        C9_old_cost()
    C9._cost = update_C9_uptime_ratio

    C10 = su.BiogenicRefineryOHX('C10', ins = C9-0, outs = ('C10_hot_gas'))
    C11 = su.BiogenicRefineryHHX('C11', ins = C10-0, outs = ('C11_hot_gas'))
    C12 = su.BiogenicRefineryHHXdryer('C12', ins = (C_ag_mixer-0, C11-0), outs = ('waste_out', 'C12_N2O', 'C12_CH4'))
    C12-0-C8

    C13 = su.Mixer('C13', ins=(C2-2, C7-1, C12-2), outs=streamC.CH4)
    C13.specification = lambda: add_fugitive_items(C13, 'CH4_item')
    C13.line = 'fugitive CH4 mixer'

    C14 = su.Mixer('C14', ins=(C2-3, C7-2, C9-1, C12-1), outs=streamC.N2O)
    C14.specification = lambda: add_fugitive_items(C14, 'N2O_item')
    C14.line = 'fugitive N2O mixer'

    ##### Simulation, TEA, and LCA #####
    sysC = System('sysC', path=(C1, C2, C3, C4, C5, C6, C7, C_ag_mixer, C8, C9, C10, C11, C12, C13, C14))

    teaC = SimpleTEA(system=sysC, discount_rate=discount_rate,
                     start_year=2020, lifetime=20, uptime_ratio=1,
                     lang_factor=None, annual_maintenance=0,
                     annual_labor=(operator_daily_wage*3*365))

    # 12 is assuming the device is running 12 hr per day (50% of the time)
    # this isn't adjusted through `uptime_ratio` because other OPEX calculation
    # in this unit needs `uptime_ratio` to be 1
    get_powerC = lambda: sum([(u.power_utility.rate*u.uptime_ratio)
                              for u in sysC.units])*(365*teaC.lifetime)*12
    LCA(system=sysC, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerC)

    return sysC


# %%

# =============================================================================
# Scenario D (sysD): pit latrine with 12,000 users without the biogenic refinery
# =============================================================================

def create_systemD(flowsheet=None):
    flowsheet = flowsheet or main_flowsheet
    streamD = flowsheet.stream
    batch_create_streams('D')

    ##### Human Inputs #####
    D1 = su.Excretion('D1', outs=('urine','feces'))

    ##### User Interface #####
    D2 = su.PitLatrine('D2', ins=(D1-0, D1-1,
                                  'toilet_paper', 'flushing_water',
                                  'cleansing_water', 'desiccant'),
                        outs=('mixed_waste', 'leachate', 'D2_CH4', 'D2_N2O'),
                        N_user=get_toilet_user(), N_toilet=get_ppl('D')/get_toilet_user(),
                        if_flushing=False, if_desiccant=False, if_toilet_paper=False,
                        OPEX_over_CAPEX=0.05, lifetime=5,
                        decay_k_COD=get_decay_k(),
                        decay_k_N=get_decay_k(),
                        max_CH4_emission=max_CH4_emission
                        )
    D2.specification = lambda: update_toilet_param(D2)

    ##### Conveyance of Waste #####
    D3 = su.Trucking('D3', ins=D2-0, outs=('transported', 'conveyance_loss'),
                      load_type='mass', distance=5, distance_unit='km',
                      interval=D2.emptying_period, interval_unit='yr',
                      loss_ratio=0.02)
    def update_D3_param():
        D3._run()
        truck = D3.single_truck
        truck.interval = D2.emptying_period*365*24
        truck.load = D3.F_mass_in*truck.interval/D2.N_toilet
        rho = D3.F_mass_in/D3.F_vol_in
        vol = truck.load/rho
        D3.fee = get_tanker_truck_fee(vol)
        D3._design()
    D3.specification = update_D3_param

    D4 = su.Lagoon('D4', ins=D3-0, outs=('anaerobic_treated', 'D4_CH4', 'D4_N2O'),
                    design_type='anaerobic',
                    decay_k_N=get_decay_k(),
                    max_CH4_emission=max_CH4_emission)

    D5 = su.DryingBed('D5', ins=D4-0, outs=('dried_sludge', 'evaporated',
                                            'D5_CH4', 'D5_N2O'),
                      design_type='unplanted',
                      decay_k_COD=get_decay_k(),
                      decay_k_N=get_decay_k(),
                      max_CH4_emission=max_CH4_emission)

    D6 = su.Mixer('D6', ins=(D2-2, D4-1, D5-2), outs=streamD.CH4)
    D6.specification = lambda: add_fugitive_items(D5, 'CH4_item')
    D6.line = 'fugitive CH4 mixer'

    D7 = su.Mixer('D7', ins=(D2-3, D4-2, D5-3), outs=streamD.N2O)
    D7.specification = lambda: add_fugitive_items(D6, 'N2O_item')
    D7.line = 'fugitive N2O mixer'

    sysD = System('sysD', path=(D1, D2, D3, D4, D5, D6, D7))
    sysD.simulate()

    teaD = SimpleTEA(system=sysD, discount_rate=discount_rate,
                     start_year=2020, lifetime=20, uptime_ratio=1,
                     lang_factor=None, annual_maintenance=0)

    # 12 is assuming the device is running 12 hr per day (50% of the time)
    # this isn't adjusted through `uptime_ratio` because other OPEX calculation
    # in this unit needs `uptime_ratio` to be 1
    get_powerD = lambda: sum([(u.power_utility.rate*u.uptime_ratio)
                              for u in sysD.units])*(365*teaD.lifetime)*12
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
        flowsheet_ID = f'br{ID}'
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