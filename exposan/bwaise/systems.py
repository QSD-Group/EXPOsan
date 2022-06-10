#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


# %%

# Filter out warnings related to solid content
import warnings
warnings.filterwarnings('ignore', message='Solid content')

import numpy as np
from sklearn.linear_model import LinearRegression as LR
from qsdsan import (
    Flowsheet, main_flowsheet,
    WasteStream,
    sanunits as su,
    ImpactItem,
    System, SimpleTEA, LCA,
    )
from exposan.bwaise._process_settings import price_dct

# =============================================================================
# Unit parameters
# =============================================================================

household_size = 4
household_per_toilet = 4
get_toilet_user = lambda: household_size * household_per_toilet

# Number of people served by the existing plant (sysA and sysC)
ppl_exist_sewer = 4e4
ppl_exist_sludge = 416667
# Number of people served by the alternative plant (sysB)
ppl_alt = 5e4
def get_ppl(kind):
    if kind.lower() in ('exist', 'existing', 'sysa', 'sysc', 'a', 'c'):
        return ppl_exist_sewer+ppl_exist_sludge
    elif kind.lower() in ('alt', 'alternative', 'sysb', 'b'):
        return ppl_alt
    else:
        raise ValueError('`kind` should be "exist" (for sysA and sysC)'
                         f'or "alt" for sysB, not {kind}.')

exchange_rate = 3700 # UGX per USD
discount_rate = 0.05

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3
# Get reduction rate constant k for COD and N, use a function so that k can be
# changed during uncertainty analysis
def get_decay_k(tau_deg=2, log_deg=3):
    k = (-1/tau_deg)*np.log(10**-log_deg)
    return k

max_CH4_emission = 0.25

# Model for tanker truck cost based on capacity (m3)
# price = a*capacity**b -> ln(price) = ln(a) + bln(capacity)
UGX_price_dct = np.array((8e4, 12e4, 20e4, 25e4))
capacities = np.array((3, 4.5, 8, 15))
emptying_fee = 0.15 # additional emptying fee, fraction of base cost
def get_tanker_truck_fee(capacity):
    price_dct = UGX_price_dct*(1+emptying_fee)/exchange_rate
    ln_p = np.log(price_dct)
    ln_cap = np.log(capacities)
    model = LR().fit(ln_cap.reshape(-1,1), ln_p.reshape(-1,1))
    predicted = model.predict(np.array((np.log(capacity))).reshape(1, -1)).item()
    cost = np.exp(predicted)
    return cost

# Flow rates for treatment plants
sewer_flow = 2750 # m3/d
sludge_flow_exist = 500 # m3/d
sludge_flow_alt = 60 # m3/d
get_sludge_flow = lambda kind: \
    sludge_flow_exist if kind.lower() in ('exist', 'sysa', 'sysc', 'a', 'c') else sludge_flow_alt

# Nutrient loss during application
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05




# =============================================================================
# Universal units and functions
# =============================================================================

def batch_create_streams(prefix):
    stream_dct = {}
    item = ImpactItem.get_item('CH4_item').copy(f'{prefix}_CH4_item', set_as_source=True)
    stream_dct['CH4'] = WasteStream(f'{prefix}_CH4', phase='g', stream_impact_item=item)

    item = ImpactItem.get_item('N2O_item').copy(f'{prefix}_N2O_item', set_as_source=True)
    stream_dct['N2O'] = WasteStream(f'{prefix}_N2O', phase='g', stream_impact_item=item)

    item = ImpactItem.get_item('N_item').copy(f'{prefix}_liq_N_item', set_as_source=True)
    stream_dct['liq_N'] = WasteStream(f'{prefix}_liq_N', phase='l', price=price_dct['N'],
                                      stream_impact_item=item)

    item = ImpactItem.get_item('N_item').copy(f'{prefix}_sol_N_item', set_as_source=True)
    stream_dct['sol_N'] = WasteStream(f'{prefix}_sol_N', phase='l', price=price_dct['N'],
                                              stream_impact_item=item)

    item = ImpactItem.get_item('P_item').copy(f'{prefix}_liq_P_item', set_as_source=True)
    stream_dct['liq_P'] = WasteStream(f'{prefix}_liq_P', phase='l', price=price_dct['P'],
                                      stream_impact_item=item)

    item = ImpactItem.get_item('P_item').copy(f'{prefix}_sol_P_item', set_as_source=True)
    stream_dct['sol_P'] = WasteStream(f'{prefix}_sol_P', phase='l', price=price_dct['P'],
                                      stream_impact_item=item)

    item = ImpactItem.get_item('K_item').copy(f'{prefix}_liq_K_item', set_as_source=True)
    stream_dct['liq_K'] = WasteStream(f'{prefix}_liq_K', phase='l', price=price_dct['K'],
                                      stream_impact_item=item)

    item = ImpactItem.get_item('K_item').copy(f'{prefix}_sol_K_item', set_as_source=True)
    stream_dct['sol_K'] = WasteStream(f'{prefix}_sol_K', phase='l', price=price_dct['K'],
                                      stream_impact_item=item)
    return stream_dct

def update_toilet_param(unit, kind):
    # Use the private attribute so that the number of users/toilets will be exactly as assigned
    # (i.e., can be fractions)
    unit._N_user = get_toilet_user()
    unit._N_toilet = get_ppl(kind)/get_toilet_user()
    unit._run()

def update_lagoon_flow_rate(unit):
    unit.flow_rate = sewer_flow + get_sludge_flow('exist')
    unit._run()

def add_fugitive_items(unit, item_ID):
    unit._run()
    for i in unit.ins:
        i.stream_impact_item = ImpactItem.get_item(item_ID).copy(set_as_source=True)

# Costs of WWTP units have been considered in the lumped unit
def clear_unit_costs(sys):
    for i in sys.units:
        if isinstance(i, su.LumpedCost):
            continue
        i.purchase_costs.clear()
        i.installed_costs.clear()


def adjust_NH3_loss(unit):
    unit._run()
    # Assume the slight higher loss of NH3 does not affect COD,
    # does not matter much since COD not considered in crop application
    unit.outs[0]._COD = unit.outs[1]._COD = unit.ins[0]._COD


# %%

# =============================================================================
# Scenario A (sysA): pit latrine with existing treatment system
# =============================================================================

def create_systemA(streamsA=None):
    streamsA = streamsA or batch_create_streams('A')

    ##### Human Inputs #####
    A1 = su.Excretion('A1', outs=('urine', 'feces'))

    ##### User Interface #####
    A2 = su.PitLatrine('A2', ins=(A1-0, A1-1,
                                  'toilet_paper', 'flushing_water',
                                  'cleansing_water', 'desiccant'),
                       outs=('mixed_waste', 'leachate', 'A2_CH4', 'A2_N2O'),
                       N_user=get_toilet_user(), N_toilet=get_ppl('exist')/get_toilet_user(),
                       OPEX_over_CAPEX=0.05,
                       decay_k_COD=get_decay_k(tau_deg, log_deg),
                       decay_k_N=get_decay_k(tau_deg, log_deg),
                       max_CH4_emission=max_CH4_emission
                       )
    A2.specification = lambda: update_toilet_param(A2, 'exist')

    ##### Conveyance #####
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
    A4 = su.LumpedCost('A4', ins=A3-0, cost_item_name='Lumped WWTP',
                       CAPEX=18606700, power=57120/(365*24), lifetime=8)
    A4.line = 'Lumped WWTP cost'
    get_A4_lifetime = lambda: A4.lifetime
    def update_A_t0():
        A5.t0 = A2.emptying_period
        A6.t0 = A5.t0 # A5.tau is for the solids
        A7.t0 = A6.t0 + A6.tau/365
        A8.t0 = A5.t0 + A5.tau/365
    A4.specification = update_A_t0
    A4.run_after_specification = True

    A5 = su.Sedimentation('A5', ins=A4-0,
                          outs=('liq', 'sol', 'A5_CH4', 'A5_N2O'),
                          decay_k_COD=get_decay_k(tau_deg, log_deg),
                          decay_k_N=get_decay_k(tau_deg, log_deg),
                          max_CH4_emission=max_CH4_emission)

    A6 = su.Lagoon('A6', ins=A5-0, outs=('anaerobic_treated', 'A6_CH4', 'A6_N2O'),
                   design_type='anaerobic',
                   flow_rate=sewer_flow+get_sludge_flow('exist'),
                   decay_k_N=get_decay_k(tau_deg, log_deg),
                   max_CH4_emission=max_CH4_emission)
    A6.specification = lambda: update_lagoon_flow_rate(A6)

    A7 = su.Lagoon('A7', ins=A6-0, outs=('facultative_treated', 'A7_CH4', 'A7_N2O'),
                   design_type='facultative',
                   flow_rate=sewer_flow+get_sludge_flow('exist'),
                   decay_k_N=get_decay_k(tau_deg, log_deg),
                   max_CH4_emission=max_CH4_emission,
                   if_N2O_emission=True)
    A7.specification = lambda: update_lagoon_flow_rate(A7)

    A8 = su.DryingBed('A8', ins=A5-1, outs=('dried_sludge', 'evaporated',
                                            'A8_CH4', 'A8_N2O'),
                      design_type='unplanted',
                      decay_k_COD=get_decay_k(tau_deg, log_deg),
                      decay_k_N=get_decay_k(tau_deg, log_deg),
                      max_CH4_emission=max_CH4_emission)
    treatA = System('treatA', path=(A4, A5, A6, A7, A8))
    A8._cost = lambda: clear_unit_costs(treatA)

    ##### Reuse or Disposal #####
    A9 = su.CropApplication('A9', ins=A7-0, outs=('liquid_fertilizer', 'reuse_loss'),
                            loss_ratio=app_loss)
    A9.specification = lambda: adjust_NH3_loss(A9)

    A10 = su.Mixer('A10', ins=(A2-2, A5-2, A6-1, A7-1, A8-2), outs=streamsA['CH4'])
    A10.specification = lambda: add_fugitive_items(A10, 'CH4_item')
    A10.line = 'fugitive CH4 mixer'

    A11 = su.Mixer('A11', ins=(A2-3, A5-3, A6-2, A7-2, A8-3), outs=streamsA['N2O'])
    A11.specification = lambda: add_fugitive_items(A11, 'N2O_item')
    A11.line = 'fugitive N2O mixer'

    A12 = su.ComponentSplitter('A12', ins=A8-0,
                               outs=(streamsA['sol_N'], streamsA['sol_P'], streamsA['sol_K'],
                                     'A_sol_non_fertilizers'),
                               split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

    A13 = su.ComponentSplitter('A13', ins=A9-0,
                               outs=(streamsA['liq_N'], streamsA['liq_P'], streamsA['liq_K'],
                                     'A_liq_non_fertilizers'),
                               split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

    ##### Simulation, TEA, and LCA #####
    sysA = System('sysA', path=(A1, A2, A3, treatA, A9, A10, A11, A12, A13))

    exist_staff_num = 12
    get_annual_labor = lambda: exist_staff_num*3e6*12/exchange_rate
    SimpleTEA(system=sysA, discount_rate=discount_rate, start_year=2018,
              lifetime=get_A4_lifetime(), uptime_ratio=1, lang_factor=None,
              annual_maintenance=0, annual_labor=get_annual_labor())

    LCA(system=sysA, lifetime=get_A4_lifetime(), lifetime_unit='yr', uptime_ratio=1,
        annualize_construction=True,
        # Assuming all additional WWTP OPEX from electricity
        E_item=lambda: A4.power_utility.rate*(365*24)*get_A4_lifetime())

    return sysA


# %%

# =============================================================================
# Scenario B (sysB): pit latrine with anaerobic treatment
# =============================================================================

def create_systemB(streamsB=None):
    streamsB = streamsB or batch_create_streams('B')
    B_biogas_item = ImpactItem.get_item('Biogas_item').copy('B_biogas_item', set_as_source=True)
    streamsB['biogas'] = WasteStream('B_biogas', phase='g', price=price_dct['Biogas'],
                                     stream_impact_item=B_biogas_item)

    ##### Human Inputs #####
    B1 = su.Excretion('B1', outs=('urine', 'feces'))

    ##### User Interface #####
    B2 = su.PitLatrine('B2', ins=(B1-0, B1-1,
                                  'toilet_paper', 'flushing_water',
                                  'cleansing_water', 'desiccant'),
                       outs=('mixed_waste', 'leachate', 'B2_CH4', 'B2_N2O'),
                       N_user=get_toilet_user(), N_toilet=get_ppl('alt')/get_toilet_user(),
                       OPEX_over_CAPEX=0.05,
                       decay_k_COD=get_decay_k(tau_deg, log_deg),
                       decay_k_N=get_decay_k(tau_deg, log_deg),
                       max_CH4_emission=max_CH4_emission)
    B2.specification = lambda: update_toilet_param(B2, 'alt')

    ##### Conveyance #####
    B3 = su.Trucking('B3', ins=B2-0, outs=('transported', 'conveyance_loss'),
                     load_type='mass', distance=5, distance_unit='km',
                     interval=B2.emptying_period, interval_unit='yr',
                     loss_ratio=0.02)
    def update_B3_param():
        B3._run()
        truck = B3.single_truck
        truck.interval = B2.emptying_period*365*24
        truck.load = B3.F_mass_in*truck.interval/B2.N_toilet
        rho = B3.F_mass_in/B3.F_vol_in
        vol = truck.load/rho
        B3.fee = get_tanker_truck_fee(vol)
        B3._design()
    B3.specification = update_B3_param

    ##### Treatment #####
    B4 = su.LumpedCost('B4', ins=B3-0, cost_item_name='Lumped WWTP',
                       CAPEX=337140, power=6854/(365*24), lifetime=10)
    B4.line = 'Lumped WWTP cost'
    get_B4_lifetime = lambda: B4.lifetime
    def update_B_t0():
        B5.t0 = B2.emptying_period
        B8.t0 = B7.t0 = B5.t0 + B5.tau/365
    B4.specification = update_B_t0
    B4.run_after_specification = True

    B5 = su.AnaerobicBaffledReactor('B5', ins=B4-0, outs=('ABR_treated', 'biogas',
                                                          'B5_CH4', 'B5_N2O'),
                                    decay_k_COD=get_decay_k(tau_deg, log_deg),
                                    max_CH4_emission=max_CH4_emission)
    def update_B5_gravel_density():
        B5.gravel_density = B2.density_dct['Gravel']
    B5.specification = update_B5_gravel_density
    B5.run_after_specification = True

    B6 = su.SludgeSeparator('B6', ins=B5-0, outs=('liq', 'sol'))

    B7 = su.LiquidTreatmentBed('B7', ins=B6-0, outs=('liquid_bed_treated', 'B7_CH4', 'B7_N2O'),
                               decay_k_COD=get_decay_k(tau_deg, log_deg),
                               decay_k_N=get_decay_k(tau_deg, log_deg),
                               max_CH4_emission=max_CH4_emission)

    B8 = su.DryingBed('B8', ins=B6-1, outs=('dried_sludge', 'evaporated',
                                            'B8_CH4', 'B8_N2O'),
                      design_type='planted',
                      decay_k_COD=get_decay_k(tau_deg, log_deg),
                      decay_k_N=get_decay_k(tau_deg, log_deg),
                      max_CH4_emission=max_CH4_emission)

    treatB = System('treatB', path=(B4, B5, B6, B7, B8))
    B8._cost = lambda: clear_unit_costs(treatB)

    ##### Reuse or Disposal #####
    B9 = su.CropApplication('B9', ins=B7-0, outs=('liquid_fertilizer', 'reuse_loss'),
                            loss_ratio=app_loss)
    B9.specification = lambda: adjust_NH3_loss(B9)

    B10 = su.Mixer('B10', ins=(B2-2, B5-2, B7-1, B8-2), outs=streamsB['CH4'])
    B10.specification = lambda: add_fugitive_items(B10, 'CH4_item')
    B10.line = 'fugitive CH4 mixer'

    B11 = su.Mixer('B11', ins=(B2-3, B5-3, B7-2, B8-3), outs=streamsB['N2O'])
    B11.specification = lambda: add_fugitive_items(B11, 'N2O_item')
    B11.line = 'fugitive N2O mixer'

    B12 = su.ComponentSplitter('B12', ins=B8-0,
                                outs=(streamsB['sol_N'], streamsB['sol_P'], streamsB['sol_K'],
                                      'B_sol_non_fertilizers'),
                                split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

    B13 = su.ComponentSplitter('B13', ins=B9-0,
                                outs=(streamsB['liq_N'], streamsB['liq_P'], streamsB['liq_K'],
                                      'B_liq_non_fertilizers'),
                                split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

    B14 = su.BiogasCombustion('B14', ins=(B5-1, 'air'),
                              outs=('used', 'lost', 'wasted'),
                              if_combustion=False,
                              biogas_loss=0.1, biogas_eff=0.55)
    B15 = su.Mixer('B15', ins=(B14-0, B14-2), outs=streamsB['biogas'])

    ##### Simulation, TEA, and LCA #####
    sysB = System('sysB', path=(B1, B2, B3, treatB, B9, B10, B11, B12, B13, B14, B15))

    skilled_num = 5
    get_skilled_num = lambda: skilled_num
    skilled_salary = 5e6 # UGX/month
    get_tot_skilled_salary = lambda: skilled_salary*get_skilled_num()

    unskilled_num = 5
    get_unskilled_num = lambda: unskilled_num
    unskilled_salary = 75e4 # UGX/month
    get_tot_unskilled_salary = lambda: unskilled_salary*get_unskilled_num()
    get_alt_salary = lambda: (get_tot_skilled_salary()+get_tot_unskilled_salary())*12/exchange_rate

    SimpleTEA(system=sysB, discount_rate=discount_rate, start_year=2018,
              lifetime=get_B4_lifetime(), uptime_ratio=1, lang_factor=None,
              annual_maintenance=0, annual_labor=get_alt_salary())

    LCA(system=sysB, lifetime=get_B4_lifetime(), lifetime_unit='yr', uptime_ratio=1,
        annualize_construction=True,
        # Assuming all additional WWTP OPEX from electricity
        E_item=lambda: B4.power_utility.rate*(365*24)*get_B4_lifetime())

    return sysB


# %%

# =============================================================================
# Scenario C (sysC): containaer-based sanitation with existing treatment system
# =============================================================================

flowsheetC = Flowsheet('sysC')
main_flowsheet.set_flowsheet(flowsheetC)
streamsC = batch_create_streams('C')

def create_systemC(streamsC=None):
    streamsC = streamsC or batch_create_streams('C')

    ##### Human Inputs #####
    C1 = su.Excretion('C1', outs=('urine', 'feces'))

    ##### User Interface #####
    C2 = su.UDDT('C2', ins=(C1-0, C1-1,
                            'toilet_paper', 'flushing_water',
                            'cleaning_water', 'desiccant'),
                 outs=('liq_waste', 'sol_waste',
                       'struvite', 'HAP', 'C2_CH4', 'C2_N2O'),
                 N_user=get_toilet_user(), N_toilet=get_ppl('exist')/get_toilet_user(),
                 OPEX_over_CAPEX=0.1,
                 decay_k_COD=get_decay_k(tau_deg, log_deg),
                 decay_k_N=get_decay_k(tau_deg, log_deg),
                 max_CH4_emission=max_CH4_emission)
    C2.specification = lambda: update_toilet_param(C2, 'exist')

    ##### Conveyance #####
    handcart_fee = 0.01 # USD/cap/d
    truck_fee = 23e3 # UGX/m3

    # Handcart fee is for both liquid/solid
    get_handcart_and_truck_fee = \
        lambda vol, ppl, include_fee: truck_fee/exchange_rate*vol \
            + int(include_fee)*handcart_fee*ppl*C2.collection_period

    # Liquid waste
    C3 = su.Trucking('C3', ins=C2-0, outs=('transported_l', 'loss_l'),
                     load_type='mass', distance=5, distance_unit='km',
                     loss_ratio=0.02)

    # Solid waste
    C4 = su.Trucking('C4', ins=C2-1, outs=('transported_s', 'loss_s'),
                     load_type='mass', load=1, load_unit='tonne',
                     distance=5, distance_unit='km',
                     loss_ratio=0.02)
    def update_C3_C4_param():
        C4._run()
        truck3, truck4 = C3.single_truck, C4.single_truck
        hr = truck3.interval = truck4.interval = C2.collection_period*24
        N_toilet = C2.N_toilet
        ppl = get_ppl('exist') / N_toilet
        truck3.load = C3.F_mass_in * hr / N_toilet
        truck4.load = C4.F_mass_in * hr / N_toilet
        rho3 = C3.F_mass_in/C3.F_vol_in
        rho4 = C4.F_mass_in/C4.F_vol_in
        C3.fee = get_handcart_and_truck_fee(truck3.load/rho3, ppl, True)
        C4.fee = get_handcart_and_truck_fee(truck4.load/rho4, ppl, False)
        C3._design()
        C4._design()
    C4.specification = update_C3_C4_param

    ##### Treatment #####
    C5 = su.LumpedCost('C5', ins=(C3-0, C4-0),
                       cost_item_name='Lumped WWTP',
                       CAPEX=18606700, power=57120/(365*24), lifetime=8)
    get_C5_lifetime = lambda: C5.lifetime
    def update_C_t0():
        C8.t0 = C6.t0 = C2.collection_period / 365
        C7.t0 = C6.t0 + C6.tau/365
    C5.specification = update_C_t0
    C5.run_after_specification = True

    C6 = su.Lagoon('C6', ins=C5-0, outs=('anaerobic_treated', 'C6_CH4', 'C6_N2O'),
                   design_type='anaerobic',
                   flow_rate=sewer_flow+get_sludge_flow('exist'),
                   decay_k_N=get_decay_k(tau_deg, log_deg),
                   max_CH4_emission=max_CH4_emission)
    C6.specification = lambda: update_lagoon_flow_rate(C6)

    C7 = su.Lagoon('C7', ins=C6-0, outs=('facultative_treated', 'C7_CH4', 'C7_N2O'),
                   design_type='facultative',
                   flow_rate=sewer_flow+get_sludge_flow('exist'),
                   decay_k_N=get_decay_k(tau_deg, log_deg),
                   max_CH4_emission=max_CH4_emission,
                   if_N2O_emission=True)
    C7.specification = lambda: update_lagoon_flow_rate(C7)

    C8 = su.DryingBed('C8', ins=C5-1, outs=('dried_sludge', 'evaporated', 'C8_CH4', 'C8_N2O'),
                     design_type='unplanted',
                     decay_k_COD=get_decay_k(tau_deg, log_deg),
                     decay_k_N=get_decay_k(tau_deg, log_deg),
                     max_CH4_emission=max_CH4_emission)
    treatC = System('treatC', path=(C5, C6, C7, C8))
    C8._cost = lambda: clear_unit_costs(treatC)

    ##### Reuse or Disposal #####
    C9 = su.CropApplication('C9', ins=C7-0, outs=('liq_fertilizer', 'liq_loss'),
                            loss_ratio=app_loss)
    C9.specification = lambda: adjust_NH3_loss(C9)

    C10 = su.CropApplication('C10', ins=C8-0, outs=('sol_fertilizer', 'sol_loss'),
                            loss_ratio=app_loss)
    def adjust_C10_loss_ratio():
        C10.loss_ratio.update(C9.loss_ratio)
        adjust_NH3_loss(C10)
    C10.specification = adjust_C10_loss_ratio

    C11 = su.Mixer('C11', ins=(C2-4, C6-1, C7-1, C8-2), outs=streamsC['CH4'])
    C11.specification = lambda: add_fugitive_items(C11, 'CH4_item')
    C11.line = 'fugitive CH4 mixer'

    C12 = su.Mixer('C12', ins=(C2-5, C6-2, C7-2, C8-3), outs=streamsC['N2O'])
    C12.specification = lambda: add_fugitive_items(C12, 'N2O_item')
    C12.line = 'fugitive N2O mixer'

    C13 = su.ComponentSplitter('C13', ins=C10-0,
                               outs=(streamsC['sol_N'], streamsC['sol_P'], streamsC['sol_K'],
                                     'C_sol_non_fertilizers'),
                               split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

    C14 = su.ComponentSplitter('C14', ins=C9-0,
                               outs=(streamsC['liq_N'], streamsC['liq_P'], streamsC['liq_K'],
                                     'C_liq_non_fertilizers'),
                               split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

    ##### Simulation, TEA, and LCA #####
    sysC = System('sysC', path=(C1, C2, C3, C4, treatC, C9, C10, C11, C12, C13, C14))

    SimpleTEA(system=sysC, discount_rate=discount_rate, start_year=2018,
              lifetime=get_C5_lifetime(), uptime_ratio=1, lang_factor=None,
              annual_maintenance=0, annual_labor=12*3e6*12/exchange_rate)

    LCA(system=sysC, lifetime=get_C5_lifetime(), lifetime_unit='yr', uptime_ratio=1,
        annualize_construction=True,
        # Assuming all additional WWTP OPEX from electricity
        E_item=lambda: C5.power_utility.rate*(365*24)*get_C5_lifetime())

    return sysC


# %%

# =============================================================================
# Wrapper function
# =============================================================================

def create_systeam(system_ID='A', flowsheet=None):
    system_ID = system_ID.lower().lstrip('sys').upper() # so that it'll work for "sysA"/"A"


    # Set flowsheet to avoid stream replacement warnings
    if flowsheet is None:
        flowsheet = Flowsheet(f'sys{system_ID}')
        main_flowsheet.set_flowsheet(flowsheet)

    streams = batch_create_streams(system_ID)

    if system_ID == 'A': system = create_systemA(streams)
    elif system_ID == 'B': system = create_systemB(streams)
    elif system_ID == 'C': system = create_systemC(streams)
    else: raise ValueError(f'`system_ID` can only be "A", "B", or "C", not "{system_ID}".')
    return system