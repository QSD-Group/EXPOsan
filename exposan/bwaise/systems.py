#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import numpy as np
import biosteam as bst
import qsdsan as qs
from collections.abc import Iterable
from sklearn.linear_model import LinearRegression as LR
from qsdsan import sanunits as su
from qsdsan import WasteStream, ImpactIndicator, ImpactItem, StreamImpactItem, SimpleTEA, LCA
from exposan.bwaise._cmps import cmps
from exposan.bwaise._lca_data import lca_data_kind, load_lca_data, _ImpactItem_LOADED

# =============================================================================
# Unit parameters
# =============================================================================

currency = qs.currency = 'USD'
qs.CEPCI = qs.CEPCI_by_year[2018]
qs.set_thermo(cmps)
bst.speed_up()

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
get_exchange_rate = lambda: exchange_rate

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
    price_dct = UGX_price_dct*(1+emptying_fee)/get_exchange_rate()
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

# Nutrient loss during applciation
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05

# Energetic content of the biogas
biogas_energy = 803 # kJ/mol CH4
LPG_energy = 50 # MJ/kg
get_biogas_factor = lambda: biogas_energy/cmps.CH4.MW/LPG_energy


# =============================================================================
# Prices and GWP CFs
# =============================================================================

# Recycled nutrients are sold at a lower price than commercial fertilizers
price_factor = 0.25
get_price_factor = lambda: price_factor

price_dct = {
    'Electricity': 0.17,
    'Concrete': 194,
    'Steel': 2.665,
    'N': 1.507*get_price_factor(),
    'P': 3.983*get_price_factor(),
    'K': 1.333*get_price_factor(),
    'Biogas': 6500/get_exchange_rate()*get_biogas_factor()
    }

GWP_dct = {
    'Electricity': 0.15,
    'CH4': 28,
    'N2O': 265,
    'N': -5.4,
    'P': -4.9,
    'K': -1.5,
    'Biogas': -3*get_biogas_factor()
    }

if not ImpactItem.get_item('Excavation'): # prevent from reloading
    load_lca_data('original')

GWP = ImpactIndicator.get_indicator('GWP')

bst.PowerUtility.price = price_dct['Electricity']
ImpactItem.get_item('Concrete').price = price_dct['Concrete']
ImpactItem.get_item('Steel').price = price_dct['Steel']

# =============================================================================
# Universal units and functions
# =============================================================================

def batch_create_stream_items(kind):
    if kind == 'original':
        for k, v in GWP_dct.items():
            if k == 'Electricity':
                ImpactItem(ID='E_item', functional_unit='kWh', GWP=v)
            else:
                StreamImpactItem(ID=f'{k}_item', GWP=v)
    elif kind == 'new':
        E_factor = {'GW2ECO': 0.000000000532, # global warming to ecosystem
                    'GW2HH': 0.0000000812, # global warming to human health
                    'OD2HH': 0.00134} # stratospheric ozone depletion to human health
        H_factor = {'GW2ECO': 0.0000000028, 'GW2HH': 0.000000928, 'OD2HH': 0.000531}
        I_factor = {'GW2ECO': 0.000000025, 'GW2HH': 0.0000125, 'OD2HH': 0.000237}

        StreamImpactItem(ID='CH4_item',
                         E_EcosystemQuality_Total=E_factor['GW2ECO']*4.8,
                         E_HumanHealth_Total=E_factor['GW2HH']*4.8,
                         H_EcosystemQuality_Total=H_factor['GW2ECO']*34,
                         H_HumanHealth_Total=H_factor['GW2HH']*34,
                         I_EcosystemQuality_Total=I_factor['GW2ECO']*84,
                         I_HumanHealth_Total=I_factor['GW2HH']*84
                         )
        StreamImpactItem(ID='N2O_item',
                         E_EcosystemQuality_Total=E_factor['GW2ECO']*78.8,
                         # From climate change + ozone depletion
                         E_HumanHealth_Total=\
                             E_factor['GW2HH']*78.8+E_factor['OD2HH']*0.017,
                         H_EcosystemQuality_Total=H_factor['GW2ECO']*298,
                         H_HumanHealth_Total=\
                             H_factor['GW2HH']*298+H_factor['OD2HH']*0.011,
                         I_EcosystemQuality_Total=I_factor['GW2ECO']*264,
                         I_HumanHealth_Total=\
                             I_factor['GW2HH']*264+I_factor['OD2HH']*0.007
                         )
    else:
        raise ValueError(f'`kind` can only be "original" or "new", not "{kind}".')

    global _ImpactItem_LOADED
    _ImpactItem_LOADED = True


def batch_create_streams(prefix):
    if not _ImpactItem_LOADED:
        batch_create_stream_items(kind=lca_data_kind)

    stream_dct = {}
    item = ImpactItem.get_item('CH4_item').copy(f'{prefix}_CH4_item', set_as_source=True)
    stream_dct['CH4'] = WasteStream(f'{prefix}_CH4', phase='g', stream_impact_item=item)
    # CH4.stream_impact_item = ImpactItem.get_item('CH4_item').copy(stream=CH4, set_as_source=True)

    item = ImpactItem.get_item('N2O_item').copy(f'{prefix}_N2O_item', set_as_source=True)
    stream_dct['N2O'] = WasteStream(f'{prefix}_N2O', phase='g', stream_impact_item=item)
    # N2O.stream_impact_item = ImpactItem.get_item('N2O_item').copy(stream=N2O, set_as_source=True)

    item = ImpactItem.get_item('N_item').copy(f'{prefix}_liq_N_item', set_as_source=True)
    stream_dct['liq_N'] = WasteStream(f'{prefix}_liq_N', phase='l', price=price_dct['N'],
                                      stream_impact_item=item)
    # liq_N.stream_impact_item = ImpactItem.get_item('N_item').copy(stream=liq_N, set_as_source=True)

    item = ImpactItem.get_item('N_item').copy(f'{prefix}_sol_N_item', set_as_source=True)
    stream_dct['sol_N'] = WasteStream(f'{prefix}_sol_N', phase='l', price=price_dct['N'],
                                              stream_impact_item=item)
    # sol_N.stream_impact_item = ImpactItem.get_item('N_item').copy(stream=sol_N, set_as_source=True)

    item = ImpactItem.get_item('P_item').copy(f'{prefix}_liq_P_item', set_as_source=True)
    stream_dct['liq_P'] = WasteStream(f'{prefix}_liq_P', phase='l', price=price_dct['P'],
                                      stream_impact_item=item)
    # liq_P.stream_impact_item=ImpactItem.get_item('P_item').copy(stream=liq_P, set_as_source=True)

    item = ImpactItem.get_item('P_item').copy(f'{prefix}_sol_P_item', set_as_source=True)
    stream_dct['sol_P'] = WasteStream(f'{prefix}_sol_P', phase='l', price=price_dct['P'],
                                      stream_impact_item=item)
    # sol_P.stream_impact_item=ImpactItem.get_item('P_item').copy(stream=sol_P, set_as_source=True)

    item = ImpactItem.get_item('K_item').copy(f'{prefix}_liq_K_item', set_as_source=True)
    stream_dct['liq_K'] = WasteStream(f'{prefix}_liq_K', phase='l', price=price_dct['K'],
                                      stream_impact_item=item)
    # liq_K.stream_impact_item = ImpactItem.get_item('K_item').copy(stream=liq_K, set_as_source=True)

    item = ImpactItem.get_item('K_item').copy(f'{prefix}_sol_K_item', set_as_source=True)
    stream_dct['sol_K'] = WasteStream(f'{prefix}_sol_K', phase='l', price=price_dct['K'],
                                      stream_impact_item=item)
    return stream_dct

def update_toilet_param(unit, kind):
    unit.N_user = get_toilet_user()
    unit.N_toilet = get_ppl(kind)/get_toilet_user()
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

flowsheetA = bst.Flowsheet('sysA')
bst.main_flowsheet.set_flowsheet(flowsheetA)
streamsA = batch_create_streams('A')

#################### Human Inputs ####################
A1 = su.Excretion('A1', outs=('urine', 'feces'))

################### User Interface ###################
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

##################### Conveyance #####################
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


###################### Treatment ######################
A4 = su.LumpedCost('A4', ins=A3-0, cost_item_name='Lumped WWTP',
                   CAPEX=18606700, power=57120/(365*24), lifetime=8)
A4.line = 'Lumped WWTP cost'
get_A4_lifetime = lambda: A4.lifetime

A5 = su.SedimentationTank('A5', ins=A4-0,
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
treatA = bst.System('treatA', path=(A4, A5, A6, A7, A8))
A8._cost = lambda: clear_unit_costs(treatA)

################## Reuse or Disposal ##################
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

############### Simulation, TEA, and LCA ###############
sysA = bst.System('sysA', path=(A1, A2, A3, treatA, A9, A10, A11, A12, A13))

exist_staff_num = 12
get_annual_labor = lambda: exist_staff_num*3e6*12/get_exchange_rate()
teaA = SimpleTEA(system=sysA, discount_rate=discount_rate, start_year=2018,
                 lifetime=get_A4_lifetime(), uptime_ratio=1, lang_factor=None,
                 annual_maintenance=0, annual_labor=get_annual_labor(),
                 construction_schedule=None)

lcaA = LCA(system=sysA, lifetime=get_A4_lifetime(), lifetime_unit='yr', uptime_ratio=1,
            # Assuming all additional WWTP OPEX from electricity
            E_item=lambda: A4.power_utility.rate*(365*24)*8)


# %%

# =============================================================================
# Scenario B (sysB): pit latrine with anaerobic treatment
# =============================================================================

flowsheetB = bst.Flowsheet('sysB')
bst.main_flowsheet.set_flowsheet(flowsheetB)
streamsB = batch_create_streams('B')
B_biogas_item = ImpactItem.get_item('Biogas_item').copy('B_biogas_item', set_as_source=True)
streamsB['biogas'] = WasteStream('B_biogas', phase='g', price=price_dct['Biogas'],
                                 stream_impact_item=B_biogas_item)

#################### Human Inputs ####################
B1 = su.Excretion('B1', outs=('urine', 'feces'))

################### User Interface ###################
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

##################### Conveyance #####################
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

###################### Treatment ######################
B4 = su.LumpedCost('B4', ins=B3-0, cost_item_name='Lumped WWTP',
                   CAPEX=337140, power=6854/(365*24), lifetime=10)
B4.line = 'Lumped WWTP cost'
get_B4_lifetime = lambda: B4.lifetime

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

treatB = bst.System('treatB', path=(B4, B5, B6, B7, B8))
B8._cost = lambda: clear_unit_costs(treatB)

################## Reuse or Disposal ##################
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

############### Simulation, TEA, and LCA ###############
sysB = bst.System('sysB', path=(B1, B2, B3, treatB, B9, B10, B11, B12, B13, B14, B15))

skilled_num = 5
unskilled_num = 5
get_unskilled_num = lambda: unskilled_num
unskilled_salary = 75e4 # UGX/month
get_unskilled_salary = lambda: unskilled_salary*get_unskilled_num()
get_alt_salary = lambda: (skilled_num*5e6+get_unskilled_salary())*12/get_exchange_rate()

teaB = SimpleTEA(system=sysB, discount_rate=discount_rate, start_year=2018,
                  lifetime=get_B4_lifetime(), uptime_ratio=1, lang_factor=None,
                  annual_maintenance=0, annual_labor=get_alt_salary(),
                  construction_schedule=None)

lcaB = LCA(system=sysB, lifetime=get_B4_lifetime(), lifetime_unit='yr', uptime_ratio=1,
           # Assuming all additional WWTP OPEX from electricity
           E_item=lambda: B4.power_utility.rate*(365*24)*10)


# %%

# =============================================================================
# Scenario C (sysC): containaer-based sanitation with existing treatment system
# =============================================================================

flowsheetC = bst.Flowsheet('sysC')
bst.main_flowsheet.set_flowsheet(flowsheetC)
streamsC = batch_create_streams('C')

#################### Human Inputs ####################
C1 = su.Excretion('C1', outs=('urine', 'feces'))

################### User Interface ###################
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

##################### Conveyance #####################
# Liquid waste
handcart_fee = 0.01 # USD/cap/d
truck_fee = 23e3 # UGX/m3

get_handcart_and_truck_fee = \
    lambda vol, ppl: truck_fee/get_exchange_rate()*vol \
        + handcart_fee*ppl*C2.collection_period
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
    ppl = get_ppl('exist') / C2.N_toilet
    truck3.load = C3.F_mass_in * hr / C2.N_toilet
    truck4.load = C4.F_mass_in * hr / C2.N_toilet
    rho3 = C3.F_mass_in/C3.F_vol_in
    rho4 = C4.F_mass_in/C4.F_vol_in
    C3.fee = get_handcart_and_truck_fee(truck3.load/rho3, ppl)
    C4.fee = get_handcart_and_truck_fee(truck4.load/rho4, ppl)
    C3._design()
    C4._design()
C4.specification = update_C3_C4_param

###################### Treatment ######################
C5 = su.LumpedCost('C5', ins=(C3-0, C4-0),
                   cost_item_name='Lumped WWTP',
                   CAPEX=18606700, power=57120/(365*24), lifetime=8)
get_C5_lifetime = lambda: C5.lifetime

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
treatC = bst.System('treatC', path=(C5, C6, C7, C8))
C8._cost = lambda: clear_unit_costs(treatC)

################## Reuse or Disposal ##################
C9 = su.CropApplication('C9', ins=C7-0, outs=('liquid_fertilizer', 'reuse_loss'),
                        loss_ratio=app_loss)
C9.specification = lambda: adjust_NH3_loss(C9)

C10 = su.Mixer('C10', ins=(C2-4, C6-1, C7-1, C8-2), outs=streamsC['CH4'])
C10.specification = lambda: add_fugitive_items(C10, 'CH4_item')
C10.line = 'fugitive CH4 mixer'

C11 = su.Mixer('C11', ins=(C2-5, C6-2, C7-2, C8-3), outs=streamsC['N2O'])
C11.specification = lambda: add_fugitive_items(C11, 'N2O_item')
C11.line = 'fugitive N2O mixer'

C12 = su.ComponentSplitter('C12', ins=C8-0,
                           outs=(streamsC['sol_N'], streamsC['sol_P'], streamsC['sol_K'],
                                 'C_sol_non_fertilizers'),
                           split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

C13 = su.ComponentSplitter('C13', ins=C9-0,
                           outs=(streamsC['liq_N'], streamsC['liq_P'], streamsC['liq_K'],
                                 'C_liq_non_fertilizers'),
                           split_keys=(('NH3', 'NonNH3'), 'P', 'K'))

############### Simulation, TEA, and LCA ###############
sysC = bst.System('sysC', path=(C1, C2, C3, C4, treatC, C9, C10, C11, C12, C13))

teaC = SimpleTEA(system=sysC, discount_rate=discount_rate, start_year=2018,
                 lifetime=get_C5_lifetime(), uptime_ratio=1, lang_factor=None,
                 annual_maintenance=0, annual_labor=12*3e6*12/get_exchange_rate(),
                 construction_schedule=None)

lcaC = LCA(system=sysC, lifetime=get_C5_lifetime(), lifetime_unit='yr', uptime_ratio=1,
           # Assuming all additional WWTP OPEX from electricity
           E_item=lambda: C5.power_utility.rate*(365*24)*8)


# %%

# =============================================================================
# Util functions
# =============================================================================

def update_lca_data(kind):
    '''
    Load impact indicator and impact item data.

    Parameters
    ----------
    kind : str
        "original" loads the data from Trimmer et al.
        (TRACI, ecoinvent v3.2),
        "new" loads the data for ReCiPe and TRACI
        (ecoinvent 3.7.1, at the point of substitution).
    '''
    global lca_data_kind

    if lca_data_kind != kind:
        load_lca_data(kind)
        batch_create_stream_items(kind)

        for lca in (lcaA, lcaB, lcaC):
            for i in lca.lca_streams:
                # To refresh the impact items
                source_ID = i.stream_impact_item.source.ID
                i.stream_impact_item.source = ImpactItem.get_item(source_ID)

        Biogas_CFs = ImpactItem.get_item('Biogas_item').CFs
        for k, v in Biogas_CFs.items():
            Biogas_CFs[k] = v * get_biogas_factor()

        for i in sysA, sysB, sysC:
            i.simulate()

        lca_data_kind = kind


def get_total_inputs(unit, ppl):
    if len(unit.ins) == 0: # Excretion units do not have ins
        ins = unit.outs
    else:
        ins = unit.ins
    inputs = {}
    inputs['COD'] = sum(i.COD*i.F_vol/1e3 for i in ins)
    inputs['N'] = sum(i.TN*i.F_vol/1e3 for i in ins)
    inputs['NH3'] = sum(i.imass['NH3'] for i in ins)
    inputs['P'] = sum(i.TP*i.F_vol/1e3 for i in ins)
    inputs['K'] = sum(i.TK*i.F_vol/1e3 for i in ins)
    for i, j in inputs.items():
        inputs[i] = j * ppl
    return inputs


def get_recovery(ins, outs, ppl):
    if not isinstance(outs, Iterable):
        outs = (outs,)

    non_g = tuple(i for i in outs if i.phase != 'g')
    recovery = {}
    recovery['COD'] = sum(i.COD*i.F_vol/1e3 for i in non_g)
    recovery['N'] = sum(i.TN*i.F_vol/1e3 for i in non_g)
    recovery['NH3'] = sum(i.imass['NH3'] for i in non_g)
    recovery['P'] = sum(i.TP*i.F_vol/1e3 for i in non_g)
    recovery['K'] = sum(i.TK*i.F_vol/1e3 for i in non_g)

    for i, j in recovery.items():
        inputs = get_total_inputs(ins, ppl)
        recovery[i] /= inputs[i]

    return recovery


sys_dct = {
    'ppl': dict(sysA=get_ppl('exist'), sysB=get_ppl('alt'), sysC=get_ppl('exist')),
    'input_unit': dict(sysA=A1, sysB=B1, sysC=C1),
    'liq_unit': dict(sysA=A13, sysB=B13, sysC=C13),
    'sol_unit': dict(sysA=A12, sysB=B12, sysC=C12),
    'gas_unit': dict(sysA=None, sysB=B15, sysC=None),
    'stream_dct': dict(sysA=streamsA, sysB=streamsB, sysC=streamsC),
    'TEA': dict(sysA=teaA, sysB=teaB, sysC=teaC),
    'LCA': dict(sysA=lcaA, sysB=lcaB, sysC=lcaC),
    'cache': dict(sysA={}, sysB={}, sysC={}),
    }


def cache_recoveries(sys):
    sys_dct['ppl'][sys.ID] = ppl = get_ppl('alt') if sys.ID=='sysB' else get_ppl('exist')
    total_COD = get_total_inputs(sys_dct['input_unit'][sys.ID], ppl)['COD']

    if sys_dct['gas_unit'][sys.ID]:
        gas_mol = sys_dct['gas_unit'][sys.ID].outs[0].imol['CH4']
        gas_COD = gas_mol*1e3*biogas_energy/14e3/total_COD
    else:
        gas_COD = 0

    sys_dct['cache'][sys.ID] = cache = {
        'liq': get_recovery(ins=sys_dct['input_unit'][sys.ID],
                            outs=sys_dct['liq_unit'][sys.ID].ins,
                            ppl=ppl),
        'sol': get_recovery(ins=sys_dct['input_unit'][sys.ID],
                            outs=sys_dct['sol_unit'][sys.ID].ins,
                            ppl=ppl),
        'gas': dict(COD=gas_COD, N=0, P=0, K=0)
        }
    return cache


sysA._set_facilities([*sysA.facilities, lambda: cache_recoveries(sysA)])
sysB._set_facilities([*sysB.facilities, lambda: cache_recoveries(sysB)])
sysC._set_facilities([*sysC.facilities, lambda: cache_recoveries(sysC)])


def get_summarizing_functions(system):
    func_dct = {}
    func_dct['get_annual_net_cost'] = lambda tea, ppl: (tea.EAC-tea.sales)/ppl
    func_dct['get_annual_cost'] = lambda tea, ppl: tea.EAC/ppl
    func_dct['get_annual_CAPEX'] = lambda tea, ppl: tea.annualized_CAPEX/ppl
    func_dct['get_annual_OPEX'] = lambda tea, ppl: tea.AOC/ppl
    func_dct['get_annual_sales'] = lambda tea, ppl: tea.sales/ppl

    for i in ('COD', 'N', 'P', 'K'):
        func_dct[f'get_liq_{i}_recovery'] = \
            lambda sys, i: sys_dct['cache'][sys.ID]['liq'][i]
        func_dct[f'get_sol_{i}_recovery'] = \
            lambda sys, i: sys_dct['cache'][sys.ID]['sol'][i]
        func_dct[f'get_gas_{i}_recovery'] = \
            lambda sys, i: sys_dct['cache'][sys.ID]['gas'][i]
        func_dct[f'get_tot_{i}_recovery'] = \
            lambda sys, i: \
                sys_dct['cache'][sys.ID]['liq'][i] + \
                sys_dct['cache'][sys.ID]['sol'][i] + \
                sys_dct['cache'][sys.ID]['gas'][i]

    return func_dct


def print_summaries(systems):
    if not isinstance(systems, Iterable):
        systems = (systems, )

    for sys in systems:
        func = get_summarizing_functions(sys)
        sys.simulate()
        ppl = sys_dct['ppl'][sys.ID]
        print(f'\n---------- Summary for {sys} ----------\n')
        for i in ('COD', 'N', 'P', 'K'):
            print(f'Total {i} recovery is {func[f"get_tot_{i}_recovery"](sys, i):.1%}, '
                  f'{func[f"get_liq_{i}_recovery"](sys, i):.1%} in liquid, '
                  f'{func[f"get_sol_{i}_recovery"](sys, i):.1%} in solid, '
                  f'{func[f"get_gas_{i}_recovery"](sys, i):.1%} in gas.')
        print('\n')

        tea = sys_dct['TEA'][sys.ID]
        tea.show()

        unit = f'{currency}/cap/yr'
        print(f'\nNet cost: {func["get_annual_net_cost"](tea, ppl):.1f} {unit}.')
        print(f'Total cost: {func["get_annual_cost"](tea, ppl):.1f} {unit}.')
        print(f'Capital: {func["get_annual_CAPEX"](tea, ppl):.1f} {unit}.')
        print(f'Operating: {func["get_annual_OPEX"](tea, ppl):.1f} {unit}.')
        print(f'Sales: {func["get_annual_sales"](tea, ppl):.1f} {unit}.')
        print('\n')

        lca = sys_dct['LCA'][sys.ID]
        lca.show()
        print('\n')

        for ind in lca.indicators:
            unit = f'{ind.unit}/cap/yr'
            print(f'\nImpact indicator {ind.ID}:')

            val = lca.total_impacts[ind.ID]/lca.lifetime/ppl
            print(f'\nNet emission: {val:.1f} {unit}.')

            val = lca.total_construction_impacts[ind.ID]/lca.lifetime/ppl
            print(f'Construction: {val:.1f} {unit}.')

            val = lca.total_transportation_impacts[ind.ID]/lca.lifetime/ppl
            print(f'Transportation: {val:.1f} {unit}.')

            val = lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='direct_emission')[ind.ID] \
                /lca.lifetime/ppl
            print(f'Direct emission: {val:.1f} {unit}.')

            val = lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='offset')[ind.ID] \
                /lca.lifetime/ppl
            print(f'Offset: {val:.1f} {unit}.')

            val = lca.total_other_impacts[ind.ID]/lca.lifetime/ppl
            print(f'Other: {val:.1} {unit}.\n')


def save_all_reports():
    import os
    path = os.path.dirname(os.path.realpath(__file__))
    path += '/results'
    if not os.path.isdir(path):
        os.path.mkdir(path)
    del os
    for i in (sysA, sysB, sysC, lcaA, lcaB, lcaC):
        if isinstance(i, bst.System):
            i.simulate()
            i.save_report(f'{path}/{i.ID}.xlsx')
        else:
            i.save_report(f'{path}/{i.system.ID}_lca.xlsx')

__all__ = ('sysA', 'sysB', 'sysC', 'teaA', 'teaB', 'teaC', 'lcaA', 'lcaB', 'lcaC',
           'print_summaries', 'save_all_reports', 'update_lca_data',
           *(i.ID for i in sysA.units),
           *(i.ID for i in sysB.units),
           *(i.ID for i in sysC.units),
           )