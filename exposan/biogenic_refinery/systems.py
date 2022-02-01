#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>

Code is based on Bwaise systems.py developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
"""
import numpy as np
import biosteam as bst
import qsdsan as qs
from sklearn.linear_model import LinearRegression as LR
from qsdsan import sanunits as su
from qsdsan import WasteStream, ImpactIndicator, ImpactItem, StreamImpactItem, SimpleTEA, LCA


from _cmps import cmps


# =============================================================================
# Unit parameters
# =============================================================================

bst.settings.set_thermo(cmps)


currency = qs.currency = 'USD'
bst.speed_up()

household_size = 4
get_household_size = lambda: household_size
household_per_toilet = 4
get_household_per_toilet = lambda: household_per_toilet
get_toilet_user = lambda: get_household_size()*get_household_per_toilet()
#get_toilet_user = lambda: 4

# Number of people served by the one Biogenic Refinery 1018
ppl_12k = 12000
ppl_10k = 10000
get_ppl = lambda kind: ppl_12k if kind=='12k' else ppl_10k

discount_rate = 0.05
get_discount_rate = lambda: discount_rate

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
get_max_CH4_emission = lambda: max_CH4_emission

# Model for tanker truck cost based on capacity (m3)
# price = a*capacity**b -> ln(price) = ln(a) + bln(capacity)
USD_price_dct = np.array((21.62, 32.43, 54.05, 67.57))
capacities = np.array((3, 4.5, 8, 15))
emptying_fee = 0.15
get_emptying_fee = lambda: emptying_fee
def get_tanker_truck_fee(capacity):
    price_dct = USD_price_dct*(1+get_emptying_fee())
    ln_p = np.log(price_dct)
    ln_cap = np.log(capacities)
    model = LR().fit(ln_cap.reshape(-1,1), ln_p.reshape(-1,1))
    [[predicted]] = model.predict(np.array((np.log(capacity))).reshape(1, -1)).tolist()
    cost = np.exp(predicted)
    #how to change cost
    #cost = 0 * capacity
    return cost

# Nutrient loss during applciation
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05

# Energetic content of the biogas
biogas_energy = 803 # kJ/mol CH4
get_biogas_energy = lambda: biogas_energy
LPG_energy = 50 # MJ/kg
get_LPG_energy = lambda: LPG_energy
get_biogas_factor = lambda: get_biogas_energy()/cmps.CH4.MW/get_LPG_energy()

# =============================================================================
# Prices and GWP CFs
# =============================================================================

# Recycled nutrients are sold at a lower price than commercial fertilizers
price_factor = 0.25
get_price_factor = lambda: price_factor

price_ratio = 1 
#changes based on country 

operator_daily_wage = 29
get_operator_daily_wage = lambda: operator_daily_wage


const_daily_wage = 17
get_const_daily_wage = lambda: const_daily_wage
const_person_days = 100

# operator_hourly_wage = 29/8
# get_operator_hr_wage = lambda: operator_hourly_wage


price_dct = {
    'Electricity': 0.13,
    'Concrete': 194*price_ratio,
    'Steel': 2.665*price_ratio,
    'N': 1.507*get_price_factor(),
    'P': 3.983*get_price_factor(),
    'K': 1.333*get_price_factor(),
    'Polymer': 1*price_ratio,
    'Resin': 3.335*price_ratio,
    'FilterBag': 4.81*price_ratio,
    'MgOH2':  0.145*price_ratio,
    'MgCO3': 0.9*price_ratio,
    'H2SO4': 0.3*price_ratio,
    'biochar': 0.014*price_ratio,
    'struvite': 3.983*(31/245)*get_price_factor(),
    'conc_NH3': 1.333*(14/17)*get_price_factor(),
    }
#filter bag is based on 1 bag that costs 4.81 and has a mass of 0.05 kg of nylon
#nylon has GWP of 9.279255342 kgCO2eq  per kg of nylon

#assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
GWP_dct = {
    'Electricity': 0.69,
    'CH4': 28,
    'N2O': 265,
    'N': -5.4,
    'P': -4.9,
    'K': -1.5,
    'Polymer':2.8,
    'Resin': 1.612,
    'FilterBag': 0.464,
    'MgOH2': 1.176277921,
    'MgCO3': 1.176277921,
    'H2SO4': 0.158899487,
    'biochar': 0, #-0.2*0.9*(44/12),
    'struvite': 0, #-4.9*(31/245),
    'conc_NH3': 0, #-5.4*(14/17),
    }

items = ImpactItem.get_all_items()

if not items.get('Excavation'): # prevent from reloading
    import os
    path = os.path.dirname(os.path.realpath(__file__)) + '/data'
    ImpactIndicator.load_from_file(path+'/impact_indicators.tsv')
    item_path = path+'/impact_items.xlsx'
    ImpactItem.load_from_file(item_path)
    del os

GWP = qs.ImpactIndicator.get_indicator('GWP')

bst.PowerUtility.price = price_dct['Electricity']
items['Concrete'].price = price_dct['Concrete']
items['Steel'].price = price_dct['Steel']

# =============================================================================
# Universal units and functions
# =============================================================================

CH4_item = items.get('CH4_item')
N2O_item = items.get('N2O_item')
if not CH4_item or not N2O_item:
    CH4_item = StreamImpactItem(ID='CH4_item', GWP=GWP_dct['CH4'])
    N2O_item = StreamImpactItem(ID='N2O_item', GWP= GWP_dct['N2O'])
    N_item = StreamImpactItem(ID='N_item', GWP= GWP_dct['N'])
    P_item = StreamImpactItem(ID='P_item', GWP= GWP_dct['P'])
    K_item = StreamImpactItem(ID='K_item', GWP= GWP_dct['K'])
    e_item = ImpactItem(ID='e_item', functional_unit='kWh', GWP= GWP_dct['Electricity'])
    polymer_item = StreamImpactItem(ID='polymer_item', GWP= GWP_dct['Polymer'])
    resin_item = StreamImpactItem(ID='resin_item', GWP= GWP_dct['Resin'])
    filter_bag_item = StreamImpactItem(ID='filter_bag_item', GWP= GWP_dct['FilterBag'])
    MgOH2_item = StreamImpactItem(ID='MgOH2_item', GWP= GWP_dct['MgOH2'])
    MgCO3_item = StreamImpactItem(ID='MgCO3_item', GWP= GWP_dct['MgCO3'])
    H2SO4_item = StreamImpactItem(ID='H2SO4_item', GWP= GWP_dct['H2SO4'])
    biochar_item = StreamImpactItem(ID='biochar_item', GWP= GWP_dct['biochar'])
    struvite_item = StreamImpactItem(ID='struvite_item', GWP= GWP_dct['struvite'])
    conc_NH3_item = StreamImpactItem(ID='conc_NH3_item', GWP= GWP_dct['conc_NH3'])
def batch_create_streams(prefix):
    
    stream_dct = {}
    item = ImpactItem.get_item('CH4_item').copy(f'{prefix}_CH4_item', set_as_source=True)
    stream_dct['CH4'] = WasteStream(f'{prefix}_CH4', phase='g', stream_impact_item=item)
    
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
    
    item = ImpactItem.get_item('polymer_item').copy(f'{prefix}_polymer_item', set_as_source=True)
    stream_dct['polymer'] = WasteStream(f'{prefix}_polymer', phase='s', price=price_dct['Polymer'],
                                      stream_impact_item=item)
    
    item = ImpactItem.get_item('resin_item').copy(f'{prefix}_resin_item', set_as_source=True)
    stream_dct['resin'] = WasteStream(f'{prefix}_resin', phase='s', price=price_dct['Resin'],
                                      stream_impact_item=item)
    
    item = ImpactItem.get_item('filter_bag_item').copy(f'{prefix}_filter_bag_item', set_as_source=True)
    stream_dct['filter_bag'] = WasteStream(f'{prefix}_filter_bag', phase='s', price=price_dct['FilterBag'],
                                      stream_impact_item=item)
    
    item = ImpactItem.get_item('MgOH2_item').copy(f'{prefix}_MgOH2_item', set_as_source=True)
    stream_dct['MgOH2'] = WasteStream(f'{prefix}_MgOH2', phase='s', price=price_dct['MgOH2'],
                                      stream_impact_item=item)
    
    item = ImpactItem.get_item('MgCO3_item').copy(f'{prefix}_MgCO3_item', set_as_source=True)
    stream_dct['MgCO3'] = WasteStream(f'{prefix}_MgCO3', phase='s', price=price_dct['MgCO3'],
                                      stream_impact_item=item)
    
    item = ImpactItem.get_item('H2SO4_item').copy(f'{prefix}_H2SO4_item', set_as_source=True)
    stream_dct['H2SO4'] = WasteStream(f'{prefix}_H2SO4', phase='s', price=price_dct['H2SO4'],
                                      stream_impact_item=item)
    
    item = ImpactItem.get_item('biochar_item').copy(f'{prefix}_biochar_item', set_as_source=True)
    stream_dct['biochar'] = WasteStream(f'{prefix}_biochar', phase='s', price=price_dct['biochar'],
                                      stream_impact_item=item)
    
    item = ImpactItem.get_item('struvite_item').copy(f'{prefix}_struvite_item', set_as_source=True)
    stream_dct['struvite'] = WasteStream(f'{prefix}_struvite', phase='s', price=price_dct['struvite'],
                                      stream_impact_item=item)
    
    item = ImpactItem.get_item('conc_NH3_item').copy(f'{prefix}_conc_NH3_item', set_as_source=True)
    stream_dct['conc_NH3'] = WasteStream(f'{prefix}_conc_NH3', phase='s', price=price_dct['conc_NH3'],
                                      stream_impact_item=item)
    
    # stream_dct['CH4'] = WasteStream(f'{prefix}_CH4', phase='g',
    #                                 impact_item=CH4_item.copy(set_as_source=True))
    # stream_dct['N2O'] = WasteStream(f'{prefix}_N2O', phase='g',
    #                                 impact_item=N2O_item.copy(set_as_source=True))
    # stream_dct['liq_N'] = WasteStream(f'{prefix}_liq_N', phase='l', price= price_dct['N'],
    #                                   impact_item=N_item.copy(set_as_source=True))
    # stream_dct['sol_N'] = WasteStream(f'{prefix}_sol_N', phase='l', price= price_dct['N'],
    #                                   impact_item=N_item.copy(set_as_source=True))
    # stream_dct['liq_P'] = WasteStream(f'{prefix}_liq_P', phase='l', price= price_dct['P'],
    #                                   impact_item=P_item.copy(set_as_source=True))
    # stream_dct['sol_P'] = WasteStream(f'{prefix}_sol_P', phase='l', price= price_dct['P'],
    #                                   impact_item=P_item.copy(set_as_source=True))
    # stream_dct['liq_K'] = WasteStream(f'{prefix}_liq_K', phase='l', price= price_dct['K'],
    #                                   impact_item=K_item.copy(set_as_source=True))
    # stream_dct['sol_K'] = WasteStream(f'{prefix}_sol_K', phase='l', price= price_dct['K'],
    #                                   impact_item=K_item.copy(set_as_source=True))
    # stream_dct['polymer'] = WasteStream(f'{prefix}_polymer', phase='s', price= price_dct['Polymer'],
    #                                   impact_item=polymer_item.copy(set_as_source=True))
    # stream_dct['resin'] = WasteStream(f'{prefix}_resin', phase='s', price=price_dct['Resin'],
    #                                   impact_item=resin_item.copy(set_as_source=True))
    # stream_dct['filter_bag'] = WasteStream(f'{prefix}_filter_bag', phase='s', price=price_dct['FilterBag'],
    #                                   impact_item=filter_bag_item.copy(set_as_source=True))
    # stream_dct['MgOH2'] = WasteStream(f'{prefix}_MgOH2', phase='s', price=price_dct['MgOH2'],
    #                                   impact_item=MgOH2_item.copy(set_as_source=True))
    # stream_dct['MgCO3'] = WasteStream(f'{prefix}_MgCO3', phase='s', price=price_dct['MgCO3'],
    #                                   impact_item=MgCO3_item.copy(set_as_source=True))
    # stream_dct['H2SO4'] = WasteStream(f'{prefix}_H2SO4', phase='s', price=price_dct['H2SO4'],
    #                                   impact_item=H2SO4_item.copy(set_as_source=True))
    # stream_dct['biochar'] = WasteStream(f'{prefix}_biochar', phase='s', price=price_dct['biochar'], 
    #                                     impact_item=biochar_item.copy(set_as_source=True))
    # stream_dct['struvite'] = WasteStream(f'{prefix}_struvite', phase='s', price=price_dct['struvite'], 
    #                                     impact_item=struvite_item.copy(set_as_source=True))
    # stream_dct['conc_NH3'] = WasteStream(f'{prefix}_conc_NH3', phase='s', price=price_dct['conc_NH3'], 
    #                                     impact_item=conc_NH3_item.copy(set_as_source=True))
    return stream_dct
    
# def add_fugitive_items(unit, item):
#     unit._run()
#     for i in unit.ins:
#         i.impact_item = item.copy(set_as_source=True)
        
        
def add_fugitive_items(unit, item_ID):
    unit._run()
    for i in unit.ins:
        i.stream_impact_item = ImpactItem.get_item(item_ID).copy(set_as_source=True)
 
carbon_dict = {
'urine': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'feces': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'toilet_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'toilet_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_sol_loss': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pretreat_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'dryer_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'dryer_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pyrolysis_biochar': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pyrolysis_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
    
#liquids
'toilet_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_liq_loss': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pretreat_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'treat_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'bed_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'bed_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
    }

nitrogen_dict = {
'urine': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'feces': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'toilet_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'toilet_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_sol_loss': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pretreat_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'dryer_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'dryer_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pyrolysis_biochar': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pyrolysis_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
    
#liquids
'toilet_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_liq_loss': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pretreat_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'treat_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'bed_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'bed_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'struvite': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'NH3': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
    }


phosphorus_dict = {
'urine': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'feces': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'toilet_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'toilet_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_sol_loss': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pretreat_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'dryer_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'dryer_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pyrolysis_biochar': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pyrolysis_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
    
#liquids
'toilet_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_liq_loss': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pretreat_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'treat_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'bed_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'bed_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'struvite': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'NH3': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
    }

potassium_dict = {
'urine': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'feces': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'toilet_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'toilet_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_sol_loss': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pretreat_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'dryer_sol': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'dryer_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pyrolysis_biochar': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pyrolysis_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
    
#liquids
'toilet_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'trans_liq_loss': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'pretreat_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'treat_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'bed_liq': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'bed_gas': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'struvite': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
'NH3': dict(sysA=0.0, sysB=0.0, sysC=0.0, sysD=0.0),
    }



# %%

# =============================================================================
# Scenario A (sysA): pit latrine with Biogenic Refinery 12k users
# =============================================================================

flowsheetA = bst.Flowsheet('sysA')
bst.main_flowsheet.set_flowsheet(flowsheetA)
streamsA = batch_create_streams('A')

#################### Human Inputs ####################
# !!! how to change excrestion based on location (e.g., calorie and protein)
A1 = su.Excretion('A1', outs=('urine','feces'))

################### User Interface ###################
# !!! how to change inputs based on location (e.g., flushing water, cleaning water, toilet paper)
A2 = su.PitLatrine('A2', ins=(A1-0, A1-1,
                              'toilet_paper', 'flushing_water',
                              'cleansing_water', 'desiccant'),
                   outs=('mixed_waste', 'leachate', 'A2_CH4', 'A2_N2O'),
                   N_user=get_toilet_user(), N_toilet=get_ppl('12k')/get_toilet_user(),
                   if_flushing=False, if_desiccant=False, if_toilet_paper=False,
                   OPEX_over_CAPEX=0.05, lifetime=10,
                  decay_k_COD=get_decay_k(tau_deg, log_deg),
                   decay_k_N=get_decay_k(tau_deg, log_deg),
                   max_CH4_emission=get_max_CH4_emission()
                   )
 
def update_A2_param():
    A2.N_user = get_toilet_user()
    A2.N_toilet = get_ppl('12k')/get_toilet_user()
    A2._run()
A2.specification = update_A2_param

##################### Conveyance of Waste #####################
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
A4 = su.IndustrialControlPanel('A4', ins=A3-0, outs='A4_in')

A5 = su.HousingBiogenicRefinery('A5', ins=A4-0, outs='A5_in', const_wage=get_const_daily_wage(), 
                                const_person_days=const_person_days)

A6 = su.ScrewPress('A6', ins=(A5-0, streamsA['polymer']), outs=('liq', 'cake_sol'))

A7 = su.LiquidTreatmentBed('A7', ins=A6-0, outs=('liquid_bed_treated', 'A7_CH4', 'A7_N2O'),
                           decay_k_COD=get_decay_k(tau_deg, log_deg),
                           decay_k_N=get_decay_k(tau_deg, log_deg),
                           max_CH4_emission=get_max_CH4_emission())

A8 = su.CarbonizerBase('A8', outs = (streamsA['biochar'], 'A8_hot_gas', 'A8_N2O'))
A9 = su.PollutionControlDevice('A9', ins = (A8-1, A8-2), outs = ('A9_hot_gas_pcd', 'A9_N2O'))

# updating uptime_ratio in all units to follow carbonizer base
A9_old_cost = A9._cost
def update_A9_uptime_ratio():
    A12.uptime_ratio = A11.uptime_ratio = A10.uptime_ratio = A9.uptime_ratio = A8.uptime_ratio 
    A9_old_cost()
A9._cost = update_A9_uptime_ratio

A10 = su.OilHeatExchanger('A10', ins = A9-0, outs = ('A10_hot_gas'))
A11 = su.HydronicHeatExchanger('A11', ins = A10-0, outs = ('A11_hot_gas'))
A12 = su.DryerFromHHX('A12', ins = (A6-1, A11-0), outs = ('waste_out', 'A12_N2O', 'A12_CH4'))

A12-0-A8


A13 = su.Mixer('A13', ins=(A2-2, A7-1, A12-2), outs=streamsA['CH4'])
A13.specification = lambda: add_fugitive_items(A13, 'CH4_item')

A13.line = 'fugitive CH4 mixer' 
        
A14 = su.Mixer('A14', ins=(A2-3, A7-2, A9-1, A12-1), outs=streamsA['N2O'])
def update_C_recoveries():
    add_fugitive_items(A14, 'N2O_item')        
    # carbon system A
    time_period = 365*24
    carbon_COD_ratio = A12.carbon_COD_ratio
    #in
    sysA_C_urine = A2.ins[0].COD*A2.ins[0].F_vol/1e3 * carbon_COD_ratio * get_ppl('12k') * time_period # kg C/yr in urine
    sysA_C_feces = A2.ins[1].COD*A2.ins[1].F_vol/1e3 * carbon_COD_ratio * get_ppl('12k') * time_period # kg C/yr in feces
    
    sysA_N_urine = A2.ins[0].TN*A2.ins[0].F_vol/1e3 * get_ppl('12k') * time_period # kg N/yr in urine
    sysA_N_feces = A2.ins[1].TN*A2.ins[1].F_vol/1e3 * get_ppl('12k') * time_period # kg N/yr in feces
    
    sysA_P_urine = A2.ins[0].TP*A2.ins[0].F_vol/1e3 * get_ppl('12k') * time_period # kg P/yr in urine
    sysA_P_feces = A2.ins[1].TP*A2.ins[1].F_vol/1e3 * get_ppl('12k') * time_period # kg P/yr in feces
    
    sysA_K_urine = A2.ins[0].TK*A2.ins[0].F_vol/1e3 * get_ppl('12k') * time_period # kg K/yr in urine
    sysA_K_feces = A2.ins[1].TK*A2.ins[1].F_vol/1e3 * get_ppl('12k') * time_period # kg K/yr in feces
    
    #toilet
    sysA_C_toilet_mixed_out = A2.outs[0].COD*A2.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr in mixed excreta
    #lechate does not contain COD or C
    # C_toilet_lechate_out = A2.outs[1].COD*A2.outs[1].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr in lechate
    sysA_C_toilet_fugitive_out = A2.outs[2].imass['CH4'] * (12/16) * time_period # kg C/yr in fugitive gas
    
    sysA_N_toilet_mixed_out = A2.outs[0].TN*A2.outs[0].F_vol/1e3 * time_period # kg N/yr in mixed excreta
    sysA_N_toilet_lechate_out = A2.outs[1].TN*A2.outs[1].F_vol/1e3 *  time_period # kg N/yr in lechate
    sysA_N_toilet_fugitive_out = A2.outs[3].imass['N2O'] * (14/44) * time_period # kg N/yr in fugitive gas
    
    sysA_P_toilet_mixed_out = A2.outs[0].TP*A2.outs[0].F_vol/1e3 * time_period # kg P/yr in mixed excreta
    sysA_P_toilet_lechate_out = A2.outs[1].TP*A2.outs[1].F_vol/1e3 *  time_period # kg P/yr in lechate
    
    sysA_K_toilet_mixed_out = A2.outs[0].TK*A2.outs[0].F_vol/1e3 * time_period # kg K/yr in mixed excreta
    sysA_K_toilet_lechate_out = A2.outs[1].TK*A2.outs[1].F_vol/1e3 *  time_period # kg K/yr in lechate
    
    #transported
    sysA_C_transport_out = A3.outs[0].COD * A3.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr transported
    sysA_C_transport_loss = A3.outs[1].COD * A3.outs[1].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr loss
    
    sysA_N_transport_out = A3.outs[0].TN * A3.outs[0].F_vol/1e3 * time_period # kg N/yr transported
    sysA_N_transport_loss = A3.outs[1].TN * A3.outs[1].F_vol/1e3 * time_period # kg N/yr loss
    
    sysA_P_transport_out = A3.outs[0].TP * A3.outs[0].F_vol/1e3 * time_period # kg P/yr transported
    sysA_P_transport_loss = A3.outs[1].TP * A3.outs[1].F_vol/1e3 * time_period # kg P/yr loss
    
    sysA_K_transport_out = A3.outs[0].TK * A3.outs[0].F_vol/1e3 * time_period # kg K/yr transported
    sysA_K_transport_loss = A3.outs[1].TK * A3.outs[1].F_vol/1e3 * time_period # kg K/yr loss
    
    #pretreatment
    sysA_C_pretreatment_liq = A6.outs[0].COD * A6.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr liq
    sysA_C_pretreatment_solid = A6.outs[1].COD * A6.outs[1].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr solid
    
    sysA_N_pretreatment_liq = A6.outs[0].TN * A6.outs[0].F_vol/1e3 * time_period # kg N/yr liq
    sysA_N_pretreatment_solid = A6.outs[1].TN * A6.outs[1].F_vol/1e3 * time_period # kg N/yr solid
    
    sysA_P_pretreatment_liq = A6.outs[0].TP * A6.outs[0].F_vol/1e3 * time_period # kg P/yr liq
    sysA_P_pretreatment_solid = A6.outs[1].TP * A6.outs[1].F_vol/1e3 * time_period # kg P/yr solid
    
    sysA_K_pretreatment_liq = A6.outs[0].TK * A6.outs[0].F_vol/1e3 * time_period # kg K/yr liq
    sysA_K_pretreatment_solid = A6.outs[1].TK * A6.outs[1].F_vol/1e3 * time_period # kg K/yr solid
    
    #liq bed
    sysA_C_liq_bed_out = A7.outs[0].COD * A7.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period
    sysA_C_liq_bed_gas = A7.outs[1].imass['CH4'] * (12/16) * time_period
    
    sysA_N_liq_bed_out = A7.outs[0].TN * A7.outs[0].F_vol/1e3 * time_period
    sysA_N_liq_bed_gas = A7.outs[2].imass['N2O'] * (14/44) * time_period
    
    sysA_P_liq_bed_out = A7.outs[0].TP * A7.outs[0].F_vol/1e3 * time_period
    
    sysA_K_liq_bed_out = A7.outs[0].TK * A7.outs[0].F_vol/1e3 * time_period
    
    #dryer
    sysA_C_dryer_out = A12.outs[0].COD * A12.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr solids to pyrolysis
    sysA_C_dryer_fugitive_out = A12.outs[2].imass['CH4'] * (12/16) * time_period # kg C/yr in fugitive gas
    
    sysA_N_dryer_out = A12.outs[0].TN * A12.outs[0].F_vol/1e3 * time_period # kg N/yr solids to pyrolysis
    sysA_N_dryer_fugitive_out = A12.outs[1].imass['N2O'] * (14/44) * time_period # kg C/yr in fugitive gas
    
    sysA_P_dryer_out = A12.outs[0].TP * A12.outs[0].F_vol/1e3 * time_period # kg P/yr solids to pyrolysis
    
    sysA_K_dryer_out = A12.outs[0].TK * A12.outs[0].F_vol/1e3 * time_period # kg K/yr solids to pyrolysis
    
    #pyrolysis
    sysA_C_pyrolysis_biochar = A8.outs[0].imass['C'] * time_period #kg C/yr in biochar
    sysA_C_pyrolysis_gas = A8.ins[0].COD * A8.carbon_COD_ratio * A8.ins[0].F_vol / 1e3 * (A8.pyrolysis_C_loss) * time_period #kg C/yr in fugitive gas
    #sysA_C_pyrolysis_gas = sysA_C_dryer_out - sysA_C_pyrolysis_biochar #kg C/yr in fugitive gas

    sysA_N_pyrolysis_biochar = A8.outs[0].imass['N'] * time_period #kg N/yr in biochar
    sysA_N_pyrolysis_gas = A8.ins[0].TN * A8.ins[0].F_vol / 1e3 * (A8.pyrolysis_N_loss) * time_period #kg N/yr in loss
    
    sysA_P_pyrolysis_biochar = A8.outs[0].imass['P'] * time_period #kg P/yr in biochar
    sysA_P_pyrolysis_gas = A8.ins[0].TP * A8.ins[0].F_vol / 1e3 * (A8.pyrolysis_P_loss) * time_period #kg P/yr in loss

    sysA_K_pyrolysis_biochar = A8.outs[0].imass['K'] * time_period #kg K/yr in biochar
    sysA_K_pyrolysis_gas = A8.ins[0].TK * A8.ins[0].F_vol / 1e3 * (A8.pyrolysis_K_loss) * time_period #kg K/yr in loss

    #add all values to dictionary
    carbon_dict['urine']['sysA']=sysA_C_urine
    carbon_dict['feces']['sysA']=sysA_C_feces
    carbon_dict['toilet_sol']['sysA']=sysA_C_toilet_mixed_out
    carbon_dict['toilet_gas']['sysA']=sysA_C_toilet_fugitive_out
    carbon_dict['trans_sol']['sysA']=sysA_C_transport_out
    carbon_dict['trans_sol_loss']['sysA']=sysA_C_transport_loss
    carbon_dict['pretreat_sol']['sysA']=sysA_C_pretreatment_solid
    carbon_dict['dryer_sol']['sysA']=sysA_C_dryer_out
    carbon_dict['dryer_gas']['sysA']=sysA_C_dryer_fugitive_out
    carbon_dict['pyrolysis_biochar']['sysA']=sysA_C_pyrolysis_biochar
    carbon_dict['pyrolysis_gas']['sysA']=sysA_C_pyrolysis_gas
    carbon_dict['toilet_liq']['sysA']=0.0
    carbon_dict['trans_liq']['sysA']=0.0
    carbon_dict['trans_liq_loss']['sysA']=0.0
    carbon_dict['pretreat_liq']['sysA']=sysA_C_pretreatment_liq
    carbon_dict['treat_liq']['sysA']=0.0
    carbon_dict['bed_liq']['sysA']=sysA_C_liq_bed_out
    carbon_dict['bed_gas']['sysA']=sysA_C_liq_bed_gas
    
    nitrogen_dict['urine']['sysA']=sysA_N_urine
    nitrogen_dict['feces']['sysA']=sysA_N_feces
    nitrogen_dict['toilet_sol']['sysA']=sysA_N_toilet_mixed_out
    nitrogen_dict['toilet_gas']['sysA']=sysA_N_toilet_fugitive_out
    nitrogen_dict['trans_sol']['sysA']=sysA_N_transport_out
    nitrogen_dict['trans_sol_loss']['sysA']=sysA_N_transport_loss
    nitrogen_dict['pretreat_sol']['sysA']=sysA_N_pretreatment_solid
    nitrogen_dict['dryer_sol']['sysA']=sysA_N_dryer_out
    nitrogen_dict['dryer_gas']['sysA']=sysA_N_dryer_fugitive_out
    nitrogen_dict['pyrolysis_biochar']['sysA']=sysA_N_pyrolysis_biochar
    nitrogen_dict['pyrolysis_gas']['sysA']=sysA_N_pyrolysis_gas
    nitrogen_dict['toilet_liq']['sysA']=0.0
    nitrogen_dict['trans_liq']['sysA']=0.0
    nitrogen_dict['trans_liq_loss']['sysA']=0.0
    nitrogen_dict['pretreat_liq']['sysA']=sysA_N_pretreatment_liq
    nitrogen_dict['treat_liq']['sysA']=0.0
    nitrogen_dict['bed_liq']['sysA']=sysA_N_liq_bed_out
    nitrogen_dict['bed_gas']['sysA']=sysA_N_liq_bed_gas
    
    phosphorus_dict['urine']['sysA']=sysA_P_urine
    phosphorus_dict['feces']['sysA']=sysA_P_feces
    phosphorus_dict['toilet_sol']['sysA']=sysA_P_toilet_mixed_out
    phosphorus_dict['trans_sol']['sysA']=sysA_P_transport_out
    phosphorus_dict['trans_sol_loss']['sysA']=sysA_P_transport_loss
    phosphorus_dict['pretreat_sol']['sysA']=sysA_P_pretreatment_solid
    phosphorus_dict['dryer_sol']['sysA']=sysA_P_dryer_out
    phosphorus_dict['pyrolysis_biochar']['sysA']=sysA_P_pyrolysis_biochar
    phosphorus_dict['pyrolysis_gas']['sysA']=sysA_P_pyrolysis_gas
    phosphorus_dict['toilet_liq']['sysA']=0.0
    phosphorus_dict['trans_liq']['sysA']=0.0
    phosphorus_dict['trans_liq_loss']['sysA']=0.0
    phosphorus_dict['pretreat_liq']['sysA']=sysA_P_pretreatment_liq
    phosphorus_dict['treat_liq']['sysA']=0.0
    phosphorus_dict['bed_liq']['sysA']=sysA_P_liq_bed_out
    
    potassium_dict['urine']['sysA']=sysA_K_urine
    potassium_dict['feces']['sysA']=sysA_K_feces
    potassium_dict['toilet_sol']['sysA']=sysA_K_toilet_mixed_out
    potassium_dict['trans_sol']['sysA']=sysA_K_transport_out
    potassium_dict['trans_sol_loss']['sysA']=sysA_K_transport_loss
    potassium_dict['pretreat_sol']['sysA']=sysA_K_pretreatment_solid
    potassium_dict['dryer_sol']['sysA']=sysA_K_dryer_out
    potassium_dict['pyrolysis_biochar']['sysA']=sysA_K_pyrolysis_biochar
    potassium_dict['pyrolysis_gas']['sysA']=sysA_K_pyrolysis_gas
    potassium_dict['toilet_liq']['sysA']=0.0
    potassium_dict['trans_liq']['sysA']=0.0
    potassium_dict['trans_liq_loss']['sysA']=0.0
    potassium_dict['pretreat_liq']['sysA']=sysA_K_pretreatment_liq
    potassium_dict['treat_liq']['sysA']=0.0
    potassium_dict['bed_liq']['sysA']=sysA_K_liq_bed_out
    
A14.specification = update_C_recoveries
A14.line = 'fugitive N2O mixer'



############### Simulation, TEA, and LCA ###############
sysA = bst.System('sysA', path= (A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14))

sysA.simulate()



teaA = SimpleTEA(system=sysA, discount_rate=get_discount_rate(), 
                  start_year=2020, lifetime=20, uptime_ratio=1, 
                  lang_factor=None, annual_maintenance=0, 
                  annual_labor=(get_operator_daily_wage()*3*365))
powerA = sum([(u.power_utility.rate*u.uptime_ratio) for u in sysA.units])
def get_powerA():
    powerA = sum([(u.power_utility.rate*u.uptime_ratio) for u in sysA.units])
    return powerA
  
lcaA = LCA(system=sysA, lifetime=20, lifetime_unit='yr', uptime_ratio=1,
            e_item=lambda: (powerA*(365*12)*20))



# %%

# =============================================================================
# Scenario B (sysB): UDDT with Biogenic Refinery 12k users
# =============================================================================

flowsheetB = bst.Flowsheet('sysB')
bst.main_flowsheet.set_flowsheet(flowsheetB)
streamsB = batch_create_streams('B')

#################### Human Inputs ####################
B1 = su.Excretion('B1', outs=('urine','feces'))

################### User Interface ###################
B2 = su.UDDT('B2', ins=(B1-0, B1-1,
                        'toilet_paper', 'flushing_water',
                        'cleaning_water', 'desiccant'),
             outs=('liq_waste', 'sol_waste',
                   'struvite', 'HAP', 'B2_CH4', 'B2_N2O'),
             N_user=get_toilet_user(), N_toilet=get_ppl('12k')/get_toilet_user(), 
             if_flushing=False, if_desiccant=False, if_toilet_paper=False, 
             OPEX_over_CAPEX=0.1, lifetime=10,
             decay_k_COD=get_decay_k(tau_deg, log_deg),
             decay_k_N=get_decay_k(tau_deg, log_deg),
             max_CH4_emission=get_max_CH4_emission())

def update_B2_param():
    B2.N_user = get_toilet_user()
    B2.N_toilet = get_ppl('12k')/get_toilet_user()
    B2._run()
B2.specification = update_B2_param

 
##################### Conveyance of Waste #####################
# Liquid waste
handcart_fee = 0.01 # USD/cap/d
get_handcart_fee = lambda: handcart_fee
truck_fee = 6.21 # USD/m3
get_truck_fee = lambda: truck_fee

get_handcart_and_truck_fee = \
    lambda vol, ppl_t: get_truck_fee()*vol \
        + get_handcart_fee()*ppl_t*B2.collection_period
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
    ppl_t = get_ppl('12k') / B2.N_toilet
    truck3.load = B3.F_mass_in * hr / B2.N_toilet
    truck4.load = B4.F_mass_in * hr / B2.N_toilet
    rho3 = B3.F_mass_in/B3.F_vol_in
    rho4 = B4.F_mass_in/B4.F_vol_in
    B3.fee = get_handcart_and_truck_fee(truck3.load/rho3, ppl_t)
    B4.fee = get_handcart_and_truck_fee(truck4.load/rho4, ppl_t)
    B3._design()
    B4._design()
B4.specification = update_B3_B4_param

###################### Treatment ######################
#B5 = su.StruvitePrecipitation('B5', ins=(B3-0, streamsB['MgOH2'], streamsB['filter_bag']), outs=('treated', streamsB['struvite']))

B5 = su.StruvitePrecipitation('B5', ins=(B3-0, streamsB['MgOH2'], streamsB['MgCO3'], 
                                      streamsB['filter_bag']), outs=('B5_treated', streamsB['struvite']))
B6 = su.IonExchangeNH3('B6', ins=(B5-0, streamsB['resin'], streamsB['H2SO4']), 
                    outs=('B6_treated', 'resin_out', streamsB['conc_NH3']))

B7 = su.LiquidTreatmentBed('B7', ins=B6-0, outs=('liquid_bed_treated', 'B7_CH4', 'B7_N2O'),
                           decay_k_COD=get_decay_k(tau_deg, log_deg),
                           decay_k_N=get_decay_k(tau_deg, log_deg),
                           max_CH4_emission=get_max_CH4_emission())

B8 = su.IndustrialControlPanel('B8', ins=B7-0, outs='B8_in')
B9 = su.HousingBiogenicRefinery('B9', ins=B8-0, outs='B9_in', const_wage=get_const_daily_wage(), 
                                const_person_days=const_person_days)

B10 = su.Grinder('B10', ins=B4-0, outs=('waste'))
B11 = su.CarbonizerBase('B11', outs = (streamsB['biochar'], 'B11_hot_gas', 'B11_N2O'))
 
B12 = su.PollutionControlDevice('B12', ins = (B11-1, B11-2), outs = ('B12_hot_gas_pcd', 'B12_N2O'))

# updating uptime_ratio in all units to follow carbonizer base
B12_old_cost = B12._cost
def update_B12_uptime_ratio():
    B15.uptime_ratio = B14.uptime_ratio = B13.uptime_ratio = B12.uptime_ratio = B11.uptime_ratio 
    B12_old_cost()
B12._cost = update_B12_uptime_ratio

B13 = su.OilHeatExchanger('B13', ins = B12-0, outs = ('B13_hot_gas'))
B14 = su.HydronicHeatExchanger('B14', ins = B13-0, outs = ('B14_hot_gas'))
B15 = su.DryerFromHHX('B15', ins = (B10-0, B14-0), outs = ('waste_out', 'B15_N2O', 'B15_CH4'))

B15-0-B11

B16 = su.Mixer('B16', ins=(B2-4, B7-1, B15-2), outs=streamsB['CH4'])
B16.specification = lambda: add_fugitive_items(B16, 'CH4_item')
B16.line = 'fugitive CH4 mixer' 
        
B17 = su.Mixer('B17', ins=(B2-5, B7-2, B12-1, B15-1), outs=streamsB['N2O'])

def update_C_recoveries():
    
    add_fugitive_items(B17, 'N2O_item')
    time_period = 365*24
    #use same carbon COD ratio for all calculations
    carbon_COD_ratio = A12.carbon_COD_ratio
    
    # carbon system B
    #in
    sysB_C_urine = B2.ins[0].COD*B2.ins[0].F_vol/1e3 * carbon_COD_ratio * get_ppl('12k') * time_period # kg C/yr in urine
    sysB_C_feces = B2.ins[1].COD*B2.ins[1].F_vol/1e3 * carbon_COD_ratio * get_ppl('12k') * time_period # kg C/yr in feces
    
    sysB_N_urine = B2.ins[0].TN*B2.ins[0].F_vol/1e3 * get_ppl('12k') * time_period # kg N/yr in urine
    sysB_N_feces = B2.ins[1].TN*B2.ins[1].F_vol/1e3 * get_ppl('12k') * time_period # kg N/yr in feces
    
    sysB_P_urine = B2.ins[0].TP*B2.ins[0].F_vol/1e3 * get_ppl('12k') * time_period # kg P/yr in urine
    sysB_P_feces = B2.ins[1].TP*B2.ins[1].F_vol/1e3 * get_ppl('12k') * time_period # kg P/yr in feces
    
    sysB_K_urine = B2.ins[0].TK*B2.ins[0].F_vol/1e3 * get_ppl('12k') * time_period # kg K/yr in urine
    sysB_K_feces = B2.ins[1].TK*B2.ins[1].F_vol/1e3 * get_ppl('12k') * time_period # kg K/yr in feces
    
    #toilet
    sysB_C_toilet_solid_out = B2.outs[1].COD*B2.outs[1].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr in solid excreta
    sysB_C_toilet_liq_out = B2.outs[0].COD*B2.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr in liq excreta
    sysB_C_toilet_fugitive_out = B2.outs[4].imass['CH4'] * (12/16) * time_period # kg C/yr in fugitive gas
    
    sysB_N_toilet_solid_out = B2.outs[1].TN*B2.outs[1].F_vol/1e3 *  time_period # kg N/yr in solid excreta
    sysB_N_toilet_liq_out = B2.outs[0].TN*B2.outs[0].F_vol/1e3 * time_period # kg N/yr in liq excreta
    sysB_N_toilet_fugitive_out = B2.outs[3].imass['N2O'] * (14/44) * time_period # kg N/yr in fugitive gas
    
    sysB_P_toilet_solid_out = B2.outs[1].TP*B2.outs[1].F_vol/1e3 *  time_period # kg P/yr in solid excreta
    sysB_P_toilet_liq_out = B2.outs[0].TP*B2.outs[0].F_vol/1e3 * time_period # kg P/yr in liq excreta
    
    sysB_K_toilet_solid_out = B2.outs[1].TK*B2.outs[1].F_vol/1e3 *  time_period # kg K/yr in solid excreta
    sysB_K_toilet_liq_out = B2.outs[0].TK*B2.outs[0].F_vol/1e3 * time_period # kg K/yr in liq excreta
    
    #solid
    #transported
    sysB_C_transport_out = B4.outs[0].COD * B4.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr transported
    sysB_C_transport_loss = B4.outs[1].COD * B4.outs[1].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr loss
    
    sysB_N_transport_out = B4.outs[0].TN * B4.outs[0].F_vol/1e3 * time_period # kg N/yr transported
    sysB_N_transport_loss = B4.outs[1].TN * B4.outs[1].F_vol/1e3 * time_period # kg N/yr loss
    
    sysB_P_transport_out = B4.outs[0].TP * B4.outs[0].F_vol/1e3 * time_period # kg P/yr transported
    sysB_P_transport_loss = B4.outs[1].TP * B4.outs[1].F_vol/1e3 * time_period # kg P/yr loss
    
    sysB_K_transport_out = B4.outs[0].TK * B4.outs[0].F_vol/1e3 * time_period # kg K/yr transported
    sysB_K_transport_loss = B4.outs[1].TK * B4.outs[1].F_vol/1e3 * time_period # kg K/yr loss
    
    #pretreatment
    sysB_C_pretreatment_solid = B10.outs[0].COD * B10.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr liq
    
    sysB_N_pretreatment_solid = B10.outs[0].TN * B10.outs[0].F_vol/1e3 * time_period # kg N/yr liq
    
    sysB_P_pretreatment_solid = B10.outs[0].TP * B10.outs[0].F_vol/1e3 * time_period # kg P/yr liq
    
    sysB_K_pretreatment_solid = B10.outs[0].TK * B10.outs[0].F_vol/1e3 * time_period # kg K/yr liq

    #dryer
    sysB_C_dryer_out = B15.outs[0].COD * B15.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr solids to pyrolysis
    sysB_C_dryer_fugitive_out = B15.outs[2].imass['CH4'] * (12/16) * time_period # kg C/yr in fugitive gas

    sysB_N_dryer_out = B15.outs[0].TN * B15.outs[0].F_vol/1e3 * time_period # kg N/yr solids to pyrolysis
    sysB_N_dryer_fugitive_out = B15.outs[1].imass['N2O'] * (14/44) * time_period # kg N/yr in fugitive gas
    
    sysB_P_dryer_out = B15.outs[0].TP * B15.outs[0].F_vol/1e3 * time_period # kg P/yr solids to pyrolysis
    
    sysB_K_dryer_out = B15.outs[0].TK * B15.outs[0].F_vol/1e3 * time_period # kg K/yr solids to pyrolysis
    
    #pyrolysis
    sysB_C_pyrolysis_biochar = B11.outs[0].imass['C'] * time_period #kg C/yr in biochar
    sysB_C_pyrolysis_gas = B11.ins[0].COD * B11.carbon_COD_ratio * B11.ins[0].F_vol / 1e3 * (B11.pyrolysis_C_loss) * time_period #kg C/yr in fugitive gas
    #sysB_C_pyrolysis_gas = sysB_C_dryer_out - sysB_C_pyrolysis_biochar #kg C/yr in fugutive gas
    
    sysB_N_pyrolysis_biochar = B11.outs[0].imass['N'] * time_period #kg N/yr in biochar
    sysB_N_pyrolysis_gas = B11.ins[0].TN *  B11.ins[0].F_vol / 1e3 * (B11.pyrolysis_N_loss) * time_period #kg N/yr in fugitive gas

    sysB_P_pyrolysis_biochar = B11.outs[0].imass['P'] * time_period #kg P/yr in biochar
    sysB_P_pyrolysis_gas = B11.ins[0].TP *  B11.ins[0].F_vol / 1e3 * (B11.pyrolysis_P_loss) * time_period #kg P/yr in loss
    
    sysB_K_pyrolysis_biochar = B11.outs[0].imass['K'] * time_period #kg K/yr in biochar
    sysB_K_pyrolysis_gas = B11.ins[0].TK *  B11.ins[0].F_vol / 1e3 * (B11.pyrolysis_K_loss) * time_period #kg K/yr in loss
    
    #liquid
    #transported
    sysB_C_liq_transport_out = B3.outs[0].COD * B3.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr transported
    sysA_C_liq_transport_loss = B3.outs[1].COD * B3.outs[1].F_vol/1e3 * carbon_COD_ratio * time_period # kg C/yr loss
    
    sysB_N_liq_transport_out = B3.outs[0].TN * B3.outs[0].F_vol/1e3 * time_period # kg N/yr transported
    sysA_N_liq_transport_loss = B3.outs[1].TN * B3.outs[1].F_vol/1e3 * time_period # kg N/yr loss
    
    sysB_P_liq_transport_out = B3.outs[0].TP * B3.outs[0].F_vol/1e3 * time_period # kg P/yr transported
    sysA_P_liq_transport_loss = B3.outs[1].TP * B3.outs[1].F_vol/1e3 * time_period # kg P/yr loss
    
    sysB_K_liq_transport_out = B3.outs[0].TK * B3.outs[0].F_vol/1e3 * time_period # kg K/yr transported
    sysA_K_liq_transport_loss = B3.outs[1].TK * B3.outs[1].F_vol/1e3 * time_period # kg K/yr loss
    
    #struvite and ion exchange
    sysB_N_liq_struvite_recovered = B5.outs[1].imass['Struvite'] * (14/245.41) * time_period # kg N as struvite/yr recovered
    sysB_N_liq_struvite_out = B5.outs[0].TN * B5.outs[1].F_vol/1e3 * time_period # kg N /yr 
    
    sysB_P_liq_struvite_recovered = B5.outs[1].imass['Struvite'] * (30/245.41) * time_period # kg P as struvite/yr recovered
    sysB_P_liq_struvite_out = B5.outs[0].TP * B5.outs[1].F_vol/1e3 * time_period # kg P /yr 

    sysB_C_liq_ionx_out = B6.outs[0].COD * B6.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period
    
    sysB_N_liq_ionx_out = B6.outs[0].TN * B6.outs[0].F_vol/1e3 * time_period
    sysB_N_liq_recovered = B6.outs[2].TN * B6.outs[2].F_vol/1e3 * time_period

    sysB_P_liq_ionx_out = B6.outs[0].TP * B6.outs[0].F_vol/1e3 * time_period
    
    sysB_K_liq_ionx_out = B6.outs[0].TK * B6.outs[0].F_vol/1e3 * time_period

    
    #liq bed
    sysB_C_liq_bed_out = B7.outs[0].COD * B7.outs[0].F_vol/1e3 * carbon_COD_ratio * time_period
    sysB_C_liq_bed_gas = B7.outs[1].imass['CH4'] * (12/16) * time_period
    
    sysB_N_liq_bed_out = B7.outs[0].TN * B7.outs[0].F_vol/1e3 * time_period
    sysB_N_liq_bed_gas = B7.outs[2].imass['N2O'] * (14/44) * time_period
    
    sysB_P_liq_bed_out = B7.outs[0].TP * B7.outs[0].F_vol/1e3 * time_period
    
    sysB_K_liq_bed_out = B7.outs[0].TK * B7.outs[0].F_vol/1e3 * time_period
    
    #add all values to dictionary
    carbon_dict['urine']['sysB']=sysB_C_urine
    carbon_dict['feces']['sysB']=sysB_C_feces
    carbon_dict['toilet_sol']['sysB']=sysB_C_toilet_solid_out
    carbon_dict['toilet_gas']['sysB']=sysB_C_toilet_fugitive_out
    carbon_dict['trans_sol']['sysB']=sysB_C_transport_out
    carbon_dict['trans_sol_loss']['sysB']=sysB_C_transport_loss
    carbon_dict['pretreat_sol']['sysB']=sysB_C_pretreatment_solid
    carbon_dict['dryer_sol']['sysB']=sysB_C_dryer_out
    carbon_dict['dryer_gas']['sysB']=sysB_C_dryer_fugitive_out
    carbon_dict['pyrolysis_biochar']['sysB']=sysB_C_pyrolysis_biochar
    carbon_dict['pyrolysis_gas']['sysB']=sysB_C_pyrolysis_gas
    carbon_dict['toilet_liq']['sysB']=sysB_C_toilet_liq_out
    carbon_dict['trans_liq']['sysB']=sysB_C_liq_transport_out
    carbon_dict['trans_liq_loss']['sysB']=sysA_C_liq_transport_loss
    carbon_dict['pretreat_liq']['sysB']=0.0
    carbon_dict['treat_liq']['sysB']=sysA_C_liq_transport_loss
    carbon_dict['bed_liq']['sysB']=sysB_C_liq_bed_out
    carbon_dict['bed_gas']['sysB']=sysB_C_liq_bed_gas
    
    nitrogen_dict['urine']['sysB']=sysB_N_urine
    nitrogen_dict['feces']['sysB']=sysB_N_feces
    nitrogen_dict['toilet_sol']['sysB']=sysB_N_toilet_solid_out
    nitrogen_dict['toilet_gas']['sysB']=sysB_N_toilet_fugitive_out
    nitrogen_dict['trans_sol']['sysB']=sysB_N_transport_out
    nitrogen_dict['trans_sol_loss']['sysB']=sysB_N_transport_loss
    nitrogen_dict['pretreat_sol']['sysB']=sysB_N_pretreatment_solid
    nitrogen_dict['dryer_sol']['sysB']=sysB_N_dryer_out
    nitrogen_dict['dryer_gas']['sysB']=sysB_N_dryer_fugitive_out
    nitrogen_dict['pyrolysis_biochar']['sysB']=sysB_N_pyrolysis_biochar
    nitrogen_dict['pyrolysis_gas']['sysB']=sysB_N_pyrolysis_gas
    nitrogen_dict['toilet_liq']['sysB']=sysB_N_toilet_liq_out
    nitrogen_dict['trans_liq']['sysB']=sysB_N_liq_transport_out
    nitrogen_dict['trans_liq_loss']['sysB']=sysA_N_liq_transport_loss
    nitrogen_dict['pretreat_liq']['sysB']=0.0
    nitrogen_dict['treat_liq']['sysB']=sysA_N_liq_transport_loss
    nitrogen_dict['bed_liq']['sysB']=sysB_N_liq_bed_out
    nitrogen_dict['bed_gas']['sysB']=sysB_N_liq_bed_gas
    nitrogen_dict['struvite']['sysB']=sysB_N_liq_struvite_recovered
    nitrogen_dict['NH3']['sysB']=sysB_N_liq_recovered
    
    phosphorus_dict['urine']['sysB']=sysB_P_urine
    phosphorus_dict['feces']['sysB']=sysB_P_feces
    phosphorus_dict['toilet_sol']['sysB']=sysB_P_toilet_solid_out
    phosphorus_dict['trans_sol']['sysB']=sysB_P_transport_out
    phosphorus_dict['trans_sol_loss']['sysB']=sysB_P_transport_loss
    phosphorus_dict['pretreat_sol']['sysB']=sysB_P_pretreatment_solid
    phosphorus_dict['dryer_sol']['sysB']=sysB_P_dryer_out
    phosphorus_dict['pyrolysis_biochar']['sysB']=sysB_P_pyrolysis_biochar
    phosphorus_dict['pyrolysis_gas']['sysB']=sysB_P_pyrolysis_gas
    phosphorus_dict['toilet_liq']['sysB']=sysB_P_toilet_liq_out
    phosphorus_dict['trans_liq']['sysB']=sysB_P_liq_transport_out
    phosphorus_dict['trans_liq_loss']['sysB']=sysA_P_liq_transport_loss
    phosphorus_dict['pretreat_liq']['sysB']=0.0
    phosphorus_dict['treat_liq']['sysB']=sysA_P_liq_transport_loss
    phosphorus_dict['bed_liq']['sysB']=sysB_P_liq_bed_out
    phosphorus_dict['struvite']['sysB']=sysB_P_liq_struvite_recovered
    
    potassium_dict['urine']['sysB']=sysB_K_urine
    potassium_dict['feces']['sysB']=sysB_K_feces
    potassium_dict['toilet_sol']['sysB']=sysB_K_toilet_solid_out
    potassium_dict['trans_sol']['sysB']=sysB_K_transport_out
    potassium_dict['trans_sol_loss']['sysB']=sysB_K_transport_loss
    potassium_dict['pretreat_sol']['sysB']=sysB_K_pretreatment_solid
    potassium_dict['dryer_sol']['sysB']=sysB_K_dryer_out
    potassium_dict['pyrolysis_biochar']['sysB']=sysB_K_pyrolysis_biochar
    potassium_dict['pyrolysis_gas']['sysB']=sysB_K_pyrolysis_gas
    potassium_dict['toilet_liq']['sysB']=sysB_K_toilet_liq_out
    potassium_dict['trans_liq']['sysB']=sysB_K_liq_transport_out
    potassium_dict['trans_liq_loss']['sysB']=sysA_K_liq_transport_loss
    potassium_dict['pretreat_liq']['sysB']=0.0
    potassium_dict['treat_liq']['sysB']=sysA_K_liq_transport_loss
    potassium_dict['bed_liq']['sysB']=sysB_K_liq_bed_out
    
B17.specification = update_C_recoveries
B17.line = 'fugitive N2O mixer'

############### Simulation, TEA, and LCA ###############
sysB = bst.System('sysB', path= (B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, B14, B15, B16, B17))

sysB.simulate()




teaB = SimpleTEA(system=sysB, discount_rate=get_discount_rate(), 
                  start_year=2020, lifetime=20, uptime_ratio=1, 
                  lang_factor=None, annual_maintenance=0, 
                  annual_labor=(get_operator_daily_wage() * 3*365))

def get_powerB():
    powerB = sum([(u.power_utility.rate*u.uptime_ratio) for u in sysB.units])
    return powerB
    

lcaB = LCA(system=sysB, lifetime=20, lifetime_unit='yr', uptime_ratio=1,
            e_item=lambda: (get_powerB()*365*12*20))



# %%

# =============================================================================
# Scenario C (sysC): pit latrine with Biogenic Refinery 10k users
# =============================================================================

flowsheetC = bst.Flowsheet('sysC')
bst.main_flowsheet.set_flowsheet(flowsheetC)

streamsC = batch_create_streams('C')

#################### Human Inputs ####################
# !!! how to change excrestion based on location (e.g., calorie and protein)
C1 = su.Excretion('C1', outs=('urine','feces'))

################### User Interface ###################
# !!! how to change inputs based on location (e.g., flushing water, cleaning water, toilet paper)
C2 = su.PitLatrine('C2', ins=(C1-0, C1-1,
                              'toilet_paper', 'flushing_water',
                              'cleansing_water', 'desiccant'),
                   outs=('mixed_waste', 'leachate', 'C2_CH4', 'C2_N2O'),
                   N_user=get_toilet_user(), N_toilet=get_ppl('10k')/get_toilet_user(),
                   if_flushing=False, if_desiccant=False, if_toilet_paper=False,
                    OPEX_over_CAPEX=0.05, lifetime=10,
                  decay_k_COD=get_decay_k(tau_deg, log_deg),
                   decay_k_N=get_decay_k(tau_deg, log_deg),
                   max_CH4_emission=get_max_CH4_emission())
 
def update_C2_param():
    C2.N_user = get_toilet_user()
    C2.N_toilet = get_ppl('10k')/get_toilet_user()
    C2._run()
C2.specification = update_C2_param

##################### Conveyance of Waste #####################

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

###################### Treatment ######################

C4 = su.IndustrialControlPanel('C4', ins=C3-0, outs='C4_in')
C5 = su.HousingBiogenicRefinery('C5', ins=C4-0, outs='C5_in', const_wage=get_const_daily_wage(), 
                                const_person_days=const_person_days)

C6 = su.ScrewPress('C6', ins=(C5-0, streamsC['polymer']), outs=('liq', 'cake_sol'))
C7 = su.LiquidTreatmentBed('C7', ins=C6-0, outs=('liquid_bed_treated', 'C7_CH4', 'C7_N2O'),
                           decay_k_COD=get_decay_k(tau_deg, log_deg),
                           decay_k_N=get_decay_k(tau_deg, log_deg),
                           max_CH4_emission=get_max_CH4_emission())

# agricultural residue flowrate is proportional to sysA effluent 
# ag_res_kg_hr = A6.outs[1].imass['OtherSS'] * (2/12)
ag_res_kg_hr = (1-(A6.outs[1].imass['H2O']/A6.outs[1].F_mass)) *(A6.outs[1].F_mass) * (2/12)*(1/0.95)
rice_husk = qs.WasteStream('RiceHusk', RiceHusk=ag_res_kg_hr, units='kg/hr', price=0) 

C_ag_mix = su.Mixer('C_ag_mix', ins=(C6-1, rice_husk), outs=('C_ag_mix_out'))

C8 = su.CarbonizerBase('C8', outs = (streamsC['biochar'], 'C8_hot_gas', 'C8_N2O'))
C9 = su.PollutionControlDevice('C9', ins = (C8-1, C8-2), outs = ('C9_hot_gas_pcd', 'C9_N2O'))

# updating uptime_ratio in all units to follow carbonizer base
C9_old_cost = C9._cost
def update_C9_uptime_ratio():
    C12.uptime_ratio = C11.uptime_ratio = C10.uptime_ratio = C9.uptime_ratio = C8.uptime_ratio 
    C9_old_cost()
C9._cost = update_C9_uptime_ratio

C10 = su.OilHeatExchanger('C10', ins = C9-0, outs = ('C10_hot_gas'))
C11 = su.HydronicHeatExchanger('C11', ins = C10-0, outs = ('C11_hot_gas'))
C12 = su.DryerFromHHX('C12', ins = (C_ag_mix-0, C11-0), outs = ('waste_out', 'C12_N2O', 'C12_CH4'))

C12-0-C8

C13 = su.Mixer('C13', ins=(C2-2, C7-1, C12-2), outs=streamsC['CH4'])
C13.specification = lambda: add_fugitive_items(C13, 'CH4_item')
C13.line = 'fugitive CH4 mixer' 
        
C14 = su.Mixer('C14', ins=(C2-3, C7-2, C9-1, C12-1), outs=streamsC['N2O'])
C14.specification = lambda: add_fugitive_items(C14, 'N2O_item')
C14.line = 'fugitive N2O mixer'

############### Simulation, TEA, and LCA ###############
sysC = bst.System('sysC', path= (C1, C2, C3, C4, C5, C6, C7, C_ag_mix, C8, C9, C10, C11, C12, C13, C14))

sysC.simulate()


teaC = SimpleTEA(system=sysC, discount_rate=get_discount_rate(), 
                  start_year=2020, lifetime=20, uptime_ratio=1, 
                  lang_factor=None, annual_maintenance=0, 
                  annual_labor=(get_operator_daily_wage() * 3*365))

def get_powerC():
    powerC = sum([(u.power_utility.rate*u.uptime_ratio) for u in sysC.units])
    return powerC

lcaC = LCA(system=sysC, lifetime=20, lifetime_unit='yr', uptime_ratio=1,
            e_item=lambda: (get_powerC()))
#%%

sysTEST = bst.System('sysTEST', path=(C6,))
sysTEST.simulate()

teaTEST = SimpleTEA(system=sysTEST, discount_rate=get_discount_rate(), 
                  start_year=2020, lifetime=20, uptime_ratio=1, 
                  lang_factor=None, annual_maintenance=0, 
                  annual_labor=(get_operator_daily_wage() * 3*365))

def get_powerTEST():
    powerTEST = sum([(u.power_utility.rate*u.uptime_ratio) for u in sysTEST.units])
    return powerTEST

lcaTEST = LCA(system=sysTEST, lifetime=20, lifetime_unit='yr', uptime_ratio=1,
            e_item=lambda: (get_powerC()*(365*18)*20))

# %%
# =============================================================================
# Scenario D (sysD): pit latrine with 12k users
# =============================================================================

flowsheetD = bst.Flowsheet('sysD')
bst.main_flowsheet.set_flowsheet(flowsheetD)
streamsD = batch_create_streams('D')

#################### Human Inputs ####################

D1 = su.Excretion('D1', outs=('urine','feces'))

################### User Interface ###################
# !!! how to change inputs based on location (e.g., flushing water, cleaning water, toilet paper)
D2 = su.PitLatrine('D2', ins=(D1-0, D1-1,
                              'toilet_paper', 'flushing_water',
                              'cleansing_water', 'desiccant'),
                    outs=('mixed_waste', 'leachate', 'A2_CH4', 'A2_N2O'),
                    N_user=get_toilet_user(), N_toilet=get_ppl('12k')/get_toilet_user(),
                    if_flushing=False, if_desiccant=False, if_toilet_paper=False,
                    OPEX_over_CAPEX=0.05, lifetime=5,
                  decay_k_COD=get_decay_k(tau_deg, log_deg),
                    decay_k_N=get_decay_k(tau_deg, log_deg),
                    max_CH4_emission=get_max_CH4_emission()
                    )
 
def update_D2_param():
    D2.N_user = get_toilet_user()
    D2.N_toilet = get_ppl('12k')/get_toilet_user()
    D2._run()
D2.specification = update_D2_param

##################### Conveyance of Waste #####################
D3 = su.Trucking('D3', ins=D2-0, outs=('transported', 'conveyance_loss'),
                  load_type='mass', distance=5, distance_unit='km',
                  interval=A2.emptying_period, interval_unit='yr',
                  loss_ratio=0.02)
def update_D3_param():
    D3._run()
    truck = D3.single_truck
    truck.interval = D2.emptying_period*365*24
    truck.load = D3.F_mass_in*truck.interval/D2.N_toilet
    rho = D3.F_mass_in/A3.F_vol_in
    vol = truck.load/rho
    D3.fee = get_tanker_truck_fee(vol)
    D3._design()
D3.specification = update_D3_param

D4 = su.Lagoon('D4', ins=D3-0, outs=('anaerobic_treated', 'D4_CH4', 'D4_N2O'),
                design_type='anaerobic',
                decay_k_N=get_decay_k(tau_deg, log_deg),
                max_CH4_emission=get_max_CH4_emission())


D5 = su.Mixer('D5', ins=(D2-2, D4-1), outs=streamsD['CH4'])
D5.specification = lambda: add_fugitive_items(D5, 'CH4_item')
D5.line = 'fugitive CH4 mixer' 
        
D6 = su.Mixer('D6', ins=(D2-3, D4-2), outs=streamsD['N2O'])
D6.specification = lambda: add_fugitive_items(D6, 'N2O_item')
D6.line = 'fugitive N2O mixer'

sysD = bst.System('sysD', path=(D1, D2, D3, D4, D5, D6))
             
sysD.simulate()


teaD = SimpleTEA(system=sysD, discount_rate=get_discount_rate(), 
                  start_year=2020, lifetime=20, uptime_ratio=1, 
                  lang_factor=None, annual_maintenance=0)

def get_powerD():
    powerD = sum([(u.power_utility.rate*u.uptime_ratio) for u in sysD.units])
    return powerD

lcaD = LCA(system=sysD, lifetime=20, lifetime_unit='yr', uptime_ratio=1,
            e_item=lambda: (get_powerD()))

# %%

# =============================================================================
# Summarizing Functions
# =============================================================================

#!!! items need to be updated in sys_dct, system_streams, learning curve assumptions
sys_dct = {
    'ppl': dict(sysA=get_ppl('12k'), sysB=get_ppl('12k'), sysC=get_ppl('10k'), sysD=get_ppl('12k')),
    'input_unit': dict(sysA=A1, sysB=B1, sysC=C1, sysD=D1),
    'liq_unit': dict(sysA=None, sysB=None, sysC=None, sysD=None),
    'sol_unit': dict(sysA=None, sysB=None, sysC=None, sysD=None),
    'gas_unit': dict(sysA=None, sysB=None, sysC=None, sysD=None),
    'stream_dct': dict(sysA=streamsA, sysB=streamsB, sysC=streamsC, sysD=streamsD),
    'TEA': dict(sysA=teaA, sysB=teaB, sysC=teaC, sysD=teaD),
    'LCA': dict(sysA=lcaA, sysB=lcaB, sysC=lcaC, sysD=lcaD),
    'cache': dict(sysA={}, sysB={}, sysC={}, sysD={}),
    }

unit_dct = {
    'front_end': dict(sysA=(A2,), sysB=(B2,), sysC=None, sysD=None),
    'transport': dict(sysA=(A3,), sysB=(B3,B4,), sysC=None, sysD=None),
    'pretreatment': dict(sysA=(A6,), sysB=(B10,), sysC=None, sysD=None),
    'liq_treatment': dict(sysA=(A7,), sysB=(B5,B6,B7), sysC=None, sysD=None),
    'controls_housing': dict(sysA=(A4,A5), sysB=(B8,B9), sysC=None, sysD=None),
    'dryer_hhx': dict(sysA=(A11,A12), sysB=(B14,B15), sysC=None, sysD=None),
    'ohx': dict(sysA=(A10,), sysB=(B13,), sysC=None, sysD=None),
    'carbonizer_base': dict(sysA=(A8,), sysB=(B11,), sysC=None, sysD=None),
    'pollution_control': dict(sysA=(A9,), sysB=(B12,), sysC=None, sysD=None)
    }

consumables_dct = {
    }

system_streams = {sysA: streamsA, sysB: streamsB, sysC: streamsC, sysD: streamsD}

#learning curve assumptions
percent_CAPEX_to_scale = 0.1
get_percent_CAPEX_to_scale = lambda: percent_CAPEX_to_scale

number_of_units = 10000

percent_limit = 0.01
get_percent_limit = lambda: percent_limit

learning_curve_percent=0.9
get_learning_curve_percent = lambda: learning_curve_percent

def get_scaled_capital(tea, sys):
    if sys == sysD:
        new_CAPEX_annualized=tea.annualized_CAPEX
    else:
        CAPEX_to_scale = tea.annualized_CAPEX * get_percent_CAPEX_to_scale()
        CAPEX_not_scaled = tea.annualized_CAPEX - CAPEX_to_scale
        scaled_limited = CAPEX_to_scale * get_percent_limit()
        b=(np.log(get_learning_curve_percent())/np.log(2))
        scaled_CAPEX_annualized  = (CAPEX_to_scale - scaled_limited)*number_of_units**b+scaled_limited    
        new_CAPEX_annualized = scaled_CAPEX_annualized + CAPEX_not_scaled    
    return new_CAPEX_annualized

def get_cost_capex(unit,tea,ppl):
    # capex = 0
    # if type(unit)== tuple:
    #     for i in unit:
    #         capex+=tea.get_unit_annualized_equipment_cost(i)/365/ppl
    return tea.get_unit_annualized_equipment_cost(unit)/365/ppl

def get_ghg_capex(unit,lca,ppl):
    capex = 0
    if type(unit)== tuple:
        for i in unit:
            capex+=lca.get_construction_impacts(i).get('GlobalWarming')/lca.lifetime/ppl
    return capex

def get_cost_opex(unit,ppl):
    opex = 0
    if type(unit)== tuple:
        for i in unit:
            opex+=sum(i._add_OPEX.values())*24/ppl
            for stream in i.ins:
                opex+=stream.imass['Polyacrylamide']*streamsA['polymer'].price*24/ppl
                opex+=stream.imass['MagnesiumHydroxide']*streamsB['MgOH2'].price*24/ppl
                opex+=stream.imass['MgCO3']*streamsB['MgCO3'].price*24/ppl
                opex+=stream.imass['H2SO4']*streamsB['H2SO4'].price*24/ppl
                opex+=stream.imass['FilterBag']*streamsB['filter_bag'].price*24/ppl
                opex+=stream.imass['Polystyrene']*streamsB['resin'].price*24/ppl
    return opex

def get_cost_electricity(unit,ppl):
    electricity = 0
    if type(unit)== tuple:
        for i in unit:
            electricity+=i.power_utility.cost*24/ppl
    return electricity
        
def get_ghg_electricity(unit,ppl):
    electricity = 0
    if type(unit)== tuple:
        for i in unit:
            electricity+=i.power_utility.consumption*i.uptime_ratio*e_item.CFs['GlobalWarming']*12*365/ppl
            electricity+=-i.power_utility.production*i.uptime_ratio*e_item.CFs['GlobalWarming']*12*365/ppl
    return electricity

def get_ghg_direct(unit,lca,ppl):
    direct = 0
    if type(unit)== tuple:
        for i in unit:
            for stream in i.outs:
                direct+=stream.imass['CH4']*CH4_item.CFs['GlobalWarming']*24*365/ppl
                direct+=stream.imass['N2O']*N2O_item.CFs['GlobalWarming']*24*365/ppl
    return direct

def get_ghg_opex(unit,lca,ppl):
    opex = 0
    if type(unit)== tuple:
        for i in unit:
            opex+=lca.get_stream_impacts(i.ins).get('GlobalWarming')/20/ppl
    return opex

def get_ghg_transport(unit,lca,ppl):
    transport = 0
    if type(unit)== tuple:
        for i in unit:
            transport+=lca.get_transportation_impacts(i).get('GlobalWarming')/20/ppl     
    return transport

def get_total_inputs(unit):
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
    hr = 365 * 24
    for i, j in inputs.items():
        inputs[i] = j * hr
    return inputs

def get_recovery(ins=None, outs=None, hr=365*24, ppl=1, if_relative=True):
    try: iter(outs)
    except: outs = (outs,)
    non_g = tuple(i for i in outs if i.phase != 'g')
    recovery = {}
    recovery['COD'] = sum(i.COD*i.F_vol/1e3 for i in non_g)
    recovery['N'] = sum(i.TN*i.F_vol/1e3 for i in non_g)
    recovery['NH3'] = sum(i.imass['NH3'] for i in non_g)
    recovery['P'] = sum(i.TP*i.F_vol/1e3 for i in non_g)
    recovery['K'] = sum(i.TK*i.F_vol/1e3 for i in non_g)
    for i, j in recovery.items():
        if if_relative:
            inputs = get_total_inputs(ins)
            recovery[i] /= inputs[i]/hr * ppl
        else:
            recovery[i] /= 1/hr * ppl
    return recovery

def get_stream_emissions(streams=None, hr=365*24, ppl=1):
    try: iter(streams)
    except: streams = (streams,)
    emission = {}
    factor = hr / ppl
    for i in streams:
        if not i.stream_impact_item: continue
        emission[f'{i.ID}'] = i.F_mass*i.stream_impact_item.CFs['GlobalWarming']*factor
    return emission


def get_fugitive_emissions(streams=None, hr=365*24*20):
    try: iter(streams)
    except: streams = (streams,)
    factor = hr
    emission=0
    for i in streams:
        if not i.stream_impact_item: continue
        emission = i.F_mass*i.stream_impact_item.CFs['GlobalWarming']*factor
    return emission

#N2O and CH4 emissions
def get_fugitive_gases(system_num=None):
    fugitive_streams = ('CH4', 'N2O')
    gases = 0
    for i in fugitive_streams:
        gases+=get_fugitive_emissions(system_streams[system_num][i])
    return gases

def cache_recoveries(sys):
    total_COD = get_total_inputs(sys_dct['input_unit'][sys.ID])['COD']
    ppl = sys_dct['ppl'][sys.ID]
    if sys_dct['gas_unit'][sys.ID]:
        gas_mol = sys_dct['gas_unit'][sys.ID].outs[0].imol['CH4']
        gas_COD = gas_mol*1e3*get_biogas_energy()*365*24/14e3/ppl/total_COD
    else:
        gas_COD = 0
    cache = {
        'liq': get_recovery(ins=sys_dct['input_unit'][sys.ID],
                            outs=sys_dct['liq_unit'][sys.ID].ins,
                            ppl=ppl),
        'sol': get_recovery(ins=sys_dct['input_unit'][sys.ID],
                            outs=sys_dct['sol_unit'][sys.ID].ins,
                            ppl=ppl),
        'gas': dict(COD=gas_COD, N=0, P=0, K=0)
        }
    return cache

def update_cache(sys):
    last_u = sys.path[-1]
    last_u._run()
    sys_dct['cache'][sys.ID] = cache_recoveries(sys)

def get_summarizing_functions():
    func_dct = {}
    func_dct['get_annual_cost'] = lambda tea, ppl, sys: (tea.EAC-tea.annualized_CAPEX
                                                    +get_scaled_capital(tea,sys))/ppl
    func_dct['get_annual_CAPEX'] = lambda tea, ppl, sys: get_scaled_capital(tea,sys)/ppl
    func_dct['get_annual_energy'] = lambda tea, ppl: sum(i.power_utility.cost for i in tea.units)*tea.operating_hours/ppl
    func_dct['get_annual_OPEX'] = lambda tea, ppl: (tea.AOC/ppl -(sum(i.power_utility.cost for i in tea.units)
                                                    *tea.operating_hours/ppl) - (tea.annual_labor/ppl))
    func_dct['get_annual_labor'] = lambda tea, ppl: tea.annual_labor/ppl
    func_dct['get_annual_sales'] = lambda tea, ppl: (tea.sales/ppl)
    ind = 'GlobalWarming'
    func_dct['get_annual_GWP'] = \
        lambda lca, ppl: lca.total_impacts[ind]/lca.lifetime/ppl
    func_dct['get_constr_GWP'] = \
        lambda lca, ppl: lca.total_construction_impacts[ind]/lca.lifetime/ppl
    func_dct['get_trans_GWP'] = \
        lambda lca, ppl: lca.total_transportation_impacts[ind]/lca.lifetime/ppl  
    func_dct['get_CH4_N2O_GWP'] = \
        lambda sys, lca, ppl: get_fugitive_gases(sys)/lca.lifetime/ppl
    func_dct['get_stream_items_emission_GWP'] = \
        lambda sys, lca, ppl: (lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='direct_emission')[ind]) \
            /lca.lifetime/ppl - (get_fugitive_gases(sys)/lca.lifetime/ppl)
    func_dct['get_offset_GWP'] = \
        lambda lca, ppl: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='offset')[ind] \
            /lca.lifetime/ppl
    func_dct['get_other_GWP'] = \
        lambda lca, ppl: lca.total_other_impacts[ind]/lca.lifetime/ppl
    
    func_dct['C_urine'] = lambda sys: carbon_dict['urine'][sys.ID]
    func_dct['C_feces'] = lambda sys: carbon_dict['feces'][sys.ID]
    func_dct['C_toilet_sol'] = lambda sys: carbon_dict['toilet_sol'][sys.ID]
    func_dct['C_toilet_gas'] = lambda sys: carbon_dict['toilet_gas'][sys.ID]
    func_dct['C_trans_sol'] = lambda sys: carbon_dict['trans_sol'][sys.ID]
    func_dct['C_trans_sol_loss'] = lambda sys: carbon_dict['trans_sol_loss'][sys.ID]
    func_dct['C_pretreat_sol'] = lambda sys: carbon_dict['pretreat_sol'][sys.ID]
    func_dct['C_dryer_sol'] = lambda sys: carbon_dict['dryer_sol'][sys.ID]
    func_dct['C_dryer_gas'] = lambda sys: carbon_dict['dryer_gas'][sys.ID]
    func_dct['C_pyrolysis_biochar'] = lambda sys: carbon_dict['pyrolysis_biochar'][sys.ID]
    func_dct['C_pyrolysis_gas'] = lambda sys: carbon_dict['pyrolysis_gas'][sys.ID]
    func_dct['C_toilet_liq'] = lambda sys: carbon_dict['toilet_liq'][sys.ID]
    func_dct['C_trans_liq'] = lambda sys: carbon_dict['trans_liq'][sys.ID]
    func_dct['C_trans_liq_loss'] = lambda sys: carbon_dict['trans_liq_loss'][sys.ID]
    func_dct['C_pretreat_liq'] = lambda sys: carbon_dict['pretreat_liq'][sys.ID]
    func_dct['C_treat_liq'] = lambda sys: carbon_dict['treat_liq'][sys.ID]
    func_dct['C_bed_liq'] = lambda sys: carbon_dict['bed_liq'][sys.ID]
    func_dct['C_bed_gas'] = lambda sys: carbon_dict['bed_gas'][sys.ID]
    
    func_dct['N_urine'] = lambda sys: nitrogen_dict['urine'][sys.ID]
    func_dct['N_feces'] = lambda sys: nitrogen_dict['feces'][sys.ID]
    func_dct['N_toilet_sol'] = lambda sys: nitrogen_dict['toilet_sol'][sys.ID]
    func_dct['N_toilet_gas'] = lambda sys: nitrogen_dict['toilet_gas'][sys.ID]
    func_dct['N_trans_sol'] = lambda sys: nitrogen_dict['trans_sol'][sys.ID]
    func_dct['N_trans_sol_loss'] = lambda sys: nitrogen_dict['trans_sol_loss'][sys.ID]
    func_dct['N_pretreat_sol'] = lambda sys: nitrogen_dict['pretreat_sol'][sys.ID]
    func_dct['N_dryer_sol'] = lambda sys: nitrogen_dict['dryer_sol'][sys.ID]
    func_dct['N_dryer_gas'] = lambda sys: nitrogen_dict['dryer_gas'][sys.ID]
    func_dct['N_pyrolysis_biochar'] = lambda sys: nitrogen_dict['pyrolysis_biochar'][sys.ID]
    func_dct['N_pyrolysis_gas'] = lambda sys: nitrogen_dict['pyrolysis_gas'][sys.ID]
    func_dct['N_toilet_liq'] = lambda sys: nitrogen_dict['toilet_liq'][sys.ID]
    func_dct['N_trans_liq'] = lambda sys: nitrogen_dict['trans_liq'][sys.ID]
    func_dct['N_trans_liq_loss'] = lambda sys: nitrogen_dict['trans_liq_loss'][sys.ID]
    func_dct['N_pretreat_liq'] = lambda sys: nitrogen_dict['pretreat_liq'][sys.ID]
    func_dct['N_treat_liq'] = lambda sys: nitrogen_dict['treat_liq'][sys.ID]
    func_dct['N_bed_liq'] = lambda sys: nitrogen_dict['bed_liq'][sys.ID]
    func_dct['N_bed_gas'] = lambda sys: nitrogen_dict['bed_gas'][sys.ID]
    func_dct['N_struvite'] = lambda sys: nitrogen_dict['struvite'][sys.ID]
    func_dct['N_NH3'] = lambda sys: nitrogen_dict['NH3'][sys.ID]
    
    func_dct['P_urine'] = lambda sys: phosphorus_dict['urine'][sys.ID]
    func_dct['P_feces'] = lambda sys: phosphorus_dict['feces'][sys.ID]
    func_dct['P_toilet_sol'] = lambda sys: phosphorus_dict['toilet_sol'][sys.ID]
    func_dct['P_toilet_gas'] = lambda sys: phosphorus_dict['toilet_gas'][sys.ID]
    func_dct['P_trans_sol'] = lambda sys: phosphorus_dict['trans_sol'][sys.ID]
    func_dct['P_trans_sol_loss'] = lambda sys: phosphorus_dict['trans_sol_loss'][sys.ID]
    func_dct['P_pretreat_sol'] = lambda sys: phosphorus_dict['pretreat_sol'][sys.ID]
    func_dct['P_dryer_sol'] = lambda sys: phosphorus_dict['dryer_sol'][sys.ID]
    func_dct['P_dryer_gas'] = lambda sys: phosphorus_dict['dryer_gas'][sys.ID]
    func_dct['P_pyrolysis_biochar'] = lambda sys: phosphorus_dict['pyrolysis_biochar'][sys.ID]
    func_dct['P_pyrolysis_gas'] = lambda sys: phosphorus_dict['pyrolysis_gas'][sys.ID]
    func_dct['P_toilet_liq'] = lambda sys: phosphorus_dict['toilet_liq'][sys.ID]
    func_dct['P_trans_liq'] = lambda sys: phosphorus_dict['trans_liq'][sys.ID]
    func_dct['P_trans_liq_loss'] = lambda sys: phosphorus_dict['trans_liq_loss'][sys.ID]
    func_dct['P_pretreat_liq'] = lambda sys: phosphorus_dict['pretreat_liq'][sys.ID]
    func_dct['P_treat_liq'] = lambda sys: phosphorus_dict['treat_liq'][sys.ID]
    func_dct['P_bed_liq'] = lambda sys: phosphorus_dict['bed_liq'][sys.ID]
    func_dct['P_bed_gas'] = lambda sys: phosphorus_dict['bed_gas'][sys.ID]
    func_dct['P_struvite'] = lambda sys: phosphorus_dict['struvite'][sys.ID]
    
    func_dct['K_urine'] = lambda sys: potassium_dict['urine'][sys.ID]
    func_dct['K_feces'] = lambda sys: potassium_dict['feces'][sys.ID]
    func_dct['K_toilet_sol'] = lambda sys: potassium_dict['toilet_sol'][sys.ID]
    func_dct['K_toilet_gas'] = lambda sys: potassium_dict['toilet_gas'][sys.ID]
    func_dct['K_trans_sol'] = lambda sys: potassium_dict['trans_sol'][sys.ID]
    func_dct['K_trans_sol_loss'] = lambda sys: potassium_dict['trans_sol_loss'][sys.ID]
    func_dct['K_pretreat_sol'] = lambda sys: potassium_dict['pretreat_sol'][sys.ID]
    func_dct['K_dryer_sol'] = lambda sys: potassium_dict['dryer_sol'][sys.ID]
    func_dct['K_dryer_gas'] = lambda sys: potassium_dict['dryer_gas'][sys.ID]
    func_dct['K_pyrolysis_biochar'] = lambda sys: potassium_dict['pyrolysis_biochar'][sys.ID]
    func_dct['K_pyrolysis_gas'] = lambda sys: potassium_dict['pyrolysis_gas'][sys.ID]
    func_dct['K_toilet_liq'] = lambda sys: potassium_dict['toilet_liq'][sys.ID]
    func_dct['K_trans_liq'] = lambda sys: potassium_dict['trans_liq'][sys.ID]
    func_dct['K_trans_liq_loss'] = lambda sys: potassium_dict['trans_liq_loss'][sys.ID]
    func_dct['K_pretreat_liq'] = lambda sys: potassium_dict['pretreat_liq'][sys.ID]
    func_dct['K_treat_liq'] = lambda sys: potassium_dict['treat_liq'][sys.ID]
    func_dct['K_bed_liq'] = lambda sys: potassium_dict['bed_liq'][sys.ID]
    func_dct['K_bed_gas'] = lambda sys: potassium_dict['bed_gas'][sys.ID]
    
    #Cost broken down by units
    #operator
    func_dct['operator_cost_opex'] = lambda tea,ppl: tea.annual_labor/365/ppl
    
    #front_end
    func_dct['front_end_cost_capex'] = lambda front_end,tea,ppl: get_cost_capex(front_end,tea,ppl)
    func_dct['front_end_cost_opex'] = lambda front_end,ppl: get_cost_opex(front_end,ppl)
    func_dct['front_end_cost_electricity'] = lambda front_end,ppl: get_cost_electricity(front_end,ppl)
    
    #transport
    func_dct['transport_cost_opex'] = lambda transport,ppl: get_cost_opex(transport,ppl)

    #pretreatment
    func_dct['pretreatment_cost_capex'] = lambda pretreatment,tea,ppl: get_cost_capex(pretreatment,tea,ppl)
    func_dct['pretreatment_cost_opex'] = lambda pretreatment,ppl: get_cost_opex(pretreatment,ppl)
    func_dct['pretreatment_cost_electricity'] = lambda pretreatment,ppl: get_cost_electricity(pretreatment,ppl)
    
    #liq_treatment
    func_dct['liq_treatment_cost_capex'] = lambda liq_treatment,tea,ppl: get_cost_capex(liq_treatment,tea,ppl)
    func_dct['liq_treatment_cost_opex'] = lambda liq_treatment,ppl: get_cost_opex(liq_treatment,ppl)
    func_dct['liq_treatment_cost_electricity'] = lambda liq_treatment,ppl: get_cost_electricity(liq_treatment,ppl)
    
    #controls_housing
    func_dct['controls_housing_cost_capex'] = lambda controls_housing,tea,ppl: get_cost_capex(controls_housing,tea,ppl)
    func_dct['controls_housing_cost_opex'] = lambda controls_housing,ppl: get_cost_opex(controls_housing,ppl)
    func_dct['controls_housing_cost_electricity'] = lambda controls_housing,ppl: get_cost_electricity(controls_housing,ppl)
    
    #dryer_hhx
    func_dct['dryer_hhx_cost_capex'] = lambda dryer_hhx,tea,ppl: get_cost_capex(dryer_hhx,tea,ppl)
    func_dct['dryer_hhx_cost_opex'] = lambda dryer_hhx,ppl: get_cost_opex(dryer_hhx,ppl)
    func_dct['dryer_hhx_cost_electricity'] = lambda dryer_hhx,ppl: get_cost_electricity(dryer_hhx,ppl)
    
    #ohx
    func_dct['ohx_cost_capex'] = lambda ohx,tea,ppl: get_cost_capex(ohx,tea,ppl)
    func_dct['ohx_cost_opex'] = lambda ohx,ppl: get_cost_opex(ohx,ppl)
    func_dct['ohx_cost_electricity'] = lambda ohx,ppl: get_cost_electricity(ohx,ppl)
    
    #carbonizer_base
    func_dct['carbonizer_base_cost_capex'] = lambda carbonizer_base,tea,ppl: get_cost_capex(carbonizer_base,tea,ppl)
    func_dct['carbonizer_base_cost_opex'] = lambda carbonizer_base,ppl: get_cost_opex(carbonizer_base,ppl)
    func_dct['carbonizer_base_cost_electricity'] = lambda carbonizer_base,ppl: get_cost_electricity(carbonizer_base,ppl)
    
    #pollution_control
    func_dct['pollution_control_cost_capex'] = lambda pollution_control,tea,ppl: get_cost_capex(pollution_control,tea,ppl)
    func_dct['pollution_control_cost_opex'] = lambda pollution_control,ppl: get_cost_opex(pollution_control,ppl)
    func_dct['pollution_control_cost_electricity'] = lambda pollution_control,ppl: get_cost_electricity(pollution_control,ppl)
    
    #GHG broken down by units
    #front_end
    func_dct['front_end_ghg_capex'] = lambda front_end,lca,ppl: get_ghg_capex(front_end,lca,ppl)
    func_dct['front_end_ghg_opex'] = lambda front_end,lca,ppl: get_ghg_opex(front_end,lca,ppl)
    func_dct['front_end_ghg_electricity'] = lambda front_end,ppl: get_ghg_electricity(front_end,ppl)
    func_dct['front_end_ghg_direct'] = lambda front_end,lca,ppl: get_ghg_direct(front_end,lca,ppl)
    
    #transport
    func_dct['transport_ghg_opex'] = lambda transport,lca,ppl: get_ghg_transport(transport,lca,ppl)
    
    #pretreatment
    func_dct['pretreatment_ghg_capex'] = lambda pretreatment,lca,ppl: get_ghg_capex(pretreatment,lca,ppl)
    func_dct['pretreatment_ghg_opex'] = lambda pretreatment,lca,ppl: get_ghg_opex(pretreatment,lca,ppl)
    func_dct['pretreatment_ghg_electricity'] = lambda pretreatment,ppl: get_ghg_electricity(pretreatment,ppl)
    func_dct['pretreatment_ghg_direct'] = lambda pretreatment,lca,ppl: get_ghg_direct(pretreatment,lca,ppl)
    
    #liq_treatment
    func_dct['liq_treatment_ghg_capex'] = lambda liq_treatment,lca,ppl: get_ghg_capex(liq_treatment,lca,ppl)
    func_dct['liq_treatment_ghg_opex'] = lambda liq_treatment,lca,ppl: get_ghg_opex(liq_treatment,lca,ppl)
    func_dct['liq_treatment_ghg_electricity'] = lambda liq_treatment,ppl: get_ghg_electricity(liq_treatment,ppl)
    func_dct['liq_treatment_ghg_direct'] = lambda liq_treatment,lca,ppl: get_ghg_direct(liq_treatment,lca,ppl)
    
    #controls_housing
    func_dct['controls_housing_ghg_capex'] = lambda controls_housing,lca,ppl: get_ghg_capex(controls_housing,lca,ppl)
    func_dct['controls_housing_ghg_opex'] = lambda controls_housing,lca,ppl: get_ghg_opex(controls_housing,lca,ppl)
    func_dct['controls_housing_ghg_electricity'] = lambda controls_housing,ppl: get_ghg_electricity(controls_housing,ppl)
    func_dct['controls_housing_ghg_direct'] = lambda controls_housing,lca,ppl: get_ghg_direct(controls_housing,lca,ppl)
    
    #dryer_hhx
    func_dct['dryer_hhx_ghg_capex'] = lambda dryer_hhx,lca,ppl: get_ghg_capex(dryer_hhx,lca,ppl)
    func_dct['dryer_hhx_ghg_opex'] = lambda dryer_hhx,lca,ppl: get_ghg_opex(dryer_hhx,lca,ppl)
    func_dct['dryer_hhx_ghg_electricity'] = lambda dryer_hhx,ppl: get_ghg_electricity(dryer_hhx,ppl)
    func_dct['dryer_hhx_ghg_direct'] = lambda dryer_hhx,lca,ppl: get_ghg_direct(dryer_hhx,lca,ppl)
    
    #ohx
    func_dct['ohx_ghg_capex'] = lambda ohx,lca,ppl: get_ghg_capex(ohx,lca,ppl)
    func_dct['ohx_ghg_opex'] = lambda ohx,lca,ppl: get_ghg_opex(ohx,lca,ppl)
    func_dct['ohx_ghg_electricity'] = lambda ohx,ppl: get_ghg_electricity(ohx,ppl)
    func_dct['ohx_ghg_direct'] = lambda ohx,lca,ppl: get_ghg_direct(ohx,lca,ppl)
    
    #carbonizer_base
    func_dct['carbonizer_base_ghg_capex'] = lambda carbonizer_base,lca,ppl: get_ghg_capex(carbonizer_base,lca,ppl)
    func_dct['carbonizer_base_ghg_opex'] = lambda carbonizer_base,lca,ppl: get_ghg_opex(carbonizer_base,lca,ppl)
    func_dct['carbonizer_base_ghg_electricity'] = lambda carbonizer_base,ppl: get_ghg_electricity(carbonizer_base,ppl)
    func_dct['carbonizer_base_ghg_direct'] = lambda carbonizer_base,lca,ppl: get_ghg_direct(carbonizer_base,lca,ppl)
    
    #pollution_control
    func_dct['pollution_control_ghg_capex'] = lambda pollution_control,lca,ppl: get_ghg_capex(pollution_control,lca,ppl)
    func_dct['pollution_control_ghg_opex'] = lambda pollution_control,lca,ppl: get_ghg_opex(pollution_control,lca,ppl)
    func_dct['pollution_control_ghg_electricity'] = lambda pollution_control,ppl: get_ghg_electricity(pollution_control,ppl)
    func_dct['pollution_control_ghg_direct'] = lambda pollution_control,lca,ppl: get_ghg_direct(pollution_control,lca,ppl)
    
    
    # func_dct['captial_front_end'] = lambda sys: cost_dict['capital']['front_end'][sys.ID]
    # func_dct['opex_transport'] = lambda sys: cost_dict['opex']['transport'][sys.ID]
    # func_dct['captial_pretreatment'] = lambda sys: cost_dict['capital']['pretreatment'][sys.ID]
    # func_dct['opex_pretreatment'] = lambda sys: cost_dict['opex']['pretreatment'][sys.ID]
    # func_dct['electricity_pretreatment'] = lambda sys: cost_dict['electricity']['pretreatment'][sys.ID]
    
    # func_dct['GHG_capital_front_end'] = lambda lca,ppl,front_end: (lca.get_construction_impacts(front_end).get('GlobalWarming')/20/ppl)
    
    
    # for i in ('COD', 'N', 'P', 'K'):
    #     func_dct[f'get_liq_{i}_recovery'] = \
    #         lambda sys, i: sys_dct['cache'][sys.ID]['liq'][i]
    #     func_dct[f'get_sol_{i}_recovery'] = \
    #         lambda sys, i: sys_dct['cache'][sys.ID]['sol'][i]
    #     func_dct[f'get_gas_{i}_recovery'] = \
    #         lambda sys, i: sys_dct['cache'][sys.ID]['gas'][i]
    #     func_dct[f'get_tot_{i}_recovery'] = \
    #         lambda sys, i: \
    #             sys_dct['cache'][sys.ID]['liq'][i] + \
    #             sys_dct['cache'][sys.ID]['sol'][i] + \
    #             sys_dct['cache'][sys.ID]['gas'][i]
    return func_dct


def print_summaries(systems):
    try: iter(systems)
    except: systems = (systems, )
    func = get_summarizing_functions()
    for sys in systems:
        sys.simulate()
        ppl = sys_dct['ppl'][sys.ID]
        print(f'\n---------- Summary for {sys} ----------\n')
        tea = sys_dct['TEA'][sys.ID]
        tea.show()
        print('\n')
        lca = sys_dct['LCA'][sys.ID]
        lca.show()
        front_end=unit_dct['front_end'][sys.ID]
        transport=unit_dct['transport'][sys.ID]
        pretreatment=unit_dct['pretreatment'][sys.ID]
        liq_treatment=unit_dct['liq_treatment'][sys.ID]
        controls_housing=unit_dct['controls_housing'][sys.ID]
        dryer_hhx=unit_dct['dryer_hhx'][sys.ID]
        ohx=unit_dct['ohx'][sys.ID]
        carbonizer_base=unit_dct['carbonizer_base'][sys.ID]
        pollution_control=unit_dct['pollution_control'][sys.ID]
        
        unit = f'{currency}/cap/yr'
        print(f'\nTotal cost: {func["get_annual_cost"](tea, ppl, sys):.2f} {unit}.')
        print(f'Capital: {func["get_annual_CAPEX"](tea, ppl, sys):.2f} {unit}.')
        print(f'Energy: {func["get_annual_energy"](tea,ppl):.2f} {unit}.')
        print(f'Operating: {func["get_annual_OPEX"](tea, ppl):.2f} {unit}.')
        print(f'Labor: {func["get_annual_labor"](tea, ppl):.2f} {unit}.')
        print(f'Sales: {func["get_annual_sales"](tea, ppl):.2f} {unit}.')
        
        unit = f'{GWP.unit}/cap/yr'
        print(f'\nNet emission: {func["get_annual_GWP"](lca, ppl):.2f} {unit}.')
        print(f'Construction: {func["get_constr_GWP"](lca, ppl):.2f} {unit}.')
        print(f'Transportation: {func["get_trans_GWP"](lca, ppl):.2f} {unit}.')
        print(f'Fugitive gas: {func["get_CH4_N2O_GWP"](sys, lca, ppl):.2f} {unit}.')
        print(f'Stream items emission: {func["get_stream_items_emission_GWP"](sys, lca, ppl):.2f} {unit}.')
        print(f'Offset: {func["get_offset_GWP"](lca, ppl):.2f} {unit}.')
        print(f'Other: {func["get_other_GWP"](lca, ppl):.2} {unit}.\n')
        
        unit = f'kg carbon/yr'
        print(f'C urine: {func["C_urine"](sys):.2} {unit}.\n')
        print(f'C feces: {func["C_feces"](sys):.2} {unit}.\n')
        print(f'C toilet sol: {func["C_toilet_sol"](sys):.2} {unit}.\n')
        print(f'C toilet gas: {func["C_toilet_gas"](sys):.2} {unit}.\n')
        print(f'C transport sol: {func["C_trans_sol"](sys):.2} {unit}.\n')
        print(f'C transport sol loss: {func["C_trans_sol_loss"](sys):.2} {unit}.\n')
        print(f'C pretreat sol: {func["C_pretreat_sol"](sys):.2} {unit}.\n')
        print(f'C dryer sol: {func["C_dryer_sol"](sys):.2} {unit}.\n')
        print(f'C dryer gas: {func["C_dryer_gas"](sys):.2} {unit}.\n')
        print(f'C pyrolysis biochar: {func["C_pyrolysis_biochar"](sys):.2} {unit}.\n')
        print(f'C pyrolysis gas: {func["C_pyrolysis_gas"](sys):.2} {unit}.\n')
        print(f'C toilet liq: {func["C_toilet_liq"](sys):.2} {unit}.\n')
        print(f'C trans liq: {func["C_trans_liq"](sys):.2} {unit}.\n')
        print(f'C trans liq loss: {func["C_trans_liq_loss"](sys):.2} {unit}.\n')
        print(f'C pretreat liq: {func["C_pretreat_liq"](sys):.2} {unit}.\n')
        print(f'C treat liq: {func["C_treat_liq"](sys):.2} {unit}.\n')
        print(f'C bed liq: {func["C_bed_liq"](sys):.2} {unit}.\n')
        print(f'C bed gas: {func["C_bed_gas"](sys):.2} {unit}.\n')
        
        unit = f'kg nitrogen/yr'
        print(f'N urine: {func["N_urine"](sys):.2} {unit}.\n')
        print(f'N feces: {func["N_feces"](sys):.2} {unit}.\n')
        print(f'N toilet sol: {func["N_toilet_sol"](sys):.2} {unit}.\n')
        print(f'N toilet gas: {func["N_toilet_gas"](sys):.2} {unit}.\n')
        print(f'N transport sol: {func["N_trans_sol"](sys):.2} {unit}.\n')
        print(f'N transport sol loss: {func["N_trans_sol_loss"](sys):.2} {unit}.\n')
        print(f'N pretreat sol: {func["N_pretreat_sol"](sys):.2} {unit}.\n')
        print(f'N dryer sol: {func["N_dryer_sol"](sys):.2} {unit}.\n')
        print(f'N dryer gas: {func["N_dryer_gas"](sys):.2} {unit}.\n')
        print(f'N pyrolysis biochar: {func["N_pyrolysis_biochar"](sys):.2} {unit}.\n')
        print(f'N pyrolysis gas: {func["N_pyrolysis_gas"](sys):.2} {unit}.\n')
        print(f'N toilet liq: {func["N_toilet_liq"](sys):.2} {unit}.\n')
        print(f'N trans liq: {func["N_trans_liq"](sys):.2} {unit}.\n')
        print(f'N trans liq loss: {func["N_trans_liq_loss"](sys):.2} {unit}.\n')
        print(f'N pretreat liq: {func["N_pretreat_liq"](sys):.2} {unit}.\n')
        print(f'N treat liq: {func["N_treat_liq"](sys):.2} {unit}.\n')
        print(f'N bed liq: {func["N_bed_liq"](sys):.2} {unit}.\n')
        print(f'N bed gas: {func["N_bed_gas"](sys):.2} {unit}.\n')
        print(f'N struvite: {func["N_struvite"](sys):.2} {unit}.\n')
        print(f'N NH3: {func["N_NH3"](sys):.2} {unit}.\n')
        
        unit = f'kg phosphorus/yr'
        print(f'P urine: {func["P_urine"](sys):.2} {unit}.\n')
        print(f'P feces: {func["P_feces"](sys):.2} {unit}.\n')
        print(f'P toilet sol: {func["P_toilet_sol"](sys):.2} {unit}.\n')
        print(f'P toilet gas: {func["P_toilet_gas"](sys):.2} {unit}.\n')
        print(f'P transport sol: {func["P_trans_sol"](sys):.2} {unit}.\n')
        print(f'P transport sol loss: {func["P_trans_sol_loss"](sys):.2} {unit}.\n')
        print(f'P pretreat sol: {func["P_pretreat_sol"](sys):.2} {unit}.\n')
        print(f'P dryer sol: {func["P_dryer_sol"](sys):.2} {unit}.\n')
        print(f'P dryer gas: {func["P_dryer_gas"](sys):.2} {unit}.\n')
        print(f'P pyrolysis biochar: {func["P_pyrolysis_biochar"](sys):.2} {unit}.\n')
        print(f'P pyrolysis gas: {func["P_pyrolysis_gas"](sys):.2} {unit}.\n')
        print(f'P toilet liq: {func["P_toilet_liq"](sys):.2} {unit}.\n')
        print(f'P trans liq: {func["P_trans_liq"](sys):.2} {unit}.\n')
        print(f'P trans liq loss: {func["P_trans_liq_loss"](sys):.2} {unit}.\n')
        print(f'P pretreat liq: {func["P_pretreat_liq"](sys):.2} {unit}.\n')
        print(f'P treat liq: {func["P_treat_liq"](sys):.2} {unit}.\n')
        print(f'P bed liq: {func["P_bed_liq"](sys):.2} {unit}.\n')
        print(f'P bed gas: {func["P_bed_gas"](sys):.2} {unit}.\n')
        print(f'P struvite: {func["P_struvite"](sys):.2} {unit}.\n')
        
        unit = f'kg potassium/yr'
        print(f'K urine: {func["K_urine"](sys):.2} {unit}.\n')
        print(f'K feces: {func["K_feces"](sys):.2} {unit}.\n')
        print(f'K toilet sol: {func["K_toilet_sol"](sys):.2} {unit}.\n')
        print(f'K toilet gas: {func["K_toilet_gas"](sys):.2} {unit}.\n')
        print(f'K transport sol: {func["K_trans_sol"](sys):.2} {unit}.\n')
        print(f'K transport sol loss: {func["K_trans_sol_loss"](sys):.2} {unit}.\n')
        print(f'K pretreat sol: {func["K_pretreat_sol"](sys):.2} {unit}.\n')
        print(f'K dryer sol: {func["K_dryer_sol"](sys):.2} {unit}.\n')
        print(f'K dryer gas: {func["K_dryer_gas"](sys):.2} {unit}.\n')
        print(f'K pyrolysis biochar: {func["K_pyrolysis_biochar"](sys):.2} {unit}.\n')
        print(f'K pyrolysis gas: {func["K_pyrolysis_gas"](sys):.2} {unit}.\n')
        print(f'K toilet liq: {func["K_toilet_liq"](sys):.2} {unit}.\n')
        print(f'K trans liq: {func["K_trans_liq"](sys):.2} {unit}.\n')
        print(f'K trans liq loss: {func["K_trans_liq_loss"](sys):.2} {unit}.\n')
        print(f'K pretreat liq: {func["K_pretreat_liq"](sys):.2} {unit}.\n')
        print(f'K treat liq: {func["K_treat_liq"](sys):.2} {unit}.\n')
        print(f'K bed liq: {func["K_bed_liq"](sys):.2} {unit}.\n')
        print(f'K bed gas: {func["K_bed_gas"](sys):.2} {unit}.\n')
    
        unit=f'{currency}/cap/d'
        print(f'operator_cost_opex: {func["operator_cost_opex"](tea,ppl):.2} {unit}.\n')

        print(f'front_end_cost_capex: {func["front_end_cost_capex"](front_end,tea,ppl):.2} {unit}.\n')
        print(f'front_end_cost_opex: {func["front_end_cost_opex"](front_end,ppl):.2} {unit}.\n')
        print(f'front_end_cost_electricity: {func["front_end_cost_electricity"](front_end,ppl):.2} {unit}.\n')
        
        print(f'transport_cost_opex: {func["transport_cost_opex"](transport,ppl):.2} {unit}.\n')
        
        print(f'pretreatment_cost_capex: {func["pretreatment_cost_capex"](pretreatment,tea,ppl):.2} {unit}.\n')
        print(f'pretreatment_cost_opex: {func["pretreatment_cost_opex"](pretreatment,ppl):.2} {unit}.\n')
        print(f'pretreatment_cost_electricity: {func["pretreatment_cost_electricity"](pretreatment,ppl):.2} {unit}.\n')
        
        print(f'liq_treatment_cost_capex: {func["liq_treatment_cost_capex"](liq_treatment,tea,ppl):.2} {unit}.\n')
        print(f'liq_treatment_cost_opex: {func["liq_treatment_cost_opex"](liq_treatment,ppl):.2} {unit}.\n')
        print(f'liq_treatment_cost_electricity: {func["liq_treatment_cost_electricity"](liq_treatment,ppl):.2} {unit}.\n')
        
        print(f'controls_housing_cost_capex: {func["controls_housing_cost_capex"](controls_housing,tea,ppl):.2} {unit}.\n')
        print(f'controls_housing_cost_opex: {func["controls_housing_cost_opex"](controls_housing,ppl):.2} {unit}.\n')
        print(f'controls_housing_cost_electricity: {func["controls_housing_cost_electricity"](controls_housing,ppl):.2} {unit}.\n')
        
        print(f'dryer_hhx_cost_capex: {func["dryer_hhx_cost_capex"](dryer_hhx,tea,ppl):.2} {unit}.\n')
        print(f'dryer_hhx_cost_opex: {func["dryer_hhx_cost_opex"](dryer_hhx,ppl):.2} {unit}.\n')
        print(f'dryer_hhx_cost_electricity: {func["dryer_hhx_cost_electricity"](dryer_hhx,ppl):.2} {unit}.\n')
        
        print(f'ohx_cost_capex: {func["ohx_cost_capex"](ohx,tea,ppl):.2} {unit}.\n')
        print(f'ohx_cost_opex: {func["ohx_cost_opex"](ohx,ppl):.2} {unit}.\n')
        print(f'ohx_cost_electricity: {func["ohx_cost_electricity"](ohx,ppl):.2} {unit}.\n')
        
        print(f'carbonizer_base_cost_capex: {func["carbonizer_base_cost_capex"](carbonizer_base,tea,ppl):.2} {unit}.\n')
        print(f'carbonizer_base_cost_opex: {func["carbonizer_base_cost_opex"](carbonizer_base,ppl):.2} {unit}.\n')
        print(f'carbonizer_base_cost_electricity: {func["carbonizer_base_cost_electricity"](carbonizer_base,ppl):.2} {unit}.\n')
        
        print(f'pollution_control_cost_capex: {func["pollution_control_cost_capex"](pollution_control,tea,ppl):.2} {unit}.\n')
        print(f'pollution_control_cost_opex: {func["pollution_control_cost_opex"](pollution_control,ppl):.2} {unit}.\n')
        print(f'pollution_control_cost_electricity: {func["pollution_control_cost_electricity"](pollution_control,ppl):.2} {unit}.\n')
        
        # print(f'opex_transport: {func["opex_transport"](sys):.2} {unit}.\n')
        # print(f'captial_pretreatment: {func["captial_pretreatment"](sys):.2} {unit}.\n')
        # print(f'opex_pretreatment: {func["opex_pretreatment"](sys):.2} {unit}.\n')
        # print(f'electricity_pretreatment: {func["electricity_pretreatment"](sys):.2} {unit}.\n')
        
        unit=f'{GWP.unit}/cap/yr'
        print(f'front_end_ghg_capex: {func["front_end_ghg_capex"](front_end,lca,ppl):.2} {unit}.\n')
        print(f'front_end_ghg_opex: {func["front_end_ghg_opex"](front_end,lca,ppl):.2} {unit}.\n')
        print(f'front_end_ghg_electricity: {func["front_end_ghg_electricity"](front_end,ppl):.2} {unit}.\n')
        print(f'front_end_ghg_direct: {func["front_end_ghg_direct"](front_end,lca,ppl):.2} {unit}.\n')
        
        print(f'transport_ghg_opex: {func["transport_ghg_opex"](transport,lca,ppl):.2} {unit}.\n')
        
        print(f'pretreatment_ghg_capex: {func["pretreatment_ghg_capex"](pretreatment,lca,ppl):.2} {unit}.\n')
        print(f'pretreatment_ghg_opex: {func["pretreatment_ghg_opex"](pretreatment,lca,ppl):.2} {unit}.\n')
        print(f'pretreatment_ghg_electricity: {func["pretreatment_ghg_electricity"](pretreatment,ppl):.2} {unit}.\n')
        print(f'pretreatment_ghg_direct: {func["pretreatment_ghg_direct"](pretreatment,lca,ppl):.2} {unit}.\n')
        
        print(f'liq_treatment_ghg_capex: {func["liq_treatment_ghg_capex"](liq_treatment,lca,ppl):.2} {unit}.\n')
        print(f'liq_treatment_ghg_opex: {func["liq_treatment_ghg_opex"](liq_treatment,lca,ppl):.2} {unit}.\n')
        print(f'liq_treatment_ghg_electricity: {func["liq_treatment_ghg_electricity"](liq_treatment,ppl):.2} {unit}.\n')
        print(f'liq_treatment_ghg_direct: {func["liq_treatment_ghg_direct"](liq_treatment,lca,ppl):.2} {unit}.\n')
        
        print(f'controls_housing_ghg_capex: {func["controls_housing_ghg_capex"](controls_housing,lca,ppl):.2} {unit}.\n')
        print(f'controls_housing_ghg_opex: {func["controls_housing_ghg_opex"](controls_housing,lca,ppl):.2} {unit}.\n')
        print(f'controls_housing_ghg_electricity: {func["controls_housing_ghg_electricity"](controls_housing,ppl):.2} {unit}.\n')
        print(f'controls_housing_ghg_direct: {func["controls_housing_ghg_direct"](controls_housing,lca,ppl):.2} {unit}.\n')
        
        print(f'dryer_hhx_ghg_capex: {func["dryer_hhx_ghg_capex"](dryer_hhx,lca,ppl):.2} {unit}.\n')
        print(f'dryer_hhx_ghg_opex: {func["dryer_hhx_ghg_opex"](dryer_hhx,lca,ppl):.2} {unit}.\n')
        print(f'dryer_hhx_ghg_electricity: {func["dryer_hhx_ghg_electricity"](dryer_hhx,ppl):.2} {unit}.\n')
        print(f'dryer_hhx_ghg_direct: {func["dryer_hhx_ghg_direct"](dryer_hhx,lca,ppl):.2} {unit}.\n')
        
        print(f'ohx_ghg_capex: {func["ohx_ghg_capex"](ohx,lca,ppl):.2} {unit}.\n')
        print(f'ohx_ghg_opex: {func["ohx_ghg_opex"](ohx,lca,ppl):.2} {unit}.\n')
        print(f'ohx_ghg_electricity: {func["ohx_ghg_electricity"](ohx,ppl):.2} {unit}.\n')
        print(f'ohx_ghg_direct: {func["ohx_ghg_direct"](ohx,lca,ppl):.2} {unit}.\n')
        
        print(f'carbonizer_base_ghg_capex: {func["carbonizer_base_ghg_capex"](carbonizer_base,lca,ppl):.2} {unit}.\n')
        print(f'carbonizer_base_ghg_opex: {func["carbonizer_base_ghg_opex"](carbonizer_base,lca,ppl):.2} {unit}.\n')
        print(f'carbonizer_base_ghg_electricity: {func["carbonizer_base_ghg_electricity"](carbonizer_base,ppl):.2} {unit}.\n')
        print(f'carbonizer_base_ghg_direct: {func["carbonizer_base_ghg_direct"](carbonizer_base,lca,ppl):.2} {unit}.\n')
        
        print(f'pollution_control_ghg_capex: {func["pollution_control_ghg_capex"](pollution_control,lca,ppl):.2} {unit}.\n')
        print(f'pollution_control_ghg_opex: {func["pollution_control_ghg_opex"](pollution_control,lca,ppl):.2} {unit}.\n')
        print(f'pollution_control_ghg_electricity: {func["pollution_control_ghg_electricity"](pollution_control,ppl):.2} {unit}.\n')
        print(f'pollution_control_ghg_direct: {func["pollution_control_ghg_direct"](pollution_control,lca,ppl):.2} {unit}.\n')
        
        # for i in ('COD', 'N', 'P', 'K'):
        #     print(f'Total {i} recovery is {func[f"get_tot_{i}_recovery"](sys, i):.1%}, '
        #           f'{func[f"get_liq_{i}_recovery"](sys, i):.1%} in liquid, '
        #           f'{func[f"get_sol_{i}_recovery"](sys, i):.1%} in solid, '
        #           f'{func[f"get_gas_{i}_recovery"](sys, i):.1%} in gas.')

def save_all_reports():
    import os
    #path = os.path.dirname(os.path.realpath(__file__))
    path = ('/Users/lane/Desktop/Coding/biogenic_refinery_20210915')
    path += '/results'
    if not os.path.isdir(path):
        os.path.mkdir(path)
    del os
    for i in (sysA, sysB, sysC, sysD, lcaA, lcaB, lcaC, lcaD):
        if isinstance(i, bst.System):
            i.simulate()
            i.save_report(f'{path}/{i.ID}.xlsx')
        else:
            i.save_report(f'{path}/{i.system.ID}_lca.xlsx')

__all__ = ('sysA', 'sysB', 'sysC', 'sysD', 'teaA', 'teaB', 'teaC', 'teaD', 'lcaA', 'lcaB', 'lcaC', 'lcaD',
            'print_summaries', 'save_all_reports',
            *(i.ID for i in sysA.units),
            *(i.ID for i in sysB.units),
            *(i.ID for i in sysC.units),
            *(i.ID for i in sysD.units),
            )



