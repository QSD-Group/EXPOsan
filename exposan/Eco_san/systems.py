#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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

#from qsdsan.systems import bwaise as bw

from _cmps import cmps


#import all units
# from _primary_ES import PrimaryES
# from _anaerobic_ES_bio import AnaerobicESBio
# from _aerobic_ES_bio import AerobicESBio
# from _anoxic_ES import AnoxicES
# from _ECR import ECR
# from _MBR import MBR
# from _murt_toilet import MURTToilet
# from _bio_ES import BioES
# from _recycling_controls import RecyclingControls
# from _solar import Solar
# from _primaryMBR import PrimaryMBR
# from _MBR_ECR import MBRECR


# =============================================================================
# Unit parameters
# =============================================================================

bst.settings.set_thermo(cmps)
#items = ImpactItem._items

currency = qs.currency = 'USD'
bst.speed_up()

household_size = 4
get_household_size = lambda: household_size
household_per_toilet = 5
get_household_per_toilet = lambda: household_per_toilet
get_toilet_user = lambda: get_household_size()*get_household_per_toilet()

# Number of people served by the Eco-san
ppl = 300

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
    return cost

# =============================================================================
# Prices and GWP CFs
# =============================================================================

# Recycled nutrients are sold at a lower price than commercial fertilizers
price_ratio = 1
price_factor = 0.25
get_price_factor = lambda: price_factor

operator_daily_wage = 29
get_operator_daily_wage = lambda: operator_daily_wage

price_dct = {
    'Electricity': 0.06,
    'Concrete': 194,
    'Steel': 2.665,
    'N': 1.507*get_price_factor(),
    'P': 3.983*get_price_factor(),
    'K': 1.333*get_price_factor(),
    'Polymer': 0.75,
    'Resin': 3.335,
    'FilterBag': 4.81,
    'MgOH2': 0,
    'MgCO3': 0.9,
    'H2SO4': 0.3,
    'struvite': 0,
    'salt': 0,
    'HCl_acid': 0
    }


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
    'struvite': 0,
    'salt': 0.266695553, 
    'HCl_acid': 0.8
    
    }



items = ImpactItem.get_all_items()

if not items.get('Excavation'): # prevent from reloading
    import os
    path = os.path.dirname(os.path.realpath(__file__)) + '/data'
    ImpactIndicator.load_from_file(path+'/impact_indicators.tsv')
    item_path = path+'/_impact_item.xlsx'
    ImpactItem.load_from_file(item_path)
    del os

GWP = qs.ImpactIndicator.get_indicator('GWP')

bst.PowerUtility.price = price_dct['Electricity']
items['Concrete'].price = price_dct['Concrete']
items['Steel'].price = price_dct['Steel']

# =============================================================================
# Universal units and functions
# =============================================================================

CH4_item = StreamImpactItem(ID='CH4_item', GWP=GWP_dct['CH4'])
N2O_item = StreamImpactItem(ID='N2O_item', GWP=GWP_dct['N2O'])
N_item = StreamImpactItem(ID='N_item', GWP=GWP_dct['N'])
P_item = StreamImpactItem(ID='P_item', GWP=GWP_dct['P'])
K_item = StreamImpactItem(ID='K_item', GWP=GWP_dct['K'])
e_item = ImpactItem(ID='e_item', functional_unit='kWh', GWP=GWP_dct['Electricity'])
polymer_item = StreamImpactItem(ID='polymer_item', GWP=GWP_dct['Polymer'])
resin_item = StreamImpactItem(ID='resin_item', GWP=GWP_dct['Resin'])
filter_bag_item = StreamImpactItem(ID='filter_bag_item', GWP=GWP_dct['FilterBag'])
MgOH2_item = StreamImpactItem(ID='MgOH2_item', GWP=GWP_dct['MgOH2'])
MgCO3_item = StreamImpactItem(ID='MgCO3_item', GWP=GWP_dct['MgCO3'])
H2SO4_item = StreamImpactItem(ID='H2SO4_item', GWP=GWP_dct['H2SO4'])
struvite_item = StreamImpactItem(ID = 'struvite_item', GWP=GWP_dct['struvite'])
salt_item = StreamImpactItem(ID = 'salt_item', GWP=GWP_dct['salt'])
HCl_acid_item = StreamImpactItem(ID = 'HCl_acid_item', GWP=GWP_dct['HCl_acid'])
#sludge_item = StreamImpactItem(ID = 'sludge_item', GWP=GWP_dct['sludge'])

def batch_create_streams(prefix):
    stream_dct = {}
    stream_dct['CH4'] = WasteStream(f'{prefix}_CH4', phase='g',
                                    stream_impact_item=CH4_item.copy(set_as_source=True))
    stream_dct['N2O'] = WasteStream(f'{prefix}_N2O', phase='g',
                                   stream_impact_item=N2O_item.copy(set_as_source=True))
    stream_dct['liq_N'] = WasteStream(f'{prefix}_liq_N', phase='l', price=price_dct['N'],
                                     stream_impact_item=N_item.copy(set_as_source=True))
    stream_dct['sol_N'] = WasteStream(f'{prefix}_sol_N', phase='l', price=price_dct['N'],
                                     stream_impact_item=N_item.copy(set_as_source=True))
    stream_dct['liq_P'] = WasteStream(f'{prefix}_liq_P', phase='l', price=price_dct['P'],
                                     stream_impact_item=P_item.copy(set_as_source=True))
    stream_dct['sol_P'] = WasteStream(f'{prefix}_sol_P', phase='l', price=price_dct['P'],
                                     stream_impact_item=P_item.copy(set_as_source=True))
    stream_dct['liq_K'] = WasteStream(f'{prefix}_liq_K', phase='l', price=price_dct['K'],
                                     stream_impact_item=K_item.copy(set_as_source=True))
    stream_dct['sol_K'] = WasteStream(f'{prefix}_sol_K', phase='l', price=price_dct['K'],
                                     stream_impact_item=K_item.copy(set_as_source=True))
    stream_dct['polymer'] = WasteStream(f'{prefix}_polymer', phase='s', price=price_dct['Polymer'],
                                     stream_impact_item=resin_item.copy(set_as_source=True))
    stream_dct['resin'] = WasteStream(f'{prefix}_resin', phase='s', price=price_dct['Resin'],
                                     stream_impact_item=resin_item.copy(set_as_source=True))
    stream_dct['filter_bag'] = WasteStream(f'{prefix}_filter_bag', phase='s', price=price_dct['FilterBag'],
                                     stream_impact_item=filter_bag_item.copy(set_as_source=True))
    stream_dct['MgOH2'] = WasteStream(f'{prefix}_MgOH2', phase='s', price=price_dct['MgOH2'],
                                     stream_impact_item=MgOH2_item.copy(set_as_source=True))
    stream_dct['MgCO3'] = WasteStream(f'{prefix}_MgCO3', phase='s', price=price_dct['MgCO3'],
                                     stream_impact_item=MgCO3_item.copy(set_as_source=True))
    stream_dct['H2SO4'] = WasteStream(f'{prefix}_H2SO4', phase='s', price=price_dct['H2SO4'],
                                     stream_impact_item=H2SO4_item.copy(set_as_source=True))
    stream_dct['struvite'] = WasteStream(f'{prefix}_struvite', phase='s', price=price_dct['struvite'],
                                     stream_impact_item=struvite_item.copy(set_as_source=True))
    stream_dct['salt'] = WasteStream(f'{prefix}_salt', phase='s', price=price_dct['salt'],
                                     stream_impact_item=salt_item.copy(set_as_source=True))
    stream_dct['HCl_acid'] = WasteStream(f'{prefix}_HCl_acid', phase='l', price=price_dct['HCl_acid'],
                                     stream_impact_item=HCl_acid_item.copy(set_as_source=True))
   # stream_dct['sludge'] = WasteStream(f'{prefix}_sludge', phase='s', price=price_dct['sludge'],
                                     #stream_impact_item=sludge_item.copy(set_as_source=True))
    
    

    return stream_dct
    

def add_fugitive_items(unit, item):
    unit._run()
    for i in unit.ins:
        i.stream_impact_item = item.copy(set_as_source=True)
        


# %%

# # =============================================================================
# # Scenario A (sysA): Original Design tested in Durban A+O+A+O+ECR
# # =============================================================================

# flowsheetA = bst.Flowsheet('sysA')
# bst.main_flowsheet.set_flowsheet(flowsheetA)

# streamsA = batch_create_streams('A')

# #################### Human Inputs ####################
# A1 = su.Excretion('A1', outs=('urine','feces'))

# ################### User Interface ###################
# A2 = su.MURTToilet('A2', ins=(A1-0, A1-1,
#                               'toilet_paper', 'flushing_water',
#                               'cleansing_water', 'desiccant'),
#                     outs=('mixed_waste', 'A2_CH4', 'A2_N2O'),
#                     decay_k_COD=get_decay_k(tau_deg, log_deg), 
#                     decay_k_N=get_decay_k(tau_deg, log_deg),
#                     max_CH4_emission=get_max_CH4_emission(),
#                     N_user=300/7, N_toilet=7,
#                     if_flushing=True, if_desiccant=False, if_toilet_paper=True,
#                     CAPEX = 0,
#                     OPEX_over_CAPEX= 0.07)

# ###################### Treatment ######################
# A3 = su.PrimaryMBR('A3', ins=(A2-0), 
#                     outs=('treated', 'A3_CH4', 'A3_N2O', 'sludge'), 
#                     decay_k_COD=get_decay_k(tau_deg, log_deg), 
#                     decay_k_N=get_decay_k(tau_deg, log_deg),
#                     max_CH4_emission=get_max_CH4_emission())

# A4 = su.AnaerobicESBio('A4', ins=(A3-0), outs=('treated', 'A4_CH4', 'A4_N2O'), 
#                     decay_k_COD=get_decay_k(tau_deg, log_deg),
#                     decay_k_N=get_decay_k(tau_deg, log_deg),
#                     max_CH4_emission=get_max_CH4_emission())
# A5 = su.AerobicESBio('A5', ins=(A4-0), outs = ('treated','A5_CH4', 'A5_N2O'), 
#                     decay_k_COD=get_decay_k(tau_deg, log_deg),
#                     decay_k_N=get_decay_k(tau_deg, log_deg),
#                     max_CH4_emission=get_max_CH4_emission())
# A6 = su.AnoxicES('A6', ins=(A5-0), outs = ('treated', 'A6_CH4', 'A6_N2O'),
#               decay_k_COD=get_decay_k(tau_deg, log_deg),
#                     decay_k_N=get_decay_k(tau_deg, log_deg),
#                     max_CH4_emission=get_max_CH4_emission()) #double check
# A7 = su.AerobicESBio('A7', ins=(A6-0), outs = ('treated','A7_CH4', 'A7_N2O'), 
#                     decay_k_COD=get_decay_k(tau_deg, log_deg),
#                     decay_k_N=get_decay_k(tau_deg, log_deg),
#                     max_CH4_emission=get_max_CH4_emission())
# A8 = su.ECR('A8', ins=(A7-0, streamsA['salt'], streamsA['HCl_acid']), 
#                     outs = ('treated'),
#                     decay_k_COD=get_decay_k(tau_deg, log_deg),)


# A10 = su.Mixer('A10', ins=(A3-1, A2-1, A4-1, A5-1, A6-1, A7-1), outs=streamsA['CH4'])
# A10.specification = lambda: add_fugitive_items(A10, CH4_item)
# A10.line = 'fugitive CH4 mixer' 
        
# A11 = su.Mixer('A11', ins=(A3-2, A2-2, A4-2, A5-2, A6-2, A7-2), outs=streamsA['N2O'])
# A11.specification = lambda: add_fugitive_items(A11, N2O_item)
# A11.line = 'fugitive N2O mixer'



# ################## Other impacts and costs ##################
# A13 = su.BioES('A13', ins=(A8-0), outs = ('A13_out'))
# A14 = su.RecyclingControls('A14', ins=(A13-0), outs = ('A14_out'))
# A15 = su.Solar('A15', ins=(A14-0), outs = ('A15_out'))

# A16 = su.Trucking('A16', ins=A3-3, outs=('transported', 'conveyance_loss'),
#                   load_type='mass', distance=5, distance_unit='km',
#                   interval=365, interval_unit='d',
#                   loss_ratio=0.02)
# def update_A16_param():
#     A16._run()
#     truck = A16.single_truck
#     truck.interval = 365*24
#     truck.load = A16.F_mass_in*truck.interval
#     rho = A16.F_mass_in/A16.F_vol_in
#     vol = truck.load/rho
#     A16.fee = get_tanker_truck_fee(vol)
#     A16._design()
# A16.specification = update_A16_param


# ############### Simulation, TEA, and LCA ###############
# sysA = bst.System('sysA', path= (A1, A2, A3, A4, A5, A6, A7, A8, A10, A11, A13, A14, A15, A16))

# sysA.simulate()

# power = sum([u.power_utility.rate for u in sysA.units])

# #!!! update labor to input country specific data and be a distribution
# teaA = SimpleTEA(system=sysA, discount_rate=get_discount_rate(), 
#                   start_year=2020, lifetime=10, uptime_ratio=1, 
#                   lang_factor=None, annual_maintenance=0, 
#                   annual_labor=(get_operator_daily_wage() * 1 * 52), construction_schedule=None)

# lcaA = LCA(system=sysA, lifetime=10, lifetime_unit='yr', uptime_ratio=1,
#             e_item=lambda: power*(365*24)*10)



# # %%

# # =============================================================================
# # Scenario B (sysB): MBR (Anoxic and Facultative Aerobic with MBR)
# # =============================================================================

# ###power sum([(u.power_utility.rate) for u in sysA.units])

# flowsheetB = bst.Flowsheet('sysB')
# bst.main_flowsheet.set_flowsheet(flowsheetB)

# streamsB = batch_create_streams('B')

# # # #################### Human Inputs #################### 

# B1 = su.Excretion('B1', outs=('urine','feces'))

# # # ################### User Interface ###################
# B2 = su.MURTToilet('B2', ins=(B1-0, B1-1,
#                                 'toilet_paper', 'flushing_water',
#                                 'cleansing_water', 'desiccant'),
#                     outs=('mixed_waste', 'B2_CH4', 'B2_N2O'),
#                     decay_k_COD=get_decay_k(tau_deg, log_deg), 
#                     decay_k_N=get_decay_k(tau_deg, log_deg),
#                     max_CH4_emission=get_max_CH4_emission(),
#                     N_user=300/7, N_toilet=7,
#                     if_flushing=True, if_desiccant=False, if_toilet_paper=True,
#                     CAPEX = 0,
#                     OPEX_over_CAPEX= 0.07)

# ###################### Treatment ######################

# B3 = su.PrimaryES('B3', ins=(B2-0), 
#                     outs=('treated', 'B3_CH4', 'B3_N2O', 'sludge'), 
#                     decay_k_COD=get_decay_k(tau_deg, log_deg), 
#                     decay_k_N=get_decay_k(tau_deg, log_deg),
#                     max_CH4_emission=get_max_CH4_emission())

# B4 = su.AnaerobicESBio('B4', ins=(B3-0), outs=('treated', 'B4_CH4', 'B4_N2O'),
#                     decay_k_COD=get_decay_k(tau_deg, log_deg), 
#                     decay_k_N=get_decay_k(tau_deg, log_deg),
#                     max_CH4_emission=get_max_CH4_emission())
# B5 = su.AerobicESBio('B5', ins=(B4-0), outs = ('treated','B5_CH4', 'B5_N2O'), 
#                     decay_k_COD=get_decay_k(tau_deg, log_deg),
#                     decay_k_N=get_decay_k(tau_deg, log_deg),
#                     max_CH4_emission=get_max_CH4_emission())

# B6 = su.MBR('B6', ins=(B5-0), outs = ('treated', 'B6_CH4', 'B6_N2O'),
#                     decay_k_COD=get_decay_k(tau_deg, log_deg), 
#                     decay_k_N=get_decay_k(tau_deg, log_deg),
#                     max_CH4_emission=get_max_CH4_emission())

# B8 = su.MBRECR('B8', ins=(B6-0, 'salt', 'HCl_acid'), 
#                     outs = ('treated'),
#                      decay_k_COD=get_decay_k(tau_deg, log_deg))


# B10 = su.Mixer('B10', ins=(B2-1, B3-1, B4-1, B5-1, B6-1), outs=streamsB['CH4'])
# B10.specification = lambda: add_fugitive_items(B10, CH4_item)
# B10.line = 'fugitive CH4 mixer' 
        
# B11 = su.Mixer('B11', ins=(B2-2, B3-2, B4-2, B5-2, B6-2), outs=streamsB['N2O'])
# B11.specification = lambda: add_fugitive_items(B11, N2O_item)
# B11.line = 'fugitive N2O mixer'

# ################## Other impacts and costs ##################
# #B13 = su.BioES('B13', ins=(B8-0), outs = ('B13_out'))
# B14 = su.RecyclingControls('B14', ins=(B8-0), outs = ('B14_out'))
# B15 = su.Solar('B15', ins=(B14-0), outs = ('B15_out'))


# B16 = su.Trucking('B16', ins=B3-3, outs=('transported', 'conveyance_loss'),
#                   load_type='mass', distance=5, distance_unit='km',
#                   interval=365, interval_unit='d',
#                   loss_ratio=0.02)
# def update_B16_param():
#     B16._run()
#     truck = B16.single_truck
#     truck.interval = 365*24
#     truck.load = B16.F_mass_in*truck.interval/7
#     rho = B16.F_mass_in/B16.F_vol_in
#     vol = truck.load/rho
#     B16.fee = get_tanker_truck_fee(vol)
#     B16._design()
# B16.specification = update_B16_param

# # ############### Simulation, TEA, and LCA ###############
# sysB = bst.System('sysB', path= (B1, B2, B3, B4, B5, B6, B8, B10, B11, B14, B15, B16))

# sysB.simulate()

# power = sum([u.power_utility.rate for u in sysB.units])


# teaB = SimpleTEA(system=sysB, discount_rate=get_discount_rate(), 
#                   start_year=2020, lifetime=10, uptime_ratio=1, 
#                   lang_factor=None, annual_maintenance=0, 
#                   annual_labor=(get_operator_daily_wage() * 1*12), construction_schedule=None)

# lcaB= LCA(system=sysB, lifetime=10, lifetime_unit='yr', uptime_ratio=1,
#             e_item=lambda: power*(365*24)*10)

# # # %%

# # # =============================================================================
# # # Scenario C (sysC): MBR + struvite precipitation (Anoxic and Facultative Aerobic with MBR)
# # # =============================================================================

flowsheetC = bst.Flowsheet('sysC')
bst.main_flowsheet.set_flowsheet(flowsheetC)

streamsC = batch_create_streams('C')


C1 = su.Excretion('C1', outs=('urine','feces'))

# ################### User Interface ###################
C2 = su.MURTToilet('C2', ins=(C1-0, C1-1,
                                'toilet_paper', 'flushing_water',
                                'cleansing_water', 'desiccant'),
                    outs=('mixed_waste', 'C2_CH4', 'C2_N2O'),
                    decay_k_COD=get_decay_k(tau_deg, log_deg), 
                    decay_k_N=get_decay_k(tau_deg, log_deg),
                    max_CH4_emission=get_max_CH4_emission(),
                    N_user=get_toilet_user(), N_toilet=ppl/get_toilet_user(),
                    if_flushing=True, if_desiccant=False, if_toilet_paper=True,
                    CAPEX = 0,
                    OPEX_over_CAPEX= 0.07)


###################### Treatment ######################
C3 = su.PrimaryMBR('C3', ins=(C2-0, streamsC['MgOH2']), 
                    outs=('treated', 'C3_CH4', 'C3_N2O', 'sludge', streamsC['struvite']), 
                    decay_k_COD=get_decay_k(tau_deg, log_deg), 
                    decay_k_N=get_decay_k(tau_deg, log_deg),
                    max_CH4_emission=get_max_CH4_emission())
C4 = su.AnaerobicESBio('C4', ins=(C3-0), outs=('treated', 'C4_CH4', 'C4_N2O'),
                    decay_k_COD=get_decay_k(tau_deg, log_deg), 
                    decay_k_N=get_decay_k(tau_deg, log_deg),
                    max_CH4_emission=get_max_CH4_emission())
C5 = su.AerobicESBio('C5', ins=(C4-0), outs = ('treated','C5_CH4', 'C5_N2O'), 
                    decay_k_COD=get_decay_k(tau_deg, log_deg),
                    decay_k_N=get_decay_k(tau_deg, log_deg),
                    max_CH4_emission=get_max_CH4_emission())
C6 = su.MBR('C6', ins=(C5-0), outs = ('treated', 'C6_CH4', 'C6_N2O'),
                    decay_k_COD=get_decay_k(tau_deg, log_deg), 
                    decay_k_N=get_decay_k(tau_deg, log_deg),
                    max_CH4_emission=get_max_CH4_emission())
C8 = su.MBRECR('C8', ins=(C6-0, 'salt', 'HCl_acid'), 
                     outs = ('treated'),
                     decay_k_COD=get_decay_k(tau_deg, log_deg))


C10 = su.Mixer('C10', ins=(C2-1, C3-1, C4-1, C5-1, C6-1), outs=streamsC['CH4'])
C10.specification = lambda: add_fugitive_items(C10, CH4_item)
C10.line = 'fugitive CH4 mixer' 
        
C11 = su.Mixer('C11', ins=(C2-2, C3-2, C4-2, C5-2, C6-2), outs=streamsC['N2O'])
C11.specification = lambda: add_fugitive_items(C1, N2O_item)
C11.line = 'fugitive N2O mixer'

################## Other impacts and costs ##################
C14 = su.RecyclingControls('C14', ins=(C8-0), outs = ('C14_out'))

# C15 = su.Solar('C15', ins=(C14-0), outs = ('C15_out'))


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

# ############### Simulation, TEA, and LCA ###############
sysC = bst.System('sysC', path= (C1, C2, C3, C4, C5, C6, C8, C10, C11, C14, C16))

sysC.simulate()

powerC = sum([u.power_utility.rate for u in sysC.units])

# power = 0
teaC = SimpleTEA(system=sysC, discount_rate=get_discount_rate(), 
                  start_year=2020, lifetime=10, uptime_ratio=1, 
                  lang_factor=None, annual_maintenance=0, 
                  annual_labor=(get_operator_daily_wage() * 1*12))

lcaC= LCA(system=sysC, lifetime=10, lifetime_unit='yr', uptime_ratio=1,
            e_item=lambda: powerC*(365*24)*10)



# %%

# =============================================================================
# Summarizing Functions
# =============================================================================

#!!! items need to be updated in sys_dct, system_streams, learning curve assumptions
sys_dct = {
    'ppl':  dict(sysC=300),
    'input_unit': dict(sysC=C1),
    'liq_unit': dict(sysA=None, sysB=None, sysC=None),
    'sol_unit': dict(sysA=None, sysB=None, sysC=None),
    'gas_unit': dict(sysA=None, sysB=None, sysC=None),
    'stream_dct': dict(sysC=streamsC),
    'TEA': dict(sysC=teaC),
    'LCA': dict(sysC=lcaC),
    'cache': dict(sysC={}),
    }

system_streams = {sysC: streamsC}

#learning curve assumptions
percent_CAPEX_to_scale = 0.1 #this number is the decimal of the fraction of scost of pecialty parts/cost of total parts
get_percent_CAPEX_to_scale = lambda: percent_CAPEX_to_scale

number_of_units = 100000

percent_limit = 0.015
get_percent_limit = lambda: percent_limit

learning_curve_percent=0.925
get_learning_curve_percent = lambda: learning_curve_percent

def get_scaled_capital(tea):
    CAPEX_to_scale = tea.annualized_CAPEX * get_percent_CAPEX_to_scale()
    CAPEX_not_scaled = tea.annualized_CAPEX - CAPEX_to_scale
    scaled_limited = CAPEX_to_scale * get_percent_limit()
    b=(np.log(get_learning_curve_percent())/np.log(2))
    scaled_CAPEX_annualized  = (CAPEX_to_scale - scaled_limited)*number_of_units**b+scaled_limited    
    new_CAPEX_annualized = scaled_CAPEX_annualized + CAPEX_not_scaled    
    return new_CAPEX_annualized
    

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
        if not i.impact_item: continue
        emission[f'{i.ID}'] = i.F_mass*i.impact_item.CFs['GlobalWarming']*factor
    return emission

#10 corresponds to year
def get_fugitive_emissions(streams=None, hr=365*24*10):
    try: iter(streams)
    except: streams = (streams,)
    factor = hr
    emission=0
    for i in streams:
        if not i.impact_item: continue
        emission = i.F_mass*i.impact_item.CFs['GlobalWarming']*factor
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

def get_summarizing_fuctions():
    func_dct = {}
    func_dct['get_annual_cost'] = lambda tea, ppl: (tea.EAC-tea.annualized_CAPEX
                                                    +get_scaled_capital(tea))/ppl
    func_dct['get_annual_CAPEX'] = lambda tea, ppl: get_scaled_capital(tea)/ppl
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
            /lca.lifetime/ppl - (get_fugitive_gases(sys)/lca.lifetime/ppl) #maybe ignore
    func_dct['get_offset_GWP'] = \
        lambda lca, ppl: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='offset')[ind] \
            /lca.lifetime/ppl
    func_dct['get_other_GWP'] = \
        lambda lca, ppl: lca.total_other_impacts[ind]/lca.lifetime/ppl
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
    func = get_summarizing_fuctions()
    for sys in systems:
        sys.simulate()
        ppl = sys_dct['ppl'][sys.ID]
        print(f'\n---------- Summary for {sys} ----------\n')
        tea = sys_dct['TEA'][sys.ID]
        tea.show()
        print('\n')
        lca = sys_dct['LCA'][sys.ID]
        lca.show()
        
        unit = f'{currency}/cap/yr'
        print(f'\nTotal cost: {func["get_annual_cost"](tea, ppl):.2f} {unit}.')
        print(f'Capital: {func["get_annual_CAPEX"](tea, ppl):.2f} {unit}.')
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

        # for i in ('COD', 'N', 'P', 'K'):
        #     print(f'Total {i} recovery is {func[f"get_tot_{i}_recovery"](sys, i):.1%}, '
        #           f'{func[f"get_liq_{i}_recovery"](sys, i):.1%} in liquid, '
        #           f'{func[f"get_sol_{i}_recovery"](sys, i):.1%} in solid, '
        #           f'{func[f"get_gas_{i}_recovery"](sys, i):.1%} in gas.')

def save_all_reports():
 # import os
    # path = os.path.dirname(os.path.realpath(__file__))
    # path += '/results'
    # if not os.path.isdir(path):
    #     os.path.mkdir(path)
    # del os
    path = '/Users/torimorgan/opt/anaconda3/lib/python3.8/site-packages/exposan/eco_san'
    for i in (sysC, lcaC):
        if isinstance(i, bst.System):
            i.simulate()
            i.save_report(f'{path}/{i.ID}.xlsx')
        else:
            i.save_report(f'{path}/{i.system.ID}_lca.xlsx')


__all__ = ('sysC', 'teaC', 'lcaC',
            'print_summaries', 'save_all_reports',
            *(i.ID for i in sysC.units),
            )
