#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np
from sklearn.linear_model import LinearRegression as LR
from qsdsan import (
    Flowsheet, main_flowsheet,
    set_thermo, WasteStream, System, PowerUtility,
    sanunits as su,
    ImpactIndicator, ImpactItem, StreamImpactItem,
    SimpleTEA, LCA,
    )
from qsdsan.utils import clear_lca_registries
from exposan.reclaimer import data_path, results_path
from exposan.reclaimer._cmps import cmps

# Clear registry to avoid replacement warning messages
clear_lca_registries()

set_thermo(cmps)
currency = 'USD'


# =============================================================================
# Unit parameters
# =============================================================================

# Number of people served by the Reclaimer
#   baseline: 30 users
#   scaled: 120 users
ppl = 120  # total population served by system
N_user = 25  # number of users per toilet

discount_rate = 0.05

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3
# Get reduction rate constant k for COD and N, use a function so that k can be
def get_decay_k(tau_deg=2, log_deg=3):
    k = (-1/tau_deg)*np.log(10**-log_deg)
    return k

max_CH4_emission = 0.25

# =============================================================================
# Prices and GWP CFs
# =============================================================================

# Recycled nutrients are sold at a lower price than commercial fertilizers
price_factor = 0.25
get_price_factor = lambda: price_factor

# Should be changed based on country
price_ratio = 1

price_dct = {
    'Electricity': 0.06,
    'wages': 3.64,
    'Concrete': 194*price_ratio,
    'Steel': 2.665*price_ratio,
    'N': 1.507*get_price_factor(),
    'P': 3.983*get_price_factor(),
    'K': 1.333*get_price_factor(),
    'Polymer': 0.75*price_ratio,
    'Resin': 3.335*price_ratio,
    'FilterBag': 4.81*price_ratio,
    'MgOH2': 0.0*price_ratio,
    'MgCO3': 0.9*price_ratio,
    'H2SO4': 0.3*price_ratio,
    'struvite': 0.0*get_price_factor(),
    'salt': 0.0*price_ratio,
    'HCl': 0.0*price_ratio,
    'KCl': 0.0*price_ratio,
    'GAC': 0.0*price_ratio,
    'Zeolite': 0.0*price_ratio,
    'Conc_NH3': 0.0*get_price_factor(),
    'sludge': 0.0
    }


GWP_dct = {
    'Electricity': 0.69,
    'CH4': 28,
    'N2O': 265,
    'N': -5.4,
    'P': -4.9,
    'K': -1.5,
    'Polymer': 2.8,
    'Resin': 1.612,
    'FilterBag': 0.464,  # based on 0.05 kg of nylon and nylon's GWP of 9.279255342 kgCO2eq per kg
    'MgOH2': 1.176277921,
    'MgCO3': 1.176277921,
    'H2SO4': 0.158899487,
    'struvite': 0,  # -4.9*(31/245),
    'salt': 0.266695553, 
    'HCl': 0.8,
    'KCl': 0.8,
    'GAC': 8.388648277,
    'Zeolite': 5.175,
    'Conc_NH3':0,  # -5.4*(14/17),
    'sludge':0
    }


ImpactIndicator.load_from_file(os.path.join(data_path, 'impact_indicators.csv'))
ImpactItem.load_from_file(os.path.join(data_path, 'impact_items.xlsx'))

GWP = ImpactIndicator.get_indicator('GWP')

PowerUtility.price = price_dct['Electricity']
ImpactItem.get_item('Concrete').price = price_dct['Concrete']
ImpactItem.get_item('Steel').price = price_dct['Steel']

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
struvite_item = StreamImpactItem(ID='struvite_item', GWP=GWP_dct['struvite'])
salt_item = StreamImpactItem(ID='salt_item', GWP=GWP_dct['salt'])
KCl_item = StreamImpactItem(ID='KCl_item', GWP=GWP_dct['KCl'])
GAC_item = StreamImpactItem(ID='GAC_item', GWP=GWP_dct['GAC'])
Zeolite_item = StreamImpactItem(ID='Zeolite_item', GWP=GWP_dct['Zeolite'])
Conc_NH3_item = StreamImpactItem(ID='Conc_NH3', GWP=GWP_dct['Conc_NH3'])
HCl_item = StreamImpactItem(ID='HCl', GWP=GWP_dct['HCl'])
sludge_item = StreamImpactItem(ID='sludge_item', GWP=GWP_dct['sludge'])


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
    stream_dct['KCl'] = WasteStream(f'{prefix}_KCl', phase='l', price=price_dct['KCl'],
                                    stream_impact_item=KCl_item.copy(set_as_source=True))
    stream_dct['GAC'] = WasteStream(f'{prefix}_GAC', phase='s', price=price_dct['GAC'],
                                    stream_impact_item=GAC_item.copy(set_as_source=True))
    stream_dct['Zeolite'] = WasteStream(f'{prefix}_Zeolite', phase='s', price=price_dct['Zeolite'],
                                        stream_impact_item=Zeolite_item.copy(set_as_source=True))
    stream_dct['Conc_NH3'] = WasteStream(f'{prefix}_Conc_NH3', phase='s', price=price_dct['Conc_NH3'],
                                         stream_impact_item=Conc_NH3_item.copy(set_as_source=True))
    stream_dct['HCl'] = WasteStream(f'{prefix}_HCl', phase='l', price=price_dct['HCl'],
                                    stream_impact_item=Conc_NH3_item.copy(set_as_source=True))
    return stream_dct
    

def add_fugitive_items(unit, item_ID):
    unit._run()
    for i in unit.ins:
        i.stream_impact_item = ImpactItem.get_item(item_ID)
        # i.stream_impact_item = ImpactItem.get_item(item_ID).copy(set_as_source=True)
        

# %%

# =============================================================================
# Scenario A (sysA): Duke Reclaimer Design 2.0 Solids Treatment Only 
# =============================================================================

# Set flowsheet to avoid stream replacement warnings
flowsheetA = Flowsheet('sysA')
main_flowsheet.set_flowsheet(flowsheetA)
streamsA = batch_create_streams('A')

#################### Human Inputs ####################
A1 = su.Excretion('A1', outs=('urine', 'feces'))

###################### Treatment ######################

A3 = su.PrimaryReclaimer('A3', ins=(A1-0, streamsA['MgOH2']),
                         outs=('A3_treated', 'A3_CH4', 'A3_N2O', 'A3_sludge', streamsA['struvite']),
                         decay_k_COD=get_decay_k(tau_deg, log_deg),
                         decay_k_N=get_decay_k(tau_deg, log_deg),
                         max_CH4_emission=max_CH4_emission,
                         ppl=ppl,
                         if_include_front_end=True
                         )

A4 = su.SludgePasteurizationReclaimer('A4', ins=(A3-3, 'air', 'lpg'), outs='treated_sludge',
                                      heat_loss=0.1, target_MC=0.1, sludge_temp=283.15,
                                      temp_pasteurization=343.15, lhv_lpg=48.5,
                                      ppl=ppl, if_sludge_service=True
                                      )

A8 = su.Mixer('A8', ins=(A3-1), outs=streamsA['CH4'])
A8.specification = lambda: add_fugitive_items(A8, 'CH4_item')
A8.line = 'fugitive CH4 mixer' 
        
A9 = su.Mixer('A9', ins=(A3-2), outs=streamsA['N2O'])
A9.specification = lambda: add_fugitive_items(A9, 'N2O_item')
A9.line = 'fugitive N2O mixer'

############### Simulation, TEA, and LCA ###############
sysA = System('sysA', path=(A1, A3, A4, A8, A9))

teaA = SimpleTEA(system=sysA, discount_rate=discount_rate,
                 start_year=2020, lifetime=20, uptime_ratio=1,
                 lang_factor=None, annual_maintenance=0,
                 annual_labor=(A4._calc_maintenance_labor_cost() * 24 * 365)
                 )

get_powerA = lambda: sum([u.power_utility.rate for u in sysA.units]) * (24 * 365 * teaA.lifetime)

lcaA = LCA(system=sysA, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerA)

# =============================================================================
# System B Duke Reclaimer 2.0 Full System on grid-tied energy
# =============================================================================

# Set flowsheet to avoid stream replacement warnings
flowsheetB = Flowsheet('sysB')
main_flowsheet.set_flowsheet(flowsheetB)
streamsB = batch_create_streams('B')

#################### Human Inputs ####################
B1 = su.Excretion('B1', outs=('urine','feces'))

################### User Interface ###################
# Reclaimer 2.0 can process ~30L/hr(net), 720L/24 hours of constant operation
# flush volume of 6L per flush determines 120 flushes in one day (~20 users if they use toilet 6 times)
B2 = su.MURTToilet('B2', ins=(B1-0, B1-1,
                              'toilet_paper', 'flushing_water',
                              'cleansing_water', 'desiccant'),
                   outs=('mixed_waste', 'B2_CH4', 'B2_N2O'),
                   decay_k_COD=get_decay_k(tau_deg, log_deg),
                   decay_k_N=get_decay_k(tau_deg, log_deg),
                   max_CH4_emission=max_CH4_emission,
                   N_user=N_user, N_tot_user=ppl, lifetime=8,
                   if_flushing=True, if_desiccant=False, if_toilet_paper=True,
                   CAPEX=0, OPEX_over_CAPEX=0.07)

##################### Treatment ######################
B3 = su.PrimaryReclaimer('B3', ins=(B2-0, streamsB['MgOH2']),
                         outs=('B3_treated', 'B3_CH4', 'B3_N2O', 'B3_sludge', streamsB['struvite']),
                         decay_k_COD=get_decay_k(tau_deg, log_deg),
                         decay_k_N=get_decay_k(tau_deg, log_deg),
                         max_CH4_emission=max_CH4_emission,
                         ppl=ppl,
                         if_include_front_end=True
                         )

B4 = su.SludgePasteurizationReclaimer('B4', ins=(B3-3, 'air', 'lpg'), outs='treated_sludge',
                                      heat_loss=0.1, target_MC=0.1, sludge_temp=283.15,
                                      temp_pasteurization=343.15, lhv_lpg=48.5,
                                      ppl=ppl, if_sludge_service=True
                                      )

B5 = su.UltrafiltrationReclaimer('B5', ins=(B3-0),
                                 outs=('B5_treated', 'retentate'),
                                 ppl=ppl,
                                 if_gridtied=True
                                 )
                        
B6 = su.IonExchangeReclaimer('B6', ins=(B5-0, streamsB['Zeolite'], streamsB['GAC'], streamsB['KCl']),
                             outs=('B6_treated', 'SpentZeolite', 'SpentGAC',streamsB['Conc_NH3']),
                             decay_k_COD=get_decay_k(tau_deg, log_deg),
                             decay_k_N=get_decay_k(tau_deg, log_deg),
                             max_CH4_emission=max_CH4_emission,
                             ppl=ppl
                             )

B7 = su.ECR_Reclaimer('B7', ins=(B6-0),
                      outs='B7_treated',
                      decay_k_COD=get_decay_k(tau_deg, log_deg),
                      ppl=ppl,
                      if_gridtied=True
                      )

B8 = su.Mixer('B8', ins=(B3-1), outs=streamsB['CH4'])
B8.specification = lambda: add_fugitive_items(B7, CH4_item)
B8.line = 'fugitive CH4 mixer' 
        
B9 = su.Mixer('B9', ins=(B3-2), outs=streamsB['N2O'])
B9.specification = lambda: add_fugitive_items(B8, N2O_item)
B9.line = 'fugitive N2O mixer'

################## Non-Treatment ##################
B10 = su.HousingReclaimer('B10', ins=(B7-0), outs='B10_out', ppl=ppl)
B11 = su.SystemReclaimer('B11', ins=(B10-0), outs='B11_out', ppl=ppl, if_gridtied=True)

############### Simulation, TEA, and LCA ###############
sysB = System('sysB', path=(B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11))

teaB = SimpleTEA(system=sysB, discount_rate=discount_rate,
                 start_year=2020, lifetime=20, uptime_ratio=1,
                 lang_factor=None, annual_maintenance=0,
                 annual_labor=((B4._calc_maintenance_labor_cost() + B6._calc_maintenance_labor_cost()) * 24 * 365)
                 )

get_powerB = lambda: sum([u.power_utility.rate for u in sysB.units]) * (24 * 365 * teaB.lifetime)

lcaB = LCA(system=sysB, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerB)

# =============================================================================
# System C Duke Reclaimer 2.0 coupled with photovoltaic system
# =============================================================================

# Set flowsheet to avoid stream replacement warnings
flowsheetC = Flowsheet('sysC')
main_flowsheet.set_flowsheet(flowsheetC)
streamsC = batch_create_streams('C')

#################### Human Inputs ####################
C1 = su.Excretion('C1', outs=('urine','feces'))

################### User Interface ###################
# Reclaimer 2.0 can process ~30L/hr(net), 720L/24 hours of constant operation
# flush volume of 6L per flush determines 120 flushes in one day (~20 users if they use toilet 6 times)
C2 = su.MURTToilet('C2', ins=(C1-0, C1-1,
                              'toilet_paper', 'flushing_water',
                              'cleansing_water', 'desiccant'),
                   outs=('mixed_waste', 'C2_CH4', 'C2_N2O'),
                   decay_k_COD=get_decay_k(tau_deg, log_deg),
                   decay_k_N=get_decay_k(tau_deg, log_deg),
                   max_CH4_emission=max_CH4_emission,
                   N_user=N_user, N_tot_user=ppl, lifetime=8,
                   if_flushing=True, if_desiccant=False, if_toilet_paper=True,
                   CAPEX=0, OPEX_over_CAPEX=0.07)

###################### Treatment ######################
C3 = su.PrimaryReclaimer('C3', ins=(C2-0, streamsC['MgOH2']),
                         outs=('C3_treated', 'C3_CH4', 'C3_N2O', 'C3_sludge', streamsC['struvite']),
                         decay_k_COD=get_decay_k(tau_deg, log_deg),
                         decay_k_N=get_decay_k(tau_deg, log_deg),
                         max_CH4_emission=max_CH4_emission,
                         ppl=ppl,
                         if_include_front_end=True
                         )

C4 = su.SludgePasteurizationReclaimer('C4', ins=(C3-3, 'air', 'lpg'), outs=('treated_sludge'),
                                      heat_loss=0.1, target_MC=0.1, sludge_temp=283.15,
                                      temp_pasteurization=343.15, lhv_lpg=48.5,
                                      ppl=ppl, if_sludge_service=True
                                      )

C5 = su.UltrafiltrationReclaimer('C5', ins=(C3-0),
                                 outs=('C5_treated', 'retentate'),
                                 ppl=ppl,
                                 if_gridtied=False
                                 )
                        
C6 = su.IonExchangeReclaimer('C6', ins=(C5-0, streamsC['Zeolite'], streamsC['GAC'], streamsC['KCl']),
                             outs=('C6_treated', 'SpentZeolite', 'SpentGAC',streamsC['Conc_NH3']),
                             decay_k_COD=get_decay_k(tau_deg, log_deg),
                             decay_k_N=get_decay_k(tau_deg, log_deg),
                             max_CH4_emission=max_CH4_emission,
                             ppl=ppl
                             )

C7 = su.ECR_Reclaimer('C7', ins=(C6-0),
                      outs='C7_treated',
                      decay_k_COD=get_decay_k(tau_deg, log_deg),
                      ppl=ppl,
                      if_gridtied=False
                      )

C8 = su.Mixer('C8', ins=(C3-1, C2-1), outs=streamsC['CH4'])
C8.specification = lambda: add_fugitive_items(C7, CH4_item)
C8.line = 'fugitive CH4 mixer' 
        
C9 = su.Mixer('B9', ins=(C3-2, C2-2), outs=streamsC['N2O'])
C9.specification = lambda: add_fugitive_items(C8, N2O_item)
C9.line = 'fugitive N2O mixer'

################## Non-Treatment ##################
C10 = su.HousingReclaimer('C10', ins=(C7-0), outs='C10_out', ppl=ppl)
C11 = su.SystemReclaimer('C11', ins=(C10-0), outs='C11_out', ppl=ppl, if_gridtied=False)
C12 = su.SolarReclaimer('C12', ins=(C11-0), outs='C12_out')

############### Simulation, TEA, and LCA ###############
sysC = System('sysC', path=(C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12))

teaC = SimpleTEA(system=sysC, discount_rate=discount_rate,
                 start_year=2020, lifetime=20, uptime_ratio=1,
                 lang_factor=None, annual_maintenance=0,
                 annual_labor=((C4._calc_maintenance_labor_cost() + C6._calc_maintenance_labor_cost() + C12._calc_maintenance_labor_cost()) * 24 * 365)
                 )

get_powerC = lambda: sum([u.power_utility.rate for u in sysC.units]) * (24 * 365 * teaC.lifetime)

lcaC = LCA(system=sysC, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerC)


# =============================================================================
# Scenario D (sysD): Targeted Nitrogen Removal Only 
# =============================================================================

# Set flowsheet to avoid stream replacement warnings
flowsheetD = Flowsheet('sysD')
main_flowsheet.set_flowsheet(flowsheetD)
streamsD = batch_create_streams('D')

#################### Human Inputs ####################
D1 = su.Excretion('D1', outs=('urine','feces')) 

###################### Treatment ######################
D3 = su.PrimaryReclaimer('D3', ins=(D1-0, streamsD['MgOH2']),
                         outs=('D3_treated', 'D3_CH4', 'D3_N2O', 'D3_sludge', streamsD['struvite']),
                         decay_k_COD=get_decay_k(tau_deg, log_deg),
                         decay_k_N=get_decay_k(tau_deg, log_deg),
                         max_CH4_emission=max_CH4_emission,
                         ppl=ppl,
                         if_include_front_end=True
                         )

D4 = su.UltrafiltrationReclaimer('D4', ins=(D3-0),
                                 outs=('D4_treated', 'retentate'),
                                 ppl=ppl,
                                 if_gridtied=True
                                 )
                        
D5 = su.IonExchangeReclaimer('D5', ins=(D4-0, streamsD['Zeolite'], streamsD['GAC'], streamsD['KCl']),
                             outs=('D5_treated', 'SpentZeolite', 'SpentGAC',streamsD['Conc_NH3']),
                             decay_k_COD=get_decay_k(tau_deg, log_deg),
                             decay_k_N=get_decay_k(tau_deg, log_deg),
                             max_CH4_emission=max_CH4_emission,
                             ppl=ppl
                             )

D6 = su.Mixer('D6', ins=(D3-1), outs=streamsD['CH4'])
D6.specification = lambda: add_fugitive_items(D6, CH4_item)
D6.line = 'fugitive CH4 mixer' 
        
D7 = su.Mixer('D7', ins=(D3-2), outs=streamsD['N2O'])
D7.specification = lambda: add_fugitive_items(D7, N2O_item)
D7.line = 'fugitive N2O mixer'

################## Non-Treatment ##################
D8 = su.HousingReclaimer('D8', ins=(D5-0), outs='D8_out', ppl=ppl)
D9 = su.SystemReclaimer('D9', ins=(D8-0), outs='D9_out', ppl=ppl, if_gridtied=True)

############### Simulation, TEA, and LCA ###############
sysD = System('sysD', path=(D1, D3, D4, D5, D6, D7, D8, D9))

teaD = SimpleTEA(system=sysD, discount_rate=discount_rate,
                 start_year=2020, lifetime=20, uptime_ratio=1,
                 lang_factor=None, annual_maintenance=0,
                 annual_labor=(D5._calc_maintenance_labor_cost() * 24 * 365)
                 )

get_powerD = lambda: sum([u.power_utility.rate for u in sysD.units]) * (24 * 365 * teaD.lifetime)

lcaD = LCA(system=sysD, lifetime=20, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerD)


# =============================================================================
# Summarizing Functions
# =============================================================================

sys_dct = {
    'ppl':  dict(sysA=ppl, sysB=ppl, sysC=ppl, sysD=ppl),
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
    'front_end':  dict(sysA=None, sysB=(B2,), sysC=(C2,), sysD=None),
    'primary': dict(sysA=(A3,), sysB=(B3,), sysC=(C3,), sysD=(D3,)),
    'sludge_pasteurization': dict(sysA=(A4,), sysB=(B4,), sysC=(C4,), sysD=None),
    'ultrafiltration': dict(sysA=None, sysB=(B5,), sysC=(C5,), sysD=(D4,)),
    'ion_exchange': dict(sysA=None, sysB=(B6,), sysC=(C6,), sysD=(D5,)),
    'ecr': dict(sysA=None, sysB=(B7,), sysC=(C7,), sysD=None),
    'housing': dict(sysA=None, sysB=(B10,), sysC=(C10,), sysD=(D8,)),
    'system': dict(sysA=None, sysB=(B11,), sysC=(C11,), sysD=(D9,)),
    'solar': dict(sysA=None, sysB=None, sysC=(C12,), sysD=None),
    }

system_streams = {sysA:streamsA, sysB:streamsB, sysC:streamsC, sysD:streamsD}

# Learning curve assumptions
percent_CAPEX_to_scale = 0.1  # ratio of cost of specialty parts divided to cost of total parts
get_percent_CAPEX_to_scale = lambda: percent_CAPEX_to_scale

number_of_units = 100000  # assume 100,000 for MURTs

percent_limit = 0.015
get_percent_limit = lambda: percent_limit

learning_curve_percent = 0.925
get_learning_curve_percent = lambda: learning_curve_percent


def get_scaled_capital(tea):
    CAPEX_to_scale = tea.annualized_CAPEX * get_percent_CAPEX_to_scale()
    CAPEX_not_scaled = tea.annualized_CAPEX - CAPEX_to_scale
    scaled_limited = CAPEX_to_scale * get_percent_limit()
    b=(np.log(get_learning_curve_percent())/np.log(2))
    scaled_CAPEX_annualized  = (CAPEX_to_scale - scaled_limited)*number_of_units**b+scaled_limited    
    new_CAPEX_annualized = scaled_CAPEX_annualized + CAPEX_not_scaled    
    return new_CAPEX_annualized

def get_cost_capex(units,tea,ppl):
    if units is None: return 0.
    return tea.get_unit_annualized_equipment_cost(units)/365/ppl

def get_ghg_capex(units,lca,ppl):
    if units is None: return 0.
    ghg = 0
    try: iter(units)
    except: units = (units,)
    for i in units:
        ghg += lca.get_construction_impacts(i).get('GlobalWarming')/lca.lifetime/ppl
    return ghg

def get_cost_opex(units,ppl):
    if units is None: return 0.
    opex = 0
    try: iter(units)
    except: units = (units,)
    for i in units:
        opex += sum(i._add_OPEX.values())*24/ppl
        for stream in i.ins:
            opex += stream.imass['Polyacrylamide']*streamsA['polymer'].price*24/ppl
            opex += stream.imass['MagnesiumHydroxide']*streamsA['MgOH2'].price*24/ppl
            opex += stream.imass['MgCO3']*streamsA['MgCO3'].price*24/ppl
            opex += stream.imass['H2SO4']*streamsA['H2SO4'].price*24/ppl
            opex += stream.imass['FilterBag']*streamsA['filter_bag'].price*24/ppl
            opex += stream.imass['Polystyrene']*streamsA['resin'].price*24/ppl
            opex += stream.imass['GAC']*streamsA['GAC'].price*24/ppl
            opex += stream.imass['Zeolite']*streamsA['Zeolite'].price*24/ppl
            opex += stream.imass['SodiumHydroxide']*streamsA['NaOH'].price*24/ppl
    return opex

def get_cost_electricity(units,ppl):
    if units is None: return 0.
    electricity = 0
    try: iter(units)
    except: units = (units,)
    for i in units:
        electricity+=i.power_utility.cost/ppl
    return electricity

def get_ghg_electricity(units,ppl):
    if units is None: return 0.
    electricity = 0
    ratio = 12 * 365 / ppl  # assuming only using electricity for 12 hr per day
    try: iter(units)
    except: units = (units,)
    for i in units:
        electricity+=i.power_utility.consumption*i.uptime_ratio*e_item.CFs['GlobalWarming']*ratio
        electricity+=-i.power_utility.production*i.uptime_ratio*e_item.CFs['GlobalWarming']*ratio
    return electricity

def get_ghg_direct(units,lca,ppl):
    if units is None: return 0.
    direct = 0
    ratio = 24*365/ppl
    try: iter(units)
    except: units = (units,)
    for i in units:
        for stream in i.outs:
            direct+=stream.imass['CH4']*CH4_item.CFs['GlobalWarming']*ratio
            direct+=stream.imass['N2O']*N2O_item.CFs['GlobalWarming']*ratio
    return direct

def get_ghg_opex(units,lca,ppl):
    if units is None: return 0.
    opex = 0
    try: iter(units)
    except: units = (units,)
    for i in units:
        opex+=lca.get_stream_impacts(i.ins).get('GlobalWarming')/lca.lifetime/ppl
    return opex

def get_total_inputs(unit):
    if unit is None: return 0.
    if len(unit.ins) == 0:  # Excretion units do not have ins
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

def get_fugitive_emissions(streams=None, hr=365*24*20):  # 20 in hr calculation corresponds to 20-year lifetime
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


def get_summarizing_functions():
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
            / lca.lifetime/ppl - (get_fugitive_gases(sys)/lca.lifetime/ppl)  # maybe ignore
    func_dct['get_offset_GWP'] = \
        lambda lca, ppl: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='offset')[ind] \
            / lca.lifetime/ppl
    func_dct['get_other_GWP'] = \
        lambda lca, ppl: lca.total_other_impacts[ind]/lca.lifetime/ppl
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


def save_all_reports():
    for i in (sysA, sysB, sysC, sysD, lcaA, lcaB, lcaC, lcaD):
        if isinstance(i, System):
            i.simulate()
            i.save_report(os.path.join(results_path, f'{i.ID}.xlsx'))
        else:
            i.save_report(os.path.join(results_path, f'{i.system.ID}_lca.xlsx'))


__all__ = ('sysA', 'sysB', 'sysC', 'sysD', 
           'teaA', 'teaB', 'teaC', 'teaD', 
           'lcaA', 'lcaB', 'lcaC', 'lcaD',
           'print_summaries', 'save_all_reports',
            *(i.ID for i in sysA.units),
            *(i.ID for i in sysB.units),
            *(i.ID for i in sysC.units),
            *(i.ID for i in sysD.units)
            )

# This prevents simulating the system when importing
if __name__ == '__main__':
    for sys in (sysA, sysB, sysC, sysD):
        sys.simulate()
