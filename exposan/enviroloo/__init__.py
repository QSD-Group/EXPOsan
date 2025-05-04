#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created by Yuyao Huang and Siqi Tang for Enviroloo Clear Toilet system

    Siqi Tang <siqit@outlook.com>
    Yuyao Huang <yuyaoh20@gmail.com>

'''

import os, qsdsan as qs
from qsdsan import ImpactItem, StreamImpactItem
from exposan.utils import (
    _init_modules,
    get_decay_k,
    get_generic_scaled_capital,
    get_generic_tanker_truck_fee as get_tanker_truck_fee,
    )

## Default settings of resource recovery in the EL system
INCLUDED_RESOURCE_RECOVERY = False

el_path = os.path.dirname(__file__)
#el_data_path = os.path.join(el_path, 'data')
#results_path = os.path.join(el_path, 'results')
module = os.path.split(el_path)[-1]
data_path, results_path = _init_modules(module, include_data_path = True)

#############################################################################################################################################
#
#                          Assign the default values of the parameters in the Enviroloo Clear (EL) system
#
############################################################################################################################################
household_size = 5  # refer to EXPOsan/exposan/pou_disinfection/__init__.py where household size meets normal distribution with mu =4,
                    # and sigma = 1.4 at global scale
household_per_toilet = 20
get_toilet_users = lambda: household_size * household_per_toilet

ppl = 1000 # the number of people served by the EL system. #TOCHANGE
baseline_ppl = 1000 # the number of people served by the EL system. #TOCHANGE
scale_factor = ppl / 100 #scale_factor for flow and dosing rates
dosing_flow = 1 #L/h base scenario

discount_rate = 0.05   # discount rate, [fraction]
#discount_rate = 0.08

max_CH4_emission = 0.25   # max CH4 emission, [g CH4/g COD]

tau_deg = 3 # time taking for full degradation, [yr]
log_deg = 4 # log reduction at full degradation, [log10]

emptying_fee = 0.15 # additional emptying fee, fraction of base cost

# recycled nutrients are sold at a lower price than commercial fertilizers
price_factor = 0.25 # conventional value, if the EL system is considering the nutrients 
                     #recovery as fertilizers sold with lower price than commercial fertilizers.

# The following parameters should be changed based on country of interest
price_ratio = 1.0
operator_daily_wage = 29  # operator daily wage, [USD/day]
const_daily_wage = 17
const_person_days = 100

# The following factors refer to  Reclaimer and Biogenic Refinery
EcosystemQuality_factor = 29320 * (2.8e-09 + 7.65e-14)   # (pt/species.yr) * (species.yr/kgCO2eq)
HumanHealth_factor = 436000 * 9.28e-07   # (pt/DALY) * (DALY/kgCO2eq)

# Ignore the nutrients loss during application

#############################################################################################################################################
#
#                             Define resource recovery settings for the Enviroloo Clear (EL) system
#
#############################################################################################################################################
def update_resource_recovery_settings():
    global INCLUDED_RESOURCE_RECOVERY
    global price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct
    RR_factor = int(bool(INCLUDED_RESOURCE_RECOVERY))
    
    price_dct = {
        'Electricity': 0.0, # $/kWh Assumption because the system is not grid-tied. 
        #'Electricity': 0.13, # $/kWh If grif-tied 
        #'Concrete': 194 * price_ratio, # $/m3
        #'Steel': 2.665 * price_ratio, # $/kg
        'N': 1.507 * price_factor * RR_factor, # $/kg, N fertilizer price if resource recovery is considered
        'P': 3.983 * price_factor * RR_factor, # $/kg, P fertilizer price if resource recovery is considered
        'K': 1.333 * price_factor * RR_factor, # $/kg, K fertilizer price if resource recovery is considered
        'ammonium': 0.8 * (14/18) * price_factor * RR_factor, # $/kg, ammonium fertilizer price if resource recovery is considered
        'struvite':3.983 * (31/245) * price_factor * RR_factor, # $/kg, as alternative P fertilizer 
                                                          #if resource recovery is considered, cited from biogenic refinery
        'NaOH': 1.2 * price_ratio, # $/kg, used for membrane cleanup
        'NaClO': 1.8 * price_ratio, # $/kg, used for membrane cleanup
        'O3': 0.9 * price_ratio, # $/kg, used for clear water tank disinfection
        'HDPE': 0.5 * price_ratio, # $/kg, used for 200M ozone connection
        'PAC': 0 / 0.2886 * price_ratio, #$/kg, kg AlOH/kg PAC conversion
        'Glucose': 0 / 1.067 * price_ratio, # $/kg, used in anoxic tank, kg COD/kg glucose conversion
        'air': 0,
        }
    
    GWP_dct = {
        'Electricity': 0.00, # kg CO2-eq/kWh Assumption because system only contains renewables and is not grid-tied
        #'Electricity': 0.69, # kg CO2-eq/kWh If system is grid-tied
        'CH4': 34.0, # kg CO2-eq/g CH4
        'N2O': 298.0, # kg CO2-eq/g N2O
        'N': -5.4 * RR_factor, # kg CO2-eq/kg N recovered, if resource recovery is considered
        'P': -4.9 * RR_factor, # kg CO2-eq/kg P recovered, if resource recovery is considered
        'K': -1.5 * RR_factor, # kg CO2-eq/kg K recovered, if resource recovery is considered
        'ammonium': -5.4 * (14/18) * RR_factor, # kg CO2-eq/kg ammonium recovered, if resource recovery is considered
        'struvite': -4.9 * (31/245) * RR_factor, # kg CO2-eq/kg, need to make sense, if resource recovery is considered
        'NaOH': 0.050981804, 
        'NaClO': 0.10287725,
        'O3': 0.39586718,
        'PAC': 0.394959 / 0.2886, #kg AlOH/kg PAC
        'Glucose': 0.033883 / 1.067, #kg COD/kg glucose conversion
        'air': 0,
        }
    
    H_Ecosystems_dct = {
        'Electricity': 0.00, #Assumption because system only contains renewables and is not grid-tied
        #'Electricity': 0.002456338, #If system is grid-tied
        'CH4': 34.0 * EcosystemQuality_factor, 
        'N2O': 298.0 * EcosystemQuality_factor, 
        'N': -0.0461961 * RR_factor, #if resource recovery is considered
        'P': -0.093269908 * RR_factor, #if resource recovery is considered
        'K': -0.01895794 * RR_factor, #if resource recovery is considered
        'ammonium': -0.0461961 * (14/18) * RR_factor, #if resource recovery is considered
        'struvite': -0.093269908 * (31/245) * RR_factor, #if resource recovery is considered
        'NaOH': 0.025732326,
        'NaClO': 0.051394080,
        'O3':0.21786132,
        'PAC': 0.1938478 / 0.2886, #kg AlOH/kg PAC
        'Glucose': 0.035664 / 1.067, #kg COD/kg glucose conversion
        'air': 0,
        }
    
    H_Health_dct = {
        'Electricity': 0.00, #Assumption because system only contains renewables and is not grid-tied
        #'Electricity': 0.040824307, If system is grid-tied
        'CH4': 34.0 * HumanHealth_factor,
        'N2O': 298.0 * HumanHealth_factor, 
        'N': -0.637826734 * RR_factor, #if resource recovery is considered
        'P': -1.774294425 * RR_factor, #if resource recovery is considered
        'K': -0.116067637 * RR_factor, #if resource recovery is considered
        'ammonium': -0.637826734 * (14/18) * RR_factor, #if resource recovery is considered
        'struvite': -1.774294425 * (31/245) * RR_factor, #if resource recovery is considered        
        'NaOH': 0.061256093,
        'NaClO': 0.12191047,
        'O3': 0.51710284,
        'PAC': 0.0173485 / 0.2886, #kg AlOH/kg PAC
        'Glucose': 0.011324 / 1.067, #kg COD/kg glucose conversion
        'air': 0,
        }
    
    H_Resources_dct = {
        'Electricity': 0.00, #Assumption because system only contains renewables and is not grid-tied
        #'Electricity': 0.027825633, #If system is grid-tied
        'CH4': 0.0, 
        'N2O': 0.0,
        'N': -0.259196888 * RR_factor, #if resource recovery is considered
        'P': -1.084191599 * RR_factor, #if resource recovery is considered
        'K': -0.054033438 * RR_factor, #if resource recovery is considered
        'ammonium': -0.259196888 * (14/18) * RR_factor, #if resource recovery is considered
        'struvite': -1.084191599 * (31/245) * RR_factor, #if resource recovery is considered
        'NaOH': 0.050981804,
        'NaClO': 0.10287725,
        'O3': 0.39586718,
        'PAC': 0.28485795 / 0.2886, #kg AlOH/kg PAC
        'Glucose': 0.02425 / 1.067, #kg COD/kg glucose conversion
        'air': 0,
        }
    
    return price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct

update_resource_recovery_settings()


#############################################################################################################################################
#
#                             Initialize data to run analysis for the Enviroloo Clear (EL) system
#
#############################################################################################################################################
from . import _components
from ._components import *  
_components_loaded = False
def _load_components(reload=False,):
    global components, _components_loaded
    if not _components_loaded or reload:
        components = create_components()
        qs.set_thermo(components)
        _components_loaded = True
          
_impact_item_loaded = False
def _load_lca_data(reload=False):
    '''
    Load impact indicators and impact item data for LCA
    '''
    global _impact_item_loaded
    if not _impact_item_loaded or reload:
        indicator_path = os.path.join(data_path, 'impact_indicators.csv')
        qs.ImpactIndicator.load_from_file(indicator_path)
                    
        item_path = os.path.join(data_path, 'impact_items.xlsx')
        qs.ImpactItem.load_from_file(item_path)

        price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct = update_resource_recovery_settings()

        # define impacts associated with streams and electricity in the EL system
        def create_stream_impact_item(item_ID, dct_key = '',):
            dct_key = dct_key or item_ID.rsplit('_item')[0] # here for example, change 'struvite_item' to 'struv'
            StreamImpactItem(ID = item_ID, 
                             GWP = GWP_dct[dct_key],
                             H_Ecosystems = H_Ecosystems_dct[dct_key],
                             H_Health = H_Health_dct[dct_key],
                             H_Resources = H_Resources_dct[dct_key],
                             )
        #TODO: stream impact items not properly initiated because 
        # lca.total_stream_impacts['GlobalWarming'] outputs 0.                
        create_stream_impact_item(item_ID='CH4_item') # link when generated in each unit
        create_stream_impact_item(item_ID='N2O_item')
        create_stream_impact_item(item_ID='N_item')
        create_stream_impact_item(item_ID='P_item')
        create_stream_impact_item(item_ID='K_item')
        create_stream_impact_item(item_ID='ammonium_item')
        create_stream_impact_item(item_ID='NaOH_item')
        create_stream_impact_item(item_ID='NaClO_item')
        create_stream_impact_item(item_ID='O3_item')
        create_stream_impact_item(item_ID='PAC_item')
        create_stream_impact_item(item_ID='Glucose_item')
        create_stream_impact_item(item_ID='air_item')
        ImpactItem(ID = 'e_item', functional_unit = 'kWh', 
                   GWP = GWP_dct['Electricity'],
                   H_Ecosystems = H_Ecosystems_dct['Electricity'],
                   H_Health = H_Health_dct['Electricity'],
                   H_Resources = H_Resources_dct['Electricity'])

        # update prices
       # ImpactItem.get_item('Steel').price = price_dct['Steel']     
                             
        _impact_item_loaded = True
                    
    return _impact_item_loaded

######################################################## Load EL system ##############################################################
from . import Enviroloo_system
from .Enviroloo_system import *
_system_loaded = False
def _load_system():
    qs.currency = 'USD'
    qs.PowerUtility.price = price_dct['Electricity']
    global sysEL, teaEL, lcaEL, _system_loaded
    sysEL = create_systemEL('EL')
    teaEL = sysEL.TEA
    lcaEL = sysEL.LCA
    _system_loaded = True

def load():
    
    _load_lca_data()
    if not _components_loaded: _load_components()
    if not _system_loaded: _load_system()
    
    
    dct = globals()
    
    for sys in (sysEL,): #dct.update(sys.flowsheet.to_dct())
        flowsheet = sys.flowsheet
        units = {unit.ID: unit for unit in flowsheet.unit}
        streams = {stream.ID: stream for stream in flowsheet.stream}
        dct.update(units)
        dct.update(streams)
                    
def __getattr__(name):
    if not _components_loaded or not _system_loaded:
        raise AttributeError(
            f'Module {__name__} does not have the attribute "{name}" '
            'and the module has not been loaded, '
            f'loading the module with `{__name__}.load()` may solve the problem.')
                                           
######################################################## Utilities functions ##############################################################
## The amount of nutrients can be recovered per year if recovery setting mode is TRUE.
hr_per_yr = 24 * 365  # full payload per day
# if considering C recovery in COD
#def get_C(stream):
#    sys = stream.source.system
#    for unit in sys.units:
#        if hasattr(unit, 'carbon_COD_ratio'):
#            carbon_COD_ratio = unit.carbon_COD_ratio
#    return stream.COD * stream.F_vol/1e3 * carbon_COD_ratio * hr_per_yr
#get_C_gas = lambda stream: stream.imol['CH4'] * 12 * hr_per_yr # break down C flow in gas
#get_N_gas = lambda stream: stream.imol['N2O'] * 28 * hr_per_yr # break down N flow in gas
get_N = lambda stream: stream.TN * stream.F_vol/1e3 * hr_per_yr
get_P = lambda stream: stream.TP * stream.F_vol/1e3 * hr_per_yr
get_K = lambda stream: stream.TK * stream.F_vol/1e3 * hr_per_yr

def get_recoveries(system, include_breakdown=False):
    # # EL = system.ID[-1]
    # EL = system.ID
    # if EL not in ('E', 'L', 'EL'):
    #     raise ValueError('This function is only available for the Enviroloo Clear system `sysEL`,'
    #                      f'not `{system.ID}`.')
    
    
    u_reg = system.flowsheet.unit
    dct = globals()
    dct['N_dct'] = N_dct = {}
    dct['P_dct'] = P_dct = {}
    dct['K_dct'] = K_dct = {}
    #dct['C_dct'] = C_dct = {}
    
    # if EL == ('E', 'L'):
#    toilet = u_reg.Toilet  # here the name 'Toilet' should be consistent with the name defined in the Enviroloo_system.py
    # CollectionTank = u_reg.CT
    # PrimaryClarifier = u_reg.PC
    AnoxicTank = u_reg.A1
    AerobicTank = u_reg.O1
    MembraneTank = u_reg.B1
    # ClearWaterTank = u_reg.CWT
    
    # # In the Toilet unit
    # #C_dct['urine'] = get_C(toilet.ins[0]) * ppl
    # #C_dct['feces'] = get_C(toilet.ins[1]) * ppl
    # #C_dct['input'] = C_dct['urine'] + C_dct['feces']
    # N_dct['urine'] = get_N(toilet.ins[0]) * ppl
    # N_dct['feces'] = get_N(toilet.ins[1]) * ppl
    # N_dct['input'] = N_dct['urine'] + N_dct['feces']
    # P_dct['urine'] = get_P(toilet.ins[0]) * ppl
    # P_dct['feces'] = get_P(toilet.ins[1]) * ppl
    # P_dct['input'] = P_dct['urine'] + P_dct['feces']
    # K_dct['urine'] = get_K(toilet.ins[0]) * ppl
    # K_dct['feces'] = get_K(toilet.ins[1]) * ppl
    # K_dct['input'] = K_dct['urine'] + K_dct['feces']

    #C_dct['toilet_gas'] = get_C_gas(toilet.outs[-1])  # CH4 in the second order of toilet outs
    #N_dct['toilet_gas'] = get_N_gas(toilet.outs[-2])  # N2O in the third order of toilet outs
    
    # consider sludge in membrane tank and aerobic tank
    # N_dct['treated_sludge'] = get_N(MembraneTank.outs[0]) + get_N(AerobicTank.outs[0]) + get_N(AnoxicTank.outs[0])
    # P_dct['treated_sludge'] = get_P(MembraneTank.outs[0]) + get_P(AerobicTank.outs[0]) + get_P(AnoxicTank.outs[0])
    # K_dct['treated_sludge'] = get_K(MembraneTank.outs[0]) + get_K(AerobicTank.outs[0]) + get_K(AnoxicTank.outs[0])

    # the code in the following domain is a sample, which will be updated later
    # functions = [
    #     lambda: N_dct['treated_sludge'] / N_dct['input'] * 100, # total N recovery percentage
    #     lambda: P_dct['treated_sludge'] / P_dct['input'] * 100, # total P recovery percentage
    #     lambda: K_dct['treated_sludge'] / K_dct['input'] * 100, # total K recovery percentage
    #     ]
    # if not include_breakdown: return functions
    
    # return [
    #     *functions,
    #     lambda: N_dct['treated_sludge'] / N_dct['input'] * 100, # N recovery percentage
    #     lambda: P_dct['treated_sludge'] / P_dct['input'] * 100, # P recovery percentage
    #     lambda: K_dct['treated_sludge'] / K_dct['input'] * 100, # K recovery percentage
    #     ]
        
######################################################## Financial and Cost parameters ##############################################################
## 
percent_CAPEX_to_scale = 0.1
number_of_units = 10000
percent_limit = 0.015
learning_curve_percent = 0.90  # assume learning curve
def get_scaled_capital(tea):
    new_CAPEX_annualized = get_generic_scaled_capital(
        tea = tea,
        percent_CAPEX_to_scale = percent_CAPEX_to_scale,
        number_of_units = number_of_units,
        percent_limit = percent_limit,
        learning_curve_percent = learning_curve_percent
        )
    return new_CAPEX_annualized

def get_TEA_metrics(system, include_breakdown=True):
    tea = system.TEA
    get_annual_electricity = lambda system: system.power_utility.cost * system.operating_hours
    functions = [lambda: (get_scaled_capital(tea) - tea.net_earnings) / ppl]
    if not include_breakdown: return functions # net cost
    return [
        *functions,
        lambda: get_scaled_capital(tea) / ppl, # means CAPEX
        lambda: get_annual_electricity(system) / ppl, # means annual electricity consumption
        lambda: tea.annual_labor / ppl, # means annual labor cost
        lambda: (tea.AOC - get_annual_electricity(system) - tea.annual_labor) / ppl, # means OPEX excluding energy and labor sectors
        lambda: tea.net_earnings / ppl, # means net earning income
        ]

def get_TEA_metrics_breakdown(system, include_breakdown=True):
    tea = system.TEA
    get_annual_electricity = lambda system: system.power_utility.cost * system.operating_hours
    metrics = {
        'CAPEX': get_scaled_capital(tea) / ppl,
        'Annual electricity': get_annual_electricity(system) / ppl,
        'Annual labor': tea.annual_labor / ppl,
        'OPEX (excl. energy & labor)': (tea.AOC - get_annual_electricity(system) - tea.annual_labor) / ppl,
        'Sales Income':  tea.net_earnings / ppl
    }
    
    net_cost = (get_scaled_capital(tea) - tea.net_earnings) / ppl
    
    if not include_breakdown:
        return net_cost
    
    total_cost = metrics['CAPEX'] + metrics['Annual electricity'] + \
                 metrics['Annual labor'] + metrics['OPEX (excl. energy & labor)']
    
    print("\nTEA Metrics Breakdown:")
    print("-" * 75)
    print(f"{'Metric':<35} {'Value':>10} {'% of Total':>14} {'Unit':>15}")
    print("-" * 75)
    
    # print cost breakdown
    for metric_name in ['CAPEX', 'Annual electricity', 'Annual labor', 'OPEX (excl. energy & labor)', 'Sales Income']:
        value = metrics[metric_name]
        percentage = (value / total_cost * 100) if total_cost != 0 else 0
        print(f"{metric_name:<35} {value:>10.2f} {percentage:>13.1f}% {'USD/cap/yr':>15}")
    
    print("-" * 75)
    print(f"{'Total Cost':<35} {total_cost:>10.2f} {100:>13.1f}% {'USD/cap/yr':>15}")
#    print(f"{'Sales Income':<35} {metrics['Sales Income']:>10.2f} {' ':>13} {'USD/cap/yr':>15}")
    print(f"{'Net Cost':<35} {net_cost:>10.2f} {' ':>13} {'USD/cap/yr':>15}")
    print("-" * 75)
    
def get_normalized_CAPEX(units): # get the CAPEX of a unit or units normalized to per capita per day
    system = units[0].system
    return system.TEA.get_unit_annualized_equipment_cost(units) / 365 / ppl

def get_normalized_electricity_cost(units): # get the electricity cost of a unit or units normalized to per capita per day
    return sum(u.power_utility.cost for u in units) / ppl

def get_normalized_OPEX(units): # get the OPEX of a unit or units normalized to per capita per day, where energy (electricity) cost is excluded.
    OPEX = sum(u.add_OPEX.values() for u in units)
    streams = sum([u.ins for u in units], [])
    OPEX += sum(s.cost for s in streams)
    return OPEX * 24 / ppl  # converted to per capita per day

def get_LCA_metrics(system, include_breakdown=False):
    lca = system.LCA
    functions = [
        lambda: lca.total_impacts['GlobalWarming'] / lca.lifetime / ppl, # annual GWP
        # ReCiPe LCA functions
        lambda: lca.total_impacts['H_Ecosystems'] / lca.lifetime / ppl,
        lambda: lca.total_impacts['H_Health'] / lca.lifetime / ppl,
        lambda: lca.total_impacts['H_Resources'] / lca.lifetime / ppl,
        ]
    if not include_breakdown: return functions
    return [
        *functions,
        lambda: lca.total_construction_impacts['GlobalWarming'] / lca.lifetime / ppl, # means construction fee
        lambda: lca.total_transportation_impacts['GlobalWarming'] / lca.lifetime / ppl, # means transportation fee
        lambda: lca.total_stream_impacts['GlobalWarming'] / lca.lifetime / ppl, # means stream impacts including fugitive gases and offsets
        lambda: lca.total_other_impacts['GlobalWarming'] / lca.lifetime / ppl, # means other impacts in GWP
        ]

def get_LCA_metrics_breakdown(system, include_breakdown=True):
    lca = system.LCA
    metrics = {
        'Annual GWP': {
            'value': lca.total_impacts['GlobalWarming'] / lca.lifetime / ppl,
            'unit': 'kg CO2-eq/cap/yr'
        },
        'H_Ecosystems': {
            'value': lca.total_impacts['H_Ecosystems'] / lca.lifetime / ppl,
            'unit': 'points/cap/yr'
        },
        'H_Health': {
            'value': lca.total_impacts['H_Health'] / lca.lifetime / ppl,
            'unit': 'points/cap/yr'
        },
        'H_Resources': {
            'value': lca.total_impacts['H_Resources'] / lca.lifetime / ppl,
            'unit': 'points/cap/yr'
        }
    }

    if not include_breakdown:
        return metrics['Annual GWP']['value']

    # GWP breakdown
    gwp_breakdown = {
        'Construction': lca.total_construction_impacts['GlobalWarming'] / lca.lifetime / ppl,
        'Transportation': lca.total_transportation_impacts['GlobalWarming'] / lca.lifetime / ppl,
        'Stream impacts': lca.total_stream_impacts['GlobalWarming'] / lca.lifetime / ppl,
        'Other impacts': lca.total_other_impacts['GlobalWarming'] / lca.lifetime / ppl
    }
    total_gwp = sum(gwp_breakdown.values())

    print("\nLCA Metrics Breakdown:")
    print("-" * 75)
    print(f"{'Metric':<35} {'Value':>10} {'% of Total':>14} {'Unit':>15}")
    print("-" * 75)

    # print GWP breakdown
    for metric_name, value in gwp_breakdown.items():
        percentage = (value / total_gwp * 100) if total_gwp != 0 else 0
        print(f"{metric_name:<35} {value:>10.2f} {percentage:>13.1f}% {'kg CO2-eq/cap/yr':>15}")

    print("-" * 75)
    print(f"{'Total GWP':<35} {total_gwp:>10.2f} {100:>13.1f}% {'kg CO2-eq/cap/yr':>15}")
    print("-" * 75)
    print("Other Impact Categories:")
    for metric_name in ['H_Ecosystems', 'H_Health', 'H_Resources']:
        value = metrics[metric_name]['value']
        unit = metrics[metric_name]['unit']
        print(f"{metric_name:<35} {value:>10.2f} {' ':>13} {unit:>15}")
    print("-" * 75)
               
def print_summaries(systems):
    try: iter(systems)
    except: systems = (systems, )

    for sys in systems:
        sys.simulate()
        print(f'\n-----------------Summary for {sys.ID}-----------------\n')
        if sys.ID in ('sysEL', ):
            recovery_functions = get_recoveries(sys)
            print(f'\nTotal N recovery: {recovery_functions[0]():.1f} %.')
            print(f'\nTotal P recovery: {recovery_functions[1]():.1f} %.')
            print(f'\nTotal K recovery: {recovery_functions[2]():.1f} %.')
            
            TEA_functions = get_TEA_metrics(sys, include_breakdown=False)
            unit = f'{qs.currency}/cap/yr'
            print(f'\nTotal cost: {TEA_functions[0]():.2f} {unit}.')
                    
            LCA_functions = get_LCA_metrics(sys, include_breakdown=False)
            print(f'\nNet emission: {LCA_functions[0]():.2f} kg CO2-eq/cap/yr.')

            unit = 'points/cap/yr' # breakdown of Impact Indicators
            print(f'\nNet ecosystems damage: {LCA_functions[1]():.2f} {unit}.')
            print(f'\nNet health damage: {LCA_functions[2]():.2f} {unit}.')
            print(f'\nNet resources damage: {LCA_functions[3]():.2f} {unit}.')
        else:
            sys.TEA.show()
            print('\n')
            sys.LCA.show()


from . import Enviroloo_model
from .Enviroloo_model import *

# This country_specific file will be done after the LCA, TEA, and uncertainty and sensitivity analysis of the EL system are all conducted.
#from exposan.enviroloo import country_specific
#from .country_specific import *

__all__ = (
    'el_path',
    'el_data_path',      
    'results_path',      
    *_components.__all__,
    *Enviroloo_system.__all__,
    *Enviroloo_model.__all__,
    #*country_specific.__all__,
    )