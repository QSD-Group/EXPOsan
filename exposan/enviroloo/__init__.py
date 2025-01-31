#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created by Yuyao Huang and Siqi Tang for Enviroloo Clear Toilet system
'''
# %%
import os, qsdsan as qs

#%%
from math import log
from qsdsan import ImpactItem, StreamImpactItem
from exposan.utils import (
          _init_modules,
          get_decay_k,
          get_generic_scaled_capital,
          get_generic_tanker_truck_fee as get_tanker_truck_fee,
          )

EL_path = os.path.dirname(__file__)
module = os.path.split(EL_path)[-1]
data_path, results_path = _init_modules(module, include_data_path = True)

# Default settings of resource recovery in the EL system
INCLUDED_RESOURCE_RECOVERY = False


#############################################################################################################################################
#
#                          Assign the default values of the parameters in the Enviroloo Clear (EL) system
#
############################################################################################################################################
household_size = 4  # refer to EXPOsan/exposan/pou_disinfection/__init__.py where household size meets normal distribution with mu =4,
                    # and sigma = 1.4 at global scale
household_per_toilet = 4
get_toilet_users = lambda: household_size * household_per_toilet

ppl = 1000; # the number of people served by the EL system. 
# Here 100-user household or 1000-user school scale will be considered.

discount_rate = 0.03   # discount rate, [fraction]

max_CH4_emission = 0.5   # max CH4 emission, [g CH4/g COD]

tau_deg = 3; # time taking for full degradation, [yr]
log_deg = 4; # log reduction at full degradation, [log10]

emptying_fee = 0.15; # additional emptying fee, fraction of base cost

# recycled nutrients are sold at a lower price than commercial fertilizers
price_factor = 0.25; # conventional value, if the EL system is considering the nutrients 
                     #recovery as fertilizers sold with lower price than commercial fertilizers.

# The following parameters should be changed based on country of interest
price_ratio = 1.0
operator_daily_wage = 29  # operator daily wage, [USD/day]
constant_daily_wage = 17
const_person_days = 100   # the working days for the constant staff, [day] (???)

# The following factors refer to  Reclaimer and Biogenic Refinery
EcosystemQuality_factor = 29320 * (2.8e-09 + 7.65e-14)   # (pt/species.yr) * (species.yr/kgCO2eq)
HumanHealth_factor = 436000 * 9.28e-07   # (pt/DALY) * (DALY/kgCO2eq)


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
        'Electricity': 0.13, # $/kWh
        'Concrete': 194, # $/m3
        'Steel': 1.8, # $/kg
        'N': 1.507 * price_factor * RR_factor, # $/kg, N fertilizer price if resource recovery is considered
        'P': 0.8 * price_factor * RR_factor, # $/kg, P fertilizer price if resource recovery is considered
        'K': 0.8 * price_factor * RR_factor, # $/kg, K fertilizer price if resource recovery is considered
        'Ammonium': 0.8 * price_factor * RR_factor, # $/kg, ammonium fertilizer price if resource recovery is considered
        'Struvite':3.983 * (31/245) * price_factor * RR_factor, # $/kg, as alternative P fertilizer 
                                                          #if resource recovery is considered, cited from biogenic refinery
        'NaOH': 1.2, # $/kg, used for membrane cleanup
        'NaClO': 1.8, # $/kg, used for membrane cleanup
        'O3': 0.9, # $/kg, used for clear water tank disinfection
        'HDPE': 0.5, # $/kg, used for 200M ozone connection
        }
    
    GWP_dct = {
        'Electricity': 0.69, # kg CO2-eq/kWh
        'Concrete': 0.0, # kg CO2-eq/m3
        'Steel': 0.0, # kg CO2-eq/kg, need to make sense
        'CH4': 34.0, # kg CO2-eq/g CH4
        'N2O': 298.0, # kg CO2-eq/g N2O
        'N': -5.4 * RR_factor, # kg CO2-eq/kg N recovered, if resource recovery is considered
        'P': -4.9 * RR_factor, # kg CO2-eq/kg P recovered, if resource recovery is considered
        'K': -1.5 * RR_factor, # kg CO2-eq/kg K recovered, if resource recovery is considered
        'Ammonium': -5.4 * RR_factor, # kg CO2-eq/kg ammonium recovered, if resource recovery is considered
        'Struvite': -0.3 * RR_factor, # kg CO2-eq/kg, need to make sense, if resource recovery is considered
        'NaOH': 0.0, # kg CO2-eq/kg, need to make sense
        'NaClO': 0.0, # kg CO2-eq/kg, need to make sense
        }
    
    H_Ecosystems_dct = {
        'Electricity': 0.0, # kg CO2-eq/kWh
        'Concrete': 0.0, # kg CO2-eq/m3
        'Steel': 0.0, # kg CO2-eq/kg, need to make sense
        'CH4': 34.0, # kg CO2-eq/g CH4
        'N2O': 298.0, # kg CO2-eq/g N2O
        'N': -5.4 * RR_factor, # kg CO2-eq/kg N recovered, if resource recovery is considered
        'P': -4.9 * RR_factor, # kg CO2-eq/kg P recovered, if resource recovery is considered
        'K': -1.5 * RR_factor, # kg CO2-eq/kg K recovered, if resource recovery is considered
        'Ammonium': -5.4 * RR_factor, # kg CO2-eq/kg ammonium recovered, if resource recovery is considered
        'Struvite': -0.3 * RR_factor, # kg CO2-eq/kg, need to make sense, if resource recovery is considered
        'NaOH': 0.0, # kg CO2-eq/kg, need to make sense
        'NaClO': 0.0, # kg CO2-eq/kg, need to make sense
        }
    
    H_Health_dct = {
        'Electricity': 0.0, # kg CO2-eq/kWh
        'Concrete': 0.0, # kg CO2-eq/m3
        'Steel': 0.0, # kg CO2-eq/kg, need to make sense
        'CH4': 34.0, # kg CO2-eq/g CH4
        'N2O': 298.0, # kg CO2-eq/g N2O
        'N': -5.4 * RR_factor, # kg CO2-eq/kg N recovered, if resource recovery is considered
        'P': -4.9 * RR_factor, # kg CO2-eq/kg P recovered, if resource recovery is considered
        'K': -1.5 * RR_factor, # kg CO2-eq/kg K recovered, if resource recovery is considered
        'Ammonium': -5.4 * RR_factor, # kg CO2-eq/kg ammonium recovered, if resource recovery is considered
        'Struvite': -0.3 * RR_factor, # kg CO2-eq/kg, need to make sense, if resource recovery is considered        
        'NaOH': 0.0, # kg CO2-eq/kg, need to make sense
        'NaClO': 0.0, # kg CO2-eq/kg, need to make sense
        }
    
    H_Resources_dct = {
        'Electricity': 0.0, # kg CO2-eq/kWh
        'Concrete': 0.0, # kg CO2-eq/m3
        'Steel': 0.0, # kg CO2-eq/kg, need to make sense
        'CH4': 34.0, # kg CO2-eq/g CH4
        'N2O': 298.0, # kg CO2-eq/g N2O
        'N': -5.4 * RR_factor, # kg CO2-eq/kg N recovered, if resource recovery is considered
        'P': -4.9 * RR_factor, # kg CO2-eq/kg P recovered, if resource recovery is considered
        'K': -1.5 * RR_factor, # kg CO2-eq/kg K recovered, if resource recovery is considered
        'Ammonium': -5.4 * RR_factor, # kg CO2-eq/kg ammonium recovered, if resource recovery is considered
        'Struvite': -0.3 * RR_factor, # kg CO2-eq/kg, need to make sense, if resource recovery is considered
        'NaOH': 0.0, # kg CO2-eq/kg, need to make sense 
        'NaClO': 0.0, # kg CO2-eq/kg, need to make sense
        }
    return price_dct, GWP_dct, H_Health_dct, H_Resources_dct

update_resource_recovery_settings()


#############################################################################################################################################
#
#                             Initialize data to run analysis for the Enviroloo Clear (EL) system
#
#############################################################################################################################################
from . import _components
from ._components import *  
_components_loaded = False
def _load_components(reload = False):
    global components, _components_loaded
    if not _components_loaded or reload:
        components = create_components()
        qs.set_thermo(components)
        _components_loaded = True
          
_impact_item_loaded = False
def _load_lca_data(reload = False):
    '''
    Load impact indicators and impact item data for LCA
    '''
    global _impact_item_loaded
    if not _impact_item_loaded or reload:
        indicator_path = os.path.join(data_path, 'impact_indicators.csv')
        qs.ImpactIndicator.load_from_file(indicator_path)
                    
        item_path = os.path.join(data_path, 'impact_items.xlsx')
        qs.ImpactItem.load_from_file(item_path)
                    
        # define impacts associated with streams and electricity in the EL system
        def create_stream_impact_item(item_ID, dct_key = ''):
            dct_key = dct_key or item_ID.rsplit('_item')[0];# change 'struvite_item' to 'struv'
            StreamImpactItem(ID = item_ID, GWP = GWP_dct[dct_key],
                             H_Ecosystems = H_Ecosystems_dct[dct_key],
                             H_Health = H_Health_dct[dct_key],
                             H_Resources = H_Resources_dct[dct_key])
                        
        create_stream_impact_item(item_ID = 'CH4_item')
        create_stream_impact_item(item_ID = 'N2O_item')
        create_stream_impact_item(item_ID = 'N_item')
        create_stream_impact_item(item_ID = 'P_item')
        create_stream_impact_item(item_ID = 'K_item')
        create_stream_impact_item(item_ID = 'Concrete_item')
        create_stream_impact_item(item_ID = 'Struvite_item')
        create_stream_impact_item(item_ID = 'Salt_item')
        create_stream_impact_item(item_ID = 'NaOH_item')
        create_stream_impact_item(item_ID = 'NaClO_item')
        ImpactItem(ID = 'e_item', functional_unit = 'kWh', 
                   GWP = GWP_dct['Electricity'],
                   H_Ecosystems = H_Ecosystems_dct['Electricity'],
                   H_Health = H_Health_dct['Electricity'],
                   H_Resources = H_Resources_dct['Electricity'])     
                         
                    
        _impact_item_loaded = True
                    
    return _impact_item_loaded

######################################################## Load EL system ##############################################################
from exposan.enviroloo import Enviroloo_system  # the name of imported module will be aligned finally
from .Enviroloo_system import create_system
_system_loaded = False
def _load_system():
    qs.currency = 'USD'
    qs.PowerUnity.price = price_dct['Electricity']
    global sysEL, teaEL, lcaEL, _system_loaded
    sysEL = create_system()
    teaEL = sysEL.TEA
    lcaEL = sysEL.LCA
    _system_loaded = True

def load():
    if not _components_loaded:
        _load_components()
    if not _system_loaded:
        _load_system()     
    dct = globals()
    for sys in (sysEL):
        dct.update(sys.flowsheet.to_dct())
                    
def _getattr_(name):
    if not _components_loaded or not _system_loaded:
        raise AttributeError(
            f'Module {__name__} does not have the attribute "{name}" '
            'and the module has not been loaded, '
            f'loading the module with `{__name__}.load()` may solve the problem.')
                                           
######################################################## Utilities functions ##############################################################
## The amount of nutrients can be recovered per year if recovery setting mode is TRUE.
hr_per_yr = 24 * 365; # full payload per day
get_N = lambda stream: stream.TN * stream.F_vol/1e3 * hr_per_yr
get_P = lambda stream: stream.TP * stream.F_vol/1e3 * hr_per_yr
get_K = lambda stream: stream.TK * stream.F_vol/1e3 * hr_per_yr
# if considering C recovery in COD
def get_C(stream):
    sys = stream.source.system
    for unit in sys.units:
        if hasattr(unit, 'carbon_COD_ratio'):
            carbon_COD_ratio = unit.carbon_COD_ratio
    return stream.COD * stream.F_vol/1e3 * carbon_COD_ratio * hr_per_yr
get_C_gas = lambda stream: stream.imol['CH4'] * 12 * hr_per_yr; # break down C flow in gas
get_N_gas = lambda stream: stream.imol['N2O'] * 44 * hr_per_yr; # break down N flow in gas
# here the code domain only serves the EL system, so a reminder is designed
def get_recoveries(system, include_breakdown = False):
    EL = system.ID[-1]
    if EL not in ('E', 'L'):
        raise ValueError('This function is only available for the Enviroloo Clear system `sysEL`,'
                         f'not `{system.ID}`.')
    
    
    u_reg = system.flowsheet.unit
    dct = globals()
    dct['N_dct'] = N_dct = {}
    dct['P_dct'] = P_dct = {}
    dct['K_dct'] = K_dct = {}
    dct['C_dct'] = C_dct = {}
    
    if EL == ('E', 'L'):
        toilet = u.reg.Toilet; # here the name 'Toilet' should be consistent with the name defined in the Enviroloo_system.py
        CollectionTank = u_reg.CT
        PrimaryClarifier = u_reg.PC
        AnoxicTank = u_reg.AnoT
        AerobicTank = u_reg.AerT
        MembraneTank = u_reg.MemT
        ClearWaterTank = u_reg.CWT
        PressureTank = u_reg.PT
    
    # In the Toilet unit
    N_dct['urine'] = get_N(toilet.ins[0]) * ppl
    N_dct['feces'] = get_N(toilet.ins[1]) * ppl
    N_dct['input'] = N_dct['urine'] + N_dct['feces']
    P_dct['urine'] = get_P(toilet.ins[0]) * ppl
    P_dct['feces'] = get_P(toilet.ins[1]) * ppl
    P_dct['input'] = P_dct['urine'] + P_dct['feces']
    K_dct['urine'] = get_K(toilet.ins[0]) * ppl
    K_dct['feces'] = get_K(toilet.ins[1]) * ppl
    K_dct['input'] = K_dct['urine'] + K_dct['feces']
    C_dct['urine'] = get_C(toilet.ins[0]) * ppl
    C_dct['feces'] = get_C(toilet.ins[1]) * ppl
    C_dct['input'] = C_dct['urine'] + C_dct['feces']
    
    N_dct['toilet_gas'] = get_N_gas(toilet.outs[-1])  # N2O in the second order of toilet outs
    C_dct['toilet_gas'] = get_C_gas(toilet.outs[-2])  # CH4 in the third order of toilet outs
    
    ###########  Add others after completing Enviroloo_system.py - - - - - - - - - - - - - ->
    
    # the code in the following domain is a sample, which will be updated later
    functions = [
        lambda: (N_dct['treated_sludge'] + N_dct['Ammonium']) / N_dct['input'] * 100, # total N recovery percentage
        lambda: (P_dct['treated_sludge'] + P_dct['Ammonium']) / P_dct['input'] * 100, # total P recovery percentage
        lambda: (K_dct['treated_sludge'] + K_dct['Ammonium']) / K_dct['input'] * 100, # total K recovery percentage
    ]
    if not included_breakdown:
        return functions
    
    return [
        *functions,
        lambda: N_dct['treated_sludge'] / N_dct['input'] * 100, # N recovery percentage
        lambda: P_dct['treated_sludge'] / P_dct['input'] * 100, # P recovery percentage
        lambda: K_dct['treated_sludge'] / K_dct['input'] * 100, # K recovery percentage
    ]
        
######################################################## Financial and Cost parameters ##############################################################
## 
percent_CAPEX_to_scale = 0.1
number_of_units = 1000  # make sense the meanings of '1000' (here as an example referring to Reclaimer and Biogenic Refinery)
percent_limit = 0.012
learning_curve_percent = 0.930  # assume learning curve
def get_scaled_capital(tea):
    return get_generic_scaled_capital(
        tea = tea,
        percent_CAPEX_to_scale = percent_CAPEX_to_scale,
        number_of_units = number_of_units,
        percent_limit = percent_limit,
        learning_curve_percent = learning_curve_percent
        )

def get_TEA_metrics(system, include_breakdown = False):
    tea = system.TEA
    get_annual_electricity = lambda system: system.power_utility.cost * system.operating_hours
    functions = [lambda: (get_scaled_capital(tea) - tea.net_earnings) / ppl]
    if not include_breakdown:
        return functions
    return [
        *functions,
        lambda: get_annual_electricity(system) / ppl, # means CAPEX
        lambda: get_annual_electrictity(system) / ppl, # means annual electricity consumption
        lambda: tea.annual_labor_cost / ppl, # means annual labor cost
        lambda: (tea.AOC - get_annual_electricity() - tea.annual_labor) / ppl, # means OPEX excluding energy and labor sectors
        lambda: tea.sales / ppl, # means sales incoming
        ]

def get_normalized_CAPEX(units): # get the CAPEX of a unit or units normalized to per capita per day
    system = units[0].system
    return system.TEA.get_unit_annualized_equipment_cost(units) / 365 / ppl

def get_normalized_electricity_cost(units): # get the electricity cost of a unit or units normalized to per capita per day
    return sum(u.power_utility.cost for u in units) / ppl

def get_normalized_OPEX(units): # get the OPEX of a unit or units normalized to per capita per day, where energy (electricity) cost is excluded.
    OPEX = sum(u.add_OPEX.values() for u in units)
    streams = sum([u.ins for u in units], [])
    OPEX += sum(s.cost for s in streams)
    return OPEX * 24 / ppl; # converted to per capita per day

def get_LCA_metrics(system, include_breakdown = False):
    lca = system.LCA
    functions = [
        lambda: lca.total_impacts['GlobalWarming'] / lca.lifetime / ppl, # annual GWP
        lambda: lca.total_impacts['H_Ecosystems'] / lca.lifetime / ppl,
        lambda: lca.total_impacts['H_HumanHealth'] / lca.lifetime / ppl,
        lambda: lca.total_impacts['H_Resources']/ lca.lifetime / ppl
        ]
    if not include_breakdown:
        return functions
    return [
        *functions,
        lambda: lca.total_construction_impacts['GlobalWarming'] / lca.lifetime / ppl, # means construction fee
        lambda: lca.total_transportation_impacts['GlobalWarming'] / lca.lifetime / ppl, # means transportation fee
        lambda: lca.total_stream_impacts['GlobalWarming'] / lca.lifetime / ppl, # means stream impacts including fugitive gases and offsets
        lambda: lca.total_other-impacts['GlobalWarming'] / lca.lifetime / ppl # means other impacts in GWP
        ]
               
def print_summaries(system):
    if sys == system:
        sys.simulate()
                   
        print(f'\n-----------------Summary for {sys.ID}-----------------\n')
        TEA_functions = get_TEA_metrics(sys);
        unit = f'{qs.currency} / cap / yr';
        print(f'\nTotal cost: {TEA_functions[0]():.2f} {unit}. ')
                    
        LCA_functions = get_LCA_metrics(sys);
        print(f'\nNet emission: {LCA_functions[0]():.2f} kg CO2-eq / cap / yr. ')


from exposan.enviroloo import Enviroloo_model
from .Enviroloo_model import *

# This country_specific file will be done after the LCA, TEA, and uncertainty and sensitivity analysis of the EL system are all conducted.
from exposan.enviroloo import country_specific
from .country_specific import *

__all__ = (
    'EL_path',
    'data_path',      
    'results_path',      
    *_components.__all__,
    *Enviroloo_system.__all__,
    *Enviroloo_model.__all__,
    #*country_specific.__all__,
    )