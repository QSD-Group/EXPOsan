#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

import os, qsdsan as qs
from qsdsan import ImpactItem, StreamImpactItem
from collections.abc import Iterable
from exposan.utils import (
    get_decay_k,
    get_generic_scaled_capital,
    get_generic_tanker_truck_fee as get_tanker_truck_fee,
    )

#%%
# Module-wise setting on whether to allow resource recovery
INCLUDE_RESOURCE_RECOVERY = False

g2rt_path = os.path.dirname(__file__)
g2rt_data_path = os.path.join(g2rt_path, 'data')
results_path = os.path.join(g2rt_path, 'results')
# To save simulation data
if not os.path.isdir(results_path): os.mkdir(results_path)
default_ppl = 6

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3

max_CH4_emission = 0.25

# =============================================================================
# Prices and GWP CFs
# =============================================================================

# Recycled nutrients & water are sold at a lower price than commercial counterparts
price_factor = 0.25

# Should be changed based on country
price_ratio = 1

EcosystemQuality_factor = 29320 * (2.8e-09+7.65e-14)  # (pt/species.yr) * (species.yr/kgCO2eq)
HumanHealth_factor = 436000 * 9.28e-07  # (pt/DALY) * (DALY/kgCO2eq)

def update_resource_recovery_settings():
    global INCLUDE_RESOURCE_RECOVERY
    global price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct
    RR_factor = int(bool(INCLUDE_RESOURCE_RECOVERY))
    price_dct = {
        'Electricity': 0.13,
        'N': 1.507*price_factor*RR_factor,
        'P': 3.983*price_factor*RR_factor,
        'K': 1.333*price_factor*RR_factor,
        'H2O': 4.01/1e3*price_ratio*RR_factor,
        #per kg of water + wastewater treated https://www.epa.gov/watersense/data-and-information-used-watersense
        'wages': 3.64,
        }
    
    GWP_dct = {
        'Electricity': 0.69,
        'CH4': 34,
        'N2O': 298,
        'N': -5.4*RR_factor,
        'P': -4.9*RR_factor,
        'K': -1.5*RR_factor,
        'H2O': -3.2756/1e4 * RR_factor, 
        # per kg of tap water production, conventional treatment
        'NO': 0,
        'SO2':0,
        'NH3': 0
        }
    
    H_Ecosystems_dct = {
        'Electricity': 0.002456338,
        'CH4': 34 * EcosystemQuality_factor,
        'N2O': 298 * EcosystemQuality_factor,
        'N': -0.0461961*RR_factor,
        'P': -0.093269908*RR_factor,
        'K': -0.01895794*RR_factor,
        'H2O': -6.648597/1e6 * RR_factor, # per kg of tap water production, conventional treatment
        'NO': 0.007178,
        'SO2':0.012817,
        'NH3': 0.0314
        }
    
    H_Health_dct = {
        'Electricity': 0.040824307,
        'CH4': 34 * HumanHealth_factor,
        'N2O': 298 * HumanHealth_factor,
        'N': -0.637826734*RR_factor,
        'P': -1.774294425*RR_factor,
        'K': -0.116067637*RR_factor,
        'H2O': -1.6207909/1e5 * RR_factor, # per kg of tap water production, conventional treatment
        'NO': 1.133,
        'SO2':1.923,
        'NH3': 1.647
        }
    
    H_Resources_dct = {
        'Electricity': 0.027825633,
        'CH4': 0,  # no GWP to Resource Depletion pathway
        'N2O': 0,  # no GWP to Resource Depletion pathway
        'N': -0.259196888*RR_factor,
        'P': -1.084191599*RR_factor,
        'K': -0.054033438*RR_factor,
        'H2O': -1.5364666/1e5 * RR_factor, # per kg of tap water production, conventional treatment
        'NO': 0,
        'SO2':0,
        'NH3': 0
        }
    
    return price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct

update_resource_recovery_settings()

# =============================================================================
# Load components and systems
# =============================================================================

from . import _components
from ._components import *
_components_loaded = False
def _load_components(reload=False):
    global components, _components_loaded
    if not _components_loaded or reload:
        components = create_components()
        qs.set_thermo(components)
        _components_loaded = True

_impact_item_loaded = False
def _load_lca_data(reload=False):
    '''
    Load impact indicator and impact item data.

    Parameters
    ----------
    reload : bool
        Whether to force reload LCA data.
    '''
    global _impact_item_loaded
    if not _impact_item_loaded or reload:
        indicator_path = os.path.join(g2rt_data_path, 'impact_indicators.csv')
        qs.ImpactIndicator.load_from_file(indicator_path)

        item_path = os.path.join(g2rt_data_path, 'impact_items.xlsx')
        qs.ImpactItem.load_from_file(item_path)
        
        price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct = update_resource_recovery_settings()

        # Impacts associated with streams and electricity
        def create_stream_impact_item(item_ID, dct_key=''):
            dct_key = dct_key or item_ID.rsplit('_item')[0] # `rstrip` will change "struvite_item" to "struv"
            StreamImpactItem(ID=item_ID, GWP=GWP_dct[dct_key],
                             H_Ecosystems=H_Ecosystems_dct[dct_key],
                             H_Health=H_Health_dct[dct_key],
                             H_Resources=H_Resources_dct[dct_key])

        create_stream_impact_item(item_ID='CH4_item')
        create_stream_impact_item(item_ID='N2O_item')
        create_stream_impact_item(item_ID='N_item')
        create_stream_impact_item(item_ID='P_item')
        create_stream_impact_item(item_ID='K_item')
        create_stream_impact_item(item_ID='H2O_item')
        create_stream_impact_item(item_ID='NO_item')
        create_stream_impact_item(item_ID='SO2_item')
        create_stream_impact_item(item_ID='NH3_item')
        ImpactItem(ID='e_item', functional_unit='kWh',
                   GWP=GWP_dct['Electricity'],
                   H_Ecosystems=H_Ecosystems_dct['Electricity'],
                   H_Health=H_Health_dct['Electricity'],
                   H_Resources=H_Resources_dct['Electricity'])

        _impact_item_loaded = True

    return _impact_item_loaded

from . import systems
from .systems import *
from exposan.g2rt.systems import create_system
_system_loaded = False
def _load_system():
    qs.currency = 'USD'
    qs.PowerUtility.price = price_dct['Electricity']
    global sysA, sysB, teaA, teaB, lcaA, lcaB, _system_loaded
    sysA = create_system('A')
    teaA = sysA.TEA
    lcaA = sysA.LCA
    # sysB = create_system('B')
    # teaB = sysB.TEA
    # lcaB = sysB.LCA
    _system_loaded = True


def load():
    if not _components_loaded: _load_components()
    if not _system_loaded: _load_system()
    dct = globals()
    for sys in (sysA, #sysB
                ): dct.update(sys.flowsheet.to_dict())


def __getattr__(name):
    if not _components_loaded or not _system_loaded:
        raise AttributeError(f'module "{__name__}" not yet loaded, '
                             f'load module with `{__name__}.load()`.')

##### Recoveries #####
# Calculate recoveries as in kg N/P/K/H2O per yr
hr_per_yr = 365 * 24
get_N = lambda stream: (stream.imass['NH3']+stream.imass['NonNH3'])*hr_per_yr
get_P = lambda stream: stream.imass['P']*hr_per_yr #TODO: check if mass in and mass out balance in each unit
get_K = lambda stream: stream.imass['K']*hr_per_yr
get_Water = lambda stream: stream.imass['H2O']*hr_per_yr 

# Only for `sysA` or `sysB`
def get_recoveries(system, ppl=default_ppl, include_breakdown=False):
    # AB = system.ID[-1]
    u_reg = system.flowsheet.unit
    dct = globals()
    dct['N_dct'] = N_dct = {}
    dct['P_dct'] = P_dct = {}
    dct['K_dct'] = K_dct = {}
    dct['Water_dct'] = Water_dct = {}
    fp = u_reg.A10
    solids_separator = u_reg.A3
    RO = u_reg.A6
    comp_splitter = u_reg.A18
    mixer = u_reg.Mixer

    # ##### Unique to A or B #####
    # if AB == 'A':
    #     toilet = u_reg.A2
    #     sludge_pasteurization = u_reg.A4
    #     ion_exchange = u_reg.A5
    # else: # unique to sysB
    #     toilet = u_reg.B2
    #     sludge_pasteurization = u_reg.B4
    #     ion_exchange = u_reg.B5

    # ##### Applicable for both A and B #####
    # # In
    # N_dct['urine'] = get_N(toilet.ins[0]) * ppl
    # N_dct['feces'] = get_N(toilet.ins[1]) * ppl
    # N_dct['input'] = N_dct['urine'] + N_dct['feces']
    
    # P_dct['urine'] = get_P(toilet.ins[0]) * ppl
    # P_dct['feces'] = get_P(toilet.ins[1]) * ppl
    # P_dct['input'] = P_dct['urine'] + P_dct['feces']
    
    # K_dct['urine'] = get_K(toilet.ins[0]) * ppl
    # K_dct['feces'] = get_K(toilet.ins[1]) * ppl
    # K_dct['input'] = K_dct['urine'] + K_dct['feces']
    
    # Water_dct['urine'] = get_Water(toilet.ins[0]) * ppl
    # Water_dct['feces'] = get_Water(toilet.ins[1]) * ppl
    # Water_dct['flushing'] = get_Water(mixer.outs[0])
    # # Reverse_osmosis
    # Water_dct['RO_permeate'] = get_Water(RO.outs[0])
    # Water_dct['input'] = (Water_dct['urine'] + Water_dct['feces']
    #                       +Water_dct['flushing'])

    # # drying tunnel produces solid cakes
    # N_dct['dried_solids'] = get_N(dry_sludge_cake.outs[0])
    # P_dct['dried_solids'] = get_P(dry_sludge_cake.outs[0])
    # K_dct['dried_solids'] = get_K(dry_sludge_cake.outs[0])

    # % N, P, and K recovered as a usable fertilizer product,
    # for model metrics and also the Resource Recovery criterion in DMsan analysis
    #TODO: check mass flow to make sure it is right after each unit
    functions = [
            lambda: get_N(comp_splitter.outs[0]) / (get_N(fp.outs[1])+get_N(RO.outs[0])+get_N(RO.outs[1])) * 100,  # total_N_recovery
            lambda: get_P(comp_splitter.outs[1]) / (get_P(fp.outs[1])+get_P(RO.outs[0])+get_P(RO.outs[1])) * 100,  # total_P_recovery
            lambda: get_K(comp_splitter.outs[2]) / (get_K(fp.outs[1])+get_K(RO.outs[0])+get_K(RO.outs[1])) * 100,  # total_K_recovery
            lambda: get_Water(RO.outs[0]) / (get_Water(RO.outs[0])+get_Water(RO.outs[1])+get_Water(fp.outs[1]))* 100 #total water recovery
            ]
    return functions

##### Costs #####
# Learning curve assumptions
percent_CAPEX_to_scale = 0.65
number_of_units = 100000
percent_limit = 0.03 #pessimistic learning curve
learning_curve_percent = 0.95 #pessimistic learning curve

def get_scaled_capital(tea):
    return get_generic_scaled_capital(
        tea=tea,
        percent_CAPEX_to_scale=percent_CAPEX_to_scale,
        number_of_units=number_of_units,
        percent_limit=percent_limit,
        learning_curve_percent=learning_curve_percent
        )

def get_TEA_metrics(system, ppl=default_ppl, include_breakdown=False):
    tea = system.TEA
    get_annual_electricity = lambda system: system.power_utility.cost*system.operating_hours 
    #TODO: check where to set 'operating_hours'
    functions = [lambda: (get_scaled_capital(tea)-tea.net_earnings) / ppl]
    if not include_breakdown: return functions # net cost
    return [
        *functions,
        lambda: get_scaled_capital(tea) / ppl, # CAPEX
        lambda: get_annual_electricity(system)/ppl, # energy (electricity)
        lambda: tea.annual_labor/ppl, # labor
        lambda: (tea.AOC-get_annual_electricity(system)-tea.annual_labor)/ppl, # OPEX (other than energy and labor)
        lambda: tea.sales / ppl, # sales
        ]


def get_normalized_CAPEX(unit, ppl=default_ppl):
    '''Get the CAPEX of a unit normalized to per capita per day.'''
    system = unit.system
    return lambda: system.TEA.get_unit_annualized_equipment_cost(unit)/365/ppl


def get_normalized_electricity_cost(unit, ppl=default_ppl):
    '''Get the energy (electricity) cost of a unit normalized to per capita per day.'''
    return lambda: unit.power_utility.cost /ppl

def get_normalized_OPEX(unit, ppl=default_ppl):
    '''
    Get the OPEX of a unit normalized to per capita per day,
    energy (electricity) and labor cost is not included.
    '''
    # # Wrap single unit in a list for consistent processing
    # if not isinstance(units, Iterable) or isinstance(units, str):
    #     units = [units]
    #   # Calculate total OPEX
    # # Calculate total OPEX by summing up individual values in `add_OPEX`
    # OPEX = sum(sum(u.add_OPEX.values()) for u in units)
    # # Sum the costs of all input streams
    # streams = sum([u.ins for u in units], [])
    # OPEX += sum(s.cost for s in streams)
    
    # # Normalize to per capita per day
    # return lambda: OPEX * 24 / ppl  # convert to per capita per day
    return lambda: unit.OPEX/ ppl # convert to per capita per day

def get_normalized_labor_cost(unit, ppl=default_ppl):
    '''
    Get the labor cost for the maintenance of a unit normalized to per capita per day.
    '''
    return lambda: unit.labor_expense / ppl # convert to per capita per day

def get_normalized_recovery_earning(unit, ppl=default_ppl):
    '''
    Calculate recovery earnings from water and nutrient (K, N, P) recovery
    for a specific unit, normalized to per capita per day.
    '''
    # # Initialize OPEX to accumulate earnings from all streams
    # OPEX = 0
    # Sum up the costs from all input and output streams
    streams = unit.ins + unit.outs
    # Return a callable lambda that computes the normalized recovery earning
    return lambda: sum(s.cost for s in streams) * 24 / ppl # convert to per capita per day

# Define a function to compute the total cost directly
def compute_unit_total_cost(u, ppl):
    system = u.system
    return lambda: (system.TEA.get_unit_annualized_equipment_cost(u)/365/ppl+
                    u.power_utility.cost/ppl+
                    u.OPEX/ppl+
                    u.labor_expense/ppl)

def get_unit_contruction_GW_impact(unit, ppl=default_ppl, time =None, time_unit ='day'):
    system = unit.system
    lca = system.LCA
    return  lambda: lca.get_construction_impacts(unit, annual=True)['GlobalWarming']/ppl
# convert to per capita per year

def get_unit_stream_GW_impact(unit, ppl=default_ppl):
    system = unit.system
    lca = system.LCA
    stream_items = {i for i in unit.ins + unit.outs if i.stream_impact_item}
    # s = lca.get_stream_impacts(stream_items=stream_items, exclude=None,
    #                              kind='all', annual=True)
    return lambda: lca.get_stream_impacts(stream_items=stream_items, exclude=None,
                                 kind='all', annual=True)['GlobalWarming']/ppl
# convert to per capita per year

def get_unit_electrcitiy_GW_impact(unit, ppl=default_ppl, time =None, time_unit ='day'):
    system = unit.system
    lca = system.LCA
    return lambda: lca.get_other_unit_impacts(unit, annual=True)['GlobalWarming']/ppl

# ['GlobalWarming']

def get_LCA_metrics(system, ppl=default_ppl, include_breakdown=False):
    lca = system.LCA
    functions = [
        lambda: lca.total_impacts['GlobalWarming']/lca.lifetime/ppl, # annual GWP
        # ReCiPe LCA functions
        lambda: lca.total_impacts['H_Ecosystems']/lca.lifetime/ppl,
        lambda: lca.total_impacts['H_Health']/lca.lifetime/ppl,
        lambda: lca.total_impacts['H_Resources']/lca.lifetime/ppl,
        ]
    if not include_breakdown: return functions
    return [
        *functions,
        lambda: lca.total_construction_impacts['GlobalWarming']/lca.lifetime/ppl, # construction
        lambda: lca.total_transportation_impacts['GlobalWarming']/lca.lifetime/ppl, # transportation
        lambda: lca.total_stream_impacts['GlobalWarming']/lca.lifetime/ppl, # stream (including fugitive gases and offsets)
        lambda: lca.total_other_impacts['GlobalWarming']/lca.lifetime/ppl, #electricity
        ]

def print_summaries(systems):
    try: iter(systems)
    except: systems = (systems, )

    for sys in systems:
        sys.simulate()

        print(f'\n---------- Summary for {sys.ID} ----------\n')
        if sys.ID in ('sysA', 
                      #'sysB'
                      ):
            recovery_functions = get_recoveries(sys)
            print(f'\nTotal N recovery: {recovery_functions[0]():.1f} %.')
            print(f'\nTotal P recovery: {recovery_functions[1]():.1f} %.')
            print(f'\nTotal K recovery: {recovery_functions[2]():.1f} %.')
            print(f'\nTotal H2O recovery: {recovery_functions[3]():.1f} %.')

            TEA_functions = get_TEA_metrics(sys)
            unit = f'{qs.currency}/cap/yr'
            print(f'\nTotal cost: {TEA_functions[0]():.2f} {unit}.')

            LCA_functions = get_LCA_metrics(sys)
            print(f'\nNet emission: {LCA_functions[0]():.2f} kg CO2-eq/cap/yr.')

            unit = 'points/cap/yr'
            print(f'\nNet ecosystems damage: {LCA_functions[1]():.2f} {unit}.')
            print(f'\nNet health damage: {LCA_functions[2]():.2f} {unit}.')
            print(f'\nNet resources damage: {LCA_functions[3]():.2f} {unit}.')

        else:
            sys.TEA.show()
            print('\n')
            sys.LCA.show()


from . import models
from .models import *

# from . import country_specific
# from .country_specific import *

__all__ = (
    'g2rt_path',
    'g2rt_data_path',
    'results_path',
    *_components.__all__,
    *systems.__all__,
    *models.__all__,
    # *country_specific.__all__,
)