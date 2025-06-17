#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>
    
    Aaron Marszewski <aaronpm3@illinois.edu>
    
    Lewis Rowles <stetsonsc@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, qsdsan as qs
from qsdsan import ImpactItem, StreamImpactItem
import numpy as np
from exposan.utils import (
    _init_modules,
    get_generic_scaled_capital,
    )







br_path = os.path.dirname(__file__)
module = os.path.split(br_path)[-1]
data_path, results_path = _init_modules(module, include_data_path=True)



# %%
# =============================================================================
# Metrics for model output 
# =============================================================================

def cost_per_ton_biochar(model):
    biochar = model.system.flowsheet.stream.biochar
    biochar_mass = biochar.F_mass * 24 * 365 / 1000  # tons/year
   
    if biochar_mass == 0:
        return np.nan
    total_cost = model.system.TEA.AOC + get_scaled_capital(model.system.TEA) + model.system.TEA.annual_labor - model.system.TEA.sales
    # print("AOC:", model.system.TEA.AOC)
    # print("scaled capital:", get_scaled_capital(model.system.TEA))
    # print("labor:", model.system.TEA.annual_labor)    
    # print("sales:", model.system.TEA.sales)
    return total_cost / biochar_mass

def gwp_per_ton_biosolids(model):
    LCA = model.system.LCA
    annual_LCA = (
    LCA.total_construction_impacts['GlobalWarming'] +
    LCA.total_transportation_impacts['GlobalWarming'] +
    LCA.total_stream_impacts['GlobalWarming']
) / LCA.lifetime #kgCO2eq/yr
    annual_C_seq = model.system.flowsheet.unit.A4.biochar_sequesterable_carbon * 1000 * 44/12 #kgCO2eq seq C/yr

    return (annual_LCA - annual_C_seq) / (model.system.flowsheet.unit.A1.biosolids_generated) # kg CO2-eq/ton-biochar

def biochar_generated(model):
    system = model.system
    return system.flowsheet.stream.biochar.F_mass * 24 * 365 / 1000  # tons biochar/year

def sequesterable_carbon(model):
    system = model.system

    return (
        0 if system.flowsheet.stream.biochar.F_mass == 0 else
        (system.flowsheet.unit.A4.biochar_sequesterable_carbon)                       #tons sequesterable C / year
    )

def drying_requirement(model):
    system = model.system
    return system.flowsheet.unit.A8.Q_drying_norm/1000  # GJ/ton-biosolids
def drying_cost(model):
    system = model.system
    biochar = model.system.flowsheet.stream.biochar
    biochar_mass = biochar.F_mass * 24 * 365 / 1000  # tons/year
    MJ_per_year = system.flowsheet.unit.A8.Q_drying/1000*24*365  #MJ/yr
    price_per_MJ = .045 #USD/MJ
    return MJ_per_year * price_per_MJ / biochar_mass #USD/ton biochar
def energy_conversion_efficiency(model):
    system = model.system
    return system.flowsheet.unit.A4.ECE             #ECE %
def char_yield(model):
    system = model.system
    return system.flowsheet.unit.A4.char_yield_db_percent
def ash_content_feedstock(model):
    system = model.system
    return system.flowsheet.unit.A4.ACf
def moisture_content_feedstock(model):
    system = model.system
    return system.flowsheet.unit.A8.MC_feedstock
def biosolids_HHV(model):
    system = model.system
    return system.flowsheet.unit.A4.HHV           #MJ/kg
# %%

# =============================================================================
# Unit parameters
# =============================================================================



max_CH4_emission = 0.25


truck_fee = 6.21 # USD/m3




# =============================================================================
# Prices and GWP CFs
# =============================================================================
# Should be changed based on location
discount_rate = .05
price_ratio = 1
operator_daily_wage = 29
const_daily_wage = 17
const_person_days = 100
INCLUDED_RESOURCE_RECOVERY = True # including resource recovery

def update_resource_recovery_settings():
    global INCLUDE_RESOURCE_RECOVERY
    global price_dct, GWP_dct
    RR_factor = int(bool(INCLUDED_RESOURCE_RECOVERY)) # RR_factor = 1 if resource recovery is included
    
    price_dct = {
        'Electricity': 0.13,
        'Concrete': 194*price_ratio,
        'Steel': 2.665*price_ratio,
        'biochar': 0,  # 0.014*price_ratio,  # assuming value of biochar is 0 for TEA - HACL
        }

    GWP_dct = {
        'Electricity': 0.69,
        'CH4': 34,
        'N2O': 298,
        # Assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
        'biochar': 0 #-0.2*0.9*(44/12), included manually, so overridden to zero for stream
        }

        
    
    return price_dct, GWP_dct

update_resource_recovery_settings()


# %%

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
        indicator_path = os.path.join(data_path, 'impact_indicators.csv')
        qs.ImpactIndicator.load_from_file(indicator_path)

        item_path = os.path.join(data_path, 'impact_items.xlsx')
        qs.ImpactItem.load_from_file(item_path)
        
        price_dct, GWP_dct = update_resource_recovery_settings()

        # Impacts associated with streams and electricity
        def create_stream_impact_item(item_ID, dct_key=''):
            dct_key = dct_key or item_ID.rsplit('_item')[0] # `rstrip` will change "struvite_item" to "struv"
            StreamImpactItem(ID=item_ID,
                             GWP=GWP_dct[dct_key],
                            )

        create_stream_impact_item(item_ID='CH4_item')
        create_stream_impact_item(item_ID='N2O_item')
        create_stream_impact_item(item_ID='biochar_item')
        ImpactItem(ID='e_item', functional_unit='kWh',
                   GWP=GWP_dct['Electricity'])

        # Update prices
        ImpactItem.get_item('Concrete').price = price_dct['Concrete']
        ImpactItem.get_item('Steel').price = price_dct['Steel']

        _impact_item_loaded = True

    return _impact_item_loaded


from . import systems
from .systems import *
_system_loaded = False
def _load_system():
    qs.currency = 'USD'
    qs.PowerUtility.price = price_dct['Electricity']
    global sysA, teaA, lcaA, _system_loaded
    sysA = create_system('A')
    teaA = sysA.TEA
    lcaA = sysA.LCA
    _system_loaded = True


def load():
    if not _components_loaded: _load_components()
    if not _system_loaded: _load_system()
    dct = globals()
    dct.update(sysA.flowsheet.to_dict())
    # for sys in (sysA): dct.update(sys.flowsheet.to_dict())


def __getattr__(name):
    if not _components_loaded or not _system_loaded:
        raise AttributeError(
            f'Module {__name__} does not have the attribute "{name}" '
            'and the module has not been loaded, '
            f'loading the module with `{__name__}.load()` may solve the issue.')


# %%

# =============================================================================
# Util functions
# =============================================================================

##### Recoveries #####
# Calculate recoveries as in kg C/N/P/K per yr
hr_per_yr = 365 * 24
def get_C(stream):
    sys = stream.source.system
    for unit in sys.units:
        if hasattr(unit, 'carbon_COD_ratio'): carbon_COD_ratio = unit.carbon_COD_ratio
    return lambda: stream.COD*stream.F_vol/1e3*carbon_COD_ratio*hr_per_yr
get_C_gas = lambda stream: stream.imol['CH4']*12*hr_per_yr
get_N_gas = lambda stream: stream.imol['N2O']*28*hr_per_yr


# Only for `sysA`
def get_sustainability_indicators(system, include_breakdown=False):
    A = system.ID[-1]
    if A not in ('A'):
        raise ValueError('This function is only applicable for `sysA`, '
                         f'not `{system.ID}`.')

    u_reg = system.flowsheet.unit
    dct = globals()


##### Costs #####
# Learning curve assumptions
percent_CAPEX_to_scale = 0.1
number_of_units = 10000
percent_limit = 0.01
learning_curve_percent = 0.9
def get_scaled_capital(tea):
    if tea.system.ID[-1] == 'D':
        new_CAPEX_annualized = tea.annualized_CAPEX
    else:
        new_CAPEX_annualized = get_generic_scaled_capital(
            tea=tea,
            percent_CAPEX_to_scale=percent_CAPEX_to_scale,
            number_of_units=number_of_units,
            percent_limit=percent_limit,
            learning_curve_percent=learning_curve_percent
        )
    return new_CAPEX_annualized


def get_TEA_metrics(system, include_breakdown=False):
    tea = system.TEA
    ppl = get_ppl(system.ID)
    get_annual_electricity = lambda system: system.power_utility.cost*system.operating_hours
    functions = [lambda: (get_scaled_capital(tea)-tea.net_earnings) / ppl]
    if not include_breakdown: return functions # net cost
    return [
        *functions,
        lambda: get_scaled_capital(system) / ppl, # CAPEX
        lambda: get_annual_electricity()/ppl, # energy (electricity)
        lambda: tea.annual_labor/ppl, # labor
        lambda: (tea.AOC-get_annual_electricity()-tea.annual_labor)/ppl, # OPEX (other than energy and labor)
        lambda: tea.sales / ppl, # sales
        ]

def get_ppl(kind):
    if kind.lower()=='10k' or kind.upper()[-1]=='C': return 10000
    return 12000

def get_normalized_CAPEX(units):
    '''Get the CAPEX of a unit/units normalized to per capita per day.'''
    system = units[0].system
    return system.TEA.get_unit_annualized_equipment_cost(units)/365/get_ppl(system.ID)


def get_noramlized_electricity_cost(units):
    '''Get the energy (electricity) cost of a unit/units normalized to per capita per day.'''
    ppl = get_ppl(units[0].ID[0])
    return sum(u.power_utility.cost for u in units)/ppl


def get_normalized_OPEX(units):
    '''
    Get the OPEX of a unit/units normalized to per capita per day,
    energy (electricity) cost is not included.
    '''
    OPEX = sum(u.add_OPEX.values() for u in units)
    streams = sum([u.ins for u in units], [])
    OPEX += sum(s.cost for s in streams)
    ppl = get_ppl(units[0].ID[0])
    return OPEX * 24 / ppl # convert to per capita per day


##### Emissions #####
def get_LCA_metrics(system, include_breakdown=False):
    lca = system.LCA
    ppl = get_ppl(system.ID)
    functions = [
        lambda: lca.total_impacts['GlobalWarming']/lca.lifetime/ppl, # annual GWP

        ]
    if not include_breakdown: return functions
    return [
        *functions,
        lambda: lca.total_construction_impacts['GlobalWarming']/lca.lifetime/ppl, # construction
        lambda: lca.total_transportation_impacts['GlobalWarming']/lca.lifetime/ppl, # transportation
        lambda: lca.total_stream_impacts['GlobalWarming']/lca.lifetime/ppl, # stream (including fugitive gases and offsets)
        lambda: lca.total_other_impacts['GlobalWarming']/lca.lifetime/ppl,
        ]


def print_summaries(systems):
    try: iter(systems)
    except: systems = (systems, )

    for sys in systems:
        sys.simulate()

        print(f'\n---------- Summary for {sys.ID} ----------\n')
        if sys.ID in ('sysA'):

            TEA_functions = get_TEA_metrics(sys)
            unit = f'{qs.currency}/cap/yr'
            print(f'\nTotal cost: {TEA_functions[0]():.2f} {unit}.')

            LCA_functions = get_LCA_metrics(sys)
            print(f'\nNet emission: {LCA_functions[0]():.2f} kg CO2-eq/cap/yr.')



        else:
            sys.TEA.show()
            print('\n')
            sys.LCA.show()



from . import models
from .models import *


__all__ = (
	'br_path',
	'data_path',
	'results_path',
	*_components.__all__,
	*systems.__all__,
    *models.__all__,
 	)