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

import os, qsdsan as qs
from qsdsan import ImpactIndicator, ImpactItem, StreamImpactItem
from exposan.utils import _init_modules
pou_path = os.path.dirname(__file__)
module = os.path.split(pou_path)[-1]
data_path, results_path = _init_modules(module, include_data_path=True)


# %%

# =============================================================================
# Unit parameters
# =============================================================================

discount_rate = 0.05
start_year = 2018
lifetime = 5

water_source = 'GW'
household_size = 6
household_per_container = 1
ppl = 1000 # 1k or 500


# %%

# =============================================================================
# Prices and GWP CFs
# =============================================================================

#!!! Might need updating
price_dct = {
    'Electricity': 0.17,
    'Concrete': 194,
    'Steel': 2.665,
    'NaClO': 1.96/0.15/1.21/0.125,
    'Polyethylene': 0,
    }

GWP_dct = {
    'Electricity': 0.1135,
    'CH4': 28,
    'N2O': 265,
    'N': -5.4,
    'NaClO': 2.6287, 
    'Polyethylene': 2.7933, 
    }


# %%

# =============================================================================
# Load components and system
# =============================================================================

from . import _components
from ._components import create_components
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
        ImpactIndicator('GWP', unit='kg CO2') # global warming potential

        item_path = os.path.join(data_path, 'impact_items.xlsx')
        qs.ImpactItem.load_from_file(item_path)

        # Impacts associated with streams and electricity
        def create_stream_impact_item(item_ID):
            StreamImpactItem(ID=item_ID,
                             GWP=GWP_dct[item_ID.rsplit('_item')[0]])

        create_stream_impact_item(item_ID='NaClO_item')
        create_stream_impact_item(item_ID='Polyethylene_item')
        ImpactItem(ID='E_item', functional_unit='kWh', GWP=GWP_dct['Electricity'])

        _impact_item_loaded = True

    return _impact_item_loaded


from . import systems
from .systems import *
_system_loaded = False
def _load_system():
    qs.currency = 'USD'
    qs.PowerUtility.price = price_dct['Electricity']
    global sysA, sysB, sysC, sysD, teaA, teaB, teaC, teaD, lcaA, lcaB, lcaC, lcaD, _system_loaded
    sysA = create_system('A')
    teaA = sysA.TEA
    lcaA = sysA.LCA
    sysB = create_system('B')
    teaB = sysB.TEA
    lcaB = sysB.LCA
    sysC = create_system('C')
    teaC = sysC.TEA
    lcaC = sysC.LCA
    sysD = create_system('D')
    teaD = sysD.TEA
    lcaD = sysD.LCA
    _system_loaded = True


def load():
    if not _components_loaded: _load_components()
    if not _system_loaded: _load_system()
    dct = globals()
    for sys in (sysA, sysB, sysC, sysD): dct.update(sys.flowsheet.to_dict())


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

def update_water_source(system, water_source='GW'):
    setattr(system.flowsheet.unit[0], 'water_source', water_source)
    
def update_number_of_householdsize(system, household_size=household_size, ppl=ppl):
    u0, u1 = system.path
    number_of_households = ppl / household_size
    setattr(u0, 'household_size', household_size)
    setattr(u0, 'number_of_households', number_of_households)
    setattr(u1, 'number_of_households', number_of_households)


def get_TEA_metrics(system, ppl=ppl, include_breakdown=False):
    tea = system.TEA
    functions = [lambda: tea.EAC/ppl] # annual cost
    
    if not include_breakdown: return functions
    
    get_AOC = lambda: tea.AOC/ppl # OPEX
    get_annual_electricity = lambda: system.power_utility.cost*system.operating_hours
    functions.extend([
        lambda: tea.annualized_CAPEX/ppl,
        get_AOC,
        get_annual_electricity,
        lambda: get_AOC()-get_annual_electricity(), # OPEX other than energy
        ])

    return functions

def get_LCA_metrics(system, ppl=ppl, include_breakdown=False):
    lca = system.LCA
    factor = lca.lifetime * ppl
    functions = [lambda: lca.total_impacts[ind.ID]/factor
                 for ind in lca.indicators]
    
    if not include_breakdown: return functions
    
    for ind in lca.indicators:
        ID = ind.ID
        functions.extend([
            lambda: lca.total_construction_impacts[ID]/factor,
            lambda: lca.total_transportation_impacts[ID]/factor,
            lambda: lca.total_stream_impacts[ID]/factor, # including fugitive gases and offsets
            lambda: lca.total_other_impacts[ID]/factor, # e.g., electricity
            ])

    return functions


def print_summaries(systems, include_breakdown=False):
    try: iter(systems)
    except: systems = (systems, )

    for system in systems:
        system.simulate()

        print(f'\n---------- Summary for {system.ID} ----------\n')

        TEA_functions = get_TEA_metrics(system, include_breakdown=include_breakdown)
        LCA_functions = get_LCA_metrics(system, include_breakdown=include_breakdown)
        indicators = system.LCA.indicators
        
        unit = f'{qs.currency}/cap/yr'
        print(f'\nTotal cost: {TEA_functions[0]():.2f} {unit}.')
        
        if include_breakdown:
            print(f'\nCAPEX: {TEA_functions[1]():.2f} {unit}.')
            print(f'\nOPEX: {TEA_functions[2]():.2f} {unit}.')
            print(f'\nElectricity: {TEA_functions[3]():.2f} {unit}.')
            print(f'\nOthers: {TEA_functions[4]():.2f} {unit}.') 
        else:
            for n, ind in enumerate(indicators):
                ID, unit = ind.ID, ind.unit+'/cap/year'
                print(f'\nTotal {ID}: {LCA_functions[n-1]():.2f} {unit}.')

        n = len(indicators)
        if include_breakdown:
            for ind in indicators:
                ID, unit = ind.ID, ind.unit+'/cap/year'
                print(f'\nFor indicator {ID}:')
                print(f'\nConstruction: {LCA_functions[n]():.2f} {unit}.')
                print(f'\nTransportation: {LCA_functions[n+1]():.2f} {unit}.')
                print(f'\nStreams: {LCA_functions[n+2]():.2f} {unit}.')
                print(f'\nOthers: {LCA_functions[n+3]():.2f} {unit}.')
                n += 4


from . import models
from .models import *

__all__ = (
    'pou_path',
    'data_path',
    'results_path',
    *_components.__all__,
    *systems.__all__,
    *models.__all__,
)