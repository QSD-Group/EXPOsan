#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>
    Tori Morgan <tvlmorgan@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, qsdsan as qs
from math import log
from qsdsan import ImpactItem, StreamImpactItem
from exposan.utils import (
    get_decay_k,
    get_generic_scaled_capital,
    get_generic_tanker_truck_fee as get_tanker_truck_fee,
    )

es_path = os.path.dirname(__file__)
data_path = os.path.join(es_path, 'data')
results_path = os.path.join(es_path, 'results')
# To save simulation data
if not os.path.isdir(results_path): os.mkdir(results_path)


# %%

# =============================================================================
# Unit parameters
# =============================================================================

household_size = 4
household_per_toilet = 5
get_toilet_user = lambda: household_size*household_per_toilet

# Number of people served by the Eco-san
ppl = 300

discount_rate = 0.05

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3

max_CH4_emission = 0.25

emptying_fee = 0.15


# %%

# =============================================================================
# Prices and GWP CFs
# =============================================================================

# Recycled nutrients are sold at a lower price than commercial fertilizers
price_ratio = 1
price_factor = 0.25

operator_daily_wage = 29

price_dct = {
    'Electricity': 0.06,
    'Concrete': 194,
    'Steel': 2.665,
    'N': 1.507*price_factor,
    'P': 3.983*price_factor,
    'K': 1.333*price_factor,
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
    'MgOH2': 1.176277921,
    'struvite': 0,
    'salt': 0.266695553,
    'HCl_acid': 0.8,
    }


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
        indicator_path = os.path.join(data_path, 'impact_indicators.tsv')
        qs.ImpactIndicator.load_from_file(indicator_path)

        item_path = os.path.join(data_path, 'impact_items.xlsx')
        qs.ImpactItem.load_from_file(item_path)

        # Impacts associated with streams and electricity
        def create_stream_impact_item(item_ID, dct_key=''):
            dct_key = dct_key or item_ID.rsplit('_item')[0] # `rstrip` will change "struvite_item" to "struv"
            StreamImpactItem(ID=item_ID, GWP=GWP_dct[dct_key])

        create_stream_impact_item(item_ID='CH4_item')
        create_stream_impact_item(item_ID='N2O_item')
        create_stream_impact_item(item_ID='N_item')
        create_stream_impact_item(item_ID='P_item')
        create_stream_impact_item(item_ID='K_item')
        create_stream_impact_item(item_ID='MgOH2_item')
        create_stream_impact_item(item_ID='struvite_item')
        create_stream_impact_item(item_ID='salt_item')
        create_stream_impact_item(item_ID='HCl_acid_item')
        ImpactItem(ID='e_item', functional_unit='kWh', GWP=GWP_dct['Electricity'])

        _impact_item_loaded = True

    return _impact_item_loaded


from . import systems
from .systems import *
_system_loaded = False
def _load_system():
    qs.currency = 'USD'
    qs.PowerUtility.price = price_dct['Electricity']
    global sysA, sysB, sysC, teaA, teaB, teaC, lcaA, lcaB, lcaC, _system_loaded
    sysA = create_system('A')
    teaA = sysA.TEA
    lcaA = sysA.LCA
    sysB = create_system('B')
    teaB = sysB.TEA
    lcaB = sysB.LCA
    sysC = create_system('C')
    teaC = sysC.TEA
    lcaC = sysC.LCA
    _system_loaded = True


def load():
    if not _components_loaded: _load_components()
    if not _system_loaded: _load_system()
    dct = globals()
    for sys in (sysA, sysB, sysC): dct.update(sys.flowsheet.to_dict())


def __getattr__(name):
    if not _components_loaded or not _system_loaded:
        raise AttributeError(f'module "{__name__}" not yet loaded, '
                             f'load module with `{__name__}.load()`.')


# %%

# =============================================================================
# Util functions
# =============================================================================

# Learning curve assumptions
percent_CAPEX_to_scale = 0.1
number_of_units = 100000
percent_limit = 0.015
learning_curve_percent = 0.925
def get_scaled_capital(tea):
    return get_generic_scaled_capital(
        tea=tea,
        percent_CAPEX_to_scale=percent_CAPEX_to_scale,
        number_of_units=number_of_units,
        percent_limit=percent_limit,
        learning_curve_percent=learning_curve_percent
        )


def get_TEA_metrics(system, include_breakdown=False):
    tea = system.TEA
    get_annual_electricity = lambda system: system.power_utility.cost*system.operating_hours
    functions = [lambda: (tea.EAC-tea.annualized_CAPEX+get_scaled_capital(tea)) / ppl]
    if not include_breakdown: return functions # net cost
    return [
        *functions,
        lambda: get_scaled_capital(system) / ppl, # CAPEX
        lambda: get_annual_electricity()/ppl, # energy (electricity)
        lambda: tea.annual_labor/ppl, # labor
        lambda: (tea.AOC-get_annual_electricity()-tea.annual_labor)/ppl, # OPEX (other than energy and labor)
        lambda: tea.sales / ppl, # sales
        ]


def get_LCA_metrics(system, include_breakdown=False):
    lca = system.LCA
    functions = [
        lambda: lca.total_impacts['GlobalWarming']/lca.lifetime/ppl,
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
        TEA_functions = get_TEA_metrics(sys)
        unit = f'{qs.currency}/cap/yr'
        print(f'\nTotal cost: {TEA_functions[0]():.2f} {unit}.')

        LCA_functions = get_LCA_metrics(sys)
        print(f'\nNet emission: {LCA_functions[0]():.2f} kg CO2-eq/cap/yr.')


#!!! Model for Eco-San not yet updated
# from . import models
# from .models import *


__all__ = (
    'es_path',
    'data_path',
    'results_path',
    *_components.__all__,
    *systems.__all__,
    # *models.__all__,
)