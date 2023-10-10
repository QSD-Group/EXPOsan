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

GW_water = {
    'Ecoli': 200000, # CFU/mg
    'turbidity': 5, # NTU
    'TOC': 5, # mg/L
    'Ca': 30, # mg/L
    'Mg': 30, # mg/L
    'UVT': 0.8,
    }

SW_water = {
    'Ecoli': 200000, # CFU/mg
    'turbidity': 20, # NTU
    'TOC': 10, # mg/L
    'Ca': 10, # mg/L
    'Mg': 10, # mg/L
    'UVT': 0.8,
    }


household_size = 6
household_per_container = 1
def get_pou_user(household_size=household_size, household_per_container=household_per_container):
    return household_size * household_per_container

ppl = 1000 # 1k or 500
def get_number_of_households(household_size=household_size, ppl=ppl):
    return ppl / household_size


# %%

# =============================================================================
# Prices and GWP CFs
# =============================================================================

#!!! Definitely not complete, need updating
price_dct = {
    'Electricity': 0.17,
    'NaClO': 1.96,
    'Polyethylene': 0,
    }

GWP_dct = {
    'Electricity': 0.1135,
    'NaClO': 1.0,
    'Polyethylene': 1.0,
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
        def create_stream_impact_item(item_ID, dct_key=''):
            StreamImpactItem(ID=item_ID, GWP=GWP_dct[dct_key])

        create_stream_impact_item(item_ID='NaClO_item')
        create_stream_impact_item(item_ID='Polyethylene_item')
        ImpactItem(ID='e_item', functional_unit='kWh', GWP=GWP_dct['Electricity'])

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

#!!! Need to add summarizing functions




#!!! Model not yet updated
# from . import models
# from .models import *


__all__ = (
    'pou_path',
    'data_path',
    'results_path',
    *_components.__all__,
    *systems.__all__,
    # *models.__all__,
)