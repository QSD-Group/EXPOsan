#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np, qsdsan as qs
from qsdsan import ImpactItem, StreamImpactItem
from exposan.utils import get_decay_k, get_generic_scaled_capital

re_path = os.path.dirname(__file__)
data_path = os.path.join(re_path, 'data')
results_path = os.path.join(re_path, 'results')
# To save simulation data
if not os.path.isdir(results_path): os.mkdir(results_path)


# %%

# =============================================================================
# Unit parameters
# =============================================================================

# Number of people served by the Reclaimer
#   baseline: 30 users
#   scaled: 120 users
ppl = 120  # total population served by system

discount_rate = 0.05

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3

max_CH4_emission = 0.25


# =============================================================================
# Prices and GWP CFs
# =============================================================================

# Recycled nutrients are sold at a lower price than commercial fertilizers
price_factor = 0.25

# Should be changed based on country
price_ratio = 1

price_dct = {
    'Electricity': 0.06,
    'wages': 3.64,
    'N': 1.507*price_factor,
    'P': 3.983*price_factor,
    'K': 1.333*price_factor,
    'conc_NH3': 1.333*(14/17)*price_factor,
    'MgOH2': 0.0*price_ratio,
    'KCl': 0.0*price_ratio,
    'GAC': 0.0*price_ratio,
    'zeolite': 0.0*price_ratio,
    }

GWP_dct = {
    'Electricity': 0.69,
    'CH4': 34,
    'N2O': 298,
    'N': -5.4,
    'P': -4.9,
    'K': -1.5,
    'conc_NH3': -5.4*(14/17),
    'MgOH2': 1.176277921,
    'KCl': 0.8,
    'GAC': 8.388648277,
    'zeolite': 5.175,
    }

EcosystemQuality_factor = 29320 * (2.8e-09+7.65e-14)  # (pt/species.yr) * (species.yr/kgCO2eq)
H_Ecosystems_dct = {
    'Electricity': 0.002456338,
    'CH4': 34 * EcosystemQuality_factor,
    'N2O': 298 * EcosystemQuality_factor,
    'N': -0.0461961,
    'P': -0.093269908,
    'K': -0.01895794,
    'conc_NH3': -0.0461961*(14/17),
    'MgOH2': 0.209556136,
    'KCl': 0.00220372,
    'GAC': 0.000437915,
    'zeolite': 0.020161669,
    }

HumanHealth_factor = 436000 * 9.28e-07  # (pt/DALY) * (DALY/kgCO2eq)
H_Health_dct = {
    'Electricity': 0.040824307,
    'CH4': 34 * HumanHealth_factor,
    'N2O': 298 * HumanHealth_factor,
    'N': -0.637826734,
    'P': -1.774294425,
    'K': -0.116067637,
    'conc_NH3': -0.637826734*(14/17),
    'MgOH2': 4.639146841,
    'KCl': 0.036770628,
    'GAC': 0.003448424,
    'zeolite': 0.36462721,
    }

H_Resources_dct = {
    'Electricity': 0.027825633,
    'CH4': 0,  # no GWP to Resource Depletion pathway
    'N2O': 0,  # no GWP to Resource Depletion pathway
    'N': -0.259196888,
    'P': -1.084191599,
    'K': -0.054033438,
    'conc_NH3': -0.259196888*(14/17),
    'MgOH2': 4.05197164,
    'KCl': 0.031653596,
    'GAC': 0.002986373,
    'zeolite': 0.224590444,
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
        indicator_path = os.path.join(data_path, 'impact_indicators.csv')
        qs.ImpactIndicator.load_from_file(indicator_path)

        item_path = os.path.join(data_path, 'impact_items.xlsx')
        qs.ImpactItem.load_from_file(item_path)

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
        create_stream_impact_item(item_ID='MgOH2_item')
        create_stream_impact_item(item_ID='KCl_item')
        create_stream_impact_item(item_ID='GAC_item')
        create_stream_impact_item(item_ID='zeolite_item')
        create_stream_impact_item(item_ID='conc_NH3_item')
        ImpactItem(ID='e_item', functional_unit='kWh',
                   GWP=GWP_dct['Electricity'],
                   H_Ecosystems=H_Ecosystems_dct['Electricity'],
                   H_Health=H_Health_dct['Electricity'],
                   H_Resources=H_Resources_dct['Electricity'])

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
        raise AttributeError(f'module "{__name__}" not yet loaded, '
                             f'load module with `{__name__}.load()`.')


# %%

# =============================================================================
# Util functions
# =============================================================================

##### Recoveries #####
# Calculate recoveries as in kg N/P/K per yr
hr_per_yr = 365 * 24
get_N = lambda stream: stream.TN*stream.F_vol/1e3*hr_per_yr
get_P = lambda stream: stream.TP*stream.F_vol/1e3*hr_per_yr
get_K = lambda stream: stream.TK*stream.F_vol/1e3*hr_per_yr

# Only for `sysB` or `sysC`
def get_recoveries(system, include_breakdown=False):
    BC = system.ID[-1]
    if BC not in ('B', 'C'):
        if BC in ('A', 'D'): return [lambda: 0]*3 # recoveries all 0 for sysA and sysD
        raise ValueError('This function is only applicable for `sysA`, `sysB`, `sysC`, or `sysD`, '
                         f'not `{system.ID}`.')

    u_reg = system.flowsheet.unit
    dct = globals()
    dct['N_dct'] = N_dct = {}
    dct['P_dct'] = P_dct = {}
    dct['K_dct'] = K_dct = {}

    ##### Unique to B or C #####
    if BC == 'B':
        toilet = u_reg.B2
        sludge_pasteurization = u_reg.B4
        ion_exchange = u_reg.B6
    else: # unique to sysC
        toilet = u_reg.C2
        sludge_pasteurization = u_reg.C4
        ion_exchange = u_reg.C6

    ##### Applicable for both B and C #####
    # In
    N_dct['urine'] = get_N(toilet.ins[0]) * ppl
    N_dct['feces'] = get_N(toilet.ins[1]) * ppl
    N_dct['input'] = N_dct['urine'] + N_dct['feces']
    P_dct['urine'] = get_P(toilet.ins[0]) * ppl
    P_dct['feces'] = get_P(toilet.ins[1]) * ppl
    P_dct['input'] = P_dct['urine'] + P_dct['feces']
    K_dct['urine'] = get_K(toilet.ins[0]) * ppl
    K_dct['feces'] = get_K(toilet.ins[1]) * ppl
    K_dct['input'] = K_dct['urine'] + K_dct['feces']

    # Sludge Pasteurization
    N_dct['treated_sludge'] = get_N(sludge_pasteurization.outs[-1])
    P_dct['treated_sludge'] = get_P(sludge_pasteurization.outs[-1])
    K_dct['treated_sludge'] = get_K(sludge_pasteurization.outs[-1])

    # Ion Exchange
    N_dct['NH3'] = ion_exchange.outs[3].imol['NH3'] * 14 * hr_per_yr

    # % N, P, and K recovered as a usable fertilizer product,
    # for model metrics and also the Resource Recovery criterion in DMsan analysis
    functions = [
            lambda: (N_dct['treated_sludge']+N_dct['NH3']) / N_dct['input'] * 100,  # total_N_recovery
            lambda: P_dct['treated_sludge'] / P_dct['input'] * 100,  # total_P_recovery
            lambda: K_dct['treated_sludge'] / K_dct['input'] * 100,  # total_K_recovery
            ]
    if not include_breakdown: return functions

    return [
        *functions,
        lambda: N_dct['treated_sludge'] / N_dct['input'] * 100,  # N_sludge_pasteurization
        lambda: N_dct['NH3'] / N_dct['input'] * 100,  # N_ion_exchange
        lambda: P_dct['treated_sludge'] / P_dct['input'] * 100,  # P_sludge_pasteurization
        lambda: K_dct['treated_sludge'] / K_dct['input'] * 100,  # K_sludge_pasteurization
        ]


##### Costs #####
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
        lambda: get_scaled_capital(tea) / ppl, # CAPEX
        lambda: get_annual_electricity()/ppl, # energy (electricity)
        lambda: tea.annual_labor/ppl, # labor
        lambda: (tea.AOC-get_annual_electricity()-tea.annual_labor)/ppl, # OPEX (other than energy and labor)
        lambda: tea.sales / ppl, # sales
        ]


def get_normalized_CAPEX(units):
    '''Get the CAPEX of a unit/units normalized to per capita per day.'''
    system = units[0].system
    return system.TEA.get_unit_annualized_equipment_cost(units)/365/ppl


def get_noramlized_electricity_cost(units):
    '''Get the energy (electricity) cost of a unit/units normalized to per capita per day.'''
    return sum(u.power_utility.cost for u in units)/ppl


def get_normalized_OPEX(units):
    '''
    Get the OPEX of a unit/units normalized to per capita per day,
    energy (electricity) cost is not included.
    '''
    OPEX = sum(u.add_OPEX.values() for u in units)
    streams = sum([u.ins for u in units], [])
    OPEX += sum(s.cost for s in streams)
    return OPEX * 24 / ppl # convert to per capita per day


##### Emissions #####
def get_LCA_metrics(system, include_breakdown=False):
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
        lambda: lca.total_other_impacts['GlobalWarming']/lca.lifetime/ppl,
        ]


def print_summaries(systems):
    try: iter(systems)
    except: systems = (systems, )

    for sys in systems:
        sys.simulate()

        print(f'\n---------- Summary for {sys.ID} ----------\n')
        if sys.ID in ('sysB', 'sysC'):
            recovery_functions = get_recoveries(sys)
            print(f'\nTotal N recovery: {recovery_functions[0]():.1f} %.')
            print(f'\nTotal P recovery: {recovery_functions[1]():.1f} %.')
            print(f'\nTotal K recovery: {recovery_functions[2]():.1f} %.')

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

from . import country_specific
from .country_specific import *

__all__ = (
    're_path',
    'data_path',
    'results_path',
    *_components.__all__,
    *systems.__all__,
    *models.__all__,
    *country_specific.__all__,
)