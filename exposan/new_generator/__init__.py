#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Hannah Lohman <hlohman94@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, qsdsan as qs
from qsdsan import ImpactItem, StreamImpactItem
from exposan.utils import (
    get_decay_k,
    get_generic_scaled_capital,
    get_generic_tanker_truck_fee as get_tanker_truck_fee,
    )

# Module-wise setting on whether to allow resource recovery
INCLUDE_RESOURCE_RECOVERY = False

# Check if has access to the private repository
from qsdsan.utils import data_path as qs_data_path
if not os.path.isdir(os.path.join(qs_data_path, 'sanunit_data/ng')):
    raise ModuleNotFoundError(
        'Files associated with the NEWgenerator system (under non-disclosure agreement) '
        'cannot be found, '
        'please set path to use the QSDsan-private repository if you have access.'
        )

ng_path = os.path.dirname(__file__)
data_path = os.path.join(ng_path, 'data')
results_path = os.path.join(ng_path, 'results')
# To save simulation data
if not os.path.isdir(results_path): os.mkdir(results_path)


# %%

# =============================================================================
# Unit parameters
# =============================================================================

# Number of people served by the one NEWgenerator100
default_ppl = 100

discount_rate = 0.05

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3

max_CH4_emission = 0.25

emptying_fee = 0.15

# Nutrient loss during application
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05

# Energetic content of the biogas
biogas_energy = 803  # kJ/mol CH4
LPG_energy = 50  # MJ/kg
get_biogas_factor = lambda: biogas_energy/cmps.CH4.MW/LPG_energy


# =============================================================================
# Prices and GWP CFs
# =============================================================================

# Recycled nutrients are sold at a lower price than commercial fertilizers
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
        'conc_NH3': 1.333*(14/17)*price_factor*RR_factor,
        'biogas': 0,
        'GAC': 1.1*price_ratio,
        'zeolite': 0.23*price_ratio,
        'NaCl': 0.2758*price_ratio,
        'NaCl1': 0.2758*price_ratio,
        'NaOH': 0.16*price_ratio,
        'LPG': 1.3916,  # USD/kg, 0.710 USD/L global average, 0.51 kg = 1 L
        'wages': 3.64,
        }
    
    GWP_dct = {
        'Electricity': 0.69,
        'CH4': 34,
        'N2O': 298,
        # # Below are for comparison with old NEWgen codes
        # 'CH4': 28,
        # 'N2O': 265,
        'N': -5.4*RR_factor,
        'P': -4.9*RR_factor,
        'K': -1.5*RR_factor,
        'conc_NH3': -5.4*(14/17)*RR_factor,
        'biogas': 0,
        'GAC': 8.388648277,
        'zeolite': 5.175,
        'NaCl': 0.266695553,
        'NaOH': 1.228848922,
        'LPG': 0.714219323022122,  # 2.93, 3.05 Bwaise
        }
    
    H_Ecosystems_dct = {
        'Electricity': 0.002456338,
        'CH4': 34 * EcosystemQuality_factor,
        'N2O': 298 * EcosystemQuality_factor,
        'N': -0.0461961*RR_factor,
        'P': -0.093269908*RR_factor,
        'K': -0.01895794*RR_factor,
        'conc_NH3': -0.0461961*(14/17)*RR_factor,
        'biogas': 0,
        'GAC': 0.000437915,
        'zeolite': 0.020161669,
        'NaCl': 0.001166777,
        'NaOH': 0.005304792,
        'LPG': 0.003610414,
        }
    
    H_Health_dct = {
        'Electricity': 0.040824307,
        'CH4': 34 * HumanHealth_factor,
        'N2O': 298 * HumanHealth_factor,
        'N': -0.637826734*RR_factor,
        'P': -1.774294425*RR_factor,
        'K': -0.116067637*RR_factor,
        'conc_NH3': -0.637826734*(14/17)*RR_factor,
        'biogas': 0,
        'GAC': 0.003448424,
        'zeolite': 0.36462721,
        'NaCl': 0.020900452,
        'NaOH': 0.092927677,
        'LPG': 0.044654087,
        }
    
    H_Resources_dct = {
        'Electricity': 0.027825633,
        'CH4': 0,  # no GWP to Resource Depletion pathway
        'N2O': 0,  # no GWP to Resource Depletion pathway
        'N': -0.259196888*RR_factor,
        'P': -1.084191599*RR_factor,
        'K': -0.054033438*RR_factor,
        'conc_NH3': -0.259196888*(14/17)*RR_factor,
        'biogas': 0,
        'GAC': 0.002986373,
        'zeolite': 0.224590444,
        'NaCl': 0.013404934,
        'NaOH': 0.052849392,
        'LPG': 0.178274445,
        }
    
    return price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct

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
        create_stream_impact_item(item_ID='conc_NH3_item')
        create_stream_impact_item(item_ID='biogas_item')
        create_stream_impact_item(item_ID='GAC_item')
        create_stream_impact_item(item_ID='zeolite_item')
        create_stream_impact_item(item_ID='NaCl_item')
        create_stream_impact_item(item_ID='NaCl1_item', dct_key='NaCl')
        create_stream_impact_item(item_ID='NaOH_item')
        create_stream_impact_item(item_ID='LPG_item')
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
    global sysA, sysB, teaA, teaB, lcaA, lcaB, _system_loaded
    sysA = create_system('A')
    teaA = sysA.TEA
    lcaA = sysA.LCA
    sysB = create_system('B')
    teaB = sysB.TEA
    lcaB = sysB.LCA
    _system_loaded = True


def load():
    if not _components_loaded: _load_components()
    if not _system_loaded: _load_system()
    dct = globals()
    for sys in (sysA, sysB): dct.update(sys.flowsheet.to_dict())


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

# Only for `sysA` or `sysB`
def get_recoveries(system, ppl=default_ppl, include_breakdown=False):
    AB = system.ID[-1]
    u_reg = system.flowsheet.unit
    dct = globals()
    dct['N_dct'] = N_dct = {}
    dct['P_dct'] = P_dct = {}
    dct['K_dct'] = K_dct = {}

    ##### Unique to A or B #####
    if AB == 'A':
        toilet = u_reg.A2
        sludge_pasteurization = u_reg.A4
        ion_exchange = u_reg.A5
    else: # unique to sysB
        toilet = u_reg.B2
        sludge_pasteurization = u_reg.B4
        ion_exchange = u_reg.B5

    ##### Applicable for both A and B #####
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
    N_dct['treated_sludge'] = get_N(sludge_pasteurization.outs[3])
    P_dct['treated_sludge'] = get_P(sludge_pasteurization.outs[3])
    K_dct['treated_sludge'] = get_K(sludge_pasteurization.outs[3])

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
percent_CAPEX_to_scale = 0.65
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


def get_TEA_metrics(system, ppl=default_ppl, include_breakdown=False):
    tea = system.TEA
    get_annual_electricity = lambda system: system.power_utility.cost*system.operating_hours
    functions = [lambda: (get_scaled_capital(tea)-tea.net_earnings) / ppl]
    if not include_breakdown: return functions # net cost
    return [
        *functions,
        lambda: get_scaled_capital(tea) / ppl, # CAPEX
        lambda: get_annual_electricity()/ppl, # energy (electricity)
        lambda: tea.annual_labor/ppl, # labor
        lambda: (tea.AOC-get_annual_electricity()-tea.annual_labor)/ppl, # OPEX (other than energy and labor)
        lambda: tea.sales / ppl, # sales
        ]


def get_normalized_CAPEX(units, ppl=default_ppl):
    '''Get the CAPEX of a unit/units normalized to per capita per day.'''
    system = units[0].system
    return system.TEA.get_unit_annualized_equipment_cost(units)/365/ppl


def get_noramlized_electricity_cost(units, ppl=default_ppl):
    '''Get the energy (electricity) cost of a unit/units normalized to per capita per day.'''
    return sum(u.power_utility.cost for u in units)/ppl


def get_normalized_OPEX(units, ppl=default_ppl):
    '''
    Get the OPEX of a unit/units normalized to per capita per day,
    energy (electricity) cost is not included.
    '''
    OPEX = sum(u.add_OPEX.values() for u in units)
    streams = sum([u.ins for u in units], [])
    OPEX += sum(s.cost for s in streams)
    return OPEX * 24 / ppl # convert to per capita per day


##### Emissions #####
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
        lambda: lca.total_other_impacts['GlobalWarming']/lca.lifetime/ppl,
        ]


def print_summaries(systems):
    try: iter(systems)
    except: systems = (systems, )

    for sys in systems:
        sys.simulate()

        print(f'\n---------- Summary for {sys.ID} ----------\n')
        if sys.ID in ('sysA', 'sysB'):
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
    'ng_path',
    'data_path',
    'results_path',
    *_components.__all__,
    *systems.__all__,
    *models.__all__,
    *country_specific.__all__,
)