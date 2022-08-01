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
from qsdsan import ImpactItem, StreamImpactItem
from exposan.utils import (
    get_decay_k,
    get_generic_scaled_capital,
    get_generic_tanker_truck_fee as get_tanker_truck_fee,
    )

br_path = os.path.dirname(__file__)
data_path = os.path.join(br_path, 'data')
results_path = os.path.join(br_path, 'results')
# To save simulation data
if not os.path.isdir(results_path): os.mkdir(results_path)


# %%

# =============================================================================
# Unit parameters
# =============================================================================

household_size = 4
household_per_toilet = 4
get_toilet_user = lambda: household_size * household_per_toilet

# Number of people served by the one Biogenic Refinery 1018
def get_ppl(kind):
    if kind.lower()=='10k' or kind.upper()[-1]=='C': return 10000
    return 12000

discount_rate = 0.05

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3

max_CH4_emission = 0.25

emptying_fee = 0.15
handcart_fee = 0.01 # USD/cap/d
truck_fee = 6.21 # USD/m3

get_handcart_and_truck_fee = \
    lambda vol, ppl, include_fee, unit: truck_fee*vol \
        + handcart_fee*ppl*unit.collection_period

# Nutrient loss during application
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05


# =============================================================================
# Prices and GWP CFs
# =============================================================================

# Recycled nutrients are sold at a lower price than commercial fertilizers
price_factor = 0.25

# Should be changed based on country
price_ratio = 1
operator_daily_wage = 29
const_daily_wage = 17
const_person_days = 100

price_dct = {
    'Electricity': 0.06,
    'Concrete': 194*price_ratio,
    'Steel': 2.665*price_ratio,
    'N': 1.507*price_factor,
    'P': 3.983*price_factor,
    'K': 1.333*price_factor,
    'Polymer': 1*price_ratio,
    'Resin': 3.335*price_ratio,
    'FilterBag': 4.81*price_ratio,
    'MgOH2':  0.145*price_ratio,
    'MgCO3': 0.9*price_ratio,
    'H2SO4': 0.3*price_ratio,
    'biochar': 0,  # 0.014*price_ratio,  # assuming value of biochar is 0 for TEA - HACL
    'struvite': 3.983*(31/245)*price_factor,
    'conc_NH3': 1.333*(14/17)*price_factor,
    }

GWP_dct = {
    'Electricity': 0.69,
    'CH4': 34,
    'N2O': 298,
    'N': -5.4,
    'P': -4.9,
    'K': -1.5,
    'Polymer': 2.8,
    'Resin': 1.612,
    'FilterBag': 0.464,  # based on 0.05 kg of nylon and nylon's GWP of 9.279255342 kgCO2eq per kg
    'MgOH2': 1.176277921,
    'MgCO3': 1.176277921,
    'H2SO4': 0.158899487,
    'biochar': -0.2*0.9*(44/12),  # assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
    'struvite': -4.9*(31/245),
    'conc_NH3': -5.4*(14/17),
    }

EcosystemQuality_factor = 29320 * (2.8e-09+7.65e-14)  # (pt/species.yr) * (species.yr/kgCO2eq)
H_Ecosystems_dct = {
    'Electricity': 0.002456338,
    'CH4': 34 * EcosystemQuality_factor,
    'N2O': 298 * EcosystemQuality_factor,
    'N': -0.0461961,
    'P': -0.093269908,
    'K': -0.01895794,
    'Polymer': 0.003527775,
    'Resin': 0.005986888,
    'FilterBag': 0.000360284,  # based on 0.05 kg of nylon and nylon's H_Ecosystems of 0.007205687 points per kg
    'MgOH2': 0.209556136,
    'MgCO3': 0.209556136,
    'H2SO4': 0.000808874,
    'biochar': -0.2*0.9*(44/12)*EcosystemQuality_factor,  # assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
    'struvite': -0.093269908*(31/245),
    'conc_NH3': -0.0461961*(14/17),
    }

HumanHealth_factor = 436000 * 9.28e-07  # (pt/DALY) * (DALY/kgCO2eq)
H_Health_dct = {
    'Electricity': 0.040824307,
    'CH4': 34 * HumanHealth_factor,
    'N2O': 298 * HumanHealth_factor,
    'N': -0.637826734,
    'P': -1.774294425,
    'K': -0.116067637,
    'Polymer': 0.054782882,
    'Resin': 0.094225663,
    'FilterBag': 0.005084823,  # based on 0.05 kg of nylon and nylon's H_Health of 0.10169646 points per kg
    'MgOH2': 4.639146841,
    'MgCO3': 4.639146841,
    'H2SO4': 0.026124187,
    'biochar': -0.2*0.9*(44/12)*HumanHealth_factor,  # assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
    'struvite': -1.774294425*(31/245),
    'conc_NH3': -0.637826734*(14/17),
    }

H_Resources_dct = {
    'Electricity': 0.027825633,
    'CH4': 0,  # no GWP to Resource Depletion pathway
    'N2O': 0,  # no GWP to Resource Depletion pathway
    'N': -0.259196888,
    'P': -1.084191599,
    'K': -0.054033438,
    'Polymer': 0.029951119,
    'Resin': 0.107890228,
    'FilterBag': 0.010811979,  # based on 0.05 kg of nylon and nylon's H_Resources of 0.216239589 points per kg
    'MgOH2': 4.05197164,
    'MgCO3': 4.05197164,
    'H2SO4': 0.025065831,
    'biochar': 0,  # no GWP to Resource Depletion pathway
    'struvite': -1.084191599*(31/245),
    'conc_NH3': -0.259196888*(14/17),
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
        create_stream_impact_item(item_ID='polymer_item', dct_key='Polymer')
        create_stream_impact_item(item_ID='resin_item', dct_key='Resin')
        create_stream_impact_item(item_ID='filter_bag_item', dct_key='FilterBag')
        create_stream_impact_item(item_ID='MgOH2_item')
        create_stream_impact_item(item_ID='MgCO3_item')
        create_stream_impact_item(item_ID='H2SO4_item')
        create_stream_impact_item(item_ID='biochar_item')
        create_stream_impact_item(item_ID='struvite_item')
        create_stream_impact_item(item_ID='conc_NH3_item')
        ImpactItem(ID='e_item', functional_unit='kWh',
                   GWP=GWP_dct['Electricity'],
                   H_Ecosystems=H_Ecosystems_dct['Electricity'],
                   H_Health=H_Health_dct['Electricity'],
                   H_Resources=H_Resources_dct['Electricity'])

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
# Calculate recoveries as in kg C/N/P/K per yr
hr_per_yr = 365 * 24
def get_C(stream):
    sys = stream.source.system
    for unit in sys.units:
        if hasattr(unit, 'carbon_COD_ratio'): carbon_COD_ratio = unit.carbon_COD_ratio
    return stream.COD*stream.F_vol/1e3*carbon_COD_ratio*hr_per_yr
get_C_gas = lambda stream: stream.imol['CH4']*12*hr_per_yr
get_N = lambda stream: stream.TN*stream.F_vol/1e3*hr_per_yr
get_N_gas = lambda stream: stream.imol['N2O']*28*hr_per_yr
get_P = lambda stream: stream.TP*stream.F_vol/1e3*hr_per_yr
get_K = lambda stream: stream.TK*stream.F_vol/1e3*hr_per_yr

# Only for `sysA` or `sysB`
def get_recoveries(system, include_breakdown=False):
    AB = system.ID[-1]
    if AB not in ('A', 'B'):
        if AB in ('C', 'D'): return [lambda: 0]*3 # recoveries all 0 for sysC and sysD
        raise ValueError('This function is only applicable for `sysA`, `sysB`, `sysC`, or `sysD`, '
                         f'not `{system.ID}`.')

    u_reg = system.flowsheet.unit
    ppl = get_ppl(AB)
    dct = globals()
    dct['C_dct'] = C_dct = {}
    dct['N_dct'] = N_dct = {}
    dct['P_dct'] = P_dct = {}
    dct['K_dct'] = K_dct = {}

    ##### Unique to A or B #####
    if AB == 'A':
        sol_num = 0
        toilet = u_reg.A2
        trans_sol = u_reg.A3
        pretreat = u_reg.A6
        liq_bed = u_reg.A7
        pyrolysis = u_reg.A8
        dryer = u_reg.A12
        # Pretreatment
        C_dct['pretreat_liq'] = get_C(pretreat.outs[0])
        N_dct['pretreat_liq'] = get_N(pretreat.outs[0])
        P_dct['pretreat_liq'] = get_P(pretreat.outs[0])
        K_dct['pretreat_liq'] = get_K(pretreat.outs[0])
        # Struvite and ion exchange all 0 for sysA
        N_dct['struvite'] = N_dct['NH3'] = P_dct['struvite'] = 0
    else: # unique to sysB
        sol_num = 1
        toilet = u_reg.B2
        trans_liq = u_reg.B3
        trans_sol = u_reg.B4
        struvite = u_reg.B5
        ix = u_reg.B6
        liq_bed = u_reg.B7
        pretreat = u_reg.B10
        pyrolysis = u_reg.B11
        dryer = u_reg.B15
        # In
        C_dct['toilet_liq'] = get_C(toilet.outs[0])
        N_dct['toilet_liq'] = get_N(toilet.outs[0])
        P_dct['toilet_liq'] = get_P(toilet.outs[0])
        K_dct['toilet_liq'] = get_K(toilet.outs[0])
        # Transported
        C_dct['trans_liq'] = get_C(trans_liq.outs[0])
        C_dct['trans_liq_loss'] = get_C(trans_liq.outs[1])
        N_dct['trans_liq'] = get_N(trans_liq.outs[0])
        N_dct['trans_liq_loss'] = get_N(trans_liq.outs[1])
        P_dct['trans_liq'] = get_P(trans_liq.outs[0])
        P_dct['trans_liq_loss'] = get_P(trans_liq.outs[1])
        K_dct['trans_liq'] = get_K(trans_liq.outs[0])
        K_dct['trans_liq_loss'] = get_K(trans_liq.outs[1])
        # Struvite and ion exchange
        N_dct['struvite'] = struvite.outs[1].imol['Struvite'] * 14 * hr_per_yr # 14 is the MW of N
        N_dct['NH3'] = ix.outs[2].imol['NH3'] * 14 * hr_per_yr
        P_dct['struvite'] = struvite.outs[1].imol['Struvite'] * 31 * hr_per_yr # 31 is the MW of P

    ##### Applicable for both A and B #####
    # In
    C_dct['urine'] = get_C(toilet.ins[0]) * ppl
    C_dct['feces'] = get_C(toilet.ins[1]) * ppl
    C_dct['input'] = C_dct['urine'] + C_dct['feces']
    N_dct['urine'] = get_N(toilet.ins[0]) * ppl
    N_dct['feces'] = get_N(toilet.ins[1]) * ppl
    N_dct['input'] = N_dct['urine'] + N_dct['feces']
    P_dct['urine'] = get_P(toilet.ins[0]) * ppl
    P_dct['feces'] = get_P(toilet.ins[1]) * ppl
    P_dct['input'] = P_dct['urine'] + P_dct['feces']
    K_dct['urine'] = get_K(toilet.ins[0]) * ppl
    K_dct['feces'] = get_K(toilet.ins[1]) * ppl
    K_dct['input'] = K_dct['urine'] + K_dct['feces']

    # Toilet, leachate does not contain COD or C
    C_dct['toilet_sol'] = get_C(toilet.outs[sol_num])
    C_dct['toilet_gas'] = get_C_gas(toilet.outs[-2]) # C in fugitive CH4
    N_dct['toilet_sol'] = get_N(toilet.outs[sol_num])
    N_dct['toilet_gas'] = get_N_gas(toilet.outs[-1]) # N in fugitive N2O
    P_dct['toilet_sol'] = get_P(toilet.outs[sol_num])
    K_dct['toilet_sol'] = get_K(toilet.outs[sol_num])

    # Transported
    C_dct['trans_sol'] = get_C(trans_sol.outs[0])
    C_dct['trans_sol_loss'] = get_C(trans_sol.outs[1])
    N_dct['trans_sol'] = get_N(trans_sol.outs[0])
    N_dct['trans_sol_loss'] = get_N(trans_sol.outs[1])
    P_dct['trans_sol'] = get_P(trans_sol.outs[0])
    P_dct['trans_sol_loss'] = get_P(trans_sol.outs[1])
    K_dct['trans_sol'] = get_K(trans_sol.outs[0])
    K_dct['trans_sol_loss'] = get_K(trans_sol.outs[1])

    # Pretreatment
    sol_num = 1 - sol_num # want 1 for sysA and 0 for sysB
    C_dct['pretreat_sol'] = get_C(pretreat.outs[sol_num])
    N_dct['pretreat_sol'] = get_N(pretreat.outs[sol_num])
    P_dct['pretreat_sol'] = get_P(pretreat.outs[sol_num])
    K_dct['pretreat_sol'] = get_K(pretreat.outs[sol_num])

    # Liq bed
    C_dct['bed_liq'] = get_C(liq_bed.outs[0])
    C_dct['bed_gas'] = get_C_gas(liq_bed.outs[1])
    N_dct['bed_liq'] = get_N(liq_bed.outs[0])
    N_dct['bed_gas'] = get_N_gas(liq_bed.outs[2])
    P_dct['bed_liq'] = get_P(liq_bed.outs[0])
    K_dct['bed_liq'] = get_K(liq_bed.outs[0])

    # Pyrolysis
    C_dct['pyrolysis_biochar'] = pyrolysis.outs[0].imass['C'] * hr_per_yr
    C_dct['pyrolysis_gas'] = get_C(pyrolysis.ins[0]) * pyrolysis.pyrolysis_C_loss
    N_dct['pyrolysis_biochar'] = pyrolysis.outs[0].imass['N'] * hr_per_yr
    N_dct['pyrolysis_gas'] = get_N(pyrolysis.ins[0]) * pyrolysis.pyrolysis_N_loss
    P_dct['pyrolysis_biochar'] = pyrolysis.outs[0].imass['P'] * hr_per_yr
    P_dct['pyrolysis_gas'] = get_P(pyrolysis.ins[0]) * pyrolysis.pyrolysis_P_loss
    K_dct['pyrolysis_biochar'] = pyrolysis.outs[0].imass['K'] * hr_per_yr
    K_dct['pyrolysis_gas'] = get_K(pyrolysis.ins[0]) * pyrolysis.pyrolysis_K_loss

    # Dryer
    C_dct['dryer_sol'] = get_C(dryer.outs[0])
    C_dct['dryer_gas'] = get_C_gas(dryer.outs[2])
    N_dct['dryer_sol'] = get_N(dryer.outs[0])
    N_dct['dryer_gas'] = get_N_gas(dryer.outs[1])
    P_dct['dryer_sol'] = get_P(dryer.outs[0])
    K_dct['dryer_sol'] = get_K(dryer.outs[0])

    # % N, P, and K recovered as a usable fertilizer product,
    # for model metrics and also the Resource Recovery criterion in DMsan analysis
    functions = [
            lambda: (N_dct['bed_liq']+N_dct['NH3']+N_dct['struvite']) / N_dct['input'] * 100, # total_N_recovery
            lambda: (P_dct['bed_liq']+P_dct['struvite']) / P_dct['input'] * 100, # total_P_recovery
            lambda: K_dct['bed_liq'] / K_dct['input'] * 100, # total_K_recovery
            ]
    if not include_breakdown: return functions

    return [
        *functions,
        lambda: N_dct['bed_liq'] / N_dct['input'] * 100, # N_liquid_treatment_bed
        lambda: N_dct['NH3'] / N_dct['input'] * 100, # N_ion_exchange
        lambda: N_dct['struvite'] / N_dct['input'] * 100, # N_struvite
        lambda: P_dct['bed_liq'] / P_dct['input'] * 100, # P_liquid_treatment_bed
        lambda: P_dct['struvite'] / P_dct['input'] * 100, # P_struvite
        lambda: K_dct['bed_liq'] / K_dct['input'] * 100, # K_liquid_treatment_bed
        ]


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
	'br_path',
	'data_path',
	'results_path',
	*_components.__all__,
	*systems.__all__,
    *models.__all__,
    *country_specific.__all__,
 	)