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
    _init_modules,
    get_decay_k,
    get_generic_scaled_capital,
    get_generic_tanker_truck_fee as get_tanker_truck_fee,
    )

# Module-wise setting on whether to allow resource recovery
INCLUDE_RESOURCE_RECOVERY = False

br_path = os.path.dirname(__file__)
module = os.path.split(br_path)[-1]
data_path, results_path = _init_modules(module, include_data_path=True)


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

EcosystemQuality_factor = 29320 * (2.8e-09+7.65e-14) # (pt/species.yr) * (species.yr/kgCO2eq)
HumanHealth_factor = 436000 * 9.28e-07 # (pt/DALY) * (DALY/kgCO2eq)

def update_resource_recovery_settings():
    global INCLUDE_RESOURCE_RECOVERY
    global price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct
    RR_factor = int(bool(INCLUDE_RESOURCE_RECOVERY))
    price_dct = {
        'Electricity': 0.13,
        'Concrete': 194*price_ratio,
        'Steel': 2.665*price_ratio,
        'N': 1.507*price_factor*RR_factor,
        'P': 3.983*price_factor*RR_factor,
        'K': 1.333*price_factor*RR_factor,
        'Polymer': 1*price_ratio,
        'Resin': 3.335*price_ratio,
        'FilterBag': 4.81*price_ratio,
        'MgOH2':  0.145*price_ratio,
        'MgCO3': 0.9*price_ratio,
        'H2SO4': 0.3*price_ratio,
        'biochar': 0*RR_factor,  # 0.014*price_ratio,  # assuming value of biochar is 0 for TEA - HACL
        'struvite': 3.983*(31/245)*price_factor*RR_factor,
        'conc_NH3': 1.333*(14/17)*price_factor*RR_factor,
        }

    GWP_dct = {
        'Electricity': 0.69,
        'CH4': 34,
        'N2O': 298,
        'N': -5.4*RR_factor,
        'P': -4.9*RR_factor,
        'K': -1.5*RR_factor,
        'Polymer': 2.8,
        'Resin': 1.612,
        'FilterBag': 0.464,  # based on 0.05 kg of nylon and nylon's GWP of 9.279255342 kgCO2eq per kg
        'MgOH2': 1.176277921,
        'MgCO3': 1.176277921,
        'H2SO4': 0.158899487,
        #!!! Assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
        'biochar': -0.2*0.9*(44/12)*RR_factor,
        'struvite': -4.9*(31/245)*RR_factor,
        'conc_NH3': -5.4*(14/17)*RR_factor,
        }

    H_Ecosystems_dct = {
        'Electricity': 0.002456338,
        'CH4': 34 * EcosystemQuality_factor,
        'N2O': 298 * EcosystemQuality_factor,
        'N': -0.0461961*RR_factor,
        'P': -0.093269908*RR_factor,
        'K': -0.01895794*RR_factor,
        'Polymer': 0.003527775,
        'Resin': 0.005986888,
        'FilterBag': 0.000360284,  # based on 0.05 kg of nylon and nylon's H_Ecosystems of 0.007205687 points per kg
        'MgOH2': 0.209556136,
        'MgCO3': 0.209556136,
        'H2SO4': 0.000808874,
        #!!! Assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
        'biochar': -0.2*0.9*(44/12)*EcosystemQuality_factor*RR_factor,
        'struvite': -0.093269908*(31/245)*RR_factor,
        'conc_NH3': -0.0461961*(14/17)*RR_factor,
        }

    H_Health_dct = {
        'Electricity': 0.040824307,
        'CH4': 34 * HumanHealth_factor,
        'N2O': 298 * HumanHealth_factor,
        'N': -0.637826734*RR_factor,
        'P': -1.774294425*RR_factor,
        'K': -0.116067637*RR_factor,
        'Polymer': 0.054782882,
        'Resin': 0.094225663,
        'FilterBag': 0.005084823,  # based on 0.05 kg of nylon and nylon's H_Health of 0.10169646 points per kg
        'MgOH2': 4.639146841,
        'MgCO3': 4.639146841,
        'H2SO4': 0.026124187,
        #!!! Assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
        'biochar': -0.2*0.9*(44/12)*HumanHealth_factor*RR_factor,
        'struvite': -1.774294425*(31/245)*RR_factor,
        'conc_NH3': -0.637826734*(14/17)*RR_factor,
        }

    H_Resources_dct = {
        'Electricity': 0.027825633,
        'CH4': 0,  # no GWP to Resource Depletion pathway
        'N2O': 0,  # no GWP to Resource Depletion pathway
        'N': -0.259196888*RR_factor,
        'P': -1.084191599*RR_factor,
        'K': -0.054033438*RR_factor,
        'Polymer': 0.029951119,
        'Resin': 0.107890228,
        'FilterBag': 0.010811979,  # based on 0.05 kg of nylon and nylon's H_Resources of 0.216239589 points per kg
        'MgOH2': 4.05197164,
        'MgCO3': 4.05197164,
        'H2SO4': 0.025065831,
        'biochar': 0*RR_factor,  # no GWP to Resource Depletion pathway
        'struvite': -1.084191599*(31/245)*RR_factor,
        'conc_NH3': -0.259196888*(14/17)*RR_factor,
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
            StreamImpactItem(ID=item_ID,
                             GWP=GWP_dct[dct_key],
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
        raise AttributeError(
            f'Module {__name__} does not have the attribute "{name}" '
            'and the module has not been loaded, '
            f'loading the module with `{__name__}.load()` may solve the issue.')


# %%

# =============================================================================
# Util functions
# =============================================================================

import numpy as np

##### Recoveries #####

# Only for sysA, sysB, sysC
def get_recoveries(system, include_breakdown=False):
    sysID = system.ID[-1]
    if sysID not in ('A', 'B', 'C'):
        if sysID in ('D'): return [lambda: 0]*3 # recoveries all 0 for sysD
        raise ValueError('This function is only applicable for `sysA`, `sysB`, `sysC`, or `sysD`, '
                         f'not `{system.ID}`.')

    u = system.flowsheet.unit
    ppl = get_ppl(sysID)
    hr_per_yr = 365 * 24

    ##### Unique to A or B #####
    if sysID == 'A':
        toilet = u.A2
        liq_bed = u.A7
        pyrolysis = u.A8        
        N_struvite = N_NH3 = P_struvite = 0 # Struvite and ion exchange all 0 for sysA
    elif sysID == 'B': # unique to sysB
        toilet = u.B2
        struvite = u.B5
        ix = u.B6
        liq_bed = u.B7
        pyrolysis = u.B11
        # Struvite and ion exchange
        N_struvite = struvite.outs[1].imol['Struvite'] * 14 * hr_per_yr # 14 is the MW of N
        N_NH3 = ix.outs[2].imol['NH3'] * 14 * hr_per_yr
        P_struvite = struvite.outs[1].imol['Struvite'] * 31 * hr_per_yr # 31 is the MW of P
    elif sysID == 'C':
        toilet = u.C2
        liq_bed = u.C7
        pyrolysis = u.C8
        N_struvite = N_NH3 = P_struvite = 0 # Struvite and ion exchange all 0 for sysC
    else:
        toilet = u.D2

    # Calculate recoveries in kg C/N/P/K per yr
    # Inputs
    COD_urine = toilet.ins[0].COD * toilet.ins[0].F_vol
    COD_feces = toilet.ins[1].COD * toilet.ins[1].F_vol
    C_input = (COD_urine + COD_feces) / 1e3 * pyrolysis.carbon_COD_ratio * hr_per_yr * ppl
    N_urine = toilet.ins[0].TN * toilet.ins[0].F_vol
    N_feces = toilet.ins[1].TN * toilet.ins[1].F_vol 
    N_input = (N_urine + N_feces) / 1e3 * hr_per_yr * ppl 
    P_urine = toilet.ins[0].TP * toilet.ins[0].F_vol
    P_feces = toilet.ins[1].TP * toilet.ins[1].F_vol
    P_input = (P_urine + P_feces) / 1e3 * hr_per_yr * ppl 
    K_urine = toilet.ins[0].TK * toilet.ins[0].F_vol
    K_feces = toilet.ins[1].TK * toilet.ins[1].F_vol
    K_input = (K_urine + K_feces) / 1e3 * hr_per_yr * ppl 

    # Liq bed
    N_bed_liq = liq_bed.outs[0].TN * liq_bed.outs[0].F_vol / 1e3 * hr_per_yr 
    P_bed_liq = liq_bed.outs[0].TP * liq_bed.outs[0].F_vol / 1e3 * hr_per_yr
    K_bed_liq = liq_bed.outs[0].TK * liq_bed.outs[0].F_vol / 1e3 * hr_per_yr

    # Pyrolysis 
    C_pyrolysis_biochar = pyrolysis.outs[0].imass['C'] * hr_per_yr


    # % N, P, and K recovered as a usable fertilizer product, %C recovered as biochar
    # for model metrics and also the Resource Recovery criterion in DMsan analysis
    functions = [
            lambda: C_pyrolysis_biochar / C_input * 100, # total_C_recovery
            lambda: (N_bed_liq + N_struvite + N_NH3) / N_input * 100, # total_N_recovery
            lambda: (P_bed_liq+P_struvite) / P_input * 100, # total_P_recovery
            lambda: K_bed_liq / K_input * 100, # total_K_recovery
            ]
    if not include_breakdown: return functions

    return [
        *functions,
        lambda: C_pyrolysis_biochar / C_input * 100, # C_carbonizer_base
        lambda: N_bed_liq / N_input * 100, # N_liquid_treatment_bed
        lambda: N_NH3 / N_input * 100, # N_ion_exchange
        lambda: N_struvite / N_input * 100, # N_struvite
        lambda: P_bed_liq / P_input * 100, # P_liquid_treatment_bed
        lambda: P_struvite / P_input * 100, # P_struvite
        lambda: K_bed_liq / K_input * 100, # K_liquid_treatment_bed
        ]

##### Biochar functions #####

# Calculate dry basis biochar yield (% mass)
def get_yield(temp, AC): 
    f_AC_dec = AC/100 #converts % ash content of feedstock to decimal
    
    # predictive equation for % biochar dry basis (db) yield - derived via Excel solver
    db_yield = 100 * (1.18 * f_AC_dec ** 0.843 + (1 - f_AC_dec) * 2.106 * np.exp(-0.0066 * temp))
    return db_yield
   
# Calculate carbon sequestration potential of biochar produced (% mass)
def get_CS(temp, AC): 
    db_yield = get_yield(temp, AC)  
    
    # predictive equation for % biochar dry ash-free (daf) yield - Neves et al. 2011
    daf_yield = 100 * (0.106 + 2.43 * np.exp(-0.0066 * temp))
    
    # predictive equation for % biochar fixed carbon - derived via Excel solver
    FC_biochar = 87.786 * daf_yield ** -0.483

    # calculate biochar volatile matter and ash content using calculated values from above eqns
    AC_biochar = (db_yield - daf_yield) * 100 / db_yield
    VM_biochar = 100 - AC_biochar - FC_biochar
    
    # calculations for carbon sequestration (CS) potential [% mass C/mass biochar] 
    C_biochar = (0.474 * VM_biochar + 0.963 * FC_biochar + 0.067 * AC_biochar) / 100 # Klasson 2017
    Cafb = (0.474 * VM_biochar + 0.963 * FC_biochar + 0.067 * AC_biochar) / (100 - AC_biochar) # Klasson 2017
    C_feedstock = -0.50 * AC + 54.51 # Krueger et al. 2021
    R50 = 0.17 * Cafb + 0.00479 # Klasson 2017
    CS = db_yield * (C_biochar*100) * R50 / C_feedstock # Zhao et al. 2013
    return CS  


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