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


import os, pickle, numpy as np, pandas as pd, qsdsan as qs

import biosteam

from qsdsan import ImpactItem, StreamImpactItem
from qsdsan.utils import time_printer
from exposan.utils import get_decay_k, get_generic_tanker_truck_fee

bw_path = os.path.dirname(__file__)
data_path = os.path.join(bw_path, 'data')
results_path = os.path.join(bw_path, 'results')
figures_path = os.path.join(bw_path, 'figures')
# To save simulation data and generated figures
if not os.path.isdir(results_path): os.mkdir(results_path)
if not os.path.isdir(figures_path): os.mkdir(figures_path)


# %%

# =============================================================================
# Unit parameters
# =============================================================================

household_size = 4
household_per_toilet = 4
get_toilet_user = lambda: household_size * household_per_toilet

# Number of people served by the existing plant (sysA and sysC)
ppl_exist_sewer = 4e4
ppl_exist_sludge = 416667
# Number of people served by the alternative plant (sysB)
ppl_alt = 5e4
def get_ppl(kind):
    if kind.lower() in ('exist', 'existing', 'sysa', 'sysc', 'a', 'c'):
        return ppl_exist_sewer+ppl_exist_sludge
    elif kind.lower() in ('alt', 'alternative', 'sysb', 'b'):
        return ppl_alt
    else:
        raise ValueError('`kind` should be "exist" (for sysA and sysC)'
                         f'or "alt" for sysB, not {kind}.')

exchange_rate = 3700 # UGX per USD
discount_rate = 0.05

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3

max_CH4_emission = 0.25

# Model for tanker truck cost based on capacity (m3)
# price = a*capacity**b -> ln(price) = ln(a) + bln(capacity)
fitting_dct = {
    3: 8e4,
    4.5: 12e4,
    8: 20e4,
    15: 25e4,
}
emptying_fee = 0.15 # additional emptying fee, fraction of base cost
def get_tanker_truck_fee(capacity):
    return get_generic_tanker_truck_fee(
        capacity=capacity,
        fitting_dct=fitting_dct,
        emptying_fee=emptying_fee,
        exchange_rate=1/exchange_rate,
        )

handcart_fee = 0.01 # USD/cap/d
truck_fee = 23e3 # UGX/m3

# Handcart fee is for both liquid/solid
get_handcart_and_truck_fee = \
    lambda vol, ppl, include_fee, unit: truck_fee/exchange_rate*vol \
        + int(include_fee)*handcart_fee*ppl*unit.collection_period

# Flow rates for treatment plants
sewer_flow = 2750 # m3/d
sludge_flow_exist = 500 # m3/d
sludge_flow_alt = 60 # m3/d
get_sludge_flow = lambda kind: \
    sludge_flow_exist if kind.lower() in ('exist', 'sysa', 'sysc', 'a', 'c') else sludge_flow_alt

# Nutrient loss during application
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05

# Recycled nutrients are sold at a lower price than commercial fertilizers
price_factor = 0.25

# Energetic content of the biogas
# biogas_energy in kJ/mol (as CH4, 16 is the MW of CH4),
# LPG_energy in MJ/kg
biogas_energy = 803
LPG_energy = 50
get_biogas_factor = lambda biogas_energy=biogas_energy, LPG_energy=LPG_energy: biogas_energy/16/LPG_energy

# Labor cost for sysB
skilled_num = 5
skilled_salary = 5e6 # UGX/month
get_tot_skilled_salary = lambda: skilled_salary*skilled_num

unskilled_num = 5
unskilled_salary = 75e4 # UGX/month
get_tot_unskilled_salary = lambda: unskilled_salary*unskilled_num
get_alt_salary = lambda: (get_tot_skilled_salary()+get_tot_unskilled_salary())*12/exchange_rate

price_dct = {
    'Electricity': 0.17,
    'Concrete': 194,
    'Steel': 2.665,
    'N': 1.507*price_factor,
    'P': 3.983*price_factor,
    'K': 1.333*price_factor,
    'Biogas': 6500/exchange_rate*get_biogas_factor()
    }

GWP_dct = {
    'Electricity': 0.1135,
    'CH4': 28,
    'N2O': 265,
    'N': -5.4,
    'P': -4.9,
    'K': -1.5,
    'Biogas': -3*get_biogas_factor(),
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


from . import _lca_data
from ._lca_data import *
_impact_item_loaded = False
def _load_lca_data(lca_kind='original', reload=False):
    '''
    Load impact indicator and impact item data.

    Parameters
    ----------
    lca_kind : str
        "original" loads the data from Trimmer et al.
        (TRACI, ecoinvent v3.2, GWP only),
        "new" loads the data for ReCiPe and TRACI
        (ecoinvent 3.7.1, at the point of substitution).
    reload : bool
        Whether to force reload LCA data.
    '''
    global _impact_item_loaded
    if _impact_item_loaded != lca_kind or reload:
        indicator_path = os.path.join(data_path, f'indicators_{lca_kind}.tsv')
        indel_col = None if lca_kind=='original' else 0
        ind_df_processed = pd.read_csv(indicator_path, sep='\t', index_col=indel_col)
        qs.ImpactIndicator.load_from_file(indicator_path)

        if lca_kind.lower() in ('original', 'traci'):
            item_path = os.path.join(data_path, 'items_original.xlsx')
            qs.ImpactItem.load_from_file(item_path)
            # Electricity and stream impact items
            for k, v in GWP_dct.items():
                if k == 'Electricity':
                    ImpactItem(ID='E_item', functional_unit='kWh', GWP=v)
                else:
                    StreamImpactItem(ID=f'{k}_item', GWP=v)
        elif lca_kind.lower() in ('new', 'recipe'):
            item_path = os.path.join(data_path, 'cf_dct.pckl')
            f = open(item_path, 'rb')
            cf_dct = pickle.load(f)
            f.close()
            create_items(ind_df_processed, cf_dct)

            # Fugitive
            EcosystemQuality_factor = 29320 # pt/species/yr
            HumanHealth_factor = 436000 # pt/DALY

            E_factor = {
                # Global warming to (terrestrial+freshwater) ecosystem
                'GW2ECO': (2.5e-08+6.82e-13)*EcosystemQuality_factor,
                'GW2HH': 1.25e-05*HumanHealth_factor, # global warming to human health
                'OD2HH': 0.00134*HumanHealth_factor, # stratospheric ozone depletion to human health
                        }
            H_factor = {
                'GW2ECO': (2.8e-09+7.65e-14)*EcosystemQuality_factor,
                'GW2HH': 9.28e-07*HumanHealth_factor,
                'OD2HH': 0.000531*HumanHealth_factor,
                }
            I_factor = {
                'GW2ECO': (5.32e-10+1.45e-14)*EcosystemQuality_factor,
                'GW2HH': 8.12e-08*HumanHealth_factor,
                'OD2HH': 0.000237*HumanHealth_factor,
                }

            StreamImpactItem(ID='CH4_item',
                             E_EcosystemQuality_Total=E_factor['GW2ECO']*4.8,
                             E_HumanHealth_Total=E_factor['GW2HH']*4.8,
                             H_EcosystemQuality_Total=H_factor['GW2ECO']*34,
                             H_HumanHealth_Total=H_factor['GW2HH']*34,
                             I_EcosystemQuality_Total=I_factor['GW2ECO']*84,
                             I_HumanHealth_Total=I_factor['GW2HH']*84
                             )
            StreamImpactItem(ID='N2O_item',
                             E_EcosystemQuality_Total=E_factor['GW2ECO']*78.8,
                             # From climate change + ozone depletion
                             E_HumanHealth_Total=\
                                 E_factor['GW2HH']*78.8+E_factor['OD2HH']*0.017,
                             H_EcosystemQuality_Total=H_factor['GW2ECO']*298,
                             H_HumanHealth_Total=\
                                 H_factor['GW2HH']*298+H_factor['OD2HH']*0.011,
                             I_EcosystemQuality_Total=I_factor['GW2ECO']*264,
                             I_HumanHealth_Total=\
                                 I_factor['GW2HH']*264+I_factor['OD2HH']*0.007
                             )
        else: raise ValueError(f'`kind` can only be "original" or "new", not "{lca_kind}".')

        # Update prices
        ImpactItem.get_item('Concrete').price = price_dct['Concrete']
        ImpactItem.get_item('Steel').price = price_dct['Steel']
        Biogas_CFs = ImpactItem.get_item('Biogas_item').CFs
        for k, v in Biogas_CFs.items(): Biogas_CFs[k] = v * get_biogas_factor()

        _impact_item_loaded = lca_kind

    return _impact_item_loaded


from . import systems
from .systems import *
_system_loaded = False
def _load_system(lca_kind='original'):
    qs.currency = 'USD'
    qs.CEPCI = qs.CEPCI_by_year[2018]
    qs.PowerUtility.price = price_dct['Electricity']
    global sysA, sysB, sysC, teaA, teaB, teaC, lcaA, lcaB, lcaC, _system_loaded
    sysA = create_system('A', lca_kind=lca_kind)
    teaA = sysA.TEA
    lcaA = sysA.LCA
    sysB = create_system('B', lca_kind=lca_kind)
    teaB = sysB.TEA
    lcaB = sysB.LCA
    sysC = create_system('C', lca_kind=lca_kind)
    teaC = sysC.TEA
    lcaC = sysC.LCA
    _system_loaded = True


def load(lca_kind='original'):
    if not _components_loaded: _load_components()
    if not _system_loaded: _load_system(lca_kind)
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

def get_recoveries(system, resource):
    sys_ID = system.ID
    if sys_ID=='sysB' and resource=='COD':
        get_gas_COD = lambda: unit.B15.outs[0].imol['CH4']*1e3*biogas_energy/14e3
    else: get_gas_COD = lambda: 0

    unit = system.flowsheet.unit
    get_outs = lambda unit: tuple(s for s in unit.outs if (s.phase!='g' and 'los' not in s.ID))
    if sys_ID == 'sysA':
        ins = unit.A1.outs # outs from Excretion
        liq = get_outs(unit.A13)
        sol = get_outs(unit.A12)
    elif sys_ID == 'sysB':
        ins = unit.B1.outs
        liq = get_outs(unit.B13)
        sol = get_outs(unit.B12)
    else: # sysC
        ins = unit.C1.outs
        liq = get_outs(unit.C14)
        sol = get_outs(unit.C13)

    attr = f'T{resource}' if resource != 'COD' else resource
    sum_attr = lambda streams: sum(getattr(s, attr)*s.F_vol/1e3 for s in streams)

    return [
        lambda: sum_attr(liq)/sum_attr(ins)/get_ppl(sys_ID),
        lambda: sum_attr(sol)/sum_attr(ins)/get_ppl(sys_ID),
        lambda: get_gas_COD()/sum_attr(ins)/get_ppl(sys_ID),
        lambda: (sum_attr(liq)+sum_attr(sol)+get_gas_COD())/sum_attr(ins)/get_ppl(sys_ID),
        ]

# added this net_earings property 7/11 as suggested by Yoel.
@property
def net_earnings(self):
    """Net earnings without accounting for annualized depreciation."""
    net_earnings = self.sales - self.AOC
    if net_earnings < 0:
        return net_earnings
    else:
        return (1 - self.income_tax) * net_earnings
qs.SimpleTEA.net_earnings = net_earnings

def get_TEA_metrics(system):
    sys_ID = system.ID
    tea = system.TEA
    return [
        lambda: (tea.annualized_equipment_cost-tea.net_earnings)/get_ppl(sys_ID),
        lambda: tea.annualized_equipment_cost/get_ppl(sys_ID),
        lambda: tea.AOC/get_ppl(sys_ID),
        lambda: tea.sales/get_ppl(sys_ID),
        ]

def get_LCA_metrics(system, indicator):
    lca = system.LCA
    sys_ID = system.ID
    ind_ID = indicator.ID
    return [
        lambda: lca.total_impacts[ind_ID]/lca.lifetime/get_ppl(sys_ID),
        lambda: lca.total_construction_impacts[ind_ID]/lca.lifetime/get_ppl(sys_ID),
        lambda: lca.total_transportation_impacts[ind_ID]/lca.lifetime/get_ppl(sys_ID),
        lambda: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='direct_emission')[ind_ID] \
            /lca.lifetime/get_ppl(sys_ID),
        lambda: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='offset')[ind_ID] \
            /lca.lifetime/get_ppl(sys_ID),
        lambda: lca.total_other_impacts[ind_ID]/lca.lifetime/get_ppl(sys_ID)
        ]


def print_summaries(systems):
    try: iter(systems)
    except: systems = (systems,)

    for sys in systems:
        sys.simulate()
        print(f'\n---------- Summary for {sys} ----------\n')
        # Recoveries
        for i in ('COD', 'N', 'P', 'K'):
            funcs = get_recoveries(sys, i)
            print(f'Total {i} recovery is {funcs[-1]():.1%}, '
                  f'{funcs[0]():.1%} in liquid, '
                  f'{funcs[1]():.1%} in solid, '
                  f'{funcs[2]():.1%} in gas.')
        print('\n')

        # TEA
        tea = sys.TEA
        tea.show()
        funcs = get_TEA_metrics(sys)
        unit = f'{qs.currency}/cap/yr'
        print(f'\nNet cost: {funcs[0]():.1f} {unit}.')
        print(f'Capital: {funcs[1]():.1f} {unit}.')
        print(f'Operating: {funcs[2]():.1f} {unit}.')
        print(f'Sales: {funcs[3]():.1f} {unit}.')
        print('\n')

        # LCA
        lca = sys.LCA
        lca.show()
        print('\n')
        for ind in lca.indicators:
            funcs = get_LCA_metrics(sys, ind)
            unit = f'{ind.unit}/cap/yr'
            print(f'\nImpact indicator {ind.ID}:')
            print(f'\nNet emission: {funcs[0]():.1f} {unit}.')
            print(f'Construction: {funcs[1]():.1f} {unit}.')
            print(f'Transportation: {funcs[2]():.1f} {unit}.')
            print(f'Direct emission: {funcs[3]():.1f} {unit}.')
            print(f'Offset: {funcs[4]():.1f} {unit}.')
            print(f'Other: {funcs[5]():.1} {unit}.\n')


# Need to be imported last
from . import models
from .models import *

@time_printer
def evaluate(model, samples=None):
    if samples is not None:
        model.load_samples(samples)
    model.evaluate()

def get_key_metrics(model, alt_names={}):
    key_metrics = [i for i in model.metrics if 'total' in i.name.lower()]
    key_metrics += [i for i in model.metrics if 'net' in i.name.lower()]
    for old, new in alt_names.items():
        for i in key_metrics:
            i.name = i.name.replace(old, new)
    return key_metrics


__all__ = (
    'bw_path',
    'data_path',
    'results_path',
    'figures_path',
    'evaluate',
    'get_key_metrics',
 	*_components.__all__,
    *_lca_data.__all__,
 	*systems.__all__,
    *models.__all__,
 	)