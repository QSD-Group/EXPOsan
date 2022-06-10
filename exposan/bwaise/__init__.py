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
from qsdsan import ImpactItem, StreamImpactItem
from qsdsan.utils import AttrGetter, FuncGetter, time_printer

bwaise_path = os.path.dirname(__file__)
data_path = os.path.join(bwaise_path, 'data')
results_path = os.path.join(bwaise_path, 'results')
figures_path = os.path.join(bwaise_path, 'figures')
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
# Get reduction rate constant k for COD and N, use a function so that k can be
# changed during uncertainty analysis
def get_decay_k(tau_deg=2, log_deg=3):
    k = (-1/tau_deg)*np.log(10**-log_deg)
    return k

max_CH4_emission = 0.25

# Model for tanker truck cost based on capacity (m3)
# price = a*capacity**b -> ln(price) = ln(a) + bln(capacity)
from sklearn.linear_model import LinearRegression
UGX_price_dct = np.array((8e4, 12e4, 20e4, 25e4))
capacities = np.array((3, 4.5, 8, 15))
emptying_fee = 0.15 # additional emptying fee, fraction of base cost
def get_tanker_truck_fee(capacity):
    price_dct = UGX_price_dct*(1+emptying_fee)/exchange_rate
    ln_p = np.log(price_dct)
    ln_cap = np.log(capacities)
    model = LinearRegression().fit(ln_cap.reshape(-1,1), ln_p.reshape(-1,1))
    predicted = model.predict(np.array((np.log(capacity))).reshape(1, -1)).item()
    cost = np.exp(predicted)
    return cost

# Flow rates for treatment plants
sewer_flow = 2750 # m3/d
sludge_flow_exist = 500 # m3/d
sludge_flow_alt = 60 # m3/d
get_sludge_flow = lambda kind: \
    sludge_flow_exist if kind.lower() in ('exist', 'sysa', 'sysc', 'a', 'c') else sludge_flow_alt

# Nutrient loss during application
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05

exchange_rate = 3700 # UGX per USD

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
    'Biogas': -3*get_biogas_factor()
    }


# %%

# =============================================================================
# Load components and systems
# =============================================================================

from . import _components, _lca_data, systems#, models
from ._components import *
from ._lca_data import *
from .systems import *

_components_loaded = False
_system_loaded = False
_impact_item_loaded = False


def _load_components():
    global components, _components_loaded
    components = create_components()
    qs.set_thermo(components)
    _components_loaded = True


def _load_lca_data(lca_kind='original'):
    '''
    Load impact indicator and impact item data.

    Parameters
    ----------
    lca_kind : str
        "original" loads the data from Trimmer et al.
        (TRACI, ecoinvent v3.2, GWP only),
        "new" loads the data for ReCiPe and TRACI
        (ecoinvent 3.7.1, at the point of substitution).
    '''
    indicator_path = os.path.join(data_path, f'indicators_{lca_kind}.tsv')
    indel_col = None if lca_kind=='original' else 0
    ind_df_processed = pd.read_csv(indicator_path, sep='\t', index_col=indel_col)
    qs.ImpactIndicator.load_from_file(indicator_path)

    global _impact_item_loaded
    _impact_item_loaded = lca_kind

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

    '''
    try: # prevent from reloading
        Excavation = ImpactItem.get_item('Excavation')
        Excavation.indicators[0]
    except:
        _load_lca_data('original')
    '''


def _load_system():
    qs.currency = 'USD'
    qs.CEPCI = qs.CEPCI_by_year[2018]
    qs.PowerUtility.price = price_dct['Electricity']
    global sysA, sysB, sysC, _system_loaded
    sysA = create_system('A')
    sysB = create_system('B')
    sysC = create_system('C')
    _system_loaded = True

def load(lca_kind='original'):
    if not _components_loaded: _load_components()
    if _impact_item_loaded is False : _load_lca_data(lca_kind)
    if not _system_loaded: _load_system()
    if _impact_item_loaded != lca_kind: # to refresh the impact items
        for sys in (sysA, sysB, sysC):
            lca = sys.LCA
            for i in lca.lca_streams:
                source_ID = i.stream_impact_item.source.ID
                i.stream_impact_item.source = ImpactItem.get_item(source_ID)
    dct = globals()
    for sys in (sysA, sysB, sysC): dct.update(sys.flowsheet.to_dict())


def __getattr__(name):
    if not _components_loaded:
        _load_components()
        if name in ('components', 'cmps'): return components
    if not _system_loaded:
        _load_system()
        dct = globals()
        for sys in (sysA, sysB, sysC): dct.update(sys.flowsheet.to_dict())
        if name in dct: return dct[name]
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")



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

    ppl = get_ppl(system.ID)

    attr = f'T{resource}' if resource != 'COD' else resource
    sum_attr = lambda streams: sum(getattr(s, attr)*s.F_vol/1e3 for s in streams)

    return (
        lambda: sum_attr(liq)/sum_attr(ins)/ppl,
        lambda: sum_attr(sol)/sum_attr(ins)/ppl,
        lambda: get_gas_COD()/sum_attr(ins)/ppl,
        lambda: (sum_attr(liq)+sum_attr(sol)+get_gas_COD())/sum_attr(ins)/ppl,
        )


def print_summaries(systems):
    try: iter(systems)
    except: systems = (systems,)

    for sys in systems:
        # func = get_summarizing_functions(sys)
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
        ppl = get_ppl(sys.ID)
        tea = sys.TEA
        tea.show()

        unit = f'{qs.currency}/cap/yr'
        val = (tea.annualized_equipment_cost-tea.net_earnings)/ppl
        print(f'\nNet cost: {val:.1f} {unit}.')

        val = tea.annualized_equipment_cost/ppl
        print(f'Capital: {val:.1f} {unit}.')

        val = tea.AOC/ppl
        print(f'Operating: {val:.1f} {unit}.')

        val = tea.sales/ppl
        print(f'Sales: {val:.1f} {unit}.')
        print('\n')

        # LCA
        lca = sys.LCA
        lca.show()

        print('\n')
        for ind in lca.indicators:
            unit = f'{ind.unit}/cap/yr'
            print(f'\nImpact indicator {ind.ID}:')

            val = lca.total_impacts[ind.ID]/lca.lifetime/ppl
            print(f'\nNet emission: {val:.1f} {unit}.')

            val = lca.total_construction_impacts[ind.ID]/lca.lifetime/ppl
            print(f'Construction: {val:.1f} {unit}.')

            val = lca.total_transportation_impacts[ind.ID]/lca.lifetime/ppl
            print(f'Transportation: {val:.1f} {unit}.')

            val = lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='direct_emission')[ind.ID] \
                /lca.lifetime/ppl
            print(f'Direct emission: {val:.1f} {unit}.')

            val = lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='offset')[ind.ID] \
                /lca.lifetime/ppl
            print(f'Offset: {val:.1f} {unit}.')

            val = lca.total_other_impacts[ind.ID]/lca.lifetime/ppl
            print(f'Other: {val:.1} {unit}.\n')


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



#!!! from . import models
# from .models import *


__all__ = (
    'evaluate',
    'get_key_metrics',
    'bwaise_path',
    'data_path',
    'results_path',
    'figures_path',
 	*_components.__all__,
    *_lca_data.__all__,
 	*systems.__all__,
  #   *models.__all__,
 	)