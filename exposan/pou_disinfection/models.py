#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import os, pickle
import numpy as np
import pandas as pd
from chaospy import distributions as shape
from thermosteam.functional import V_to_rho, rho_to_V
from biosteam import PowerUtility
from biosteam.evaluation import Model, Metric
from qsdsan import currency, ImpactItem
from qsdsan.utils import (
    ospath, load_data, data_path, dct_from_str,
    AttrSetter, AttrFuncSetter, DictAttrSetter,
    FuncGetter,
    time_printer
    )
from exposan import POU_dis as pou

c_path = pou._lca_data.c_path
lca_data_kind = pou.systems.lca_data_kind

__all__ = ('modelD','modelE')


# # %%

# # =============================================================================
# # Functions for batch-making metrics and -setting parameters
# # =============================================================================

systems = pou.systems
sys_dct = systems.sys_dct
price_dct = systems.price_dct
GWP_dct = systems.GWP_dct
get_summarizing_functions = systems.get_summarizing_functions

def add_LCA_metrics(system, metrics, kind):
    systems.update_lca_data(kind)
    lca = sys_dct['LCA'][system.ID]
    ppl = sys_dct['ppl'][system.ID]
    funcs = [
        lambda ID: lca.total_impacts[ID]/lca.lifetime/ppl,
        lambda ID: lca.total_construction_impacts[ID]/lca.lifetime/ppl,
        lambda ID: lca.total_transportation_impacts[ID]/lca.lifetime/ppl,
        lambda ID: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='direct_emission')[ID] \
            /lca.lifetime/ppl,
        lambda ID: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='offset')[ID] \
            /lca.lifetime/ppl,
        lambda ID: lca.total_other_impacts[ID]/lca.lifetime/ppl
        ]

    for ind in lca.indicators:
        unit = f'{ind.unit}/cap/yr'
        cat = 'LCA results'
        metrics.extend([
            Metric(f'Net emission {ind.ID}', FuncGetter(funcs[0], (ind.ID,)), unit, cat),
            Metric(f'Construction {ind.ID}', FuncGetter(funcs[1], (ind.ID,)),unit, cat),
            Metric(f'Transportation {ind.ID}', FuncGetter(funcs[2], (ind.ID,)),unit, cat),
            Metric(f'Direct emission {ind.ID}', FuncGetter(funcs[3], (ind.ID,)),unit, cat),
            Metric(f'Offset {ind.ID}', FuncGetter(funcs[4], (ind.ID,)),unit, cat),
            Metric(f'Other {ind.ID}', FuncGetter(funcs[5], (ind.ID,)),unit, cat),
            ])

    return metrics

def add_metrics(system, kind):
    sys_ID = system.ID
    tea = sys_dct['TEA'][sys_ID]
    ppl = sys_dct['ppl'][sys_ID]
    func = get_summarizing_functions(system)

    metrics = []

    unit = f'{currency}/cap/yr'
    cat = 'TEA results'
    metrics.extend([
        Metric('Annual net cost', lambda: func['get_annual_net_cost'](tea, ppl), unit, cat),
        Metric('Annual CAPEX', lambda: func['get_annual_CAPEX'](tea, ppl), unit, cat),
        Metric('Annual OPEX', lambda: func['get_annual_OPEX'](tea, ppl), unit, cat),
        Metric('Annual sales', lambda: func['get_annual_sales'](tea, ppl), unit, cat)
        ])

    metrics = add_LCA_metrics(system, metrics, kind)

    return metrics

#!!! leave out update_metrics for now, may need to add in later
# def update_metrics(model, kind):
#     metrics = [i for i in model.metrics if i.element_name!='LCA results']
#     model.metrics = add_LCA_metrics(model.system, metrics, kind)
#     return model


def batch_setting_unit_params(df, model, unit, exclude=()):
    for para in df.index:
        if para in exclude: continue
        b = getattr(unit, para)
        lower = float(df.loc[para]['low'])
        upper = float(df.loc[para]['high'])
        dist = df.loc[para]['distribution']
        if dist == 'uniform':
            D = shape.Uniform(lower=lower, upper=upper)
        elif dist == 'triangular':
            D = shape.Triangle(lower=lower, midpoint=b, upper=upper)
        elif dist == 'constant': continue
        else:
            raise ValueError(f'Distribution {dist} not recognized for unit {unit}.')

        su_type = type(unit).__name__
        if su_type.lower() == 'lagoon':
            su_type = f'{unit.design_type.capitalize()} lagoon'
        name = f'{su_type} {para}'
        model.parameter(setter=AttrSetter(unit, para),
                        name=name, element=unit,
                        kind='coupled', units=df.loc[para]['unit'],
                        baseline=b, distribution=D)


# %%

# =============================================================================
# Shared by all systems
# =============================================================================
su_data_path = ospath.join(data_path, 'sanunit_data/')

def add_shared_parameters(model, water):
    ########## Related to multiple units ##########
    sys = model.system
###########     sys.path[0], sys.path[1] ?????????????
    param = model.parameter
    streams = sys_dct['stream_dct'][sys.ID]
    tea = sys_dct['TEA'][sys.ID]



    ########## Related to raw water ##########

    
    
    # Household size
    b = systems.household_size
    D = shape.Trunc(shape.Normal(mu=b, sigma=1.8), lower=1)
    @param(name='Household size', element=water, kind='coupled', units='cap/household',
            baseline=b, distribution=D)
    def set_household_size(i):
        systems.household_size = i



    ######## General TEA settings ########


    # Money discount rate
    b = systems.discount_rate
    D = shape.Uniform(lower=0.03, upper=0.06)
    @param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
            baseline=b, distribution=D)
    def set_discount_rate(i):
        systems.discount_rate = tea.discount_rate = i

    # Electricity price
    b = price_dct['Electricity']
    D = shape.Triangle(lower=0.08, midpoint=b, upper=0.21)
    @param(name='Electricity price', element='TEA', kind='isolated',
            units='$/kWh', baseline=b, distribution=D)
    def set_electricity_price(i):
        PowerUtility.price = i

    return model


def add_LCA_CF_parameters(model, kind=pou._lca_data.lca_data_kind):
    param = model.parameter
    sys = model.system
    lca = sys_dct['LCA'][sys.ID]

    ######## LCA CF ########
    if kind == 'original':
        b = GWP_dct['CH4']
        D = shape.Uniform(lower=28, upper=34)
        @param(name='CH4 CF', element='LCA', kind='isolated', units='kg CO2-eq/kg CH4',
                baseline=b, distribution=D)
        def set_CH4_CF(i):
            GWP_dct['CH4'] = ImpactItem.get_item('CH4_item').CFs['GlobalWarming'] = i

        b = GWP_dct['N2O']
        D = shape.Uniform(lower=265, upper=298)
        @param(name='N2O CF', element='LCA', kind='isolated', units='kg CO2-eq/kg N2O',
                baseline=b, distribution=D)
        def set_N2O_CF(i):
            GWP_dct['N2O'] = ImpactItem.get_item('N2O_item').CFs['GlobalWarming'] = i

        b = GWP_dct['Electricity']
        D = shape.Uniform(lower=0.106, upper=0.121)
        @param(name='Electricity CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kWh', baseline=b, distribution=D)
        def set_electricity_CF(i):
            GWP_dct['Electricity'] = ImpactItem.get_item('E_item').CFs['GlobalWarming'] = i

        b = GWP_dct['NaClO']
        D = shape.Triangle(lower=1.0*0.75, midpoint=b, upper=1.0*1.25)
        @param(name='NaClO CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kg NaClO', baseline=b, distribution=D)
        def set_NaClO_CF(i):
            GWP_dct['NaClO'] = ImpactItem.get_item('NaClO_item').CFs['GlobalWarming'] = i

        b = GWP_dct['Polyethylene']
        D = shape.Triangle(lower=1.0*0.75, midpoint=b, upper=1.0*1.25)
        @param(name='Polyethylene CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kg Polyethylene', baseline=b, distribution=D)
        def set_Polyethylene_CF(i):
            GWP_dct['Polyethylene'] = ImpactItem.get_item('Polyethylene_item').CFs['GlobalWarming'] = i

       
        item_path = ospath.join(pou._lca_data.data_path, 'items_original.xlsx')
        data = load_data(item_path, sheet='GWP')
        for p in data.index:
            item = ImpactItem.get_item(p)

            b = item.CFs['GlobalWarming']
            lower = float(data.loc[p]['low'])
            upper = float(data.loc[p]['high'])
            dist = data.loc[p]['distribution']
            if dist == 'uniform':
                D = shape.Uniform(lower=lower, upper=upper)
            elif dist == 'triangular':
                D = shape.Triangle(lower=lower, midpoint=b, upper=upper)
            elif dist == 'constant': continue
            else:
                raise ValueError(f'Distribution {dist} not recognized.')
            model.parameter(name=p+'CF',
                            setter=DictAttrSetter(item, 'CFs', 'GlobalWarming'),
                            element='LCA', kind='isolated',
                            units=f'kg CO2-eq/{item.functional_unit}',
                            baseline=b, distribution=D)

       
    return model


# %%

# =============================================================================
# Scenario D (sysD)
# =============================================================================

sysD = systems.sysD
sysD.simulate()
modelD = Model(sysD, add_metrics(sysD, lca_data_kind))
paramD = modelD.parameter

# Shared parameters
modelD = add_shared_parameters(modelD, systems.D1)
modelD = add_LCA_CF_parameters(modelD)

# RawWater
D1 = systems.D1
path = ospath.join(su_data_path, '_raw_water.tsv')
data = load_data(path)
batch_setting_unit_params(data, modelD, D1)
    
    
# Chlorination
D2 = systems.D2
path = ospath.join(su_data_path, '_pou_chlorination.csv')
data = load_data(path)
batch_setting_unit_params(data, modelD, D2)



# # %%

# # =============================================================================
# # Scenario E (sysE)
# # =============================================================================

sysE = systems.sysE
sysE.simulate()
modelE = Model(sysE, add_metrics(sysE, lca_data_kind))
paramE = modelE.parameter

# Shared parameters
modelE = add_shared_parameters(modelE, systems.E1)
modelE = add_LCA_CF_parameters(modelE)

# RawWater
E1 = systems.E1
path = ospath.join(su_data_path, '_raw_water.tsv')
data = load_data(path)
batch_setting_unit_params(data, modelE, E1)
    
    
# AgNP CWF
E2 = systems.E2
path = ospath.join(su_data_path, '_AgNP_CWF_2.csv')
data = load_data(path)
batch_setting_unit_params(data, modelE, E2)
# # %%

# # =============================================================================
# # Functions to run simulation and generate plots
# # =============================================================================

result_dct = {
        'sysD': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        'sysE': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        }

@time_printer
def run_uncertainty(model, seed=None, N=1000, rule='L',
                    percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                    spearman_metrics='default'):
    if seed:
        np.random.seed(seed)

    samples = model.sample(N, rule)

    model.load_samples(samples)
    model.evaluate()

    # Spearman's rank correlation,
    # metrics default to net cost, net emission, and total recoveries
    spearman_results = None
    if spearman_metrics:
        if spearman_metrics.lower() == 'default':
            spearman_metrics = [i for i in model.metrics
                                if 'net' in i.name.lower() or 'total' in i.name.lower()]

        # Different versions of BioSTEAM
        try: spearman_results = model.spearman_r(model.parameters, spearman_metrics)[0]
        except: spearman_results = model.spearman_r(model.parameters, spearman_metrics)

        spearman_results.columns = pd.Index([i.name_with_units for i in spearman_metrics])

    dct = organize_uncertainty_results(model, spearman_results, percentiles)
    return dct


# Data organization
def organize_uncertainty_results(model, spearman_results,
                                  percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)):
    global result_dct
    dct = result_dct[model._system.ID]
    index_p = len(model.parameters)
    dct['parameters'] = model.table.iloc[:, :index_p].copy()
    dct['data'] = model.table.iloc[:, index_p:].copy()

    if percentiles is not None:
        dct['percentiles'] = dct['data'].quantile(q=percentiles)

    if spearman_results is not None:
        dct['spearman'] = spearman_results
    return dct


def save_uncertainty_results(model, dct=None, path=''):
    if not path:
        path = ospath.join(c_path, 'results')

        if not ospath.isdir(path):
            os.mkdir(path)
        path = ospath.join(path, f'sys{model._system.ID[-1]}_model.xlsx')

    elif not (path.endswith('xlsx') or path.endswith('xls')):
        extension = path.split('.')[-1]
        raise ValueError(f'Only "xlsx" and "xls" are supported, not {extension}.')

    dct = dct or result_dct[model._system.ID]
    if dct['parameters'] is None:
        raise ValueError('No cached result, run model first.')
    with pd.ExcelWriter(path) as writer:
        dct['parameters'].to_excel(writer, sheet_name='Parameters')
        dct['data'].to_excel(writer, sheet_name='Uncertainty results')
        if 'percentiles' in dct.keys():
            dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
        dct['spearman'].to_excel(writer, sheet_name='Spearman')
        model.table.to_excel(writer, sheet_name='Raw data')