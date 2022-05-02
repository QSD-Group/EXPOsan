#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txtt
for license details.
'''


# %%

import os, numpy as np, pandas as pd
from chaospy import distributions as shape
from biosteam.evaluation import Metric
from qsdsan import Model, PowerUtility, ImpactItem
from qsdsan.utils import (
    load_data, data_path, AttrSetter, DictAttrSetter, time_printer, dct_from_str
    )
from exposan import reclaimer as re
from exposan.reclaimer import data_path as re_data_path, results_path

__all__ = (
    'create_model', 'result_dct',
    'run_uncertainty', 'save_uncertainty_results',
    )

# %%

# =============================================================================
# Functions for batch-making metrics and -setting parameters
# =============================================================================

systems = re.systems
currency = systems.currency
sys_dct = systems.sys_dct
# unit_dct = systems.unit_dct
price_dct = systems.price_dct
GWP_dct = systems.GWP_dct
GWP = systems.GWP
get_summarizing_functions = systems.get_summarizing_functions
func = get_summarizing_functions()


def add_metrics(system):
    sys_ID = system.ID
    tea = sys_dct['TEA'][sys_ID]
    lca = sys_dct['LCA'][sys_ID]
    ppl = sys_dct['ppl'][sys_ID]
    unit = f'{currency}/cap/yr'
    cat = 'TEA results'
    metrics = [
        Metric('Net cost', lambda: func['get_annual_cost'](tea, ppl), unit, cat),
        Metric('Annual CAPEX', lambda: func['get_annual_CAPEX'](tea, ppl), unit, cat),
        Metric('Energy', lambda: func['get_annual_energy'](tea, ppl), unit, cat),
        Metric('Annual OPEX', lambda: func['get_annual_OPEX'](tea, ppl), unit, cat),
        Metric('Labor', lambda: func['get_annual_labor'](tea, ppl), unit, cat),
        Metric('Sales', lambda: func['get_annual_sales'](tea, ppl), unit, cat),
        ]
    unit = f'{GWP.unit}/cap/yr'
    cat = 'LCA results'
    metrics.extend([
        Metric('Net emission', lambda: func['get_annual_GWP'](lca, ppl), unit, cat),
        Metric('Construction', lambda: func['get_constr_GWP'](lca, ppl), unit, cat),
        Metric('Transportation', lambda: func['get_trans_GWP'](lca, ppl), unit, cat),
        Metric('Fugitive gas', lambda: func['get_CH4_N2O_GWP'](system, lca, ppl), unit, cat),
        Metric('Stream items', lambda: func['get_stream_items_emission_GWP'](system, lca, ppl), unit, cat),
        Metric('Offset', lambda: func['get_offset_GWP'](lca, ppl), unit, cat),
        Metric('Other', lambda: func['get_other_GWP'](lca, ppl), unit, cat),
        ])
    return metrics


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
        model.parameter(setter=AttrSetter(unit, para),
                        name=para, element=unit, kind='coupled', units=df.loc[para]['unit'],
                        baseline=b, distribution=D)


# %%

# =============================================================================
# Pre-load data sheets as they will be used in multiple systems
# =============================================================================

join_path = lambda prefix, file_name: os.path.join(prefix, file_name)
su_data_path = join_path(data_path, 'sanunit_data')
excretion_data = load_data(join_path(su_data_path, '_excretion.tsv'))
murt_toilet_data = load_data(join_path(su_data_path, '_murt_toilet.tsv'))
primary_data = load_data(join_path(su_data_path, '_primary_reclaimer.csv'))
sludge_pasteurization_data = load_data(join_path(su_data_path, '_sludge_pasteurization_reclaimer.tsv'))
ultrafiltration_data = load_data(join_path(su_data_path, '_ultrafiltration_reclaimer.csv'))
ion_exchange_data = load_data(join_path(su_data_path, '_ion_exchange_reclaimer.csv'))
ecr_data = load_data(join_path(su_data_path, '_ECR_reclaimer.csv'))
housing_data = load_data(join_path(su_data_path, '_housing_reclaimer.csv'))
system_data = load_data(join_path(su_data_path, '_system_reclaimer.csv'))
solar_data = load_data(join_path(su_data_path, '_solar_reclaimer.csv'))


# %%

# =============================================================================
# Shared by all systems
# =============================================================================

def add_shared_parameters(sys, model, unit_dct, country_specific=False):
    param = model.parameter
    streams = sys_dct['stream_dct'][sys.ID]

    # Add these parameters if not running country-specific analysis,
    # in which they would be updated separately
    if not country_specific:      

        # Labor wage
        b = price_dct['wages']
        D = shape.Triangle(lower=1.82, midpoint=b, upper=5.46)  # wage in USD/hour
        @param(name='Labor wages', element='TEA', kind='cost', units='USD/h',
               baseline=b, distribution=D)
        def set_labor_wages(i):
            labor_cost = 0
            for u in sys.units:
                if hasattr(u, '_calc_maintenance_labor_cost'):
                    u.wages = i
                    labor_cost += u._calc_maintenance_labor_cost()
            sys.TEA.annual_labor = labor_cost * 365 * 24  # converting labor_cost (USD/hr) to annual_labor (USD/yr)

        # Electricity price
        b = price_dct['Electricity']
        D = shape.Triangle(lower=0.04, midpoint=b, upper=0.1)
        @param(name='Electricity price', element='TEA', kind='isolated',
           units='$/kWh', baseline=b, distribution=D)
        def set_electricity_price(i):
            PowerUtility.price = i

        # Electricity GWP
        b = GWP_dct['Electricity']
        D = shape.Triangle(lower=0.6212, midpoint=b, upper=0.7592)
        @param(name='Electricity CF', element='LCA', kind='isolated',
                   units='kg CO2-eq/kWh', baseline=b, distribution=D)
        def set_electricity_CF(i):
            GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i

    # ########## Specific units ##########
    # Diet and excretion
    excretion_unit = unit_dct['Excretion']
    exclude = ('e_cal', 'p_anim', 'p_veg') if country_specific else ()
    batch_setting_unit_params(excretion_data, model, excretion_unit, exclude)

    # Septic tank (primary)
    primary_unit = unit_dct['Primary']
    batch_setting_unit_params(primary_data, model, primary_unit)

    ##### Universal degradation parameters #####
    # Max methane emission
    unit = sys.path[1] # the first unit that involves degradation
    b = systems.max_CH4_emission
    D = shape.Triangle(lower=0.175, midpoint=b, upper=0.325)
    @param(name='Max CH4 emission', element=unit, kind='coupled', units='g CH4/g COD',
           baseline=b, distribution=D)
    def set_max_CH4_emission(i):
        systems.max_CH4_emission = i

    # Time to full degradation
    b = systems.tau_deg
    D = shape.Uniform(lower=1, upper=3)
    @param(name='Full degradation time', element=unit, kind='coupled', units='yr',
           baseline=b, distribution=D)
    def set_tau_deg(i):
        systems.tau_deg = i
    
    # Reduction at full degradation
    b = systems.log_deg
    D = shape.Uniform(lower=2, upper=4)
    @param(name='Log degradation', element=unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_log_deg(i):
        systems.log_deg = i

    ######## General TEA settings ########
    # # Keeping discount rate constant
    # b = systems.discount_rate
    # D = shape.Uniform(lower=0.03, upper=0.06)
    # @param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
    #         baseline=b, distribution=D)
    # def set_discount_rate(i):
    #     systems.discount_rate = i

    ######## General LCA settings ########
    b = GWP_dct['CH4']
    D = shape.Uniform(lower=28, upper=34)
    @param(name='CH4 CF', element='LCA', kind='isolated', units='kg CO2-eq/kg CH4',
           baseline=b, distribution=D)
    def set_CH4_CF(i):
        GWP_dct['CH4'] = systems.CH4_item.CFs['GlobalWarming'] = i    
    
    b = GWP_dct['N2O']
    D = shape.Uniform(lower=265, upper=298)
    @param(name='N2O CF', element='LCA', kind='isolated', units='kg CO2-eq/kg N2O',
           baseline=b, distribution=D)
    def set_N2O_CF(i):
        GWP_dct['N2O'] = systems.N2O_item.CFs['GlobalWarming'] = i

    item_path = join_path(re_data_path, 'impact_items.xlsx')
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
        model.parameter(name=p,
                        setter=DictAttrSetter(item, 'CFs', 'GlobalWarming'),
                        element='LCA', kind='isolated',
                        units=f'kg CO2-eq/{item.functional_unit}',
                        baseline=b, distribution=D)
    
    return model


# %%

# =============================================================================
# Functions to create models
# =============================================================================

# System A: Solids removal only
def create_modelA(country_specific=False, **model_kwargs):
    sysA = systems.sysA
    modelA = Model(sysA, add_metrics(sysA), **model_kwargs)

    # Shared parameters
    unit_dctA = {
        'Excretion': systems.A1,
        'Primary': systems.A3,
    }
    modelA = add_shared_parameters(sysA, modelA, unit_dctA, country_specific)

    # Sludge pasteurization
    sludge_pasteurization_data = load_data(join_path(su_data_path, '_sludge_pasteurization_reclaimer.tsv'))
    exclude = 'wages' if country_specific else ()
    batch_setting_unit_params(sludge_pasteurization_data, modelA, systems.A4, exclude)

    return modelA


# System B: Full Duke system with grid electricity source
def create_modelB(country_specific=False, **model_kwargs):
    sysB = systems.sysB
    modelB = Model(sysB, add_metrics(sysB), **model_kwargs)
    paramB = modelB.parameter

    # Shared parameters
    unit_dctB = {
        'Excretion': systems.B1,
        'Primary': systems.B3,
    }
    modelB = add_shared_parameters(sysB, modelB, unit_dctB, country_specific)

    # MURT toilet
    B2 = systems.B2
    murt_toilet_data = load_data(join_path(su_data_path, '_murt_toilet.tsv'))
    exclude = ('MCF_decay', 'N2O_EF_decay', 'OPEX_over_CAPEX')
    batch_setting_unit_params(murt_toilet_data, modelB, systems.B2, exclude)

    b = B2.OPEX_over_CAPEX
    D = shape.Uniform(lower=0.02, upper=0.08)
    @paramB(name='MURT Toilet operating cost', element=B2, kind='coupled', units='cost',
            baseline=b, distribution=D)
    def set_OPEX_over_CAPEX(i):
        B2.OPEX_over_CAPEX = i

    b = B2.MCF_decay
    D = shape.Triangle(lower=0.05, midpoint=b, upper=0.15)
    @paramB(name='MCF_decay', element=B2, kind='coupled',
            units='fraction of anaerobic conversion of degraded COD',
            baseline=b, distribution=D)
    def set_MCF_decay(i):
        B2.MCF_decay = i

    b = B2.N2O_EF_decay
    D = shape.Triangle(lower=0, midpoint=b, upper=0.001)
    @paramB(name='N2O_EF_decay', element=B2, kind='coupled',
            units='fraction of N emitted as N2O',
            baseline=b, distribution=D)
    def set_N2O_EF_decay(i):
        B2.N2O_EF_decay = i

    # Sludge pasteurization
    sludge_pasteurization_data = load_data(join_path(su_data_path, '_sludge_pasteurization_reclaimer.tsv'))
    exclude = 'wages' if country_specific else ()
    batch_setting_unit_params(sludge_pasteurization_data, modelB, systems.B4, exclude)

    # Ultrafiltration
    ultrafiltration_data = load_data(join_path(su_data_path, '_ultrafiltration_reclaimer.csv'))
    batch_setting_unit_params(ultrafiltration_data, modelB, systems.B5)

    # Ion exchange
    ion_exchange_data = load_data(join_path(su_data_path, '_ion_exchange_reclaimer.csv'))
    exclude = 'wages' if country_specific else ()
    batch_setting_unit_params(ion_exchange_data, modelB, systems.B6, exclude)

    # ECR
    ecr_data = load_data(join_path(su_data_path, '_ECR_reclaimer.csv'))
    batch_setting_unit_params(ecr_data, modelB, systems.B7)

    # Housing
    housing_data = load_data(join_path(su_data_path, '_housing_reclaimer.csv'))
    batch_setting_unit_params(housing_data, modelB, systems.B10)

    # System
    system_data = load_data(join_path(su_data_path, '_system_reclaimer.csv'))
    batch_setting_unit_params(system_data, modelB, systems.B11)

    return modelB


# System C: Full Duke system with solar electricity source
def create_modelC(country_specific=False, **model_kwargs):
    sysC = systems.sysC
    modelC = Model(sysC, add_metrics(sysC), **model_kwargs)
    paramC = modelC.parameter

    # Shared parameters
    unit_dctC = {
        'Excretion': systems.C1,
        'Primary': systems.C3,
    }
    modelC = add_shared_parameters(sysC, modelC, unit_dctC, country_specific)

    # MURT toilet
    C2 = systems.C2
    murt_toilet_data = load_data(join_path(su_data_path, '_murt_toilet.tsv'))
    exclude = ('MCF_decay', 'N2O_EF_decay', 'OPEX_over_CAPEX')
    batch_setting_unit_params(murt_toilet_data, modelC, systems.C2, exclude)

    b = C2.OPEX_over_CAPEX
    D = shape.Uniform(lower=0.02, upper=0.08)
    @paramC(name='MURT Toilet operating cost', element=C2, kind='coupled', units='cost',
            baseline=b, distribution=D)
    def set_OPEX_over_CAPEX(i):
        C2.OPEX_over_CAPEX = i

    b = C2.MCF_decay
    D = shape.Triangle(lower=0.05, midpoint=b, upper=0.15)
    @paramC(name='MCF_decay', element=C2, kind='coupled',
            units='fraction of anaerobic conversion of degraded COD',
            baseline=b, distribution=D)
    def set_MCF_decay(i):
        C2.MCF_decay = i

    b = C2.N2O_EF_decay
    D = shape.Triangle(lower=0, midpoint=b, upper=0.001)
    @paramC(name='N2O_EF_decay', element=C2, kind='coupled',
            units='fraction of N emitted as N2O',
            baseline=b, distribution=D)
    def set_N2O_EF_decay(i):
        C2.N2O_EF_decay = i

    # Sludge pasteurization
    sludge_pasteurization_data = load_data(join_path(su_data_path, '_sludge_pasteurization_reclaimer.tsv'))
    exclude = 'wages' if country_specific else ()
    batch_setting_unit_params(sludge_pasteurization_data, modelC, systems.C4, exclude)

    # Ultrafiltration
    ultrafiltration_data = load_data(join_path(su_data_path, '_ultrafiltration_reclaimer.csv'))
    batch_setting_unit_params(ultrafiltration_data, modelC, systems.C5)

    # Ion exchange
    ion_exchange_data = load_data(join_path(su_data_path, '_ion_exchange_reclaimer.csv'))
    exclude = 'wages' if country_specific else ()
    batch_setting_unit_params(ion_exchange_data, modelC, systems.C6, exclude)

    # ECR
    ecr_data = load_data(join_path(su_data_path, '_ECR_reclaimer.csv'))
    batch_setting_unit_params(ecr_data, modelC, systems.C7)

    # Housing
    housing_data = load_data(join_path(su_data_path, '_housing_reclaimer.csv'))
    batch_setting_unit_params(housing_data, modelC, systems.C10)

    # System
    system_data = load_data(join_path(su_data_path, '_system_reclaimer.csv'))
    batch_setting_unit_params(system_data, modelC, systems.C11)

    # Solar
    solar_data = load_data(join_path(su_data_path, '_solar_reclaimer.csv'))
    exclude = 'wages' if country_specific else ()
    batch_setting_unit_params(solar_data, modelC, systems.C12, exclude)

    return modelC


# System D: Targeted nitrogen removal (created for NSS preliminary analysis)
def create_modelD(country_specific=False, **model_kwargs):
    sysD = systems.sysD
    modelD = Model(sysD, add_metrics(sysD), **model_kwargs)
    paramD = modelD.parameter

    # Shared parameters
    unit_dctD = {
        'Excretion': systems.D1,
        'Primary': systems.D3,
    }
    modelD = add_shared_parameters(sysD, modelD, unit_dctD, country_specific)

    # Ultrafiltration
    ultrafiltration_data = load_data(join_path(su_data_path, '_ultrafiltration_reclaimer.csv'))
    batch_setting_unit_params(ultrafiltration_data, modelD, systems.D4)

    # Ion exchange
    ion_exchange_data = load_data(join_path(su_data_path, '_ion_exchange_reclaimer.csv'))
    batch_setting_unit_params(ion_exchange_data, modelD, systems.D5)

    # Housing
    housing_data = load_data(join_path(su_data_path, '_housing_reclaimer.csv'))
    batch_setting_unit_params(housing_data, modelD, systems.D8)

    # System
    system_data = load_data(join_path(su_data_path, '_system_reclaimer.csv'))
    batch_setting_unit_params(system_data, modelD, systems.D9)

    return modelD


# Wrapper function so that it'd work for all
def create_model(model_ID='A', country_specific=False, **model_kwargs):
    model_ID = model_ID.lstrip('model').lstrip('sys') # so that it'll work for "modelA"/"sysA"/"A"
    if model_ID == 'A': model = create_modelA(country_specific, **model_kwargs)
    elif model_ID == 'B': model = create_modelB(country_specific, **model_kwargs)
    elif model_ID == 'C': model = create_modelC(country_specific, **model_kwargs)
    else: model = create_modelD(country_specific, **model_kwargs)
    return model


# %%

# =============================================================================
# Functions to run simulation and generate plots
# =============================================================================

result_dct = {
        'sysA': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        'sysB': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        'sysC': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        'sysD':  dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman'))
        }

@time_printer
def run_uncertainty(model, seed=None, N=10000, rule='L',
                    percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                    print_time=False):
    global result_dct
    if seed:
        np.random.seed(seed)

    samples = model.sample(N, rule)
    model.load_samples(samples)
    model.evaluate()

    # Data organization
    dct = result_dct[model.system.ID]
    index_p = len(model.get_parameters())
    dct['parameters'] = model.table.iloc[:, :index_p].copy()
    dct['data'] = model.table.iloc[:, index_p:].copy()
    if percentiles:
        dct['percentiles'] = dct['data'].quantile(q=percentiles)
        dct['percentiles_parameters'] = dct['parameters'].quantile(q=percentiles)

    # Spearman's rank correlation
    spearman_metrics = model.metrics[:13]
    spearman_results = model.spearman(model.get_parameters(), spearman_metrics)
    spearman_results.columns = pd.Index([i.name_with_units for i in spearman_metrics])
    dct['spearman'] = spearman_results
    return dct


def save_uncertainty_results(model, dct={}, path=''):
    sys_ID = model.system.ID
    population = systems.ppl
    path = join_path(results_path, f'uncertainty{sys_ID[-1]}_{population}users.xlsx') if path=='' else path
    dct = dct or result_dct[sys_ID]
    if dct['parameters'] is None:
        raise ValueError('No cached result, run model first.')
    with pd.ExcelWriter(path) as writer:
        dct['parameters'].to_excel(writer, sheet_name='Parameters')
        dct['data'].to_excel(writer, sheet_name='Uncertainty results')
        if 'percentiles' in dct.keys():
            dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
        dct['spearman'].to_excel(writer, sheet_name='Spearman')
        model.table.to_excel(writer, sheet_name='Raw data')
