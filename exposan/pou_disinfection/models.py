#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Bright Elijah <be05055@georgiasouthern.edu & brightcarlelijah@gmail.com>
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import os, qsdsan as qs
from chaospy import distributions as shape
from qsdsan import Model, Metric, PowerUtility, ImpactItem
from qsdsan.utils import load_data
from exposan.utils import (
    batch_setting_LCA_params,
    batch_setting_unit_params,
    run_uncertainty as run
    )
from exposan import pou_disinfection as pou
from exposan.pou_disinfection import (
    create_system,
    data_path,
    get_LCA_metrics,
    get_TEA_metrics,
    results_path,
    update_number_of_householdsize,
    )

__all__ = ('create_model', 'run_uncertainty',)


# %%

# =============================================================================
# Functions for batch-making metrics and -setting parameters
# =============================================================================

def add_metrics(model, ppl=pou.ppl, include_breakdown=False):
    pou._load_lca_data()
    system = model.system
    TEA_functions = get_TEA_metrics(system, ppl=ppl, include_breakdown=include_breakdown)
    LCA_functions = get_LCA_metrics(system, ppl=ppl, include_breakdown=include_breakdown)
    indicators = system.LCA.indicators
    
    unit = f'{qs.currency}/cap/yr'
    metrics = [
        Metric('Total cost', TEA_functions[0], unit, 'TEA results'),
        ]
    
    if include_breakdown:
        metrics.extend([
            Metric('CAPEX', TEA_functions[1], unit, 'TEA results'),
            Metric('OPEX', TEA_functions[2], unit, 'TEA results'),
            Metric('Electricity', TEA_functions[3], unit, 'TEA results'),
            Metric('Others', TEA_functions[4], unit, 'TEA results'),
            ])
    else:
        for n, ind in enumerate(indicators):
            ID, unit = ind.ID, ind.unit+'/cap/year'
            metrics.append(
                Metric(f'Total {ID}', LCA_functions[n-1], unit, 'LCA results'),
                )
        model.metrics = metrics
        return

    n = len(indicators)
    if include_breakdown:
        for ind in indicators:
            ID, unit = ind.ID, ind.unit+'/cap/year'
            metrics.extend([
                Metric(f'Construction {ID}', LCA_functions[n], unit, 'LCA results'),
                Metric(f'Transportation {ID}', LCA_functions[n+1], unit, 'LCA results'),
                Metric(f'Streams {ID}', LCA_functions[n+2], unit, 'LCA results'),
                Metric(f'Others {ID}', LCA_functions[n+3], unit, 'LCA results'),
                ])
            n += 4
    model.metrics = metrics


def add_shared_parameters(model, ppl=pou.ppl):
    sys = model.system
    param = model.parameter

    ########## Related to multiple units ##########
    # Household size, Stetson compiled 10/19/2023
    # https://www.un.org/development/desa/pd/data/household-size-and-composition
    raw_water = sys.path[0]
    b = pou.household_size
    # This use of sigma should be correct since we are not using it to 
    # calculate the number of unit needed in a community
    # More details: https://github.com/QSD-Group/EXPOsan/issues/39
    D = shape.Trunc(shape.Normal(mu=b, sigma=1.4), lower=1)
    @param(name='Household size', element=raw_water, kind='coupled', units='cap/household',
            baseline=b, distribution=D)
    def set_household_size(i):
        pou.household_size = max(1, i)
        update_number_of_householdsize(sys, household_size=i, ppl=ppl)


    ######## General TEA and LCA settings ########
    # Money discount rate
    tea = sys.TEA
    b = pou.discount_rate
    D = shape.Uniform(lower=0.03, upper=0.06)
    @param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
           baseline=b, distribution=D)
    def set_discount_rate(i):
        pou.discount_rate = tea.discount_rate = i

    # Electricity price
    b = PowerUtility.price
    D = shape.Triangle(lower=0.08, midpoint=b, upper=0.21)
    @param(name='Electricity price', element='TEA', kind='isolated',
            units='$/kWh', baseline=b, distribution=D)
    def set_electricity_price(i):
        PowerUtility.price = i

    # CF of PE_stream should be the same as PE,
    # will be added separately for POU chlorination (the only system that uses it)
    item_path = os.path.join(pou.data_path, 'impact_items.xlsx')
    batch_setting_LCA_params(item_path, model, exclude=('PE_stream'))


def import_water_data(water_source):
    water_path = os.path.join(data_path, '_raw_water.xlsx')
    sheet = 'water1_gw' if water_source.lower() in ('gw', 'groundwater') \
        else 'water2_sw'
    return load_data(water_path, sheet=sheet)


# %%

# =============================================================================
# Functions to create and run models
# =============================================================================

# System A: POU Chlorination
def create_modelA(system_kwargs={}, **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysA = create_system('A', flowsheet=flowsheet, **system_kwargs)
    unitA = sysA.flowsheet.unit
    
    # Shared metrics/parameters
    modelA = Model(sysA, **model_kwargs)
    ppl = system_kwargs.get('ppl', pou.ppl)
    add_metrics(modelA, ppl=ppl)
    add_shared_parameters(modelA, ppl=ppl)
    
    # NaClO price
    A_naclo = sysA.flowsheet.stream['A_naclo']
    b = A_naclo.price
    D = shape.Uniform(lower=b*0.75, upper=b*1.25)
    @modelA.parameter(
        name='NaClO price', element='TEA', kind='isolated',
        units='$/kg', baseline=b, distribution=D,
        )
    def set_NaClO_price(i):
        A_naclo.price = i
    
    # Since the chlorine bottle also uses PE as the container,
    # adjust so they will be using the same value
    A_naclo = sysA.flowsheet.stream['A_cl_bottle']
    PE, PE_stream = ImpactItem.get_item('PE'), ImpactItem.get_item('PE_stream')
    # The distribution of this dummy parameter actually doesn't matter
    # (since it's not actually used anywhere)
    b = PE_stream.CFs['GWP']
    D = shape.Uniform(lower=b*0.75, upper=b*1.25)
    @modelA.parameter(
        name='PE_stream_dummy', element='LCA', kind='isolated',
        units='kg CO2-eq/kg', baseline=b, distribution=D,
        )
    def set_PE_stream_CF(i):
        PE_stream.CFs['GWP'] = PE.CFs['GWP']
    
    # RawWater
    water_data = import_water_data(system_kwargs.get('water_source', pou.water_source))
    batch_setting_unit_params(water_data, modelA, unitA.A1)

    # Chlorination
    chlorination_data_path = os.path.join(data_path, '_pou_chlorination.csv')
    chlorination_data = load_data(chlorination_data_path)
    batch_setting_unit_params(chlorination_data, modelA, unitA.A2)

    return modelA


# System B: AgNP CWF
def create_modelB(system_kwargs={}, **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysB = create_system('B', flowsheet=flowsheet, **system_kwargs)
    unitB = sysB.flowsheet.unit
    
    # Shared metrics/parameters
    modelB = Model(sysB, **model_kwargs)
    ppl = system_kwargs.get('ppl', pou.ppl)
    add_metrics(modelB, ppl=ppl)
    add_shared_parameters(modelB, ppl=ppl)
    
    # RawWater
    water_data = import_water_data(system_kwargs.get('water_source', pou.water_source))
    batch_setting_unit_params(water_data, modelB, unitB.B1)
    
    # AgNP CWF
    agnp_cwf_path = os.path.join(data_path, '_agnp_cwf.csv')
    agnp_cwf_data = load_data(agnp_cwf_path)
    batch_setting_unit_params(agnp_cwf_data, modelB, unitB.B2)

    return modelB

# System C: POU UV
def create_modelC(system_kwargs={}, **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysC = create_system('C', flowsheet=flowsheet, **system_kwargs)
    unitC = sysC.flowsheet.unit
    
    # Shared metrics/parameters
    modelC = Model(sysC, **model_kwargs)
    ppl = system_kwargs.get('ppl', pou.ppl)
    add_metrics(modelC, ppl=ppl)
    add_shared_parameters(modelC, ppl=ppl)
    
    # RawWater
    water_data = import_water_data(system_kwargs.get('water_source', pou.water_source))
    batch_setting_unit_params(water_data, modelC, unitC.C1)
    
    # UV lamp
    pou_uv_path = os.path.join(data_path, '_pou_uv.csv')
    pou_uv_data = load_data(pou_uv_path)
    batch_setting_unit_params(pou_uv_data, modelC, unitC.C2)

    return modelC


# System D: UV LED
def create_modelD(system_kwargs={}, **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysD = create_system('D', flowsheet=flowsheet, **system_kwargs)
    unitD = sysD.flowsheet.unit
    
    # Shared metrics/parameters
    modelD = Model(sysD, **model_kwargs)
    ppl = system_kwargs.get('ppl', pou.ppl)
    add_metrics(modelD, ppl=ppl)
    add_shared_parameters(modelD, ppl=ppl)
    
    # RawWater
    water_data = import_water_data(system_kwargs.get('water_source', pou.water_source))
    batch_setting_unit_params(water_data, modelD, unitD.D1)
    
    # UV LED
    uv_led_path = os.path.join(data_path, '_uv_led.csv')
    uv_led_data = load_data(uv_led_path)
    batch_setting_unit_params(uv_led_data, modelD, unitD.D2)

    return modelD


# Wrapper function so that it'd work for all
def create_model(model_ID='A',
                 water_source=pou.water_source,
                 household_size=pou.household_size,
                 ppl=pou.ppl,
                 **model_kwargs):
    model_ID = model_ID.lower().rsplit('model')[-1].rsplit('sys')[-1].upper() # works for "modelA"/"sysA"/"A"
    system_kwargs = {
        'water_source': water_source,
        'household_size': household_size,
        'ppl': ppl,
        }
    if model_ID == 'A': model = create_modelA(system_kwargs=system_kwargs, **model_kwargs)
    elif model_ID == 'B': model = create_modelB(system_kwargs=system_kwargs, **model_kwargs)
    elif model_ID == 'C': model = create_modelC(system_kwargs=system_kwargs, **model_kwargs)
    elif model_ID == 'D': model = create_modelD(system_kwargs=system_kwargs, **model_kwargs)
    else: raise ValueError(f'`model_ID` can only be "A", "B", "C", or "D", not "{model_ID}".')
    return model


def run_uncertainty(model, path='', **kwargs):
    kwargs['path'] = os.path.join(results_path, f'sys{model.system.ID[-1]}_model.xlsx') if path=='' else path
    run(model=model, **kwargs)
    return
