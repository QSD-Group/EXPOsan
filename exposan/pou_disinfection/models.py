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
from qsdsan.utils import (
    DictAttrSetter,
    load_data,
    )
from exposan.utils import batch_setting_unit_params, run_uncertainty as run
from exposan import pou_disinfection as pou
from exposan.pou_disinfection import (
    create_system,
    data_path,
    get_LCA_metrics,
    get_TEA_metrics,
    GWP_dct,
    price_dct,
    results_path,
    )

__all__ = ('create_model', 'run_uncertainty',)


# %%

# =============================================================================
# Functions for batch-making metrics and -setting parameters
# =============================================================================

def add_metrics(model, include_breakdown=False):
    pou._load_lca_data()
    system = model.system
    TEA_functions = get_TEA_metrics(system, include_breakdown=include_breakdown)
    LCA_functions = get_LCA_metrics(system, include_breakdown=include_breakdown)
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
                Metric('Total {ID}', LCA_functions[n-1], unit, 'LCA results'),
                )
        return metrics

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

    return metrics


def add_shared_parameters(model, water):
    sys = model.system
    sys_stream = sys.flowsheet.stream
    param = model.parameter

    ########## Related to multiple units ##########
    # Household size
    raw_water = sys.path[0]
    b = pou.household_size
    #!!! Need to update this sigma
    D = shape.Trunc(shape.Normal(mu=b, sigma=1.8), lower=1)
    @param(name='Household size', element=raw_water, kind='coupled', units='cap/household',
            baseline=b, distribution=D)
    def set_household_size(i):
        pou.household_size = max(1, i)

    ######## General TEA settings ########
    # Money discount rate
    tea = sys.TEA
    b = pou.discount_rate
    D = shape.Uniform(lower=0.03, upper=0.06)
    @param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
            baseline=b, distribution=D)
    def set_discount_rate(i):
        pou.discount_rate = tea.discount_rate = i

    # Electricity price
    b = price_dct['Electricity']
    D = shape.Triangle(lower=0.08, midpoint=b, upper=0.21)
    @param(name='Electricity price', element='TEA', kind='isolated',
            units='$/kWh', baseline=b, distribution=D)
    def set_electricity_price(i):
        PowerUtility.price = i
        
    # NaClO price
    b = price_dct['NaClO']
    D = shape.Uniform(lower=b*0.75, upper=b*1.25)
    @param(name='NaClO price', element='TEA', kind='isolated',
            units='$/kg', baseline=b, distribution=D)
    def set_NaClO_price(i):
        price_dct['NaClO'] = sys_stream['NaClO'].price = i
        
    ######## General LCA settings ########
    lca = sys.LCA
    #!!! Why not add electricity, etc., through spreadsheet?
    b = GWP_dct['Electricity']
    D = shape.Uniform(lower=0.106, upper=0.121)
    @param(name='Electricity CF', element='LCA', kind='isolated',
            units='kg CO2-eq/kWh', baseline=b, distribution=D)
    def set_electricity_CF(i):
        GWP_dct['Electricity'] = ImpactItem.get_item('E_item').CFs['GlobalWarming'] = i

    b = GWP_dct['NaClO']
    D = shape.Triangle(lower=2.6287*0.75, midpoint=b, upper=2.6287*1.25)
    @param(name='NaClO CF', element='LCA', kind='isolated',
            units='kg CO2-eq/kg NaClO', baseline=b, distribution=D)
    def set_NaClO_CF(i):
        GWP_dct['NaClO'] = ImpactItem.get_item('NaClO_item').CFs['GlobalWarming'] = i

    b = GWP_dct['Polyethylene']
    D = shape.Triangle(lower=2.7933*0.75, midpoint=b, upper=2.7933*1.25)
    @param(name='Polyethylene CF', element='LCA', kind='isolated',
            units='kg CO2-eq/kg Polyethylene', baseline=b, distribution=D)
    def set_Polyethylene_CF(i):
        GWP_dct['Polyethylene'] = ImpactItem.get_item('Polyethylene_item').CFs['GlobalWarming'] = i

    # b = GWP_dct['PVC']
    # D = shape.Triangle(lower=1.0*0.75, midpoint=b, upper=1.0*1.25)
    # @param(name='PVC CF', element='LCA', kind='isolated',
    #         units='kg CO2-eq/kg PVC', baseline=b, distribution=D)
    # def set_PVC_CF(i):
    #     GWP_dct['PVC'] = ImpactItem.get_item('PVC_item').CFs['GlobalWarming'] = i
        
    # b = GWP_dct['Mecury']
    # D = shape.Triangle(lower=1.0*0.75, midpoint=b, upper=1.0*1.25)
    # @param(name='Mecury CF', element='LCA', kind='isolated',
    #         units='kg CO2-eq/kg Mecury', baseline=b, distribution=D)
    # def set_Mecury_CF(i):
    #     GWP_dct['Mecury'] = ImpactItem.get_item('Mecury_item').CFs['GlobalWarming'] = i
        
    # b = GWP_dct['Aluminum']
    # D = shape.Triangle(lower=1.0*0.75, midpoint=b, upper=1.0*1.25)
    # @param(name='Aluminum CF', element='LCA', kind='isolated',
    #         units='kg CO2-eq/kg Mecury', baseline=b, distribution=D)
    # def set_Aluminum_CF(i):
    #     GWP_dct['Aluminum'] = ImpactItem.get_item('Aluminum_item').CFs['GlobalWarming'] = i
        
    item_path = os.path.join(pou.data_path, 'impact_items.xlsx')
    for ind in lca.indicators:
        data = load_data(item_path, sheet=ind.ID)
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


# =============================================================================
# Datasheets
# =============================================================================

# Groundwater
gw_path = os.path.join(data_path, '_raw_water1_gw.tsv')
gw_data = load_data(gw_path)
# Surface water
sw_path = os.path.join(data_path, '_raw_water2_sw.tsv')
sw_data = load_data(sw_path)
chlorination_data_path = os.path.join(data_path, '_pou_chlorination.csv')
chlorination_data = load_data(chlorination_data_path)
cwf_path = os.path.join(data_path, '_AgNP_CWF_2.csv')
cwf_data = load_data(cwf_path)
pou_uv_path = os.path.join(data_path, '_pou_uv.csv')
pou_uv_data = load_data(pou_uv_path)
uv_led_path = os.path.join(data_path, '_uv_led.csv')
uv_led_data = load_data(uv_led_path)


# %%

# =============================================================================
# Functions to create and run models
# =============================================================================

# System A: POU Chlorination
def create_modelA(**model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysA = create_system('A', flowsheet=flowsheet)
    unitA = sysA.flowsheet.unit
    
    # Shared metrics/parameters
    modelA = Model(sysA, **model_kwargs)
    add_metrics(modelA)
    add_shared_parameters(modelA)
    
    # RawWater
    batch_setting_unit_params(gw_data, modelA, unitA.A1)

    # Chlorination
    batch_setting_unit_params(chlorination_data, modelA, unitA.A2)

    return modelA


# System B: AgNP CWF
def create_modelB(**model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysB = create_system('B', flowsheet=flowsheet)
    unitB = sysB.flowsheet.unit
    
    # Shared metrics/parameters
    modelB = Model(sysB, **model_kwargs)
    add_metrics(modelB)
    add_shared_parameters(modelB)
    
    # RawWater
    batch_setting_unit_params(gw_data, modelB, unitB.B1)
    
    # AgNP CWF
    batch_setting_unit_params(cwf_data, modelB, unitB.B2)

    return modelB

# System C: POU UV
def create_modelC(**model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysC = create_system('C', flowsheet=flowsheet)
    unitC = sysC.flowsheet.unit
    
    # Shared metrics/parameters
    modelC = Model(sysC, **model_kwargs)
    add_metrics(modelC)
    add_shared_parameters(modelC)
    
    # RawWater
    batch_setting_unit_params(gw_data, modelC, unitC.C1)
    
    # UV lamp
    batch_setting_unit_params(pou_uv_data, modelC, unitC.C2)

    return modelC


# System D: UV LED
def create_modelD(**model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysD = create_system('D', flowsheet=flowsheet)
    unitD = sysD.flowsheet.unit
    
    # Shared metrics/parameters
    modelD = Model(sysD, **model_kwargs)
    add_metrics(modelD)
    add_shared_parameters(modelD)
    
    # RawWater
    batch_setting_unit_params(gw_data, modelD, unitD.D1)
    
    # UV LED
    batch_setting_unit_params(uv_led_data, modelD, unitD.D2)

    return modelD


# Wrapper function so that it'd work for all
def create_model(model_ID='A', country_specific=False, **model_kwargs):
    model_ID = model_ID.lower().rsplit('model')[-1].rsplit('sys')[-1].upper() # works for "modelA"/"sysA"/"A"
    if model_ID == 'A': model = create_modelA(country_specific, **model_kwargs)
    elif model_ID == 'B': model = create_modelB(country_specific, **model_kwargs)
    elif model_ID == 'C': model = create_modelC(country_specific, **model_kwargs)
    elif model_ID == 'D': model = create_modelD(country_specific, **model_kwargs)
    else: raise ValueError(f'`model_ID` can only be "A", "B", "C", or "D", not "{model_ID}".')
    return model


def run_uncertainty(model, path='', **kwargs):
    kwargs['path'] = os.path.join(results_path, f'sys{model.system.ID[-1]}_model.xlsx') if path=='' else path
    run(model=model, **kwargs)
    return
