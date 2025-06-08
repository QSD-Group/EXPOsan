#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Lewis Rowles <stetsonsc@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

    Hannah Lohman <hlohman94@gmail.com>

    Lane To <lane20@illinois.edu>
    
    Aaron Marszewski <aaronpm3@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import os, pandas as pd, qsdsan as qs
from chaospy import distributions as shape
from qsdsan import Model, Metric, PowerUtility, ImpactItem, TEA
from qsdsan.utils import (
    AttrSetter,
    data_path,
    DictAttrSetter,
    dct_from_str,
    load_data,
    )
from exposan.utils import batch_setting_unit_params, run_uncertainty as run
from exposan import biogenic_refinery_context as br
from exposan.biogenic_refinery_context import (
    create_system,
    data_path, #as br_data_path,
    br_path,
    results_path,
    #get_decay_k, not found in _init, commented for now #TODO
    get_LCA_metrics,
    get_TEA_metrics,
    get_sustainability_indicators,
    results_path,
    update_resource_recovery_settings,
    )

from . import (
    cost_per_ton_biochar,
    gwp_per_ton_biochar,
    biochar_generated,
    sequesterable_carbon,
    drying_requirement,
    drying_cost,
    energy_conversion_efficiency,
    char_yield,
    ash_content_feedstock,
    moisture_content_feedstock,
    biosolids_HHV
)

__all__ = ('create_model', 'run_uncertainty',)


# %%

# =============================================================================
# Functions for batch-making metrics and -setting parameters
# =============================================================================

def add_metrics(model):
    br._load_lca_data()
    system = model.system
    # Sustainability indicators
    ##TODO need to define this function in _init_ and add other indicators of interest here. 
    funcs = get_sustainability_indicators(system)
    
    # model.TEA = TEA(
    # system=model.system,
    # discount_rate = 0.05,
    # lifetime = 20,
# )
    
    metrics = [
        Metric('Cost per ton biochar', lambda: cost_per_ton_biochar(model), 'USD/ton biochar'),
        Metric('GWP per ton biochar', lambda: gwp_per_ton_biochar(model), 'kg CO2-eq/ton biochar'),
        Metric('Biochar generated', lambda: biochar_generated(model), 'ton biochar/yr'),
        Metric('Sequesterable carbon', lambda: sequesterable_carbon(model), 'ton C/yr'),
        Metric('Drying requirement', lambda: drying_requirement(model), 'MJ/ton biosolids'),
        Metric('Drying cost normalized', lambda: drying_cost(model), 'USD/ton biochar'),
        Metric('Energy Conversion Efficiency', lambda: energy_conversion_efficiency(model), '%'),
        Metric('Char Yield (db %)', lambda: char_yield(model), '%'),
        Metric('Ash Content Feedstock (db %)', lambda: ash_content_feedstock(model), '%'),
        Metric('Moisture Content Feedstock', lambda: moisture_content_feedstock(model), '%'),
        Metric('Biosolids HHV', lambda: biosolids_HHV(model), 'MJ/kg'),
    ]

    # Net cost
    # metrics.append(
    #     Metric('Annual net cost', get_TEA_metrics(system)[0], f'{qs.currency}/cap/yr', 'TEA results'),
    #     )
    # # Net emissions
    # funcs = get_LCA_metrics(system)
    # cat = 'LCA results'
    # metrics.extend([
    #     Metric('GlobalWarming', funcs[0], 'kg CO2-eq/cap/yr', cat),
    #     ])
    model.metrics = metrics


# %%

# =============================================================================
# Data sheets
# =============================================================================

su_data_path = os.path.join(data_path)
br_su_data_path = os.path.join(su_data_path, 'br')

def load_br_su_data(file_name):
    if file_name.startswith('_br'):
        return load_data(os.path.join(br_su_data_path, file_name))
    return load_data(os.path.join(su_data_path, file_name))

biosolids_data = load_br_su_data('_br_biosolids.tsv')
controls_data = load_br_su_data('_br_controls.tsv')
housing_data = load_br_su_data('_br_housing.tsv')
carbonizer_data = load_br_su_data('_br_carbonizer_base.tsv')
pollution_control_data = load_br_su_data('_br_pollution_control.tsv')
ohx_data = load_br_su_data('_br_ohx.tsv')
hhx_data = load_br_su_data('_br_hhx.tsv')
hhx_dryer_data = load_br_su_data('_br_hhx_dryer.tsv')



# %%

# =============================================================================
# Shared by systems A-C
# =============================================================================

def add_shared_parameters(model, unit_dct, location_specific=False):
    sys = model.system
    sys_stream = sys.flowsheet.stream
    param = model.parameter
    price_dct, GWP_dct = update_resource_recovery_settings()

    # Add these parameters if not running country-specific analysis,
    # in which they would be updated separately
    if not location_specific:
        # Price ratio
        # Just want to have this parameter so that can be used in other analyses,
        # set the distribution to be a really tight one
        old_price_dct = price_dct.copy()

        b = 1
        D = shape.Uniform(lower=b-(10**(-6)), upper=b+(10**(-6)))
        item_ref = {
            'Concrete': 'Concrete',
            'Steel': 'Steel',
        }
        stream_ref = {
            'biochar': 'biochar',
        }

        def set_price_ratio(i):
            br.price_ratio = i
            for obj_name in (*item_ref.keys(), *stream_ref.keys()):
                old_price = old_price_dct[obj_name]
                new_price = old_price * i
                if obj_name in item_ref.keys():
                    ImpactItem.get_item(item_ref[obj_name]).price = new_price
                else:
                    getattr(sys_stream, stream_ref[obj_name]).price = new_price
            for u in sys.units:
                if hasattr(u, 'price_ratio'):
                    u.price_ratio = i


        # Operator labor wage TODO: Scale to scale_factor?
        b = br.operator_daily_wage
        D = shape.Triangle(lower=(14.55), midpoint=b, upper=(43.68))
        @param(name='Operator daily wages', element='TEA', kind='cost', units='USD/d',
              baseline=b, distribution=D)
        def set_operator_daily_wage(i):
            sys._TEA.annual_labor = i*3*365

        # Construction labor wage
        b = br.const_daily_wage
        D = shape.Triangle(lower=(b*0.5), midpoint=b, upper=(b*1.5))
        @param(name='Construction daily wages', element='TEA', kind='cost', units='USD/d',
              baseline=b, distribution=D)
        def set_const_daily_wage(i):
            for u in sys.units:
                if isinstance(u, qs.sanunits.BiogenicRefineryHousing): break
                u.const_daily_wage = i
        
        # if br.INCLUDE_RESOURCE_RECOVERY:
        #     # Commented out because not taking into account economic value of biochar
        #     # D = shape.Uniform(lower=(0.014*0.8), upper=(0.014*1.2))
        #     # @param(name='Biochar  price', element='TEA', kind='isolated', units='USD/kg biochar',
        #     #         baseline=(0.014), distribution=D)
        #     # def set_biochar_price(i):
        #     #     price_dct['biochar'] = sys_stream.biochar.price = i * br.price_factor

        # Electricity price TODO: update contextual parameters to U.S. values
        b = price_dct['Electricity']
        D = shape.Triangle(lower=0.08, midpoint=b, upper=0.14)
        @param(name='Electricity price', element='TEA', kind='isolated',
               units='$/kWh', baseline=b, distribution=D)
        def set_electricity_price(i):
            PowerUtility.price = i

        # Electricity GWP
        b = GWP_dct['Electricity']
        D = shape.Triangle(lower=b*0.9, midpoint=b, upper=b*1.1)
        @param(name='Electricity CF', element='LCA', kind='isolated',
                   units='kg CO2-eq/kWh', baseline=b, distribution=D)
        def set_electricity_CF(i):
            GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i



    ##### Specific units #####
    # biosolids_unit = unit_dct['Biosolids']          #overriding simulation values when used
    # batch_setting_unit_params(biosolids_data, model, biosolids_unit)
    
    # Control box (industrial control panel)
    control_unit = unit_dct['ControlBox']
    exclude = ('certified_electrician_wages', 'service_team_wages', 'facility_manager_wages', 'biomass_controls_wages',) if location_specific else ()
    batch_setting_unit_params(controls_data, model, control_unit, exclude)

    # Housing biogenic refinery
    housing_unit = unit_dct['Housing']
    batch_setting_unit_params(housing_data, model, housing_unit)

    # Carbonizer base
    carbonizer_unit = unit_dct['Carbonizer']
    exclude = ('service_team_wages',) if location_specific else ()
    batch_setting_unit_params(carbonizer_data, model, carbonizer_unit, exclude)

    # Pollution control device
    pollution_control_unit = unit_dct['PCD']
    exclude = ('service_team_wages',) if location_specific else ()
    batch_setting_unit_params(pollution_control_data, model, pollution_control_unit, exclude)

    # Oil heat exchanger
    ohx_unit = unit_dct['OilHX']
    batch_setting_unit_params(ohx_data, model, ohx_unit)

    # Hydronic heat exchanger
    hhx_unit = unit_dct['HHX']
    exclude = ('service_team_wages',) if location_specific else ()
    batch_setting_unit_params(hhx_data, model, hhx_unit, exclude)

    # Dryer from HHX
    hhx_dryer_unit = unit_dct['HHXdryer']
    batch_setting_unit_params(hhx_dryer_data, model, hhx_dryer_unit)

    ##### Universal degradation parameters #####
    # Max methane emission
    HHX_unit = sys.path[6] # the first unit that involves degradation
    b = br.max_CH4_emission
    D = shape.Triangle(lower=0.175, midpoint=b, upper=0.325)
    @param(name='Max CH4 emission', element=HHX_unit, kind='coupled', units='g CH4/g COD',
           baseline=b, distribution=D)
    def set_max_CH4_emission(i):
        br.max_CH4_emission = i
        for unit in sys.units:
            if hasattr(unit, 'max_CH4_emission'):
                setattr(unit, 'max_CH4_emission', i)
    ##### General TEA settings #####
    # # Keeping discount rate constant
    # b = br.discount_rate
    # D = shape.Uniform(lower=0.03, upper=0.06)
    # @param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
    #         baseline=b, distribution=D)
    # def set_discount_rate(i):
    #     br.discount_rate = i

    # Discount factor for the excreta-derived fertilizers
    b = br.price_factor
    D = shape.Uniform(lower=0.1, upper=0.4)
    @param(name='Price factor', element='TEA', kind='isolated', units='-',
           baseline=b, distribution=D)
    def set_price_factor(i):
        br.price_factor = i


    ##### General LCA settings #####
    # TODO: there's probably a more elegant way of adding these CFs
    # CH4
    b = GWP_dct['CH4']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='CH4 CF', element='LCA', kind='isolated', units='kg CO2-eq/kg CH4',
           baseline=b, distribution=D)
    def set_CH4_CF(i):
        GWP_dct['CH4'] = ImpactItem.get_item('CH4_item').CFs['GlobalWarming'] = i

    # N20
    b = GWP_dct['N2O']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='N2O CF', element='LCA', kind='isolated', units='kg CO2-eq/kg N2O',
           baseline=b, distribution=D)
    def set_N2O_CF(i):
        GWP_dct['N2O'] = ImpactItem.get_item('N2O_item').CFs['GlobalWarming'] = i

    
        # Recovered biochar
        b = 0.2 * 0.9 * (44 / 12)  # assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
        D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
        @param(name='biochar CF', element='LCA', kind='isolated',
               units='kg CO2-eq/kg biochar', baseline=b, distribution=D)
        def set_biochar_CF(i):
            GWP_dct['biochar'] = ImpactItem.get_item('biochar_item').CFs['GlobalWarming'] = -i
    

    # Other CFs
    #item_path = os.path.join(br_data_path, 'impact_items.xlsx')
    item_path = os.path.join(data_path, 'impact_items.xlsx')

    for indicator in ('GlobalWarming', 'H_Ecosystems', 'H_Health', 'H_Resources'):
        sheet_name = indicator if indicator!='GlobalWarming' else 'GWP'
        data = load_data(item_path, sheet=sheet_name)
        for p in data.index:
            item = ImpactItem.get_item(p)
            b = item.CFs[indicator]
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
            model.parameter(name=p+f'-{indicator}',
                            setter=DictAttrSetter(item, 'CFs', indicator),
                            element='LCA', kind='isolated',
                            units=f'kg CO2-eq/{item.functional_unit}',
                            baseline=b, distribution=D)


# %%

# =============================================================================
# Functions to create models
# =============================================================================

def create_modelA(location_specific=False, **model_kwargs):
    br.load()
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysA = create_system('A', flowsheet=flowsheet)
    unitA = sysA.flowsheet.unit
   
    # Shared metrics/parameters
    modelA = Model(sysA, **model_kwargs)
    add_metrics(modelA)
    unit_dctA = {
        'Biosolids': unitA.A1,
        'ControlBox': unitA.A2, # industrial control panel
        'Housing': unitA.A3,
        'Carbonizer': unitA.A4,
        'PCD': unitA.A5, # pollution control device
        'OilHX': unitA.A6,
        'HHX': unitA.A7,
        'HHXdryer': unitA.A8,
        }
    add_shared_parameters(modelA, unit_dctA, location_specific)

    return modelA



# Wrapper function so that it'd work for all
def create_model(model_ID='A', location_specific=False, **model_kwargs):
    model_ID = model_ID.lower().rsplit('model')[-1].rsplit('sys')[-1].upper() # works for "modelA"/"sysA"/"A"
    if model_ID == 'A': model = create_modelA(location_specific, **model_kwargs)
    # can add if aother models are used
    #elif model_ID == 'B': model = create_modelB(location_specific, **model_kwargs) 
    else: raise ValueError(f'`model_ID` can only be "A" not "{model_ID}".')
    return model


def run_uncertainty(model, path='', **kwargs):
    kwargs['path'] = os.path.join(results_path, f'sys{model.system.ID[-1]}_model.xlsx') if path=='' else path
    run(model=model, **kwargs)
    return