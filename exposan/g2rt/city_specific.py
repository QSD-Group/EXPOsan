#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

#%%

from chaospy import distributions as shape
from qsdsan import ImpactItem, PowerUtility
from exposan.utils import general_city_specific_inputs, run_module_city_specific
from exposan import g2rt
from exposan.g2rt import (
    results_path,
    update_resource_recovery_settings,
    )
from exposan.g2rt.models import (
    create_model,
    run_uncertainty
    )
__all__ = ('create_city_specific_model',
           'run_city')

# Filter out warnings related to uptime ratio
import warnings
warnings.filterwarnings('ignore', message='uptime_ratio')

from copy import copy
old_price_ratio = copy(g2rt.price_ratio)

#%%

# =============================================================================
# Create model for city-specific analysis
# =============================================================================

def get_default_uniform(b, ratio, lb=None, ub=None): # lb/ub for upper/lower bounds
    lower = max(b*(1-ratio), lb) if lb else b*(1-ratio)
    upper = min(b*(1+ratio), ub) if ub else b*(1+ratio)
    return shape.Uniform(lower=lower, upper=upper)

format_key = lambda key: (' '.join(key.split('_'))).capitalize()

def create_city_specific_model(ID, city, model=None, city_data=None,
                                  include_non_contextual_params=True):
    city_data = city_data or general_city_specific_inputs[city]
    g2rt.set_default_ppl(city_data['household_size'])
    model = create_model(model_ID=ID, city_specific=True)
    param = model.parameter
    sys = model.system
    sys_stream = sys.flowsheet.stream
    # print(f'current default number is {g2rt.get_default_ppl()},mixed waste flow rate is {sys.flowsheet.unit.A2.outs[0].F_mass} kg/hr')
    price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct = update_resource_recovery_settings()

    if include_non_contextual_params:
        param_dct = {p.name: p for p in model.parameters}
        def get_param_name_b_D(key, ratio, lb=None, ub=None):
            name = format_key(key)
            p = param_dct.get(name)
            # Throw out this parameter (so that a new one can be added with updated values)
            # if it's already there
            if p: model.parameters = [i for i in model.parameters if i is not p]
            b = city_data[key]
            D = get_default_uniform(b, ratio, lb=lb, ub=ub)
            return name, p, b, D
    else:
        model.parameters = [] # throw out all parameters
        def get_param_name_b_D(key, ratio, lb=None, ub=None):
            name = format_key(key)
            p = '' # not used, doesn't matter what the value is
            b = city_data[key]
            D = get_default_uniform(b, ratio, lb=lb, ub=ub)
            return name, p, b, D

    # Diet and excretion
    if sys.ID == 'sysA':
        excretion_unit = sys.flowsheet.unit.A1
    elif sys.ID == 'sysB':
        excretion_unit = sys.flowsheet.unit.B1
    else:
        raise ValueError(f"Unexpected system ID: {sys.ID}")
    dietary_D_ratio = 0.1
    key = 'p_anim'
    name, p, b, D = get_param_name_b_D(key, dietary_D_ratio)
    @param(name=name, element=excretion_unit, kind='coupled', units='g/cap/d',
           baseline=b, distribution=D)
    def set_animal_protein(i):
        excretion_unit.p_anim = i
        
    key = 'p_veg'
    name, p, b, D = get_param_name_b_D(key, dietary_D_ratio)
    @param(name=name, element=excretion_unit, kind='coupled', units='g/cap/d',
           baseline=b, distribution=D)
    def set_vegetal_protein(i):
        excretion_unit.p_veg = i
        
    key = 'e_cal'
    name, p, b, D = get_param_name_b_D(key, dietary_D_ratio)
    @param(name=name, element=excretion_unit, kind='coupled', units='kcal/cap/d',
           baseline=b, distribution=D)
    def set_caloric_intake(i):
        excretion_unit.e_cal = i
        
    key = 'food_waste_ratio'
    name, p, b, D = get_param_name_b_D(key, dietary_D_ratio, lb=0, ub=1)
    @param(name=name, element=excretion_unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_food_waste_ratio(i):
        excretion_unit.waste_ratio = i

    # Price ratio
    price_ratio_D_ratio = 0.2
    key = 'price_ratio'
    name, p, b, D = get_param_name_b_D(key, price_ratio_D_ratio, lb=0)
    stream_ref = {
        'H2O': 'H2O',
    }

    @param(name=name, element=excretion_unit, kind='cost', units='-',
           baseline=b, distribution=D)
    def set_price_ratio(i):
        g2rt.price_ratio = i
        for obj_name in stream_ref.keys():
            old_price = price_dct[obj_name]
            new_price = old_price * i/old_price_ratio
            getattr(sys_stream, stream_ref[obj_name]).price = new_price
        for u in sys.units:
            if hasattr(u, 'price_ratio'):
                u.price_ratio = i

    # Wages
    key = 'wages'
    wage_D_ratio = 0.5
    name, p, b, D = get_param_name_b_D(key, wage_D_ratio)
    @param(name=name,
           element=excretion_unit,  # just want to add to the 0th unit of the system
           kind='coupled', units='USD/h',
           baseline=b, distribution=D)
    def set_labor_wages(i):
        for u in sys.units:
            if hasattr(u, '_calc_maintenance_labor_cost'):
                u.wages = i
        # print(f'current wage  is {i} $/hr, annual OPEX is {sys.TEA.AOC} $/yr')
    
        
    # # User number
    # household_size_D_ratio = 0.1
    # key = 'household_size'
    # name, p, b, D = get_param_name_b_D(key, household_size_D_ratio)
    # @param(name=name, 
    #        element=excretion_unit, 
    #        kind='coupled', 
    #        units='-',
    #        baseline=b, distribution=D)
    # def set_household_size(i):
    #     g2rt.set_default_ppl(i)
    #     # global model
    #     # model = create_model(model_ID=ID, city_specific=True)
    #     print(f'current default number is {g2rt.get_default_ppl()},mixed waste flow rate is {sys.flowsheet.unit.A2.outs[0].F_mass} kg/hr')

    try:
        # Energy price
        energy_price_D_ratio = 0.1
        key = 'energy_price'
        name, p, b, D = get_param_name_b_D(key, energy_price_D_ratio)
        @param(name=name, element='TEA', kind='cost', units='USD/kWh',
               baseline=b, distribution=D)
        def set_energy_price(i):
            PowerUtility.price = i
    except AssertionError:
        warnings.warn("Energy price is zero!!!")


    # Energy GWP
    key = 'energy_GWP'
    GWP_D_ratio = 0.3
    name, p, b, D = get_param_name_b_D(key, GWP_D_ratio)
    @param(name=name, element='LCA', kind='isolated',
           units='kg CO2-eq/kWh', baseline=b, distribution=D)
    def set_electricity_CF(i):
        GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i

    # Energy H_Ecosystems
    key = 'energy_H_Ecosystems'
    name, p, b, D = get_param_name_b_D(key, GWP_D_ratio)
    @param(name=name, element='LCA', kind='isolated',
           units='points/kWh', baseline=b, distribution=D)
    def set_electricity_ecosystems_CF(i):
        H_Ecosystems_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Ecosystems'] = i

    # Energy H_Health
    key = 'energy_H_Health'
    name, p, b, D = get_param_name_b_D(key, GWP_D_ratio)
    @param(name=name, element='LCA', kind='isolated',
           units='points/kWh', baseline=b, distribution=D)
    def set_electricity_health_CF(i):
        H_Health_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Health'] = i

    # Energy H_Resources
    key = 'energy_H_Resources'
    name, p, b, D = get_param_name_b_D(key, GWP_D_ratio)
    @param(name=name, element='LCA', kind='isolated',
           units='points/kWh', baseline=b, distribution=D)
    def set_electricity_resources_CF(i):
        H_Resources_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Resources'] = i

    # ##### Commented out because I am not using country-specific NaCl and am just converting cost
    # ##### with price level ratio - HACL 9/20/22

    # # NaCl price
    # key = 'NaCl_price'
    # price_D_ratio = 0.2
    # name, p, b, D = get_param_name_b_D(key, price_D_ratio)
    # @param(name=name, element='TEA', kind='isolated', units='USD/kg',
    #        baseline=b, distribution=D)
    # def set_NaCl_price(i):
    #     price_dct['NaCl'] = sys_stream.NaCl.price = sys_stream.NaCl1.price = i

    # # LPG
    # key = 'LPG_price'
    # price_D_ratio = 0.2
    # name, p, b, D = get_param_name_b_D(key, price_D_ratio)
    # @param(name=name, element='TEA', kind='isolated', units='USD/kg',
    #        baseline=b, distribution=D)
    # def set_LPG_price(i):
    #     price_dct['LPG'] = sys_stream.LPG.price = i
    price_D_ratio = 0.3
    if g2rt.INCLUDE_RESOURCE_RECOVERY:
        # N fertilizer price
        key = 'N_fertilizer_price'
        name, p, b, D = get_param_name_b_D(key, price_D_ratio)
        @param(name=name, element='TEA', kind='isolated', units='USD/kg N',
               baseline=b, distribution=D)
        def set_N_price(i):
            price_dct['N'] = sys_stream.liq_N.price = sys_stream.sol_N.price = i * g2rt.price_factor
    
        # P fertilizer price
        key = 'P_fertilizer_price'
        name, p, b, D = get_param_name_b_D(key, price_D_ratio)
        @param(name=name, element='TEA', kind='isolated', units='USD/kg P',
               baseline=b, distribution=D)
        def set_P_price(i):
            price_dct['P'] = sys_stream.liq_P.price = sys_stream.sol_P.price = i * g2rt.price_factor
    
        # K fertilizer price
        key = 'K_fertilizer_price'
        name, p, b, D = get_param_name_b_D(key, price_D_ratio)
        @param(name=name, element='TEA', kind='isolated', units='USD/kg K',
               baseline=b, distribution=D)
        def set_K_price(i):
            price_dct['K'] = sys_stream.liq_K.price = sys_stream.sol_K.price = i * g2rt.price_factor

    return model

#%%
def run_city(system_IDs, seed=None, N=1000, city=None, note=''):
    return run_module_city_specific(
        create_city_specific_model_func=create_city_specific_model,
        run_uncertainty_func=run_uncertainty,
        folder_path=results_path,
        system_IDs=system_IDs, #e.g., ['A','B']
        seed=seed,
        N=N,
        city=city, 
        note=note
    )

if __name__ == '__main__':
    results = run_city(
        system_IDs=('sysA', 'sysB'),
        seed=5, N=1000,
        )

