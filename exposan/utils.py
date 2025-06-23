#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Hannah Lohman <hlohman94@gmail.com>

This module is under the University of Illinois/NCSA Open Source License. Please refer to
https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np, pandas as pd
from math import log
from datetime import datetime
from sklearn.linear_model import LinearRegression as LR
from chaospy import distributions as shape
from thermosteam.functional import rho_to_V
from qsdsan import ImpactItem, sanunits as su
from qsdsan.utils import (
    AttrSetter,
    DictAttrSetter,
    load_data,
    time_printer, 
    )
from . import es_path



__all__ = (
    '_init_modules',
    'add_fugitive_items',
    'add_V_from_rho',
    'batch_setting_LCA_params',
    'batch_setting_unit_params',
    'clear_unit_costs',
    'general_country_specific_inputs',
    'get_decay_k',
    'get_generic_scaled_capital',
    'get_generic_tanker_truck_fee',
    'organize_and_save_results',
    'run_module_country_specific',
    'run_uncertainty',
    )


def _init_modules(module_name, include_data_path=False, include_figures_path=False):
    module_path = os.path.join(es_path, module_name)
    dirnames = ['results']
    if include_data_path: dirnames.insert(0, 'data')
    if include_figures_path: dirnames.append('figures')
    paths = []
    for dirname in dirnames:
        p = os.path.join(module_path, dirname)
        paths.append(p)
        if not os.path.isdir(p): os.mkdir(p)
    return paths


def add_fugitive_items(unit, item_ID):
    unit._run()
    for i in unit.ins:
        i.stream_impact_item = ImpactItem.get_item(item_ID).copy(set_as_source=True)


def add_V_from_rho(component, rho, phase=''):
    '''
    Add constant molar volume model for a component by providing density.

    Parameters
    ----------
    component : obj
        The component object.
    rho : float
        Density in kg/m3.
    phase : str
        For which phase the constant molar volume model should be add,
        must be provided if the component is not locked at a certain phase.

    Examples
    --------    
    >>> from qsdsan import Component
    >>> from exposan.utils import add_V_from_rho
    >>> H2O = Component('H2O', particle_size='Soluble',
    ...                 degradability='Undegradable', organic=False)
    >>> # Add a constant model when H2O is not locked at a state
    >>> add_V_from_rho(H2O, rho=1000, phase='l')
    >>> H2O.V.l(300, 101325) # doctest: +ELLIPSIS
    1.80...
    >>> # Lock H2O at liquid phase then add the model
    >>> H2O.at_state('l')
    >>> add_V_from_rho(H2O, rho=1500)
    >>> H2O.V(300, 101325) # doctest: +ELLIPSIS
    1.20...
    '''
    V_model = rho_to_V(rho, component.MW)
    state = component.locked_state

    if state:
        if phase and state != 'phase':
            raise ValueError(f'The component "{component.ID}" is locked at phase '
                             f'"{state}", cannot add model for phase "{phase}".')
        handle = component.V
    else:
        if not phase:
            raise ValueError('`phase` must be provided for component that is not locked '
                             f'at a certain phase ("{component.ID}").')
        handle = getattr(component.V, phase)
    handle.add_model(V_model)


def batch_setting_LCA_params(path, model, exclude=()):
    '''
    Adding LCA uncertain parameters (characterization factors of different impact items)
    through an Excel spreadsheet.
    
    Parameters
    ----------
    path : str
        Path of the Excel spreadsheet containing the CF distributions.
    model : obj
        The model object where the parameters will be added to.
    exclude : iterable
        IDs of the indicators to be excluded from being added to the model.
        
    See Also
    --------
    Refer to the `pou_disinfection` module (models.py) for an example usage.
    '''
    lca = model.system.LCA
    if isinstance(exclude, str): exclude = (exclude,)
    for ind in lca.indicators:
        ind_ID = ind.ID
        ind_unit = ind.unit
        if ind_ID in exclude: continue
        data = load_data(path, sheet=ind_ID)
        for p in data.index:
            item = ImpactItem.get_item(p)
            b = item.CFs[ind_ID]
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
                            setter=DictAttrSetter(item, 'CFs', ind_ID),
                            element='LCA', kind='isolated',
                            units=f'{ind_unit}/{item.functional_unit}',
                            baseline=b, distribution=D)


def batch_setting_unit_params(df, model, unit, exclude=()):
    '''
    Adding unit uncertain parameters through an Excel spreadsheet.
    
    Parameters
    ----------
    df : pandas DataFrame
        Datasheet containing the uncertain parameters and distributions.
    model : obj
        The model object where the parameters will be added to.
    unit : obj
        The unit object that the parameters are related to.
    exclude : iterable
        IDs of the parameters (i.e., unit attributes) to be excluded from being added to the model.
        
    See Also
    --------
    Refer to the `pou_disinfection` module (models.py) for an example usage.
    '''
    if isinstance(exclude, str): exclude = (exclude,)
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
            raise ValueError(f'Distribution {dist} not recognized for unit {unit} with parameter {para}.')

        su_type = type(unit).__name__
        if su_type.lower() == 'lagoon':
            su_type = f'{unit.design_type.capitalize()} lagoon'
        name = f'{su_type} {para}'
        model.parameter(setter=AttrSetter(unit, para),
                        name=name, element=unit,
                        kind='coupled', units=df.loc[para]['unit'],
                        baseline=b, distribution=D)


# Costs of WWTP units have been considered in the lumped unit
def clear_unit_costs(sys):
    for i in sys.units:
        if isinstance(i, su.LumpedCost): continue
        i.purchase_costs.clear()
        i.installed_costs.clear()


# Get reduction rate constant k for COD and N, use a function so that k can be
# changed during uncertainty analysis
def get_decay_k(tau_deg=2, log_deg=3):
    k = (-1/tau_deg)*np.log(10**-log_deg)
    return k


def get_generic_scaled_capital(tea, percent_CAPEX_to_scale, number_of_units,
                               percent_limit, learning_curve_percent):
    '''
    Scale capital cost based for the Nth system
    (would be lower than the cost for a single system due to scaling effect).

    Parameters
    ----------
    tea : obj
        TEA obj for the system of interest.
    percent_CAPEX_to_scale : float
        The fraction of the cost of specialty parts/cost of total parts.
    number_of_units : int
        Number of units to be constructed.
    percent_limit : float
        Percent of the lowest cost of the normal cost of a single system.
    learning_curve_percent : float
        The percent factor of the learning curve.
    '''
    CAPEX_to_scale = tea.annualized_CAPEX * percent_CAPEX_to_scale
    CAPEX_not_scaled = tea.annualized_CAPEX - CAPEX_to_scale
    scaled_limited = CAPEX_to_scale * percent_limit
    b = log(learning_curve_percent)/log(2)
    scaled_CAPEX_annualized  = (CAPEX_to_scale - scaled_limited)*number_of_units**b + scaled_limited
    new_CAPEX_annualized = scaled_CAPEX_annualized + CAPEX_not_scaled
    return new_CAPEX_annualized


fitting_dct = {
    3: 21.62,
    4.5: 32.43,
    8: 54.05,
    15: 67.57,
}
def get_generic_tanker_truck_fee(capacity,
                                 fitting_dct=fitting_dct,
                                 emptying_fee=0.15,
                                 exchange_rate=1):
    '''
    Exponential fitting to get the tanker truck fee based on capacity.

    cost = a*capacity**b -> ln(price) = ln(a) + bln(capacity)

    Parameters
    ----------
    capacity : float
        The capacity at which the tanker truck fee will be calculated.
    fitting_dct : dict(float, float)
        Capacity-based cost to develop the exponential fitting correlation,
        keys should be the capacities and values should be the corresponding costs.
        Capacities for fitting.
    emptying_fee : float
        Additional fraction of fee that will be added on top of the given prices.
    exchange_rate : float
        Exchange that will be multiplied to the prices.
    '''
    capacities = np.array(tuple(fitting_dct.keys()))
    costs = np.array(tuple(fitting_dct.values()))
    costs *= (1+emptying_fee)*exchange_rate
    ln_p = np.log(costs)
    ln_cap = np.log(np.array(capacities))
    model = LR().fit(ln_cap.reshape(-1,1), ln_p.reshape(-1,1))
    predicted = model.predict(np.array((np.log(capacity))).reshape(1, -1)).item()
    fee = np.exp(predicted)
    return fee


def organize_and_save_results(
        model, percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
        spearman_results=None, path=''):
    '''
    Organize model simulation results and save as an Excel spreadsheet.
    
    Parameters
    ----------
    model : obj
        Model object (with `model.table` containing simulation results).
    percentiles : iterable
        Characteristic percentiles of the parameters values/results to be generated.
    spearman_results : df
        Either None (if no Spearman rank correlation has been run),
        or one df containing the rho values of the Spearman results,
        or two dfs containing the rho and p-values of the Spearman results.
    path : str
        Path where the output Excel to be saved to.
        Will be saved to the "/results" folder if not provided.
        
    Note
    ----
    The `model` object needs to be evaluated (i.e., simulated) first.
    '''
    dct = {}
    index_p = len(model.get_parameters())
    dct['parameters'] = model.table.iloc[:, :index_p].copy()
    dct['data'] = model.table.iloc[:, index_p:].copy()
    if percentiles:
        dct['percentiles_parameters'] = dct['parameters'].quantile(q=percentiles)
        dct['percentiles_results'] = dct['data'].quantile(q=percentiles)
    try:
        iter(spearman_results)
        dct['spearman_rho'] = spearman_results[0]
        dct['spearman_p'] = spearman_results[1]
    except: dct['spearman_rho'] = spearman_results

    path = os.path.join(es_path, f'sys{model.system.ID[-1]}_model.xlsx') if path=='' else path
    with pd.ExcelWriter(path) as writer:
        dct['parameters'].to_excel(writer, sheet_name='Parameters')
        dct['data'].to_excel(writer, sheet_name='Uncertainty results')
        if percentiles:
            dct['percentiles_parameters'].to_excel(writer, sheet_name='Parameter percentiles')
            dct['percentiles_results'].to_excel(writer, sheet_name='Result percentiles')
        if spearman_results is not None:
            dct['spearman_rho'].to_excel(writer, sheet_name='Spearman_rho')
            if dct.get('spearman_p') is not None: dct['spearman_p'].to_excel(writer, sheet_name='Spearman_p')
        model.table.to_excel(writer, sheet_name='Raw data')


@time_printer
def run_uncertainty(model, seed=None, N=1000, rule='L',
                    percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                    path='', print_time=False):
    if seed: np.random.seed(seed)

    samples = model.sample(N, rule)
    model.load_samples(samples)
    model.evaluate()
    # Spearman's rank correlation
    rho, p = model.spearman_r()
    rho.columns = p.columns = pd.Index([i.name_with_units for i in model.metrics])
    organize_and_save_results(model=model, percentiles=percentiles,
                              spearman_results=(rho, p), path=path)


# Example input dict for country-specific analysis
general_country_specific_inputs = {
    'China': {
        'energy_GWP': 0.745,
        'energy_H_Ecosystems': 0.002336342,
        'energy_H_Health': 0.037590269,
        'energy_H_Resources': 0.02714852,
        'energy_price': 0.084,
        'wages': 6.5736875,  # MURT labor wage USD/hour
        'operator_daily_wage': 52.5895,  # USD/day
        'const_wage': 32.66,  # USD/day
        'certified_electrician_wages': 6.251324074,  # USD/hour
        'service_team_wages': 6.5736875,  # USD/hour
        'facility_manager_wages': 6.018625,  # USD/hour
        'biomass_controls_wages': 6.018625,  # USD/hour
        'e_cal': 3191,
        'p_anim': 40,
        'p_veg': 60.63,
        'food_waste_ratio': 0.15,
        'price_ratio': 0.610,
        'household_size': 3,
        'N_fertilizer_price': 0.939,
        'P_fertilizer_price': 1.744,
        'K_fertilizer_price': 1.119,
        'NH3_fertilizer_price': 0.939*(14/17),
        'struvite_fertilizer_price': 1.744*(31/245),
        'NaCl_price': 0.35,  # for NEWgen
        'LPG_price': 0.954  # for NEWgen
        },
    'India': {
        'energy_GWP': 0.852,
        'energy_H_Ecosystems': 0.002616438,
        'energy_H_Health': 0.043194571,
        'energy_H_Resources': 0.030496415,
        'energy_price': 0.081,
        'wages': 1.315625,  # MURT labor wage USD/hour
        'operator_daily_wage': 10.525,  # USD/day
        'const_wage': 10.3285,  # USD/day
        'certified_electrician_wages': 2.1189375,  # USD/hour
        'service_team_wages': 1.315625,  # USD/hour
        'facility_manager_wages': 2.1026875,  # USD/hour
        'biomass_controls_wages': 2.1026875,  # USD/hour
        'e_cal': 2533,
        'p_anim': 15,
        'p_veg': 48.35,
        'food_waste_ratio': 0.03,
        'price_ratio': 0.300,
        'household_size': 5,
        'N_fertilizer_price': 0.158,
        'P_fertilizer_price': 0.567,
        'K_fertilizer_price': 0.445,
        'NH3_fertilizer_price': 0.158 * (14 / 17),
        'struvite_fertilizer_price': 0.567 * (31 / 245),
        'NaCl_price': 0.47,  # for NEWgen
        'LPG_price': 1.488  # for NEWgen
        },
    'South Africa': {
        'energy_GWP': 0.955,
        'energy_H_Ecosystems': 0.002734378,
        'energy_H_Health': 0.042074692,
        'energy_H_Resources': 0.033168799,
        'energy_price': 0.14,
        'wages': 1.695565104,  # MURT labor wage USD/hour (Africa average)
        'operator_daily_wage': 14.79638636,  # USD/day
        'const_wage': 14.06925,  # USD/day
        'certified_electrician_wages': 2.669459239,  # USD/hour
        'service_team_wages': 1.849548295,  # USD/hour
        'facility_manager_wages': 3.122466346,  # USD/hour
        'biomass_controls_wages': 3.122466346,  # USD/hour
        'e_cal': 2899,
        'p_anim': 36.03,
        'p_veg': 48.33,
        'food_waste_ratio': 0.02,
        'price_ratio': 0.460,
        'household_size': 3,
        'N_fertilizer_price': 0.807,
        'P_fertilizer_price': 5.062,
        'K_fertilizer_price': 0.872,
        'NH3_fertilizer_price': 0.807 * (14 / 17),
        'struvite_fertilizer_price': 5.062 * (31 / 245),
        'NaCl_price': 0.225,  # for NEWgen
        'LPG_price': 2.257  # for NEWgen
        },
    'Senegal': {
        'energy_GWP': 0.939,
        'energy_H_Ecosystems': 0.002819172,
        'energy_H_Health': 0.04520618,
        'energy_H_Resources': 0.033261027,
        'energy_price': 0.186,
        'wages': 1.695565104,  # MURT labor wage USD/hour (Africa average)
        'operator_daily_wage': 14.79638636,  # USD/day
        'const_wage': 14.06925,  # USD/day
        'certified_electrician_wages': 2.669459239,  # USD/hour
        'service_team_wages': 1.849548295,  # USD/hour
        'facility_manager_wages': 3.122466346,  # USD/hour
        'biomass_controls_wages': 3.122466346,  # USD/hour
        'e_cal': 2545,
        'p_anim': 13.69,
        'p_veg': 48.67,
        'food_waste_ratio': 0.02,
        'price_ratio': 0.408,
        'household_size': 9,
        'N_fertilizer_price': 1.400,
        'P_fertilizer_price': 14.049,
        'K_fertilizer_price': 1.506,  # Africa average
        'NH3_fertilizer_price': 1.400 * (14 / 17),
        'struvite_fertilizer_price': 14.049 * (31 / 245),
        'NaCl_price': 0.05,  # for NEWgen
        'LPG_price': 1.422  # for NEWgen (Africa average)
        },
    'Uganda': {
        'energy_GWP': 0.159,
        'energy_H_Ecosystems': 0.001594625,
        'energy_H_Health': 0.036965578,
        'energy_H_Resources': 0.012043899,
        'energy_price': 0.184,
        'wages': 1.3920625,  # MURT labor wage (USD/hour)
        'operator_daily_wage': 11.1365,  # USD/day
        'const_wage': 5.0795,  # USD/day
        'certified_electrician_wages': 1.2930625,  # USD/hour
        'service_team_wages': 1.3920625,  # USD/hour
        'facility_manager_wages': 1.6025625,  # USD/hour
        'biomass_controls_wages': 1.6025625,  # USD/hour
        'e_cal': 1981,
        'p_anim': 12.25,
        'p_veg': 34.69,
        'food_waste_ratio': 0.02,
        'price_ratio': 0.348,
        'household_size': 5,
        'N_fertilizer_price': 1.790,
        'P_fertilizer_price': 3.965,
        'K_fertilizer_price': 1.329,
        'NH3_fertilizer_price': 1.790 * (14 / 17),
        'struvite_fertilizer_price': 3.965 * (31 / 245),
        'NaCl_price': 0.284,  # for NEWgen
        'LPG_price': 1.700  # for NEWgen
        },
    'Median': {
        'energy_GWP': 0.686,
        'energy_H_Ecosystems': 0.002456338,
        'energy_H_Health': 0.040824307,
        'energy_H_Resources': 0.027825633,
        'energy_price': 0.129,
        'wages': 3.6228125,  # MURT labor wage (USD/hour)
        'operator_daily_wage': 29.654,  # USD/day
        'const_wage': 28.153,  # USD/day
        'certified_electrician_wages': 4.5995,  # USD/hour
        'service_team_wages': 3.70675,  # USD/hour
        'facility_manager_wages': 5.1823125,  # USD/hour
        'biomass_controls_wages': 5.1823125,  # USD/hour
        'e_cal': 2864,
        'p_anim': 36.04,
        'p_veg': 43.75,
        'food_waste_ratio': 0.06,
        'price_ratio': 0.479,
        'household_size': 4,
        'N_fertilizer_price': 1.465,
        'P_fertilizer_price': 3.965,
        'K_fertilizer_price': 1.268,
        'NH3_fertilizer_price': 1.465 * (14 / 17),
        'struvite_fertilizer_price': 3.965 * (31 / 245),
        'NaCl_price': 0.284,
        'LPG_price': 1.3916
        },
    'Worst_ECON': {
        'energy_GWP': 1.046968,
        'energy_H_Ecosystems': 0.004594958,
        'energy_H_Health': 0.136580714,
        'energy_H_Resources': 0.036204436,
        'energy_price': 0.378,
        'wages': 42.1214375,  # MURT labor wage (USD/hour)
        'operator_daily_wage': 336.9715,  # USD/day
        'const_wage': 332.76,  # USD/day
        'certified_electrician_wages': 56.512875,  # USD/hour
        'service_team_wages': 42.1214375,  # USD/hour
        'facility_manager_wages': 58.8279375,  # USD/hour
        'biomass_controls_wages': 58.8279375,  # USD/hour
        'e_cal': 1786,
        'p_anim': 6.55,
        'p_veg': 24.81,
        'food_waste_ratio': 0.22,
        'price_ratio': 1.370785956,
        'household_size': 2,
        'N_fertilizer_price': 0.158,
        'P_fertilizer_price': 0.567,
        'K_fertilizer_price': 0.315,
        'NH3_fertilizer_price': 0.158 * (14 / 17),
        'struvite_fertilizer_price': 0.567 * (31 / 245),
        'NaCl_price': 0.47,
        'LPG_price': 2.68128  # Maximum because LPG is an input to the system
        },
    'Best_ECON': {
        'energy_GWP': 0.012,
        'energy_H_Ecosystems': 0.000516572,
        'energy_H_Health': 0.011264692,
        'energy_H_Resources': 0.005237881,
        'energy_price': 0,
        'wages': 0.1251875,  # MURT labor wage (USD/hour)
        'operator_daily_wage': 1.0015,  # USD/day
        'const_wage': 3.25,  # USD/day
        'certified_electrician_wages': 0.1875,  # USD/hour
        'service_team_wages': 0.1251875,  # USD/hour
        'facility_manager_wages': 0.692625,  # USD/hour
        'biomass_controls_wages': 0.692625,  # USD/hour
        'e_cal': 3885,
        'p_anim': 104.98,
        'p_veg': 73.29,
        'food_waste_ratio': 0.02,
        'price_ratio': 0.174,
        'household_size': 9,
        'N_fertilizer_price': 3.283,
        'P_fertilizer_price': 15.244,
        'K_fertilizer_price': 2.560,
        'NH3_fertilizer_price': 3.283 * (14 / 17),
        'struvite_fertilizer_price': 15.244 * (31 / 245),
        'NaCl_price': 0.050,
        'LPG_price': 0.13132  # Minimum because LPG is an input into the system
        },
    'Worst_ENV': {
        'energy_GWP': 1.046968,
        'energy_H_Ecosystems': 0.004594958,
        'energy_H_Health': 0.136580714,
        'energy_H_Resources': 0.036204436,
        'energy_price': 0.378,
        'wages': 42.1214375,  # MURT labor wage (USD/hour)
        'operator_daily_wage': 336.9715,  # USD/day
        'const_wage': 332.76,  # USD/day
        'certified_electrician_wages': 56.512875,  # USD/hour
        'service_team_wages': 42.1214375,  # USD/hour
        'facility_manager_wages': 58.8279375,  # USD/hour
        'biomass_controls_wages': 58.8279375,  # USD/hour
        'e_cal': 3885,
        'p_anim': 104.98,
        'p_veg': 73.29,
        'food_waste_ratio': 0.02,
        'price_ratio': 1.370785956,
        'household_size': 2,
        'N_fertilizer_price': 0.158,
        'P_fertilizer_price': 0.567,
        'K_fertilizer_price': 0.315,
        'NH3_fertilizer_price': 0.158 * (14 / 17),
        'struvite_fertilizer_price': 0.567 * (31 / 245),
        'NaCl_price': 0.47,
        'LPG_price': 2.68128  # Maximum because LPG is an input to the system
        },
    'Best_ENV': {
        'energy_GWP': 0.012,
        'energy_H_Ecosystems': 0.000516572,
        'energy_H_Health': 0.011264692,
        'energy_H_Resources': 0.005237881,
        'energy_price': 0,
        'wages': 0.1251875,  # MURT labor wage (USD/hour)
        'operator_daily_wage': 1.0015,  # USD/day
        'const_wage': 3.25,  # USD/day
        'certified_electrician_wages': 0.1875,  # USD/hour
        'service_team_wages': 0.1251875,  # USD/hour
        'facility_manager_wages': 0.692625,  # USD/hour
        'biomass_controls_wages': 0.692625,  # USD/hour
        'e_cal': 1786,
        'p_anim': 6.55,
        'p_veg': 24.81,
        'price_ratio': 0.174,
        'food_waste_ratio': 0.22,
        'household_size': 9,
        'N_fertilizer_price': 3.283,
        'P_fertilizer_price': 15.244,
        'K_fertilizer_price': 2.560,
        'NH3_fertilizer_price': 3.283 * (14 / 17),
        'struvite_fertilizer_price': 15.244 * (31 / 245),
        'NaCl_price': 0.050,
        'LPG_price': 0.13132  # Minimum because LPG is an input into the system
        },
    }

general_city_specific_inputs = {
    # 'Hongkong': {
    #     'price_ratio': 0.7,
    #     'p_anim': 83.79,
    #     'p_veg': 42.99,
    #     'e_cal': 2996,
    #     'food_waste_ratio': 0.15,
    #     'wages': 10,  # average labor wages USD/hour
    #     'energy_price': 0.1857,  # USD/kWh
    #     'energy_GWP': 0.753928903,  # kg CO2-eq/kWh
    #     'N_fertilizer_price': 3.08695652,  # USD/kg N
    #     'P_fertilizer_price': 4.98866213,  # USD/kg P
    #     'K_fertilizer_price': 1.666666667,  # USD/kg K
    #     'household_size': 2.9,  # cap/household
    #     'energy_H_Ecosystems': 0.013938093,  # points/kWh
    #     'energy_H_Health': 0.027525855,  # points/kWh
    #     'energy_H_Resources': 0.029968687,  # points/kWh
    #     'certified_electrician_wages': 16.17,  # USD/hour
    #     'service_team_wages': 10,  # USD/hour
    #     'facility_manager_wages': 12.21  # USD/hour
    # },
    # 'Addis Ababa': {
    #     'price_ratio': 0.4,
    #     'p_anim': 6,
    #     'p_veg': 59,
    #     'e_cal': 2111,
    #     'food_waste_ratio': 0.02,
    #     'wages': 0.9,  # average labor wages USD/hour
    #     'energy_price': 0.01413,  # USD/kWh
    #     'energy_GWP': 0.006710502,  # kg CO2-eq/kWh
    #     'N_fertilizer_price': 0.529,  # USD/kg N
    #     'P_fertilizer_price': 0.695,  # USD/kg P
    #     'K_fertilizer_price': 0,  # USD/kg K
    #     'household_size': 4.68,  # cap/household
    #     'energy_H_Ecosystems': 0.000652665,  # points/kWh
    #     'energy_H_Health': 0.000511132,  # points/kWh
    #     'energy_H_Resources': 0.000819693  # points/kWh
    # },
    
    # 'Paris': {
    # 'price_ratio': 0.8,
    # 'p_anim': 40.2,
    # 'p_veg': 40.2,
    # 'e_cal': 1986,
    # 'food_waste_ratio': 0.17,
    # 'wages': 17.948,  # average labor wages USD/hour
    # 'energy_price': 0.2016,  # USD/kWh
    # 'energy_GWP': 0.054,  # kg CO2-eq/kWh
    # 'N_fertilizer_price': 9.11,  # USD/kg N
    # 'P_fertilizer_price': 33.8,  # USD/kg P
    # 'K_fertilizer_price': 18.42,  # USD/kg K
    # 'household_size': 1.86,  # cap/household
    # 'energy_H_Ecosystems': 0.001046026,  # points/kWh
    # 'energy_H_Health': 0.002296732,  # points/kWh
    # 'energy_H_Resources': 0.003274977,  # points/kWh
    # 'certified_electrician_wages': 25.10998263,  # USD/hour
    # 'service_team_wages': 17.948,  # USD/hour
    # 'facility_manager_wages': 17.948  # USD/hour
    # },
    
    # 'New York City': {
    # 'price_ratio': 1.77,
    # 'p_anim': 51.96721357,
    # 'p_veg': 26.31278643,
    # 'e_cal': 2093.14,
    # 'food_waste_ratio': 0.3,
    # 'wages': 30.84,  # average labor wages USD/hour
    # 'energy_price': 0.2469,  # USD/kWh
    # 'energy_GWP': 0.596453139,  # kg CO2-eq/kWh
    # 'N_fertilizer_price': 4.485913197,  # USD/kg N
    # 'P_fertilizer_price': 29.28801254,  # USD/kg P
    # 'K_fertilizer_price': 9.487666034,  # USD/kg K
    # 'household_size': 2.67,  # cap/household
    # 'energy_H_Ecosystems': 0.010927451,  # points/kWh
    # 'energy_H_Health': 0.017937517,  # points/kWh
    # 'energy_H_Resources': 0.031353155,  # points/kWh
    # 'certified_electrician_wages': 48.83,  # USD/hour
    # 'service_team_wages': 30.84,  # USD/hour
    # 'facility_manager_wages': 87.15  # USD/hour
    # },
    "Albania": {
       "price_ratio": 0.365,
       "p_anim": 61.75,
       "p_veg": 54.0,
       "e_cal": 3360.0,
       "food_waste_ratio": 0.17,
       "wages": 2.35,
       "energy_price": 0.11,
       "energy_GWP": 0.06,
       "household_size": 3.3
   },
   "Armenia": {
       "price_ratio": 0.324,
       "p_anim": 45.34,
       "p_veg": 49.02,
       "e_cal": 2997.0,
       "food_waste_ratio": 0.17,
       "wages": 1.37,
       "energy_price": 0.08,
       "energy_GWP": 0.47,
       "household_size": 3.54
   },
   "Austria": {
       "price_ratio": 0.83,
       "p_anim": 65.83,
       "p_veg": 43.29,
       "e_cal": 3695.0,
       "food_waste_ratio": 0.17,
       "wages": 21.28,
       "energy_price": 0.24,
       "energy_GWP": 0.31,
       "household_size": 2.27
   },
   "Bangladesh": {
       "price_ratio": 0.374,
       "p_anim": 12.56,
       "p_veg": 48.17,
       "e_cal": 2563.0,
       "food_waste_ratio": 0.03,
       "wages": 1.15,
       "energy_price": 0.07,
       "energy_GWP": 1.03,
       "household_size": 4.47
   },
   "Barbados": {
       "price_ratio": 1.111,
       "p_anim": 51.02,
       "p_veg": 37.67,
       "e_cal": 2956.0,
       "food_waste_ratio": 0.06,
       "wages": 9.19,
       "energy_price": 0.21,
       "energy_GWP": 1.01,
       "household_size": 2.85
   },
   "Belarus": {
       "price_ratio": 0.333,
       "p_anim": 52.2,
       "p_veg": 39.75,
       "e_cal": 3270.0,
       "food_waste_ratio": 0.17,
       "wages": 3.62,
       "energy_price": 0.07,
       "energy_GWP": 1.03,
       "household_size": 2.48
   },
   "Belgium": {
       "price_ratio": 0.824,
       "p_anim": 58.58,
       "p_veg": 41.29,
       "e_cal": 3769.0,
       "food_waste_ratio": 0.17,
       "wages": 26.18,
       "energy_price": 0.3,
       "energy_GWP": 0.43,
       "household_size": 2.36
   },
   "Belize": {
       "price_ratio": 0.658,
       "p_anim": 30.82,
       "p_veg": 41.11,
       "e_cal": 2775.0,
       "food_waste_ratio": 0.06,
       "wages": 3.19,
       "energy_price": 0.23,
       "energy_GWP": 0.22,
       "household_size": 4.3
   },
   "Bolivia": {
       "price_ratio": 0.39,
       "p_anim": 35.46,
       "p_veg": 37.78,
       "e_cal": 2412.0,
       "food_waste_ratio": 0.06,
       "wages": 4.65,
       "energy_price": 0.12,
       "energy_GWP": 0.74,
       "household_size": 3.53
   },
   "Botswana": {
       "price_ratio": 0.429,
       "p_anim": 28.02,
       "p_veg": 37.87,
       "e_cal": 2342.0,
       "food_waste_ratio": 0.02,
       "wages": 4.71,
       "energy_price": 0.11,
       "energy_GWP": 1.05,
       "household_size": 3.52
   },
   "Bulgaria": {
       "price_ratio": 0.388,
       "p_anim": 41.96,
       "p_veg": 41.3,
       "e_cal": 2854.0,
       "food_waste_ratio": 0.17,
       "wages": 3.59,
       "energy_price": 0.14,
       "energy_GWP": 0.5,
       "household_size": 2.34
   },
   "Cambodia": {
       "price_ratio": 0.359,
       "p_anim": 19.02,
       "p_veg": 46.93,
       "e_cal": 2492.0,
       "food_waste_ratio": 0.03,
       "wages": 1.06,
       "energy_price": 0.15,
       "energy_GWP": 0.45,
       "household_size": 4.61
   },
   "Cameroon": {
       "price_ratio": 0.396,
       "p_anim": 11.26,
       "p_veg": 60.63,
       "e_cal": 2733.0,
       "food_waste_ratio": 0.02,
       "wages": 1.45,
       "energy_price": 0.09,
       "energy_GWP": 0.45,
       "household_size": 4.99
   },
   "Chile": {
       "price_ratio": 0.552,
       "p_anim": 49.6,
       "p_veg": 42.04,
       "e_cal": 3029.0,
       "food_waste_ratio": 0.06,
       "wages": 7.58,
       "energy_price": 0.19,
       "energy_GWP": 0.6,
       "household_size": 3.58
   },
   "China": {
       "price_ratio": 0.61,
       "p_anim": 40.0,
       "p_veg": 60.63,
       "e_cal": 3191.0,
       "food_waste_ratio": 0.15,
       "wages": 6.57,
       "energy_price": 0.08,
       "energy_GWP": 0.74,
       "household_size": 3.38
   },
   "Croatia": {
       "price_ratio": 0.48,
       "p_anim": 55.7,
       "p_veg": 34.9,
       "e_cal": 3074.0,
       "food_waste_ratio": 0.17,
       "wages": 7.99,
       "energy_price": 0.16,
       "energy_GWP": 0.34,
       "household_size": 2.8
   },
   "Cyprus": {
       "price_ratio": 0.662,
       "p_anim": 45.6,
       "p_veg": 46.07,
       "e_cal": 3019.0,
       "food_waste_ratio": 0.17,
       "wages": 13.12,
       "energy_price": 0.28,
       "energy_GWP": 0.95,
       "household_size": 2.75
   },
   "Czech Republic": {
       "price_ratio": 0.53,
       "p_anim": 52.83,
       "p_veg": 34.15,
       "e_cal": 3277.0,
       "food_waste_ratio": 0.17,
       "wages": 9.11,
       "energy_price": 0.24,
       "energy_GWP": 0.6,
       "household_size": 2.4
   },
   "Denmark": {
       "price_ratio": 0.969,
       "p_anim": 79.73,
       "p_veg": 37.35,
       "e_cal": 3401.0,
       "food_waste_ratio": 0.17,
       "wages": 35.74,
       "energy_price": 0.32,
       "energy_GWP": 0.4,
       "household_size": 2.1
   },
   "Dominican Republic": {
       "price_ratio": 0.431,
       "p_anim": 33.48,
       "p_veg": 32.82,
       "e_cal": 2892.0,
       "food_waste_ratio": 0.06,
       "wages": 1.33,
       "energy_price": 0.09,
       "energy_GWP": 0.91,
       "household_size": 3.48
   },
   "Ecuador": {
       "price_ratio": 0.521,
       "p_anim": 30.73,
       "p_veg": 35.24,
       "e_cal": 2606.0,
       "food_waste_ratio": 0.06,
       "wages": 3.71,
       "energy_price": 0.1,
       "energy_GWP": 0.32,
       "household_size": 3.78
   },
   "Egypt": {
       "price_ratio": 0.246,
       "p_anim": 26.42,
       "p_veg": 71.1,
       "e_cal": 3292.0,
       "food_waste_ratio": 0.09,
       "wages": 1.08,
       "energy_price": 0.04,
       "energy_GWP": 0.96,
       "household_size": 4.13
   },
   "El Salvador": {
       "price_ratio": 0.457,
       "p_anim": 26.91,
       "p_veg": 50.12,
       "e_cal": 2696.0,
       "food_waste_ratio": 0.06,
       "wages": 2.42,
       "energy_price": 0.19,
       "energy_GWP": 0.34,
       "household_size": 4.07
   },
   "Estonia": {
       "price_ratio": 0.593,
       "p_anim": 65.42,
       "p_veg": 39.96,
       "e_cal": 3247.0,
       "food_waste_ratio": 0.17,
       "wages": 10.23,
       "energy_price": 0.19,
       "energy_GWP": 0.91,
       "household_size": 2.3
   },
   "Finland": {
       "price_ratio": 0.917,
       "p_anim": 75.43,
       "p_veg": 42.57,
       "e_cal": 3343.0,
       "food_waste_ratio": 0.17,
       "wages": 23.83,
       "energy_price": 0.18,
       "energy_GWP": 0.31,
       "household_size": 2.07
   },
   "France": {
       "price_ratio": 0.794,
       "p_anim": 64.45,
       "p_veg": 43.8,
       "e_cal": 3502.0,
       "food_waste_ratio": 0.17,
       "wages": 21.63,
       "energy_price": 0.21,
       "energy_GWP": 0.11,
       "household_size": 2.22
   },
   "Georgia": {
       "price_ratio": 0.3,
       "p_anim": 31.24,
       "p_veg": 45.53,
       "e_cal": 2835.0,
       "food_waste_ratio": 0.17,
       "wages": 1.81,
       "energy_price": 0.06,
       "energy_GWP": 0.22,
       "household_size": 3.34
   },
   "Germany": {
       "price_ratio": 0.807,
       "p_anim": 64.03,
       "p_veg": 41.37,
       "e_cal": 3554.0,
       "food_waste_ratio": 0.17,
       "wages": 27.51,
       "energy_price": 0.38,
       "energy_GWP": 0.58,
       "household_size": 2.05
   },
   "Ghana": {
       "price_ratio": 0.39,
       "p_anim": 15.37,
       "p_veg": 46.23,
       "e_cal": 3035.0,
       "food_waste_ratio": 0.02,
       "wages": 0.13,
       "energy_price": 0.06,
       "energy_GWP": 0.55,
       "household_size": 3.49
   },
   "Greece": {
       "price_ratio": 0.602,
       "p_anim": 60.07,
       "p_veg": 48.09,
       "e_cal": 3382.0,
       "food_waste_ratio": 0.17,
       "wages": 6.15,
       "energy_price": 0.22,
       "energy_GWP": 0.74,
       "household_size": 2.44
   },
   "Guatemala": {
       "price_ratio": 0.512,
       "p_anim": 21.75,
       "p_veg": 47.24,
       "e_cal": 2551.0,
       "food_waste_ratio": 0.06,
       "wages": 2.09,
       "energy_price": 0.25,
       "energy_GWP": 0.49,
       "household_size": 4.81
   },
   "Guyana": {
       "price_ratio": 0.484,
       "p_anim": 35.54,
       "p_veg": 50.97,
       "e_cal": 2913.0,
       "food_waste_ratio": 0.06,
       "wages": 2.46,
       "energy_price": 0.26,
       "energy_GWP": 1.03,
       "household_size": 3.8
   },
   "Honduras": {
       "price_ratio": 0.43,
       "p_anim": 21.94,
       "p_veg": 41.01,
       "e_cal": 2673.0,
       "food_waste_ratio": 0.06,
       "wages": 1.53,
       "energy_price": 0.19,
       "energy_GWP": 0.4,
       "household_size": 4.47
   },
   "Hungary": {
       "price_ratio": 0.479,
       "p_anim": 51.23,
       "p_veg": 38.41,
       "e_cal": 3316.0,
       "food_waste_ratio": 0.17,
       "wages": 6.71,
       "energy_price": 0.12,
       "energy_GWP": 0.44,
       "household_size": 2.6
   },
   "India": {
       "price_ratio": 0.3,
       "p_anim": 15.0,
       "p_veg": 48.35,
       "e_cal": 2533.0,
       "food_waste_ratio": 0.03,
       "wages": 1.32,
       "energy_price": 0.08,
       "energy_GWP": 0.85,
       "household_size": 4.57
   },
   "Israel": {
       "price_ratio": 1.016,
       "p_anim": 74.69,
       "p_veg": 52.26,
       "e_cal": 3528.0,
       "food_waste_ratio": 0.09,
       "wages": 21.94,
       "energy_price": 0.17,
       "energy_GWP": 1.02,
       "household_size": 3.14
   },
   "Italy": {
       "price_ratio": 0.727,
       "p_anim": 56.94,
       "p_veg": 49.62,
       "e_cal": 3503.0,
       "food_waste_ratio": 0.17,
       "wages": 16.71,
       "energy_price": 0.25,
       "energy_GWP": 0.66,
       "household_size": 2.4
   },
   "Jordan": {
       "price_ratio": 0.419,
       "p_anim": 23.22,
       "p_veg": 45.99,
       "e_cal": 2732.0,
       "food_waste_ratio": 0.09,
       "wages": 2.82,
       "energy_price": 0.1,
       "energy_GWP": 0.94,
       "household_size": 4.72
   },
   "Kazakhstan": {
       "price_ratio": 0.357,
       "p_anim": 54.7,
       "p_veg": 36.79,
       "e_cal": 3067.0,
       "food_waste_ratio": 0.09,
       "wages": 1.98,
       "energy_price": 0.04,
       "energy_GWP": 0.94,
       "household_size": 3.5
   },
   "Kenya": {
       "price_ratio": 0.402,
       "p_anim": 14.92,
       "p_veg": 46.83,
       "e_cal": 2197.0,
       "food_waste_ratio": 0.02,
       "wages": 1.03,
       "energy_price": 0.2,
       "energy_GWP": 0.2,
       "household_size": 3.64
   },
   "Korea, Republic of": {
       "price_ratio": 0.724,
       "p_anim": 51.73,
       "p_veg": 48.45,
       "e_cal": 3420.0,
       "food_waste_ratio": 0.15,
       "wages": 20.45,
       "energy_price": 0.11,
       "energy_GWP": 0.77,
       "household_size": 2.53
   },
   "Kyrgyzstan": {
       "price_ratio": 0.239,
       "p_anim": 33.91,
       "p_veg": 50.5,
       "e_cal": 2729.0,
       "food_waste_ratio": 0.09,
       "wages": 1.31,
       "energy_price": 0.01,
       "energy_GWP": 0.13,
       "household_size": 4.21
   },
   "Latvia": {
       "price_ratio": 0.54,
       "p_anim": 63.08,
       "p_veg": 39.64,
       "e_cal": 3229.0,
       "food_waste_ratio": 0.17,
       "wages": 8.3,
       "energy_price": 0.19,
       "energy_GWP": 0.57,
       "household_size": 2.58
   },
   "Lebanon": {
       "price_ratio": 0.499,
       "p_anim": 22.47,
       "p_veg": 46.68,
       "e_cal": 2857.0,
       "food_waste_ratio": 0.09,
       "wages": 3.78,
       "energy_price": 0.08,
       "energy_GWP": 1.03,
       "household_size": 3.8
   },
   "Lithuania": {
       "price_ratio": 0.49,
       "p_anim": 78.85,
       "p_veg": 47.8,
       "e_cal": 3411.0,
       "food_waste_ratio": 0.17,
       "wages": 9.54,
       "energy_price": 0.18,
       "energy_GWP": 0.31,
       "household_size": 2.32
   },
   "Luxembourg": {
       "price_ratio": 0.921,
       "p_anim": 67.47,
       "p_veg": 41.18,
       "e_cal": 3465.0,
       "food_waste_ratio": 0.17,
       "wages": 33.97,
       "energy_price": 0.25,
       "energy_GWP": 0.71,
       "household_size": 2.41
   },
   "Malaysia": {
       "price_ratio": 0.385,
       "p_anim": 43.18,
       "p_veg": 34.28,
       "e_cal": 2845.0,
       "food_waste_ratio": 0.03,
       "wages": 3.82,
       "energy_price": 0.06,
       "energy_GWP": 0.88,
       "household_size": 4.56
   },
   "Mauritius": {
       "price_ratio": 0.465,
       "p_anim": 41.21,
       "p_veg": 47.72,
       "e_cal": 3051.0,
       "food_waste_ratio": 0.02,
       "wages": 3.93,
       "energy_price": 0.15,
       "energy_GWP": 0.89,
       "household_size": 3.48
   },
   "Moldova, Republic of": {
       "price_ratio": 0.33,
       "p_anim": 30.27,
       "p_veg": 31.37,
       "e_cal": 2383.0,
       "food_waste_ratio": 0.17,
       "wages": 2.61,
       "energy_price": 0.11,
       "energy_GWP": 0.98,
       "household_size": 2.89
   },
   "Mongolia": {
       "price_ratio": 0.337,
       "p_anim": 56.85,
       "p_veg": 30.47,
       "e_cal": 2579.0,
       "food_waste_ratio": 0.09,
       "wages": 2.17,
       "energy_price": 0.04,
       "energy_GWP": 0.97,
       "household_size": 4.32
   },
   "Montenegro": {
       "price_ratio": 0.371,
       "p_anim": 70.27,
       "p_veg": 44.26,
       "e_cal": 3500.0,
       "food_waste_ratio": 0.17,
       "wages": 5.03,
       "energy_price": 0.1,
       "energy_GWP": 0.45,
       "household_size": 3.21
   },
   "Netherlands": {
       "price_ratio": 0.854,
       "p_anim": 69.81,
       "p_veg": 36.88,
       "e_cal": 3297.0,
       "food_waste_ratio": 0.17,
       "wages": 24.71,
       "energy_price": 0.18,
       "energy_GWP": 0.84,
       "household_size": 2.23
   },
   "New Zealand": {
       "price_ratio": 0.927,
       "p_anim": 51.57,
       "p_veg": 42.44,
       "e_cal": 3191.0,
       "food_waste_ratio": 0.22,
       "wages": 23.02,
       "energy_price": 0.24,
       "energy_GWP": 0.22,
       "household_size": 2.67
   },
   "Nigeria": {
       "price_ratio": 0.416,
       "p_anim": 7.34,
       "p_veg": 51.28,
       "e_cal": 2572.0,
       "food_waste_ratio": 0.02,
       "wages": 0.72,
       "energy_price": 0.06,
       "energy_GWP": 0.86,
       "household_size": 4.9
   },
   "Norway": {
       "price_ratio": 1.077,
       "p_anim": 67.91,
       "p_veg": 45.52,
       "e_cal": 3371.0,
       "food_waste_ratio": 0.17,
       "wages": 29.96,
       "energy_price": 0.1,
       "energy_GWP": 0.08,
       "household_size": 2.22
   },
   "Pakistan": {
       "price_ratio": 0.262,
       "p_anim": 26.49,
       "p_veg": 40.51,
       "e_cal": 2486.0,
       "food_waste_ratio": 0.03,
       "wages": 1.07,
       "energy_price": 0.06,
       "energy_GWP": 0.68,
       "household_size": 6.8
   },
   "Philippines": {
       "price_ratio": 0.375,
       "p_anim": 26.32,
       "p_veg": 36.22,
       "e_cal": 2662.0,
       "food_waste_ratio": 0.03,
       "wages": 1.91,
       "energy_price": 0.18,
       "energy_GWP": 0.81,
       "household_size": 4.23
   },
   "Poland": {
       "price_ratio": 0.446,
       "p_anim": 57.93,
       "p_veg": 47.66,
       "e_cal": 3537.0,
       "food_waste_ratio": 0.17,
       "wages": 7.76,
       "energy_price": 0.19,
       "energy_GWP": 0.93,
       "household_size": 2.81
   },
   "Portugal": {
       "price_ratio": 0.613,
       "p_anim": 74.37,
       "p_veg": 42.58,
       "e_cal": 3480.0,
       "food_waste_ratio": 0.17,
       "wages": 6.95,
       "energy_price": 0.31,
       "energy_GWP": 0.55,
       "household_size": 2.66
   },
   "Romania": {
       "price_ratio": 0.388,
       "p_anim": 53.58,
       "p_veg": 55.66,
       "e_cal": 3581.0,
       "food_waste_ratio": 0.17,
       "wages": 6.44,
       "energy_price": 0.18,
       "energy_GWP": 0.46,
       "household_size": 2.88
   },
   "Russian Federation": {
       "price_ratio": 0.383,
       "p_anim": 55.59,
       "p_veg": 45.81,
       "e_cal": 3345.0,
       "food_waste_ratio": 0.17,
       "wages": 3.45,
       "energy_price": 0.06,
       "energy_GWP": 0.69,
       "household_size": 2.58
   },
   "Rwanda": {
       "price_ratio": 0.353,
       "p_anim": 8.32,
       "p_veg": 50.8,
       "e_cal": 2188.0,
       "food_waste_ratio": 0.02,
       "wages": 0.63,
       "energy_price": 0.26,
       "energy_GWP": 0.52,
       "household_size": 4.26
   },
   "Saudi Arabia": {
       "price_ratio": 0.472,
       "p_anim": 32.74,
       "p_veg": 54.2,
       "e_cal": 3307.0,
       "food_waste_ratio": 0.09,
       "wages": 11.54,
       "energy_price": 0.05,
       "energy_GWP": 1.05,
       "household_size": 5.6
   },
   "Slovakia": {
       "price_ratio": 0.575,
       "p_anim": 37.41,
       "p_veg": 33.16,
       "e_cal": 2871.0,
       "food_waste_ratio": 0.17,
       "wages": 7.67,
       "energy_price": 0.21,
       "energy_GWP": 0.28,
       "household_size": 2.8
   },
   "Slovenia": {
       "price_ratio": 0.611,
       "p_anim": 51.23,
       "p_veg": 45.16,
       "e_cal": 3191.0,
       "food_waste_ratio": 0.17,
       "wages": 11.98,
       "energy_price": 0.16,
       "energy_GWP": 0.36,
       "household_size": 2.47
   },
   "Spain": {
       "price_ratio": 0.681,
       "p_anim": 66.73,
       "p_veg": 40.82,
       "e_cal": 3322.0,
       "food_waste_ratio": 0.17,
       "wages": 16.78,
       "energy_price": 0.22,
       "energy_GWP": 0.45,
       "household_size": 2.69
   },
   "Sri Lanka": {
       "price_ratio": 0.282,
       "p_anim": 17.74,
       "p_veg": 47.48,
       "e_cal": 2737.0,
       "food_waste_ratio": 0.03,
       "wages": 1.39,
       "energy_price": 0.08,
       "energy_GWP": 0.59,
       "household_size": 3.8
   },
   "Sweden": {
       "price_ratio": 0.911,
       "p_anim": 67.55,
       "p_veg": 38.01,
       "e_cal": 3184.0,
       "food_waste_ratio": 0.17,
       "wages": 23.75,
       "energy_price": 0.18,
       "energy_GWP": 0.08,
       "household_size": 2.2
   },
   "Switzerland": {
       "price_ratio": 1.133,
       "p_anim": 59.86,
       "p_veg": 35.61,
       "e_cal": 3354.0,
       "food_waste_ratio": 0.17,
       "wages": 42.12,
       "energy_price": 0.22,
       "energy_GWP": 0.07,
       "household_size": 2.21
   },
   "Tanzania": {
       "price_ratio": 0.404,
       "p_anim": 12.14,
       "p_veg": 46.79,
       "e_cal": 2373.0,
       "food_waste_ratio": 0.02,
       "wages": 1.6,
       "energy_price": 0.1,
       "energy_GWP": 0.71,
       "household_size": 4.85
   },
   "Thailand": {
       "price_ratio": 0.405,
       "p_anim": 27.06,
       "p_veg": 36.3,
       "e_cal": 2804.0,
       "food_waste_ratio": 0.03,
       "wages": 3.95,
       "energy_price": 0.12,
       "energy_GWP": 0.91,
       "household_size": 3.69
   },
   "Turkey": {
       "price_ratio": 0.324,
       "p_anim": 39.5,
       "p_veg": 71.18,
       "e_cal": 3711.0,
       "food_waste_ratio": 0.09,
       "wages": 5.04,
       "energy_price": 0.09,
       "energy_GWP": 0.72,
       "household_size": 4.07
   },
   "Uganda": {
       "price_ratio": 0.348,
       "p_anim": 12.25,
       "p_veg": 34.69,
       "e_cal": 1981.0,
       "food_waste_ratio": 0.02,
       "wages": 1.39,
       "energy_price": 0.18,
       "energy_GWP": 0.16,
       "household_size": 4.53
   },
   "Ukraine": {
       "price_ratio": 0.274,
       "p_anim": 38.26,
       "p_veg": 48.18,
       "e_cal": 3102.0,
       "food_waste_ratio": 0.17,
       "wages": 2.27,
       "energy_price": 0.05,
       "energy_GWP": 0.41,
       "household_size": 2.46
   },
   "United Kingdom": {
       "price_ratio": 0.848,
       "p_anim": 57.96,
       "p_veg": 45.96,
       "e_cal": 3344.0,
       "food_waste_ratio": 0.17,
       "wages": 21.24,
       "energy_price": 0.26,
       "energy_GWP": 0.53,
       "household_size": 2.27
   },
   "United States": {
       "price_ratio": 1.0,
       "p_anim": 73.48,
       "p_veg": 40.27,
       "e_cal": 3782.0,
       "food_waste_ratio": 0.22,
       "wages": 29.32,
       "energy_price": 0.15,
       "energy_GWP": 0.68,
       "household_size": 2.49
   },
   "Uruguay": {
       "price_ratio": 0.719,
       "p_anim": 46.22,
       "p_veg": 40.28,
       "e_cal": 3202.0,
       "food_waste_ratio": 0.06,
       "wages": 5.85,
       "energy_price": 0.2,
       "energy_GWP": 0.13,
       "household_size": 2.78
   }
}

def run_module_country_specific(
        create_country_specific_model_func,
        run_uncertainty_func,
        folder_path,
        system_IDs,
        country_specific_inputs=None,
        seed=None,
        N=1000
        ):
    country_specific_inputs = country_specific_inputs or general_country_specific_inputs
    models = dict.fromkeys(system_IDs)
    for sys_ID in system_IDs:
        sys_dct = {}
        for n, country in enumerate(country_specific_inputs.keys()):
            kwargs = {
                'ID': sys_ID,
                'country': country,
                'country_data': country_specific_inputs[country],
                }
            if n == 0: # create the model
                model = create_country_specific_model_func(**kwargs)
            else: # reuse the model, just update parameters
                kwargs['model'] = model
                model = create_country_specific_model_func(**kwargs)
            # Run analysis and save results
            path = os.path.join(folder_path, f'{sys_ID}_{country}.xlsx')
            run_uncertainty_func(model=model, path=path, seed=seed, N=N)
            sys_dct[country] = model
        models[sys_ID] = sys_dct
    return models

def run_module_city_specific(
        create_city_specific_model_func,
        run_uncertainty_func,
        folder_path,
        system_IDs,
        note = '',
        city = None,
        city_specific_inputs=None,
        seed=None,
        N=1000,
        **model_kwargs
        ):
    city_specific_inputs = city_specific_inputs or general_city_specific_inputs
    if city is not None:
        for c in city:
            if c not in city_specific_inputs:
                print(f"Warning: {c} is not in the database.")
        city_specific_inputs = {k: v for k, v in city_specific_inputs.items() if k in city}
    models = dict.fromkeys(system_IDs)
    for sys_ID in system_IDs:
        sys_dct = {}
        for n, city in enumerate(city_specific_inputs.keys()):
            kwargs = {
                'ID': sys_ID,
                'city': city,
                'city_data': city_specific_inputs[city],
                **model_kwargs
                }
            # if n == 0: # create the model
            #     model = create_city_specific_model_func(**kwargs)
            # else: # reuse the model, just update parameters
            #     kwargs['model'] = model
            model = create_city_specific_model_func(**kwargs)
            # Run analysis and save results
            path = os.path.join(folder_path, f'sys{sys_ID}_model_{city}_')
            date = datetime.now().strftime('%Y-%m-%d')
            path += date 
            if note:
                path += f'_{note}'
            path += f'_{N}.xlsx'
            run_uncertainty_func(model=model, path=path, seed=seed, N=N)
            sys_dct[city] = model
        models[sys_ID] = sys_dct
    return models