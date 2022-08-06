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
from sklearn.linear_model import LinearRegression as LR
from chaospy import distributions as shape
from thermosteam.functional import rho_to_V
from qsdsan import ImpactItem, sanunits as su
from qsdsan.utils import time_printer, AttrSetter
from . import es_path



__all__ = (
    'add_fugitive_items',
    'add_V_from_rho',
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


def add_fugitive_items(unit, item_ID):
    unit._run()
    for i in unit.ins:
        i.stream_impact_item = ImpactItem.get_item(item_ID).copy(set_as_source=True)


def add_V_from_rho(cmp, rho):
    V_model = rho_to_V(rho, cmp.MW)
    try: cmp.V.add_model(V_model)
    except:
        handle = getattr(cmp.V, cmp.locked_state)
        handle.add_model(V_model)


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
    dct = {}
    index_p = len(model.get_parameters())
    dct['parameters'] = model.table.iloc[:, :index_p].copy()
    dct['data'] = model.table.iloc[:, index_p:].copy()
    if percentiles:
        dct['percentiles_parameters'] = dct['parameters'].quantile(q=percentiles)
        dct['percentiles_results'] = dct['data'].quantile(q=percentiles)
    dct['spearman'] = spearman_results

    path = os.path.join(es_path, f'sys{model.system.ID[-1]}_model.xlsx') if path=='' else path
    with pd.ExcelWriter(path) as writer:
        dct['parameters'].to_excel(writer, sheet_name='Parameters')
        dct['data'].to_excel(writer, sheet_name='Uncertainty results')
        if percentiles:
            dct['percentiles_parameters'].to_excel(writer, sheet_name='Parameter percentiles')
            dct['percentiles_results'].to_excel(writer, sheet_name='Result percentiles')
        if spearman_results is not None: dct['spearman'].to_excel(writer, sheet_name='Spearman')
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
    spearman_results = model.spearman()
    spearman_results.columns = pd.Index([i.name_with_units for i in model.metrics])
    organize_and_save_results(model=model, percentiles=percentiles,
                              spearman_results=spearman_results, path=path)


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