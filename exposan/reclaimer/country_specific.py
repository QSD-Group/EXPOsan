#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


import os, qsdsan as qs
from chaospy import distributions as shape
from qsdsan import ImpactItem, PowerUtility
from exposan import reclaimer as re
from exposan.reclaimer import (
    create_model,
    GWP_dct,
    H_Ecosystems_dct,
    H_Health_dct,
    H_Resources_dct,
    price_dct,
    run_uncertainty,
    results_path,
    )

__all__ = ('create_country_specific_model',)

# Filter out warnings related to uptime ratio
import warnings
warnings.filterwarnings('ignore', message='uptime_ratio')

from copy import copy
old_price_ratio = copy(re.price_ratio)


# %%

# =============================================================================
# Country-specific inputs
# =============================================================================

input_dct = {
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
        'NaCl': 0.35,  # for NEWgen
        'LPG': 0.954  # for NEWgen
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
        'NaCl': 0.47,  # for NEWgen
        'LPG': 1.488  # for NEWgen
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
        'NaCl': 0.225,  # for NEWgen
        'LPG': 2.257  # for NEWgen
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
        'NaCl': 0.05,  # for NEWgen
        'LPG': 1.422  # for NEWgen (Africa average)
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
        'NaCl': 0.284,  # for NEWgen
        'LPG': 1.700  # for NEWgen
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
        'NaCl': 0.284,
        'LPG': 1.3916
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
        'NaCl': 0.47,
        'LPG': 2.68128  # Maximum because LPG is an input to the system
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
        'NaCl': 0.050,
        'LPG': 0.13132  # Minimum because LPG is an input into the system
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
        'NaCl': 0.47,
        'LPG': 2.68128  # Maximum because LPG is an input to the system
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
        'NaCl': 0.050,
        'LPG': 0.13132  # Minimum because LPG is an input into the system
        },
    }


# %%

# =============================================================================
# Create model for country-specific analysis
# =============================================================================

def create_country_specific_model(ID, country):
    model = create_model(model_ID=ID, country_specific=True)
    param = model.parameter
    sys = model.system
    sys_stream = sys.flowsheet.stream
    country_dct = input_dct[country]

    # Diet and excretion
    excretion_unit = sys.path[0]
    excretion_unit.p_anim = country_dct['p_anim']
    excretion_unit.p_veg = country_dct['p_veg']
    excretion_unit.e_cal = country_dct['e_cal']
    excretion_unit.waste_ratio = country_dct['food_waste_ratio']

    # Wages
    b = country_dct['wages']
    D = shape.Uniform(lower=(b * 0.5), upper=(b * 1.5))
    @param(name='Labor wages',
           element=excretion_unit,  # just want to add to the 0th unit of the system
           kind='coupled', units='USD/h',
           baseline=b, distribution=D)
    def set_labor_wages(i):
        for u in sys.units:
            if hasattr(u, '_calc_maintenance_labor_cost'):
                u.wages = i

    # Price ratio
    i = country_dct['price_ratio']
    re.price_ratio = i
    for u in sys.units:
        if hasattr(u, 'price_ratio'):
            u.price_ratio = i

    # Energy GWP
    b = country_dct['energy_GWP']
    D = shape.Uniform(lower=(b * 0.9), upper=(b * 1.1))
    @param(name='Electricity CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kWh', baseline=b, distribution=D)
    def set_electricity_CF(i):
        GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i

    # Energy H_Ecosystems
    b = country_dct['energy_H_Ecosystems']
    D = shape.Uniform(lower=(b * 0.9), upper=(b * 1.1))
    @param(name='Electricity Ecosystems CF', element='LCA', kind='isolated',
           units='points/kWh', baseline=b, distribution=D)
    def set_electricity_ecosystems_CF(i):
        H_Ecosystems_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Ecosystems'] = i

    # Energy H_Health
    b = country_dct['energy_H_Health']
    D = shape.Uniform(lower=(b * 0.9), upper=(b * 1.1))
    @param(name='Electricity Health CF', element='LCA', kind='isolated',
           units='points/kWh', baseline=b, distribution=D)
    def set_electricity_health_CF(i):
        H_Health_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Health'] = i

    # Energy H_Resources
    b = country_dct['energy_H_Resources']
    D = shape.Uniform(lower=(b * 0.9), upper=(b * 1.1))
    @param(name='Electricity Resources CF', element='LCA', kind='isolated',
           units='points/kWh', baseline=b, distribution=D)
    def set_electricity_resources_CF(i):
        H_Resources_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Resources'] = i

    # Energy price
    PowerUtility.price = country_dct['energy_price']

    # N fertilizer price
    b = country_dct['N_fertilizer_price']
    D = shape.Uniform(lower=(b * 0.8), upper=(b * 1.2))
    @param(name='N fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
           baseline=b, distribution=D)
    def set_N_price(i):
        price_dct['N'] = sys_stream.liq_N.price = sys_stream.sol_N.price = i * re.price_factor

    # P fertilizer price
    b = country_dct['P_fertilizer_price']
    D = shape.Uniform(lower=(b * 0.8), upper=(b * 1.2))
    @param(name='P fertilizer price', element='TEA', kind='isolated', units='USD/kg P',
           baseline=b, distribution=D)
    def set_P_price(i):
        price_dct['P'] = sys_stream.liq_P.price = sys_stream.sol_P.price = i * re.price_factor

    # K fertilizer price
    b = country_dct['K_fertilizer_price']
    D = shape.Uniform(lower=(b * 0.8), upper=(b * 1.2))
    @param(name='K fertilizer price', element='TEA', kind='isolated', units='USD/kg K',
           baseline=b, distribution=D)
    def set_K_price(i):
        price_dct['K'] = sys_stream.liq_K.price = sys_stream.sol_K.price = i * re.price_factor

    # Concentrated NH3 fertilizer price
    b = country_dct['NH3_fertilizer_price']
    D = shape.Uniform(lower=(b * 0.8), upper=(b * 1.2))
    @param(name='NH3 fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
           baseline=b, distribution=D)
    def set_con_NH3_price(i):
        price_dct['conc_NH3'] = sys_stream.conc_NH3.price = i * re.price_factor

    return model

def run_country(input_dct, system_IDs=(), folder_path='', seed=None, N=1000):
    global models
    models = {}
    for country, country_dct in input_dct.items():
        for sys_ID in system_IDs:
            # Create model for country-specific analysis
            model = create_country_specific_model(ID=sys_ID, country=country)
            # Run analysis and save results
            folder_path = folder_path or results_path
            path = os.path.join(folder_path, f're_{country}_{sys_ID}.xlsx')
            run_uncertainty(model=model, path=path, seed=seed, N=N)
            models[country] = model
    return models


if __name__ == '__main__':
    results = run_country(
        input_dct, seed=5, N=10,
        # system_IDs=('sysA', 'sysB', 'sysC', 'sysD'),
        system_IDs=('sysB', 'sysC'),
        )
