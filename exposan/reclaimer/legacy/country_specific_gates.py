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


# %%

import os
from chaospy import distributions as shape
from qsdsan import ImpactItem, PowerUtility
from exposan import reclaimer as re
from exposan.reclaimer import (
    results_path, create_model, run_uncertainty, save_uncertainty_results,
    systems as re_systems,
    )

# These codes are written so that it'll be easier/possible/less confusing
# when running multiple systems at the same time
sys_dct = re_systems.sys_dct
price_dct = re_systems.price_dct
GWP_dct = re_systems.GWP_dct
GWP = re_systems.GWP
old_price_ratio = re_systems.price_ratio


def update_price_ratio(ratio):
    re_systems.price_ratio


# Country-specific inputs
input_dct = {
    'China': {
        'energy_GWP': ('uniform', 0.712133646, 0.745, 0.848557635),
        # add energy ReCiPe CFs (one line of code for each - 3 total)
        'energy_price': 0.084,
        'wages': ('uniform', 4.695, 6.26, 7.825),
        'e_cal': 3191,
        'p_anim': 40,
        'p_veg': 60.63,
        'price_ratio': 0.610,
        },
    'India': {
        'energy_GWP': ('uniform', 0.805105241, 0.852, 0.958663233),
        'energy_price': 0.081,
        'wages': ('uniform', 3.094, 3.64, 4.186),
        'e_cal': 2533,
        'p_anim': 15,
        'p_veg': 48.35,
        'price_ratio': 0.300,
        },
    'South Africa': {
       'energy_GWP': ('uniform', 0.909088322, 0.955, 1.086711278),
       'energy_price': 0.14,
       'wages': ('uniform', 2.2125, 2.95, 3.6875),
       'e_cal': 2899,
       'p_anim': 36.03,
       'p_veg': 48.33,
       'price_ratio': 0.460,
       },
    'Senegal': {
        'energy_GWP': ('uniform', 0.850772468, 0.939, 1.016179461),
        'energy_price': 0.186,
        'wages': ('uniform', 3.094, 3.64, 4.186),
        'e_cal': 2545,
        'p_anim': 13.69,
        'p_veg': 48.67,
        'price_ratio': 0.408,
        },
    'Uganda': {
        'energy_GWP': ('uniform', 0.1113, 0.159, 0.19875),
        'energy_price': 0.184,
        'wages': ('uniform', 0.9975, 1.33, 1.6625),
        'e_cal': 1981,
        'p_anim': 12.25,
        'p_veg': 34.69,
        'price_ratio': 0.348,
        },
    'worst': {
        'energy_GWP': ('uniform', 0.785, 1.046968, 1.3087),
        'energy_price': 0.378,
        'wages': ('uniform', 30.087, 40.11565476, 50.143),
        'e_cal': 1786,
        'p_anim': 6.55,
        'p_veg': 24.81,
        'price_ratio':1.370785956,
        },
    'median': {
        'energy_GWP': ('uniform', 0.5145, 0.686, 0.8575),
        'energy_price': 0.129,
        'wages': ('uniform', 2.775, 3.70, 4.625),
        'e_cal': 2864,
        'p_anim': 36.04,
        'p_veg': 43.75,
        'price_ratio':0.479,
        },
    'best': {
        'energy_GWP': ('uniform', 0.009, 0.012, 0.015),
        'energy_price': 0,
        'wages': ('uniform', 0.0375, 0.05, 0.0625),
        'e_cal': 3885,
        'p_anim': 104.98,
        'p_veg': 73.29,
        'price_ratio':0.174,
        },
    'general': {
        'energy_GWP': ('uniform', 0.5175, 0.69, 0.8625),
        'energy_price': 0.13,
        'wages': ('uniform', 2.71875, 3.625, 4.53125),
        'e_cal': 2864,
        'p_anim': 36.04,
        'p_veg': 43.75,
        'price_ratio': 1,
         },
    'worst2_GHG': {
        'energy_GWP': ('uniform', 0.785, 1.046968, 1.3087),
        'energy_price': 0.378,
        'wages': ('uniform', 30.087, 40.11565476, 50.143),
        'e_cal': 3885,
        'p_anim': 104.98,
        'p_veg': 73.29,
        'price_ratio': 1.370785956,
        },
    'best2_GHG': {
        'energy_GWP': ('uniform', 0.009, 0.012, 0.015),
        'energy_price': 0,
        'wages': ('uniform',  0.0375, 0.05, 0.0625),
        'e_cal': 1786,
        'p_anim': 6.55,
        'p_veg': 24.81,
        'price_ratio': 0.174,
        },
    }

def run_country(input_dct, systems=(), seed=None, N=1000, save_results=True):
    results = {}
    for country, country_dct in input_dct.items():
        for sys in systems:
            streams = sys_dct['stream_dct'][sys.ID]

            # Create model for country-specific analysis
            model = create_model(model_ID=sys.ID, country_specific=True)
            param = model.parameter

            # Diet and excretion
            excretion_unit = sys.path[0]
            excretion_unit.p_anim = country_dct['p_anim']
            excretion_unit.p_veg = country_dct['p_veg']
            excretion_unit.e_cal = country_dct['e_cal']

            # Wages
            kind, low_val, peak_val, max_val = country_dct['wages']
            b = peak_val
            if kind == 'triangle':
                D = shape.Triangle(lower=low_val, midpoint=peak_val, upper=max_val)
            else:
                D = shape.Uniform(lower=low_val,upper=max_val)
            @param(name='Labor wages',
                   element=excretion_unit,  # just want to add to the 0th unit of the system
                   kind='coupled', units='USD/h',
                   baseline=b, distribution=D)
            def set_labor_wages(i):
                labor_cost = 0
                for u in sys.units:
                    if hasattr(u, '_calc_maintenance_labor_cost'):
                        u.wages = i
                        labor_cost += u._calc_maintenance_labor_cost()
                sys.TEA.annual_labor = labor_cost

            # Price ratio
            i = country_dct['price_ratio']
            concrete_item = ImpactItem.get_item('Concrete')
            steel_item = ImpactItem.get_item('Steel')

            price_ref = {
            # 'Concrete': concrete_item,
            # 'Steel': steel_item,
            # add items here that amount of material will influence the cost
            'Polymer': streams['polymer'],
            'Resin': streams['resin'],
            'FilterBag': streams['filter_bag'],
            'MgOH2': streams['MgOH2'],
            'MgCO3': streams['MgCO3'],
            'H2SO4': streams['H2SO4'],
            }
            for obj_name, obj in price_ref.items():
                old_price = price_dct[obj_name]
                new_price = old_price * i/old_price_ratio
                obj.price = new_price
            update_price_ratio(i)
            for u in sys.units:
                if hasattr(u, 'price_ratio'):
                    u.price_ratio = i

            # Energy GWP
            kind, low_val, peak_val, max_val = country_dct['energy_GWP']
            b = peak_val
            if kind == 'triangle':
                D = shape.Triangle(lower=low_val, midpoint=peak_val, upper=max_val)
            else:
                D = shape.Uniform(lower=low_val,upper=max_val)

            @param(name='Electricity CF', element='LCA', kind='isolated',
                    units='kg CO2-eq/kWh', baseline=b, distribution=D)
            def set_electricity_CF(i):
                GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i

            # Create this block of code for each of the 3 ReCiPe indicators (can create a loop to do the 4 indicators)
            # # Energy GWP
            # kind, low_val, peak_val, max_val = country_dct['energy_GWP']
            # b = peak_val
            # if kind == 'triangle':
            #     D = shape.Triangle(lower=low_val, midpoint=peak_val, upper=max_val)
            # else:
            #     D = shape.Uniform(lower=low_val,upper=max_val)
            #
            # @param(name='Electricity CF', element='LCA', kind='isolated',
            #         units='kg CO2-eq/kWh', baseline=b, distribution=D)
            # def set_electricity_CF(i):
            #     GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i

            # Energy price
            PowerUtility.price = country_dct['energy_price']

            # Run analysis and save results
            results[country] = dct = run_uncertainty(model=model, seed=seed, N=N)
            if save_results:
                file_name = os.path.join(results_path, f're_{country}_{sys.ID}.xlsx')
                save_uncertainty_results(model, dct, file_name)

    return results

if __name__ == '__main__':
    results = run_country(input_dct, seed=5, N=1000,
                          # systems=(re.sysA, re.sysB, re.sysC, re.sysD),
                          systems=(re.sysB, re.sysC),
                          save_results=True)
