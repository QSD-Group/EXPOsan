#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import models as m, systems
import pandas as pd
from chaospy import distributions as shape
from qsdsan.utils import load_data, data_path
from qsdsan import ImpactItem, PowerUtility, Model

sys_dct = systems.sys_dct
price_dct = systems.price_dct
GWP_dct = systems.GWP_dct
GWP = systems.GWP

import os
su_data_path = os.path.join(data_path, 'sanunit_data/')

# add energy_price as a parameter to the new_model

input_dct = {
#     'China': {'energy_GWP': ('uniform', 0.745, 0.712133646, 0.848557635), 
#               'energy_price': 0.084,
# #             'operator': 50.08,
#               'wages': ('uniform', 6.26, 4.695, 7.825), 
#               'construction': 31.12,
#               'e_cal': 3191,
#               'p_anim': 40,
#               'p_veg': 60.63,
#               'price_ratio': 0.610,
#               'household_size': 3,
#               'N_fertilizer_price': 0.94,
#               'P_fertilizer_price': 1.74,
#               'K_fertilizer_price': 1.12
#                },
      # 'India': {'energy_GWP': ('uniform', 0.852, 0.805105241, 0.958663233), 
      #           'energy_price': 0.081,
      #           # 'operator': 29.12,
      #           'wages': ('uniform', 3.64, 3.094, 4.186), 
      #           'construction': 16.85,
      #           'e_cal': 2533,
      #           'p_anim': 15,
      #           'p_veg': 48.35,
      #           'price_ratio': 0.300,
      #           'household_size': 5,
      #           'N_fertilizer_price': 0.16,
      #           'P_fertilizer_price': 0.57,
      #           'K_fertilizer_price': 0.44
      #           },
      # 'South Africa': {'energy_GWP': ('uniform', 0.955, 0.909088322, 1.086711278),
      #           'energy_price': 0.14,
      #           # 'operator': 23.6,
      #           'wages': ('uniform', 2.95, 2.2125, 3.6875), 
      #           'construction': 11.76,
      #           'e_cal': 2899,
      #           'p_anim': 36.03,
      #           'p_veg': 48.33,
      #           'price_ratio': 0.460,
      #           'household_size': 3,
      #           'N_fertilizer_price': 0.81,
      #           'P_fertilizer_price': 0.78, #source: https://www.namc.co.za/wp-content/uploads/2021/07/Trends-in-selected-Agricultural-input-prices-June-2021.pdf
      #           'K_fertilizer_price': 0.87
      #           },
      # 'Senegal': {'energy_GWP': ('uniform', 0.939, 0.850772468, 1.016179461),  
      #           'energy_price': 0.186,
      #           # 'operator': 29.12,
      #           'wages': ('uniform', 3.64, 3.094, 4.186), 
      #           'construction': 16.85,
      #           'e_cal': 2545,
      #           'p_anim': 13.69,
      #           'p_veg': 48.67,
      #           'price_ratio': 0.408,
      #           'household_size': 9,
      #           'N_fertilizer_price': 1.40,
      #           'P_fertilizer_price': 2.50, 
      #           'K_fertilizer_price': 1.27
       #           },
       # 'Uganda': {'energy_GWP': ('uniform', 0.159, 0.1113, 0.19875), 
       #           'energy_price': 0.184,
       #           # 'operator': 10.64,
       #           'wages': ('uniform', 1.33, 0.9975, 1.6625), 
       #           'construction': 4.80,
       #           'e_cal': 1981,
       #           'p_anim': 12.25,
       #           'p_veg': 34.69,
       #           'price_ratio': 0.348,
       #           'household_size': 5,
       #           'N_fertilizer_price': 1.79,
       #           'P_fertilizer_price': 3.97,
       #           'K_fertilizer_price': 1.33
       #           },
      # 'worst': {'energy_GWP': ('uniform', 1.046968, 0.785, 1.3087), 
      #         'energy_price': 0.378,
      #         'wages': ('uniform', 40.11565476, 30.087, 50.143), 
      #         'construction': 40.7322619,
      #         'e_cal': 1786,
      #         'p_anim': 6.55,
      #         'p_veg': 24.81,
      #         'price_ratio':1.370785956,
      #         'household_size': 2,
      #         'N_fertilizer_price': 0.158446749,
      #         'P_fertilizer_price': 0.567104405,
      #         'K_fertilizer_price': 0.444586022
      #           },
    # 'median': {'energy_GWP': ('uniform', 0.686, 0.5145, 0.8575), 
    #           'energy_price': 0.129,
    #           'wages': ('uniform', 3.70, 2.775, 4.625), 
    #           'construction': 3.64,
    #           'e_cal': 2864,
    #           'p_anim': 36.04,
    #           'p_veg': 43.75,
    #           'price_ratio':0.479,
    #           'household_size': 4,
    #           'N_fertilizer_price': 1.53,
    #           'P_fertilizer_price': 2.50,
    #           'K_fertilizer_price': 1.27
    #             },
    # 'best': {'energy_GWP': ('uniform', 0.012, 0.009, 0.015), 
    #           'energy_price': 0,
    #           'wages': ('uniform', 0.05, 0.0375, 0.0625), 
    #           'construction': 0.14,
    #           'e_cal': 3885,
    #           'p_anim': 104.98,
    #           'p_veg': 73.29,
    #           'price_ratio':0.174,
    #           'household_size': 9,
    #           'N_fertilizer_price': 1.79,
    #           'P_fertilizer_price': 11.20,
    #           'K_fertilizer_price': 2.56
    #             }
    # 'general': {'energy_GWP': ('uniform', 0.69, 0.5175, 0.8625), 
    #           'energy_price': 0.13,
    #           'wages': ('uniform', 3.625, 2.71875, 4.53125), 
    #           'construction': 17,
    #           'e_cal': 2864,
    #           'p_anim': 36.04,
    #           'p_veg': 43.75,
    #           'price_ratio': 1,
    #           'household_size': 4,
    #           'N_fertilizer_price': 1.507,
    #           'P_fertilizer_price': 3.983,
    #           'K_fertilizer_price': 1.333
    #             }
    # 'Worst2_GHG': {'energy_GWP': ('uniform', 1.046968, 0.785, 1.3087),
    #       'energy_price': 0.378,
    #       'wages': ('uniform', 40.11565476, 30.087, 50.143),
    #       'construction': 40.7322619,
    #       'e_cal': 3885,
    #       'p_anim': 104.98,
    #       'p_veg': 73.29,
    #       'price_ratio': 1.370785956,
    #       'household_size': 2,
    #       'N_fertilizer_price': 0.158446749,
    #       'P_fertilizer_price': 0.567104405,
    #       'K_fertilizer_price': 0.444586022
    #       },
    'Best2_GHG': {'energy_GWP': ('uniform', 0.009, 0.012, 0.015),
          'energy_price': 0,
          'wages': ('uniform',  0.0375, 0.05, 0.0625),
          'construction': 0.14,
          'e_cal': 1786,
          'p_anim': 6.55,
          'p_veg': 24.81,
          'price_ratio': 0.174,
          'household_size': 9,
          'N_fertilizer_price': 1.79,
          'P_fertilizer_price': 11.20,
          'K_fertilizer_price': 2.56},
    }

b = price_dct['Electricity']

def run_country(dct):
    results = {}
    

    
    # all_paramsA = modelA.get_parameters()

    
    for country, country_dct in dct.items():
        sysA = systems.sysA
        sysA.simulate()
        streams = sys_dct['stream_dct'][sysA.ID]
        # operator
        #sysA._TEA.annual_labor = country_dct['operator']* 3*365
        
        modelA = Model(sysA, m.add_metrics(sysA))
        # paramA = modelA.parameter
        
        # Shared parameters
        modelA = m.add_shared_parameters(sysA, modelA, country_specific=True)
        param = modelA.parameter
        
        
        #price ratio
        i = country_dct['price_ratio']
        concrete_item = ImpactItem.get_item('Concrete')
        steel_item = ImpactItem.get_item('Steel')
        price_dct = systems.price_dct
        old_price_ratio = systems.price_ratio
        # a much better way would be to have consistent naming in `price_dct` and `streams` dict
        price_ref = {
        'Concrete': concrete_item,
        'Steel': steel_item,
        'Polymer': streams['polymer'],
        'Resin': streams['resin'],
        'FilterBag': streams['filter_bag'],
        'MgOH2': streams['MgOH2'],
        'MgCO3': streams['MgCO3'],
        'H2SO4': streams['H2SO4'],
        'biochar': streams['biochar'],
        }
        for obj_name, obj in price_ref.items():
            old_price = price_dct[obj_name]
            new_price = old_price * i/old_price_ratio
            obj.price = new_price
        systems.price_ratio = i
        for u in sysA.units:
            if hasattr(u, 'price_ratio'):
                u.price_ratio = i
        
        # wages
        kind, low_val, peak_val, max_val = country_dct['wages']
        b=peak_val
        if kind == 'triangle':
            D = shape.Triangle(lower=low_val, midpoint=peak_val, upper=max_val)
        else:
            D = shape.Uniform(lower=low_val,upper=max_val)
        @param(name='Labor wages', element=systems.A1, kind='coupled', units='USD/h',
               baseline=b, distribution=D)
        def set_labor_wages(i):
            labor_cost = 0 
            sysA = systems.sysA
            for u in sysA.units:
                if hasattr(u, '_calc_maintenance_labor_cost'):
                    u.wages = i
                    labor_cost += u._calc_maintenance_labor_cost()
            sysA.TEA.annual_labor = labor_cost
        
        #energy_GWP
        kind, low_val, peak_val, max_val = country_dct['energy_GWP']
        b=peak_val
        if kind == 'triangle':
            D = shape.Triangle(lower=low_val, midpoint=peak_val, upper=max_val)
        else:
            D = shape.Uniform(lower=low_val,upper=max_val)
            
        @param(name='Electricity CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kWh', baseline=b, distribution=D)
        def set_electricity_CF(i):
            GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i
            
        # energy_price
        PowerUtility.price = country_dct['energy_price']
        
        # # household size
        unit = sysA.path[0]
        b = country_dct['household_size']
        D = shape.Normal(mu=b, sigma=1.8)
        @param(name='Household size', element=unit, kind='coupled', units='cap/household',
                baseline=b, distribution=D)
        def set_household_size(i):
            systems.household_size = max(1, i)

        #fertilizer prices
        #N fertilizer
        get_price_factor = systems.get_price_factor
        b = country_dct['N_fertilizer_price']
        D = shape.Uniform(lower=0.16, upper=1.79)
        @param(name='N fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
                baseline=b, distribution=D)
        def set_N_price(i):
            price_dct['N'] = streams['liq_N'] = streams['sol_N'] = i * get_price_factor()
        
        #P fertilizer    
        b = country_dct['P_fertilizer_price']
        D = shape.Uniform(lower=0.57, upper=11.20)
        @param(name='P fertilizer price', element='TEA', kind='isolated', units='USD/kg P',
                baseline=b, distribution=D)
        def set_P_price(i):
            price_dct['P'] = streams['liq_P'] = streams['sol_P'] = i * get_price_factor()
        
        #K fertilizer
        b = country_dct['K_fertilizer_price']
        D = shape.Uniform(lower=0.44, upper=.56)
        @param(name='K fertilizer price', element='TEA', kind='isolated', units='USD/kg K',
                baseline=b, distribution=D)
        def set_K_price(i):
            price_dct['K'] = streams['liq_K'] = streams['sol_K'] = i * get_price_factor()
        
        # Diet and excretion
        A1 = systems.A1
        # p_anim
        A1.p_anim = country_dct['p_anim']
        
        # p_veg
        A1.p_veg = country_dct['p_veg']
        
        # e_cal
        A1.e_cal = country_dct['e_cal']
        
        path = su_data_path + '_excretion.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A1, exclude=('e_cal','p_anim','p_veg')) #adds distributions to tsv file
        
        # Pit latrine and conveyance
        #modelA = m.add_pit_latrine_parameters(sysA, modelA)
        
        # Industrial control panel
        A4 = systems.A4
        path = su_data_path + '_industrial_control_panel.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A4)
        
        # Housing biogenic refinery
        A5 = systems.A5
        
        # construction
        A5.const_wage = country_dct['construction']
        
        path = su_data_path + '_housing_biogenic_refinery.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A5)
        
        # Screw press
        A6 = systems.A6
        path = su_data_path + '_screw_press.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A6)
        
        # Liquid treatment bed
        A7 = systems.A7
        path = su_data_path + '_liquid_treatment_bed.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A7)
        
        # Carbonizer base
        A8 = systems.A8
        path = su_data_path + '_carbonizer_base.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A8)
        
        # Pollution control device
        A9 = systems.A9
        path = su_data_path + '_pollution_control_device.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A9)
        
        # Oil heat exchanger
        A10 = systems.A10
        path = su_data_path + '_oil_heat_exchanger.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A10)
        
        # Hydronic heat exchanger
        A11 = systems.A11
        path = su_data_path + '_hydronic_heat_exchanger.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A11)
        
        # Dryer from HHx
        A12 = systems.A12
        path = su_data_path + '_dryer_from_hhx.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelA, A12)

        

        results[country] = m.run_uncertainty(model=modelA, seed=5, N=1000) 
        del sysA, modelA
        
    return results


results = run_country(dct=input_dct)


#need results folder in the file where this is located
def save_uncertainty_results(results):
    import os
    path = os.path.dirname(os.path.realpath(__file__))
    path = '/Users/lane/Desktop/Coding/GitHub/EXPOsan/exposan/biogenic_refinery'
    path += '/results'
    if not os.path.isdir(path):
         os.mkdir(path)
    del os


    for country, dct in results.items():
        file_name = path+'/'+country+'_A'+'.xlsx'
        if dct['parameters'] is None:
            raise ValueError('No cached result, run model first.')
        with pd.ExcelWriter(file_name) as writer:
            dct['parameters'].to_excel(writer, sheet_name='Parameters')
            dct['data'].to_excel(writer, sheet_name='Uncertainty results')
            if 'percentiles' in dct.keys():
                dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
            dct['spearman'].to_excel(writer, sheet_name='Spearman')
            # model.table.to_excel(writer, sheet_name='Raw data')

save_uncertainty_results(results)

        
        
