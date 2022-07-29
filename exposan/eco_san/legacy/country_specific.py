#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import models as m
import pandas as pd
from systems import price_dct
from chaospy import distributions as shape
import biosteam as bst
from biosteam.evaluation import Model
from exposan import eco_san as es
from qsdsan.utils import (
    load_data, data_path,
    AttrSetter, AttrFuncSetter, DictAttrSetter,
    FuncGetter,
    time_printer
    )
from qsdsan import currency, ImpactItem

systems = es.systems
sys_dct = systems.sys_dct
price_dct = systems.price_dct
GWP_dct = systems.GWP_dct
GWP = systems.GWP


su_data_path = data_path + 'sanunit_data/'

# add energy_price as a parameter to the new_model

input_dct = {
    # 'China': {'energy_GWP': ('triangle', 0.575709657, 0.712133646, 0.848557635), 
    #           'energy_price': 0.08,
    #           'operator': 52.59,
    #           'construction': 32.66,
    #           'e_cal': 3108,
    #           'p_anim':39.28,
    #           'p_veg':58.74,
    #           'price_ratio':0.54},
    # 'India': {'energy_GWP': ('triangle', 0.651547249, 0.805105241, 0.958663233), 
    #           'energy_price': 0.08,
    #           'operator': 29.12,
    #           'construction': 16.85,
    #           'e_cal': 2459,
    #           'p_anim':11.99,
    #           'p_veg': 48.26,
    #           'price_ratio':0.26},
    # 'South Africa': {'energy_GWP': ('triangle', 0.731465366, 0.909088322, 1.086711278), 
    #           'energy_price': 0.11,
    #           'operator': 29.39,
    #           'construction': 11.76,
    #           'e_cal': 3022,
    #           'p_anim':36.39,
    #           'p_veg': 48.94,
    #           'price_ratio':0.46},
    # 'Senegal': {'energy_GWP': ('triangle', 0.685365474, 0.850772468, 1.016179461), 
    #           'energy_price': 0.17,
    #           'operator': 29.12,
    #           'construction': 16.8,
    #           'e_cal': 2456,
    #           'p_anim': 15.25,
    #           'p_veg': 43.24,
    #           'price_ratio':0.4},
    'Uganda': {'energy_GWP': ('triangle', 0.111731123, 0.133335061, 0.154938999), 
              'energy_price': 0.18,
              'operator': 5.37,
              'construction': 6.12,
              'e_cal': 2130,
              'p_anim': 12.39,
              'p_veg': 40.29,
              'price_ratio':0.46},
    }

b = price_dct['Electricity']

def run_country(dct):
    results = {}
    

    
    # all_paramsA = modelA.get_parameters()

    
    for country, country_dct in dct.items():
        sysC = systems.sysC
        sysC.simulate()
        
        # operator
        sysC._TEA.annual_labor = country_dct['operator']* 1*12
        
        modelC = Model(sysC, m.add_metrics(sysC))
        # paramA = modelA.parameter

        
        # Shared parameters
        modelC = m.add_shared_parameters(sysC, modelC, country_specific=True)
        param = modelC.parameter
        
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
        bst.PowerUtility.price = country_dct['energy_price']
        
        
        # # Diet and excretion
        # A1 = systems.A1
        # # p_anim
        # A1.p_anim = country_dct['p_anim']
        
        # # p_veg
        # A1.p_veg = country_dct['p_veg']
        
        # # e_cal
        # A1.e_cal = country_dct['e_cal']
        
        # path = su_data_path + '_excretion.tsv'
        # data = load_data(path)
        # m.batch_setting_unit_params(data, modelA, A1, exclude=('e_cal','p_anim','p_veg'))
        
        
        # # #MURT TOILET
        # A2 = systems.A2
        # path = su_data_path + '_murt_toilet.tsv'
        # data = load_data(path)
        # m.batch_setting_unit_params(data, modelA, A2)
        
        # #primary treatment without struvite
        # A3 = systems.A3
        # path = su_data_path + '_primaryES.tsv'
        # data = load_data(path)
        # m.batch_setting_unit_params(data, modelA, A3)
        
        # ##SysA Bio treatment: A + O + A + O 
        # #Anaerboic 
        # A4 = systems.A4
        # path = su_data_path + '_anaerobic_ES_bio.tsv'
        # data = load_data(path)
        # m.batch_setting_unit_params(data, modelA, A4)


        # # Aerobic 1st Cycle 
        # A5 = systems.A5
        # path = su_data_path + '_aerobic_ES_bio.tsv'
        # data = load_data(path)
        # m.batch_setting_unit_params(data, modelA, A5)
        
        # # Anoxic 
        # A6 = systems.A6
        # path = su_data_path + '_anoxic_ES.tsv'
        # data = load_data(path)
        # m.batch_setting_unit_params(data, modelA, A6)

        # # Aerobic 2nd cycle
        # A7 = systems.A7
        # path = su_data_path + '_aerobic_ES_bio.tsv'
        # data = load_data(path)
        # m.batch_setting_unit_params(data, modelA, A7)

        # ## Chemical treatment
        #     #ECR
        # A8 = systems.A8
        # path = su_data_path + '_ECR.tsv'
        # data = load_data(path)
        # m.batch_setting_unit_params(data, modelA, A8)


        # # Biological Capital Costs and Impacts
        # A13 = systems.A13
        # path = su_data_path + '_bio_ES_1.tsv'
        # data = load_data(path)
        # m.batch_setting_unit_params(data, modelA, A13)
        
        # #Mischelaneous costs and impacts
        # A14 = systems.A14
        # path = su_data_path + '_recycling_controls.tsv'
        # data = load_data(path)
        # m.batch_setting_unit_params(data, modelA, A14)

        # #Solar costs and impacts
        # A15 = systems.A15
        # path = su_data_path + '_solar_ES.tsv'
        # data = load_data(path)
        # m.batch_setting_unit_params(data, modelA, A15)
        
        # # Diet and excretion
        C1 = systems.C1
        #p_anim
        C1.p_anim = country_dct['p_anim']
        
        #p_veg
        C1.p_veg = country_dct['p_veg']
        
        #e_cal
        C1.e_cal = country_dct['e_cal']
        
        path = su_data_path + '_excretion.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelC, C1, exclude=('e_cal','p_anim','p_veg'))
        
        #MURT TOILET
        C2 = systems.C2
        path = su_data_path + '_murt_toilet.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelC, C2)

        #primary treatment without struvite
        C3 = systems.C3
        path = su_data_path + '_primary_MBR.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelC, C3)

        ##SysA Bio treatment: A + O + MBR
        #Anaerboic 
        C4 = systems.C4
        path = su_data_path + '_anaerobic_ES_bio.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelC, C4)

        # Facultative Aerobic 
        C5 = systems.C5
        path = su_data_path + '_aerobic_ES_bio.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelC, C5)

        # MBR
        C6 = systems.C6
        path = su_data_path + '_MBR.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelC, C6)

        # MBR with ECR
        C8 = systems.C8
        path = su_data_path + '_MBR_ECR.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelC, C8)


        #Mischelaneous costs and impacts
        C14 = systems.C14
        path = su_data_path + '_recycling_controls.tsv'
        data = load_data(path)
        m.batch_setting_unit_params(data, modelC, C14)

        







        results[country] = m.run_uncertainty(model=modelC)
        del sysC, modelC
        
    return results


results = run_country(dct=input_dct)


#need results folder in the file where this is located
def save_uncertainty_results(results):
    import os
    #path = os.path.dirname(os.path.realpath(__file__))
    path = ('/Users/torimorgan/opt/anaconda3/lib/python3.8/site-packages/exposan')
    path += '/results'
    if not os.path.isdir(path):
         os.mkdir(path)
    del os


    for country, dct in results.items():
        file_name = path+'/'+country+'.xlsx'
        if dct['parameters'] is None:
            raise ValueError('No cached result, run model first.')
        with pd.ExcelWriter(file_name) as writer:
            dct['parameters'].to_excel(writer, sheet_name='Parameters')
            dct['data'].to_excel(writer, sheet_name='Uncertainty results')
            if 'percentiles' in dct.keys():
                dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
            dct['spearman'].to_excel(writer, sheet_name='Spearman')
            # model.table.to_excel(writer, sheet_name='Raw data')



 # b = systems.get_operator_daily_wage()
 #    D = shape.Triangle(lower=14.55, midpoint=b, upper=43.68)
 #    @param(name='Operator daily wages', element='TEA', kind='cost', units='USD/d',
 #          baseline=b, distribution=D)
 #    def set_operator_daily_wage(i):
 #        sys._TEA.annual_labor = i* 3*365 
        
        
 #        b = price_dct['Electricity']
 #    D = shape.Triangle(lower=0.04, midpoint=b, upper=0.1)
 #    @param(name='Electricity price', element='TEA', kind='isolated',
 #           units='$/kWh', baseline=b, distribution=D)
 #    def set_electricity_price(i):
 #        PowerUtility.price = i
        
        
 #        b = GWP_dct['Electricity']
 #    D = shape.Uniform(lower=0.106, upper=0.121)
 #    @param(name='Electricity CF', element='LCA', kind='isolated',
 #               units='kg CO2-eq/kWh', baseline=b, distribution=D)
 #    def set_electricity_CF(i):
 #        GWP_dct['Electricity'] = ImpactItem.get_item('E_item').CFs['GlobalWarming'] = i
        
        
