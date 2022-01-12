#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is modified for Biogenic Refinery by:
    Lewis Rowles <stetsonsc@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import numpy as np
import pandas as pd
from chaospy import distributions as shape
from thermosteam.functional import V_to_rho, rho_to_V
from biosteam import PowerUtility
from biosteam.evaluation import Model, Metric
from qsdsan import currency, ImpactItem
from qsdsan.utils import (
    load_data, data_path,
    AttrSetter, AttrFuncSetter, DictAttrSetter,
    FuncGetter,
    time_printer
    )
from exposan import reclaimer as R

getattr = getattr
eval = eval

item_path = R.systems.item_path

__all__ = ( 'modelB', 'modelC', 'result_dct', 
           'run_uncertainty', 'save_uncertainty_results', 'add_metrics',
           'batch_setting_unit_params', 'add_shared_parameters') 
           


# %%

# =============================================================================
# Functions for batch-making metrics and -setting parameters
# =============================================================================

systems = R.systems
sys_dct = systems.sys_dct
#unit_dct = systems.unit_dct
price_dct = systems.price_dct
GWP_dct = systems.GWP_dct
GWP = systems.GWP
get_summarizing_fuctions = systems.get_summarizing_fuctions
func = get_summarizing_fuctions()


def add_metrics(system):
    sys_ID = system.ID
    tea = sys_dct['TEA'][sys_ID]
    lca = sys_dct['LCA'][sys_ID]
    ppl = sys_dct['ppl'][sys_ID]
    unit = f'{currency}/cap/yr'
    cat = 'TEA results'
    metrics = [
        Metric('Net cost', lambda: func['get_annual_cost'](tea, ppl), unit, cat),
        Metric('Annual CAPEX', lambda: func['get_annual_CAPEX'](tea, ppl), unit, cat),
        Metric('Energy', lambda: func['get_annual_energy'](tea, ppl), unit, cat),
        Metric('Annual OPEX', lambda: func['get_annual_OPEX'](tea, ppl), unit, cat),
        Metric('Labor', lambda: func['get_annual_labor'](tea, ppl), unit, cat),
        Metric('Sales', lambda: func['get_annual_sales'](tea, ppl), unit, cat),
        ]
    unit = f'{GWP.unit}/cap/yr'
    cat = 'LCA results'
    metrics.extend([
        Metric('Net emission', lambda: func['get_annual_GWP'](lca, ppl), unit, cat),
        Metric('Construction', lambda: func['get_constr_GWP'](lca, ppl), unit, cat),
        Metric('Transportation', lambda: func['get_trans_GWP'](lca, ppl), unit, cat),
        Metric('Fugitive gas', lambda: func['get_CH4_N2O_GWP'](system, lca, ppl), unit, cat),
        Metric('Stream items', lambda: func['get_stream_items_emission_GWP'](system, lca, ppl), unit, cat),
        Metric('Offset', lambda: func['get_offset_GWP'](lca, ppl), unit, cat),
        Metric('Other', lambda: func['get_other_GWP'](lca, ppl), unit, cat),
        ])
    # for i in ('COD', 'N', 'P', 'K'):
    #     cat = f'{i} recovery'
    #     metrics.extend([
    #         Metric(f'Liquid {i}', FuncGetter(func[f'get_liq_{i}_recovery'], (system, i)), '', cat),
    #         Metric(f'Solid {i}', FuncGetter(func[f'get_sol_{i}_recovery'], (system, i)), '', cat),
    #         Metric(f'Gas {i}', FuncGetter(func[f'get_gas_{i}_recovery'], (system, i)), '', cat),
    #         Metric(f'Total {i}', FuncGetter(func[f'get_tot_{i}_recovery'], (system, i)), '', cat)
    #         ])
    return metrics


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
        model.parameter(setter=AttrSetter(unit, para),
                        name=para, element=unit, kind='coupled', units=df.loc[para]['unit'],
                        baseline=b, distribution=D)  



# %%

# =============================================================================
# Shared by all three systems
# =============================================================================

su_data_path = data_path + 'sanunit_data/'
path = su_data_path + '_drying_bed.tsv'
drying_bed_data = load_data(path)

def add_shared_parameters(sys, model, country_specific=False):
    ########## Related to multiple units ##########
    unit = sys.path[0]
    param = model.parameter
    streams = sys_dct['stream_dct'][sys.ID]
            
        # Electricity price
    b = price_dct['Electricity']
    D = shape.Triangle(lower=0.04, midpoint=b, upper=0.1)
    @param(name='Electricity price', element='TEA', kind='isolated',
        units='$/kWh', baseline=b, distribution=D)
        
    def set_electricity_price(i):
            PowerUtility.price = i

    b = GWP_dct['Electricity']
    D = shape.Triangle(lower=0.6212, midpoint=b, upper=0.7592)
    @param(name='Electricity CF', element='LCA', kind='isolated',
                   units='kg CO2-eq/kWh', baseline=b, distribution=D)
    def set_electricity_CF(i):
            GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i
        
            # GWP_dct['Electricity'] = systems.e_item.CFs['GlobalWarming'] = i''
            
            
    #Paramter for flushing water to avoid default assumptions found in the _toilet.tsv
    #Amount of water used for flushing [kg/cap/hr]
    unit = sys.path[1] #this is my second san unit so teh path is 1 (starts at 0)
    b = unit.flushing_water #unit in toilet
    D = shape.Uniform(lower=0.1667, upper=0.25) #Processes 20-30 L/h -> for 120 users = 0.16667-0.25 kg/cap/hr
    @param(name='Flushing Water Processed',
        element = unit, 
        kind='coupled',
        units='kg/cap/hr',
        baseline=b, distribution=D)
    def set_flushing_water(i):
        unit.flushing_water = i          
        
    # # Comment out to see sensitvity to number of users
    # unit = sys.path[1]
    # b = unit.ppl
    # D = shape.Uniform(lower=50, upper=1000) 
    # @param(name='Number of users for the system',
    #     element = unit, 
    #     kind='coupled',
    #     units='cap/system',
    #     baseline=b, distribution=D)
    # def set_number_users_system(i):
    #     unit.ppl = i      
     
    #Sludge Pasteurization Unit
    unit = sys.path[3]
    tea = sys.TEA
    sysID = sys.ID
    b = unit.sludge_labor_maintenance 
    D = shape.Uniform(lower=2.1, upper=3.9) # or whatever the distribution should be
    @param(name='Sludge Pasteurization for all Units',
        element='Sludge Pasteurization', 
        kind='isolated',
        units='hr/year',
        baseline=b, distribution=D)
    def set_SludgePasteurizationHours(i):
        unit.sludge_labor_maintenance = i
        tea.annual_labor = systems.update_labor_cost(sysID)
        
    unit = sys.path[3]
    tea = sys.TEA
    sysID = sys.ID
    b = unit.wages
    D = shape.Triangle(lower=14.55, midpoint = 29.11, upper=43.68)
    @param(name='Wages that affect all labor costs in the system',
        element='Sludge Pasteurization Wages', 
        kind='isolated',
        units='USD/cap/day',
        baseline=b, distribution=D)
    def set_SludgePasteurizationWages(i):
        unit.wages = i
        tea.annual_labor = systems.update_labor_cost(sysID)
    
    #Ion Exchange Unit for Labor Uncertainty Costs
    unit = sys.path[5]
    tea = sys.TEA
    sysID = sys.ID
    b = unit.labor_maintenance_zeolite_regeneration # I'm using a fake parameter as an example, you'll need to update all the `XXX`
    D = shape.Triangle(lower=10, midpoint = 20, upper=30) # or whatever the distribution should be
    @param(name='Ion Exchnage Labor Hours for all systems',
        element='Ion Exchange', 
        kind='isolated',
        units='hr/year',
        baseline=b, distribution=D)
    def set_IonExchangeHours(i):
        unit.labor_maintenance_zeolite_regeneration = i
        tea.annual_labor = systems.update_labor_cost(sysID)
        
    unit = sys.path[5]
    tea = sys.TEA
    sysID = sys.ID
    b = unit.wages
    D = shape.Triangle(lower=14.55, midpoint = 29.11, upper=43.68)
    @param(name='Wages that affect all labor costs in the system',
        element='IonExchange Wages', 
        kind='isolated',
        units='USD/cap/day',
        baseline=b, distribution=D)
    def set_IonExchangeWages(i):
        unit.wages = i
        tea.annual_labor = systems.update_labor_cost(sysID)
    

    
    ########## Related to human input ##########
    # Diet and excretion
    path = data_path + 'sanunit_data/_excretion.tsv'
    data = load_data(path)
    batch_setting_unit_params(data, model, unit=sys.path[0])
    
    # Household size
    b = systems.get_household_size()
    D = shape.Normal(mu=b, sigma=1.8)
    @param(name='Household size', element=unit, kind='coupled', units='cap/household',
            baseline=b, distribution=D)
    def set_household_size(i):
        systems.household_size = max(1, i)
    
    # Toilet density
    b = systems.get_household_per_toilet()
    D = shape.Uniform(lower=3, upper=5)
    @param(name='Toilet density', element=unit, kind='coupled', units='household/toilet',
            baseline=b, distribution=D)
    def set_toilet_density(i):
        systems.household_per_toilet = i

    ##### Universal degradation parameters #####
    # Max methane emission
    unit = sys.path[1] # the first unit that involves degradation
    b = systems.get_max_CH4_emission()
    D = shape.Triangle(lower=0.175, midpoint=b, upper=0.325)
    @param(name='Max CH4 emission', element=unit, kind='coupled', units='g CH4/g COD',
           baseline=b, distribution=D)
    def set_max_CH4_emission(i):
        systems.max_CH4_emission = i
        
    
    # Time to full degradation
    b = systems.tau_deg
    D = shape.Uniform(lower=1, upper=3)
    @param(name='Full degradation time', element=unit, kind='coupled', units='yr',
           baseline=b, distribution=D)
    def set_tau_deg(i):
        systems.tau_deg = i
    
    #!!! Reduction at full degradation
    b = systems.log_deg
    D = shape.Uniform(lower=2, upper=4)
    @param(name='Log degradation', element=unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_log_deg(i):
        systems.log_deg = i

    ######## General TEA settings ########
    
    
    # Money discount rate
    # keep discount rate constant
    b = systems.get_discount_rate()
    D = shape.Uniform(lower=0.03, upper=0.06)
    @param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
            baseline=b, distribution=D)
    def set_discount_rate(i):
        systems.discount_rate = i
 
    
    ######## General LCA settings ########
    b = GWP_dct['CH4']
    D = shape.Uniform(lower=28, upper=34)
    @param(name='CH4 CF', element='LCA', kind='isolated', units='kg CO2-eq/kg CH4',
           baseline=b, distribution=D)
    def set_CH4_CF(i):
        GWP_dct['CH4'] = systems.CH4_item.CFs['GlobalWarming'] = i    
   
    
    b = GWP_dct['N2O']
    D = shape.Uniform(lower=265, upper=298)
    @param(name='N2O CF', element='LCA', kind='isolated', units='kg CO2-eq/kg N2O',
           baseline=b, distribution=D)
    def set_N2O_CF(i):
        GWP_dct['N2O'] = systems.N2O_item.CFs['GlobalWarming'] = i
       
    b = -GWP_dct['N']
    D = shape.Triangle(lower=1.8, midpoint=b, upper=8.9)
    @param(name='N fertilizer CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg N', baseline=b, distribution=D)
    def set_N_fertilizer_CF(i):
        GWP_dct['N'] = systems.N_item.CFs['GlobalWarming'] = -i
        
    b = -GWP_dct['P']
    D = shape.Triangle(lower=4.3, midpoint=b, upper=5.4)
    @param(name='P fertilizer CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg P', baseline=b, distribution=D)
    def set_P_fertilizer_CF(i):
        GWP_dct['P'] = systems.P_item.CFs['GlobalWarming'] = -i
        
    b = -GWP_dct['K']
    D = shape.Triangle(lower=1.1, midpoint=b, upper=2)
    @param(name='K fertilizer CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg K', baseline=b, distribution=D)
    def set_K_fertilizer_CF(i):
        GWP_dct['K'] = systems.K_item.CFs['GlobalWarming'] = -i

    data = load_data(item_path, sheet='GWP')    
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
        model.parameter(name=p,
                        setter=DictAttrSetter(item, 'CFs', 'GlobalWarming'),
                        element='LCA', kind='isolated',
                        units=f'kg CO2-eq/{item.functional_unit}',
                        baseline=b, distribution=D)
    
    return model


# %%

# =============================================================================
# Scenario A (sysA): Trucking Scenario (O&M sludge?)
# =============================================================================

#!!! System is still under progress (might be a possible sludge pasteruziation 
#!!! O&M scenario or link with the zyclone)

# sysA = systems.sysA
# sysA.simulate()
# modelA = Model(sysA, add_metrics(sysA))
# paramA = modelA.parameter

# # Shared parameters
# modelA = add_shared_parameters(sysA, modelA)



# # # Toilet and conveyance
# # modelA = add_toilet_parameters(sysA, modelA)

# # #MURT TOILET
# A2 = systems.A2
# path = su_data_path + '_murt_toilet.tsv'
# data = load_data(path)
# batch_setting_unit_params(data, modelA, A2)

# #primary treatment without struvite
# A3 = systems.A3
# path = su_data_path + '_primary_reclaimer.csv'
# data = load_data(path)
# batch_setting_unit_params(data, modelA, A3)

# ##Ultrafiltrtion
# A4 = systems.A4
# path = su_data_path + '_ultrafiltration_reclaimer.csv'
# data = load_data(path)
# batch_setting_unit_params(data, modelA, A4)

# # Ion exchange
# A5 = systems.A5
# path = su_data_path + '_ion_exchange_reclaimer.csv'
# data = load_data(path)
# batch_setting_unit_params(data, modelA, A5)

# # ECR
# A6 = systems.A6
# path = su_data_path + '_ECR_reclaimer.csv'
# data = load_data(path)
# batch_setting_unit_params(data, modelA, A6)

# #Housing
# A9 = systems.A9
# path = su_data_path + '_housing_reclaimer.csv'
# data = load_data(path)
# batch_setting_unit_params(data, modelA, A9)

# #System Controls
# A10 = systems.A10
# path = su_data_path + '_system_reclaimer.csv'
# data = load_data(path)
# batch_setting_unit_params(data, modelA, A10)



# # Conveyance
# A11 = systems.A11
# b = A11.loss_ratio
# D = shape.Uniform(lower=0.02, upper=0.05)
# @paramA(name='Transportation loss', element=A11, kind='coupled', units='fraction',
#         baseline=b, distribution=D)
# def set_trans_loss(i):
#     A11.loss_ratio = A11.loss_ratio = i

# b = A11.single_truck.distance
# D = shape.Uniform(lower=2, upper=10)
# @paramA(name='Transportation distance', element=A11, kind='coupled', units='km',
#         baseline=b, distribution=D)
# def set_trans_distance(i):
#     A11.single_truck.distance = i

# all_paramsA = modelA.get_parameters()


# %%

# =============================================================================
# Scenario B (sysB): Sludge Pasteurization - Baseline Scenario
# =============================================================================

sysB = systems.sysB
sysB.simulate()
modelB = Model(sysB, add_metrics(sysB))
paramB = modelB.parameter


# Shared parameters
modelB = add_shared_parameters(sysB, modelB)

#primary treatment without struvite
B3 = systems.B3
path = su_data_path + '_primary_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelB, B3)

#SysB Sludge Pasteurization
B4 = systems.B4
path = su_data_path + '_sludge_pasteurization.tsv'
data = load_data(path)
batch_setting_unit_params(data, modelB, B4)

# Ultrafiltration
B5 = systems.B5
path = su_data_path + '_ultrafiltration_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelB, B5)

# Ion exchange
B6 = systems.B6
path = su_data_path + '_ion_exchange_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelB, B6)

# ECR
B7 = systems.B7
path = su_data_path + '_ECR_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelB, B7)


#Housing
B10 = systems.B10
path = su_data_path + '_housing_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelB, B10)

#Mischelaneous
B11 = systems.B11
path = su_data_path + '_system_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelB, B11)

all_paramsB = modelB.get_parameters()


# =============================================================================
# Scenario C (sysC): Primary with struvite + MBR
# =============================================================================
sysC = systems.sysC
sysC.simulate()
modelC = Model(sysC, add_metrics(sysC))
paramC = modelC.parameter


# Model C shared parameters
modelC = add_shared_parameters(sysC, modelC) 


#MURT TOILET
C2 = systems.C2
path = su_data_path + '_murt_toilet.tsv'
data = load_data(path)
batch_setting_unit_params(data, modelC, C2)

#primary treatment without struvite
C3 = systems.C3
path = su_data_path + '_primary_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelC, C3)

# Sludge Pasteurization
C4 = systems.C4
path = su_data_path + '_sludge_pasteurization.tsv'
data = load_data(path)
batch_setting_unit_params(data, modelC, C4)

# Ultrafiltration
C5 = systems.C5
path = su_data_path + '_ultrafiltration_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelC, C5)

# Ion exchange
C6 = systems.C6
path = su_data_path + '_ion_exchange_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelC, C6)

# ECR
C7 = systems.C7
path = su_data_path + '_ECR_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelC, C7)


#Housing
C10 = systems.C10
path = su_data_path + '_housing_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelC, C10)

#Mischelaneous
C11 = systems.C11
path = su_data_path + '_system_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelC, C11)

#Solar costs and impacts
C12 = systems.C12
teaC = systems.sysC.TEA
b = C12.pannel_cleaning 
D = shape.Uniform(lower=10, upper=15) 
@paramC(name='Solar Labor Hours for all systems',
    element='Solar', 
    kind='isolated',
    units='hr/year',
    baseline=b, distribution=D)
def set_SolarLaborHours(i):
    C12.pannel_cleaning = i
    teaC.annual_labor = systems.update_labor_cost(sysC)

path = su_data_path + '_solar_reclaimer.csv'
data = load_data(path)
batch_setting_unit_params(data, modelC, C12)

all_paramsC = modelC.get_parameters()

# %%

# =============================================================================
# Functions to run simulation and generate plots
# =============================================================================

result_dct = {
        'sysB': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        }
models=modelB
@time_printer
def run_uncertainty(model, seed=None, N=10000, rule='L',
                    percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                    print_time=False):
    global result_dct
    if seed:
        np.random.seed(seed)

    samples = model.sample(N, rule)
    model.load_samples(samples)
    model.evaluate()

    # Data organization
    dct = result_dct[model._system.ID]
    index_p = len(model.get_parameters())
    dct['parameters'] = model.table.iloc[:, :index_p].copy()
    dct['data'] = model.table.iloc[:, index_p:].copy()
    if percentiles:
        dct['percentiles'] = dct['data'].quantile(q=percentiles)
        dct['percentiles_parameters'] = dct['parameters'].quantile(q=percentiles)
        

    # Spearman's rank correlation
    spearman_metrics = [model.metrics[i] for i in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
    #spearman_metrics = [model.metrics[i] for i in (0, 3, 12, 16, 20, 24)]
    spearman_results, pvalues = model.spearman_r(model.get_parameters(), spearman_metrics)
    cols = [i.name_with_units for i in spearman_metrics]
    spearman_results.columns = pd.Index(cols)
    pvalues.columns = pd.Index(cols)
    dct['spearman'] = spearman_results
    dct['spearman_p'] = pvalues

# 

def save_uncertainty_results(model, path=''):
    if not path:
        import os
        path = os.path.dirname(os.path.realpath(__file__))
        path += '/results'
        if not os.path.isdir(path):
            os.mkdir(path)
        path += f'/model{model._system.ID[-1]}.xlsx'
        del os
    elif not (path.endswith('xlsx') or path.endswith('xls')):
        extension = path.split('.')[-1]
        raise ValueError(f'Only "xlsx" and "xls" are supported, not {extension}.')
    
    dct = result_dct[model._system.ID]
    if dct['parameters'] is None:
        raise ValueError('No cached result, run model first.')
    with pd.ExcelWriter(path) as writer:
        dct['parameters'].to_excel(writer, sheet_name='Parameters')
        dct['data'].to_excel(writer, sheet_name='Uncertainty results')
        if 'percentiles' in dct.keys():
            dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
        dct['spearman'].to_excel(writer, sheet_name='Spearman')
        dct['spearman_p'].to_excel(writer, sheet_name='Spearman_pvalues')
        model.table.to_excel(writer, sheet_name='Raw data')



