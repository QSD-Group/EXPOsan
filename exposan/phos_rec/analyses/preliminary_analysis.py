#!/usr/bin/env python3 for the CI and Cost heatmap of food_sludge_ratio and HRT
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Xuan Wang <easonbiubiu99@gmail.com>
    
    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

#%% initialization

import os, pandas as pd, qsdsan as qs
from exposan.phos_rec import create_system, create_model
from datetime import date

folder = os.path.dirname(os.path.dirname(__file__))

#%% heatmaps

FePO4_MSP_5th = []
FePO4_MSP_50th = []
FePO4_MSP_95th = []
FePO4_MSP_mean = []

FePO4_GWP_5th = []
FePO4_GWP_50th = []
FePO4_GWP_95th = []
FePO4_GWP_mean = []

sludge_cost_5th = []
sludge_cost_50th = []
sludge_cost_95th = []
sludge_cost_mean = []

sludge_GWP_5th = []
sludge_GWP_50th = []
sludge_GWP_95th = []
sludge_GWP_mean = []

for perspective in ['FePO4','sludge']:
    ratio_list = []
    HRT_list = []
    for ratio in [0, 1/3, 2/3, 1, 4/3]:
        if ratio == 4/3:
            HRTs = [132,]
        else:
            HRTs = range(0, 144, 12)

        for HRT in HRTs:
            try:
                sys = create_system(temp_ratio=10, food_sludge_ratio=ratio, HRT=HRT)
            except RuntimeError:
                continue
            model = create_model(sys, perspective=perspective)
            
            # TODO: increase N for formal analyses
            kwargs = {'N':10,'rule':'L','seed':3221}
            samples = model.sample(**kwargs)
            model.load_samples(samples)
            model.evaluate()
            
            idx = len(model.parameters)
            parameters = model.table.iloc[:, :idx]
            results = model.table.iloc[:, idx:]
            
            results_no_nan = results.dropna()
            
            ratio_list.append(ratio)
            HRT_list.append(HRT)
            
            if perspective == 'FePO4':
                FePO4_MSP_5th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.05))
                FePO4_MSP_50th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.5))
                FePO4_MSP_95th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.95))
                FePO4_MSP_mean.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].mean())
                
                FePO4_GWP_5th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.05))
                FePO4_GWP_50th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.5))
                FePO4_GWP_95th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.95))
                FePO4_GWP_mean.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].mean())
            
            if perspective == 'sludge':
                sludge_cost_5th.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].quantile(0.05))
                sludge_cost_50th.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].quantile(0.5))
                sludge_cost_95th.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].quantile(0.95))
                sludge_cost_mean.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].mean())
                
                sludge_GWP_5th.append(results_no_nan[('LCA','Sludge management GWP [kg_CO2_eq/tonne]')].quantile(0.05))
                sludge_GWP_50th.append(results_no_nan[('LCA','Sludge management GWP [kg_CO2_eq/tonne]')].quantile(0.5))
                sludge_GWP_95th.append(results_no_nan[('LCA','Sludge management GWP [kg_CO2_eq/tonne]')].quantile(0.95))
                sludge_GWP_mean.append(results_no_nan[('LCA','Sludge management GWP [kg_CO2_eq/tonne]')].mean())
    
    if perspective == 'FePO4':
        FePO4_result = pd.DataFrame()
        FePO4_result['ratio'] = ratio_list
        FePO4_result['HRT'] = HRT_list
        FePO4_result['FePO4_MSP_5th'] = FePO4_MSP_5th
        FePO4_result['FePO4_MSP_50th'] = FePO4_MSP_50th
        FePO4_result['FePO4_MSP_95th'] = FePO4_MSP_95th
        FePO4_result['FePO4_MSP_mean'] = FePO4_MSP_mean
        FePO4_result['FePO4_GWP_5th'] = FePO4_GWP_5th
        FePO4_result['FePO4_GWP_50th'] = FePO4_GWP_50th
        FePO4_result['FePO4_GWP_95th'] = FePO4_GWP_95th
        FePO4_result['FePO4_GWP_mean'] = FePO4_GWP_mean
    
    if perspective == 'sludge':
        sludge_result = pd.DataFrame()
        sludge_result['ratio'] = ratio_list
        sludge_result['HRT'] = HRT_list
        sludge_result['sludge_cost_5th'] = sludge_cost_5th
        sludge_result['sludge_cost_50th'] = sludge_cost_50th
        sludge_result['sludge_cost_95th'] = sludge_cost_95th
        sludge_result['sludge_cost_mean'] = sludge_cost_mean
        sludge_result['sludge_GWP_5th'] = sludge_GWP_5th
        sludge_result['sludge_GWP_50th'] = sludge_GWP_50th
        sludge_result['sludge_GWP_95th'] = sludge_GWP_95th
        sludge_result['sludge_GWP_mean'] = sludge_GWP_mean
        
FePO4_result.to_excel(os.path.join(folder, f'results/FePO4_result_{date.today()}.xlsx'))
sludge_result.to_excel(os.path.join(folder, f'results/sludge_result_{date.today()}.xlsx'))

#%% for the simulation of N=1000 and HRT = 132 and food_sludge = 1 to the sensitivity analysis

sys = create_system()
model = create_model(sys)
model.parameters
kwargs = {'N':1000,'rule':'L','seed':3221}
samples = model.sample(**kwargs)
model.load_samples(samples)
model.evaluate()
spearman_kwargs={'nan_policy': 'omit'}
r_df, p_df = qs.stats.get_correlations(model, kind='Spearman', **spearman_kwargs)

r_df.to_excel(os.path.join(folder, f'results/r_df_{date.today()}.xlsx'))
p_df.to_excel(os.path.join(folder, f'results/p_df_{date.today()}.xlsx'))

#%% for the simulation of N=1000 and HRT = 132 and food_sludge_ratio = 1 to the avoided waste management cost and CI analysis (sludge_credit is a parameter)

credit_list = []

FePO4_MSP_5th = []
FePO4_MSP_50th = []
FePO4_MSP_95th = []
FePO4_MSP_mean = []

FePO4_GWP_5th = []
FePO4_GWP_50th = []
FePO4_GWP_95th = []
FePO4_GWP_mean = []

for credit in range(0, 1100, 100):
    sys = create_system(temp_ratio=10, food_sludge_ratio=1, HRT=132, sludge_credit=credit)
    model = create_model(sys, perspective='FePO4')
    
    kwargs = {'N':1000,'rule':'L','seed':3221}
    samples = model.sample(**kwargs)
    model.load_samples(samples)
    model.evaluate()
    
    idx = len(model.parameters)
    parameters = model.table.iloc[:, :idx]
    results = model.table.iloc[:, idx:]
    
    results_no_nan = results.dropna()
    
    credit_list.append(credit)
    
    FePO4_MSP_5th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.05))
    FePO4_MSP_50th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.5))
    FePO4_MSP_95th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.95))
    FePO4_MSP_mean.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].mean())
    
    FePO4_GWP_5th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.05))
    FePO4_GWP_50th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.5))
    FePO4_GWP_95th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.95))
    FePO4_GWP_mean.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].mean())
    
credit_result = pd.DataFrame()
credit_result['credit'] = credit_list
credit_result['FePO4_MSP_5th'] = FePO4_MSP_5th
credit_result['FePO4_MSP_50th'] = FePO4_MSP_50th
credit_result['FePO4_MSP_95th'] = FePO4_MSP_95th
credit_result['FePO4_MSP_mean'] = FePO4_MSP_mean
credit_result['FePO4_GWP_5th'] = FePO4_GWP_5th
credit_result['FePO4_GWP_50th'] = FePO4_GWP_50th
credit_result['FePO4_GWP_95th'] = FePO4_GWP_95th
credit_result['FePO4_GWP_mean'] = FePO4_GWP_mean

credit_result.to_excel(os.path.join(folder, f'results/FePO4_result_credit_{date.today()}.xlsx'),)

#%% for the simulation of N=1000 and HRT = 132 and food_sludge = 1 to the cost and CI (size of WRRF is a parameter)


size_list = []

FePO4_MSP_5th = []
FePO4_MSP_50th = []
FePO4_MSP_95th = []
FePO4_MSP_mean = []

FePO4_GWP_5th = []
FePO4_GWP_50th = []
FePO4_GWP_95th = []
FePO4_GWP_mean = []

sludge_cost_5th = []
sludge_cost_50th = []
sludge_cost_95th = []
sludge_cost_mean = []

sludge_GWP_5th = []
sludge_GWP_50th = []
sludge_GWP_95th = []
sludge_GWP_mean = []

for perspective in ['FePO4','sludge']:
    size_list = []

    for size in range(20, 320, 20):
        sys = create_system(temp_ratio=size, food_sludge_ratio=1, HRT=132, sludge_credit=300)
        sys.simulate()
        
        if perspective == 'FePO4':
            model = create_model(sys, perspective='FePO4')
        else:
            model = create_model(sys, perspective='sludge')
        
        kwargs = {'N':1000,'rule':'L','seed':3221}
        samples = model.sample(**kwargs)
        model.load_samples(samples)
        model.evaluate()
        
        idx = len(model.parameters)
        parameters = model.table.iloc[:, :idx]
        results = model.table.iloc[:, idx:]
        
        results_no_nan = results.dropna()
        
        size_list.append(size)
        
        if perspective == 'FePO4':
            FePO4_MSP_5th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.05))
            FePO4_MSP_50th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.5))
            FePO4_MSP_95th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.95))
            FePO4_MSP_mean.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].mean())
            
            FePO4_GWP_5th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.05))
            FePO4_GWP_50th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.5))
            FePO4_GWP_95th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.95))
            FePO4_GWP_mean.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].mean())
        
        if perspective == 'sludge':
            sludge_cost_5th.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].quantile(0.05))
            sludge_cost_50th.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].quantile(0.5))
            sludge_cost_95th.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].quantile(0.95))
            sludge_cost_mean.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].mean())
            
            sludge_GWP_5th.append(results_no_nan[('LCA','Sludge management GWP [kg_CO2_eq/tonne]')].quantile(0.05))
            sludge_GWP_50th.append(results_no_nan[('LCA','Sludge management GWP [kg_CO2_eq/tonne]')].quantile(0.5))
            sludge_GWP_95th.append(results_no_nan[('LCA','Sludge management GWP [kg_CO2_eq/tonne]')].quantile(0.95))
            sludge_GWP_mean.append(results_no_nan[('LCA','Sludge management GWP [kg_CO2_eq/tonne]')].mean())

    if perspective == 'FePO4':
        FePO4_result = pd.DataFrame()
        FePO4_result['size'] = size_list
        FePO4_result['FePO4_MSP_5th'] = FePO4_MSP_5th
        FePO4_result['FePO4_MSP_50th'] = FePO4_MSP_50th
        FePO4_result['FePO4_MSP_95th'] = FePO4_MSP_95th
        FePO4_result['FePO4_MSP_mean'] = FePO4_MSP_mean
        FePO4_result['FePO4_GWP_5th'] = FePO4_GWP_5th
        FePO4_result['FePO4_GWP_50th'] = FePO4_GWP_50th
        FePO4_result['FePO4_GWP_95th'] = FePO4_GWP_95th
        FePO4_result['FePO4_GWP_mean'] = FePO4_GWP_mean
    
    if perspective == 'sludge':
        sludge_result = pd.DataFrame()
        sludge_result['size'] = size_list
        sludge_result['sludge_cost_5th'] = sludge_cost_5th
        sludge_result['sludge_cost_50th'] = sludge_cost_50th
        sludge_result['sludge_cost_95th'] = sludge_cost_95th
        sludge_result['sludge_cost_mean'] = sludge_cost_mean
        sludge_result['sludge_GWP_5th'] = sludge_GWP_5th
        sludge_result['sludge_GWP_50th'] = sludge_GWP_50th
        sludge_result['sludge_GWP_95th'] = sludge_GWP_95th
        sludge_result['sludge_GWP_mean'] = sludge_GWP_mean

FePO4_result.to_excel(os.path.join(folder, f'results/FePO4_result_ratio_credit_{date.today()}.xlsx'))
sludge_result.to_excel(os.path.join(folder, f'results/sludge_result_ratio_credit_{date.today()}.xlsx'))


# model = create_model(sys)
# result = []
# result.to_excel(os.path.join(folder, f'results/FePO4_result_{date.today()}.xlsx'))

#%% for the simulation of N=1000 and HRT = 132 and food_sludge = 1 to the cost and CI (size of WRRF=10)
sys = create_system(temp_ratio=10, food_sludge_ratio=1, HRT=132, sludge_credit=300)
sys.simulate()
model = create_model(sys, perspective='sludge')
kwargs = {'N':1000,'rule':'L','seed':3221}
samples = model.sample(**kwargs)
model.load_samples(samples)
model.evaluate()
idx = len(model.parameters)
parameters = model.table.iloc[:, :idx]
results = model.table.iloc[:, idx:]
results.to_excel(os.path.join(folder, f'results/sludge_management_cost_CI_basline_{date.today()}.xlsx'))


#%% for the simulation of N=1000 and HRT = 132 and food_sludge_ratio = 1 to the avoided waste management CI analysis (sludge_CI_credit is a parameter)

credit_list = []

FePO4_MSP_5th = []
FePO4_MSP_50th = []
FePO4_MSP_95th = []
FePO4_MSP_mean = []

FePO4_GWP_5th = []
FePO4_GWP_50th = []
FePO4_GWP_95th = []
FePO4_GWP_mean = []

for credit in range(-300, 2100, 200):
    sys = create_system(temp_ratio=10, food_sludge_ratio=1, HRT=132, sludge_CI_credit=credit)
    model = create_model(sys, perspective='FePO4')
    
    kwargs = {'N':1000,'rule':'L','seed':3221}
    samples = model.sample(**kwargs)
    model.load_samples(samples)
    model.evaluate()
    
    idx = len(model.parameters)
    parameters = model.table.iloc[:, :idx]
    results = model.table.iloc[:, idx:]
    
    results_no_nan = results.dropna()
    
    credit_list.append(credit)
    
    FePO4_MSP_5th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.05))
    FePO4_MSP_50th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.5))
    FePO4_MSP_95th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.95))
    FePO4_MSP_mean.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].mean())
    
    FePO4_GWP_5th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.05))
    FePO4_GWP_50th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.5))
    FePO4_GWP_95th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.95))
    FePO4_GWP_mean.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].mean())
    
credit_result = pd.DataFrame()
credit_result['credit'] = credit_list
credit_result['FePO4_MSP_5th'] = FePO4_MSP_5th
credit_result['FePO4_MSP_50th'] = FePO4_MSP_50th
credit_result['FePO4_MSP_95th'] = FePO4_MSP_95th
credit_result['FePO4_MSP_mean'] = FePO4_MSP_mean
credit_result['FePO4_GWP_5th'] = FePO4_GWP_5th
credit_result['FePO4_GWP_50th'] = FePO4_GWP_50th
credit_result['FePO4_GWP_95th'] = FePO4_GWP_95th
credit_result['FePO4_GWP_mean'] = FePO4_GWP_mean

credit_result.to_excel(os.path.join(folder, f'results/FePO4_result_CI_credit_{date.today()}.xlsx'),)


#%% for the simulation of N=1000 and HRT = 132 and food_sludge = 1 to the cost and CI (size of WRRF is a parameter)


sys = create_system(temp_ratio=10, food_sludge_ratio=1, HRT=132, sludge_credit=300  )
model = create_model(sys, perspective='FePO4')
model.parameters
kwargs = {'N':1000,'rule':'L','seed':3221}
samples = model.sample(**kwargs)
model.load_samples(samples)
model.evaluate()
spearman_kwargs={'nan_policy': 'omit'}
r_df, p_df = qs.stats.get_correlations(model, kind='Spearman', **spearman_kwargs)

target_parameters = [
    'Food waste moisture [-]',
    'Org to gas [-]',
    'Org to vfa [-]',
    'Org to ethanol [-]',
    'Org to residue [-]',
    'Fe reduction [-]',
    'Acid dose [-]',
    'Oxidant excess [-]',
]

target_metrics = [
    'Fe recovery [-]',
    'P recovery [-]',
    'Fe po4 MSP [$/kg]',
    'Fe po4 GWP [kg_CO2_eq/kg]',
]

r_df_selected = r_df.loc[
    r_df.index.get_level_values(-1).isin(target_parameters),
    r_df.columns.get_level_values(-1).isin(target_metrics)
]

p_df_selected = p_df.loc[
    p_df.index.get_level_values(-1).isin(target_parameters),
    p_df.columns.get_level_values(-1).isin(target_metrics)
]

print(r_df_selected.shape)

r_df_selected.to_excel(
    os.path.join(folder, f'results/r_df_selected_FePO4_{date.today()}.xlsx')
)
p_df_selected.to_excel(
    os.path.join(folder, f'results/p_df_selected_FePO4_{date.today()}.xlsx')
)

#%% for the simulation of N=1000 and HRT = 132 and food_sludge = 1 to the cost and CI (size of WRRF is a parameter)

sys = create_system(temp_ratio=10, food_sludge_ratio=1, HRT=132, sludge_credit=300)
model = create_model(sys, perspective='sludge')
model.parameters
kwargs = {'N':1000,'rule':'L','seed':3221}
samples = model.sample(**kwargs)
model.load_samples(samples)
model.evaluate()
spearman_kwargs={'nan_policy':'omit'}
r_df, p_df = qs.stats.get_correlations(model, kind='Spearman', **spearman_kwargs)

target_parameters = [
    'Food waste moisture [-]',
    'Org to gas [-]',
    'Org to vfa [-]',
    'Org to ethanol [-]',
    'Org to residue [-]',
    'Fe reduction [-]',
    'Acid dose [-]',
    'Oxidant excess [-]',
]

target_metrics = [
    'Fe recovery [-]',
    'P recovery [-]',
    'Sludge management cost [$/tonne]',
    'Sludge management GWP [kg_CO2_eq/tonne]',
]

r_df_selected = r_df.loc[
    r_df.index.get_level_values(-1).isin(target_parameters),
    r_df.columns.get_level_values(-1).isin(target_metrics)
]

p_df_selected = p_df.loc[
    p_df.index.get_level_values(-1).isin(target_parameters),
    p_df.columns.get_level_values(-1).isin(target_metrics)
]

print(r_df_selected.shape)

r_df_selected.to_excel(
    os.path.join(folder, f'results/r_df_selected_sludge_{date.today()}.xlsx')
)
p_df_selected.to_excel(
    os.path.join(folder, f'results/p_df_selected_sludge_{date.today()}.xlsx')
)