#!/usr/bin/env python3
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

# TODO: rerun all analyses and remake all figures

# TODO: add cost and CI breakdowns

#%% initialization

import os, pandas as pd, qsdsan as qs, numpy as np
from exposan.phos_rec import create_system, create_model
from datetime import date

folder = os.path.dirname(os.path.dirname(__file__))

#%% Fig. 3 - sludge management perspective, fermentation_time = 132, food_sludge = 1

sys = create_system(dry_solids_tonne_per_day=100, food_sludge_ratio=1, fermentation_time=132, perspective='sludge')
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

#%% Fig. 3 - FePO4 perspective, fermentation_time = 132, food_sludge_ratio = 1, avoided sludge management cost as credits

credit_list = []

FePO4_MSP_5th = []
FePO4_MSP_50th = []
FePO4_MSP_95th = []
FePO4_MSP_mean = []

# TODO: use 120 h or 132 h?
for credit in range(0, 1100, 100):
    sys = create_system(dry_solids_tonne_per_day=100, food_sludge_ratio=1, fermentation_time=132, sludge_cost_credit=credit, perspective='FePO4')
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
    
credit_result = pd.DataFrame()
credit_result['credit'] = credit_list
credit_result['FePO4_MSP_5th'] = FePO4_MSP_5th
credit_result['FePO4_MSP_50th'] = FePO4_MSP_50th
credit_result['FePO4_MSP_95th'] = FePO4_MSP_95th
credit_result['FePO4_MSP_mean'] = FePO4_MSP_mean

credit_result.to_excel(os.path.join(folder, f'results/FePO4_result_cost_credit_{date.today()}.xlsx'))

#%% Fig. 3 - FePO4 perspective, fermentation_time = 132, food_sludge_ratio = 1, avoided sludge management CI as credits

credit_list = []

FePO4_GWP_5th = []
FePO4_GWP_50th = []
FePO4_GWP_95th = []
FePO4_GWP_mean = []

for credit in range(0, 4500, 500):
    sys = create_system(dry_solids_tonne_per_day=100, food_sludge_ratio=1, fermentation_time=132, sludge_CI_credit=credit, perspective='FePO4')
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
    
    FePO4_GWP_5th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.05))
    FePO4_GWP_50th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.5))
    FePO4_GWP_95th.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].quantile(0.95))
    FePO4_GWP_mean.append(results_no_nan[('LCA','Fe po4 GWP [kg_CO2_eq/kg]')].mean())
    
credit_result = pd.DataFrame()
credit_result['credit'] = credit_list
credit_result['FePO4_GWP_5th'] = FePO4_GWP_5th
credit_result['FePO4_GWP_50th'] = FePO4_GWP_50th
credit_result['FePO4_GWP_95th'] = FePO4_GWP_95th
credit_result['FePO4_GWP_mean'] = FePO4_GWP_mean

credit_result.to_excel(os.path.join(folder, f'results/FePO4_result_CI_credit_{date.today()}.xlsx'),)

#%% Fig. 4 - heatmaps - food waste:sludge ratio and fermentation time

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
    fermentation_time_list = []
    for ratio in [0, 1/3, 2/3, 1, 4/3]:
        if ratio == 4/3:
            times = [132,]
        else:
            times = range(0, 144, 12)

        for time in times:
            # use try, since some ratio and reaction time combinations not feasible
            try:
                sys = create_system(dry_solids_tonne_per_day=100, food_sludge_ratio=ratio, fermentation_time=time, perspective=perspective)
            except RuntimeError:
                continue
            model = create_model(sys, perspective=perspective)
            
            kwargs = {'N':1000,'rule':'L','seed':3221}
            samples = model.sample(**kwargs)
            model.load_samples(samples)
            model.evaluate()
            
            idx = len(model.parameters)
            parameters = model.table.iloc[:, :idx]
            results = model.table.iloc[:, idx:]
            
            results_no_nan = results.dropna()
            
            ratio_list.append(ratio)
            fermentation_time_list.append(time)
            
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
        FePO4_result['fermentation_time'] = fermentation_time_list
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
        sludge_result['fermentation_time'] = fermentation_time_list
        sludge_result['sludge_cost_5th'] = sludge_cost_5th
        sludge_result['sludge_cost_50th'] = sludge_cost_50th
        sludge_result['sludge_cost_95th'] = sludge_cost_95th
        sludge_result['sludge_cost_mean'] = sludge_cost_mean
        sludge_result['sludge_GWP_5th'] = sludge_GWP_5th
        sludge_result['sludge_GWP_50th'] = sludge_GWP_50th
        sludge_result['sludge_GWP_95th'] = sludge_GWP_95th
        sludge_result['sludge_GWP_mean'] = sludge_GWP_mean
        
FePO4_result.to_excel(os.path.join(folder, f'results/decision_heatmap_FePO4_{date.today()}.xlsx'))
sludge_result.to_excel(os.path.join(folder, f'results/decision_heatmap_sludge_{date.today()}.xlsx'))

#%% sensitivity analysis, FePO4 perspective

sys = create_system(dry_solids_tonne_per_day=100, food_sludge_ratio=1, fermentation_time=132, perspective='FePO4')
model = create_model(sys, perspective='FePO4')
model.parameters
kwargs = {'N':1000,'rule':'L','seed':3221}
samples = model.sample(**kwargs)
model.load_samples(samples)

model.evaluate()
spearman_kwargs={'nan_policy': 'omit'}
r_df, p_df = qs.stats.get_correlations(model, kind='Spearman', **spearman_kwargs)

r_df.to_excel(os.path.join(folder, f'results/r_df_{date.today()}.xlsx'))
p_df.to_excel(os.path.join(folder, f'results/p_df_{date.today()}.xlsx'))

#%% Fig. 5 - heatmaps - VFA price and residual disposal cost

FePO4_MSP_5th = []
FePO4_MSP_50th = []
FePO4_MSP_95th = []
FePO4_MSP_mean = []

VFA_price_list = []
residue_cost_list = []

for VFA_price in np.linspace(0, 3/1000, 10):
    for residue_cost in np.linspace(0, 120/1000, 10):
        sys = create_system(dry_solids_tonne_per_day=100, food_sludge_ratio=1, fermentation_time=132, perspective='FePO4')
        sys.flowsheet.stream.precipitation_supernatant.price = VFA_price
        sys.flowsheet.stream.residue.price = -residue_cost
        
        model = create_model(sys, perspective='FePO4', exclude_context=True)
        
        kwargs = {'N':1000,'rule':'L','seed':3221}
        samples = model.sample(**kwargs)
        model.load_samples(samples)
        model.evaluate()
        
        idx = len(model.parameters)
        results = model.table.iloc[:, idx:].dropna()
                   
        VFA_price_list.append(VFA_price)
        residue_cost_list.append(residue_cost)
        
        FePO4_MSP_5th.append(results[('TEA','Fe po4 MSP [$/kg]')].quantile(0.05))
        FePO4_MSP_50th.append(results[('TEA','Fe po4 MSP [$/kg]')].quantile(0.5))
        FePO4_MSP_95th.append(results[('TEA','Fe po4 MSP [$/kg]')].quantile(0.95))
        FePO4_MSP_mean.append(results[('TEA','Fe po4 MSP [$/kg]')].mean())
        
heatmap_result = pd.DataFrame()
heatmap_result['VFA_price'] = VFA_price_list
heatmap_result['residue_cost'] = residue_cost_list
heatmap_result['FePO4_MSP_5th'] = FePO4_MSP_5th
heatmap_result['FePO4_MSP_50th'] = FePO4_MSP_50th
heatmap_result['FePO4_MSP_95th'] = FePO4_MSP_95th
heatmap_result['FePO4_MSP_mean'] = FePO4_MSP_mean

heatmap_result.to_excel(os.path.join(folder, f'results/context_heatmap_{date.today()}.xlsx'))

#%% SI - fermentation_time = 132, food_sludge = 1, different sizes

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

    for size in [1, 2, 4, 8, 16, 32, 64, 128, 256]:
        sys = create_system(dry_solids_tonne_per_day=size, food_sludge_ratio=1, fermentation_time=132, perspective=perspective)
        sys.simulate()
        model = create_model(sys, perspective=perspective)
        
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

FePO4_result.to_excel(os.path.join(folder, f'results/FePO4_result_size_{date.today()}.xlsx'))
sludge_result.to_excel(os.path.join(folder, f'results/sludge_result_size_{date.today()}.xlsx'))

#%% SI - different IRRs

IRR_list = []

FePO4_MSP_5th = []
FePO4_MSP_50th = []
FePO4_MSP_95th = []
FePO4_MSP_mean = []

sludge_cost_5th = []
sludge_cost_50th = []
sludge_cost_95th = []
sludge_cost_mean = []

for perspective in ['FePO4','sludge']:
    IRR_list = []
    
    if perspective == 'FePO4':
        IRR_range = np.arange(0.05, 0.2, 0.05)
    else:
        IRR_range = np.arange(0, 0.06, 0.01)

    for IRR in IRR_range:
        sys = create_system(dry_solids_tonne_per_day=100, food_sludge_ratio=1, fermentation_time=132, perspective=perspective)
        sys.TEA.IRR = IRR
        sys.simulate()
        
        model = create_model(sys, perspective=perspective, exclude_IRR=True)
        
        kwargs = {'N':1000,'rule':'L','seed':3221}
        samples = model.sample(**kwargs)
        model.load_samples(samples)
        model.evaluate()
        
        idx = len(model.parameters)
        parameters = model.table.iloc[:, :idx]
        results = model.table.iloc[:, idx:]
        
        results_no_nan = results.dropna()
        
        IRR_list.append(IRR)
        
        if perspective == 'FePO4':
            FePO4_MSP_5th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.05))
            FePO4_MSP_50th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.5))
            FePO4_MSP_95th.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].quantile(0.95))
            FePO4_MSP_mean.append(results_no_nan[('TEA','Fe po4 MSP [$/kg]')].mean())
        
        if perspective == 'sludge':
            sludge_cost_5th.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].quantile(0.05))
            sludge_cost_50th.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].quantile(0.5))
            sludge_cost_95th.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].quantile(0.95))
            sludge_cost_mean.append(results_no_nan[('TEA','Sludge management cost [$/tonne]')].mean())

    if perspective == 'FePO4':
        FePO4_result = pd.DataFrame()
        FePO4_result['IRR'] = IRR_list
        FePO4_result['FePO4_MSP_5th'] = FePO4_MSP_5th
        FePO4_result['FePO4_MSP_50th'] = FePO4_MSP_50th
        FePO4_result['FePO4_MSP_95th'] = FePO4_MSP_95th
        FePO4_result['FePO4_MSP_mean'] = FePO4_MSP_mean

    if perspective == 'sludge':
        sludge_result = pd.DataFrame()
        sludge_result['IRR'] = IRR_list
        sludge_result['sludge_cost_5th'] = sludge_cost_5th
        sludge_result['sludge_cost_50th'] = sludge_cost_50th
        sludge_result['sludge_cost_95th'] = sludge_cost_95th
        sludge_result['sludge_cost_mean'] = sludge_cost_mean


FePO4_result.to_excel(os.path.join(folder, f'results/FePO4_cost_IRR_{date.today()}.xlsx'))
sludge_result.to_excel(os.path.join(folder, f'results/sludge_cost_IRR_{date.today()}.xlsx'))