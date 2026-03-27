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

import os, pandas as pd
from exposan.phos_rec import create_system, create_model
from datetime import date

folder = os.path.dirname(os.path.dirname(__file__))

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