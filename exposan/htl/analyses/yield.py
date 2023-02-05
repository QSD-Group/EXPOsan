#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Jianan Feng <jiananf2@illinois.edu>
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References:

(1) Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
    Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
    Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
    Process Design and Economics for the Conversion of Algal Biomass to
    Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
    PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    
(2) Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels via
    Liquefaction - Hydrothermal Liquefaction Reactor Design: April 5, 2013;
    NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462, 1111191.
    https://doi.org/10.2172/1111191.
'''

import os, pandas as pd
from exposan.htl import results_path, create_model

yield_results_dict = {'CHG_yield':[], 'HT_yield':[],
                        'TEA_5th':[], 'TEA_50th':[], 'TEA_95th':[],
                        'sludge_price_5th':[], 'sludge_price_50th':[], 'sludge_price_95th':[],
                        'LCA_diesel_5th':[], 'LCA_diesel_50th':[], 'LCA_diesel_95th':[],
                        'LCA_sludge_5th':[], 'LCA_sludge_50th':[], 'LCA_sludge_95th':[]}

yield_results = pd.DataFrame(yield_results_dict)
CHG_yields = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
HT_yields = (0.78, 0.82, 0.86, 0.90, 0.94, 0.98)

model = create_model('baseline', exclude_yield=True, key_metrics_only=True)
N = 50
samples = model.sample(N=N, rule='L', seed=3221)
model.load_samples(samples)

sys = model.system
tea = sys.TEA
CHG = sys.flowsheet.unit.CHG
HT = sys.flowsheet.unit.HT

get_quantiles = lambda data, quantiles=(0.05, 0.5, 0.95): [data.quantile(q) for q in quantiles]

for CHG_yield in CHG_yields:
    for HT_yield in HT_yields:
        CHG.gas_C_2_total_C = CHG_yield
        HT.hydrocarbon_ratio = HT_yield
        print('\n\n', f'CHG_yield: {CHG_yield}\n', f'HT_yield: {HT_yield}', '\n\n')
        model.evaluate()            
        MFSP = model.table['TEA']['MFSP [$/GGE]'].dropna()
        sludge_price = model.table['TEA']['sludge_management_price [$/ton dry sludge]'].dropna()
        LCA_diesel = model.table['LCA']['GWP_diesel [g CO2/MMBTU diesel]'].dropna()
        LCA_sludge = model.table['LCA']['GWP_sludge [kg CO2/ton dry sludge]'].dropna()
            
        yield_results.loc[len(yield_results.index)] = (
            [CHG_yield, HT_yield,] +
            get_quantiles(MFSP) +
            get_quantiles(sludge_price) +
            get_quantiles(LCA_diesel) +
            get_quantiles(LCA_sludge)
            )

yield_results.to_excel(os.path.join(results_path, f'_yield_{N}.xlsx'))