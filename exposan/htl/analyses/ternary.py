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
from exposan.htl import results_path, create_system, create_model
from datetime import date

ternary_results_dict = {'lipid':[], 'protein':[], 'carbohydrate':[],
                        'MDSP_5th':[], 'MDSP_50th':[], 'MDSP_95th':[],
                        'sludge_price_5th':[], 'sludge_price_50th':[], 'sludge_price_95th':[],
                        'GWP_diesel_5th':[], 'GWP_diesel_50th':[], 'GWP_diesel_95th':[],
                        'GWP_sludge_5th':[], 'GWP_sludge_50th':[], 'GWP_sludge_95th':[]}

ternary_results = pd.DataFrame(ternary_results_dict)
lipids = (0, )
proteins = (0, )

get_quantiles = lambda data, quantiles=(0.05, 0.5, 0.95): [data.quantile(q) for q in quantiles]

for lipid in lipids:
    for protein in proteins:
        if lipid + protein <= 1:
            sys = create_system()
            WWTP = sys.flowsheet.unit.WWTP
            WWTP.sludge_afdw_lipid = lipid
            WWTP.sludge_afdw_protein = protein
            sys.simulate()
            print('\n\n', f'Lipid: {WWTP.sludge_afdw_lipid}\n', f'Protein: {WWTP.sludge_afdw_protein}\n', f'Carbohydrate: {WWTP.sludge_afdw_carbo}', '\n\n')
            model = create_model(sys, exclude_sludge_compositions=True,
                                 include_HTL_yield_as_metrics=False,
                                 include_other_metrics=False,
                                 include_other_CFs_as_metrics=False,
                                 include_check=False)
            N = 1000
            samples = model.sample(N=N, rule='L', seed=3221)
            model.load_samples(samples)
            model.evaluate()            
            MDSP = model.table['TEA']['MDSP [$/gal diesel]'].dropna()
            sludge_price = model.table['TEA']['Sludge management price [$/ton dry sludge]'].dropna()
            LCA_diesel = model.table['LCA']['GWP diesel [kg CO2/MMBTU diesel]'].dropna()
            LCA_sludge = model.table['LCA']['GWP sludge [kg CO2/ton dry sludge]'].dropna()
            
            ternary_results.loc[len(ternary_results.index)] = (
                [100*lipid, 100*protein, round(100-100*lipid-100*protein),] +
                get_quantiles(MDSP) +
                get_quantiles(sludge_price) +
                get_quantiles(LCA_diesel) +
                get_quantiles(LCA_sludge)
                )

ternary_results.to_excel(os.path.join(results_path, f'{lipids}_{proteins}_{date.today()}_ternary_{N}.xlsx'))