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
'''

import numpy as np, pandas as pd, qsdsan as qs
from exposan.htl import create_model

model = create_model()
 
np.random.seed(3221)
samples = model.sample(N=1000, rule='L')
model.load_samples(samples)
model.evaluate()
model.table

#%%
def organize_results(model, path):
    idx = len(model.parameters)
    parameters = model.table.iloc[:, :idx]
    results = model.table.iloc[:, idx:]
    percentiles = results.quantile([0, 0.05, 0.25, 0.5, 0.75, 0.95, 1])
    with pd.ExcelWriter(path) as writer:
        parameters.to_excel(writer, sheet_name='Parameters')
        results.to_excel(writer, sheet_name='Results')
        percentiles.to_excel(writer, sheet_name='Percentiles')
organize_results(model, 'example_model.xlsx')

#%%
fig, ax = qs.stats.plot_uncertainties(model)
fig

#%%
fig, ax = qs.stats.plot_uncertainties(model, x_axis=model.metrics[0], y_axis=model.metrics[1],
                                      kind='kde-kde', center_kws={'fill': True})
fig

#%%
r_df, p_df = qs.stats.get_correlations(model, kind='Spearman', nan_policy='omit')
fig, ax = qs.stats.plot_correlations(r_df)
fig