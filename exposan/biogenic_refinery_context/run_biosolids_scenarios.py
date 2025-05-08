#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Lewis Rowles <stetsonsc@gmail.com>

    Aaron Marszewski <aaronpm3@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


import sys
import os
import pandas as pd
import numpy as np
import qsdsan as qs
from exposan import biogenic_refinery_context as brc


sys.path.insert(0, '/Users/stetsonrowles/Dropbox/Mac (3)/Documents/GitHub/EXPOsan/exposan/')

from exposan.biogenic_refinery_context import models
model = models.create_modelA()



df = pd.read_excel('Flow + Biosolids Combined Data (Major-Facilities) (1).xlsx')

biosolids_data = df[['NPDES ID2', 'Amount of Biosolids Generated']]

# Prepare a list to collect results
results = []

# ===== Loop through facilities =====
for index, row in biosolids_data.iterrows():
    facility = row['NPDES ID2']
    biosolids_tpy = row['Amount of Biosolids Generated']


    # Update biosolids stream
    biosolids_stream = model.system.get_stream('biosolids')
    biosolids_stream.F_mass = biosolids_tpy * 1000  # Convert tons to kg

    # Simulate system with updated input
    model.system.simulate()

    # ===== Monte Carlo sampling =====
    N = 1000  # Number of Latin Hypercube samples
    model.sample(N, rule='LHS')
    model.evaluate()

    # ===== Collect metrics =====
    metric_values = {}
    for metric in model.metrics:
        values = model.table[metric.name]
        mean = np.mean(values)
        lower = np.percentile(values, 5)
        upper = np.percentile(values, 95)
        metric_values[metric.name] = (mean, lower, upper)

    # Store results for this facility
    result = {'Facility': facility, 'Biosolids_tpy': biosolids_tpy}
    for name, (mean, lower, upper) in metric_values.items():
        result[f'{name} (mean)'] = mean
        result[f'{name} (5th pct)'] = lower
        result[f'{name} (95th pct)'] = upper

    results.append(result)

# ===== Save results =====
results_df = pd.DataFrame(results)
results_df.to_csv('biosolids_scenario_results.csv', index=False)

print('\nAll scenarios completed. Results saved to biosolids_scenario_results.csv.')
