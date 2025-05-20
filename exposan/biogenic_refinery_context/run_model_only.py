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

folder = os.path.dirname(__file__)
# sys.path.insert(0, '/Users/stetsonrowles/Dropbox/Mac (3)/Documents/GitHub/EXPOsan/exposan/')

from exposan.biogenic_refinery_context import models, results_path
model = models.create_modelA()

path = os.path.join(folder, 'Flow + Biosolids Combined Data (Major-Facilities) (1).xlsx')
df = pd.read_excel(path)

biosolids_data = df[['NPDES ID2', 'Amount of Biosolids Generated']]

# Prepare a list to collect results
results = []

# ===== Loop through facilities =====



    # Simulate system with updated input
model.system.simulate()
biochar = model.system.flowsheet.stream.biochar
print('Biochar flow [kg/hr]:', biochar.F_mass)
print('Carbon in biochar:', biochar.imass['C'])

    # ===== Monte Carlo sampling =====
mpath = os.path.join (results_path,'biosolids_base_scenario_results.xlsx')
N = 100  # Number of Latin Hypercube samples
samples = model.sample(N, rule='L')
model.load_samples(samples)
model.evaluate()
for metric in model.metrics:
    print(f"{metric.name}: {metric()}")  # Will call the lambda
print(model.system.flowsheet.stream.biochar)
model.table.to_excel(mpath)

    # ===== Collect metrics =====
metric_values = {}
for metric in model.metrics:
    for col in model.table.columns:
        if metric.name in col[-1]:  # Match just the name part
            values = model.table[col]
            break
    else:
        raise KeyError(f"Metric '{metric.name}' not found in model.table.")
    
    mean = np.mean(values)
    lower = np.percentile(values, 5)
    upper = np.percentile(values, 95)
    metric_values[metric.name] = (mean, lower, upper)

#     # Store results for this facility
# result = {'Facility': facility, 'Biosolids_tpy': biosolids_tpy}
# for name, (mean, lower, upper) in metric_values.items():
#     result[f'{name} (mean)'] = mean
#     result[f'{name} (5th pct)'] = lower
#     result[f'{name} (95th pct)'] = upper

# results.append(result)

# ===== Save results =====
mpath = os.path.join (results_path,'biosolids_scenario_results.csv')
results_df = pd.DataFrame(results)
results_df.to_csv(mpath, index=False)

print('\nAll scenarios completed. Results saved to biosolids_scenario_results.csv.')