#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Lewis Rowles <stetsonsc@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# Run uncertainty analysis and Spearman without country-specific settings
from exposan.biogenic_refinery import create_model, run_uncertainty, save_uncertainty_results

def run(model_IDs, seed=None, N=1000, country_specific=False):
    # Make it possible to run with one or more models
    if isinstance(model_IDs, str): model_IDs = (model_IDs, )
    for ID in model_IDs:
        model = create_model(ID)
        run_uncertainty(model, seed=seed, N=N, country_specific=country_specific)
        save_uncertainty_results(model)


if __name__ == '__main__':
    run(('A', 'B'), seed=5, N=1000)