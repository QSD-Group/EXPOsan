#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# Just to run the uncertainty analysis for sysA and sysB
from exposan import biogenic_refinery as br

m = br.models
modelA = m.modelA
modelB = m.modelB

def run(models, seed=None, N=1000):
    try: iter(models)
    except: models = (models,) # make it possible to run with one or more models
    for model in models:
        m.run_uncertainty(model, seed=seed, N=N)
        m.save_uncertainty_results(model)


if __name__ == '__main__':
    run((modelA, modelB), seed=5, N=1000)