# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from qsdsan.utils import ospath, load_data
from exposan.metab import create_system, create_model, run_model, data_path
# import numpy as np, pandas as pd

#%%

def run_discrete_DVs(samples_path):
    dct = load_data(samples_path, sheet=None)
    for n in (1, 2):
        for i in ('UASB', 'FB', 'PB'):
        # for i in ('UASB',):
            for j in 'PHM':
                sys = create_system(n_stages=n, reactor_type=i, gas_extraction=j)
                mdl = create_model(sys, kind='DV')
                sample = dct[i+'_'+j].to_numpy()
                run_model(mdl, sample)
        
#%%
if __name__ == '__main__':
    path = ospath.join(data_path, 'analysis_framework.xlsx')
    run_discrete_DVs(path)
