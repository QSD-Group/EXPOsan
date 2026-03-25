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

import matplotlib.pyplot as plt
from exposan.phos_rec import create_system, create_model, simulate_and_save

for perspective in ['FePO4','sludge']:
    for ratio in [0, 1/3, 2/3, 1, 4/3]:
        if ratio == 4/3:
            sys = create_system(temp_ratio=10, food_sludge_ratio=ratio,HRT=132)
        
            model = create_model(sys, perspective=perspective)

            simulate_and_save(model, notes=perspective+'_'+str(round(ratio, 2))+'_132')
        else:
            for HRT in [0,12,24,36,48,60,72,84,96,108,120,132]:
                sys = create_system(temp_ratio=10, food_sludge_ratio=ratio,HRT=HRT)
            
                model = create_model(sys, perspective=perspective)

                simulate_and_save(model, notes=perspective+'_'+str(round(ratio, 2))+'_'+str(HRT))