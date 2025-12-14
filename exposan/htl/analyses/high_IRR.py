#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

from exposan.htl import create_system, create_model, simulate_and_save

for feedstock in ['sludge','food','fogs','green','manure']:
    waste_cost_value=0
    waste_GWP_value=0
    print('\n\nfeedstock: ' + feedstock + '\ncost: ' + str(waste_cost_value) +\
          '\nGWP: ' + str(waste_GWP_value) + '\n\n')
    sys = create_system(waste_cost=waste_cost_value, waste_GWP=waste_GWP_value)
    model = create_model(sys, feedstock=feedstock)
    
    # if the independent variable is cost, uncomment the next three lines
    # simulate_and_save(model,
    #                   samples_kwargs={'N':1000, 'rule':'L', 'seed':3221},
    #                   notes=f'{feedstock}_high_IRR_cost_{waste_cost_value}')
    
    # if the independent variable is GWP, uncomment the next three lines
    # simulate_and_save(model,
    #                   samples_kwargs={'N':1000, 'rule':'L', 'seed':3221},
    #                   notes=f'{feedstock}_high_IRR_GWP_{waste_GWP_value}')