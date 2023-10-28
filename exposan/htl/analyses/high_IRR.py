#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 08:45:48 2023

@author: jiananfeng
"""

from exposan.htl import create_system, create_model, simulate_and_save

for feedstock in ['sludge','food','fogs','green','manure']:
    waste_price_value=0
    waste_GWP_value=0
    print('\n\nfeedstock: ' + feedstock + '\ncost: ' + str(waste_price_value) +\
          '\nGWP: ' + str(waste_GWP_value) + '\n\n')
    sys = create_system(waste_price=waste_price_value, waste_GWP=waste_GWP_value)
    model = create_model(sys, feedstock=feedstock)
    
    # if the independent variable is cost, uncomment the next three lines
    # simulate_and_save(model,
    #                   samples_kwargs={'N':1000, 'rule':'L', 'seed':3221},
    #                   notes=f'{feedstock}_high_IRR_cost_{waste_price_value}')
    
    # if the independent variable is GWP, uncomment the next three lines
    # simulate_and_save(model,
    #                   samples_kwargs={'N':1000, 'rule':'L', 'seed':3221},
    #                   notes=f'{feedstock}_high_IRR_GWP_{waste_GWP_value}')