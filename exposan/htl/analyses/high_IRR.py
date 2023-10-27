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
    sys = create_system(waste_price=waste_price_value, waste_GWP=waste_GWP_value)
    model = create_model(sys, feedstock=feedstock)
    simulate_and_save(model, samples_kwargs={'N':1000, 'rule':'L', 'seed':3221}, notes=f'high_IRR_{feedstock}_cost_{waste_price_value}')
    # simulate_and_save(model, samples_kwargs={'N':1000, 'rule':'L', 'seed':3221}, notes=f'high_IRR_{feedstock}_GWP_{waste_GWP_value}')