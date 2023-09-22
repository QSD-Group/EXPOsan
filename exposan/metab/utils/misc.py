# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 12:55:58 2023

@author: joy_c
"""

from math import ceil

__all__ = ('add_prefix',
           '_construct_water_pump',
           '_construct_vacuum_pump',
           '_F_mass',
           '_F_vol')

add_prefix = lambda dct, prefix: {f'{prefix} - {k}':v for k,v in dct.items()}

def _construct_water_pump(pump):
    hp = pump.design_results['Power'] * pump.parallel['self']
    if hp > 29.5: 
        q22 = hp/29.5
        q40 = 0    
    else:
        q22 = 0
        q40 = ceil(hp/0.05364)
    return q22, q40

def _construct_vacuum_pump(pump):
    # kW = pump.design_results['Power'] * pump.parallel['self'] * 0.7457
    kW = pump.power_utility.consumption
    return (kW/4)**0.6

def _F_mass(ws):
    cmps = ws.components
    return sum(ws.mass * cmps.i_mass)

def _F_vol(ws):
    cmps = ws.components
    mol = sum(ws.mass * cmps.i_mass / cmps.MW) * 1e3 # mol/hr
    Q = mol * 8.314 * ws.T / ws.P  # m3/hr
    return Q