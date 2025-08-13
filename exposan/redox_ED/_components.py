#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  9 10:09:19 2025

QSDSan

This module is developed by:
    Buai Shi <buaishi2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

#%%

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_components', )

def create_components(set_thermo=True):
    
    Na = Component('Na', phase='l', particle_size='Soluble',
                  degradability='Undegradable', organic=False)
    K = Component('K', phase='l', particle_size='Soluble',
                  degradability='Undegradable', organic=False)
    H2O = Component('H2O', phase='l', particle_size='Soluble',
                    degradability='Undegradable', organic=False)
    Cl = Component('Cl', phase='l', particle_size='Soluble',
                    degradability='Undegrdable', organic=False)
    Propionate = Component('Propionate', formula='C3H6O2', phase='l', particle_size='Soluble',
                           degradability='Readily', organic=True)
    Butyrate = Component('Butyrate', phase='l', particle_size='Soluble',
                         degradability='Readily', organic=True)
    Hexanoic = Component('Hexanoic', phase='l', particle_size='Soluble',
                         degradability='Readily', organic=True)
    cmps = Components((Na, K, H2O,Cl,Propionate,Butyrate))
    
    if set_thermo: qs_set_thermo(cmps)
    
    return cmps
