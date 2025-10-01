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
    
    NaCl = Component('NaCl', phase='l', particle_size='Soluble',
                  degradability='Undegradable', organic=False)
    KCl = Component('KCl', phase='l', particle_size='Soluble',
                  degradability='Undegradable', organic=False)
    H2O = Component('H2O', phase='l', particle_size='Soluble',
                    degradability='Undegradable', organic=False)
# !!! model ions separately in the future
    #Cl = Component('Cl', phase='l', particle_size='Soluble',
                    #degradability='Undegrdable', organic=False)
    Propionate = Component('Propionate', formula='C3H6O2', phase='l', particle_size='Soluble',
                           degradability='Readily', organic=True)
    Butyrate = Component('Butyrate', search_ID = '461-55-2', phase='l', particle_size='Soluble',
                         degradability='Readily', organic=True)
    add_V_from_rho(Butyrate, 960)
    Hexanoate = Component('Hexanoate', phase='l', particle_size='Soluble',
                          degradability='Readily', organic=True)
    add_V_from_rho(Hexanoate, 930)
    
    Digestate_solids = Component('Digestate_solids', formula='C5H7O2N', phase='s', particle_size='Particulate',
                                 degradability='Degradable', organic=True)
    # !!! Assuming digestate solids are undegradable because will leave system as sludge
    add_V_from_rho(Digestate_solids, 1400)
    # !!! Assuming average formula and density for digestate_solids, density unit: kg/m3
    Digestate_solids.Cn.add_model(1.25*10**3*Digestate_solids.MW/1000)
    # adding heat capacity of digestate_solids, reference: shijie leow et al. 2015 Green Chemistry
    cmps = Components((NaCl, KCl, H2O, Propionate, Butyrate, Hexanoate, Digestate_solids))
    
    
    if set_thermo: qs_set_thermo(cmps)
    
    return cmps
