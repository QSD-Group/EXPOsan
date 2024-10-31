#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import Components, set_thermo as qs_set_thermo
from exposan.saf import create_components as create_saf_components

__all__ = ('create_components',)

def create_components(set_thermo=True):
    saf_cmps = create_saf_components(set_thermo=False)
    biobinder_cmps = Components([i for i in saf_cmps])

    biobinder_cmps.compile()
    biobinder_cmps.set_alias('H2O', 'Water')
    biobinder_cmps.set_alias('H2O', '7732-18-5')
    biobinder_cmps.set_alias('Carbohydrates', 'Carbs')
    biobinder_cmps.set_alias('C', 'Carbon')
    biobinder_cmps.set_alias('N', 'Nitrogen')
    biobinder_cmps.set_alias('P', 'Phosphorus')
    biobinder_cmps.set_alias('K', 'Potassium')
    biobinder_cmps.set_alias('C16H34', 'Biofuel') # Tb = 559 K
    biobinder_cmps.set_alias('TRICOSANE', 'Biobinder') # Tb = 654 K

    if set_thermo: qs_set_thermo(biobinder_cmps)

    return biobinder_cmps