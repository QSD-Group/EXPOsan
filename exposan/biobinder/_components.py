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

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
# from exposan.utils import add_V_from_rho
from exposan import htl



__all__ = ('create_components',)

def create_components(set_thermo=True):
    htl_cmps = htl.create_components()
    biobinder_cmps = list(htl_cmps)
    
    replace_dct = {
        'Lipids': 'Sludge_lipid',
        'Proteins': 'Sludge_protein',
        'Carbohydrates': 'Sludge_carbo',
        'Ash': 'Sludge_ash',
        }
    
    for new_ID, old_ID in replace_dct.items():
        old_cmp = htl_cmps[old_ID]
        new_cmp = old_cmp.copy(new_ID)
        biobinder_cmps.remove(old_cmp)
        biobinder_cmps.append(new_cmp)

    
    biobinder_cmps = Components(biobinder_cmps)
    
    # for i in cmps:
    #     for attr in ('HHV', 'LHV', 'Hf'):
    #         if getattr(i, attr) is None: setattr(i, attr, 0)

    biobinder_cmps.compile()
    biobinder_cmps.set_alias('H2O', 'Water')
    biobinder_cmps.set_alias('Carbohydrates', 'Carbs')
    biobinder_cmps.set_alias('C', 'Carbon')
    biobinder_cmps.set_alias('N', 'Nitrogen')
    biobinder_cmps.set_alias('P', 'Phosphorus')

    if set_thermo: qs_set_thermo(biobinder_cmps)

    return biobinder_cmps