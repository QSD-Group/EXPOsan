#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import Components, set_thermo as qs_set_thermo
from exposan.reclaimer import create_components as create_re_components

__all__ = ('create_components',)

def create_components(set_thermo=True):
    #Reuse components in the bwaise and reclaimer for consistency
    re_cmps = create_re_components(set_thermo = False)
    #remove components that are not used
    not_used_cmps = ('KCl','GAC','Zeolite')
    cmps = Components(cmp for cmp in re_cmps if cmp.id not in not_used_cmps)
    
    #assign value 0 to heat value that does not exist
    for cmp in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(cmp, attr) is None: 
                setattr(cmp, attr, 0)

    cmps.compile()
    
    cmps.set_alias('H2O','Water')
    
    if set_thermo: qs_set_thermo(cmps)
    
    return cmps
    
    
    

