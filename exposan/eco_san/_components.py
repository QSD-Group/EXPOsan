#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>
    Lewis Rowles <stetsonsc@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import Components, set_thermo as qs_set_thermo
from exposan.reclaimer import create_components as create_re_components

__all__ = ('create_components', )

def create_components(set_thermo=True):
    re_cmps = create_re_components(set_thermo=False)

    # Reuse components in the bwaise module for consistency
    not_used_cmps = ('KCl', 'GAC', 'LPG', 'Zeolite')
    cmps = Components((cmp for cmp in re_cmps if cmp.ID not in not_used_cmps))

    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_synonym('H2O', 'Water')
    cmps.set_synonym('MgOH2', 'MagnesiumHydroxide')
    cmps.set_synonym('NaCl', 'SodiumChloride')
    cmps.set_synonym('HCl', 'HydrogenChloride')

    if set_thermo: qs_set_thermo(cmps)

    return cmps