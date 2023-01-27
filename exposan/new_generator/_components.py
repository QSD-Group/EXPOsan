#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>
    Shion Watabe <shionwatabe@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''


# %%

from qsdsan import Component, Components, set_thermo as qs_set_thermo
from exposan.reclaimer import create_components as create_re_components

__all__ = ('create_components',)

def create_components(set_thermo=True):
    re_cmps = create_re_components(set_thermo=False)

    C = Component('C', phase='l', particle_size='Soluble',
                  degradability='Undegradable', organic=False)

    SolubleCH4 = Component('SolubleCH4', search_ID='CH4', phase='l', particle_size='Soluble',
                           degradability='Slowly', organic=True)

    SO2 = Component('SO2', phase='g', particle_size='Dissolved',
                    degradability='Undegradable', organic=False)

    NaOH = Component('SodiumHydroxide', search_ID='SodiumHydroxide', formula='NaOH',
                     phase='s', particle_size='Particulate',
                     degradability='Slowly', organic=False)

    # Reuse components in the bwaise module for consistency
    cmps = Components(
        (*(cmp for cmp in re_cmps if cmp not in (re_cmps.KCl, re_cmps.MgOH2)),
        C, SolubleCH4, SO2, NaOH,)
        )

    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_alias('H2O', 'Water')
    if set_thermo: qs_set_thermo(cmps)

    return cmps