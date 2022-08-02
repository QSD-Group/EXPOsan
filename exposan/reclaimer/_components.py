#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''


# %%

from qsdsan import Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho
from exposan.bwaise import create_components as create_bw_components

__all__ = ('create_components', )


def create_components(set_thermo=True):
    bw_cmps = create_bw_components(set_thermo=False)

    KCl = Component('KCl', phase='s', particle_size='Soluble',
                    degradability='Undegradable', organic=False)

    MgOH2 = Component('MgOH2', phase='s', particle_size='Particulate',
                      degradability='Undegradable', organic=False)
    MgOH2.default()

    NaCl = Component('NaCl', phase='s', particle_size='Particulate', 
                     degradability='Undegradable', organic=False)

    HCl = Component('HCl', phase='s', particle_size='Particulate',
                    degradability='Undegradable', organic=False)


    GAC = Component('GAC', search_ID='C', phase='s', particle_size='Particulate',
                    degradability='Undegradable', organic=True)
    # 485 kg/m3 is average from:
    # https://www.systematixusa.com/products/media/active_media/gac.htm (accessed 2022-07-25)
    add_V_from_rho(GAC, 485)

    LPG = Component('LPG', search_ID='PubChem=6334', formula='CH3CH2CH3',
                    phase='g', particle_size='Dissolved gas',
                    degradability='Slowly', organic=True)

    Zeolite = Component('Zeolite', search_ID='PubChem=9942228', formula='Na2Al2Si2O8',
                        phase='s', particle_size='Particulate',
                        degradability='Undegradable', organic=False)
    # Assume the same density of water, may or may not be good,
    # densities can vary from ~0.5->2 kg/m3 depending on the Zeolite type
    Zeolite.copy_models_from(bw_cmps.H2O, ('V',))

    # Reuse components in the bwaise module for consistency
    cmps = Components((*bw_cmps, KCl, MgOH2, NaCl, HCl, GAC, LPG, Zeolite))

    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)
    
    cmps.compile()
    
    cmps.set_alias('H2O', 'Water')
    cmps.set_alias('KCl', 'PotassiumChloride')
    cmps.set_alias('MgOH2', 'MagnesiumHydroxide')
    cmps.set_alias('NaCl', 'SodiumChloride')
    cmps.set_alias('HCl', 'HydrogenChloride')    
    
    if set_thermo: qs_set_thermo(cmps)

    return cmps
