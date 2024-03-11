#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References:

'''

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_components',)

def create_components(set_thermo=True):
    
    # bauxite
    Bauxite = Component(ID='Bauxite',
                        phase='s',
                        particle_size='Particulate',
                        formula='',
                        degradability='Undegradable',
                        organic=False)
    
    # Al(OH)3
    AlOH3 = Component(ID='AlOH3',
                      phase='s',
                      particle_size='Particulate',
                      degradability='Undegradable',
                      organic=False)
    
    # formatic acid
    HCOOH = Component(ID='HCOOH',
                      phase='l',
                      particle_size='Soluble',
                      degradability='Readily',
                      organic=True)
    
    # ALF 7360-53-4
    ALF = Component(ID='C3H3AlO6',
                    search_ID='7360-53-4',
                    particle_size='Particulate',
                    degradability='Undegradable',
                    organic=False)
    
    # Al(HCOO)2(OH) 
    AlHCOO2OH = Component(ID='AlHCOO2OH',
                          search_ID='',
                          particle_size='Particulate',
                          degradability='Undegradable',
                          organic=False)
    
    H2O = Component(ID='H2O',
                    particle_size='Soluble',
                    degradability='Undegradable',
                    organic=False)
    
    O2 = Component(ID='O2',
                   phase='g',
                   particle_size='Dissolved gas',
                   degradability='Undegradable',
                   organic=False)
    
    N2 = Component(ID='N2',
                   phase='g',
                   particle_size='Dissolved gas',
                   degradability='Undegradable',
                   organic=False)
    
    CH4 = Component(ID='CH4',
                    phase='g',
                    particle_size='Dissolved gas',
                    degradability='Slowly',
                    organic=True)
    
    CO2 = Component(ID='CO2',
                    phase='g',
                    particle_size='Dissolved gas',
                    degradability='Undegradable',
                    organic=False)
    
    # values for Membrane are made up since we just need to calculate cost for membrane
    Membrane = Component('Membrane', phase='s', particle_size='Particulate',
                             degradability='Undegradable', organic=False)
    add_V_from_rho(Membrane, 1500)
    Membrane.copy_models_from(Chemical('CaCO3'),('Cn',))
    
    cmps = Components([Bauxite, AlOH3, HCOOH, ALF,
                       AlHCOO2OH, H2O, O2, N2,
                       CH4, CO2])
    
    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_alias('H2O', 'Water')
    if set_thermo: qs_set_thermo(cmps)

    return cmps
    
    # examples:
    Sludge_lipid = Component('Sludge_lipid', phase='s',
                              particle_size='Particulate',
                              formula='C56H95O24N9P',
                              degradability='Undegradable',
                              organic=False)
    add_V_from_rho(Sludge_lipid, 1400)
    Sludge_lipid.HHV = 22.0*10**6*Sludge_lipid.MW/1000
    Sludge_lipid.Cn.add_model(1.25*10**3*Sludge_lipid.MW/1000) 
    Sludge_lipid.mu.add_model(6000)
    
    H2SO4 = Component('H2SO4', phase='l', particle_size='Soluble',
                      degradability='Undegradable', organic=False)
    

    
    NaOH = Component('NaOH', phase='l', particle_size='Soluble',
                     degradability='Undegradable', organic=False)