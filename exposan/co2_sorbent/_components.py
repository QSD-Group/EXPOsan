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

import qsdsan as qs
from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_components',)

def create_components(set_thermo=True):
    
    # bauxite
    Bauxite = Component(ID='Bauxite',
                        search_ID='Bauxite',
                        phase='s',
                        particle_size='Particulate',
                        degradability='Undegradable',
                        organic=False)
    # https://www.standardboksit.com.tr/en/Arge/Arge/8747#:~:text=Bauxite%20is%20a%20mixture%20of,2.5%2D3.5%20g%2Fcm3.
    add_V_from_rho(Bauxite, 3000)

    
    # Al(OH)3
    AlH3O3 = Component(ID='AlH3O3',
                      search_ID='21645-51-2',
                      phase='s',
                      particle_size='Particulate',
                      degradability='Undegradable',
                      organic=False)
    
    # ALF
    ALF = Component(ID='C3H3AlO6',
                    search_ID='7360-53-4',
                    phase='s',
                    particle_size='Particulate',
                    degradability='Undegradable',
                    organic=False)
    add_V_from_rho(ALF, 1441) # Evans et al. Science Advances SI. Table S1
    
    # https://www.chemsrc.com/en/cas/7360-53-4_311690.html
    ALF.Tm = 8.4+273.15
    # ALF.Tb = 100.6+273.15
    
    # assume similar to Al
    # https://link.springer.com/article/10.1023/B:JMSC.0000048735.50256.96
    ALF.mu.add_model(1000)
    
    # formatic acid
    HCOOH = Component(ID='HCOOH',
                      search_ID='64-18-6',
                      particle_size='Soluble',
                      degradability='Readily',
                      organic=True)
    
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
    
    H2 = Component(ID='H2',
                   phase='g',
                   particle_size='Dissolved gas',
                   degradability='Undegradable',
                   organic=False)
    
    CH4 = Component(ID='CH4',
                    phase='g',
                    particle_size='Dissolved gas',
                    degradability='Slowly',
                    organic=True)
    
    C2H4 = Component(ID='C2H4',
                     phase='g',
                     particle_size='Dissolved gas',
                     degradability='Slowly',
                     organic=True)
    
    CO2 = Component(ID='CO2',
                    phase='g',
                    particle_size='Dissolved gas',
                    degradability='Undegradable',
                    organic=False)
    
    CO = Component(ID='CO',
                   phase='g',
                   particle_size='Dissolved gas',
                   degradability='Undegradable',
                   organic=False)
    
    # methanol
    CH4O = Component(ID='CH4O',
                     search_ID='67-56-1',
                     particle_size='Soluble',
                     degradability='Readily',
                     organic=True)
    
    # ethanol
    C2H6O = Component(ID='C2H6O',
                      search_ID='64-17-5',
                      particle_size='Soluble',
                      degradability='Readily',
                      organic=True)
    
    # n-propanol
    C3H8O = Component(ID='C3H8O',
                      search_ID='71-23-8',
                      particle_size='Soluble',
                      degradability='Readily',
                      organic=True)
    
    cmps = Components([Bauxite, AlH3O3, ALF, HCOOH, H2O,
                       O2, N2, H2, CH4, CO2, CO, C2H6O,
                       C2H4, CH4O, C3H8O])
    
    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_alias('H2O', 'Water')
    if set_thermo: qs_set_thermo(cmps)

    return cmps