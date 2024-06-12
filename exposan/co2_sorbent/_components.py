#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_components',)

def create_components(set_thermo=True):
    
    # Al2O3
    Al2O3 = Component(ID='Al2O3',
                      search_ID='1344-28-1',
                      particle_size='Particulate',
                      degradability='Undegradable',
                      organic=False)
    
    # SiO2
    SiO2 = Component(ID='SiO2',
                     search_ID='7631-86-9',
                     particle_size='Particulate',
                     degradability='Undegradable',
                     organic=False)
    
    # Fe2O3
    Fe2O3 = Component(ID='Fe2O3',
                      search_ID='1309-37-1',
                      particle_size='Particulate',
                      degradability='Undegradable',
                      organic=False)
    
    # Fe
    Fe = Component(ID='Fe',
                   search_ID='7439-89-6',
                   particle_size='Particulate',
                   degradability='Undegradable',
                   organic=False)
    
    # Al(OH)3
    AlH3O3 = Component(ID='AlH3O3',
                       search_ID='21645-51-2',
                       particle_size='Particulate',
                       degradability='Undegradable',
                       organic=False)
    
    # ALF
    # ALF is typical solid at relevant (temperature and pressure) conditions
    # set the phase of ALF as 's', otherwise, it is missing key thermodynamic properties
    # also tested: changing the phase to 'l' does not change the results too much for all systems
    ALF = Component(ID='C3H3AlO6',
                    search_ID='7360-53-4',
                    phase='s',
                    particle_size='Particulate',
                    degradability='Undegradable',
                    organic=False)
    add_V_from_rho(ALF, 1441) # Evans et al. Science Advances SI. Table S1
    # https://www.harrellindustries.com/wp-content/uploads/2015/10/0200-Aluminum-Formate-SDS-US.pdf (accessed 05-30-2024)
    ALF.Tm = 660+273.15 # normal melting temperature [K]
    # no information on ALF.mu
    # assume it as 1 (changing it has minimal effect on TEA and LCA)
    ALF.mu.add_model(1)
    
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
                   particle_size='Dissolved gas',
                   degradability='Undegradable',
                   organic=False)
    
    N2 = Component(ID='N2',
                   particle_size='Dissolved gas',
                   degradability='Undegradable',
                   organic=False)
    
    H2 = Component(ID='H2',
                   particle_size='Dissolved gas',
                   degradability='Undegradable',
                   organic=False)
    
    CH4 = Component(ID='CH4',
                    particle_size='Dissolved gas',
                    degradability='Slowly',
                    organic=True)
    
    C2H4 = Component(ID='C2H4',
                     particle_size='Dissolved gas',
                     degradability='Slowly',
                     organic=True)
    
    CO2 = Component(ID='CO2',
                    particle_size='Dissolved gas',
                    degradability='Undegradable',
                    organic=False)
    
    CO = Component(ID='CO',
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
    
    cmps = Components([Al2O3, SiO2, Fe2O3, Fe, AlH3O3, ALF,
                       HCOOH, H2O, O2, N2, H2, CH4, CO2, CO,
                       C2H6O, C2H4, CH4O, C3H8O])
    
    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_alias('H2O', 'Water')
    if set_thermo: qs_set_thermo(cmps)

    return cmps