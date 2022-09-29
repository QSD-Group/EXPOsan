# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Jianan Feng <jiananf2@illinois.edu>
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('creat_components',)

def create_components(set_thermo=True):
    
    # Solids components
    
    sludge = Component('sludge', phase='s', formula='C56H95O27N6P',
                       particle_size='Particulate',
                       degradability='Undegradable',organic=False)
    add_V_from_rho(sludge,721)   
    #https://www.aqua-calc.com/page/density-table/substance/sewage-coma-and-blank-sludge
    sludge.HHV = 22.0  #Li et al 2018
    sludge.copy_models_from(Chemical('glucose'),('Cn',)) # glucose will be replaced
    
    struvite = Component('struvite',search_ID='MagnesiumAmmoniumPhosphate',
                         formula='NH4MgPO4·H12O6',phase='s',
                         particle_size='Particulate', 
                         degradability='Undegradable',
                         organic=False)
    add_V_from_rho(struvite, 1710)
    #http://webmineral.com/data/Struvite.shtml#.YzYvqOzMIiM
    
    # Oil components
    
    # Gas components
    
    CH4 = Component('CH4',phase='g',particle_size='Dissolved gas',
                    degradability='Slowly',organic=True)
    
    C2H6 = Component('C2H6',phase='g',particle_size='Dissolved gas',
                   degradability='Slowly',organic=True)
    
    C3H8 = Component('C3H8',phase='g',particle_size='Dissolved gas',
                   degradability='Slowly',organic=True)
    
    C5H12 = Component('C5H12',phase='g',particle_size='Dissolved gas',
                   degradability='Slowly',organic=True)
    
    CO2 = Component('CO2',phase='g',particle_size='Dissolved gas',
                    degradability='Undegradable',organic=False)
    
    CO = Component('CO',phase='g',particle_size='Dissolved gas',
                   degradability='Undegradable',organic=False)
    
    H2 = Component('H2',phase='g',particle_size='Dissolved gas',
                   degradability='Undegradable',organic=False)
    
    # Other components
    
    H2SO4 = Component('H2SO4',phase='l',particle_size='Soluble',
                      degradability='Undegradable',organic=False)
    
    MgCl2 = Component('MgCl2',phase='l',particle_size='Soluble',
                      degradability='Undegradable',organic=False)
    
    NaOH = Component('NaOH',phase='l',particle_size='Soluble',
                     degradability='Undegradable',organic=False)
    
    NH42SO4 = Component('NH42SO4',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False)
    add_V_from_rho(NH42SO4, 1770)
    #https://en.wikipedia.org/wiki/Ammonium_sulfate
    
    H2O = Component('H2O',phase='l',particle_size='Soluble',
                    degradability='Undegradable',organic=False)
         
    NH3 = Component('NH3',phase='g',particle_size='Dissolved gas',
                    degradability='Undegradable',organic=False)


    cmps = Components([sludge,struvite,
                       CH4,C2H6,C3H8,C5H12,CO2,CO,H2,
                       H2SO4,MgCl2,NaOH,NH42SO4,H2O,NH3])
    
    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()

    cmps.set_alias('H2O', 'Water')
