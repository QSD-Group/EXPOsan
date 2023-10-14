#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Bright Elijah <be05055@georgiasouthern.edu & brightcarlelijah@gmail.com>
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

from qsdsan import Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_components', )

def create_components(set_thermo=True):
    Mg = Component('Mg', phase='l', particle_size='Soluble',
                    degradability='Undegradable', organic=False)
    
    Ca = Component('Ca', phase='l', particle_size='Soluble',
                   degradability='Undegradable', organic=False)
    
    H2O = Component('H2O', phase='l', particle_size='Soluble',
                    degradability='Undegradable', organic=False)
    
    NaClO = Component('NaClO', phase='l', particle_size = 'Soluble', 
                      degradability = 'Undegradable', organic=False)
    
    Ecoli = Component('Ecoli', MW = 1, phase = 's', particle_size = 'Particluate', 
                      degradability = 'Undegradable', organic = False)
    
    add_V_from_rho(Ecoli, 500)
    # Ecoli.copy_models_from(tmo.Chemical('Glucose'), ('Cn', 'mu'))
    
    Polyethylene = Component('Polyethylene', formula='C2H4',
                             phase='s', particle_size='Particulate',
                             degradability='Slowly', organic=True)
    
    # PVC = Component('PVC', formula='H2C–CHCl',
    #                 phase='s', particle_size='Particulate',
    #                 degradability='Undegradable', organic=True)
    
    # Quartz = Component('Quartz', formula='SiO2',
    #                   phase='s', particle_size='Particulate', 
    #                   degradability='Undegradable', organic=False)
    
    AgNp = Component('AgNp', search_ID='Ag',
                     phase='l', particle_size='Particulate',
                     degradability='Undegradable', organic=False)
    
    components = Components((Mg, Ca, H2O, NaClO, Ecoli, Polyethylene, AgNp))
    
    components.default_compile(soluble_ref='H2O') # default properties to those of water
    components.set_alias('H2O', 'Water')

    if set_thermo: qs_set_thermo(components)

    return components