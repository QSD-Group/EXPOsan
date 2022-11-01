# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

'''
TODO:
    - Add cost
    - Some of the design may refer to Table 10-3 on page 1068 of Metcalf and Eddy
        - Process (a) low loaded anaerobic lagoon system or (n) plug flow anaerobic system
'''

import qsdsan as qs
from qsdsan import Component, Components

__all__ = ('create_components',)

def create_components(set_thermo=True):

    X = Component('X', phase='s', measured_as='COD', i_COD=0, description='Biomass',
                  organic=True, particle_size='Particulate', degradability='Readily')
    X_inert = Component('X_inert', phase='s', description='Inert biomass', i_COD=0,
                        organic=True, particle_size='Particulate', degradability='Undegradable')
    Substrate = Component('Substrate', phase='l', i_COD=1,
                          organic=True, particle_size='Soluble', degradability='Readily')
    CH4 = Component('CH4', phase='g', organic=True, particle_size='Dissolved gas',
                    degradability='Readily')
    cmps = Components([X, X_inert, Substrate, CH4])
    
    # Default-adding components needed for combustion
    cmps = Components.append_combustion_components(cmps, lock_state_at='')
    
    # Define component groups
    cmps.define_group('active_biomass', IDs=('X',))
    cmps.define_group('inert_biomass', IDs=('X_inert',))
    cmps.define_group('substrates', IDs=('Substrate',))
    
    if set_thermo: qs.set_thermo(cmps)
    
    return cmps