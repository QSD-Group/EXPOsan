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
from exposan.bwaise import create_components as create_HTL_components

__all__ = ('creat_components',)

def create_components(set_thermo=True):
    HTL_cmps=create_HTL_components(set_thermo=False)
    
    H2SO4 = Component('H2SO4',phase='l',particle_size='Soluble',degradability='Undegradable',organic=False)
    MgCl2 = Component('MgCl2',phase='l',particle_size='Soluble',degradability='Undegradable',organic=False)
    NaOH = Component('NaOH',phase='l',particle_size='Soluble',degradability='Undegradable',organic=False)
    NH42SO4 = Component('NH42SO4',phase='l',particle_size='Soluble',degradability='Undegradable',organic=False)
    add_V_from_rho(NH42SO4, 1770)
    H2O = Component('H2O',phase='l',particle_size='Soluble',degradability='Undegradable',organic=False)
    CH4 = Component('CH4',phase='g',particle_size='Dissolved gas',degradability='Slowly',organic=True)
    CO2 = Component('CO2',phase='g',particle_size='Dissolved gas',degradability='Undegradable',organic=False)
    CO = Component('CO',phase='g',particle_size='Dissolved gas',degradability='Undegradable',organic=False)
    H2 = Component('H2',phase='g',particle_size='Dissolved gas',degradability='Undegradable',organic=False)
    NH3 = Component('NH3',phase='g',particle_size='Dissolved gas',degradability='Undegradable',organic=False)
    struvite = Component('struvite',search_ID='MagnesiumAmmoniumPhosphate',formula='NH4MgPO4·H12O6',phase='s',
                         particle_size='Particulate', degradability='Undegradable',organic=False)
    HT_Flue_gas_hydrocarbon = Component('HTL_flue_gas',phase='g',formula='C3H13O8',particle_size='Dissolved gas',degradability='Undegradable',organic=True)
    HT_Flue_gas_hydrocarbon.copy_models_from(Chemical('Propane'),('Cn',))
    CHG_Flue_gas_hydrocarbon = Component('CHG_flue_gas',phase='g',formula='C4H9O5',particle_size='Dissolved gas',degradability='Undegradable',organic=True)
    CHG_Flue_gas_hydrocarbon.copy_models_from(Chemical('butane'),('Cn',))

