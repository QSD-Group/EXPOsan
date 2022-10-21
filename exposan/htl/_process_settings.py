#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Jianan Feng <jiananf2@illinois.edu>
    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import biosteam as bst, qsdsan as qs

def load_process_settings():
    #Use p-Terphenyl as the heating agent
    
    Terphenyl = qs.Component('Terphenyl',CAS='92-94-4',phase='l',
                          particle_size='Soluble',
                          degradability='Slowly',organic=True)
    heating_oil_thermo = bst.Thermo((Terphenyl,))
    
    heating_oil = bst.UtilityAgent(
                'heating_oil',
                Terphenyl=1, T=693.15, P=1516847.2, phase='l',
                thermo=heating_oil_thermo,T_limit=500,
                regeneration_price = 0.2378,
                heat_transfer_efficiency = 0.9,
            ) #values are not true
    
    '''
    T_limit is the temperature limit of outlet utility streams [K].
    If no limit is given, phase change is assumed.
    If utility agent heats up, `T_limit` is the maximum temperature.
    If utility agent cools down, `T_limit` is the minimum temperature.
    
    
    So we should give a `T_limit` as the heating oil is not changing phase,
    and we should give a `T_limit` that is smaller than `T` since the heating oil
    cools down when the sludge feed heats up.
    
    500 is made up
    '''
    bst.HeatUtility.heating_agents.append(heating_oil)