# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Saumitra Rai <raisaumitra9@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np, qsdsan as qs
from qsdsan import processes as pc, sanunits as su, WasteStream, System

__all__ = ('create_components',)

def create_components():
     asm2d_cmps = pc.create_asm2d_cmps()
     asm2d_cmps.X_S.f_BOD5_COD = 0.54
     CO2 = qs.Component('S_CO2', search_ID='CO2', particle_size='Dissolved gas',
                        degradability='Undegradable', organic=False)
     CH4 = qs.Component('S_CH4', search_ID='CH4', particle_size='Dissolved gas',
                        degradability='Undegradable', organic=False)
     H2 = qs.Component('S_H2', search_ID='H2', particle_size='Dissolved gas',
                       degradability='Undegradable', organic=False)
     cmps1 = qs.Components.load_default()
     ash = cmps1.X_Ig_ISS.copy('ash')
     cmps = qs.Components([*asm2d_cmps, CO2, CH4, H2, ash])
     cmps.compile()
     return cmps
