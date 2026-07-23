# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = (
    'create_system',
    'SYSTEM_CREATORS',
	)

from . import (
    B1, B2, B3, 
    B1E,
    C1, C2, C3, 
    C1E,
    E2, E2P, 
    F1, F1E,
    G1, G2, G3, 
    G1E,
    H1, H1E,
    I1, I2, I3, 
    I1E,
    N1, N2,
    N1E
    )

from .B1 import *
from .B1E import *
from .B2 import *
from .B3 import *
from .C1 import *
from .C1E import *
from .C2 import *
from .C3 import *
from .E2 import *
from .E2P import *
from .F1 import *
from .F1E import *
from .G1 import *
from .G1E import *
from .G2 import *
from .G3 import *
from .H1 import *
from .H1E import *
from .I1 import *
from .I1E import *
from .I2 import *
from .I3 import *
from .N1 import *
from .N1E import *
from .N2 import *

SYSTEM_CREATORS = {
    'B1': create_b1_system,
    'B1E': create_b1e_system,
    'B2': create_b2_system,
    'B3': create_b3_system,
    'C1': create_c1_system,
    'C1E': create_c1e_system,
    'C2': create_c2_system,
    'C3': create_c3_system,
    'E2': create_e2_system,
    'E2P': create_e2p_system,
    'F1': create_f1_system,
    'F1E': create_f1e_system,
    'G1': create_g1_system,
    'G1E': create_g1e_system,
    'G2': create_g2_system,
    'G3': create_g3_system,
    'H1': create_h1_system,
    'H1E': create_h1e_system,
    'I1': create_i1_system,
    'I1E': create_i1e_system,
    'I2': create_i2_system,
    'I3': create_i3_system,
    'N1': create_n1_system,
    'N1E': create_n1e_system,
    'N2': create_n2_system,
    }

def create_system(ID):
    return SYSTEM_CREATORS[ID.upper()]()
