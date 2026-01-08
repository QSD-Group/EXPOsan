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
	)

from . import (
    B1, B2, B3, 
    C1, C2, C3, 
    E2, E2P, 
    F1, 
    G1, G2, G3, 
    H1, 
    I1, I2, I3, 
    N1, N2
    )

from .B1 import *
from .B2 import *
from .B3 import *
from .C1 import *
from .C2 import *
from .C3 import *
from .E2 import *
from .E2P import *
from .F1 import *
from .G1 import *
from .G2 import *
from .G3 import *
from .H1 import *
from .I1 import *
from .I2 import *
from .I3 import *
from .N1 import *
from .N2 import *

def create_system(ID):
    f = globals()[f'create_{ID.lower()}_system']
    return f()