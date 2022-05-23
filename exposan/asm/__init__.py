#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os
asm_path = os.path.dirname(__file__)
data_path = os.path.join(asm_path, 'data')
results_path = os.path.join(asm_path, 'results')
# To save simulation data
if not os.path.isdir(results_path): os.mkdir(results_path)
del os

from . import systems
from .systems import *