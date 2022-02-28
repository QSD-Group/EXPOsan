#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os
bsm1_path = os.path.dirname(__file__)
data_path = os.path.join(bsm1_path, 'data')
results_path = os.path.join(bsm1_path, 'results')
figures_path = os.path.join(bsm1_path, 'figures')
del os

from . import system
from .system import *

from . import model
from .model import *

__all__ = (
    'bsm1_path',
    'data_path',
    'results_path',
    'figures_path',
    *system.__all__,
    *model.__all__,
	)