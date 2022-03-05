#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os
br_path = os.path.dirname(__file__)
data_path = os.path.join(br_path, 'data')
results_path = os.path.join(br_path, 'results')
# To save simulation data and generated figures
if not os.path.isdir(results_path): os.mkdir(results_path)
del os

from . import _cmps, systems, models

from ._cmps import *
from .systems import *
from .models import *

__all__ = (
	'br_path',
	'data_path',
	*_cmps.__all__,
	*systems.__all__,
    *models.__all__,
	)
