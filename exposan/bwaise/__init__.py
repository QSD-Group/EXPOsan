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
bwaise_path = os.path.dirname(__file__)
data_path = os.path.join(bwaise_path, 'data')
results_path = os.path.join(bwaise_path, 'results')
figures_path = os.path.join(bwaise_path, 'figures')
del os


from . import _cmps, _lca_data, systems, models

from ._cmps import *
from ._lca_data import *
from .systems import *
from .models import *

__all__ = (
	*_cmps.__all__,
    *_lca_data.__all__,
	*systems.__all__,
    *models.__all__,
	)