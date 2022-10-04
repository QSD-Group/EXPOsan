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
from exposan.utils import _init_modules
adm_path = os.path.dirname(__file__)
module = os.path.split(adm_path)[-1]
data_path, results_path, figures_path = \
    _init_modules(module, include_data_path=True, include_figures_path=True)
del os


from . import system
from .system import *

from . import model
from .model import *

__all__ = (
    'adm_path',
    'data_path',
    'results_path',
    'figures_path',
    *system.__all__,
    *model.__all__,
	)