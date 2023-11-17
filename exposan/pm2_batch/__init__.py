# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os
from exposan.utils import _init_modules
pm2_path = os.path.dirname(__file__)
module = os.path.split(pm2_path)[-1]
data_path, results_path, figures_path = \
    _init_modules(module, include_data_path=True, include_figures_path=True)
del os

from . import system

from .system import *

_system_loaded = False
def load(reload=False, pm2_kwargs={}, init_conds={}):

    global _system_loaded
    if not _system_loaded: reload = True
    if reload:
        global cmps, components, pm2, sys
        sys = create_system(
            pm2_kwargs=pm2_kwargs,
            init_conds=init_conds,
            )
        PBR = sys.flowsheet.unit.PBR
        cmps = components = PBR.components
        pm2 = PBR.model
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    _system_loaded = True

def __getattr__(name):
    if not _system_loaded:
        raise AttributeError(f'module "{__name__}" not yet loaded, '
                             f'load module with `{__name__}.load()`.')

from . import model
from .model import *

__all__ = (
    'pm2_path',
    'data_path',
    'results_path',
    'figures_path',
    *system.__all__,
    *model.__all__,
	)