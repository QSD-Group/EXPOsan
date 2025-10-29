# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os
from exposan.utils import _init_modules
edu_path = os.path.dirname(__file__)
module = os.path.split(edu_path)[-1]
data_path, results_path, figures_path = \
    _init_modules(module, include_data_path=True, include_figures_path=True)
del os

from . import system
from .system import *

_system_loaded = False
def load(reload=False, inf_kwargs={}, asm_kwargs={}, init_conds={}):
    global _system_loaded
    if not _system_loaded: reload = True
    if reload:
        global cmps, components, asm, sys
        sys = create_system(
                inf_kwargs=inf_kwargs,
                asm_kwargs=asm_kwargs,
                init_conds=init_conds,
                )    
        O1 = sys.flowsheet.unit.O1
        cmps = components = O1.components
        asm = O1.suspended_growth_model
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    _system_loaded = True

def __getattr__(name):
    if not _system_loaded:
        raise AttributeError(f'module "{__name__}" not yet loaded, '
                             f'load module with `{__name__}.load()`.')

__all__ = (
    'edu_path',
    'data_path',
    'results_path',
    'figures_path',
    *system.__all__,
	)