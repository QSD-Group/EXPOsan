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
from exposan.utils import _init_modules
bsm1_path = os.path.dirname(__file__)
module = os.path.split(bsm1_path)[-1]
data_path, results_path, figures_path = \
    _init_modules(module, include_data_path=True, include_figures_path=True)
del os

from . import system
from .system import *

_system_loaded = False
def load(reload=False, suspended_growth_model='ASM1', reactor_model='CSTR', 
         inf_kwargs={}, asm_kwargs={}, settler_kwargs={}, 
         init_conds={}, aeration_processes=()):
    global _system_loaded
    if not _system_loaded: reload = True
    if reload:
        global cmps, components, asm, sys
        sys = create_system(
            suspended_growth_model=suspended_growth_model,
            reactor_model=reactor_model,
            inf_kwargs=inf_kwargs,
            asm_kwargs=asm_kwargs,
            settler_kwargs=settler_kwargs,
            init_conds=init_conds,
            aeration_processes=aeration_processes,
            )
        # O1 = sys.flowsheet.unit.O1
        # cmps = components = O1.components
        # asm = O1.suspended_growth_model
        # Legacy names
        global PE, SE#, RE
        stream = sys.flowsheet.stream
        PE = stream.wastewater
        SE = stream.effluent
        # RE = stream.RWW
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    _system_loaded = True


def __getattr__(name):
    if not _system_loaded:
        raise AttributeError(
            f'Module {__name__} does not have the attribute "{name}" '
            'and the module has not been loaded, '
            f'loading the module with `{__name__}.load()` may solve the issue.')


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