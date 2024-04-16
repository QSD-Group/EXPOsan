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

_system_loaded = False
def load(reload=False, inf_kwargs={}, adm_kwargs={}, init_conds={}):
    global _system_loaded
    if not _system_loaded: reload = True
    if reload:
        global cmps, components, adm, sys
        sys = create_system(
            inf_kwargs=inf_kwargs,
            adm_kwargs=adm_kwargs,
            init_conds=init_conds,
            )
        AD = sys.flowsheet.unit.AD
        cmps = components = AD.components
        adm = AD.model
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
    'adm_path',
    'data_path',
    'results_path',
    'figures_path',
    *system.__all__,
    *model.__all__,
	)