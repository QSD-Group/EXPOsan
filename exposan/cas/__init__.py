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

import os, qsdsan as qs
from exposan.utils import _init_modules
cas_path = os.path.dirname(__file__)
module = os.path.split(cas_path)[-1]
data_path, results_path = _init_modules(module, include_data_path=True)
del os


# %%

# =============================================================================
# Load components and system
# =============================================================================

from . import _components
from ._components import create_components

from . import system
from .system import create_system
_loaded = False
def load():
    # `create_system` creates the components and sets the thermo itself,
    # so `components` is sourced from the built system to stay consistent.
    global sys, components, _loaded
    sys = create_system()
    sys.simulate()
    components = qs.get_thermo().chemicals
    _loaded = True
    dct = globals()
    dct.update(sys.flowsheet.to_dict())


def __getattr__(name):
    if not _loaded:
        raise AttributeError(
            f'Module {__name__} does not have the attribute "{name}" '
            'and the module has not been loaded, '
            f'loading the module with `{__name__}.load()` may solve the issue.')


__all__ = (
    'cas_path',
    *_components.__all__,
    *system.__all__,
 	)