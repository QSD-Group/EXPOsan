# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, qsdsan as qs
from exposan.utils import _init_modules

bsm2_path = os.path.dirname(__file__)
module = os.path.split(bsm2_path)[-1]
data_path, results_path, figures_path = \
    _init_modules(module, include_data_path=True, include_figures_path=True)


# %%

# =============================================================================
# Load the system
# =============================================================================

from . import system
from .system import *
_system_loaded = False
def _load_system():
    global sys, _system_loaded
    sys = create_system()
    _system_loaded = True


def load():
    if not _system_loaded: _load_system()
    dct = globals()
    dct.update(sys.flowsheet.to_dict())


def __getattr__(name):
    if not _system_loaded:
        raise AttributeError(
            f'Module {__name__} does not have the attribute "{name}" '
            'and the module has not been loaded, '
            f'loading the module with `{__name__}.load()` may solve the issue.')


__all__ = (
    'bsm2_path',
    'figures_path',
    'results_path',
    *system.__all__,
)