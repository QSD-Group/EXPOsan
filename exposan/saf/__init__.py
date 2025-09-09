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

saf_path = os.path.dirname(__file__)
module = os.path.split(saf_path)[-1]
data_path, results_path = _init_modules(module, include_data_path=True)


# %%

# =============================================================================
# Load components and systems
# =============================================================================

from . import utils
from .utils import *

# Default settings for consistency across the module
from . import _process_settings
from ._process_settings import *

from . import _components
from ._components import *
_components_loaded = False
def _load_components(reload=False):
    global components, _components_loaded
    if not _components_loaded or reload:
        components = create_components()
        qs.set_thermo(components)
        _components_loaded = True

from . import _units
from ._units import *

from . import systems
from .systems import *

_system_loaded = False
def load(configuration='baseline',):
    global sys, tea, lca, flowsheet, _system_loaded
    configuration = configuration.lower()   
    if configuration == 'baseline':
        kwargs = config_baseline
    elif configuration == 'ec':
        kwargs = config_EC
    elif configuration in ('ec-future', 'ec_future', 'ec future'):
        kwargs = config_EC_future
    elif configuration == 'no_psa':
        kwargs = config_no_PSA
    else:
        raise ValueError('Configuration only be "baseline", "EC", "EC future", or "no PSA", '
                         f'not {kwargs}.')

    sys = create_system(**kwargs)
    tea = sys.TEA
    lca = sys.LCA
    flowsheet = sys.flowsheet
    _system_loaded = True
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    simulate_and_print(sys)

def __getattr__(name):
    if not _components_loaded or not _system_loaded:
        raise AttributeError(
            f'Module {__name__} does not have the attribute "{name}" '
            'and the module has not been loaded, '
            f'loading the module with `{__name__}.load()` may solve the issue.')

__all__ = (
    'saf_path',
    'data_path',
    'results_path',
    *_components.__all__,
    *_process_settings.__all__,
    *_units.__all__,
    *systems.__all__,
    *utils.__all__,
)