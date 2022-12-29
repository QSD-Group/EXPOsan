#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Jianan Feng <jiananf2@illinois.edu>
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, qsdsan as qs
from qsdsan.utils import auom
from exposan.utils import _init_modules

htl_path = os.path.dirname(__file__)
module = os.path.split(htl_path)[-1]
data_path, results_path = _init_modules(module, include_data_path=True)

_kg_to_g = auom('kg').conversion_factor('g')
_m3perh_to_MGD = auom('m3/h').conversion_factor('MGD')
_MJ_to_MMBTU = auom('MJ').conversion_factor('MMBTU')
_MMgal_to_L = auom('gal').conversion_factor('L')*1000000


# %%

# =============================================================================
# Load components and systems
# =============================================================================

from . import _components
from ._components import *
_components_loaded = False
def _load_components(reload=False):
    global components, _components_loaded
    if not _components_loaded or reload:
        components = create_components()
        qs.set_thermo(components)
        _components_loaded = True

# from . import _process_settings
# from ._process_settings import *

# from . import _tea
# from ._tea import *

# from . import systems
# from .systems import *
# _system_loaded = False
# def _load_system(configuration='baseline'):
#     global sys, tea, lca, _system_loaded
#     sys = create_system(configuration)
#     tea = sys.TEA
#     lca = sys.LCA
#     _system_loaded = True


from . import _process_settings
from ._process_settings import *

from . import _tea
from ._tea import *

from . import systems
from .systems import *

# def _load_system(configuration='baseline'):
#     global sys, tea, lca, _system_loaded
#     sys = create_system(configuration)
#     tea = sys.TEA
#     lca = sys.LCA
#     _system_loaded = True

_system_loaded = False
def load(configuration='baseline'):
    global sys, tea, lca, flowsheet, _system_loaded
    sys = create_system(configuration)
    tea = sys.TEA
    lca = sys.LCA
    flowsheet = sys.flowsheet
    _system_loaded = True
    dct = globals()
    dct.update(sys.flowsheet.to_dict())


def __getattr__(name):
    if not _components_loaded or not _system_loaded:
        raise AttributeError(
            f'Module {__name__} does not have the attribute "{name}" '
            'and the module has not been loaded, '
            f'loading the module with `{__name__}.load()` may solve the issue.')
        
from . import models
from .models import *

__all__ = (
    'htl_path',
    'data_path',
    'results_path',
    *_components.__all__,
    *_process_settings.__all__,
    *_tea.__all__,
    *systems.__all__,
    *models.__all__,
)