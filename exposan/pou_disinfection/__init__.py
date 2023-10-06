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

import os
from qsdsan.utils import time_printer
from exposan.utils import _init_modules
pou_path = os.path.dirname(__file__)
module = os.path.split(pou_path)[-1]
data_path, results_path = _init_modules(module, include_data_path=True)
del os



@time_printer
def evaluate(model, samples=None):
    if samples is not None:
        model.load_samples(samples)
    model.evaluate()

def get_key_metrics(model, alt_names={}):
    key_metrics = [i for i in model.metrics if 'total' in i.name.lower()]
    key_metrics += [i for i in model.metrics if 'net' in i.name.lower()]
    for old, new in alt_names.items():
        for i in key_metrics:
            i.name = i.name.replace(old, new)
    return key_metrics


from . import systems, models

from .systems import *
from .models import *

__all__ = (
	*systems.__all__,
    *models.__all__,
    'evaluate',
    'get_key_metrics',
	)

# %%

# # =============================================================================
# # Load components and system
# # =============================================================================

# from . import _components
# from ._components import create_components
# _components_loaded = False
# def _load_components(reload=False):
#     global components, _components_loaded
#     if not _components_loaded or reload:
#         components = create_components()
#         qs.set_thermo(components)
#         _components_loaded = True


# from . import system
# from .system import create_system
# _system_loaded = False
# def _load_system():
#     global sys, _system_loaded
#     sys = create_system()
#     _system_loaded = True


# def load():
#     if not _components_loaded: _load_components()
#     if not _system_loaded: _load_system()
#     sys.simulate()
#     dct = globals()
#     dct.update(sys.flowsheet.to_dict())


# def __getattr__(name):
#     if not _components_loaded or not _system_loaded:
#         raise AttributeError(
#             f'Module {__name__} does not have the attribute "{name}" '
#             'and the module has not been loaded, '
#             f'loading the module with `{__name__}.load()` may solve the issue.')


# __all__ = (
#     'cas_path',
#     *_components.__all__,
#     *system.__all__,
#  	)