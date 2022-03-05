#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os
bwaise_path = os.path.dirname(__file__)
data_path = os.path.join(bwaise_path, 'data')
results_path = os.path.join(bwaise_path, 'results')
figures_path = os.path.join(bwaise_path, 'figures')
# To save simulation data and generated figures
if not os.path.isdir(results_path): os.mkdir(results_path)
if not os.path.isdir(figures_path): os.mkdir(figures_path)
del os

from qsdsan.utils import time_printer

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


from . import _cmps, _lca_data, systems, models

from ._cmps import *
from ._lca_data import *
from .systems import *
from .models import *

__all__ = (
    'evaluate',
    'get_key_metrics',
    'bwaise_path',
    'data_path',
    'results_path',
    'figures_path',
	*_cmps.__all__,
    *_lca_data.__all__,
	*systems.__all__,
    *models.__all__,
	)