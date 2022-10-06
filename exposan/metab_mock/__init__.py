# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 13:50:10 2022

@author: joy_c
"""

import os
folder = os.path.dirname(__file__)
data_path = os.path.join(folder, 'data')
results_path = os.path.join(folder, 'results')
figures_path = os.path.join(folder, 'figures')
# To save simulation results and generated figures
if not os.path.isdir(results_path): os.mkdir(results_path)
if not os.path.isdir(figures_path): os.mkdir(figures_path)
del os

from . import units
from .units import *

from . import system
from .system import *

from . import model
from .model import *

__all__ = (
    'folder',
    'data_path',
    'results_path',
    'figures_path',
    *units.__all__,
    *system.__all__,
    *model.__all__,
	)