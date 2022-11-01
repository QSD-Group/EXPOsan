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

from . import systems
from .systems import *

from . import models
from .models import *

__all__ = (
    'folder',
    'data_path',
    'results_path',
    'figures_path',
    *units.__all__,
    *systems.__all__,
    *models.__all__,
	)