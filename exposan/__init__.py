#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License. Please refer to
https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import importlib.metadata as impmeta
try:
    __version__ = impmeta.version('exposan')
except impmeta.PackageNotFoundError:
    __version__ = None
    
import os
es_path = os.path.dirname(__file__)
del os, impmeta

from . import utils