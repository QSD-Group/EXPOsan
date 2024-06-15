#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Jianan Feng <jiananf2@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs
from exposan.co2_sorbent import create_system_A, create_system_B, create_system_C
from chaospy import distributions as shape
from qsdsan.utils import DictAttrSetter

# TODO: add all parameters first with distributions for some parameters which are already determined
# TODO: later we can add distributions for other parameters

__all__ = (
    'create_model_A', # ALF production using Al(OH)3
    'create_model_B', # ALF production using Bauxite
    'create_model_C' # CO2 capture and utilization
    )

# =============================================================================
# ALF production: Al(OH)3 + HCOOH
# =============================================================================
def create_model_A(): pass
    
# =============================================================================
# ALF production: Bauxite + HCOOH
# =============================================================================
def create_model_B(): pass

# =============================================================================
# CO2 capture and utilization
# =============================================================================
def create_model_C(): pass