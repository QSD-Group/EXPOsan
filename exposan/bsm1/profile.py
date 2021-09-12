#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 20:19:58 2021

@author: Yalin Li
"""

# Add path for the local thermosteam and biosteam,
# shouldn't need to do this if using Spyder
# and having bst/tmo path added to sys.path
import os, sys
dir_path = os.path.dirname('systems.py')
tmo_path = os.path.join(dir_path, '../../../tmo')
bst_path = os.path.join(dir_path, '../../../bst')
sys.path.extend([tmo_path, bst_path])

# Run cProfile for systems.py
import cProfile
from systems import bsm1
cProfile.run('bsm1.simulate(t_span = (0, 0.1))', 'systems_01.prof')