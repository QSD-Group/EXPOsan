#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

'''
#TODO
Other things that I think we might want to include:
    - A figure in the SI (or just on GitHub) showing that we can converge to
    similar steady-states conditions
'''

from qsdsan.utils import ords
from exposan.bsm1 import model_bsm1

# Why not use something meaningful as the seed?
# guess what this is (hint: check out the built-in function `ord`)
seed = ords('bsm1')
N = 100 # number of simulations
samples = model_bsm1.sample(N=N, rule='L', seed=seed)
model_bsm1.load_samples(samples)