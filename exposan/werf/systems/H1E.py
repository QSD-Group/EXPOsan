# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from .H1 import create_h1_system
from exposan.werf.utils._chp import create_chp_system

__all__ = ('create_h1e_system',)

def create_h1e_system(flowsheet=None, default_init_conds=True):
    return create_chp_system('H1E', create_h1_system, flowsheet, default_init_conds)
