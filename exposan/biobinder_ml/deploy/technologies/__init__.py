# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Mature waste-treatment alternatives with explicit screening assumptions."""
from ._common import Feedstock, Economics, Equipment, Result, FEEDSTOCKS
from .ad import ADConfig, create_system as create_ad_system, adbc_dairy_costs
from .composting import CompostConfig, create_system as create_composting_system
from ._pactea import PaCTEAConfig

__all__ = ('Feedstock','Economics','Equipment','Result','FEEDSTOCKS','ADConfig',
           'CompostConfig','PaCTEAConfig','create_ad_system','create_composting_system','adbc_dairy_costs')
