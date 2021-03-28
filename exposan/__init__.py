#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License. Please refer to
https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import pkg_resources
try:
    __version__ = pkg_resources.get_distribution('exposan').version
except pkg_resources.DistributionNotFound:
    __version__ = None