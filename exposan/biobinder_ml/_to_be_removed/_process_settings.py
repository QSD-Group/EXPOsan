#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import biosteam as bst, qsdsan as qs
# from biosteam.units.design_tools import CEPCI_by_year
from exposan import htl

__all__ = (
    '_load_process_settings',
    'CEPCI_by_year',
    )

#!!! Update the numbers in QSDsan
CEPCI_by_year = {
    'Seider': 567,
    1990: 357.6,
    1991: 361.3,
    1992: 358.2,
    1993: 359.2,
    1994: 368.1,
    1995: 381.1,
    1996: 381.7,
    1997: 386.5,
    1998: 389.5,
    1999: 390.6,
    2000: 394.1,
    2001: 394.3,
    2002: 395.6,
    2003: 402,
    2004: 444.2,
    2005: 468.2,
    2006: 499.6,
    2007: 525.4,
    2008: 575.4,
    2009: 521.9,
    2010: 550.8,
    2011: 585.7,
    2012: 584.6,
    2013: 567.3,
    2014: 576.1,
    2015: 556.8,
    2016: 541.7,
    2017: 567.5,
    2018: 603.1,
    2019: 607.5,
    2020: 596.2,
    2021: 708.8,
    2022: 816,
    2023: 798,
    }


#!!! Need to update process settings such as utility price
def _load_process_settings():
    htl._load_process_settings()
    bst.CE = 2023
    # bst.PowerUtility().price = 