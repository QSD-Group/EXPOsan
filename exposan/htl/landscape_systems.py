#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs, biosteam as bst
from qsdsan import sanunits as qsu
from qsdsan.utils import auom, clear_lca_registries
from exposan.htl import _load_components, _sanunits as su, create_tea

# TODO: add composting
__all__ = ('create_C1_system',
           'create_C1_plus_system',
           'create_C2_system',
           'create_C2_system',
           'create_C3_system',
           'create_C3_plus_system',
           'create_C4_system',
           'create_C4_plus_system',
           'create_C5_system',
           'create_C5_plus_system',
           'create_C6_system',
           'create_C6_plus_system',
           'create_C7_system',
           'create_C7_plus_system',
           'create_H1_system',
           'create_H1_plus_system',
           'create_H2_system',
           'create_H2_plus_system',
           'create_H3_system',
           'create_T1_system',
           'create_T2_system',)

#%% system C1 (landfilling)

def create_C1_system():
    pass

#%% system C1+ (landfilling + leachate treatment for PFAS removal)

def create_C1_plus_system():
    pass

#%% system C2 (heat drying + landfilling)

def create_C2_system():
    pass

#%% system C2+ (heat drying + landfilling + leachate treatment for PFAS removal)

def create_C2_plus_system():
    pass

#%% system C3 (lime stabilization + heat drying + landfilling)

def create_C3_system():
    pass

#%% system C3+ (lime stabilization + heat drying + landfilling + leachate treatment for PFAS removal)

def create_C3_plus_system():
    pass

#%% system C4 (land application)

def create_C4_system():
    pass

#%% system C4+ (land application + biosolids remediation)

def create_C4_plus_system():
    pass

#%% system C5 (heat drying + land application)

def create_C5_system():
    pass

#%% system C5+ (heat drying + land application + biosolids remediation)

def create_C5_plus_system():
    pass

#%% system C6 (lime stabilization + heat drying + land application)

def create_C6_system():
    pass

#%% system C6+ (lime stabilization + heat drying + land application + biosolids remediation)

def create_C6_plus_system():
    pass

#%% system C7 (heat drying + incineration)

def create_C7_system():
    pass

#%% system C7+ (heat drying + high temperature incineration)

def create_C7_plus_system():
    pass

#%% system H1 (hydrothermal liquefaction, biocrude as the main product)

def create_H1_system():
    pass

#%% system H1+ (hydrothermal alkaline treatment, biocrude as the main product)

def create_H1_plus_system():
    pass

#%% system H2 (hydrothermal liquefaction, renewable diesel as the main product)

def create_H2_system():
    pass

#%% system H2+ (hydrothermal alkaline treatment, renewable diesel as the main product)

def create_H2_plus_system():
    pass

#%% system H3 (supercritical water oxidation)

def create_H3_system():
    pass

#%% system T1 (pyrolysis)

def create_T1_system():
    pass

#%% system T2 (gasification)

def create_T2_system():
    pass