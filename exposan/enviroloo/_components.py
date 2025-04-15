#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''

This module is developed by:
    Siqi Tang <siqit@outlook.com>
    Yuyao Huang <yuyaoh20@gmail.com>
               for Enviroloo Clear Toilet system
               
'''
import qsdsan as qs
import os
from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho
from exposan.bwaise import create_components as create_bw_components
import numpy as np
from qsdsan import processes as pc, sanunits as su, Model as mod
from chaospy import distributions as shape

__all__ = ('create_components', )

## Pre-define the involving components in the toilet system of concern such as Enviroloo Clear
def create_components(set_thermo = True
                      #adjust_MW_to_measured_as=False
                      ):
    # bw_cmps = create_bw_components(set_thermo=False)
    masm2d_cmps = pc.create_masm2d_cmps(set_thermo=False)
    cmps = Components((*masm2d_cmps, 
                       # Tissue, WoodAsh, H2O,
                       # C, SolubleCH4, 
                       # #H2O, CO2, 
                       # CH4, 
                       # N2O, NH3
                       # Glucose, 
                       # O3, air, PAC, NaOH, NaClO
                       ))
    
    # for i in cmps:
    #     for attr in ('HHV', 'LHV', 'Hf'):
    #         if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()    
    # cmps.compile(ignore_inaccurate_molar_weight=True) #override for runtime error where N2_S molecular weight was not found

    # cmps.set_alias('H2O', 'Water')
    # #cmps.set_alias('CO2', 'Carbon Dioxide')
    # cmps.set_alias('CH4', 'Methane')
    if set_thermo: qs_set_thermo(cmps)

    return cmps

# cmps = create_components()
# print("List of Components:")
# for cmp in cmps:
#     print(cmp.ID)
#     print(f"Description: {cmp.description}")
