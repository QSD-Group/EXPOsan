#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Jianan Feng <jiananf2@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References:

'''

import os, qsdsan as qs, biosteam as bst
from qsdsan import sanunits as qsu
from biosteam.units import IsenthalpicValve
from qsdsan.utils import clear_lca_registries
from exposan.co2_sorbent import (
    _load_components,
    create_tea,
    )
from exposan.htl import _sanunits as su
from biosteam import settings

__all__ = (
    'create_system_A',
    'create_system_B'
    )

# =============================================================================
# Al(OH)3 + HCOOH
# =============================================================================
def create_system_A():

    flowsheet_ID = f'ALF_A'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    bst.CE = qs.CEPCI_by_year[2020]

Mixer

H1

ALF_reactor

ALF_dryer

CO2_absorber

CO2_stripper

Electrochemical_cell

HXN

Cooling_tower

CHP

    sys = qs.System.from_units(
        f'sys_ALF_A',
        units=list(flowsheet.unit), 
        operating_hours=,
        )
    sys.register_alias('sys')

    qs.StreamImpactItem(ID='',
                        linked_stream=stream.,
                        Acidification=,
                        Ecotoxicity=,
                        Eutrophication=,
                        GlobalWarming=,
                        OzoneDepletion=,
                        PhotochemicalOxidation=,
                        Carcinogenics=,
                        NonCarcinogenics=,
                        RespiratoryEffects=)
    
    create_tea(system=sys,
               IRR_value=,
               income_tax_value=,
               finance_interest_value=,
               labor_cost_value=)
    
    qs.LCA(system=sys,
           lifetime=,
           lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    sys.simulate()
    
    return sys

# =============================================================================
# Bauxite + HCOOH
# =============================================================================
def create_system_B():

    flowsheet_ID = f'ALF_B'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    bst.CE = qs.CEPCI_by_year[2020]

Mixer

H1

ALF_reactor

ALF_dryer

CO2_absorber

CO2_stripper

Electrochemical_cell

HXN

Cooling_tower

CHP


    sys = qs.System.from_units(
        f'sys_ALF_B',
        units=list(flowsheet.unit), 
        operating_hours=,
        )
    sys.register_alias('sys')

    qs.StreamImpactItem(ID='',
                        linked_stream=stream.,
                        Acidification=,
                        Ecotoxicity=,
                        Eutrophication=,
                        GlobalWarming=,
                        OzoneDepletion=,
                        PhotochemicalOxidation=,
                        Carcinogenics=,
                        NonCarcinogenics=,
                        RespiratoryEffects=)
    
    create_tea(system=sys,
               IRR_value=,
               income_tax_value=,
               finance_interest_value=,
               labor_cost_value=)
    
    qs.LCA(system=sys,
           lifetime=,
           lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    sys.simulate()
    
    return sys