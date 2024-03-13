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
    # create_tea,
    )
from exposan.co2_sorbent import _sanunits as su
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
    # TODO: update CEPCI year
    bst.CE = qs.CEPCI_by_year[2020]
    
    # TODO: change the flow rate of AlOH3
    # TODO: add price for Al(OH)3
    aluminum_hydroxide = qs.WasteStream('aluminum_hydroxide', phase='s', AlH3O3=100, units='kg/h', T=25+273.15)
    
    # TODO: change the flow rate of DI H2O
    # TODO: add price for H2O
    water = qs.WasteStream('water', H2O=100, units='kg/h', T=25+273.15)

    M1 = qsu.Mixer('feed_mixer_1', ins=(aluminum_hydroxide, water), outs='AlOH3_H2O', init_with='Stream', rigorous=True, conserve_phases=True)
    M1.register_alias('M1')
    
    # TODO: change the flow rate of formic acid
    # TODO: add price for formic acid
    HCOOH_H2O = qs.WasteStream('formic_acid', HCOOH=350, H2O=100, units='kg/h', T=25+273.15)
    
    M2 = qsu.Mixer('feed_mixer_2', ins=(M1-0, HCOOH_H2O), outs='AlOH3_H2O_HCOOH', init_with='Stream', rigorous=True, conserve_phases=True)
    M2.register_alias('M2')
    
    R1 = su.ALFProduction('ALF_production', M2-0, outs='ALF_solution')
    R1.register_alias('R1')
        
    # update operating hours if necessary
    sys = qs.System.from_units(f'sys_ALF_A', units=list(flowsheet.unit), operating_hours=7920)
    sys.register_alias('sys')

    sys.simulate()
    sys.diagram()

# =============================================================================
# ALF_crystallizer
# 
# ALF_dryer
# 
# CO2_absorber
# 
# CO2_stripper
# 
# Electrochemical_cell
# 
# HXN
# 
# Cooling_tower
# 
# CHP
# 
#     sys = qs.System.from_units(
#         f'sys_ALF_A',
#         units=list(flowsheet.unit), 
#         operating_hours=,
#         )
#     
# 
#     qs.StreamImpactItem(ID='',
#                         linked_stream=stream.,
#                         Acidification=,
#                         Ecotoxicity=,
#                         Eutrophication=,
#                         GlobalWarming=,
#                         OzoneDepletion=,
#                         PhotochemicalOxidation=,
#                         Carcinogenics=,
#                         NonCarcinogenics=,
#                         RespiratoryEffects=)
#     
#     create_tea(system=sys,
#                IRR_value=,
#                income_tax_value=,
#                finance_interest_value=,
#                labor_cost_value=)
#     
#     qs.LCA(system=sys,
#            lifetime=,
#            lifetime_unit='yr',
#            Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
#            Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
#     
#     sys.simulate()
# =============================================================================
    
    return sys

# =============================================================================
# Bauxite + HCOOH
# =============================================================================
def create_system_B():

    flowsheet_ID = f'ALF_B'