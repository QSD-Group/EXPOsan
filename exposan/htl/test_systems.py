#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 08:49:12 2025

@author: jiananfeng
"""

import qsdsan as qs, biosteam as bst
from qsdsan.utils import clear_lca_registries
from exposan.htl import _load_components, create_tea, _sanunits as su, test_sanunits as tsu

__all__ = ('create_test_system',)

def create_test_system(size=10,
                             sludge_distance=100,
                             ww_2_dry_sludge_ratio=1,
                             # use balancing-area-level in the analysis
                             # kg CO2 eq/kWh
                             elec_GHG=0.40,
                             operation_hours=8760
                             ):
    
    
    flowsheet_ID = 'test_system'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # =========================================================================
    # pretreatment
    # =========================================================================
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)
    

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.8,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    Landfilling = tsu.Landfilling(ID='Landfilling',
                                  ins=WWTP-0,
                                  outs=('landfilled_solids','landfilling_emission'),
                                  MSW_wet_density=1040,
                                  hauling_distance=sludge_distance,
                                  methane_density=0.68,
                                  emission_time=1,
                                  other_pollutants=[['C2H6',1,28],],
                                  US_state='Illinois',
                                  tipping_fee=None,
                                  operation_hours=operation_hours)
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)

    
    
    sys.simulate()
        
        
        
        
    qs.LCA(system=sys, lifetime=30, lifetime_unit='yr',)
           # 0.48748859 is the GHG level with the Electricity item from ecoinvent,
           # adjust the electricity amount to reflect different GHG of electricity at different states
           # Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30/0.48748859*elec_GHG)
    
    # =========================================================================
    # TEA
    # =========================================================================

    
    create_tea(sys, IRR_value=0.03,
               income_tax_value=0.21,
               finance_interest_value=0.03)
    
    return sys