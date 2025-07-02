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
from exposan.htl import _load_components, _sanunits as su, landscape_sanunits as lsu, create_tea

_mile_to_km = auom('mile').conversion_factor('km')
_ton_to_tonne = auom('ton').conversion_factor('tonne')
_lb_to_kg = auom('lb').conversion_factor('kg')

# GDPCTPI (Gross Domestic Product: Chain-type Price Index)
# https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20)
GDPCTPI = {1978: 33.339,
           2007: 86.352,
           2008: 87.977,
           2009: 88.557,
           2010: 89.619,
           2011: 91.466,
           2012: 93.176,
           2013: 94.786,
           2014: 96.436,
           2015: 97.277,
           2016: 98.208,
           2017: 100.000,
           2018: 102.290,
           2019: 104.008,
           2020: 105.407,
           2021: 110.220,
           2022: 117.995,
           2023: 122.284}
# TODO: different systems should have different labor costs (and others?)

# TODO: address all warnings
__all__ = (
    'create_C1_system',
    'create_C2_system',
    'create_C3_system',
    'create_C4_system',
    'create_C5_system',
    'create_C6_system',
    'create_C7_system',
    'create_C8_system',
    'create_C9_system',
    'create_C10_system',
    'create_C11_system',
    'create_C12_system',
    'create_C13_system',
    'create_C14_system',
    'create_C15_system',
    'create_C16_system',
    'create_C17_system',
    'create_C18_system',
    'create_C19_system',
    'create_C20_system',
    'create_C21_system',
    'create_C22_system',
    'create_C23_system',
    'create_C24_system',
    'create_C25_system',
    # 'create_T1_system',
    # 'create_T2_system',
    # 'create_T3_system',
    # 'create_T4_system',
    # 'create_T5_system',
    # 'create_T6_system',
    # 'create_T7_system',
    # 'create_T8_system',
    # 'create_T9_system',
    # 'create_T10_system',
    # 'create_T11_system',
    # 'create_T12_system',
    # 'create_T13_system',
    # 'create_T14_system',
    # 'create_T15_system'
    )

#%% system C1

def create_C1_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100):
    flowsheet_ID = 'C1'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C2

def create_C2_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100):
    flowsheet_ID = 'C2'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    
    AlkalineStabilization = lsu.AlkalineStabilization(ID='AlkalineStabilization',
                                                      ins=(Dewatering-0, 'lime'),
                                                      outs=('stabilized_solids','NG_as_CO2'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=AlkalineStabilization-0,
                                  outs=('fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C3

def create_C3_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100):
    flowsheet_ID = 'C3'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    AlkalineStabilization = lsu.AlkalineStabilization(ID='AlkalineStabilization',
                                                      ins=(Dewatering-0, 'lime'),
                                                      outs=('stabilized_solids','NG_as_CO2'))
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(AlkalineStabilization-0, 'diesel'),
                                          outs=('fugitive_methane','fugitive_nitrous_oxide','carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C4

def create_C4_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100):
    flowsheet_ID = 'C4'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    # TODO: make sure the Composting unit is composting + land application
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel'),
                                outs=('fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C5

def create_C5_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100):
    flowsheet_ID = 'C5'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=HeatDrying-0,
                                  outs=('fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C6

def create_C6_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100):
    flowsheet_ID = 'C6'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(HeatDrying-0, 'diesel'),
                                          outs=('fugitive_methane','fugitive_nitrous_oxide','carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C7

def create_C7_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100):
    flowsheet_ID = 'C7'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor_drying'))
    
    Incineration = lsu.Incineration(ID='Incineration',
                                    ins=(HeatDrying-0, 'natural_gas'),
                                    outs=('ash','vapor_incineration','fugitive_methane','fugitive_nitrous_oxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C8

def create_C8_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100):
    flowsheet_ID = 'C8'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('fugitive_methane_LF','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C9

def create_C9_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100):
    flowsheet_ID = 'C9'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(Dewatering-0, 'diesel'),
                                          outs=('fugitive_methane_LA','fugitive_nitrous_oxide','carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C10

def create_C10_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C10'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    # TODO: make sure the Composting unit is composting + land application
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel'),
                                outs=('fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C11

def create_C11_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C11'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=HeatDrying-0,
                                  outs=('fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C12

def create_C12_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C12'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(HeatDrying-0, 'diesel'),
                                          outs=('fugitive_methane','fugitive_nitrous_oxide','carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C13

def create_C13_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C13'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    Incineration = lsu.Incineration(ID='Incineration',
                                    ins=(HeatDrying-0, 'natural_gas_incineration'),
                                    outs=('ash','vapor_incineration','fugitive_methane','fugitive_nitrous_oxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C14

def create_C14_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C14'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('fugitive_methane_LF','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C15

def create_C15_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C15'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(Dewatering-0, 'diesel'),
                                          outs=('fugitive_methane_LA','fugitive_nitrous_oxide','carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C16

def create_C16_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C16'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    # TODO: make sure the Composting unit is composting + land application
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel'),
                                outs=('fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C17

def create_C17_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C17'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=HeatDrying-0,
                                  outs=('fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C18

def create_C18_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C18'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(HeatDrying-0, 'diesel'),
                                          outs=('fugitive_methane','fugitive_nitrous_oxide','carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C19

def create_C19_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C19'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    Incineration = lsu.Incineration(ID='Incineration',
                                    ins=(HeatDrying-0, 'natural_gas_incineration'),
                                    outs=('ash','vapor_incineration','fugitive_methane','fugitive_nitrous_oxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C20

def create_C20_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C20'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('fugitive_methane_LF','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP ', 'air'),
                                outs=('emission','solid_ash'),
                                supplement_power_utility=False)
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C21

def create_C21_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C21'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(Dewatering-0, 'diesel'),
                                          outs=('fugitive_methane_LA','fugitive_nitrous_oxide','carbon_dioxide'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', include_construction=True,
                                ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C22

def create_C22_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C22'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    # TODO: make sure the Composting unit is composting + land application
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel'),
                                outs=('fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', include_construction=True,
                                ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C23

def create_C23_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C23'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=HeatDrying-0,
                                  outs=('fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', include_construction=True,
                                ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C24

def create_C24_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C24'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(HeatDrying-0, 'diesel'),
                                          outs=('fugitive_methane','fugitive_nitrous_oxide','carbon_dioxide'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', include_construction=True,
                                ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C25

def create_C25_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100):
    flowsheet_ID = 'C25'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WWTP = su.WWTP(ID='WWTP',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WWTP-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    Incineration = lsu.Incineration(ID='Incineration',
                                    ins=(HeatDrying-0, 'natural_gas_incineration'),
                                    outs=('ash','vapor_incineration','fugitive_methane','fugitive_nitrous_oxide'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', include_construction=True,
                                ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys