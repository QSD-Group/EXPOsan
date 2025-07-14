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
from biosteam import settings
from qsdsan.utils import auom, clear_lca_registries, tea_indices
from exposan.htl import _load_components, landscape_sanunits as lsu, create_tea
from biosteam.units import IsenthalpicValve

_mile_to_km = auom('mile').conversion_factor('km')
_ton_to_tonne = auom('ton').conversion_factor('tonne')
_lb_to_kg = auom('lb').conversion_factor('kg')

# TODO: add in writing: use 2023$ (since CEPCI is not available after that, need a fee to access the full 2024 data)

# GDPCTPI (Gross Domestic Product: Chain-type Price Index)
# https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20)
GDPCTPI = {2014: 96.436,
           2015: 97.277,
           2016: 98.208,
           2017: 100.000,
           2018: 102.290,
           2019: 103.982,
           2020: 105.380,
           2021: 110.173,
           2022: 118.042,
           2023: 122.272,
           2024: 125.231}	

labor_index = tea_indices['labor']

GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                   method='TRACI',
                                   category='environmental impact',
                                   unit='kg CO2-eq',
                                   description='Global Warming Potential')

# TODO: no HXN in all systems

# TODO: address all warnings (especially in T pathways)

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
    'create_T1_system',
    'create_T2_system',
    'create_T3_system',
    'create_T4_system',
    'create_T5_system',
    'create_T6_system',
    'create_T7_system',
    'create_T8_system',
    'create_T9_system',
    'create_T10_system',
    'create_T11_system',
    'create_T12_system',
    'create_T13_system',
    'create_T14_system',
    'create_T15_system'
    )

#%% system C1

# TODO: keep ww_2_dry_sludge_ratio = 1
def create_C1_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.15):
    flowsheet_ID = 'C1'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    # TODO: 0.155 is the average of 77 countries in Lohman et al. may need update later, and need to ensure using the industrial electricity price (may use https://www.globalpetrolprices.com/electricity_prices/)
    bst.PowerUtility.price = 0.155
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.11
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.11
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('landfilled_solids','fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    # TODO: 0.594 is the average of 77 countries in Lohman et al. may need update later
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.594)

    qs.StreamImpactItem(ID='polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='fugitive_methane', linked_stream=stream.fugitive_methane, GlobalWarming=29.8)
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size*ww_2_dry_sludge_ratio/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys, IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C2

def create_C2_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.2):
    flowsheet_ID = 'C2'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    
    AlkalineStabilization = lsu.AlkalineStabilization(ID='AlkalineStabilization',
                                                      ins=(Dewatering-0, 'lime'),
                                                      outs=('stabilized_solids','NG_as_CO2'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=AlkalineStabilization-0,
                                  outs=('landfilled_solids','fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C3

def create_C3_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.2):
    flowsheet_ID = 'C3'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
                     solids_distance=100, FTE=0.2):
    flowsheet_ID = 'C4'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
                     solids_distance=100, FTE=0.3):
    flowsheet_ID = 'C5'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=HeatDrying-0,
                                  outs=('landfilled_solids','fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C6

def create_C6_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.3):
    flowsheet_ID = 'C6'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
                     solids_distance=100, FTE=0.4):
    flowsheet_ID = 'C7'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
                     solids_distance=100, FTE=0.3):
    flowsheet_ID = 'C8'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('landfilled_solids','fugitive_methane_LF','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C9

def create_C9_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.3):
    flowsheet_ID = 'C9'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
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
                      solids_distance=100, FTE=0.35):
    flowsheet_ID = 'C10'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
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
                      solids_distance=100, FTE=0.45):
    flowsheet_ID = 'C11'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=HeatDrying-0,
                                  outs=('landfilled_solids','fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C12

def create_C12_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.45):
    flowsheet_ID = 'C12'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
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
                      solids_distance=100, FTE=0.55):
    flowsheet_ID = 'C13'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
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
                      solids_distance=100, FTE=0.3):
    flowsheet_ID = 'C14'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('landfilled_solids','fugitive_methane_LF','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C15

def create_C15_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.3):
    flowsheet_ID = 'C15'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
                      solids_distance=100, FTE=0.35):
    flowsheet_ID = 'C16'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
                      solids_distance=100, FTE=0.45):
    flowsheet_ID = 'C17'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
                                  outs=('landfilled_solids','fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C18

def create_C18_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.45):
    flowsheet_ID = 'C18'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
                      solids_distance=100, FTE=0.55):
    flowsheet_ID = 'C19'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
                      solids_distance=100, FTE=0.45):
    flowsheet_ID = 'C20'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('landfilled_solids','fugitive_methane_LF','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP ', 'air'),
                                outs=('emission','solid_ash'),
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C21

def create_C21_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.45):
    flowsheet_ID = 'C21'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(Dewatering-0, 'diesel'),
                                          outs=('fugitive_methane_LA','fugitive_nitrous_oxide','carbon_dioxide'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C22

def create_C22_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.5):
    flowsheet_ID = 'C22'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    # TODO: make sure the Composting unit is composting + land application
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel'),
                                outs=('fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C23

def create_C23_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.6):
    flowsheet_ID = 'C23'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
                                  outs=('landfilled_solids','fugitive_methane','fugitive_nitrous_oxide','sequestered_carbon_dioxide'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C24

def create_C24_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.6):
    flowsheet_ID = 'C24'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system C25

def create_C25_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.7):
    flowsheet_ID = 'C25'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
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
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T1

def create_T1_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.55):
    flowsheet_ID = 'T1'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HTL = lsu.HydrothermalLiquefaction(ID='HTL', ins=Dewatering-0,
                                       outs=('hydrochar','HTLaqueous','biocrude','offgas_HTL'))
    
    H2SO4_Tank = qsu.StorageTank(ID='H2SO4_Tank', ins='H2SO4', outs=('H2SO4_out'),
                                 init_with='WasteStream', tau=24, vessel_material='Stainless steel')
    # 0.5 M H2SO4: ~5%
    H2SO4_Tank.ins[0].price = 0.00658 # based on 93% H2SO4 and fresh water (dilute to 5%) price found in Davis 2020$/kg
    H2SO4_Tank.lifetime = 20

    SP1 = qsu.ReversedSplitter(ID='SP1', ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                               init_with='Stream')
    # must put after AcidEx and MemDis in path during simulation to ensure input
    # not empty
    SP1.lifetime = 20
    
    AcidEx = lsu.AcidExtraction(ID='AcidEx', ins=(HTL-0, SP1-0),
                               outs=('residual','extracted'))
    # AcidEx.outs[0].price = -0.055 # SS 2021 SOT PNNL report page 24 Table 9
    # not include residual for TEA and LCA for now
    
    M1_outs1 = AcidEx.outs[1]
    M1 = lsu.HTLmixer(ID='M1', ins=(HTL-1, M1_outs1), outs=('mixture',))
    
    StruPre = lsu.StruvitePrecipitation(ID='StruPre', ins=(M1-0,'MgCl2','NH4Cl','MgO'),
                                        outs=('struvite','CHG_feed'))
    StruPre.ins[1].price = 0.5452
    StruPre.ins[2].price = 0.13
    StruPre.ins[3].price = 0.2
    StruPre.outs[0].price = 0.661
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(StruPre-1, '7.8%_Ru/C'),
                                                outs=('CHG_out','7.8%_Ru/C_out'))
    CHG.ins[1].price = 134.53
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    F1.lifetime = 20
    
    MemDis = lsu.MembraneDistillation(ID='MemDis', ins=(F1-1, SP1-1, 'NaOH'),
                                      outs=('ammonium_sulfate','MemDis_ww','solution'))
    MemDis.ins[2].price = 0.5256
    MemDis.outs[0].price = 0.3236
    
    # TODO: just CHG fuel gas to CHP, not HTL; need to add fugitive emissions for both
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(F1-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T2

def create_T2_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.55):
    flowsheet_ID = 'T2'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HALT = lsu.HydrothermalAlkalineTreatment(ID='HALT', ins=(Dewatering-0, 'sodium_hydroxide', 'hydrochloric_acid'),
                                             outs=('hydrochar','HTLaqueous','biocrude','offgas_HALT'))
    
    H2SO4_Tank = qsu.StorageTank(ID='H2SO4_Tank', ins='H2SO4', outs=('H2SO4_out'),
                                 init_with='WasteStream', tau=24, vessel_material='Stainless steel')
    # 0.5 M H2SO4: ~5%
    H2SO4_Tank.ins[0].price = 0.00658 # based on 93% H2SO4 and fresh water (dilute to 5%) price found in Davis 2020$/kg
    H2SO4_Tank.lifetime = 20

    SP1 = qsu.ReversedSplitter(ID='SP1', ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                               init_with='Stream')
    # must put after AcidEx and MemDis in path during simulation to ensure input
    # not empty
    SP1.lifetime = 20
    
    AcidEx = lsu.AcidExtraction(ID='AcidEx', ins=(HALT-0, SP1-0),
                               outs=('residual','extracted'))
    # AcidEx.outs[0].price = -0.055 # SS 2021 SOT PNNL report page 24 Table 9
    # not include residual for TEA and LCA for now
    
    M1_outs1 = AcidEx.outs[1]
    M1 = lsu.HTLmixer(ID='M1', ins=(HALT-1, M1_outs1), outs=('mixture',))
    
    StruPre = lsu.StruvitePrecipitation(ID='StruPre', ins=(M1-0,'MgCl2','NH4Cl','MgO'),
                                        outs=('struvite','CHG_feed'))
    StruPre.ins[1].price = 0.5452
    StruPre.ins[2].price = 0.13
    StruPre.ins[3].price = 0.2
    StruPre.outs[0].price = 0.661
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(StruPre-1, '7.8%_Ru/C'),
                                                outs=('CHG_out','7.8%_Ru/C_out'))
    CHG.ins[1].price = 134.53
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    F1.lifetime = 20
    
    MemDis = lsu.MembraneDistillation(ID='MemDis', ins=(F1-1, SP1-1, 'NaOH'),
                                      outs=('ammonium_sulfate','MemDis_ww','solution'))
    MemDis.ins[2].price = 0.5256
    MemDis.outs[0].price = 0.3236
    
    # TODO: just CHG fuel gas to CHP, not HTL; need to add fugitive emissions for both
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(F1-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T3

def create_T3_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.4):
    flowsheet_ID = 'T3'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    SCWO = lsu.SupercriticalWaterOxidation(ID='SCWO', ins=Dewatering-0, outs=('ash','offgas_SCWO'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T4

def create_T4_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.7):
    flowsheet_ID = 'T4'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))

    Pyrolysis = lsu.Pyrolysis(ID='Pyrolysis', ins=HeatDrying-0,
                              outs=('biooil','biochar','pyrogas','fugitive_methane','fugitive_nitrous_oxide'))
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(Pyrolysis-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T5

def create_T5_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.7):
    flowsheet_ID = 'T5'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))

    Gasification = lsu.Gasification(ID='Gasification', ins=HeatDrying-0,
                                    outs=('tar','biochar','syngas','fugitive_methane','fugitive_nitrous_oxide'))
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(Gasification-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T6

def create_T6_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.7):
    flowsheet_ID = 'T6'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HTL = lsu.HydrothermalLiquefaction(ID='HTL', ins=Dewatering-0,
                                       outs=('hydrochar','HTLaqueous','biocrude','offgas_HTL'))
    
    H2SO4_Tank = qsu.StorageTank(ID='H2SO4_Tank', ins='H2SO4', outs=('H2SO4_out'),
                                 init_with='WasteStream', tau=24, vessel_material='Stainless steel')
    # 0.5 M H2SO4: ~5%
    H2SO4_Tank.ins[0].price = 0.00658 # based on 93% H2SO4 and fresh water (dilute to 5%) price found in Davis 2020$/kg

    SP1 = qsu.ReversedSplitter(ID='SP1', ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                               init_with='Stream')
    # must put after AcidEx and MemDis in path during simulation to ensure input
    # not empty
    
    AcidEx = lsu.AcidExtraction(ID='AcidEx', ins=(HTL-0, SP1-0),
                               outs=('residual','extracted'))
    # AcidEx.outs[0].price = -0.055 # SS 2021 SOT PNNL report page 24 Table 9
    # not include residual for TEA and LCA for now
    
    M1_outs1 = AcidEx.outs[1]
    M1 = lsu.HTLmixer(ID='M1', ins=(HTL-1, M1_outs1), outs=('mixture',))
    
    StruPre = lsu.StruvitePrecipitation(ID='StruPre', ins=(M1-0,'MgCl2','NH4Cl','MgO'),
                                        outs=('struvite','CHG_feed'))
    StruPre.ins[1].price = 0.5452
    StruPre.ins[2].price = 0.13
    StruPre.ins[3].price = 0.2
    StruPre.outs[0].price = 0.661
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(StruPre-1, '7.8%_Ru/C'),
                                                outs=('CHG_out','7.8%_Ru/C_out'))
    CHG.ins[1].price = 134.53
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    
    MemDis = lsu.MembraneDistillation(ID='MemDis', ins=(F1-1, SP1-1, 'NaOH'),
                                      outs=('ammonium_sulfate','MemDis_ww','solution'))
    MemDis.ins[2].price = 0.5256
    MemDis.outs[0].price = 0.3236
    
    # TODO: just CHG fuel gas to CHP, not HTL; need to add fugitive emissions for both
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(F1-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T7

def create_T7_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.7):
    flowsheet_ID = 'T7'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HALT = lsu.HydrothermalAlkalineTreatment(ID='HALT', ins=(Dewatering-0, 'sodium_hydroxide', 'hydrochloric_acid'),
                                             outs=('hydrochar','HTLaqueous','biocrude','offgas_HALT'))
    
    H2SO4_Tank = qsu.StorageTank(ID='H2SO4_Tank', ins='H2SO4', outs=('H2SO4_out'),
                                 init_with='WasteStream', tau=24, vessel_material='Stainless steel')
    # 0.5 M H2SO4: ~5%
    H2SO4_Tank.ins[0].price = 0.00658 # based on 93% H2SO4 and fresh water (dilute to 5%) price found in Davis 2020$/kg

    SP1 = qsu.ReversedSplitter(ID='SP1', ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                               init_with='Stream')
    # must put after AcidEx and MemDis in path during simulation to ensure input
    # not empty
    
    AcidEx = lsu.AcidExtraction(ID='AcidEx', ins=(HALT-0, SP1-0),
                               outs=('residual','extracted'))
    # AcidEx.outs[0].price = -0.055 # SS 2021 SOT PNNL report page 24 Table 9
    # not include residual for TEA and LCA for now
    
    M1_outs1 = AcidEx.outs[1]
    M1 = lsu.HTLmixer(ID='M1', ins=(HALT-1, M1_outs1), outs=('mixture',))
    
    StruPre = lsu.StruvitePrecipitation(ID='StruPre', ins=(M1-0,'MgCl2','NH4Cl','MgO'),
                                        outs=('struvite','CHG_feed'))
    StruPre.ins[1].price = 0.5452
    StruPre.ins[2].price = 0.13
    StruPre.ins[3].price = 0.2
    StruPre.outs[0].price = 0.661
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(StruPre-1, '7.8%_Ru/C'),
                                                outs=('CHG_out','7.8%_Ru/C_out'))
    CHG.ins[1].price = 134.53
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    
    MemDis = lsu.MembraneDistillation(ID='MemDis', ins=(F1-1, SP1-1, 'NaOH'),
                                      outs=('ammonium_sulfate','MemDis_ww','solution'))
    MemDis.ins[2].price = 0.5256
    MemDis.outs[0].price = 0.3236
    
    # TODO: just CHG fuel gas to CHP, not HTL; need to add fugitive emissions for both
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(F1-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T8

def create_T8_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.55):
    flowsheet_ID = 'T8'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    SCWO = lsu.SupercriticalWaterOxidation(ID='SCWO', ins=Dewatering-0, outs=('ash','offgas_SCWO'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T9

def create_T9_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.85):
    flowsheet_ID = 'T9'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))

    Pyrolysis = lsu.Pyrolysis(ID='Pyrolysis', ins=HeatDrying-0,
                              outs=('biooil','biochar','pyrogas','fugitive_methane','fugitive_nitrous_oxide'))
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(Pyrolysis-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T10

def create_T10_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                     solids_distance=100, FTE=0.85):
    flowsheet_ID = 'T10'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))

    Gasification = lsu.Gasification(ID='Gasification', ins=HeatDrying-0,
                                    outs=('tar','biochar','syngas','fugitive_methane','fugitive_nitrous_oxide'))
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(Gasification-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T11

def create_T11_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.7):
    flowsheet_ID = 'T11'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HTL = lsu.HydrothermalLiquefaction(ID='HTL', ins=Dewatering-0,
                                       outs=('hydrochar','HTLaqueous','biocrude','offgas_HTL'))
    
    H2SO4_Tank = qsu.StorageTank(ID='H2SO4_Tank', ins='H2SO4', outs=('H2SO4_out'),
                                 init_with='WasteStream', tau=24, vessel_material='Stainless steel')
    # 0.5 M H2SO4: ~5%
    H2SO4_Tank.ins[0].price = 0.00658 # based on 93% H2SO4 and fresh water (dilute to 5%) price found in Davis 2020$/kg

    SP1 = qsu.ReversedSplitter(ID='SP1', ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                               init_with='Stream')
    # must put after AcidEx and MemDis in path during simulation to ensure input
    # not empty
    
    AcidEx = lsu.AcidExtraction(ID='AcidEx', ins=(HTL-0, SP1-0),
                               outs=('residual','extracted'))
    # AcidEx.outs[0].price = -0.055 # SS 2021 SOT PNNL report page 24 Table 9
    # not include residual for TEA and LCA for now
    
    M1_outs1 = AcidEx.outs[1]
    M1 = lsu.HTLmixer(ID='M1', ins=(HTL-1, M1_outs1), outs=('mixture',))
    
    StruPre = lsu.StruvitePrecipitation(ID='StruPre', ins=(M1-0,'MgCl2','NH4Cl','MgO'),
                                        outs=('struvite','CHG_feed'))
    StruPre.ins[1].price = 0.5452
    StruPre.ins[2].price = 0.13
    StruPre.ins[3].price = 0.2
    StruPre.outs[0].price = 0.661
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(StruPre-1, '7.8%_Ru/C'),
                                                outs=('CHG_out','7.8%_Ru/C_out'))
    CHG.ins[1].price = 134.53
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    
    MemDis = lsu.MembraneDistillation(ID='MemDis', ins=(F1-1, SP1-1, 'NaOH'),
                                      outs=('ammonium_sulfate','MemDis_ww','solution'))
    MemDis.ins[2].price = 0.5256
    MemDis.outs[0].price = 0.3236
    
    GasMixer = qsu.Mixer('GasMixer', ins=(AnaerobicDigestion-1, F1-0),
                         outs=('fuel_gas'), init_with='WasteStream')
    
    # TODO: just CHG fuel gas to CHP, not HTL; need to add fugitive emissions for both
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T12

def create_T12_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.7):
    flowsheet_ID = 'T12'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HALT = lsu.HydrothermalAlkalineTreatment(ID='HALT', ins=(Dewatering-0, 'sodium_hydroxide', 'hydrochloric_acid'),
                                             outs=('hydrochar','HTLaqueous','biocrude','offgas_HALT'))
    
    H2SO4_Tank = qsu.StorageTank(ID='H2SO4_Tank', ins='H2SO4', outs=('H2SO4_out'),
                                 init_with='WasteStream', tau=24, vessel_material='Stainless steel')
    # 0.5 M H2SO4: ~5%
    H2SO4_Tank.ins[0].price = 0.00658 # based on 93% H2SO4 and fresh water (dilute to 5%) price found in Davis 2020$/kg

    SP1 = qsu.ReversedSplitter(ID='SP1', ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                               init_with='Stream')
    # must put after AcidEx and MemDis in path during simulation to ensure input
    # not empty
    
    AcidEx = lsu.AcidExtraction(ID='AcidEx', ins=(HALT-0, SP1-0),
                               outs=('residual','extracted'))
    # AcidEx.outs[0].price = -0.055 # SS 2021 SOT PNNL report page 24 Table 9
    # not include residual for TEA and LCA for now
    
    M1_outs1 = AcidEx.outs[1]
    M1 = lsu.HTLmixer(ID='M1', ins=(HALT-1, M1_outs1), outs=('mixture',))
    
    StruPre = lsu.StruvitePrecipitation(ID='StruPre', ins=(M1-0,'MgCl2','NH4Cl','MgO'),
                                        outs=('struvite','CHG_feed'))
    StruPre.ins[1].price = 0.5452
    StruPre.ins[2].price = 0.13
    StruPre.ins[3].price = 0.2
    StruPre.outs[0].price = 0.661
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(StruPre-1, '7.8%_Ru/C'),
                                                outs=('CHG_out','7.8%_Ru/C_out'))
    CHG.ins[1].price = 134.53
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    
    MemDis = lsu.MembraneDistillation(ID='MemDis', ins=(F1-1, SP1-1, 'NaOH'),
                                      outs=('ammonium_sulfate','MemDis_ww','solution'))
    MemDis.ins[2].price = 0.5256
    MemDis.outs[0].price = 0.3236
    
    GasMixer = qsu.Mixer('GasMixer', ins=(AnaerobicDigestion-1, F1-0),
                         outs=('fuel_gas'), init_with='WasteStream')
    
    # TODO: just CHG fuel gas to CHP, not HTL; need to add fugitive emissions for both
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T13

def create_T13_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.55):
    flowsheet_ID = 'T13'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    SCWO = lsu.SupercriticalWaterOxidation(ID='SCWO', ins=Dewatering-0, outs=('ash','offgas_SCWO'))
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T14

def create_T14_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.85):
    flowsheet_ID = 'T14'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))

    Pyrolysis = lsu.Pyrolysis(ID='Pyrolysis', ins=HeatDrying-0,
                              outs=('biooil','biochar','pyrogas','fugitive_methane','fugitive_nitrous_oxide'))
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(Pyrolysis-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys

#%% system T15

def create_T15_system(size=10, ww_2_dry_sludge_ratio=1, operation_hours=8760,
                      solids_distance=100, FTE=0.85):
    flowsheet_ID = 'T15'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                   ins=raw_wastewater,
                   outs=('sludge','treated_water'),
                   ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                   sludge_moisture=0.99,
                   sludge_dw_ash=0.436,
                   sludge_afdw_lipid=0.193,
                   sludge_afdw_protein=0.510,
                   sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas','fugitive_methane_AD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=Dewatering-0,
                                outs=('dried_solids','vapor'))

    Gasification = lsu.Gasification(ID='Gasification', ins=HeatDrying-0,
                                    outs=('tar','biochar','syngas','fugitive_methane','fugitive_nitrous_oxide'))
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(Gasification-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.lifetime = 20
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    return sys