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
from qsdsan.utils import auom, clear_lca_registries
from exposan.htl import _load_components, create_tea, _sanunits as su, test_sanunits as tsu

_mile_to_km = auom('mile').conversion_factor('km')

# GDPCTPI (Gross Domestic Product: Chain-type Price Index)
# https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20)
GDPCTPI = {2007: 86.352,
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

__all__ = ('create_landfilling_system','create_land_application_system',)

def create_landfilling_system(size=10,
                              sludge_distance=13,
                              ww_2_dry_sludge_ratio=1,
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
                                  methane_density=0.68,
                                  emission_time=1,
                                  # TODO: remove this later
                                  other_pollutants=[['C2H6',1,28],],
                                  US_state='Illinois',
                                  tipping_fee=None,
                                  operation_hours=operation_hours)
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)

    sys.simulate()

    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    wastewater_solids_hauling = qs.ImpactItem('wastewater_solids_hauling', functional_unit='tonne*km')
    wastewater_solids_hauling.add_indicator(GlobalWarming, 0.13004958)
    # price for per functional unit
    # 4.56 $/m3, 0.072 $/m3/mile (likely 2015$, [2])
    wastewater_solids_hauling.price = (4.56/1040 + 0.072/_mile_to_km/1040*sludge_distance)/\
        GDPCTPI[2015]*GDPCTPI[2022]/sludge_distance*1000
    
    wastewater_solids_transportation = qs.Transportation('wastewater_solids_transportation',
                                                         linked_unit=Landfilling,
                                                         item=wastewater_solids_hauling,
                                                         load_type='mass',
                                                         load=stream.sludge.F_mass/1000,
                                                         load_unit='tonne',
                                                         distance=sludge_distance,
                                                         distance_unit='km',
                                                         # set to 1 h since load = tonne/h
                                                         interval='1',
                                                         interval_unit='h')
    Landfilling.transportation = wastewater_solids_transportation
    
    qs.StreamImpactItem(ID='landfill_methane',
                        linked_stream=stream.landfilling_emission,
                        GlobalWarming=29.8*stream.landfilling_emission.imass['CH4']/stream.landfilling_emission.F_mass)
    
    qs.LCA(system=sys, lifetime=30, lifetime_unit='yr')
    
    create_tea(sys, IRR_value=0.03,
               income_tax_value=0.21,
               finance_interest_value=0.03,
               labor_cost_value=0)
    
    return sys

def create_land_application_system(size=10,
                              sludge_distance=13,
                              ww_2_dry_sludge_ratio=1,
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
    
    LandApplication = tsu.LandApplication(ID='LandApplication',
                                          ins=WWTP-0,
                                          outs=('applied_biosolids','land_application_emission'),
                                          # TODO: remove this later
                                          other_pollutants=[['C2H6',1,28],],
                                          biosolids_grade='B')
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
        
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    biosolids_hauling = qs.ImpactItem('biosolids_hauling', functional_unit='tonne*km')
    biosolids_hauling.add_indicator(GlobalWarming, 0.13004958)
    # price for per functional unit
    # 4.56 $/m3, 0.072 $/m3/mile (likely 2015$, [2])
    biosolids_hauling.price = (4.56/1040 + 0.072/_mile_to_km/1040*sludge_distance)/\
        GDPCTPI[2015]*GDPCTPI[2022]/sludge_distance*1000
    
    biosolids_transportation = qs.Transportation('biosolids_transportation',
                                                 linked_unit=LandApplication,
                                                 item=biosolids_hauling,
                                                 load_type='mass',
                                                 load=stream.sludge.F_mass/1000,
                                                 load_unit='tonne',
                                                 distance=sludge_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = tonne/h
                                                 interval='1',
                                                 interval_unit='h')
    LandApplication.transportation = biosolids_transportation
    
    qs.StreamImpactItem(ID='land_application_nitrous_oxide',
                        linked_stream=stream.land_application_emission,
                        GlobalWarming=273*stream.land_application_emission.imass['N2O']/stream.land_application_emission.F_mass)
    
    qs.LCA(system=sys, lifetime=30, lifetime_unit='yr')
    
    create_tea(sys, IRR_value=0.03,
               income_tax_value=0.21,
               finance_interest_value=0.03,
               labor_cost_value=0)
    
    return sys