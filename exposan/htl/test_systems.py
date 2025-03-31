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

__all__ = ('create_test_system',)

def create_test_system(size=10,
                       solids_distance=13,
                       ww_2_dry_sludge_ratio=1,
                       lime_type='quick_lime',
                       disposal='landfilling',
                       operation_hours=8760):
        
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
    
    LimeStabilization = tsu.LimeStabilization(ID='LimeStabilization',
                                              ins=(WWTP-0,'lime','water','carbon_dioxide'),
                                              outs=('stablized_solids','vapor'),
                                              lime_type=lime_type, lime_purity=None,
                                              water_ratio=None, vapor_ratio=0.2,
                                              CaOH2_ratio=0.2, unit_electricity=5)
    # TODO: update here if needed
    # from Process Design Manual for Sludge Treatment and Disposal; United States
    # Environmental Protection Agency, 1979:
    # quick lime with a purity of 85% CaO is 40 1978$/ton
    # hydrated lime with a purity of 47% CaO is 44.5 1978$/ton
    if lime_type == 'quick_lime':
        LimeStabilization.ins[1].price = 40/_ton_to_tonne/1000/GDPCTPI[1978]*GDPCTPI[2022]
    else:
        LimeStabilization.ins[1].price = 44.5/_ton_to_tonne/1000/GDPCTPI[1978]*GDPCTPI[2022]
    
    LimeStabilization.ins[2].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    
    HeatDrying = tsu.HeatDrying(ID='HeatDrying',
                                ins=LimeStabilization-0,
                                outs='heated_sludge',
                                target_moisture=0.2,
                                rigorous=True)
        
    Demoisturizer = tsu.Demoisturizer(ID='Demoisturizer',
                                      ins=HeatDrying-0,
                                      outs=('dried_sludge','vapor'))
    
    if disposal == 'landfilling':
        Landfilling = tsu.Landfilling(ID='Landfilling',
                                      ins=Demoisturizer-0,
                                      outs=('landfilled_solids','landfilling_emission'),
                                      methane_density=0.68,
                                      emission_time=1,
                                      US_state='Illinois',
                                      tipping_fee=None,
                                      operation_hours=operation_hours)
        disposal_process = Landfilling
    elif disposal == 'land_application':
        LandApplication = tsu.LandApplication(ID='LandApplication',
                                              ins=Demoisturizer-0,
                                              outs=('applied_biosolids','land_application_emission'),
                                              biosolids_grade='B')
        disposal_process = LandApplication
    
    sys = qs.System.from_units(ID='sys_test',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
        
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    hauling = qs.ImpactItem('hauling', functional_unit='tonne*km')
    hauling.add_indicator(GlobalWarming, 0.13004958)
    # price for per functional unit
    # 4.56 $/m3, 0.072 $/m3/mile (likely 2015$, [2])
    hauling.price = (4.56/1040 + 0.072/_mile_to_km/1040*solids_distance)/\
        GDPCTPI[2015]*GDPCTPI[2022]/solids_distance*1000
    
    solids_hauling = qs.Transportation('solids_hauling',
                                       linked_unit=disposal_process,
                                       item=hauling,
                                       load_type='mass',
                                       load=stream.dried_sludge.F_mass/1000,
                                       load_unit='tonne',
                                       distance=solids_distance,
                                       distance_unit='km',
                                       # set to 1 h since load = tonne/h
                                       interval='1',
                                       interval_unit='h')
    disposal_process.transportation = solids_hauling
    
    # TODO: add LCA data, pay attention to purity, should be consistent with price
    # ecoinvent lacks purity data (market for lime, market for lime, hydrated, packed, both using RER)
    # if lime_type == 'quick_lime':
        # qs.StreamImpactItem(ID='quick_lime',
        #                     linked_stream=stream.lime,
        #                     GlobalWarming=)
    # else:
        # qs.StreamImpactItem(ID='hydrated_lime',
        #                     linked_stream=stream.lime,
        #                     GlobalWarming=)
    
    if disposal == 'landfilling':
        qs.StreamImpactItem(ID='landfill_methane',
                            linked_stream=stream.landfilling_emission,
                            GlobalWarming=29.8*stream.landfilling_emission.imass['CH4']/stream.landfilling_emission.F_mass)
        
    elif disposal == 'land_application':
        qs.StreamImpactItem(ID='land_application_nitrous_oxide',
                            linked_stream=stream.land_application_emission,
                            GlobalWarming=273*stream.land_application_emission.imass['N2O']/stream.land_application_emission.F_mass)
    
    qs.LCA(system=sys, lifetime=30, lifetime_unit='yr')
    
    create_tea(sys, IRR_value=0.03,
               income_tax_value=0.21,
               finance_interest_value=0.03,
               labor_cost_value=0)
    
    return sys