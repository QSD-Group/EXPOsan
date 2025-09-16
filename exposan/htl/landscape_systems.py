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

#%% housekeeping

# TODO: LCA data update
# TODO 1: update lime LCA data from ecoinvent 3.8 to 3.11
# TODO 2: check why the previous code has different CI for transporting sludge and biooils
# TODO 3: country-level electricty CI data need more digits (use the code below to process ecoinvent data)
# =============================================================================
# data = pd.read_excel('/Users/jiananfeng/Desktop/XXX.xlsx')
# data['country_code'] = data['country'].str.extract(r'\((.*?)\)')
# data['country'] = data['country'].str.replace(r"\s*\(.*?\)", "", regex=True)
# data.to_excel('/Users/jiananfeng/Desktop/XXX.xlsx')
# =============================================================================

# TODO: TODOs in _components.py (and maybe others as well)

# TODO: need to decide how to integrate the following TODO into the roadmap analysis
# TODO: consider update TEA parameters for pioneer T systems based on the suggestions in Poddar, T. K.; Scown, C. D. Technoeconomic Analysis for Near-Term Scale-up of Bioprocesses. Current Opinion in Biotechnology 2025, 92, 103258. https://doi.org/10.1016/j.copbio.2025.103258

# TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)

# TODO: https://www.globalpetrolprices.com/electricity_prices/ also has prices for electricty, diesel, natural gas, etc.
# TODO: make sure to use the price for the industrial sector
# TODO: 0.16 is the median electricity of of 77 countries in Lohman et al. may need update later, and need to ensure using the industrial electricity price (may use https://www.globalpetrolprices.com/electricity_prices/)

# TODO: check all parameters in TEA, especially the ones listed here
# TODO: consider adding IRR as a WRRF typology parameter?
# TODO: income tax may be found in Steward et al. ES&T, 2023

#%% initialization

import numpy as np, qsdsan as qs, biosteam as bst
from qsdsan import sanunits as qsu
from biosteam import settings
from qsdsan.utils import auom, clear_lca_registries, tea_indices
from exposan.htl import _load_components, landscape_sanunits as lsu, create_tea
from biosteam.units import IsenthalpicValve

# use methane density, probably consistent with BioSTEAM, kg/m3
# https://en.wikipedia.org/wiki/Methane (accessed 2025-02-10)
natural_gas_density = 0.657

# diesel density, kg/m3
diesel_density = 850

_mile_to_km = auom('mile').conversion_factor('km')
_ton_to_tonne = auom('ton').conversion_factor('tonne')
_lb_to_kg = auom('lb').conversion_factor('kg')
_gal_to_liter = auom('gallon').conversion_factor('liter')
_oil_barrel_to_m3 = auom('oil_barrel').conversion_factor('m3')

# GDPCTPI (Gross Domestic Product: Chain-type Price Index)
# https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20)
GDPCTPI = {2005: 81.537,
           2006: 84.075,
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
           2019: 103.982,
           2020: 105.380,
           2021: 110.173,
           2022: 118.042,
           2023: 122.272,
           2024: 125.231}	

labor_index = tea_indices['labor']

# TODO: try addressing all warnings (especially in T pathways)

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

def create_C1_system(size=10, operation_hours=7884, LF_distance=100, tipping_fee=0.0626, FTE=0.15, lifetime=20):
    flowsheet_ID = 'C1'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: landfilling tipping fees can be an uncertainty parameter or a typology parameter
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('landfilled_solids','methane_LF','nitrous_oxide_LF','sequestered_carbon_dioxide_LF'),
                                  solids_distance=LF_distance)
    # see exposan/htl/data/landfilling_tipping_fee.xlsx
    # U.S. average: 56.8 $·wet ton-1 ~ 62.6 $·wet tonne-1
    Landfilling.outs[0].price = -tipping_fee
    
    sys = qs.System.from_units(ID='system_C1',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    sludge_trucking.price = (0.00551 + 0.0000541*Landfilling.solids_distance)/Landfilling.solids_distance
    
    sludge_transportation = qs.Transportation(ID='Sludge_transportation',
                                              linked_unit=Landfilling,
                                              item=sludge_trucking,
                                              load_type='mass',
                                              load=stream.landfilled_solids.F_mass,
                                              load_unit='kg',
                                              distance=Landfilling.solids_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Landfilling.transportation.append(sludge_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C2

def create_C2_system(size=10, operation_hours=7884, LF_distance=100, tipping_fee=0.0626, FTE=0.2, lifetime=20):
    flowsheet_ID = 'C2'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AlkalineStabilization = lsu.AlkalineStabilization(ID='AlkalineStabilization',
                                                      ins=(Dewatering-0, 'lime'),
                                                      outs='stabilized_solids')
    AlkalineStabilization.ins[1].price = 0.1189/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=AlkalineStabilization-0,
                                  outs=('landfilled_solids','methane_LF','nitrous_oxide_LF','sequestered_carbon_dioxide_LF'),
                                  solids_distance=LF_distance)
    # see exposan/htl/data/landfilling_tipping_fee.xlsx
    # U.S. average: 56.8 $·wet ton-1 ~ 62.6 $·wet tonne-1
    Landfilling.outs[0].price = -tipping_fee
    
    sys = qs.System.from_units(ID='system_C2',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # TODO: discuss with J.S.G.: whehter the CO2 emission from limestone decomposition is fossil-origin
    qs.StreamImpactItem(ID='Lime', linked_stream=stream.lime, GlobalWarming=0.04261602+44/56)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    sludge_trucking.price = (0.00551 + 0.0000541*Landfilling.solids_distance)/Landfilling.solids_distance
    
    sludge_transportation = qs.Transportation(ID='Sludge_transportation',
                                              linked_unit=Landfilling,
                                              item=sludge_trucking,
                                              load_type='mass',
                                              load=stream.landfilled_solids.F_mass,
                                              load_unit='kg',
                                              distance=Landfilling.solids_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Landfilling.transportation.append(sludge_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C3

def create_C3_system(size=10, operation_hours=7884, LA_distance=100, biosolids_price=0, FTE=0.2, lifetime=20):
    flowsheet_ID = 'C3'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AlkalineStabilization = lsu.AlkalineStabilization(ID='AlkalineStabilization',
                                                      ins=(Dewatering-0, 'lime'),
                                                      outs='stabilized_solids')
    AlkalineStabilization.ins[1].price = 0.1189/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(AlkalineStabilization-0, 'diesel_LA'),
                                          outs=('biosolids_cost','methane_LA','nitrous_oxide_LA','carbon_dioxide_LA'),
                                          solids_distance=LA_distance)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    LandApplication.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    LandApplication.outs[0].price = biosolids_price
    
    sys = qs.System.from_units(ID='system_C3',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # TODO: discuss with J.S.G.: whehter the CO2 emission from limestone decomposition is fossil-origin
    qs.StreamImpactItem(ID='Lime', linked_stream=stream.lime, GlobalWarming=0.04261602+44/56)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    biosolids_trucking.price = (0.00551 + 0.0000541*LandApplication.solids_distance)/LandApplication.solids_distance
    
    biosolids_transportation = qs.Transportation(ID='Biosolids_transportation',
                                                 linked_unit=LandApplication,
                                                 item=biosolids_trucking,
                                                 load_type='mass',
                                                 load=stream.biosolids_cost.F_mass,
                                                 load_unit='kg',
                                                 distance=LandApplication.solids_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = kg/h
                                                 interval='1',
                                                 interval_unit='h')
    LandApplication.transportation.append(biosolids_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C4

def create_C4_system(size=10, operation_hours=7884, LA_distance=100, compost_price=0.05, FTE=0.2, lifetime=20):
    flowsheet_ID = 'C4'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel_composting'),
                                outs=('compost_cost','methane_composting','nitrous_oxide_composting','sequestered_carbon_dioxide_composting'),
                                feedstock_digested=False, solids_distance=LA_distance)
    # TODO: uncertainty range (uniform) 18-36 2005$/tonne
    Composting.ins[1].price = 27/1000/GDPCTPI[2005]*GDPCTPI[2023]
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    Composting.ins[2].price = 4.224/_gal_to_liter*1000/diesel_density
    Composting.outs[0].price = compost_price
    
    sys = qs.System.from_units(ID='system_C4',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Bulking_agent = qs.ImpactItem('Bulking_agent', functional_unit='kg')
    Bulking_agent.add_indicator(GlobalWarming, 0.04291794596775965)
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_composting', linked_stream=stream.diesel_composting, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_composting', linked_stream=stream.methane_composting, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_composting', linked_stream=stream.nitrous_oxide_composting, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_composting', linked_stream=stream.sequestered_carbon_dioxide_composting, GlobalWarming=1)
    
    compost_trucking = qs.ImpactItem(ID='Compost_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    compost_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    compost_trucking.price = (0.00551 + 0.0000541*Composting.solids_distance)/Composting.solids_distance
    
    compost_transportation = qs.Transportation(ID='Compost_transportation',
                                               linked_unit=Composting,
                                               item=compost_trucking,
                                               load_type='mass',
                                               load=stream.compost_cost.F_mass,
                                               load_unit='kg',
                                               distance=Composting.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    Composting.transportation.append(compost_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Bulking_agent=lambda:sys.flowsheet.Composting.ins[1].imass['Sawdust']*operation_hours*lifetime,
           N_fertilizer=lambda:sys.flowsheet.Composting.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.Composting.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C5

def create_C5_system(size=10, operation_hours=7884, LF_distance=100, tipping_fee=0.0626, FTE=0.3, lifetime=20):
    flowsheet_ID = 'C5'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=HeatDrying-0,
                                  outs=('landfilled_solids','methane_LF','nitrous_oxide_LF','sequestered_carbon_dioxide_LF'),
                                  solids_distance=LF_distance)
    # see exposan/htl/data/landfilling_tipping_fee.xlsx
    # U.S. average: 56.8 $·wet ton-1 ~ 62.6 $·wet tonne-1
    Landfilling.outs[0].price = -tipping_fee
    
    sys = qs.System.from_units(ID='system_C5',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    sludge_trucking.price = (0.00551 + 0.0000541*Landfilling.solids_distance)/Landfilling.solids_distance
    
    sludge_transportation = qs.Transportation(ID='Sludge_transportation',
                                              linked_unit=Landfilling,
                                              item=sludge_trucking,
                                              load_type='mass',
                                              load=stream.landfilled_solids.F_mass,
                                              load_unit='kg',
                                              distance=Landfilling.solids_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Landfilling.transportation.append(sludge_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:sys.flowsheet.HeatDrying.ins[1].F_vol*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C6

def create_C6_system(size=10, operation_hours=7884, LA_distance=100, biosolids_price=0, FTE=0.3, lifetime=20):
    flowsheet_ID = 'C6'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(HeatDrying-0, 'diesel_LA'),
                                          outs=('biosolids_cost','methane_LA','nitrous_oxide_LA','carbon_dioxide_LA'),
                                          solids_distance=LA_distance)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    LandApplication.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    LandApplication.outs[0].price = biosolids_price
    
    sys = qs.System.from_units(ID='system_C6',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    biosolids_trucking.price = (0.00551 + 0.0000541*LandApplication.solids_distance)/LandApplication.solids_distance
    
    biosolids_transportation = qs.Transportation(ID='Biosolids_transportation',
                                                 linked_unit=LandApplication,
                                                 item=biosolids_trucking,
                                                 load_type='mass',
                                                 load=stream.biosolids_cost.F_mass,
                                                 load_unit='kg',
                                                 distance=LandApplication.solids_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = kg/h
                                                 interval='1',
                                                 interval_unit='h')
    LandApplication.transportation.append(biosolids_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:sys.flowsheet.HeatDrying.ins[1].F_vol*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C7

def create_C7_system(size=10, operation_hours=7884, FTE=0.4, lifetime=20):
    flowsheet_ID = 'C7'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor_drying'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Incineration = lsu.Incineration(ID='Incineration',
                                    ins=(HeatDrying-0, 'natural_gas_incineration'),
                                    outs=('ash_incineration','vapor_incineration','methane_IN','nitrous_oxide_IN'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    Incineration.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    Incineration.outs[0].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    sys = qs.System.from_units(ID='system_C7',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Ash_incineration', linked_stream=stream.ash_incineration, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_IN', linked_stream=stream.methane_IN, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_IN', linked_stream=stream.nitrous_oxide_IN, GlobalWarming=273)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.Incineration.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C8

def create_C8_system(size=10, operation_hours=7884, LF_distance=100, tipping_fee=0.0626, FTE=0.3, lifetime=20):
    flowsheet_ID = 'C8'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('landfilled_solids','methane_LF','nitrous_oxide_LF','sequestered_carbon_dioxide_LF'),
                                  solids_distance=LF_distance)
    # see exposan/htl/data/landfilling_tipping_fee.xlsx
    # U.S. average: 56.8 $·wet ton-1 ~ 62.6 $·wet tonne-1
    Landfilling.outs[0].price = -tipping_fee
    
    sys = qs.System.from_units(ID='system_C8',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    sludge_trucking.price = (0.00551 + 0.0000541*Landfilling.solids_distance)/Landfilling.solids_distance
    
    sludge_transportation = qs.Transportation(ID='Sludge_transportation',
                                              linked_unit=Landfilling,
                                              item=sludge_trucking,
                                              load_type='mass',
                                              load=stream.landfilled_solids.F_mass,
                                              load_unit='kg',
                                              distance=Landfilling.solids_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Landfilling.transportation.append(sludge_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C9

def create_C9_system(size=10, operation_hours=7884, LA_distance=100, biosolids_price=0, FTE=0.3, lifetime=20):
    flowsheet_ID = 'C9'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(Dewatering-0, 'diesel_LA'),
                                          outs=('biosolids_cost','methane_LA','nitrous_oxide_LA','carbon_dioxide_LA'),
                                          solids_distance=LA_distance)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    LandApplication.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    LandApplication.outs[0].price = biosolids_price
    
    sys = qs.System.from_units(ID='system_C9',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    biosolids_trucking.price = (0.00551 + 0.0000541*LandApplication.solids_distance)/LandApplication.solids_distance
    
    biosolids_transportation = qs.Transportation(ID='Biosolids_transportation',
                                                 linked_unit=LandApplication,
                                                 item=biosolids_trucking,
                                                 load_type='mass',
                                                 load=stream.biosolids_cost.F_mass,
                                                 load_unit='kg',
                                                 distance=LandApplication.solids_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = kg/h
                                                 interval='1',
                                                 interval_unit='h')
    LandApplication.transportation.append(biosolids_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C10

def create_C10_system(size=10, operation_hours=7884, LA_distance=100, compost_price=0.05, FTE=0.35, lifetime=20):
    flowsheet_ID = 'C10'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel_composting'),
                                outs=('compost_cost','methane_composting','nitrous_oxide_composting','sequestered_carbon_dioxide_composting'),
                                feedstock_digested=True, solids_distance=LA_distance)
    # TODO: uncertainty range (uniform) 18-36 2005$/tonne
    Composting.ins[1].price = 27/1000/GDPCTPI[2005]*GDPCTPI[2023]
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    Composting.ins[2].price = 4.224/_gal_to_liter*1000/diesel_density
    Composting.outs[0].price = compost_price
    
    sys = qs.System.from_units(ID='system_C10',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Bulking_agent = qs.ImpactItem('Bulking_agent', functional_unit='kg')
    Bulking_agent.add_indicator(GlobalWarming, 0.04291794596775965)
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_composting', linked_stream=stream.diesel_composting, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_composting', linked_stream=stream.methane_composting, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_composting', linked_stream=stream.nitrous_oxide_composting, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_composting', linked_stream=stream.sequestered_carbon_dioxide_composting, GlobalWarming=1)
    
    compost_trucking = qs.ImpactItem(ID='Compost_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    compost_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    compost_trucking.price = (0.00551 + 0.0000541*Composting.solids_distance)/Composting.solids_distance
    
    compost_transportation = qs.Transportation(ID='Compost_transportation',
                                               linked_unit=Composting,
                                               item=compost_trucking,
                                               load_type='mass',
                                               load=stream.compost_cost.F_mass,
                                               load_unit='kg',
                                               distance=Composting.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    Composting.transportation.append(compost_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Bulking_agent=lambda:sys.flowsheet.Composting.ins[1].imass['Sawdust']*operation_hours*lifetime,
           N_fertilizer=lambda:sys.flowsheet.Composting.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.Composting.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C11

def create_C11_system(size=10, operation_hours=7884, LF_distance=100, tipping_fee=0.0626, FTE=0.45, lifetime=20):
    flowsheet_ID = 'C11'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=HeatDrying-0,
                                  outs=('landfilled_solids','methane_LF','nitrous_oxide_LF','sequestered_carbon_dioxide_LF'),
                                  solids_distance=LF_distance)
    # see exposan/htl/data/landfilling_tipping_fee.xlsx
    # U.S. average: 56.8 $·wet ton-1 ~ 62.6 $·wet tonne-1
    Landfilling.outs[0].price = -tipping_fee
    
    sys = qs.System.from_units(ID='system_C11',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    sludge_trucking.price = (0.00551 + 0.0000541*Landfilling.solids_distance)/Landfilling.solids_distance
    
    sludge_transportation = qs.Transportation(ID='Sludge_transportation',
                                              linked_unit=Landfilling,
                                              item=sludge_trucking,
                                              load_type='mass',
                                              load=stream.landfilled_solids.F_mass,
                                              load_unit='kg',
                                              distance=Landfilling.solids_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Landfilling.transportation.append(sludge_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:sys.flowsheet.HeatDrying.ins[1].F_vol*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C12

def create_C12_system(size=10, operation_hours=7884, LA_distance=100, biosolids_price=0, FTE=0.45, lifetime=20):
    flowsheet_ID = 'C12'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(HeatDrying-0, 'diesel_LA'),
                                          outs=('biosolids_cost','methane_LA','nitrous_oxide_LA','carbon_dioxide_LA'),
                                          solids_distance=LA_distance)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    LandApplication.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    LandApplication.outs[0].price = biosolids_price
    
    sys = qs.System.from_units(ID='system_C12',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    biosolids_trucking.price = (0.00551 + 0.0000541*LandApplication.solids_distance)/LandApplication.solids_distance
    
    biosolids_transportation = qs.Transportation(ID='Biosolids_transportation',
                                                 linked_unit=LandApplication,
                                                 item=biosolids_trucking,
                                                 load_type='mass',
                                                 load=stream.biosolids_cost.F_mass,
                                                 load_unit='kg',
                                                 distance=LandApplication.solids_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = kg/h
                                                 interval='1',
                                                 interval_unit='h')
    LandApplication.transportation.append(biosolids_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:sys.flowsheet.HeatDrying.ins[1].F_vol*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C13

def create_C13_system(size=10, operation_hours=7884, FTE=0.55, lifetime=20):
    flowsheet_ID = 'C13'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Incineration = lsu.Incineration(ID='Incineration',
                                    ins=(HeatDrying-0, 'natural_gas_incineration'),
                                    outs=('ash_incineration','vapor_incineration','methane_IN','nitrous_oxide_IN'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    Incineration.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    Incineration.outs[0].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    sys = qs.System.from_units(ID='system_C13',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Ash_incineration', linked_stream=stream.ash_incineration, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_IN', linked_stream=stream.methane_IN, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_IN', linked_stream=stream.nitrous_oxide_IN, GlobalWarming=273)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.Incineration.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C14

def create_C14_system(size=10, operation_hours=7884, LF_distance=100, tipping_fee=0.0626, FTE=0.3, lifetime=20):
    flowsheet_ID = 'C14'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0, biogas_RNG_ratio=0.9)
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    AnaerobicDigestion.outs[1].price = 0.218
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('landfilled_solids','methane_LF','nitrous_oxide_LF','sequestered_carbon_dioxide_LF'),
                                  solids_distance=LF_distance)
    # see exposan/htl/data/landfilling_tipping_fee.xlsx
    # U.S. average: 56.8 $·wet ton-1 ~ 62.6 $·wet tonne-1
    Landfilling.outs[0].price = -tipping_fee
    
    sys = qs.System.from_units(ID='system_C14',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    sludge_trucking.price = (0.00551 + 0.0000541*Landfilling.solids_distance)/Landfilling.solids_distance
    
    sludge_transportation = qs.Transportation(ID='Sludge_transportation',
                                              linked_unit=Landfilling,
                                              item=sludge_trucking,
                                              load_type='mass',
                                              load=stream.landfilled_solids.F_mass,
                                              load_unit='kg',
                                              distance=Landfilling.solids_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Landfilling.transportation.append(sludge_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C15

def create_C15_system(size=10, operation_hours=7884, LA_distance=100, biosolids_price=0, FTE=0.3, lifetime=20):
    flowsheet_ID = 'C15'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0, biogas_RNG_ratio=0.9)
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    AnaerobicDigestion.outs[1].price = 0.218
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(Dewatering-0, 'diesel_LA'),
                                          outs=('biosolids_cost','methane_LA','nitrous_oxide_LA','carbon_dioxide_LA'),
                                          solids_distance=LA_distance)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    LandApplication.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    LandApplication.outs[0].price = biosolids_price
    
    sys = qs.System.from_units(ID='system_C15',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    biosolids_trucking.price = (0.00551 + 0.0000541*LandApplication.solids_distance)/LandApplication.solids_distance
    
    biosolids_transportation = qs.Transportation(ID='Biosolids_transportation',
                                                 linked_unit=LandApplication,
                                                 item=biosolids_trucking,
                                                 load_type='mass',
                                                 load=stream.biosolids_cost.F_mass,
                                                 load_unit='kg',
                                                 distance=LandApplication.solids_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = kg/h
                                                 interval='1',
                                                 interval_unit='h')
    LandApplication.transportation.append(biosolids_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C16

def create_C16_system(size=10, operation_hours=7884, LA_distance=100, compost_price=0.05, FTE=0.35, lifetime=20):
    flowsheet_ID = 'C16'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0, biogas_RNG_ratio=0.9)
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    AnaerobicDigestion.outs[1].price = 0.218
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel_composting'),
                                outs=('compost_cost','methane_composting','nitrous_oxide_composting','sequestered_carbon_dioxide_composting'),
                                feedstock_digested=True, solids_distance=LA_distance)
    # TODO: uncertainty range (uniform) 18-36 2005$/tonne
    Composting.ins[1].price = 27/1000/GDPCTPI[2005]*GDPCTPI[2023]
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    Composting.ins[2].price = 4.224/_gal_to_liter*1000/diesel_density
    Composting.outs[0].price = compost_price
    
    sys = qs.System.from_units(ID='system_C16',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Bulking_agent = qs.ImpactItem('Bulking_agent', functional_unit='kg')
    Bulking_agent.add_indicator(GlobalWarming, 0.04291794596775965)
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_composting', linked_stream=stream.diesel_composting, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_composting', linked_stream=stream.methane_composting, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_composting', linked_stream=stream.nitrous_oxide_composting, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_composting', linked_stream=stream.sequestered_carbon_dioxide_composting, GlobalWarming=1)
    
    compost_trucking = qs.ImpactItem(ID='Compost_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    compost_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    compost_trucking.price = (0.00551 + 0.0000541*Composting.solids_distance)/Composting.solids_distance
    
    compost_transportation = qs.Transportation(ID='Compost_transportation',
                                               linked_unit=Composting,
                                               item=compost_trucking,
                                               load_type='mass',
                                               load=stream.compost_cost.F_mass,
                                               load_unit='kg',
                                               distance=Composting.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    Composting.transportation.append(compost_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Bulking_agent=lambda:sys.flowsheet.Composting.ins[1].imass['Sawdust']*operation_hours*lifetime,
           N_fertilizer=lambda:sys.flowsheet.Composting.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.Composting.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C17

def create_C17_system(size=10, operation_hours=7884, LF_distance=100, tipping_fee=0.0626, FTE=0.45, lifetime=20):
    flowsheet_ID = 'C17'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0, biogas_RNG_ratio=0.9)
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    AnaerobicDigestion.outs[1].price = 0.218
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=HeatDrying-0,
                                  outs=('landfilled_solids','methane_LF','nitrous_oxide_LF','sequestered_carbon_dioxide_LF'),
                                  solids_distance=LF_distance)
    # see exposan/htl/data/landfilling_tipping_fee.xlsx
    # U.S. average: 56.8 $·wet ton-1 ~ 62.6 $·wet tonne-1
    Landfilling.outs[0].price = -tipping_fee
    
    sys = qs.System.from_units(ID='system_C17',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    sludge_trucking.price = (0.00551 + 0.0000541*Landfilling.solids_distance)/Landfilling.solids_distance
    
    sludge_transportation = qs.Transportation(ID='Sludge_transportation',
                                              linked_unit=Landfilling,
                                              item=sludge_trucking,
                                              load_type='mass',
                                              load=stream.landfilled_solids.F_mass,
                                              load_unit='kg',
                                              distance=Landfilling.solids_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Landfilling.transportation.append(sludge_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol+sys.flowsheet.HeatDrying.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C18

def create_C18_system(size=10, operation_hours=7884, LA_distance=100, biosolids_price=0, FTE=0.45, lifetime=20):
    flowsheet_ID = 'C18'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0, biogas_RNG_ratio=0.9)
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    AnaerobicDigestion.outs[1].price = 0.218
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(HeatDrying-0, 'diesel_LA'),
                                          outs=('biosolids_cost','methane_LA','nitrous_oxide_LA','carbon_dioxide_LA'),
                                          solids_distance=LA_distance)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    LandApplication.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    LandApplication.outs[0].price = biosolids_price
    
    sys = qs.System.from_units(ID='system_C18',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    biosolids_trucking.price = (0.00551 + 0.0000541*LandApplication.solids_distance)/LandApplication.solids_distance
    
    biosolids_transportation = qs.Transportation(ID='Biosolids_transportation',
                                                 linked_unit=LandApplication,
                                                 item=biosolids_trucking,
                                                 load_type='mass',
                                                 load=stream.biosolids_cost.F_mass,
                                                 load_unit='kg',
                                                 distance=LandApplication.solids_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = kg/h
                                                 interval='1',
                                                 interval_unit='h')
    LandApplication.transportation.append(biosolids_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol+sys.flowsheet.HeatDrying.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C19

def create_C19_system(size=10, operation_hours=7884, FTE=0.55, lifetime=20):
    flowsheet_ID = 'C19'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0, biogas_RNG_ratio=0.9)
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    AnaerobicDigestion.outs[1].price = 0.218
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Incineration = lsu.Incineration(ID='Incineration',
                                    ins=(HeatDrying-0, 'natural_gas_incineration'),
                                    outs=('ash_incineration','vapor_incineration','methane_IN','nitrous_oxide_IN'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    Incineration.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    Incineration.outs[0].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    sys = qs.System.from_units(ID='system_C19',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Ash_incineration', linked_stream=stream.ash_incineration, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_IN', linked_stream=stream.methane_IN, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_IN', linked_stream=stream.nitrous_oxide_IN, GlobalWarming=273)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol+sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.Incineration.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C20

def create_C20_system(size=10, operation_hours=7884, LF_distance=100, tipping_fee=0.0626, FTE=0.45, lifetime=20):
    flowsheet_ID = 'C20'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('landfilled_solids','methane_LF','nitrous_oxide_LF','sequestered_carbon_dioxide_LF'),
                                  solids_distance=LF_distance)
    # see exposan/htl/data/landfilling_tipping_fee.xlsx
    # U.S. average: 56.8 $·wet ton-1 ~ 62.6 $·wet tonne-1
    Landfilling.outs[0].price = -tipping_fee
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    sys = qs.System.from_units(ID='system_C20',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    sludge_trucking.price = (0.00551 + 0.0000541*Landfilling.solids_distance)/Landfilling.solids_distance
    
    sludge_transportation = qs.Transportation(ID='Sludge_transportation',
                                              linked_unit=Landfilling,
                                              item=sludge_trucking,
                                              load_type='mass',
                                              load=stream.landfilled_solids.F_mass,
                                              load_unit='kg',
                                              distance=Landfilling.solids_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Landfilling.transportation.append(sludge_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C21

def create_C21_system(size=10, operation_hours=7884, LA_distance=100, biosolids_price=0, FTE=0.45, lifetime=20):
    flowsheet_ID = 'C21'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(Dewatering-0, 'diesel_LA'),
                                          outs=('biosolids_cost','methane_LA','nitrous_oxide_LA','carbon_dioxide_LA'),
                                          solids_distance=LA_distance)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    LandApplication.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    LandApplication.outs[0].price = biosolids_price
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    sys = qs.System.from_units(ID='system_C21',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    biosolids_trucking.price = (0.00551 + 0.0000541*LandApplication.solids_distance)/LandApplication.solids_distance
    
    biosolids_transportation = qs.Transportation(ID='Biosolids_transportation',
                                                 linked_unit=LandApplication,
                                                 item=biosolids_trucking,
                                                 load_type='mass',
                                                 load=stream.biosolids_cost.F_mass,
                                                 load_unit='kg',
                                                 distance=LandApplication.solids_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = kg/h
                                                 interval='1',
                                                 interval_unit='h')
    LandApplication.transportation.append(biosolids_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C22

def create_C22_system(size=10, operation_hours=7884, LA_distance=100, compost_price=0.05, FTE=0.5, lifetime=20):
    flowsheet_ID = 'C22'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel_composting'),
                                outs=('compost_cost','methane_composting','nitrous_oxide_composting','sequestered_carbon_dioxide_composting'),
                                feedstock_digested=True, solids_distance=LA_distance)
    # TODO: uncertainty range (uniform) 18-36 2005$/tonne
    Composting.ins[1].price = 27/1000/GDPCTPI[2005]*GDPCTPI[2023]
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    Composting.ins[2].price = 4.224/_gal_to_liter*1000/diesel_density
    Composting.outs[0].price = compost_price
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    sys = qs.System.from_units(ID='system_C22',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Bulking_agent = qs.ImpactItem('Bulking_agent', functional_unit='kg')
    Bulking_agent.add_indicator(GlobalWarming, 0.04291794596775965)
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_composting', linked_stream=stream.diesel_composting, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_composting', linked_stream=stream.methane_composting, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_composting', linked_stream=stream.nitrous_oxide_composting, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_composting', linked_stream=stream.sequestered_carbon_dioxide_composting, GlobalWarming=1)
    
    compost_trucking = qs.ImpactItem(ID='Compost_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    compost_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    compost_trucking.price = (0.00551 + 0.0000541*Composting.solids_distance)/Composting.solids_distance
    
    compost_transportation = qs.Transportation(ID='Compost_transportation',
                                               linked_unit=Composting,
                                               item=compost_trucking,
                                               load_type='mass',
                                               load=stream.compost_cost.F_mass,
                                               load_unit='kg',
                                               distance=Composting.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    Composting.transportation.append(compost_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Bulking_agent=lambda:sys.flowsheet.Composting.ins[1].imass['Sawdust']*operation_hours*lifetime,
           N_fertilizer=lambda:sys.flowsheet.Composting.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.Composting.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C23

def create_C23_system(size=10, operation_hours=7884, LF_distance=100, tipping_fee=0.0626, FTE=0.6, lifetime=20):
    flowsheet_ID = 'C23'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=HeatDrying-0,
                                  outs=('landfilled_solids','methane_LF','nitrous_oxide_LF','sequestered_carbon_dioxide_LF'),
                                  solids_distance=LF_distance)
    # see exposan/htl/data/landfilling_tipping_fee.xlsx
    # U.S. average: 56.8 $·wet ton-1 ~ 62.6 $·wet tonne-1
    Landfilling.outs[0].price = -tipping_fee
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    sys = qs.System.from_units(ID='system_C23',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    sludge_trucking.price = (0.00551 + 0.0000541*Landfilling.solids_distance)/Landfilling.solids_distance
    
    sludge_transportation = qs.Transportation(ID='Sludge_transportation',
                                              linked_unit=Landfilling,
                                              item=sludge_trucking,
                                              load_type='mass',
                                              load=stream.landfilled_solids.F_mass,
                                              load_unit='kg',
                                              distance=Landfilling.solids_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Landfilling.transportation.append(sludge_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C24

def create_C24_system(size=10, operation_hours=7884, LA_distance=100, biosolids_price=0, FTE=0.6, lifetime=20):
    flowsheet_ID = 'C24'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(HeatDrying-0, 'diesel_LA'),
                                          outs=('biosolids_cost','methane_LA','nitrous_oxide_LA','carbon_dioxide_LA'),
                                          solids_distance=LA_distance)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    LandApplication.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    LandApplication.outs[0].price = biosolids_price
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    sys = qs.System.from_units(ID='system_C24',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    biosolids_trucking.price = (0.00551 + 0.0000541*LandApplication.solids_distance)/LandApplication.solids_distance
    
    biosolids_transportation = qs.Transportation(ID='Biosolids_transportation',
                                                 linked_unit=LandApplication,
                                                 item=biosolids_trucking,
                                                 load_type='mass',
                                                 load=stream.biosolids_cost.F_mass,
                                                 load_unit='kg',
                                                 distance=LandApplication.solids_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = kg/h
                                                 interval='1',
                                                 interval_unit='h')
    LandApplication.transportation.append(biosolids_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*lifetime,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*lifetime,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C25

def create_C25_system(size=10, operation_hours=7884, FTE=0.7, lifetime=20):
    flowsheet_ID = 'C25'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Incineration = lsu.Incineration(ID='Incineration',
                                    ins=(HeatDrying-0, 'natural_gas_incineration'),
                                    outs=('ash_incineration','vapor_incineration','methane_IN','nitrous_oxide_IN'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    Incineration.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    Incineration.outs[0].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    sys = qs.System.from_units(ID='system_C25',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Ash_incineration', linked_stream=stream.ash_incineration, GlobalWarming=0.018281422578429424)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_IN', linked_stream=stream.methane_IN, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_IN', linked_stream=stream.nitrous_oxide_IN, GlobalWarming=273)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.Incineration.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T1

def create_T1_system(size=10, operation_hours=7884, refinery_distance=100, FOAK=True, FTE=0.55, lifetime=20):
    flowsheet_ID = 'T1'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HTL = lsu.HydrothermalLiquefaction(ID='HTL', ins=Dewatering-0,
                                       outs=('biocrude','HTL_aqueous_undefined','hydrochar','offgas_HTL'),
                                       biocrude_distance=refinery_distance, FOAK=FOAK)
    # assume hydrochar from HTL is disposed of by sending it to landfills
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    HTL.outs[2].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    first_year_factor = HTL.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    Analyzer = lsu.Analyzer(ID='Analyzer',
                            ins=HTL-1,
                            outs='HTL_aqueous_defined')
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(Analyzer-0, 'virgin_CHG_catalyst'),
                                                outs=('CHG_out','used_CHG_catalyst'))
    # CHG catalyst price, https://doi.org/10.2172/1126336
    CHG.ins[1].price = 60/_lb_to_kg/GDPCTPI[2011]*GDPCTPI[2023]
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    F1.lifetime = 20
    
    GasMixer = qsu.Mixer(ID='GasMixer', ins=(HTL-3, F1-0), outs=('fuel_gas'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T1',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T1_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    # biocrude replacing crude oil of the same amount of energy
    # TODO: may need update the price to a more general number
    # 76.1 $/barrel crude oil, U.S. EIA, 2023 monthly average
    HTL.outs[0].price = 76.1/_oil_barrel_to_m3/HTL.crude_oil_density/HTL.crude_oil_HHV*HTL.biocrude_HHV
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # use market for petroleum to offset transportation and then add the transportation part
    # 0.6254770020033417 kg CO2 eq/kg petroleum ('market for petroleum')
    # assume biocrude is used to produce biofuel for combustion
    # pertoleum is 84% C, https://en.wikipedia.org/wiki/Petroleum
    qs.StreamImpactItem(ID='Biocrude', linked_stream=stream.biocrude, GlobalWarming=-(0.6254770020033417 + 0.84/12*44)/HTL.crude_oil_HHV*HTL.biocrude_HHV)
    qs.StreamImpactItem(ID='Hydrochar', linked_stream=stream.hydrochar, GlobalWarming=0.018281422578429424)
    qs.StreamImpactItem(ID='CHG_catalyst', linked_stream=stream.used_CHG_catalyst, GlobalWarming=1269.60621783954)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    biocrude_trucking = qs.ImpactItem('Biocrude_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biocrude_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # transportation cost: 5.67 2008$/m3 (fixed cost) and 0.07 2008$/m3/km (variable cost), https://doi.org/10.1016/j.biortech.2010.03.136
    biocrude_trucking.price = (5.67 + 0.07*HTL.biocrude_distance)/HTL.biocrude_density/HTL.biocrude_distance/GDPCTPI[2008]*GDPCTPI[2023]
    
    biocrude_transportation = qs.Transportation('Biocrude_transportation',
                                                linked_unit=HTL,
                                                item=biocrude_trucking,
                                                load_type='mass',
                                                load=stream.biocrude.F_mass,
                                                load_unit='kg',
                                                distance=HTL.biocrude_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    HTL.transportation.append(biocrude_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T2

def create_T2_system(size=10, operation_hours=7884, refinery_distance=100, LA_distance=100, FOAK=True, FTE=0.55, lifetime=20):
    flowsheet_ID = 'T2'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HALT = lsu.HydrothermalAlkalineTreatment(ID='HALT', ins=(Dewatering-0, 'sodium_hydroxide', 'hydrochloric_acid', 'diesel_HALT'),
                                             outs=('biocrude','HALT_aqueous_undefined','hydrochar','offgas_HALT'),
                                             biocrude_distance=refinery_distance, hydrochar_distance=LA_distance, FOAK=FOAK)
    # 0.2384 2016$/lb, https://doi.org/10.2172/1483234
    HALT.ins[1].price = 0.2384/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # 0.49 2020$/lb, from A.J.K. HALT model
    HALT.ins[2].price = 0.49/_lb_to_kg/GDPCTPI[2020]*GDPCTPI[2023]
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    HALT.ins[3].price = 4.224/_gal_to_liter*1000/diesel_density
    # 0 - 0.093 2018$, Figure 6 in https://www.sciencedirect.com/science/article/pii/S0306261919318021?via%3Dihub
    HALT.outs[2].price = 0.0465/GDPCTPI[2018]*GDPCTPI[2023]
    
    first_year_factor = HALT.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    Analyzer = lsu.Analyzer(ID='Analyzer',
                            ins=HALT-1,
                            outs='HALT_aqueous_defined')
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(Analyzer-0, 'virgin_CHG_catalyst'),
                                                outs=('CHG_out','used_CHG_catalyst'))
    # CHG catalyst price, https://doi.org/10.2172/1126336
    CHG.ins[1].price = 60/_lb_to_kg/GDPCTPI[2011]*GDPCTPI[2023]
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    F1.lifetime = 20
    
    GasMixer = qsu.Mixer(ID='GasMixer', ins=(HALT-3, F1-0), outs=('fuel_gas'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T2',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T2_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    # biocrude replacing crude oil of the same amount of energy
    # TODO: may need update the price to a more general number
    # 76.1 $/barrel crude oil, U.S. EIA, 2023 monthly average
    HALT.outs[0].price = 76.1/_oil_barrel_to_m3/HALT.crude_oil_density/HALT.crude_oil_HHV*HALT.biocrude_HHV
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Sodium_hydroxide', linked_stream=stream.sodium_hydroxide, GlobalWarming=1.4080631253939493)
    qs.StreamImpactItem(ID='Hydrochloric_acid', linked_stream=stream.hydrochloric_acid, GlobalWarming=0.8701969784483516)
    # use market for petroleum to offset transportation and then add the transportation part
    # 0.6254770020033417 kg CO2 eq/kg petroleum ('market for petroleum')
    # assume biocrude is used to produce biofuel for combustion
    # pertoleum is 84% C, https://en.wikipedia.org/wiki/Petroleum
    qs.StreamImpactItem(ID='Biocrude', linked_stream=stream.biocrude, GlobalWarming=-(0.6254770020033417 + 0.84/12*44)/HALT.crude_oil_HHV*HALT.biocrude_HHV)
    qs.StreamImpactItem(ID='CHG_catalyst', linked_stream=stream.used_CHG_catalyst, GlobalWarming=1269.60621783954)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_HALT', linked_stream=stream.diesel_HALT, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    # carbon sequestration
    # directly using hydrochar, unlike landfilling, land application, and composting which have carbon dioxide as fake streams
    # assume the carbon sequestration credit is the higher-end value of biosolids carbon sequestration (0.745 kg CO2 eq/kg dry solids, BEAM*2024)
    qs.StreamImpactItem(ID='Hydrochar', linked_stream=stream.hydrochar, GlobalWarming=-0.745)
    
    biocrude_trucking = qs.ImpactItem('Biocrude_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biocrude_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # transportation cost: 5.67 2008$/m3 (fixed cost) and 0.07 2008$/m3/km (variable cost), https://doi.org/10.1016/j.biortech.2010.03.136
    biocrude_trucking.price = (5.67 + 0.07*HALT.biocrude_distance)/HALT.biocrude_density/HALT.biocrude_distance/GDPCTPI[2008]*GDPCTPI[2023]
    
    biocrude_transportation = qs.Transportation('Biocrude_transportation',
                                                linked_unit=HALT,
                                                item=biocrude_trucking,
                                                load_type='mass',
                                                load=stream.biocrude.F_mass,
                                                load_unit='kg',
                                                distance=HALT.biocrude_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    HALT.transportation.append(biocrude_transportation)
    
    hydrochar_trucking = qs.ImpactItem('Hydrochar_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    hydrochar_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    hydrochar_trucking.price = (0.00551 + 0.0000541*HALT.hydrochar_distance)/HALT.hydrochar_distance
    
    hydrochar_transportation = qs.Transportation('Hydrochar_transportation',
                                                 linked_unit=HALT,
                                                 item=hydrochar_trucking,
                                                 load_type='mass',
                                                 load=stream.hydrochar.F_mass,
                                                 load_unit='kg',
                                                 distance=HALT.hydrochar_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = kg/h
                                                 interval='1',
                                                 interval_unit='h')
    HALT.transportation.append(hydrochar_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T3

def create_T3_system(size=10, operation_hours=7884, FOAK=True, FTE=0.4, lifetime=20):
    flowsheet_ID = 'T3'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    SCWO = lsu.SupercriticalWaterOxidation(ID='SCWO', ins=Dewatering-0, outs=('ash_SCWO','offgas_SCWO'), FOAK=FOAK)
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    SCWO.outs[0].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    first_year_factor = SCWO.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T3',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T3_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Ash_SCWO', linked_stream=stream.ash_SCWO, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T4

def create_T4_system(size=10, operation_hours=7884, refinery_distance=100, LA_distance=100, FOAK=True, FTE=0.7, lifetime=20):
    flowsheet_ID = 'T4'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Pyrolysis = lsu.Pyrolysis(ID='Pyrolysis', ins=(HeatDrying-0, 'diesel_pyrolysis'),
                              outs=('biooil','biochar','pyrogas'),
                              biooil_distance=refinery_distance, biochar_distance=LA_distance, FOAK=FOAK)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    Pyrolysis.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    # biooil replacing crude oil of the same amount of energy
    # TODO: may need update the price to a more general number
    # 76.1 $/barrel crude oil, U.S. EIA, 2023 monthly average
    Pyrolysis.outs[0].price = 76.1/_oil_barrel_to_m3/Pyrolysis.crude_oil_density/Pyrolysis.crude_oil_HHV*Pyrolysis.biooil_HHV
    # https://cloverly.com/blog/the-ultimate-business-guide-to-biochar-everything-you-need-to-know
    Pyrolysis.outs[1].price = 0.131
    
    first_year_factor = Pyrolysis.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(Pyrolysis-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T4',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T4_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # use market for petroleum to offset transportation and then add the transportation part
    # 0.6254770020033417 kg CO2 eq/kg petroleum ('market for petroleum')
    # assume biooil is used to produce biofuel for combustion
    # pertoleum is 84% C, https://en.wikipedia.org/wiki/Petroleum
    qs.StreamImpactItem(ID='Biooil', linked_stream=stream.biooil, GlobalWarming=-(0.6254770020033417 + 0.84/12*44)/Pyrolysis.crude_oil_HHV*Pyrolysis.biooil_HHV)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_pyrolysis', linked_stream=stream.diesel_pyrolysis, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    # carbon sequestration
    # the following number is consistent with the best guess from Table SI-12 in the SI of https://pubs.acs.org/doi/full/10.1021/acs.est.2c06083 (assuming 60% biochar is C)
    # directly using biochar, unlike landfilling, land application, and composting which have carbon dioxide as fake streams
    # carbon sequestration credit: 2 kg CO2 eq/kg dry solids, https://www.bioflux.earth/blog/what-is-biochar-carbon-removal
    qs.StreamImpactItem(ID='Biochar', linked_stream=stream.biochar, GlobalWarming=-2)
    
    biooil_trucking = qs.ImpactItem('Biooil_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biooil_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # transportation cost: 5.67 2008$/m3 (fixed cost) and 0.07 2008$/m3/km (variable cost), https://doi.org/10.1016/j.biortech.2010.03.136
    biooil_trucking.price = (5.67 + 0.07*Pyrolysis.biooil_distance)/Pyrolysis.biooil_density/Pyrolysis.biooil_distance/GDPCTPI[2008]*GDPCTPI[2023]
    
    biooil_transportation = qs.Transportation('Biooil_transportation',
                                              linked_unit=Pyrolysis,
                                              item=biooil_trucking,
                                              load_type='mass',
                                              load=stream.biooil.F_mass,
                                              load_unit='kg',
                                              distance=Pyrolysis.biooil_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Pyrolysis.transportation.append(biooil_transportation)
    
    biochar_trucking = qs.ImpactItem('Biochar_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biochar_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    biochar_trucking.price = (0.00551 + 0.0000541*Pyrolysis.biochar_distance)/Pyrolysis.biochar_distance
    
    bioochar_transportation = qs.Transportation('Bioochar_transportation',
                                                linked_unit=Pyrolysis,
                                                item=biochar_trucking,
                                                load_type='mass',
                                                load=stream.biochar.F_mass,
                                                load_unit='kg',
                                                distance=Pyrolysis.biochar_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    Pyrolysis.transportation.append(bioochar_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T5

def create_T5_system(size=10, operation_hours=7884, FOAK=True, FTE=0.7, lifetime=20):
    flowsheet_ID = 'T5'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(Thickening-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Gasification = lsu.Gasification(ID='Gasification', ins=HeatDrying-0, outs=('tar','ash_gasification','syngas'), FOAK=FOAK)
    # 300 to 510 $·tonne-1, including removal, transport, and disposal as a Resource Conservation and Recovery Act (RCRA) permitted facility
    # https://www.profitableventure.com/cost-dispose-hazardous-waste-per-ton/
    Gasification.outs[0].price = -0.405
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    Gasification.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    first_year_factor = Gasification.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(Gasification-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T5',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T5_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Tar', linked_stream=stream.tar, GlobalWarming=1.3479189232049629)
    qs.StreamImpactItem(ID='Ash_gasification', linked_stream=stream.ash_gasification, GlobalWarming=0.018281422578429424)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T6

def create_T6_system(size=10, operation_hours=7884, refinery_distance=100, FOAK=True, FTE=0.7, lifetime=20):
    flowsheet_ID = 'T6'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HTL = lsu.HydrothermalLiquefaction(ID='HTL', ins=Dewatering-0,
                                       outs=('biocrude','HTL_aqueous_undefined','hydrochar','offgas_HTL'),
                                       biocrude_distance=refinery_distance, FOAK=FOAK)
    # assume hydrochar from HTL is disposed of by sending it to landfills
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    HTL.outs[2].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    first_year_factor = HTL.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    Analyzer = lsu.Analyzer(ID='Analyzer',
                            ins=HTL-1,
                            outs='HTL_aqueous_defined')
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(Analyzer-0, 'virgin_CHG_catalyst'),
                                                outs=('CHG_out','used_CHG_catalyst'))
    # CHG catalyst price, https://doi.org/10.2172/1126336
    CHG.ins[1].price = 60/_lb_to_kg/GDPCTPI[2011]*GDPCTPI[2023]
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    
    GasMixer = qsu.Mixer(ID='GasMixer', ins=(HTL-3, F1-0), outs=('fuel_gas'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T6',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T6_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    # biocrude replacing crude oil of the same amount of energy
    # TODO: may need update the price to a more general number
    # 76.1 $/barrel crude oil, U.S. EIA, 2023 monthly average
    HTL.outs[0].price = 76.1/_oil_barrel_to_m3/HTL.crude_oil_density/HTL.crude_oil_HHV*HTL.biocrude_HHV
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # use market for petroleum to offset transportation and then add the transportation part
    # 0.6254770020033417 kg CO2 eq/kg petroleum ('market for petroleum')
    # assume biocrude is used to produce biofuel for combustion
    # pertoleum is 84% C, https://en.wikipedia.org/wiki/Petroleum
    qs.StreamImpactItem(ID='Biocrude', linked_stream=stream.biocrude, GlobalWarming=-(0.6254770020033417 + 0.84/12*44)/HTL.crude_oil_HHV*HTL.biocrude_HHV)
    qs.StreamImpactItem(ID='Hydrochar', linked_stream=stream.hydrochar, GlobalWarming=0.018281422578429424)
    qs.StreamImpactItem(ID='CHG_catalyst', linked_stream=stream.used_CHG_catalyst, GlobalWarming=1269.60621783954)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    biocrude_trucking = qs.ImpactItem('Biocrude_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biocrude_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # transportation cost: 5.67 2008$/m3 (fixed cost) and 0.07 2008$/m3/km (variable cost), https://doi.org/10.1016/j.biortech.2010.03.136
    biocrude_trucking.price = (5.67 + 0.07*HTL.biocrude_distance)/HTL.biocrude_density/HTL.biocrude_distance/GDPCTPI[2008]*GDPCTPI[2023]
    
    biocrude_transportation = qs.Transportation('Biocrude_transportation',
                                                linked_unit=HTL,
                                                item=biocrude_trucking,
                                                load_type='mass',
                                                load=stream.biocrude.F_mass,
                                                load_unit='kg',
                                                distance=HTL.biocrude_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    HTL.transportation.append(biocrude_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T7

def create_T7_system(size=10, operation_hours=7884, refinery_distance=100, LA_distance=100, FOAK=True, FTE=0.7, lifetime=20):
    flowsheet_ID = 'T7'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HALT = lsu.HydrothermalAlkalineTreatment(ID='HALT', ins=(Dewatering-0, 'sodium_hydroxide', 'hydrochloric_acid', 'diesel_HALT'),
                                             outs=('biocrude','HALT_aqueous_undefined','hydrochar','offgas_HALT'),
                                             biocrude_distance=refinery_distance, hydrochar_distance=LA_distance, FOAK=FOAK)
    # 0.2384 2016$/lb, https://doi.org/10.2172/1483234
    HALT.ins[1].price = 0.2384/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # 0.49 2020$/lb, from A.J.K. HALT model
    HALT.ins[2].price = 0.49/_lb_to_kg/GDPCTPI[2020]*GDPCTPI[2023]
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    HALT.ins[3].price = 4.224/_gal_to_liter*1000/diesel_density
    # 0 - 0.093 2018$, Figure 6 in https://www.sciencedirect.com/science/article/pii/S0306261919318021?via%3Dihub
    HALT.outs[2].price = 0.0465/GDPCTPI[2018]*GDPCTPI[2023]
    
    first_year_factor = HALT.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    Analyzer = lsu.Analyzer(ID='Analyzer',
                            ins=HALT-1,
                            outs='HALT_aqueous_defined')
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(Analyzer-0, 'virgin_CHG_catalyst'),
                                                outs=('CHG_out','used_CHG_catalyst'))
    # CHG catalyst price, https://doi.org/10.2172/1126336
    CHG.ins[1].price = 60/_lb_to_kg/GDPCTPI[2011]*GDPCTPI[2023]
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    
    GasMixer = qsu.Mixer(ID='GasMixer', ins=(HALT-3, F1-0), outs=('fuel_gas'))
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T7',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T7_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    # biocrude replacing crude oil of the same amount of energy
    # TODO: may need update the price to a more general number
    # 76.1 $/barrel crude oil, U.S. EIA, 2023 monthly average
    HALT.outs[0].price = 76.1/_oil_barrel_to_m3/HALT.crude_oil_density/HALT.crude_oil_HHV*HALT.biocrude_HHV
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Sodium_hydroxide', linked_stream=stream.sodium_hydroxide, GlobalWarming=1.4080631253939493)
    qs.StreamImpactItem(ID='Hydrochloric_acid', linked_stream=stream.hydrochloric_acid, GlobalWarming=0.8701969784483516)
    # use market for petroleum to offset transportation and then add the transportation part
    # 0.6254770020033417 kg CO2 eq/kg petroleum ('market for petroleum')
    # assume biocrude is used to produce biofuel for combustion
    # pertoleum is 84% C, https://en.wikipedia.org/wiki/Petroleum
    qs.StreamImpactItem(ID='Biocrude', linked_stream=stream.biocrude, GlobalWarming=-(0.6254770020033417 + 0.84/12*44)/HALT.crude_oil_HHV*HALT.biocrude_HHV)
    qs.StreamImpactItem(ID='CHG_catalyst', linked_stream=stream.used_CHG_catalyst, GlobalWarming=1269.60621783954)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_HALT', linked_stream=stream.diesel_HALT, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    # carbon sequestration
    # directly using hydrochar, unlike landfilling, land application, and composting which have carbon dioxide as fake streams
    # assume the carbon sequestration credit is the higher-end value of biosolids carbon sequestration (0.745 kg CO2 eq/kg dry solids, BEAM*2024)
    qs.StreamImpactItem(ID='Hydrochar', linked_stream=stream.hydrochar, GlobalWarming=-0.745)
    
    biocrude_trucking = qs.ImpactItem('Biocrude_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biocrude_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # transportation cost: 5.67 2008$/m3 (fixed cost) and 0.07 2008$/m3/km (variable cost), https://doi.org/10.1016/j.biortech.2010.03.136
    biocrude_trucking.price = (5.67 + 0.07*HALT.biocrude_distance)/HALT.biocrude_density/HALT.biocrude_distance/GDPCTPI[2008]*GDPCTPI[2023]
    
    biocrude_transportation = qs.Transportation('Biocrude_transportation',
                                                linked_unit=HALT,
                                                item=biocrude_trucking,
                                                load_type='mass',
                                                load=stream.biocrude.F_mass,
                                                load_unit='kg',
                                                distance=HALT.biocrude_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    HALT.transportation.append(biocrude_transportation)
    
    hydrochar_trucking = qs.ImpactItem('Hydrochar_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    hydrochar_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    hydrochar_trucking.price = (0.00551 + 0.0000541*HALT.hydrochar_distance)/HALT.hydrochar_distance
    
    hydrochar_transportation = qs.Transportation('Hydrochar_transportation',
                                                 linked_unit=HALT,
                                                 item=hydrochar_trucking,
                                                 load_type='mass',
                                                 load=stream.hydrochar.F_mass,
                                                 load_unit='kg',
                                                 distance=HALT.hydrochar_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = kg/h
                                                 interval='1',
                                                 interval_unit='h')
    HALT.transportation.append(hydrochar_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T8

def create_T8_system(size=10, operation_hours=7884, FOAK=True, FTE=0.55, lifetime=20):
    flowsheet_ID = 'T8'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    SCWO = lsu.SupercriticalWaterOxidation(ID='SCWO', ins=Dewatering-0, outs=('ash_SCWO','offgas_SCWO'), FOAK=FOAK)
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    SCWO.outs[0].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    first_year_factor = SCWO.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T8',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T8_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Ash_SCWO', linked_stream=stream.ash_SCWO, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T9

def create_T9_system(size=10, operation_hours=7884, refinery_distance=100, LA_distance=100, FOAK=True, FTE=0.85, lifetime=20):
    flowsheet_ID = 'T9'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Pyrolysis = lsu.Pyrolysis(ID='Pyrolysis', ins=(HeatDrying-0, 'diesel_pyrolysis'),
                              outs=('biooil','biochar','pyrogas'),
                              biooil_distance=refinery_distance, biochar_distance=LA_distance, FOAK=FOAK)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    Pyrolysis.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    # biooil replacing crude oil of the same amount of energy
    # TODO: may need update the price to a more general number
    # 76.1 $/barrel crude oil, U.S. EIA, 2023 monthly average
    Pyrolysis.outs[0].price = 76.1/_oil_barrel_to_m3/Pyrolysis.crude_oil_density/Pyrolysis.crude_oil_HHV*Pyrolysis.biooil_HHV
    # https://cloverly.com/blog/the-ultimate-business-guide-to-biochar-everything-you-need-to-know
    Pyrolysis.outs[1].price = 0.131
    
    first_year_factor = Pyrolysis.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(Pyrolysis-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T9',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T9_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # use market for petroleum to offset transportation and then add the transportation part
    # 0.6254770020033417 kg CO2 eq/kg petroleum ('market for petroleum')
    # assume biooil is used to produce biofuel for combustion
    # pertoleum is 84% C, https://en.wikipedia.org/wiki/Petroleum
    qs.StreamImpactItem(ID='Biooil', linked_stream=stream.biooil, GlobalWarming=-(0.6254770020033417 + 0.84/12*44)/Pyrolysis.crude_oil_HHV*Pyrolysis.biooil_HHV)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_pyrolysis', linked_stream=stream.diesel_pyrolysis, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    # carbon sequestration
    # the following number is consistent with the best guess from Table SI-12 in the SI of https://pubs.acs.org/doi/full/10.1021/acs.est.2c06083 (assuming 60% biochar is C)
    # directly using biochar, unlike landfilling, land application, and composting which have carbon dioxide as fake streams
    # carbon sequestration credit: 2 kg CO2 eq/kg dry solids, https://www.bioflux.earth/blog/what-is-biochar-carbon-removal
    qs.StreamImpactItem(ID='Biochar', linked_stream=stream.biochar, GlobalWarming=-2)
    
    biooil_trucking = qs.ImpactItem('Biooil_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biooil_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # transportation cost: 5.67 2008$/m3 (fixed cost) and 0.07 2008$/m3/km (variable cost), https://doi.org/10.1016/j.biortech.2010.03.136
    biooil_trucking.price = (5.67 + 0.07*Pyrolysis.biooil_distance)/Pyrolysis.biooil_density/Pyrolysis.biooil_distance/GDPCTPI[2008]*GDPCTPI[2023]
    
    biooil_transportation = qs.Transportation('Biooil_transportation',
                                              linked_unit=Pyrolysis,
                                              item=biooil_trucking,
                                              load_type='mass',
                                              load=stream.biooil.F_mass,
                                              load_unit='kg',
                                              distance=Pyrolysis.biooil_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Pyrolysis.transportation.append(biooil_transportation)
    
    biochar_trucking = qs.ImpactItem('Biochar_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biochar_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    biochar_trucking.price = (0.00551 + 0.0000541*Pyrolysis.biochar_distance)/Pyrolysis.biochar_distance
    
    bioochar_transportation = qs.Transportation('Bioochar_transportation',
                                                linked_unit=Pyrolysis,
                                                item=biochar_trucking,
                                                load_type='mass',
                                                load=stream.biochar.F_mass,
                                                load_unit='kg',
                                                distance=Pyrolysis.biochar_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    Pyrolysis.transportation.append(bioochar_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T10

def create_T10_system(size=10, operation_hours=7884, FOAK=True, FTE=0.85, lifetime=20):
    flowsheet_ID = 'T10'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AerobicDigestion = lsu.AerobicDigestion(ID='AerobicDigestion', ins=(Thickening-0, 'air'),
                                            outs=('digested_sludge','offgas_AeD'))
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Gasification = lsu.Gasification(ID='Gasification', ins=HeatDrying-0, outs=('tar','ash_gasification','syngas'), FOAK=FOAK)
    # 300 to 510 $·tonne-1, including removal, transport, and disposal as a Resource Conservation and Recovery Act (RCRA) permitted facility
    # https://www.profitableventure.com/cost-dispose-hazardous-waste-per-ton/
    Gasification.outs[0].price = -0.405
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    Gasification.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    first_year_factor = Gasification.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(Gasification-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T10',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T10_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Tar', linked_stream=stream.tar, GlobalWarming=1.3479189232049629)
    qs.StreamImpactItem(ID='Ash_gasification', linked_stream=stream.ash_gasification, GlobalWarming=0.018281422578429424)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T11

def create_T11_system(size=10, operation_hours=7884, refinery_distance=100, FOAK=True, FTE=0.7, lifetime=20):
    flowsheet_ID = 'T11'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HTL = lsu.HydrothermalLiquefaction(ID='HTL', ins=Dewatering-0,
                                       outs=('biocrude','HTL_aqueous_undefined','hydrochar','offgas_HTL'),
                                       biocrude_distance=refinery_distance, FOAK=FOAK)
    # assume hydrochar from HTL is disposed of by sending it to landfills
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    HTL.outs[2].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    first_year_factor = HTL.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    Analyzer = lsu.Analyzer(ID='Analyzer',
                            ins=HTL-1,
                            outs='HTL_aqueous_defined')
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(Analyzer-0, 'virgin_CHG_catalyst'),
                                                outs=('CHG_out','used_CHG_catalyst'))
    # CHG catalyst price, https://doi.org/10.2172/1126336
    CHG.ins[1].price = 60/_lb_to_kg/GDPCTPI[2011]*GDPCTPI[2023]
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    
    GasMixer = qsu.Mixer('GasMixer', ins=(AnaerobicDigestion-1, HTL-3, F1-0),
                         outs=('fuel_gas'), init_with='WasteStream')
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T11',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T11_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    # biocrude replacing crude oil of the same amount of energy
    # TODO: may need update the price to a more general number
    # 76.1 $/barrel crude oil, U.S. EIA, 2023 monthly average, may need to update this to a more general number
    HTL.outs[0].price = 76.1/_oil_barrel_to_m3/HTL.crude_oil_density/HTL.crude_oil_HHV*HTL.biocrude_HHV
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # use market for petroleum to offset transportation and then add the transportation part
    # 0.6254770020033417 kg CO2 eq/kg petroleum ('market for petroleum')
    # assume biocrude is used to produce biofuel for combustion
    # pertoleum is 84% C, https://en.wikipedia.org/wiki/Petroleum
    qs.StreamImpactItem(ID='Biocrude', linked_stream=stream.biocrude, GlobalWarming=-(0.6254770020033417 + 0.84/12*44)/HTL.crude_oil_HHV*HTL.biocrude_HHV)
    qs.StreamImpactItem(ID='Hydrochar', linked_stream=stream.hydrochar, GlobalWarming=0.018281422578429424)
    qs.StreamImpactItem(ID='CHG_catalyst', linked_stream=stream.used_CHG_catalyst, GlobalWarming=1269.60621783954)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    biocrude_trucking = qs.ImpactItem('Biocrude_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biocrude_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # transportation cost: 5.67 2008$/m3 (fixed cost) and 0.07 2008$/m3/km (variable cost), https://doi.org/10.1016/j.biortech.2010.03.136
    biocrude_trucking.price = (5.67 + 0.07*HTL.biocrude_distance)/HTL.biocrude_density/HTL.biocrude_distance/GDPCTPI[2008]*GDPCTPI[2023]
    
    biocrude_transportation = qs.Transportation('Biocrude_transportation',
                                                linked_unit=HTL,
                                                item=biocrude_trucking,
                                                load_type='mass',
                                                load=stream.biocrude.F_mass,
                                                load_unit='kg',
                                                distance=HTL.biocrude_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    HTL.transportation.append(biocrude_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T12

def create_T12_system(size=10, operation_hours=7884, refinery_distance=100, LA_distance=100, FOAK=True, FTE=0.7, lifetime=20):
    flowsheet_ID = 'T12'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HALT = lsu.HydrothermalAlkalineTreatment(ID='HALT', ins=(Dewatering-0, 'sodium_hydroxide', 'hydrochloric_acid', 'diesel_HALT'),
                                             outs=('biocrude','HALT_aqueous_undefined','hydrochar','offgas_HALT'),
                                             biocrude_distance=refinery_distance, hydrochar_distance=LA_distance, FOAK=FOAK)
    # 0.2384 2016$/lb, https://doi.org/10.2172/1483234
    HALT.ins[1].price = 0.2384/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # 0.49 2020$/lb, from A.J.K. HALT model
    HALT.ins[2].price = 0.49/_lb_to_kg/GDPCTPI[2020]*GDPCTPI[2023]
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    HALT.ins[3].price = 4.224/_gal_to_liter*1000/diesel_density
    # 0 - 0.093 2018$, Figure 6 in https://www.sciencedirect.com/science/article/pii/S0306261919318021?via%3Dihub
    HALT.outs[2].price = 0.0465/GDPCTPI[2018]*GDPCTPI[2023]
    
    first_year_factor = HALT.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    Analyzer = lsu.Analyzer(ID='Analyzer',
                            ins=HALT-1,
                            outs='HALT_aqueous_defined')
    
    CHG = lsu.CatalyticHydrothermalGasification(ID='CHG', ins=(Analyzer-0, 'virgin_CHG_catalyst'),
                                                outs=('CHG_out','used_CHG_catalyst'))
    # CHG catalyst price, https://doi.org/10.2172/1126336
    CHG.ins[1].price = 60/_lb_to_kg/GDPCTPI[2011]*GDPCTPI[2023]
    
    V1 = IsenthalpicValve(ID='V1', ins=CHG-0, outs='depressurized_cooled_CHG', P=50*6894.76, vle=True)
    
    F1 = qsu.Flash(ID='F1', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    
    GasMixer = qsu.Mixer('GasMixer', ins=(AnaerobicDigestion-1, HALT-3, F1-0),
                         outs=('fuel_gas'), init_with='WasteStream')
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T12',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T12_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    # biocrude replacing crude oil of the same amount of energy
    # TODO: may need update the price to a more general number
    # 76.1 $/barrel crude oil, U.S. EIA, 2023 monthly average
    HALT.outs[0].price = 76.1/_oil_barrel_to_m3/HALT.crude_oil_density/HALT.crude_oil_HHV*HALT.biocrude_HHV
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Sodium_hydroxide', linked_stream=stream.sodium_hydroxide, GlobalWarming=1.4080631253939493)
    qs.StreamImpactItem(ID='Hydrochloric_acid', linked_stream=stream.hydrochloric_acid, GlobalWarming=0.8701969784483516)
    # use market for petroleum to offset transportation and then add the transportation part
    # 0.6254770020033417 kg CO2 eq/kg petroleum ('market for petroleum')
    # assume biocrude is used to produce biofuel for combustion
    # pertoleum is 84% C, https://en.wikipedia.org/wiki/Petroleum
    qs.StreamImpactItem(ID='Biocrude', linked_stream=stream.biocrude, GlobalWarming=-(0.6254770020033417 + 0.84/12*44)/HALT.crude_oil_HHV*HALT.biocrude_HHV)
    qs.StreamImpactItem(ID='CHG_catalyst', linked_stream=stream.used_CHG_catalyst, GlobalWarming=1269.60621783954)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_HALT', linked_stream=stream.diesel_HALT, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    # carbon sequestration
    # directly using hydrochar, unlike landfilling, land application, and composting which have carbon dioxide as fake streams
    # assume the carbon sequestration credit is the higher-end value of biosolids carbon sequestration (0.745 kg CO2 eq/kg dry solids, BEAM*2024)
    qs.StreamImpactItem(ID='Hydrochar', linked_stream=stream.hydrochar, GlobalWarming=-0.745)
    
    biocrude_trucking = qs.ImpactItem('Biocrude_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biocrude_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # transportation cost: 5.67 2008$/m3 (fixed cost) and 0.07 2008$/m3/km (variable cost), https://doi.org/10.1016/j.biortech.2010.03.136
    biocrude_trucking.price = (5.67 + 0.07*HALT.biocrude_distance)/HALT.biocrude_density/HALT.biocrude_distance/GDPCTPI[2008]*GDPCTPI[2023]
    
    biocrude_transportation = qs.Transportation('Biocrude_transportation',
                                                linked_unit=HALT,
                                                item=biocrude_trucking,
                                                load_type='mass',
                                                load=stream.biocrude.F_mass,
                                                load_unit='kg',
                                                distance=HALT.biocrude_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    HALT.transportation.append(biocrude_transportation)
    
    hydrochar_trucking = qs.ImpactItem('Hydrochar_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    hydrochar_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    hydrochar_trucking.price = (0.00551 + 0.0000541*HALT.hydrochar_distance)/HALT.hydrochar_distance
    
    hydrochar_transportation = qs.Transportation('Hydrochar_transportation',
                                                 linked_unit=HALT,
                                                 item=hydrochar_trucking,
                                                 load_type='mass',
                                                 load=stream.hydrochar.F_mass,
                                                 load_unit='kg',
                                                 distance=HALT.hydrochar_distance,
                                                 distance_unit='km',
                                                 # set to 1 h since load = kg/h
                                                 interval='1',
                                                 interval_unit='h')
    HALT.transportation.append(hydrochar_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T13

def create_T13_system(size=10, operation_hours=7884, FOAK=True, FTE=0.55, lifetime=20):
    flowsheet_ID = 'T13'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    SCWO = lsu.SupercriticalWaterOxidation(ID='SCWO', ins=Dewatering-0, outs=('ash_SCWO','offgas_SCWO'), FOAK=FOAK)
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    SCWO.outs[0].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    first_year_factor = SCWO.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T13',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T13_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Ash_SCWO', linked_stream=stream.ash_SCWO, GlobalWarming=0.018281422578429424)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T14

def create_T14_system(size=10, operation_hours=7884, refinery_distance=100, LA_distance=100, FOAK=True, FTE=0.85, lifetime=20):
    flowsheet_ID = 'T14'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Pyrolysis = lsu.Pyrolysis(ID='Pyrolysis', ins=(HeatDrying-0, 'diesel_pyrolysis'),
                              outs=('biooil','biochar','pyrogas'),
                              biooil_distance=refinery_distance, biochar_distance=LA_distance, FOAK=FOAK)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    Pyrolysis.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    # biooil replacing crude oil of the same amount of energy
    # TODO: may need update the price to a more general number
    # 76.1 $/barrel crude oil, U.S. EIA, 2023 monthly average
    Pyrolysis.outs[0].price = 76.1/_oil_barrel_to_m3/Pyrolysis.crude_oil_density/Pyrolysis.crude_oil_HHV*Pyrolysis.biooil_HHV
    # https://cloverly.com/blog/the-ultimate-business-guide-to-biochar-everything-you-need-to-know
    Pyrolysis.outs[1].price = 0.131
    
    first_year_factor = Pyrolysis.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    GasMixer = qsu.Mixer('GasMixer', ins=(AnaerobicDigestion-1, Pyrolysis-2),
                         outs=('fuel_gas'), init_with='WasteStream')
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T14',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T14_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    # use market for petroleum to offset transportation and then add the transportation part
    # 0.6254770020033417 kg CO2 eq/kg petroleum ('market for petroleum')
    # assume biooil is used to produce biofuel for combustion
    # pertoleum is 84% C, https://en.wikipedia.org/wiki/Petroleum
    qs.StreamImpactItem(ID='Biooil', linked_stream=stream.biooil, GlobalWarming=-(0.6254770020033417 + 0.84/12*44)/Pyrolysis.crude_oil_HHV*Pyrolysis.biooil_HHV)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_pyrolysis', linked_stream=stream.diesel_pyrolysis, GlobalWarming=0.8608804649420178 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    # carbon sequestration
    # the following number is consistent with the best guess from Table SI-12 in the SI of https://pubs.acs.org/doi/full/10.1021/acs.est.2c06083 (assuming 60% biochar is C)
    # directly using biochar, unlike landfilling, land application, and composting which have carbon dioxide as fake streams
    # carbon sequestration credit: 2 kg CO2 eq/kg dry solids, https://www.bioflux.earth/blog/what-is-biochar-carbon-removal
    qs.StreamImpactItem(ID='Biochar', linked_stream=stream.biochar, GlobalWarming=-2)
    
    biooil_trucking = qs.ImpactItem('Biooil_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biooil_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # transportation cost: 5.67 2008$/m3 (fixed cost) and 0.07 2008$/m3/km (variable cost), https://doi.org/10.1016/j.biortech.2010.03.136
    biooil_trucking.price = (5.67 + 0.07*Pyrolysis.biooil_distance)/Pyrolysis.biooil_density/Pyrolysis.biooil_distance/GDPCTPI[2008]*GDPCTPI[2023]
    
    biooil_transportation = qs.Transportation('Biooil_transportation',
                                              linked_unit=Pyrolysis,
                                              item=biooil_trucking,
                                              load_type='mass',
                                              load=stream.biooil.F_mass,
                                              load_unit='kg',
                                              distance=Pyrolysis.biooil_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    Pyrolysis.transportation.append(biooil_transportation)
    
    biochar_trucking = qs.ImpactItem('Biochar_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biochar_trucking.add_indicator(GlobalWarming, 0.1611858456466717/1000)
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    biochar_trucking.price = (0.00551 + 0.0000541*Pyrolysis.biochar_distance)/Pyrolysis.biochar_distance
    
    bioochar_transportation = qs.Transportation('Bioochar_transportation',
                                                linked_unit=Pyrolysis,
                                                item=biochar_trucking,
                                                load_type='mass',
                                                load=stream.biochar.F_mass,
                                                load_unit='kg',
                                                distance=Pyrolysis.biochar_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    Pyrolysis.transportation.append(bioochar_transportation)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T15

def create_T15_system(size=10, operation_hours=7884, FOAK=True, FTE=0.85, lifetime=20):
    flowsheet_ID = 'T15'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        
    clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]
    
    bst.PowerUtility.price = 0.16
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)

    WRRF = lsu.WRRF(ID='WRRF',
                    ins=raw_wastewater,
                    outs=('sludge','treated_water'),
                    ww_2_dry_sludge=1,
                    sludge_moisture=0.99,
                    sludge_dw_ash=0.436,
                    sludge_afdw_lipid=0.193,
                    sludge_afdw_protein=0.510,
                    sludge_wet_density=1040)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Thickening = lsu.Thickening(ID='Thickening', ins=(WRRF-0, 'polymer_thickening'),
                                outs=('thickened_sludge','reject_thickening'))
    Thickening.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD','carbon_dioxide_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    HeatDrying = lsu.HeatDrying(ID='HeatDrying', ins=(Dewatering-0, 'natural_gas_heat_drying'),
                                outs=('dried_solids','vapor'))
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    HeatDrying.ins[1].price = 0.218
    
    Gasification = lsu.Gasification(ID='Gasification', ins=HeatDrying-0, outs=('tar','ash_gasification','syngas'), FOAK=FOAK)
    # 300 to 510 $·tonne-1, including removal, transport, and disposal as a Resource Conservation and Recovery Act (RCRA) permitted facility
    # https://www.profitableventure.com/cost-dispose-hazardous-waste-per-ton/
    Gasification.outs[0].price = -0.405
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    Gasification.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    first_year_factor = Gasification.plant_performance_factor/100
    
    performance_factor_list = [min(first_year_factor*(1.2)**i, 1) for i in range(0, lifetime)]
    
    operation_hours *= np.mean(performance_factor_list)
    
    GasMixer = qsu.Mixer('GasMixer', ins=(AnaerobicDigestion-1, Gasification-2),
                         outs=('fuel_gas'), init_with='WasteStream')
    
    CHP = qsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # treatment of CT residual is not considerd (since it is not in the model)
    # but ecoinvent has it and indicates it might be disposed of through landfilling
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_makeup_water
    CT.ins[1].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # cooling_tower_chemicals: 1.7842 2016$/lb, https://doi.org/10.2172/1483234
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    CT.lifetime = 20
    
    sys = qs.System.from_units(ID='system_T15',
                               units=list(flowsheet.unit),
                               operating_hours=operation_hours)
    
    try:
        sys.simulate()
    except AttributeError as e:
        if 'CoolingTower: CT' in str(e) and "'NoneType' object has no attribute 'T'" in str(e):
            flowsheet.remove_unit_and_associated_streams(ID='CT')
            
            sys = qs.System.from_units(ID='system_T15_no_CT',
                                       units=list(flowsheet.unit),
                                       operating_hours=operation_hours)
            
            sys.simulate()
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.691007559959689)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.12677990083093105)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.03882149971451173)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.5780189368532676 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.0680678777230173)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00047332694442793645)
    
    def deionized_water_quantity():
        try:
            # the cooling_tower_chemicals is a mix of chemicals, therefore, the following approach is used to approximate the CI of the cooling_tower_chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*lifetime
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.670217713628895)
    qs.StreamImpactItem(ID='Tar', linked_stream=stream.tar, GlobalWarming=1.3479189232049629)
    qs.StreamImpactItem(ID='Ash_gasification', linked_stream=stream.ash_gasification, GlobalWarming=0.018281422578429424)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.018281422578429424)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    qs.LCA(system=sys, lifetime=lifetime, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys