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


#%% initialization

import qsdsan as qs, biosteam as bst
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

# TODO: add in writing: use 2023$ (since CEPCI is not available after that, need a fee to access the full 2024 data)

# GDPCTPI (Gross Domestic Product: Chain-type Price Index)
# https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20)
GDPCTPI = {2005: 81.537,
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

def create_C1_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.15):
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Landfilling.transportation = sludge_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C2

def create_C2_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.2):
    flowsheet_ID = 'C2'
    
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
                                                      ins=(Dewatering-0, 'lime', 'natural_gas_lime'),
                                                      outs='stabilized_solids')
    AlkalineStabilization.ins[1].price = 0.1189/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    AlkalineStabilization.ins[2].price = 0.218
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=AlkalineStabilization-0,
                                  outs=('landfilled_solids','methane_LF','nitrous_oxide_LF','sequestered_carbon_dioxide_LF'),
                                  solids_distance=LF_distance)
    
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Lime', linked_stream=stream.lime, GlobalWarming=0.04261602)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Landfilling.transportation = sludge_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:sys.flowsheet.AlkalineStabilization.ins[2].F_vol*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C3

def create_C3_system(size=10, operation_hours=8760, LA_distance=100, FTE=0.2):
    flowsheet_ID = 'C3'
    
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
                                                      ins=(Dewatering-0, 'lime', 'natural_gas_lime'),
                                                      outs='stabilized_solids')
    AlkalineStabilization.ins[1].price = 0.1189/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    AlkalineStabilization.ins[2].price = 0.218
    
    LandApplication = lsu.LandApplication(ID='LandApplication',
                                          ins=(AlkalineStabilization-0, 'diesel_LA'),
                                          outs=('biosolids_cost','methane_LA','nitrous_oxide_LA','carbon_dioxide_LA'),
                                          solids_distance=LA_distance)
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    LandApplication.ins[1].price = 4.224/_gal_to_liter*1000/diesel_density
    # TODO: update cost for biosolids_cost with references (should be based on dry weight since heat drying in some systems can reduce moisture)
    LandApplication.outs[0].price = 10/1000
    
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Lime', linked_stream=stream.lime, GlobalWarming=0.04261602)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    LandApplication.transportation = biosolids_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:sys.flowsheet.AlkalineStabilization.ins[2].F_vol*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C4

def create_C4_system(size=10, operation_hours=8760, LA_distance=100, FTE=0.2):
    flowsheet_ID = 'C4'
    
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
    
    # TODO: make sure the Composting unit is composting + land application
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel_composting'),
                                outs=('compost_cost','methane_composting','nitrous_oxide_composting','sequestered_carbon_dioxide_composting'),
                                solids_distance=LA_distance)
    # TODO: need to decide how to deal with the bulking_agent in composting and its costing
    # TODO: uncertainty range (uniform) 18-36 2005$/tonne
    Composting.ins[1].price = 27/1000/GDPCTPI[2005]*GDPCTPI[2023]
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    Composting.ins[2].price = 4.224/_gal_to_liter*1000/diesel_density
    # TODO: update cost for compost_cost with references
    Composting.outs[0].price = 10/1000
    
    sys = qs.System.from_units(ID='system_C4',
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Bulking_agent', linked_stream=stream.bulking_agent, GlobalWarming=0.041056332)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_composting', linked_stream=stream.diesel_composting, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_composting', linked_stream=stream.methane_composting, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_composting', linked_stream=stream.nitrous_oxide_composting, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_composting', linked_stream=stream.sequestered_carbon_dioxide_composting, GlobalWarming=1)
    
    compost_trucking = qs.ImpactItem(ID='Compost_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    compost_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Composting.transportation = compost_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.Composting.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.Composting.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C5

def create_C5_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.3):
    flowsheet_ID = 'C5'
    
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Landfilling.transportation = sludge_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:sys.flowsheet.HeatDrying.ins[1].F_vol*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C6

def create_C6_system(size=10, operation_hours=8760, LA_distance=100, FTE=0.3):
    flowsheet_ID = 'C6'
    
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    LandApplication.transportation = biosolids_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:sys.flowsheet.HeatDrying.ins[1].F_vol*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C7

def create_C7_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.4):
    flowsheet_ID = 'C7'
    
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
                                    outs=('ash_incineration','vapor_incineration','methane_IN','nitrous_oxide_IN'),
                                    solids_distance=LF_distance)
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_incineration', linked_stream=stream.ash_incineration, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_IN', linked_stream=stream.methane_IN, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_IN', linked_stream=stream.nitrous_oxide_IN, GlobalWarming=273)
    
    ash_incineration_trucking = qs.ImpactItem(ID='Ash_incineration_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_incineration_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_incineration_trucking.price = (0.00551 + 0.0000541*Incineration.solids_distance)/Incineration.solids_distance
    
    ash_incineration_transportation = qs.Transportation(ID='Ash_incineration_transportation',
                                                        linked_unit=Incineration,
                                                        item=ash_incineration_trucking,
                                                        load_type='mass',
                                                        load=stream.ash_incineration.F_mass,
                                                        load_unit='kg',
                                                        distance=Incineration.solids_distance,
                                                        distance_unit='km',
                                                        # set to 1 h since load = kg/h
                                                        interval='1',
                                                        interval_unit='h')
    Incineration.transportation = ash_incineration_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.Incineration.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C8

def create_C8_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.3):
    flowsheet_ID = 'C8'
    
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Landfilling.transportation = sludge_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C9

def create_C9_system(size=10, operation_hours=8760, LA_distance=100, FTE=0.3):
    flowsheet_ID = 'C9'
    
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    LandApplication.transportation = biosolids_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C10

def create_C10_system(size=10, operation_hours=8760, LA_distance=100, FTE=0.35):
    flowsheet_ID = 'C10'
    
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
    
    # TODO: make sure the Composting unit is composting + land application
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel_composting'),
                                outs=('compost_cost','methane_composting','nitrous_oxide_composting','sequestered_carbon_dioxide_composting'),
                                solids_distance=LA_distance)
    # TODO: need to decide how to deal with the bulking_agent in composting and its costing
    # TODO: uncertainty range (uniform) 18-36 2005$/tonne
    Composting.ins[1].price = 27/1000/GDPCTPI[2005]*GDPCTPI[2023]
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    Composting.ins[2].price = 4.224/_gal_to_liter*1000/diesel_density
    # TODO: update cost for compost_cost with references
    Composting.outs[0].price = 10/1000
    
    sys = qs.System.from_units(ID='system_C10',
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Bulking_agent', linked_stream=stream.bulking_agent, GlobalWarming=0.041056332)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_composting', linked_stream=stream.diesel_composting, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_composting', linked_stream=stream.methane_composting, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_composting', linked_stream=stream.nitrous_oxide_composting, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_composting', linked_stream=stream.sequestered_carbon_dioxide_composting, GlobalWarming=1)
    
    compost_trucking = qs.ImpactItem(ID='Compost_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    compost_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Composting.transportation = compost_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.Composting.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.Composting.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C11

def create_C11_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.45):
    flowsheet_ID = 'C11'
    
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Landfilling.transportation = sludge_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:sys.flowsheet.HeatDrying.ins[1].F_vol*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C12

def create_C12_system(size=10, operation_hours=8760, LA_distance=100, FTE=0.45):
    flowsheet_ID = 'C12'
    
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    LandApplication.transportation = biosolids_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:sys.flowsheet.HeatDrying.ins[1].F_vol*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C13

def create_C13_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.55):
    flowsheet_ID = 'C13'
    
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
                                    outs=('ash_incineration','vapor_incineration','methane_IN','nitrous_oxide_IN'),
                                    solids_distance=LF_distance)
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_incineration', linked_stream=stream.ash_incineration, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_IN', linked_stream=stream.methane_IN, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_IN', linked_stream=stream.nitrous_oxide_IN, GlobalWarming=273)
    
    ash_incineration_trucking = qs.ImpactItem(ID='Ash_incineration_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_incineration_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_incineration_trucking.price = (0.00551 + 0.0000541*Incineration.solids_distance)/Incineration.solids_distance
    
    ash_incineration_transportation = qs.Transportation(ID='Ash_incineration_transportation',
                                                        linked_unit=Incineration,
                                                        item=ash_incineration_trucking,
                                                        load_type='mass',
                                                        load=stream.ash_incineration.F_mass,
                                                        load_unit='kg',
                                                        distance=Incineration.solids_distance,
                                                        distance_unit='km',
                                                        # set to 1 h since load = kg/h
                                                        interval='1',
                                                        interval_unit='h')
    Incineration.transportation = ash_incineration_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.Incineration.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C14

def create_C14_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.3):
    flowsheet_ID = 'C14'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Landfilling.transportation = sludge_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C15

def create_C15_system(size=10, operation_hours=8760, LA_distance=100, FTE=0.3):
    flowsheet_ID = 'C15'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    LandApplication.transportation = biosolids_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C16

def create_C16_system(size=10, operation_hours=8760, LA_distance=100, FTE=0.35):
    flowsheet_ID = 'C16'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
                                                biogas_CHP_ratio=0, biogas_RNG_ratio=0.9)
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    AnaerobicDigestion.outs[1].price = 0.218
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: make sure the Composting unit is composting + land application
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel_composting'),
                                outs=('compost_cost','methane_composting','nitrous_oxide_composting','sequestered_carbon_dioxide_composting'),
                                solids_distance=LA_distance)
    # TODO: need to decide how to deal with the bulking_agent in composting and its costing
    # TODO: uncertainty range (uniform) 18-36 2005$/tonne
    Composting.ins[1].price = 27/1000/GDPCTPI[2005]*GDPCTPI[2023]
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    Composting.ins[2].price = 4.224/_gal_to_liter*1000/diesel_density
    # TODO: update cost for compost_cost with references
    Composting.outs[0].price = 10/1000
    
    sys = qs.System.from_units(ID='system_C16',
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Bulking_agent', linked_stream=stream.bulking_agent, GlobalWarming=0.041056332)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_composting', linked_stream=stream.diesel_composting, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_composting', linked_stream=stream.methane_composting, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_composting', linked_stream=stream.nitrous_oxide_composting, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_composting', linked_stream=stream.sequestered_carbon_dioxide_composting, GlobalWarming=1)
    
    compost_trucking = qs.ImpactItem(ID='Compost_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    compost_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Composting.transportation = compost_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.Composting.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.Composting.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C17

def create_C17_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.45):
    flowsheet_ID = 'C17'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Landfilling.transportation = sludge_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol+sys.flowsheet.HeatDrying.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C18

def create_C18_system(size=10, operation_hours=8760, LA_distance=100, FTE=0.45):
    flowsheet_ID = 'C18'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    LandApplication.transportation = biosolids_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol+sys.flowsheet.HeatDrying.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C19

def create_C19_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.55):
    flowsheet_ID = 'C19'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
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
    
    # TODO: add ash disposal cost
    Incineration = lsu.Incineration(ID='Incineration',
                                    ins=(HeatDrying-0, 'natural_gas_incineration'),
                                    outs=('ash_incineration','vapor_incineration','methane_IN','nitrous_oxide_IN'),
                                    solids_distance=LF_distance)
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_incineration', linked_stream=stream.ash_incineration, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_IN', linked_stream=stream.methane_IN, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_IN', linked_stream=stream.nitrous_oxide_IN, GlobalWarming=273)
    
    ash_incineration_trucking = qs.ImpactItem(ID='Ash_incineration_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_incineration_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_incineration_trucking.price = (0.00551 + 0.0000541*Incineration.solids_distance)/Incineration.solids_distance
    
    ash_incineration_transportation = qs.Transportation(ID='Ash_incineration_transportation',
                                                        linked_unit=Incineration,
                                                        item=ash_incineration_trucking,
                                                        load_type='mass',
                                                        load=stream.ash_incineration.F_mass,
                                                        load_unit='kg',
                                                        distance=Incineration.solids_distance,
                                                        distance_unit='km',
                                                        # set to 1 h since load = kg/h
                                                        interval='1',
                                                        interval_unit='h')
    Incineration.transportation = ash_incineration_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(-sys.flowsheet.AnaerobicDigestion.outs[1].F_vol+sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.Incineration.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C20

def create_C20_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.45):
    flowsheet_ID = 'C20'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    Landfilling = lsu.Landfilling(ID='Landfilling',
                                  ins=Dewatering-0,
                                  outs=('landfilled_solids','methane_LF','nitrous_oxide_LF','sequestered_carbon_dioxide_LF'),
                                  solids_distance=LF_distance)
    
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Landfilling.transportation = sludge_transportation
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C21

def create_C21_system(size=10, operation_hours=8760, LA_distance=100, LF_distance=100, FTE=0.45):
    flowsheet_ID = 'C21'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
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
    
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    LandApplication.transportation = biosolids_transportation
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C22

def create_C22_system(size=10, operation_hours=8760, LA_distance=100, LF_distance=100, FTE=0.5):
    flowsheet_ID = 'C22'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: make sure the Composting unit is composting + land application
    Composting = lsu.Composting(ID='Composting', ins=(Dewatering-0, 'bulking_agent', 'diesel_composting'),
                                outs=('compost_cost','methane_composting','nitrous_oxide_composting','sequestered_carbon_dioxide_composting'),
                                solids_distance=LA_distance)
    # TODO: need to decide how to deal with the bulking_agent in composting and its costing
    # TODO: uncertainty range (uniform) 18-36 2005$/tonne
    Composting.ins[1].price = 27/1000/GDPCTPI[2005]*GDPCTPI[2023]
    # 2023 weekly average from U.S. EIA: 4.224 $/gallon
    Composting.ins[2].price = 4.224/_gal_to_liter*1000/diesel_density
    # TODO: update cost for compost_cost with references
    Composting.outs[0].price = 10/1000
    
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
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
    
    # BEAM
    N_fertilizer = qs.ImpactItem('N_fertilizer', functional_unit='kg')
    N_fertilizer.add_indicator(GlobalWarming, -3.0)
    
    # BEAM
    P_fertilizer = qs.ImpactItem('P_fertilizer', functional_unit='kg')
    P_fertilizer.add_indicator(GlobalWarming, -2.0)
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Bulking_agent', linked_stream=stream.bulking_agent, GlobalWarming=0.041056332)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_composting', linked_stream=stream.diesel_composting, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_composting', linked_stream=stream.methane_composting, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_composting', linked_stream=stream.nitrous_oxide_composting, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_composting', linked_stream=stream.sequestered_carbon_dioxide_composting, GlobalWarming=1)
    
    compost_trucking = qs.ImpactItem(ID='Compost_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    compost_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Composting.transportation = compost_transportation
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.Composting.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.Composting.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C23

def create_C23_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.6):
    flowsheet_ID = 'C23'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
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
    
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LF', linked_stream=stream.methane_LF, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LF', linked_stream=stream.nitrous_oxide_LF, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Sequestered_carbon_dioxide_LF', linked_stream=stream.sequestered_carbon_dioxide_LF, GlobalWarming=1)
    
    sludge_trucking = qs.ImpactItem(ID='Sludge_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    sludge_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    Landfilling.transportation = sludge_transportation
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production()-sys.flowsheet.Landfilling.electricity_kW*operation_hours)*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C24

def create_C24_system(size=10, operation_hours=8760, LA_distance=100, LF_distance=100, FTE=0.6):
    flowsheet_ID = 'C24'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
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
    
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    # diesel average chemical formula: C12H23
    # https://en.wikipedia.org/wiki/Diesel_fuel (accessed 2025-08-15)
    qs.StreamImpactItem(ID='Diesel_LA', linked_stream=stream.diesel_LA, GlobalWarming=0.4776041 + 44*12/(12*12 + 23*1))
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_LA', linked_stream=stream.methane_LA, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_LA', linked_stream=stream.nitrous_oxide_LA, GlobalWarming=273)
    
    # carbon sequestration
    qs.StreamImpactItem(ID='Carbon_dioxide_LA', linked_stream=stream.carbon_dioxide_LA, GlobalWarming=1)
    
    biosolids_trucking = qs.ImpactItem(ID='Biosolids_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    biosolids_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
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
    LandApplication.transportation = biosolids_transportation
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           N_fertilizer=lambda:sys.flowsheet.LandApplication.nitrogen_mass_flow*operation_hours*20,
           P_fertilizer=lambda:sys.flowsheet.LandApplication.phosphorus_mass_flow*operation_hours*20,
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system C25

def create_C25_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.7):
    flowsheet_ID = 'C25'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
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
                                    outs=('ash_incineration','vapor_incineration','methane_IN','nitrous_oxide_IN'),
                                    solids_distance=LF_distance)
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    Incineration.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    Incineration.outs[0].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_incineration', linked_stream=stream.ash_incineration, GlobalWarming=0.0082744841)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_IN', linked_stream=stream.methane_IN, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_IN', linked_stream=stream.nitrous_oxide_IN, GlobalWarming=273)
    
    ash_incineration_trucking = qs.ImpactItem(ID='Ash_incineration_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_incineration_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_incineration_trucking.price = (0.00551 + 0.0000541*Incineration.solids_distance)/Incineration.solids_distance
    
    ash_incineration_transportation = qs.Transportation(ID='Ash_incineration_transportation',
                                                        linked_unit=Incineration,
                                                        item=ash_incineration_trucking,
                                                        load_type='mass',
                                                        load=stream.ash_incineration.F_mass,
                                                        load_unit='kg',
                                                        distance=Incineration.solids_distance,
                                                        distance_unit='km',
                                                        # set to 1 h since load = kg/h
                                                        interval='1',
                                                        interval_unit='h')
    Incineration.transportation = ash_incineration_transportation
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.Incineration.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           )
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T1

# TODO: need to decide whether to recover N and P or not
def create_T1_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.55):
    flowsheet_ID = 'T1'
    
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
    
    # TODO: add carbon sequestration benefits for hydrochar and biochar in LCA of T systems
    # TODO: add transportation of products (and waste disposal, e.g., in landfilling) in all relevant systems
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
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(F1-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T2

def create_T2_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.55):
    flowsheet_ID = 'T2'
    
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
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(F1-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T3

def create_T3_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.4):
    flowsheet_ID = 'T3'
    
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
    
    # TODO: add ash disposal
    SCWO = lsu.SupercriticalWaterOxidation(ID='SCWO', ins=Dewatering-0, outs=('ash_SCWO','offgas_SCWO'),
                                           solids_distance=LF_distance)
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    SCWO.outs[0].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_SCWO', linked_stream=stream.ash_SCWO, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    ash_SCWO_trucking = qs.ImpactItem(ID='Ash_SCWO_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_SCWO_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_SCWO_trucking.price = (0.00551 + 0.0000541*SCWO.solids_distance)/SCWO.solids_distance
    
    ash_SCWO_transportation = qs.Transportation(ID='Ash_SCWO_transportation',
                                                linked_unit=SCWO,
                                                item=ash_SCWO_trucking,
                                                load_type='mass',
                                                load=stream.ash_SCWO.F_mass,
                                                load_unit='kg',
                                                distance=SCWO.solids_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    SCWO.transportation = ash_SCWO_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T4

def create_T4_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.7):
    flowsheet_ID = 'T4'
    
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
    
    Pyrolysis = lsu.Pyrolysis(ID='Pyrolysis', ins=HeatDrying-0,
                              outs=('biooil','biochar','pyrogas','methane_pyrolysis','nitrous_oxide_pyrolysis'))
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(Pyrolysis-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_pyrolysis', linked_stream=stream.methane_pyrolysis, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_pyrolysis', linked_stream=stream.nitrous_oxide_pyrolysis, GlobalWarming=273)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T5

def create_T5_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.7):
    flowsheet_ID = 'T5'
    
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
    
    Gasification = lsu.Gasification(ID='Gasification', ins=HeatDrying-0,
                                    outs=('tar','biochar','syngas','methane_gasification','nitrous_oxide_gasification'))
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(Gasification-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_gasification', linked_stream=stream.methane_gasification, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_gasification', linked_stream=stream.nitrous_oxide_gasification, GlobalWarming=273)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T6

def create_T6_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.7):
    flowsheet_ID = 'T6'
    
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
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(F1-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T7

def create_T7_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.7):
    flowsheet_ID = 'T7'
    
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
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(F1-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T8

def create_T8_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.55):
    flowsheet_ID = 'T8'
    
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
    
    SCWO = lsu.SupercriticalWaterOxidation(ID='SCWO', ins=Dewatering-0, outs=('ash_SCWO','offgas_SCWO'),
                                           solids_distance=LF_distance)
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    SCWO.outs[0].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_SCWO', linked_stream=stream.ash_SCWO, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    ash_SCWO_trucking = qs.ImpactItem(ID='Ash_SCWO_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_SCWO_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_SCWO_trucking.price = (0.00551 + 0.0000541*SCWO.solids_distance)/SCWO.solids_distance
    
    ash_SCWO_transportation = qs.Transportation(ID='Ash_SCWO_transportation',
                                                linked_unit=SCWO,
                                                item=ash_SCWO_trucking,
                                                load_type='mass',
                                                load=stream.ash_SCWO.F_mass,
                                                load_unit='kg',
                                                distance=SCWO.solids_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    SCWO.transportation = ash_SCWO_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T9

def create_T9_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.85):
    flowsheet_ID = 'T9'
    
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
    
    Pyrolysis = lsu.Pyrolysis(ID='Pyrolysis', ins=HeatDrying-0,
                              outs=('biooil','biochar','pyrogas','methane_pyrolysis','nitrous_oxide_pyrolysis'))
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(Pyrolysis-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_pyrolysis', linked_stream=stream.methane_pyrolysis, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_pyrolysis', linked_stream=stream.nitrous_oxide_pyrolysis, GlobalWarming=273)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T10

def create_T10_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.85):
    flowsheet_ID = 'T10'
    
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
    
    Gasification = lsu.Gasification(ID='Gasification', ins=HeatDrying-0,
                                    outs=('tar','biochar','syngas','methane_gasification','nitrous_oxide_gasification'))
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(Gasification-2, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_gasification', linked_stream=stream.methane_gasification, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_gasification', linked_stream=stream.nitrous_oxide_gasification, GlobalWarming=273)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T11

def create_T11_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.7):
    flowsheet_ID = 'T11'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
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
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T12

def create_T12_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.7):
    flowsheet_ID = 'T12'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
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
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T13

def create_T13_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.55):
    flowsheet_ID = 'T13'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
                                                biogas_CHP_ratio=0.9, biogas_RNG_ratio=0)
    
    # note disposal_cost (add_OPEX here, and other similar funcions) does not work since TEA is from BioSTEAM, but not QSDsan
    Dewatering = lsu.Dewatering(ID='Dewatering', ins=(AnaerobicDigestion-0, 'polymer_dewatering'),
                                outs=('dewatered_solids','reject_dewatering','methane_dewatering'))
    Dewatering.ins[1].price = 2.6282/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2023]
    
    SCWO = lsu.SupercriticalWaterOxidation(ID='SCWO', ins=Dewatering-0, outs=('ash_SCWO','offgas_SCWO'),
                                           solids_distance=LF_distance)
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    SCWO.outs[0].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(AnaerobicDigestion-1, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_SCWO', linked_stream=stream.ash_SCWO, GlobalWarming=0.0082744841)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    
    ash_SCWO_trucking = qs.ImpactItem(ID='Ash_SCWO_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_SCWO_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_SCWO_trucking.price = (0.00551 + 0.0000541*SCWO.solids_distance)/SCWO.solids_distance
    
    ash_SCWO_transportation = qs.Transportation(ID='Ash_SCWO_transportation',
                                                linked_unit=SCWO,
                                                item=ash_SCWO_trucking,
                                                load_type='mass',
                                                load=stream.ash_SCWO.F_mass,
                                                load_unit='kg',
                                                distance=SCWO.solids_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    SCWO.transportation = ash_SCWO_transportation
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T14

def create_T14_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.85):
    flowsheet_ID = 'T14'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
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
    
    Pyrolysis = lsu.Pyrolysis(ID='Pyrolysis', ins=HeatDrying-0,
                              outs=('biooil','biochar','pyrogas','methane_pyrolysis','nitrous_oxide_pyrolysis'))
    
    GasMixer = qsu.Mixer('GasMixer', ins=(AnaerobicDigestion-1, Pyrolysis-2),
                         outs=('fuel_gas'), init_with='WasteStream')
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_pyrolysis', linked_stream=stream.methane_pyrolysis, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_pyrolysis', linked_stream=stream.nitrous_oxide_pyrolysis, GlobalWarming=273)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys

#%% system T15

def create_T15_system(size=10, operation_hours=8760, LF_distance=100, FTE=0.85):
    flowsheet_ID = 'T15'
    
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
    
    # TODO: consider replacing natural_gas with biogas (CH4 + CO2) or create a stream for CO2
    AnaerobicDigestion = lsu.AnaerobicDigestion(ID='AnaerobicDigestion', ins=Thickening-0,
                                                outs=('digested_sludge','natural_gas_AD','methane_AD'),
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
    
    Gasification = lsu.Gasification(ID='Gasification', ins=HeatDrying-0,
                                    outs=('tar','biochar','syngas','methane_gasification','nitrous_oxide_gasification'))
    
    GasMixer = qsu.Mixer('GasMixer', ins=(AnaerobicDigestion-1, Gasification-2),
                         outs=('fuel_gas'), init_with='WasteStream')
    
    # TODO: consider adding pyrogas to CHP (for other unit sending gas to CHP, may add fugitive emissions together in the CHP unit)
    CHP = lsu.CombinedHeatPower(ID='CHP', ins=(GasMixer-0, 'natural_gas_CHP', 'air_CHP'),
                                outs=('emission','ash_CHP'), supplement_power_utility=False,
                                solids_distance=LF_distance)
    CHP.lifetime = 20
    # from _heat_utility.py (BioSTEAM): 3.49672 $/kmol
    # assume the MW of natural gas is 16.04 g/mol (same as CH4, probably consistent with BioSTEAM)
    CHP.ins[1].price = 0.218
    # 1.41 MM 2016$/year for 4270/4279 kg/h ash, 7880 annual operating hours, https://doi.org/10.2172/1483234
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2023]
    
    # TODO: check: construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
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
    Electricity.add_indicator(GlobalWarming, 0.67877501)
    
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.1194374)
    
    Natural_gas_E = qs.ImpactItem('Natural_gas_E', functional_unit='MJ')
    Natural_gas_E.add_indicator(GlobalWarming, 0.036990763)
    
    Natural_gas_V = qs.ImpactItem('Natural_gas_V', functional_unit='m3')
    Natural_gas_V.add_indicator(GlobalWarming, 0.47016123 + natural_gas_density/16*44)
        
    Cooling = qs.ImpactItem('Cooling', functional_unit='MJ')
    Cooling.add_indicator(GlobalWarming, 0.065877932)
    
    Deionized_water = qs.ImpactItem('Deionized_water', functional_unit='kg')
    Deionized_water.add_indicator(GlobalWarming, 0.00045239247)
    
    # TODO: will any other places have deionized water as well
    def deionized_water_quantity():
        try:
            # TODO: is there a direct number or a better way for cooling tower chemicals
            # assume the ratio between the CI of the cooling_tower_chemicals and the CI of cooling_tower_makeup_water is the same as the ratio of the prices of these two streams
            # note LCA for water_steam, cooling_tower_makeup_water, and cooling_tower_chemicals were included in the 'Other' category while it should be in the 'Stream' category
            # the effect is minimal since (i) this part of LCA is small and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
            result = (sys.flowsheet.CT.ins[1].F_mass+sys.flowsheet.CT.ins[2].F_mass*1.7842/0.0002)*operation_hours*20
        except AttributeError as e:
            if 'no registered item' in str(e) and 'CT' in str(e):
                result = 0
        return result
    
    qs.StreamImpactItem(ID='Polymer_thickening', linked_stream=stream.polymer_thickening, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Polymer_dewatering', linked_stream=stream.polymer_dewatering, GlobalWarming=3.1940311)
    qs.StreamImpactItem(ID='Ash_CHP', linked_stream=stream.ash_CHP, GlobalWarming=0.0082744841)
    
    # fugitive emissions
    qs.StreamImpactItem(ID='Methane_AD', linked_stream=stream.methane_AD, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_dewatering', linked_stream=stream.methane_dewatering, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Methane_gasification', linked_stream=stream.methane_gasification, GlobalWarming=29.8)
    qs.StreamImpactItem(ID='Nitrous_oxide_gasification', linked_stream=stream.nitrous_oxide_gasification, GlobalWarming=273)
    
    ash_CHP_trucking = qs.ImpactItem(ID='Ash_CHP_trucking', functional_unit='kg*km')
    # based on one-way distance, empty return trips included
    ash_CHP_trucking.add_indicator(GlobalWarming, 0.13673337/1000)
    
    # for sludge (with an assumed density of 1040 kg/m3): 4.56 $/m3, 0.072 $/m3/mile (likely 2015$)
    # https://doi.org/10.1016/j.tra.2015.02.001
    # converted to 2023$/kg/km
    ash_CHP_trucking.price = (0.00551 + 0.0000541*CHP.solids_distance)/CHP.solids_distance
    
    ash_CHP_transportation = qs.Transportation(ID='Ash_CHP_transportation',
                                               linked_unit=CHP,
                                               item=ash_CHP_trucking,
                                               load_type='mass',
                                               load=stream.ash_CHP.F_mass,
                                               load_unit='kg',
                                               distance=CHP.solids_distance,
                                               distance_unit='km',
                                               # set to 1 h since load = kg/h
                                               interval='1',
                                               interval_unit='h')
    CHP.transportation = ash_CHP_transportation
    
    # TODO: for both TEA and LCA, consider use a lifetime of 50 years (a duration from 2023 to 2073) or maybe not since need to avoid the replacement of thermochemical units in the greenfield construction strategy (if this strategy is kept)
    # TODO: when calculate the total quantity of LCA items, make sure the lifetime is consistent with TEA and LCA
    qs.LCA(system=sys, lifetime=20, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*20,
           Steam=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if 'steam' in i.ID)),
           Natural_gas_E=lambda:abs(sum(i.duty/1000*operation_hours*20 for i in sys.heat_utilities if i.ID == 'natural_gas')),
           Natural_gas_V=lambda:(sys.flowsheet.HeatDrying.ins[1].F_vol+sys.flowsheet.CHP.ins[1].F_vol)*operation_hours*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Deionized_water=lambda:deionized_water_quantity())
    
    FTE_labor_cost = (0.34/labor_index[2014]*labor_index[2023]+\
                      0.48/labor_index[2014]*labor_index[2023]*size/100)*10**6
    
    # TODO: check all parameters in TEA, especially the ones listed here
    # TODO: consider adding IRR as a WRRF typology parameter?
    # TODO: income tax may be found in Steward et al. ES&T, 2023
    create_tea(sys,
               duration=(2023, 2043),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=FTE*FTE_labor_cost)
    
    return sys