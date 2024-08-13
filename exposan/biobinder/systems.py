# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References
[1] Snowden-Swan et al., Wet Waste Hydrothermal Liquefaction and Biocrude Upgrading to Hydrocarbon Fuels:
    2021 State of Technology; PNNL-32731; Pacific Northwest National Lab. (PNNL), Richland, WA (United States), 2022.
    https://doi.org/10.2172/1863608.
'''

# !!! Temporarily ignoring warnings
import warnings
warnings.filterwarnings('ignore')

import os, biosteam as bst, qsdsan as qs
# from biosteam.units import IsenthalpicValve
# from biosteam import settings
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import data_path as htl_data_path
from exposan.biobinder import (
    data_path,
    results_path,
    _load_components,
    _load_process_settings,
    create_tea,
    _units as u
    )


__all__ = ('create_system',)

#!!! Placeholder function for now, update when flowsheet ready
def create_system():
    pass


# %%

# Create and set flowsheet
configuration = 'DC' # decentralized HTL, centralized upgrading
flowsheet_ID = f'biobinder_{configuration}'
flowsheet = qs.Flowsheet(flowsheet_ID)
qs.main_flowsheet.set_flowsheet(flowsheet)

_load_components()
_load_process_settings()

# Desired feedstock flowrate, in dry kg/hr
decentralized_dry_flowrate = 11.46 # feedstock mass flowrate, dry kg/hr
N_decentralized_HTL = 1000 # number of parallel HTL reactor, PNNL is about 1900x of UIUC pilot reactor
target_HTL_solid_loading = 0.2

# %%

# =============================================================================
# Hydrothermal Liquefaction
# =============================================================================

scaled_feedstock = qs.WasteStream('scaled_feedstock')
# fresh_process_water = qs.WasteStream('fresh_process_water')

# Adjust feedstock composition
FeedstockScaler = u.Scaler(
    'FeedstockScaler', ins=scaled_feedstock, outs='feedstock',
    scaling_factor=N_decentralized_HTL, reverse=True,
    )

ProcessWaterScaler = u.Scaler(
    'ProcessWaterScaler', ins='scaled_process_water', outs='htl_process_water',
    scaling_factor=N_decentralized_HTL, reverse=True,
    )

FeedstockTrans = u.Transportation(
    'FeedstockTrans',
    ins=(FeedstockScaler-0, 'feedstock_trans_surrogate'),
    outs=('transported_feedstock',),
    N_unit=N_decentralized_HTL,
    copy_ins_from_outs=True,
    transportation_distance=78, # km ref [1]
    )

FeedstockCond = u.Conditioning(
    'FeedstockCond', ins=(FeedstockTrans-0, ProcessWaterScaler-0),
    outs='conditioned_feedstock',
    feedstock_composition=u.salad_dressing_waste_composition,
    feedstock_dry_flowrate=decentralized_dry_flowrate,
    N_unit=N_decentralized_HTL,
    )

HTL = u.PilotHTL(
    'HTL', ins=FeedstockCond-0, outs=('hydrochar','HTL_aqueous','biocrude','HTL_offgas'),
    afdw_yields=u.salad_dressing_waste_yields,
    N_unit=N_decentralized_HTL,
    )
HTL.register_alias('PilotHTL')


# %%

# =============================================================================
# Biocrude Upgrading
# =============================================================================

BiocrudeDeashing = u.BiocrudeDeashing(
    'BiocrudeDeashing', ins=HTL-2, outs=('deashed_biocrude', 'biocrude_ash'),
    N_unit=N_decentralized_HTL,)

BiocrudeAshScaler = u.Scaler(
    'BiocrudeAshScaler', ins=BiocrudeDeashing-1, outs='scaled_biocrude_ash',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )

BiocrudeDewatering = u.BiocrudeDewatering(
    'BiocrudeDewatering', ins=BiocrudeDeashing-0, outs=('dewatered_biocrude', 'biocrude_water'),
    N_unit=N_decentralized_HTL,)

BiocrudeWaterScaler = u.Scaler(
    'BiocrudeWaterScaler', ins=BiocrudeDewatering-1, outs='scaled_biocrude_water',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )

BiocrudeTrans = u.Transportation(
    'BiocrudeTrans',
    ins=(BiocrudeDewatering-0, 'biocrude_trans_surrogate'),
    outs=('transported_biocrude',),
    N_unit=N_decentralized_HTL,
    transportation_distance=78, # km ref [1]
    )

BiocrudeScaler = u.Scaler(
    'BiocrudeScaler', ins=BiocrudeTrans-0, outs='scaled_biocrude',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )

BiocrudeSplitter = u.BiocrudeSplitter(
    'BiocrudeSplitter', ins=BiocrudeScaler-0, outs='splitted_biocrude',
    cutoff_Tb=343+273.15, light_frac=0.5316)

# Shortcut column uses the Fenske-Underwood-Gilliland method,
# better for hydrocarbons according to the tutorial
# https://biosteam.readthedocs.io/en/latest/API/units/distillation.html
FracDist = u.ShortcutColumn(
    'FracDist', ins=BiocrudeSplitter-0,
    outs=('biocrude_light','biocrude_heavy'),
    LHK=('Biofuel', 'Biobinder'), # will be updated later
    P=50*6894.76, # outflow P, 50 psig
    # Lr=0.1, Hr=0.5,
    y_top=188/253, x_bot=53/162,
    k=2, is_divided=True)
@FracDist.add_specification
def adjust_LHK():
    FracDist.LHK = (BiocrudeSplitter.light_key, BiocrudeSplitter.heavy_key)
    FracDist._run()

LightFracStorage = qsu.StorageTank(
    'LightFracStorage',
    FracDist-0, outs='biofuel_additives',
    tau=24*7, vessel_material='Stainless steel')
HeavyFracStorage = qsu.StorageTank(
    'HeavyFracStorage', FracDist-1, outs='biobinder',
    tau=24*7, vessel_material='Stainless steel')


# %%

# =============================================================================
# Aqueous Product Treatment
# =============================================================================

AqueousFiltration = u.AqueousFiltration(
    'AqueousFiltration',
    ins=(HTL-1,),
    outs=('fertilizer', 'recycled_water', 'filtered_solids'),
    N_unit=N_decentralized_HTL,)

FertilizerScaler = u.Scaler(
    'FertilizerScaler', ins=AqueousFiltration-0, outs='scaled_fertilizer',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )

RecycledWaterScaler = u.Scaler(
    'RecycledWaterScaler', ins=AqueousFiltration-1, outs='scaled_recycled_water',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )

FilteredSolidsScaler = u.Scaler(
    'FilteredSolidsScaler', ins=AqueousFiltration-2, outs='filterd_solids',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )


# %%

# =============================================================================
# Facilities and waste disposal
# =============================================================================

# Scale flows
HydrocharScaler = u.Scaler(
    'HydrocharScaler', ins=HTL-0, outs='scaled_hydrochar',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )
@HydrocharScaler.add_specification
def scale_feedstock_flows():
    FeedstockTrans._run()
    FeedstockScaler._run()
    ProcessWaterScaler._run()

GasScaler = u.Scaler(
    'GasScaler', ins=HTL-3, outs='scaled_gas',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )

# Potentially recycle the water from aqueous filtration (will be ins[2])
ProcessWaterCenter = u.ProcessWaterCenter(
    'ProcessWaterCenter',
    process_water_streams=[ProcessWaterScaler.ins[0]],
    )


# No need to consider transportation as priced are based on mass
AshDisposal = u.Disposal('AshDisposal', ins=(BiocrudeAshScaler-0, FilteredSolidsScaler-0),
                         outs=('ash_disposal', 'ash_others'),
                         exclude_components=('Water',))

WWDisposal = u.Disposal('WWDisposal', ins=BiocrudeWaterScaler-0,
                        outs=('ww_disposal', 'ww_others'),
                        exclude_components=('Water',))

# Heat exchanger network
# 86 K: Jones et al. PNNL, 2014
HXN = qsu.HeatExchangerNetwork('HXN', T_min_app=86, force_ideal_thermo=True)


# %%

# =============================================================================
# Assemble System
# =============================================================================

sys = qs.System.from_units(
    f'sys_{configuration}',
    units=list(flowsheet.unit), 
    operating_hours=7920, # same as the HTL module, about 90% uptime
    )
sys.register_alias('sys')
stream = sys.flowsheet.stream

# =============================================================================
# TEA
# =============================================================================

cost_year = 2020

# U.S. Energy Information Administration (EIA) Annual Energy Outlook (AEO)
# GDP_indices = {
#   2003: 0.808,
#   2005: 0.867,
#   2007: 0.913,
#   2008: 0.941,
#   2009: 0.951,
#   2010: 0.962,
#   2011: 0.983,
#   2012: 1.000,
#   2013: 1.014,
#   2014: 1.033,
#   2015: 1.046,
#   2016: 1.059,
#   2017: 1.078,
#   2018: 1.100,
#   2019: 1.123,
#   2020: 1.133,
#   2021: 1.181,
#   2022: 1.269,
#   2023: 1.322,
#   2024: 1.354,
#     }

#Federal Reserve Economic Data, Personal Consumption Expenditures: Chain-type Price Index, Index 2017=1.00, Annual, Seasonally Adjusted

PCE_indices = {
    2000: 0.738,
    2001: 0.753,
    2002: 0.763,
    2003: 0.779,
    2004: 0.798,
    2005: 0.821,
    2006: 0.844,
    2007: 0.866,
    2008: 0.892,
    2009: 0.889,
    2010: 0.905,
    2011: 0.928,
    2012: 0.945,
    2013: 0.958,
    2014: 0.971,
    2015: 0.973,
    2016: 0.983,
    2017: 1.000,
    2018: 1.020,
    2019: 1.035,
    2020: 1.046,
    2021: 1.090,
    2022: 1.160,
    2023: 1.204,
    2024: 1.220,
    }

# Inputs
scaled_feedstock.price = -69.14/907.185 # tipping fee 69.14±21.14 for IL, https://erefdn.org/analyzing-municipal-solid-waste-landfill-tipping-fees/

# Utilities, price from Table 17.1 in Seider et al., 2016$
# Use bst.HeatUtility.cooling_agents/heating_agents to see all the heat utilities
Seider_factor = PCE_indices[cost_year]/PCE_indices[2016]

transport_cost = 50/1e3 * PCE_indices[cost_year]/PCE_indices[2016] # $/kg ref [1]
FeedstockTrans.transportation_cost = BiocrudeTrans.transportation_cost = transport_cost

ProcessWaterCenter.process_water_price = 0.8/1e3/3.785*Seider_factor # process water for moisture adjustment

hps = bst.HeatUtility.get_agent('high_pressure_steam') # 450 psig
hps.regeneration_price = 17.6/(1000/18)*Seider_factor

mps = bst.HeatUtility.get_agent('medium_pressure_steam') # 150 psig
mps.regeneration_price = 15.3/(1000/18)*Seider_factor

lps = bst.HeatUtility.get_agent('low_pressure_steam') # 50 psig
lps.regeneration_price = 13.2/(1000/18)*Seider_factor

heating_oil = bst.HeatUtility.get_agent('HTF') # heat transfer fluids, added in the HTL module
crude_oil_density = 3.205 # kg/gal, GREET1 2023, "Fuel_Specs", US conventional diesel
heating_oil.regeneration_price = 3.5/crude_oil_density*Seider_factor

cw = bst.HeatUtility.get_agent('cooling_water')
cw.regeneration_price = 0.1*3.785/(1000/18)*Seider_factor # $0.1/1000 gal to $/kmol

for i in (hps, mps, lps, heating_oil, cw):
    i.heat_transfer_price = 0

# Annual Energy Outlook 2023 https://www.eia.gov/outlooks/aeo/data/browser/
# Table 8. Electricity Supply, Disposition, Prices, and Emissions
# End-Use Prices, Industrial, nominal 2024 value in $/kWh
bst.PowerUtility.price = 0.07*Seider_factor

# Waste disposal
AshDisposal.disposal_price = 0.17*Seider_factor # deashing, landfill price
WWDisposal.disposal_price = 0.33*Seider_factor # dewater, for organics removed

# Products
diesel_density = 3.167 # kg/gal, GREET1 2023, "Fuel_Specs", US conventional diesel
biofuel_additives = stream.biofuel_additives
biofuel_additives.price = 4.07/diesel_density # diesel, https://afdc.energy.gov/fuels/prices.html

hydrochar = stream.hydrochar
hydrochar.price = 0

biobinder = stream.biobinder
biobinder.price = 0.67 # bitumnous, https://idot.illinois.gov/doing-business/procurements/construction-services/transportation-bulletin/price-indices.html


# Other TEA assumptions
bst.CE = qs.CEPCI_by_year[cost_year]
lifetime = 30

base_labor = 338256 # for 1000 kg/hr

tea = create_tea(
    sys,
    labor_cost=lambda: (scaled_feedstock.F_mass-scaled_feedstock.imass['Water'])/1000*base_labor,
    land=0, #!!! need to be updated
    )

# To see out-of-boundary-limits units
# tea.OSBL_units

# =============================================================================
# LCA
# =============================================================================

# Load impact indicators, TRACI
clear_lca_registries()
qs.ImpactIndicator.load_from_file(os.path.join(data_path, 'impact_indicators.csv'))
qs.ImpactItem.load_from_file(os.path.join(data_path, 'impact_items.xlsx'))

# Add impact for streams
streams_with_impacts = [i for i in sys.feeds+sys.products if (
    i.isempty() is False and 
    i.imass['Water']!=i.F_mass and
    'surrogate' not in i.ID
    )]
for i in streams_with_impacts: print (i.ID)

# scaled_feedstock
# biofuel_additives
# biobinder
# scaled_gas
feedstock_item = qs.StreamImpactItem(
    ID='feedstock_item',
    linked_stream=scaled_feedstock,
    Acidification=0,
    Ecotoxicity=0,
    Eutrophication=0,
    GlobalWarming=0,
    OzoneDepletion=0,
    PhotochemicalOxidation=0,
    Carcinogenics=0,
    NonCarcinogenics=0,
    RespiratoryEffects=0
    )
qs.ImpactItem.get_item('Diesel').linked_stream = biofuel_additives

#!!! Need to get heating duty
lca = qs.LCA(
    system=sys,
    lifetime=lifetime,
    uptime_ratio=sys.operating_hours/(365*24),
    Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
    # Heating=lambda:sys.get_heating_duty()/1000*lifetime,
    Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
    )


# %%

def simulate_and_print(save_report=False):
    sys.simulate()
    
    biobinder.price = biobinder_price = tea.solve_price(biobinder)
    print(f'Minimum selling price of the biobinder is ${biobinder_price:.2f}/kg.')
    c = qs.currency
    for attr in ('NPV','AOC', 'sales', 'net_earnings'):
        uom = c if attr in ('NPV', 'CAPEX') else (c+('/yr'))
        print(f'{attr} is {getattr(tea, attr):,.0f} {uom}')
        
    all_impacts = lca.get_allocated_impacts(streams=(biobinder,), operation_only=True, annual=True)
    GWP = all_impacts['GlobalWarming']/biobinder.F_mass
    print(f'Global warming potential of the biobinder is {GWP:.2f} kg CO2e/kg.')
    if save_report:
        # Use `results_path` and the `join` func can make sure the path works for all users
        sys.save_report(file=os.path.join(results_path, 'sys.xlsx'))

if __name__ == '__main__':
    simulate_and_print()
    