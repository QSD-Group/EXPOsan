#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>
      
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import biosteam as bst, qsdsan as qs

__all__ = (
    '_HHV_per_GGE',
    '_load_process_settings',
    'dry_flowrate',
    'feedstock_composition',
    'gwp_dct',
    'HTL_yields',
    'price_dct',
    'tea_kwargs',
    'uptime_ratio',
    )

_ton_to_kg = 907.185
_MJ_to_MMBtu = 0.000948
_HHV_per_GGE = 46.52*2.82 # MJ/gal
# DOE properties
# https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels
# Conventional Gasoline: HHV=46.52 MJ/kg, rho=2.82 kg/gal
# U.S. Conventional Diesel: HHV=45.76 MJ/kg, rho=3.17 kg/gal
# diesel_density = 3.167 # kg/gal, GREET1 2023, "Fuel_Specs", US conventional diesel
# Fuel properties
# https://afdc.energy.gov/fuels/properties

ratio = 1
tpd = 110*ratio # dry mass basis
uptime_ratio = 0.9

dry_flowrate = tpd*907.185/(24*uptime_ratio) # 110 dry sludge tpd [1]

moisture = 0.7580
ash = (1-moisture)*0.0614
feedstock_composition = {
    'Water': moisture,
    'Lipids': (1-moisture)*0.5315,
    'Proteins': (1-moisture)*0.0255,
    'Carbohydrates': (1-moisture)*0.3816,
    'Ash': ash,
    }

# Salad dressing waste, char is separated out through distillation
HTL_yields = {
    'gas': 0.006,
    'aqueous': 0.192,
    'biocrude': 0.802,
    'char': 0,
    }

# All in 2020 $/kg unless otherwise noted, needs to do a thorough check to update values
tea_indices = qs.utils.indices.tea_indices
cost_year = 2020
PCE_indices = tea_indices['PCEPI']

AEO_year = 2022
AEO_factor = PCE_indices[cost_year]/PCE_indices[AEO_year]

Seider_year = 2016
Seider_factor = PCE_indices[cost_year]/PCE_indices[Seider_year]

SnowdenSwan_year = 2016
SnowdenSwan_factor = PCE_indices[cost_year]/PCE_indices[SnowdenSwan_year]

bst_utility_price = bst.stream_utility_prices
price_dct = {
    'tipping': -39.7/1e3*SnowdenSwan_factor, # PNNL 2022, -$39.7/wet tonne is the weighted average
    'trans_feedstock': 50/1e3*SnowdenSwan_factor, # $50 dry tonne for 78 km, Snowden-Swan PNNL 32731 
    'trans_biocrude': 0.092*SnowdenSwan_factor, # $0.092/GGE of biocrude, 100 miles, Snowden-Swan PNNL 32731
    # $6.77/MMBtu in 2018 from petroleum and coal products in https://www.eia.gov/todayinenergy/detail.php?id=61763
    # 141.88 MJ/kg from https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels
    # This is too high, ~$50/kg H2, but on par with CLEAN CITIES and COMMUNITIES Alternative Fuel Price Report
    # $33.37/GGE, or $33.37/kg (1 kg H2 is 1 GGE, https://afdc.energy.gov/fuels/properties)
    # 'H2': 6.77/(141.88*_MJ_to_MMBtu),
    'H2': 1.61, # Feng et al., 2024, in 2020$
    'HCcatalyst': 3.52, # Fe-ZSM5, CatCost modified from ZSM5, in 2020$
    'HTcatalyst': 75.18, # Pd/Al2O3, CatCost modified from 2% Pt/TiO2, in 2020$
    'natural_gas': 0.213/0.76**Seider_factor, # $0.213/SCM, $0.76 kg/SCM per https://www.henergy.com/conversion
    'process_water': 0.27/1e3*Seider_factor, # $0.27/m3, higher than $0.80/1,000 gal
    'gasoline': 3.32*AEO_factor, # EIA AEO 2023, Table 12, Transportation Sector, 2024 price in 2022$
    'jet': 2.92*AEO_factor, # EIA AEO 2023, Table 12, Transportation Sector, 2024 price in 2022$
    'diesel': 4.29*AEO_factor, # EIA AEO 2023, Table 12, Transportation Sector, 2024 price in 2022$
    # Fertilizers from USDA, for the week ending 10/4/2024, https://www.ams.usda.gov/mnreports/ams_3195.pdf
    # Not good that it's just for one week, but it has negligible impacts on the results
    'N': 0.90, # recovered N in $/kg N
    'P': 1.14, # recovered P in $/kg P
    'K': 0.81, # recovered K in $/kg K
    'solids': -0.17*Seider_factor,
    'COD': -0.3676, # $/kg, Li et al., 2023; Seider has 0.33 for organics removed
    }

# GREET 2023, unless otherwise noted

# Ecoinvent 3.10, cutoff, market group for transport, freight, lorry, unspecified, GLO
# https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/17617/impact_assessment
gwp_trans = 0.152/1e3 # 0.152 kg CO2e/tonne/km

gwp_dct = {
    'feedstock': 0,
    'landfill': 400/1e3, # nearly 400 kg CO2e/tonne, Nordahl et al., 2020
    'composting': -41/1e3, # -41 kg CO2e/tonne, Nordahl et al., 2020
    'anaerobic_digestion': (-36-2)/2, # -36 to -2 kg CO2e/tonne, Nordahl et al., 2020
    # Ecoinvent 3.10, cutoff, market group for transport, freight, lorry, unspecified, GLO
    # https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/17617/impact_assessment
    'trans_feedstock': gwp_trans*78, # 78 km, Snowden-Swan PNNL 32731
    'trans_biocrude': gwp_trans*100*1.6, # 100 miles, Snowden-Swan PNNL 32731
    'H2': 11.0469, # Mix: Central Plants: Compressed G.H2 production (100% from Natural Gas)
    'H2_electrolysis': 0.9514, # Central Plants: Compressed G.H2 production from Electrolysis with HGTR
    'HCcatalyst': 6.1901, # Feng et al., 2024
    'natural_gas': 0.3877+44/16, # NA NG from Shale and Conventional Recovery, with CO2 after combustion
    'process_water': 0,
    'electricity': 0.4181, # kg CO2e/kWh Non Distributed - U.S. Mix
    'steam': 86.3928/1e3, # 86.3928 g CO2e/MJ, Mix: Natural Gas and Still Gas
    'cooling': 0.066033, # kg CO2e/MJ, Feng et al., 2024
    'gasoline': 2.3722, # kg CO2e/gal, 0.8415 kg CO2e/kg, no combustion emission, Gasoline Blendstock from Crude oil for Use in US Refineries
    'jet': 1.4599, # kg CO2e/gal, 0.4809 kg CO2e/kg, no combustion emission, Conventional Jet Fuel from Crude Oil
    'diesel': 2.0696, # kg CO2e/gal, 0.6535 kg CO2e/kg, no combustion emission, Conventional Diesel from Crude Oil for US Refineries
    'N': -3.46, # 3.46 kg CO2e/kg N, Mix: Nitrogen Average
    'P': -1.6379*142/(31*2), # 1.6379 kg CO2e/kg P2O5, Mix: Phosphate (P2O5) from MAP and DAP
    'K': -0.4830*94/(39*2), # 0.4830 kg CO2e/kg K2O, Potassium Oxide Production
    'COD': 1.7, # Li et al., 2023
    'wastewater': 0.2851/1e3, # Industrial Wastewater Treatment
    }
gwp_dct['HTcatalyst'] = gwp_dct['HCcatalyst']
gwp_dct['solids'] = gwp_dct['trans_feedstock'] # only account for transportation

labor_indices = tea_indices['labor']
size_ratio = tpd/1339
tea_kwargs = dict(
    IRR=0.1,
    duration=(2020, 2050),
    income_tax=0.21,
    finance_interest=0.08,
    warehouse=0.04,
    site_development=0.1, # Snowden-Swan et al. 2022
    additional_piping=0.045,
    labor_cost=2.36e6*size_ratio*labor_indices[cost_year]/labor_indices[2011], # PNNL 2014
    land=0., #!!! need to update
    )

def _load_process_settings():
    bst.CE = tea_indices['CEPCI'][cost_year]
    
    # Utilities, price from Table 17.1 in Seider et al., 2016$
    # Use bst.HeatUtility.cooling_agents/heating_agents to see all the heat utilities
    
    # Steams are provided by CHP thus cost already considered
    # setting the regeneration price to 0 or not will not affect the final results
    # as the utility cost will be positive for the unit that consumes it
    # but negative for HXN/CHP as they produce it
    hps = bst.HeatUtility.get_agent('high_pressure_steam') # 450 psig, $17.6/1000 kg to $/kmol
    hps.regeneration_price = 17.6/1000*18*Seider_factor
    
    mps = bst.HeatUtility.get_agent('medium_pressure_steam') # 150 psig, $15.3/1000 kg to $/kmol
    mps.regeneration_price = 15.3/1000*18*Seider_factor
    
    lps = bst.HeatUtility.get_agent('low_pressure_steam') # 50 psig, $13.2/1000 kg to $/kmol
    lps.regeneration_price = 13.2/1000*18*Seider_factor
    
    cw = bst.HeatUtility.get_agent('cooling_water')
    cw.regeneration_price = 0.1*3.785/1000*18*Seider_factor # $0.1/1000 gal to $/kmol (higher than the SI unit option)
    
    chilled = bst.HeatUtility.get_agent('chilled_water')
    chilled.heat_transfer_price = 5/1e6*Seider_factor # $5/GJ to $/kJ
    chilled.regeneration_price = 0
    
    for i in (hps, mps, lps, cw):
        i.heat_transfer_price = 0

    bst.PowerUtility.price = 0.076*AEO_factor # EIA AEO 2023, Table 8, End-Use Industrial Sector, 2024 price in 2022$