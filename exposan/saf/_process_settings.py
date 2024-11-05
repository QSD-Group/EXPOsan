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
    'HTL_yields',
    'price_dct',
    'tea_kwargs',
    'uptime_ratio',
    )

_ton_to_kg = 907.185
_HHV_per_GGE = 46.52*2.82 # MJ/gal
# DOE properties
# https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels
# Conventional Gasoline: HHV=46.52 MJ/kg, rho=2.82 kg/gal
# U.S. Conventional Diesel: HHV=45.76 MJ/kg, rho=3.17 kg/gal
# diesel_density = 3.167 # kg/gal, GREET1 2023, "Fuel_Specs", US conventional diesel

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

# Salad dressing waste
HTL_yields = {
    'gas': 0.006,
    'aqueous': 0.192,
    'biocrude': 0.802-ash,
    'char': ash,
    }

# All in 2020 $/kg unless otherwise noted, needs to do a thorough check to update values
tea_indices = qs.utils.indices.tea_indices
cost_year = 2020
PCE_indices = tea_indices['PCEPI']
SnowdenSwan_year = 2016
SnowdenSwan_factor = PCE_indices[cost_year]/PCE_indices[SnowdenSwan_year]

Seider_year = 2016
Seider_factor = PCE_indices[cost_year]/PCE_indices[Seider_year]

bst_utility_price = bst.stream_utility_prices
price_dct = {
    'tipping': -39.7/1e3*SnowdenSwan_factor, # PNNL 2022, -$39.7/wet tonne is the weighted average
    'trans_feedstock': 50/1e3*SnowdenSwan_factor, # $50 dry tonne for 78 km, PNNL 32731
    'trans_biocrude': 0.092*SnowdenSwan_factor, # $0.092/GGE of biocrude, 100 miles, PNNL 32731
    'H2': 1.61, # Feng et al., 2024
    'HCcatalyst': 3.52, # Fe-ZSM5, CatCost modified from ZSM5
    'HTcatalyst': 75.18, # Pd/Al2O3, CatCost modified from 2% Pt/TiO2
    'natural_gas': 0.213/0.76**Seider_factor, # $0.213/SCM, $0.76 kg/SCM per https://www.henergy.com/conversion
    'process_water': 0.27/1e3*Seider_factor, # $0.27/m3, higher than $0.80/1,000 gal
    'gasoline': 2.5, # target $/gal
    'jet': 3.53, # 2024$/gal
    'diesel': 3.45, # 2024$/gal, https://afdc.energy.gov/fuels/prices.html
    'N': 0.90, # recovered N in $/kg N
    'P': 1.14, # recovered P in $/kg P
    'K': 0.81, # recovered K in $/kg K
    'solids': -0.17*Seider_factor,
    #!!! Should look at the PNNL report source, at least compare
    'COD': -0.3676, # $/kg; Seider has 0.33 for organics removed
    'wastewater': -0.03/1e3, # $0.03/m3
    }

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

    bst.PowerUtility.price = 0.07*Seider_factor

    #!!! Should set production vs. consumption price    
    # Annual Energy Outlook 2023 https://www.eia.gov/outlooks/aeo/data/browser/
    # Table 8. Electricity Supply, Disposition, Prices, and Emissions
    # End-Use Prices, Industrial, nominal 2024 value in $/kWh