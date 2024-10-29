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
    '_load_process_settings',
    'uptime_ratio',
    'dry_flowrate',
    'feedstock_composition',
    'wet_flowrate',
    'HTL_yields',
    'annual_hours',
    'price_dct',
    )


tpd = 110 # dry mass basis
uptime_ratio = 0.9

dry_flowrate = tpd*907.185/(24*uptime_ratio) # 110 dry sludge tpd [1]

#!!! Need to update the composition (moisture/ash)
moisture = 0.7566
feedstock_composition = {
    'Water': moisture,
    'Lipids': (1-moisture)*0.5315,
    'Proteins': (1-moisture)*0.0255,
    'Carbohydrates': (1-moisture)*0.3816,
    'Ash': (1-moisture)*0.0614,
    }
wet_flowrate = dry_flowrate / (1-moisture)

HTL_yields = {
    'gas': 0.006,
    'aqueous': 0.192,
    'biocrude': 0.802,
    'char': 0,
    }

annual_hours = 365*24*uptime_ratio

# All in 2020 $/kg unless otherwise noted, needs to do a thorough check to update values
bst_utility_price = bst.stream_utility_prices
price_dct = {
    'tipping': -69.14/907.185, # tipping fee 69.14±21.14 for IL, https://erefdn.org/analyzing-municipal-solid-waste-landfill-tipping-fees/
    'transportation': 50/1e3, # $50 kg/tonne for 78 km, 2016$ ref [1]
    'H2': 1.61, # Feng et al.
    'HCcatalyst': 3.52, # Fe-ZSM5, CatCost modified from ZSM5
    'HTcatalyst': 75.18, # Pd/Al2O3, CatCost modified from 2% Pt/TiO2
    'natural_gas': 0.1685,
    'process_water': bst_utility_price['Process water'],
    'gasoline': 2.5, # target $/gal
    'jet': 3.53, # 2024$/gal
    'diesel': 3.45, # 2024$/gal
    'N': 0.90, # recovered N in $/kg N
    'P': 1.14, # recovered P in $/kg P
    'K': 0.81, # recovered K in $/kg K
    'solids': bst_utility_price['Ash disposal'],
    'COD': -0.3676, # $/kg
    'wastewater': -0.03/1e3, # $0.03/m3
    }

def _load_process_settings():
    bst.CE = qs.CEPCI_by_year[2020]
    bst.PowerUtility.price = 0.06879 
    
    # # These utilities are provided by CHP thus cost already considered
    # # setting the regeneration price to 0 or not will not affect the final results
    # # as the utility cost will be positive for the unit that consumes it
    # # but negative for HXN/CHP as they produce it
    # for adj in ('low', 'medium', 'high'):
    #     steam = bst.HeatUtility.get_agent(f'{adj}_pressure_steam')
    #     steam.heat_transfer_price = steam.regeneration_price = 0.
