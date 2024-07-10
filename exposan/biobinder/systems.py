# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# !!! Temporarily ignoring warnings
import warnings
warnings.filterwarnings('ignore')

import os, math, biosteam as bst, qsdsan as qs
# from biosteam.units import IsenthalpicValve
# from biosteam import settings
from qsdsan import sanunits as qsu
# from qsdsan.utils import clear_lca_registries
from exposan.biobinder import (
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

#!!! PAUSED
# Have two separate ReverseSplitters to fix feedstock tipping fee/process water
# also need to consider recycling for the process water


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
centralized_dry_flowrate = decentralized_dry_flowrate*1000 # PNNL is about 1900x of UIUC pilot reactor

# Salad dressing waste, all on weight basis
feedstock_composition = {
    'Water': 0.7566,
    'Lipids': 0.2434*0.6245,
    'Proteins': 0.2434*0.0238,
    'Carbohydrates': 0.2434*0.2946,
    'Ash': 0.2434*0.0571,
    }

target_HTL_solid_loading = 0.2

# %%

# =============================================================================
# Area 100 Hydrothermal Liquefaction
# =============================================================================

feedstock = qs.WasteStream('feedstock')
htl_process_water = qs.WasteStream('htl_process_water')

# Adjust feedstock moisture
FeedstockPrep = T101 = u.PreProcessing(
    'T101', ins=(feedstock, htl_process_water),
    decentralized_dry_flowrate=decentralized_dry_flowrate,
    centralized_dry_flowrate=centralized_dry_flowrate,    
    )

HTL = u.PilotHTL(
    'R102', ins=T101-0, outs=('hydrochar','HTL_aqueous','biocrude','HTL_offgas'),
    feedstock_composition=u.salad_dressing_composition,
    decentralized_dry_flowrate=decentralized_dry_flowrate,
    centralized_dry_flowrate=centralized_dry_flowrate,
    ) 
HTL.register_alias('HTL')


# %%

# =============================================================================
# Area 200 Aqueous Product Treatment
# =============================================================================

AqueousFiltration = u.AqueousFiltration(
    'S201', ins=HTL-1, outs=('treated_aq'), init_with='WasteStream')

AqStorage = qsu.StorageTank(
    'T202', ins=AqueousFiltration-0, outs=('stored_aq'),
    init_with='WasteStream', tau=24*7, vessel_material='Stainless steel')


# %%

# =============================================================================
# Area 300 Biocrude Upgrading
# =============================================================================

Deashing = u.BiocrudeDeashing('A301', ins=HTL-2, outs=('deashed', 'biocrude_ash'))
Deashing.register_alias('Deashing')
Dewatering = u.BiocrudeDewatering('A302', ins=Deashing-0, outs=('dewatered', 'biocrude_water'))
Dewatering.register_alias('Dewatering')
BiocrudeTrans = u.Transportation('U301', ins=Dewatering-0, outs='transported_biocrude')

BiocrudeSplitter = u.BiocrudeSplitter('S303', ins=Dewatering-0,
                                      cutoff_Tb=343+273.15, light_frac=0.5316)
BiocrudeSplitter.register_alias('BiocrudeSplitter')

# Shortcut column uses the Fenske-Underwood-Gilliland method,
# better for hydrocarbons according to the tutorial
# https://biosteam.readthedocs.io/en/latest/API/units/distillation.html
FracDist = u.ShortcutColumn('D304', ins=BiocrudeSplitter-0,
                        outs=('biocrude_light','biocrude_heavy'),
                        LHK=('Biofuel', 'Biobinder'), # will be updated later
                        P=50*6894.76, # outflow P, 50 psig
                        y_top=188/253, x_bot=53/162, k=2, is_divided=True)
FracDist.register_alias('FracDist')
@FracDist.add_specification
def adjust_LHK():
    FracDist.LHK = (BiocrudeSplitter.light_key, BiocrudeSplitter.heavy_key)
    FracDist._run()


# FracDist = qsu.BinaryDistillation('D304', ins=BiocrudeSplitter-0,
#                         outs=('biocrude_light','biocrude_heavy'),
#                         LHK=('7MINDOLE', 'C16:0FA'), # will be updated later
#                         # P=50*6894.76, # outflow P
#                         # P=101325*5, # outflow P
#                         # Lr=0.1, Hr=0.5,
#                         y_top=0.1134, x_bot=0.0136,
#                         # y_top=188/253, x_bot=53/162,
#                         k=2, is_divided=True)
# FracDist.register_alias('FracDist')
# @FracDist.add_specification
# def adjust_LHK():
#     FracDist.LHK = (BiocrudeSplitter.light_key, BiocrudeSplitter.heavy_key)
#     FracDist._run()

LightFracStorage = qsu.StorageTank('T305', FracDist-0, outs='biofuel_additives',
                                   tau=24*7, vessel_material='Stainless steel')
LightFracStorage.register_alias('LightFracStorage')
HeavyFracStorage = qsu.StorageTank('T306', FracDist-1, outs='biobinder',
                                   tau=24*7, vessel_material='Stainless steel')
HeavyFracStorage.register_alias('HeavyFracStorage')

# %%

# =============================================================================
# Area 400 Waste disposal
# =============================================================================

AshDisposal = u.Disposal('U401', ins=Deashing-1,
                         outs=('ash_disposal', 'ash_others'),
                         exclude_components=('Water',))
WaterDisposal = u.Disposal('U402', ins=Dewatering-1,
                           outs=('water_disposal', 'water_others'),
                           exclude_components=('Water', 'Ash'))


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
GDP_indices = {
    2003: 0.808,
    2005: 0.867,
    2007: 0.913,
    2008: 0.941,
    2009: 0.951,
    2010: 0.962,
    2011: 0.983,
    2012: 1.000,
    2013: 1.014,
    2014: 1.033,
    2015: 1.046,
    2016: 1.059,
    2017: 1.078,
    2018: 1.100,
    2019: 1.123,
    2020: 1.133,
    2021: 1.181,
    2022: 1.269,
    2023: 1.322,
    2024: 1.354,
    }

# Inputs
feedstock.price = -69.14/907.185 # tipping fee 69.14±21.14 for IL, https://erefdn.org/analyzing-municipal-solid-waste-landfill-tipping-fees/

# Utilities, price from Table 17.1 in Seider et al., 2016$
# Use bst.HeatUtility.cooling_agents/heating_agents to see all the heat utilities
Seider_factor = GDP_indices[cost_year]/GDP_indices[2016]

htl_process_water.price = 0.8/1e3/3.758*Seider_factor # process water for moisture adjustment

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
WaterDisposal.disposal_price = 0.17*Seider_factor # dewater, for organics removed

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

base_labor = 338256 # for 1000 kg/hr

tea = create_tea(
    sys,
    labor_cost=lambda: (feedstock.F_mass-feedstock.imass['Water'])/1000*base_labor,
    )

# To see out-of-boundary-limits units
# tea.OSBL_units

# =============================================================================
# LCA
# =============================================================================


# %%

def simulate_and_print(save_report=False):
    sys.simulate()
    
    biobinder.price = biobinder_price = tea.solve_price(biobinder)
    print(f'Minimum selling price of the biobinder is ${biobinder_price:.2f}/kg.')
    c = qs.currency
    for attr in ('NPV','AOC', 'sales', 'net_earnings'):
        uom = c if attr in ('NPV', 'CAPEX') else (c+('/yr'))
        print(f'{attr} is {getattr(tea, attr):,.0f} {uom}')
    if save_report:
        # Use `results_path` and the `join` func can make sure the path works for all users
        sys.save_report(file=os.path.join(results_path, 'sys.xlsx'))

if __name__ == '__main__':
    simulate_and_print()
    