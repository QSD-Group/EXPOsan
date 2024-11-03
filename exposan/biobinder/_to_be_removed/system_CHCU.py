# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 14:23:30 2024

@author: aliah
"""

import os
import biosteam as bst
import qsdsan as qs
import numpy as np
from qsdsan import sanunits as qsu
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

# Placeholder function for now
def create_system():
    pass

# Create and set flowsheet
configuration = 'CC'  # Centralized HTL, centralized upgrading
flowsheet_ID = f'biobinder_{configuration}'
flowsheet = qs.Flowsheet(flowsheet_ID)
qs.main_flowsheet.set_flowsheet(flowsheet)

_load_components()
_load_process_settings()

# Desired feedstock flowrate, in dry kg/hr
centralized_dry_flowrate = 1000  # can be updated

# Centralized Hydrothermal Liquefaction
scaled_feedstock = qs.WasteStream('scaled_feedstock')

FeedstockCond = u.Conditioning(
    'FeedstockCond', ins=(scaled_feedstock, 'fresh_process_water'),
    outs='conditioned_feedstock',
    feedstock_composition=u.salad_dressing_waste_composition,
    feedstock_dry_flowrate=centralized_dry_flowrate,
)

HTL = u.CentralizedHTL(
    'HTL', ins=FeedstockCond-0, outs=('hydrochar', 'HTL_aqueous', 'biocrude', 'HTL_offgas'),
    afdw_yields=u.salad_dressing_waste_yields,
    N_unit=1,  # Single centralized reactor
)

# Centralized Biocrude Upgrading
BiocrudeDeashing = u.BiocrudeDeashing(
    'BiocrudeDeashing', ins=HTL-2, outs=('deashed_biocrude', 'biocrude_ash'),
    N_unit=1,
)

BiocrudeDewatering = u.BiocrudeDewatering(
    'BiocrudeDewatering', ins=BiocrudeDeashing-0, outs=('dewatered_biocrude', 'biocrude_water'),
    N_unit=1,
)

FracDist = u.ShortcutColumn(
    'FracDist', ins=BiocrudeDewatering-0,
    outs=('biocrude_light', 'biocrude_heavy'),
    LHK=('Biofuel', 'Biobinder'),
    P=50 * 6894.76,
    y_top=188/253, x_bot=53/162,
    k=2, is_divided=True
)

@FracDist.add_specification
def adjust_LHK():
    FracDist.LHK = (BiocrudeDeashing.light_key, BiocrudeDeashing.heavy_key)
    FracDist._run()

LightFracStorage = qsu.StorageTank(
    'LightFracStorage',
    FracDist-0, outs='biofuel_additives',
    tau=24*7, vessel_material='Stainless steel'
)

HeavyFracStorage = qsu.StorageTank(
    'HeavyFracStorage', FracDist-1, outs='biobinder',
    tau=24*7, vessel_material='Stainless steel'
)

# Aqueous Product Treatment
ElectrochemicalOxidation = u.ElectrochemicalOxidation(
    'MicrobialFuelCell',
    ins=(HTL-1,),
    outs=('fertilizer', 'recycled_water', 'filtered_solids'),
    N_unit=1,
)

# Facilities and waste disposal
AshDisposal = u.Disposal(
    'AshDisposal', 
    ins=(BiocrudeDeashing-1, 'filtered_solids'),
    outs=('ash_disposal', 'ash_others'),
    exclude_components=('Water',)
)

WWDisposal = u.Disposal(
    'WWDisposal', 
    ins='biocrude_water',
    outs=('ww_disposal', 'ww_others'),
    exclude_components=('Water',)
)

# Heat exchanger network
HXN = qsu.HeatExchangerNetwork('HXN', T_min_app=86, force_ideal_thermo=True)

# Assemble System
sys = qs.System.from_units(
    f'sys_{configuration}',
    units=list(flowsheet.unit),
    operating_hours=7920,  # 90% uptime
)

sys.register_alias('sys')
stream = sys.flowsheet.stream


# ...

def simulate_and_print(save_report=False):
    sys.simulate()
    
    biobinder.price = biobinder_price = tea.solve_price(biobinder)
    print(f'Minimum selling price of the biobinder is ${biobinder_price:.2f}/kg.')
    c = qs.currency
    for attr in ('NPV', 'AOC', 'sales', 'net_earnings'):
        uom = c if attr in ('NPV', 'CAPEX') else (c + ('/yr'))
        print(f'{attr} is {getattr(tea, attr):,.0f} {uom}')
        
    all_impacts = lca.get_allocated_impacts(streams=(biobinder,), operation_only=True, annual=True)
    GWP = all_impacts['GlobalWarming'] / (biobinder.F_mass * lca.system.operating_hours)
    print(f'Global warming potential of the biobinder is {GWP:.4f} kg CO2e/kg.')
    
    if save_report:
        sys.save_report(file=os.path.join(results_path, 'centralized_sys.xlsx'))
