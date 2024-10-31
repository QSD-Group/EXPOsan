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
import numpy as np
import matplotlib.pyplot as plt
# from biosteam.units import IsenthalpicValve
# from biosteam import settings
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import create_tea
from exposan.saf import (
    _units as safu,
    price_dct,
    tea_kwargs,
    )
from exposan.biobinder import (
    data_path,
    results_path,
    _load_components,
    _load_process_settings,
    _units as u,
    )


__all__ = ('create_system',)

#!!! Placeholder function for now, update when flowsheet ready
def create_system():
    pass


# %%

# Create and set flowsheet
# (De)centralized HTL, centralized upgrading
configuration = 'CHCU'
flowsheet_ID = f'bb_{configuration}'
flowsheet = qs.Flowsheet(flowsheet_ID)
qs.main_flowsheet.set_flowsheet(flowsheet)

_load_components()
_load_process_settings()


# Desired feedstock flowrate, in dry kg/hr
decentralized_dry_flowrate = 11.46 # feedstock mass flowrate, dry kg/hr
N_decentralized_HTL = 1300 # number of parallel HTL reactor, PNNL is about 1900x of UIUC pilot reactor
# target_HTL_solid_loading = 0.2 # not adjusted


# %%

# =============================================================================
# Hydrothermal Liquefaction
# =============================================================================

scaled_feedstock = qs.WasteStream('scaled_feedstock', price=price_dct['tipping'])
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

FeedstockTrans = safu.Transportation(
    'FeedstockTrans',
    ins=(FeedstockScaler-0, 'feedstock_trans_surrogate'),
    outs=('transported_feedstock',),
    N_unit=N_decentralized_HTL,
    copy_ins_from_outs=True,
    transportation_distance=25, # km ref [1]
    )
#!!! where was 64.1 from?
# trans_unit_cost = 64.1/1e3 * PCE_indices[cost_year]/PCE_indices[2016] # $/kg/km PNNL 32731
# FeedstockTrans.transportation_unit_cost = BiocrudeTrans.transportation_unit_cost = trans_unit_cost

FeedstockCond = safu.Conditioning(
    'FeedstockCond', ins=(FeedstockTrans-0, ProcessWaterScaler-0),
    outs='conditioned_feedstock',
    feedstock_composition=safu.salad_dressing_waste_composition,
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
    target_moisture=0.001, #!!! so that FracDist can work
    N_unit=N_decentralized_HTL,)

BiocrudeWaterScaler = u.Scaler(
    'BiocrudeWaterScaler', ins=BiocrudeDewatering-1, outs='scaled_biocrude_water',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )

BiocrudeTrans = safu.Transportation(
    'BiocrudeTrans',
    ins=(BiocrudeDewatering-0, 'biocrude_trans_surrogate'),
    outs=('transported_biocrude',),
    N_unit=N_decentralized_HTL,
    transportation_distance=1, # cost considered average transportation distance
    # transportation_distance=biocrude_transportation_distance, # km ref [1]
    )

BiocrudeScaler = u.Scaler(
    'BiocrudeScaler', ins=BiocrudeTrans-0, outs='scaled_biocrude',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )

BiocrudeSplitter = safu.BiocrudeSplitter(
    'BiocrudeSplitter', ins=BiocrudeScaler-0, outs='splitted_biocrude',
    cutoff_Tb=343+273.15, light_frac=0.5316)

# Shortcut column uses the Fenske-Underwood-Gilliland method,
# better for hydrocarbons according to the tutorial
# https://biosteam.readthedocs.io/en/latest/API/units/distillation.html
FracDist = qsu.ShortcutColumn(
    'FracDist', ins=BiocrudeSplitter-0,
    outs=('biocrude_light','biocrude_heavy'),
    LHK=('Biofuel', 'Biobinder'), # will be updated later
    P=50*6894.76, # outflow P, 50 psig
    # Lr=0.1, Hr=0.5,
    y_top=188/253, x_bot=53/162,
    k=2, is_divided=True)
@FracDist.add_specification
def adjust_LHK():
    FracDist.LHK = BiocrudeSplitter.keys[0]
    FracDist._run()
    
# Lr_range = Hr_range = np.linspace(0.05, 0.95, 19)
# Lr_range = Hr_range = np.linspace(0.01, 0.2, 20)
# results = find_Lr_Hr(FracDist, Lr_trial_range=Lr_range, Hr_trial_range=Hr_range)
# results_df, Lr, Hr = results

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

ElectrochemicalOxidation = u.ElectrochemicalOxidation(
    'ElectrochemicalCell',
    ins=(HTL-1,),
    outs=('fertilizer', 'recycled_water', 'filtered_solids'),
    N_unit=N_decentralized_HTL,)

FertilizerScaler = u.Scaler(
    'FertilizerScaler', ins=ElectrochemicalOxidation-0, outs='scaled_fertilizer',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )

RecycledWaterScaler = u.Scaler(
    'RecycledWaterScaler', ins=ElectrochemicalOxidation-1, outs='scaled_recycled_water',
    scaling_factor=N_decentralized_HTL, reverse=False,
    )

FilteredSolidsScaler = u.Scaler(
    'FilteredSolidsScaler', ins=ElectrochemicalOxidation-2, outs='filterd_solids',
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
PWC = safu.ProcessWaterCenter(
    'ProcessWaterCenter',
    process_water_streams=[ProcessWaterScaler.ins[0]],
    process_water_price=price_dct['process_water']
    )
PWC.register_alias('PWC')

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
# System, TEA, and LCA
# =============================================================================

sys = qs.System.from_units(
    f'sys_{configuration}',
    units=list(flowsheet.unit), 
    operating_hours=7920, # same as the HTL module, about 90% uptime
    )
sys.register_alias('sys')
stream = sys.flowsheet.stream

#!!! burn hydrochar
hydrochar = stream.hydrochar
hydrochar.price = 0

base_labor = 338256 # for 1000 kg/hr #!!! where was this from?
tea_kwargs['labor_cost'] = lambda: (scaled_feedstock.F_mass-scaled_feedstock.imass['Water'])/1000*base_labor


biofuel_additives = stream.biofuel_additives
biofuel_additives.price = price_dct['diesel']

# https://idot.illinois.gov/doing-business/procurements/construction-services/transportation-bulletin/price-indices.html
biobinder_price = 0.67 # bitumnous, IL
tea = create_tea(sys, **tea_kwargs)


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
biobinder = stream.biobinder
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

lifetime = tea.duration[1]

#!!! Need to get heating duty
lca = qs.LCA(
    system=sys,
    lifetime=lifetime,
    uptime_ratio=sys.operating_hours/(365*24),
    Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
    # Heating=lambda:sys.get_heating_duty()/1000*lifetime,
    Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
    )



def simulate_and_print(save_report=False):
    sys.simulate()
    
    
    biobinder.price = biobinder_price = tea.solve_price(biobinder)
    print(f'Minimum selling price of the biobinder is ${biobinder_price:.2f}/kg.')
    c = qs.currency
    for attr in ('NPV','AOC', 'sales', 'net_earnings'):
        uom = c if attr in ('NPV', 'CAPEX') else (c+('/yr'))
        print(f'{attr} is {getattr(tea, attr):,.0f} {uom}')
        
    all_impacts = lca.get_allocated_impacts(streams=(biobinder,), operation_only=True, annual=True)
    GWP = all_impacts['GlobalWarming']/(biobinder.F_mass*lca.system.operating_hours)
    print(f'Global warming potential of the biobinder is {GWP:.4f} kg CO2e/kg.')
    if save_report:
        # Use `results_path` and the `join` func can make sure the path works for all users
        sys.save_report(file=os.path.join(results_path, 'sys.xlsx'))

if __name__ == '__main__':
    simulate_and_print()

