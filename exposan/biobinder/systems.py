# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

#!!! Temporarily ignoring warnings
# import warnings
# warnings.filterwarnings('ignore')

import os, qsdsan as qs
from biosteam.units import IsenthalpicValve
from biosteam import settings
from exposan.htl import create_tea
# from biorefineries.tea import create_cellulosic_ethanol_tea
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.biobinder import (
    _load_components,
    _load_process_settings,
    create_tea,
    _units as u
    )


# Create and set flowsheet
configuration = 'pilot'
flowsheet_ID = f'biobinder_{configuration}'
flowsheet = qs.Flowsheet(flowsheet_ID)
qs.main_flowsheet.set_flowsheet(flowsheet)

_load_components()
_load_process_settings()

__all__ = ('create_system',)

#!!! Placeholder function for now, update when flowsheet ready
def create_system():
    pass

# %%

# =============================================================================
# Area 100 Hydrothermal Liquefaction
# =============================================================================

feedstock = qs.WasteStream(
    'feedstock', Lipids=62.45*24.34, Proteins=2.38*24.34, Carbohydrates=29.46*24.34, Ash=5.71*24.34,
    Water=75.66*100,
    )
feed_factor= 0.93*0.2 #SDW feedstock desnity 0.93 g/ml, solid content=20 wt%
feedstock.F_mass = 537.63*feed_factor  #dry feedstock flow rate kg/hr
#all component flowrates will be adjsuted accordingly

# Adjust feedstock moisture
feedstock_water = qs.WasteStream('feedstock_water')
T101 = qsu.MixTank('T101', ins=(feedstock, feedstock_water))
@T101.add_specification
def adjust_feedstock_water():
    feedstock_water.imass['Water'] = max(0, (feedstock.F_mass-feedstock.imass['Water'])/0.2-feedstock.imass['Water'])
    T101._run()

HTL = u.PilotHTL(
    'R102', ins=T101-0, outs=('hydrochar','HTL_aqueous','biocrude','HTL_offgas'))
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


Deashing = u.BiocrudeDeashing('A301', HTL-2, outs=('deashed', 'excess_ash'))
Deashing.register_alias('Deashing')
Dewatering = u.BiocrudeDewatering('A302', Deashing-0, outs=('dewatered', 'excess_water'))
Dewatering.register_alias('Dewatering')
BiocrudeSplitter = u.BiocrudeSplitter('S303', ins=Dewatering-0,
                                      cutoff_Tb=343+273.15, light_frac=0.5316)
BiocrudeSplitter.register_alias('BiocrudeSplitter')

# Shortcut column uses the Fenske-Underwood-Gilliland method,
# better for hydrocarbons according to the tutorial
# https://biosteam.readthedocs.io/en/latest/API/units/distillation.html
FracDist = u.ShortcutColumn('D304', ins=BiocrudeSplitter-0,
                        outs=('biocrude_light','biocrude_heavy'),
                        LHK=('Biofuel', 'Biobinder'), # will be updated later
                        P=50*6894.76, # outflow P
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
# Assemble System, TEA, LCA
# =============================================================================

sys = qs.System.from_units(
    f'sys_{configuration}',
    units=list(flowsheet.unit), 
    operating_hours=7920, # same as the HTL module, about 90% updtime
    )
sys.register_alias('sys')

stream = sys.flowsheet.stream
biofuel_additives = stream.biofuel_additives
biofuel_additives.price = 1.4 # $/kg

biobinder = stream.biobinder

base_labor = 338256 # for 1000 kg/hr


tea = create_tea(
    sys,
    labor_cost_value=(feedstock.F_mass-feedstock.imass['Water'])/1000*base_labor,
    )
sys.simulate()

#%%

biobinder.price = biobinder_price = tea.solve_price(biobinder)
print(f'Minimum selling price of the biobinder is ${biobinder_price:.2f}/kg.')
c = qs.currency
for attr in ('NPV','AOC', 'sales', 'net_earnings'):
    uom = c if attr in ('NPV', 'CAPEX') else (c+('/yr'))
    print(f'{attr} is {getattr(tea, attr):,.0f} {uom}')
#sys.save_report(file= 'C:\Work\Rutgers\Biobinder\sys.xlsx')
