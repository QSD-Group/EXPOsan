# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, qsdsan as qs
from qsdsan import sanunits as qsu
from biosteam.units import IsenthalpicValve
from qsdsan.utils import clear_lca_registries
from exposan.biobinder import (
    _load_components,
    _load_process_settings,
    create_tea,
    _units as u
    )
from biosteam import settings

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
    'feedstock', Lipids=62.45, Proteins=2.38, Carbohydrates=29.46, Ash=5.71,
    Water=100/24.34*75.66,
    )
feedstock.F_mass = 11.46 # all component flowrates will be adjsuted accordingly

# Adjust feedstock moisture
feedstock_water = qs.WasteStream('feedstock_water')
T101 = qsu.MixTank('T101', ins=(feedstock, feedstock_water))
@T101.add_specification
def adjust_feedstock_water():
    feedstock_water.imass['Water'] = max(0, (feedstock.F_mass-feedstock.imass['Water'])/0.2-feedstock.imass['Water'])
    T101._run()

HTL = u.PilotHTL(
    'A102', ins=T101-0, outs=('hydrochar','HTL_aqueous','biocrude','offgas_HTL'))
HTL.register_alias('HTL')

# %%

# =============================================================================
# Area 200 Aqueous Product Treatment
# =============================================================================

HTLaq_Tank = qsu.StorageTank(
    'T200', ins=HTL-1, outs=('treated_aq'),
    init_with='WasteStream', tau=24*7, vessel_material='Stainless steel')


# %%

# =============================================================================
# Area 300 Biocrude Upgrading
# =============================================================================

#!!! Need to connect to the biocrude product
D1 = qsu.BinaryDistillation('A370', ins=H3-0,
                        outs=('HT_light','HT_heavy'),
                        LHK=('C4H10','TWOMBUTAN'), P=50*6894.76, # outflow P
                        y_top=188/253, x_bot=53/162, k=2, is_divided=True)
D1.register_alias('D1')


# %%

# =============================================================================
# Assemble System, TEA, LCA
# =============================================================================

sys = qs.System.from_units(
    f'sys_{configuration}',
    units=list(flowsheet.unit), 
    )
sys.register_alias('sys')

# sys.diagram() # see a diagram of the system
# Currently won't work since there are many conversion factors not included
sys.simulate()
