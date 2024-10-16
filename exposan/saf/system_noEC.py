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
# import warnings
# warnings.filterwarnings('ignore')

import os, biosteam as bst, qsdsan as qs
# from biosteam.units import IsenthalpicValve
# from biosteam import settings
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import (
    create_tea,    
    )
from exposan.biobinder import _units as bbu
from exposan.saf import (
    create_components
    # data_path,
    # results_path,
    # _load_components,
    # _load_process_settings,
    # create_tea,
    # _units as u
    )

flowsheet_ID = f'saf'
flowsheet = qs.Flowsheet(flowsheet_ID)
qs.main_flowsheet.set_flowsheet(flowsheet)
saf_cmps = create_components(set_thermo=True)

feedstock = qs.WasteStream('feedstock')
feedstock_water = qs.WasteStream('feedstock_water', Water=1)

FeedstockTrans = bbu.Transportation(
    'FeedstockTrans',
    ins=(feedstock,),
    outs=('transported_feedstock',),
    N_unit=1,
    copy_ins_from_outs=False,
    transportation_distance=78, # km ref [1]
    )

#!!! Need to update the composition (moisture/ash)
moisture = 0.7566
feedstock_composition = {
    'Water': moisture,
    'Lipids': (1-moisture)*0.5315,
    'Proteins': (1-moisture)*0.0255,
    'Carbohydrates': (1-moisture)*0.3816,
    'Ash': (1-moisture)*0.0614,
    }
FeedstockCond = bbu.Conditioning(
    'FeedstockCond', ins=(FeedstockTrans-0, feedstock_water),
    outs='conditioned_feedstock',
    feedstock_composition=feedstock_composition,
    feedstock_dry_flowrate=110*907.185/(24*0.9), # 110 dry sludge tpd ref [1]; 90% upfactor
    N_unit=1,
    )


# =============================================================================
# Hydrothermal Liquefaction (HTL)
# =============================================================================

#!!! Consider adding a pump for feedstock

HX1 = qsu.HXutility('HX1', include_construction=True,
                   ins=FeedstockCond-0, outs='heated_feedstock', T=280+273.15,
                   U=0.0198739, init_with='Stream', rigorous=True)

#!!! Need to update the HTL unit so that it can use actual yields
HTL = qsu.HydrothermalLiquefaction('HTL', ins=HX1-0, outs=('hydrochar','HTL_aqueous','biocrude','offgas_HTL'),
                                   mositure_adjustment_exist_in_the_system=True)

CrudePump = qsu.Pump('CrudePump', ins=HTL-2, outs='press_biocrude', P=1530.0*6894.76,
          init_with='Stream')
# Jones 2014: 1530.0 psia


# Separate water from organics
CrudeDis = qsu.BinaryDistillation('CrudeDis', ins=CrudePump-0,
                        outs=('HT_light','HT_heavy'),
                        LHK=('C4H10','TWOMBUTAN'), P=50*6894.76, # outflow P
                        y_top=188/253, x_bot=53/162, k=2, is_divided=True)

# Separate biocrude from char
CrudeDis = qsu.BinaryDistillation('CrudeDis', ins=CrudePump-0,
                        outs=('HT_light','HT_heavy'),
                        LHK=('C4H10','TWOMBUTAN'), P=50*6894.76, # outflow P
                        y_top=188/253, x_bot=53/162, k=2, is_divided=True)

# =============================================================================
# Hydrocracking
# =============================================================================

include_PSA = False # want to compare with vs. w/o PSA

# External H2, if needed
RSP1 = qsu.ReversedSplitter('RSP1', ins='H2', outs=('HC_H2', 'HT_H2'),
                            init_with='WasteStream')
# reversed splitter, write before HT and HC, simulate after HT and HC
RSP1.ins[0].price = 1.61
RSP1.register_alias('RSP1')

#!!! Need to update the catalyst and price
HC = qsu.Hydrocracking('Hydrocracking', ins=(P3-0, RSP1-1, 'CoMo_alumina_HC'),
                   outs=('HC_out','CoMo_alumina_HC_out'))
HC.ins[2].price = 38.79
HC.register_alias('HC')

CrackedOilPump = qsu.Pump('CrackedOilPump', ins=HTL-2, outs='press_biocrude', P=1530.0*6894.76,
          init_with='Stream')


# =============================================================================
# Hydrotreating
# =============================================================================

#!!! Need to update the catalyst and price
HT = qsu.Hydrotreating('Hydrotreating', ins=(HTL-0, RSP1-0, 'CoMo_alumina_HT'),
           outs=('HTout','CoMo_alumina_HT_out'), include_PSA=include_PSA)
HT.ins[2].price = 38.79
HT.register_alias('HT')

TreatedOilPump = qsu.Pump('TreatedOilPump', ins=HTL-2, outs='press_biocrude', P=1530.0*6894.76,
          init_with='Stream')

# =============================================================================
# Electrochemical Units
# =============================================================================


# =============================================================================
# Products and Wastes
# =============================================================================

# Storage time assumed to be 3 days per [1]
GasolineTank = qsu.StorageTank('T500', ins=PC1-0, outs=('gasoline'),
                                tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')
GasolineTank.outs[0].price = 0.9388

SAFTank = qsu.StorageTank('T510', ins=PC2-0, outs=('diesel'),
                              tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')
# store for 3 days based on Jones 2014
SAFTank.outs[0].price = 0.9722

DieselTank = qsu.StorageTank('T510', ins=PC2-0, outs=('diesel'),
                              tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')
# store for 3 days based on Jones 2014
DieselTank.register_alias('DieselTank')
DieselTank.outs[0].price = 0.9722

# All fuel gases sent to CHP for heat generation
GasMixer = qsu.Mixer('S580', ins=(HTL-3, F1-0, F2-0, D1-0, F3-0),
                      outs=('fuel_gas'), init_with='Stream')
GasMixer.register_alias('GasMixer')

# All wastewater, assumed to be sent to municipal wastewater treatment plant
WWmixer = qsu.Mixer('S590', ins=(SluC-0, MemDis-1, SP2-0),
                    outs='wastewater', init_with='Stream')

# All solids, assumed to be disposed to landfill
SolidsMixer = qsu.Mixer('S590', ins=(SluC-0, MemDis-1, SP2-0),
                    outs='wastewater', init_with='Stream')

# =============================================================================
# Facilities
# =============================================================================

qsu.HeatExchangerNetwork('HXN', T_min_app=86, force_ideal_thermo=True)
# 86 K: Jones et al. PNNL, 2014

CHP = qsu.CombinedHeatPower('CHP', include_construction=True,
                            ins=(GasMixer-0, 'natural_gas', 'air'),
                            outs=('emission','solid_ash'), init_with='WasteStream',
                            supplement_power_utility=False)
CHP.ins[1].price = 0.1685

sys = qs.System.from_units(
    'sys_noEC',
    units=list(flowsheet.unit),
    operating_hours=7920, # 90% uptime
    )

tea = create_tea(sys, IRR_value=0.1, income_tax_value=0.21, finance_interest_value=0.08,
                 labor_cost_value=1.81*10**6)