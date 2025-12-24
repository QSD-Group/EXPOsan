# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Xuan Wang <easonbiubiu99@gmail.com>
    
    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

import qsdsan as qs, qsdsan.sanunits as qsu, exposan.phos_rec._sanunits as su
from exposan.phos_rec._components import create_components

_C_to_K = 273.15

create_components()

# Inert represents other inorganics; now VSS/TSS = 0.74, which is reasonable
fe_sludge = qs.WasteStream(ID='fe_sludge', Fe3=180, Org=5000, PO4=300, Water=1000000,
                           Ca2=150, Mg2=100, Inert=1000, units='kg/d')

P1 = qsu.SludgePump(ID='P1', ins=fe_sludge, outs='sludge_to_ST')

# the same as AF: 4 days of reaction+ 3 hr for cleaning and unloading
ST = qsu.StorageTank(ID='ST', ins=P1-0, outs='sludge_to_P2', tau=4*24+3)

P2 = qsu.SludgePump(ID='P2', ins=ST-0, outs='sludge_to_AF')

AF = su.AcidogenicFermenter(ID='AF', ins=(P2-0, 'food_waste'), outs=('fermentate', 'fermentation_gas'),
                            food_sludge_ratio=1)

FC = qsu.SludgeCentrifuge(ID='FC', ins=AF-0, outs=('fermentation_supernatant', 'residue'),
                          sludge_moisture=0.85, solids=('Inert','Residue'))

SP = su.SelectivePrecipitation(ID='SP', ins=(FC-0, 'acid', 'oxidant'), outs='slurry',
                               P_recovery=0.82)

# TODO: Org seems low, but conservative: less organics for heat generation, and more go back to WRRF headworks
PC = qsu.SludgeCentrifuge(ID='PC', ins=SP-0, outs=('precipitation_supernatant', 'precipitate'), 
                          sludge_moisture=0.92, solids=('FePO4_2H2O',))

HD = su.HeatDrying(ID='HD', ins=(PC-1, 'heat_drying_natural_gas'), outs=('dried_precipitate', 'heat_drying_vapor'),
                   T=105+_C_to_K)

SI = su.Sintering(ID='ST', ins=(HD-0, 'sintering_natural_gas', 'air'), outs=('product', 'sintering_vapor'))

sys = qs.System.from_units('phos_rec', units=[P1, ST, P2, AF, FC, SP, PC, HD, SI])
sys.simulate()
sys.diagram()