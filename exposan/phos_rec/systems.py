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

fe_sludge = qs.WasteStream('sludge', Fe3=180, Org=5000, PO4=300, Water=1000000,
                           Ca2=150, Mg2=100, Inert=1000, units='kg/d')

# TODO: check mixing power
AF = su.AcidogenicFermenter(ID='AF', ins=(fe_sludge,'food_waste'), outs=('gas', 'fermentate'),
                            food_sludge_ratio=1)

FC = qsu.SludgeCentrifuge(ID='FC', ins=AF-1, outs=('fermentation_supernatant', 'residue'), 
                          sludge_moisture=0.85, solids=('Inert','Residue'))

# TODO: check mixing power
SP = su.SelectivePrecipitation(ID='SP', ins=(FC-0, 'acid', 'oxidant'), outs='slurry',
                               P_recovery=0.82)

PC = qsu.SludgeCentrifuge(ID='PC', ins=SP-0, outs=('precipitation_supernatant', 'precipitate'), 
                          sludge_moisture=0.9, solids=('FePO4_2H2O',))

HD = su.HeatDrying(ID='HD', ins=(PC-0, 'natural_gas'), outs=('dried_precipitate', 'vapor'),
                   target_moisture=0.2, T_out=90 + _C_to_K)

sys = qs.System.from_units('phos_rec', units=[AF, FC, SP, PC, HD])
sys.simulate()
sys.diagram()