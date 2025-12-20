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
import qsdsan as qs
import qsdsan.sanunits as su
from exposan.phos_rec._components import create_components
from exposan.phos_rec._sanunits import AcidogenicFermenter, SelectivePrecipitation

create_components()

fe_sludge = qs.WasteStream('sludge', Fe3=180, Org=5000, PO4=300, Water=1000000,
                           Ca2=150, Mg2=100, Inert=1000, units='kg/d')

AF = AcidogenicFermenter(ID='AF', ins=(fe_sludge,'food_waste'), outs=('gas', 'fermentate'),
                         food_sludge_ratio=1)

FC = su.SludgeCentrifuge(ID='FC', ins=AF-1, outs=('fermentation_supernatant', 'residue'), 
                         sludge_moisture=0.85, solids=('Inert','Residue'))

SP = SelectivePrecipitation(ID='SP', ins=(FC-0, 'acid', 'oxidant'), outs='slurry',
                            P_recovery=0.82)

PC = su.SludgeCentrifuge(ID='PC', ins=SP-0, outs=('precipitation_supernatant', 'precipitate'), 
                         sludge_moisture=0.9, solids=('FePO4_2H2O',))

sys = qs.System.from_units('phos_rec',units=[AF, FC, SP, PC])
sys.simulate()
sys.diagram()

