#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Jianan Feng <jiananf2@illinois.edu>
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs
import exposan.htl._sanunits as su
from qsdsan import sanunits as suu
from exposan.htl._components import create_components


cmps = create_components()

sludge = qs.WasteStream('sludge',Sludge_lipid=308,Sludge_protein=464,Sludge_carbo=228,H2O=4000, units='kg/hr',T=25+273.15)
acidforP = qs.WasteStream('acidforP',H2SO4=33,units='kg/hr')
supply_mgcl2 = qs.WasteStream('MgCl2',MgCl2=66,units='kg/hr')
acidforN = qs.WasteStream('acidforN', H2SO4=114,units='kg/hr')

HX1 = suu.HXutility('HX1',ins=sludge,outs=('HX1out',),T=350+273.15)
HX2 = suu.HXprocess

HTL = su.HTL('HTL',ins=HX1-0,outs=('biocrude','aqueous','offgas','biochar'))
HT = su.HT('HT',ins=HTL-0,outs=('biooil', 'fuelgas_HT', 'char'))
Acidex = su.AcidExtraction('Acidex',ins=(HTL-3,acidforP),outs=('residual','H3PO4'))
M1 = su.HTLmixer('M1',ins=(HTL-1,Acidex-1),outs=('mixture',))
Strupre = su.StruvitePrecipitation('Strupre',ins=(M1-0,supply_mgcl2),outs=('struvite','chgfeed'))
CHG = su.CHG('CHG',ins=Strupre-1,outs=('fuelgas_CHG','chgeffluent'))
MemDis = su.MembraneDistillation('MemDis',ins=(CHG-1,acidforN),outs=('AmmoniaSulfate','ww'))

sys = qs.System('sys',path=(HX1,HTL,HT,Acidex,M1,Strupre,CHG,MemDis))
    
sys.simulate()
sys.diagram()
