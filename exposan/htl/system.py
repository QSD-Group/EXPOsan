#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 15:31:22 2022

@author: yalinli_cabbi
"""

import qsdsan as qs
import exposan.htl._sanunits as su
from exposan.htl._components import create_components


cmps = create_components()

sludge = qs.WasteStream('sludge',Sludge_lipid=308,Sludge_protein=464,Sludge_carbo=228,units='kg/hr')
acidforP = qs.WasteStream('acidforP',H2SO4=33,units='kg/hr')
supply_mgcl2 = qs.WasteStream('MgCl2',MgCl2=66,units='kg/hr')
acidforN = qs.WasteStream('acidforN', H2SO4=114,units='kg/hr')

htl = su.HTL('htl',ins=sludge,outs=('biocrude','aqueous','offgas','biochar'))
ht = su.HT('ht',ins=htl-0,outs=('biooil', 'fuelgas_HT', 'char'))
acidex = su.AcidExtraction('acidex',ins=(htl-3,acidforP),outs=('residual','H3PO4'))
M1 = su.HTLmixer('M1',ins=(htl-1,acidex-1),outs=('mixture',))
strupre = su.StruvitePrecipitation('strupre',ins=(M1-0,supply_mgcl2),outs=('struvite','chgfeed'))
chg = su.CHG('chg',ins=strupre-1,outs=('fuelgas_CHG','chgeffluent'))
memdis = su.MembraneDistillation('memdis',ins=(chg-1,acidforN),outs=('AmmoniaSulfate','ww'))

sys = qs.System('sys',path=(htl,ht,acidex,M1,strupre,chg,memdis))
    
sys.simulate()
sys.diagram()
