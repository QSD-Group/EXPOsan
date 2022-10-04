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

sludge = qs.WasteStream('sludge',Sludge_lipid=308,Sludge_protein=464,Sludge_carbo=228,units='kg/h')
acidforP = qs.WasteStream('acidforP',H2SO4=33,units='kg/h')

htl = su.HTL('htl',ins=sludge,outs=('biocrude','aqueous','offgas','biochar'))
ht = su.HT('ht',ins=htl-0,outs=('biooil', 'fuelgas_HT', 'char'))
acidex = su.AcidExtraction('acidex',ins=(acidforP,htl-3),outs=('residual','H3PO4'))


sys = qs.System('sys',path=(htl,ht,acidex))
    
sys.simulate()
sys.diagram()
