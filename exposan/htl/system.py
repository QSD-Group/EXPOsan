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
acidforP = qs.WasteStream('acidforP')
supply_mgcl2 = qs.WasteStream('MgCl2')
acidforN = qs.WasteStream('acidforN')


HTL = su.HTL('A120',ins=sludge,outs=('biochar','others'))
S1 = suu.Splitter('A130',ins=HTL-1,outs=('aqueous_biocrude','offgas'),
                  split={
                      'CO2':0,
                      'C_c':1,'N_c':1,'Others_c':1,'C_l':1,'N_l':1,'P_l':1,'Others_l':1,'H2O':1
                      })
S2 = suu.Splitter('A140',ins=S1-0,outs=('aqueous','biocrude'),
                  split={
                      'C_c':0,'N_c':0,'Others_c':0,'C_l':1,'N_l':1,'P_l':1,'Others_l':1,'H2O':1
                      })
HT = su.HT('A330',ins=S2-1,outs=('char', 'others'))
S3 = suu.Splitter('A340',ins=HT-1,outs=('biooil','fuel_gas'),
                  split={
                      'C_o':1,'N_o':1,'Others_o':1,
                      'CO':0,'CO2':0,'CH4':0,'C2H6':0,'C3H8':0,'C5H12':0
                      })
Acidex = su.AcidExtraction('A210',ins=(HTL-0,acidforP),outs=('residual','extracted'))
M1 = su.HTLmixer('A220',ins=(S2-0,Acidex-1),outs=('mixture'))
StruPre = su.StruvitePrecipitation('A230',ins=(M1-0,supply_mgcl2),outs=('struvite','CHGfeed'))
CHG = su.CHG('A250',ins=StruPre-1,outs=('fuelgas','effluent'))
MemDis = su.MembraneDistillation('A260',ins=(CHG-1,acidforN),outs=('AmmoniaSulfate','ww'))

sys = qs.System('sys',path=(HTL,S1,S2,HT,S3,Acidex,M1,StruPre,CHG,MemDis))
    
sys.simulate()
sys.diagram()
