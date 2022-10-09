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

sludge = qs.Stream('sludge',Sludge_lipid=308,Sludge_protein=464,Sludge_carbo=228,
                        H2O=99000, units='kg/hr',T=25+273.15)

sludge_dry_weight = sludge.F_mass-sludge.imass['H2O']
sludge_C_ratio = 0.508
sludge_P_ratio = 0.0235
sludge_H_ratio = 0.0717
sludge_S_ratio = 0.0124
sludge_N_ratio = sludge.imass['Sludge_protein']/sludge_dry_weight/4.78
sludge_O_ratio = 1-sludge_C_ratio-sludge_P_ratio-sludge_H_ratio-sludge_S_ratio-sludge_N_ratio
AOSc = (3*sludge_N_ratio/14+2*sludge_O_ratio/16-sludge_H_ratio/1)/(sludge_C_ratio/12)
sludge_carbo_ratio = sludge.imass['Sludge_carbo']/sludge_dry_weight
sludge_protein_ratio = sludge.imass['Sludge_protein']/sludge_dry_weight


acidforP = qs.Stream('H2SO4_1')
supply_mgcl2 = qs.Stream('MgCl2')
acidforN = qs.Stream('H2S04_2')


SludgeThickener = suu.SludgeThickening('A000',init_with='Stream',ins=sludge,
                                       outs=('Supernatant_1','Compressed_sludge_1'),
                                       solids=('Sludge_lipid','Sludge_protein','Sludge_carbo'))

SludgeCentrifuge = suu.SludgeCentrifuge('A010',init_with='Stream',ins=SludgeThickener-1,
                                        outs=('Supernatant_2','Compressed_sludge_2'),
                                        solids=('Sludge_lipid','Sludge_protein','Sludge_carbo'))

HTL = su.HTL('A110',ins=SludgeCentrifuge-1,outs=('biochar','HTLaqueous','biocrude','offgas'))
HTL.simulate()
biochar_C_ratio = min(1.75*sludge_carbo_ratio,0.65)
biochar_C_N_ratio=15.5
biochar_C_P_ratio=2.163
biochar_N_ratio = biochar_C_ratio/biochar_C_N_ratio
biochar_P_ratio = biochar_C_ratio/biochar_C_P_ratio
biochar_C = HTL.outs[0].F_mass*biochar_C_ratio
biochar_N = HTL.outs[0].F_mass*biochar_N_ratio
biochar_P = min((sludge.F_mass-sludge.imass['H2O'])*sludge_P_ratio,HTL.outs[0].F_mass*biochar_P_ratio)

biocrude_C_ratio = (AOSc*(-8.37)+68.55)/100
biocrude_N_ratio = 0.133*sludge_protein_ratio
biocrude_C = HTL.outs[2].imass['Biocrude']*biocrude_C_ratio
biocrude_N = HTL.outs[2].imass['Biocrude']*biocrude_N_ratio

offgas_C = HTL.outs[3].F_mass*12/44

HTLaqueous_C = sludge_dry_weight*sludge_C_ratio - biochar_C - biocrude_C - offgas_C
HTLaqueous_N = sludge_dry_weight*sludge_N_ratio - biochar_N - biocrude_N
HTLaqueous_P = sludge_dry_weight*sludge_P_ratio - biochar_P

HT = su.HT('A410',ins=HTL-2,outs=('HTaqueous','fuel_gas','biooil'))
HT.simulate()
HTfuel_gas_C = 0
fuelgas_carbo_ratio = {
    'C5H12':60/72,
    'CO':12/28,
    'CO2':12/44,
    'C2H6':24/30,
    'C3H8':36/44,
    'CH4':12/16
     }
for name, ratio in fuelgas_carbo_ratio.items():
    HTfuel_gas_C+=HT.outs[1].imass[name]*ratio
    
biooil_C_ratio = 0.855 
biooil_N_ratio = 0.01
biooil_C = HT.outs[2].F_mass*biooil_C_ratio
biooil_N = HT.outs[2].F_mass*biooil_N_ratio

HTaqueous_C = biocrude_C - HTfuel_gas_C - biooil_C
HTaqueous_N = biocrude_N - biooil_N


Acidex = su.AcidExtraction('A200',ins=(HTL-0,acidforP),outs=('residual','extracted'),
                           biochar_C_ratio=biochar_C_ratio,biochar_C_P_ratio=biochar_C_P_ratio)
Acidex.simulate()
residual_C = biochar_C
redidual_N = biochar_N
residual_P = biochar_P-Acidex.outs[1].imass['P']

extracted_P = Acidex.outs[1].imass['P']


M1 = su.HTLmixer('A300',ins=(HTL-1,Acidex-1),outs=('mixture'),HTLaqueous_C=HTLaqueous_C,
                 HTLaqueous_N=HTLaqueous_N,HTLaqueous_P=HTLaqueous_P)

P_in_struvite=0.127
StruPre = su.StruvitePrecipitation('A310',ins=(M1-0,supply_mgcl2),outs=('struvite','CHGfeed'),P_in_struvite=P_in_struvite)

CHG = su.CHG('A330',ins=StruPre-1,outs=('fuelgas','effluent'))

MemDis = su.MembraneDistillation('A260',ins=(CHG-1,acidforN),outs=('AmmoniaSulfate','ww'))

sys=qs.System('sys',path=(SludgeThickener,SludgeCentrifuge,HTL,HT,Acidex,M1,StruPre,CHG,MemDis))

sys.simulate()

struvite_P=StruPre.outs[0].imass['Struvite']*P_in_struvite
struvite_N=struvite_P*14/31

CHG_C = CHG.outs[1].imass['C']
CHG_N = CHG.outs[1].imass['N']
CHG_P = CHG.outs[1].imass['P']

ww_C = MemDis.outs[1].imass['C']
ww_N = MemDis.outs[1].imass['N']
ww_P = MemDis.outs[1].imass['P']

sys.diagram()

#%%
from qsdsan import Model

model = Model(sys)

from chaospy import distributions as shape

param=model.parameter

dist = shape.Uniform(250,350)
@param(name='sludge_lipid_content',
       element=sludge,
       kind='coupled',
       units='kg/hr',
       baseline=308,
       distribution=dist)
def set_sludge_lipid_content(i):
    sludge.imass['Sludge_lipid']=i
    
dist = shape.Uniform(0.02,0.03)
@param(name='sludge_p_ratio',
       element=sludge,
       kind='coupled',
       units='kg/hr',
       baseline=0.0235,
       distribution=dist)
def set_sludge_p_ratio(i):
    M1.HTLaqueous_P = i
    
metric = model.metric
@metric(name='Struvite_production',units='kg/hr',element='TEA')
def get_struvite_production():
    return StruPre.outs[0].imass['Struvite']

#%%
import numpy as np
np.random.seed(3221)
samples = model.sample(N=100, rule='L')
model.load_samples(samples)
model.evaluate()
model.table