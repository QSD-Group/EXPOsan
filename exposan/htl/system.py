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

HT = su.HT('A410',ins=HTL-2,outs=('HTaqueous','fuel_gas','biooil'))

Acidex = su.AcidExtraction('A200',ins=(HTL-0,acidforP),outs=('residual','extracted'))

M1 = su.HTLmixer('A300',ins=(HTL-1,Acidex-1),outs=('mixture'))

StruPre = su.StruvitePrecipitation('A310',ins=(M1-0,supply_mgcl2),outs=('struvite','CHGfeed'))

CHG = su.CHG('A330',ins=StruPre-1,outs=('fuelgas','effluent'))

MemDis = su.MembraneDistillation('A340',ins=(CHG-1,acidforN),outs=('AmmoniaSulfate','ww'))

sys=qs.System('sys',path=(SludgeThickener,SludgeCentrifuge,HTL,HT,Acidex,M1,StruPre,CHG,MemDis))

sys.simulate()

sys.diagram()

#%%
from qsdsan import Model

model = Model(sys)

from chaospy import distributions as shape

param=model.parameter

dist = shape.Triangle(0.9,0.95,0.98)
@param(name='p_acid_recovery_ratio',
       element=Acidex,
       kind='coupled',
       units='-',
       baseline=0.95,
       distribution=dist)
def set_p_acid_recovery_ratio(i):
    Acidex.P_acid_recovery_rate=i
    
dist = shape.Triangle(0.9,0.95,0.98)
@param(name='p_struvite_recovery_ratio',
       element=StruPre,
       kind='coupled',
       units='-',
       baseline=0.95,
       distribution=dist)
def set_p_struvite_recovery_ratio(i):
    StruPre.P_pre_recovery_rate=i
    
dist = shape.Uniform(0.0235,0.025)
@param(name='sludge_p_ratio',
       element=sludge,
       kind='coupled',
       units='-',
       baseline=0.0235,
       distribution=dist)
def set_sludge_p_ratio(i):
    HTL.sludge_P_ratio=i
    
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
 #%%
fig, ax = qs.stats.plot_uncertainties(model)
fig
#%%
r_df, p_df = qs.stats.get_correlations(model, kind='Spearman')
fig, ax = qs.stats.plot_correlations(r_df)
fig