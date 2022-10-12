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

# __all__ = ('create_system',)

# def create_system():

cmps = create_components()

fake_sludge = qs.Stream('fake_sludge',H2O=100000,units='kg/hr',T=298.15)
#set H2O equal to the total sludge input flow
acidforP = qs.Stream('H2SO4_1')
supply_mgcl2 = qs.Stream('MgCl2')
acidforN = qs.Stream('H2SO4_2')

SluL = su.SludgeLab('S000',ins=fake_sludge,outs='real_sludge')

SluT = suu.SludgeThickening('A000',ins=SluL-0,outs=('Supernatant_1','Compressed_sludge_1'),
                            init_with='Stream',solids=('Sludge_lipid','Sludge_protein','Sludge_carbo',
                                                       'Sludge_ash'))

SluC = suu.SludgeCentrifuge('A010',ins=SluT-1,outs=('Supernatant_2','Compressed_sludge_2'),
                            init_with='Stream',solids=('Sludge_lipid','Sludge_protein','Sludge_carbo',
                                                       'Sludge_ash'))

HTL = su.HTL('A110',ins=SluC-1,outs=('biochar','HTLaqueous','biocrude','offgas'))

HT = su.HT('A410',ins=HTL-2,outs=('HTaqueous','fuel_gas','biooil'))

Acidex = su.AcidExtraction('A200',ins=(HTL-0,acidforP),outs=('residual','extracted'))

M1 = su.HTLmixer('A300',ins=(HTL-1,Acidex-1),outs=('mixture'))

StruPre = su.StruvitePrecipitation('A310',ins=(M1-0,supply_mgcl2),outs=('struvite','CHGfeed'))

CHG = su.CHG('A330',ins=StruPre-1,outs=('fuelgas','effluent'))

MemDis = su.MembraneDistillation('A340',ins=(CHG-1,acidforN),outs=('AmmoniaSulfate','ww'))

sys=qs.System('sys',path=(SluL,SluT,SluC,HTL,HT,Acidex,M1,StruPre,CHG,MemDis))

sys.simulate()

sys.diagram()

# return sys

#%%
from qsdsan import Model

model = Model(sys)

from chaospy import distributions as shape

param=model.parameter

dist = shape.Triangle(0.9895,0.99,0.9905)
@param(name='sludge_moisture',
        element=SluL,
        kind='coupled',
        units='-',
        baseline=0.99,
        distribution=dist)
def set_sludge_moisture(i):
    SluL.sludge_moisture=i

dist = shape.Triangle(0.29,0.32575,0.376)
@param(name='sludge_dw_protein',
        element=SluL,
        kind='coupled',
        units='-',
        baseline=0.341,
        distribution=dist)
def set_sludge_dw_protein(i):
    SluL.sludge_dw_protein=i

dist = shape.Triangle(0.167,0.22925,0.308)
@param(name='sludge_dw_carbo',
        element=SluL,
        kind='coupled',
        units='-',
        baseline=0.167,
        distribution=dist)
def set_sludge_dw_carbo(i):
    SluL.sludge_dw_carbo=i

dist = shape.Triangle(0.116,0.1725,0.226)
@param(name='sludge_dw_lipid',
        element=SluL,
        kind='coupled',
        units='-',
        baseline=0.226,
        distribution=dist)
def set_sludge_dw_lipid(i):
    SluL.sludge_dw_lipid=i

dist = shape.Triangle(0.343,0.4206,0.478)
@param(name='sludge_C_ratio',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.411,
        distribution=dist)
def set_sludge_C_ratio(i):
    HTL.sludge_C_ratio=i

dist = shape.Triangle(0.007,0.0168,0.02)
@param(name='sludge_P_ratio',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.019,
        distribution=dist)
def set_sludge_P_ratio(i):
    HTL.sludge_P_ratio=i

dist = shape.Triangle(0.047,0.0586,0.065)
@param(name='sludge_H_ratio',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.058,
        distribution=dist)
def set_sludge_H_ratio(i):
    HTL.sludge_H_ratio=i

dist = shape.Triangle(0.005,0.0096,0.016)
@param(name='sludge_S_ratio',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.01,
        distribution=dist)
def set_sludge_S_ratio(i):
    HTL.sludge_S_ratio=i

dist = shape.Triangle(0.036,0.0484,0.061)
@param(name='sludge_N_ratio',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.05,
        distribution=dist)
def set_sludge_N_ratio(i):
    HTL.sludge_N_ratio=i

dist = shape.Triangle(0.261,0.2808,0.336)
@param(name='sludge_O_ratio',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.261,
        distribution=dist)
def set_sludge_O_ratio(i):
    HTL.sludge_O_ratio=i

dist = shape.Triangle(9.57,15.5,23.8)
@param(name='biochar_C_N_ratio',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=15.5,
        distribution=dist)
def set_biochar_C_N_ratio(i):
    HTL.biochar_C_N_ratio=i

dist = shape.Triangle(1.49,2.16,2.90)
@param(name='biochar_C_P_ratio',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=2.16,
        distribution=dist)
def set_biochar_C_P_ratio(i):
    HTL.biochar_C_P_ratio=i

dist = shape.Triangle(0.035,0.0648,0.102)
@param(name='biocrude_moisture_content',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.044,
        distribution=dist)
def set_biocrude_moisture_content(i):
    HTL.biocrude_moisture_content=i

dist = shape.Triangle(0.75,0.78,0.82)
@param(name='biooil_ratio',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.77,
        distribution=dist)
def set_biooil_ratio(i):
    HT.biooil_ratio=i

dist = shape.Triangle(0.04,0.073,0.1)
@param(name='gas_ratio',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.07,
        distribution=dist)
def set_gas_ratio(i):
    HT.gas_ratio=i

dist = shape.Normal(0.128,0.00064)
@param(name='co_ratio',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.128,
        distribution=dist)
def set_co_ratio_HT(i):
    HT.co_ratio=i

dist = shape.Normal(0.007,0.000035)
@param(name='co2_ratio',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.007,
        distribution=dist)
def set_co2_ratio_HT(i):
    HT.co2_ratio=i

dist = shape.Normal(0.188,0.00094)
@param(name='c2h6_ratio',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.188,
        distribution=dist)
def set_c2h6_ratio_HT(i):
    HT.c2h6_ratio=i

dist = shape.Normal(0.107,0.000535)
@param(name='c3h8_ratio',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.107,
        distribution=dist)
def set_c3h8_ratio_HT(i):
    HT.c3h8_ratio=i

dist = shape.Normal(0.09,0.00045)
@param(name='c5h12_ratio',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.09,
        distribution=dist)
def set_c5h12_ratio_HT(i):
    HT.c5h12_ratio=i

dist = shape.Triangle(0.846,0.854,0.86)
@param(name='biooil_C_ratio',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.855,
        distribution=dist)
def set_biooil_C_ratio(i):
    HT.biooil_C_ratio=i

dist = shape.Triangle(0.0004,0.004111,0.016)
@param(name='biooil_N_ratio',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.01,
        distribution=dist)
def set_biooil_N_ratio(i):
    HT.biooil_N_ratio=i

dist = shape.Triangle(0.92625,0.95,0.97375)
@param(name='P_acid_recovery_ratio',
        element=Acidex,
        kind='coupled',
        units='-',
        baseline=0.95,
        distribution=dist)
def set_P_acid_recovery_ratio(i):
    Acidex.P_acid_recovery_ratio=i

dist = shape.Triangle(0.92625,0.95,0.97375)
@param(name='P_pre_recovery_ratio',
        element=StruPre,
        kind='coupled',
        units='-',
        baseline=0.95,
        distribution=dist)
def set_P_pre_recovery_ratio(i):
    StruPre.P_pre_recovery_ratio=i

dist = shape.Triangle(0.097,0.11,0.127)
@param(name='P_in_struvite',
        element=StruPre,
        kind='coupled',
        units='-',
        baseline=0.127,
        distribution=dist)
def set_P_in_struvite(i):
    StruPre.P_in_struvite=i

dist = shape.Normal(0.244,0.00122)
@param(name='ch4_ratio',
        element=CHG,
        kind='coupled',
        units='-',
        baseline=0.244,
        distribution=dist)
def set_ch4_ratio_CHG(i):
    CHG.ch4_ratio=i

dist = shape.Normal(0.029,0.000145)
@param(name='co_ratio',
        element=CHG,
        kind='coupled',
        units='-',
        baseline=0.029,
        distribution=dist)
def set_co_ratio_CHG(i):
    CHG.co_ratio=i

dist = shape.Normal(0.15,0.00075)
@param(name='co2_ratio',
        element=CHG,
        kind='coupled',
        units='-',
        baseline=0.15,
        distribution=dist)
def set_co2_ratio_CHG(i):
    CHG.co2_ratio=i

dist = shape.Normal(0.043,0.000215)
@param(name='c2h6_ratio',
        element=CHG,
        kind='coupled',
        units='-',
        baseline=0.043,
        distribution=dist)
def set_c2h6_ratio_CHG(i):
    CHG.c2h6_ratio=i

dist = shape.Normal(0.764,0.049)
@param(name='toc_tc_ratio',
        element=CHG,
        kind='coupled',
        units='-',
        baseline=0.764,
        distribution=dist)
def set_toc_tc_ratio(i):
    CHG.toc_tc_ratio=i

dist = shape.Normal(0.262,0.06)
@param(name='toc_to_gas_c_ratio',
        element=CHG,
        kind='coupled',
        units='-',
        baseline=0.262,
        distribution=dist)
def set_toc_to_gas_c_ratio(i):
    CHG.toc_to_gas_c_ratio=i

dist = shape.Triangle(0.8,0.825,0.85)
@param(name='N_recovery_rate',
        element=MemDis,
        kind='coupled',
        units='-',
        baseline=0.825,
        distribution=dist)
def set_N_recovery_rate(i):
    MemDis.N_recovery_rate=i   

# metric = model.metric
# @metric(name='P_recovery_rate',units='-',element='Production')
# def get_P_recovery_rate():
#     return StruPre.struvite_P/((HTL.ins[0].F_mass-HTL.ins[0].imass['H2O'])*HTL.sludge_P_ratio)

metric = model.metric
@metric(name='Struvite',units='kg/hr',element='Production')
def get_struvite_production():
    return StruPre.outs[0].imass['Struvite']

@metric(name='(NH4)2SO4',units='kg/hr',element='Production')
def get_nh42so4_production():
    return MemDis.outs[0].imass['NH42SO4']

@metric(name='Biooil',units='kg/hr',element='Production')
def get_bioil_production():
    return HT.outs[2].imass['Biooil']

@metric(name='H2',units='kg/hr',element='Production')
def get_h2_production():
    return CHG.outs[0].imass['H2']

@metric(name='CO',units='kg/hr',element='Production')
def get_co_production():
    return HT.outs[1].imass['CO']+CHG.outs[0].imass['CO']

@metric(name='CO2',units='kg/hr',element='Production')
def get_co2_production():
    return HTL.outs[3].imass['CO2']+HT.outs[1].imass['CO']+CHG.outs[0].imass['CO']

@metric(name='CH4',units='kg/hr',element='Production')
def get_ch4_production():
    return HT.outs[1].imass['CH4']+CHG.outs[0].imass['CH4']

#%%
import numpy as np
np.random.seed(3221)
samples = model.sample(N=1000, rule='L')
model.load_samples(samples)
model.evaluate()
model.table
 #%%
fig, ax = qs.stats.plot_uncertainties(model)
fig
#%%
fig, ax = qs.stats.plot_uncertainties(model, x_axis=model.metrics[0], y_axis=model.metrics[1],
                                      kind='kde-kde', center_kws={'fill': True})
fig
#%%
r_df, p_df = qs.stats.get_correlations(model, kind='Spearman')
fig, ax = qs.stats.plot_correlations(r_df)
fig