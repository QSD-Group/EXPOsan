#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Jianan Feng <jiananf2@illinois.edu>
    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References:
(1) Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels via
    Liquefaction - Hydrothermal Liquefaction Reactor Design: April 5, 2013;
    NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462, 1111191.
    https://doi.org/10.2172/1111191.

'''

import qsdsan as qs
import biosteam as bst
import exposan.htl._sanunits as su
from qsdsan import sanunits as suu
from exposan.htl._process_settings import load_process_settings
from exposan.htl._components import create_components


# __all__ = ('create_system',)

# def create_system():
    
load_process_settings()
cmps = create_components()

fake_sludge = qs.Stream('fake_sludge',H2O=1814370,units='kg/hr',T=298.15)
#set H2O equal to the total sludge input flow
#assume 99% moisture, 20 us tons of dw sludge per h


SluL = su.SludgeLab('S000',ins=fake_sludge,outs='real_sludge')

SluT = suu.SludgeThickening('A000',ins=SluL-0,outs=('supernatant_1','compressed_sludge_1'),
                            init_with='Stream',solids=('Sludge_lipid','Sludge_protein','Sludge_carbo',
                                                       'Sludge_ash'))

SluC = suu.SludgeCentrifuge('A010',ins=SluT-1,outs=('supernatant_2','compressed_sludge_2'),
                            init_with='Stream',solids=('Sludge_lipid','Sludge_protein','Sludge_carbo',
                                                       'Sludge_ash'))

P1 = suu.Pump('A100',ins=SluC-1,outs='press_sludge',P=3049.7*6894.76) #Jones 2014: 3049.7 psia

H1 = suu.HXutility('A110',ins=P1-0,outs='heated_sludge',T=350+273.15,U=0.874,init_with='Stream')
# H1: from NREL 2013: 154 (40-446) Btu/hr/ft2/F ~ U = 0.874 (0.227-2.531) kW/m2/k
# U is just needed for H1? Right? I think high viscosity of sludge is just here but not in other pumps
# unit conversion: http://www.unitconversion.org/heat-transfer-coefficient/watts-
# per-square-meter-per-k-to-btus-th--per-hour-per-square-foot-per-f-conversion.html

#!!! currently the power/heat utilities of HTL_hx aren't added to HTL,
# it'll be added in after I pulled in some updates from biosteam,
# so we don't need to worry about it (I'm leaving this note to reminder myself)

HTL = su.HTL('A120',ins=H1-0,outs=('biochar','HTLaqueous','biocrude','offgas_HTL'))
HTL_hx = HTL.heat_exchanger #include this into path?

AcidEx = su.AcidExtraction('A200',ins=(HTL-0,'H2SO4_P'),outs=('residual','extracted'))

M1 = su.HTLmixer('A210',ins=(HTL-1,AcidEx-1),outs=('mixture'))

StruPre = su.StruvitePrecipitation('A220',ins=(M1-0,'MgCl2'),outs=('struvite','CHG_feed'))

P2 = suu.Pump('A230',ins=StruPre-1,outs='press_aqueous',P=3089.7*6894.76) #Jones 2014: 3089.7 psia

H2 = suu.HXutility('A240',ins=P2-0,outs='heated_aqueous',T=350+273.15,init_with='Stream')

CHG = su.CHG('A250',ins=H2-0,outs=('CHG_fuel_gas','effluent'))
CHG_hx = CHG.heat_exchanger

MemDis = su.MembraneDistillation('A260',ins=(CHG-1,'H2SO4_N'),outs=('AmmoniaSulfate','ww'))

P3 = suu.Pump('A300',ins=HTL-2,outs='press_biocrude',P=1530.0*6894.76) #Jones 2014: 1530.0 psia

H3 = suu.HXutility('A310',ins=P3-0,outs='heated_biocrude',T=405+273.15,init_with='Stream')

HT = su.HT('A320',ins=(H3-0,'H2_HT'),outs=('HTaqueous','HT_fuel_gas','gasoline_HT','diesel_HT','heavy_oil'))
HT_hx = HT.heat_exchanger

P4 = suu.Pump('A330',ins=HT-4,outs='press_heavy_oil',P=1034.7*6894.76) #Jones 2014: 1034.7 psia

H4 = suu.HXutility('A340',ins=P4-0,outs='heated_heavy_oil',T=395+273.15,init_with='Stream')
#All temperatures and pressures are from Jones et al., 2014

HC = su.HC('A350',ins=(H4-0,'H2_HC'),outs=('gasoline_HC', 'diesel_HC', 'offgas_HC'))
HC_hx = HC.heat_exchanger

GasolineMixer = su.FuelMixer('S100',ins=(HT-2,HC-0),outs='mixed_gasoline')

DieselMixer = su.FuelMixer('S110',ins=(HT-3,HC-1),outs='mixed_diesel')

H5 = suu.HXutility('A360',ins=GasolineMixer-0,outs='cooled_gasoline',T=60+273.15,init_with='Stream')

H6 = suu.HXutility('A370',ins=DieselMixer-0,outs='cooled_diesel',T=60+273.15,init_with='Stream')

GasolineTank = suu.StorageTank('A380',ins=H5-0,outs=('gasoline_out'),tau=3*24,init_with='Stream')
#store for 3 days based on Jones 2014

DieselTank = suu.StorageTank('A390',ins=H6-0,outs=('diesel_out'),tau=3*24,init_with='Stream')
#store for 3 days based on Jones 2014

GasMixer = su.GasMixer('S200',ins=(HTL-3,CHG-0,HT-1,HC-2),outs=('fuel_gas'))

CHP = suu.CHP('A400',ins=(GasMixer-0,'natural_gas','air'),outs=('emission','solid_ash'))

# HXN = suu.HeatExchangerNetwork('HXN')

# for unit in (SluL,SluT,SluC,P1,H1,HTL,AcidEx,M1,StruPre,P2,H2,CHG,MemDis,
#              P3,H3,HT,P4,H4,HC,GasolineMixer,DieselMixer,H5,H6,GasolineTank,
#              DieselTank,GasMixer,CHP):
#     unit.register_alias(f'{unit=}'.split('=')[0].split('.')[-1]) # so that qs.main_flowsheet.H1 works as well

sys=qs.System('sys',path=(SluL,SluT,SluC,
                          P1,H1,HTL,
                          AcidEx,M1,StruPre,
                          P2,H2,CHG,MemDis,
                          P3,H3,HT,
                          P4,H4,HC,
                          GasolineMixer,DieselMixer,
                          H5,H6,
                          GasolineTank,DieselTank,
                          GasMixer,CHP))#,facilities=(HXN,))


sys.operating_hours=7884 # NRES 2013

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
@param(name='c4h10_ratio',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.09,
        distribution=dist)
def set_c4h10_ratio_HT(i):
    HT.c4h10_ratio=i

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
        element=AcidEx,
        kind='coupled',
        units='-',
        baseline=0.95,
        distribution=dist)
def set_P_acid_recovery_ratio(i):
    AcidEx.P_acid_recovery_ratio=i

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

@metric(name='CO2',units='kg/hr',element='Production')
def get_co2_production():
    return CHP.outs[0].imass['CO2']


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
fig, ax = qs.stats.plot_uncertainties(model, x_axis=model.metrics[0], y_axis=model.metrics[1],
                                      kind='kde-kde', center_kws={'fill': True})
fig
#%%
r_df, p_df = qs.stats.get_correlations(model, kind='Spearman')
fig, ax = qs.stats.plot_correlations(r_df)
fig