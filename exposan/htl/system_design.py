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

(1) Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
    Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
    Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
    Process Design and Economics for the Conversion of Algal Biomass to
    Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
    PNNL--23227, 1126336; 2014; 
    https://doi.org/10.2172/1126336.
    
(2) Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels via
    Liquefaction - Hydrothermal Liquefaction Reactor Design: April 5, 2013;
    NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462, 1111191.
    https://doi.org/10.2172/1111191.
'''

import qsdsan as qs
import exposan.htl._sanunits_design as su
from qsdsan import sanunits as qsu
from biosteam.units import Flash, IsenthalpicValve, BinaryDistillation
from exposan.htl._process_settings_design import load_process_settings
from exposan.htl._components_design import create_components

# __all__ = ('create_system',)

# def create_system():
    
load_process_settings()
cmps = create_components()

fake_sludge = qs.Stream('sludge', H2O=100000, units='kg/hr', T=25+273.15)
# set H2O equal to the total sludge input flow
# assume 99% moisture, 50 metric tons of dw sludge per h

# =============================================================================
# pretreatment (Area 000)
# =============================================================================

SluL = su.SludgeLab('S000', ins=fake_sludge, outs='real_sludge',
                    sludge_moisture=0.99, sludge_dw_protein=0.341,
                    sludge_dw_lipid=0.226, sludge_dw_carbo=0.167)

SluT = qsu.SludgeThickening('A000', ins=SluL-0, 
                            outs=('supernatant_1','compressed_sludge_1'),
                            init_with='Stream', 
                            solids=('Sludge_lipid','Sludge_protein',
                                    'Sludge_carbo','Sludge_ash'))

SluC = qsu.SludgeCentrifuge('A010', ins=SluT-1,
                            outs=('supernatant_2','compressed_sludge_2'),
                            init_with='Stream',
                            solids=('Sludge_lipid','Sludge_protein',
                                    'Sludge_carbo','Sludge_ash'))

# =============================================================================
# HTL (Area 100)
# =============================================================================

P1 = qsu.Pump('A100', ins=SluC-1, outs='press_sludge', P=3049.7*6894.76,
              init_with='Stream')
# Jones 2014: 3049.7 psia

H1 = qsu.HXutility('A110', ins=P1-0, outs='heated_sludge', T=351+273.15,
                   U=0.874, init_with='Stream')
# H1: NREL 2013: 154 (40-446) Btu/hr/ft2/F ~ U = 0.874 (0.227-2.531) kW/m2/k
# U is just needed for H1? Right? I think high viscosity of sludge is just here
# but not in other pumps
# unit conversion: http://www.unitconversion.org/heat-transfer-coefficient/
# watts-per-square-meter-per-k-to-btus-th--per-hour-per-square-foot-per-f-
# conversion.html

HTL = su.HTL('A120', ins=(H1-0), outs=('biochar','HTL_aqueous',
             'biocrude','offgas_HTL'))
HTL_hx = HTL.heat_exchanger

# =============================================================================
# CHG (Area 200)
# =============================================================================

H2SO4_Tank = qsu.StorageTank('T200', ins='H2SO4', outs=('H2SO4_out'),
                             init_with='Stream', tau=3*24)

SP1 = su.HTLsplitter('S200',ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                     init_with='Stream')
# must put after AcidEx and MemDis in path during simulation to ensure input
# not empty

AcidEx = su.AcidExtraction('A200', ins=(HTL-0,SP1-0),
                           outs=('residual','extracted'))

M1 = su.HTLmixer('A210', ins=(HTL-1,AcidEx-1), outs=('mixture'))

StruPre = su.StruvitePrecipitation('A220', ins=(M1-0,'MgCl2','NH4Cl','NaOH'),
                                   outs=('struvite','CHG_feed'))

P2 = qsu.Pump('A230', ins=StruPre-1, outs='press_aqueous',
              P=3089.7*6894.76, init_with='Stream') # Jones 2014: 3089.7 psia

H2 = qsu.HXutility('A240', ins=P2-0, outs='heated_aqueous',
                   T=350+273.15, init_with='Stream')

CHG = su.CHG('A250', ins=H2-0, outs='CHG_out')

H3 = qsu.HXutility('A260', ins=CHG-0, outs='cooled_CHG', T=60+273.15,
                   init_with='Stream')


F1 = Flash('A270', ins=H3-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
           T=60+273.15, P=50*6894.76)

MemDis = su.MembraneDistillation('A280', ins=(F1-1, SP1-1),
                                 outs=('ammonium_sulfate','MemDis_ww'))
#ammonium sulfate fate? crystallation or transport out?

# =============================================================================
# HT (Area 300)
# =============================================================================

P3 = qsu.Pump('A300', ins=HTL-2, outs='press_biocrude', P=1530.0*6894.76,
              init_with='Stream')
# Jones 2014: 1530.0 psia

# Tin = 174 C (345 F) based on Jones PNNL report. However, the reaction
# releases a lot of heat and increase the temperature of effluent to 402 C
# (755.5 F).

RSP1 = qsu.ReversedSplitter('S600', ins='H2', outs=('HT_H2','HC_H2'),
                            init_with='Stream')
# reversed splitter, write before HT and HC, simulate after HT and HC

HT = su.HT('A310', ins=(P3-0, RSP1-0), outs='HTout')
HT_compressor = HT.compressor
HT_hx = HT.heat_exchanger

V1 = IsenthalpicValve('A320', ins=HT-0, outs='depressed_HT', P=717.4*6894.76)

H4 = qsu.HXutility('A330', ins=V1-0, outs='cooled_HT', T=60+273.15,
                   init_with='Stream')

F2 = Flash('A340', ins=H4-0, outs=('HT_fuel_gas','HT_aqueous'), T=43+273.15,
           P=717.4*6894.76) # outflow P

V2 = IsenthalpicValve('A350', ins=F2-1, outs='depressed_flash_effluent', P=55*6894.76)

SP2 = qsu.Splitter('S300', ins=V2-0, outs=('HT_ww','HT_oil'),
                   split={'H2O':1}, init_with='Stream')
# separate water and oil based on gravity

H5 = qsu.HXutility('A360', ins=SP2-1, outs='heated_oil', T=104+273.15)
# temperature: Jones stream #334 (we remove the first distillation column)

C1 = BinaryDistillation('A370', ins=H5-0,
                        outs=('HT_light','HT_heavy'),
                        LHK=('C4H10','TWOMBUTAN'), P=50*6894.76, # outflow P
                        y_top=188/253, x_bot=53/162, k=2, is_divided=True)

C2 = BinaryDistillation('A380', ins=C1-1,
                        outs=('HT_Gasoline','HT_other_oil'),
                        LHK=('C10H22','C4BENZ'), P=25*6894.76, # outflow P
                        y_top=116/122, x_bot=114/732, k=2, is_divided=True)

C3 = BinaryDistillation('A390', ins=C2-1,
                        outs=('HT_Diesel','HT_heavy_oil'),
                        LHK=('C19H40','C21H44'),P=18.7*6894.76, # outflow P
                        y_top=2421/2448, x_bot=158/2448, k=2, is_divided=True)

# =============================================================================
# HC (Area 400)
# =============================================================================

P4 = qsu.Pump('A400', ins=C3-1, outs='press_heavy_oil', P=1034.7*6894.76,
              init_with='Stream')
# Jones 2014: 1034.7 psia

# Tin = 394 C (741.2 F) based on Jones PNNL report. However, the reaction
# releases a lot of heat and increase the temperature of effluent to 451 C
# (844.6 F).

HC = su.HC('A410', ins=(P4-0, RSP1-1), outs='HC_out')
HC_compressor = HC.compressor
HC_hx = HC.heat_exchanger

H6 = qsu.HXutility('A420', ins=HC-0, outs='cooled_HC', T=60+273.15,
                   init_with='Stream')

F3 = Flash('A430', ins=H6-0, outs=('HC_fuel_gas','HC_aqueous'), T=60.2+273,
           P=30*6894.76) # outflow P

C4 = BinaryDistillation('A440', ins=F3-1, outs=('HC_Gasoline','HC_Diesel'),
                        LHK=('C9H20','C10H22'), P=20*6894.76, # outflow P
                        y_top=360/546, x_bot=7/708, k=2, is_divided=True)

# =============================================================================
# CHP and storage (Area 500)
# =============================================================================

GasolineMixer = qsu.Mixer('S500', ins=(C2-0, C4-0), outs='mixed_gasoline',
                          init_with='Stream')

DieselMixer = qsu.Mixer('S510', ins=(C3-0, C4-1), outs='mixed_diesel',
                        init_with='Stream')

H7 = qsu.HXutility('A500', ins=GasolineMixer-0, outs='cooled_gasoline',
                   T=60+273.15, init_with='Stream', rigorous=True)

H8 = qsu.HXutility('A510', ins=DieselMixer-0, outs='cooled_diesel',
                   T=60+273.15, init_with='Stream', rigorous=True)

PC1 = su.PhaseChanger('S520', ins=H7-0, outs='cooled_gasoline_liquid')

PC2 = su.PhaseChanger('S530', ins=H8-0, outs='cooled_diesel_liquid')

GasolineTank = qsu.StorageTank('T500', ins=PC1-0, outs=('gasoline'),
                               tau=3*24, init_with='Stream')
# store for 3 days based on Jones 2014

DieselTank = qsu.StorageTank('T510', ins=PC2-0, outs=('diesel'),
                             tau=3*24, init_with='Stream')
# store for 3 days based on Jones 2014

GasMixer = qsu.Mixer('S540', ins=(HTL-3, F1-0, F2-0, C1-0, F3-0),
                     outs=('fuel_gas'), init_with='Stream')
# F3-0 and F5-0 are empty, may change remove/modify F3 and/or F5 later

CHP = qsu.CHP('A520', ins=(GasMixer-0,'natural_gas','air'),
              outs=('emission','solid_ash'), init_with='Stream')

WWmixer = su.WWmixer('S550', ins=(SluT-0, SluC-0, MemDis-1, SP2-0),
                    outs='wastewater', init_with='Stream')
# effluent of WWmixer goes back to WWTP

# =============================================================================
# facilities
# =============================================================================

# HXN = qsu.HeatExchangerNetwork('HXN')

# for unit in (SluL, SluT, SluC, P1, H1, HTL, H2SO4_Tank, AcidEx,
#              M1, StruPre, P2, H2, CHG, F1, MemDis, SP1, P3, H3,
#              HT, F2, F3, SP3, HX_biocrude, C1, C2, P4, H4, HC,
#              F4, F5, C3, GasolineMixer, DieselMixer, H5, H6,
#              GasolineTank, DieselTank, GasMixer, CHP, SP2,
#              HX_H2_HT, IC_H2_HT, HX_H2_HC, IC_H2_HC):
#     unit.register_alias(f'{unit=}'.split('=')[0].split('.')[-1])
# so that qs.main_flowsheet.H1 works as well

sys = qs.System('sys', path=(SluL, SluT, SluC, P1, H1, HTL, H2SO4_Tank, AcidEx,
                             M1, StruPre, P2, H2, CHG, H3, F1, MemDis, SP1,
                             P3, HT, V1, H4, F2, V2, SP2, H5, C1, C2, C3, P4,
                             HC, H6, F3, C4, GasolineMixer, DieselMixer,
                             H7, H8, PC1, PC2, GasolineTank, DieselTank,
                             GasMixer, CHP, WWmixer, RSP1)) 
                # , facilities=(HXN,))

sys.operating_hours = 7884 # NRES 2013

sys.simulate()

sys.diagram()

# return sys

#%%
from qsdsan import Model

model = Model(sys)

from chaospy import distributions as shape

param = model.parameter

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

dist = shape.Triangle(0.007,0.0168,0.02)
@param(name='sludge_P_ratio',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.019,
        distribution=dist)
def set_sludge_P_ratio(i):
    HTL.sludge_P_ratio=i

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

dist = shape.Triangle(0.5,0.7,0.9)
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

metric = model.metric
@metric(name='Struvite',units='kg/hr',element='Production')
def get_struvite_production():
    return StruPre.outs[0].imass['Struvite']

@metric(name='(NH4)2SO4',units='kg/hr',element='Production')
def get_nh42so4_production():
    return MemDis.outs[0].imass['NH42SO4']

@metric(name='Gasoline',units='kg/hr',element='Production')
def get_gasoline_production():
    return GasolineTank.outs[0].F_mass

@metric(name='Diesel',units='kg/hr',element='Production')
def get_diesel_production():
    return DieselTank.outs[0].F_mass

#%%
import numpy as np
np.random.seed(3221)
samples = model.sample(N=10, rule='L')
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