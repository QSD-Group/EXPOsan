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
    PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    
(2) Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels via
    Liquefaction - Hydrothermal Liquefaction Reactor Design: April 5, 2013;
    NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462, 1111191.
    https://doi.org/10.2172/1111191.
'''

import qsdsan as qs
import exposan.htl._sanunits as su
from qsdsan import sanunits as qsu
from biosteam.units import Flash, IsenthalpicValve, BinaryDistillation
from exposan.htl._process_settings import load_process_settings
from exposan.htl._components import create_components
from exposan.htl._TEA import *

# __all__ = ('create_system',)

# def create_system():

load_process_settings()
cmps = create_components()

fake_sludge = qs.Stream('sludge', H2O=1000000*365/7920, units='kg/hr', T=25+273.15)
# set H2O equal to the total sludge input flow
# assume 99% moisture, 10 metric tons of dw sludge per d: (1000000*365/7920) kg/hr

# =============================================================================
# pretreatment (Area 000)
# =============================================================================

SluL = su.SludgeLab('S000', ins=fake_sludge, outs='real_sludge',
                    sludge_moisture=0.99, sludge_dw_ash=0.266, 
                    sludge_afdw_protein=0.465, sludge_afdw_lipid=0.308, yearly_operation_hour=7920)

SluT = qsu.SludgeThickening('A000', ins=SluL-0,
                            outs=('supernatant_1','compressed_sludge_1'),
                            init_with='Stream', 
                            solids=('Sludge_lipid','Sludge_protein',
                                    'Sludge_carbo','Sludge_ash'),
                            sludge_moisture=0.96)

SluC = qsu.SludgeCentrifuge('A010', ins=SluT-1,
                            outs=('supernatant_2','compressed_sludge_2'),
                            init_with='Stream',
                            solids=('Sludge_lipid','Sludge_protein',
                                    'Sludge_carbo','Sludge_ash'),
                            sludge_moisture=0.8)

# =============================================================================
# HTL (Area 100)
# =============================================================================

P1 = qsu.Pump('A100', ins=SluC-1, outs='press_sludge', P=3049.7*6894.76,
              init_with='Stream')
# Jones 2014: 3049.7 psia

H1 = qsu.HXutility('A110', ins=P1-0, outs='heated_sludge', T=351+273.15,
                   U=0.284, init_with='Stream')
# H1: SS PNNL 2020: 50 (17-76) Btu/hr/ft2/F ~ U = 0.284 (0.096-0.4313) kW/m2/k
# U is just needed for H1? Right? I think high viscosity of sludge is just here
# but not in other pumps
# unit conversion: http://www.unitconversion.org/heat-transfer-coefficient/
# watts-per-square-meter-per-k-to-btus-th--per-hour-per-square-foot-per-f-
# conversion.html

HTL = su.HTL('A120', ins=H1-0, outs=('biochar','HTL_aqueous',
             'biocrude','offgas_HTL'))
HTL_hx = HTL.heat_exchanger
HTL_drum = HTL.kodrum

# =============================================================================
# CHG (Area 200)
# =============================================================================

H2SO4_Tank = qsu.StorageTank('T200', ins='H2SO4', outs=('H2SO4_out'),
                             init_with='Stream', tau=3*24)
H2SO4_Tank.ins[0].price = 0.0055 # based on 93% H2SO4 and fresh water (dilute to 5%) price found in Davis 2016$/kg

SP1 = su.HTLsplitter('S200',ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                     init_with='Stream')
# must put after AcidEx and MemDis in path during simulation to ensure input
# not empty

AcidEx = su.AcidExtraction('A200', ins=(HTL-0, SP1-0),
                           outs=('residual','extracted'))

M1 = su.HTLmixer('A210', ins=(HTL-1, AcidEx-1), outs=('mixture'))

StruPre = su.StruvitePrecipitation('A220', ins=(M1-0,'MgCl2','NH4Cl','MgO'),
                                   outs=('struvite','CHG_feed'))

CHG = su.CHG('A250', ins=(StruPre-1, '7.8% Ru/C'), outs=('CHG_out', '7.8% Ru/C_out'))
CHG_pump = CHG.pump
CHG_heating = CHG.heat_ex_heating
CHG_cooling = CHG.heat_ex_cooling

V1 = IsenthalpicValve('A270', ins=CHG-0, outs='depressed_cooled_CHG', P=50*6894.76)

F1 = Flash('A280', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
            T=60+273.15, P=50*6894.76)

MemDis = su.MembraneDistillation('A290', ins=(F1-1, SP1-1, 'NaOH', 'Membrane_in'),
                                  outs=('ammonium_sulfate','MemDis_ww', 'Membrane_out'))

# =============================================================================
# HT (Area 300)
# =============================================================================

P3 = qsu.Pump('A300', ins=HTL-2, outs='press_biocrude', P=1530.0*6894.76,
              init_with='Stream')
# Jones 2014: 1530.0 psia

# Tin = 174 C (345 F) based on Jones PNNL report. However, the reaction
# releases a lot of heat and increase the temperature of effluent to 402 C
# (755.5 F).

RSP1 = qsu.ReversedSplitter('S300', ins='H2', outs=('HT_H2','HC_H2'),
                            init_with='Stream')
# reversed splitter, write before HT and HC, simulate after HT and HC

HT = su.HT('A310', ins=(P3-0, RSP1-0, 'CoMo_alumina_HT'), outs=('HTout', 'CoMo_alumina_HT_out'))
HT_compressor = HT.compressor
HT_hx = HT.heat_exchanger

V2 = IsenthalpicValve('A320', ins=HT-0, outs='depressed_HT', P=717.4*6894.76)

H4 = qsu.HXutility('A330', ins=V2-0, outs='cooled_HT', T=60+273.15,
                    init_with='Stream')

F2 = Flash('A340', ins=H4-0, outs=('HT_fuel_gas','HT_aqueous'), T=43+273.15,
            P=717.4*6894.76) # outflow P

V3 = IsenthalpicValve('A350', ins=F2-1, outs='depressed_flash_effluent', P=55*6894.76)

SP2 = qsu.Splitter('S310', ins=V3-0, outs=('HT_ww','HT_oil'),
                    split={'H2O':1}, init_with='Stream')
# separate water and oil based on gravity

H5 = qsu.HXutility('A360', ins=SP2-1, outs='heated_oil', T=104+273.15)
# temperature: Jones stream #334 (we remove the first distillation column)

D1 = BinaryDistillation('A370', ins=H5-0,
                        outs=('HT_light','HT_heavy'),
                        LHK=('C4H10','TWOMBUTAN'), P=50*6894.76, # outflow P
                        y_top=188/253, x_bot=53/162, k=2, is_divided=True)

D2 = BinaryDistillation('A380', ins=D1-1,
                        outs=('HT_Gasoline','HT_other_oil'),
                        LHK=('C10H22','C4BENZ'), P=25*6894.76, # outflow P
                        y_top=116/122, x_bot=114/732, k=2, is_divided=True)

D3 = BinaryDistillation('A390', ins=D2-1,
                        outs=('HT_Diesel','HT_heavy_oil'),
                        LHK=('C19H40','C21H44'),P=18.7*6894.76, # outflow P
                        y_top=2421/2448, x_bot=158/2448, k=2, is_divided=True)

# =============================================================================
# HC (Area 400)
# =============================================================================

P4 = qsu.Pump('A400', ins=D3-1, outs='press_heavy_oil', P=1034.7*6894.76,
              init_with='Stream')
# Jones 2014: 1034.7 psia

# Tin = 394 C (741.2 F) based on Jones PNNL report. However, the reaction
# releases a lot of heat and increase the temperature of effluent to 451 C
# (844.6 F).

HC = su.HC('A410', ins=(P4-0, RSP1-1, 'CoMo_alumina_HC'), outs=('HC_out', 'CoMo_alumina_HC_out'))
HC_compressor = HC.compressor
HC_hx = HC.heat_exchanger

H6 = qsu.HXutility('A420', ins=HC-0, outs='cooled_HC', T=60+273.15,
                    init_with='Stream')

V4 = IsenthalpicValve('A430', ins=H6-0, outs='cooled_depressed_HC', P=30*6894.76)


F3 = Flash('A440', ins=V4-0, outs=('HC_fuel_gas','HC_aqueous'), T=60.2+273,
            P=30*6894.76) # outflow P

D4 = BinaryDistillation('A450', ins=F3-1, outs=('HC_Gasoline','HC_Diesel'),
                        LHK=('C9H20','C10H22'), P=20*6894.76, # outflow P
                        y_top=360/546, x_bot=7/708, k=2, is_divided=True)

# =============================================================================
# CHP and storage (Area 500)
# =============================================================================

GasolineMixer = qsu.Mixer('S500', ins=(D2-0, D4-0), outs='mixed_gasoline',
                          init_with='Stream')

DieselMixer = qsu.Mixer('S510', ins=(D3-0, D4-1), outs='mixed_diesel',
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

FuelMixer = su.FuelMixer('S540', ins=(GasolineTank-0, DieselTank-0),\
                         outs='fuel', target='diesel')
# integrate gasoline and diesel based on their LHV for MFSP calculation

GasMixer = qsu.Mixer('S550', ins=(HTL-3, F1-0, F2-0, D1-0, F3-0),
                      outs=('fuel_gas'), init_with='Stream')

# The system produces more energy than needed (heating+power)
CHP = qsu.CHP('A520', ins=(GasMixer-0,'natural_gas','air'),
              outs=('emission','solid_ash'), init_with='Stream', supplement_power_utility=True)

CHP.ins[1].price = 5.1/1000/0.02391792567 # from $/scf to $/kg, values from Jones

WWmixer = su.WWmixer('S560', ins=(SluT-0, SluC-0, MemDis-1, SP2-0),
                    outs='wastewater', init_with='Stream')
# effluent of WWmixer goes back to WWTP

# =============================================================================
# facilities
# =============================================================================

HXN = qsu.HeatExchangerNetwork('HXN')

for unit in (SluL, SluT, SluC, P1, H1, HTL, HTL_hx, HTL_drum, H2SO4_Tank, AcidEx,
             M1, StruPre, CHG, CHG_pump, CHG_heating, CHG_cooling, V1, F1, MemDis, SP1,
             P3, HT, HT_compressor, HT_hx, V2, H4, F2, V3, SP2, H5, D1, D2, D3, P4,
             HC, HC_compressor, HC_hx, H6, V4, F3, D4, GasolineMixer, DieselMixer,
             H7, H8, PC1, PC2, GasolineTank, DieselTank, FuelMixer,
             GasMixer, CHP, WWmixer, RSP1, HXN):
    unit.register_alias(f'{unit=}'.split('=')[0].split('.')[-1])
# so that qs.main_flowsheet.H1 works as well

sys = qs.System('sys', path=(SluL, SluT, SluC, P1, H1, HTL, H2SO4_Tank, AcidEx,
                             M1, StruPre, CHG, V1, F1, MemDis, SP1,
                             P3, HT, V2, H4, F2, V3, SP2, H5, D1, D2, D3, P4,
                             HC, H6, V4, F3, D4, GasolineMixer, DieselMixer,
                             H7, H8, PC1, PC2, GasolineTank, DieselTank, FuelMixer,
                             GasMixer, CHP, WWmixer, RSP1
                             ), facilities=(HXN,))

sys.operating_hours = SluL.operation_hour # 7920 hr Jones

sys.simulate()

sys.diagram()

tea = create_tea(sys)

table = capex_table(tea)

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
# @metric(name='Struvite',units='kg/hr',element='Production')
# def get_struvite_production():
#     return StruPre.outs[0].imass['Struvite']

# @metric(name='(NH4)2SO4',units='kg/hr',element='Production')
# def get_nh42so4_production():
#     return MemDis.outs[0].imass['NH42SO4']

# @metric(name='Gasoline',units='kg/hr',element='Production')
# def get_gasoline_production():
#     return GasolineTank.outs[0].F_mass

# @metric(name='Diesel',units='kg/hr',element='Production')
# def get_diesel_production():
#     return DieselTank.outs[0].F_mass

# @metric(name='Total installed cost',units='$',element='TEA')
# def get_total_installed_cost():
#     return sys.installed_equipment_cost

@metric(name='MFSP',units='$',element='TEA')
def get_MFSP():
    return tea.solve_price(FuelMixer.outs[0])*3.2245 # from $/kg to $/gal diesel

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