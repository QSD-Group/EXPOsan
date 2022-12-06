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
from biosteam.units import IsenthalpicValve
from exposan.htl._process_settings import load_process_settings
from exposan.htl._components import create_components
from exposan.htl._TEA import *
from qsdsan import PowerUtility

# __all__ = ('create_system',)

# def create_system():

load_process_settings()
cmps = create_components()

# Construction here, StreamImpactItem after TEA
qs.ImpactIndicator.load_from_file('/Users/jiananfeng/Desktop/PhD CEE/coding/Cloned packages/EXPOsan/exposan/htl/data/impact_indicators.csv')

qs.ImpactItem.load_from_file('/Users/jiananfeng/Desktop/PhD CEE/coding/Cloned packages/EXPOsan/exposan/htl/data/impact_items.xlsx')

raw_wastewater = qs.Stream('raw_wastewater', H2O=20, units='MGD', T=25+273.15)
# Jones baseline: 1276.6 MGD, 1.066e-4 $/kg ww
# set H2O equal to the total raw wastewater into the WWTP

# =============================================================================
# pretreatment (Area 000)
# =============================================================================

WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                    ww_2_dry_sludge=0.94,
                    # how much metric ton/day sludge can be produced by 1 MGD of ww
                    sludge_moisture=0.99, sludge_dw_ash=0.257, 
                    sludge_afdw_lipid=0.204, sludge_afdw_protein=0.463, yearly_operation_hour=7920)

SluC = su.HTL_sludge_centrifuge('A000', ins=WWTP-0,
                            outs=('supernatant','compressed_sludge'),
                            init_with='Stream',
                            solids=('Sludge_lipid','Sludge_protein',
                                    'Sludge_carbo','Sludge_ash'),
                            sludge_moisture=0.8)

# =============================================================================
# HTL (Area 100)
# =============================================================================

P1 = su.HTLpump('A100', ins=SluC-1, outs='press_sludge', P=3049.7*6894.76,
              init_with='Stream')
# Jones 2014: 3049.7 psia

H1 = su.HTLHX('A110', ins=P1-0, outs='heated_sludge', T=351+273.15,
                   U=0.0794957, init_with='Stream')
# feed T is low, thus high viscosity and low U (case B in Knorr 2013)
# U: 3, 14, 15 BTU/hr/ft2/F as minimum, baseline, and maximum
# U: 0.0170348, 0.0794957, 0.085174 kW/m2/K
# H1: SS PNNL 2020: 50 (17-76) Btu/hr/ft2/F ~ U = 0.284 (0.096-0.4313) kW/m2/K
# but not in other pumps (low viscosity, don't need U to enforce total heat transfer efficiency)
# unit conversion: https://www.unitsconverters.com/en/Btu(It)/Hmft2mdegf-To-W/M2mk/Utu-4404-4398

HTL = su.HTL('A120', ins=H1-0, outs=('biochar','HTL_aqueous',
             'biocrude','offgas_HTL'))
HTL_hx = HTL.heat_exchanger
HTL_drum = HTL.kodrum

# =============================================================================
# CHG (Area 200)
# =============================================================================

H2SO4_Tank = su.HTL_storage_tank('T200', ins='H2SO4', outs=('H2SO4_out'),
                             init_with='WasteStream', tau=24, vessel_material='Carbon steel')
H2SO4_Tank.ins[0].price = 0.00658 # based on 93% H2SO4 and fresh water (dilute to 5%) price found in Davis 2020$/kg

SP1 = su.HTLsplitter('S200',ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                     init_with='Stream')
# must put after AcidEx and MemDis in path during simulation to ensure input
# not empty

AcidEx = su.AcidExtraction('A200', ins=(HTL-0, SP1-0),
                           outs=('residual','extracted'))

AcidEx.outs[0].price = -0.055 # SS 2021 SOT PNNL report page 24 Table 9

M1 = su.HTLmixer('A210', ins=(HTL-1, AcidEx-1), outs=('mixture'))

StruPre = su.StruvitePrecipitation('A220', ins=(M1-0,'MgCl2','NH4Cl','MgO'),
                                   outs=('struvite','CHG_feed'))
StruPre.ins[1].price = 0.5452
StruPre.ins[2].price = 0.13
StruPre.ins[3].price = 0.2
StruPre.outs[0].price = 0.661

CHG = su.CHG('A230', ins=(StruPre-1, '7.8% Ru/C'), outs=('CHG_out', '7.8% Ru/C_out'))
CHG_pump = CHG.pump
CHG_heating = CHG.heat_ex_heating
CHG_cooling = CHG.heat_ex_cooling
CHG.ins[1].price = 134.53

V1 = IsenthalpicValve('A240', ins=CHG-0, outs='depressed_cooled_CHG', P=50*6894.76)

F1 = su.HTLflash('A250', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                 T=60+273.15, P=50*6894.76)

MemDis = su.MembraneDistillation('A260', ins=(F1-1, SP1-1, 'NaOH', 'Membrane_in'),
                                  outs=('ammonium_sulfate','MemDis_ww', 'Membrane_out','solution'), init_with='WasteStream')
MemDis.ins[2].price = 0.5256
MemDis.outs[0].price = 0.3236

# =============================================================================
# HT (Area 300)
# =============================================================================

P2 = su.HTLpump('A300', ins=HTL-2, outs='press_biocrude', P=1530.0*6894.76,
              init_with='Stream')
# Jones 2014: 1530.0 psia

# Tin = 174 C (345 F) based on Jones PNNL report. However, the reaction
# releases a lot of heat and increase the temperature of effluent to 402 C
# (755.5 F).

RSP1 = qsu.ReversedSplitter('S300', ins='H2', outs=('HT_H2','HC_H2'),
                            init_with='WasteStream')
# reversed splitter, write before HT and HC, simulate after HT and HC
RSP1.ins[0].price = 1.61

HT = su.HT('A310', ins=(P2-0, RSP1-0, 'CoMo_alumina_HT'), outs=('HTout', 'CoMo_alumina_HT_out'))
HT_compressor = HT.compressor
HT_hx = HT.heat_exchanger
HT.ins[2].price = 38.79

V2 = IsenthalpicValve('A320', ins=HT-0, outs='depressed_HT', P=717.4*6894.76)

H2 = su.HTLHX('A330', ins=V2-0, outs='cooled_HT', T=60+273.15,
                    init_with='Stream')

F2 = su.HTLflash('A340', ins=H2-0, outs=('HT_fuel_gas','HT_aqueous'), T=43+273.15,
                 P=717.4*6894.76) # outflow P

V3 = IsenthalpicValve('A350', ins=F2-1, outs='depressed_flash_effluent', P=55*6894.76)

SP2 = qsu.Splitter('S310', ins=V3-0, outs=('HT_ww','HT_oil'),
                    split={'H2O':1}, init_with='Stream')
# separate water and oil based on gravity

H3 = su.HTLHX('A360', ins=SP2-1, outs='heated_oil', T=104+273.15)
# temperature: Jones stream #334 (we remove the first distillation column)

D1 = su.HTLdistillation('A370', ins=H3-0,
                        outs=('HT_light','HT_heavy'),
                        LHK=('C4H10','TWOMBUTAN'), P=50*6894.76, # outflow P
                        y_top=188/253, x_bot=53/162, k=2, is_divided=True)

D2 = su.HTLdistillation('A380', ins=D1-1,
                        outs=('HT_Gasoline','HT_other_oil'),
                        LHK=('C10H22','C4BENZ'), P=25*6894.76, # outflow P
                        y_top=116/122, x_bot=114/732, k=2, is_divided=True)

D3 = su.HTLdistillation('A390', ins=D2-1,
                        outs=('HT_Diesel','HT_heavy_oil'),
                        LHK=('C19H40','C21H44'),P=18.7*6894.76, # outflow P
                        y_top=2421/2448, x_bot=158/2448, k=2, is_divided=True)

# =============================================================================
# HC (Area 400)
# =============================================================================

P3 = su.HTLpump('A400', ins=D3-1, outs='press_heavy_oil', P=1034.7*6894.76,
              init_with='Stream')
# Jones 2014: 1034.7 psia

# Tin = 394 C (741.2 F) based on Jones PNNL report. However, the reaction
# releases a lot of heat and increase the temperature of effluent to 451 C
# (844.6 F).

HC = su.HC('A410', ins=(P3-0, RSP1-1, 'CoMo_alumina_HC'), outs=('HC_out', 'CoMo_alumina_HC_out'))
HC_compressor = HC.compressor
HC_hx = HC.heat_exchanger
HC.ins[2].price = 38.79

H4 = su.HTLHX('A420', ins=HC-0, outs='cooled_HC', T=60+273.15,
                    init_with='Stream')

V4 = IsenthalpicValve('A430', ins=H4-0, outs='cooled_depressed_HC', P=30*6894.76)


F3 = su.HTLflash('A440', ins=V4-0, outs=('HC_fuel_gas','HC_aqueous'), T=60.2+273,
                 P=30*6894.76) # outflow P

D4 = su.HTLdistillation('A450', ins=F3-1, outs=('HC_Gasoline','HC_Diesel'),
                        LHK=('C9H20','C10H22'), P=20*6894.76, # outflow P
                        y_top=360/546, x_bot=7/708, k=2, is_divided=True)

# =============================================================================
# CHP and storage (Area 500)
# =============================================================================

GasolineMixer = qsu.Mixer('S500', ins=(D2-0, D4-0), outs='mixed_gasoline',
                          init_with='Stream')

DieselMixer = qsu.Mixer('S510', ins=(D3-0, D4-1), outs='mixed_diesel',
                        init_with='Stream')

H5 = su.HTLHX('A500', ins=GasolineMixer-0, outs='cooled_gasoline',
                    T=60+273.15, init_with='Stream', rigorous=True)

H6 = su.HTLHX('A510', ins=DieselMixer-0, outs='cooled_diesel',
                    T=60+273.15, init_with='Stream', rigorous=True)

PC1 = su.PhaseChanger('S520', ins=H5-0, outs='cooled_gasoline_liquid')

PC2 = su.PhaseChanger('S530', ins=H6-0, outs='cooled_diesel_liquid')

GasolineTank = su.HTL_storage_tank('T500', ins=PC1-0, outs=('gasoline'),
                                tau=3*24, init_with='Stream', vessel_material='Carbon steel')
# store for 3 days based on Jones 2014

DieselTank = su.HTL_storage_tank('T510', ins=PC2-0, outs=('diesel'),
                              tau=3*24, init_with='Stream', vessel_material='Carbon steel')
# store for 3 days based on Jones 2014

FuelMixer = su.FuelMixer('S540', ins=(GasolineTank-0, DieselTank-0),\
                         outs='fuel', target='diesel')
# integrate gasoline and diesel based on their LHV for MFSP calculation

GasMixer = qsu.Mixer('S550', ins=(HTL-3, F1-0, F2-0, D1-0, F3-0),
                      outs=('fuel_gas'), init_with='Stream')

# The system produces more energy than needed (heating+power)
CHP = su.HTLCHP('A520', ins=(GasMixer-0,'natural_gas','air'),
              outs=('emission','solid_ash'), init_with='WasteStream', supplement_power_utility=False)

CHP.ins[1].price = 0.1685

WWmixer = su.WWmixer('S560', ins=(SluC-0, MemDis-1, SP2-0),
                    outs='wastewater', init_with='Stream')
# effluent of WWmixer goes back to WWTP

# =============================================================================
# facilities
# =============================================================================

HXN = qsu.HeatExchangerNetwork('HXN')

for unit in (WWTP, SluC, P1, H1, HTL, HTL_hx, HTL_drum, H2SO4_Tank, AcidEx,
             M1, StruPre, CHG, CHG_pump, CHG_heating, CHG_cooling, V1, F1, MemDis, SP1,
             P2, HT, HT_compressor, HT_hx, V2, H2, F2, V3, SP2, H3, D1, D2, D3, P3,
             HC, HC_compressor, HC_hx, H4, V4, F3, D4, GasolineMixer, DieselMixer,
             H5, H6, PC1, PC2, GasolineTank, DieselTank, FuelMixer,
             GasMixer, CHP, WWmixer, RSP1, HXN):
    unit.register_alias(f'{unit=}'.split('=')[0].split('.')[-1])
# so that qs.main_flowsheet.H1 works as well

sys = qs.System('sys', path=(WWTP, SluC, P1, H1, HTL, H2SO4_Tank, AcidEx,
                             M1, StruPre, CHG, V1, F1, MemDis, SP1,
                             P2, HT, V2, H2, F2, V3, SP2, H3, D1, D2, D3, P3,
                             HC, H4, V4, F3, D4, GasolineMixer, DieselMixer,
                             H5, H6, PC1, PC2, GasolineTank, DieselTank, FuelMixer,
                             GasMixer, CHP, WWmixer, RSP1
                              ), facilities=(HXN,))

sys.operating_hours = WWTP.operation_hour # 7920 hr Jones

sys.simulate()

sys.diagram()

tea = create_tea(sys)

table = capex_table(tea)

RO_item = qs.StreamImpactItem(linked_stream=MemDis.ins[3],
                              Acidification=0.53533,
                              Ecotoxicity=0.90848,
                              Eutrophication=0.0028322,
                              GlobalWarming=2.2663,
                              OzoneDepletion=0.00000025541,
                              PhotochemicalOxidation=0.0089068,
                              Carcinogenics=0.034791,
                              NonCarcinogenics=31.8,
                              RespiratoryEffects=0.0028778)
                              
H2SO4_item = qs.StreamImpactItem(linked_stream=H2SO4_Tank.ins[0],
                                 Acidification=0.019678922,
                                 Ecotoxicity=0.069909345,
                                 Eutrophication=4.05E-06,
                                 GlobalWarming=0.008205666,
                                 OzoneDepletion=8.94E-10,
                                 PhotochemicalOxidation=5.04E-05,
                                 Carcinogenics=1.74E-03,
                                 NonCarcinogenics=1.68237815,
                                 RespiratoryEffects=9.41E-05)

MgCl2_item = qs.StreamImpactItem(linked_stream=StruPre.ins[1],
                                 Acidification=0.77016,
                                 Ecotoxicity=0.97878,
                                 Eutrophication=0.00039767,
                                 GlobalWarming=2.8779,
                                 OzoneDepletion=4.94E-08,
                                 PhotochemicalOxidation=0.0072306,
                                 Carcinogenics=0.0050938,
                                 NonCarcinogenics=8.6916,
                                 RespiratoryEffects=0.004385)

H2_item = qs.StreamImpactItem(linked_stream=RSP1.ins[0],
                              Acidification=0.81014,
                              Ecotoxicity=0.42747,
                              Eutrophication=0.0029415,
                              GlobalWarming=1.5624,
                              OzoneDepletion=1.80E-06,
                              PhotochemicalOxidation=0.0052545,
                              Carcinogenics=0.0026274,
                              NonCarcinogenics=8.5687,
                              RespiratoryEffects=0.0036698)

MgO_item = qs.StreamImpactItem(linked_stream=StruPre.ins[3],
                               Acidification=0.12584,
                               Ecotoxicity=2.7949,
                               Eutrophication=0.00063607,
                               GlobalWarming=1.1606,
                               OzoneDepletion=1.54E-08,
                               PhotochemicalOxidation=0.0017137,
                               Carcinogenics=0.018607,
                               NonCarcinogenics=461.54,
                               RespiratoryEffects=0.0008755)

NaOH_item = qs.StreamImpactItem(linked_stream=MemDis.ins[2],
                                Acidification=0.33656,
                                Ecotoxicity=0.77272,
                                Eutrophication=0.00032908,
                                GlobalWarming=1.2514,
                                OzoneDepletion=7.89E-07,
                                PhotochemicalOxidation=0.0033971,
                                Carcinogenics=0.0070044,
                                NonCarcinogenics=13.228,
                                RespiratoryEffects=0.0024543)

NH4Cl_item = qs.StreamImpactItem(linked_stream=StruPre.ins[2],
                                 Acidification=0.34682,
                                 Ecotoxicity=0.90305, 
                                 Eutrophication=0.0047381,
                                 GlobalWarming=1.525,
                                 OzoneDepletion=9.22E-08,
                                 PhotochemicalOxidation=0.0030017,
                                 Carcinogenics=0.010029,
                                 NonCarcinogenics=14.85,
                                 RespiratoryEffects=0.0018387)

struvite_item = qs.StreamImpactItem(linked_stream=StruPre.outs[0],
                                    Acidification=-0.122829597,
                                    Ecotoxicity=-0.269606396,
                                    Eutrophication=-0.000174952,
                                    GlobalWarming=-0.420850152,
                                    OzoneDepletion=-2.29549E-08,
                                    PhotochemicalOxidation=-0.001044087,
                                    Carcinogenics=-0.002983018,
                                    NonCarcinogenics=-4.496533528,
                                    RespiratoryEffects=-0.00061764)

NH42SO4_item = qs.StreamImpactItem(linked_stream=MemDis.outs[0],
                                   Acidification=-0.72917,
                                   Ecotoxicity=-3.4746,
                                   Eutrophication=-0.0024633,
                                   GlobalWarming=-1.2499,
                                   OzoneDepletion=-6.12E-08,
                                   PhotochemicalOxidation=-0.0044519,
                                   Carcinogenics=-0.036742,
                                   NonCarcinogenics=-62.932,
                                   RespiratoryEffects=-0.0031315)

Natural_gas_item = qs.StreamImpactItem(linked_stream=CHP.ins[1],
                                       Acidification=0.1032,
                                       Ecotoxicity=0.1071,
                                       Eutrophication=7.87E-05,
                                       GlobalWarming=2.754919474,
                                       OzoneDepletion=2.17E-07,
                                       PhotochemicalOxidation=0.00092588,
                                       Carcinogenics=0.00088329,
                                       NonCarcinogenics=4.3413,
                                       RespiratoryEffects=0.00048979)

CHG_catalyst_item = qs.StreamImpactItem(linked_stream=CHG.ins[1],
                                        Acidification=991.6544196,
                                        Ecotoxicity=15371.08292,
                                        Eutrophication=0.45019348,
                                        GlobalWarming=484.7862509,
                                        OzoneDepletion=2.23437E-05,
                                        PhotochemicalOxidation=6.735405072,
                                        Carcinogenics=1.616793132,
                                        NonCarcinogenics=27306.37232,
                                        RespiratoryEffects=3.517184526)

HT_catalyst_item = qs.StreamImpactItem(linked_stream=HT.ins[2],
                                        Acidification=4.056401283,
                                        Ecotoxicity=50.26926274,
                                        Eutrophication=0.005759274,
                                        GlobalWarming=6.375878231,
                                        OzoneDepletion=1.39248E-06,
                                        PhotochemicalOxidation=0.029648759,
                                        Carcinogenics=0.287516945,
                                        NonCarcinogenics=369.791688,
                                        RespiratoryEffects=0.020809293)

HC_catalyst_item = qs.StreamImpactItem(linked_stream=HC.ins[2],
                                        Acidification=4.056401283,
                                        Ecotoxicity=50.26926274,
                                        Eutrophication=0.005759274,
                                        GlobalWarming=6.375878231,
                                        OzoneDepletion=1.39248E-06,
                                        PhotochemicalOxidation=0.029648759,
                                        Carcinogenics=0.287516945,
                                        NonCarcinogenics=369.791688,
                                        RespiratoryEffects=0.020809293)

lca = qs.LCA(system=sys, lifetime=30, lifetime_unit='yr', Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30)

# return sys

#%%
from qsdsan import Model

model = Model(sys)

from chaospy import distributions as shape

param = model.parameter

# =============================================================================
# WWTP
# =============================================================================
dist = shape.Uniform(0.846,1.034)
@param(name='ww_2_dry_sludge',
        element=WWTP,
        kind='coupled',
        units='ton/d/MGD',
        baseline=0.94,
        distribution=dist)
def set_ww_2_dry_sludge(i):
    WWTP.ww_2_dry_sludge=i

dist = shape.Uniform(0.97,0.995)
@param(name='sludge_moisture',
        element=WWTP,
        kind='coupled',
        units='-',
        baseline=0.99,
        distribution=dist)
def set_WWTP_sludge_moisture(i):
    WWTP.sludge_moisture=i

dist = shape.Triangle(0.174,0.2567,0.414)
@param(name='sludge_dw_ash',
        element=WWTP,
        kind='coupled',
        units='-',
        baseline=0.2567,
        distribution=dist)
def set_sludge_dw_ash(i):
    WWTP.sludge_dw_ash=i

dist = shape.Triangle(0.38,0.4634,0.51)
@param(name='sludge_afdw_protein',
        element=WWTP,
        kind='coupled',
        units='-',
        baseline=0.4634,
        distribution=dist)
def set_sludge_afdw_protein(i):
    WWTP.sludge_afdw_protein=i

dist = shape.Triangle(0.08,0.2043,0.308)
@param(name='sludge_afdw_lipid',
        element=WWTP,
        kind='coupled',
        units='-',
        baseline=0.2043,
        distribution=dist)
def set_sludge_afdw_lipid(i):
    WWTP.sludge_afdw_lipid=i

dist = shape.Uniform(0.675,0.825)
@param(name='lipid_2_C',
        element=WWTP,
        kind='coupled',
        units='-',
        baseline=0.75,
        distribution=dist)
def set_lipid_2_C(i):
    WWTP.lipid_2_C=i

dist = shape.Uniform(0.4905,0.5995)
@param(name='protein_2_C',
        element=WWTP,
        kind='coupled',
        units='-',
        baseline=0.545,
        distribution=dist)
def set_protein_2_C(i):
    WWTP.protein_2_C=i

dist = shape.Uniform(0.36,0.44)
@param(name='carbo_2_C',
        element=WWTP,
        kind='coupled',
        units='-',
        baseline=0.4,
        distribution=dist)
def set_carbo_2_C(i):
    WWTP.carbo_2_C=i

dist = shape.Triangle(0.1348,0.1427,0.1647)
@param(name='C_2_H',
        element=WWTP,
        kind='coupled',
        units='-',
        baseline=0.1427,
        distribution=dist)
def set_C_2_H(i):
    WWTP.C_2_H=i

dist = shape.Uniform(0.1431,0.1749)
@param(name='protein_2_N',
        element=WWTP,
        kind='coupled',
        units='-',
        baseline=0.159,
        distribution=dist)
def set_protein_2_N(i):
    WWTP.protein_2_N=i
    
dist = shape.Triangle(0.1944,0.3927,0.5556)
@param(name='N_2_P',
        element=WWTP,
        kind='coupled',
        units='-',
        baseline=0.3927,
        distribution=dist)
def set_N_2_P(i):
    WWTP.N_2_P=i

dist = shape.Triangle(7392,7920,8448)
@param(name='operation_hour',
        element=WWTP,
        kind='coupled',
        units='hr/yr',
        baseline=7920,
        distribution=dist)
def set_operation_hour(i):
    WWTP.operation_hour=i

# =============================================================================
# HTL
# =============================================================================
dist = shape.Triangle(0.017035,0.0795,0.085174)
@param(name='enforced heating transfer coefficient',
        element=H1,
        kind='coupled',
        units='kW/m2/K',
        baseline=0.0795,
        distribution=dist)
def set_U(i):
    H1.U=i

dist = shape.Triangle(0.692,0.846,1)
@param(name='lipid_2_biocrude',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.846,
        distribution=dist)
def set_lipid_2_biocrude(i):
    HTL.lipid_2_biocrude=i

dist = shape.Normal(0.445,0.030)
@param(name='protein_2_biocrude',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.445,
        distribution=dist)
def set_protein_2_biocrude(i):
    HTL.protein_2_biocrude=i

dist = shape.Normal(0.205,0.050)
@param(name='carbo_2_biocrude',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.205,
        distribution=dist)
def set_carbo_2_biocrude(i):
    HTL.carbo_2_biocrude=i

dist = shape.Normal(0.074,0.020)
@param(name='protein_2_gas',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.074,
        distribution=dist)
def set_protein_2_gas(i):
    HTL.protein_2_gas=i

dist = shape.Normal(0.418,0.030)
@param(name='carbo_2_gas',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.418,
        distribution=dist)
def set_carbo_2_gas(i):
    HTL.carbo_2_gas=i
    
dist = shape.Normal(-8.370,0.939)
@param(name='biocrude_C_slope',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=-8.370,
        distribution=dist)
def set_biocrude_C_slope(i):
    HTL.biocrude_C_slope=i
    
dist = shape.Normal(68.55,0.367)
@param(name='biocrude_C_intercept',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=68.55,
        distribution=dist)
def set_biocrude_C_intercept(i):
    HTL.biocrude_C_intercept=i
    
dist = shape.Normal(0.133,0.005)
@param(name='biocrude_N_slope',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.133,
        distribution=dist)
def set_biocrude_N_slope(i):
    HTL.biocrude_N_slope=i
    
dist = shape.Normal(-2.610,0.352)
@param(name='biocrude_H_slope',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=-2.610,
        distribution=dist)
def set_biocrude_H_slope(i):
    HTL.biocrude_H_slope=i

dist = shape.Normal(8.200,0.138)
@param(name='biocrude_H_intercept',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=8.200,
        distribution=dist)
def set_biocrude_H_intercept(i):
    HTL.biocrude_H_intercept=i
    
dist = shape.Normal(478,18.878)
@param(name='HTLaqueous_C_slope',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=478,
        distribution=dist)
def set_HTLaqueous_C_slope(i):
    HTL.HTLaqueous_C_slope=i
    
dist = shape.Triangle(0.715,0.764,0.813)
@param(name='TOC_TC',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.764,
        distribution=dist)
def set_TOC_TC(i):
    HTL.TOC_TC=i

dist = shape.Normal(1.750,0.122)
@param(name='biochar_C_slope',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=1.750,
        distribution=dist)
def set_biochar_C_slope(i):
    HTL.biochar_C_slope=i

dist = shape.Triangle(0.035,0.063125,0.102)
@param(name='biocrude_moisture_content',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.063125,
        distribution=dist)
def set_biocrude_moisture_content(i):
    HTL.biocrude_moisture_content=i

dist = shape.Uniform(0.84,0.88)
@param(name='biochar_P_recovery_ratio',
        element=HTL,
        kind='coupled',
        units='-',
        baseline=0.86,
        distribution=dist)
def set_biochar_P_recovery_ratio(i):
    HTL.biochar_P_recovery_ratio=i

# =============================================================================
# AcidEx
# =============================================================================
dist = shape.Uniform(4,10)
@param(name='acid_vol',
        element=AcidEx,
        kind='coupled',
        units='-',
        baseline=7,
        distribution=dist)
def set_acid_vol(i):
    AcidEx.acid_vol=i

dist = shape.Uniform(0.7,0.9)
@param(name='P_recovery_ratio',
        element=AcidEx,
        kind='coupled',
        units='-',
        baseline=0.8,
        distribution=dist)
def set_P_recovery_ratio(i):
    AcidEx.P_recovery_ratio=i

# =============================================================================
# StruPre
# =============================================================================
dist = shape.Uniform(8.5,9.5)
@param(name='target_pH',
        element=StruPre,
        kind='coupled',
        units='-',
        baseline=9,
        distribution=dist)
def set_StruPre_target_pH(i):
    StruPre.target_pH=i

dist = shape.Uniform(1,4)
@param(name='Mg_P_ratio',
        element=StruPre,
        kind='coupled',
        units='-',
        baseline=2.5,
        distribution=dist)
def set_Mg_P_ratio(i):
    StruPre.Mg_P_ratio=i
    
dist = shape.Triangle(0.7,0.828,0.95)
@param(name='P_pre_recovery_ratio',
        element=StruPre,
        kind='coupled',
        units='-',
        baseline=0.828,
        distribution=dist)
def set_P_pre_recovery_ratio(i):
    StruPre.P_pre_recovery_ratio=i

# =============================================================================
# CHG
# =============================================================================
dist = shape.Triangle(2.86,3.562,3.99)
@param(name='WHSV',
        element=CHG,
        kind='coupled',
        units='kg/hr/kg',
        baseline=3.562,
        distribution=dist)
def set_CHG_WHSV(i):
    CHG.WHSV=i

dist = shape.Triangle(3960,7920,15840)
@param(name='catalyst_lifetime',
        element=CHG,
        kind='coupled',
        units='hr',
        baseline=7920,
        distribution=dist)
def set_CHG_catalyst_lifetime(i):
    CHG.catalyst_lifetime=i

dist = shape.Triangle(0.1893,0.5981,0.7798)
@param(name='gas_C_to_total_C',
        element=CHG,
        kind='coupled',
        units='-',
        baseline=0.5981,
        distribution=dist)
def set_gas_C_to_total_C(i):
    CHG.gas_C_to_total_C=i

# =============================================================================
# MemDis
# =============================================================================
dist = shape.Uniform(7.91,8.41)
@param(name='influent_pH',
        element=MemDis,
        kind='coupled',
        units='-',
        baseline=8.16,
        distribution=dist)
def set_influent_pH(i):
    MemDis.influent_pH=i

dist = shape.Uniform(10,11.8)
@param(name='target_pH',
        element=MemDis,
        kind='coupled',
        units='-',
        baseline=10,
        distribution=dist)
def set_MemDis_target_pH(i):
    MemDis.target_pH=i

dist = shape.Uniform(0.00075,0.000917)
@param(name='m2_2_m3',
        element=MemDis,
        kind='coupled',
        units='-',
        baseline=0.0008333,
        distribution=dist)
def set_m2_2_m3(i):
    MemDis.m2_2_m3=i

dist = shape.Uniform(0.0000205,0.0000251)
@param(name='Dm',
        element=MemDis,
        kind='coupled',
        units='m2/s',
        baseline=0.0000228,
        distribution=dist)
def set_Dm(i):
    MemDis.Dm=i

dist = shape.Uniform(0.81,0.99)
@param(name='porosity',
        element=MemDis,
        kind='coupled',
        units='-',
        baseline=0.9,
        distribution=dist)
def set_porosity(i):
    MemDis.porosity=i

dist = shape.Uniform(0.000063,0.000077)
@param(name='thickness',
        element=MemDis,
        kind='coupled',
        units='m',
        baseline=0.00007,
        distribution=dist)
def set_thickness(i):
    MemDis.thickness=i

dist = shape.Uniform(1.08,1.32)
@param(name='tortuosity',
        element=MemDis,
        kind='coupled',
        units='-',
        baseline=1.2,
        distribution=dist)
def set_tortuosity(i):
    MemDis.tortuosity=i

dist = shape.Uniform(0.0000158,0.0000193)
@param(name='Ka',
        element=MemDis,
        kind='coupled',
        units='-',
        baseline=0.0000175,
        distribution=dist)
def set_Ka(i):
    MemDis.Ka=i

dist = shape.Uniform(5.409,6.611)
@param(name='capacity',
        element=MemDis,
        kind='coupled',
        units='-',
        baseline=6.01,
        distribution=dist)
def set_capacity(i):
    MemDis.capacity=i

# =============================================================================
# HT
# =============================================================================
dist = shape.Uniform(0.5625,0.6875)
@param(name='WHSV',
        element=HT,
        kind='coupled',
        units='kg/hr/kg',
        baseline=0.625,
        distribution=dist)
def set_HT_WHSV(i):
    HT.WHSV=i

dist = shape.Triangle(7920,15840,39600)
@param(name='catalyst_lifetime',
        element=HT,
        kind='coupled',
        units='hr',
        baseline=15840,
        distribution=dist)
def set_HT_catalyst_lifetime(i):
    HT.catalyst_lifetime=i

dist = shape.Uniform(0.0414,0.0506)
@param(name='hydrogen_rxned_to_biocrude',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.046,
        distribution=dist)
def set_HT_hydrogen_rxned_to_biocrude(i):
    HT.hydrogen_rxned_to_biocrude=i

dist = shape.Uniform(2.4,3.6)
@param(name='hydrogen_excess',
        element=HT,
        kind='coupled',
        units='-',
        baseline=3,
        distribution=dist)
def set_HT_hydrogen_excess(i):
    HT.hydrogen_excess=i

dist = shape.Uniform(0.7875,0.9625)
@param(name='hydrocarbon_ratio',
        element=HT,
        kind='coupled',
        units='-',
        baseline=0.875,
        distribution=dist)
def set_HT_hydrocarbon_ratio(i):
    HT.hydrocarbon_ratio=i

# =============================================================================
# HC
# =============================================================================
dist = shape.Uniform(0.5625,0.6875)
@param(name='WHSV',
        element=HC,
        kind='coupled',
        units='kg/hr/kg',
        baseline=0.625,
        distribution=dist)
def set_HC_WHSV(i):
    HC.WHSV=i

dist = shape.Uniform(35640,43560)
@param(name='catalyst_lifetime',
        element=HC,
        kind='coupled',
        units='hr',
        baseline=39600,
        distribution=dist)
def set_HC_catalyst_lifetime(i):
    HC.catalyst_lifetime=i

dist = shape.Uniform(0.010125,0.012375)
@param(name='hydrogen_rxned_to_heavy_oil',
        element=HC,
        kind='coupled',
        units='-',
        baseline=0.01125,
        distribution=dist)
def set_HC_hydrogen_rxned_to_heavy_oil(i):
    HC.hydrogen_rxned_to_heavy_oil=i

dist = shape.Uniform(4.4448,6.6672)
@param(name='hydrogen_excess',
        element=HC,
        kind='coupled',
        units='-',
        baseline=5.556,
        distribution=dist)
def set_HC_hydrogen_excess(i):
    HC.hydrogen_excess=i

dist = shape.Uniform(0.9,1)
@param(name='hydrocarbon_ratio',
        element=HC,
        kind='coupled',
        units='-',
        baseline=1,
        distribution=dist)
def set_HC_hydrocarbon_ratio(i):
    HC.hydrocarbon_ratio=i

# =============================================================================
# TEA
# =============================================================================

dist = shape.Triangle(0.6,1,1.4)
@param(name='CAPEX_factor',
        element='TEA',
        kind='isolated',
        units='-',
        baseline=1,
        distribution=dist)
def set_CAPEX_factor(i):
    HTL.CAPEX_factor=CHG.CAPEX_factor=HT.CAPEX_factor=i


dist = shape.Uniform(1102.5,1347.5)
@param(name='unit_CAPEX',
        element='TEA',
        kind='isolated',
        units='-',
        baseline=1225,
        distribution=dist)
def set_unit_CAPEX(i):
    CHP.unit_CAPEX=i

dist = shape.Triangle(0,0.1,0.2)
@param(name='IRR',
        element='TEA',
        kind='isolated',
        units='-',
        baseline=0.1,
        distribution=dist)
def set_IRR(i):
    tea.IRR=i

dist = shape.Triangle(0.005994,0.00658,0.014497)
@param(name='5% H2SO4 price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=0.00658,
        distribution=dist)
def set_H2SO4_price(i):
    H2SO4_Tank.ins[0].price=i

dist = shape.Triangle(0.525,0.5452,0.57)
@param(name='MgCl2 price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=0.5452,
        distribution=dist)
def set_MgCl2_price(i):
    StruPre.ins[1].price=i

dist = shape.Uniform(0.12,0.14)
@param(name='NH4Cl price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=0.13,
        distribution=dist)
def set_NH4Cl_price(i):
    StruPre.ins[2].price=i

dist = shape.Uniform(0.1,0.3)
@param(name='MgO price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=0.2,
        distribution=dist)
def set_MgO_price(i):
    StruPre.ins[3].price=i
        
        
dist = shape.Triangle(0.419,0.661,1.213)
@param(name='struvite price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=0.661,
        distribution=dist)
def set_struvite_price(i):
    StruPre.outs[0].price=i        
        
dist = shape.Uniform(1.450,1.772)
@param(name='H2 price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=1.611,
        distribution=dist)
def set_H2_price(i):
    RSP1.ins[0].price=i

dist = shape.Uniform(0.473,0.578)
@param(name='NaOH price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=0.5256,
        distribution=dist)
def set_NaOH_price(i):
    MemDis.ins[2].price=i

dist = shape.Triangle(0.1636,0.3236,0.463)
@param(name='ammonium sulfate price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=0.3236,
        distribution=dist)
def set_ammonium_sulfate_price(i):
    MemDis.outs[0].price=i

dist = shape.Uniform(83.96,102.62)
@param(name='membrane price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=93.29,
        distribution=dist)
def set_membrane_price(i):
    MemDis.membrane_price=i

dist = shape.Triangle(67.27,134.53,269.07)
@param(name='CHG catalyst price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=134.53,
        distribution=dist)
def set_catalyst_price(i):
    CHG.ins[1].price=i

dist = shape.Triangle(0.121,0.1685,0.3608)
@param(name='CH4 price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=0.1685,
        distribution=dist)
def set_CH4_price(i):
    CHP.ins[1].price=i

dist = shape.Uniform(34.91,42.67)
@param(name='HT & HC catalyst price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=38.79,
        distribution=dist)
def set_HT_HC_catalyst_price(i):
    HT.ins[2].price=HC.ins[2].price=i
    
dist = shape.Triangle(0.7085,0.9388,1.4493)
@param(name='gasoline price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=0.9388,
        distribution=dist)
def set_gasoline_price(i):
    FuelMixer.gasoline_price=i
    
dist = shape.Triangle(0.7458,0.9722,1.6579)
@param(name='diesel price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=0.9722,
        distribution=dist)
def set_diesel_price(i):
    FuelMixer.diesel_price=i
    
dist = shape.Uniform(-0.0605,-0.0495)
@param(name='residual disposal',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=-0.055,
        distribution=dist)
def set_residual_disposal(i):
    AcidEx.outs[0].price =i

dist = shape.Triangle(0.0667,0.06879,0.07180)
@param(name='electrivity price',
        element='TEA',
        kind='isolated',
        units='$/kg',
        baseline=0.06879,
        distribution=dist)
def set_electrivity_price(i):
    PowerUtility.price=i

metric = model.metric
@metric(name='MFSP',units='$',element='TEA')
def get_MFSP():
    return tea.solve_price(FuelMixer.outs[0])*3.220628346 # from $/kg to $/gal diesel

#%%
import numpy as np
np.random.seed(3221)
samples = model.sample(N=5, rule='L')
model.load_samples(samples)
model.evaluate()
model.table
#%%
import pandas as pd
def organize_results(model, path):
    idx = len(model.parameters)
    parameters = model.table.iloc[:, :idx]
    results = model.table.iloc[:, idx:]
    percentiles = results.quantile([0, 0.05, 0.25, 0.5, 0.75, 0.95, 1])
    with pd.ExcelWriter(path) as writer:
        parameters.to_excel(writer, sheet_name='Parameters')
        results.to_excel(writer, sheet_name='Results')
        percentiles.to_excel(writer, sheet_name='Percentiles')
organize_results(model, 'example_model.xlsx')
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