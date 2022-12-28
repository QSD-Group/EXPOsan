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

import os, qsdsan as qs
import exposan.htl._sanunits as su
from qsdsan import sanunits as qsu
from biosteam.units import IsenthalpicValve
from exposan.htl._process_settings import load_process_settings
from exposan.htl._components import create_components
from exposan.htl._tea import *
from qsdsan import PowerUtility, Model
from qsdsan.utils import auom, DictAttrSetter
from chaospy import distributions as shape

_m3perh_to_MGD = auom('m3/h').conversion_factor('MGD')

_kg_to_g = auom('kg').conversion_factor('g')

_MJ_to_MMBTU = auom('MJ').conversion_factor('MMBTU')

_MMgal_to_L = auom('gal').conversion_factor('L')*1000000

# __all__ = ('create_system',)

# def create_system():

load_process_settings()
cmps = create_components()

# Construction here, StreamImpactItem after TEA
folder = os.path.dirname(__file__)
qs.ImpactIndicator.load_from_file(os.path.join(folder, 'data/impact_indicators.csv'))
qs.ImpactItem.load_from_file(os.path.join(folder, 'data/impact_items.xlsx'))

raw_wastewater = qs.Stream('raw_wastewater', H2O=100, units='MGD', T=25+273.15)
# Jones baseline: 1276.6 MGD, 1.066e-4 $/kg ww
# set H2O equal to the total raw wastewater into the WWTP

# =============================================================================
# pretreatment (Area 000)
# =============================================================================
            
WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                    ww_2_dry_sludge=0.94,
                    # how much metric ton/day sludge can be produced by 1 MGD of ww
                    sludge_moisture=0.99, sludge_dw_ash=0.257, 
                    sludge_afdw_lipid=0.204, sludge_afdw_protein=0.463, operation_hours=7920)

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
                   U=0.0795, init_with='Stream', rigorous=True)
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
                             init_with='WasteStream', tau=24, vessel_material='Stainless steel')
H2SO4_Tank.ins[0].price = 0.00658 # based on 93% H2SO4 and fresh water (dilute to 5%) price found in Davis 2020$/kg

SP1 = su.HTLsplitter('S200',ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                     init_with='Stream')
# must put after AcidEx and MemDis in path during simulation to ensure input
# not empty

AcidEx = su.AcidExtraction('A200', ins=(HTL-0, SP1-0),
                           outs=('residual','extracted'))

# AcidEx.outs[0].price = -0.055 # SS 2021 SOT PNNL report page 24 Table 9
# not include residual for TEA and LCA for now

# baseline:
M1 = su.HTLmixer('A210', ins=(HTL-1, AcidEx-1), outs=('mixture'))

# if remove AcidEx:
# M1 = su.HTLmixer('A210', ins=(HTL-1, ''), outs=('mixture'))

StruPre = su.StruvitePrecipitation('A220', ins=(M1-0,'MgCl2','NH4Cl','MgO'),
                                   outs=('struvite','CHG_feed'))
StruPre.ins[1].price = 0.5452
StruPre.ins[2].price = 0.13
StruPre.ins[3].price = 0.2
StruPre.outs[0].price = 0.661

CHG = su.CHG('A230', ins=(StruPre-1, '7.8%_Ru/C'), outs=('CHG_out', '7.8%_Ru/C_out'))
CHG_pump = CHG.pump
CHG_heating = CHG.heat_ex_heating
CHG_cooling = CHG.heat_ex_cooling
CHG.ins[1].price = 134.53

V1 = IsenthalpicValve('A240', ins=CHG-0, outs='depressed_cooled_CHG', P=50*6894.76, vle=True)

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
HT_hx_H2 = HT.heat_exchanger_H2
HT_hx_oil = HT.heat_exchanger_oil
HT.ins[2].price = 38.79

V2 = IsenthalpicValve('A320', ins=HT-0, outs='depressed_HT', P=717.4*6894.76, vle=True)

H2 = su.HTLHX('A330', ins=V2-0, outs='cooled_HT', T=60+273.15,
                    init_with='Stream', rigorous=True)

F2 = su.HTLflash('A340', ins=H2-0, outs=('HT_fuel_gas','HT_aqueous'), T=43+273.15,
                 P=717.4*6894.76) # outflow P

V3 = IsenthalpicValve('A350', ins=F2-1, outs='depressed_flash_effluent', P=55*6894.76, vle=True)

SP2 = qsu.Splitter('S310', ins=V3-0, outs=('HT_ww','HT_oil'),
                    split={'H2O':1}, init_with='Stream')
# separate water and oil based on gravity

H3 = su.HTLHX('A360', ins=SP2-1, outs='heated_oil', T=104+273.15, rigorous=True)
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
HC_hx_H2 = HC.heat_exchanger_H2
HC_hx_oil = HC.heat_exchanger_oil
HC.ins[2].price = 38.79

H4 = su.HTLHX('A420', ins=HC-0, outs='cooled_HC', T=60+273.15,
                    init_with='Stream', rigorous=True)

V4 = IsenthalpicValve('A430', ins=H4-0, outs='cooled_depressed_HC', P=30*6894.76, vle=True)


F3 = su.HTLflash('A440', ins=V4-0, outs=('HC_fuel_gas','HC_aqueous'), T=60.2+273,
                 P=30*6894.76) # outflow P

D4 = su.HTLdistillation('A450', ins=F3-1, outs=('HC_Gasoline','HC_Diesel'),
                        LHK=('C9H20','C10H22'), P=20*6894.76, # outflow P
                        y_top=360/546, x_bot=7/708, k=2, is_divided=True)

# =============================================================================
# Storage, and disposal (Area 500)
# =============================================================================

GasolineMixer = qsu.Mixer('S500', ins=(D2-0, D4-0), outs='mixed_gasoline',
                          init_with='Stream', rigorous=True)

DieselMixer = qsu.Mixer('S510', ins=(D3-0, D4-1), outs='mixed_diesel',
                        init_with='Stream', rigorous=True)

H5 = su.HTLHX('A500', ins=GasolineMixer-0, outs='cooled_gasoline',
                    T=60+273.15, init_with='Stream', rigorous=True)

H6 = su.HTLHX('A510', ins=DieselMixer-0, outs='cooled_diesel',
                    T=60+273.15, init_with='Stream', rigorous=True)

PC1 = su.PhaseChanger('S520', ins=H5-0, outs='cooled_gasoline_liquid')

PC2 = su.PhaseChanger('S530', ins=H6-0, outs='cooled_diesel_liquid')

PC3 = su.PhaseChanger('S540', ins=CHG-1, outs='CHG_catalyst_out', phase='s')

PC4 = su.PhaseChanger('S550', ins=HT-1, outs='HT_catalyst_out', phase='s')

PC5 = su.PhaseChanger('S560', ins=HC-1, outs='HC_catalyst_out', phase='s')

GasolineTank = su.HTL_storage_tank('T500', ins=PC1-0, outs=('gasoline'),
                                tau=3*24, init_with='Stream', vessel_material='Carbon steel')
# store for 3 days based on Jones 2014

DieselTank = su.HTL_storage_tank('T510', ins=PC2-0, outs=('diesel'),
                              tau=3*24, init_with='Stream', vessel_material='Carbon steel')
# store for 3 days based on Jones 2014

FuelMixer = su.FuelMixer('S570', ins=(GasolineTank-0, DieselTank-0),
                         outs='fuel', target='diesel')
# integrate gasoline and diesel based on their LHV for MFSP calculation

GasMixer = qsu.Mixer('S580', ins=(HTL-3, F1-0, F2-0, D1-0, F3-0),
                      outs=('fuel_gas'), init_with='Stream')

WWmixer = su.WWmixer('S590', ins=(SluC-0, MemDis-1, SP2-0),
                    outs='wastewater', init_with='Stream')
# effluent of WWmixer goes back to WWTP

# =============================================================================
# facilities
# =============================================================================

HXN = su.HTLHXN('HXN')

CHP = su.HTLCHP('CHP', ins=(GasMixer-0, 'natural_gas', 'air'),
              outs=('emission','solid_ash'), init_with='WasteStream', supplement_power_utility=False)

CHP.ins[1].price = 0.1685

RO_item = qs.StreamImpactItem(ID='RO_item',
                              linked_stream=MemDis.ins[3],
                              Acidification=0.53533,
                              Ecotoxicity=0.90848,
                              Eutrophication=0.0028322,
                              GlobalWarming=2.2663,
                              OzoneDepletion=0.00000025541,
                              PhotochemicalOxidation=0.0089068,
                              Carcinogenics=0.034791,
                              NonCarcinogenics=31.8,
                              RespiratoryEffects=0.0028778)
                              
H2SO4_item = qs.StreamImpactItem(ID='H2SO4_item',
                                 linked_stream=H2SO4_Tank.ins[0],
                                 Acidification=0.019678922,
                                 Ecotoxicity=0.069909345,
                                 Eutrophication=4.05E-06,
                                 GlobalWarming=0.008205666,
                                 OzoneDepletion=8.94E-10,
                                 PhotochemicalOxidation=5.04E-05,
                                 Carcinogenics=1.74E-03,
                                 NonCarcinogenics=1.68237815,
                                 RespiratoryEffects=9.41E-05)

MgCl2_item = qs.StreamImpactItem(ID='MgCl2_item',
                                 linked_stream=StruPre.ins[1],
                                 Acidification=0.77016,
                                 Ecotoxicity=0.97878,
                                 Eutrophication=0.00039767,
                                 GlobalWarming=2.8779,
                                 OzoneDepletion=4.94E-08,
                                 PhotochemicalOxidation=0.0072306,
                                 Carcinogenics=0.0050938,
                                 NonCarcinogenics=8.6916,
                                 RespiratoryEffects=0.004385)

H2_item = qs.StreamImpactItem(ID='H2_item',
                              linked_stream=RSP1.ins[0],
                              Acidification=0.81014,
                              Ecotoxicity=0.42747,
                              Eutrophication=0.0029415,
                              GlobalWarming=1.5624,
                              OzoneDepletion=1.80E-06,
                              PhotochemicalOxidation=0.0052545,
                              Carcinogenics=0.0026274,
                              NonCarcinogenics=8.5687,
                              RespiratoryEffects=0.0036698)

MgO_item = qs.StreamImpactItem(ID='MgO_item',
                               linked_stream=StruPre.ins[3],
                               Acidification=0.12584,
                               Ecotoxicity=2.7949,
                               Eutrophication=0.00063607,
                               GlobalWarming=1.1606,
                               OzoneDepletion=1.54E-08,
                               PhotochemicalOxidation=0.0017137,
                               Carcinogenics=0.018607,
                               NonCarcinogenics=461.54,
                               RespiratoryEffects=0.0008755)

NaOH_item = qs.StreamImpactItem(ID='NaOH_item',
                                linked_stream=MemDis.ins[2],
                                Acidification=0.33656,
                                Ecotoxicity=0.77272,
                                Eutrophication=0.00032908,
                                GlobalWarming=1.2514,
                                OzoneDepletion=7.89E-07,
                                PhotochemicalOxidation=0.0033971,
                                Carcinogenics=0.0070044,
                                NonCarcinogenics=13.228,
                                RespiratoryEffects=0.0024543)

NH4Cl_item = qs.StreamImpactItem(ID='NH4Cl_item',
                                 linked_stream=StruPre.ins[2],
                                 Acidification=0.34682,
                                 Ecotoxicity=0.90305, 
                                 Eutrophication=0.0047381,
                                 GlobalWarming=1.525,
                                 OzoneDepletion=9.22E-08,
                                 PhotochemicalOxidation=0.0030017,
                                 Carcinogenics=0.010029,
                                 NonCarcinogenics=14.85,
                                 RespiratoryEffects=0.0018387)

struvite_item = qs.StreamImpactItem(ID='struvite_item',
                                    linked_stream=StruPre.outs[0],
                                    Acidification=-0.122829597,
                                    Ecotoxicity=-0.269606396,
                                    Eutrophication=-0.000174952,
                                    GlobalWarming=-0.420850152,
                                    OzoneDepletion=-2.29549E-08,
                                    PhotochemicalOxidation=-0.001044087,
                                    Carcinogenics=-0.002983018,
                                    NonCarcinogenics=-4.496533528,
                                    RespiratoryEffects=-0.00061764)

NH42SO4_item = qs.StreamImpactItem(ID='NH42SO4_item',
                                   linked_stream=MemDis.outs[0],
                                   Acidification=-0.72917,
                                   Ecotoxicity=-3.4746,
                                   Eutrophication=-0.0024633,
                                   GlobalWarming=-1.2499,
                                   OzoneDepletion=-6.12E-08,
                                   PhotochemicalOxidation=-0.0044519,
                                   Carcinogenics=-0.036742,
                                   NonCarcinogenics=-62.932,
                                   RespiratoryEffects=-0.0031315)

natural_gas_item = qs.StreamImpactItem(ID='natural_gas_item',
                                       linked_stream=CHP.ins[1],
                                       Acidification=0.083822558,
                                       Ecotoxicity=0.063446198,
                                       Eutrophication=7.25E-05,
                                       GlobalWarming=1.584234288,
                                       OzoneDepletion=1.23383E-07,
                                       PhotochemicalOxidation=0.000973731,
                                       Carcinogenics=0.000666424,
                                       NonCarcinogenics=3.63204,
                                       RespiratoryEffects=0.000350917)

CHG_catalyst_item = qs.StreamImpactItem(ID='CHG_catalyst_item',
                                        linked_stream=PC3.outs[0],
                                        Acidification=991.6544196,
                                        Ecotoxicity=15371.08292,
                                        Eutrophication=0.45019348,
                                        GlobalWarming=484.7862509,
                                        OzoneDepletion=2.23437E-05,
                                        PhotochemicalOxidation=6.735405072,
                                        Carcinogenics=1.616793132,
                                        NonCarcinogenics=27306.37232,
                                        RespiratoryEffects=3.517184526)

HT_catalyst_item = qs.StreamImpactItem(ID='HT_catalyst_item',
                                       linked_stream=PC4.outs[0],
                                       Acidification=4.056401283,
                                       Ecotoxicity=50.26926274,
                                       Eutrophication=0.005759274,
                                       GlobalWarming=6.375878231,
                                       OzoneDepletion=1.39248E-06,
                                       PhotochemicalOxidation=0.029648759,
                                       Carcinogenics=0.287516945,
                                       NonCarcinogenics=369.791688,
                                       RespiratoryEffects=0.020809293)

HC_catalyst_item = qs.StreamImpactItem(ID='HC_catalyst_item',
                                       linked_stream=PC5.outs[0],
                                       Acidification=4.056401283,
                                       Ecotoxicity=50.26926274,
                                       Eutrophication=0.005759274,
                                       GlobalWarming=6.375878231,
                                       OzoneDepletion=1.39248E-06,
                                       PhotochemicalOxidation=0.029648759,
                                       Carcinogenics=0.287516945,
                                       NonCarcinogenics=369.791688,
                                       RespiratoryEffects=0.020809293)

dct = globals()
for unit_alias in (
        'WWTP', 'SluC', 'P1', 'H1', 'HTL', 'HTL_hx', 'HTL_drum', 'H2SO4_Tank', 'AcidEx',
        'M1', 'StruPre', 'CHG', 'CHG_pump', 'CHG_heating', 'CHG_cooling', 'V1', 'F1', 'MemDis', 'SP1',
        'P2', 'HT', 'HT_compressor', 'HT_hx_H2', 'HT_hx_oil', 'V2', 'H2', 'F2', 'V3', 'SP2', 'H3', 'D1', 'D2', 'D3', 'P3',
        'HC', 'HC_compressor', 'HC_hx_H2', 'HC_hx_oil', 'H4', 'V4', 'F3', 'D4', 'GasolineMixer', 'DieselMixer',
        'H5', 'H6', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'GasolineTank', 'DieselTank', 'FuelMixer',
        'GasMixer', 'WWmixer', 'RSP1', 'HXN', 'CHP'):
    dct[unit_alias].register_alias(unit_alias)
# so that qs.main_flowsheet.H1 works as well

sys = qs.System('sys', path=(WWTP, SluC, P1, H1, HTL, H2SO4_Tank, AcidEx,
                             M1, StruPre, CHG, V1, F1, MemDis, SP1,
                             P2, HT, V2, H2, F2, V3, SP2, H3, D1, D2, D3, P3,
                             HC, H4, V4, F3, D4, GasolineMixer, DieselMixer,
                             H5, H6, PC1, PC2, PC3, PC4, PC5, GasolineTank, DieselTank, FuelMixer,
                             GasMixer, WWmixer, RSP1), facilities=(HXN, CHP,))

sys.operating_hours = WWTP.operation_hours # 7920 hr Jones

lca_diesel = qs.LCA(system=sys, lifetime=30, lifetime_unit='yr',
                    Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30,
                    Cooling=lambda:sys.get_cooling_duty()/1000*30)

diesel_item = qs.StreamImpactItem(ID='diesel_item',
                                  linked_stream=FuelMixer.outs[0],
                                  Acidification=-0.25164,
                                  Ecotoxicity=-0.18748,
                                  Eutrophication=-0.0010547,
                                  GlobalWarming=-0.47694,
                                  OzoneDepletion=-6.42E-07,
                                  PhotochemicalOxidation=-0.0019456,
                                  Carcinogenics=-0.00069252,
                                  NonCarcinogenics=-2.9281,
                                  RespiratoryEffects=-0.0011096)

lca_sludge = qs.LCA(system=sys, lifetime=30, lifetime_unit='yr',
              Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30,
              Cooling=lambda:sys.get_cooling_duty()/1000*30)

tea = create_tea(sys)