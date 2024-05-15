#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Jianan Feng <jiananf2@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References:

(1) David, E.; Stanciu, V.; Sandru, C.; Armeanu, A.; Niculescu, V.
    Exhaust Gas Treatment Technologies for Pollutant Emission Abatement from
    Fossil Fuel Power Plants. In Sustainable Development and Planning III;
    WIT Press: Algarve, Portugal, 2007; Vol. II, pp 923–932.
    https://doi.org/10.2495/SDP070882.
(2) Huang, Z.; Grim, R. G.; Schaidle, J. A.; Tao, L. The Economic Outlook for
    Converting CO2 and Electrons to Molecules. Energy Environ. Sci.
    2021, 14 (7), 3664–3678. https://doi.org/10.1039/D0EE03525D.
(3) Hu, L.; Wrubel, J. A.; Baez-Cotto, C. M.; Intia, F.; Park, J. H.;
    Kropf, A. J.; Kariuki, N.; Huang, Z.; Farghaly, A.; Amichi, L.;
    Saha, P.; Tao, L.; Cullen, D. A.; Myers, D. J.; Ferrandon, M. S.;
    Neyerlin, K. C. A Scalable Membrane Electrode Assembly Architecture
    for Efficient Electrochemical Conversion of CO2 to Formic Acid.
    Nat Commun 2023, 14 (1), 7605. https://doi.org/10.1038/s41467-023-43409-6.

'''

import os, qsdsan as qs, biosteam as bst
from qsdsan import sanunits as qsu
from biosteam.units import PressureFilter, DrumDryer, AdsorptionColumnTSA
from qsdsan.utils import clear_lca_registries
from exposan.co2_sorbent import (
    _load_components,
    create_tea,
    )
from exposan.co2_sorbent import _sanunits as su
from biosteam import settings

__all__ = (
    'create_system_A',
    'create_system_B'
    )

# =============================================================================
# Al(OH)3 + HCOOH
# =============================================================================
def create_system_A(product='formic acid'):

    flowsheet_ID = 'ALF_A'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    # TODO: update CEPCI year
    bst.CE = qs.CEPCI_by_year[2020]
    
    # TODO: change the flow rate of AlOH3
    # TODO: add price for Al(OH)3
    aluminum_hydroxide = qs.WasteStream(ID='aluminum_hydroxide',
                                        phase='s',
                                        AlH3O3=100,
                                        units='kg/h',
                                        T=25+273.15)
    
    # TODO: change the flow rate of DI H2O
    # TODO: add price for H2O
    water = qs.WasteStream(ID='water',
                           H2O=100,
                           units='kg/h',
                           T=25+273.15)

    M1 = qsu.Mixer(ID='feed_mixer_1',
                   ins=(aluminum_hydroxide, water),
                   outs='AlOH3_H2O',
                   init_with='Stream',
                   rigorous=True,
                   conserve_phases=True)
    M1.register_alias('M1')
    
    # TODO: change the flow rate of formic acid (98% w/w, Xue et al. 2018)
    # TODO: add price for formic acid
    HCOOH_H2O = qs.WasteStream(ID='formic_acid',
                               HCOOH=350,
                               H2O=100,
                               units='kg/h',
                               T=25+273.15)
    
    M2 = qsu.Mixer(ID='feed_mixer_2',
                   ins=(M1-0, HCOOH_H2O),
                   outs='AlOH3_H2O_HCOOH',
                   init_with='Stream',
                   rigorous=True,
                   conserve_phases=True)
    M2.register_alias('M2')
    
    R1 = su.ALFProduction(ID='ALF_production',
                          ins=M2-0,
                          outs='ALF_solution')
    R1.register_alias('R1')
    
    # TODO: determine crystallization temperature
    # TODO: confirm no cooling utility is needed (it said atmospheric batch crystallization, see _batch_crystallizer.py)
    C1 = su.ALFCrystallizer(ID='ALF_crystallizer',
                            ins=R1-0,
                            outs='ALF_mixed',
                            T=273.15,
                            crystal_ALF_yield=1)
    C1.register_alias('C1')   
    
    # TODO: determine split ratio for ALF and HCOOH
    F1 = PressureFilter(ID='ALF_filter',
                        ins=C1-0,
                        outs=('retentate','permeate'),
                        moisture_content=0.35,
                        split={'C3H3AlO6':0.83,'HCOOH':0.05})
    F1.register_alias('F1')
    
    # TODO: determine moisture content for ALF
    D1 = DrumDryer(ID='ALF_dryer',
                   ins=(F1-0,'dryer_air','natural_gas'),
                   outs=('dryed_ALF','hot_air','emissions'),
                   moisture_content=0.0,
                   split={'HCOOH':1})
    D1.register_alias('D1')
    
    # coal-fired power plant flue gas composition: 13% CO2, 5% O2, and 82% N2 (David et al. 2007)
    # for a typical 1000 MWh plant, assume CO2=780000 kg/h (Huang et al. 2021)
    flue_gas = qs.WasteStream(ID='flue_gas',
                              CO2=780000,
                              H2O=300000,
                              N2=4920000,
                              phase='g',
                              units='kg/h',
                              T=25+273.15)
    
    TSA = su.ALFTemperatureSwingAdsorption(ID='ALF_TSA',
                                           ins=(flue_gas, D1-0, 'regenerated_ALF_in'),
                                           outs=('offgas','CO2','used_ALF','regenerated_ALF_out'))
    TSA.register_alias('TSA')
    
    # TODO: replace LowTemperatureElectrolysis with CO2ElectrolyzerSystem
    E1 = su.CO2ElectrolyzerSystem(ID='CO2_electrolyzer',
                                  ins=(TSA-1, 'process_water'),
                                  outs=('product','mixed_offgas'),
                                  target_product=product,
                                  current_density=0.2,
                                  cell_voltage=2.3,
                                  cathodic_overpotential=0.454,
                                  product_selectivity=0.9,
                                  converstion=0.5)
    E1.register_alias('E1')
    
    # TODO: determine the offgas fate for E1 (maybe sending it to a CHP unit?)
    
    # TODO: need to match up temperatures
    # HXN = qsu.HeatExchangerNetwork('HXN', T_min_app=5, force_ideal_thermo=True)
    # HXN.register_alias('HXN')
    
    # TODO: only add a cooling tower if there is a big finanical benefit
    # CT = bst.facilities.CoolingTower('CT')
    # CT.register_alias('CT')
    # from biorefineries.lactic import price
    # CT.ins[-1].price = price['Cooling tower chems']

    # opearting hours: 7884 h/year (Hu et al. 2023 SI)
    sys = qs.System.from_units('sys_ALF_A', units=list(flowsheet.unit), operating_hours=7884)
    sys.register_alias('sys')

    create_tea(sys)

    sys.simulate()
    sys.diagram()

# =============================================================================
# Electrochemical_cell
# 
# distillation
# 
# HXN
# 
# Cooling_tower
# 
#     sys = qs.System.from_units(
#         f'sys_ALF_A',
#         units=list(flowsheet.unit), 
#         operating_hours=,
#         )
# 
#     qs.StreamImpactItem(ID='',
#                         linked_stream=stream.,
#                         Acidification=,
#                         Ecotoxicity=,
#                         Eutrophication=,
#                         GlobalWarming=,
#                         OzoneDepletion=,
#                         PhotochemicalOxidation=,
#                         Carcinogenics=,
#                         NonCarcinogenics=,
#                         RespiratoryEffects=)
#     
#     create_tea(system=sys,
#                IRR_value=,
#                income_tax_value=,
#                finance_interest_value=,
#                labor_cost_value=)
#     
#     qs.LCA(system=sys,
#            lifetime=,
#            lifetime_unit='yr',
#            Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
#            Cooling=lambda:sys.get_cooling_duty()/1000*lifetime)
#     
#     sys.simulate()
# =============================================================================
    
    return sys

# =============================================================================
# Bauxite + HCOOH
# =============================================================================
def create_system_B():

    flowsheet_ID = 'ALF_B'