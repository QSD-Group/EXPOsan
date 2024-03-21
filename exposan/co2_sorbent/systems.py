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

(1) Hu, L.; Wrubel, J. A.; Baez-Cotto, C. M.; Intia, F.; Park, J. H.;
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
def create_system_A():

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
    aluminum_hydroxide = qs.WasteStream('aluminum_hydroxide', phase='s', AlH3O3=100, units='kg/h', T=25+273.15)
    
    # TODO: change the flow rate of DI H2O
    # TODO: add price for H2O
    water = qs.WasteStream('water', H2O=100, units='kg/h', T=25+273.15)

    M1 = qsu.Mixer('feed_mixer_1', ins=(aluminum_hydroxide, water), outs='AlOH3_H2O', init_with='Stream', rigorous=True, conserve_phases=True)
    M1.register_alias('M1')
    
    # TODO: change the flow rate of formic acid
    # TODO: add price for formic acid
    HCOOH_H2O = qs.WasteStream('formic_acid', HCOOH=350, H2O=100, units='kg/h', T=25+273.15)
    
    M2 = qsu.Mixer('feed_mixer_2', ins=(M1-0, HCOOH_H2O), outs='AlOH3_H2O_HCOOH', init_with='Stream', rigorous=True, conserve_phases=True)
    M2.register_alias('M2')
    
    R1 = su.ALFProduction('ALF_production', M2-0, outs='ALF_solution')
    R1.register_alias('R1')
    
    # TODO: determine crystallization temperature
    # TODO: confirm no cooling utility is needed (it said atmospheric batch crystallization, see _batch_crystallizer.py)
    C1 = su.ALFCrystallizer('ALF_crystallizer', R1-0, outs='ALF_mixed', T=273.15, crystal_ALF_yield=1)
    C1.register_alias('C1')   
    
    # TODO: determine split ratio for ALF and HCOOH
    F1 = PressureFilter('ALF_filter', ins=C1-0, outs=('retentate','permeate'), moisture_content=0.35, split={'C3H3AlO6':0.83, 'HCOOH':0.05})
    F1.register_alias('F1')
    
    # TODO: determine moisture content for ALF
    D1 = DrumDryer('ALF_dryer', (F1-0,'dryer_air','natural_gas'), ('dryed_ALF','hot_air','emissions'), moisture_content=0.0, split={'HCOOH':1})
    D1.register_alias('D1')
    
    flue_gas = qs.WasteStream('flue_gas', CO2=780000, H2O=300000, N2=4920000, phase='g', units='kg/h', T=25+273.15)
    
    TSA = su.ALFTSA('ALF_TSA', ins=(flue_gas, D1-0, 'regenerated_ALF_in'), outs=('offgas','CO2','used_ALF','regenerated_ALF_in'))
    TSA.register_alias('TSA')
    
    # TODO: need to match up temperatures
    # HXN = qsu.HeatExchangerNetwork('HXN', T_min_app=5, force_ideal_thermo=True)
    # HXN.register_alias('HXN')
    
    # TODO: only add a cooling tower if there is a big finanical benefit
    # CT = bst.facilities.CoolingTower('CT')
    # CT.register_alias('CT')
    # from biorefineries.lactic import price
    # CT.ins[-1].price = price['Cooling tower chems']

    # opearting hours: Hu et al. 2023 SI
    sys = qs.System.from_units('sys_ALF_A', units=list(flowsheet.unit), operating_hours=7884)
    sys.register_alias('sys')

    create_tea(sys)

    sys.simulate()
    sys.diagram()

# =============================================================================
# Electrochemical_cell
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