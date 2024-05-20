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

[1] Huang, Z.; Grim, R. G.; Schaidle, J. A.; Tao, L. The Economic Outlook for
    Converting CO2 and Electrons to Molecules. Energy Environ. Sci.
    2021, 14 (7), 3664–3678. https://doi.org/10.1039/D0EE03525D.
[2] Evans, H. A.; Mullangi, D.; Deng, Z.; Wang, Y.; Peh, S. B.; Wei, F.;
    Wang, J.; Brown, C. M.; Zhao, D.; Canepa, P.; Cheetham, A. K.
    Aluminum Formate, Al(HCOO)3: An Earth-Abundant, Scalable, and Highly
    Selective Material for CO2 Capture. Science Advances 2022, 8 (44),
    eade1473. https://doi.org/10.1126/sciadv.ade1473.
[3] Mikola, M.; Rämö, J.; Sarpola, A.; Tanskanen, J. Removal of Different NOM
    Fractions from Surface Water with Aluminium Formate. Separation and
    Purification Technology 2013, 118, 842–846.
    https://doi.org/10.1016/j.seppur.2013.08.037.
[4] David, E.; Stanciu, V.; Sandru, C.; Armeanu, A.; Niculescu, V.
    Exhaust Gas Treatment Technologies for Pollutant Emission Abatement from
    Fossil Fuel Power Plants. In Sustainable Development and Planning III;
    WIT Press: Algarve, Portugal, 2007; Vol. II, pp 923–932.
    https://doi.org/10.2495/SDP070882.
[5] Hu, L.; Wrubel, J. A.; Baez-Cotto, C. M.; Intia, F.; Park, J. H.;
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

# GDPCTPI (Gross Domestic Product: Chain-type Price Index): https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20)
GDPCTPI = {2016: 98.208,
           2017: 100.000,
           2018: 102.290,
           2019: 104.008,
           2020: 105.407,
           2021: 110.220,
           2022: 117.995,
           2023: 122.284}

# =============================================================================
# Al(OH)3 + HCOOH
# =============================================================================
# TODO: determine the relationship between ALF production and CO2 capture (one systems or two systems)
# TODO: if one system, set waste_ratio in TSA? (the current waste_ratio = 0.7 seems not correct)
# TODO: if two systems, calculate ALF price, and set it as an argument in TSA?
def create_system_A(product='formic acid',
                    flue_gas=6000000, # Huang et al. 2021, for a 1000 MWh coal-fired power plant
                    purity=0.975, # Evans et al. 2022
                    recovery=0.945, # Evans et al. 2022
                    electricity_price=0.0832,
                    yearly_operating_days=350):

    flowsheet_ID = 'ALF_A'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    # 2022 industrial electricity price is 0.0832 $/kWh, https://www.eia.gov/electricity/annual/html/epa_02_04.html (accessed 2024-05-20)
    # electricity price in the future and theoretical scenarios are based on Huang et al. 2021:
    # 0.030 $/kWh (future scenario), 0.020 $/kWh (theoretical scenario)
    bst.PowerUtility.price = electricity_price
    
    _load_components()
    # all costs in this analysis are 2022$
    bst.CE = qs.CEPCI_by_year[2022]
    # TODO: update all price/cost to the baseline year (2022)
    
    # TODO: determine the flow rate of AlOH3 based on the flue gas amount
    aluminum_hydroxide = qs.WasteStream(ID='aluminum_hydroxide',
                                        phase='s',
                                        AlH3O3=100,
                                        units='kg/h',
                                        T=25+273.15)
    aluminum_hydroxide.price = 0.424 # 2022 average: https://businessanalytiq.com/procurementanalytics/index/aluminum-hydroxide-price-index/ (accessed 2024-05-20)
    
    # TODO: determine the flow rate of water needed to generate a sludge that contains 5.1% w/w of aluminum, Mikola et al. 2013
    water = qs.WasteStream(ID='water',
                           H2O=[mass of AlH3O3]/78*27/0.051-[mass of AlH3O3],
                           units='kg/h',
                           T=25+273.15)
    water.price = 0.00044/GDPCTPI[2016]*GDPCTPI[2022] # Davis et al. 2018 NREL ($0.0002/lb, 2016$)

    M1 = qsu.Mixer(ID='feed_mixer_1',
                   ins=(aluminum_hydroxide, water),
                   outs='AlH3O3_H2O',
                   init_with='Stream',
                   rigorous=True,
                   conserve_phases=True)
    M1.register_alias('M1')
    
    # TODO: determine the flow rate of 80% w/w HCOOH needed to reach a 1:3.5 of Al/HCOOH ratio, Mikola et al. 2013
    HCOOH_H2O = qs.WasteStream(ID='formic_acid',
                               HCOOH=[mass of AlH3O3]/78*3.5*46,
                               H2O=[mass of AlH3O3]/78*3.5*46/0.8*0.2,
                               units='kg/h',
                               T=25+273.15)
    HCOOH_H2O.price = 0.82*0.8 + 0.00044/GDPCTPI[2016]*GDPCTPI[2022]*0.2 # formic acid, 2022 average: https://businessanalytiq.com/procurementanalytics/index/aluminum-hydroxide-price-index/ (accessed 2024-05-20)
    
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
    # TODO: decide if it takes cost to treat 'permeate'
    F1 = PressureFilter(ID='ALF_filter',
                        ins=C1-0,
                        outs=('retentate','permeate'),
                        moisture_content=0.35,
                        split={'C3H3AlO6':0.83,'HCOOH':0.05})
    F1.register_alias('F1')
    
    # TODO: determine moisture content for ALF
    # TODO: note the price of natural gas is included, but is it 2022$?
    D1 = DrumDryer(ID='ALF_dryer',
                   ins=(F1-0,'dryer_air','natural_gas'),
                   outs=('dryed_ALF','hot_air','emissions'),
                   moisture_content=0.0,
                   split={'HCOOH':1})
    D1.register_alias('D1')
    
    # coal-fired power plant flue gas composition: 13% CO2, 5% O2, and 82% N2 (David et al. 2007)
    # for a typical 1000 MWh plant, assume CO2=780000 kg/h (Huang et al. 2021)
    flue_gas = qs.WasteStream(ID='flue_gas',
                              CO2=flue_gas*0.13,
                              O2=flue_gas*0.05,
                              N2=flue_gas*0.82,
                              phase='g',
                              units='kg/h',
                              T=25+273.15)
    
    TSA = su.ALFTemperatureSwingAdsorption(ID='ALF_TSA',
                                           split=dict(O2=(0.13*recovery/purity-0.13*recovery)/0.87,
                                                      N2=(0.13*recovery/purity-0.13*recovery)/0.87,
                                                      CO2=recovery),
                                           ins=(flue_gas, D1-0, 'regenerated_ALF_in'),
                                           outs=('offgas','CO2','used_ALF','regenerated_ALF_out'))
    TSA.register_alias('TSA')
    
    E1 = su.CO2ElectrolyzerSystem(ID='CO2_electrolyzer',
                                  ins=(TSA-1, 'process_water'),
                                  outs=('product','mixed_offgas'),
                                  target_product=product,
                                  current_density=0.2,
                                  cell_voltage=2.3,
                                  cathodic_overpotential=0.454,
                                  product_selectivity=0.9,
                                  converstion=0.5,
                                  PSA_operating_cost=0.25,
                                  operating_days_per_year=yearly_operating_days)
    E1.register_alias('E1')
    E1.ins[1].price = 0.00044/GDPCTPI[2016]*GDPCTPI[2022] # Davis et al. 2018 NREL ($0.0002/lb, 2016$)
    
    # TODO: determine the offgas fate for E1 (maybe sending it to a CHP unit?)
    # TODO: or separating and selling H2? See Huang et al. 2021 SI, Table S1 for price
    
    # TODO: need to match up temperatures
    # HXN = qsu.HeatExchangerNetwork('HXN', T_min_app=5, force_ideal_thermo=True)
    # HXN.register_alias('HXN')
    
    # TODO: only add a cooling tower if there is a big finanical benefit
    # CT = bst.facilities.CoolingTower('CT')
    # CT.register_alias('CT')
    # from biorefineries.lactic import price
    # CT.ins[-1].price = price['Cooling tower chems']

    # opearting hours: 7884 h/year (Hu et al. 2023 SI)
    sys = qs.System.from_units('sys_ALF_A', units=list(flowsheet.unit), operating_hours=yearly_operating_days*24)
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