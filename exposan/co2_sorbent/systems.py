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
[5] Xue, M.; Gao, B.; Li, R.; Sun, J. Aluminum Formate (AF): Synthesis,
    Characterization and Application in Dye Wastewater Treatment.
    Journal of Environmental Sciences 2018, 74, 95–106.
[6] Davis, R. E.; Grundl, N. J.; Tao, L.; Biddy, M. J.; Tan, E. C.;
    Beckham, G. T.; Humbird, D.; Thompson, D. N.; Roni, M. S. Process Design
    and Economics for the Conversion of Lignocellulosic Biomass to Hydrocarbon
    Fuels and Coproducts: 2018 Biochemical Design Case Update;
    Biochemical Deconstruction and Conversion of Biomass to Fuels and Products
    via Integrated Biorefinery Pathways; NREL/TP--5100-71949, 1483234; 2018;
    p NREL/TP--5100-71949, 1483234. https://doi.org/10.2172/1483234.
'''

import os, qsdsan as qs, biosteam as bst
from qsdsan import sanunits as qsu
from biosteam.units import DrumDryer
from qsdsan.utils import clear_lca_registries
from exposan.co2_sorbent import (
    _load_components,
    create_tea,
    )
from exposan.co2_sorbent import _sanunits as su
from qsdsan.utils import auom
from math import ceil

# HXN (heat exchanger network), CWP (chilled water package)
# and BT (BoilerTurbogenerator) don't help in these systems

# TODO: from Jeremy: check if LCA studies (e.g., reports from DOE) for
# DAC (direct air capture) includes construction
# if not, we don't need to include LCA for constructions
# but we can add LCA for constructions if that's of interest later

# TODO: change all TEA and LCA data to the U.S.-based values
# if impossible, especially for LCA: use US, then RER, then RoW, then GLO
# also, remember to use LCA data collected for the 'cutoff' model and 'IPCC' method
# TODO: confirm IPCC 2013 (no LT) is ok

__all__ = (
    'create_system_A', # ALF production using Al(OH)3
    'create_system_B', # ALF production using Bauxite
    'create_system_C' # CO2 capture and utilization
    )

_lb_to_kg = auom('lb').conversion_factor('kg')
_MMBTU_to_MJ = auom('MMBTU').conversion_factor('MJ')

# GDPCTPI (Gross Domestic Product: Chain-type Price Index)
# https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20)
GDPCTPI = {2010: 89.619,
           2011: 91.466,
           2012: 93.176,
           2013: 94.786,
           2014: 96.436,
           2015: 97.277,
           2016: 98.208,
           2017: 100.000,
           2018: 102.290,
           2019: 104.008,
           2020: 105.407,
           2021: 110.220,
           2022: 117.995,
           2023: 122.284}

# 5 2016$/MMBTU (Davis et al. 2018)
# natural gas heating value 39 MJ/m3 (https://ecoquery.ecoinvent.org/3.8/apos/dataset/14395/documentation)
# natural gas density: 0.678 kg/m3 (https://en.wikipedia.org/wiki/Natural_gas, accessed 2024-06-02)
bst.stream_prices['Natural gas'] = 5/_MMBTU_to_MJ*39/0.678/GDPCTPI[2016]*GDPCTPI[2022]

# 0.0054 2010$/gal (Jouny et al. 2018)
# 1 gal water = 3.7854 kg water (https://www.inchcalculator.com/convert/gallon-to-kilogram/, accessed 2024-06-02)
bst.stream_prices['Reverse osmosis water'] = 0.0054/3.7854/GDPCTPI[2010]*GDPCTPI[2022]

# =============================================================================
# ALF production: Al(OH)3 + HCOOH
# =============================================================================
def create_system_A(AlH3O3=2416.7, # to produce 100 metric ton of ALF per day
                    electricity_price=0.0832,
                    yearly_operating_days=350,
                    lifetime=20):

    flowsheet_ID = 'ALF_production_A'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    # 2022 industrial electricity price is 0.0832 $/kWh
    # https://www.eia.gov/electricity/annual/html/epa_02_04.html (accessed 2024-05-20)
    # electricity price in the future and theoretical scenarios are based on Huang et al. 2021:
    # 0.030 $/kWh (future scenario), 0.020 $/kWh (theoretical scenario)
    bst.PowerUtility.price = electricity_price
    
    _load_components()
    bst.CE = qs.CEPCI_by_year[2022]
    
    folder = os.path.dirname(__file__)
    qs.ImpactIndicator.load_from_file(os.path.join(folder, 'data/impact_indicators.csv'))
    qs.ImpactItem.load_from_file(os.path.join(folder, 'data/impact_items.xlsx'))
    
    aluminum_hydroxide = qs.WasteStream(ID='aluminum_hydroxide',
                                        phase='s',
                                        AlH3O3=AlH3O3,
                                        units='kg/h',
                                        T=25+273.15)
    # TODO: check if there are other sources for the price since this source seems not that reliable
    # https://businessanalytiq.com/procurementanalytics/index/aluminum-hydroxide-price-index/
    # (accessed 2024-05-20)
    # 2022 average:
    aluminum_hydroxide.price = 0.424
    
    # add water to generate a sludge that contains 5.1% w/w of aluminum, Mikola et al. 2013
    water = qs.WasteStream(ID='water',
                           H2O=AlH3O3/78*27/0.051-AlH3O3,
                           units='kg/h',
                           T=25+273.15)
    water.price = bst.stream_prices['Reverse osmosis water']

    M1 = qsu.Mixer(ID='feed_mixer_1',
                   ins=(aluminum_hydroxide, water),
                   outs='AlH3O3_H2O',
                   init_with='Stream',
                   rigorous=True,
                   conserve_phases=True)
    M1.register_alias('M1')
    
    # add 80% w/w HCOOH to reach a 1:3.5 of Al/HCOOH ratio, Mikola et al. 2013
    formic_acid = qs.WasteStream(ID='formic_acid',
                                 HCOOH=AlH3O3/78*3.5*46,
                                 H2O=AlH3O3/78*3.5*46/0.8*0.2,
                                 units='kg/h',
                                 T=25+273.15)
    # TODO: check if there are other sources for the price since this source seems not that reliable
    # https://businessanalytiq.com/procurementanalytics/index/formic-acid-price-index/
    # (accessed 2024-05-20)
    # 2022 average:
    formic_acid.price = 0.82*0.8 + bst.stream_prices['Reverse osmosis water']*0.2
    
    M2 = qsu.Mixer(ID='feed_mixer_2',
                   ins=(M1-0, formic_acid),
                   outs='AlOH3_H2O_HCOOH',
                   init_with='Stream',
                   rigorous=True,
                   conserve_phases=True)
    M2.register_alias('M2')
    
    R1 = su.ALFProduction(ID='ALF_production',
                          ins=M2-0,
                          outs='ALF_solution',
                          T=333.15, # [K], Xue et al. 2018
                          P=101325, # [Pa]
                          tau=2) # [h], 2 h in Xue et al. 2018, 1 h in Mikola et al. 2013
    R1.register_alias('R1')
    
    C1 = su.ALFCrystallizer(ID='ALF_crystallizer',
                            ins=R1-0,
                            outs='ALF_mixed',
                            tau=5,
                            N=3,
                            T=298.15, # room temperature, Evans et al. 2022
                            crystal_ALF_yield=0.83) # crystal_ALF_yield = 0.83, Evans et al. 2022
    C1.register_alias('C1')
    
    # split ratio can be found in biosteam/units/solids_separation.py (soluble chemicals~0.036)
    F1 = su.ALFPressureFilter(ID='ALF_filter',
                              ins=C1-0,
                              outs=('retentate','permeate'),
                              moisture_content=0.35,
                              split={'C3H3AlO6_s':1,
                                     'C3H3AlO6_l':0,
                                     'HCOOH':0.036})
    F1.register_alias('F1')
    
    # assume the RO cake (produced after evaporator) can be potentially reused
    # therefore, no cost and environmental impact are given to it
    RO = su.ReverseOsmosis(ID='Reverse_osmosis',
                            ins=F1-1,
                            outs=('RO_water','brine'),
                            water_recovery=0.987)
    RO.outs[0].price = bst.stream_prices['Reverse osmosis water']
    RO.register_alias('RO')
    
    # note the price of natural price is already included (natural gas as a utility instead of a stream)
    # note hot_air contains a small amount of HCOOH, which is not considered as a greenhouse gas
    D1 = DrumDryer(ID='ALF_dryer',
                   ins=(F1-0,'dryer_air','natural_gas'),
                   outs=('dryed_ALF','hot_air','emissions'),
                   moisture_content=0,
                   split={'HCOOH':1})
    D1.register_alias('D1')
    
    SP1 = qsu.Splitter('SP1',
                       ins=D1-2,
                       outs=('natural_gas_LCA','H2O_discharge'),
                       split={'CO2':1,'H2O':0})
    SP1.register_alias('SP1')
    
    # TODO: do we need to activate ALF here? see Evans et al. 2022 and Ben's slides
    
    S1 = bst.StorageTank('ALF_storage_tank', ins=D1-0, outs='ALF',
                         tau=24*7, vessel_material='Stainless steel')
    S1.register_alias('S1')
    
    sys = qs.System.from_units(ID='sys_ALF_A',
                               units=list(flowsheet.unit),
                               operating_hours=yearly_operating_days*24)
    sys.register_alias('sys')
    
    # natural gas density: 0.678 kg/m3 (https://en.wikipedia.org/wiki/Natural_gas, accessed 06-02-2024)
    # natural gas MW: 19 (16.8-22.8 from https://en.wikipedia.org/wiki/Natural_gas, accessed 06-02-2024)
    # market group for natural gas, high pressure, GLO
    # does not include emissions
    # use produced CO2 amount to calculte the GlobalWarming for natural gas
    qs.StreamImpactItem(ID='natural_gas',
                        linked_stream=stream.natural_gas_LCA,
                        GlobalWarming=1+0.33658049/0.678*19/44)
    
    # market for aluminium hydroxide, GLO
    qs.StreamImpactItem(ID='aluminum_hydroxide',
                        linked_stream=stream.aluminum_hydroxide,
                        GlobalWarming=0.98811449)
    
    # water production, ultrapure, RER
    qs.StreamImpactItem(ID='water',
                        linked_stream=stream.water,
                        GlobalWarming=0.0028140858)
    
    # market for formic acid, RoW
    # for 80% HCOOH
    qs.StreamImpactItem(ID='formic_acid',
                        linked_stream=stream.formic_acid,
                        GlobalWarming=2.8683493*0.8+0.0028140858*0.2)
    
    # water production, ultrapure, RER
    qs.StreamImpactItem(ID='RO_water',
                        linked_stream=stream.RO_water,
                        GlobalWarming=-0.0028140858)
    
    create_tea(sys, lifetime=lifetime)
    
    qs.LCA(system=sys,
           lifetime=lifetime,
           lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-\
                               sys.get_electricity_production())*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Heating=lambda:sys.get_heating_duty()/1000*lifetime)
    
    # note simulating twice in this system may cause errors
    # TODO: does this affect creating models?
    
    sys.diagram()
    
    print(f"production: {sys.flowsheet.ALF.F_mass/1000*24} metric ton ALF/day")
    print(f"TEA: {sys.TEA.solve_price(sys.flowsheet.ALF)} $/kg ALF")
    print(f"LCA: {sys.LCA.get_total_impacts()['GlobalWarming']/sys.flowsheet.ALF.F_mass/sys.operating_hours/sys.LCA.lifetime} kg CO2 eq/kg ALF")
    
    return sys

# =============================================================================
# ALF production: Bauxite + HCOOH
# =============================================================================
def create_system_B(bauxite=2730.8, # to produce 100 metric ton of ALF per day
                    bauxite_price=0.075, # $/kg, related to bauxite_purity
                    bauxite_Al2O3=0.6, # Al2O3 weight ratio in bauxite, related to bauxite_price
                    bauxite_SiO2=0.11, # SiO2 weight ratio in bauxite, related to bauxite_price
                    electricity_price=0.0832,
                    yearly_operating_days=350,
                    lifetime=20):

    flowsheet_ID = 'ALF_production_B'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    # 2022 industrial electricity price is 0.0832 $/kWh
    # https://www.eia.gov/electricity/annual/html/epa_02_04.html (accessed 2024-05-20)
    # electricity price in the future and theoretical scenarios are based on Huang et al. 2021:
    # 0.030 $/kWh (future scenario), 0.020 $/kWh (theoretical scenario)
    bst.PowerUtility.price = electricity_price
    
    _load_components()
    bst.CE = qs.CEPCI_by_year[2022]
    
    folder = os.path.dirname(__file__)
    qs.ImpactIndicator.load_from_file(os.path.join(folder, 'data/impact_indicators.csv'))
    qs.ImpactItem.load_from_file(os.path.join(folder, 'data/impact_items.xlsx'))
    
    # TODO: confirm bauxite composition and make updates in _components.py
    bauxite_ore = qs.WasteStream(ID='bauxite_ore',
                                 phase='s',
                                 Al2O3=bauxite*bauxite_Al2O3,
                                 SiO2=bauxite*bauxite_SiO2,
                                 Fe2O3=bauxite*(1-bauxite_Al2O3-bauxite_SiO2),
                                 units='kg/h',
                                 T=25+273.15)
    # TODO: check if there are other sources for the price since this source seems not that reliable
    # https://en.institut-seltene-erden.de/aktuelle-preise-von-basismetallen/ (accessed 2024-05-23)
    # Al2O3 60% min
    bauxite_ore.price = bauxite_price
    
    G1 = su.BauxiteHammerMill(ID='bauxite_hammer_mill',
                              ins=bauxite_ore,
                              outs='crushed_bauxite')
    G1.register_alias('G1')
    
    # add water to generate a sludge that contains 5.1% w/w of aluminum, Mikola et al. 2013
    water = qs.WasteStream(ID='water',
                           H2O=bauxite*bauxite_Al2O3/102*54/0.051-bauxite,
                           units='kg/h',
                           T=25+273.15)
    water.price = bst.stream_prices['Reverse osmosis water']

    M1 = qsu.Mixer(ID='feed_mixer_1',
                    ins=(G1-0, water),
                    outs='bauxite_H2O',
                    init_with='Stream',
                    rigorous=True,
                    conserve_phases=True)
    M1.register_alias('M1')
    
    Al_mol = bauxite*bauxite_Al2O3/102*2
    Fe_mol = bauxite*(1-bauxite_Al2O3-bauxite_SiO2)/160*2
    
    # add 80% w/w HCOOH to reach a 1:3.5 of Al/HCOOH ratio, Mikola et al. 2013
    formic_acid = qs.WasteStream(ID='formic_acid',
                                 HCOOH=(Al_mol+Fe_mol)*3.5*46,
                                 H2O=(Al_mol+Fe_mol)*3.5*46/0.8*0.2,
                                 units='kg/h',
                                 T=25+273.15)
    # TODO: check if there are other sources for the price since this source seems not that reliable
    # https://businessanalytiq.com/procurementanalytics/index/aluminum-hydroxide-price-index/
    # (accessed 2024-05-20)
    # 2022 average:
    formic_acid.price = 0.82*0.8 + bst.stream_prices['Reverse osmosis water']*0.2
    
    M2 = qsu.Mixer(ID='feed_mixer_2',
                    ins=(M1-0, formic_acid),
                    outs='bauxite_H2O_HCOOH',
                    init_with='Stream',
                    rigorous=True,
                    conserve_phases=True)
    M2.register_alias('M2')    
    
    # TODO: does Fe2O3 react with HCOOH? Seems yes. Then how to separate ALF with FEF? And if yes, we need to further increase HCOOH amount
    # TODO: adjust reaction conditions if necessary
    R1 = su.ALFProduction(ID='ALF_production',
                          ins=M2-0,
                          outs='ALF_solution',
                          T=333.15, # [K], Xue et al. 2018
                          P=101325, # [Pa]
                          tau=2) # [h], 2 h in Xue et al. 2018, 1 h in Mikola et al. 2013
    R1.register_alias('R1')    
    
    # TODO: confirm the reaction between Fe2O3 and HCOOH product CO2
    # TODO: if not, no need to include PC1 and S2WS1 here
    PC1 = su.PhaseChanger(ID='PhaseChanger',
                          ins=R1-0,
                          outs=('slurry','carbon_dioxide'))
    PC1.register_alias('PC1')
    
    S2WS1 = su.S2WS(ID='S2WS1',
                    ins=PC1-1,
                    outs='carbon_dioxide_LCA')
    S2WS1.register_alias('S2WS1')
    
    # split ratio can be found in biosteam/units/solids_separation.py (soluble chemicals~0.036)
    F1 = su.SolidPressureFilter(ID='solid_filter',
                                ins=PC1-0,
                                outs=('solid_waste','ALF_permeate'),
                                moisture_content=0.35,
                                split={'SiO2':1,
                                       'Fe':1,
                                       'C3H3AlO6':0.036,
                                       'HCOOH':0.036})
    F1.register_alias('F1')
    
    S2WS2 = su.S2WS(ID='S2WS2',
                    ins=F1-0,
                    outs='solid_waste_LCA')
    S2WS2.register_alias('S2WS2')
    # TODO: update the price here
    # 5-9 $/ton (1 ton = 907.185 kg) from https://www.sciencedirect.com/science/article/pii/S0892687521003137 # TODO: add this as a reference once deciding to use this source
    S2WS2.outs[0].price=-0.0077
    
    C1 = su.ALFCrystallizer(ID='ALF_crystallizer',
                            ins=F1-1,
                            outs='ALF_mixed',
                            tau=5,
                            N=3,
                            T=298.15, # room temperature, Evans et al. 2022
                            crystal_ALF_yield=0.83) # crystal_ALF_yield = 0.83, Evans et al. 2022
    C1.register_alias('C1')
    
    # split ratio can be found in biosteam/units/solids_separation.py (soluble chemicals~0.036)
    F2 = su.ALFPressureFilter(ID='ALF_filter',
                              ins=C1-0,
                              outs=('retentate','permeate'),
                              moisture_content=0.35,
                              split={'C3H3AlO6_s':1,
                                     'C3H3AlO6_l':0,
                                     'HCOOH':0.036})
    F2.register_alias('F2')
    
    # assume the RO cake (produced after evaporator) can be potentially reused
    # therefore, no cost and environmental impact are given to it
    RO = su.ReverseOsmosis(ID='Reverse_osmosis',
                           ins=F2-1,
                           outs=('RO_water','brine'),
                           water_recovery=0.987)
    RO.outs[0].price = bst.stream_prices['Reverse osmosis water']
    RO.register_alias('RO')

    # note the price of natural price is already included (natural gas as a utility instead of a stream)
    # note hot_air contains a small amount of HCOOH, which is not considered as a greenhouse gas
    D1 = DrumDryer(ID='ALF_dryer',
                    ins=(F2-0,'dryer_air','natural_gas'),
                    outs=('dryed_ALF','hot_air','emissions'),
                    moisture_content=0,
                    split={'HCOOH':1})
    D1.register_alias('D1')
    
    SP1 = qsu.Splitter('SP1',
                       ins=D1-2,
                       outs=('natural_gas_LCA','H2O_discharge'),
                       split={'CO2':1,'H2O':0})
    SP1.register_alias('SP1')
    
    # TODO: do we need to activate ALF here? see Evans et al. 2022 and Ben's slides
    
    S1 = bst.StorageTank('ALF_storage_tank', ins=D1-0, outs='ALF',
                          tau=24*7, vessel_material='Stainless steel')
    S1.register_alias('S1')
    
    sys = qs.System.from_units(ID='sys_ALF_B',
                               units=list(flowsheet.unit),
                               operating_hours=yearly_operating_days*24)
    sys.register_alias('sys')
    
    # natural gas density: 0.678 kg/m3 (https://en.wikipedia.org/wiki/Natural_gas, accessed 06-02-2024)
    # natural gas MW: 19 (16.8-22.8 from https://en.wikipedia.org/wiki/Natural_gas, accessed 06-02-2024)
    # market group for natural gas, high pressure, GLO
    # does not include emissions
    # use produced CO2 amount to calculte the GlobalWarming for natural gas
    qs.StreamImpactItem(ID='natural_gas',
                        linked_stream=stream.natural_gas_LCA,
                        GlobalWarming=1+0.33658049/0.678*19/44)
    
    # market for bauxite, GLO
    qs.StreamImpactItem(ID='bauxite_ore',
                        linked_stream=stream.bauxite_ore,
                        GlobalWarming=0.026592271)
    
    # water production, ultrapure, RER
    qs.StreamImpactItem(ID='water',
                        linked_stream=stream.water,
                        GlobalWarming=0.0028140858)
    
    # market for formic acid, RoW
    # for 80% HCOOH
    qs.StreamImpactItem(ID='formic_acid',
                        linked_stream=stream.formic_acid,
                        GlobalWarming=2.8683493*0.8+0.0028140858*0.2)
    
    # market for redmud from bauxite digestion, GLO
    qs.StreamImpactItem(ID='solid_waste',
                        linked_stream=stream.solid_waste_LCA,
                        GlobalWarming=0.011066165)
    
    # water production, ultrapure, RER
    qs.StreamImpactItem(ID='RO_water',
                        linked_stream=stream.RO_water,
                        GlobalWarming=-0.0028140858)
    
    # CO2 direct emission
    qs.StreamImpactItem(ID='carbon_dioxide',
                        linked_stream=stream.carbon_dioxide_LCA,
                        GlobalWarming=1)
    
    create_tea(sys, lifetime=lifetime)
    
    qs.LCA(system=sys,
           lifetime=lifetime,
           lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-\
                               sys.get_electricity_production())*lifetime,
           Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
           Heating=lambda:sys.get_heating_duty()/1000*lifetime)
    
    # note simulating twice in this system may cause errors
    # TODO: does this affect creating models?
    
    sys.diagram()
    
    print(f"production: {sys.flowsheet.ALF.F_mass/1000*24} metric ton ALF/day")
    print(f"TEA: {sys.TEA.solve_price(sys.flowsheet.ALF)} $/kg ALF")
    print(f"LCA: {sys.LCA.get_total_impacts()['GlobalWarming']/sys.flowsheet.ALF.F_mass/sys.operating_hours/sys.LCA.lifetime} kg CO2 eq/kg ALF")
    
    return sys    

# =============================================================================
# CO2 capture and utilization
# =============================================================================
def create_system_C(product='formic acid',
                    ALF_system='A',
                    flue_gas_flow_rate=6000000, # Huang et al. 2021, for a typical 1000 MWh coal-fired power plant
                    flue_gas_CO2=0.13,
                    flue_gas_N2=0.82,
                    purity=0.975, # Evans et al. 2022
                    recovery=0.945, # Evans et al. 2022
                    electricity_price=0.0832,
                    yearly_operating_days=350,
                    lifetime=20,
                    upgrade=True):

    flowsheet_ID = 'CCU'
    
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
    bst.CE = qs.CEPCI_by_year[2022]
    
    folder = os.path.dirname(__file__)
    qs.ImpactIndicator.load_from_file(os.path.join(folder, 'data/impact_indicators.csv'))
    qs.ImpactItem.load_from_file(os.path.join(folder, 'data/impact_items.xlsx'))
    
    # coal-fired power plant flue gas composition: 13% CO2, 5% O2, and 82% N2 (David et al. 2007)
    # temperature of the flue gas: 160 C (David et al. 2007)
    # for a typical 1000 MWh plant, assume CO2=780000 kg/h (Huang et al. 2021)
    flue_gas = qs.WasteStream(ID='flue_gas',
                              CO2=flue_gas_flow_rate*flue_gas_CO2,
                              N2=flue_gas_flow_rate*flue_gas_N2,
                              O2=flue_gas_flow_rate*(1-flue_gas_CO2-flue_gas_N2),
                              phase='g',
                              units='kg/h',
                              T=160+273.15)
    
    # TODO: do we need to cool down flue gas first so ALF can absorb CO2?
    
    # TODO: update ALF price and CI if necessary
    # for system A, producing 100 metric ton ALF per day needs 2416.7 kg Al)OH)3 per hour
    if ALF_system == 'A':
        adsorbent_cost = 1.441*1441/35.3147 # 1.437 $/kg ALF to $/ft3 (1441 kg/m3, 35.3147 ft3/m3)
        adsorbent_CI = 4.907 # kg CO2 eq/kg ALF
        
    # for system B, producing 100 metric ton ALF per day needs 2730.8 kg bauxite per hour
    if ALF_system == 'B':
        adsorbent_cost = 1.943*1441/35.3147 # 1.935 $/kg ALF to $/ft3 (1441 kg/m3, 35.3147 ft3/m3)
        adsorbent_CI = 7.703 # kg CO2 eq/kg ALF
    
    TSA = su.ALFTSA(ID='ALF_TSA',
                    ins=(flue_gas,'air'),
                    outs=('captured_carbon_dioxide','TSA_offgas'),
                    adsorbent='ALF',
                    adsorbent_cost=adsorbent_cost,
                    adsorbent_lifetime=10,
                    adsorbate_ID='CO2',
                    split=dict(O2=(flue_gas_CO2*recovery/purity-flue_gas_CO2*recovery)/(1-flue_gas_CO2),
                               N2=(flue_gas_CO2*recovery/purity-flue_gas_CO2*recovery)/(1-flue_gas_CO2),
                               CO2=recovery),
                    superficial_velocity=2160,
                    regeneration_velocity=1332,
                    cycle_time=5,
                    rho_adsorbent=1441,
                    adsorbent_capacity=0.1188,
                    T_regeneration=418,
                    vessel_material='Stainless steel 316',
                    vessel_type='Vertical',
                    length_unused=1.219,    
                    treatment_capacity=10000,
                    regen_CO2=purity,
                    regen_N2=(1-purity)*flue_gas_N2/(1-flue_gas_CO2),
                    regen_O2=(1-purity)*(1-flue_gas_CO2-flue_gas_N2)/(1-flue_gas_CO2))
    TSA.register_alias('TSA')
    
    if not upgrade:
        sys = qs.System.from_units('sys_ALF_A', units=list(flowsheet.unit), operating_hours=yearly_operating_days*24)
        sys.register_alias('sys')
    
        create_tea(sys, lifetime=lifetime)
        
        # simulate here before LCA to enable the calculation of outlet streams
        # simulate twice (here and later in qs.LCA) in this system does not affect the results
        sys.simulate()
    
        # note LCA for ALF and CO2 are calculated as 'Other', not 'Construction' or 'Stream'
        # use ceil to calculate the needed ALF to be conservative (assume no salvage value)
        qs.LCA(system=sys,
                lifetime=lifetime,
                lifetime_unit='yr',
                Electricity=lambda:(sys.get_electricity_consumption()-\
                                    sys.get_electricity_production())*lifetime,
                Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
                Heating=lambda:sys.get_heating_duty()/1000*lifetime,
                CO2=(-stream.flue_gas.imass['CO2']+stream.TSA_offgas.imass['CO2'])*24*yearly_operating_days*lifetime,
                ALF=TSA.design_results['Number of sets']*\
                    TSA.design_results['Number of reactors']*\
                    TSA.vessel_volume*1441*\
                    ceil(lifetime/TSA.equipment_lifetime['ALF'])*\
                    adsorbent_CI) # 1441 kg/m3
        
        print(f"TEA: {sys.TEA.solve_price(sys.flowsheet.captured_carbon_dioxide)} $/kg CO2 captured")
        print(f"LCA: {sys.LCA.get_total_impacts()['GlobalWarming']/sys.flowsheet.flue_gas.imass['CO2']/sys.operating_hours/sys.LCA.lifetime} kg CO2 eq reduction/kg CO2 captured")
    
    else:
        # TODO: do we need to cool down CO2 first?
        
        E1 = su.CO2ElectrolyzerSystem(ID='CO2_electrolyzer',
                                      ins=(TSA-0,'process_water'),
                                      outs=('product','hydrogen','electro_offgas'),
                                      target_product=product,
                                      current_density=0.2,
                                      cell_voltage=2.3,
                                      cathodic_overpotential=0.454,
                                      product_selectivity=0.9,
                                      converstion=0.5,
                                      operating_days_per_year=yearly_operating_days)
        E1.register_alias('E1')
        E1.ins[1].price = bst.stream_prices['Reverse osmosis water']
        # use the price of gray hydrogen since the LCA data for hydrogen only has gray hydrogen
        # and it makes more sense to replace gray hydrogen first
        # 0.7306 2016$/lb (Davis et al. 2018)
        E1.outs[1].price = 0.7306/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
        
        # TODO: if gas products: do we want to add gas compressors and storage vessels for gas products and hydrogen?
        # TODO: if liquid products: do we want to add storage tanks for liquid products and gas compressors and storage vessels for hydrogen?
        # Jouny et al. 2018 ignored these
        
        sys = qs.System.from_units('sys_ALF_A', units=list(flowsheet.unit), operating_hours=yearly_operating_days*24)
        sys.register_alias('sys')
        
        # water production, ultrapure, RER
        qs.StreamImpactItem(ID='water',
                            linked_stream=stream.process_water,
                            GlobalWarming=0.0028140858)
        
        # market for hydrogen, gaseous, GLO
        qs.StreamImpactItem(ID='hydrogen',
                            linked_stream=stream.hydrogen,
                            GlobalWarming=-1.5621537)
        
        create_tea(sys, lifetime=lifetime)
        
        # simulate here before LCA to enable the calculation of outlet streams
        # simulate twice (here and later in qs.LCA) in this system does not affect the results
        sys.simulate()
        
        # note LCA for ALF and CO2 are calculated as 'Other', not 'Construction' or 'Stream'
        # use ceil to calculate the needed ALF to be conservative (assume no salvage value)
        qs.LCA(system=sys,
                lifetime=lifetime,
                lifetime_unit='yr',
                Electricity=lambda:(sys.get_electricity_consumption()-\
                                    sys.get_electricity_production())*lifetime,
                Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
                Heating=lambda:sys.get_heating_duty()/1000*lifetime,
                CO2=(-stream.flue_gas.imass['CO2']+stream.TSA_offgas.imass['CO2'])*24*yearly_operating_days*lifetime,
                ALF=TSA.design_results['Number of sets']*\
                    TSA.design_results['Number of reactors']*\
                    TSA.vessel_volume*1441*\
                    ceil(lifetime/TSA.equipment_lifetime['ALF'])*\
                    adsorbent_CI) # 1441 kg/m3
        
        print(f"TEA: {sys.TEA.solve_price(sys.flowsheet.product)} $/kg {product}")
        print(f"LCA: {sys.LCA.get_total_impacts()['GlobalWarming']/sys.flowsheet.product.F_mass/sys.operating_hours/sys.LCA.lifetime} kg CO2 eq/kg {product}")
    
    sys.diagram()
    
    return sys