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

# TODO: add LCA for all systems
# TODO: from Jeremy: check if LCA studies (e.g., reports from DOE) for DAC (direct air capture) includes constructin
# TODO: if not, we don't need to include LCA for constructions
# TODO: but we can add LCA for constructions if that's of interest later
# TODO: update all price/cost to the baseline year (2022)
__all__ = (
    'create_system_A', # ALF production using Al(OH)3
    'create_system_B', # ALF production using Bauxite
    # 'create_system_C' # CO2 capture and utilization
    )

# GDPCTPI (Gross Domestic Product: Chain-type Price Index)
# https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20)
GDPCTPI = {2016: 98.208,
           2017: 100.000,
           2018: 102.290,
           2019: 104.008,
           2020: 105.407,
           2021: 110.220,
           2022: 117.995,
           2023: 122.284}

# =============================================================================
# ALF production: Al(OH)3 + HCOOH
# =============================================================================
def create_system_A(AlH3O3=1000,
                    electricity_price=0.0832,
                    yearly_operating_days=350):

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
    # all costs in this analysis are 2022$
    bst.CE = qs.CEPCI_by_year[2022]
    
    folder = os.path.dirname(__file__)
    qs.ImpactIndicator.load_from_file(os.path.join(folder, 'data/impact_indicators.csv'))
    qs.ImpactItem.load_from_file(os.path.join(folder, 'data/impact_items.xlsx'))
    
    aluminum_hydroxide = qs.WasteStream(ID='aluminum_hydroxide',
                                        phase='s',
                                        AlH3O3=AlH3O3,
                                        units='kg/h',
                                        T=25+273.15)
    # 2022 average:
    # https://businessanalytiq.com/procurementanalytics/index/aluminum-hydroxide-price-index/
    # (accessed 2024-05-20)
    # TODO: check if there are other sources for the price since this source seems not that reliable
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
    # formic acid, 2022 average:
    # https://businessanalytiq.com/procurementanalytics/index/aluminum-hydroxide-price-index/
    # (accessed 2024-05-20)
    # TODO: check if there are other sources for the price since this source seems not that reliable
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
                            tau=5, # TODO: test in the uncertainty analysis if this is important, if important, then decide a more appropriate value
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
    
    # TODO: check if the RO cake (produced after evaporator) can be potentially reused
    # asssume no cost and environmental impact associated with it
    RO1 = su.ReverseOsmosis(ID='Reverse_osmosis',
                            ins=F1-1,
                            outs=('RO_water','brine'),
                            water_recovery=0.987)
    RO1.outs[0].price = bst.stream_prices['Reverse osmosis water']
    RO1.register_alias('RO1')

    # natural price is already included (bst.stream_prices['Natural gas'] = 0.218 $/kg)
    D1 = DrumDryer(ID='ALF_dryer',
                   ins=(F1-0,'dryer_air','natural_gas'),
                   outs=('dryed_ALF','hot_air','emissions'),
                   moisture_content=0,
                   split={'HCOOH':1})
    D1.register_alias('D1')
    
    S2WS = su.S2WS(ID='S2WS',
                   ins=D1-2,
                   outs='natural_gas_LCA')
    S2WS.register_alias('S2WS')
    
    S1 = bst.StorageTank('ALF_storage_tank', ins=D1-0, outs='ALF',
                         tau=24*7, vessel_material='Stainless steel')
    S1.register_alias('S1')
    
    # HXN (heat exchanger network), CWP (chilled water package)
    # and BT (BoilerTurbogenerator) don't help in this system
    
    sys = qs.System.from_units(ID='sys_ALF_A',
                               units=list(flowsheet.unit),
                               operating_hours=yearly_operating_days*24)
    sys.register_alias('sys')
    
    # for LCA, using ecoinvent database 3.8 apos (default) and TRACI method for now
    # TODO: can change to ecoinvent database 3.8 cutoff and IPCC method if necessary
    # TODO: confirm HCOOH in the stream.hot_air does not contribute to GWP
    
    # market group for natural gas, high pressure, GLO
    # does not include emissions
    # 39 MJ/m3 (https://ecoquery.ecoinvent.org/3.8/apos/dataset/14395/documentation)
    # 53.6 MJ/kg (https://www.sciencedirect.com/science/article/pii/B9780128095973003357)
    qs.StreamImpactItem(ID='natural_gas',
                        linked_stream=stream.natural_gas_LCA,
                        GlobalWarming=1+0.33690797/39*53.6/44*16)
    
    # market for aluminium hydroxide, GLO
    qs.StreamImpactItem(ID='aluminum_hydroxide',
                        linked_stream=stream.aluminum_hydroxide,
                        GlobalWarming=0.98853282)
    
    # water production, deionised, RoW
    qs.StreamImpactItem(ID='water',
                        linked_stream=stream.water,
                        GlobalWarming=0.00030228066)
    
    # market for formic acid, RoW
    # for 80% HCOOH
    qs.StreamImpactItem(ID='formic_acid',
                        linked_stream=stream.formic_acid,
                        GlobalWarming=2.8612976*0.8+0.00030228066*0.2)
    
    # water production, deionised, RoW
    qs.StreamImpactItem(ID='RO_water',
                        linked_stream=stream.RO_water,
                        GlobalWarming=-0.00030228066)
    
    create_tea(sys)
    
    # from Jeremy: LCA for construction may not be needed
    # TODO: postpone doing LCA for construction
    # TODO: make sure the lifetimes for TEA and LCA are the same and decide if 20 years is appropriate
    qs.LCA(system=sys,
           lifetime=20,
           lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-\
                               sys.get_electricity_production())*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Heating=lambda:sys.get_heating_duty()/1000*20)
    
    # !!! errors like 'UndefinedPhase: <DrumDryer: ALF_dryer> 'G'' is because the system being simulated for multiple times (note qs.LCA simulate the system as well)
    # TODO: decide whether the above glitch affects creating models
    
    sys.diagram()
    
    print(f"production: {sys.flowsheet.ALF.F_mass/1000*24:.2f} metric ton ALF/day")
    print(f"TEA: {sys.TEA.solve_price(sys.flowsheet.ALF):.2f} $/kg ALF")
    print(f"LCA: {sys.LCA.get_total_impacts()['GlobalWarming']/sys.flowsheet.ALF.F_mass/sys.operating_hours/sys.LCA.lifetime:.2f} kg CO2 eq/kg ALF")
    
    return sys

# =============================================================================
# ALF production: Bauxite + HCOOH
# =============================================================================
def create_system_B(bauxite=1000,
                    bauxite_price=0.075, # $/kg, related to bauxite_purity
                    bauxite_Al2O3=0.6, # Al2O3 weight ratio in bauxite, related to bauxite_price
                    bauxite_SiO2=0.11, # SiO2 weight ratio in bauxite, related to bauxite_price
                    electricity_price=0.0832,
                    yearly_operating_days=350):

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
    # all costs in this analysis are 2022$
    bst.CE = qs.CEPCI_by_year[2022]
    
    folder = os.path.dirname(__file__)
    qs.ImpactIndicator.load_from_file(os.path.join(folder, 'data/impact_indicators.csv'))
    qs.ImpactItem.load_from_file(os.path.join(folder, 'data/impact_items.xlsx'))
    
    bauxite_ore = qs.WasteStream(ID='bauxite_ore',
                                 phase='s',
                                 Al2O3=bauxite*bauxite_Al2O3,
                                 SiO2=bauxite*bauxite_SiO2,
                                 Fe2O3=bauxite*(1-bauxite_Al2O3-bauxite_SiO2),
                                 units='kg/h',
                                 T=25+273.15)
    # https://en.institut-seltene-erden.de/aktuelle-preise-von-basismetallen/ (accessed 2024-05-23)
    # Al2O3 60% min
    # TODO: check if there are other sources for the price since this source seems not that reliable
    bauxite_ore.price = bauxite_price
    
    G1 = su.BauxiteHammerMill(ID='bauxite_hammer_mill',
                              ins=bauxite_ore,
                              outs='crushed_bauxite',
                              Al2O3_ratio=0.6)
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
    # formic acid, 2022 average:
    # https://businessanalytiq.com/procurementanalytics/index/aluminum-hydroxide-price-index/
    # (accessed 2024-05-20)
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
    
    # TODO: adding CO2 in the 'gas' to LCA?
    PC1 = su.PhaseChanger(ID='PhaseChanger',
                          ins=R1-0,
                          outs=('slurry','carbon_dioxide'))
    PC1.register_alias('PC1')
    
    S2WS1 = su.S2WS(ID='S2WS1',
                    ins=PC1-1,
                    outs='carbon_dioxide_LCA')
    S2WS1.register_alias('S2WS1')
    
    # TODO: how to dispose the 'retentate' here?
    # split ratio can be found in biosteam/units/solids_separation.py (soluble chemicals~0.036)
    F1 = su.SolidPressureFilter(ID='solid_filter',
                                ins=PC1-0,
                                outs=('solid_waste','ALF_permeate'),
                                moisture_content=0.35,
                                split={'SiO2':1,
                                       'Fe':1,
                                       'C3H3AlO6':0.036, # TODO: part of C3H3AlO6 may become solid during filtering. Increase the C3H3AlO6 split ratio here?
                                       'HCOOH':0.036})
    F1.register_alias('F1')
    
    # TODO: update the price here
    # 5-9 $/ton from https://www.sciencedirect.com/science/article/pii/S0892687521003137 # TODO: add this as a reference once deciding to use this source
    F1.outs[0].price=-0.0077
    
    S2WS2 = su.S2WS(ID='S2WS2',
                    ins=F1-0,
                    outs='solid_waste_LCA')
    S2WS2.register_alias('S2WS2')
    
    C1 = su.ALFCrystallizer(ID='ALF_crystallizer',
                            ins=F1-1,
                            outs='ALF_mixed',
                            tau=5, # TODO: test in the uncertainty analysis if this is important, if important, then decide a more appropriate value
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
    
    # TODO: check if the RO cake (produced after evaporator) can be potentially reused
    # asssume no cost and environmental impact associated with it
    RO1 = su.ReverseOsmosis(ID='Reverse_osmosis',
                            ins=F2-1,
                            outs=('RO_water','brine'),
                            water_recovery=0.987)
    RO1.outs[0].price = bst.stream_prices['Reverse osmosis water']
    RO1.register_alias('RO1')

    # natural price is already included (bst.stream_prices['Natural gas'] = 0.218 $/kg)
    D1 = DrumDryer(ID='ALF_dryer',
                    ins=(F2-0,'dryer_air','natural_gas'),
                    outs=('dryed_ALF','hot_air','emissions'),
                    moisture_content=0,
                    split={'HCOOH':1})
    D1.register_alias('D1')
    
    S2WS3 = su.S2WS(ID='S2WS3',
                    ins=D1-2,
                    outs='natural_gas_LCA')
    S2WS3.register_alias('S2WS3')
    
    S1 = bst.StorageTank('ALF_storage_tank', ins=D1-0, outs='ALF',
                          tau=24*7, vessel_material='Stainless steel')
    S1.register_alias('S1')
    
    # # HXN (heat exchanger network), CWP (chilled water package)
    # # and BT (BoilerTurbogenerator) don't help in this system
    
    sys = qs.System.from_units(ID='sys_ALF_B',
                               units=list(flowsheet.unit),
                               operating_hours=yearly_operating_days*24)
    sys.register_alias('sys')
    
    # for LCA, using ecoinvent database 3.8 apos (default) and TRACI method for now
    # TODO: can change to ecoinvent database 3.8 cutoff and IPCC method if necessary
    # TODO: confirm HCOOH in the stream.hot_air does not contribute to GWP
    
    # market group for natural gas, high pressure, GLO
    # does not include emissions
    # 39 MJ/m3 (https://ecoquery.ecoinvent.org/3.8/apos/dataset/14395/documentation)
    # 53.6 MJ/kg (https://www.sciencedirect.com/science/article/pii/B9780128095973003357)
    qs.StreamImpactItem(ID='natural_gas',
                        linked_stream=stream.natural_gas_LCA,
                        GlobalWarming=1+0.33690797/39*53.6/44*16)
    
    # market for bauxite, GLO
    qs.StreamImpactItem(ID='bauxite_ore',
                        linked_stream=stream.bauxite_ore,
                        GlobalWarming=0.026571501)
    
    # water production, deionised, RoW
    qs.StreamImpactItem(ID='water',
                        linked_stream=stream.water,
                        GlobalWarming=0.00030228066)
    
    # market for formic acid, RoW
    # for 80% HCOOH
    qs.StreamImpactItem(ID='formic_acid',
                        linked_stream=stream.formic_acid,
                        GlobalWarming=2.8612976*0.8+0.00030228066*0.2)
    
    # market for redmud from bauxite digestion
    qs.StreamImpactItem(ID='solid_waste',
                        linked_stream=stream.solid_waste_LCA,
                        GlobalWarming=0.011072384)
    
    # water production, deionised, RoW
    qs.StreamImpactItem(ID='RO_water',
                        linked_stream=stream.RO_water,
                        GlobalWarming=-0.00030228066)
    
    # TODO: confirm this
    # CO2 direct emission
    qs.StreamImpactItem(ID='carbon_dioxide',
                        linked_stream=stream.carbon_dioxide_LCA,
                        GlobalWarming=1)
    
    create_tea(sys)
    
    # from Jeremy: LCA for construction may not be needed
    # TODO: postpone doing LCA for construction
    # TODO: make sure the lifetimes for TEA and LCA are the same and decide if 20 years is appropriate
    qs.LCA(system=sys,
           lifetime=20,
           lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-\
                               sys.get_electricity_production())*20,
           Cooling=lambda:sys.get_cooling_duty()/1000*20,
           Heating=lambda:sys.get_heating_duty()/1000*20)
    
    # !!! errors like 'UndefinedPhase: <DrumDryer: ALF_dryer> 'G'' is because the system being simulated for multiple times (note qs.LCA simulate the system as well)
    # TODO: decide whether the above glitch affects creating models
    
    sys.diagram()
    
    print(f"production: {sys.flowsheet.ALF.F_mass/1000*24:.2f} metric ton ALF/day")
    print(f"TEA: {sys.TEA.solve_price(sys.flowsheet.ALF):.2f} $/kg ALF")
    print(f"LCA: {sys.LCA.get_total_impacts()['GlobalWarming']/sys.flowsheet.ALF.F_mass/sys.operating_hours/sys.LCA.lifetime:.2f} kg CO2 eq/kg ALF")
    
    return sys    

# =============================================================================
# CO2 capture and utilization
# =============================================================================
def create_system_C(product='formic acid',
                    ALF_system='A',
                    flue_gas=6000000, # Huang et al. 2021, for a 1000 MWh coal-fired power plant
                    purity=0.975, # Evans et al. 2022
                    recovery=0.945, # Evans et al. 2022
                    electricity_price=0.0832,
                    yearly_operating_days=350):

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
    # all costs in this analysis are 2022$
    bst.CE = qs.CEPCI_by_year[2022]
    # TODO: update all price/cost to the baseline year (2022)
    
    # coal-fired power plant flue gas composition: 13% CO2, 5% O2, and 82% N2 (David et al. 2007)
    # for a typical 1000 MWh plant, assume CO2=780000 kg/h (Huang et al. 2021)
    flue_gas = qs.WasteStream(ID='flue_gas',
                              CO2=flue_gas*0.13,
                              O2=flue_gas*0.05,
                              N2=flue_gas*0.82,
                              phase='g',
                              units='kg/h',
                              T=25+273.15)
    
    # TODO: update ALF price and CI if necessary
    if ALF_system == 'A':
        adsorbent_cost = 1.40 # $/kg ALF
        adsorbent_CI = 4.98 # kg CO2 eq/kg ALF
    if ALF_system == 'B':
        adsorbent_cost = 1.95 # $/kg ALF
        adsorbent_CI = 7.77 # kg CO2 eq/kg ALF
    
    TSA = su.ALFTemperatureSwingAdsorption(ID='ALF_TSA',
                                            split=dict(O2=(0.13*recovery/purity-0.13*recovery)/0.87,
                                                        N2=(0.13*recovery/purity-0.13*recovery)/0.87,
                                                        CO2=recovery),
                                            ins=(flue_gas, D1-0, 'regenerated_ALF_in'),
                                            outs=('offgas','CO2','used_ALF','regenerated_ALF_out'),
                                            adsorbent_cost=adsorbent_cost)
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
    E1.ins[1].price = bst.stream_prices['Reverse osmosis water']
    
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
    
    sys = qs.System.from_units('sys_ALF_A', units=list(flowsheet.unit), operating_hours=yearly_operating_days*24)
    sys.register_alias('sys')
    
    
    create_tea(sys)
    
    # TODO: add LCA for ALF in others
    qs.LCA(system=sys,
           lifetime=20,
           lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-\
                               sys.get_electricity_production())*20,
               Cooling=lambda:sys.get_cooling_duty()/1000*20,
               Heating=lambda:sys.get_heating_duty()/1000*20,
               ALF=)
    
    
    



    sys.simulate()
    sys.diagram()

# =============================================================================
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