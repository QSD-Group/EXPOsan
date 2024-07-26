#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Jianan Feng <jiananf2@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import biosteam as bst
from qsdsan.utils import auom
from exposan.co2_sorbent import create_system_A, create_system_B, create_system_C

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
#%% analysis

# =============================================================================
# # system A cost and CI breakdown
# =============================================================================
sys = create_system_A(AlH3O3=2416.7, electricity_price=0.0832, clean_electricity=False)

# construction cost [$]
sys.installed_cost
sys.flowsheet.R1.installed_cost
sys.flowsheet.C1.installed_cost
sys.flowsheet.F1.installed_cost
sys.flowsheet.RO.installed_cost
sys.flowsheet.D1.installed_cost
sys.flowsheet.S1.installed_cost

# material cost [$/year]
sys.flowsheet.aluminum_hydroxide.F_mass*sys.flowsheet.aluminum_hydroxide.price*24*350
sys.flowsheet.water.F_mass*sys.flowsheet.water.price*24*350
sys.flowsheet.formic_acid.F_mass*sys.flowsheet.formic_acid.price*24*350
sys.flowsheet.RO_water.F_mass*(-sys.flowsheet.RO_water.price)*24*350

# utility [$/year]
 # duty > 0 is heating utility, duty < 0 is cooling utility
sys.heat_utilities[0].duty
sys.heat_utilities[0].cost*24*350
sys.heat_utilities[1].duty
sys.heat_utilities[1].cost*24*350
sys.power_utility.cost*24*350
sys.flowsheet.natural_gas.F_mass*bst.stream_prices['Natural gas']*24*350

# material CI [kg CO2 eq/year]
table = sys.LCA.get_impact_table('Stream')
print(table['GlobalWarming [kg CO2-eq]'])

# utility CI [kg CO2 eq/year]
table = sys.LCA.get_impact_table('Other')
print(table['GlobalWarming [kg CO2-eq]'])

# =============================================================================
# # system B cost and CI breakdown
# =============================================================================
sys = create_system_B(bauxite=2730.8, electricity_price=0.0832, clean_electricity=False)

# construction cost [$]
sys.installed_cost
sys.flowsheet.G1.installed_cost
sys.flowsheet.R1.installed_cost
sys.flowsheet.F1.installed_cost
sys.flowsheet.C1.installed_cost
sys.flowsheet.F2.installed_cost
sys.flowsheet.RO.installed_cost
sys.flowsheet.D1.installed_cost
sys.flowsheet.S1.installed_cost

# material cost [$/year]
sys.flowsheet.bauxite_ore.F_mass*sys.flowsheet.bauxite_ore.price*24*350
sys.flowsheet.water.F_mass*sys.flowsheet.water.price*24*350
sys.flowsheet.formic_acid.F_mass*sys.flowsheet.formic_acid.price*24*350
sys.flowsheet.solid_waste_LCA.F_mass*(-sys.flowsheet.solid_waste_LCA.price)*24*350
sys.flowsheet.RO_water.F_mass*(-sys.flowsheet.RO_water.price)*24*350

# utility [$/year]
 # duty > 0 is heating utility, duty < 0 is cooling utility
sys.heat_utilities[0].duty
sys.heat_utilities[0].cost*24*350
sys.heat_utilities[1].duty
sys.heat_utilities[1].cost*24*350
sys.power_utility.cost*24*350
sys.flowsheet.natural_gas.F_mass*bst.stream_prices['Natural gas']*24*350

# material CI [kg CO2 eq/year]
table = sys.LCA.get_impact_table('Stream')
print(table['GlobalWarming [kg CO2-eq]'])

# utility CI [kg CO2 eq/year]
table = sys.LCA.get_impact_table('Other')
print(table['GlobalWarming [kg CO2-eq]'])

# =============================================================================
# # system C (ALF from systems A and B) CO2 capture
# =============================================================================
sys = create_system_C(product='formic acid', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.0832, clean_electricity=False, upgrade=False)

sys = create_system_C(product='formic acid', ALF_system='B', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.0832, clean_electricity=False, upgrade=False)

# =============================================================================
# # system C (ALF from system A) cost and CI breakdown
# =============================================================================
sys = create_system_C(product='formic acid', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.0832, clean_electricity=False, upgrade=True)

# construction cost [$]
sys.installed_cost
sys.flowsheet.ALF_TSA.installed_cost # incluse ALF
sys.flowsheet.CO2_electrolyzer.installed_cost

# material cost [$/year]
sys.flowsheet.process_water.F_mass*sys.flowsheet.process_water.price*24*350
sys.flowsheet.hydrogen.F_mass*(-sys.flowsheet.hydrogen.price)*24*350

# utility [$/year]
 # duty > 0 is heating utility, duty < 0 is cooling utility
sys.heat_utilities[0].duty
sys.heat_utilities[0].cost*24*350
sys.power_utility.cost*24*350

# material CI [kg CO2 eq/year]
table = sys.LCA.get_impact_table('Stream')
print(table['GlobalWarming [kg CO2-eq]'])

# utility CI [kg CO2 eq/year]
table = sys.LCA.get_impact_table('Other')
print(table['GlobalWarming [kg CO2-eq]'])

# =============================================================================
# # system C (ALF from system B) cost and CI breakdown
# =============================================================================
sys = create_system_C(product='formic acid', ALF_system='B', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.0832, clean_electricity=False, upgrade=True)

# construction cost [$]
sys.installed_cost
sys.flowsheet.ALF_TSA.installed_cost # incluse ALF
sys.flowsheet.CO2_electrolyzer.installed_cost

# material cost [$/year]
sys.flowsheet.process_water.F_mass*sys.flowsheet.process_water.price*24*350
sys.flowsheet.hydrogen.F_mass*(-sys.flowsheet.hydrogen.price)*24*350

# utility [$/year]
 # duty > 0 is heating utility, duty < 0 is cooling utility
sys.heat_utilities[0].duty
sys.heat_utilities[0].cost*24*350
sys.power_utility.cost*24*350

# material CI [kg CO2 eq/year]
table = sys.LCA.get_impact_table('Stream')
print(table['GlobalWarming [kg CO2-eq]'])

# utility CI [kg CO2 eq/year]
table = sys.LCA.get_impact_table('Other')
print(table['GlobalWarming [kg CO2-eq]'])

# =============================================================================
# # system C (ALF from system A) CO2 valorization
# =============================================================================
sys = create_system_C(product='carbon monoxide', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.0832, clean_electricity=False, upgrade=True)
sys = create_system_C(product='methane', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.0832, clean_electricity=False, upgrade=True)
sys = create_system_C(product='methanol', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.0832, clean_electricity=False, upgrade=True)
sys = create_system_C(product='formic acid', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.0832, clean_electricity=False, upgrade=True)
sys = create_system_C(product='ethylene', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.0832, clean_electricity=False, upgrade=True)
sys = create_system_C(product='ethanol', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.0832, clean_electricity=False, upgrade=True)
sys = create_system_C(product='propanol', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.0832, clean_electricity=False, upgrade=True)

sys = create_system_C(product='carbon monoxide', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.03, clean_electricity=True, upgrade=True)
sys = create_system_C(product='methane', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.03, clean_electricity=True, upgrade=True)
sys = create_system_C(product='methanol', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.03, clean_electricity=True, upgrade=True)
sys = create_system_C(product='formic acid', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.03, clean_electricity=True, upgrade=True)
sys = create_system_C(product='ethylene', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.03, clean_electricity=True, upgrade=True)
sys = create_system_C(product='ethanol', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.03, clean_electricity=True, upgrade=True)
sys = create_system_C(product='propanol', ALF_system='A', flue_gas_flow_rate=6000000, adsorbent_lifetime=10, electricity_price=0.03, clean_electricity=True, upgrade=True)