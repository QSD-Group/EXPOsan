#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

This is a temporary script to translate BEAM in to QSDsan.

BEAM*2024 model:
North East Biosolids and Residuals Association (NEBRA), Northern Tilth LLC, and 
orthwest Biosolids. Estimating Greenhouse Gas Emissions from Biosolids Management.
BEAM*2024 Spreadsheet Model and Supporting Information, 2024.
https://www.BiosolidsGHGs.org (accessed 2025-06-16).

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# TODO: some parameters were rounded in the BEAM*2024, update after unlocking the spreadsheet

from qsdsan.utils import auom

_MGD_to_MLD = auom('gallon').conversion_factor('liter')
_BTU_to_kWh = auom('BTU').conversion_factor('kWh')
_GJ_to_BTU = auom('GJ').conversion_factor('BTU')
_MJ_to_BTU = auom('MJ').conversion_factor('BTU')
_C_to_K = 273.15
_N_to_N2O = 44/28
_C_to_CH4 = 16/12
_C_to_CO2 = 44/12

# raw wastewater residual solids density (before dewatering), tonne/m3
raw_solids_density = 1

# TODO: double check this, why this is less than 1
# dewatered wastewater residual solids density (after dewatering), tonne/m3
dewatered_solids_density = 1.1

# sawdust density, tonne/m3
bulking_agent_density = 0.25

# CH4 density, kg/m3
CH4_density = 0.707

# CH4 heat content, BTU/m3
CH4_heat = 35830

# natural gas heat content, BTU/m3
NG_heat = 36263

# TODO: CF / CI values based on BEAM, but may need to update this to use ecoinvent data
# TODO: this only include NG combustion (scope 1), consider adding the production of NG (scope 2)
# g/m3
NG_CI = 1901
# g/kWh
elec_CI = 226
# kg CO2 eq/kg
CH4_CF = 25
N2O_CF = 298
polymer_CI = 1.6
lime_CI = {'waste': 0,
           'virgin': 0.9}
# TODO: consider discounting N and P credit
N_fertilizer_CI = 3
P_fertilizer_CI = 2
# g CO2 eq/L
diesel_CI = 2697

#%% storage

# m3/day
raw_solids_volume = 100

# MGD
influent_flow_rate = 100

# mg/L
influent_BOD = 350

# -
BOD_removal = 0.9

# kW/1000 m3 sludge
elec_lagoon = 5.6

# see BEAM user guide, -
percentage_days_per_year_above_15_C = 0.5

# m
lagoon_depth = 3.5

# kg CH4/kg BOD (shallow: <= 2 m)
CH4_EF_shallow_lagoon = 0.12

# kg CH4/kg BOD (deep: > 2 m)
CH4_EF_deep_lagoon = 0.48

# kg/day
BOD_to_storage = influent_flow_rate*1000000*_MGD_to_MLD*influent_BOD/1000000*BOD_removal

aerated = True

if aerated:
    # kWh/day
    elec = raw_solids_volume*elec_lagoon*24/1000
else:
    # kWh/day
    elec = 0

# tonne CO2 eq/day
elec_emission = elec*elec_CI/1000000

if aerated:
    # tonne CO2 eq/day
    CH4_emission = 0
elif lagoon_depth <= 2:
    # tonne CO2 eq/day
    CH4_emission = BOD_to_storage*CH4_EF_shallow_lagoon*percentage_days_per_year_above_15_C/1000*CH4_CF
    
# output, tonne CO2 eq/year
# scope 1
CH4_emission*365
# scope 2
elec_emission*365

#%% conditioning thickening

# m3/day
raw_solids_volume = 10000

# -
solids_ratio_before = 0.01

# -
solids_ratio_after = 0.03

# polymer, kg/dry tonne
polymer_usage = 5

# dry tonne/day
solids_dry = raw_solids_volume*raw_solids_density*solids_ratio_before

# kg/day
polymer = solids_dry*polymer_usage

# tonne CO2 eq/day
polymer_emission = polymer*polymer_CI/1000

# electricity, kWh/dry tonne
elec_thickener = {'centrifuge': 33,
                  'other': 4.9}

# kWh/day
elec = solids_dry*elec_thickener['centrifuge']

# tonne CO2 eq/day
elec_emission = elec*elec_CI/1000000

# output, tonne CO2 eq/year
# scope 2
elec_emission*365
# scope 3
polymer_emission*365

#%% aerobic digestion

# thickened solids, m3/day
thickened_solids_volume = raw_solids_volume*solids_ratio_before/solids_ratio_after

# day
AeD_SRT = 40

# VS before digestion, -
VS_before_digestion_ratio = 0.7

# dry kg/day
VS_before_digestion = thickened_solids_volume*raw_solids_density*1000*solids_ratio_after*VS_before_digestion_ratio

# -
solids_volume_reduction = 0.2

# output solids, m3/day
aerobic_digested_solids_volume = thickened_solids_volume*(1 - solids_volume_reduction)

# -
VS_reduction = 0.45

# kg/day
VS_destroyed = VS_before_digestion*VS_reduction

# only not 0 if AeD is heated (not common), m3 natural gas/day
NG = 0

# kW/m3 of digestor volume
AeD_unit_elec = 0.03

# kWh/day
AeD_elec_requirement = thickened_solids_volume*AeD_SRT*AeD_unit_elec*24

# tonne CO2 eq/day
NG_combustion_emission = NG*NG_CI/1000000

# tonne CO2 eq/day
elec_emission = elec*elec_CI/1000000

# output, tonne CO2 eq/year
# scope 1
NG_combustion_emission*365
# scope 2
elec_emission*365

#%% anaerobic digestion

# thickened solids, m3/day
thickened_solids_volume = raw_solids_volume*solids_ratio_before/solids_ratio_after

# dry tonne/day
dry_solids = thickened_solids_volume*solids_ratio_after

# VS before digestion, -
VS_before_digestion_ratio = 0.7

# dry tonne/day
VS_before_digestion = dry_solids*VS_before_digestion_ratio

# nitrogen before digestion, -
N_before_digestion_ratio = 0.05

# dry tonne/day
N_before_digestion = dry_solids*N_before_digestion_ratio

# phosphorus before digestion, -
P_before_digestion_ratio = 0.01

# dry tonne/day
P_before_digestion = dry_solids*P_before_digestion_ratio

# day
AD_SRT = 22

# -
solids_volume_reduction = 0.2

# output solids, m3/day
anaerobic_digested_solids_volume = thickened_solids_volume*(1 - solids_volume_reduction)

# -
VS_reduction = 0.45

# kg/day
VS_destroyed = VS_before_digestion*1000*VS_reduction

# output solids, dry tonne/day
anaerobic_digested_solids_mass = dry_solids - VS_destroyed/1000

# VS after digestion, -
VS_after_digestion_ratio = (VS_before_digestion - VS_destroyed/1000)/anaerobic_digested_solids_mass

# nitrogen after digestion, -
N_after_digestion_ratio = N_before_digestion/anaerobic_digested_solids_mass

# phosphorus after digestion, -
P_after_digestion_ratio = P_before_digestion/anaerobic_digested_solids_mass

# biogas yield from destroyed VS, m3/kg VS
biogas_yield = 0.9

# m3/day
biogas_volume = VS_destroyed*biogas_yield

# CH4 in biogas, -
CH4_biogas = 0.65

# m3/day
CH4_volume = biogas_volume*CH4_biogas

# combined heat and power
CHP = True

# renewable natural gas
RNG = True

# - 
biogas_CHP_ratio = 0.55

# - 
biogas_RNG_ratio = 0.35

# - 
biogas_flare_ratio = 0.09

# -
biogas_fugitive_ratio = 0.01

assert biogas_CHP_ratio + biogas_RNG_ratio + biogas_flare_ratio + biogas_fugitive_ratio == 1

# -
CHP_heat_ratio = 0.42

# -
CHP_elec_ratio = 0.35

# -
CHP_heat = biogas_CHP_ratio*CHP_heat_ratio

# - 
CHP_elec = biogas_CHP_ratio*CHP_elec_ratio

# inefficient data is used by standard GHG accounting protocols including IPCC and EPA, - 
CHP_fugitive_CH4_ratio = {'normal': 0.003,
                          'inefficient': 0.01}

# tonne CO2 eq/day
CHP_fugitive_CH4 = CH4_volume*biogas_CHP_ratio*CHP_fugitive_CH4_ratio['inefficient']*CH4_density/1000*CH4_CF

# default data is used by standard GHG accounting protocols including IPCC and EPA, - 
flare_fugitive_CH4_ratio = {'default': 0.01,
                            'enclosed': 0,
                            'candlestick': 0.05}

# tonne CO2 eq/day
flare_fugitive_CH4 = CH4_volume*biogas_flare_ratio*flare_fugitive_CH4_ratio['default']*CH4_density/1000*CH4_CF

# CH4 heat content, BTU/m3
CH4_heat = 35830

# natural gas heat content, BTU/m3
NG_heat = 36263

# -
RNG_parasitic = 0.169

# NG equivalent generated, m3 natural gas/day
NG_generated = CH4_volume*biogas_CHP_ratio*CHP_heat_ratio*CH4_heat/NG_heat + CH4_volume*biogas_RNG_ratio*CH4_heat/NG_heat*(1 - RNG_parasitic)

# CH4 to electricity factor, -
net_capacity_factor = 0.85

# elec equivalent generated, kWh/day
elec_generated = CH4_volume*biogas_CHP_ratio*CHP_elec_ratio*CH4_heat*_BTU_to_kWh*net_capacity_factor

# m3 natural gas/m3 solids
AD_unit_heat_requiement = 3.3

# m3 natural gas/day
AD_heat_requirement = thickened_solids_volume*AD_unit_heat_requiement

# kW/m3 of digestor volume
AD_unit_elec_requiement = 0.0065

# kWh/day
AD_elec_requirement = thickened_solids_volume*AD_SRT*AD_unit_elec_requiement*24

# m3 natural gas/day
net_NG = AD_heat_requirement - NG_generated

# kWh/day
net_elec = AD_elec_requirement - elec_generated

# tonne CO2 eq/day
NG_combustion_emission = net_NG*NG_CI/1000000

# tonne CO2 eq/day
elec_emission = net_elec*elec_CI/1000000

# tonne CO2 eq/day
fugitive_emission = CH4_volume*biogas_fugitive_ratio*CH4_density/1000*CH4_CF + CHP_fugitive_CH4 + flare_fugitive_CH4

# output, tonne CO2 eq/year
# scope 1
(NG_combustion_emission + fugitive_emission)*365
# scope 2
elec_emission*365

#%% dewatering

# -
volume_reduction_before_dewatering = 0.2

# m3/day
solids_volume_before_dewatering = thickened_solids_volume*(1 - volume_reduction_before_dewatering)

# -
solids_content_after_dewatering = 0.2

# dry tonne/day
solids_dry = thickened_solids_volume*solids_ratio_after*(1 - VS_before_digestion_ratio*VS_reduction)

# polymer, kg/dry tonne
polymer_usage = 5

# kg/day
polymer = solids_dry*polymer_usage

# tonne CO2 eq/day
polymer_emission = polymer*polymer_CI/1000

# electricity, kWh/dry tonne
elec_dewatering = {'centrifuge': 107,
                   'other': 11.3}

# kWh/day
elec = solids_dry*elec_dewatering['centrifuge']

# tonne CO2 eq/day
elec_emission = elec*elec_CI/1000000

AD_before = True

# CH4 in anaerobically digested biosolids, mg CH4/L
CH4_digested = 6

# CH4 lost during dewatering, -
CH4_lost = 0.9

if AD_before:
    # tonne CO2 eq/day
    CH4_emission = solids_volume_before_dewatering*1000*CH4_digested*CH4_lost/1000000000*CH4_CF

# output, tonne CO2 eq/year
# scope 1
CH4_emission*365
# scope 2
elec_emission*365
# scope 3
polymer_emission*365

#%% thermal drying

# tonne/day
dewatered_solids_mass = 100

# -
solids_content_before_drying = 0.2

# dry tonne/day
dry_solids = dewatered_solids_mass*solids_content_before_drying

# -
solids_content_after_drying = 0.8

# tonne/day
water_removal = dewatered_solids_mass*(1 - solids_content_before_drying) - dry_solids/solids_content_after_drying*(1 - solids_content_after_drying)

# GJ/tonne
energy_unit_water_removal = 4.5

# BTU/day
energy_water_removal = water_removal*energy_unit_water_removal*_GJ_to_BTU

# m3 natural gas/day
NG = energy_water_removal/NG_heat

# tonne CO2 eq/day
NG_combustion_emission = NG*NG_CI/1000000

# kWh/dry tonne
elec_drying = 214

# kWh/day
elec = dry_solids*elec_drying

# tonne CO2 eq/day
elec_emission = elec*elec_CI/1000000

# output, tonne CO2 eq/year
# scope 1
NG_combustion_emission*365
# scope 2
elec_emission*365

#%% alkaline stabilization

# tonne/day
dewatered_solids_mass = 100

# -
solids_content_before_lime = 0.2

# dry tonne/day
dry_solids = dewatered_solids_mass*solids_content_before_lime

# tonne/dry tonne
lime_dose = {'class_A': 0.3,
             'class_B': 0.2}

# tonne/day
lime = dry_solids*lime_dose['class_A']

# tonne CO2 eq/day
lime_emission = lime*lime_CI['virgin']

# TODO: confirm this is just for combustion
# kg CO2 eq/dry tonne
NG_unit_emission = {'class_A': 15.6,
                    'class_B': 0}

# TODO: confirm this is just for combustion
# tonne CO2 eq/day
NG_emission = dry_solids*NG_unit_emission['class_A']/1000

# TODO: check, B>A?
# kWh/wet tonne
elec_unit = {'class_A': 3.7,
             'class_B': 4.9}

# kWh/day
elec = dewatered_solids_mass*elec['class_A']

# tonne CO2 eq/day
elec_emission = elec*elec_CI/1000000

# output, tonne CO2 eq/year
# scope 1
NG_emission*365
# scope 2
elec_emission*365
# scope 3
lime_emission*365

#%% composting

# TODO: the current code assume all bulking agents go to compost without reuse, which may not be true, can add a parameter to consider reuse

# tonne/day
dewatered_solids_mass = 100

# -
solids_content_before_composting = 0.2

# dry tonne/day
dry_solids = dewatered_solids_mass*solids_content_before_composting

# m3/day
solids_volume = dewatered_solids_mass/dewatered_solids_density

# organic carbon to VS ratio, -
OC_VS = 0.56

# organic carbon before digestion, -
OC_before_digestion_ratio = VS_before_digestion_ratio*OC_VS

# organic carbon after digestion, -
OC_after_digestion_ratio = VS_after_digestion_ratio*OC_VS

# volumetric ratio of bulking agent (sawdust) to solids
bulking_agent_solids = 3

# tonne/day
bulking_agent = solids_volume*bulking_agent_solids*bulking_agent_density

# -
bulking_agent_solids_ratio = 0.61

# based on dry weight, -
bulking_agent_OC_ratio = 0.518

# based on dry weight, -
bulking_agent_nitrogen_ratio = 0.0024

# blended feedstock C:N ratio, -
blended_C_N = (dry_solids*OC_after_digestion_ratio + bulking_agent*bulking_agent_solids_ratio*bulking_agent_OC_ratio)/\
              (dry_solids*N_after_digestion_ratio + bulking_agent*bulking_agent_solids_ratio*bulking_agent_nitrogen_ratio)

# blended feedstock solids content, -
blended_solids_ratio = (dry_solids + bulking_agent*bulking_agent_solids_ratio)/(dewatered_solids_mass + bulking_agent)

# -
bulking_agent_grinding_onsite = True

# L diesel/wet tonne
grind_unit_diesel = 3.3

if bulking_agent_grinding_onsite:
    # L/day
    grinding_diesel = bulking_agent*grind_unit_diesel

# for setting up and breaking down piles (i.e., other), L diesel/wet tonne
other_unit_diesel = {'windrow': 5,
                     'ASP': 2.5,
                     'in_vessel': 2.5}

# L/day
other_composting_diesel = (dewatered_solids_mass + bulking_agent)*other_unit_diesel['ASP']

# TODO: this is not included in the BEAM*2024 model
# compost volume to the dewatered solids volume, -
compost_solids_ratio = 2

# size of loads per truck, m3/load
load_size = 13

# load frequency, load/hr
load_frequency = 3

# tractor fuel, L diesel/hr
tractor_fuel = 25

# L/day
land_application_diesel = (dewatered_solids_mass/dewatered_solids_density*compost_solids_ratio/load_size/load_frequency*tractor_fuel)

# total diesel, L/day
diesel = grinding_diesel + other_composting_diesel + land_application_diesel

# tonne CO2 eq/day
diesel_emission = diesel*diesel_CI/1000000

# kWh/dry tonne solids (not including the bulking agent)
elec_unit = {'ASP': 180,
             'in_vessel': 291}

# kWh/day
elec = dry_solids*elec_unit['ASP']

# tonne CO2 eq/day
elec_emission = elec*elec_CI/1000000

# TODO: confirm these ratios are CH4:C, not CH4-C:C
# -
CH4_fugitive_ratio = {'well_aerated': 0.0001,
                      'inadequately_aerated': 0.017}

# fugitive CH4, tonne/day
CH4_fugitive = dry_solids*OC_after_digestion_ratio*CH4_fugitive_ratio['inadequately_aerated']

# tonne CO2 eq/day
CH4_fugitive_emission = CH4_fugitive*CH4_CF

# TODO: confirm these ratios are N2O:N, not N2O-N:N
# -
N2O_N = {'digested_solids': 0.00076,
         'undigested_solids': 0.018}

# TODO: use N_before_digestion_ratio for undigested_solids
# fugitive N2O, tonne/day
N2O_fugitive = dry_solids*N_after_digestion_ratio*N2O_N['digested_solids']

# tonne CO2 eq/day
N2O_fugitive_emission = N2O_fugitive*N2O_CF

# carbon sequestration, tonne CO2 eq/dry tonne solids
carbon_unit_sequestration = {'low': 0.15,
                             'mid': 0.4475,
                             'high': 0.745}

# TODO: note the BEAM*2024 model only considered carbon sequestratred by solids (not including the bulking agent)
# tonne CO2 eq/day
carbon_sequestration = -(dry_solids + bulking_agent*bulking_agent_solids_ratio)*carbon_unit_sequestration['mid']

# N fertilizer credit, tonne CO2 eq/day
N_credit = -dry_solids*N_after_digestion_ratio*N_fertilizer_CI

# P fertilizer credit, tonne CO2 eq/day
P_credit = -dry_solids*P_after_digestion_ratio*P_fertilizer_CI

# output, tonne CO2 eq/year
# scope 1
(diesel_emission + CH4_fugitive_emission + N2O_fugitive_emission + carbon_sequestration)*365
# scope 2
elec_emission*365
# scope 3
(N_credit + P_credit)*365

#%% landfilling (typical)

# tonne/day
OC = dry_solids*OC_after_digestion_ratio

# landfilling uncertainty factor, -
landfilling_uncertainty = 0.75

# CH4 in landfilling gas, -
CH4_landfilling_gas = 0.5

# CH4 correction factor, -
MCF = 1

# fraction of degradable organic carbon that can decompose, -
OC_decomposition_ratio = {'complete_digested': 0.5,
                          'partial_digested': 0.65,
                          'undigested': 0.8}

# TODO: check the unit
# dacay rate, -
k_decay = {'cool_dry': 0.06,
           'cool_wet': 0.185,
           'warm_dry': 0.085,
           'warm_wet': 0.4}

# decayed OC schedule
decayed_OC_schedule = lambda x: k_decay['cool_dry']*(1 - k_decay['cool_dry'])**x

# CH4 lost schedule
CH4_lost_schedule = {'0_1': 1,
                     '2_4': 0.5,
                     '5_14': 0.25,
                     'after_capping': 0.1}

# oxidized CH4 lost schedule
oxidized_CH4_lost_schedule = {'0_1': 0.1,
                              '2_4': 0.25,
                              '5_14': 0.25,
                              'after_capping': 0.35}

# tonne CH4/day
CH4_fugitive_years_0_1 = OC*landfilling_uncertainty*OC_decomposition_ratio['complete_digested']*\
                         sum(decayed_OC_schedule(i) for i in range(0, 2))*_C_to_CH4*CH4_landfilling_gas*\
                         CH4_lost_schedule['0_1']*MCF*(1-oxidized_CH4_lost_schedule['0_1'])

# tonne CH4/day
CH4_fugitive_years_2_4 = OC*landfilling_uncertainty*OC_decomposition_ratio['complete_digested']*\
                         sum(decayed_OC_schedule(i) for i in range(2, 5))*_C_to_CH4*CH4_landfilling_gas*\
                         CH4_lost_schedule['2_4']*MCF*(1-oxidized_CH4_lost_schedule['2_4'])

# tonne CH4/day
CH4_fugitive_years_5_14 = OC*landfilling_uncertainty*OC_decomposition_ratio['complete_digested']*\
                          sum(decayed_OC_schedule(i) for i in range(5, 15))*_C_to_CH4*CH4_landfilling_gas*\
                          CH4_lost_schedule['5_14']*MCF*(1-oxidized_CH4_lost_schedule['5_14'])

# tonne CH4/day
CH4_fugitive_years_after_capping = OC*landfilling_uncertainty*OC_decomposition_ratio['complete_digested']*\
                                   sum(decayed_OC_schedule(i) for i in range(15, 30))*_C_to_CH4*CH4_landfilling_gas*\
                                   CH4_lost_schedule['after_capping']*MCF*(1-oxidized_CH4_lost_schedule['after_capping'])

# tonne CH4/day
CH4_combustion_0_1 = OC*landfilling_uncertainty*OC_decomposition_ratio['complete_digested']*\
                     sum(decayed_OC_schedule(i) for i in range(0, 2))*_C_to_CH4*CH4_landfilling_gas*\
                     (1 - CH4_lost_schedule['0_1'])*MCF

# tonne CH4/day
CH4_combustion_2_4 = OC*landfilling_uncertainty*OC_decomposition_ratio['complete_digested']*\
                     sum(decayed_OC_schedule(i) for i in range(2, 5))*_C_to_CH4*CH4_landfilling_gas*\
                     (1 - CH4_lost_schedule['2_4'])*MCF

# tonne CH4/day
CH4_combustion_5_14 = OC*landfilling_uncertainty*OC_decomposition_ratio['complete_digested']*\
                      sum(decayed_OC_schedule(i) for i in range(5, 15))*_C_to_CH4*CH4_landfilling_gas*\
                      (1 - CH4_lost_schedule['5_14'])*MCF

# tonne CH4/day
CH4_combustion_after_capping = OC*landfilling_uncertainty*OC_decomposition_ratio['complete_digested']*\
                               sum(decayed_OC_schedule(i) for i in range(15, 30))*_C_to_CH4*CH4_landfilling_gas*\
                               (1 - CH4_lost_schedule['after_capping'])*MCF

# tonne/day
CH4_combustion = CH4_combustion_0_1 + CH4_combustion_2_4 +\
                 CH4_combustion_5_14 + CH4_combustion_after_capping

# tonne CH4/day
CH4_fugitive = CH4_fugitive_years_0_1 + CH4_fugitive_years_2_4 + CH4_fugitive_years_5_14 +\
               CH4_fugitive_years_after_capping + CH4_combustion*flare_fugitive_CH4_ratio['default']

# tonne CO2 eq/day
CH4_fugitive_emission = CH4_fugitive*CH4_CF

# -
C_N_cutoff = 30

if OC_after_digestion_ratio/N_after_digestion_ratio < C_N_cutoff:
    # -
    N2O_N_landfilling = 0.015
    # tonne N2O/day
    N2O_fugitive = dry_solids*N_after_digestion_ratio*N2O_N_landfilling*_N_to_N2O
else:
    # tonne N2O/day
    N2O_fugitive = 0

# tonne CO2 eq/day
N2O_fugitive_emission = N2O_fugitive*N2O_CF

# carbon sequestration, tonne CO2 eq/day
carbon_sequestration = -OC*(1 - OC_decomposition_ratio['complete_digested'])*_C_to_CO2

# percentage of captured CH4 used to generate electricity, -
CH4_elec = 0.5

# BTU to kWh, not at 100% efficiency
_BTU_to_kWh_partial_efficiency = 0.0000854

# kWh/day
elec = CH4_combustion*CH4_elec*1000/CH4_density*CH4_heat*_BTU_to_kWh_partial_efficiency*net_capacity_factor

# tonne CO2 eq/day
elec_emission = -elec*elec_CI/1000000

# output, tonne CO2 eq/year
# scope 1
(CH4_fugitive_emission + N2O_fugitive_emission + carbon_sequestration)*365
# scope 2
elec_emission*365

#%% combustion

# wet tonne/day
dried_solids_mass = dry_solids/solids_content_after_drying

# dry tonne/day
ash_dry = dry_solids*(1 - VS_after_digestion_ratio)

# -
ash_solids_content = 0.99

# dry tonne/day
ash_wet = ash_dry/ash_solids_content

# energy for evaporating water, BTU/day
energy_water_removal = dried_solids_mass*(1 - solids_content_after_drying)*energy_unit_water_removal*_GJ_to_BTU

# heating value of unit solids, MJ/dry tonne
heating_unit_solids = {'digested_solids': 12000,
                       'undigested_solids': 23000}

# energy potential of solids, BTU/day
energy_solids = dry_solids*heating_unit_solids['digested_solids']*_MJ_to_BTU

# multiple hearth incinerator (MHI) or fluidized bed incinerator (FBI)
incinerator = 'MHI'

if incinerator == 'MHI':
    # -
    MHI_additional_fuel = 0.2
    
    # m3 natural gas/day
    NG_use = (1 + MHI_additional_fuel)*energy_water_removal/NG_heat
else:
    # m3 natural gas/day
    NG_use = energy_water_removal/NG_heat

# -
heat_recovered_ratio = 0.5

# -
heat_recovered_efficiency = 0.8

# m3/day
NG_avoided = energy_solids*heat_recovered_ratio/NG_heat*heat_recovered_efficiency

# m3/day
NG = NG_use - NG_avoided

# tonne CO2 eq/day
NG_combustion_emission = NG*NG_CI/1000000

# kWh/dry tonne
elec_unit = {'MHI': 285,
             'FBI': 200}

# kWh/day
elec_use = dry_solids*elec_unit['MHI']

# -
elec_recovered = 0

# kWh/day
elec_avoided = energy_solids*elec_recovered*_BTU_to_kWh_partial_efficiency*net_capacity_factor

# kWh/day
elec = elec_use - elec_avoided

# tonne CO2 eq/day
elec_emission = elec*elec_CI/1000000

# TODO: note this value assumes 20% solids
# tonne CH4/dry tonne solids
CH4_combustion = 0.0000097

# TODO: note BEAM*2024 use wet tonne
# tonne CH4/day
CH4_fugitive = dry_solids*CH4_combustion

# tonne CO2 eq/day
CH4_fugitive_emission = CH4_fugitive*CH4_CF

# N, dry tonne/day
N = dry_solids*N_after_digestion_ratio

# -
suzuki_constant_1 = 161.3

# -
suzuki_constant_2 = 0.14

# °C
combustion_temperature = 850

# °C
suzuki_lowest_temperature = 750

if N*(suzuki_constant_1 - (suzuki_constant_2*(combustion_temperature + _C_to_K)))/100 < 0:
    # tonne N2O/day
    N2O_before_adjustment = 0
else:
    # tonne N2O/day
    N2O_before_adjustment = (N*(suzuki_constant_1 - (suzuki_constant_2*(max(combustion_temperature, suzuki_lowest_temperature) + _C_to_K)))/100*_N_to_N2O)

# urea catalyst
urea_catalyst = True

if urea_catalyst:
    # N20 increase ratio if urea catalyst is used, -
    N2O_adjustment_factor_urea = 0.2
    
    # tonne N2O/day
    N2O_urea_catalyst = N2O_before_adjustment*N2O_adjustment_factor_urea

if solids_content_after_drying < 0.24:
    # - 
    N2O_reduction_ratio = 0
elif solids_content_after_drying < 0.87:
    # - 
    N2O_reduction_ratio = 0.5
else:
    # - 
    N2O_reduction_ratio = 0.6

# tonne N2O/day
N2O_reduction = -N2O_before_adjustment*N2O_reduction_ratio

# tonne N2O/day
N2O_fugitive = N2O_before_adjustment + N2O_urea_catalyst + N2O_reduction

# tonne CO2 eq/day
N2O_fugitive_emission = N2O_fugitive*N2O_CF

# TODO: is it possible to recover both (by adding split ratio parameter)
# cement or phosphorus
resource_recovered = 'cement'
# cement replacement

if resource_recovered == 'cement':
    # cement unit replacement, kg CO2 eq/dry tonne solids
    cement_unit_credit = 1.2675
    
    # cement replacement, tonne CO2 eq/day
    resource_credit = -dry_solids*cement_unit_credit/1000
else:
    # phosphorus fertilizer replacement, tonne CO2 eq/day
    resource_credit = -dry_solids*P_after_digestion_ratio*P_fertilizer_CI

# output, tonne CO2 eq/year
# scope 1
(NG_combustion_emission + CH4_fugitive_emission + N2O_fugitive_emission)*365
# scope 2
elec_emission*365
# scope 3
resource_credit*365

#%% pyrolysis

# wet tonne/day
dried_solids_mass = dry_solids/solids_content_after_drying

# solids mass loss, - 
mass_loss = 0.5

# dry tonne/day
biochar_dry = dry_solids*(1 - mass_loss)

# -
biochar_solids_content = 0.99

# wet tonne/day
biochar_wet = biochar_dry/biochar_solids_content

# pyrolysis can be autothermal, so NG = 0
# m3 natural gas/day
NG = 0

# tonne CO2 eq/day
NG_combustion_emission = NG*NG_CI/1000000

# assume no electricity generated in pyrolysis
# electricity for anciliary equipment, kWh/dry tonne
elec_unit = 123.424

# kWh/day
elec = dry_solids*elec_unit

# tonne CO2 eq/day
elec_emission = elec*elec_CI/1000000

# g CH4/dry tonne
CH4_fugitive_unit = 2.65

# fugitive CH4, tonne/day
CH4_fugitive = dry_solids*CH4_fugitive_unit/1000000

# tonne CO2 eq/day
CH4_fugitive_emission = CH4_fugitive*CH4_CF

# g N2O/dry tonne
N2O_fugitive_unit = 5.23

# fugitive N2O
N2O_fugitive = dry_solids*N2O_fugitive_unit/1000000

# tonne CO2 eq/day
N2O_fugitive_emission = N2O_fugitive*N2O_CF

# output, tonne CO2 eq/year
# scope 1
(NG_combustion_emission + CH4_fugitive_emission + N2O_fugitive_emission)*365
# scope 2
elec_emission*365

#%% land application

# L/day
land_application_diesel = (dewatered_solids_mass/dewatered_solids_density/load_size/load_frequency*tractor_fuel)

# tonne CO2 eq/day
diesel_emission = land_application_diesel*diesel_CI/1000000

# minimum solids ratio above which no fugitive emission from storage
min_solids_content_no_fugitive = 0.55

# days
storage_time = 10

# CH4 during storage before land application, kg CH4/m3/day
CH4_unit_before_land_application = 0.0091

if solids_content_after_dewatering < min_solids_content_no_fugitive:
    # tonne CH4/day
    CH4_storage = dewatered_solids_mass/dewatered_solids_density*storage_time*CH4_unit_before_land_application/1000
else:
    # tonne CH4/day
    CH4_storage = 0

# tonne CH4/day
CH4_fugitive = CH4_storage

# tonne CO2 eq/day
CH4_fugitive_emission = CH4_fugitive*CH4_CF

# fine textured soils ratio, -
fine_textured_ratio = 0.5

# N2O EF from land application (fine textured soils), -
N2O_EF_fine_textured_soils = 0.0275

if OC_after_digestion_ratio/N_after_digestion_ratio < C_N_cutoff:
    # tonne N2O/day
    N2O_fine_textured_soils = dry_solids*N_after_digestion_ratio*fine_textured_ratio*N2O_EF_fine_textured_soils*_N_to_N2O
else:
    # tonne N2O/day
    N2O_fine_textured_soils = 0

# minimum solids content for N2O reduction on fine textured soils
min_solids_content_N2O_reduction = 0.8

if solids_content_after_dewatering >= min_solids_content_N2O_reduction:
    # -
    N2O_reduction_fine_textured_soils_ratio = 1
    
    # tonne N2O/day
    N2O_reduction_fine_textured_soils = -N2O_fine_textured_soils*N2O_reduction_fine_textured_soils_ratio
elif solids_content_after_dewatering > min_solids_content_no_fugitive:
    # -
    N2O_reduction_slope = 0.276
    
    # -
    N2O_reduction_intercept = 0.1518
    
    # TODO: double check the calculation here and in the BEAM*2024 model
    # tonne N2O/day
    N2O_reduction_fine_textured_soils = -N2O_fine_textured_soils*(N2O_reduction_slope*solids_content_after_dewatering - N2O_reduction_intercept)
else:
    # tonne N2O/day
    N2O_reduction_fine_textured_soils = 0

# N2O during storage before land application, kg/m3/day
N2O_unit_before_land_application = 0.00043

if solids_content_after_dewatering < min_solids_content_no_fugitive:
    # tonne N2O/day
    N2O_storage = dewatered_solids_mass/dewatered_solids_density*storage_time*N2O_unit_before_land_application/1000
else:
    # tonne N2O/day
    N2O_storage = 0

N2O_fugitive = N2O_fine_textured_soils + N2O_reduction_fine_textured_soils + N2O_storage

# climate at land application sites, can be 'humid' and 'arid'
climate_land_application = 'humid'

if climate_land_application == 'humid':
    N2O_fugitive = N2O_fine_textured_soils + N2O_reduction_fine_textured_soils + N2O_storage
else:
    N2O_fugitive = N2O_storage

# tonne CO2 eq/day
N2O_fugitive_emission = N2O_fugitive*N2O_CF

# TODO: need to match up with the solids content
# can be 'thermal_drying', 'alkaline_stabilization', 'pyrolysis', and 'other'
previous_step = 'alkaline_stablization'

if previous_step == 'pyrolysis':
    # OC sequestered ratio, -
    OC_sequestered_ratio = 0.75
    
    # tonne CO2 eq/day
    carbon_sequestration = -dry_solids*OC_after_digestion_ratio*OC_sequestered_ratio*_C_to_CO2
else:
    # tonne CO2 eq/day
    carbon_sequestration = -dry_solids*carbon_unit_sequestration['mid']

# N fertilizer credit, tonne CO2 eq/day
N_credit = -dry_solids*N_after_digestion_ratio*N_fertilizer_CI

# P fertilizer credit, tonne CO2 eq/day
P_credit = -dry_solids*P_after_digestion_ratio*P_fertilizer_CI

# TODO: update CaCO3 equivalence (%-dry weight) value (this might be related to lime stablization)
# -
CaCO3_ratio = 0.1

# -
lime_from_waste = True

# -
lime_replace_lime = True

# -
C_CaCO3 = 12/100

# debit from CaCO3 applied to soil
if lime_from_waste and lime_replace_lime:
    # CaCO3 debit, tonne CO2 eq/day
    CaCO3_debit = 0
else:
    # CaCO3 debit, tonne CO2 eq/day
    CaCO3_debit = dry_solids*CaCO3_ratio*C_CaCO3*_C_to_CO2

# output, tonne CO2 eq/year
# scope 1
(diesel_emission + CH4_fugitive_emission + N2O_fugitive_emission + carbon_sequestration + CaCO3_debit)*365
# scope 3
(N_credit + P_credit)*365