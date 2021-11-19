#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 17:06:29 2018

@author: John Trimmer; Yalin Li (minor modification)
"""

# These functions track resources, costs, and emissions for various centralized treatment options
#  (i) anaerobic_digestion
#  (ii) sedimentation (used at the existing Lubigi plant)
#  (iii) sludge_separator (used to estimate solid/liquid separation after anaerobic treatment in the alternative plant)
#  (iv) anaerobic_lagoon (used at the existing Lubigi plant)
#  (v) facultative_lagoon (used at the existing Lubigi plant)
#  (vi) unplanted_drying_bed (used at the existing Lubigi plant)
#  (vii) drying_beds_alt (drying beds of a different design for the anaerobic treatment plant)
#  (viii) ABR (anaerobic  baffled reactor used in the alternative plant)\
#  (ix) secondary_liquid_bed (planted bed used for treatment of ABR liquid effluent in alternative plant)

import numpy as np
import pandas as pd
import copy
import lhs
import math

def first_order_decay(k, t0, t, max_decay, tot=1):
    tf = t0 + t
    Cdeg = tot * max_decay
    Cavg = Cdeg/(k*t) * (np.exp(-k*t0)-np.exp(-k*tf))
    loss = Cdeg - Cavg
    return loss


#%% anaerobic digestion function

def anaerobic_digestion(inputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, biogas, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, CH4_energy, previous_storage_time, additional_storage, sludge_flow_alt, concrete_thickness, concrete_IF_GHG, excavation_IF_GHG, plant_lifetime, sludge_pop_alt, discount_rate):
    outputs = copy.deepcopy(inputs)

    mass = np.reshape(inputs[:,0], (-1,1))
    mass_dry = np.reshape(inputs[:,1], (-1,1))
    N_total = np.reshape(inputs[:,2], (-1,1))
    P_total = np.reshape(inputs[:,3], (-1,1))
    K_total = np.reshape(inputs[:,4], (-1,1))
    Mg_total = np.reshape(inputs[:,5], (-1,1))
    Ca_total = np.reshape(inputs[:,6], (-1,1))
    energy = np.reshape(inputs[:,7], (-1,1))
    N_ammonia = np.reshape(inputs[:,8], (-1,1))

    # biogas production
    # parameters
    MCF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_AD, correlation_distributions, correlation_parameters, n_samples)
    COD_removal_AD, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_removal_AD, correlation_distributions, correlation_parameters, n_samples)
    residence_time_AD, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.residence_time_AD, correlation_distributions, correlation_parameters, n_samples)

    # calculations
    # carbon (kg COD/yr; assume 14 kJ/g COD in wastewater)
    COD_total = (energy/14/1000)
    COD_degrade = COD_total*(COD_removal_AD/100)
    CH4_production = COD_degrade*(MCF/100)*maximum_methane_emission
    COD_total = COD_total - COD_degrade
    energy = COD_total*14*1000
    if parameters.CH4_captured_AD.expected == 'yes':
        CH4eq = np.full([n_samples, 1], 0)
    else:
        CH4eq = CH4_production*CH4_GWP
        CH4_production = 0

    # nitrogen emissions (kg N/yr; N2O expressed as kg N2O/yr)
    if parameters.N_emission_in_biogas_AD.expected == 'yes':
        N2O_EF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_AD, correlation_distributions, correlation_parameters, n_samples)
        N_denitrification, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_max_denitrification_AD, correlation_distributions, correlation_parameters, n_samples)

        #!!! Changed algorithm
        t0 = previous_storage_time + additional_storage
        t = residence_time_AD / 365
        N_loss = first_order_decay(rate_constant, t0, t,
                                    max_decay=N_denitrification/100,
                                    tot=N_total)

        # N_degradable = N_total * (N_denitrification/100)
        # N_initial = (N_degradable*rate_constant*(previous_storage_time))/(np.exp(-rate_constant*(additional_storage)) - np.exp(-rate_constant*(additional_storage+previous_storage_time)))
        # N_after = (N_initial/(rate_constant*(previous_storage_time)))*(np.exp(-rate_constant*(additional_storage+residence_time_AD/365)) - np.exp(-rate_constant*(previous_storage_time + additional_storage + residence_time_AD/365)))
        # N_loss = N_degradable - N_after

        N2O_emission = N_loss*(N2O_EF/100)*(44/28)
        N2Oeq = N2O_emission*N2O_GWP
        N_total = N_total - N_loss
        N_ammonia = N_ammonia - N_loss
        for i in range(0, len(N_ammonia)):
            if N_ammonia[i] < 0:
                N_ammonia[i] = 0
    else:
        N2Oeq = np.full([n_samples, 1], 0)


    # solids loss (based on COD loss)
    mass = mass - mass_dry * (COD_removal_AD/100)
    mass_dry = mass_dry - mass_dry * (COD_removal_AD/100)

    # energy in biogas (kJ/yr)
    biogas_new = CH4_production / 16 * CH4_energy * 1000
    biogas = biogas + biogas_new

    outputs[:,0:9] = np.concatenate((mass, mass_dry, N_total, P_total, K_total,
                                     Mg_total, Ca_total, energy, N_ammonia), 1)


    # construction costs and emissions
    number, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.number_AD, correlation_distributions, correlation_parameters, n_samples)
    aspect_ratio, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.aspect_ratio_AD, correlation_distributions, correlation_parameters, n_samples)
    headspace, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.headspace_AD, correlation_distributions, correlation_parameters, n_samples)

    reactor_volume = ((residence_time_AD*sludge_flow_alt)/(number-1))/(1 - headspace/100)
    reactor_diameter = ((4*reactor_volume*aspect_ratio)/math.pi)**(1/3)
    reactor_height = reactor_diameter/aspect_ratio
    concrete_volume = number*(concrete_thickness*((2*(math.pi/4)*(reactor_diameter**2))+(reactor_height*math.pi*reactor_diameter)))
    concrete_emissions = concrete_volume * concrete_IF_GHG

    excavation_volume = reactor_volume * number
    excavation_emissions = excavation_volume * excavation_IF_GHG

    construction_emissions_annual = (concrete_emissions + excavation_emissions)/plant_lifetime/sludge_pop_alt

    tech_construction_emissions[:,3:4] = tech_construction_emissions[:,3:4] + construction_emissions_annual


    direct_emissions[:,3:4] = direct_emissions[:,3:4] + CH4eq + N2Oeq

    # cost
    if parameters.use_total_price.expected == 'no':
        concrete_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.concrete_cost, correlation_distributions, correlation_parameters, n_samples)
        steel_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.steel_cost, correlation_distributions, correlation_parameters, n_samples)
        construction_cost = (concrete_volume + concrete_volume)*concrete_cost

        opex_percent, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.alternative_plant_opex, correlation_distributions, correlation_parameters, n_samples)
        annual_opex = (construction_cost*(opex_percent/100)) / sludge_pop_alt
        capex_annualized = (construction_cost * ((discount_rate*(1 + discount_rate)**plant_lifetime)/(((1 + discount_rate)**plant_lifetime) - 1))) / sludge_pop_alt

        construction_cost[:,3:4] = construction_cost[:,3:4] + capex_annualized
        operating_cost[:,3:4] = operating_cost[:,3:4] + annual_opex

        electricity_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_cost, correlation_distributions, correlation_parameters, n_samples)
        electricity_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_GHG, correlation_distributions, correlation_parameters, n_samples)
        # assume all operating expenses come from electricity
        annual_kWh = annual_opex / electricity_cost
        operating_emissions = (annual_kWh * electricity_GHG)
        tech_operating_emissions[:,3:4] = tech_operating_emissions[:,3:4] + operating_emissions


    additional_storage = additional_storage + residence_time_AD/365

    return outputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, biogas, correlation_distributions, correlation_parameters, additional_storage

#%% sedimentation function
def sedimentation(inputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage, concrete_thickness, roof_slope, roof_mass, concrete_IF_GHG, steel_IF_GHG, plant_lifetime, existing_population, discount_rate):
    solids_outputs = copy.deepcopy(inputs)
    liquid_outputs = copy.deepcopy(inputs)

    mass = np.reshape(inputs[:,0], (-1,1))
    mass_dry = np.reshape(inputs[:,1], (-1,1))
    N_total = np.reshape(inputs[:,2], (-1,1))
    P_total = np.reshape(inputs[:,3], (-1,1))
    K_total = np.reshape(inputs[:,4], (-1,1))
    Mg_total = np.reshape(inputs[:,5], (-1,1))
    Ca_total = np.reshape(inputs[:,6], (-1,1))
    energy = np.reshape(inputs[:,7], (-1,1))
    N_ammonia = np.reshape(inputs[:,8], (-1,1))

    # calculate retention in settled solids (before degradation)
    COD_retention_sedimentation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_retention_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    TS_retention_sedimentation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.TS_retention_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    N_retention_sedimentation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_retention_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    P_retention_sedimentation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.P_retention_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    K_retention_sedimentation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.K_retention_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    Mg_retention_sedimentation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.Mg_retention_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    Ca_retention_sedimentation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.Ca_retention_sedimentation, correlation_distributions, correlation_parameters, n_samples)

    solids_COD = (energy * (COD_retention_sedimentation/100))/14/1000
    solids_TS = mass_dry * (TS_retention_sedimentation/100)
    solids_N = N_total * (N_retention_sedimentation/100)
    solids_P = P_total * (P_retention_sedimentation/100)
    solids_K = K_total * (K_retention_sedimentation/100)
    solids_Mg = Mg_total * (Mg_retention_sedimentation/100)
    solids_Ca = Ca_total * (Ca_retention_sedimentation/100)

    liquid_COD = (energy * ((100 - COD_retention_sedimentation)/100))/14/1000
    liquid_TS = mass_dry * ((100 - TS_retention_sedimentation)/100)
    liquid_N = N_total * ((100 - N_retention_sedimentation)/100)
    liquid_P = P_total * ((100 - P_retention_sedimentation)/100)
    liquid_K = K_total * ((100 - K_retention_sedimentation)/100)
    liquid_Mg = Mg_total * ((100 - Mg_retention_sedimentation)/100)
    liquid_Ca = Ca_total * ((100 - Ca_retention_sedimentation)/100)

    # assume as much ammonia as possible drains with liquid
    liquid_N_ammonia = np.full(np.shape(N_ammonia), np.nan)
    solids_N_ammonia = np.full(np.shape(N_ammonia), np.nan)
    N_ammonia_temp = copy.deepcopy(N_ammonia)
    liquid_N_temp = copy.deepcopy(liquid_N)
    for i in range (0, len(N_ammonia_temp)):
        if N_ammonia_temp[i] <= liquid_N[i]:
            liquid_N_ammonia[i] = N_ammonia_temp[i]
            solids_N_ammonia[i] = 0
        else:
            liquid_N_ammonia[i] = liquid_N_temp[i]
            solids_N_ammonia[i] = N_ammonia_temp[i] - liquid_N_temp[i]

    # calculate degradation in settled solids
    MCF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    retention_time_sedimentation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.retention_time_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    COD_degradation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_degradation_sedimentation, correlation_distributions, correlation_parameters, n_samples)

    #!!! Changed algorithm
    t0 = previous_storage_time + additional_storage
    t = retention_time_sedimentation / 365
    COD_loss = first_order_decay(rate_constant, t0, t,
                                 max_decay=COD_degradation/100,
                                 tot=solids_COD)

    # COD_degradable = solids_COD * COD_degradation/100
    # COD_initial = (COD_degradable*rate_constant*previous_storage_time)/(np.exp(-rate_constant*(additional_storage)) - np.exp(-rate_constant*(additional_storage+previous_storage_time)))
    # COD_after = (COD_initial/(rate_constant*(previous_storage_time)))*(np.exp(-rate_constant*(additional_storage+retention_time_sedimentation/365)) - np.exp(-rate_constant*(previous_storage_time + additional_storage + retention_time_sedimentation/365)))
    # COD_loss = COD_degradable - COD_after

    COD_reduction = COD_loss/solids_COD
    CH4_emission = COD_loss*(MCF/100)*maximum_methane_emission
    CH4eq = CH4_emission*CH4_GWP
    solids_COD = solids_COD - COD_loss
    solids_energy = solids_COD*14*1000
    liquid_energy = liquid_COD*14*1000
    #solids loss from COD degradation
    mass = mass - solids_TS*COD_reduction
    solids_TS = solids_TS - solids_TS*COD_reduction

    if parameters.N_emission_from_sedimentation.expected == 'yes':
        N2O_EF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_sedimentation, correlation_distributions, correlation_parameters, n_samples)
        N_denitrification, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_max_denitrification_sedimentation, correlation_distributions, correlation_parameters, n_samples)

        #!!! Changed algorithm
        t0 = previous_storage_time + additional_storage
        t = retention_time_sedimentation / 365
        N_loss = first_order_decay(rate_constant, t0, t,
                                    max_decay=N_denitrification/100,
                                    tot=N_total)
        N2O_emission = N_loss*(N2O_EF/100)*(44/28)

        # N_degradable = solids_N * (N_denitrification/100)
        # N_initial = (N_degradable*rate_constant*(previous_storage_time))/(np.exp(-rate_constant*(additional_storage)) - np.exp(-rate_constant*(additional_storage+previous_storage_time)))
        # N_after = (N_initial/(rate_constant*(previous_storage_time)))*(np.exp(-rate_constant*(additional_storage+retention_time_sedimentation/365)) - np.exp(-rate_constant*(previous_storage_time + additional_storage + retention_time_sedimentation/365)))
        # N_loss = N_degradable - N_after
        # N2O_emission = solids_N*(N_loss/N_degradable)*(N2O_EF/100)*(44/28)

        for i in range(0, len(N2O_emission)):
            if N2O_emission[i] > N_loss[i]*(44/28):
                N2O_emission[i] = N_loss[i]*(44/28)
        N2Oeq = N2O_emission*N2O_GWP
        solids_N = solids_N - N_loss
        solids_N_ammonia = solids_N_ammonia - N_loss
        for i in range(0, len(solids_N_ammonia)):
            if solids_N_ammonia[i] < 0:
                solids_N_ammonia[i] = 0
    else:
        N2Oeq = np.full([n_samples, 1], 0)

    # calculate total mass of settled solids based on final solids content
    final_solids_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.final_solids_content_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    liquid = np.full(np.shape(final_solids_content), np.nan)
    solids = np.full(np.shape(final_solids_content), np.nan)
    for i in range(0, len(solids_TS)):
        if final_solids_content[i] > ((solids_TS[i]/mass[i])*100):
            solids[i] = solids_TS[i] / (final_solids_content[i]/100)
            drained_water = mass[i] - solids[i]
            liquid[i] = drained_water + liquid_TS[i]
        else:
            solids[i] = mass[i]
            liquid[i] = liquid_TS[i]

    # construction costs and emissions
    volume, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.volume_sedimentation_tank, correlation_distributions, correlation_parameters, n_samples)
    length_width_ratio, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.length_width_ratio_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    width_height_ratio, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.width_height_ratio_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    number, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.number_sedimentation_tanks, correlation_distributions, correlation_parameters, n_samples)
    columns_per_side, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.columns_per_side_sedimentation, correlation_distributions, correlation_parameters, n_samples)

    average_height = (volume/(length_width_ratio*(width_height_ratio**2)))**(1/3)
    width = average_height * width_height_ratio
    length = width * length_width_ratio
    concrete_volume = (concrete_thickness*(width*length + 2*width*average_height + 2*length*average_height))*number
    concrete_columns = (concrete_thickness**2)*(average_height)*(columns_per_side*2)*number
    concrete_emissions = (concrete_volume + concrete_columns)*concrete_IF_GHG

    roof_area = np.full(np.shape(roof_slope), np.nan)
    for i in range(0,len(roof_slope)):
        roof_area[i] = (number[i]*length[i]*width[i])/math.cos(roof_slope[i]*math.pi/180)
    siding_area = number*(2*length*average_height + 2*width*average_height)
    mass_steel = (roof_area + siding_area)*roof_mass
    steel_emissions = mass_steel*steel_IF_GHG
    construction_emissions_annual = (concrete_emissions + steel_emissions)/plant_lifetime/existing_population

    tech_construction_emissions[:,3:4] = tech_construction_emissions[:,3:4] + construction_emissions_annual

    solids_outputs[:,0:9] = np.concatenate((solids, solids_TS, solids_N, solids_P, solids_K,
                                     solids_Mg, solids_Ca, solids_energy, solids_N_ammonia), 1)
    liquid_outputs[:,0:9] = np.concatenate((liquid, liquid_TS, liquid_N, liquid_P, liquid_K,
                                     liquid_Mg, liquid_Ca, liquid_energy, liquid_N_ammonia), 1)

    direct_emissions[:,3:4] = direct_emissions[:,3:4] + CH4eq + N2Oeq

    # cost
    if parameters.use_total_price.expected == 'no':
        concrete_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.concrete_cost, correlation_distributions, correlation_parameters, n_samples)
        steel_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.steel_cost, correlation_distributions, correlation_parameters, n_samples)
        cap_cost = (concrete_volume + concrete_columns)*concrete_cost + mass_steel*steel_cost

        opex_percent, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.existing_plant_opex, correlation_distributions, correlation_parameters, n_samples)
        annual_opex = (cap_cost*(opex_percent/100)) / existing_population
        capex_annualized = (cap_cost * ((discount_rate*(1 + discount_rate)**plant_lifetime)/(((1 + discount_rate)**plant_lifetime) - 1))) / existing_population

        construction_cost[:,3:4] = construction_cost[:,3:4] + capex_annualized
        operating_cost[:,3:4] = operating_cost[:,3:4] + annual_opex

        electricity_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_cost, correlation_distributions, correlation_parameters, n_samples)
        electricity_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_GHG, correlation_distributions, correlation_parameters, n_samples)
        # assume all operating expenses come from electricity
        annual_kWh = annual_opex / electricity_cost
        operating_emissions = (annual_kWh * electricity_GHG)
        tech_operating_emissions[:,3:4] = tech_operating_emissions[:,3:4] + operating_emissions


    solid_additional_storage = additional_storage + retention_time_sedimentation/365
    liquid_additional_storage = additional_storage

    return liquid_outputs, solids_outputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, liquid_additional_storage, solid_additional_storage

#%% sludge separator function
def sludge_separator(inputs, direct_emissions, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage):
    solids_outputs = copy.deepcopy(inputs)
    liquid_outputs = copy.deepcopy(inputs)

    mass = np.reshape(inputs[:,0], (-1,1))
    mass_dry = np.reshape(inputs[:,1], (-1,1))
    N_total = np.reshape(inputs[:,2], (-1,1))
    P_total = np.reshape(inputs[:,3], (-1,1))
    K_total = np.reshape(inputs[:,4], (-1,1))
    Mg_total = np.reshape(inputs[:,5], (-1,1))
    Ca_total = np.reshape(inputs[:,6], (-1,1))
    energy = np.reshape(inputs[:,7], (-1,1))
    N_ammonia = np.reshape(inputs[:,8], (-1,1))

    # calculate retention in settled solids (before degradation)
    COD_retention, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_retention_separator, correlation_distributions, correlation_parameters, n_samples)
    TS_retention, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.TS_retention_separator, correlation_distributions, correlation_parameters, n_samples)
    N_retention, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_retention_separator, correlation_distributions, correlation_parameters, n_samples)
    P_retention, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.P_retention_separator, correlation_distributions, correlation_parameters, n_samples)
    K_retention, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.K_retention_separator, correlation_distributions, correlation_parameters, n_samples)
    Mg_retention, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.Mg_retention_separator, correlation_distributions, correlation_parameters, n_samples)
    Ca_retention, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.Ca_retention_separator, correlation_distributions, correlation_parameters, n_samples)

    solids_COD = (energy * (COD_retention/100))/14/1000
    solids_TS = mass_dry * (TS_retention/100)
    solids_N = N_total * (N_retention/100)
    solids_P = P_total * (P_retention/100)
    solids_K = K_total * (K_retention/100)
    solids_Mg = Mg_total * (Mg_retention/100)
    solids_Ca = Ca_total * (Ca_retention/100)

    liquid_COD = (energy * ((100 - COD_retention)/100))/14/1000
    liquid_TS = mass_dry * ((100 - TS_retention)/100)
    liquid_N = N_total * ((100 - N_retention)/100)
    liquid_P = P_total * ((100 - P_retention)/100)
    liquid_K = K_total * ((100 - K_retention)/100)
    liquid_Mg = Mg_total * ((100 - Mg_retention)/100)
    liquid_Ca = Ca_total * ((100 - Ca_retention)/100)

    # assume as much ammonia as possible drains with liquid
    liquid_N_ammonia = np.full(np.shape(N_ammonia), np.nan)
    solids_N_ammonia = np.full(np.shape(N_ammonia), np.nan)
    N_ammonia_temp = copy.deepcopy(N_ammonia)
    liquid_N_temp = copy.deepcopy(liquid_N)
    for i in range (0, len(N_ammonia_temp)):
        if N_ammonia_temp[i] <= liquid_N[i]:
            liquid_N_ammonia[i] = N_ammonia_temp[i]
            solids_N_ammonia[i] = 0
        else:
            liquid_N_ammonia[i] = liquid_N_temp[i]
            solids_N_ammonia[i] = N_ammonia_temp[i] - liquid_N_temp[i]

    solids_energy = solids_COD*14*1000
    liquid_energy = liquid_COD*14*1000

    # calculate total mass of settled solids based on final solids content
    final_solids_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.final_solids_content_sedimentation, correlation_distributions, correlation_parameters, n_samples)
    liquid = np.full(np.shape(final_solids_content), np.nan)
    solids = np.full(np.shape(final_solids_content), np.nan)
    for i in range(0, len(solids_TS)):
        if final_solids_content[i] > ((solids_TS[i]/mass[i])*100):
            solids[i] = solids_TS[i] / (final_solids_content[i]/100)
            drained_water = mass[i] - solids[i]
            liquid[i] = drained_water + liquid_TS[i]
        else:
            solids[i] = mass[i]
            liquid[i] = liquid_TS[i]

    solids_outputs[:,0:9] = np.concatenate((solids, solids_TS, solids_N, solids_P, solids_K,
                                     solids_Mg, solids_Ca, solids_energy, solids_N_ammonia), 1)
    liquid_outputs[:,0:9] = np.concatenate((liquid, liquid_TS, liquid_N, liquid_P, liquid_K,
                                     liquid_Mg, liquid_Ca, liquid_energy, liquid_N_ammonia), 1)

    liquid_additional_storage = copy.deepcopy(additional_storage)
    solid_additional_storage = copy.deepcopy(additional_storage)

    return liquid_outputs, solids_outputs, direct_emissions, correlation_distributions, correlation_parameters, liquid_additional_storage, solid_additional_storage

#%% anaerobic lagoon function
def anaerobic_lagoon(inputs, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage, flow_rate, excavation_IF_GHG, liner_mass, liner_IF_GHG, plant_lifetime, existing_population):
    outputs = copy.deepcopy(inputs)

    mass = np.reshape(inputs[:,0], (-1,1))
    mass_dry = np.reshape(inputs[:,1], (-1,1))
    N_total = np.reshape(inputs[:,2], (-1,1))
    P_total = np.reshape(inputs[:,3], (-1,1))
    K_total = np.reshape(inputs[:,4], (-1,1))
    Mg_total = np.reshape(inputs[:,5], (-1,1))
    Ca_total = np.reshape(inputs[:,6], (-1,1))
    energy = np.reshape(inputs[:,7], (-1,1))
    N_ammonia = np.reshape(inputs[:,8], (-1,1))

    # parameters
    COD_removal, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_removal_anaerobic_lagoon, correlation_distributions, correlation_parameters, n_samples)
    COD_degradation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_degradation_anaerobic_lagoon, correlation_distributions, correlation_parameters, n_samples)
    MCF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_anaerobic_lagoon, correlation_distributions, correlation_parameters, n_samples)
    volume, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.volume_anaerobic_lagoon, correlation_distributions, correlation_parameters, n_samples)
    number, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.number_anaerobic_lagoons, correlation_distributions, correlation_parameters, n_samples)
    length, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.length_anaerobic_lagoon, correlation_distributions, correlation_parameters, n_samples)
    width, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.width_anaerobic_lagoon, correlation_distributions, correlation_parameters, n_samples)

    # calculations
    retention_time = (volume*number)/flow_rate
    COD_influent = energy/14/1000
    COD_removed = COD_influent*(COD_removal/100)
    mass = mass - mass_dry*(COD_removal/100)
    mass_dry = mass_dry - mass_dry*(COD_removal/100)   # assume solids removal is similar to COD removal
    COD_loss = COD_removed * COD_degradation/100
    CH4_emission = COD_loss*(MCF/100)*maximum_methane_emission
    CH4eq = CH4_emission*CH4_GWP
    COD_effluent = COD_influent - COD_removed
    energy = COD_effluent*14*1000

    if parameters.N_emission_from_anaerobic_lagoon.expected == 'yes':
        N_denitrification, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_max_denitrification_anaerobic_lagoon, correlation_distributions, correlation_parameters, n_samples)
        N2O_EF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_anaerobic_lagoon, correlation_distributions, correlation_parameters, n_samples)

        #!!! Changed algorithm
        t0 = previous_storage_time + additional_storage
        t = retention_time / 365
        N_loss = first_order_decay(rate_constant, t0, t,
                                    max_decay=N_denitrification/100,
                                    tot=N_total)
        N2O_emission = N_loss*(N2O_EF/100)*(44/28)

        # N_degradable = N_total * (N_denitrification/100)
        # N_initial = (N_degradable*rate_constant*previous_storage_time)/(np.exp(-rate_constant*additional_storage) - np.exp(-rate_constant*(additional_storage+previous_storage_time)))
        # N_after = (N_initial/(rate_constant*(previous_storage_time)))*(np.exp(-rate_constant*(additional_storage+retention_time/365)) - np.exp(-rate_constant*(previous_storage_time + additional_storage + retention_time/365)))
        # N_loss = N_degradable - N_after
        # N2O_emission = N_total*(N_loss/N_degradable)*(N2O_EF/100)*(44/28)

        for i in range(0, len(N2O_emission)):
            if N2O_emission[i] > N_loss[i]*(44/28):
                N2O_emission[i] = N_loss[i]*(44/28)
        N2Oeq = N2O_emission*N2O_GWP
        N_total = N_total - N_loss
        N_ammonia = N_ammonia - N_loss
        for i in range(0, len(N_ammonia)):
            if N_ammonia[i] < 0:
                N_ammonia[i] = 0
    else:
        N2Oeq = np.full((n_samples,1), 0)

    # construction costs and emissions
    average_depth = volume/(length*width)
    liner_area = ((length*width)+(average_depth*(2*length+2*width)))*number

    excavation_emissions = (volume*number)*excavation_IF_GHG
    liner_emissions = liner_area*liner_mass*liner_IF_GHG

    construction_emissions_annual = (excavation_emissions+liner_emissions)/plant_lifetime/existing_population

    tech_construction_emissions[:,3:4] = tech_construction_emissions[:,3:4] + construction_emissions_annual

    outputs[:,0:9] = np.concatenate((mass, mass_dry, N_total, P_total, K_total,
                                     Mg_total, Ca_total, energy, N_ammonia), 1)

    direct_emissions[:,3:4] = direct_emissions[:,3:4] + CH4eq + N2Oeq

    additional_storage = additional_storage + retention_time/365

    return outputs, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, additional_storage

#%% facultative lagoon function
def facultative_lagoon(inputs, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage, flow_rate, excavation_IF_GHG, liner_mass, liner_IF_GHG, plant_lifetime, existing_population):
    outputs = copy.deepcopy(inputs)

    mass = np.reshape(inputs[:,0], (-1,1))
    mass_dry = np.reshape(inputs[:,1], (-1,1))
    N_total = np.reshape(inputs[:,2], (-1,1))
    P_total = np.reshape(inputs[:,3], (-1,1))
    K_total = np.reshape(inputs[:,4], (-1,1))
    Mg_total = np.reshape(inputs[:,5], (-1,1))
    Ca_total = np.reshape(inputs[:,6], (-1,1))
    energy = np.reshape(inputs[:,7], (-1,1))
    N_ammonia = np.reshape(inputs[:,8], (-1,1))

    # parameters
    COD_removal, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_removal_facultative_lagoon, correlation_distributions, correlation_parameters, n_samples)
    COD_degradation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_degradation_facultative_lagoon, correlation_distributions, correlation_parameters, n_samples)
    MCF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_facultative_lagoon, correlation_distributions, correlation_parameters, n_samples)
    N2O_EF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_facultative_lagoon, correlation_distributions, correlation_parameters, n_samples)
    N_denitrification, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_max_denitrification_facultative_lagoon, correlation_distributions, correlation_parameters, n_samples)
    P_removal, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.P_removal_facultative_lagoon, correlation_distributions, correlation_parameters, n_samples)
    volume, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.volume_facultative_lagoon, correlation_distributions, correlation_parameters, n_samples)
    number, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.number_facultative_lagoons, correlation_distributions, correlation_parameters, n_samples)
    length, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.length_facultative_lagoon, correlation_distributions, correlation_parameters, n_samples)
    width, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.width_facultative_lagoon, correlation_distributions, correlation_parameters, n_samples)

    # calculations
    retention_time = (volume*number)/flow_rate
    COD_influent = energy/14/1000
    COD_removed = COD_influent*(COD_removal/100)
    mass = mass - mass_dry*(COD_removal/100)
    mass_dry = mass_dry - mass_dry*(COD_removal/100)   # assume solids removal is similar to COD removal
    COD_loss = COD_removed * COD_degradation/100
    CH4_emission = COD_loss*(MCF/100)*maximum_methane_emission
    CH4eq = CH4_emission*CH4_GWP
    COD_effluent = COD_influent - COD_removed
    energy = COD_effluent*14*1000

    #!!! Changed algorithm
    t0 = previous_storage_time + additional_storage
    t = retention_time / 365
    N_loss = first_order_decay(rate_constant, t0, t,
                               max_decay=N_denitrification/100,
                               tot=N_total)
    N2O_emission = N_loss*(N2O_EF/100)*(44/28)

    # N_degradable = N_total * (N_denitrification/100)
    # N_initial = (N_degradable*rate_constant*previous_storage_time)/(np.exp(-rate_constant*additional_storage) - np.exp(-rate_constant*(additional_storage+previous_storage_time)))
    # N_after = (N_initial/(rate_constant*(previous_storage_time)))*(np.exp(-rate_constant*(additional_storage+retention_time/365)) - np.exp(-rate_constant*(previous_storage_time + additional_storage + retention_time/365)))
    # N_loss = N_degradable - N_after
    # N2O_emission = N_total*(N_loss/N_degradable)*(N2O_EF/100)*(44/28)

    for i in range(0, len(N2O_emission)):
        if N2O_emission[i] > N_loss[i]*(44/28):
            N2O_emission[i] = N_loss[i]*(44/28)
    N2Oeq = N2O_emission*N2O_GWP
    N_total = N_total - N_loss
    N_ammonia = N_ammonia - N_loss
    for i in range(0, len(N_ammonia)):
        if N_ammonia[i] < 0:
            N_ammonia[i] = 0

    P_loss = P_total * (P_removal/100)
    P_total = P_total - P_loss

    # construction costs and emissions
    average_depth = volume/(length*width)
    liner_area = ((length*width)+(average_depth*(2*length+2*width)))*number

    excavation_emissions = (volume*number)*excavation_IF_GHG
    liner_emissions = liner_area*liner_mass*liner_IF_GHG

    construction_emissions_annual = (excavation_emissions+liner_emissions)/plant_lifetime/existing_population

    tech_construction_emissions[:,3:4] = tech_construction_emissions[:,3:4] + construction_emissions_annual

    outputs[:,0:9] = np.concatenate((mass, mass_dry, N_total, P_total, K_total,
                                     Mg_total, Ca_total, energy, N_ammonia), 1)

    additional_storage = additional_storage + retention_time/365

    direct_emissions[:,3:4] = direct_emissions[:,3:4] + CH4eq + N2Oeq

    return outputs, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, additional_storage

#%% unplanted drying bed function
def unplanted_drying_bed(inputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage, concrete_thickness, roof_slope, roof_mass, concrete_IF_GHG, steel_IF_GHG, plant_lifetime, existing_population, discount_rate):
    outputs = copy.deepcopy(inputs)

    mass = np.reshape(inputs[:,0], (-1,1))
    mass_dry = np.reshape(inputs[:,1], (-1,1))
    N_total = np.reshape(inputs[:,2], (-1,1))
    P_total = np.reshape(inputs[:,3], (-1,1))
    K_total = np.reshape(inputs[:,4], (-1,1))
    Mg_total = np.reshape(inputs[:,5], (-1,1))
    Ca_total = np.reshape(inputs[:,6], (-1,1))
    energy = np.reshape(inputs[:,7], (-1,1))
    N_ammonia = np.reshape(inputs[:,8], (-1,1))

    # parameters
    COD_degradation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_degradation_drying_bed, correlation_distributions, correlation_parameters, n_samples)
    retention_time, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.retention_time_drying_bed, correlation_distributions, correlation_parameters, n_samples)
    MCF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_drying_bed, correlation_distributions, correlation_parameters, n_samples)
    N2O_EF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_drying_bed, correlation_distributions, correlation_parameters, n_samples)
    N_denitrification, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_max_denitrification_drying_bed, correlation_distributions, correlation_parameters, n_samples)
    final_solids_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.final_solids_content_drying_bed, correlation_distributions, correlation_parameters, n_samples)

    # calculations
    COD_influent = energy/14/1000

    #!!! Changed algorithm
    t0 = previous_storage_time + additional_storage
    t = retention_time / 365
    COD_loss = first_order_decay(rate_constant, t0, t,
                                 max_decay=COD_degradation/100,
                                 tot=COD_influent)

    # COD_degradable = COD_influent * COD_degradation/100
    # COD_initial = (COD_degradable*rate_constant*previous_storage_time)/(np.exp(-rate_constant*additional_storage) - np.exp(-rate_constant*(additional_storage+previous_storage_time)))
    # COD_after = (COD_initial/(rate_constant*(previous_storage_time)))*(np.exp(-rate_constant*(additional_storage+retention_time/365)) - np.exp(-rate_constant*(previous_storage_time + additional_storage + retention_time/365)))
    # COD_loss = COD_degradable - COD_after

    COD_reduction = COD_loss/COD_influent
    CH4_emission = COD_loss*(MCF/100)*maximum_methane_emission
    CH4eq = CH4_emission*CH4_GWP
    COD_effluent = COD_influent - COD_loss
    energy = COD_effluent*14*1000
    mass = mass - mass_dry*COD_reduction
    mass_dry = mass_dry - mass_dry*COD_reduction

    #!!! Changed algorithm
    N_loss = first_order_decay(rate_constant, t0, t,
                               max_decay=N_denitrification/100,
                               tot=N_total)
    N2O_emission = N_loss*(N2O_EF/100)*(44/28)

    # N_degradable = N_total * (N_denitrification/100)
    # N_initial = (N_degradable*rate_constant*previous_storage_time)/(np.exp(-rate_constant*additional_storage) - np.exp(-rate_constant*(additional_storage+previous_storage_time)))
    # N_after = (N_initial/(rate_constant*(previous_storage_time)))*(np.exp(-rate_constant*(additional_storage+retention_time/365)) - np.exp(-rate_constant*(previous_storage_time + additional_storage + retention_time/365)))
    # N_loss = N_degradable - N_after
    # N2O_emission = N_total*(N_loss/N_degradable)*(N2O_EF/100)*(44/28)

    for i in range(0, len(N2O_emission)):
        if N2O_emission[i] > N_loss[i]*(44/28):
            N2O_emission[i] = N_loss[i]*(44/28)
    N2Oeq = N2O_emission*N2O_GWP
    N_total = N_total - N_loss
    N_ammonia = N_ammonia - N_loss
    for i in range(0, len(N_ammonia)):
        if N_ammonia[i] < 0:
            N_ammonia[i] = 0

    for i in range(0, len(mass_dry)):
        if final_solids_content[i] > ((mass_dry[i]/mass[i])*100):
            mass[i] = mass_dry[i] / (final_solids_content[i]/100)

    # construction costs and emissions
    number_storage, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.number_storage_beds, correlation_distributions, correlation_parameters, n_samples)
    number_covered, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.number_covered_drying_beds, correlation_distributions, correlation_parameters, n_samples)
    number_uncovered, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.number_uncovered_drying_beds, correlation_distributions, correlation_parameters, n_samples)
    height_storage, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.storage_wall_height, correlation_distributions, correlation_parameters, n_samples)
    width_covered, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.covered_bed_width, correlation_distributions, correlation_parameters, n_samples)
    length_covered, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.covered_bed_length, correlation_distributions, correlation_parameters, n_samples)
    height_drying, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.drying_bed_wall_height, correlation_distributions, correlation_parameters, n_samples)
    width_uncovered, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.uncovered_bed_width, correlation_distributions, correlation_parameters, n_samples)
    length_uncovered, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.uncovered_bed_length, correlation_distributions, correlation_parameters, n_samples)
    column_mass_per_meter, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.column_mass_per_meter, correlation_distributions, correlation_parameters, n_samples)
    covered_columns_per_side, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.covered_columns_per_side, correlation_distributions, correlation_parameters, n_samples)
    covered_column_height, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.covered_column_height, correlation_distributions, correlation_parameters, n_samples)

    concrete_volume_covered = number_covered*(concrete_thickness*(length_covered*width_covered + 2*length_covered*height_drying + 2*width_covered*height_drying))
    concrete_volume_uncovered = number_uncovered*(concrete_thickness*(length_uncovered*width_uncovered + 2*length_uncovered*height_drying + 2*width_uncovered*height_drying))
    concrete_volume_storage = number_storage*(concrete_thickness*(length_covered*width_covered + 2*length_covered*height_storage + 2*width_covered*height_storage))
    concrete_volume = concrete_volume_covered + concrete_volume_uncovered + concrete_volume_storage
    concrete_emissions = concrete_IF_GHG*concrete_volume

    roof_area = np.full(np.shape(roof_slope), np.nan)
    for i in range(0,len(roof_slope)):
        roof_area[i] = ((number_covered[i]+number_storage[i])*length_covered[i]*width_covered[i])/math.cos(roof_slope[i]*math.pi/180)
    roof_mass_steel = roof_area*roof_mass
    column_mass = (covered_column_height*column_mass_per_meter)*(covered_columns_per_side*2)*(number_storage+number_covered)
    steel_emissions = steel_IF_GHG*(roof_mass_steel+column_mass)

    construction_emissions_annual = (concrete_emissions+steel_emissions)/plant_lifetime/existing_population

    tech_construction_emissions[:,3:4] = tech_construction_emissions[:,3:4] + construction_emissions_annual

    # cost
    if parameters.use_total_price.expected == 'no':
        concrete_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.concrete_cost, correlation_distributions, correlation_parameters, n_samples)
        steel_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.steel_cost, correlation_distributions, correlation_parameters, n_samples)
        cap_cost = (concrete_volume)*concrete_cost + (roof_mass_steel+column_mass)*steel_cost

        opex_percent, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.existing_plant_opex, correlation_distributions, correlation_parameters, n_samples)
        annual_opex = (cap_cost*(opex_percent/100)) / existing_population
        capex_annualized = (cap_cost * ((discount_rate*(1 + discount_rate)**plant_lifetime)/(((1 + discount_rate)**plant_lifetime) - 1))) / existing_population

        construction_cost[:,3:4] = construction_cost[:,3:4] + capex_annualized
        operating_cost[:,3:4] = operating_cost[:,3:4] + annual_opex

        electricity_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_cost, correlation_distributions, correlation_parameters, n_samples)
        electricity_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_GHG, correlation_distributions, correlation_parameters, n_samples)
        # assume all operating expenses come from electricity
        annual_kWh = annual_opex / electricity_cost
        operating_emissions = (annual_kWh * electricity_GHG)
        tech_operating_emissions[:,3:4] = tech_operating_emissions[:,3:4] + operating_emissions


    outputs[:,0:9] = np.concatenate((mass, mass_dry, N_total, P_total, K_total,
                                     Mg_total, Ca_total, energy, N_ammonia), 1)

    direct_emissions[:,3:4] = direct_emissions[:,3:4] + CH4eq + N2Oeq

    additional_storage = additional_storage + retention_time/365

    return outputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, additional_storage

#%% drying beds alternate function

def drying_beds_alt(inputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage, concrete_thickness, roof_slope, roof_mass, concrete_IF_GHG, steel_IF_GHG, plant_lifetime, existing_population, discount_rate):
    outputs = copy.deepcopy(inputs)

    mass = np.reshape(inputs[:,0], (-1,1))
    mass_dry = np.reshape(inputs[:,1], (-1,1))
    N_total = np.reshape(inputs[:,2], (-1,1))
    P_total = np.reshape(inputs[:,3], (-1,1))
    K_total = np.reshape(inputs[:,4], (-1,1))
    Mg_total = np.reshape(inputs[:,5], (-1,1))
    Ca_total = np.reshape(inputs[:,6], (-1,1))
    energy = np.reshape(inputs[:,7], (-1,1))
    N_ammonia = np.reshape(inputs[:,8], (-1,1))

    # parameters
    COD_degradation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_degradation_drying_bed_alt, correlation_distributions, correlation_parameters, n_samples)
    retention_time, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.retention_time_drying_bed_alt, correlation_distributions, correlation_parameters, n_samples)
    MCF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_drying_bed_alt, correlation_distributions, correlation_parameters, n_samples)
    N2O_EF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_drying_bed_alt, correlation_distributions, correlation_parameters, n_samples)
    N_denitrification, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_max_denitrification_drying_bed_alt, correlation_distributions, correlation_parameters, n_samples)
    final_solids_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.final_solids_content_drying_bed_alt, correlation_distributions, correlation_parameters, n_samples)

    # calculations
    COD_influent = energy/14/1000

    #!!! Changed algorithm
    t0 = previous_storage_time + additional_storage
    t = retention_time / 365
    COD_loss = first_order_decay(rate_constant, t0, t,
                                 max_decay=COD_degradation/100,
                                 tot=COD_influent)

    # COD_degradable = COD_influent * COD_degradation/100
    # COD_initial = (COD_degradable*rate_constant*previous_storage_time)/(np.exp(-rate_constant*additional_storage) - np.exp(-rate_constant*(additional_storage+previous_storage_time)))
    # COD_after = (COD_initial/(rate_constant*(previous_storage_time)))*(np.exp(-rate_constant*(additional_storage+retention_time/365)) - np.exp(-rate_constant*(previous_storage_time + additional_storage + retention_time/365)))
    # COD_loss = COD_degradable - COD_after

    COD_reduction = COD_loss/COD_influent
    CH4_emission = COD_loss*(MCF/100)*maximum_methane_emission
    CH4eq = CH4_emission*CH4_GWP
    COD_effluent = COD_influent - COD_loss
    energy = COD_effluent*14*1000
    mass = mass - mass_dry*COD_reduction
    mass_dry = mass_dry - mass_dry*COD_reduction

    N_degradable = N_total * (N_denitrification/100)
    N_initial = (N_degradable*rate_constant*previous_storage_time)/(np.exp(-rate_constant*additional_storage) - np.exp(-rate_constant*(additional_storage+previous_storage_time)))
    N_after = (N_initial/(rate_constant*(previous_storage_time)))*(np.exp(-rate_constant*(additional_storage+retention_time/365)) - np.exp(-rate_constant*(previous_storage_time + additional_storage + retention_time/365)))
    N_loss = N_degradable - N_after
    N2O_emission = N_total*(N_loss/N_degradable)*(N2O_EF/100)*(44/28)
    for i in range(0, len(N2O_emission)):
        if N2O_emission[i] > N_loss[i]*(44/28):
            N2O_emission[i] = N_loss[i]*(44/28)
    N2Oeq = N2O_emission*N2O_GWP
    N_total = N_total - N_loss
    N_ammonia = N_ammonia - N_loss
    for i in range(0, len(N_ammonia)):
        if N_ammonia[i] < 0:
            N_ammonia[i] = 0

    for i in range(0, len(mass_dry)):
        if final_solids_content[i] > ((mass_dry[i]/mass[i])*100):
            mass[i] = mass_dry[i] / (final_solids_content[i]/100)

    # construction costs and emissions
    number_unplanted, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.number_unplanted_beds_alt, correlation_distributions, correlation_parameters, n_samples)
    unplanted_height, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.unplanted_bed_wall_height_alt, correlation_distributions, correlation_parameters, n_samples)
    unplanted_width, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.unplanted_bed_width_alt, correlation_distributions, correlation_parameters, n_samples)
    unplanted_length, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.unplanted_bed_length_alt, correlation_distributions, correlation_parameters, n_samples)

    number_planted, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.number_planted_beds_alt, correlation_distributions, correlation_parameters, n_samples)
    planted_height, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.planted_bed_wall_height_alt, correlation_distributions, correlation_parameters, n_samples)
    planted_width, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.planted_bed_width_alt, correlation_distributions, correlation_parameters, n_samples)
    planted_length, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.planted_bed_length_alt, correlation_distributions, correlation_parameters, n_samples)

    concrete_volume_unplanted = number_unplanted*(concrete_thickness*(unplanted_length*unplanted_width + 2*unplanted_length*unplanted_height + 2*unplanted_width*unplanted_height))
    concrete_volume_planted = number_planted*(concrete_thickness*(planted_length*planted_width + 2*planted_length*planted_height + 2*planted_width*planted_height))
    concrete_volume = concrete_volume_unplanted + concrete_volume_planted
    concrete_emissions = concrete_IF_GHG*concrete_volume

    construction_emissions_annual = (concrete_emissions)/plant_lifetime/existing_population

    tech_construction_emissions[:,3:4] = tech_construction_emissions[:,3:4] + construction_emissions_annual

    # cost
    if parameters.use_total_price.expected == 'no':
        concrete_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.concrete_cost, correlation_distributions, correlation_parameters, n_samples)
        cap_cost = (concrete_volume)*concrete_cost

        opex_percent, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.existing_plant_opex, correlation_distributions, correlation_parameters, n_samples)
        annual_opex = (cap_cost*(opex_percent/100)) / existing_population
        capex_annualized = (cap_cost * ((discount_rate*(1 + discount_rate)**plant_lifetime)/(((1 + discount_rate)**plant_lifetime) - 1))) / existing_population

        construction_cost[:,3:4] = construction_cost[:,3:4] + capex_annualized
        operating_cost[:,3:4] = operating_cost[:,3:4] + annual_opex

        electricity_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_cost, correlation_distributions, correlation_parameters, n_samples)
        electricity_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_GHG, correlation_distributions, correlation_parameters, n_samples)
        # assume all operating expenses come from electricity
        annual_kWh = annual_opex / electricity_cost
        operating_emissions = (annual_kWh * electricity_GHG)
        tech_operating_emissions[:,3:4] = tech_operating_emissions[:,3:4] + operating_emissions


    outputs[:,0:9] = np.concatenate((mass, mass_dry, N_total, P_total, K_total,
                                     Mg_total, Ca_total, energy, N_ammonia), 1)

    direct_emissions[:,3:4] = direct_emissions[:,3:4] + CH4eq + N2Oeq

    additional_storage = additional_storage + retention_time/365

    return outputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, additional_storage

#%% anaerobic baffled reactor function

def ABR(inputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, biogas, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, CH4_energy, previous_storage_time, additional_storage, sludge_flow_alt, concrete_thickness, concrete_IF_GHG, gravel_IF_GHG, gravel_bulk_density, excavation_IF_GHG, plant_lifetime, sludge_pop_alt, discount_rate):
    outputs = copy.deepcopy(inputs)

    mass = np.reshape(inputs[:,0], (-1,1))
    mass_dry = np.reshape(inputs[:,1], (-1,1))
    N_total = np.reshape(inputs[:,2], (-1,1))
    P_total = np.reshape(inputs[:,3], (-1,1))
    K_total = np.reshape(inputs[:,4], (-1,1))
    Mg_total = np.reshape(inputs[:,5], (-1,1))
    Ca_total = np.reshape(inputs[:,6], (-1,1))
    energy = np.reshape(inputs[:,7], (-1,1))
    N_ammonia = np.reshape(inputs[:,8], (-1,1))

    # biogas production
    # parameters
    MCF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_ABR, correlation_distributions, correlation_parameters, n_samples)
    COD_removal, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_removal_ABR, correlation_distributions, correlation_parameters, n_samples)
    residence_time, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.retention_time_ABR, correlation_distributions, correlation_parameters, n_samples)

    # calculations
    # carbon (kg COD/yr; assume 14 kJ/g COD in wastewater)
    COD_total = (energy/14/1000)
    COD_degrade = COD_total*(COD_removal/100)
    CH4_production = COD_degrade*(MCF/100)*maximum_methane_emission
    COD_total = COD_total - COD_degrade
    energy = COD_total*14*1000
    if parameters.CH4_captured_ABR.expected == 'yes':
        CH4eq = np.full([n_samples, 1], 0)
    else:
        CH4eq = CH4_production*CH4_GWP
        CH4_production = 0

    # nitrogen emissions (kg N/yr; N2O expressed as kg N2O/yr)
    N_removal, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_removal_ABR, correlation_distributions, correlation_parameters, n_samples)
    N_removed = N_total * (N_removal/100)
    N_total = N_total - N_removed
    N_ammonia = N_ammonia- N_removed
    for i in range(0, len(N_ammonia)):
        if N_ammonia[i] < 0:
            N_ammonia[i] = 0

    if parameters.N_emission_in_biogas_ABR.expected == 'yes':
        N2O_EF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_ABR, correlation_distributions, correlation_parameters, n_samples)
        N_denitrification, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_max_denitrification_ABR, correlation_distributions, correlation_parameters, n_samples)
        N_degradable = N_removed * (N_denitrification/100)
        N2O_emission = N_removed*(N2O_EF/100)*(44/28)
        for i in range(0, len(N2O_emission)):
            if N2O_emission[i] > N_degradable[i]*(44/28):
                N2O_emission[i] = N_degradable[i]*(44/28)
        N2Oeq = N2O_emission*N2O_GWP
    else:
        N2Oeq = np.full([n_samples, 1], 0)


    # solids loss (based on COD loss)
    mass = mass - mass_dry * (COD_removal/100)
    mass_dry = mass_dry - mass_dry * (COD_removal/100)

    # energy in biogas (kJ/yr)
    biogas_new = CH4_production / 16 * CH4_energy * 1000
    biogas = biogas + biogas_new

    outputs[:,0:9] = np.concatenate((mass, mass_dry, N_total, P_total, K_total,
                                     Mg_total, Ca_total, energy, N_ammonia), 1)

    direct_emissions[:,3:4] = direct_emissions[:,3:4] + CH4eq + N2Oeq

    # construction costs and emissions
    number, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.number_ABR, correlation_distributions, correlation_parameters, n_samples)
    length, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.length_ABR, correlation_distributions, correlation_parameters, n_samples)
    width, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.width_ABR, correlation_distributions, correlation_parameters, n_samples)
    height, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.height_ABR, correlation_distributions, correlation_parameters, n_samples)
    baffles, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.baffles_ABR, correlation_distributions, correlation_parameters, n_samples)
    add_concrete, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.additional_concrete_ABR, correlation_distributions, correlation_parameters, n_samples)

    volume = length*width*height

    #!!! Original code didn't divide add_concrete by 100
    concrete_volume = number*concrete_thickness*(2*length*width+2*length*height+(2+baffles)*width*height)*(add_concrete/100)
    # concrete_volume = number*concrete_thickness*(2*length*width+2*length*height+(2+baffles)*width*height)*add_concrete

    concrete_emissions = concrete_volume * concrete_IF_GHG

    gravel_volume = number*length*width*height/(baffles+1)
    gravel_mass = gravel_volume * gravel_bulk_density
    gravel_emissions = gravel_mass * gravel_IF_GHG

    excavation_volume = number*volume
    excavation_emissions = excavation_volume * excavation_IF_GHG

    construction_emissions_annual = (concrete_emissions+gravel_emissions+excavation_emissions)/plant_lifetime/sludge_pop_alt

    tech_construction_emissions[:,3:4] = tech_construction_emissions[:,3:4] + construction_emissions_annual

    # cost
    if parameters.use_total_price.expected == 'no':
        concrete_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.concrete_cost, correlation_distributions, correlation_parameters, n_samples)
        steel_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.steel_cost, correlation_distributions, correlation_parameters, n_samples)
        cap_cost = (concrete_volume)*concrete_cost

        opex_percent, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.alternative_plant_opex, correlation_distributions, correlation_parameters, n_samples)
        annual_opex = (cap_cost*(opex_percent/100)) / sludge_pop_alt
        capex_annualized = (cap_cost * ((discount_rate*(1 + discount_rate)**plant_lifetime)/(((1 + discount_rate)**plant_lifetime) - 1))) / sludge_pop_alt

        construction_cost[:,3:4] = construction_cost[:,3:4] + capex_annualized
        operating_cost[:,3:4] = operating_cost[:,3:4] + annual_opex

        electricity_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_cost, correlation_distributions, correlation_parameters, n_samples)
        electricity_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_GHG, correlation_distributions, correlation_parameters, n_samples)
        # assume all operating expenses come from electricity
        annual_kWh = annual_opex / electricity_cost
        operating_emissions = (annual_kWh * electricity_GHG)
        tech_operating_emissions[:,3:4] = tech_operating_emissions[:,3:4] + operating_emissions


    additional_storage = additional_storage + residence_time/365

    return outputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, biogas, correlation_distributions, correlation_parameters, additional_storage

#%% secondary liquid bed function
def secondary_liquid_bed(inputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage, concrete_thickness, roof_slope, roof_mass, concrete_IF_GHG, steel_IF_GHG, plant_lifetime, existing_population, discount_rate):
    outputs = copy.deepcopy(inputs)

    mass = np.reshape(inputs[:,0], (-1,1))
    mass_dry = np.reshape(inputs[:,1], (-1,1))
    N_total = np.reshape(inputs[:,2], (-1,1))
    P_total = np.reshape(inputs[:,3], (-1,1))
    K_total = np.reshape(inputs[:,4], (-1,1))
    Mg_total = np.reshape(inputs[:,5], (-1,1))
    Ca_total = np.reshape(inputs[:,6], (-1,1))
    energy = np.reshape(inputs[:,7], (-1,1))
    N_ammonia = np.reshape(inputs[:,8], (-1,1))

    # parameters
    COD_degradation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.COD_degradation_sec, correlation_distributions, correlation_parameters, n_samples)
    retention_time, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.retention_time_sec, correlation_distributions, correlation_parameters, n_samples)
    MCF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_sec, correlation_distributions, correlation_parameters, n_samples)

    # calculations
    COD_influent = energy/14/1000

    #!!! Changed algorithm
    t0 = previous_storage_time + additional_storage
    t = retention_time / 365
    COD_loss = first_order_decay(rate_constant, t0, t,
                                 max_decay=COD_degradation/100,
                                 tot=COD_influent)

    # COD_degradable = COD_influent * COD_degradation/100
    # COD_initial = (COD_degradable*rate_constant*previous_storage_time)/(np.exp(-rate_constant*additional_storage) - np.exp(-rate_constant*(additional_storage+previous_storage_time)))
    # COD_after = (COD_initial/(rate_constant*(previous_storage_time)))*(np.exp(-rate_constant*(additional_storage+retention_time/365)) - np.exp(-rate_constant*(previous_storage_time + additional_storage + retention_time/365)))
    # COD_loss = COD_degradable - COD_after

    COD_reduction = COD_loss/COD_influent
    CH4_emission = COD_loss*(MCF/100)*maximum_methane_emission
    CH4eq = CH4_emission*CH4_GWP
    COD_effluent = COD_influent - COD_loss
    energy = COD_effluent*14*1000
    mass = mass - mass_dry*COD_reduction
    mass_dry = mass_dry - mass_dry*COD_reduction

    if parameters.N_emission_from_anaerobic_lagoon.expected == 'yes':
        N2O_EF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_sec, correlation_distributions, correlation_parameters, n_samples)
        N_denitrification, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_max_denitrification_sec, correlation_distributions, correlation_parameters, n_samples)

        #!!! Changed algorithm
        t0 = previous_storage_time + additional_storage
        t = retention_time / 365
        N_loss = first_order_decay(rate_constant, t0, t,
                                   max_decay=N_denitrification/100,
                                   tot=N_total)
        N2O_emission = N_loss*(N2O_EF/100)*(44/28)

        # N_degradable = N_total * (N_denitrification/100)
        # N_initial = (N_degradable*rate_constant*previous_storage_time)/(np.exp(-rate_constant*additional_storage) - np.exp(-rate_constant*(additional_storage+previous_storage_time)))
        # N_after = (N_initial/(rate_constant*(previous_storage_time)))*(np.exp(-rate_constant*(additional_storage+retention_time/365)) - np.exp(-rate_constant*(previous_storage_time + additional_storage + retention_time/365)))
        # N_loss = N_degradable - N_after
        # N2O_emission = N_total*(N_loss/N_degradable)*(N2O_EF/100)*(44/28)

        for i in range(0, len(N2O_emission)):
            if N2O_emission[i] > N_loss[i]*(44/28):
                N2O_emission[i] = N_loss[i]*(44/28)
        N2Oeq = N2O_emission*N2O_GWP
        N_total = N_total - N_loss
        N_ammonia = N_ammonia - N_loss
        for i in range(0, len(N_ammonia)):
            if N_ammonia[i] < 0:
                N_ammonia[i] = 0
    else:
        N2Oeq = np.full([n_samples, 1], 0)

    # construction costs and emissions
    number, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.number_sec, correlation_distributions, correlation_parameters, n_samples)
    height, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.height_sec, correlation_distributions, correlation_parameters, n_samples)
    width, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.width_sec, correlation_distributions, correlation_parameters, n_samples)
    length, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.length_sec, correlation_distributions, correlation_parameters, n_samples)

    concrete_volume = number*(concrete_thickness*(length*width + 2*length*height + 2*width*height))
    concrete_emissions = concrete_IF_GHG*concrete_volume

    construction_emissions_annual = (concrete_emissions)/plant_lifetime/existing_population

    tech_construction_emissions[:,3:4] = tech_construction_emissions[:,3:4] + construction_emissions_annual

    # cost
    if parameters.use_total_price.expected == 'no':
        concrete_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.concrete_cost, correlation_distributions, correlation_parameters, n_samples)
        cap_cost = (concrete_volume)*concrete_cost

        opex_percent, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.existing_plant_opex, correlation_distributions, correlation_parameters, n_samples)
        annual_opex = (cap_cost*(opex_percent/100)) / existing_population
        capex_annualized = (cap_cost * ((discount_rate*(1 + discount_rate)**plant_lifetime)/(((1 + discount_rate)**plant_lifetime) - 1))) / existing_population

        construction_cost[:,3:4] = construction_cost[:,3:4] + capex_annualized
        operating_cost[:,3:4] = operating_cost[:,3:4] + annual_opex

        electricity_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_cost, correlation_distributions, correlation_parameters, n_samples)
        electricity_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_GHG, correlation_distributions, correlation_parameters, n_samples)
        # assume all operating expenses come from electricity
        annual_kWh = annual_opex / electricity_cost
        operating_emissions = (annual_kWh * electricity_GHG)
        tech_operating_emissions[:,3:4] = tech_operating_emissions[:,3:4] + operating_emissions


    outputs[:,0:9] = np.concatenate((mass, mass_dry, N_total, P_total, K_total,
                                     Mg_total, Ca_total, energy, N_ammonia), 1)

    direct_emissions[:,3:4] = direct_emissions[:,3:4] + CH4eq + N2Oeq

    additional_storage = additional_storage + retention_time/365

    return outputs, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, additional_storage


#%% main function

def main(input_excel_name, excreta_inputs, liquid_inputs, solid_inputs, direct_emissions, tech_construction_emissions, tech_operating_emissions, construction_cost, operating_cost, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, discount_rate, exchange_rate):
    # import module parameters from input spreadsheet
    parameters = pd.DataFrame.transpose(pd.read_excel(input_excel_name, sheet_name = 'treatment').set_index('parameters'))

    CH4_energy, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.CH4_energy, correlation_distributions, correlation_parameters, n_samples)
    sewer_flow_rate, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sewer_flow_existing, correlation_distributions, correlation_parameters, n_samples)
    sludge_flow_rate, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sludge_flow_existing, correlation_distributions, correlation_parameters, n_samples)
    sludge_flow_alt, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sludge_flow_alternative, correlation_distributions, correlation_parameters, n_samples)
    concrete_thickness, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.concrete_thickness, correlation_distributions, correlation_parameters, n_samples)
    roof_slope, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.roof_slope, correlation_distributions, correlation_parameters, n_samples)
    roof_mass, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.roof_mass, correlation_distributions, correlation_parameters, n_samples)
    concrete_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.concrete_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    stainless_steel_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.stainless_steel_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    stainless_steel_sheet_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.stainless_steel_sheet_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    excavation_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.excavation_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    liner_mass, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.liner_mass, correlation_distributions, correlation_parameters, n_samples)
    liner_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.liner_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    gravel_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.gravel_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    gravel_bulk_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.gravel_bulk_density, correlation_distributions, correlation_parameters, n_samples)
    existing_plant_lifetime, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.existing_plant_lifetime, correlation_distributions, correlation_parameters, n_samples)
    alternative_plant_lifetime, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.alternative_plant_lifetime, correlation_distributions, correlation_parameters, n_samples)
    existing_sewer_population, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.existing_sewer_population_served, correlation_distributions, correlation_parameters, n_samples)
    existing_sludge_population, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.existing_sludge_population_served, correlation_distributions, correlation_parameters, n_samples)
    sludge_pop_alt, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sludge_population_alternative, correlation_distributions, correlation_parameters, n_samples)

    existing_population = existing_sewer_population + existing_sludge_population

    steel_IF_GHG = stainless_steel_IF_GHG + stainless_steel_sheet_IF_GHG
    flow_rate = sewer_flow_rate + sludge_flow_rate

    additional_storage_excreta = 0
    additional_storage_liquid = 0
    additional_storage_solid = 0

    # define the module(s)
    excreta_module = [parameters.mixed_excreta_module_1.expected, parameters.mixed_excreta_module_2.expected, parameters.mixed_excreta_module_3.expected]
    liquid_module = [parameters.liquid_module_1.expected, parameters.liquid_module_2.expected, parameters.liquid_module_3.expected]
    solid_module = [parameters.solid_module_1.expected, parameters.solid_module_2.expected, parameters.solid_module_3.expected]

    # create temporary variables to track excreta, solids, and liquids
    excreta_temp = copy.deepcopy(excreta_inputs)
    liquid_temp = copy.deepcopy(liquid_inputs)
    solid_temp = copy.deepcopy(solid_inputs)

    biogas = np.full([len(excreta_temp), 1], 0)

    for i in range(0, len(excreta_module)):
        if (type(excreta_module[i]) is float) and (type(liquid_module[i]) is float) and (type(solid_module[i]) is float):
            # other numerical inputs are not valid
            if (not np.isnan(excreta_module[i])):
                raise ValueError('The specified excreta treatment module is not valid.')
            if (not np.isnan(liquid_module[i])):
                raise ValueError('The specified liquid treatment module is not valid.')
            if (not np.isnan(solid_module[i])):
                raise ValueError('The specified solid treatment module is not valid.')

        # otherwise, are both mixed and split stream options entered?
        elif (type(excreta_module[i]) is str) and ((type(liquid_module[i]) is str) or (type(solid_module[i]) is str)):
            raise ValueError('Modules for both the mixed and separated cases should not be evaluated simultaneously.')

        # check mixed stream options first
        if type(excreta_module[i]) is str:

            # anaerobic baffled reactor module
            if excreta_module[i] == 'ABR':
                (excreta_temp, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, biogas, correlation_distributions,
                 correlation_parameters, additional_storage_liquid) = ABR(excreta_temp, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, biogas, parameters,
                                       correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, CH4_energy, previous_storage_time, additional_storage_liquid, sludge_flow_alt, concrete_thickness, concrete_IF_GHG, gravel_IF_GHG, gravel_bulk_density, excavation_IF_GHG, alternative_plant_lifetime, sludge_pop_alt, discount_rate)

            # sedimentation module (mixed input, separate outputs)
            elif excreta_module[i] == 'sedimentation':
                (liquid_temp, solid_temp, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions,
                 correlation_parameters, additional_storage_liquid, additional_storage_solid) = sedimentation(excreta_temp, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters,
                                       correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage_excreta, concrete_thickness, roof_slope, roof_mass, concrete_IF_GHG, steel_IF_GHG, existing_plant_lifetime, existing_sludge_population, discount_rate)

                excreta_temp = np.full(np.shape(excreta_temp), np.nan)

            # sludge separator module (mixed input, separate outputs)
            elif excreta_module[i] == 'sludge_separator':
                (liquid_temp, solid_temp, direct_emissions, correlation_distributions,
                 correlation_parameters, additional_storage_liquid, additional_storage_solid) = sludge_separator(excreta_temp, direct_emissions, parameters,
                                       correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage_excreta)

                excreta_temp = np.full(np.shape(excreta_temp), np.nan)

            # if the excreta module input is not supported/valid
            else:
                raise ValueError('The treatment module specified for excreta is not valid.')

        # check liquid stream
        if (type(liquid_module[i]) is str):

            # anaerobic lagoon module
            if liquid_module[i] == 'anaerobic_lagoon':
                (liquid_temp, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions,
                 correlation_parameters, additional_storage_liquid) = anaerobic_lagoon(liquid_temp, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters,
                                       correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage_liquid, flow_rate, excavation_IF_GHG, liner_mass, liner_IF_GHG, existing_plant_lifetime, existing_population)

            elif liquid_module[i] == 'facultative_lagoon':
                (liquid_temp, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions,
                 correlation_parameters, additional_storage_liquid) = facultative_lagoon(liquid_temp, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters,
                                       correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage_liquid, flow_rate, excavation_IF_GHG, liner_mass, liner_IF_GHG, existing_plant_lifetime, existing_population)

            elif liquid_module[i] == 'secondary_liquid_bed':
                (liquid_temp, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions,
                 correlation_parameters, additional_storage_solid) = secondary_liquid_bed(liquid_temp, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters,
                                       correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage_solid, concrete_thickness, roof_slope, roof_mass, concrete_IF_GHG, steel_IF_GHG, existing_plant_lifetime, existing_sludge_population, discount_rate)
            # if the liquid module input is not supported/valid
            else:
               raise ValueError('The treatment module specified for liquid is not valid.')

        # check solid stream
        if (type(solid_module[i]) is str):

            # unplanted drying bed module
            if solid_module[i] == 'unplanted_drying_bed':
                (solid_temp, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions,
                 correlation_parameters, additional_storage_solid) = unplanted_drying_bed(solid_temp, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters,
                                       correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage_solid, concrete_thickness, roof_slope, roof_mass, concrete_IF_GHG, steel_IF_GHG, existing_plant_lifetime, existing_sludge_population, discount_rate)

            # drying beds alternate module
            elif solid_module[i] == 'drying_beds_alt':
                (solid_temp, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, correlation_distributions,
                 correlation_parameters, additional_storage_solid) = drying_beds_alt(solid_temp, construction_cost, operating_cost, direct_emissions, tech_construction_emissions, tech_operating_emissions, parameters,
                                       correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, additional_storage_solid, concrete_thickness, roof_slope, roof_mass, concrete_IF_GHG, steel_IF_GHG, existing_plant_lifetime, existing_sludge_population, discount_rate)

            # if the solid module input is not supported/valid
            else:
                raise ValueError('The treatment module specified for solid is not valid.')

    # after iteration, set outputs equal to current values of temporary variables
    excreta_outputs = excreta_temp
    liquid_outputs = liquid_temp
    solid_outputs = solid_temp

    electricity_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_cost, correlation_distributions, correlation_parameters, n_samples)
    electricity_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.electricity_GHG, correlation_distributions, correlation_parameters, n_samples)

    if parameters.use_total_price.expected == 'yes':
        if parameters.use_existing_plant.expected == 'yes':
            cap_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.existing_plant_capex, correlation_distributions, correlation_parameters, n_samples)
            annual_electricity, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.existing_plant_electricity, correlation_distributions, correlation_parameters, n_samples)
            staff, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.existing_plant_staff, correlation_distributions, correlation_parameters, n_samples)
            monthly_salary, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.existing_plant_salary, correlation_distributions, correlation_parameters, n_samples)
            annual_labor = staff * (monthly_salary/exchange_rate) * 12
            annual_opex = (annual_electricity*electricity_cost + annual_labor) / existing_population
            capex_annualized = (cap_cost * ((discount_rate*(1 + discount_rate)**existing_plant_lifetime)/(((1 + discount_rate)**existing_plant_lifetime) - 1))) / existing_population
            operating_emissions = (annual_electricity * electricity_GHG) / existing_population
        elif parameters.use_existing_plant.expected == 'no':
            cap_cost, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.alternative_plant_capex, correlation_distributions, correlation_parameters, n_samples)
            annual_electricity, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.alternative_plant_electricity, correlation_distributions, correlation_parameters, n_samples)
            skilled_staff, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.alternative_plant_skilled_staff, correlation_distributions, correlation_parameters, n_samples)
            unskilled_staff, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.alternative_plant_unskilled_staff, correlation_distributions, correlation_parameters, n_samples)
            skilled_monthly_salary, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.alternative_plant_skilled_salary, correlation_distributions, correlation_parameters, n_samples)
            unskilled_monthly_salary, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.alternative_plant_unskilled_salary, correlation_distributions, correlation_parameters, n_samples)
            annual_labor = (skilled_staff*(skilled_monthly_salary/exchange_rate) + unskilled_staff*(unskilled_monthly_salary/exchange_rate)) * 12
            annual_opex = (annual_electricity*electricity_cost + annual_labor) / sludge_pop_alt
            capex_annualized = (cap_cost * ((discount_rate*(1 + discount_rate)**alternative_plant_lifetime)/(((1 + discount_rate)**alternative_plant_lifetime) - 1))) / sludge_pop_alt
            operating_emissions = (annual_electricity * electricity_GHG) / sludge_pop_alt
        construction_cost[:,3:4] = construction_cost[:,3:4] + capex_annualized
        operating_cost[:,3:4] = operating_cost[:,3:4] + annual_opex
        # assume all operating expenses come from electricity
        tech_operating_emissions[:,3:4] = tech_operating_emissions[:,3:4] + operating_emissions

    return excreta_outputs, liquid_outputs, solid_outputs, direct_emissions, tech_construction_emissions, tech_operating_emissions, construction_cost, operating_cost, biogas, correlation_distributions, correlation_parameters