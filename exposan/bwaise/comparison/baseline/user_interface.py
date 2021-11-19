#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 14:00:56 2017

@author: John Trimmer; Yalin Li (minor modification)
"""

# These functions track resources, costs, and emissions through various user interface options
#  (i) dry toilet (pit latrine receiving mixed excreta)
#  (ii) UDDT (urine-diverting dry toilet that separates urine and feces)

import numpy as np
import pandas as pd
import lhs



#%% Dry Toilet function
def dry_toilet(urine_inputs, feces_inputs, parameters, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, n_samples, discount_rate):
    # separate out input variables from concatenated arrays
    urine_input = np.reshape(urine_inputs[:,0], (-1,1))
    urine_dry_input = np.reshape(urine_inputs[:,1], (-1,1))
    N_urine_input = np.reshape(urine_inputs[:,2], (-1,1))
    P_urine_input = np.reshape(urine_inputs[:,3], (-1,1))
    K_urine_input = np.reshape(urine_inputs[:,4], (-1,1))
    Mg_urine_input = np.reshape(urine_inputs[:,5], (-1,1))
    Ca_urine_input = np.reshape(urine_inputs[:,6], (-1,1))
    energy_urine_input = np.reshape(urine_inputs[:,7], (-1,1))
    N_ammonia_urine_input = np.reshape(urine_inputs[:,8], (-1,1))

    feces_input = np.reshape(feces_inputs[:,0], (-1,1))
    feces_dry_input = np.reshape(feces_inputs[:,1], (-1,1))
    N_feces_input = np.reshape(feces_inputs[:,2], (-1,1))
    P_feces_input = np.reshape(feces_inputs[:,3], (-1,1))
    K_feces_input = np.reshape(feces_inputs[:,4], (-1,1))
    Mg_feces_input = np.reshape(feces_inputs[:,5], (-1,1))
    Ca_feces_input = np.reshape(feces_inputs[:,6], (-1,1))
    energy_feces_input = np.reshape(feces_inputs[:,7], (-1,1))
    N_ammonia_feces_input = np.reshape(feces_inputs[:,8], (-1,1))

    # generate distributions of parameters
    # cleansing material: toilet paper
    if parameters.toilet_paper.expected == 'yes':
        toilet_paper_addition, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_addition, correlation_distributions, correlation_parameters, n_samples)
        toilet_paper_COD_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_COD_content, correlation_distributions, correlation_parameters, n_samples)
        toilet_paper_N_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_N_content, correlation_distributions, correlation_parameters, n_samples)
        toilet_paper_P_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_P_content, correlation_distributions, correlation_parameters, n_samples)
        toilet_paper_K_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_K_content, correlation_distributions, correlation_parameters, n_samples)
        toilet_paper_TS_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_TS_content, correlation_distributions, correlation_parameters, n_samples)
    elif parameters.toilet_paper.expected == 'no':
        toilet_paper_addition = 0
        toilet_paper_COD_content = 0
        toilet_paper_N_content = 0
        toilet_paper_P_content = 0
        toilet_paper_K_content = 0
        toilet_paper_TS_content = 0
    # flushing water
    if parameters.flushing_water.expected == 'yes':
        flushing_water_addition, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.flushing_water_use, correlation_distributions, correlation_parameters, n_samples)
    elif parameters.flushing_water.expected == 'no':
        flushing_water_addition = 0
    # cleansing material: water
    if (parameters.cleansing_water.expected == 'yes') and (parameters.cleansing_water_hole.expected == 'no'):
        cleansing_water_addition, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.cleansing_water_use, correlation_distributions, correlation_parameters, n_samples)
    elif parameters.cleansing_water.expected == 'no':
        cleansing_water_addition = 0
    # desiccant
    if parameters.desiccant.expected == 'yes':
        desiccant_volume, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_volume, correlation_distributions, correlation_parameters, n_samples)
        desiccant_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_density, correlation_distributions, correlation_parameters, n_samples)
        desiccant_C_N_ratio, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_C_N_ratio, correlation_distributions, correlation_parameters, n_samples)
        desiccant_N_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_N_content, correlation_distributions, correlation_parameters, n_samples)
        desiccant_P_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_P_content, correlation_distributions, correlation_parameters, n_samples)
        desiccant_K_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_K_content, correlation_distributions, correlation_parameters, n_samples)
        desiccant_Mg_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_Mg_content, correlation_distributions, correlation_parameters, n_samples)
        desiccant_Ca_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_Ca_content, correlation_distributions, correlation_parameters, n_samples)
    elif parameters.desiccant.expected == 'no':
        desiccant_volume = 0
        desiccant_density = 0
        desiccant_N_content = 0
        desiccant_P_content = 0
        desiccant_K_content = 0
        desiccant_Mg_content = 0
        desiccant_Ca_content = 0

    # household and toilet characteristics
    household_size, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.household_size, correlation_distributions, correlation_parameters, n_samples)
    household_use_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.household_use_density, correlation_distributions, correlation_parameters, n_samples)

    # calculate resources in materials added by user (per capita per year)
    # N (kg N/cap/yr)
    additives_N = ((toilet_paper_addition * toilet_paper_N_content * 365 / 1000000)
                    + (desiccant_volume * desiccant_density * (desiccant_N_content / 100) * 365 / 1000000))
    # P (kg P/cap/yr)
    additives_P = ((toilet_paper_addition * toilet_paper_P_content * 365 / 1000000)
                    + (desiccant_volume * desiccant_density * (desiccant_P_content / 100) * 365 / 1000000))
    # K (kg K/cap/yr)
    additives_K = ((toilet_paper_addition * toilet_paper_K_content * 365 / 1000000)
                    + (desiccant_volume * desiccant_density * (desiccant_K_content / 100) * 365 / 1000000))
    # Mg (kg Mg/cap/yr)
    additives_Mg = (desiccant_volume * desiccant_density * (desiccant_Mg_content / 100) * 365 / 1000000)
    # Ca (kg Ca/cap/yr)
    additives_Ca = (desiccant_volume * desiccant_density * (desiccant_Ca_content / 100) * 365 / 1000000)
    # energy (kJ/cap/yr) - assume desiccant contains no additional energy (wood ash), so only include TP
    #   - also assume 3.86 kWh energy production per kg COD oxidized to CO2 and H20
    additives_energy = toilet_paper_addition * toilet_paper_COD_content * 365 / 1000000 * 3.86 * 3600
    # total per capita addition (kg/cap/yr)
    additives_input = ((toilet_paper_addition * toilet_paper_TS_content * 365 / 1000000)
                        + (desiccant_volume * desiccant_density * 365 / 1000000))
    # total dry addition (kg/cap/yr) - assume toilet paper and desiccant have moisture content ~ 0%
    additives_dry_input = ((toilet_paper_addition * toilet_paper_TS_content * 365 / 1000000)
                            + (desiccant_volume * desiccant_density * 365 / 1000000))
    # water input
    water_input = (flushing_water_addition + cleansing_water_addition) * 365

    # mix materials together to form single stream from all HHs using toilet - total input into toilet per year
    # nitrogen (kg N/cap/yr)
    N_excreta = (N_urine_input + N_feces_input + additives_N)
    # phosphorus (kg P/cap/yr)
    P_excreta = (P_urine_input + P_feces_input + additives_P)
    # potassium (kg K/cap/yr)
    K_excreta = (K_urine_input + K_feces_input + additives_K)
    # magnesium (kg Mg/cap/yr)
    Mg_excreta = (Mg_urine_input + Mg_feces_input + additives_Mg)
    # calcium (kg Ca/cap/yr)
    Ca_excreta = (Ca_urine_input + Ca_feces_input + additives_Ca)
    # energy (kJ/cap/yr)
    energy_excreta = (energy_urine_input + energy_feces_input + additives_energy)
    # total excreta input (kg/cap/yr)
    excreta = (urine_input + feces_input + additives_input + water_input)
    # total dry excreta input (kg/cap/yr)
    excreta_dry = (urine_dry_input + feces_dry_input + additives_dry_input)
    # ammonia (kg N/cap/yr)
    N_ammonia_excreta = (N_ammonia_urine_input + N_ammonia_feces_input)

    # concatenate excreta outputs to simplify function output
    excreta_outputs = np.concatenate((excreta, excreta_dry, N_excreta, P_excreta, K_excreta,
                                      Mg_excreta, Ca_excreta, energy_excreta, N_ammonia_excreta), 1)

    number_users = household_use_density * household_size

    # cost for one facility
    capex, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.dry_toilet_capex, correlation_distributions, correlation_parameters, n_samples)
    opex_percent, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.dry_toilet_opex, correlation_distributions, correlation_parameters, n_samples)
    lifetime, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_lifetime, correlation_distributions, correlation_parameters, n_samples)
    annual_opex = (capex*(opex_percent/100)) / number_users
    capex_annualized = (capex * ((discount_rate*(1 + discount_rate)**lifetime)/(((1 + discount_rate)**lifetime) - 1))) / number_users
    construction_cost[:,:1] = construction_cost[:,:1] + capex_annualized
    operating_cost[:,:1] = operating_cost[:,:1] + annual_opex

    # GHG emissions from construction
    cement, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.cement_dry_toilet, correlation_distributions, correlation_parameters, n_samples)
    gravel, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.gravel_dry_toilet, correlation_distributions, correlation_parameters, n_samples)
    sand, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sand_dry_toilet, correlation_distributions, correlation_parameters, n_samples)
    bricks, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.bricks_dry_toilet, correlation_distributions, correlation_parameters, n_samples)
    plastic, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.plastic_dry_toilet, correlation_distributions, correlation_parameters, n_samples)
    steel, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.steel_dry_toilet, correlation_distributions, correlation_parameters, n_samples)
    excavation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.excavation_dry_toilet, correlation_distributions, correlation_parameters, n_samples)
    wood, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.wood_dry_toilet, correlation_distributions, correlation_parameters, n_samples)

    plastic_mass, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.plastic_mass, correlation_distributions, correlation_parameters, n_samples)
    brick_volume, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.brick_volume, correlation_distributions, correlation_parameters, n_samples)
    brick_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.brick_density, correlation_distributions, correlation_parameters, n_samples)
    gravel_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.gravel_bulk_density, correlation_distributions, correlation_parameters, n_samples)
    sand_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sand_bulk_density, correlation_distributions, correlation_parameters, n_samples)
    steel_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.steel_density, correlation_distributions, correlation_parameters, n_samples)

    steel_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.steel_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    stainless_steel_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.stainless_steel_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    excavation_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.excavation_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    plastic_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.plastic_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    gravel_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.gravel_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    sand_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sand_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    cement_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.cement_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    bricks_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.bricks_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    wood_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.wood_IF_GHG, correlation_distributions, correlation_parameters, n_samples)

    cement_emissions = cement * cement_IF_GHG
    gravel_emissions = gravel*gravel_density*gravel_IF_GHG
    sand_emissions = sand*sand_density*sand_IF_GHG
    brick_emissions = bricks*brick_volume*brick_density*bricks_IF_GHG
    plastic_emissions = plastic*plastic_mass*plastic_IF_GHG
    steel_emissions = steel*steel_density*steel_IF_GHG
    excavation_emissions = excavation * excavation_IF_GHG
    wood_emissions = wood * wood_IF_GHG

    construction_emissions_annual = ((cement_emissions+gravel_emissions+sand_emissions+brick_emissions+plastic_emissions+steel_emissions
                                     +excavation_emissions+wood_emissions)/lifetime/number_users)

    tech_construction_emissions[:,:1] = tech_construction_emissions[:,:1] + construction_emissions_annual


    return excreta_outputs, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, number_users

#%% Urine-Diverting Dry Toilet function
def UDDT(urine_inputs, feces_inputs, parameters, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, n_samples, discount_rate):
    # separate out input variables from concatenated arrays
    urine_input = np.reshape(urine_inputs[:,0], (-1,1))
    urine_dry_input = np.reshape(urine_inputs[:,1], (-1,1))
    N_urine_input = np.reshape(urine_inputs[:,2], (-1,1))
    P_urine_input = np.reshape(urine_inputs[:,3], (-1,1))
    K_urine_input = np.reshape(urine_inputs[:,4], (-1,1))
    Mg_urine_input = np.reshape(urine_inputs[:,5], (-1,1))
    Ca_urine_input = np.reshape(urine_inputs[:,6], (-1,1))
    energy_urine_input = np.reshape(urine_inputs[:,7], (-1,1))
    N_ammonia_urine_input = np.reshape(urine_inputs[:,8], (-1,1))

    feces_input = np.reshape(feces_inputs[:,0], (-1,1))
    feces_dry_input = np.reshape(feces_inputs[:,1], (-1,1))
    N_feces_input = np.reshape(feces_inputs[:,2], (-1,1))
    P_feces_input = np.reshape(feces_inputs[:,3], (-1,1))
    K_feces_input = np.reshape(feces_inputs[:,4], (-1,1))
    Mg_feces_input = np.reshape(feces_inputs[:,5], (-1,1))
    Ca_feces_input = np.reshape(feces_inputs[:,6], (-1,1))
    energy_feces_input = np.reshape(feces_inputs[:,7], (-1,1))
    N_ammonia_feces_input = np.reshape(feces_inputs[:,8], (-1,1))

    # generate distributions of parameters
    # cleansing material: toilet paper
    if parameters.toilet_paper.expected == 'yes':
        toilet_paper_addition, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_addition, correlation_distributions, correlation_parameters, n_samples)
        toilet_paper_COD_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_COD_content, correlation_distributions, correlation_parameters, n_samples)
        toilet_paper_N_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_N_content, correlation_distributions, correlation_parameters, n_samples)
        toilet_paper_P_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_P_content, correlation_distributions, correlation_parameters, n_samples)
        toilet_paper_K_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_K_content, correlation_distributions, correlation_parameters, n_samples)
        toilet_paper_TS_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_paper_TS_content, correlation_distributions, correlation_parameters, n_samples)
    elif parameters.toilet_paper.expected == 'no':
        toilet_paper_addition = 0
        toilet_paper_COD_content = 0
        toilet_paper_N_content = 0
        toilet_paper_P_content = 0
        toilet_paper_K_content = 0
        toilet_paper_TS_content = 0
    # flushing water
    if parameters.flushing_water.expected == 'yes':
        flushing_water_addition, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.flushing_water_use, correlation_distributions, correlation_parameters, n_samples)
    elif parameters.flushing_water.expected == 'no':
        flushing_water_addition = 0
    # cleansing material: water
    if (parameters.cleansing_water.expected == 'yes') and (parameters.cleansing_water_hole.expected == 'no'):
        cleansing_water_addition, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.cleansing_water_use, correlation_distributions, correlation_parameters, n_samples)
    else:
        cleansing_water_addition = 0
    # desiccant
    if parameters.desiccant.expected == 'yes':
        desiccant_volume, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_volume, correlation_distributions, correlation_parameters, n_samples)
        desiccant_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_density, correlation_distributions, correlation_parameters, n_samples)
        desiccant_C_N_ratio, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_C_N_ratio, correlation_distributions, correlation_parameters, n_samples)
        desiccant_N_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_N_content, correlation_distributions, correlation_parameters, n_samples)
        desiccant_P_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_P_content, correlation_distributions, correlation_parameters, n_samples)
        desiccant_K_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_K_content, correlation_distributions, correlation_parameters, n_samples)
        desiccant_Mg_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_Mg_content, correlation_distributions, correlation_parameters, n_samples)
        desiccant_Ca_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desiccant_Ca_content, correlation_distributions, correlation_parameters, n_samples)
    elif parameters.desiccant.expected == 'no':
        desiccant_volume = 0
        desiccant_density = 0
        desiccant_N_content = 0
        desiccant_P_content = 0
        desiccant_K_content = 0
        desiccant_Mg_content = 0
        desiccant_Ca_content = 0
    # household and toilet characteristics
    household_size, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.household_size, correlation_distributions, correlation_parameters, n_samples)
    household_use_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.household_use_density, correlation_distributions, correlation_parameters, n_samples)

    # calculate resources in materials added by user (per capita per year)
    # N (kg N/cap/yr)
    additives_N = ((toilet_paper_addition * toilet_paper_N_content * 365 / 1000000)
                    + (desiccant_volume * desiccant_density * (desiccant_N_content / 100) * 365 / 1000000))
    # P (kg P/cap/yr)
    additives_P = ((toilet_paper_addition * toilet_paper_P_content * 365 / 1000000)
                    + (desiccant_volume * desiccant_density * (desiccant_P_content / 100) * 365 / 1000000))
    # K (kg K/cap/yr)
    additives_K = ((toilet_paper_addition * toilet_paper_K_content * 365 / 1000000)
                    + (desiccant_volume * desiccant_density * (desiccant_K_content / 100) * 365 / 1000000))
    # Mg (kg Mg/cap/yr)
    additives_Mg = (desiccant_volume * desiccant_density * (desiccant_Mg_content / 100) * 365 / 1000000)
    # Ca (kg Ca/cap/yr)
    additives_Ca = (desiccant_volume * desiccant_density * (desiccant_Ca_content / 100) * 365 / 1000000)
    # energy (kJ/cap/yr) - assume desiccant contains no additional energy (wood ash), so only include TP
    #   - also assume 3.86 kWh energy production per kg COD oxidized to CO2 and H20
    additives_energy = toilet_paper_addition * toilet_paper_COD_content * 365 / 1000000 * 3.86 * 3600
    # total per capita addition (kg/cap/yr)
    additives_input = ((toilet_paper_addition * toilet_paper_TS_content * 365 / 1000000)
                        + (desiccant_volume * desiccant_density * 365 / 1000000))
    # total water input (L/cap/yr)
    water_input = (flushing_water_addition + cleansing_water_addition) * 365
    # total dry addition (kg/cap/yr) - assume toilet paper and desiccant have moisture content ~ 0%
    additives_dry_input = ((toilet_paper_addition * toilet_paper_TS_content * 365 / 1000000)
                            + (desiccant_volume * desiccant_density * 365 / 1000000))


    # mix materials together to form two streams (urine, feces) from all HHs using toilet - total input per year
    # nitrogen (kg N/cap/yr)
    N_urine = (N_urine_input)
    N_feces = (N_feces_input + additives_N)
    # phosphorus (kg P/cap/yr)
    P_urine = (P_urine_input)
    P_feces = (P_feces_input + additives_P)
    # potassium (kg K/cap/yr)
    K_urine = (K_urine_input)
    K_feces = (K_feces_input + additives_K)
    # magnesium (kg Mg/cap/yr)
    Mg_urine = (Mg_urine_input)
    Mg_feces = (Mg_feces_input + additives_Mg)
    # calcium (kg Ca/cap/yr)
    Ca_urine = (Ca_urine_input)
    Ca_feces = (Ca_feces_input + additives_Ca)
    # energy (kJ/cap/yr)
    energy_urine = (energy_urine_input)
    energy_feces = (energy_feces_input + additives_energy)
    # ammonia (kg N/cap/yr)
    N_ammonia_urine = N_ammonia_urine_input
    N_ammonia_feces = N_ammonia_feces_input

    # total excreta input (kg/cap/yr)
    urine_total = (urine_input)
    # if design includes a separate hole for cleansing water, do not include it with fecal stream
    feces_total = (feces_input + additives_input + water_input)

    # total dry excreta input (kg/cap/yr)
    urine_dry = (urine_dry_input)
    feces_dry = (feces_dry_input + additives_dry_input)

    # concatenate excreta outputs to simplify function output
    urine_outputs = np.concatenate((urine_total, urine_dry, N_urine, P_urine, K_urine,
                                    Mg_urine, Ca_urine, energy_urine, N_ammonia_urine), 1)
    feces_outputs = np.concatenate((feces_total, feces_dry, N_feces, P_feces, K_feces,
                                    Mg_feces, Ca_feces, energy_feces, N_ammonia_feces), 1)

    # users
    number_users = household_use_density * household_size

    # cost for one facility
    capex, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.UDDT_capex, correlation_distributions, correlation_parameters, n_samples)
    opex_percent, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.UDDT_opex, correlation_distributions, correlation_parameters, n_samples)
    lifetime, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.toilet_lifetime, correlation_distributions, correlation_parameters, n_samples)
    annual_opex = (capex*(opex_percent/100)) / number_users
    capex_annualized = (capex * ((discount_rate*(1 + discount_rate)**lifetime)/(((1 + discount_rate)**lifetime) - 1))) / number_users
    construction_cost[:,:1] = construction_cost[:,:1] + capex_annualized
    operating_cost[:,:1] = operating_cost[:,:1] + annual_opex

    # GHG emissions from construction
    cement, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.cement_UDDT, correlation_distributions, correlation_parameters, n_samples)
    gravel, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.gravel_UDDT, correlation_distributions, correlation_parameters, n_samples)
    sand, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sand_UDDT, correlation_distributions, correlation_parameters, n_samples)
    bricks, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.bricks_UDDT, correlation_distributions, correlation_parameters, n_samples)
    plastic, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.plastic_UDDT, correlation_distributions, correlation_parameters, n_samples)
    steel, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.steel_UDDT, correlation_distributions, correlation_parameters, n_samples)
    stainless_steel_sheet, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.stainless_steel_sheet_UDDT, correlation_distributions, correlation_parameters, n_samples)
    wood, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.wood_UDDT, correlation_distributions, correlation_parameters, n_samples)

    plastic_mass, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.plastic_mass, correlation_distributions, correlation_parameters, n_samples)
    brick_volume, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.brick_volume, correlation_distributions, correlation_parameters, n_samples)
    brick_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.brick_density, correlation_distributions, correlation_parameters, n_samples)
    steel_sheet_mass, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.steel_sheet_mass, correlation_distributions, correlation_parameters, n_samples)
    gravel_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.gravel_bulk_density, correlation_distributions, correlation_parameters, n_samples)
    sand_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sand_bulk_density, correlation_distributions, correlation_parameters, n_samples)
    steel_density, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.steel_density, correlation_distributions, correlation_parameters, n_samples)

    steel_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.steel_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    stainless_steel_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.stainless_steel_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    stainless_steel_sheet_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.stainless_steel_sheet_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    plastic_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.plastic_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    gravel_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.gravel_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    sand_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sand_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    cement_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.cement_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    bricks_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.bricks_IF_GHG, correlation_distributions, correlation_parameters, n_samples)
    wood_IF_GHG, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.wood_IF_GHG, correlation_distributions, correlation_parameters, n_samples)

    cement_emissions = cement * cement_IF_GHG
    gravel_emissions = gravel*gravel_density*gravel_IF_GHG
    sand_emissions = sand*sand_density*sand_IF_GHG
    brick_emissions = bricks*brick_volume*brick_density*bricks_IF_GHG
    plastic_emissions = plastic*plastic_mass*plastic_IF_GHG
    steel_emissions = steel*steel_density*steel_IF_GHG
    stainless_steel_emissions = stainless_steel_sheet*steel_sheet_mass*(stainless_steel_IF_GHG+stainless_steel_sheet_IF_GHG)
    wood_emissions = wood * wood_IF_GHG

    construction_emissions_annual = ((cement_emissions+gravel_emissions+sand_emissions+brick_emissions+plastic_emissions+steel_emissions
                                     +stainless_steel_emissions+wood_emissions)/lifetime/number_users)

    tech_construction_emissions[:,:1] = tech_construction_emissions[:,:1] + construction_emissions_annual


    return urine_outputs, feces_outputs, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, number_users


#%% User Interface function - main function
def main(input_excel_name, urine_inputs, feces_inputs, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, n_samples, discount_rate):
    # import user interface parameters from input spreadsheet
    parameters = pd.DataFrame.transpose(pd.read_excel(input_excel_name, sheet_name = 'user_interface').set_index('parameters'))

    # define the module
    module = parameters.module.expected

    # modules are defined by text strings
    if type(module) is str:
        # dry toilet module
        if module == 'dry_toilet':
            (excreta_outputs, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions, correlation_distributions,
             correlation_parameters, number_users) = dry_toilet(urine_inputs, feces_inputs, parameters, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions,
                                   correlation_distributions, correlation_parameters, n_samples, discount_rate)
            # no separate urine/feces streams
            urine_outputs = np.full(np.shape(urine_inputs), np.nan)
            feces_outputs = np.full(np.shape(feces_inputs), np.nan)

        # UDDT module
        elif module == 'UDDT':
            (urine_outputs, feces_outputs, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions, correlation_distributions,
             correlation_parameters, number_users) = UDDT(urine_inputs, feces_inputs, parameters, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions,
                                   correlation_distributions, correlation_parameters, n_samples, discount_rate)
            # no mixed excreta stream
            excreta_outputs = np.full(np.shape(urine_inputs), np.nan)

        # if the user interface module input is not supported/valid
        else:
            raise ValueError('The specified user interface module is not valid.')

    elif type(module) is float:
        # no module specified: pass through (inputs = outputs)
        if np.isnan(module):
            # urine and feces pass through with no processing
            # they remain as per capita excretion, with no additives
            urine_outputs = urine_inputs
            feces_outputs = feces_inputs
            # no mixed excreta stream
            excreta_outputs = np.full(np.shape(urine_inputs), np.nan)

        # any other numerical input is not supported
        else:
            raise ValueError('The specified user interface module is not valid.')

    # any other data types are not supported
    else:
        raise ValueError('The specified user interface module is not valid.')

    return excreta_outputs, urine_outputs, feces_outputs, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, number_users