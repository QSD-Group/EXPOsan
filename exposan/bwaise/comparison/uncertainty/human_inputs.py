#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 17:44:38 2017

@author: John Trimmer; Yalin Li (minor modification)
"""

# This function defines the excretion inputs entering the sanitation system

import numpy as np  # import NumPy library to use mathematical functions and preserve index information
import lhs  # import lhs, which contains functions for Latin Hypercube Sampling

#%% initial inputs function
def human_inputs(initial_inputs, correlation_distributions, correlation_parameters, n_samples):
    # create uncertainty distributions using Latin Hypercube Sampling
    # use the LHS distribution function to identify and generate desired distributions
    caloric_intake, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.caloric_intake, correlation_distributions, correlation_parameters, n_samples)
    protein_vegetal_intake, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.protein_vegetal_intake, correlation_distributions, correlation_parameters, n_samples)
    protein_animal_intake, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.protein_animal_intake, correlation_distributions, correlation_parameters, n_samples)
    N_content_protein, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.N_content_protein, correlation_distributions, correlation_parameters, n_samples)
    P_content_protein_vegetal, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.P_content_protein_vegetal, correlation_distributions, correlation_parameters, n_samples)
    P_content_protein_animal, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.P_content_protein_animal, correlation_distributions, correlation_parameters, n_samples)
    K_content_caloric_intake, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.K_content_caloric_intake, correlation_distributions, correlation_parameters, n_samples)
    N_excretion, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.N_excretion, correlation_distributions, correlation_parameters, n_samples)
    P_excretion, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.P_excretion, correlation_distributions, correlation_parameters, n_samples)
    K_excretion, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.K_excretion, correlation_distributions, correlation_parameters, n_samples)
    energy_excretion, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.energy_excretion, correlation_distributions, correlation_parameters, n_samples)
    N_in_urine, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.N_in_urine, correlation_distributions, correlation_parameters, n_samples)
    P_in_urine, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.P_in_urine, correlation_distributions, correlation_parameters, n_samples)
    K_in_urine, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.K_in_urine, correlation_distributions, correlation_parameters, n_samples)
    energy_in_feces, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.energy_in_feces, correlation_distributions, correlation_parameters, n_samples)
    urine_excretion, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.urine_excretion, correlation_distributions, correlation_parameters, n_samples)
    feces_excretion, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.feces_excretion, correlation_distributions, correlation_parameters, n_samples)
    urine_moisture_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.urine_moisture_content, correlation_distributions, correlation_parameters, n_samples)
    feces_moisture_content, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.feces_moisture_content, correlation_distributions, correlation_parameters, n_samples)
    Mg_in_urine, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.Mg_in_urine, correlation_distributions, correlation_parameters, n_samples)
    Mg_in_feces, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.Mg_in_feces, correlation_distributions, correlation_parameters, n_samples)
    Ca_in_urine, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.Ca_in_urine, correlation_distributions, correlation_parameters, n_samples)
    Ca_in_feces, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.Ca_in_feces, correlation_distributions, correlation_parameters, n_samples)
    N_ammonia_in_urine, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.N_reduced_inorganic_in_urine, correlation_distributions, correlation_parameters, n_samples)
    N_ammonia_in_feces, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.N_reduced_inorganic_in_feces, correlation_distributions, correlation_parameters, n_samples)

    # calculate resources in urine and feces excreted by a person
    # urine
    # nitrogen (kg N/yr)
    N_urine_input = ((protein_vegetal_intake + protein_animal_intake) * (N_content_protein / 100)
                    * (N_excretion / 100) * (N_in_urine / 100) * 365 / 1000)
    # phosphorus (kg P/yr)
    P_urine_input = (((protein_vegetal_intake * (P_content_protein_vegetal / 100))
                    + (protein_animal_intake * (P_content_protein_animal / 100)))
                    * (P_excretion / 100) * (P_in_urine / 100) * 365 / 1000)
    # potassium (kg K/yr)
    K_urine_input = (((caloric_intake / 1000) * K_content_caloric_intake)
                    * (K_excretion / 100) * (K_in_urine / 100) * 365 / 1000)
    # magnesium (kg Mg/yr)
    Mg_urine_input = Mg_in_urine * 365 / 1000
    # calcium (kg Ca/yr)
    Ca_urine_input = Ca_in_urine * 365 / 1000
    # energy (kJ/yr)
    energy_urine_input = (caloric_intake * (energy_excretion / 100)
                            * ((100 - energy_in_feces) / 100) * 365 * 4.184)
    # total household excretion (kg/yr)
    urine_input = urine_excretion * 365 / 1000
    # total dry matter household excretion (kg/yr)
    urine_dry_input = urine_input * ((100 - urine_moisture_content) / 100)
    # total ammonia (kg N/yr)
    N_ammonia_urine = N_urine_input * (N_ammonia_in_urine / 100)

    # feces
    # nitrogen (kg N/yr)
    N_feces_input = ((protein_vegetal_intake + protein_animal_intake) * (N_content_protein / 100)
                    * (N_excretion / 100) * ((100 - N_in_urine) / 100) * 365 / 1000)
    # phosphorus (kg P/yr)
    P_feces_input = (((protein_vegetal_intake * (P_content_protein_vegetal / 100))
                    + (protein_animal_intake * (P_content_protein_animal / 100)))
                    * (P_excretion / 100) * ((100 - P_in_urine) / 100) * 365 / 1000)
    # potassium (kg K/yr)
    K_feces_input = (((caloric_intake / 1000) * K_content_caloric_intake)
                    * (K_excretion / 100) * ((100 - K_in_urine) / 100) * 365 / 1000)
    # magnesium (kg Mg/yr)
    Mg_feces_input = Mg_in_feces * 365 / 1000
    # calcium (kg Ca/yr)
    Ca_feces_input = Ca_in_feces * 365 / 1000
    # energy (kJ/yr)
    energy_feces_input = (caloric_intake * (energy_excretion / 100)
                            * (energy_in_feces / 100) * 365 * 4.184)
    # total per capita excretion (kg/yr)
    feces_input = feces_excretion * 365 / 1000
    # total dry matter per capita excretion (kg/yr)
    feces_dry_input = feces_input * ((100 - feces_moisture_content) / 100)
    # total ammonia (kg N/yr)
    N_ammonia_feces = N_feces_input * (N_ammonia_in_feces / 100)

    # concatenate all urine inputs and all feces inputs to simplify function variables
    system_urine_inputs = np.concatenate((urine_input, urine_dry_input, N_urine_input,
                                                  P_urine_input, K_urine_input, Mg_urine_input, Ca_urine_input,
                                                  energy_urine_input, N_ammonia_urine), 1)
    system_feces_inputs = np.concatenate((feces_input, feces_dry_input, N_feces_input,
                                                  P_feces_input, K_feces_input, Mg_feces_input, Ca_feces_input,
                                                  energy_feces_input, N_ammonia_feces), 1)

    return system_urine_inputs, system_feces_inputs, correlation_distributions, correlation_parameters