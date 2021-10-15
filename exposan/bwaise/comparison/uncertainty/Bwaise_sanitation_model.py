#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 11:18:42 2017

@author: John Trimmer; Yalin Li (minor modification)
"""

# This Python script evaluates resource recovery potential, economics, and GHG emissions
# in a given scenario. Inputs from the Bwaise_sanitation_inputs Excel spreadsheet define
# the scenario and all parameter values. Any changes to inputs should be implemented in
# the spreadsheet. The script outputs a summary spreadsheet, showing results for the
# overall system and for each stage in the sanitation chain. The output spreadsheet also
# provides results of sensitivity analysis (Spearman coefficients).


#%% SECTION 0 - IMPORT DATA AND CALCULATE EXCRETA INPUTS TO SANITATION SYSTEM

import pandas as pd  # import pandas for matrix data manipulation
import numpy as np  # import NumPy library to use mathematical functions and preserve index information
import copy
from scipy import stats
import lhs  # import lhs, which contains functions for Latin Hypercube Sampling
from human_inputs import human_inputs
import user_interface
import decentralized_storage
import conveyance
import treatment
import reuse_disposal

for system in ('A', 'B', 'C'):
    # import input spreadsheet
    input_excel_name = f'Bwaise_sanitation_inputs{system}.xlsx'
    initial_inputs = pd.DataFrame.transpose(pd.read_excel(input_excel_name, sheet_name = 'initial_inputs').set_index('parameters'))

    # define the number of scenarios to be generated for the uncertainty analysis
    n_samples = int(initial_inputs.n_samples.expected)

    # create a matrix to hold distributions of input parameters for which correlations will be assessed
    correlation_distributions = np.full((n_samples, n_samples), np.nan)

    correlation_parameters = np.full((n_samples, 1), np.nan)
    correlation_parameters = correlation_parameters.tolist()

    # create matrices for global parameters
    maximum_methane_emission, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.maximum_methane_emission, correlation_distributions, correlation_parameters, n_samples)
    time_full_degradation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.time_full_degradation, correlation_distributions, correlation_parameters, n_samples)
    reduction_full_degradation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.reduction_full_degradation, correlation_distributions, correlation_parameters, n_samples)

    reduction_full_degradation = reduction_full_degradation.astype('float')

    reduction_full_degradation = (1 - 10**(-reduction_full_degradation))*100

    rate_constant = (-1/time_full_degradation)*np.log((100-reduction_full_degradation)/100)

    N2O_GWP, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.N2O_GWP, correlation_distributions, correlation_parameters, n_samples)
    CH4_GWP, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.CH4_GWP, correlation_distributions, correlation_parameters, n_samples)

    exchange_rate, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.exchange_rate, correlation_distributions, correlation_parameters, n_samples)
    discount_rate, correlation_distributions, correlation_parameters = lhs.lhs_distribution(initial_inputs.discount_rate, correlation_distributions, correlation_parameters, n_samples)
    discount_rate = discount_rate/100

    # create a matrix to hold cost data (USD/cap/yr)
    construction_cost = np.full((n_samples,5), 0.0)
    operating_cost = np.full((n_samples,5), 0.0)
    income = np.full((n_samples,5), 0.0)

    # emissions matrices (kg CO2eq/cap/yr)
    direct_emissions = np.full((n_samples,5), 0.0)
    tech_construction_emissions = np.full((n_samples,5), 0.0)
    tech_operating_emissions = np.full((n_samples,5), 0.0)
    emission_offsets = np.full((n_samples,5), 0.0)

    # run human inputs function to estimate per capita input to sanitation system
    (system_liquid_inputs, system_solid_inputs, correlation_distributions,
     correlation_parameters) = human_inputs(initial_inputs, correlation_distributions,
                           correlation_parameters, n_samples)

    #%% SECTION 1 - MODULE 1: USER INTERFACE

    # run main user interface function
    (user_interface_excreta_outputs, user_interface_liquid_outputs, user_interface_solid_outputs, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions,
     correlation_distributions, correlation_parameters, number_users) = user_interface.main(input_excel_name, system_liquid_inputs,
                                 system_solid_inputs, construction_cost, operating_cost, tech_construction_emissions, tech_operating_emissions, correlation_distributions, correlation_parameters, n_samples, discount_rate)


    #%% SECTION 2 - MODULE 2: DECENTRALIZED COLLECTION AND STORAGE

    # run main decentralized collection and storage function
    (decentralized_storage_excreta_outputs, decentralized_storage_liquid_outputs, decentralized_storage_solid_outputs, direct_emissions,
     correlation_distributions, correlation_parameters, previous_storage_time) = decentralized_storage.main(input_excel_name, user_interface_excreta_outputs,
                                        user_interface_liquid_outputs, user_interface_solid_outputs, direct_emissions, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, number_users)


    #%% SECTION 3 - MODULE 3: CONVEYANCE

    # run main conveyance function
    (conveyance_excreta_outputs, conveyance_liquid_outputs, conveyance_solid_outputs, tech_operating_emissions, operating_cost,
     correlation_distributions, correlation_parameters) = conveyance.main(input_excel_name, decentralized_storage_excreta_outputs,
                                        decentralized_storage_liquid_outputs, decentralized_storage_solid_outputs, tech_operating_emissions, operating_cost, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, exchange_rate, discount_rate, number_users, previous_storage_time)


    #%% SECTION 4 - MODULE 4: TREATMENT

    # run main treatment function
    (treatment_excreta_outputs, treatment_liquid_outputs, treatment_solid_outputs, direct_emissions, tech_construction_emissions, tech_operating_emissions, construction_cost, operating_cost, biogas,
     correlation_distributions, correlation_parameters) = treatment.main(input_excel_name, conveyance_excreta_outputs,
                                        conveyance_liquid_outputs, conveyance_solid_outputs, direct_emissions, tech_construction_emissions, tech_operating_emissions, construction_cost, operating_cost, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, previous_storage_time, discount_rate, exchange_rate)


    #%% SECTION 5 - MODULE 5: REUSE OR DISPOSAL

    # run main reuse or disposal function
    (reuse_disposal_excreta_outputs, reuse_disposal_liquid_outputs, reuse_disposal_solid_outputs, emission_offsets,
     biogas_output, income, correlation_distributions, correlation_parameters) = reuse_disposal.main(input_excel_name, treatment_excreta_outputs,
                                        treatment_liquid_outputs, treatment_solid_outputs, emission_offsets, biogas, income, correlation_distributions, correlation_parameters, exchange_rate, n_samples)


    #%% FINAL SECTION - SUMMARIZE RESULTS

    # simplify correlation matrix
    i = 0
    while not (np.isnan(correlation_distributions[0,i])):
        i = i + 1

    correlation_distributions = correlation_distributions[:,:i]
    correlation_parameters = correlation_parameters[:i]

    # write results to Excel file
    writer = pd.ExcelWriter(f'Bwaise_sanitation_outputs{system}.xlsx', engine='xlsxwriter')

    # define the summary percentile value to report
    summary_percentile = initial_inputs.summary_percentile.expected

    # define the percentile values to report
    i = 1
    percentiles = np.full((5), np.nan)
    while (i <= 5) and (not np.isnan(initial_inputs.distribution_percentiles[i])):
        percentiles[i-1] = initial_inputs.distribution_percentiles[i]
        i = i + 1

    percentiles = percentiles[:i-1]

    # create strings for each percentile value, to be used when labeling output spreadsheet
    percentile_str = percentiles.tolist()
    for i in range(0, len(percentiles)):
        if (percentiles[i] % 10) == 1:
            percentile_str[i] = (str(int(percentiles[i]))+'st percentile')
        elif (percentiles[i] % 10) == 2:
            percentile_str[i] = (str(int(percentiles[i]))+'nd percentile')
        elif (percentiles[i] % 10) == 3:
            percentile_str[i] = (str(int(percentiles[i]))+'rd percentile')
        else:
            percentile_str[i] = (str(int(percentiles[i]))+'th percentile')

    # summarize results at each stage
    # mixed
    summary_excreta = np.array([np.percentile((system_liquid_inputs[:,0] + system_solid_inputs[:,0]), summary_percentile),
                              np.percentile(user_interface_excreta_outputs[:,0], summary_percentile),
                              np.percentile(decentralized_storage_excreta_outputs[:,0], summary_percentile),
                              np.percentile(conveyance_excreta_outputs[:,0], summary_percentile),
                              np.percentile(treatment_excreta_outputs[:,0], summary_percentile),
                              np.percentile(reuse_disposal_excreta_outputs[:,0], summary_percentile)])
    summary_excreta_dry = np.array([np.percentile((system_liquid_inputs[:,1] + system_solid_inputs[:,1]), summary_percentile),
                                  np.percentile(user_interface_excreta_outputs[:,1], summary_percentile),
                                  np.percentile(decentralized_storage_excreta_outputs[:,1], summary_percentile),
                                  np.percentile(conveyance_excreta_outputs[:,1], summary_percentile),
                                  np.percentile(treatment_excreta_outputs[:,1], summary_percentile),
                                  np.percentile(reuse_disposal_excreta_outputs[:,1], summary_percentile)])
    summary_excreta_N = np.array([np.percentile((system_liquid_inputs[:,2] + system_solid_inputs[:,2]), summary_percentile),
                                np.percentile(user_interface_excreta_outputs[:,2], summary_percentile),
                                np.percentile(decentralized_storage_excreta_outputs[:,2], summary_percentile),
                                np.percentile(conveyance_excreta_outputs[:,2], summary_percentile),
                                np.percentile(treatment_excreta_outputs[:,2], summary_percentile),
                                np.percentile(reuse_disposal_excreta_outputs[:,2], summary_percentile)])
    summary_excreta_P = np.array([np.percentile((system_liquid_inputs[:,3] + system_solid_inputs[:,3]), summary_percentile),
                                np.percentile(user_interface_excreta_outputs[:,3], summary_percentile),
                                np.percentile(decentralized_storage_excreta_outputs[:,3], summary_percentile),
                                np.percentile(conveyance_excreta_outputs[:,3], summary_percentile),
                                np.percentile(treatment_excreta_outputs[:,3], summary_percentile),
                                np.percentile(reuse_disposal_excreta_outputs[:,3], summary_percentile)])
    summary_excreta_K = np.array([np.percentile((system_liquid_inputs[:,4] + system_solid_inputs[:,4]), summary_percentile),
                                np.percentile(user_interface_excreta_outputs[:,4], summary_percentile),
                                np.percentile(decentralized_storage_excreta_outputs[:,4], summary_percentile),
                                np.percentile(conveyance_excreta_outputs[:,4], summary_percentile),
                                np.percentile(treatment_excreta_outputs[:,4], summary_percentile),
                                np.percentile(reuse_disposal_excreta_outputs[:,4], summary_percentile)])
    summary_excreta_Mg = np.array([np.percentile((system_liquid_inputs[:,5] + system_solid_inputs[:,5]), summary_percentile),
                                 np.percentile(user_interface_excreta_outputs[:,5], summary_percentile),
                                 np.percentile(decentralized_storage_excreta_outputs[:,5], summary_percentile),
                                 np.percentile(conveyance_excreta_outputs[:,5], summary_percentile),
                                 np.percentile(treatment_excreta_outputs[:,5], summary_percentile),
                                 np.percentile(reuse_disposal_excreta_outputs[:,5], summary_percentile)])
    summary_excreta_Ca = np.array([np.percentile((system_liquid_inputs[:,6] + system_solid_inputs[:,6]), summary_percentile),
                                 np.percentile(user_interface_excreta_outputs[:,6], summary_percentile),
                                 np.percentile(decentralized_storage_excreta_outputs[:,6], summary_percentile),
                                 np.percentile(conveyance_excreta_outputs[:,6], summary_percentile),
                                 np.percentile(treatment_excreta_outputs[:,6], summary_percentile),
                                 np.percentile(reuse_disposal_excreta_outputs[:,6], summary_percentile)])
    summary_excreta_energy = np.array([np.percentile((system_liquid_inputs[:,7] + system_solid_inputs[:,7]), summary_percentile),
                                     np.percentile(user_interface_excreta_outputs[:,7], summary_percentile),
                                     np.percentile(decentralized_storage_excreta_outputs[:,7], summary_percentile),
                                     np.percentile(conveyance_excreta_outputs[:,7], summary_percentile),
                                     np.percentile(treatment_excreta_outputs[:,7], summary_percentile),
                                     np.percentile(reuse_disposal_excreta_outputs[:,7], summary_percentile)])
    summary_excreta_ammonia_N = np.array([np.percentile((system_liquid_inputs[:,8] + system_solid_inputs[:,8]), summary_percentile),
                                        np.percentile(user_interface_excreta_outputs[:,8], summary_percentile),
                                        np.percentile(decentralized_storage_excreta_outputs[:,8], summary_percentile),
                                        np.percentile(conveyance_excreta_outputs[:,8], summary_percentile),
                                        np.percentile(treatment_excreta_outputs[:,8], summary_percentile),
                                        np.percentile(reuse_disposal_excreta_outputs[:,8], summary_percentile)])
    if initial_inputs.excreta_decentralized_storage.expected == 'single_pit':
        summary_excreta_filling_time = np.array([np.nan, np.nan, np.percentile(decentralized_storage_excreta_outputs[:,9], summary_percentile), np.nan, np.nan, np.nan])
    else:
        summary_excreta_filling_time = np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

    if (initial_inputs.excreta_reuse_disposal_1.expected == 'biogas_combustion') or (initial_inputs.excreta_reuse_disposal_2.expected == 'biogas_combustion'):
        summary_excreta_biogas = np.array([np.nan, np.nan, np.nan, np.nan, np.percentile(biogas, summary_percentile), np.percentile(biogas_output, summary_percentile)])
    else:
        summary_excreta_biogas = np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

    summary_direct_emissions = np.array([0, np.percentile(direct_emissions[:,0], summary_percentile),
                                        np.percentile(direct_emissions[:,1], summary_percentile),
                                        np.percentile(direct_emissions[:,2], summary_percentile),
                                        np.percentile(direct_emissions[:,3], summary_percentile),
                                        np.percentile(direct_emissions[:,4], summary_percentile)])
    summary_tech_construction_emissions = np.array([0, np.percentile(tech_construction_emissions[:,0], summary_percentile),
                                        np.percentile(tech_construction_emissions[:,1], summary_percentile),
                                        np.percentile(tech_construction_emissions[:,2], summary_percentile),
                                        np.percentile(tech_construction_emissions[:,3], summary_percentile),
                                        np.percentile(tech_construction_emissions[:,4], summary_percentile)])
    summary_tech_operating_emissions = np.array([0, np.percentile(tech_operating_emissions[:,0], summary_percentile),
                                        np.percentile(tech_operating_emissions[:,1], summary_percentile),
                                        np.percentile(tech_operating_emissions[:,2], summary_percentile),
                                        np.percentile(tech_operating_emissions[:,3], summary_percentile),
                                        np.percentile(tech_operating_emissions[:,4], summary_percentile)])
    summary_emission_offsets = np.array([0, np.percentile(emission_offsets[:,0], summary_percentile),
                                        np.percentile(emission_offsets[:,1], summary_percentile),
                                        np.percentile(emission_offsets[:,2], summary_percentile),
                                        np.percentile(emission_offsets[:,3], summary_percentile),
                                        np.percentile(emission_offsets[:,4], summary_percentile)])
    summary_construction_cost = np.array([0, np.percentile(construction_cost[:,0], summary_percentile),
                                        np.percentile(construction_cost[:,1], summary_percentile),
                                        np.percentile(construction_cost[:,2], summary_percentile),
                                        np.percentile(construction_cost[:,3], summary_percentile),
                                        np.percentile(construction_cost[:,4], summary_percentile)])
    summary_operating_cost = np.array([0, np.percentile(operating_cost[:,0], summary_percentile),
                                        np.percentile(operating_cost[:,1], summary_percentile),
                                        np.percentile(operating_cost[:,2], summary_percentile),
                                        np.percentile(operating_cost[:,3], summary_percentile),
                                        np.percentile(operating_cost[:,4], summary_percentile)])
    summary_income = np.array([0, np.percentile(income[:,0], summary_percentile),
                                        np.percentile(income[:,1], summary_percentile),
                                        np.percentile(income[:,2], summary_percentile),
                                        np.percentile(income[:,3], summary_percentile),
                                        np.percentile(income[:,4], summary_percentile)])

    summary_excreta_list = [summary_excreta, summary_excreta_dry, summary_excreta_N,
                                        summary_excreta_P, summary_excreta_K, summary_excreta_Mg,
                                        summary_excreta_Ca, summary_excreta_energy, summary_excreta_ammonia_N,
                                        summary_excreta_filling_time, summary_excreta_biogas, summary_direct_emissions, summary_tech_construction_emissions, summary_tech_operating_emissions, summary_emission_offsets, summary_construction_cost, summary_operating_cost, summary_income]
    excreta_index_list = ['Excreta (kg/cap/yr)', 'Excreta dry matter (kg/cap/yr)', 'Excreta N (kg/cap/yr)',
                                                'Excreta P (kg/cap/yr)', 'Excreta K (kg/cap/yr)', 'Excreta Mg (kg/cap/yr)',
                                                'Excreta Ca (kg/cap/yr)', 'Excreta energy (kJ/cap/yr)', 'Excreta ammonia as N (kg/cap/yr)',
                                                'Pit filling time (years)', 'Biogas energy (kJ/cap/yr)', 'Direct emissions (kg CO2eq/cap/yr)', 'Tech construction emissions (kg CO2eq/cap/yr)', 'Tech operating emissions (kg CO2eq/cap/yr)', 'Emission offsets (kg CO2eq/cap/yr)', 'Construction cost (USD/cap/yr)', 'Operating cost (USD/cap/yr)', 'Income (USD/cap/yr)']

    summary_excreta_dataframe = pd.DataFrame(summary_excreta_list,
                                            index=excreta_index_list,
                                            columns=['Single user inputs', 'User interface outputs', 'Decentralized storage outputs',
                                                     'Conveyance outputs', 'Treatment outputs', 'Reuse and disposal outputs'])
    summary_excreta_dataframe.to_excel(writer, sheet_name='summary', startcol=0)

    recovery_excreta = (reuse_disposal_excreta_outputs[:,0:9] / user_interface_excreta_outputs)
    excreta_biogas_recovery = biogas_output / user_interface_excreta_outputs[:,7:8]
    total_direct_emissions = np.reshape(np.nansum(direct_emissions, axis=1), (-1,1))
    total_tech_construction_emissions = np.reshape(np.nansum(tech_construction_emissions, axis=1), (-1,1))
    total_tech_operating_emissions = np.reshape(np.nansum(tech_operating_emissions, axis=1), (-1,1))
    total_emission_offsets = np.reshape(np.nansum(emission_offsets, axis=1), (-1,1))
    total_construction_cost = np.reshape(np.nansum(construction_cost, axis=1), (-1,1))
    total_operating_cost = np.reshape(np.nansum(operating_cost, axis=1), (-1,1))
    total_income = np.reshape(np.nansum(income, axis=1), (-1,1))
    recovery_excreta = np.concatenate((recovery_excreta, np.full(np.shape(excreta_biogas_recovery), np.nan), excreta_biogas_recovery, total_direct_emissions, total_tech_construction_emissions, total_tech_operating_emissions, total_emission_offsets, total_construction_cost, total_operating_cost, total_income), 1)

    summary_recovery_excreta_dataframe = pd.DataFrame(np.transpose(np.percentile(recovery_excreta, percentiles, axis=0)),
                                                index=['Recovered excreta mass (%)', 'Recovered excreta dry matter (%)', 'Recovered excreta N (%)',
                                                   'Recovered excreta P (%)', 'Recovered excreta K (%)', 'Recovered excreta Mg (%)',
                                                   'Recovered excreta Ca (%)', 'Recovered energy in sludge (%)', 'Recovered excreta ammonia as N (%)', '', 'Recovered energy in biogas', 'Total direct emissions (kg CO2eq/cap/yr)', 'Total tech construction emissions (kg CO2eq/cap/yr)', 'Total tech operating emissions (kg CO2eq/cap/yr)', 'Total emission offsets (kg CO2eq/cap/yr)', 'Total construction cost (USD/cap/yr)', 'Total operating cost (USD/cap/yr)', 'Total income (USD/cap/yr)'],
                                                columns=percentile_str)
    summary_recovery_excreta_dataframe.to_excel(writer, sheet_name='summary', startcol=len(summary_excreta)+2)

    # liquid
    summary_liquid = np.array([np.percentile(system_liquid_inputs[:,0], summary_percentile),
                              np.percentile(user_interface_liquid_outputs[:,0], summary_percentile),
                              np.percentile(decentralized_storage_liquid_outputs[:,0], summary_percentile),
                              np.percentile(conveyance_liquid_outputs[:,0], summary_percentile),
                              np.percentile(treatment_liquid_outputs[:,0], summary_percentile),
                              np.percentile(reuse_disposal_liquid_outputs[:,0], summary_percentile)])
    summary_liquid_dry = np.array([np.percentile(system_liquid_inputs[:,1], summary_percentile),
                                  np.percentile(user_interface_liquid_outputs[:,1], summary_percentile),
                                  np.percentile(decentralized_storage_liquid_outputs[:,1], summary_percentile),
                                  np.percentile(conveyance_liquid_outputs[:,1], summary_percentile),
                                  np.percentile(treatment_liquid_outputs[:,1], summary_percentile),
                                  np.percentile(reuse_disposal_liquid_outputs[:,1], summary_percentile)])
    summary_liquid_N = np.array([np.percentile(system_liquid_inputs[:,2], summary_percentile),
                                np.percentile(user_interface_liquid_outputs[:,2], summary_percentile),
                                np.percentile(decentralized_storage_liquid_outputs[:,2], summary_percentile),
                                np.percentile(conveyance_liquid_outputs[:,2], summary_percentile),
                                np.percentile(treatment_liquid_outputs[:,2], summary_percentile),
                                np.percentile(reuse_disposal_liquid_outputs[:,2], summary_percentile)])
    summary_liquid_P = np.array([np.percentile(system_liquid_inputs[:,3], summary_percentile),
                                np.percentile(user_interface_liquid_outputs[:,3], summary_percentile),
                                np.percentile(decentralized_storage_liquid_outputs[:,3], summary_percentile),
                                np.percentile(conveyance_liquid_outputs[:,3], summary_percentile),
                                np.percentile(treatment_liquid_outputs[:,3], summary_percentile),
                                np.percentile(reuse_disposal_liquid_outputs[:,3], summary_percentile)])
    summary_liquid_K = np.array([np.percentile(system_liquid_inputs[:,4], summary_percentile),
                                np.percentile(user_interface_liquid_outputs[:,4], summary_percentile),
                                np.percentile(decentralized_storage_liquid_outputs[:,4], summary_percentile),
                                np.percentile(conveyance_liquid_outputs[:,4], summary_percentile),
                                np.percentile(treatment_liquid_outputs[:,4], summary_percentile),
                                np.percentile(reuse_disposal_liquid_outputs[:,4], summary_percentile)])
    summary_liquid_Mg = np.array([np.percentile(system_liquid_inputs[:,5], summary_percentile),
                                 np.percentile(user_interface_liquid_outputs[:,5], summary_percentile),
                                 np.percentile(decentralized_storage_liquid_outputs[:,5], summary_percentile),
                                 np.percentile(conveyance_liquid_outputs[:,5], summary_percentile),
                                 np.percentile(treatment_liquid_outputs[:,5], summary_percentile),
                                 np.percentile(reuse_disposal_liquid_outputs[:,5], summary_percentile)])
    summary_liquid_Ca = np.array([np.percentile(system_liquid_inputs[:,6], summary_percentile),
                                 np.percentile(user_interface_liquid_outputs[:,6], summary_percentile),
                                 np.percentile(decentralized_storage_liquid_outputs[:,6], summary_percentile),
                                 np.percentile(conveyance_liquid_outputs[:,6], summary_percentile),
                                 np.percentile(treatment_liquid_outputs[:,6], summary_percentile),
                                 np.percentile(reuse_disposal_liquid_outputs[:,6], summary_percentile)])
    summary_liquid_energy = np.array([np.percentile(system_liquid_inputs[:,7], summary_percentile),
                                     np.percentile(user_interface_liquid_outputs[:,7], summary_percentile),
                                     np.percentile(decentralized_storage_liquid_outputs[:,7], summary_percentile),
                                     np.percentile(conveyance_liquid_outputs[:,7], summary_percentile),
                                     np.percentile(treatment_liquid_outputs[:,7], summary_percentile),
                                     np.percentile(reuse_disposal_liquid_outputs[:,7], summary_percentile)])
    summary_liquid_ammonia_N = np.array([np.percentile(system_liquid_inputs[:,8], summary_percentile),
                                        np.percentile(user_interface_liquid_outputs[:,8], summary_percentile),
                                        np.percentile(decentralized_storage_liquid_outputs[:,8], summary_percentile),
                                        np.percentile(conveyance_liquid_outputs[:,8], summary_percentile),
                                        np.percentile(treatment_liquid_outputs[:,8], summary_percentile),
                                        np.percentile(reuse_disposal_liquid_outputs[:,8], summary_percentile)])

    if initial_inputs.liquid_decentralized_storage.expected == 'storage_tank':
        summary_liquid_filling_time = np.array([np.nan, np.nan,
                                                np.percentile(decentralized_storage_liquid_outputs[:,9], summary_percentile),
                                                np.nan, np.nan, np.nan])
        summary_liquid_treatment_time = np.array([np.nan, np.nan,
                                                  np.percentile(decentralized_storage_liquid_outputs[:,10], summary_percentile),
                                                  np.nan, np.nan, np.nan])
        summary_liquid_treatment_volume = np.array([np.nan, np.nan,
                                                    np.percentile(decentralized_storage_liquid_outputs[:,11], summary_percentile),
                                                    np.nan, np.nan, np.nan])
    else:
        summary_liquid_filling_time = np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
        summary_liquid_treatment_time = np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
        summary_liquid_treatment_volume = np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

    summary_liquid_list = [summary_liquid, summary_liquid_dry, summary_liquid_N,
                                        summary_liquid_P, summary_liquid_K, summary_liquid_Mg,
                                        summary_liquid_Ca, summary_liquid_energy, summary_liquid_ammonia_N,
                                        summary_liquid_filling_time, summary_liquid_treatment_time, summary_liquid_treatment_volume]
    liquid_index_list = ['liquid (kg/yr)', 'liquid dry matter (kg/yr)', 'liquid N (kg/yr)',
                                                'liquid P (kg/yr)', 'liquid K (kg/yr)', 'liquid Mg (kg/yr)',
                                                'liquid Ca (kg/yr)', 'liquid energy (kJ/yr)', 'liquid ammonia as N (kg/yr)',
                                                'Liquid container filling time (days)', 'liquid treatment time (days)', 'liquid treatment volume (L)']

    summary_liquid_dataframe = pd.DataFrame(summary_liquid_list,
                                            index=liquid_index_list,
                                            columns=['Single user inputs', 'User interface outputs', 'Decentralized storage outputs',
                                                     'Conveyance outputs', 'Treatment outputs', 'Reuse and disposal outputs'])
    summary_liquid_dataframe.to_excel(writer, sheet_name='summary', startcol=0, startrow=len(excreta_index_list)+2)

    # solids
    summary_solid = np.array([np.percentile(system_solid_inputs[:,0], summary_percentile),
                              np.percentile(user_interface_solid_outputs[:,0], summary_percentile),
                              np.percentile(decentralized_storage_solid_outputs[:,0], summary_percentile),
                              np.percentile(conveyance_solid_outputs[:,0], summary_percentile),
                              np.percentile(treatment_solid_outputs[:,0], summary_percentile),
                              np.percentile(reuse_disposal_solid_outputs[:,0], summary_percentile)])
    summary_solid_dry = np.array([np.percentile(system_solid_inputs[:,1], summary_percentile),
                                  np.percentile(user_interface_solid_outputs[:,1], summary_percentile),
                                  np.percentile(decentralized_storage_solid_outputs[:,1], summary_percentile),
                                  np.percentile(conveyance_solid_outputs[:,1], summary_percentile),
                                  np.percentile(treatment_solid_outputs[:,1], summary_percentile),
                                  np.percentile(reuse_disposal_solid_outputs[:,1], summary_percentile)])
    summary_solid_N = np.array([np.percentile(system_solid_inputs[:,2], summary_percentile),
                                np.percentile(user_interface_solid_outputs[:,2], summary_percentile),
                                np.percentile(decentralized_storage_solid_outputs[:,2], summary_percentile),
                                np.percentile(conveyance_solid_outputs[:,2], summary_percentile),
                                np.percentile(treatment_solid_outputs[:,2], summary_percentile),
                                np.percentile(reuse_disposal_solid_outputs[:,2], summary_percentile)])
    summary_solid_P = np.array([np.percentile(system_solid_inputs[:,3], summary_percentile),
                                np.percentile(user_interface_solid_outputs[:,3], summary_percentile),
                                np.percentile(decentralized_storage_solid_outputs[:,3], summary_percentile),
                                np.percentile(conveyance_solid_outputs[:,3], summary_percentile),
                                np.percentile(treatment_solid_outputs[:,3], summary_percentile),
                                np.percentile(reuse_disposal_solid_outputs[:,3], summary_percentile)])
    summary_solid_K = np.array([np.percentile(system_solid_inputs[:,4], summary_percentile),
                                np.percentile(user_interface_solid_outputs[:,4], summary_percentile),
                                np.percentile(decentralized_storage_solid_outputs[:,4], summary_percentile),
                                np.percentile(conveyance_solid_outputs[:,4], summary_percentile),
                                np.percentile(treatment_solid_outputs[:,4], summary_percentile),
                                np.percentile(reuse_disposal_solid_outputs[:,4], summary_percentile)])
    summary_solid_Mg = np.array([np.percentile(system_solid_inputs[:,5], summary_percentile),
                                 np.percentile(user_interface_solid_outputs[:,5], summary_percentile),
                                 np.percentile(decentralized_storage_solid_outputs[:,5], summary_percentile),
                                 np.percentile(conveyance_solid_outputs[:,5], summary_percentile),
                                 np.percentile(treatment_solid_outputs[:,5], summary_percentile),
                                 np.percentile(reuse_disposal_solid_outputs[:,5], summary_percentile)])
    summary_solid_Ca = np.array([np.percentile(system_solid_inputs[:,6], summary_percentile),
                                 np.percentile(user_interface_solid_outputs[:,6], summary_percentile),
                                 np.percentile(decentralized_storage_solid_outputs[:,6], summary_percentile),
                                 np.percentile(conveyance_solid_outputs[:,6], summary_percentile),
                                 np.percentile(treatment_solid_outputs[:,6], summary_percentile),
                                 np.percentile(reuse_disposal_solid_outputs[:,6], summary_percentile)])
    summary_solid_energy = np.array([np.percentile(system_solid_inputs[:,7], summary_percentile),
                                     np.percentile(user_interface_solid_outputs[:,7], summary_percentile),
                                     np.percentile(decentralized_storage_solid_outputs[:,7], summary_percentile),
                                     np.percentile(conveyance_solid_outputs[:,7], summary_percentile),
                                     np.percentile(treatment_solid_outputs[:,7], summary_percentile),
                                     np.percentile(reuse_disposal_solid_outputs[:,7], summary_percentile)])
    summary_solid_ammonia_N = np.array([np.percentile(system_solid_inputs[:,8], summary_percentile),
                                        np.percentile(user_interface_solid_outputs[:,8], summary_percentile),
                                        np.percentile(decentralized_storage_solid_outputs[:,8], summary_percentile),
                                        np.percentile(conveyance_solid_outputs[:,8], summary_percentile),
                                        np.percentile(treatment_solid_outputs[:,8], summary_percentile),
                                        np.percentile(reuse_disposal_solid_outputs[:,8], summary_percentile)])
    if initial_inputs.solid_decentralized_storage.expected == 'dehydration_vault':
        summary_solid_vault_volume = np.array([np.nan, np.nan,
                                           np.percentile(decentralized_storage_solid_outputs[:,9], summary_percentile),
                                           np.nan, np.nan, np.nan,])
    else:
        summary_solid_vault_volume = np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

    if (initial_inputs.solid_reuse_disposal_1.expected == 'biogas_combustion') or (initial_inputs.solid_reuse_disposal_2.expected == 'biogas_combustion'):
        summary_solid_biogas = np.array([np.nan, np.nan, np.nan, np.nan,
                                           np.percentile(biogas, summary_percentile), np.percentile (biogas_output, summary_percentile)])
    else:
        summary_solid_biogas = np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

    summary_solid_list = [summary_solid, summary_solid_dry, summary_solid_N,
                                        summary_solid_P, summary_solid_K, summary_solid_Mg,
                                        summary_solid_Ca, summary_solid_energy, summary_solid_ammonia_N,
                                        summary_solid_vault_volume, summary_solid_biogas]
    solid_index_list = ['solid (kg/yr)', 'solid dry matter (kg/yr)', 'solid N (kg/yr)',
                                               'solid P (kg/yr)', 'solid K (kg/yr)', 'solid Mg (kg/yr)',
                                               'solid Ca (kg/yr)', 'solid energy (kJ/yr)', 'solid ammonia as N (kg/yr)',
                                               'solid vault volume (m3)', 'Biogas energy (kJ/yr)']

    summary_solid_dataframe = pd.DataFrame(summary_solid_list,
                                            index=solid_index_list,
                                            columns=['Single user inputs', 'User interface outputs', 'Decentralized storage outputs',
                                                     'Conveyance outputs', 'Treatment outputs', 'Reuse and disposal outputs'])
    summary_solid_dataframe.to_excel(writer, sheet_name='summary', startcol=0, startrow=len(liquid_index_list)+len(excreta_index_list)+4)

    if initial_inputs.user_interface.expected == 'dry_toilet':
        recovery_liquid = (reuse_disposal_liquid_outputs[:,0:9] / user_interface_excreta_outputs)
        recovery_solid = (reuse_disposal_solid_outputs[:,0:9] / user_interface_excreta_outputs)
        solid_biogas_recovery = (biogas_output / user_interface_excreta_outputs[:,7:8])
        recovery_solid = np.concatenate((recovery_solid, np.full(np.shape(solid_biogas_recovery), np.nan), solid_biogas_recovery), 1)
        recovery = ((reuse_disposal_solid_outputs[:,0:9] + reuse_disposal_liquid_outputs[:,0:9])
                      / user_interface_excreta_outputs)
    elif initial_inputs.user_interface.expected == 'UDDT':
        recovery_liquid = (reuse_disposal_liquid_outputs[:,0:9] / (user_interface_solid_outputs + user_interface_liquid_outputs))
        recovery_solid = (reuse_disposal_solid_outputs[:,0:9] / (user_interface_solid_outputs + user_interface_liquid_outputs))
        solid_biogas_recovery = (biogas_output / (user_interface_solid_outputs[:,7:8] + user_interface_liquid_outputs[:,7:8]))
        recovery_solid = np.concatenate((recovery_solid, np.full(np.shape(solid_biogas_recovery), np.nan), solid_biogas_recovery), 1)
        recovery = ((reuse_disposal_solid_outputs[:,0:9] + reuse_disposal_liquid_outputs[:,0:9])
                      / (user_interface_solid_outputs + user_interface_liquid_outputs))

    recovery = np.concatenate((recovery, np.full(np.shape(solid_biogas_recovery), np.nan), solid_biogas_recovery, (total_direct_emissions+total_tech_construction_emissions+total_tech_operating_emissions-total_emission_offsets), (total_construction_cost+total_operating_cost-total_income), np.full(np.shape(solid_biogas_recovery), np.nan),
                               (total_direct_emissions/(total_direct_emissions+total_tech_construction_emissions+total_tech_operating_emissions)), (total_tech_construction_emissions/(total_direct_emissions+total_tech_construction_emissions+total_tech_operating_emissions)), (total_tech_operating_emissions/(total_direct_emissions+total_tech_construction_emissions+total_tech_operating_emissions)), (total_emission_offsets/(total_direct_emissions+total_tech_construction_emissions+total_tech_operating_emissions)),
                               np.full(np.shape(solid_biogas_recovery), np.nan), (total_construction_cost/(total_construction_cost+total_operating_cost)), (total_operating_cost/(total_construction_cost+total_operating_cost)), (total_income/(total_construction_cost+total_operating_cost))), 1)

    summary_recovery_liquid_dataframe = pd.DataFrame(np.transpose(np.percentile(recovery_liquid, percentiles, axis=0)),
                                                index=['Recovered liquid mass (%)', 'Recovered liquid dry matter (%)', 'Recovered liquid N (%)',
                                                   'Recovered liquid P (%)', 'Recovered liquid K (%)', 'Recovered liquid Mg (%)',
                                                   'Recovered liquid Ca (%)', 'Recovered liquid energy (%)', 'Recovered liquid ammonia as N (%)'],
                                                columns=percentile_str)
    summary_recovery_liquid_dataframe.to_excel(writer, sheet_name='summary', startcol=len(summary_liquid)+2, startrow=len(excreta_index_list)+2)

    summary_recovery_solid_dataframe = pd.DataFrame(np.transpose(np.percentile(recovery_solid, percentiles, axis=0)),
                                                index=['Recovered solid mass (%)', 'Recovered solid dry matter (%)', 'Recovered solid N (%)',
                                                   'Recovered solid P (%)', 'Recovered solid K (%)', 'Recovered solid Mg (%)',
                                                   'Recovered solid Ca (%)', 'Recovered solid energy (%)', 'Recovered solid ammonia as N (%)', '', 'Recovered biogas (%)'],
                                                columns=percentile_str)
    summary_recovery_solid_dataframe.to_excel(writer, sheet_name='summary', startcol=len(summary_solid)+2, startrow=len(liquid_index_list)+len(excreta_index_list)+4)

    summary_recovery_dataframe = pd.DataFrame(np.transpose(np.percentile(recovery, percentiles, axis=0)),
                                                index=['Recovered mass (%)', 'Recovered dry matter (%)', 'Recovered N (%)',
                                                   'Recovered P (%)', 'Recovered K (%)', 'Recovered Mg (%)',
                                                   'Recovered Ca (%)', 'Recovered energy in sludge (%)', 'Recovered ammonia as N (%)', '', 'Recovered energy in biogas (%)', 'Total net emissions (kg CO2eq/cap/yr)', 'Total net cost (USD/cap/yr)', '',
                                                   'Direct emissions (% of total)', 'Tech construction emissions (% of total)', 'Tech operating emissions (% of total)', 'Emission offsets (% of total)', '', 'Construction cost (% of total)', 'Operating cost (% of total)', 'Income (% of total cost)'],
                                                columns=percentile_str)
    summary_recovery_dataframe.to_excel(writer, sheet_name='summary', startcol=len(summary_solid)+2, startrow=len(liquid_index_list)+len(excreta_index_list)+len(solid_index_list)+6)

    # initial system inputs
    system_liquid_in = np.percentile(system_liquid_inputs[:,0], percentiles)
    system_liquid_dry_in = np.percentile(system_liquid_inputs[:,1], percentiles)
    system_liquid_N_in = np.percentile(system_liquid_inputs[:,2], percentiles)
    system_liquid_P_in = np.percentile(system_liquid_inputs[:,3], percentiles)
    system_liquid_K_in = np.percentile(system_liquid_inputs[:,4], percentiles)
    system_liquid_Mg_in = np.percentile(system_liquid_inputs[:,5], percentiles)
    system_liquid_Ca_in = np.percentile(system_liquid_inputs[:,6], percentiles)
    system_liquid_energy_in = np.percentile(system_liquid_inputs[:,7], percentiles)
    system_liquid_ammonia_N_in = np.percentile(system_liquid_inputs[:,8], percentiles)
    system_liquid_dataframe = pd.DataFrame([system_liquid_in, system_liquid_dry_in, system_liquid_N_in,
                                           system_liquid_P_in, system_liquid_K_in, system_liquid_Mg_in,
                                           system_liquid_Ca_in, system_liquid_energy_in, system_liquid_ammonia_N_in],
                                           index=['liquid (kg/yr)', 'liquid dry matter (kg/yr)', 'liquid N (kg/yr)',
                                                    'liquid P (kg/yr)', 'liquid K (kg/yr)', 'liquid Mg (kg/yr)',
                                                    'liquid Ca (kg/yr)', 'liquid energy (kJ/yr)', 'liquid ammonia as N (kg/yr)'],
                                           columns=percentile_str)
    system_liquid_dataframe.to_excel(writer, sheet_name='system_inputs', startcol=0)

    system_solid_in = np.percentile(system_solid_inputs[:,0], percentiles)
    system_solid_dry_in = np.percentile(system_solid_inputs[:,1], percentiles)
    system_solid_N_in = np.percentile(system_solid_inputs[:,2], percentiles)
    system_solid_P_in = np.percentile(system_solid_inputs[:,3], percentiles)
    system_solid_K_in = np.percentile(system_solid_inputs[:,4], percentiles)
    system_solid_Mg_in = np.percentile(system_solid_inputs[:,5], percentiles)
    system_solid_Ca_in = np.percentile(system_solid_inputs[:,6], percentiles)
    system_solid_energy_in = np.percentile(system_solid_inputs[:,7], percentiles)
    system_solid_ammonia_N_in = np.percentile(system_solid_inputs[:,8], percentiles)
    system_solid_dataframe = pd.DataFrame([system_solid_in, system_solid_dry_in, system_solid_N_in,
                                           system_solid_P_in, system_solid_K_in, system_solid_Mg_in,
                                           system_solid_Ca_in, system_solid_energy_in, system_solid_ammonia_N_in],
                                           index=['solid (kg/yr)', 'solid dry matter (kg/yr)', 'solid N (kg/yr)',
                                                    'solid P (kg/yr)', 'solid K (kg/yr)', 'solid Mg (kg/yr)',
                                                    'solid Ca (kg/yr)', 'solid energy (kJ/yr)', 'solid ammonia as N (kg/yr)'],
                                           columns=percentile_str)
    system_solid_dataframe.to_excel(writer, sheet_name='system_inputs', startcol=len(percentiles)+2)

    # user interface outputs
    if not np.isnan(user_interface_excreta_outputs[0,0]):
        user_interface_excreta_out = np.percentile(user_interface_excreta_outputs[:,0], percentiles)
        user_interface_excreta_dry_out = np.percentile(user_interface_excreta_outputs[:,1], percentiles)
        user_interface_excreta_N_out = np.percentile(user_interface_excreta_outputs[:,2], percentiles)
        user_interface_excreta_P_out = np.percentile(user_interface_excreta_outputs[:,3], percentiles)
        user_interface_excreta_K_out = np.percentile(user_interface_excreta_outputs[:,4], percentiles)
        user_interface_excreta_Mg_out = np.percentile(user_interface_excreta_outputs[:,5], percentiles)
        user_interface_excreta_Ca_out = np.percentile(user_interface_excreta_outputs[:,6], percentiles)
        user_interface_excreta_energy_out = np.percentile(user_interface_excreta_outputs[:,7], percentiles)
        user_interface_excreta_ammonia_N_out = np.percentile(user_interface_excreta_outputs[:,8], percentiles)

        user_interface_excreta_construction_cost = np.percentile(construction_cost[:,0], percentiles)
        user_interface_excreta_operating_cost = np.percentile(operating_cost[:,0], percentiles)

        user_interface_excreta_tech_construction_emissions = np.percentile(tech_construction_emissions[:,0], percentiles)
        user_interface_excreta_tech_operating_emissions = np.percentile(tech_operating_emissions[:,0], percentiles)

        user_interface_excreta_construction_cost_percent = np.percentile(construction_cost[:,0]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)
        user_interface_excreta_operating_cost_percent = np.percentile(operating_cost[:,0]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)

        user_interface_excreta_tech_construction_emissions_percent = np.percentile(tech_construction_emissions[:,0]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)
        user_interface_excreta_tech_operating_emissions_percent = np.percentile(tech_operating_emissions[:,0]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)


        user_interface_excreta_dataframe = pd.DataFrame([user_interface_excreta_out, user_interface_excreta_dry_out, user_interface_excreta_N_out,
                                           user_interface_excreta_P_out, user_interface_excreta_K_out, user_interface_excreta_Mg_out,
                                           user_interface_excreta_Ca_out, user_interface_excreta_energy_out, user_interface_excreta_ammonia_N_out, user_interface_excreta_construction_cost, user_interface_excreta_operating_cost,
                                           user_interface_excreta_tech_construction_emissions, user_interface_excreta_tech_operating_emissions, user_interface_excreta_construction_cost_percent, user_interface_excreta_operating_cost_percent,
                                           user_interface_excreta_tech_construction_emissions_percent, user_interface_excreta_tech_operating_emissions_percent],
                                           index=['Excreta (kg/yr)', 'Excreta dry matter (kg/yr)', 'Excreta N (kg/yr)',
                                                    'Excreta P (kg/yr)', 'Excreta K (kg/yr)', 'Excreta Mg (kg/yr)',
                                                    'Excreta Ca (kg/yr)', 'Excreta energy (kJ/yr)', 'Excreta ammonia as N (kg/yr)', 'Construction cost (USD/cap/yr)', 'Operating cost (USD/cap/yr)', 'Tech construction emissions (kg CO2eq/yr)', 'Tech operating emissions (kg CO2eq/yr)',
                                                    'Construction cost (% of total)', 'Operating cost (% of total)', 'Tech construction emissions (% of total)', 'Tech operating emissions (% of total)'],
                                           columns=percentile_str)
        user_interface_excreta_dataframe.to_excel(writer, sheet_name='user_interface', startcol=0)

    elif not (np.isnan(user_interface_liquid_outputs[0,0]) or np.isnan(user_interface_solid_outputs[0,0])):
        user_interface_liquid_out = np.percentile(user_interface_liquid_outputs[:,0], percentiles)
        user_interface_liquid_dry_out = np.percentile(user_interface_liquid_outputs[:,1], percentiles)
        user_interface_liquid_N_out = np.percentile(user_interface_liquid_outputs[:,2], percentiles)
        user_interface_liquid_P_out = np.percentile(user_interface_liquid_outputs[:,3], percentiles)
        user_interface_liquid_K_out = np.percentile(user_interface_liquid_outputs[:,4], percentiles)
        user_interface_liquid_Mg_out = np.percentile(user_interface_liquid_outputs[:,5], percentiles)
        user_interface_liquid_Ca_out = np.percentile(user_interface_liquid_outputs[:,6], percentiles)
        user_interface_liquid_energy_out = np.percentile(user_interface_liquid_outputs[:,7], percentiles)
        user_interface_liquid_ammonia_N_out = np.percentile(user_interface_liquid_outputs[:,8], percentiles)
        user_interface_liquid_dataframe = pd.DataFrame([user_interface_liquid_out, user_interface_liquid_dry_out, user_interface_liquid_N_out,
                                           user_interface_liquid_P_out, user_interface_liquid_K_out, user_interface_liquid_Mg_out,
                                           user_interface_liquid_Ca_out, user_interface_liquid_energy_out, user_interface_liquid_ammonia_N_out],
                                           index=['liquid (kg/yr)', 'liquid dry matter (kg/yr)', 'liquid N (kg/yr)',
                                                    'liquid P (kg/yr)', 'liquid K (kg/yr)', 'liquid Mg (kg/yr)',
                                                    'liquid Ca (kg/yr)', 'liquid energy (kJ/yr)', 'liquid ammonia as N (kg/yr)'],
                                           columns=percentile_str)
        user_interface_liquid_dataframe.to_excel(writer, sheet_name='user_interface', startcol=0)

        user_interface_solid_out = np.percentile(user_interface_solid_outputs[:,0], percentiles)
        user_interface_solid_dry_out = np.percentile(user_interface_solid_outputs[:,1], percentiles)
        user_interface_solid_N_out = np.percentile(user_interface_solid_outputs[:,2], percentiles)
        user_interface_solid_P_out = np.percentile(user_interface_solid_outputs[:,3], percentiles)
        user_interface_solid_K_out = np.percentile(user_interface_solid_outputs[:,4], percentiles)
        user_interface_solid_Mg_out = np.percentile(user_interface_solid_outputs[:,5], percentiles)
        user_interface_solid_Ca_out = np.percentile(user_interface_solid_outputs[:,6], percentiles)
        user_interface_solid_energy_out = np.percentile(user_interface_solid_outputs[:,7], percentiles)
        user_interface_solid_ammonia_N_out = np.percentile(user_interface_solid_outputs[:,8], percentiles)

        user_interface_solid_construction_cost = np.percentile(construction_cost[:,0], percentiles)
        user_interface_solid_operating_cost = np.percentile(operating_cost[:,0], percentiles)

        user_interface_solid_tech_construction_emissions = np.percentile(tech_construction_emissions[:,0], percentiles)
        user_interface_solid_tech_operating_emissions = np.percentile(tech_operating_emissions[:,0], percentiles)

        user_interface_solid_construction_cost_percent = np.percentile(construction_cost[:,0]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)
        user_interface_solid_operating_cost_percent = np.percentile(operating_cost[:,0]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)

        user_interface_solid_tech_construction_emissions_percent = np.percentile(tech_construction_emissions[:,0]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)
        user_interface_solid_tech_operating_emissions_percent = np.percentile(tech_operating_emissions[:,0]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)

        user_interface_solid_dataframe = pd.DataFrame([user_interface_solid_out, user_interface_solid_dry_out, user_interface_solid_N_out,
                                           user_interface_solid_P_out, user_interface_solid_K_out, user_interface_solid_Mg_out,
                                           user_interface_solid_Ca_out, user_interface_solid_energy_out, user_interface_solid_ammonia_N_out, user_interface_solid_construction_cost, user_interface_solid_operating_cost, user_interface_solid_tech_construction_emissions, user_interface_solid_tech_operating_emissions,
                                           user_interface_solid_construction_cost_percent, user_interface_solid_operating_cost_percent, user_interface_solid_tech_construction_emissions_percent, user_interface_solid_tech_operating_emissions_percent],
                                           index=['solid (kg/yr)', 'solid dry matter (kg/yr)', 'solid N (kg/yr)',
                                                    'solid P (kg/yr)', 'solid K (kg/yr)', 'solid Mg (kg/yr)',
                                                    'solid Ca (kg/yr)', 'solid energy (kJ/yr)', 'solid ammonia as N (kg/yr)', 'Construction cost (USD/cap/yr)', 'Operating cost (USD/cap/yr)', 'Tech construction emissions (kg CO2eq/yr)', 'Tech operating emissions (kg CO2eq/yr)',
                                                    'Construction cost (% of total)', 'Operating cost (% of total)', 'Tech construction emissions (% of total)', 'Tech operating emissions (% of total)'],
                                           columns=percentile_str)
        user_interface_solid_dataframe.to_excel(writer, sheet_name='user_interface', startcol=len(percentiles)+2)

    # decentralized storage outputs
    if not np.isnan(decentralized_storage_excreta_outputs[0,0]):
        decentralized_storage_excreta_out = np.percentile(decentralized_storage_excreta_outputs[:,0], percentiles)
        decentralized_storage_excreta_dry_out = np.percentile(decentralized_storage_excreta_outputs[:,1], percentiles)
        decentralized_storage_excreta_N_out = np.percentile(decentralized_storage_excreta_outputs[:,2], percentiles)
        decentralized_storage_excreta_P_out = np.percentile(decentralized_storage_excreta_outputs[:,3], percentiles)
        decentralized_storage_excreta_K_out = np.percentile(decentralized_storage_excreta_outputs[:,4], percentiles)
        decentralized_storage_excreta_Mg_out = np.percentile(decentralized_storage_excreta_outputs[:,5], percentiles)
        decentralized_storage_excreta_Ca_out = np.percentile(decentralized_storage_excreta_outputs[:,6], percentiles)
        decentralized_storage_excreta_energy_out = np.percentile(decentralized_storage_excreta_outputs[:,7], percentiles)
        decentralized_storage_excreta_ammonia_N_out = np.percentile(decentralized_storage_excreta_outputs[:,8], percentiles)
        decentralized_storage_excreta_filling_time = np.percentile(decentralized_storage_excreta_outputs[:,9], percentiles)

        decentralized_storage_excreta_direct_emissions = np.percentile(direct_emissions[:,1], percentiles)

        decentralized_storage_excreta_direct_emissions_percent = np.percentile(direct_emissions[:,1]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)

        decentralized_storage_excreta_dataframe = pd.DataFrame([decentralized_storage_excreta_out, decentralized_storage_excreta_dry_out, decentralized_storage_excreta_N_out,
                                           decentralized_storage_excreta_P_out, decentralized_storage_excreta_K_out, decentralized_storage_excreta_Mg_out,
                                           decentralized_storage_excreta_Ca_out, decentralized_storage_excreta_energy_out, decentralized_storage_excreta_ammonia_N_out,
                                           decentralized_storage_excreta_filling_time, decentralized_storage_excreta_direct_emissions, decentralized_storage_excreta_direct_emissions_percent],
                                           index=['Excreta (kg/yr)', 'Excreta dry matter (kg/yr)', 'Excreta N (kg/yr)',
                                                    'Excreta P (kg/yr)', 'Excreta K (kg/yr)', 'Excreta Mg (kg/yr)',
                                                    'Excreta Ca (kg/yr)', 'Excreta energy (kJ/yr)', 'Excreta ammonia as N (kg/yr)', 'Pit filling time (years)', 'Direct emissions (kg CO2eq/yr)', 'Direct emissions (% of total)'],
                                           columns=percentile_str)
        decentralized_storage_excreta_dataframe.to_excel(writer, sheet_name='decentralized_storage', startcol=0)

    elif not (np.isnan(decentralized_storage_liquid_outputs[0,0]) or np.isnan(decentralized_storage_solid_outputs[0,0])):
        decentralized_storage_liquid_out = np.percentile(decentralized_storage_liquid_outputs[:,0], percentiles)
        decentralized_storage_liquid_dry_out = np.percentile(decentralized_storage_liquid_outputs[:,1], percentiles)
        decentralized_storage_liquid_N_out = np.percentile(decentralized_storage_liquid_outputs[:,2], percentiles)
        decentralized_storage_liquid_P_out = np.percentile(decentralized_storage_liquid_outputs[:,3], percentiles)
        decentralized_storage_liquid_K_out = np.percentile(decentralized_storage_liquid_outputs[:,4], percentiles)
        decentralized_storage_liquid_Mg_out = np.percentile(decentralized_storage_liquid_outputs[:,5], percentiles)
        decentralized_storage_liquid_Ca_out = np.percentile(decentralized_storage_liquid_outputs[:,6], percentiles)
        decentralized_storage_liquid_energy_out = np.percentile(decentralized_storage_liquid_outputs[:,7], percentiles)
        decentralized_storage_liquid_ammonia_N_out = np.percentile(decentralized_storage_liquid_outputs[:,8], percentiles)
        decentralized_storage_liquid_filling_time = np.percentile(decentralized_storage_liquid_outputs[:,9], percentiles)
        decentralized_storage_liquid_treatment_time = np.percentile(decentralized_storage_liquid_outputs[:,10], percentiles)
        decentralized_storage_liquid_treatment_volume = np.percentile(decentralized_storage_liquid_outputs[:,11], percentiles)
        decentralized_storage_liquid_dataframe = pd.DataFrame([decentralized_storage_liquid_out, decentralized_storage_liquid_dry_out, decentralized_storage_liquid_N_out,
                                           decentralized_storage_liquid_P_out, decentralized_storage_liquid_K_out, decentralized_storage_liquid_Mg_out,
                                           decentralized_storage_liquid_Ca_out, decentralized_storage_liquid_energy_out, decentralized_storage_liquid_ammonia_N_out,
                                           decentralized_storage_liquid_filling_time, decentralized_storage_liquid_treatment_time, decentralized_storage_liquid_treatment_volume],
                                           index=['liquid (kg/yr)', 'liquid dry matter (kg/yr)', 'liquid N (kg/yr)',
                                                    'liquid P (kg/yr)', 'liquid K (kg/yr)', 'liquid Mg (kg/yr)',
                                                    'liquid Ca (kg/yr)', 'liquid energy (kJ/yr)', 'liquid ammonia as N (kg/yr)',
                                                    'liquid container filling time (days)', 'liquid treatment time (days)', 'liquid treatment volume (L)'],
                                           columns=percentile_str)
        decentralized_storage_liquid_dataframe.to_excel(writer, sheet_name='decentralized_storage', startcol=0)

        decentralized_storage_solid_out = np.percentile(decentralized_storage_solid_outputs[:,0], percentiles)
        decentralized_storage_solid_dry_out = np.percentile(decentralized_storage_solid_outputs[:,1], percentiles)
        decentralized_storage_solid_N_out = np.percentile(decentralized_storage_solid_outputs[:,2], percentiles)
        decentralized_storage_solid_P_out = np.percentile(decentralized_storage_solid_outputs[:,3], percentiles)
        decentralized_storage_solid_K_out = np.percentile(decentralized_storage_solid_outputs[:,4], percentiles)
        decentralized_storage_solid_Mg_out = np.percentile(decentralized_storage_solid_outputs[:,5], percentiles)
        decentralized_storage_solid_Ca_out = np.percentile(decentralized_storage_solid_outputs[:,6], percentiles)
        decentralized_storage_solid_energy_out = np.percentile(decentralized_storage_solid_outputs[:,7], percentiles)
        decentralized_storage_solid_ammonia_N_out = np.percentile(decentralized_storage_solid_outputs[:,8], percentiles)
        decentralized_storage_solid_vault_volume = np.percentile(decentralized_storage_solid_outputs[:,9], percentiles)

        decentralized_storage_solid_direct_emissions = np.percentile(direct_emissions[:,1], percentiles)

        decentralized_storage_solid_direct_emissions_percent = np.percentile(direct_emissions[:,1]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)

        decentralized_storage_solid_dataframe = pd.DataFrame([decentralized_storage_solid_out, decentralized_storage_solid_dry_out, decentralized_storage_solid_N_out,
                                           decentralized_storage_solid_P_out, decentralized_storage_solid_K_out, decentralized_storage_solid_Mg_out,
                                           decentralized_storage_solid_Ca_out, decentralized_storage_solid_energy_out, decentralized_storage_solid_ammonia_N_out,
                                           decentralized_storage_solid_vault_volume, decentralized_storage_solid_direct_emissions, decentralized_storage_solid_direct_emissions_percent],
                                           index=['solid (kg/cap/yr)', 'solid dry matter (kg/cap/yr)', 'solid N (kg/cap/yr)',
                                                    'solid P (kg/cap/yr)', 'solid K (kg/cap/yr)', 'solid Mg (kg/cap/yr)',
                                                    'solid Ca (kg/cap/yr)', 'solid energy (kJ/cap/yr)', 'solid ammonia as N (kg/cap/yr)',
                                                    'solid vault volume (m3)', 'Direct emissions (kg CO2eq/yr)', 'Direct emissions (% of total)'],
                                           columns=percentile_str)
        decentralized_storage_solid_dataframe.to_excel(writer, sheet_name='decentralized_storage', startcol=len(percentiles)+2)

    # conveyance outputs
    if not np.isnan(conveyance_excreta_outputs[0,0]):
        conveyance_excreta_out = np.percentile(conveyance_excreta_outputs[:,0], percentiles)
        conveyance_excreta_dry_out = np.percentile(conveyance_excreta_outputs[:,1], percentiles)
        conveyance_excreta_N_out = np.percentile(conveyance_excreta_outputs[:,2], percentiles)
        conveyance_excreta_P_out = np.percentile(conveyance_excreta_outputs[:,3], percentiles)
        conveyance_excreta_K_out = np.percentile(conveyance_excreta_outputs[:,4], percentiles)
        conveyance_excreta_Mg_out = np.percentile(conveyance_excreta_outputs[:,5], percentiles)
        conveyance_excreta_Ca_out = np.percentile(conveyance_excreta_outputs[:,6], percentiles)
        conveyance_excreta_energy_out = np.percentile(conveyance_excreta_outputs[:,7], percentiles)
        conveyance_excreta_ammonia_N_out = np.percentile(conveyance_excreta_outputs[:,8], percentiles)

        conveyance_excreta_tech_operating_emissions = np.percentile(tech_operating_emissions[:,2], percentiles)

        conveyance_excreta_operating_cost = np.percentile(operating_cost[:,2], percentiles)

        conveyance_excreta_tech_operating_emissions_percent = np.percentile(tech_operating_emissions[:,2]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)

        conveyance_excreta_operating_cost_percent = np.percentile(operating_cost[:,2]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)

        conveyance_excreta_dataframe = pd.DataFrame([conveyance_excreta_out, conveyance_excreta_dry_out, conveyance_excreta_N_out,
                                           conveyance_excreta_P_out, conveyance_excreta_K_out, conveyance_excreta_Mg_out,
                                           conveyance_excreta_Ca_out, conveyance_excreta_energy_out, conveyance_excreta_ammonia_N_out, conveyance_excreta_tech_operating_emissions, conveyance_excreta_operating_cost, conveyance_excreta_tech_operating_emissions_percent, conveyance_excreta_operating_cost_percent],
                                           index=['Excreta (kg/yr)', 'Excreta dry matter (kg/yr)', 'Excreta N (kg/yr)',
                                                    'Excreta P (kg/yr)', 'Excreta K (kg/yr)', 'Excreta Mg (kg/yr)',
                                                    'Excreta Ca (kg/yr)', 'Excreta energy (kJ/yr)', 'Excreta ammonia as N (kg/yr)', 'Tech operating emissions (kg CO2eq/yr)', 'Operating cost (USD/cap/yr)', 'Tech operating emissions (% of total)', 'Operating cost (% of total)'],
                                           columns=percentile_str)
        conveyance_excreta_dataframe.to_excel(writer, sheet_name='conveyance', startcol=0)

    elif not (np.isnan(conveyance_liquid_outputs[0,0]) or np.isnan(conveyance_solid_outputs[0,0])):
        conveyance_liquid_out = np.percentile(conveyance_liquid_outputs[:,0], percentiles)
        conveyance_liquid_dry_out = np.percentile(conveyance_liquid_outputs[:,1], percentiles)
        conveyance_liquid_N_out = np.percentile(conveyance_liquid_outputs[:,2], percentiles)
        conveyance_liquid_P_out = np.percentile(conveyance_liquid_outputs[:,3], percentiles)
        conveyance_liquid_K_out = np.percentile(conveyance_liquid_outputs[:,4], percentiles)
        conveyance_liquid_Mg_out = np.percentile(conveyance_liquid_outputs[:,5], percentiles)
        conveyance_liquid_Ca_out = np.percentile(conveyance_liquid_outputs[:,6], percentiles)
        conveyance_liquid_energy_out = np.percentile(conveyance_liquid_outputs[:,7], percentiles)
        conveyance_liquid_ammonia_N_out = np.percentile(conveyance_liquid_outputs[:,8], percentiles)
        conveyance_liquid_dataframe = pd.DataFrame([conveyance_liquid_out, conveyance_liquid_dry_out, conveyance_liquid_N_out,
                                           conveyance_liquid_P_out, conveyance_liquid_K_out, conveyance_liquid_Mg_out,
                                           conveyance_liquid_Ca_out, conveyance_liquid_energy_out, conveyance_liquid_ammonia_N_out],
                                           index=['liquid (kg/yr)', 'liquid dry matter (kg/yr)', 'liquid N (kg/yr)',
                                                    'liquid P (kg/yr)', 'liquid K (kg/yr)', 'liquid Mg (kg/yr)',
                                                    'liquid Ca (kg/yr)', 'liquid energy (kJ/yr)', 'liquid ammonia as N (kg/yr)'],
                                           columns=percentile_str)
        conveyance_liquid_dataframe.to_excel(writer, sheet_name='conveyance', startcol=0)

        conveyance_solid_out = np.percentile(conveyance_solid_outputs[:,0], percentiles)
        conveyance_solid_dry_out = np.percentile(conveyance_solid_outputs[:,1], percentiles)
        conveyance_solid_N_out = np.percentile(conveyance_solid_outputs[:,2], percentiles)
        conveyance_solid_P_out = np.percentile(conveyance_solid_outputs[:,3], percentiles)
        conveyance_solid_K_out = np.percentile(conveyance_solid_outputs[:,4], percentiles)
        conveyance_solid_Mg_out = np.percentile(conveyance_solid_outputs[:,5], percentiles)
        conveyance_solid_Ca_out = np.percentile(conveyance_solid_outputs[:,6], percentiles)
        conveyance_solid_energy_out = np.percentile(conveyance_solid_outputs[:,7], percentiles)
        conveyance_solid_ammonia_N_out = np.percentile(conveyance_solid_outputs[:,8], percentiles)

        conveyance_solid_tech_operating_emissions = np.percentile(tech_operating_emissions[:,2], percentiles)

        conveyance_solid_operating_cost = np.percentile(operating_cost[:,2], percentiles)

        conveyance_solid_tech_operating_emissions_percent = np.percentile(tech_operating_emissions[:,2]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)

        conveyance_solid_operating_cost_percent = np.percentile(operating_cost[:,2]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)

        conveyance_solid_dataframe = pd.DataFrame([conveyance_solid_out, conveyance_solid_dry_out, conveyance_solid_N_out,
                                           conveyance_solid_P_out, conveyance_solid_K_out, conveyance_solid_Mg_out,
                                           conveyance_solid_Ca_out, conveyance_solid_energy_out, conveyance_solid_ammonia_N_out, conveyance_solid_tech_operating_emissions, conveyance_solid_operating_cost, conveyance_solid_tech_operating_emissions_percent, conveyance_solid_operating_cost_percent],
                                           index=['solid (kg/yr)', 'solid dry matter (kg/yr)', 'solid N (kg/yr)',
                                                    'solid P (kg/yr)', 'solid K (kg/yr)', 'solid Mg (kg/yr)',
                                                    'solid Ca (kg/yr)', 'solid energy (kJ/yr)', 'solid ammonia as N (kg/yr)', 'Tech operating emissions (kg CO2eq/yr)', 'Total operating cost (USD/cap/yr)', 'Tech operating emissions (% of total)', 'Operating cost (% of total)'],
                                           columns=percentile_str)
        conveyance_solid_dataframe.to_excel(writer, sheet_name='conveyance', startcol=len(percentiles)+2)

    # treatment outputs
    if not np.isnan(treatment_excreta_outputs[0,0]):
        treatment_excreta_out = np.percentile(treatment_excreta_outputs[:,0], percentiles)
        treatment_excreta_dry_out = np.percentile(treatment_excreta_outputs[:,1], percentiles)
        treatment_excreta_N_out = np.percentile(treatment_excreta_outputs[:,2], percentiles)
        treatment_excreta_P_out = np.percentile(treatment_excreta_outputs[:,3], percentiles)
        treatment_excreta_K_out = np.percentile(treatment_excreta_outputs[:,4], percentiles)
        treatment_excreta_Mg_out = np.percentile(treatment_excreta_outputs[:,5], percentiles)
        treatment_excreta_Ca_out = np.percentile(treatment_excreta_outputs[:,6], percentiles)
        treatment_excreta_energy_out = np.percentile(treatment_excreta_outputs[:,7], percentiles)
        treatment_excreta_ammonia_N_out = np.percentile(treatment_excreta_outputs[:,8], percentiles)
        treatment_excreta_biogas_out = np.percentile(biogas, percentiles)

        treatment_excreta_direct_emissions = np.percentile(direct_emissions[:,3], percentiles)
        treatment_excreta_tech_construction_emissions = np.percentile(tech_construction_emissions[:,3], percentiles)
        treatment_excreta_tech_operating_emissions = np.percentile(tech_operating_emissions[:,3], percentiles)

        treatment_excreta_construction_cost = np.percentile(construction_cost[:,3], percentiles)
        treatment_excreta_operating_cost = np.percentile(operating_cost[:,3], percentiles)

        treatment_excreta_direct_emissions_percent = np.percentile(direct_emissions[:,3]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)
        treatment_excreta_tech_construction_emissions_percent = np.percentile(tech_construction_emissions[:,3]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)
        treatment_excreta_tech_operating_emissions_percent = np.percentile(tech_operating_emissions[:,3]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)

        treatment_excreta_construction_cost_percent = np.percentile(construction_cost[:,3]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)
        treatment_excreta_operating_cost_percent = np.percentile(operating_cost[:,3]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)

        treatment_excreta_dataframe = pd.DataFrame([treatment_excreta_out, treatment_excreta_dry_out, treatment_excreta_N_out,
                                           treatment_excreta_P_out, treatment_excreta_K_out, treatment_excreta_Mg_out,
                                           treatment_excreta_Ca_out, treatment_excreta_energy_out, treatment_excreta_ammonia_N_out, treatment_excreta_biogas_out, treatment_excreta_direct_emissions, treatment_excreta_tech_construction_emissions, treatment_excreta_tech_operating_emissions, treatment_excreta_construction_cost, treatment_excreta_operating_cost,
                                           treatment_excreta_direct_emissions_percent, treatment_excreta_tech_construction_emissions_percent, treatment_excreta_tech_operating_emissions_percent, treatment_excreta_construction_cost_percent, treatment_excreta_operating_cost_percent],
                                           index=['Excreta (kg/yr)', 'Excreta dry matter (kg/yr)', 'Excreta N (kg/yr)',
                                                    'Excreta P (kg/yr)', 'Excreta K (kg/yr)', 'Excreta Mg (kg/yr)',
                                                    'Excreta Ca (kg/yr)', 'Excreta energy (kJ/yr)', 'Excreta ammonia as N (kg/yr)', 'Excreta biogas (kJ/yr)', 'Direct emissions (kg CO2eq/yr)', 'Tech construction emissions (kg CO2eq/yr)', 'Tech operating emissions (kg CO2eq/yr)', 'Construction cost (USD/cap/yr)', 'Operating cost (USD/cap/yr)',
                                                    'Direct emissions (% of total)', 'Tech construction emissions (% of total)', 'Tech operating emissions (% of total)', 'Construction cost (% of total)', 'Operating cost (% of total)'],
                                           columns=percentile_str)
        treatment_excreta_dataframe.to_excel(writer, sheet_name='treatment', startcol=0)

    elif not (np.isnan(treatment_liquid_outputs[0,0]) or np.isnan(treatment_solid_outputs[0,0])):
        treatment_liquid_out = np.percentile(treatment_liquid_outputs[:,0], percentiles)
        treatment_liquid_dry_out = np.percentile(treatment_liquid_outputs[:,1], percentiles)
        treatment_liquid_N_out = np.percentile(treatment_liquid_outputs[:,2], percentiles)
        treatment_liquid_P_out = np.percentile(treatment_liquid_outputs[:,3], percentiles)
        treatment_liquid_K_out = np.percentile(treatment_liquid_outputs[:,4], percentiles)
        treatment_liquid_Mg_out = np.percentile(treatment_liquid_outputs[:,5], percentiles)
        treatment_liquid_Ca_out = np.percentile(treatment_liquid_outputs[:,6], percentiles)
        treatment_liquid_energy_out = np.percentile(treatment_liquid_outputs[:,7], percentiles)
        treatment_liquid_ammonia_N_out = np.percentile(treatment_liquid_outputs[:,8], percentiles)
        treatment_liquid_dataframe = pd.DataFrame([treatment_liquid_out, treatment_liquid_dry_out, treatment_liquid_N_out,
                                           treatment_liquid_P_out, treatment_liquid_K_out, treatment_liquid_Mg_out,
                                           treatment_liquid_Ca_out, treatment_liquid_energy_out, treatment_liquid_ammonia_N_out],
                                           index=['liquid (kg/yr)', 'liquid dry matter (kg/yr)', 'liquid N (kg/yr)',
                                                    'liquid P (kg/yr)', 'liquid K (kg/yr)', 'liquid Mg (kg/yr)',
                                                    'liquid Ca (kg/yr)', 'liquid energy (kJ/yr)', 'liquid ammonia as N (kg/yr)'],
                                           columns=percentile_str)
        treatment_liquid_dataframe.to_excel(writer, sheet_name='treatment', startcol=0)

        treatment_solid_out = np.percentile(treatment_solid_outputs[:,0], percentiles)
        treatment_solid_dry_out = np.percentile(treatment_solid_outputs[:,1], percentiles)
        treatment_solid_N_out = np.percentile(treatment_solid_outputs[:,2], percentiles)
        treatment_solid_P_out = np.percentile(treatment_solid_outputs[:,3], percentiles)
        treatment_solid_K_out = np.percentile(treatment_solid_outputs[:,4], percentiles)
        treatment_solid_Mg_out = np.percentile(treatment_solid_outputs[:,5], percentiles)
        treatment_solid_Ca_out = np.percentile(treatment_solid_outputs[:,6], percentiles)
        treatment_solid_energy_out = np.percentile(treatment_solid_outputs[:,7], percentiles)
        treatment_solid_ammonia_N_out = np.percentile(treatment_solid_outputs[:,8], percentiles)
        treatment_solid_biogas_out = np.percentile(biogas, percentiles)

        treatment_solid_direct_emissions = np.percentile(direct_emissions[:,3], percentiles)
        treatment_solid_tech_construction_emissions = np.percentile(tech_construction_emissions[:,3], percentiles)
        treatment_solid_tech_operating_emissions = np.percentile(tech_operating_emissions[:,3], percentiles)

        treatment_solid_construction_cost = np.percentile(construction_cost[:,3], percentiles)
        treatment_solid_operating_cost = np.percentile(operating_cost[:,3], percentiles)

        treatment_solid_direct_emissions_percent = np.percentile(direct_emissions[:,3]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)
        treatment_solid_tech_construction_emissions_percent = np.percentile(tech_construction_emissions[:,3]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)
        treatment_solid_tech_operating_emissions_percent = np.percentile(tech_operating_emissions[:,3]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)

        treatment_solid_construction_cost_percent = np.percentile(construction_cost[:,3]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)
        treatment_solid_operating_cost_percent = np.percentile(operating_cost[:,3]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)

        treatment_solid_dataframe = pd.DataFrame([treatment_solid_out, treatment_solid_dry_out, treatment_solid_N_out,
                                           treatment_solid_P_out, treatment_solid_K_out, treatment_solid_Mg_out,
                                           treatment_solid_Ca_out, treatment_solid_energy_out, treatment_solid_ammonia_N_out, treatment_solid_biogas_out, treatment_solid_direct_emissions, treatment_solid_tech_construction_emissions, treatment_solid_tech_operating_emissions, treatment_solid_construction_cost, treatment_solid_operating_cost,
                                           treatment_solid_direct_emissions_percent, treatment_solid_tech_construction_emissions_percent, treatment_solid_tech_operating_emissions_percent, treatment_solid_construction_cost_percent, treatment_solid_operating_cost_percent],
                                           index=['solid (kg/yr)', 'solid dry matter (kg/yr)', 'solid N (kg/yr)',
                                                    'solid P (kg/yr)', 'solid K (kg/yr)', 'solid Mg (kg/yr)',
                                                    'solid Ca (kg/yr)', 'solid energy (kJ/yr)', 'solid ammonia as N (kg/yr)', 'solid biogas (kJ/yr)', 'Direct emissions (kg CO2eq/yr)', 'Tech construction emissions (kg CO2eq/yr)', 'Tech operating emissions (kg CO2eq/yr)', 'Construction cost (USD/cap/yr)', 'Operating cost (USD/cap/yr)',
                                                    'Direct emissions (% of total)', 'Tech construction emissions (% of total)', 'Tech operating emissions (% of total)', 'Construction cost (% of total)', 'Operating cost (% of total)'],
                                           columns=percentile_str)
        treatment_solid_dataframe.to_excel(writer, sheet_name='treatment', startcol=len(percentiles)+2)

    # reuse disposal outputs
    if not np.isnan(reuse_disposal_excreta_outputs[0,0]):
        reuse_disposal_excreta_out = np.percentile(reuse_disposal_excreta_outputs[:,0], percentiles)
        reuse_disposal_excreta_dry_out = np.percentile(reuse_disposal_excreta_outputs[:,1], percentiles)
        reuse_disposal_excreta_N_out = np.percentile(reuse_disposal_excreta_outputs[:,2], percentiles)
        reuse_disposal_excreta_P_out = np.percentile(reuse_disposal_excreta_outputs[:,3], percentiles)
        reuse_disposal_excreta_K_out = np.percentile(reuse_disposal_excreta_outputs[:,4], percentiles)
        reuse_disposal_excreta_Mg_out = np.percentile(reuse_disposal_excreta_outputs[:,5], percentiles)
        reuse_disposal_excreta_Ca_out = np.percentile(reuse_disposal_excreta_outputs[:,6], percentiles)
        reuse_disposal_excreta_energy_out = np.percentile(reuse_disposal_excreta_outputs[:,7], percentiles)
        reuse_disposal_excreta_ammonia_N_out = np.percentile(reuse_disposal_excreta_outputs[:,8], percentiles)
        reuse_disposal_excreta_biogas_out = np.percentile(biogas_output, percentiles)

        reuse_disposal_excreta_emission_offsets = np.percentile(emission_offsets[:,4], percentiles)

        reuse_disposal_excreta_income = np.percentile(income[:,4], percentiles)

        reuse_disposal_excreta_emission_offsets_percent = np.percentile(emission_offsets[:,4]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)

        reuse_disposal_excreta_income_percent = np.percentile(income[:,4]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)

        reuse_disposal_excreta_dataframe = pd.DataFrame([reuse_disposal_excreta_out, reuse_disposal_excreta_dry_out, reuse_disposal_excreta_N_out,
                                           reuse_disposal_excreta_P_out, reuse_disposal_excreta_K_out, reuse_disposal_excreta_Mg_out,
                                           reuse_disposal_excreta_Ca_out, reuse_disposal_excreta_energy_out, reuse_disposal_excreta_ammonia_N_out, reuse_disposal_excreta_biogas_out, reuse_disposal_excreta_emission_offsets, reuse_disposal_excreta_income,
                                           reuse_disposal_excreta_emission_offsets_percent, reuse_disposal_excreta_income_percent],
                                           index=['Excreta (kg/yr)', 'Excreta dry matter (kg/yr)', 'Excreta N (kg/yr)',
                                                    'Excreta P (kg/yr)', 'Excreta K (kg/yr)', 'Excreta Mg (kg/yr)',
                                                    'Excreta Ca (kg/yr)', 'Excreta energy (kJ/yr)', 'Excreta ammonia as N (kg/yr)', 'Excreta biogas (kJ/yr)', 'Emission offsets (kg CO2eq/yr)', 'Income (USD/cap/yr)',
                                                    'Emission offsets (% of total emissions)', 'Income (% of total cost)'],
                                           columns=percentile_str)
        reuse_disposal_excreta_dataframe.to_excel(writer, sheet_name='reuse_disposal', startcol=0)

        # spearman correlations with final outputs
        excreta_spearman_input = np.concatenate((reuse_disposal_excreta_outputs[:,0:9], biogas_output, (total_direct_emissions+total_tech_construction_emissions+total_tech_operating_emissions-total_emission_offsets), (total_construction_cost+total_operating_cost-total_income)), 1)
        excreta_spearman_columns = ['Excreta (kg/yr)', 'Excreta dry matter (kg/yr)', 'Excreta N (kg/yr)',
                                     'Excreta P (kg/yr)', 'Excreta K (kg/yr)', 'Excreta Mg (kg/yr)',
                                     'Excreta Ca (kg/yr)', 'Excreta energy (kJ/yr)', 'Excreta ammonia as N (kg/yr)', 'Excreta biogas (kJ/yr)', 'Total net emissions (kg CO2eq/yr)', 'Total net cost (USD/cap/yr)']

        spearman_excreta_rho = np.full((np.size(correlation_distributions, 1), np.size(excreta_spearman_input, 1)), np.nan)
        spearman_excreta_p = copy.deepcopy(spearman_excreta_rho)
        for i in range(0, np.size(correlation_distributions, 1)):
            for j in range(0, np.size(excreta_spearman_input, 1)):
                if (np.std(correlation_distributions[:,i]) > 0) and (np.std(excreta_spearman_input[:,j]) > 0):
                    spearman_excreta_rho[i,j], spearman_excreta_p[i,j] = stats.spearmanr(correlation_distributions[:,i],
                                      excreta_spearman_input[:,j])
        spearman_excreta_rho = pd.DataFrame(spearman_excreta_rho, index=correlation_parameters,
                                          columns=excreta_spearman_columns)
        spearman_excreta_rho.to_excel(writer, sheet_name='correlations_rho', startcol=0)
        spearman_excreta_p = pd.DataFrame(spearman_excreta_p, index=correlation_parameters,
                                          columns=excreta_spearman_columns)
        spearman_excreta_p.to_excel(writer, sheet_name='correlations_p', startcol=0)

    elif not (np.isnan(reuse_disposal_liquid_outputs[0,0]) or np.isnan(reuse_disposal_solid_outputs[0,0])):
        reuse_disposal_liquid_out = np.percentile(reuse_disposal_liquid_outputs[:,0], percentiles)
        reuse_disposal_liquid_dry_out = np.percentile(reuse_disposal_liquid_outputs[:,1], percentiles)
        reuse_disposal_liquid_N_out = np.percentile(reuse_disposal_liquid_outputs[:,2], percentiles)
        reuse_disposal_liquid_P_out = np.percentile(reuse_disposal_liquid_outputs[:,3], percentiles)
        reuse_disposal_liquid_K_out = np.percentile(reuse_disposal_liquid_outputs[:,4], percentiles)
        reuse_disposal_liquid_Mg_out = np.percentile(reuse_disposal_liquid_outputs[:,5], percentiles)
        reuse_disposal_liquid_Ca_out = np.percentile(reuse_disposal_liquid_outputs[:,6], percentiles)
        reuse_disposal_liquid_energy_out = np.percentile(reuse_disposal_liquid_outputs[:,7], percentiles)
        reuse_disposal_liquid_ammonia_N_out = np.percentile(reuse_disposal_liquid_outputs[:,8], percentiles)
        reuse_disposal_liquid_dataframe = pd.DataFrame([reuse_disposal_liquid_out, reuse_disposal_liquid_dry_out, reuse_disposal_liquid_N_out,
                                           reuse_disposal_liquid_P_out, reuse_disposal_liquid_K_out, reuse_disposal_liquid_Mg_out,
                                           reuse_disposal_liquid_Ca_out, reuse_disposal_liquid_energy_out, reuse_disposal_liquid_ammonia_N_out],
                                           index=['liquid (kg/yr)', 'liquid dry matter (kg/yr)', 'liquid N (kg/yr)',
                                                    'liquid P (kg/yr)', 'liquid K (kg/yr)', 'liquid Mg (kg/yr)',
                                                    'liquid Ca (kg/yr)', 'liquid energy (kJ/yr)', 'liquid ammonia as N (kg/yr)'],
                                           columns=percentile_str)
        reuse_disposal_liquid_dataframe.to_excel(writer, sheet_name='reuse_disposal', startcol=0)

        reuse_disposal_solid_out = np.percentile(reuse_disposal_solid_outputs[:,0], percentiles)
        reuse_disposal_solid_dry_out = np.percentile(reuse_disposal_solid_outputs[:,1], percentiles)
        reuse_disposal_solid_N_out = np.percentile(reuse_disposal_solid_outputs[:,2], percentiles)
        reuse_disposal_solid_P_out = np.percentile(reuse_disposal_solid_outputs[:,3], percentiles)
        reuse_disposal_solid_K_out = np.percentile(reuse_disposal_solid_outputs[:,4], percentiles)
        reuse_disposal_solid_Mg_out = np.percentile(reuse_disposal_solid_outputs[:,5], percentiles)
        reuse_disposal_solid_Ca_out = np.percentile(reuse_disposal_solid_outputs[:,6], percentiles)
        reuse_disposal_solid_energy_out = np.percentile(reuse_disposal_solid_outputs[:,7], percentiles)
        reuse_disposal_solid_ammonia_N_out = np.percentile(reuse_disposal_solid_outputs[:,8], percentiles)
        reuse_disposal_solid_biogas_out = np.percentile(biogas_output, percentiles)

        reuse_disposal_solid_emission_offsets = np.percentile(emission_offsets[:,4], percentiles)

        reuse_disposal_solid_income = np.percentile(income[:,4], percentiles)

        reuse_disposal_solid_emission_offsets_percent = np.percentile(emission_offsets[:,4]/(total_direct_emissions[:,0]+total_tech_construction_emissions[:,0]+total_tech_operating_emissions[:,0]), percentiles)

        reuse_disposal_solid_income_percent = np.percentile(income[:,4]/(total_construction_cost[:,0]+total_operating_cost[:,0]), percentiles)

        reuse_disposal_solid_dataframe = pd.DataFrame([reuse_disposal_solid_out, reuse_disposal_solid_dry_out, reuse_disposal_solid_N_out,
                                           reuse_disposal_solid_P_out, reuse_disposal_solid_K_out, reuse_disposal_solid_Mg_out,
                                           reuse_disposal_solid_Ca_out, reuse_disposal_solid_energy_out, reuse_disposal_solid_ammonia_N_out, reuse_disposal_solid_biogas_out, reuse_disposal_solid_emission_offsets, reuse_disposal_solid_income,
                                           reuse_disposal_solid_emission_offsets_percent, reuse_disposal_solid_income_percent],
                                           index=['solid (kg/yr)', 'solid dry matter (kg/yr)', 'solid N (kg/yr)',
                                                    'solid P (kg/yr)', 'solid K (kg/yr)', 'solid Mg (kg/yr)',
                                                    'solid Ca (kg/yr)', 'solid energy (kJ/yr)', 'solid ammonia as N (kg/yr)', 'solid biogas (kJ/yr)', 'Emission offsets (kg CO2eq/yr)', 'Income (USD/cap/yr)',
                                                    'Emission offsets (% of total emissions)', 'Income (% of total cost)'],
                                           columns=percentile_str)
        reuse_disposal_solid_dataframe.to_excel(writer, sheet_name='reuse_disposal', startcol=len(percentiles)+2)

        # spearman correlations with final outputs
        spearman_input = np.concatenate((reuse_disposal_solid_outputs[:,0:9]+reuse_disposal_liquid_outputs[:,0:9], biogas_output, (total_direct_emissions+total_tech_construction_emissions+total_tech_operating_emissions-total_emission_offsets), (total_construction_cost+total_operating_cost-total_income)),1)
        spearman_columns = ['mass (kg/yr)', 'mass dry (kg/yr)', 'N (kg/yr)',
                                     'P (kg/yr)', 'K (kg/yr)', 'Mg (kg/yr)',
                                     'Ca (kg/yr)', 'energy (kJ/yr)', 'ammonia as N (kg/yr)', 'biogas (kJ/yr)', 'Total net emissions (kg CO2eq/cap/yr)', 'Total net cost (USD/cap/yr)']

        spearman_rho = np.full((np.size(correlation_distributions, 1), np.size(spearman_input, 1)), np.nan)
        spearman_p = copy.deepcopy(spearman_rho)
        for i in range(0, np.size(correlation_distributions, 1)):
            for j in range(0, np.size(spearman_input, 1)):
                if (np.std(correlation_distributions[:,i]) > 0) and (np.std(spearman_input[:,j]) > 0):
                    spearman_rho[i,j], spearman_p[i,j] = stats.spearmanr(correlation_distributions[:,i],
                                      spearman_input[:,j])
        spearman_rho = pd.DataFrame(spearman_rho, index=correlation_parameters,
                                          columns=spearman_columns)
        spearman_rho.to_excel(writer, sheet_name='correlations_rho')
        spearman_p = pd.DataFrame(spearman_p, index=correlation_parameters,
                                          columns=spearman_columns)
        spearman_p.to_excel(writer, sheet_name='correlations_p')


    writer.save()
    # writer.close()