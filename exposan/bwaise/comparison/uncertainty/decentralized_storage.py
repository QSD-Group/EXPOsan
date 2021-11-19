#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 15:25:13 2017

@author: John Trimmer; Yalin Li (minor modification)
"""

# These functions track resources, costs, and emissions from various decentralized storage options
#  (i) single_pit (semi-permeable latrine pit)
#  (ii) storage_tank (urine storage tank for UDDT)
#  (iii) dehydration_vault (feces container/vault for UDDT)

import numpy as np
import pandas as pd
import lhs
import copy




#%% Single Pit function
def single_pit(excreta_inputs, direct_emissions, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, number_users):
    # separate out input variables from concatenated arrays
    excreta = np.reshape(excreta_inputs[:,0], (-1,1))
    excreta_dry = np.reshape(excreta_inputs[:,1], (-1,1))
    N_total = np.reshape(excreta_inputs[:,2], (-1,1))
    P_total = np.reshape(excreta_inputs[:,3], (-1,1))
    K_total = np.reshape(excreta_inputs[:,4], (-1,1))
    Mg_total = np.reshape(excreta_inputs[:,5], (-1,1))
    Ca_total = np.reshape(excreta_inputs[:,6], (-1,1))
    energy_total = np.reshape(excreta_inputs[:,7], (-1,1))
    N_ammonia = np.reshape(excreta_inputs[:,8], (-1,1))

    # 1. Losses
    # does infiltration occur?
    if parameters.infiltration.expected == 'yes':
        # define parameters
        N_leaching, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_leaching, correlation_distributions, correlation_parameters, n_samples)
        P_leaching, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.P_leaching, correlation_distributions, correlation_parameters, n_samples)
        K_leaching, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.K_leaching, correlation_distributions, correlation_parameters, n_samples)
    elif parameters.infiltration.expected == 'no':
        N_leaching = 0
        P_leaching = 0
        K_leaching = 0

    # compute losses
    N_infiltrated = N_total * (N_leaching / 100)
    P_infiltrated = P_total * (P_leaching / 100)
    K_infiltrated = K_total * (K_leaching / 100)

    N_total = N_total - N_infiltrated
    N_ammonia = N_ammonia - N_infiltrated
    for i in range(0, len(N_ammonia)):
        if N_ammonia[i] < 0:
            N_ammonia[i] = 0
    P_total = P_total - P_infiltrated
    K_total = K_total - K_infiltrated

    # do air emissions occur?
    pit_emptying_period, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.pit_emptying_period, correlation_distributions, correlation_parameters, n_samples)
    if parameters.air_emissions_pit.expected == 'yes':
        # parameters
        N_volatilization, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_pit_volatilization, correlation_distributions, correlation_parameters, n_samples)
        N_vol = N_total * (N_volatilization)/100
        N_total = N_total - N_vol
        N_ammonia = N_ammonia - N_vol
        for i in range(0, len(N_ammonia)):
            if N_ammonia[i] < 0:
                N_ammonia[i] = 0

        if parameters.pit_above_water_table.expected == 'yes':
            if parameters.shared.expected == 'no':
                MCF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_single_above_water, correlation_distributions, correlation_parameters, n_samples)
                N2O_EF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_single_above_water, correlation_distributions, correlation_parameters, n_samples)
            elif parameters.shared.expected == 'yes':
                MCF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_communal_above_water, correlation_distributions, correlation_parameters, n_samples)
                N2O_EF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_communal_above_water, correlation_distributions, correlation_parameters, n_samples)
        elif parameters.pit_above_water_table.expected == 'no':
            MCF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_below_water, correlation_distributions, correlation_parameters, n_samples)
            N2O_EF, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_below_water, correlation_distributions, correlation_parameters, n_samples)
        OD_max_removal_storage, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.OD_max_removal_storage, correlation_distributions, correlation_parameters, n_samples)
        # calculations
        # carbon (kg COD/yr; assume 14 kJ/g COD in wastewater)
        COD_total = (energy_total/14/1000)
        COD_degrade = COD_total*(OD_max_removal_storage/100)
        COD_after = (COD_degrade/(rate_constant*pit_emptying_period))*(1 - np.exp(-rate_constant*pit_emptying_period))
        COD_loss = COD_degrade - COD_after
        CH4_emission = COD_loss*(MCF/100)*maximum_methane_emission
        CH4eq = CH4_emission*CH4_GWP
        COD_reduction = COD_loss/COD_total
        COD_total = COD_total - COD_loss
        energy_total = COD_total*14*1000
        # nitrogen (kg N/yr; N2O expressed as kg N2O/yr)
        N_max_denit, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_max_denitrification_storage, correlation_distributions, correlation_parameters, n_samples)
        N_denit = N_total * (N_max_denit/100)
        N_after = (N_denit/(rate_constant*pit_emptying_period))*(1 - np.exp(-rate_constant*pit_emptying_period))
        N_loss = N_denit - N_after
        N2O_emission = N_total*(N_loss/N_denit)*(N2O_EF/100)*(44/28)
        for i in range(0, len(N2O_emission)):
            if N2O_emission[i] > N_loss[i]*(44/28):
                N2O_emission[i] = N_loss[i]*(44/28)
        N2Oeq = N2O_emission*N2O_GWP
        N_total = N_total - N_loss
        N_ammonia = N_ammonia - N_loss
        for i in range(0, len(N_ammonia)):
            if N_ammonia[i] < 0:
                N_ammonia[i] = 0
        # solids loss (based on COD loss)
        excreta = excreta - excreta_dry*COD_reduction
        excreta_dry = excreta_dry - (excreta_dry * COD_reduction)
    else:
        CH4eq = np.full((n_samples,1), 0)
        N2Oeq = np.full((n_samples,1), 0)

    # total accumulation (kg/yr)
    sludge_accumulation_rate, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sludge_accumulation_rate, correlation_distributions, correlation_parameters, n_samples)
    total_accumulation = sludge_accumulation_rate * number_users
    excreta_out = copy.deepcopy(excreta)
    for i in range(0, len(total_accumulation)):
        # if total accumulation is between total solids and total influent mass
        if (excreta_dry[i]*number_users[i] < total_accumulation[i]) and (total_accumulation[i] < excreta_out[i]*number_users[i]):
            excreta_out[i] = total_accumulation[i]/number_users[i]
        # if total solids is greater than total accumulation, assume all water has drained
        elif (excreta_dry[i]*number_users[i] > total_accumulation[i]):
            excreta_out[i] = excreta_dry[i]
        # otherwise (total influent mass < expected accumulation), then do nothing (no drainage)

    # 2. Pit filling time
    # define pit dimensions
    pit_depth, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.pit_depth, correlation_distributions, correlation_parameters, n_samples)
    pit_area, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.pit_area, correlation_distributions, correlation_parameters, n_samples)
    # calculate pit volume (m3)
    pit_volume = pit_depth * pit_area
    # calculate filling time (years), assuming density of contents is approximately 1000 kg/m3
    filling_time = pit_volume / (excreta_out*number_users / 1000)

    # non-ideal emptying (some latrines are emptied inappropriately, discharged directly to drainage channels or aquatic environments)
    if parameters.ideal_emptying.expected == 'no':
        # parameters
        appropriate_emptying, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.appropriate_emptying, correlation_distributions, correlation_parameters, n_samples)
        MCF_aquatic_discharge, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_aquatic_discharge, correlation_distributions, correlation_parameters, n_samples)
        N2O_EF_aquatic_discharge, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_aquatic_discharge, correlation_distributions, correlation_parameters, n_samples)
        # carbon
        CH4_emission = COD_total*(1-(appropriate_emptying/100))*(OD_max_removal_storage/100)*(MCF_aquatic_discharge/100)*maximum_methane_emission
        CH4eq = CH4eq + CH4_emission*CH4_GWP
        COD_total = COD_total*(appropriate_emptying/100)
        energy_total = COD_total*14*1000
        # nitrogen
        N2O_emission = N_total*(1-(appropriate_emptying/100))*(N2O_EF_aquatic_discharge/100)*(44/28)
        N2Oeq = N2Oeq + N2O_emission*N2O_GWP
        N_total = N_total*(appropriate_emptying/100)
        # others
        N_ammonia = N_ammonia*(appropriate_emptying/100)
        P_total = P_total*(appropriate_emptying/100)
        K_total = K_total*(appropriate_emptying/100)
        Mg_total = Mg_total*(appropriate_emptying/100)
        Ca_total = Ca_total*(appropriate_emptying/100)
        excreta_out = excreta_out*(appropriate_emptying/100)
        excreta_dry = excreta_dry*(appropriate_emptying/100)

    direct_emissions[:,1:2] = direct_emissions[:,1:2] + CH4eq + N2Oeq

    # concatenate outputs
    excreta_outputs = np.concatenate((excreta_out, excreta_dry, N_total, P_total, K_total,
                                     Mg_total, Ca_total, energy_total, N_ammonia,
                                     filling_time), 1)

    previous_storage_time = pit_emptying_period

    return excreta_outputs, direct_emissions, correlation_distributions, correlation_parameters, previous_storage_time

#%% Storage Tank function
def storage_tank(urine_inputs, direct_emissions, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, number_users, collection_period):
    # separate out input variables from concatenated arrays
    urine = np.reshape(urine_inputs[:,0], (-1,1))
    urine_dry = np.reshape(urine_inputs[:,1], (-1,1))
    N_total = np.reshape(urine_inputs[:,2], (-1,1))
    P_total = np.reshape(urine_inputs[:,3], (-1,1))
    K_total = np.reshape(urine_inputs[:,4], (-1,1))
    Mg_total = np.reshape(urine_inputs[:,5], (-1,1))
    Ca_total = np.reshape(urine_inputs[:,6], (-1,1))
    energy_total = np.reshape(urine_inputs[:,7], (-1,1))
    N_ammonia = np.reshape(urine_inputs[:,8], (-1,1))

    # 1. N losses due to ammonia volatilization (kg N/yr)
    N_volatilization, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_urine_volatilization, correlation_distributions, correlation_parameters, n_samples)
    # calculate mass of N volatilized (kg N/yr)
    N_volatilized = N_total * (N_volatilization / 100)
    # remaining total N and ammonia (kg N/yr)
    N_total = N_total - N_volatilized
    N_ammonia = N_ammonia - N_volatilized
    for i in range(0, len(N_ammonia)):
        if N_ammonia[i] < 0:
            N_ammonia[i] = 0
    # remaining dry matter and total mass
    urine_dry = urine_dry - (N_volatilized * 17 / 14)
    urine = urine - (N_volatilized * 17 / 14)

    # 2. N and P losses due to precipitation of struvite and HAP (assume as much Mg and Ca consumed as possible)
    if parameters.precipitation_losses.expected == 'yes':
        # initialize matrices to hold precipitate quantities
        struvite = np.full(np.shape(P_total), np.nan)
        HAP = np.full(np.shape(P_total), np.nan)
        # molar concentrations
        amm_M = ((N_ammonia/urine)*1000)/14.01
        P_M = ((P_total/urine)*1000)/30.97
        Mg_M = ((Mg_total/urine)*1000)/24.31
        # conditional Ksp
        pKsp, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.struvite_cond_pKsp, correlation_distributions, correlation_parameters, n_samples)
        Ksp = 10.0**(-pKsp)

        # use cond Ksp to calculate equilibrium conditions
        for i in range(0, len(P_total)):
            # coefficients of cubic equation (Ksp = (initial N - struvite)(initial P - struvite)(initial Mg - struvite))
            coeff = [1, -(Mg_M[i] + amm_M[i] + P_M[i]),
                     (Mg_M[i]*amm_M[i] + amm_M[i]*P_M[i] + Mg_M[i]*P_M[i]),
                     (Ksp[i] - Mg_M[i]*amm_M[i]*P_M[i])]
            # calculate roots of cubic
            r = np.roots(coeff)
            # identify true struvite production
            for j in range(0, len(r)):
                if np.min([Mg_M[i], amm_M[i], P_M[i]]) > r[j]:
                    struvite[i] = (r[j]*245.41/1000)*urine[i]
            if struvite[i] < 0:
                struvite[i] = 0

        # calculate nutrients in unrecoverable struvite (scaling; sludge is recoverable)
        precipitate_sludge, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.precipitate_sludge, correlation_distributions, correlation_parameters, n_samples)
        struvite = struvite * ((100 - precipitate_sludge)/100)
        P_struvite = struvite * 30.97 / 245.41
        N_struvite = struvite * 14 / 245.41
        Mg_struvite = struvite * 24.31 / 245.41

        # remaining after precipitation (and accounting for sludge fraction, which can still be recovered)
        N_total = N_total - N_struvite
        N_ammonia = N_ammonia - N_struvite
        P_total = P_total - P_struvite
        Mg_total = Mg_total - Mg_struvite

        for i in range(0, len(P_total)):
            # HAP losses (precipitates second according to Udert et al., 2003) - expect all precipitates, due to very high pKsp (57.5)
            # compare molar ratios to determine which element is completely consumed (3 P : 5 Ca)
            if (P_total[i] / (3 * 30.97)) > (Ca_total[i] / (5 * 40.08)):
                # all Ca consumed
                HAP[i] = Ca_total[i] * 502.31 / (5 * 40.08)
            else:
                # all P consumed
                HAP[i] = P_total[i] * 502.31 / (3 * 30.97)

        # calculate precipitated nutrients in HAP (in scale; sludge is recoverable)
        HAP = HAP * ((100 - precipitate_sludge)/100)
        P_HAP = HAP * (3 * 30.97) / 502.31
        Ca_HAP = HAP * (5 * 40.08) / 502.31

        # total precipitated
        precipitate = struvite + HAP

        # remaining after precipitation
        P_total = P_total - P_HAP
        Ca_total = Ca_total - Ca_HAP

        urine = urine - precipitate
        # do not include complexed water or hydroxide ions as coming from solids
        urine_dry = urine_dry - (struvite * 137.29 / 245.41) - (HAP * 485.3 / 502.31)

        # correct for any rounding errors that create negative numbers
        N_total[N_total < 0] = 0
        P_total[P_total < 0] = 0
        Mg_total[Mg_total < 0] = 0
        Ca_total[Ca_total < 0] = 0
        N_ammonia[N_ammonia < 0] = 0
    else:
        precipitate = 0

    # 3. pathogen inactivation (using Fidjeland et al., 2015 model for Ascaris egg inactivation through ammonia)
    if parameters.in_situ_treatment.expected == 'yes':
        # define parameters
        log_inactivation, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.desired_pathogen_inactivation, correlation_distributions, correlation_parameters, n_samples)
        safety_factor, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.safety_factor, correlation_distributions, correlation_parameters, n_samples)
        temperature, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.temperature, correlation_distributions, correlation_parameters, n_samples)
        pH, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.stored_urine_pH, correlation_distributions, correlation_parameters, n_samples)

        # define total ammonia concentration (mM N)
        TAN_conc = (N_ammonia * 1000000 / 14) / urine
        # define dry matter content (%)
        DM = (urine_dry / urine) * 100

        # calculate ammonia pKa
        pKa = 0.09018 + (2729.92 / (273.15 + temperature))

        # calculate fraction of total ammonia present as free ammonia (using equation from Emerson et al., 1975)
        f_NH3_Emerson = 1 / (10 ** (pKa - pH) + 1)

        # convert Emerson fraction to Pitzer fraction as required by Fidjeland et al. model
        alpha = 0.82 - 0.011 * np.sqrt(TAN_conc + 1700 * (DM / 100))
        beta = 1.17 + 0.02 * np.sqrt(TAN_conc + 1100 * (DM / 100))
        f_NH3_Pitzer = f_NH3_Emerson * (alpha + ((1 - alpha) * (f_NH3_Emerson ** beta)))

        # calculate free ammonia concentration
        NH3_conc = TAN_conc * f_NH3_Pitzer

        # calculate time (in days) to reach desired inactivation level (Fidjeland et al., 2015)
        treatment_time = (((3.2 + log_inactivation)
                            / (10 ** (-3.7 + 0.062 * temperature) * (NH3_conc ** 0.7)))
                            * 1.14 * safety_factor)

        # total required treatment volume (L)
        treatment_volume = treatment_time * (urine * number_users / 365)

    elif parameters.in_situ_treatment.expected == 'no':
        treatment_time = np.full(np.shape(urine), np.nan)
        treatment_volume = np.full(np.shape(urine), np.nan)

    # 4. storage tank filling time
    tank_volume, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.tank_volume, correlation_distributions, correlation_parameters, n_samples)
    # calculate filling time (days), assuming urine density is approximately 1 kg/L and including precipitates
    filling_time = (tank_volume / ((urine + precipitate) * number_users)) * 365

    CH4eq = np.full([n_samples, 1], 0)
    N2Oeq = np.full([n_samples, 1], 0)

    # non-ideal emptying (some latrines are emptied inappropriately, discharged directly to drainage channels or aquatic environments)
    if parameters.ideal_emptying.expected == 'no':
        # parameters
        OD_max_removal_storage, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.OD_max_removal_storage, correlation_distributions, correlation_parameters, n_samples)
        appropriate_emptying, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.appropriate_emptying, correlation_distributions, correlation_parameters, n_samples)
        MCF_aquatic_discharge, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_aquatic_discharge, correlation_distributions, correlation_parameters, n_samples)
        N2O_EF_aquatic_discharge, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_aquatic_discharge, correlation_distributions, correlation_parameters, n_samples)
        # carbon
        COD_total = (energy_total/14/1000)
        CH4_emission = COD_total*(1-(appropriate_emptying/100))*(OD_max_removal_storage/100)*(MCF_aquatic_discharge/100)*maximum_methane_emission
        CH4eq = CH4eq + CH4_emission*CH4_GWP
        COD_total = COD_total*(appropriate_emptying/100)
        energy_total = COD_total*14*1000
        # nitrogen
        N2O_emission = N_total*(1-(appropriate_emptying/100))*(N2O_EF_aquatic_discharge/100)*(44/28)
        N2Oeq = N2Oeq + N2O_emission*N2O_GWP
        N_total = N_total*(appropriate_emptying/100)
        # others
        N_ammonia = N_ammonia*(appropriate_emptying/100)
        P_total = P_total*(appropriate_emptying/100)
        K_total = K_total*(appropriate_emptying/100)
        Mg_total = Mg_total*(appropriate_emptying/100)
        Ca_total = Ca_total*(appropriate_emptying/100)
        urine = urine*(appropriate_emptying/100)
        urine_dry = urine_dry*(appropriate_emptying/100)

    direct_emissions[:,1:2] = direct_emissions[:,1:2] + CH4eq + N2Oeq

    # concatenate outputs
    urine_outputs = np.concatenate((urine, urine_dry, N_total, P_total, K_total,
                                    Mg_total, Ca_total, energy_total, N_ammonia,
                                    filling_time, treatment_time, treatment_volume), 1)

    previous_storage_time = collection_period / 365

    return urine_outputs, direct_emissions, correlation_distributions, correlation_parameters, previous_storage_time

#%% Dehydration Vault function
def dehydration_vault(feces_inputs, direct_emissions, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, number_users, collection_period):
    # separate out input variables from concatenated arrays
    feces = np.reshape(feces_inputs[:,0], (-1,1))
    feces_dry = np.reshape(feces_inputs[:,1], (-1,1))
    N_total = np.reshape(feces_inputs[:,2], (-1,1))
    P_total = np.reshape(feces_inputs[:,3], (-1,1))
    K_total = np.reshape(feces_inputs[:,4], (-1,1))
    Mg_total = np.reshape(feces_inputs[:,5], (-1,1))
    Ca_total = np.reshape(feces_inputs[:,6], (-1,1))
    energy_total = np.reshape(feces_inputs[:,7], (-1,1))
    N_ammonia = np.reshape(feces_inputs[:,8], (-1,1))

    # 1. emissions through degradation
    # convert storage time from days to years
    storage_time = collection_period / 365
    # initial moisture content
    MC_initial = (1 - feces_dry/feces) * 100
    # do air emissions occur?
    if parameters.air_emissions_vault.expected == 'yes':
        MCF_vault, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_vault, correlation_distributions, correlation_parameters, n_samples)
        N2O_EF_vault, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_vault, correlation_distributions, correlation_parameters, n_samples)
        OD_max_removal_storage, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.OD_max_removal_storage, correlation_distributions, correlation_parameters, n_samples)
        # carbon (kg COD/yr; assume 14 kJ/g COD in wastewater)
        COD_total = (energy_total/14/1000)
        COD_degrade = COD_total*(OD_max_removal_storage/100)

        #!!! Storage time shouldn't be counted twice
        # COD_after = (COD_degrade/(rate_constant*storage_time))*(np.exp(-rate_constant*storage_time) - np.exp(-rate_constant*(storage_time+storage_time)))
        COD_after = (COD_degrade/(rate_constant*storage_time))*(np.exp(-rate_constant*0) - np.exp(-rate_constant*(storage_time)))

        COD_loss = COD_degrade - COD_after
        COD_reduction = COD_loss/COD_total
        CH4_emission = COD_loss*(MCF_vault/100)*maximum_methane_emission
        CH4eq = CH4_emission*CH4_GWP
        COD_total = COD_total - COD_loss
        energy_total = COD_total*14*1000
        # nitrogen (kg N/yr; N2O expressed as kg N2O/yr)
        N_max_denit, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_max_denitrification_storage, correlation_distributions, correlation_parameters, n_samples)
        N_denit = N_total * (N_max_denit/100)
        N_after = (N_denit/(rate_constant*storage_time))*(np.exp(-rate_constant*storage_time) - np.exp(-rate_constant*(storage_time+storage_time)))
        N_loss = N_denit - N_after
        N2O_emission = N_total*(N_loss/N_denit)*(N2O_EF_vault/100)*(44/28)
        N2Oeq = N2O_emission*N2O_GWP
        N_total = N_total - N_loss
        N_ammonia = N_ammonia - N_loss
        for i in range(0, len(N_ammonia)):
            if N_ammonia[i] < 0:
                N_ammonia[i] = 0
        # solids loss (based on COD loss)
        feces = feces - feces_dry*COD_reduction
        feces_dry = feces_dry - feces_dry * COD_reduction
    else:
        CH4eq = np.full((n_samples,1), 0)
        N2Oeq = np.full((n_samples,1), 0)

    # water losses (evaporation/drying) if desiccant added
    if parameters.desiccant_added.expected == 'yes':
        MC_minimum, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.minimum_moisture_content, correlation_distributions, correlation_parameters, n_samples)
        decay_rate, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.exponential_MC_decay, correlation_distributions, correlation_parameters, n_samples)
        # express time in days to match decay rate
        time_days = storage_time * 365
        # calculate final MC as average of each addition over time
        MC_final = ((MC_initial - MC_minimum)/(decay_rate*time_days))*(1 - np.exp(-decay_rate*time_days)) + MC_minimum
        # calculate final total mass of feces
        feces = feces_dry / (1 - (MC_final / 100))

    # dehydration vault required volume (m3): assume feces density ~ 1000 kg/m3
    vault_volume = (feces * storage_time * number_users) / 1000

    # non-ideal emptying (some latrines are emptied inappropriately, discharged directly to drainage channels or aquatic environments)
    if parameters.ideal_emptying.expected == 'no':
        # parameters
        appropriate_emptying, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.appropriate_emptying, correlation_distributions, correlation_parameters, n_samples)
        MCF_aquatic_discharge, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.MCF_aquatic_discharge, correlation_distributions, correlation_parameters, n_samples)
        N2O_EF_aquatic_discharge, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N2O_EF_aquatic_discharge, correlation_distributions, correlation_parameters, n_samples)
        # carbon
        CH4_emission = COD_total*(1-(appropriate_emptying/100))*(OD_max_removal_storage/100)*(MCF_aquatic_discharge/100)*maximum_methane_emission
        CH4eq = CH4eq + CH4_emission*CH4_GWP
        COD_total = COD_total*(appropriate_emptying/100)
        energy_total = COD_total*14*1000
        # nitrogen
        N2O_emission = N_total*(1-(appropriate_emptying/100))*(N2O_EF_aquatic_discharge/100)*(44/28)
        N2Oeq = N2Oeq + N2O_emission*N2O_GWP
        N_total = N_total*(appropriate_emptying/100)
        # others
        N_ammonia = N_ammonia*(appropriate_emptying/100)
        P_total = P_total*(appropriate_emptying/100)
        K_total = K_total*(appropriate_emptying/100)
        Mg_total = Mg_total*(appropriate_emptying/100)
        Ca_total = Ca_total*(appropriate_emptying/100)
        feces = feces*(appropriate_emptying/100)
        feces_dry = feces_dry*(appropriate_emptying/100)

    direct_emissions[:,1:2] = direct_emissions[:,1:2] + CH4eq + N2Oeq

    # concatenate outputs
    feces_outputs = np.concatenate((feces, feces_dry, N_total, P_total, K_total,
                                    Mg_total, Ca_total, energy_total, N_ammonia,
                                    vault_volume), 1)

    previous_storage_time = storage_time

    return feces_outputs, direct_emissions, correlation_distributions, correlation_parameters, previous_storage_time

#%% Decentralized Collection and Storage function - main function
def main(input_excel_name, excreta_inputs, urine_inputs, feces_inputs, direct_emissions, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, number_users):
    # import module parameters from input spreadsheet
    parameters = pd.DataFrame.transpose(pd.read_excel(input_excel_name, sheet_name = 'decentralized_storage').set_index('parameters'))

    # define the module(s)
    excreta_module = parameters.mixed_excreta_module.expected
    urine_module = parameters.urine_module.expected
    feces_module = parameters.feces_module.expected

    if (type(excreta_module) is float) and (type(urine_module) is float) and (type(feces_module) is float):
        # if no modules specified, pass through (inputs = outputs)
        if np.isnan(excreta_module) and np.isnan(urine_module) and np.isnan(feces_module):
            excreta_outputs = excreta_inputs
            urine_outputs = urine_inputs
            feces_outputs = feces_inputs

        # other numerical inputs are not valid
        elif (not np.isnan(excreta_module)):
            raise ValueError('The decentralized collection and storage module specified for excreta is not valid.')
        elif (not np.isnan(urine_module)):
            raise ValueError('The decentralized collection and storage module specified for urine is not valid.')
        elif (not np.isnan(feces_module)):
            raise ValueError('The decentralized collection and storage module specified for feces is not valid.')

    # otherwise, are both mixed and split stream options entered?
    elif (type(excreta_module) is str) and ((type(urine_module) is str) or (type(feces_module) is str)):
        raise ValueError('Modules for both the mixed and separated cases should not be evaluated simultaneously.')

    # otherwise, check mixed stream options first
    elif type(excreta_module) is str:
        # single pit module
        if excreta_module == 'single_pit':
            (excreta_outputs, direct_emissions, correlation_distributions,
             correlation_parameters, previous_storage_time) = single_pit(excreta_inputs, direct_emissions, parameters,
                                   correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, number_users)

        # if the excreta module input is not supported/valid
        else:
            raise ValueError('The decentralized collection and storage module specified for excreta is not valid.')

        # no separate urine/feces streams
        urine_outputs = np.full(np.shape(urine_inputs), np.nan)
        feces_outputs = np.full(np.shape(feces_inputs), np.nan)

    elif (type(urine_module) is str) and (type(feces_module) is str):
        collection_period, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.CBS_collection_period, correlation_distributions, correlation_parameters, n_samples)
        # storage tank module
        if urine_module == 'storage_tank':
            (urine_outputs, direct_emissions, correlation_distributions,
             correlation_parameters, previous_storage_time) = storage_tank(urine_inputs, direct_emissions, parameters,
                                   correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, number_users, collection_period)

        # if the urine module input is not supported/valid
        else:
            raise ValueError('The decentralized collection and storage module specified for urine is not valid.')

        # dehydration vault module
        if feces_module == 'dehydration_vault':
            (feces_outputs, direct_emissions, correlation_distributions,
             correlation_parameters, previous_storage_time) = dehydration_vault(feces_inputs, direct_emissions, parameters,
                                   correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, number_users, collection_period)

        # if the feces module input is not supported/valid
        else:
            raise ValueError('The decentralized collection and storage module specified for feces is not valid.')

        # no mixed excreta stream
        excreta_outputs = np.full(np.shape(excreta_inputs), np.nan)

    elif (type(urine_module) is str) or (type(feces_module) is str):
        raise ValueError('For systems with separated urine and fecal streams, both urine and feces modules must be specified for decentralized collection and storage.')

    else:
        raise ValueError('The specified decentralized collection and storage modules are not valid.')

    return excreta_outputs, urine_outputs, feces_outputs, direct_emissions, correlation_distributions, correlation_parameters, previous_storage_time