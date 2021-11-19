#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 14:57:22 2018

@author: John Trimmer; Yalin Li (minor modification)
"""

# THese functions track resources, costs, and emissions for various reuse/disposal options
#  (i) fill_cover (covering over a filled pit latrine, such that a planted tree benefits from the nutrients)
#  (ii) discharge (environmental discharge, offering no recovery)
#  (iii) crop_application (sale of liquid/solid nutrient products to farmers)
#  (iv) biogas_combustion (sale of biogas from anaerobic treatment to households as cooking fuel)

import numpy as np
import pandas as pd
import copy
import lhs

#%% Fill and Cover function
def fill_cover(inputs, parameters, correlation_distributions, correlation_parameters, n_samples):
    if parameters.tree_planted.expected == 'yes':
        # mass and nutrients are recovered; and keep extra outputs (filling time)
        outputs = copy.deepcopy(inputs)
    else:
        # nothing is recovered, but keep extra outputs (filling time)
        outputs = copy.deepcopy(inputs)
        outputs[:,0:9] = 0

    return outputs, correlation_distributions, correlation_parameters

#%% Discharge function
def discharge(inputs, parameters, correlation_distributions, correlation_parameters, n_samples):
    outputs = copy.deepcopy(inputs)

    outputs[:,0:9] = 0

    return outputs, correlation_distributions, correlation_parameters

#%% Crop Application function
def crop_application(inputs, emission_offsets, income, parameters, correlation_distributions, correlation_parameters, n_samples):
    # mass and nutrients are recovered, but not energy; this also maintains extra outputs
    outputs = copy.deepcopy(inputs)

    mass = np.reshape(inputs[:,0], (-1,1))
    mass_dry = np.reshape(inputs[:,1], (-1,1))
    N_total = np.reshape(inputs[:,2], (-1,1))
    P_total = np.reshape(inputs[:,3], (-1,1))
    K_total = np.reshape(inputs[:,4], (-1,1))
    Mg_total = np.reshape(inputs[:,5], (-1,1))
    Ca_total = np.reshape(inputs[:,6], (-1,1))
    energy = np.reshape(inputs[:,7], (-1,1))
    N_amm = np.reshape(inputs[:,8], (-1,1))

    # income from sale of sludge or liquid
    if parameters.fertilizer_sale.expected == 'yes':
        N_price, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_fertilizer_price, correlation_distributions, correlation_parameters, n_samples)
        P_price, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.P_fertilizer_price, correlation_distributions, correlation_parameters, n_samples)
        K_price, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.K_fertilizer_price, correlation_distributions, correlation_parameters, n_samples)
        discount_factor, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.sludge_fertilizer_discount_factor, correlation_distributions, correlation_parameters, n_samples)
        income_fertilizer = (((N_total/1000)*N_price) + ((P_total/1000)*P_price) + ((K_total/1000)*K_price)) * discount_factor
        income[:,4:] = income[:,4:] + income_fertilizer

    # nutrient losses during transfer to cropland
    if parameters.transfer_losses_application.expected == 'yes':
        N_amm_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_amm_loss_application, correlation_distributions, correlation_parameters, n_samples)
        N_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_loss_application, correlation_distributions, correlation_parameters, n_samples)
        P_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.P_loss_application, correlation_distributions, correlation_parameters, n_samples)
        K_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.K_loss_application, correlation_distributions, correlation_parameters, n_samples)
        Mg_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.Mg_loss_application, correlation_distributions, correlation_parameters, n_samples)
        Ca_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.Ca_loss_application, correlation_distributions, correlation_parameters, n_samples)
        C_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.C_loss_application, correlation_distributions, correlation_parameters, n_samples)
    else:
        N_amm_loss = 0
        N_loss = 0
        P_loss = 0
        K_loss = 0
        Mg_loss = 0
        Ca_loss = 0
        C_loss = 0

    # compute remaining nutrients after losses
    N_other_lost = (N_total - N_amm) * (N_loss/100)
    N_amm_lost = N_amm * ((N_amm_loss/100))
    N_total_lost = N_other_lost + N_amm_lost

    P_total_lost = P_total * (P_loss/100)
    K_total_lost = K_total * (K_loss/100)
    Mg_total_lost = Mg_total * (Mg_loss/100)
    Ca_total_lost = Ca_total * (Ca_loss/100)

    # correct if N amm loss + other N loss is > 100%
    for i in range(0, len(N_total)):
        if N_total[i] < N_total_lost[i]:
            N_total_lost[i] = N_total[i]

    N_total = N_total - N_total_lost
    N_amm = N_amm - N_amm_lost
    P_total = P_total - P_total_lost
    K_total = K_total - K_total_lost
    Mg_total = Mg_total - Mg_total_lost
    Ca_total = Ca_total - Ca_total_lost
    energy = energy * ((100 - C_loss)/100)

    # fertilizer GHG offsets
    if parameters.fertilizer_offsets.expected == 'yes':
        N_emissions, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_fertilizer_emissions, correlation_distributions, correlation_parameters, n_samples)
        P_emissions, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.P_fertilizer_emissions, correlation_distributions, correlation_parameters, n_samples)
        K_emissions, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.K_fertilizer_emissions, correlation_distributions, correlation_parameters, n_samples)
        nutrient_offsets = N_total*N_emissions + P_total*P_emissions + K_total*K_emissions
        emission_offsets[:,4:] = emission_offsets[:,4:] + nutrient_offsets

    # return outputs to output matrix
    outputs[:,0:9] = np.concatenate((mass, mass_dry, N_total, P_total, K_total,
                                     Mg_total, Ca_total, energy, N_amm), 1)

    return outputs, emission_offsets, income, correlation_distributions, correlation_parameters

#%% biogas combustion modulte
def biogas_combustion(biogas, emission_offsets, income, parameters, correlation_distributions, correlation_parameters, exchange_rate, n_samples):
    # biogas losses from fittings, etc.
    loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.biogas_loss, correlation_distributions, correlation_parameters, n_samples)
    biogas_delivered = biogas * ((100 - loss)/100)

    # income from biogas sale (USD/cap/yr)
    if parameters.biogas_sale.expected == 'yes':
        selling_price, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.LPG_selling_price, correlation_distributions, correlation_parameters, n_samples)
        specific_energy, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.LPG_specific_energy, correlation_distributions, correlation_parameters, n_samples)
        income_biogas = (biogas_delivered/1000) * ((selling_price/specific_energy)/exchange_rate)
        income[:,4:] = income[:,4:] + income_biogas

    # emissions offsets
    if parameters.biogas_offsets.expected == 'yes':
        emission_factor, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.LPG_emissions, correlation_distributions, correlation_parameters, n_samples)
        offsets_biogas = (biogas_delivered/1000) * (emission_factor/specific_energy)
        emission_offsets[:,4:] = emission_offsets[:,4:] + offsets_biogas

    # combustion efficiency
    if parameters.consider_combustion_efficiency.expected == 'yes':
        efficiency, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.efficiency_biogas, correlation_distributions, correlation_parameters, n_samples)
        biogas_output = biogas_delivered * (efficiency/100)
    else:
        biogas_output = copy.deepcopy(biogas_delivered)

    return biogas_output, emission_offsets, income, correlation_distributions, correlation_parameters

#%% Reuse or disposal function - main function
def main(input_excel_name, excreta_inputs, liquid_inputs, solid_inputs, emission_offsets, biogas, income, correlation_distributions, correlation_parameters, exchange_rate, n_samples):
    # import module parameters from input spreadsheet
    parameters = pd.DataFrame.transpose(pd.read_excel(input_excel_name, sheet_name = 'reuse_disposal').set_index('parameters'))

    # define the module(s)
    excreta_module = [parameters.excreta_module_1.expected, parameters.excreta_module_2.expected]
    liquid_module = [parameters.liquid_module_1.expected, parameters.liquid_module_2.expected]
    solid_module = [parameters.solid_module_1.expected, parameters.solid_module_2.expected]

    # create temporary variables to track excreta, solids, and liquids
    excreta_temp = copy.deepcopy(excreta_inputs)
    liquid_temp = copy.deepcopy(liquid_inputs)
    solid_temp = copy.deepcopy(solid_inputs)

    biogas_output = np.full([len(excreta_temp), 1], np.nan)

    for i in range(0, len(excreta_module)):
        if (type(excreta_module[i]) is float) and (type(liquid_module[i]) is float) and (type(solid_module[i]) is float):
            # other numerical inputs are not valid
            if (not np.isnan(excreta_module[i])):
                raise ValueError('The reuse or disposal module specified for excreta is not valid.')
            if (not np.isnan(liquid_module[i])):
                raise ValueError('The reuse or disposal module specified for liquid is not valid.')
            if (not np.isnan(solid_module[i])):
                raise ValueError('The reuse or disposal module specified for solid is not valid.')

        # otherwise, are both mixed and split stream options entered?
        elif (type(excreta_module[i]) is str) and ((type(liquid_module[i]) is str) or (type(solid_module[i]) is str)):
            raise ValueError('Modules for both the mixed and separated cases should not be evaluated simultaneously.')

        # otherwise, check mixed stream options first
        if type(excreta_module[i]) is str:
            # single pit module
            if excreta_module[i] == 'fill_cover':
                (excreta_temp, correlation_distributions,
                 correlation_parameters) = fill_cover(excreta_temp, parameters,
                                       correlation_distributions, correlation_parameters, n_samples)

            elif excreta_module[i] == 'crop_application':
                (excreta_temp, emission_offsets, income, correlation_distributions,
                 correlation_parameters) = crop_application(excreta_temp, emission_offsets, income, parameters,
                                       correlation_distributions, correlation_parameters, n_samples)

            elif excreta_module[i] == 'biogas_combustion':
                (biogas_output, emission_offsets, income, correlation_distributions,
                 correlation_parameters) = biogas_combustion(biogas, emission_offsets, income, parameters,
                                       correlation_distributions, correlation_parameters, exchange_rate, n_samples)

            elif excreta_module[i] == 'discharge':
                (excreta_temp, correlation_distributions,
                 correlation_parameters) = discharge(excreta_temp, parameters,
                                       correlation_distributions, correlation_parameters, n_samples)

            # if the excreta module input is not supported/valid
            else:
                raise ValueError('The reuse or disposal module specified for excreta is not valid.')

        if (type(liquid_module[i]) is str):
            # storage tank module
            if liquid_module[i] == 'crop_application':
                (liquid_temp, emission_offsets, income, correlation_distributions,
                 correlation_parameters) = crop_application(liquid_temp, emission_offsets, income, parameters,
                                       correlation_distributions, correlation_parameters, n_samples)

            elif liquid_module[i] == 'discharge':
                (liquid_temp, correlation_distributions,
                 correlation_parameters) = discharge(liquid_temp, parameters,
                                       correlation_distributions, correlation_parameters, n_samples)
            # if the liquid module input is not supported/valid
            else:
                raise ValueError('The reuse or disposal module specified for liquid is not valid.')

        if (type(solid_module[i]) is str):
            # dehydration vault module
            if solid_module[i] == 'crop_application':
                (solid_temp, emission_offsets, income, correlation_distributions,
                 correlation_parameters) = crop_application(solid_temp, emission_offsets, income, parameters,
                                       correlation_distributions, correlation_parameters, n_samples)

            elif solid_module[i] == 'biogas_combustion':
                (biogas_output, emission_offsets, income, correlation_distributions,
                 correlation_parameters) = biogas_combustion(biogas, emission_offsets, income, parameters,
                                       correlation_distributions, correlation_parameters, exchange_rate, n_samples)

            elif solid_module[i] == 'discharge':
                (solid_temp, correlation_distributions,
                 correlation_parameters) = discharge(solid_temp, parameters,
                                       correlation_distributions, correlation_parameters, n_samples)

            # if the solid module input is not supported/valid
            else:
                raise ValueError('The reuse or disposal module specified for solid is not valid.')

    # after iteration, set outputs equal to current values of temporary variables
    excreta_outputs = excreta_temp
    liquid_outputs = liquid_temp
    solid_outputs = solid_temp

    return excreta_outputs, liquid_outputs, solid_outputs, emission_offsets, biogas_output, income, correlation_distributions, correlation_parameters