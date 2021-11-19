#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 16:34:53 2019

@author: John Trimmer; Yalin Li (minor modification)
"""

# These functions track resources, costs, and emissions for various conveyance options
#  (i) tanker_truck (vacuum collection trucks for latrine emptying and transport)
#  (ii) handcart_and_truck (manual pushcarts for collecting contains from CBS systems, which are then loaded on a truck for transport to treatment plant)

import numpy as np
import pandas as pd
import lhs
import copy
from sklearn.linear_model import LinearRegression


#%% Tanker truck function
def tanker_truck(inputs, tech_operating_emissions, operating_cost, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, exchange_rate, discount_rate, number_users, previous_storage_time):
    # mass and nutrients are recovered; this also maintains extra outputs
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

    # nutrient losses during conveyance
    if parameters.losses_tanker_truck.expected == 'yes':
        N_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_loss_tanker_truck, correlation_distributions, correlation_parameters, n_samples)
        P_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.P_loss_tanker_truck, correlation_distributions, correlation_parameters, n_samples)
        K_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.K_loss_tanker_truck, correlation_distributions, correlation_parameters, n_samples)
        Mg_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.Mg_loss_tanker_truck, correlation_distributions, correlation_parameters, n_samples)
        Ca_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.Ca_loss_tanker_truck, correlation_distributions, correlation_parameters, n_samples)
        C_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.C_loss_tanker_truck, correlation_distributions, correlation_parameters, n_samples)
    else:
        N_loss = 0
        P_loss = 0
        K_loss = 0
        Mg_loss = 0
        Ca_loss = 0
        C_loss = 0

    # compute remaining nutrients after transfer losses
    N_total_lost = N_total * (N_loss/100)
    P_total_lost = P_total * (P_loss/100)
    K_total_lost = K_total * (K_loss/100)
    Mg_total_lost = Mg_total * (Mg_loss/100)
    Ca_total_lost = Ca_total * (Ca_loss/100)

    N_total = N_total - N_total_lost
    N_amm = N_amm - N_total_lost  # assume first N losses are ammonia
    for i in range(0, len(N_amm)):
        if N_amm[i] < 0:
            N_amm[i] = 0
    P_total = P_total - P_total_lost
    K_total = K_total - K_total_lost
    Mg_total = Mg_total - Mg_total_lost
    Ca_total = Ca_total - Ca_total_lost
    energy = energy * ((100 - C_loss)/100)

    # emissions
    transport_distance, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.transport_distance_tanker_truck, correlation_distributions, correlation_parameters, n_samples)
    emissions_factor, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.transport_emissions_factor_tanker_truck, correlation_distributions, correlation_parameters, n_samples)
    truck_emissions = (mass/1000)*transport_distance*emissions_factor
    tech_operating_emissions[:,2:3] = tech_operating_emissions[:,2:3] + truck_emissions

    # cost: pit emptying is a future cost, to be annualized and normalized per capita
    capacity = np.array([parameters.truck_emptying_capacity_1.expected,
                         parameters.truck_emptying_capacity_2.expected,
                         parameters.truck_emptying_capacity_3.expected,
                         parameters.truck_emptying_capacity_4.expected]).reshape(1,-1)
    base_price = np.array([parameters.truck_emptying_cost_1.expected,
                           parameters.truck_emptying_cost_2.expected,
                           parameters.truck_emptying_cost_3.expected,
                           parameters.truck_emptying_cost_4.expected]).reshape(1,-1)
    add_fee_percentage, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.truck_emptying_add_fee_percentage, correlation_distributions, correlation_parameters, n_samples)

    emptying_price = (np.repeat(base_price, n_samples, axis=0) *
                      (1 + np.repeat(add_fee_percentage/100, np.size(base_price), axis=1)))
    truck_capacity = np.repeat(capacity, n_samples, axis=0)

    ln_capacity = np.log(truck_capacity)
    ln_price = np.log(emptying_price)

    a = np.full((n_samples,1), np.nan)
    b = np.full((n_samples,1), np.nan)
    R2 = np.full((n_samples,1), np.nan)
    for i in range(0, n_samples):
        model = LinearRegression().fit(ln_capacity[i,:].reshape(-1,1), ln_price[i,:].reshape(-1,1))
        a[i] = np.exp(model.intercept_)
        b[i] = model.coef_
        R2[i] = model.score(ln_capacity[i,:].reshape(-1,1), ln_price[i,:].reshape(-1,1))

    cost_per_trip = a * ((mass*number_users*previous_storage_time/1000)**b)
    annualized_emptying_cost = (cost_per_trip/number_users/exchange_rate) * (discount_rate/(((1 + discount_rate)**previous_storage_time) - 1))

    operating_cost[:,2:3] = operating_cost[:,2:3] + annualized_emptying_cost

    # return outputs to output matrix
    outputs[:,0:9] = np.concatenate((mass, mass_dry, N_total, P_total, K_total,
                                     Mg_total, Ca_total, energy, N_amm), 1)

    return outputs, tech_operating_emissions, operating_cost, correlation_distributions, correlation_parameters

#%% Tanker truck function
def handcart_and_truck(inputs, tech_operating_emissions, operating_cost, parameters, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, exchange_rate, discount_rate, number_users, previous_storage_time):
    # mass and nutrients are recovered; this also maintains extra outputs
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

    # nutrient losses during conveyance
    if parameters.losses_handcart_and_truck.expected == 'yes':
        N_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.N_loss_handcart_and_truck, correlation_distributions, correlation_parameters, n_samples)
        P_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.P_loss_handcart_and_truck, correlation_distributions, correlation_parameters, n_samples)
        K_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.K_loss_handcart_and_truck, correlation_distributions, correlation_parameters, n_samples)
        Mg_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.Mg_loss_handcart_and_truck, correlation_distributions, correlation_parameters, n_samples)
        Ca_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.Ca_loss_handcart_and_truck, correlation_distributions, correlation_parameters, n_samples)
        C_loss, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.C_loss_handcart_and_truck, correlation_distributions, correlation_parameters, n_samples)
    else:
        N_loss = 0
        P_loss = 0
        K_loss = 0
        Mg_loss = 0
        Ca_loss = 0
        C_loss = 0

    # compute remaining nutrients after transfer losses
    N_total_lost = N_total * (N_loss/100)
    P_total_lost = P_total * (P_loss/100)
    K_total_lost = K_total * (K_loss/100)
    Mg_total_lost = Mg_total * (Mg_loss/100)
    Ca_total_lost = Ca_total * (Ca_loss/100)

    N_total = N_total - N_total_lost
    N_amm = N_amm - N_total_lost  # assume first N losses are ammonia
    for i in range(0, len(N_amm)):
        if N_amm[i] < 0:
            N_amm[i] = 0
    P_total = P_total - P_total_lost
    K_total = K_total - K_total_lost
    Mg_total = Mg_total - Mg_total_lost
    Ca_total = Ca_total - Ca_total_lost
    energy = energy * ((100 - C_loss)/100)

    # emissions
    transport_distance, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.transport_distance_CBS_truck, correlation_distributions, correlation_parameters, n_samples)
    emissions_factor, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.transport_emissions_factor_CBS_truck, correlation_distributions, correlation_parameters, n_samples)
    truck_emissions = (mass/1000)*transport_distance*emissions_factor
    tech_operating_emissions[:,2:3] = tech_operating_emissions[:,2:3] + truck_emissions

    # cost: pit emptying is a future cost, to be annualized and normalized per capita
    truck_cost_UGX, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.transport_cost_CBS_truck, correlation_distributions, correlation_parameters, n_samples)

    annualized_truck_cost = ((truck_cost_UGX)/exchange_rate)*(mass*previous_storage_time/1000) * (discount_rate/(((1 + discount_rate)**previous_storage_time) - 1))
    # annualized_truck_cost = (truck_cost_UGX/exchange_rate)*(mass/1000)

    operating_cost[:,2:3] = operating_cost[:,2:3] + annualized_truck_cost

    # return outputs to output matrix
    outputs[:,0:9] = np.concatenate((mass, mass_dry, N_total, P_total, K_total,
                                     Mg_total, Ca_total, energy, N_amm), 1)

    return outputs, tech_operating_emissions, operating_cost, correlation_distributions, correlation_parameters

#%% Conveyance function - main function
def main(input_excel_name, excreta_inputs, urine_inputs, feces_inputs, tech_operating_emissions, operating_cost, correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, exchange_rate, discount_rate, number_users, previous_storage_time):
    # import module parameters from input spreadsheet
    parameters = pd.DataFrame.transpose(pd.read_excel(input_excel_name, sheet_name = 'conveyance').set_index('parameters'))

    # define the module(s)
    excreta_module = parameters.mixed_excreta_module.expected
    urine_module = parameters.urine_module.expected
    feces_module = parameters.feces_module.expected

    excreta_outputs = np.full(np.shape(excreta_inputs), np.nan)
    feces_outputs = np.full(np.shape(feces_inputs), np.nan)
    urine_outputs = np.full(np.shape(urine_inputs), np.nan)

    if (type(excreta_module) is float) and (type(urine_module) is float) and (type(feces_module) is float):
        # if no modules specified, pass through (inputs = outputs)
        if np.isnan(excreta_module) and np.isnan(urine_module) and np.isnan(feces_module):
            excreta_outputs = excreta_inputs
            urine_outputs = urine_inputs
            feces_outputs = feces_inputs

        # other numerical inputs are not valid
        elif (not np.isnan(excreta_module)):
            raise ValueError('The conveyance module specified for excreta is not valid.')
        elif (not np.isnan(urine_module)):
            raise ValueError('The conveyance module specified for urine is not valid.')
        elif (not np.isnan(feces_module)):
            raise ValueError('The conveyance module specified for feces is not valid.')

    # otherwise, are both mixed and split stream options entered?
    elif (type(excreta_module) is str) and ((type(urine_module) is str) or (type(feces_module) is str)):
        raise ValueError('Modules for both the mixed and separated cases should not be evaluated simultaneously.')

    # otherwise, check mixed stream options first
    if type(excreta_module) is str:
        # tanker truck module
        if excreta_module == 'tanker_truck':
            (excreta_outputs, tech_operating_emissions, operating_cost, correlation_distributions,
             correlation_parameters) = tanker_truck(excreta_inputs, tech_operating_emissions, operating_cost, parameters,
                                   correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, exchange_rate, discount_rate, number_users, previous_storage_time)

        # if the excreta module input is not supported/valid
        else:
            raise ValueError('The conveyance module specified for excreta is not valid.')

    # separated streams
    if (type(urine_module) is str):
        # container based module
        if urine_module == 'handcart_and_truck':
            (urine_outputs, tech_operating_emissions, operating_cost, correlation_distributions,
             correlation_parameters) = handcart_and_truck(urine_inputs,tech_operating_emissions, operating_cost, parameters,
                                   correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, exchange_rate, discount_rate, number_users, previous_storage_time)

        # if the urine module input is not supported/valid
        else:
            raise ValueError('The conveyance module specified for urine is not valid.')

    if (type(feces_module) is str):
        # container based module
        if feces_module == 'handcart_and_truck':
            (feces_outputs, tech_operating_emissions, operating_cost, correlation_distributions,
             correlation_parameters) = handcart_and_truck(feces_inputs, tech_operating_emissions, operating_cost, parameters,
                                   correlation_distributions, correlation_parameters, n_samples, rate_constant, maximum_methane_emission, CH4_GWP, N2O_GWP, exchange_rate, discount_rate, number_users, previous_storage_time)

        # if the feces module input is not supported/valid
        else:
            raise ValueError('The conveyance module specified for feces is not valid.')

    if (urine_module == 'handcart_and_truck') or (feces_module == 'handcart_and_truck'):
        # add cost for CBS emptying to cover both urine and feces collection
        emptying_cost_CBS, correlation_distributions, correlation_parameters = lhs.lhs_distribution(parameters.emptying_cost_CBS_handcart, correlation_distributions, correlation_parameters, n_samples)

        emptying_cost_annualized = emptying_cost_CBS * (discount_rate/((1+discount_rate)**(1/365)-1))
        # emptying_cost_annualized = emptying_cost_CBS * 365

        operating_cost[:,2:3] = operating_cost[:,2:3] + emptying_cost_annualized

    return excreta_outputs, urine_outputs, feces_outputs, tech_operating_emissions, operating_cost, correlation_distributions, correlation_parameters