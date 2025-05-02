#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''

This module is developed by:
    Siqi Tang <siqit@outlook.com>
    Yuyao Huang <yuyaoh2@illinois.edu>
    
This python file is used to perform uncertainty and sensitivity analysis for Enviroloo Clear Reinvented Toilet system.

'''

import os, qsdsan as qs 
from chaospy import distributions as shape
from qsdsan import Model, Metric, PowerUtility, ImpactItem
from qsdsan.utils import (
    AttrSetter,
    DictAttrSetter,
    dct_from_str,
    load_data,
    )
from exposan.utils import batch_setting_unit_params, run_uncertainty as run
from exposan import enviroloo as el
from exposan.enviroloo import (
    create_systemEL,
    #el_data_path,
    data_path,
    get_decay_k,
    get_LCA_metrics_breakdown,
    get_TEA_metrics_breakdown,
    get_LCA_metrics,
    get_TEA_metrics,
    get_normalized_CAPEX,
    get_recoveries,
    results_path,
    update_resource_recovery_settings,
    )
import numpy as np
__all__ = ('create_model', 'run_uncertainty',)


# ####################################################################################################################################################
# 
#                                         Functions for batch-making metrics and -setting parameters
#
# #####################################################################################################################################################
######################################### define function enabling the addition of metrics to the model ###############################################
# --- GWP emission conversion functions ---
GWP_CH4 = 27.2  # or your model's default #CH4_CO2_EQ
GWP_N2O = 273.0  # update these if different in your setup #N2O_CO2_eq
 

# def calc_CH4_emissions_from_unit(unit, CH4_EF):
#     COD_in = unit.ins[0].COD * unit.ins[0].F_vol / 1e3 + unit.ins[1].COD * unit.ins[1].F_vol / 1e3 # kg/day
#     COD_out = unit.outs[0].COD * unit.outs[0].F_vol / 1e3
#     COD_removed = COD_in - COD_out
#     COD_removed = max(COD_in - COD_out, 0)
#     return COD_removed * CH4_EF * GWP_CH4  # kg CO2-eq/day

# def calc_N2O_emissions_from_unit(unit, N2O_EF):
#     N_in = unit.ins[0].TN * unit.ins[0].F_vol / 1e3 + unit.ins[1].TN * unit.ins[0].F_vol / 1e3 # kg/day
#     N_out = unit.outs[0].TN * unit.outs[0].F_vol / 1e3
#     N_removed = max(N_in - N_out, 0)
#     return N_removed * N2O_EF * GWP_N2O  # kg CO2-eq/day

def calc_CH4_emissions_from_unit(unit, CH4_EF):
    COD_in = sum([s.COD * s.F_vol for s in unit.ins if s and hasattr(s, 'COD')]) / 1e3  # kg/day
    COD_out = unit.outs[0].COD * unit.outs[0].F_vol / 1e3
    COD_removed = max(COD_in - COD_out, 0)
    return COD_removed * CH4_EF * GWP_CH4  # kg CO2-eq/day

def calc_N2O_emissions_from_unit(unit, N2O_EF):
    N_in = sum([s.TN * s.F_vol for s in unit.ins if s and hasattr(s, 'TN')]) / 1e3  # kg/day
    N_out = unit.outs[0].TN * unit.outs[0].F_vol / 1e3
    N_removed = max(N_in - N_out, 0)
    return N_removed * N2O_EF * GWP_N2O  # kg CO2-eq/day


def add_metrics(model):
    el._load_lca_data() # load LCA data for EL system
    system = model.system
    # Recoveries TODO: add recoveries to system?
    # funcs = get_recoveries(system)
    # metrics = []
    #     Metric('Total N', funcs[0], '% N', 'N recovery'), #here the order of funcs in line with the function get_recoveries in _init_.py
    #     Metric('Total P', funcs[1], '% P', 'P recovery'),
    #     Metric('Total K', funcs[2], '% K', 'K recovery'),
    # ]
    # Net cost of the EL system in TEA
    # metrics = ([
    #     Metric('Annualized CAPEX', get_TEA_metrics(system)[1], f'{qs.currency}/cap/yr', 'TEA results'),
    #     Metric('Annual electricity consumption', get_TEA_metrics(system)[2], 'kWh/cap/yr', 'TEA results'),
    #     Metric('Annual labor cost', get_TEA_metrics(system)[3], f'{qs.currency}/cap/yr', 'TEA results'),
    #     Metric('OPEX excluding labor/energy', get_TEA_metrics(system)[4], f'{qs.currency}/cap/yr', 'TEA results'),
    #     Metric('Revenue', get_TEA_metrics(system)[5], f'{qs.currency}/cap/yr', 'TEA results')
        
    #     #Metric('Annual net cost', get_TEA_metrics_breakdown(system), f'{qs.currency}/cap/yr', 'TEA results'),
    #     ])
    # Net emissions of the EL system in LCA
    funcs = get_LCA_metrics(system)  # extract LCA metrics from the EL system's LCA results
    cat = 'LCA results'  # assign the same index to all LCA metrics
    metrics = [
        Metric('Annualized CAPEX', get_TEA_metrics(system)[1], f'{qs.currency}/cap/yr', 'TEA results'),
        Metric('Annual electricity consumption', get_TEA_metrics(system)[2], 'kWh/cap/yr', 'TEA results'),
        Metric('Annual labor cost', get_TEA_metrics(system)[3], f'{qs.currency}/cap/yr', 'TEA results'),
        Metric('OPEX excluding labor/energy', get_TEA_metrics(system)[4], f'{qs.currency}/cap/yr', 'TEA results'),
        Metric('Revenue', get_TEA_metrics(system)[5], f'{qs.currency}/cap/yr', 'TEA results'),
        Metric('GlobalWarming', funcs[0], 'kg CO2-eq/cap/yr', 'LCA results'),
        Metric('H_Ecosystems', funcs[1], 'points/cap/yr', 'LCA results'),
        Metric('H_Health', funcs[2], 'points/cap/yr', 'LCA results'),
        Metric('H_Resources', funcs[3], 'points/cap/yr', 'LCA results'),
    ]

        # Load CH4 and N2O emission factors from unit TSVs
        # Strip all index whitespace first
    CT_data.index = CT_data.index.str.strip()
    AnoxicTank_data.index = AnoxicTank_data.index.str.strip()
    AerobicTank_data.index = AerobicTank_data.index.str.strip()
    MembTank_data.index = MembTank_data.index.str.strip()

    
    # Now try fetching values robustly
    try:
        CH4_EF_CT = CT_data.loc['EL_CT_methane_yield', 'expected']
    except KeyError:
        print('[ERROR] Could not find CH4 EF for CT.')
        CH4_EF_CT = 0  # or raise / default fallback
    
    try:
        CH4_EF_A1 = AnoxicTank_data.loc['EL_anoT_methane_yield', 'expected']
    except KeyError:
        print('[ERROR] Could not find CH4 EF for A1.')
        CH4_EF_A1 = 0
    
    try:
        N2O_EF_O1 = AerobicTank_data.loc['N2O_EF_decay', 'expected']
        N2O_EF_B1 = MembTank_data.loc['N2O_EF_decay', 'expected']
    except KeyError as e:
        print(f'[ERROR] Could not find N2O EF:', e)
        N2O_EF_O1 = N2O_EF_B1 = N2O_EF_CWT = 0


    metrics.extend([
        Metric('CH4 from CT', lambda: calc_CH4_emissions_from_unit(system.flowsheet.unit.CT, CH4_EF=CH4_EF_CT)*365/el.ppl, 'kg CO2-eq/cap/yr', 'CH4 emissions'),
        Metric('CH4 from A1', lambda: calc_CH4_emissions_from_unit(system.flowsheet.unit.A1, CH4_EF=CH4_EF_A1)*365/el.ppl, 'kg CO2-eq/cap/yr', 'CH4 emissions'),
        Metric('N2O from O1', lambda: calc_N2O_emissions_from_unit(system.flowsheet.unit.O1, N2O_EF=N2O_EF_O1)*365/el.ppl, 'kg CO2-eq/cap/yr', 'N2O emissions'),
        Metric('N2O from B1', lambda: calc_N2O_emissions_from_unit(system.flowsheet.unit.B1, N2O_EF=N2O_EF_B1)*365/el.ppl, 'kg CO2-eq/cap/yr', 'N2O emissions'),
        
    ])


    model.metrics = metrics

    
# %%
# load data for the EL model
su_data_path = os.path.join(data_path, 'sanunit_data')
el_su_data_path = os.path.join(data_path, 'units_data')

######################## load parameters related to separate units of the EL system for uncertainty and sensitivity analysis ############################
def load_el_su_data(file_name):
    if file_name.startswith('_EL'):
        return load_data(os.path.join(el_su_data_path, file_name))
    return load_data(os.path.join(su_data_path, file_name))

#TODO: standardize naming conventions
CT_data = load_el_su_data('_EL_CT.tsv')
PC_data = load_el_su_data('_EL_PC.tsv')
AnoxicTank_data = load_el_su_data('_EL_Anoxic.tsv')
AerobicTank_data = load_el_su_data('_EL_Aerobic.tsv')
MembTank_data = load_el_su_data('_EL_MBR.tsv')
ClearWaterTank_data = load_el_su_data('_EL_CWT.tsv')
system_data = load_el_su_data('_EL_system.tsv')
photovoltaic_wind_data = load_el_su_data('_EL_photovoltaic_wind.tsv')
housing_data = load_el_su_data('_EL_housing.tsv')

############################################## define parameters of interest for EL system ################################################################
def add_parameters(model, unit_dct, country_specific=False):
    sys = model.system
    sys_stream = sys.flowsheet.stream
    param = model.parameter
    price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct  = update_resource_recovery_settings()
    

    # if specific country is not selected, add these parameters below
    if not country_specific:
        # Price ratio
        # Just want to have this parameter so that can be used in other analyses,
        # set the distribution to be a really tight one
        # b = 1
        # D = shape.Uniform(lower=b-(10**(-6)), upper=b+(10**(-6)))
        # @param(name='Price ratio', element=excretion_unit, kind='cost', units='-',
        #        baseline=b, distribution=D)
        # def set_price_ratio(i):
        #     el.price_ratio = i
        #     for u in sys.units:
        #         if hasattr(u, 'price_ratio'):
        #             u.price_ratio = i

        # price ratio
        #old_price_dct = price_dct.copy()
         # need change (??)
        #b = 1
        #D = shape.Uniform(lower = b - (10**(-6)), upper = b + (10**(-6)))
        #item_ref = {
            #'Concrete': 'Concrete',
            #'Steel': 'Steel',
        #}
        #stream_ref = {
            #'ClearWater': 'ClearWater',
            #'NaOH': 'NaOH',
            #'NaClO': 'NaClO',
            #'O3': 'O3',
            #'PAC': 'PAC',
            #'HDPE': 'HDPE',
        #}
        #@param(name = 'Price ratio', element = excretion_unit, kind = 'cost', units = '-',
               #baseline = b, distribution = D)
        #def set_price_ratio(i):
            #el.price_ratio = i
            #for obj_name in (*item_ref.keys(), *stream_ref.keys()):
                #old_price = old_price_dct[obj_name]
                #new_price = old_price * i
                #if obj_name in item_ref.keys():
                    #ImpactItem.get_item(item_ref[obj_name]).price = new_price
                #else:
                    #getattr(sys_stream, stream_ref[obj_name]).price = new_price
            #for u in sys.units:
                #if hasattr(u, 'price_ratio'):
                    #u.price_ratio = i
        

        
        b = el.operator_daily_wage # operator_daily_wage is constant, needing initialization in _init_.py
        D = shape.Triangle(lower = (14.55), midpoint = b, upper = (43.68))
        @param(name = 'Operator daily wage', element = 'TEA', kind = 'cost', units = 'USD/d',
            baseline = b, distribution = D)
        def set_operator_daily_wage(i):
            sys._TEA.annual_labor = i * 3 * 365 # need change TODO

        # Construction labor wage
        # b = el.const_daily_wage   # const_daily_wage is constant, needing initialization in _init_.py
        # D = shape.Triangle(lower = (b * 0.5), midpoint = b, upper = (b * 1.5))
        # @param(name = 'Construction daily wage', element = 'TEA', kind = 'cost', units = 'USD/d',
        #     baseline = b, distribution = D)
        # def set_const_daily_wage(i):
        #     for u in sys.units:
        #         if isinstance(u, qs.sanunits.EL_Housing): break
        #         u.const_daily_wage = i
        
        if el.INCLUDED_RESOURCE_RECOVERY:
            # N fertilizer price
            b = 1.507
            D = shape.Uniform(lower = b * 0.8, upper = b * 1.2)
            @param(name = 'N fertilizer price', element = 'TEA', kind = 'isolated', units = 'USD/kg N',
                   baseline = b, distribution = D)
            def set_N_fert_price(i):
                price_dct['N'] = sys_stream.liq_N.price = sys_stream.sol_N.price = i * el.price_factor
            
            # P fertilizer price
            b = 3.983
            D = shape.Uniform(lower = b * 0.8, upper = b * 1.2)
            @param(name = 'P fertilizer price', element = 'TEA', kind = 'isolated', units = 'USD/kg P',
                   baseline = b, distribution = D)
            def set_P_fert_price(i):
                price_dct['P'] = sys_stream.liq_P.price = sys_stream.sol_P.price = i * el.price_factor
            
            # K fertilizer price
            b = 1.333
            D = shape.Uniform(lower = b * 0.8, upper = b * 1.2)
            @param(name = 'K fertilizer price', element = 'TEA', kind = 'isolated', units = 'USD/kg K',
                   baseline = b, distribution = D)
            def set_K_fert_price(i):
                price_dct['K'] = sys_stream.liq_K.price = sys_stream.sol_K.price = i * el.price_factor
            
            # Ammonium fertilizer price
            b = 1.507 * (14/18)
            D = shape.Uniform(lower = b * 0.8, upper = b * 1.2)
            @param(name = 'Ammonium fertilizer price', element = 'TEA', kind = 'isolated', units = 'USD/kg N',
                   baseline = b, distribution = D)
            def set_NH4_fert_price(i):
                price_dct['ammonium'] = sys_stream.liq_NH4.price = sys_stream.sol_NH4.price = i * el.price_factor
            
            # Struvite fertilizer price
            b = 3.983 * (31 / 245) # need check up to now
            D = shape.Uniform(lower = b * 0.8, upper = b * 1.2)
            @param(name = 'Struvite fertilizer price', element = 'TEA', kind = 'isolated', units = 'USD/kg P',
                   baseline = b, distribution = D)
            def set_struvite_fert_price(i):
                price_dct['struvite'] = sys_stream.liq_struvite.price = sys_stream.sol_struvite.price = i * el.price_factor
        
        # # Electricity price, set to zero as renewable assumption.
        # b = price_dct['Electricity']
        # D = shape.Uniform(lower = 0.00, upper = 0.00)
        # @param(name = 'Electricity price', element = 'TEA', kind = 'isolated', units = 'USD/kWh',
        #        baseline = b, distribution = D)
        # def set_electricity_price(i):
        #     PowerUtility.price = i

        # # Electricity GWP, set to zero as renewable assumption
        # b = GWP_dct['Electricity']
        # D = shape.Uniform(lower = b * 0.9, upper = b * 1.1)
        # @param(name = 'Electricity CF', element = 'LCA', kind = 'isolated', units = 'kg CO2-eq/kWh',
        #        baseline = b, distribution = D)
        # def set_electricity_CF(i):
        #     GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i
        
        # # Electricity H_Ecosystems, if applicable
        # b = H_Ecosystems_dct['Electricity']
        # D = shape.Uniform(lower = b * 0.9, upper = b * 1.1)
        # @param(name = 'Electricity Ecosystems CF', element = 'LCA', kind = 'isolated', units = 'points/kWh',
        #        baseline = b, distribution = D)
        # def set_electricity_ecosystems_CF(i):
        #     H_Ecosystems_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Ecosystems'] = i
        
        # # Electricity H_Health, if applicable
        # b = H_Health_dct['Electricity']
        # D = shape.Uniform(lower = b * 0.9, upper = b * 1.1)
        # @param(name = 'Electricity Health CF', element = 'LCA', kind = 'isolated', units = 'points/kWh',
        #        baseline = b, distribution = D)
        # def set_electricity_health_CF(i):
        #     H_Health_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Health'] = i
        
        # # Electricity H_Resources, if applicable
        # b = H_Resources_dct['Electricity']
        # D = shape.Uniform(lower = b * 0.9, upper = b * 1.1)
        # @param(name = 'Electricity Resources CF', element = 'LCA', kind = 'isolated', units = 'points/kWh',
        #        baseline = b, distribution = D)
        # def set_electricity_resources_CF(i):
        #     H_Resources_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Resources'] = i
    
    ############################# Specific Units having parameters engaged with uncertainty and sensitivity analysis ###############################

    CT_unit = unit_dct['Collection_Tank']
    if CT_unit: 
        batch_setting_unit_params(CT_data, model, CT_unit)
    
    # In Primary Clarifier Tank
    PC_unit = unit_dct['PrimaryClarifierTank']
    if PC_unit: 
        batch_setting_unit_params(PC_data, model, PC_unit)
    
    # In Anoxic Tank
    AnoxT_unit = unit_dct['AnoxicTank']
    if AnoxT_unit: 
        batch_setting_unit_params(AnoxicTank_data, model, AnoxT_unit)
    
    # In Aerobic Tank
    AeroT_unit = unit_dct['AerobicTank']
    if AeroT_unit: 
        batch_setting_unit_params(AerobicTank_data, model, AeroT_unit)
    
    # In Membrane Tank
    MembT_unit = unit_dct['MembraneTank']
    if MembT_unit: 
        batch_setting_unit_params(MembTank_data, model, MembT_unit)
    
    # In Clear Water Tank
    ClearWaterT_unit = unit_dct['ClearWaterTank']
    if ClearWaterT_unit: 
        batch_setting_unit_params(ClearWaterTank_data, model, ClearWaterT_unit)
    
    # EL PV system
    PV_unit = unit_dct['PhotovoltaicWind']
    if PV_unit: 
        batch_setting_unit_params(photovoltaic_wind_data, model, PV_unit)
   
    # EL housing
    ELH_unit = unit_dct['Housing']
    if ELH_unit:
        batch_setting_unit_params(housing_data, model, ELH_unit)


    ################################################ Universal degradation parameters ##########################################################
    #TODO: system does not have emissions from run processes
    # # Max methane emission
    # toilet_unit = sys.path[1] # the first unit involving degradation
    # # toilet_unit = toilet_unit  # the first unit involving degradation
    # b = el.max_CH4_emission
    # D = shape.Triangle(lower = 0.175, midpoint = b, upper = 0.325)
    # @param(name = 'Max CH4 emission', element = toilet_unit, kind = 'coupled', units = 'g CH4/g COD',
    #        baseline = b, distribution = D)
    # def set_max_CH4_emission(i):
    #     el.max_CH4_emission = i
    #     for unit in sys.units:
    #         if hasattr(unit, 'max_CH4_emission'):
    #             setattr(unit, 'max_CH4_emission', i)
    
    # # time to full degradation
    # b =el.tau_deg
    # D = shape.Uniform(lower = 1, upper = 3)
    # @param(name = 'Time to full degradation', element = toilet_unit, kind = 'coupled', units = 'yr',
    #        baseline = b, distribution = D)
    # def set_tau_deg(i):
    #     el.tau_deg = i
    #     k = get_decay_k(i, el.log_deg)
    #     for unit in sys.units:
    #         if hasattr(unit, 'decay_k_COD'):
    #             setattr(unit, 'decay_k_COD', k)
    #         if hasattr(unit, 'decay_k_N'):
    #             setattr(unit, 'decay_k_N', k)

    # # reduction at full degradation
    # b = el.log_deg
    # D = shape.Uniform(lower = 2, upper = 4)
    # @param(name = 'Log of degradation rate', element = toilet_unit, kind = 'coupled', units = '-',
    #        baseline = b, distribution = D)
    # def set_log_deg(i):
    #     el.log_deg = i
    #     k = get_decay_k(el.tau_deg, i)
    #     for unit in sys.units:
    #         if hasattr(unit, 'decay_k_COD'):
    #             setattr(unit, 'decay_k_COD', k)
    #         if hasattr(unit, 'decay_k_N'):
    #             setattr(unit, 'decay_k_N', k)
    
    ################################################# General TEA Settings ##################################################
    # Discount rate if changing in TEA
    b = el.discount_rate
    D = shape.Uniform(lower = 0.02, upper = 0.06)
    @param(name = 'Discount rate', element = 'TEA', kind = 'isolated', units = 'fraction',
          baseline = b, distribution = D)
    def set_discount_rate(i):
        el.discount_rate = i
    
    # if resource recovery settings are updated, change price_factor for uncertainty analysis
    b = el.price_factor
    D = shape.Uniform(lower = 0.1, upper = 0.4)
    @param(name = 'Price factor', element = 'TEA', kind = 'isolated', units = '-',
          baseline = b, distribution = D)
    def set_price_factor(i):
        el.price_factor = i

    # if emptying operation is scheduled
    b = el.emptying_fee
    D = shape.Uniform(lower=0, upper=0.3)
    @param(name = 'Emptying fee', element = 'TEA', kind = 'isolated', units = 'USD',
            baseline = b, distribution = D)
    def set_emptying_fee(i):
        el.emptying_fee = i    
    
    # #PAC
    # b = price_dct['PAC']
    # D = shape.Uniform(lower = (b * 0.9), upper = (b * 1.1))
    # @param(name = 'Price of NaClO', element = 'TEA', kind = 'isolated', units = 'USD/kg PAC',
    #       baseline = (b), distribution = D)
    # def set_PAC_price(i):
    #     price_dct['PAC'] = sys_stream.PAC.price = i
    
    ###################################################### General LCA Settings ###############################################################
    # CH4 in GWP, H_Ecosystems, H_Health
    b = GWP_dct['CH4']
    D = shape.Uniform(lower = b * 0.9, upper = b * 1.1)
    @param(name = 'CH4 CF', element = 'LCA', kind = 'isolated', units = 'kg CO2-eq/kg CH4',
          baseline = b, distribution = D)
    def set_CH4_CF(i):
        GWP_dct['CH4'] = ImpactItem.get_item('CH4_item').CFs['GlobalWarming'] = i
        

    b = H_Ecosystems_dct['CH4']
    D = shape.Uniform(lower = b * 0.9, upper = b * 1.1)
    @param(name = 'CH4 Ecosystems CF', element = 'LCA', kind = 'isolated', units = 'points/kg CH4',
          baseline = b, distribution = D)
    def set_CH4_Ecosystems_CF(i):
        H_Ecosystems_dct['CH4'] = ImpactItem.get_item('CH4_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['CH4']
    D = shape.Uniform(lower = b * 0.9, upper = b * 1.1)
    @param(name = 'CH4 Health CF', element = 'LCA', kind = 'isolated', units = 'points/kg CH4',
          baseline = b, distribution = D)
    def set_CH4_Health_CF(i):
        H_Health_dct['CH4'] = ImpactItem.get_item('CH4_item').CFs['H_Health'] = i
    
    # N2O in GWP, H_Ecosystems, H_Health
    b = GWP_dct['N2O']
    D = shape.Uniform(lower = b * 0.9, upper = b * 1.1)
    @param(name = 'N2O CF', element = 'LCA', kind = 'isolated', units = 'kg CO2-eq/kg N2O',
          baseline = b, distribution = D)
    def set_N2O_CF(i):
        GWP_dct['N2O'] = ImpactItem.get_item('N2O_item').CFs['GlobalWarming'] = i
        
    
    b = H_Ecosystems_dct['N2O']
    D = shape.Uniform(lower = b * 0.9, upper = b * 1.1)
    @param(name = 'N2O Ecosystems CF', element = 'LCA', kind = 'isolated', units = 'points/kg N2O',
          baseline = b, distribution = D)
    def set_N2O_Ecosystems_CF(i):
        H_Ecosystems_dct['N2O'] = ImpactItem.get_item('N2O_item').CFs['H_Ecosystems'] = i
    
    b = H_Health_dct['N2O']
    D = shape.Uniform(lower = b * 0.9, upper = b * 1.1)
    @param(name = 'N2O Health CF', element = 'LCA', kind = 'isolated', units = 'points/kg N2O',
          baseline = b, distribution = D)
    def set_N2O_Health_CF(i):
        H_Health_dct['N2O'] = ImpactItem.get_item('N2O_item').CFs['H_Health'] = i
    
    # # PAC in GWP, H_Ecosystems, H_Health, H_Resources
    b = GWP_dct['PAC']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'PAC CF', element = 'LCA', kind = 'isolated', units = 'kg CO2-eq/kg PAC',
          baseline = b, distribution = D)
    def set_PAC_CF(i):
        GWP_dct['PAC'] = ImpactItem.get_item('PAC_item').CFs['GlobalWarming'] = i
    
    b = H_Ecosystems_dct['PAC']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'PAC Ecosystems CF', element = 'LCA', kind = 'isolated', units = 'points/kg PAC',
          baseline = b, distribution = D)
    def set_PAC_Ecosystems_CF(i):
        H_Ecosystems_dct['PAC'] = ImpactItem.get_item('PAC_item').CFs['H_Ecosystems'] = i
    
    b = H_Health_dct['PAC']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'PAC Health CF', element = 'LCA', kind = 'isolated', units = 'points/kg PAC',
          baseline = b, distribution = D)
    def set_PAC_Health_CF(i):
        H_Health_dct['PAC'] = ImpactItem.get_item('PAC_item').CFs['H_Health'] = i
    
    b = H_Resources_dct['PAC']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'PAC Resources CF', element = 'LCA', kind = 'isolated', units = 'points/kg PAC',
          baseline = b, distribution = D)
    def set_PAC_Resources_CF(i):
        H_Resources_dct['PAC'] = ImpactItem.get_item('PAC_item').CFs['H_Resources'] = i

    if el.INCLUDED_RESOURCE_RECOVERY:
        # Recovered N as alternative fertilizer, in GWP, H_Ecosystems, H_Health, H_Resources
        b = -GWP_dct['N']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'N fertilizer CF', element = 'LCA', kind = 'isolated',
               units = 'kg CO2-eq/kg N', baseline = b, distribution = D)
        def set_N_fertilizer_CF(i):
            GWP_dct['N'] = ImpactItem.get_item('N_item').CFs['GlobalWarming'] = -i
        
        b = -H_Ecosystems_dct['N']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'N fertilizer Ecosystems CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg N', baseline = b, distribution = D)
        def set_N_fertilizer_Ecosystems_CF(i):
            H_Ecosystems_dct['N'] = ImpactItem.get_item('N_item').CFs['H_Ecosystems'] = -i
        
        b = -H_Health_dct['N']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'N fertilizer Health CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg N', baseline = b, distribution = D)
        def set_N_fertilizer_Health_CF(i):
            H_Health_dct['N'] = ImpactItem.get_item('N_item').CFs['H_Health'] = -i
        
        b = -H_Resources_dct['N']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'N fertilizer Resources CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg N', baseline = b, distribution = D)
        def set_N_fertilizer_Resources_CF(i):
            H_Resources_dct['N'] = ImpactItem.get_item('N_item').CFs['H_Resources'] = -i
        
        # Recovered P as alternative fertilizer, in GWP, H_Ecosystems, H_Health, H_Resources
        b = -GWP_dct['P']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'P fertilizer CF', element = 'LCA', kind = 'isolated',
               units = 'kg CO2-eq/kg P', baseline = b, distribution = D)
        def set_P_fertilizer_CF(i):
            GWP_dct['P'] = ImpactItem.get_item('P_item').CFs['GlobalWarming'] = -i
        
        b = -H_Ecosystems_dct['P']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'P fertilizer Ecosystems CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg P', baseline = b, distribution = D)
        def set_P_fertilizer_Ecosystems_CF(i):
            H_Ecosystems_dct['P'] = ImpactItem.get_item('P_item').CFs['H_Ecosystems'] = -i
        
        b = -H_Health_dct['P']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'P fertilizer Health CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg P', baseline = b, distribution = D)
        def set_P_fertilizer_Health_CF(i):
            H_Health_dct['P'] = ImpactItem.get_item('P_item').CFs['H_Health'] = -i
        
        b = -H_Resources_dct['P']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'P fertilizer Resources CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg P', baseline = b, distribution = D)
        def set_P_fertilizer_Resources_CF(i):
            H_Resources_dct['P'] = ImpactItem.get_item('P_item').CFs['H_Resources'] = -i
        
        # Recovered K as alternative fertilizer, in GWP, H_Ecosystems, H_Health, H_Resources
        b = -GWP_dct['K']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'K fertilizer CF', element = 'LCA', kind = 'isolated',
               units = 'kg CO2-eq/kg K', baseline = b, distribution = D)
        def set_K_fertilizer_CF(i):
            GWP_dct['K'] = ImpactItem.get_item('K_item').CFs['GlobalWarming'] = -i
        
        b = -H_Ecosystems_dct['K']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'K fertilizer Ecosystems CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg K', baseline = b, distribution = D)
        def set_K_fertilizer_Ecosystems_CF(i):
            H_Ecosystems_dct['K'] = ImpactItem.get_item('K_item').CFs['H_Ecosystems'] = -i

        b = -H_Health_dct['K']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'K fertilizer Health CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg K', baseline = b, distribution = D)
        def set_K_fertilizer_Health_CF(i):
            H_Health_dct['K'] = ImpactItem.get_item('K_item').CFs['H_Health'] = -i
        
        b = -H_Resources_dct['K']
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'K fertilizer Resources CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg K', baseline = b, distribution = D)
        def set_K_fertilizer_Resources_CF(i):
            H_Resources_dct['K'] = ImpactItem.get_item('K_item').CFs['H_Resources'] = -i
        
        # Recovered Ammonium as alternative N fertilizer, in GWP, H_Ecosystems, H_Health, H_Resources
        b = 5.4 * (14 / 18) # refer to Reclaimer and Biogenic Refinery, as below
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Ammonium fertilizer CF', element = 'LCA', kind = 'isolated',
               units = 'kg CO2-eq/kg Ammonium', baseline = b, distribution = D)
        def set_Ammonium_fertilizer_CF(i):
            GWP_dct['ammonium'] = ImpactItem.get_item('ammonium_item').CFs['GlobalWarming'] = -i
        
        b = 0.0461961 * (14 / 18)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Ammonium fertilizer Ecosystems CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Ammonium', baseline = b, distribution = D)
        def set_Ammonium_fertilizer_Ecosystems_CF(i):
            H_Ecosystems_dct['ammonium'] = ImpactItem.get_item('ammonium_item').CFs['H_Ecosystems'] = -i

        b = 0.637826734 * (14 / 18)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Ammonium fertilizer Health CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Ammonium', baseline = b, distribution = D)
        def set_Ammonium_fertilizer_Health_CF(i):
            H_Health_dct['ammonium'] = ImpactItem.get_item('ammonium_item').CFs['H_Health'] = -i
        
        b = 0.259196888 * (14 / 18)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Ammonium fertilizer Resources CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Ammonium', baseline = b, distribution = D)
        def set_Ammonium_fertilizer_Resources_CF(i):
            H_Resources_dct['ammonium'] = ImpactItem.get_item('ammonium_item').CFs['H_Resources'] = -i
        
        # Recovered Struvite as alternative fertilizer, in GWP, H_Ecosystems, H_Health, H_Resources
        b = 4.9 * (30.97 / 245.41)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Struvite fertilizer CF', element = 'LCA', kind = 'isolated',
               units = 'kg CO2-eq/kg Struvite', baseline = b, distribution = D)
        def set_Struvite_fertilizer_CF(i):
            GWP_dct['struvite'] = ImpactItem.get_item('struvite_item').CFs['GlobalWarming'] = -i
        
        b = 0.093269908 * (30.97 / 245.41)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Struvite fertilizer Ecosystems CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Struvite', baseline = b, distribution = D)
        def set_Struvite_fertilizer_Ecosystems_CF(i):
            H_Ecosystems_dct['struvite'] = ImpactItem.get_item('struvite_item').CFs['H_Ecosystems'] = -i
        
        b = 1.774294425 * (30.97 / 245.41)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Struvite fertilizer Health CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Struvite', baseline = b, distribution = D)
        def set_Struvite_fertilizer_Health_CF(i):
            H_Health_dct['struvite'] = ImpactItem.get_item('struvite_item').CFs['H_Health'] = -i
        
        b = 1.084191599 * (30.97 / 245.41)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Struvite fertilizer Resources CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Struvite', baseline = b, distribution = D)
        def set_Struvite_fertilizer_Resources_CF(i):
            H_Resources_dct['struvite'] = ImpactItem.get_item('struvite_item').CFs['H_Resources'] = -i
    
    ########################################################### Add other CFs ###############################################################
    item_path = os.path.join(data_path, 'impact_items.xlsx')

    for indicator in ('GlobalWarming', 'H_Ecosystems', 'H_Health', 'H_Resources'):
        sheet_name = indicator if indicator != 'GlobalWarming' else 'GWP'   # GWP is the default sheet name
        data = load_data(item_path, sheet=sheet_name)   # load data from Excel file
        for para in data.index:
            item = ImpactItem.get_item(para)
            b = item.CFs[indicator]
            lower = float(data.loc[para]['low'])
            upper = float(data.loc[para]['high'])
            dist = data.loc[para]['distribution']
            
            if dist == 'uniform':
                D = shape.Uniform(lower=lower, upper=upper)
            elif dist == 'triangular':
                D = shape.Triangle(lower=lower, midpoint=b, upper=upper)
            elif dist == 'constant': 
                continue
            else:
                raise ValueError(f'Distribution {dist} not recognized.')
            model.parameter(name=para + f'-{indicator}',
                            setter=DictAttrSetter(item, 'CFs', indicator),
                            element='LCA',
                            kind='isolated',
                            units=f'kg CO2-eq/{item.functional_unit}',
                            baseline=b, distribution=D)

    # === Uncertainty for CH4/N2O emission factors ===
    b = CT_data.loc['EL_CT_methane_yield', 'expected']
    low = CT_data.loc['EL_CT_methane_yield', 'low']
    high = CT_data.loc['EL_CT_methane_yield', 'high']
    if high > low:
        D = shape.Uniform(lower=low, upper=high)
        @param(name='CH4 EF - CT', element='Collection_Tank', kind='isolated',
               units='m3 CH4/kg COD removed', baseline=b, distribution=D)
        def set_CH4_EF_CT(i):
            global CH4_EF_CT
            CH4_EF_CT = i
    else:
        CH4_EF_CT = b  # just assign the value if no variation

    
    b = AnoxicTank_data.loc['EL_anoT_methane_yield', 'expected']
    low = AnoxicTank_data.loc['EL_anoT_methane_yield', 'low']
    high = AnoxicTank_data.loc['EL_anoT_methane_yield', 'high']
    # D = shape.Constant(lower=low, upper=high)
    # @param(name='CH4 EF - A1', element='AnoxicTank', kind='isolated',
    #        units='m3 CH4/kg COD removed', baseline=b, distribution=D)
    # def set_CH4_EF_A1(i):
    #     global CH4_EF_A1
    #     CH4_EF_A1 = i
    if high > low:
        D = shape.Uniform(lower=low, upper=high)
        @param(name='CH4 EF - A1', element='AnoxicTank', kind='isolated',
               units='m3 CH4/kg COD removed', baseline=b, distribution=D)
        def set_CH4_EF_A1(i):
            global CH4_EF_A1
            CH4_EF_A1 = i
    else:
        CH4_EF_A1 = b  # just assign the value if no variation
    
    b = AerobicTank_data.loc['N2O_EF_decay', 'expected']
    low = AerobicTank_data.loc['N2O_EF_decay', 'low']
    high = AerobicTank_data.loc['N2O_EF_decay', 'high']
    D = shape.Uniform(lower=low, upper=high)
    @param(name='N2O EF - O1', element='AerobicTank', kind='isolated',
           units='kg N2O/kg N removed', baseline=b, distribution=D)
    def set_N2O_EF_O1(i):
        global N2O_EF_O1
        N2O_EF_O1 = i
    
    b = MembTank_data.loc['N2O_EF_decay', 'expected']
    low = MembTank_data.loc['N2O_EF_decay', 'low']
    high = MembTank_data.loc['N2O_EF_decay', 'high']
    D = shape.Uniform(lower=low, upper=high)
    @param(name='N2O EF - B1', element='MembraneTank', kind='isolated',
           units='kg N2O/kg N removed', baseline=b, distribution=D)
    def set_N2O_EF_B1(i):
        global N2O_EF_B1
        N2O_EF_B1 = i

            

#Create Model for EL system
def create_modelEL(country_specific=False, **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysEL = create_systemEL(flowsheet = flowsheet)
    unitEL = sysEL.flowsheet.unit

    modelEL = Model(sysEL, **model_kwargs)
    add_metrics(modelEL)
    
   
    unit_dctEL = { # name here needs to be aligned with those in Enviroloo_system.py (!!!)
        'Collection_Tank': unitEL.CT,
        'PrimaryClarifierTank': unitEL.PC,
        'AnoxicTank': unitEL.A1,
        'AerobicTank': unitEL.O1,
        'MembraneTank': unitEL.B1,
        'ClearWaterTank': unitEL.CWT,
        'PhotovoltaicWind': unitEL.PV,
        'Housing': unitEL.ELH
        }
    add_parameters(modelEL, unit_dctEL, country_specific)
    
    return modelEL

def create_model(model_ID='EL', country_specific=False, **model_kwargs):
    model_ID = model_ID.lower().rsplit('model')[-1].rsplit('sys')[-1].upper() # works for "modelEL"/"sysEL"/"E"
    if model_ID == 'EL': f = create_modelEL
    elif model_ID == 'E': f = create_modelEL
    elif model_ID == 'L': f = create_modelEL
    else: raise ValueError(f'`model_ID` can only be "EL", "E", or "L", not "{model_ID}".')
    model = f(country_specific, **model_kwargs)
    return model

# define runing function intializing uncertainty and sensitivity analysis, TODO: make kwargs changeable
def run_uncertainty(model, N=10, rule='L', T=2, t_step=2, mpath='', **kwargs):
    #generate sample
    sample = model.sample(N = N, rule = rule)
    
    run_model(model, sample, T=T, t_step=t_step, method='RK23', mpath=mpath)

    return

def run_model(model, sample, T=2, t_step=.1, method='BDF', 
              mpath='', tpath='', seed=None):
    model.load_samples(sample)
    t_span = (0, T)
    t_eval = np.arange(0, T+t_step, t_step)
    
    #mpath for metrics data
    mpath = mpath or os.path.join(results_path, 'enviroloo_results.xlsx')
    
    #tpath for time-series data
    if not tpath:
        folder = os.path.join(results_path, 'enviroloo_time_series')
        if not os.path.isdir(folder): os.mkdir(folder)
        tpath = os.path.join(folder, 'state.npy')
    
    model.evaluate(
        state_reset_hook='reset_cache',
        t_span=t_span,
        t_eval=t_eval,
        method=method,
        export_state_to=tpath
        )
    model.table.to_excel(mpath)