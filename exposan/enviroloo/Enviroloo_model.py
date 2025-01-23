#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Created by Yuyao Huang and Siqi Tang

siqit@outlook.com

This python file is used to perform uncertainty and sensitivity analysis for Enviroloo Clear Reinvented Toilet system.
'''

import os, pandas as pd, qsdsan as qs
from chaospy import distributions as shape
from qsdsan import Model, Metric, PowerUtility, ImpactItem
from qsdsan.utils import (
    AttrSetter,
    data_path,
    DictAttrSetter,
    dct_from_str,
    load_data,
    )
from exposan.utils import batch_setting_unit_params, run_uncertainty as run
from exposan import enviroloo as EL
from exposan.enviroloo import (
    create_system,
    data_path as EL_data_path,
    get_decay_k,
    get_LCA_metrics,
    get_TEA_metrics,
    get_recoveries,
    results_path,
    update_resource_recovery_settings,
    )

__all__ = ('create_model', 'run_uncertainty',)


# ####################################################################################################################################################
# 
#                                         Functions for batch-making metrics and -setting parameters
#
# #####################################################################################################################################################
######################################### define function enabling the addition of metrics to the model ###############################################
def add_metrics(model):
    EL._load_lca_data(); # load LCA data for EL system
    system = model.system
    # Recoveries
    funcs = get_recoveries(system);
    metrics = [
        Metric('Total N', funcs[0], '% N', 'N recovery'), #here the order of funcs in line with the function get_recoveries in _init_.py
        Metric('Total P', funcs[1], '% P', 'P recovery'),
        Metric('Total K', funcs[2], '% K', 'K recovery'),
    ];
    # Net cost of the EL system in TEA
    metrics.append(
        Metric('Annual net cost', get_TEA_metrics(system)[0], f'{qs.currency} / cap / yr', 'EL_TEA results'),
        );
    # Net emissions of the EL system in LCA
    funcs = get_LCA_metrics(system); # extract LCA metrics from the EL system's LCA results
    cat = 'EL_LCA results'; # assign the same index to all LCA metrics
    metrics.extend([
        Metric('GlobalWarming', funcs[0], 'kg CO2-eq/cap/yr', cat),
        Metric('H_Ecosystems', funcs[1], 'points/cap/yr', cat), 
                          # will be considered later but still placed here
        Metric('H_Health', funcs[2], 'points/cap/yr', cat),
        Metric('H_Resources', funcs[3], 'points/cap/yr', cat),
        ]);
    model.metrics = metrics
    
# load data for the EL model
su_data_path = os.path.join(data_path, 'sanunit_data'); # here 'su' has not be defined yet in the header section,
                                                        # could this calling impact data loading???
EL_su_data_path = os.path.join(su_data_path, 'el');

# locate path of the EL system data used in _EnvirolooUnits.py

######################## load parameters related to separate units of the EL system for uncertainty and sensitivity analysis ############################
def load_EL_su_data(file_name):
    if file_name.startswith('_EL'):
        return load_data(os.path.join(EL_su_data_path, file_name))
    return load_data(os.path.join(su_data_path, file_name))

''' here names in .tsv files should be consistent with those in _enviroloo_clear.py, which will be 
    modified finally.
'''
excretion_data = load_EL_su_data('_excretion.tsv'); # prepare excretion data, directly imported from qsdsan/sanunit_data/_excretion.tsv
toilet_data = load_EL_su_data('_toilet.tsv'); # prepare toilet data, directly imported from qsdsan/sanunit_data/_toilet.tsv
CT_data = load_EL_su_data('_EL_CT.tsv'); # prepare collection tank (CT) data
PC_data = load_EL_su_data('_EL_PC.tsv'); # prepare primary clarifier (PC) data
AnoxicTank_data = load_EL_su_data('_EL_anoxic.tsv'); # prepare anoxic tank data
AerobicTank_data = load_EL_su_data('_EL_aerobic.tsv'); # prepare aerobic tank data
MembTank_data = load_EL_su_data('_EL_MBR.tsv'); # prepare membrane tank data
Blower = load_EL_su_data('_EL_blower.tsv'); # prepare data of blower in membrane tank
ClearWaterTank_data = load_EL_su_data('_EL_CWT.tsv'); # prepare clear water tank data
PressureTank_data = load_EL_su_data('_EL_PT.tsv'); # prepare pressure tank data

############################################## define parameters of interest for EL system ################################################################
def add_parameters(model, unit_dct, country_specific = False):
    sys = model.system
    sys_stream = sys.flowsheet.stream
    param = model.parameter
    price_dct, GWP_dct, H_Ecosystems, H_Health, H_Resources  = update_resource_recovery_settings()
    
    # if specific country is not selected, add these parameters below
    if not country_specific:
        # price ratio
        old_price_dct = price_dct.copy()
        excretion_unit = unit_dct['Excretion']; # need change (??)
        b = 1
        D = shape.Uniform(lower = b - (10**(-6)), upper = b + (10**(-6)))
        item_ref = {
            'Concrete': 'Concrete',
            'Steel': 'Steel',
            }
        stream_ref = {
            'Struvite': 'struvite',
            'Salt': 'salt',
            'ClearWater': 'ClearWater',
            'NaOH': 'NaOH',
            'NaClO': 'NaClO',
            }
        @param(name = 'Price ratio', element = excretion_unit,# how to use element (???)
               kind = 'cost', units = '-',
               baseline = b, distribution = D)
        def set_price_ratio():
            EL.price_ratio = i
            for obj_name in (*item_ref.keys(), *stream_ref.keys()):
                old_price = old_price_dct[obj_name]
                new_price = old_price * i
                if obj_name in item_ref.keys():
                    ImpactItem.get_item(item_ref[obj_name]).price = new_price
                else:
                    getattr(sys_stream, stream_ref[obj_name]).price = new_price
            for u in sys.units:
                if hasattr(u, 'price_ratio'):
                    u.price_ratio = i
        
        # Household size
        b = EL.household_size
        D = shape.Normal(mu = b, sigma = 1.4)
        @param(name = 'Household size', element = excretion_unit, kind = 'coupled', units = 'cap / household',
            baseline = b, distribution = D)
        def set_household_size():
            EL.household_size = max(1, i)
        
        # Operator labor wage
        b = EL.operator_daily_wage; # operator_daily_wage is constant, needing initialization in _init_.py
        D = shape.Triangle(lower = (14.55), midpoint = b, upper = (43.68))
        D = shape.Triangle(lower = (14.55), midpoint = b, upper = (43.68))
        @param(name = 'Operator daily wage', element = 'TEA', kind = 'cost', units = 'USD / day',
            baseline = b, distribution = D)
        def set_operator_daily_wage(i):
            sys._TEA.annual_labor = i * 3 * 365 # need change

        # Construction labor wage
        b = EL.const_daily_wage   # const_daily_wage is constant, needing initialization in _init_.py
        D = shape.Triangle(lower = (b * 0.5), midpoint = b, upper = (b * 1.5))
        @param(name = 'Construction daily wage', element = 'TEA', kind = 'cost', units = 'USD / day',
            baseline = b, distribution = D)
        def set_const_daily_wage(i):
            for u in sys.units:
                if isinstance(u, qs.sanunits.EL_Housing): # Here EL_housing needs change
                    break
                u.const_daily_wage = i
        
        if EL.INCLUDE_RESOURCE_RECOVERY:
            # N fertilizer price
            b = 1.51; # need check up to now
            D = shape.Uniform(lower = b * 0.8, upper = b * 1.2)
            @param(name = 'N fertilizer price', element = 'TEA', kind = 'isolated', units = 'USD/kg N',
                   baseline = b, distribution = D)
            def set_N_fert_price(i):
                price_dct['N'] = sys_stream.liq_N.price = sys_stream.sol_N.price = i * EL.price_factor
            
            # P fertilizer price
            b = 3.9; # need check up to now
            D = shape.Uniform(lower = b * 0.8, upper = b * 1.2)
            @param(name = 'P fertilizer price', element = 'TEA', kind = 'isolated', units = 'USD/kg P',
                   baseline = b, distribution = D)
            def set_P_fert_price(i):
                price_dct['P'] = sys_stream.liq_P.price = sys_stream.sol_P.price = i * EL.price_factor
            
            # K fertilizer price
            b = 1.5; # need check up to now
            D = shape.Uniform(lower = b * 0.8, upper = b * 1.2)
            @param(name = 'K fertilizer price', element = 'TEA', kind = 'isolated', units = 'USD/kg K',
                   baseline = b, distribution = D)
            def set_K_fert_price(i):
                price_dct['K'] = sys_stream.liq_K.price = sys_stream.sol_K.price = i * EL.price_factor
            
            # Ammonium fertilizer price
            b = 0.8; # need check up to now
            D = shape.Uniform(lower = b * 0.8, upper = b * 1.2)
            @param(name = 'Ammonium fertilizer price', element = 'TEA', kind = 'isolated', units = 'USD/kg NH4',
                   baseline = b, distribution = D)
            def set_NH4_fert_price(i):
                price_dct['Ammonium'] = sys_stream.liq_NH4.price = sys_stream.sol_NH4.price = i * EL.price_factor
            
            # Struvite fertilizer price
            b = 3.983 * (31 / 245); # need check up to now
            D = shape.Uniform(lower = b * 0.8, upper = b * 1.2)
            @param(name = 'Struvite fertilizer price', element = 'TEA', kind = 'isolated', units = 'USD/kg P',
                   baseline = b, distribution = D)
            def set_struvite_fert_price(i):
                price_dct['Struvite'] = sys_stream.liq_struvite.price = sys_stream.sol_struvite.price = i * EL.price_factor
        
        # Electricity price
        b = price_dct['Electricity']
        D = shape.Triangle(lower = 0.08, midpoint = b, upper = 0.4);
        @param(name = 'Electricity price', element = 'TEA', kind = 'isolated', units = 'USD/kWh',
               baseline = b, distribution = D)
        def set_electricity_price(i):
            PowerUtility.price = i
        # Electricity GWP
        b = GWP_dct['Electricity']
        D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
        @param(name = 'Electricity CF', element = 'LCA', kind = 'isolated', units = 'kg CO2-eq/kWh',
               baseline = b, distribution = D)
        def set_electricity_CF(i):
            GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i
        # Electricity H_Ecosystems, if applicable
        b = H_Ecosystems_dct['Electricity']
        D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
        @param(name = 'Electricity Ecosystems CF', element = 'LCA', kind = 'isolated', units = 'points/kWh',
               baseline = b, distribution = D)
        def set_electricity_ecosystems_CF(i):
            H_Ecosystems_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Ecosystems'] = i
        # Electricity H_Health, if applicable
        b = H_Health_dct['Electricity']
        D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
        @param(name = 'Electricity Health CF', element = 'LCA', kind = 'isolated', units = 'points/kWh',
               baseline = b, distribution = D)
        def set_electricity_health_CF(i):
            H_Health_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Health'] = i
        # Electricity H_Resources, if applicable
        b = H_Resources_dct['Electricity']
        D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
        @param(name = 'Electricity Resources CF', element = 'LCA', kind = 'isolated', units = 'points/kWh',
               baseline = b, distribution = D)
        def set_electricity_resources_CF(i):
            H_Resources_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Resources'] = i
    
    ############################# Specific Units having parameters engaged with uncertainty and sensitivity analysis ###############################
    # In diet and excretion section
    excretion_unit = unit_dct['Excretion']
    exclude = ('e_cal', 'p_anim', 'p_veg') if country_specific else (); # e_cal: caloric_intake, p_anim: protein_animal_intake, 
                                                                        # p_veg: protein_vegetal_intake, all defined in _Excretion.tsv.
    batch_setting_unit_params(excretion_data, model, excretion_unit, exclude)
    #here [batch_setting_unit_params] means adding unit uncertain parameters through an Excel spreadsheet.
    
    # In toilet section
    toilet_unit = unit_dct['Toilet']
    b = EL.household_per_toilet
    D = shape.Uniform(lower = 6, upper = 12)
    @param(name = 'Toilet density', element = toilet_unit, kind = 'coupled', units = 'household/toilet',
           baseline = b, distribution = D)
    def set_toilet_density(i):
        EL.household_per_toilet = i
    
    batch_setting_unit_params(toilet_data, model, toilet_unit)
    
    # In Collection Tank
    CT_unit = unit_dct['Collection_Tank']
    batch_setting_unit_params(CT_data, model, CT_unit)
    
    # In Primary Clarifier Tank
    PC_unit = unit_dct['PrimaryClarifierTank']
    batch_setting_unit_params(PC_data, model, PC_unit)
    
    # In Anoxic Tank
    AnoxT_unit = unit_dct['AnoxicTank']
    batch_setting_unit_params(AnoxicTank_data, model, AnoxT_unit)
    
    # In Aerobic Tank
    AeroT_unit = unit_dct['AerobicTank']
    batch_setting_unit_params(AerobicTank_data, model, AeroT_unit)

    # In Blower for Aerobic Tank
    Blower_AeroT_unit = unit_dct['AerobicTankBlower']
    batch_setting_unit_params(Blower_AerobicTank, model, Blower_AeroT_unit)
    
    # In Membrane Tank
    MembT_unit = unit_dct['MembraneTank']
    batch_setting_unit_params(MembTank_data, model, MembT_unit)
    
    # In Membrane Tank Blower
    Blower_MembT_unit = unit_dct['MembraneTankBlower']
    batch_setting_unit_params(Blower_MembTank, model, Blower_MembT_unit)
    
    # In Clear Water Tank
    ClearWaterT_unit = unit_dct['ClearWaterTank']
    batch_setting_unit_params(ClearWaterTank_data, model, ClearWaterT_unit)
    
    # In Pressure Tank
    PressureT_unit = unit_dct['PressureTank']
    batch_setting_unit_params(PressureTank_data, model, PressureT_unit)
    

    ################################################ Universal degradation parameters ##########################################################
    # max methane emission
    toilet_unit = sys.path[1] # the first unit involving degradation
    b = EL.max_CH4_emission
    D = shape.Triangle(lower = 0.175, midpoint = b, upper = 0.325)
    @param(name = 'Max CH4 emission', element = toilet_unit, kind = 'coupled', units = 'g CH4/g COD',
           baseline = b, distribution = D)
    def set_max_CH4_emission(i):
        EL.max_CH4_emission = i
        for unit in sys.units:
            if hasattr(unit, 'max_CH4_emission'):
                setattr(unit, 'max_CH4_emission', i)
    
    # time to full degradation
    b =EL.tau_deg
    D = shape.Uniform(lower = 1, upper = 3)
    @param(name = 'Time to full degradation', element = toilet_unit, kind = 'coupled', units = 'yr',
           baseline = b, distribution = D)
    def set_tau_deg(i):
        EL.tau_deg = i
        k = get_decay_k(i, EL.log_deg)
        for unit in sys.units:
            if hasattr(unit, 'decay_k_COD'):
                setattr(unit, 'decay_k_COD', k)
            if hasattr(unit, 'decay_k_N'):
                setattr(unit, 'decay_k_N', k)

    # reduction at full degradation
    b = EL.log_deg
    D = shape.Uniform(lower = 2, upper = 4)
    @param(name = 'Log of degradation rate', element = toilet_unit, kind = 'coupled', units = '-',
           baseline = b, distribution = D)
    def set_log_deg(i):
        EL.log_deg = i
        k = get_decay_k(EL.tau_deg, i)
        for unit in sys.units:
            if hasattr(unit, 'decay_k_COD'):
                setattr(unit, 'decay_k_COD', k)
            if hasattr(unit, 'decay_k_N'):
                setattr(unit, 'decay_k_N', k)
    
    # toilet material properties
    density = toilet_unit.density_dct
    #########################################################################################################################
    b = density['Plastic']
    D = shape.Uniform(lower = 0.31, upper = 1.24)
    param(setter = DictAttrSetter(toilet_unit, 'density_dct', 'Plastic'),
          name = 'Density of plastic', element = toilet_unit, kind = 'isolated', units = 'kg/m3',
          baseline = b, distribution = D)
    ########################################################################################################################
    b = density['Brick']
    D = shape.Uniform(lower = 1500, upper = 2000);
    param(setter = DictAttrSetter(toilet_unit, 'density_dct', 'Brick'),
          name = 'Density of brick', element = toilet_unit, kind = 'isolated', units = 'kg/m3',
          baseline = b, distribution = D)
    ########################################################################################################################
    b = density['StainlessSteelSheet']
    D = shape.Uniform(lower = 2.26, upper = 3.58)
    param(setter = DictAttrSetter(toilet_unit, 'density_dct', 'StainlessSteelSheet'),
          name = 'Density of stainless steel sheet', element = toilet_unit, kind = 'isolated', units = 'kg/m3',
          baseline = b, distribution = D)
    ########################################################################################################################
    b = density['Gravel']
    D = shape.Uniform(lower = 1520, upper = 1680)
    param(setter = DictAttrSetter(toilet_unit, 'density_dct', 'Gravel'),
          name = 'Density of gravel', element = toilet_unit, kind = 'isolated', units = 'kg/m3',
          baseline = b, distribution = D)
    ########################################################################################################################      
    b = density['Sand']
    D = shape.Uniform(lower = 1281, upper = 1602)
    param(setter = DictAttrSetter(toilet_unit, 'density_dct', 'Sand'),
          name = 'Density of sand', element = toilet_unit, kind = 'isolated', units = 'kg/m3',
          baseline = b, distribution = D)
    #########################################################################################################################
    b = density['Steel']
    D = shape.Uniform(lower = 7750, upper = 8050)
    param(setter = DictAttrSetter(toilet_unit, 'density_dct', 'Steel'),
          name = 'Density of steel', element = toilet_unit, kind = 'isolated', units = 'kg/m3',
          baseline = b, distribution = D)
    
    ################################################# General TEA Settings ##################################################
    # Discount rate if changing in TEA
    b = EL.discount_rate
    D = shape.Uniform(lower = 0.02, upper = 0.06)
    @param(name = 'Discount rate', element = 'TEA', kind = 'isolated', units = 'fraction',
          baseline = b, distribution = D)
    def set_discount_rate(i):
        EL.discount_rate = i
    
    # if resource recovery settings are updated, change price_factor for uncertainty analysis
    b = EL.price_factor
    D = shape.Uniform(lower = 0.1, upper = 0.4)
    @param(name = 'Price factor', element = 'TEA', kind = 'isolated', units = '-',
          baseline = b, distribution = D)
    def set_price_factor(i):
        EL.price_factor = i
    
    # Some chemicals defined in price_dct of _init_.py need update, the below format can be used to add it for uncertainty analysis
    # NaOH
    b = price_dct['NaOH']
    D = shape.Uniform(lower = (b * 0.9), upper = (b * 1.1)); # here (b * 0.9), for example, 
                                                             #can enable b firstly to be calculated if b was present in terms of complicate,
                                                             # explicit formulas.
    @param(name = 'Price of NaOH', element = 'TEA', kind = 'isolated', units = 'USD/kg NaOH',
          baseline = (b), distribution = D)
    def set_NaOH_price(i):
        price_dct['NaOH'] = sys_stream.NaOH.price = i
    
    # NaClO   
    b = price_dct['NaClO']
    D = shape.Uniform(lower = (b * 0.9), upper = (b * 1.1));
    @param(name = 'Price of NaClO', element = 'TEA', kind = 'isolated', units = 'USD/kg NaClO',
          baseline = (b), distribution = D)
    def set_NaClO_price(i):
        price_dct['NaClO'] = sys_stream.NaClO.price = i

    
    ###################################################### General LCA Settings ###############################################################
    # CH4 in GWP, H_Ecosystems, H_Health
    b = GWP_dct['CH4']
    D = shape.Uniform(lower = b * 0.9, upper = b * 1.1)
    @param(name = 'CH4 CF', element = 'LCA', kind = 'isolated', units = 'kg CO2-eq/kg CH4',
          baseline = b, distribution = D)
    def set_CH4_CF(i):
        GWP_dct['CH4'] = ImpactItem.get_item('CH4_item').CFs['GlobalWarming'] = i
    # if H_Ecosystems, H_Health, and H_Resources indicators will be scheduled in LCA, the below format can be used to define CH4 in them
    # for uncertainty analysis.
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
    
    # NaOH in GWP, H_Ecosystems, H_Health, H_Resources
    b = GWP_dct['NaOH']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'NaOH CF', element = 'LCA', kind = 'isolated', units = 'kg CO2-eq/kg NaOH',
          baseline = b, distribution = D)
    def set_NaOH_CF(i):
        GWP_dct['NaOH'] = ImpactItem.get_item('NaOH_item').CFs['GlobalWarming'] = i
    
    b = H_Ecosystems_dct['NaOH']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'NaOH Ecosystems CF', element = 'LCA', kind = 'isolated', units = 'points/kg NaOH',
          baseline = b, distribution = D)
    def set_NaOH_Ecosystems_CF(i):
        H_Ecosystems_dct['NaOH'] = ImpactItem.get_item('NaOH_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['NaOH']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'NaOH Health CF', element = 'LCA', kind = 'isolated', units = 'points/kg NaOH',
          baseline = b, distribution = D)
    def set_NaOH_Health_CF(i):
        H_Health_dct['NaOH'] = ImpactItem.get_item('NaOH_item').CFs['H_Health'] = i
    
    b = H_Resources_dct['NaOH']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'NaOH Resources CF', element = 'LCA', kind = 'isolated', units = 'points/kg NaOH',
          baseline = b, distribution = D)
    def set_NaOH_Resources_CF(i):
        H_Resources_dct['NaOH'] = ImpactItem.get_item('NaOH_item').CFs['H_Resources'] = i
    
    # NaClO in GWP, H_Ecosystems, H_Health, H_Resources
    b = H_GWP_dct['NaClO']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'NaClO CF', element = 'LCA', kind = 'isolated', units = 'kg CO2-eq/kg NaClO',
          baseline = b, distribution = D)
    def set_NaClO_CF(i):
        H_GWP_dct['NaClO'] = ImpactItem.get_item('NaClO_item').CFs['GlobalWarming'] = i
    
    b = H_Ecosystems_dct['NaClO']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'NaClO Ecosystems CF', element = 'LCA', kind = 'isolated', units = 'points/kg NaClO',
          baseline = b, distribution = D)
    def set_NaClO_Ecosystems_CF(i):
        H_Ecosystems_dct['NaClO'] = ImpactItem.get_item('NaClO_item').CFs['H_Ecosystems'] = i
    
    b = H_Health_dct['NaClO']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'NaClO Health CF', element = 'LCA', kind = 'isolated', units = 'points/kg NaClO',
          baseline = b, distribution = D)
    def set_NaClO_Health_CF(i):
        H_Health_dct['NaClO'] = ImpactItem.get_item('NaClO_item').CFs['H_Health'] = i
    
    b = H_Resources_dct['NaClO']
    D = shape.Triangle(lower = b * 0.9, midpoint = b, upper = b * 1.1)
    @param(name = 'NaClO Resources CF', element = 'LCA', kind = 'isolated', units = 'points/kg NaClO',
          baseline = b, distribution = D)
    def set_NaClO_Resources_CF(i):
        H_Resources_dct['NaClO'] = ImpactItem.get_item('NaClO_item').CFs['H_Resources'] = i
    
    if EL.INCLUDED_RESOURCE_RECOVERY:
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
        b = 5.4 * (14 / 18); # refer to Reclaimer and Biogenic Refinery, as below
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Ammonium fertilizer CF', element = 'LCA', kind = 'isolated',
               units = 'kg CO2-eq/kg Ammonium', baseline = b, distribution = D)
        def set_Ammonium_fertilizer_CF(i):
            GWP_dct['Ammonium'] = ImpactItem.get_item('Ammonium_item').CFs['GlobalWarming'] = -i
        
        b = 0.0461961 * (14 / 18)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Ammonium fertilizer Ecosystems CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Ammonium', baseline = b, distribution = D)
        def set_Ammonium_fertilizer_Ecosystems_CF(i):
            H_Ecosystems_dct['Ammonium'] = ImpactItem.get_item('Ammonium_item').CFs['H_Ecosystems'] = -i

        b = 0.637826734 * (14 / 18)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Ammonium fertilizer Health CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Ammonium', baseline = b, distribution = D)
        def set_Ammonium_fertilizer_Health_CF(i):
            H_Health_dct['Ammonium'] = ImpactItem.get_item('Ammonium_item').CFs['H_Health'] = -i
        
        b = 0.259196888 * (14 / 18)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Ammonium fertilizer Resources CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Ammonium', baseline = b, distribution = D)
        def set_Ammonium_fertilizer_Resources_CF(i):
            H_Resources_dct['Ammonium'] = ImpactItem.get_item('Ammonium_item').CFs['H_Resources'] = -i
        
        # Recovered Struvite as alternative fertilizer, in GWP, H_Ecosystems, H_Health, H_Resources
        b = 4.9 * (30.97 / 245.41)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Struvite fertilizer CF', element = 'LCA', kind = 'isolated',
               units = 'kg CO2-eq/kg Struvite', baseline = b, distribution = D)
        def set_Struvite_fertilizer_CF(i):
            GWP_dct['Struvite'] = ImpactItem.get_item('struvite_item').CFs['GlobalWarming'] = -i
        
        b = 0.093269908 * (30.97 / 245.41)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Struvite fertilizer Ecosystems CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Struvite', baseline = b, distribution = D)
        def set_Struvite_fertilizer_Ecosystems_CF(i):
            H_Ecosystems_dct['Struvite'] = ImpactItem.get_item('struvite_item').CFs['H_Ecosystems'] = -i
        
        b = 1.774294425 * (30.97 / 245.41)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Struvite fertilizer Health CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Struvite', baseline = b, distribution = D)
        def set_Struvite_fertilizer_Health_CF(i):
            H_Health_dct['Struvite'] = ImpactItem.get_item('struvite_item').CFs['H_Health'] = -i
        
        b = 1.084191599 * (30.97 / 245.41)
        D = shape.Triangle(lower = b * 0.90, midpoint = b, upper = b * 1.1)
        @param(name = 'Struvite fertilizer Resources CF', element = 'LCA', kind = 'isolated',
               units = 'points/kg Struvite', baseline = b, distribution = D)
        def set_Struvite_fertilizer_Resources_CF(i):
            H_Resources_dct['Struvite'] = ImpactItem.get_item('struvite_item').CFs['H_Resources'] = -i
        

########################################################### Add other CFs ###############################################################
item_path = os.path.join(EL_data_path, 'impact_items.xlsx')

for indicator in ('GlobalWarming',
                  'H_Ecosystems',
                  'H_Health',
                  'H_Resources'):
    sheet_name = indicator if indicator != 'GlobalWarming' else 'GWP'   # GWP is the default sheet name
    data = load_data(item_path, sheet = sheet_name)   # load data from Excel file
    for p in data.index:
        item = ImpactItem.get_item(p)
        b = item.CFs[indicator]
        lower = float(data.loc[p]['low'])   # need aligning with that in prepared files
        upper = float(data.loc[p]['high'])    # need aligning with that in prepared files
        dist = data.loc[p]['distribution']    # need aligning with that in prepared files
        if dist == 'uniform':
            D = shape.Uniform(lower = lower, upper = upper)
        elif dist == 'triangular':
            D = shape.Triangle(lower = lower, midpoint = b, upper = upper)
        elif dist == 'constant':
            continue
        else:
            raise ValueError(f'Distribution {dist} not recognized.')
        model.parameter(name = p + f'-{indicator}',
                         setter = DictAttrSetter(item, 'CFs', indicator),
                         element = 'LCA',
                         kind = 'isolated',
                         units = f'kg CO2-eq/{item.functional_unit}', # here make sense function_unit
                         baseline = b, distribution = D)
                        
                        

#Import data of separate units prepared in the _EnvirolooUnits.py file for uncertainty and sensitivity analysis'''
# for example, _EL_CT.tsv is called
#def add_CollectionTank_parameters(model, unit_dct):
    #CT_unit = unit_dct['CollectionTank'];
    #param = model.parameter;
    # need continuing

#Create Model for EL system
def create_model(country_specific = False, **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)    # figure out the content of **model_kwargs
    sysEL = create_system('sysEL', flowsheet=flowsheet)
    unitEL = sysEL.flowsheet.unit
    
    # add parameters for the EL system
    modelEL = Model(sysEL, **model_kwargs)
    add_metrics(modelEL)
    unit_dct = { # name here needs to be aligned with those in Enviroloo_system.py (!!!)
        'Excretion': unitEL.WasteWaterGenerator,
        'Toilet': unitEL.Toilet,
        'Collection_Tank': unitEL.CT,
        'Lifting_Pump': unitEL.P_CT_lift,
        'PrimaryClarifierTank': unitEL.PC,
        'PrimaryClarifierReturnPump': unitEL.P_PC_return,
        'AnoxicTank': unitEL.AnoxT,
        'GlucoseAgitationPump': unitEL.P_Glu_agitation,
        'GlucoseDosingPump': unitEL.P_Glu_dosing,
        'AnoxicTankAgitationPump': unitEL.P_AnoxT_agitation,
        'AerobicTank': unitEL.AerT,
        'PACAgitationPump': unitEL.P_PAC_agitation,
        'PACDosingPump': unitEL.P_PAC_dosing,
        'AerobicTankBlower': unitEL.Aerobic_blower,
        'MembraneTank': unitEL.MembT,
        'NitrateReturnPumpToPC': unitEL.P_NitrateReturn_PC,
        'NitrateReturnPumpToAnoxicTank': unitEL.P_NitrateReturn_AnoxT,
        'MembraneTankBlower': unitEL.MembT_blower,
        'SelfPrimingPump': unitEL.P_MT_selfpriming,
        'ClearWaterTank': unitEL.CWT,
        'O3Generator': unitEL.O3_gen,
        'O3DosingPump': unitEL.P_O3_dosing,
        'AirDissolvedPump': unitEL.P_AirDissolvedP,
        'PressurePumptoClearWater': unitEL.P_CWT,
        'PressureTank': unitEL.PT,
        }
    
    add_parameters(modelEL, unit_dct, country_specific);
    
    return modelEL

# define runing function intializing uncertainty and sensitivity analysis
def run_uncertainty(model, path='', **kwargs):
    kwargs['path'] = os.path.join(results_path, f'sysEL_model.xlsx') if path=='' else path
    run(model = model, **kwargs)
    return