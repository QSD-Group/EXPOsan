#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
DMsan: Decision-making for sanitation and resource recovery systems

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>
'''


# %%

import os, qsdsan as qs
from chaospy import distributions as shape
from qsdsan import Model, Metric, PowerUtility, ImpactItem
from qsdsan.utils import (
    data_path,
    DictAttrSetter,
    load_data,
    )
from exposan.utils import batch_setting_unit_params, run_uncertainty as run
from exposan import reclaimer as re
from exposan.reclaimer import (
    _load_components,
    create_system,
    data_path as re_data_path,
    get_decay_k,
    get_LCA_metrics,
    get_TEA_metrics,
    get_recoveries,
    GWP_dct,
    H_Ecosystems_dct,
    H_Health_dct,
    H_Resources_dct,
    price_dct,
    results_path,
    )

__all__ = ('create_model', 'run_uncertainty',)


# %%

# =============================================================================
# Functions for batch-making metrics and setting parameters
# =============================================================================

def add_metrics(model):
    re._load_lca_data()
    system = model.system
    # Recoveries
    funcs = get_recoveries(system)
    metrics = [
        Metric('Total N', funcs[0], '% N', 'N recovery'),
        Metric('Total P', funcs[1], '% P', 'P recovery'),
        Metric('Total K', funcs[2], '% K', 'K recovery'),
    ]
    # Net cost
    metrics.append(
        Metric('Annual net cost', get_TEA_metrics(system)[0], f'{qs.currency}/cap/yr', 'TEA results'),
        )
    # Net emissions
    funcs = get_LCA_metrics(system)
    cat = 'LCA results'
    metrics.extend([
        Metric('GlobalWarming', funcs[0], 'kg CO2-eq/cap/yr', cat),
        Metric('H_Ecosystems', funcs[1], 'points/cap/yr', cat),
        Metric('H_Health', funcs[2], 'points/cap/yr', cat),
        Metric('H_Resources', funcs[3], 'points/cap/yr', cat),
        ])
    model.metrics = metrics


# %%

# =============================================================================
# Data sheets
# =============================================================================

su_data_path = os.path.join(data_path, 'sanunit_data')
re_su_data_path = os.path.join(su_data_path, 're')

def load_re_su_data(file_name):
    if file_name.startswith('_re'):
        return load_data(os.path.join(re_su_data_path, file_name))
    return load_data(os.path.join(su_data_path, file_name))

excretion_data = load_re_su_data('_excretion.tsv')
murt_data = load_re_su_data('_murt.tsv')
pasteurization_data = load_re_su_data('_sludge_pasteurization.tsv')
septic_tank_data = load_re_su_data('_septic_tank.csv')

ultrafiltration_data = load_re_su_data('_re_ultrafiltration.csv')
ix_data = load_re_su_data('_re_ion_exchange.csv')
ecr_data = load_re_su_data('_re_ecr.csv')
housing_data = load_re_su_data('_re_housing.csv')
system_data = load_re_su_data('_re_system.csv')
solar_data = load_re_su_data('_re_solar.csv')


# %%

# =============================================================================
# Shared by all systems
# =============================================================================

def add_shared_parameters(model, unit_dct, country_specific=False):
    sys = model.system
    sys_stream = sys.flowsheet.stream
    param = model.parameter

    # Add these parameters if not running country-specific analysis,
    # in which they would be updated separately

    if not country_specific:

        # Labor wage
        b = price_dct['wages']
        D = shape.Triangle(lower=1.82, midpoint=b, upper=5.46)  # wage in USD/hour
        @param(name='Labor wages', element='TEA', kind='cost', units='USD/h',
               baseline=b, distribution=D)
        def set_labor_wages(i):
            for u in sys.units:
                if hasattr(u, '_calc_maintenance_labor_cost'):
                    u.wages = i

        # Electricity price
        b = price_dct['Electricity']
        D = shape.Uniform(lower=0.045, upper=0.075)
        @param(name='Electricity price', element='TEA', kind='isolated',
               units='$/kWh', baseline=b, distribution=D)
        def set_electricity_price(i):
            PowerUtility.price = i

        # Electricity GWP
        b = GWP_dct['Electricity']
        D = shape.Triangle(lower=b*0.9, midpoint=b, upper=b*1.1)
        @param(name='Electricity CF', element='LCA', kind='isolated',
                   units='kg CO2-eq/kWh', baseline=b, distribution=D)
        def set_electricity_CF(i):
            GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i

        # Electricity H_Ecosystems
        b = H_Ecosystems_dct['Electricity']
        D = shape.Triangle(lower=b*0.9, midpoint=b, upper=b*1.1)
        @param(name='Electricity Ecosystems CF', element='LCA', kind='isolated',
                   units='points/kWh', baseline=b, distribution=D)
        def set_electricity_ecosystems_CF(i):
            H_Ecosystems_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Ecosystems'] = i

        # Electricity H_Health
        b = H_Health_dct['Electricity']
        D = shape.Triangle(lower=b*0.9, midpoint=b, upper=b*1.1)
        @param(name='Electricity Health CF', element='LCA', kind='isolated',
                   units='points/kWh', baseline=b, distribution=D)
        def set_electricity_health_CF(i):
            H_Health_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Health'] = i

        # Electricity H_Resources
        b = H_Resources_dct['Electricity']
        D = shape.Triangle(lower=b*0.9, midpoint=b, upper=b*1.1)
        @param(name='Electricity Resources CF', element='LCA', kind='isolated',
                   units='points/kWh', baseline=b, distribution=D)
        def set_electricity_resources_CF(i):
            H_Resources_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['H_Resources'] = i

        # N fertilizer price
        b = 1.507
        D = shape.Uniform(lower=b * 0.8, upper=b * 1.2)

        @param(name='N fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
               baseline=b, distribution=D)
        def set_N_price(i):
            price_dct['N'] = sys_stream.liq_N.price = sys_stream.sol_N.price = i * re.price_factor

        # P fertilizer price
        b = 3.983
        D = shape.Uniform(lower=b * 0.8, upper=b * 1.2)

        @param(name='P fertilizer price', element='TEA', kind='isolated', units='USD/kg P',
               baseline=b, distribution=D)
        def set_P_price(i):
            price_dct['P'] = sys_stream.liq_P.price = sys_stream.sol_P.price = i * re.price_factor

        # K fertilizer price
        b = 1.333
        D = shape.Uniform(lower=b * 0.8, upper=b * 1.2)

        @param(name='K fertilizer price', element='TEA', kind='isolated', units='USD/kg K',
               baseline=b, distribution=D)
        def set_K_price(i):
            price_dct['K'] = sys_stream.liq_K.price = sys_stream.sol_K.price = i * re.price_factor

        # NH3 fertilizer price
        D = shape.Uniform(lower=(1.507 * (14 / 17) * 0.8), upper=(1.507 * (14 / 17) * 1.2))

        @param(name='NH3 fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
               baseline=(1.507 * (14 / 17)), distribution=D)
        def set_con_NH3_price(i):
            price_dct['conc_NH3'] = sys_stream.conc_NH3.price = i * re.price_factor

    ##### Specific units #####
    # Diet and excretion
    excretion_unit = unit_dct['Excretion']
    exclude = ('e_cal', 'p_anim', 'p_veg') if country_specific else ()
    batch_setting_unit_params(excretion_data, model, excretion_unit, exclude)

    # MURT
    murt_unit = unit_dct.get('Toilet')
    if murt_unit:
        exclude = ('MCF_decay', 'N2O_EF_decay', 'OPEX_over_CAPEX')
        batch_setting_unit_params(murt_data, model, murt_unit, exclude)

        b = murt_unit.OPEX_over_CAPEX
        D = shape.Uniform(lower=0.02, upper=0.08)
        @param(name='MURT operating cost', element=murt_unit, kind='coupled', units='cost',
               baseline=b, distribution=D)
        def set_OPEX_over_CAPEX(i):
            murt_unit.OPEX_over_CAPEX = i

        b = murt_unit.MCF_decay
        D = shape.Triangle(lower=0.05, midpoint=b, upper=0.15)
        @param(name='MCF_decay', element=murt_unit, kind='coupled',
               units='fraction of anaerobic conversion of degraded COD',
               baseline=b, distribution=D)
        def set_MCF_decay(i):
            murt_unit.MCF_decay = i

        b = murt_unit.N2O_EF_decay
        D = shape.Triangle(lower=0, midpoint=b, upper=0.001)
        @param(name='N2O_EF_decay', element=murt_unit, kind='coupled',
               units='fraction of N emitted as N2O',
               baseline=b, distribution=D)
        def set_N2O_EF_decay(i):
            murt_unit.N2O_EF_decay = i

    # Septic tank (primary treatment)
    primary_unit = unit_dct['Primary']
    batch_setting_unit_params(septic_tank_data, model, primary_unit)

    b = primary_unit.sludge_moisture_content
    D = shape.Uniform(lower=0.9, upper=0.98)
    param(name='Sludge moisture content', element=primary_unit, kind='coupled', units='fraction',
          baseline=b, distribution=D)
    def set_primary_sludge_moisture_content(i):
        primary_unit.sludge_moisture_content = i

    # Sludge pasteurization
    # `wages` should be excluded for all units this this attribute
    exclude = 'wages' if country_specific else ()
    pasteurization_unit = unit_dct.get('Pasteurization')
    if pasteurization_unit: batch_setting_unit_params(pasteurization_data, model, pasteurization_unit, exclude)

    # Ultrafiltration
    ultrafiltration_unit = unit_dct.get('Ultrafiltration')
    if ultrafiltration_unit: batch_setting_unit_params(ultrafiltration_data, model, ultrafiltration_unit)

    # Ion exchange
    ix_unit = unit_dct.get('Ion exchange')
    if ix_unit: batch_setting_unit_params(ix_data, model, ix_unit, exclude)

    # Electrochemical reactor
    ecr_unit = unit_dct.get('ECR')
    if ecr_unit: batch_setting_unit_params(ecr_data, model, ecr_unit, exclude)

    # Housing
    housing_unit = unit_dct.get('Housing')
    if housing_unit: batch_setting_unit_params(housing_data, model, housing_unit)

    # System
    system_unit = unit_dct.get('System')
    if system_unit: batch_setting_unit_params(system_data, model, system_unit)

    # Solar
    solar_unit = unit_dct.get('Solar')
    if solar_unit: batch_setting_unit_params(solar_data, model, solar_unit, exclude)

    ##### Universal degradation parameters #####
    # Max methane emission
    unit = sys.path[1]  # the first unit that involves degradation
    b = re.max_CH4_emission
    D = shape.Triangle(lower=0.175, midpoint=b, upper=0.325)
    @param(name='Max CH4 emission', element=unit, kind='coupled', units='g CH4/g COD',
           baseline=b, distribution=D)
    def set_max_CH4_emission(i):
        re.max_CH4_emission = i
        for unit in sys.units:
            if hasattr(unit, 'max_CH4_emission'):
                setattr(unit, 'max_CH4_emission', i)

    # Time to full degradation
    b = re.tau_deg
    D = shape.Uniform(lower=1, upper=3)
    @param(name='Full degradation time', element=unit, kind='coupled', units='yr',
           baseline=b, distribution=D)
    def set_tau_deg(i):
        re.tau_deg = i
        k = get_decay_k(i, re.log_deg)
        for unit in sys.units:
            if hasattr(unit, 'decay_k_COD'):
                setattr(unit, 'decay_k_COD', k)
            if hasattr(unit, 'decay_k_N'):
                setattr(unit, 'decay_k_N', k)

    # Reduction at full degradation
    b = re.log_deg
    D = shape.Uniform(lower=2, upper=4)
    @param(name='Log degradation', element=unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_log_deg(i):
        re.log_deg = i
        k = get_decay_k(re.tau_deg, i)
        for unit in sys.units:
            if hasattr(unit, 'decay_k_COD'):
                setattr(unit, 'decay_k_COD', k)
            if hasattr(unit, 'decay_k_N'):
                setattr(unit, 'decay_k_N', k)

    ##### General TEA settings #####

    # # Keeping discount rate constant
    # b = systems.discount_rate
    # D = shape.Uniform(lower=0.03, upper=0.06)
    # @param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
    #         baseline=b, distribution=D)
    # def set_discount_rate(i):
    #     systems.discount_rate = i

    # Discount factor for the excreta-derived fertilizers
    b = re.price_factor
    D = shape.Uniform(lower=0.1, upper=0.4)
    @param(name='Price factor', element='TEA', kind='isolated', units='-',
           baseline=b, distribution=D)
    def set_price_factor(i):
        re.price_factor = i

    ##### General LCA settings #####

    # CH4
    b = GWP_dct['CH4']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='CH4 CF', element='LCA', kind='isolated', units='kg CO2-eq/kg CH4',
           baseline=b, distribution=D)
    def set_CH4_CF(i):
        GWP_dct['CH4'] = ImpactItem.get_item('CH4_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['CH4']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='CH4 ecosystems CF', element='LCA', kind='isolated', units='points/kg CH4',
           baseline=b, distribution=D)
    def set_CH4_ecosystems_CF(i):
        H_Ecosystems_dct['CH4'] = ImpactItem.get_item('CH4_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['CH4']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='CH4 health CF', element='LCA', kind='isolated', units='points/kg CH4',
           baseline=b, distribution=D)
    def set_CH4_health_CF(i):
        H_Health_dct['CH4'] = ImpactItem.get_item('CH4_item').CFs['H_Health'] = i

    # N2O
    b = GWP_dct['N2O']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='N2O CF', element='LCA', kind='isolated', units='kg CO2-eq/kg N2O',
           baseline=b, distribution=D)
    def set_N2O_CF(i):
        GWP_dct['N2O'] = ImpactItem.get_item('N2O_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['N2O']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='N2O ecosystems CF', element='LCA', kind='isolated', units='points/kg N2O',
           baseline=b, distribution=D)
    def set_N2O_ecosystems_CF(i):
        H_Ecosystems_dct['N2O'] = ImpactItem.get_item('N2O_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['N2O']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='N2O health CF', element='LCA', kind='isolated', units='points/kg N2O',
           baseline=b, distribution=D)
    def set_N2O_health_CF(i):
        H_Health_dct['N2O'] = ImpactItem.get_item('N2O_item').CFs['H_Health'] = i

    # Recovered N fertilizer
    b = -GWP_dct['N']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='N fertilizer CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg N', baseline=b, distribution=D)
    def set_N_fertilizer_CF(i):
        GWP_dct['N'] = ImpactItem.get_item('N_item').CFs['GlobalWarming'] = -i

    b = -H_Ecosystems_dct['N']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='N fertilizer ecosystems CF', element='LCA', kind='isolated',
           units='points/kg N', baseline=b, distribution=D)
    def set_N_fertilizer_ecosystems_CF(i):
        H_Ecosystems_dct['N'] = ImpactItem.get_item('N_item').CFs['H_Ecosystems'] = -i

    b = -H_Health_dct['N']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='N fertilizer health CF', element='LCA', kind='isolated',
           units='points/kg N', baseline=b, distribution=D)
    def set_N_fertilizer_health_CF(i):
        H_Health_dct['N'] = ImpactItem.get_item('N_item').CFs['H_Health'] = -i

    b = -H_Resources_dct['N']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='N fertilizer resources CF', element='LCA', kind='isolated',
           units='points/kg N', baseline=b, distribution=D)
    def set_N_fertilizer_resources_CF(i):
        H_Resources_dct['N'] = ImpactItem.get_item('N_item').CFs['H_Resources'] = -i

    # Recovered P fertilizer
    b = -GWP_dct['P']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='P fertilizer CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg P', baseline=b, distribution=D)
    def set_P_fertilizer_CF(i):
        GWP_dct['P'] = ImpactItem.get_item('P_item').CFs['GlobalWarming'] = -i

    b = -H_Ecosystems_dct['P']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='P fertilizer ecosystems CF', element='LCA', kind='isolated',
           units='points/kg P', baseline=b, distribution=D)
    def set_P_fertilizer_ecosystems_CF(i):
        H_Ecosystems_dct['P'] = ImpactItem.get_item('P_item').CFs['H_Ecosystems'] = -i

    b = -H_Health_dct['P']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='P fertilizer health CF', element='LCA', kind='isolated',
           units='points/kg P', baseline=b, distribution=D)
    def set_P_fertilizer_health_CF(i):
        H_Health_dct['P'] = ImpactItem.get_item('P_item').CFs['H_Health'] = -i

    b = -H_Resources_dct['P']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='P fertilizer resources CF', element='LCA', kind='isolated',
           units='points/kg P', baseline=b, distribution=D)
    def set_P_fertilizer_resources_CF(i):
        H_Resources_dct['P'] = ImpactItem.get_item('P_item').CFs['H_Resources'] = -i

    # Recovered K fertilizer
    b = -GWP_dct['K']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='K fertilizer CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg K', baseline=b, distribution=D)
    def set_K_fertilizer_CF(i):
        GWP_dct['K'] = ImpactItem.get_item('K_item').CFs['GlobalWarming'] = -i

    b = -H_Ecosystems_dct['K']
    D = shape.Triangle(lower=b * 0.90, midpoint=b, upper=b * 1.1)
    @param(name='K fertilizer ecosystems CF', element='LCA', kind='isolated',
           units='points/kg K', baseline=b, distribution=D)
    def set_K_fertilizer_ecosystems_CF(i):
        H_Ecosystems_dct['K'] = ImpactItem.get_item('K_item').CFs['H_Ecosystems'] = -i

    b = -H_Health_dct['K']
    D = shape.Triangle(lower=b * 0.90, midpoint=b, upper=b * 1.1)
    @param(name='K fertilizer health CF', element='LCA', kind='isolated',
           units='points/kg K', baseline=b, distribution=D)
    def set_K_fertilizer_health_CF(i):
        H_Health_dct['K'] = ImpactItem.get_item('K_item').CFs['H_Health'] = -i

    b = -H_Resources_dct['K']
    D = shape.Triangle(lower=b * 0.90, midpoint=b, upper=b * 1.1)
    @param(name='K fertilizer resources CF', element='LCA', kind='isolated',
           units='points/kg K', baseline=b, distribution=D)
    def set_K_fertilizer_resources_CF(i):
        H_Resources_dct['K'] = ImpactItem.get_item('K_item').CFs['H_Resources'] = -i

    # Recovered concentrated NH3 fertilizer
    b = -GWP_dct['conc_NH3']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Conc NH3 CF', element='LCA', kind='isolated',
            units='kg CO2-eq/kg', baseline=b, distribution=D)
    def set_conc_NH3_CF(i):
        GWP_dct['conc_NH3'] = ImpactItem.get_item('conc_NH3_item').CFs['GlobalWarming'] = -i

    b = -H_Ecosystems_dct['conc_NH3']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Conc NH3 ecosystems CF', element='LCA', kind='isolated',
            units='points/kg', baseline=b, distribution=D)
    def set_conc_NH3_ecosystems_CF(i):
        H_Ecosystems_dct['conc_NH3'] = ImpactItem.get_item('conc_NH3_item').CFs['H_Ecosystems'] = -i

    b = -H_Health_dct['conc_NH3']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Conc NH3 health CF', element='LCA', kind='isolated',
            units='points/kg', baseline=b, distribution=D)
    def set_conc_NH3_health_CF(i):
        H_Health_dct['conc_NH3'] = ImpactItem.get_item('conc_NH3_item').CFs['H_Health'] = -i

    b = -H_Resources_dct['conc_NH3']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Conc NH3 resources CF', element='LCA', kind='isolated',
            units='points/kg', baseline=b, distribution=D)
    def set_conc_NH3_resources_CF(i):
        H_Resources_dct['conc_NH3'] = ImpactItem.get_item('conc_NH3_item').CFs['H_Resources'] = -i

    # GAC
    b = GWP_dct['GAC']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='GAC CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg', baseline=b, distribution=D)
    def set_GAC_CF(i):
        GWP_dct['GAC'] = ImpactItem.get_item('GAC_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['GAC']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='GAC ecosystems CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_GAC_ecosystems_CF(i):
        H_Ecosystems_dct['GAC'] = ImpactItem.get_item('GAC_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['GAC']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='GAC health CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_GAC_health_CF(i):
        H_Health_dct['GAC'] = ImpactItem.get_item('GAC_item').CFs['H_Health'] = i

    b = H_Resources_dct['GAC']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='GAC resources CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_GAC_resources_CF(i):
        H_Resources_dct['GAC'] = ImpactItem.get_item('GAC_item').CFs['H_Resources'] = i

    # Zeolite
    b = GWP_dct['zeolite']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='zeolite CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg', baseline=b, distribution=D)
    def set_zeolite_CF(i):
        GWP_dct['zeolite'] = ImpactItem.get_item('zeolite_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['zeolite']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='Zeolite ecosystems CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_zeolite_ecosystems_CF(i):
        H_Ecosystems_dct['zeolite'] = ImpactItem.get_item('zeolite_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['zeolite']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='Zeolite health CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_zeolite_health_CF(i):
        H_Health_dct['zeolite'] = ImpactItem.get_item('zeolite_item').CFs['H_Health'] = i

    b = H_Resources_dct['zeolite']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='Zeolite resources CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_zeolite_resources_CF(i):
        H_Resources_dct['zeolite'] = ImpactItem.get_item('zeolite_item').CFs['H_Resources'] = i

    # KCl
    b = GWP_dct['KCl']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='KCl CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg', baseline=b, distribution=D)
    def set_KCl_CF(i):
        GWP_dct['KCl'] = ImpactItem.get_item('KCl_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['KCl']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='KCl ecosystems CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_KCl_ecosystems_CF(i):
        H_Ecosystems_dct['KCl'] = ImpactItem.get_item('KCl_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['KCl']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='KCl health CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_KCl_health_CF(i):
        H_Health_dct['KCl'] = ImpactItem.get_item('KCl_item').CFs['H_Health'] = i

    b = H_Resources_dct['KCl']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='KCl resources CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_KCl_resources_CF(i):
        H_Resources_dct['KCl'] = ImpactItem.get_item('KCl_item').CFs['H_Resources'] = i

    # MgOH2
    b = GWP_dct['MgOH2']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='MgOH2 CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg', baseline=b, distribution=D)
    def set_MgOH2_CF(i):
        GWP_dct['MgOH2'] = ImpactItem.get_item('MgOH2_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['MgOH2']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='MgOH2 ecosystems CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_MgOH2_ecosystems_CF(i):
        H_Ecosystems_dct['MgOH2'] = ImpactItem.get_item('MgOH2_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['MgOH2']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='MgOH2 health CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_MgOH2_health_CF(i):
        H_Health_dct['MgOH2'] = ImpactItem.get_item('MgOH2_item').CFs['H_Health'] = i

    b = H_Resources_dct['MgOH2']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='MgOH2 resources CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_MgOH2_resources_CF(i):
        H_Resources_dct['MgOH2'] = ImpactItem.get_item('MgOH2_item').CFs['H_Resources'] = i

    # Other CFs
    item_path = os.path.join(re_data_path, 'impact_items.xlsx')

    for indicator in ('GlobalWarming', 'H_Ecosystems', 'H_Health', 'H_Resources'):
        sheet_name = indicator if indicator != 'GlobalWarming' else 'GWP'
        data = load_data(item_path, sheet=sheet_name)
        for p in data.index:
            item = ImpactItem.get_item(p)
            b = item.CFs[indicator]
            lower = float(data.loc[p]['low'])
            upper = float(data.loc[p]['high'])
            dist = data.loc[p]['distribution']
            if dist == 'uniform':
                D = shape.Uniform(lower=lower, upper=upper)
            elif dist == 'triangular':
                D = shape.Triangle(lower=lower, midpoint=b, upper=upper)
            elif dist == 'constant':
                continue
            else:
                raise ValueError(f'Distribution {dist} not recognized.')
            model.parameter(name=p + f'-{indicator}',
                            setter=DictAttrSetter(item, 'CFs', indicator),
                            element='LCA', kind='isolated',
                            units=f'kg CO2-eq/{item.functional_unit}',
                            baseline=b, distribution=D)


# %%

# =============================================================================
# Functions to create models
# =============================================================================

# System A: Solids removal only
def create_modelA(country_specific=False, **model_kwargs):
    sysA = create_system('A')
    unitA = sysA.flowsheet.unit

    # Shared metrics/parameters
    modelA = Model(sysA, **model_kwargs)
    add_metrics(modelA)
    unit_dctA = {
        'Excretion': unitA.A1,
        'Primary': unitA.A3,
        'Pasteurization': unitA.A4,
    }
    add_shared_parameters(modelA, unit_dctA, country_specific)

    return modelA


# System B: Full Duke system with grid electricity source
def create_modelB(country_specific=False, **model_kwargs):
    sysB = create_system('B')
    unitB = sysB.flowsheet.unit

    # Shared parameters
    modelB = Model(sysB, **model_kwargs)
    add_metrics(modelB)
    unit_dctB = {
        'Excretion': unitB.B1,
        'Toilet': unitB.B2,
        'Primary': unitB.B3,
        'Pasteurization': unitB.B4,
        'Ultrafiltration': unitB.B5,
        'Ion exchange': unitB.B6,
        'ECR': unitB.B7,
        'Housing': unitB.B10,
        'System': unitB.B11,
    }
    add_shared_parameters(modelB, unit_dctB, country_specific)

    return modelB


# System C: Full Duke system with solar electricity source
def create_modelC(country_specific=False, **model_kwargs):
    sysC = create_system('C')
    unitC = sysC.flowsheet.unit

    # Shared parameters
    modelC = Model(sysC, **model_kwargs)
    add_metrics(modelC)
    unit_dctC = {
        'Excretion': unitC.C1,
        'Primary': unitC.C3,
        'Pasteurization': unitC.C4,
        'Ultrafiltration': unitC.C5,
        'Ion exchange': unitC.C6,
        'ECR': unitC.C7,
        'Housing': unitC.C10,
        'System': unitC.C11,
        'Solar': unitC.C12,
    }
    add_shared_parameters(modelC, unit_dctC, country_specific)

    return modelC


# System D: Targeted nitrogen removal (created for NSS preliminary analysis)
def create_modelD(country_specific=False, **model_kwargs):
    sysD = create_system('D')
    unitD = sysD.flowsheet.unit

    # Shared parameters
    modelD = Model(sysD, **model_kwargs)
    add_metrics(modelD)
    unit_dctD = {
        'Excretion': unitD.D1,
        'Primary': unitD.D3,
        'Ultrafiltration': unitD.D4,
        'Ion exchange': unitD.D5,
        'Housing': unitD.D8,
        'System': unitD.D9,
    }
    add_shared_parameters(modelD, unit_dctD, country_specific)

    return modelD


# Wrapper function so that it'd work for all
def create_model(model_ID='A', country_specific=False, **model_kwargs):
    model_ID = model_ID.lower().rsplit('model')[-1].rsplit('sys')[-1].upper() # works for "modelA"/"sysA"/"A"
    if model_ID == 'A': model = create_modelA(country_specific, **model_kwargs)
    elif model_ID == 'B': model = create_modelB(country_specific, **model_kwargs)
    elif model_ID == 'C': model = create_modelC(country_specific, **model_kwargs)
    elif model_ID == 'D': model = create_modelD(country_specific, **model_kwargs)
    else: raise ValueError(f'`model_ID` can only be "A", "B", "C", or "D", not "{model_ID}".')
    return model


def run_uncertainty(model, path='', **kwargs):
    kwargs['path'] = os.path.join(results_path, f'sys{model.system.ID[-1]}_model.xlsx') if path=='' else path
    run(model=model, **kwargs)
    return