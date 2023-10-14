#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Shion Watabe <shionwatabe@gmail.com>

    Hannah Lohman <hlohman94@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import os, qsdsan as qs
from chaospy import distributions as shape
from qsdsan import Model, Metric, PowerUtility, ImpactItem
from qsdsan.utils import (
    AttrSetter,
    data_path,
    DictAttrSetter,
    load_data,
    )
from exposan.utils import batch_setting_unit_params, run_uncertainty as run
from exposan import new_generator as ng
from exposan.new_generator import (
    create_system,
    data_path as ng_data_path,
    default_ppl,
    get_decay_k,
    get_LCA_metrics,
    get_TEA_metrics,
    get_recoveries,
    results_path,
    update_resource_recovery_settings,
    )

__all__ = ('create_model', 'run_uncertainty',)


# %%

# =============================================================================
# Functions for batch-making metrics and setting parameters
# =============================================================================

def add_metrics(model, ppl=default_ppl):
    ng._load_lca_data()
    system = model.system

    # Recoveries
    funcs = get_recoveries(system, ppl)
    metrics = [
        Metric('Total N', funcs[0], '% N', 'N recovery'),
        Metric('Total P', funcs[1], '% P', 'P recovery'),
        Metric('Total K', funcs[2], '% K', 'K recovery'),
    ]
    # Net cost
    metrics.append(
        Metric('Annual net cost', get_TEA_metrics(system, ppl)[0], f'{qs.currency}/cap/yr', 'TEA results'),
        )
    # Net emissions
    funcs = get_LCA_metrics(system, ppl)
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
ng_su_data_path = os.path.join(su_data_path, 'ng')

def load_ng_su_data(file_name):
    if file_name.startswith('_ng'):
        return load_data(os.path.join(ng_su_data_path, file_name))
    return load_data(os.path.join(su_data_path, file_name))

excretion_data = load_ng_su_data('_excretion.tsv')
murt_data = load_ng_su_data('_murt.tsv')
pasteurization_data = load_ng_su_data('_sludge_pasteurization.tsv')

anmbr_data = load_ng_su_data('_ng_anmbr.csv')
ix_data = load_ng_su_data('_ng_ion_exchange.tsv')
chlorination_data = load_ng_su_data('_ng_chlorination.tsv')
pv_data = load_ng_su_data('_ng_photovoltaic.tsv')
grid_data = load_ng_su_data('_ng_gridtied.tsv')
controls_data = load_ng_su_data('_ng_controls.tsv')
foundation_data = load_ng_su_data('_ng_foundation.csv')
housing_data = load_ng_su_data('_ng_housing.csv')
pretreatment_data = load_ng_su_data('_ng_pretreatment.csv')


# %%

# =============================================================================
# Shared by all systems
# =============================================================================

def add_shared_parameters(model, unit_dct, country_specific=False):
    sys = model.system
    sys_stream = sys.flowsheet.stream
    param = model.parameter
    price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct = update_resource_recovery_settings()

    # Add these parameters if not running country-specific analysis,
    # in which they would be updated separately

    excretion_unit = unit_dct['excretion']

    if not country_specific:
        # Price ratio
        old_price_dct = price_dct.copy()
        b = 1
        D = shape.Uniform(lower=b-(10**(-6)), upper=b+(10**(-6)))
        stream_ref = {
            'GAC': 'GAC',
            'zeolite': 'zeolite',
            'NaOH': 'NaOH',
            'NaCl': 'NaCl',
            'NaCl1': 'NaCl1',
        }
        @param(name='Price ratio', element=excretion_unit, kind='cost', units='-',
               baseline=b, distribution=D)
        def set_price_ratio(i):
            ng.price_ratio = i
            for obj_name in stream_ref.keys():
                old_price = old_price_dct[obj_name]
                new_price = old_price * i
                getattr(sys_stream, stream_ref[obj_name]).price = new_price
            for u in sys.units:
                if hasattr(u, 'price_ratio'):
                    u.price_ratio = i

        # Labor wage
        b = price_dct['wages']
        D = shape.Triangle(lower=b*0.5, midpoint=b, upper=b*1.5)
        @param(name='Labor wages', element='TEA', kind='cost', units='USD/h',
               baseline=b, distribution=D)
        def set_labor_wages(i):
            for u in sys.units:
                if hasattr(u, '_calc_maintenance_labor_cost'):
                    u.wages = i

        # Electricity price
        b = price_dct['Electricity']
        D = shape.Triangle(lower=0.08, midpoint=b, upper=0.14)
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

        if ng.INCLUDE_RESOURCE_RECOVERY:
            # N fertilizer price
            b = 1.507
            D = shape.Uniform(lower=b*0.8, upper=b*1.2)
            @param(name='N fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
                    baseline=b, distribution=D)
            def set_N_price(i):
                price_dct['N'] = sys_stream.liq_N.price = sys_stream.sol_N.price = i * ng.price_factor
    
            # P fertilizer price
            b = 3.983
            D = shape.Uniform(lower=b*0.8, upper=b*1.2)
            @param(name='P fertilizer price', element='TEA', kind='isolated', units='USD/kg P',
                   baseline=b, distribution=D)
            def set_P_price(i):
                price_dct['P'] = sys_stream.liq_P.price = sys_stream.sol_P.price = i * ng.price_factor
    
            # K fertilizer price
            b = 1.333
            D = shape.Uniform(lower=b*0.8, upper=b*1.2)
            @param(name='K fertilizer price', element='TEA', kind='isolated', units='USD/kg K',
                   baseline=b, distribution=D)
            def set_K_price(i):
                price_dct['K'] = sys_stream.liq_K.price = sys_stream.sol_K.price = i * ng.price_factor
    
            # NH3 fertilizer price
            D = shape.Uniform(lower=(1.507*(14/17)*0.8), upper=(1.507*(14/17)*1.2))
            @param(name='NH3 fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
                   baseline=(1.507*(14/17)), distribution=D)
            def set_con_NH3_price(i):
                price_dct['conc_NH3'] = sys_stream.conc_NH3.price = i * ng.price_factor

        # NaCl price
        b = price_dct['NaCl']
        D = shape.Uniform(lower=0.20685, upper=0.34475)
        @param(name='NaCl price', element='TEA', kind='isolated', units='USD/kg',
               baseline=b, distribution=D)
        def set_NaCl_price(i):
            price_dct['NaCl'] = sys_stream.NaCl.price = sys_stream.NaCl1.price = i

        # LPG price
        b = price_dct['LPG']
        D = shape.Uniform(lower=1.0661, upper=1.9799)
        @param(name='LPG price', element='TEA', kind='isolated', units='USD/kg',
               baseline=b, distribution=D)
        def set_LPG_price(i):
            price_dct['LPG'] = sys_stream.LPG.price = i

    ##### Specific units #####
    param = model.parameter

    # Diet and Excretion
    exclude = ('e_cal', 'p_anim', 'p_veg') if country_specific else ()
    batch_setting_unit_params(excretion_data, model, excretion_unit, exclude)

    # MURT Toilet
    toilet_unit = unit_dct['toilet']
    exclude = ('MCF_decay', 'N2O_EF_decay', 'OPEX_over_CAPEX')
    batch_setting_unit_params(murt_data, model, toilet_unit, exclude)

    b = toilet_unit.OPEX_over_CAPEX
    D = shape.Uniform(lower=0.02, upper=0.08)
    param(name='MURT operating cost', element=toilet_unit, kind='coupled', units='cost',
          baseline=b, distribution=D)
    def set_OPEX_over_CAPEX(i):
        toilet_unit.OPEX_over_CAPEX = i

    b = toilet_unit.MCF_decay
    D = shape.Triangle(lower=0.05, midpoint=b, upper=0.15)
    param(name='MCF_decay', element=toilet_unit, kind='coupled',
          units='fraction of anaerobic conversion of degraded COD',
          baseline=b, distribution=D)
    def set_toilet_MCF_decay(i):
        toilet_unit.MCF_decay = i

    b = toilet_unit.N2O_EF_decay
    D = shape.Triangle(lower=0, midpoint=b, upper=0.001)
    param(name='N2O_EF_decay', element=toilet_unit, kind='coupled',
          units='fraction of N emitted as N2O',
          baseline=b, distribution=D)
    def set_toilet_N2O_EF_decay(i):
        toilet_unit.N2O_EF_decay = i

    # AnMBR
    anmbr_unit = unit_dct['anmbr']
    exclude = ('wages', 'MCF_decay', 'N2O_EF_decay')
    batch_setting_unit_params(anmbr_data, model, anmbr_unit, exclude)

    b = anmbr_unit.MCF_decay
    D = shape.Triangle(lower=0.05, midpoint=b, upper=0.15)
    param(name='MCF_decay', element=anmbr_unit, kind='coupled',
          units='fraction of anaerobic conversion of degraded COD',
          baseline=b, distribution=D)
    def set_anmbr_MCF_decay(i):
        anmbr_unit.MCF_decay = i

    b = anmbr_unit.N2O_EF_decay
    D = shape.Triangle(lower=0, midpoint=b, upper=0.001)
    param(name='N2O_EF_decay', element=anmbr_unit, kind='coupled',
          units='fraction of N emitted as N2O',
          baseline=b, distribution=D)
    def set_anmbr_N2O_EF_decay(i):
        anmbr_unit.N2O_EF_decay = i

    b = anmbr_unit.sludge_moisture_content
    D = shape.Uniform(lower=0.9, upper=0.96)
    param(name='Sludge moisture content', element=anmbr_unit, kind='coupled', units='fraction',
          baseline=b, distribution=D)
    def set_anmbr_sludge_moisture_content(i):
        anmbr_unit.sludge_moisture_content = i

    # Sludge pasteurization
    pasteurization_unit = unit_dct['pasteurization']
    exclude = ('biogas_loss', 'biogas_eff', 'heat_loss', 'sludge_temp', 'wages')
    batch_setting_unit_params(pasteurization_data, model, pasteurization_unit, exclude)

    b = pasteurization_unit.biogas_loss
    D = shape.Uniform(lower=0.07, upper=0.13)
    param(setter=AttrSetter(pasteurization_unit, 'biogas_loss'),
          name='Fraction biogas lost', element=pasteurization_unit, kind='coupled', units='fraction',
          baseline=b, distribution=D)

    b = pasteurization_unit.biogas_eff
    D = shape.Uniform(lower=0.385, upper=0.715)
    param(setter=AttrSetter(pasteurization_unit, 'biogas_eff'),
          name='Biogas efficiency', element=pasteurization_unit, kind='coupled', units='fraction',
          baseline=b, distribution=D)

    b = pasteurization_unit.heat_loss
    D = shape.Uniform(lower=0.075, upper=0.125)
    param(setter=AttrSetter(pasteurization_unit, 'heat_loss'),
          name='Heat loss', element=pasteurization_unit, kind='coupled', units='fraction',
          baseline=b, distribution=D)

    b = pasteurization_unit.sludge_temp
    D = shape.Uniform(lower=212.36, upper=353.94)
    param(setter=AttrSetter(pasteurization_unit, 'sludge_temp'),
          name='Sludge temperature', element=pasteurization_unit, kind='coupled', units='fraction',
          baseline=b, distribution=D)

    # Ion exchange
    ix_unit = unit_dct['ion_exchange']
    exclude = ('wages',)
    batch_setting_unit_params(ix_data, model, ix_unit, exclude)

    # Chlorination
    chlorination_unit = unit_dct['chlorination']
    batch_setting_unit_params(chlorination_data, model, chlorination_unit, exclude)

    # Controls
    controls_unit = unit_dct['control_system']
    batch_setting_unit_params(controls_data, model, controls_unit, exclude)

    # Housing
    housing_unit = unit_dct['housing']
    batch_setting_unit_params(housing_data, model, housing_unit)

    # Pretreatment
    pretreatment_unit = unit_dct['pretreatment']
    batch_setting_unit_params(pretreatment_data, model, pretreatment_unit, exclude)

    # Foundation
    foundation_unit = unit_dct['foundation']
    batch_setting_unit_params(foundation_data, model, foundation_unit)

    ##### Universal degradation parameters #####
    # Max methane emission
    unit = sys.path[1]  # the first unit that involves degradation
    b = ng.max_CH4_emission
    D = shape.Triangle(lower=0.175, midpoint=b, upper=0.325)
    @param(name='Max CH4 emission', element=unit, kind='coupled', units='g CH4/g COD',
           baseline=b, distribution=D)
    def set_max_CH4_emission(i):
        ng.max_CH4_emission = i
        for unit in sys.units:
            if hasattr(unit, 'max_CH4_emission'):
                setattr(unit, 'max_CH4_emission', i)

    # Time to full degradation
    b = ng.tau_deg
    D = shape.Uniform(lower=1, upper=3)
    @param(name='Full degradation time', element=unit, kind='coupled', units='yr',
           baseline=b, distribution=D)
    def set_tau_deg(i):
        ng.tau_deg = i
        k = get_decay_k(i, ng.log_deg)
        for unit in sys.units:
            if hasattr(unit, 'decay_k_COD'):
                setattr(unit, 'decay_k_COD', k)
            if hasattr(unit, 'decay_k_N'):
                setattr(unit, 'decay_k_N', k)

    # Reduction at full degradation
    b = ng.log_deg
    D = shape.Uniform(lower=2, upper=4)
    @param(name='Log degradation', element=unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_log_deg(i):
        ng.log_deg = i
        k = get_decay_k(ng.tau_deg, i)
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
    b = ng.price_factor
    D = shape.Uniform(lower=0.1, upper=0.4)
    @param(name='Price factor', element='TEA', kind='isolated', units='-',
           baseline=b, distribution=D)
    def set_price_factor(i):
        ng.price_factor = i

    b = price_dct['zeolite']
    D = shape.Uniform(lower=0.20, upper=0.27)
    @param(name='Zeolite price', element='TEA', kind='isolated', units='USD/kg',
           baseline=b, distribution=D)
    def set_zeolite_price(i):
        price_dct['zeolite'] = sys_stream.zeolite.price = i

    b = price_dct['GAC']
    D = shape.Uniform(lower=0.83, upper=1.38)
    @param(name='GAC price', element='TEA', kind='isolated', units='USD/kg',
           baseline=b, distribution=D)
    def set_GAC_price(i):
        price_dct['GAC'] = sys_stream.GAC.price = i

    b = price_dct['NaOH']
    D = shape.Uniform(lower=0.12, upper=0.20)
    @param(name='NaOH price', element='TEA', kind='isolated', units='USD/kg',
           baseline=b, distribution=D)
    def set_NaOH_price(i):
        price_dct['NaOH'] = sys_stream.NaOH.price = i

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

    if ng.INCLUDE_RESOURCE_RECOVERY:
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
        D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
        @param(name='K fertilizer ecosystems CF', element='LCA', kind='isolated',
               units='points/kg K', baseline=b, distribution=D)
        def set_K_fertilizer_ecosystems_CF(i):
            H_Ecosystems_dct['K'] = ImpactItem.get_item('K_item').CFs['H_Ecosystems'] = -i

        b = -H_Health_dct['K']
        D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
        @param(name='K fertilizer health CF', element='LCA', kind='isolated',
               units='points/kg K', baseline=b, distribution=D)
        def set_K_fertilizer_health_CF(i):
            H_Health_dct['K'] = ImpactItem.get_item('K_item').CFs['H_Health'] = -i

        b = -H_Resources_dct['K']
        D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
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
    @param(name='Zeolite CF', element='LCA', kind='isolated',
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
    @param(name='zeolite resources CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_zeolite_resources_CF(i):
        H_Resources_dct['zeolite'] = ImpactItem.get_item('zeolite_item').CFs['H_Resources'] = i

    # NaCl
    b = GWP_dct['NaCl']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='NaCl CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg', baseline=b, distribution=D)
    def set_NaCl_CF(i):
        GWP_dct['NaCl'] = ImpactItem.get_item('NaCl_item').CFs['GlobalWarming'] = \
            ImpactItem.get_item('NaCl1_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['NaCl']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='NaCl ecosystems CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_NaCl_ecosystems_CF(i):
        H_Ecosystems_dct['NaCl'] = ImpactItem.get_item('NaCl_item').CFs['H_Ecosystems'] = \
            ImpactItem.get_item('NaCl1_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['NaCl']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='NaCl health CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_NaCl_health_CF(i):
        H_Health_dct['NaCl'] = ImpactItem.get_item('NaCl_item').CFs['H_Health'] = \
            ImpactItem.get_item('NaCl1_item').CFs['H_Health'] = i

    b = H_Resources_dct['NaCl']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='NaCl resources CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_NaCl_resources_CF(i):
        H_Resources_dct['NaCl'] = ImpactItem.get_item('NaCl_item').CFs['H_Resources'] = \
            ImpactItem.get_item('NaCl1_item').CFs['H_Resources'] = i

    # NaOH
    b = GWP_dct['NaOH']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='NaOH CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg', baseline=b, distribution=D)
    def set_NaOH_CF(i):
        GWP_dct['NaOH'] = ImpactItem.get_item('NaOH_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['NaOH']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='NaOH ecosystems CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_NaOH_ecosystems_CF(i):
        H_Ecosystems_dct['NaOH'] = ImpactItem.get_item('NaOH_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['NaOH']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='NaOH health CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_NaOH_health_CF(i):
        H_Health_dct['NaOH'] = ImpactItem.get_item('NaOH_item').CFs['H_Health'] = i

    b = H_Resources_dct['NaOH']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='NaOH resources CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_NaOH_resources_CF(i):
        H_Resources_dct['NaOH'] = ImpactItem.get_item('NaOH_item').CFs['H_Resources'] = i

    # LPG
    b = GWP_dct['LPG']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='LPG CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg', baseline=b, distribution=D)
    def set_LPG_CF(i):
        GWP_dct['LPG'] = ImpactItem.get_item('LPG_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['LPG']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='LPG ecosystems CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_LPG_ecosystems_CF(i):
        H_Ecosystems_dct['LPG'] = ImpactItem.get_item('LPG_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['LPG']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='LPG health CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_LPG_health_CF(i):
        H_Health_dct['LPG'] = ImpactItem.get_item('LPG_item').CFs['H_Health'] = i

    b = H_Resources_dct['LPG']
    D = shape.Uniform(lower=b*0.9, upper=b*1.1)
    @param(name='LPG resources CF', element='LCA', kind='isolated',
           units='points/kg', baseline=b, distribution=D)
    def set_LPG_resources_CF(i):
        H_Resources_dct['LPG'] = ImpactItem.get_item('LPG_item').CFs['H_Resources'] = i

    # Other CFs
    item_path = os.path.join(ng_data_path, 'impact_items.xlsx')

    for indicator in ('GlobalWarming', 'H_Ecosystems', 'H_Health', 'H_Resources'):
        sheet_name = indicator if indicator!='GlobalWarming' else 'GWP'
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
            elif dist == 'constant': continue
            else:
                raise ValueError(f'Distribution {dist} not recognized.')
            model.parameter(name=p+f'-{indicator}',
                            setter=DictAttrSetter(item, 'CFs', indicator),
                            element='LCA', kind='isolated',
                            units=f'kg CO2-eq/{item.functional_unit}',
                            baseline=b, distribution=D)


# %%

# =============================================================================
# Functions to create models
# =============================================================================

# System A: Baseline Photovoltaic System
def create_modelA(country_specific=False, ppl=default_ppl, **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysA = create_system('A', ppl=ppl, flowsheet=flowsheet)
    unitA = sysA.flowsheet.unit

    # Shared metrics/parameters
    modelA = Model(sysA, **model_kwargs)
    add_metrics(modelA, ppl=ppl)
    unit_dctA = {
        'excretion': unitA.A1,
        'toilet': unitA.A2,
        'anmbr': unitA.A3,
        'pasteurization': unitA.A4,
        'ion_exchange': unitA.A5,
        'chlorination': unitA.A6,
        'control_system': unitA.A8,
        'housing': unitA.A9,
        'pretreatment': unitA.A10,
        'foundation': unitA.A11,
        }
    add_shared_parameters(modelA, unit_dctA, country_specific)

    # Photovoltaic energy
    exclude = ('wages',)
    batch_setting_unit_params(pv_data, modelA, unitA.A7, exclude)

    return modelA


# System B: Alternative Grid-tied System
def create_modelB(country_specific=False, ppl=default_ppl, **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysB = create_system('B', ppl=ppl, flowsheet=flowsheet)
    unitB = sysB.flowsheet.unit

    # Shared parameters
    modelB = Model(sysB, **model_kwargs)
    add_metrics(modelB, ppl=ppl)
    unit_dctB = {
        'excretion': unitB.B1,
        'toilet': unitB.B2,
        'anmbr': unitB.B3,
        'pasteurization': unitB.B4,
        'ion_exchange': unitB.B5,
        'chlorination': unitB.B6,
        'control_system': unitB.B8,
        'housing': unitB.B9,
        'pretreatment': unitB.B10,
        'foundation': unitB.B11,
        }
    add_shared_parameters(modelB, unit_dctB, country_specific)

    # Grid-tied energy
    exclude = ('wages',)
    batch_setting_unit_params(grid_data, modelB, unitB.B7, exclude)

    return modelB


# Wrapper function so that it'd work for all
def create_model(model_ID='A', country_specific=False, ppl=default_ppl, **model_kwargs):
    model_ID = model_ID.lower().rsplit('model')[-1].rsplit('sys')[-1].upper() # works for "modelA"/"sysA"/"A"
    if model_ID == 'A': model = create_modelA(country_specific, ppl=ppl, **model_kwargs)
    elif model_ID == 'B': model = create_modelB(country_specific, ppl=ppl, **model_kwargs)
    else: raise ValueError(f'`model_ID` can only be "A" or "B", not "{model_ID}".')
    return model


def run_uncertainty(model, path='', **kwargs):
    kwargs['path'] = os.path.join(results_path, f'sys{model.system.ID[-1]}_model.xlsx') if path=='' else path
    run(model=model, **kwargs)
    return