#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

#%%
import os, qsdsan as qs
from datetime import datetime
from chaospy import distributions as shape
from qsdsan import Model, Metric, PowerUtility, ImpactItem, StreamImpactItem
from qsdsan.utils import (
    AttrSetter,
    data_path,
    DictAttrSetter,
    load_data,
    )
from exposan.utils import (
    batch_setting_unit_params, 
    run_uncertainty as run,
    get_generic_tanker_truck_fee as get_tanker_truck_fee,
    get_decay_k,
    get_generic_scaled_capital,
    )
from exposan import g2rt
from exposan.g2rt import (
    create_system,
    g2rt_data_path,
    default_ppl,
    default_lifetime,
    get_decay_k,
    get_LCA_metrics,
    get_TEA_metrics,
    get_normalized_CAPEX,
    get_normalized_electricity_cost,
    get_normalized_OPEX,
    get_recoveries,
    results_path,
    update_resource_recovery_settings,
    get_normalized_labor_cost,
    get_unit_contruction_GW_impact,
    get_unit_stream_GW_impact,
    get_unit_electrcitiy_GW_impact,
    get_unit_total_impact,
    compute_unit_total_cost,
    get_normalized_recovery_earning
    )
import warnings

__all__ = ('create_model', 'run_uncertainty',)


#%%

# =============================================================================
# Functions for batch-making metrics and setting parameters
# =============================================================================

def add_metrics(model, ppl=default_ppl):
    g2rt._load_lca_data()
    system = model.system

    # Recoveries
    funcs = get_recoveries(system, ppl)
    metrics = [
        # Metric('Total N', funcs[0], '% N', 'N recovery'),
        # Metric('Total P', funcs[1], '% P', 'P recovery'),
        # Metric('Total K', funcs[2], '% K', 'K recovery'),
        Metric('Total Water', funcs[3], '% Water', 'Water recovery'),
    ]
    
    # Net cost
    metrics.append(
        Metric('Annual net cost', get_TEA_metrics(system, ppl)[0], f'{qs.currency}/cap/yr', 'TEA results'),
        )
    metrics.extend([
        Metric('Total CAPEX', get_TEA_metrics(system, ppl,include_breakdown=True)[1], f'{qs.currency}', 'TEA results'),
        Metric('Annual OPEX', get_TEA_metrics(system, ppl,include_breakdown=True)[7], f'{qs.currency}/yr','TEA results'),
        Metric('Construction emissions', get_LCA_metrics(system, ppl,include_breakdown=True)[4], 'kg CO2-eq', 'LCA results'),
        Metric('Annual operating emissions', get_LCA_metrics(system, ppl,include_breakdown=True)[8], 'kg CO2-eq/yr', 'LCA results')
        ])
    
    # metrics.extend([
    #     Metric('Energy consumption', get_TEA_metrics(system, ppl,include_breakdown=True)[2], 'kWh /cap/day', 'TEA results')
    #     ])
    # for u in system.TEA.units:
    #     class_name = u.__class__.__name__
    #     metrics.extend([
    #         Metric(f'{class_name} CAPEX', get_normalized_CAPEX(u, ppl), f'{qs.currency}/cap/day', f'{class_name} TEA'),
    #         Metric(f'{class_name} Electricity cost', get_normalized_electricity_cost(u, ppl), f'{qs.currency}/cap/day', f'{class_name} TEA'),
    #         Metric(f'{class_name} OPEX_no_labor_electricity', 
    #                get_normalized_OPEX(u,ppl), f'{qs.currency}/cap/day', f'{class_name} TEA'),
    #         Metric(f'{class_name} Labor', 
    #                get_normalized_labor_cost(u, ppl), f'{qs.currency}/cap/day', f'{class_name} TEA'),
    #         Metric(f'{class_name} total cost', 
    #                compute_unit_total_cost(u, ppl), f'{qs.currency}/cap/day', f'{class_name} TEA'),])

    # for u in system.TEA.units:
    #     class_name = u.__class__.__name__
    #     metrics.extend([
    #         Metric(f'{class_name} capital GW', get_unit_contruction_GW_impact(u,ppl),'kg CO2-eq/cap/yr', f'{class_name} LCA'),
    #         Metric(f'{class_name} stream GW', get_unit_stream_GW_impact(u,ppl),'kg CO2-eq/cap/yr', f'{class_name} LCA'),
    #         Metric(f'{class_name} electricity GW', get_unit_electrcitiy_GW_impact(u,ppl),'kg CO2-eq/cap/yr', f'{class_name} LCA'),
    #         Metric(f'{class_name} total GW', get_unit_total_impact(u,ppl),'kg CO2-eq/cap/yr', f'{class_name} LCA'),
    #     ])

    # direct_emission_units = {'A15','A16','A17','B14','B15','B16'}
    # for u in system.flowsheet.unit:
    #     if u.ID in direct_emission_units:
    #         class_name = u.__class__.__name__
    #         metrics.append(
    #             Metric(f'{class_name} {u.ID} direct emission GW', get_unit_stream_GW_impact(u,ppl),'kg CO2-eq/cap/yr', f'{class_name} LCA'),
    #         )
    # if g2rt.INCLUDE_RESOURCE_RECOVERY:
    #     #TODO: separate LCA because Mixer, A14, A15, A16 need to be included...
    #     recovery_units = {'A6','A18'}
    #     for u in system.flowsheet.unit:
    #         if u.ID in recovery_units:
    #             class_name = u.__class__.__name__
    #             metrics.extend([
    #                 Metric(f'{class_name} {u.ID} recovery earning', get_normalized_recovery_earning(u, ppl), f'{qs.currency}/cap/day', f'Recovery {class_name} TEA'),
    #                 Metric(f'{class_name} {u.ID} reduced GW', get_unit_stream_GW_impact(u,ppl),'kg CO2-eq/cap/yr', f'Recovery {class_name} LCA'),
    #             ])
    
    # Net emissions
    funcs = get_LCA_metrics(system, ppl)
    metrics.extend([
        Metric('GlobalWarming', funcs[0], 'kg CO2-eq/cap/yr', 'LCA results'),
        # Metric('H_Ecosystems', funcs[1], 'points/cap/yr', cat),
        # Metric('H_Health', funcs[2], 'points/cap/yr', cat),
        # Metric('H_Resources', funcs[3], 'points/cap/yr', cat),
        ])
    model.metrics = metrics

# %%

# =============================================================================
# Data sheets
# =============================================================================

su_data_path = os.path.join(data_path, 'sanunit_data')
g2rt_su_data_path = os.path.join(su_data_path, 'g2rt')

def load_g2rt_su_data(file_name):
    if file_name.startswith(('_g2rt_', '_vr_','_mscwo_')):
        return load_data(os.path.join(g2rt_su_data_path, file_name))
    else: 
        return load_data(os.path.join(su_data_path, file_name))

excretion_path = load_g2rt_su_data('_excretion.tsv')
murt_path = load_g2rt_su_data('_murt.tsv')
mixer_toilet_path = load_g2rt_su_data('_g2rt_toilet.csv')
vr_pasteurization_path = load_g2rt_su_data('_vr_pasteurization.csv')
G2RT_homogenizer_path = load_g2rt_su_data('_g2rt_homogenizer.csv')
vr_filter_press_path = load_g2rt_su_data('_vr_filter_press.csv')
vr_concentrator_path = load_g2rt_su_data('_vr_concentrator.csv')
g2rt_liquids_tank_path = load_g2rt_su_data('_g2rt_liquids_tank.csv')
g2rt_solids_tank_path = load_g2rt_su_data('_g2rt_solids_tank.csv')
vr_drying_tunnel_path = load_g2rt_su_data('_vr_dryingtunnel.csv')
g2rt_solids_separation_path = load_g2rt_su_data('_g2rt_solids_separation.csv')
g2rt_belt_separation_path = load_g2rt_su_data('_g2rt_belt_separation.csv')
ultrafiltration_path = load_g2rt_su_data('_g2rt_ultrafiltration.csv')
reverse_osmosis_path = load_g2rt_su_data('_g2rt_reverse_osmosis.csv')
combustor_path = load_g2rt_su_data('_vr_combustor.csv')
g2rt_housing_path = load_g2rt_su_data('_g2rt_housing.csv')
g2rt_controls_path = load_g2rt_su_data('_g2rt_controls.csv')
mscwo_concentrator_path = load_g2rt_su_data('_mscwo_concentrator_module.csv')
mscwo_gas_handling_module_path = load_g2rt_su_data('_mscwo_gas_handling_module.csv')
mscwo_reactor_module_path = load_g2rt_su_data('_mscwo_reactor_module.csv')

# %%

# =============================================================================
# Shared by all systems
# =============================================================================

def add_shared_parameters(model, unit_dct, e_CF=None, flush_water=None, city_specific=False):
    sys = model.system
    sys_stream = sys.flowsheet.stream
    param = model.parameter
    price_dct, GWP_dct, H_Ecosystems_dct, H_Health_dct, H_Resources_dct = update_resource_recovery_settings()

    # Add these parameters if not running context-specific analysis,
    # in which they would be updated separately
    
    excretion_unit = unit_dct['excretion']

    if not city_specific:
        # Price ratio
        old_price_dct = price_dct.copy()
        b = 1
        D = shape.Uniform(lower=b-(10**(-6)), upper=b+(10**(-6)))
        stream_ref = {
            'H2O': 'H2O',
        }
        @param(name='Price ratio', element=excretion_unit, kind='cost', units='-', #TODO: why not element = 'TEA'? 
               baseline=b, distribution=D)
        
        def set_price_ratio(i):
            price_ratio = i
            for obj_name, dic_key in stream_ref.items():
                old_price = old_price_dct[dic_key]
                new_price = old_price * i
                getattr(sys_stream, obj_name).price = new_price
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
        D = shape.Triangle(lower=0.06, midpoint=b, upper=0.20)
        @param(name='Electricity price', element='TEA', kind='isolated',
               units='$/kWh', baseline=b, distribution=D)
        def set_electricity_price(i):
            PowerUtility.price = i

        # Electricity GWP
        if e_CF is not None:
            GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = e_CF
        else:
            b = GWP_dct['Electricity']
            D = shape.Triangle(lower=b*0.7, midpoint=b, upper=b*1.3)
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

        if g2rt.INCLUDE_RESOURCE_RECOVERY:
            # N fertilizer price
            b = 1.507
            D = shape.Uniform(lower=b*0.8, upper=b*1.2)
            @param(name='N fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
                    baseline=b, distribution=D)
            def set_N_price(i):
                price_dct['N'] = sys_stream.liq_N.price = sys_stream.sol_N.price = i * g2rt.price_factor
    
            # P fertilizer price
            b = 3.983
            D = shape.Uniform(lower=b*0.8, upper=b*1.2)
            @param(name='P fertilizer price', element='TEA', kind='isolated', units='USD/kg P',
                   baseline=b, distribution=D)
            def set_P_price(i):
                price_dct['P'] = sys_stream.liq_P.price = sys_stream.sol_P.price = i * g2rt.price_factor
    
            # K fertilizer price
            b = 1.333
            D = shape.Uniform(lower=b*0.8, upper=b*1.2)
            @param(name='K fertilizer price', element='TEA', kind='isolated', units='USD/kg K',
                   baseline=b, distribution=D)
            def set_K_price(i):
                price_dct['K'] = sys_stream.liq_K.price = sys_stream.sol_K.price = i * g2rt.price_factor
                
            # H2O price
            b = 4.01/1e3
            D = shape.Uniform(lower=1.59/1e3, upper=8.84/1e3) #consider inflation and data from 2018 https://webassets.bv.com/2019-10/50_Largest_Cities_Rate_Survey_2018_2019_Report.pdf
            @param(name='H2O price', element='TEA', kind='isolated', units='USD/kg H2O',
                   baseline= b, distribution=D)
            def set_H2O_price(i):
                price_dct['H2O'] = sys_stream.H2O.price = i * g2rt.price_factor
    ##### Specific units #####
    param = model.parameter

    # Diet and Excretion
    exclude = ('e_cal', 'p_anim', 'p_veg') if city_specific else ()
    batch_setting_unit_params(excretion_path, model, excretion_unit, exclude)
    
    # MURT Toilet
    toilet_unit = unit_dct['toilet']
    exclude = ('OPEX_over_CAPEX')
    batch_setting_unit_params(murt_path, model, toilet_unit, exclude)
    
    b = toilet_unit.OPEX_over_CAPEX
    D = shape.Uniform(lower=0.02, upper=0.08)
    param(name='MURT operating cost', element=toilet_unit, kind='coupled', units='cost',
          baseline=b, distribution=D)
    def set_OPEX_over_CAPEX(i):
        toilet_unit.OPEX_over_CAPEX = i
    
    #Mixer for flushing   
    mixer_unit = unit_dct['mixer']
    if flush_water is not None:
        exclude = ('flushing_water')
        batch_setting_unit_params(mixer_toilet_path, model, mixer_unit,exclude)
    else: 
        batch_setting_unit_params(mixer_toilet_path, model, mixer_unit)
    
    #Solids separation
    solids_separation_unit = unit_dct['solids_separation']
    exclude = ('wages') #exclude because it is a context or city-specific parameter
    batch_setting_unit_params(g2rt_solids_separation_path, model, solids_separation_unit, exclude)
    
    #Belt separation
    belt_separation_unit = unit_dct['belt_separation']
    exclude = ('wages') 
    batch_setting_unit_params(g2rt_belt_separation_path, model, belt_separation_unit, exclude)
    
    #Solids tank
    solids_tank_unit = unit_dct['solids_tank']
    exclude = ('wages')
    batch_setting_unit_params(g2rt_solids_tank_path, model, solids_tank_unit, exclude)
    
    #Liquids tank
    liquids_tank_unit = unit_dct['liquids_tank']
    exclude = ('wages')
    batch_setting_unit_params(g2rt_liquids_tank_path, model, liquids_tank_unit, exclude)
    
    #Ultrafiltration
    ultrafiltration_unit = unit_dct['ultrafiltration']
    exclude = ('wages')
    batch_setting_unit_params(ultrafiltration_path, model, ultrafiltration_unit, exclude)
    
    #Reverse osmosis
    reverse_osmosis_unit = unit_dct['reverse_osmosis']
    exclude = ('wages')
    batch_setting_unit_params(reverse_osmosis_path, model, reverse_osmosis_unit, exclude)
    
    #Homogenizer
    homogenizer_unit = unit_dct['homogenizer']
    exclude = ('wages')
    batch_setting_unit_params(G2RT_homogenizer_path, model, homogenizer_unit, exclude)
    
    #Drying tunnel
    drying_tunnel_unit = unit_dct['drying_tunnel']
    exclude = ('wages')
    batch_setting_unit_params(vr_drying_tunnel_path, model, drying_tunnel_unit, exclude)
    
    #Controls
    controls_unit = unit_dct['controls']
    exclude = ('wages')
    batch_setting_unit_params(g2rt_controls_path, model, controls_unit, exclude)
    
    #Housing
    housing_unit = unit_dct['housing']
    exclude = ('wages')
    batch_setting_unit_params(g2rt_housing_path, model, housing_unit, exclude)

    ##### Universal degradation parameters #####
    # Max methane emission
    unit = sys.path[1]  # the first unit that involves degradation
    b = g2rt.max_CH4_emission
    D = shape.Triangle(lower=0.175, midpoint=b, upper=0.325)
    @param(name='Max CH4 emission', element=unit, kind='coupled', units='g CH4/g COD',
           baseline=b, distribution=D)
    def set_max_CH4_emission(i):
        g2rt.max_CH4_emission = i
        for unit in sys.units:
            if hasattr(unit, 'max_CH4_emission'):
                setattr(unit, 'max_CH4_emission', i)

    # Time to full degradation
    b = g2rt.tau_deg
    D = shape.Uniform(lower=1, upper=3)
    @param(name='Full degradation time', element=unit, kind='coupled', units='yr',
           baseline=b, distribution=D)
    def set_tau_deg(i):
        g2rt.tau_deg = i
        k = get_decay_k(i, g2rt.log_deg)
        for unit in sys.units:
            if hasattr(unit, 'decay_k_COD'):
                setattr(unit, 'decay_k_COD', k)
            if hasattr(unit, 'decay_k_N'):
                setattr(unit, 'decay_k_N', k)

    # Reduction at full degradation
    b = g2rt.log_deg
    D = shape.Uniform(lower=2, upper=4)
    @param(name='Log degradation', element=unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_log_deg(i):
        g2rt.log_deg = i
        k = get_decay_k(g2rt.tau_deg, i)
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
    b = g2rt.price_factor
    D = shape.Uniform(lower=0.1, upper=0.4)
    @param(name='Price factor', element='TEA', kind='isolated', units='-',
           baseline=b, distribution=D)
    def set_price_factor(i):
        g2rt.price_factor = i

    if g2rt.INCLUDE_RESOURCE_RECOVERY:
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

        # Recovered H2O
        b = -GWP_dct['H2O']
        D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
        @param(name='H2O CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kg', baseline=b, distribution=D)
        def set_H2O_CF(i):
            GWP_dct['H2O'] = ImpactItem.get_item('H2O_item').CFs['GlobalWarming'] = -i

        b = -H_Ecosystems_dct['H2O']
        D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
        @param(name='H2O ecosystems CF', element='LCA', kind='isolated',
                units='points/kg', baseline=b, distribution=D)
        def set_H2O_ecosystems_CF(i):
            H_Ecosystems_dct['H2O'] = ImpactItem.get_item('H2O_item').CFs['H_Ecosystems'] = -i

        b = -H_Health_dct['H2O']
        D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
        @param(name='H2O health CF', element='LCA', kind='isolated',
                units='points/kg', baseline=b, distribution=D)
        def set_H2O_health_CF(i):
            H_Health_dct['H2O'] = ImpactItem.get_item('H2O_item').CFs['H_Health'] = -i

        b = -H_Resources_dct['H2O']
        D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
        @param(name='H2O resources CF', element='LCA', kind='isolated',
                units='points/kg', baseline=b, distribution=D)
        def set_H2O_resources_CF(i):
            H_Resources_dct['H2O'] = ImpactItem.get_item('H2O_item').CFs['H_Resources'] = -i

    # Other CFs
    item_path = os.path.join(g2rt_data_path, 'impact_items.xlsx')
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
# System A: volume reduction toilet
def create_modelA(city_specific=False, ppl=default_ppl, lifetime=default_lifetime,
                  flush_water= None, combustion_CH4_EF= None,e_CF=None,
                  **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysA = create_system('A', ppl=ppl, lifetime=lifetime, flowsheet=flowsheet, 
                         flush_water =flush_water, combustion_CH4_EF=combustion_CH4_EF)
    unitA = sysA.flowsheet.unit
    if flush_water is not None:
        if not (isinstance(flush_water, float) or 0 <= flush_water <= 4.2):
            warnings.warn(
                "Warning: flush_water should be a float between 0 and 4.2 kg/cap/hr "
                "will be ignored. Please provide another proper value if intended.",
                UserWarning
            )
            flush_water = None
    if combustion_CH4_EF is not None:
        if not (isinstance(combustion_CH4_EF, float) or 0 <= combustion_CH4_EF <= 0.1):
            warnings.warn(
                "Warning: combustion_CH4_EF should be a float between 0 and 0.1 g-CH4/g-feces solids"
                "will be ignored. Please provide another proper value if intended.",
                UserWarning
            )
            combustion_CH4_EF = None
    if e_CF is not None:
        if not (isinstance(e_CF, float) or 0 <= e_CF <= 2):
            warnings.warn(
                "Warning: e_CF should be a float (between 0 and 2) "
                "will be ignored. Please provide another proper value if intended.",
                UserWarning
            )
            e_CF = None
    # Shared metrics/parameters
    modelA = Model(sysA, **model_kwargs)
    add_metrics(modelA, ppl=ppl)
    unit_dctA = {
        'excretion': unitA.A1,
        'mixer':unitA.Mixer,
        'toilet': unitA.A2,
        'solids_separation': unitA.A3,
        'belt_separation': unitA.A4,
        'solids_tank': unitA.A12,
        'liquids_tank': unitA.A13,
        'ultrafiltration': unitA.A5,
        'reverse_osmosis': unitA.A6,
        'homogenizer': unitA.A8,
        'drying_tunnel': unitA.A11,
        'controls': unitA.A19,
        'housing': unitA.A20
        }
    add_shared_parameters(model=modelA, unit_dct=unit_dctA, e_CF=e_CF, 
                          flush_water=flush_water, city_specific=city_specific)
    
    #Add parameters unique to System A: volume reduction toilet
    #Concentrator
    exclude = ('wages')
    batch_setting_unit_params(vr_concentrator_path, modelA, unitA.A7, exclude)
    
    #Pasteurization
    exclude = ('wages')
    batch_setting_unit_params(vr_pasteurization_path, modelA, unitA.A9, exclude)
    
    #Filter press
    exclude = ('wages')
    batch_setting_unit_params(vr_filter_press_path, modelA, unitA.A10, exclude)

    #Combustor
    if combustion_CH4_EF is not None:
        exclude = ('wages','CH4_emission_factor') #remove CH4_emission_factor and to compute with a range of fixed values
    else:
        exclude = ('wages')
    batch_setting_unit_params(combustor_path, modelA, unitA.A14, exclude)
    
    return modelA

#System B: micro supercritical water oxidation toilet
def create_modelB(city_specific=False, ppl=default_ppl,lifetime=default_lifetime, 
                  mscwo_replacement_cost =None, e_CF=None, flush_water= None, **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysB = create_system('B', ppl=ppl,lifetime=lifetime, flowsheet=flowsheet, 
                         mscwo_replacement_cost=mscwo_replacement_cost, flush_water =flush_water)
    unitB = sysB.flowsheet.unit
    if mscwo_replacement_cost is not None:
        if not (isinstance(mscwo_replacement_cost, float) or 0 <= mscwo_replacement_cost <= 1):
            warnings.warn(
                "Warning: mscwo_replacement_cost should be a fraction (between 0 and 1) "
                "will be ignored. Please provide another proper value if intended.",
                UserWarning
            )
            mscwo_replacement_cost = None
    if e_CF is not None:
        if not (isinstance(e_CF, float) or  0 <= e_CF <= 2):
            warnings.warn(
                "Warning: e_CF should be a float (between 0 and 2) "
                "will be ignored. Please provide another proper value if intended.",
                UserWarning
            )
            e_CF = None
    if flush_water is not None:
        if not (isinstance(flush_water, float) or  0 <= flush_water <= 4.2):
            warnings.warn(
                "Warning: flush_water should be a float between 0 and 4.2 kg/cap/hr "
                "will be ignored. Please provide another proper value if intended.",
                UserWarning
            )
            flush_water = None
    # Shared parameters
    modelB = Model(sysB, **model_kwargs)
    add_metrics(modelB, ppl=ppl)
    unit_dctB = {
        'excretion': unitB.B1,
        'mixer':unitB.Mixer,
        'toilet': unitB.B2,
        'solids_separation': unitB.B3,
        'belt_separation': unitB.B4,
        'solids_tank': unitB.B5,
        'homogenizer': unitB.B6,
        'liquids_tank': unitB.B7,
        'ultrafiltration': unitB.B8,
        'reverse_osmosis': unitB.B9,
        'drying_tunnel': unitB.B13,
        'controls': unitB.B18,
        'housing': unitB.B19
        }
    add_shared_parameters(model=modelB, unit_dct=unit_dctB, e_CF=e_CF, 
                          flush_water=flush_water, city_specific=city_specific)
    #Add parameters unique to System B: micro supercritical water oxidation toilet
    #Gas handling module
    exclude = ('wages')
    batch_setting_unit_params(mscwo_gas_handling_module_path, modelB, unitB.B10, exclude)
    
    #mSCWO reactor module
    if mscwo_replacement_cost is not None:
        exclude = ('wages','material_replacement_cost') #remove material_replacement_cost and to compute with a range of fixed values
    else:
        exclude = ('wages')
    batch_setting_unit_params(mscwo_reactor_module_path, modelB, unitB.B11, exclude)

    #mSCWO concentrator module
    exclude = ('wages')
    batch_setting_unit_params(mscwo_concentrator_path, modelB, unitB.B12, exclude)

    return modelB

# Wrapper function so that it'd work for all
def create_model(model_ID='A', city_specific=False, ppl=default_ppl,lifetime=default_lifetime, **model_kwargs):
    model_ID = model_ID.lower().rsplit('model')[-1].rsplit('sys')[-1].upper() # works for "modelA"/"sysA"/"A"
    if model_ID == 'A': model = create_modelA(city_specific, ppl=ppl,lifetime= lifetime, **model_kwargs)
    elif model_ID == 'B': model = create_modelB(city_specific, ppl=ppl,lifetime= lifetime, **model_kwargs)
    else: raise ValueError(f'`model_ID` can only be "A" or "B", not "{model_ID}".')
    return model

def run_uncertainty(model, path='', date = None, note = '', **kwargs):
    # Use current date if date not provided
    date = date or datetime.now().strftime('%Y-%m-%d')
    file_name = f'sys{model.system.ID[-1]}_model_{date}'
    if note:
        file_name += f'_{note}'
    file_name += '.xlsx'
    kwargs['path'] = os.path.join(results_path,file_name) if path=='' else path
    run(model=model, **kwargs)
    return