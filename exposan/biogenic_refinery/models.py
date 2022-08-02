#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Lane To <lane20@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


# %%

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
from exposan import biogenic_refinery as br
from exposan.biogenic_refinery import (
    _load_components,
    create_system,
    data_path as br_data_path,
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
# Functions for batch-making metrics and -setting parameters
# =============================================================================

def add_metrics(model):
    br._load_lca_data()
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
br_su_data_path = os.path.join(su_data_path, 'br')

def load_br_su_data(file_name):
    if file_name.startswith('_br'):
        return load_data(os.path.join(br_su_data_path, file_name))
    return load_data(os.path.join(su_data_path, file_name))

excretion_data = load_br_su_data('_excretion.tsv')
toilet_data = load_br_su_data('_toilet.tsv')
pit_latrine_data = load_br_su_data('_pit_latrine.tsv')
uddt_data = load_br_su_data('_uddt.tsv')
all_uddt_data = pd.concat((toilet_data, uddt_data))
anaerobic_lagoon_data = load_br_su_data('_anaerobic_lagoon.tsv')
drying_bed_data = load_br_su_data('_drying_bed.tsv')
liquid_bed_data = load_br_su_data('_liquid_treatment_bed.tsv')

controls_data = load_br_su_data('_br_controls.tsv')
housing_data = load_br_su_data('_br_housing.tsv')
# Screw press for sysA/sysC, grinder for sysB
screw_data = load_br_su_data('_br_screw_press.tsv')
grinder_data = load_br_su_data('_br_grinder.tsv')
carbonizer_data = load_br_su_data('_br_carbonizer_base.tsv')
pollution_control_data = load_br_su_data('_br_pollution_control.tsv')
ohx_data = load_br_su_data('_br_ohx.tsv')
hhx_data = load_br_su_data('_br_hhx.tsv')
hhx_dryer_data = load_br_su_data('_br_hhx_dryer.tsv')
struvite_data = load_br_su_data('_br_struvite_precipitation.tsv')
ix_data = load_br_su_data('_br_ion_exchange.tsv')


# %%

# =============================================================================
# Shared by systems A-C
# =============================================================================

def add_shared_parameters(model, unit_dct, country_specific=False):
    sys = model.system
    sys_stream = sys.flowsheet.stream
    param = model.parameter

    # Add these parameters if not running country-specific analysis,
    # in which they would be updated separately
    if not country_specific:
        # Household size
        excretion_unit = unit_dct['Excretion']
        b = br.household_size
        D = shape.Normal(mu=b, sigma=1.8)
        @param(name='Household size', element=excretion_unit, kind='coupled', units='cap/household',
               baseline=b, distribution=D)
        def set_household_size(i):
            br.household_size = max(1, i)

        # Operator labor wage
        b = br.operator_daily_wage
        D = shape.Triangle(lower=(14.55), midpoint=b, upper=(43.68))
        @param(name='Operator daily wages', element='TEA', kind='cost', units='USD/d',
              baseline=b, distribution=D)
        def set_operator_daily_wage(i):
            sys._TEA.annual_labor = i*3*365

        # Construction labor wage
        b = br.const_daily_wage
        D = shape.Triangle(lower=(b*0.5), midpoint=b, upper=(b*1.5))
        @param(name='Construction daily wages', element='TEA', kind='cost', units='USD/d',
              baseline=b, distribution=D)
        def set_const_daily_wage(i):
            for u in sys.units:
                if isinstance(u, qs.sanunits.BiogenicRefineryHousing): break
                u.const_daily_wage = i

        # N fertilizer price
        b = 1.507
        D = shape.Uniform(lower=b*0.8, upper=b*1.2)
        @param(name='N fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
                baseline=b, distribution=D)
        def set_N_price(i):
            price_dct['N'] = sys_stream.liq_N.price = sys_stream.sol_N.price = i * br.price_factor

        # P fertilizer price
        b = 3.983
        D = shape.Uniform(lower=b*0.8, upper=b*1.2)
        @param(name='P fertilizer price', element='TEA', kind='isolated', units='USD/kg P',
               baseline=b, distribution=D)
        def set_P_price(i):
            price_dct['P'] = sys_stream.liq_P.price = sys_stream.sol_P.price = i * br.price_factor

        # K fertilizer price
        b = 1.333
        D = shape.Uniform(lower=b*0.8, upper=b*1.2)
        @param(name='K fertilizer price', element='TEA', kind='isolated', units='USD/kg K',
               baseline=b, distribution=D)
        def set_K_price(i):
            price_dct['K'] = sys_stream.liq_K.price = sys_stream.sol_K.price = i * br.price_factor

        # NH3 fertilizer price
        D = shape.Uniform(lower=(1.507*(14/17)*0.8), upper=(1.507*(14/17)*1.2))
        @param(name='NH3 fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
               baseline=(1.507*(14/17)), distribution=D)
        def set_con_NH3_price(i):
            price_dct['conc_NH3'] = sys_stream.conc_NH3.price = i * br.price_factor

        # Struvite fertilizer price
        D = shape.Uniform(lower=(3.983*(31/245)*0.8), upper=(3.983*(31/245)*1.2))
        @param(name='Struvite fertilizer price', element='TEA', kind='isolated', units='USD/kg P',
               baseline=(3.983*(31/245)), distribution=D)
        def set_struvite_price(i):
            price_dct['struvite'] = sys_stream.struvite.price = i * br.price_factor

        # Commented out because not taking into account economic value of biochar
        # D = shape.Uniform(lower=(0.014*0.8), upper=(0.014*1.2))
        # @param(name='Biochar  price', element='TEA', kind='isolated', units='USD/kg biochar',
        #         baseline=(0.014), distribution=D)
        # def set_biochar_price(i):
        #     price_dct['biochar'] = sys_stream.biochar.price = i * br.price_factor

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

    ##### Specific units #####
    # Diet and excretion
    excretion_unit = unit_dct['Excretion']
    exclude = ('e_cal','p_anim','p_veg',) if country_specific else ()
    batch_setting_unit_params(excretion_data, model, excretion_unit, exclude)

    # Toilet
    toilet_unit = unit_dct['Toilet']
    b = br.household_per_toilet
    D = shape.Uniform(lower=3, upper=5)
    @param(name='Toilet density', element=toilet_unit, kind='coupled', units='household/toilet',
            baseline=b, distribution=D)
    def set_toilet_density(i):
        br.household_per_toilet = i

    # Control box (industrial control panel)
    control_unit = unit_dct['ControlBox']
    exclude = ('certified_electrician_wages', 'service_team_wages', 'facility_manager_wages', 'biomass_controls_wages',) if country_specific else ()
    batch_setting_unit_params(controls_data, model, control_unit, exclude)

    # Housing biogenic refinery
    housing_unit = unit_dct['Housing']
    batch_setting_unit_params(housing_data, model, housing_unit)

    # Screw press (sysA/sysC) or grinder (sysB)
    solid_decomp_unit = unit_dct.get('Screw', None)
    if solid_decomp_unit:
        batch_setting_unit_params(screw_data, model, solid_decomp_unit)
    else:
        batch_setting_unit_params(grinder_data, model, unit_dct['Grinder'])

    # Liquid treatment bed
    liquid_bed_unit = unit_dct['LiquidBed']
    batch_setting_unit_params(liquid_bed_data, model, liquid_bed_unit)

    # Carbonizer base
    carbonizer_unit = unit_dct['Carbonizer']
    exclude = ('service_team_wages',) if country_specific else ()
    batch_setting_unit_params(carbonizer_data, model, carbonizer_unit, exclude)

    # Pollution control device
    pollution_control_unit = unit_dct['PCD']
    exclude = ('service_team_wages',) if country_specific else ()
    batch_setting_unit_params(pollution_control_data, model, pollution_control_unit, exclude)

    # Oil heat exchanger
    ohx_unit = unit_dct['OilHX']
    batch_setting_unit_params(ohx_data, model, ohx_unit)

    # Hydronic heat exchanger
    hhx_unit = unit_dct['HHX']
    exclude = ('service_team_wages',) if country_specific else ()
    batch_setting_unit_params(hhx_data, model, hhx_unit, exclude)

    # Dryer from HHX
    hhx_dryer_unit = unit_dct['HHXdryer']
    batch_setting_unit_params(hhx_dryer_data, model, hhx_dryer_unit)

    ##### Universal degradation parameters #####
    # Max methane emission
    toilet_unit = sys.path[1] # the first unit that involves degradation
    b = br.max_CH4_emission
    D = shape.Triangle(lower=0.175, midpoint=b, upper=0.325)
    @param(name='Max CH4 emission', element=toilet_unit, kind='coupled', units='g CH4/g COD',
           baseline=b, distribution=D)
    def set_max_CH4_emission(i):
        br.max_CH4_emission = i
        for unit in sys.units:
            if hasattr(unit, 'max_CH4_emission'):
                setattr(unit, 'max_CH4_emission', i)

    # Time to full degradation
    b = br.tau_deg
    D = shape.Uniform(lower=1, upper=3)
    @param(name='Full degradation time', element=toilet_unit, kind='coupled', units='yr',
           baseline=b, distribution=D)
    def set_tau_deg(i):
        br.tau_deg = i
        k = get_decay_k(i, br.log_deg)
        for unit in sys.units:
            if hasattr(unit, 'decay_k_COD'):
                setattr(unit, 'decay_k_COD', k)
            if hasattr(unit, 'decay_k_N'):
                setattr(unit, 'decay_k_N', k)

    # Reduction at full degradation
    b = br.log_deg
    D = shape.Uniform(lower=2, upper=4)
    @param(name='Log degradation', element=toilet_unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_log_deg(i):
        br.log_deg = i
        k = get_decay_k(br.tau_deg, i)
        for unit in sys.units:
            if hasattr(unit, 'decay_k_COD'):
                setattr(unit, 'decay_k_COD', k)
            if hasattr(unit, 'decay_k_N'):
                setattr(unit, 'decay_k_N', k)

    ##### Toilet material properties #####
    density = toilet_unit.density_dct
    b = density['Plastic']
    D = shape.Uniform(lower=0.31, upper=1.24)
    param(setter=DictAttrSetter(toilet_unit, 'density_dct', 'Plastic'),
          name='Plastic density', element=toilet_unit, kind='isolated', units='kg/m2',
          baseline=b, distribution=D)

    b = density['Brick']
    D = shape.Uniform(lower=1500, upper=2000)
    param(setter=DictAttrSetter(toilet_unit, 'density_dct', 'Brick'),
          name='Brick density', element=toilet_unit, kind='isolated', units='kg/m3',
          baseline=b, distribution=D)

    b = density['StainlessSteelSheet']
    D = shape.Uniform(lower=2.26, upper=3.58)
    param(setter=DictAttrSetter(toilet_unit, 'density_dct', 'StainlessSteelSheet'),
          name='SS sheet density', element=toilet_unit, kind='isolated', units='kg/m2',
          baseline=b, distribution=D)

    b = density['Gravel']
    D = shape.Uniform(lower=1520, upper=1680)
    param(setter=DictAttrSetter(toilet_unit, 'density_dct', 'Gravel'),
          name='Gravel density', element=toilet_unit, kind='isolated', units='kg/m3',
          baseline=b, distribution=D)

    b = density['Sand']
    D = shape.Uniform(lower=1281, upper=1602)
    param(setter=DictAttrSetter(toilet_unit, 'density_dct', 'Sand'),
          name='Sand density', element=toilet_unit, kind='isolated', units='kg/m3',
          baseline=b, distribution=D)

    b = density['Steel']
    D = shape.Uniform(lower=7750, upper=8050)
    param(setter=DictAttrSetter(toilet_unit, 'density_dct', 'Steel'),
          name='Steel density', element=toilet_unit, kind='isolated', units='kg/m3',
          baseline=b, distribution=D)

    ##### General TEA settings #####
    # # Keeping discount rate constant
    # b = br.discount_rate
    # D = shape.Uniform(lower=0.03, upper=0.06)
    # @param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
    #         baseline=b, distribution=D)
    # def set_discount_rate(i):
    #     br.discount_rate = i

    # Discount factor for the excreta-derived fertilizers
    b = br.price_factor
    D = shape.Uniform(lower=0.1, upper=0.4)
    @param(name='Price factor', element='TEA', kind='isolated', units='-',
           baseline=b, distribution=D)
    def set_price_factor(i):
        br.price_factor = i

    b = price_dct['Polymer']
    D = shape.Uniform(lower=(b*0.9), upper=(b*1.1))
    @param(name='Polymer price', element='TEA', kind='isolated', units='USD/kg polymer',
            baseline=(b), distribution=D)
    def set_polymber_price(i):
        price_dct['Polymer'] = sys_stream.polymer.price = i

    b = price_dct['Resin']
    D = shape.Uniform(lower=(b*0.9), upper=(b*1.1))
    @param(name='Resin price', element='TEA', kind='isolated', units='USD/kg resin',
            baseline=(b), distribution=D)
    def set_resin_price(i):
        price_dct['Resin'] = sys_stream.resin.price = i

    b = price_dct['FilterBag']
    D = shape.Uniform(lower=(b*0.9), upper=(b*1.1))
    @param(name='Filter bag price', element='TEA', kind='isolated', units='USD/kg filter bag',
            baseline=(b), distribution=D)
    def set_filterbag_price(i):
        price_dct['FilterBag'] = sys_stream.filter_bag.price = i

    ##### General LCA settings #####
    # TODO: there's probably a more elegant way of adding these CFs
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

    # N20
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

    # Polymer
    b = GWP_dct['Polymer']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Polymer CF', element='LCA', kind='isolated',
            units='kg CO2-eq/kg N', baseline=b, distribution=D)
    def set_polymer_CF(i):
        GWP_dct['Polymer'] = ImpactItem.get_item('polymer_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['Polymer']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Polymer ecosystems CF', element='LCA', kind='isolated',
            units='points/kg N', baseline=b, distribution=D)
    def set_polymer_ecosystems_CF(i):
        H_Ecosystems_dct['Polymer'] = ImpactItem.get_item('polymer_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['Polymer']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Polymer health CF', element='LCA', kind='isolated',
            units='points/kg N', baseline=b, distribution=D)
    def set_polymer_health_CF(i):
        H_Health_dct['Polymer'] = ImpactItem.get_item('polymer_item').CFs['H_Health'] = i

    b = H_Resources_dct['Polymer']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Polymer resources CF', element='LCA', kind='isolated',
            units='points/kg N', baseline=b, distribution=D)
    def set_polymer_resources_CF(i):
        H_Resources_dct['Polymer'] = ImpactItem.get_item('polymer_item').CFs['H_Resources'] = i

    # Resin
    b = GWP_dct['Resin']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Resin CF', element='LCA', kind='isolated',
            units='kg CO2-eq/kg N', baseline=b, distribution=D)
    def set_resin_CF(i):
        GWP_dct['Resin'] = ImpactItem.get_item('resin_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['Resin']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Resin ecosystems CF', element='LCA', kind='isolated',
            units='points/kg N', baseline=b, distribution=D)
    def set_resin_ecosystems_CF(i):
        H_Ecosystems_dct['Resin'] = ImpactItem.get_item('resin_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['Resin']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Resin health CF', element='LCA', kind='isolated',
            units='points/kg N', baseline=b, distribution=D)
    def set_resin_health_CF(i):
        H_Health_dct['Resin'] = ImpactItem.get_item('resin_item').CFs['H_Health'] = i

    b = H_Resources_dct['Resin']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Resin resources CF', element='LCA', kind='isolated',
            units='points/kg N', baseline=b, distribution=D)
    def set_resin_resources_CF(i):
        H_Resources_dct['Resin'] = ImpactItem.get_item('resin_item').CFs['H_Resources'] = i

    # Filter Bag
    b = GWP_dct['FilterBag']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Filter bag CF', element='LCA', kind='isolated',
            units='kg CO2-eq/kg N', baseline=b, distribution=D)
    def set_filter_bag_CF(i):
        GWP_dct['FilterBag'] = ImpactItem.get_item('filter_bag_item').CFs['GlobalWarming'] = i

    b = H_Ecosystems_dct['FilterBag']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Filter bag ecosystems CF', element='LCA', kind='isolated',
            units='points/kg N', baseline=b, distribution=D)
    def set_filter_bag_ecosystems_CF(i):
        H_Ecosystems_dct['FilterBag'] = ImpactItem.get_item('filter_bag_item').CFs['H_Ecosystems'] = i

    b = H_Health_dct['FilterBag']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Filter bag health CF', element='LCA', kind='isolated',
            units='points/kg N', baseline=b, distribution=D)
    def set_filter_bag_health_CF(i):
        H_Health_dct['FilterBag'] = ImpactItem.get_item('filter_bag_item').CFs['H_Health'] = i

    b = H_Resources_dct['FilterBag']
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='Filter bag resources CF', element='LCA', kind='isolated',
            units='points/kg N', baseline=b, distribution=D)
    def set_filter_bag_resources_CF(i):
        H_Resources_dct['FilterBag'] = ImpactItem.get_item('filter_bag_item').CFs['H_Resources'] = i

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
    b = 5.4 * (14 / 17)
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='conc_NH3 fertilizer CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg conc_NH3', baseline=b, distribution=D)
    def set_conc_NH3_fertilizer_CF(i):
        GWP_dct['conc_NH3'] = ImpactItem.get_item('conc_NH3_item').CFs['GlobalWarming'] = -i

    b = 0.0461961 * (14 / 17)
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='conc_NH3 fertilizer ecosystems CF', element='LCA', kind='isolated',
           units='points/kg conc_NH3', baseline=b, distribution=D)
    def set_conc_NH3_fertilizer_ecosystems_CF(i):
        H_Ecosystems_dct['conc_NH3'] = ImpactItem.get_item('conc_NH3_item').CFs['H_Ecosystems'] = -i

    b = 0.637826734 * (14 / 17)
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='conc_NH3 fertilizer health CF', element='LCA', kind='isolated',
           units='points/kg conc_NH3', baseline=b, distribution=D)
    def set_conc_NH3_fertilizer_health_CF(i):
        H_Health_dct['conc_NH3'] = ImpactItem.get_item('conc_NH3_item').CFs['H_Health'] = -i

    b = 0.259196888 * (14 / 17)
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='conc_NH3 fertilizer resources CF', element='LCA', kind='isolated',
           units='points/kg conc_NH3', baseline=b, distribution=D)
    def set_conc_NH3_fertilizer_resources_CF(i):
        H_Resources_dct['conc_NH3'] = ImpactItem.get_item('conc_NH3_item').CFs['H_Resources'] = -i

    # Recovered struvite fertilizer
    b = 4.9 * (31 / 245)
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='struvite fertilizer CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg struvite', baseline=b, distribution=D)
    def set_struvite_fertilizer_CF(i):
        GWP_dct['struvite'] = ImpactItem.get_item('struvite_item').CFs['GlobalWarming'] = -i

    b = 0.093269908 * (31 / 245)
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='struvite fertilizer ecosystems CF', element='LCA', kind='isolated',
           units='points/kg struvite', baseline=b, distribution=D)
    def set_struvite_fertilizer_ecosystems_CF(i):
        H_Ecosystems_dct['struvite'] = ImpactItem.get_item('struvite_item').CFs['H_Ecosystems'] = -i

    b = 1.774294425 * (31 / 245)
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='struvite fertilizer health CF', element='LCA', kind='isolated',
           units='points/kg struvite', baseline=b, distribution=D)
    def set_struvite_fertilizer_health_CF(i):
        H_Health_dct['struvite'] = ImpactItem.get_item('struvite_item').CFs['H_Health'] = -i

    b = 1.084191599 * (31 / 245)
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='struvite fertilizer resources CF', element='LCA', kind='isolated',
           units='points/kg struvite', baseline=b, distribution=D)
    def set_struvite_fertilizer_resources_CF(i):
        H_Resources_dct['struvite'] = ImpactItem.get_item('struvite_item').CFs['H_Resources'] = -i

    # Recovered biochar
    b = 0.2 * 0.9 * (44 / 12)  # assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='biochar CF', element='LCA', kind='isolated',
           units='kg CO2-eq/kg biochar', baseline=b, distribution=D)
    def set_biochar_CF(i):
        GWP_dct['biochar'] = ImpactItem.get_item('biochar_item').CFs['GlobalWarming'] = -i

    b = 0.2 * 0.9 * (44 / 12) * br.EcosystemQuality_factor  # assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='biochar ecosystems CF', element='LCA', kind='isolated',
           units='points/kg biochar', baseline=b, distribution=D)
    def set_biochar_ecosystems_CF(i):
        H_Ecosystems_dct['biochar'] = ImpactItem.get_item('biochar_item').CFs['H_Ecosystems'] = -i

    b = 0.2 * 0.9 * (44 / 12) * br.HumanHealth_factor  # assume biochar 20% by mass is fixed C with 90% of that being stable (44/12) carbon to CO2
    D = shape.Triangle(lower=b*0.90, midpoint=b, upper=b*1.1)
    @param(name='biochar health CF', element='LCA', kind='isolated',
           units='points/kg biochar', baseline=b, distribution=D)
    def set_biochar_health_CF(i):
        H_Health_dct['biochar'] = ImpactItem.get_item('biochar_item').CFs['H_Health'] = -i

    # Other CFs
    item_path = os.path.join(br_data_path, 'impact_items.xlsx')

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
# Shared by systems A, C, and D
# =============================================================================

MCF_lower_dct = dct_from_str(pit_latrine_data.loc['MCF_decay']['low'])
MCF_upper_dct = dct_from_str(pit_latrine_data.loc['MCF_decay']['high'])
N2O_EF_lower_dct = dct_from_str(pit_latrine_data.loc['N2O_EF_decay']['low'])
N2O_EF_upper_dct = dct_from_str(pit_latrine_data.loc['N2O_EF_decay']['high'])

def add_pit_latrine_parameters(model, unit_dct):
    pit_unit = unit_dct['Toilet']
    param = model.parameter
    ##### Related to the toilet #####
    all_pit_data = pd.concat((toilet_data, pit_latrine_data))
    batch_setting_unit_params(all_pit_data, model, pit_unit, exclude=('MCF_decay', 'N2O_EF_decay',))

    kind = pit_unit._return_MCF_EF()
    if model.system.ID != 'sysD': # Parameters not applicable to sysD
        ##### Related to conveyance #####
        convey_unit = unit_dct['Conveyance']
        b = convey_unit.loss_ratio
        D = shape.Uniform(lower=0.02, upper=0.05)
        param(setter=AttrSetter(convey_unit, 'loss_ratio'),
              name='Transportation loss', element=convey_unit, kind='coupled', units='fraction',
              baseline=b, distribution=D)

        b = convey_unit.single_truck.distance
        D = shape.Uniform(lower=2, upper=10)
        param(setter=AttrSetter(convey_unit.single_truck, 'distance'),
              name='Transportation distance', element=convey_unit, kind='coupled', units='km',
              baseline=b, distribution=D)

        b = br.emptying_fee
        D = shape.Uniform(lower=0, upper=0.3)
        @param(name='Emptying fee', element=convey_unit, kind='coupled', units='USD',
                baseline=b, distribution=D)
        def set_emptying_fee(i):
            br.emptying_fee = i

    # Parameters applicable to all other systems with the pit latrine
    b = pit_unit.MCF_decay
    D = shape.Triangle(lower=MCF_lower_dct[kind], midpoint=b, upper=MCF_upper_dct[kind])
    param(setter=DictAttrSetter(pit_unit, '_MCF_decay', kind),
          name='MCF_decay', element=pit_unit, kind='coupled',
          units='fraction of anaerobic conversion of degraded COD',
          baseline=b, distribution=D)

    b = pit_unit.N2O_EF_decay
    D = shape.Triangle(lower=N2O_EF_lower_dct[kind], midpoint=b, upper=N2O_EF_upper_dct[kind])
    param(setter=DictAttrSetter(pit_unit, '_N2O_EF_decay', kind),
          name='N2O_EF_decay', element=pit_unit, kind='coupled',
          units='fraction of N emitted as N2O',
          baseline=b, distribution=D)

    # Costs
    b = pit_unit.CAPEX
    D = shape.Uniform(lower=523, upper=673)
    param(setter=AttrSetter(pit_unit, 'CAPEX'),
          name='Pit latrine capital cost', element=pit_unit, kind='cost',
          units='USD', baseline=b, distribution=D)

    b = pit_unit.OPEX_over_CAPEX
    D = shape.Uniform(lower=0.02, upper=0.08)
    param(setter=AttrSetter(pit_unit, 'OPEX_over_CAPEX'),
          name='Pit latrine operating cost', element=pit_unit, kind='cost',
          units='fraction of capital cost', baseline=b, distribution=D)


# %%

# =============================================================================
# Functions to create models
# =============================================================================

def create_modelA(country_specific=False, **model_kwargs):
    sysA = create_system('A')
    unitA = sysA.flowsheet.unit

    # Shared metrics/parameters
    modelA = Model(sysA, **model_kwargs)
    add_metrics(modelA)
    unit_dctA = {
        'Excretion': unitA.A1,
        'Toilet': unitA.A2,
        'Conveyance': unitA.A3,
        'ControlBox': unitA.A4, # industrial control panel
        'Housing': unitA.A5,
        'Screw': unitA.A6,
        'LiquidBed': unitA.A7,
        'Carbonizer': unitA.A8,
        'PCD': unitA.A9, # pollution control device
        'OilHX': unitA.A10,
        'HHX': unitA.A11,
        'HHXdryer': unitA.A12,
        }
    add_shared_parameters(modelA, unit_dctA, country_specific)

    # Pit latrine and conveyance
    add_pit_latrine_parameters(modelA, unit_dctA)

    return modelA


def create_modelB(country_specific=False, **model_kwargs):
    sysB = create_system('B')
    unitB = sysB.flowsheet.unit

    # Shared parameters
    modelB = Model(sysB, **model_kwargs)
    add_metrics(modelB)
    unit_dctB = {
        'Excretion': unitB.B1,
        'Toilet': unitB.B2,
        'Conveyance': (unitB.B3, unitB.B4),
        'LiquidBed': unitB.B7,
        'ControlBox': unitB.B8, # industrial control panel
        'Housing': unitB.B9,
        'Grinder': unitB.B10,
        'Carbonizer': unitB.B11,
        'PCD': unitB.B12, # pollution control device
        'OilHX': unitB.B13,
        'HHX': unitB.B14,
        'HHXdryer': unitB.B15,
        }
    add_shared_parameters(modelB, unit_dctB, country_specific)

    # UDDT
    B2 = unitB.B2
    batch_setting_unit_params(all_uddt_data, modelB, B2)

    paramB = modelB.parameter
    b = B2.CAPEX
    D = shape.Uniform(lower=571, upper=756)
    @paramB(name='UDDT capital cost', element=B2, kind='cost',
            units='USD', baseline=b, distribution=D)
    def set_UDDT_CAPEX(i):
        B2.CAPEX = i

    b = B2.OPEX_over_CAPEX
    D = shape.Uniform(lower=0.05, upper=0.1)
    @paramB(name='UDDT operating cost', element=B2, kind='cost',
            units='fraction of capital cost', baseline=b, distribution=D)
    def set_UDDT_OPEX(i):
        B2.OPEX_over_CAPEX = i

    # Conveyance
    B3 = unitB.B3
    B4 = unitB.B4
    b = B3.loss_ratio
    D = shape.Uniform(lower=0.02, upper=0.05)
    @paramB(name='Transportation loss', element=B3, kind='coupled', units='fraction',
           baseline=b, distribution=D)
    def set_trans_loss(i):
        B3.loss_ratio = B4.loss_ratio = i

    b = B3.single_truck.distance
    D = shape.Uniform(lower=2, upper=10)
    @paramB(name='Transportation distance', element=B3, kind='coupled', units='km',
           baseline=b, distribution=D)
    def set_trans_distance(i):
        B3.single_truck.distance = B4.single_truck.distance = i

    b = br.handcart_fee
    D = shape.Uniform(lower=0.004, upper=0.015)
    @paramB(name='Handcart fee', element=B3, kind='cost', units='USD',
           baseline=b, distribution=D)
    def set_handcart_fee(i):
        br.handcart_fee = i

    b = br.truck_fee
    D = shape.Uniform(lower=4.82, upper=8.5)
    @paramB(name='Truck fee', element=B3, kind='cost', units='USD',
           baseline=b, distribution=D)
    def set_truck_fee(i):
        br.truck_fee = i

    # Struvite precipitation
    batch_setting_unit_params(struvite_data, modelB, unitB.B5)

    # Ion exchange
    batch_setting_unit_params(ix_data, modelB, unitB.B6)

    return modelB


def create_modelC(country_specific=False, **model_kwargs):
    sysC = create_system('C')
    unitC = sysC.flowsheet.unit

    # Shared parameters
    modelC = Model(sysC, **model_kwargs)
    add_metrics(modelC)
    unit_dctC = {
        'Excretion': unitC.C1,
        'Toilet': unitC.C2,
        'Conveyance': unitC.C3,
        'ControlBox': unitC.C4, # industrial control panel
        'Housing': unitC.C5,
        'Screw': unitC.C6,
        'LiquidBed': unitC.C7,
        'Carbonizer': unitC.C8,
        'PCD': unitC.C9, # pollution control device
        'OilHX': unitC.C10,
        'HHX': unitC.C11,
        'HHXdryer': unitC.C12,
        }
    add_shared_parameters(modelC, unit_dctC, country_specific)

    # Pit latrine and conveyance
    add_pit_latrine_parameters(modelC, unit_dctC)

    return modelC


# The `country_specific` kwarg is just a placeholder to be consistent
# with other systems, isn't actually being used
def create_modelD(country_specific=False, **model_kwargs):
    sysD = create_system('D')
    unitD = sysD.flowsheet.unit
    modelD = Model(sysD, **model_kwargs)
    add_metrics(modelD)
    add_pit_latrine_parameters(modelD, unit_dct={'Toilet': unitD.D2})
    paramD = modelD.parameter

    batch_setting_unit_params(anaerobic_lagoon_data, modelD, unitD.D4)

    D5 = unitD.D5
    batch_setting_unit_params(drying_bed_data, modelD, D5, exclude=('sol_frac', 'bed_H',))

    b = D5.sol_frac
    if D5.design_type == 'unplanted':
        D = shape.Uniform(lower=0.3, upper=0.4)
    elif D5.design_type == 'planted':
        D = shape.Uniform(lower=0.4, upper=0.7)
    paramD(setter=DictAttrSetter(D5, '_sol_frac', getattr(D5, 'design_type')),
           name='sol_frac', element=D5, kind='coupled', units='fraction',
           baseline=b, distribution=D)

    b = D5.bed_H['covered']
    D = shape.Uniform(lower=0.45, upper=0.75)
    paramD(setter=DictAttrSetter(D5, 'bed_H', ('covered', 'uncovered')),
           name='non_storage_bed_H', element=D5, kind='coupled', units='m',
           baseline=b, distribution=D)

    b = D5.bed_H['storage']
    D = shape.Uniform(lower=1.2, upper=1.8)
    paramD(DictAttrSetter(D5, 'bed_H', 'storage'),
           name='storage_bed_H', element=D5, kind='coupled', units='m',
           baseline=b, distribution=D)

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