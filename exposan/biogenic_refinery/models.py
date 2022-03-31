#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import os, numpy as np, pandas as pd
from chaospy import distributions as shape
from biosteam.evaluation import Metric
from qsdsan import Model, PowerUtility, ImpactItem
from qsdsan.utils import (
    load_data, data_path, AttrSetter, DictAttrSetter, time_printer, dct_from_str
    )
from exposan import biogenic_refinery as br
from exposan.biogenic_refinery import data_path as br_data_path, results_path

__all__ = (
    'create_model', 'result_dct',
    'run_uncertainty', 'save_uncertainty_results',
    )

# %%

# =============================================================================
# Functions for batch-making metrics and -setting parameters
# =============================================================================

systems = br.systems
currency = systems.currency
sys_dct = systems.sys_dct
unit_dct = systems.unit_dct
price_dct = systems.price_dct
GWP_dct = systems.GWP_dct
GWP = systems.GWP
get_summarizing_functions = systems.get_summarizing_functions
func = get_summarizing_functions()

def add_metrics(system):
    sys_ID = system.ID
    tea = sys_dct['TEA'][sys_ID]
    lca = sys_dct['LCA'][sys_ID]
    ppl = sys_dct['ppl'][sys_ID]
    front_end = unit_dct['front_end'][sys_ID]
    transport = unit_dct['transport'][sys_ID]
    pretreatment= unit_dct['pretreatment'][sys_ID]
    liq_treatment = unit_dct['liq_treatment'][sys_ID]
    controls_housing = unit_dct['controls_housing'][sys_ID]
    dryer_hhx = unit_dct['dryer_hhx'][sys_ID]
    ohx = unit_dct['ohx'][sys_ID]
    carbonizer_base = unit_dct['carbonizer_base'][sys_ID]
    pollution_control = unit_dct['pollution_control'][sys_ID]
    unit = f'{currency}/cap/yr'
    cat = 'TEA results'
    metrics = [
        Metric('Net cost', lambda: func['get_annual_cost'](tea, ppl, system), unit, cat),
        Metric('Annual CAPEX', lambda: func['get_annual_CAPEX'](tea, ppl, system), unit, cat),
        Metric('Energy', lambda: func['get_annual_energy'](tea, ppl), unit, cat),
        Metric('Annual OPEX', lambda: func['get_annual_OPEX'](tea, ppl), unit, cat),
        Metric('Labor', lambda: func['get_annual_labor'](tea, ppl), unit, cat),
        Metric('Sales', lambda: func['get_annual_sales'](tea, ppl), unit, cat),
        ]
    unit = f'{GWP.unit}/cap/yr'
    cat = 'LCA results'
    metrics.extend([
        Metric('Net emission', lambda: func['get_annual_GWP'](lca, ppl), unit, cat),
        Metric('Construction', lambda: func['get_constr_GWP'](lca, ppl), unit, cat),
        Metric('Transportation', lambda: func['get_trans_GWP'](lca, ppl), unit, cat),
        Metric('Fugitive gas', lambda: func['get_CH4_N2O_GWP'](system, lca, ppl), unit, cat),
        Metric('Stream items', lambda: func['get_stream_items_emission_GWP'](system, lca, ppl), unit, cat),
        Metric('Offset', lambda: func['get_offset_GWP'](lca, ppl), unit, cat),
        Metric('Other', lambda: func['get_other_GWP'](lca, ppl), unit, cat),
        ])
    unit = 'kg carbon/yr'
    cat = 'Carbon recovery'
    metrics.extend([
        Metric('C urine', lambda: func['C_urine'](system), unit, cat),
        Metric('C feces', lambda: func['C_feces'](system), unit, cat),
        Metric('C toilet sol', lambda: func['C_toilet_sol'](system), unit, cat),
        Metric('C toilet gas', lambda: func['C_toilet_gas'](system), unit, cat),
        Metric('C transport sol', lambda: func['C_trans_sol'](system), unit, cat),
        Metric('C transport sol loss', lambda: func['C_trans_sol_loss'](system), unit, cat),
        Metric('C pretreat sol', lambda: func['C_pretreat_sol'](system), unit, cat),
        Metric('C dryer sol', lambda: func['C_dryer_sol'](system), unit, cat),
        Metric('C dryer gas', lambda: func['C_dryer_gas'](system), unit, cat),
        Metric('C pyrolysis biochar', lambda: func['C_pyrolysis_biochar'](system), unit, cat),
        Metric('C pyrolysis gas', lambda: func['C_pyrolysis_gas'](system), unit, cat),
        Metric('C toilet liq', lambda: func['C_toilet_liq'](system), unit, cat),
        Metric('C trans liq', lambda: func['C_trans_liq'](system), unit, cat),
        Metric('C trans liq loss', lambda: func['C_trans_liq_loss'](system), unit, cat),
        Metric('C pretreat liq', lambda: func['C_pretreat_liq'](system), unit, cat),
        Metric('C treat liq', lambda: func['C_treat_liq'](system), unit, cat),
        Metric('C bed liq', lambda: func['C_bed_liq'](system), unit, cat),
        Metric('C bed gas', lambda: func['C_bed_gas'](system), unit, cat),
        ])
    unit = 'kg nitrogen/yr'
    cat = 'Nitrogen recovery'
    metrics.extend([
        Metric('N urine', lambda: func['N_urine'](system), unit, cat),
        Metric('N feces', lambda: func['N_feces'](system), unit, cat),
        Metric('N toilet sol', lambda: func['N_toilet_sol'](system), unit, cat),
        Metric('N toilet gas', lambda: func['N_toilet_gas'](system), unit, cat),
        Metric('N transport sol', lambda: func['N_trans_sol'](system), unit, cat),
        Metric('N transport sol loss', lambda: func['N_trans_sol_loss'](system), unit, cat),
        Metric('N pretreat sol', lambda: func['N_pretreat_sol'](system), unit, cat),
        Metric('N dryer sol', lambda: func['N_dryer_sol'](system), unit, cat),
        Metric('N dryer gas', lambda: func['N_dryer_gas'](system), unit, cat),
        Metric('N pyrolysis biochar', lambda: func['N_pyrolysis_biochar'](system), unit, cat),
        Metric('N pyrolysis gas', lambda: func['N_pyrolysis_gas'](system), unit, cat),
        Metric('N toilet liq', lambda: func['N_toilet_liq'](system), unit, cat),
        Metric('N trans liq', lambda: func['N_trans_liq'](system), unit, cat),
        Metric('N trans liq loss', lambda: func['N_trans_liq_loss'](system), unit, cat),
        Metric('N pretreat liq', lambda: func['N_pretreat_liq'](system), unit, cat),
        Metric('N treat liq', lambda: func['N_treat_liq'](system), unit, cat),
        Metric('N bed liq', lambda: func['N_bed_liq'](system), unit, cat),
        Metric('N bed gas', lambda: func['N_bed_gas'](system), unit, cat),
        Metric('N struvite', lambda: func['N_struvite'](system), unit, cat),
        Metric('N NH3', lambda: func['N_NH3'](system), unit, cat),
        ])
    unit = 'kg phosphorus/yr'
    cat = 'phosphorus recovery'
    metrics.extend([
        Metric('P urine', lambda: func['P_urine'](system), unit, cat),
        Metric('P feces', lambda: func['P_feces'](system), unit, cat),
        Metric('P toilet sol', lambda: func['P_toilet_sol'](system), unit, cat),
        Metric('P toilet gas', lambda: func['P_toilet_gas'](system), unit, cat),
        Metric('P transport sol', lambda: func['P_trans_sol'](system), unit, cat),
        Metric('P transport sol loss', lambda: func['P_trans_sol_loss'](system), unit, cat),
        Metric('P pretreat sol', lambda: func['P_pretreat_sol'](system), unit, cat),
        Metric('P dryer sol', lambda: func['P_dryer_sol'](system), unit, cat),
        Metric('P dryer gas', lambda: func['P_dryer_gas'](system), unit, cat),
        Metric('P pyrolysis biochar', lambda: func['P_pyrolysis_biochar'](system), unit, cat),
        Metric('P pyrolysis gas', lambda: func['P_pyrolysis_gas'](system), unit, cat),
        Metric('P toilet liq', lambda: func['P_toilet_liq'](system), unit, cat),
        Metric('P trans liq', lambda: func['P_trans_liq'](system), unit, cat),
        Metric('P trans liq loss', lambda: func['P_trans_liq_loss'](system), unit, cat),
        Metric('P pretreat liq', lambda: func['P_pretreat_liq'](system), unit, cat),
        Metric('P treat liq', lambda: func['P_treat_liq'](system), unit, cat),
        Metric('P bed liq', lambda: func['P_bed_liq'](system), unit, cat),
        Metric('P bed gas', lambda: func['P_bed_gas'](system), unit, cat),
        Metric('P struvite', lambda: func['P_struvite'](system), unit, cat),
        ])
    unit = 'kg potassium/yr'
    cat = 'potassium recovery'
    metrics.extend([
        Metric('K urine', lambda: func['K_urine'](system), unit, cat),
        Metric('K feces', lambda: func['K_feces'](system), unit, cat),
        Metric('K toilet sol', lambda: func['K_toilet_sol'](system), unit, cat),
        Metric('K toilet gas', lambda: func['K_toilet_gas'](system), unit, cat),
        Metric('K transport sol', lambda: func['K_trans_sol'](system), unit, cat),
        Metric('K transport sol loss', lambda: func['K_trans_sol_loss'](system), unit, cat),
        Metric('K pretreat sol', lambda: func['K_pretreat_sol'](system), unit, cat),
        Metric('K dryer sol', lambda: func['K_dryer_sol'](system), unit, cat),
        Metric('K dryer gas', lambda: func['K_dryer_gas'](system), unit, cat),
        Metric('K pyrolysis biochar', lambda: func['K_pyrolysis_biochar'](system), unit, cat),
        Metric('K pyrolysis gas', lambda: func['K_pyrolysis_gas'](system), unit, cat),
        Metric('K toilet liq', lambda: func['K_toilet_liq'](system), unit, cat),
        Metric('K trans liq', lambda: func['K_trans_liq'](system), unit, cat),
        Metric('K trans liq loss', lambda: func['K_trans_liq_loss'](system), unit, cat),
        Metric('K pretreat liq', lambda: func['K_pretreat_liq'](system), unit, cat),
        Metric('K treat liq', lambda: func['K_treat_liq'](system), unit, cat),
        Metric('K bed liq', lambda: func['K_bed_liq'](system), unit, cat),
        Metric('K bed gas', lambda: func['K_bed_gas'](system), unit, cat),
        ])

    unit = f'{currency}/cap/d'
    cat = 'cost breakdown'
    metrics.extend([
        Metric('operator_cost_opex', lambda: func['operator_cost_opex'](tea,ppl), unit, cat),

        Metric('cost_capex_front_end', lambda: func['front_end_cost_capex'](front_end,tea,ppl), unit, cat),
        Metric('front_end_cost_opex', lambda: func['front_end_cost_opex'](front_end,ppl), unit, cat),
        Metric('front_end_cost_electricity', lambda: func['front_end_cost_electricity'](front_end,ppl), unit, cat),

        Metric('transport_cost_opex', lambda: func['transport_cost_opex'](transport,ppl), unit, cat),

        Metric('pretreatment_cost_capex', lambda: func['pretreatment_cost_capex'](pretreatment,tea,ppl), unit, cat),
        Metric('pretreatment_cost_opex', lambda: func['pretreatment_cost_opex'](pretreatment,ppl), unit, cat),
        Metric('pretreatment_cost_electricity', lambda: func['pretreatment_cost_electricity'](pretreatment,ppl), unit, cat),

        Metric('liq_treatment_cost_capex', lambda: func['liq_treatment_cost_capex'](liq_treatment,tea,ppl), unit, cat),
        Metric('liq_treatment_cost_opex', lambda: func['liq_treatment_cost_opex'](liq_treatment,ppl), unit, cat),
        Metric('liq_treatment_cost_electricity', lambda: func['liq_treatment_cost_electricity'](liq_treatment,ppl), unit, cat),

        Metric('controls_housing_cost_capex', lambda: func['controls_housing_cost_capex'](controls_housing,tea,ppl), unit, cat),
        Metric('controls_housing_cost_opex', lambda: func['controls_housing_cost_opex'](controls_housing,ppl), unit, cat),
        Metric('controls_housing_cost_electricity', lambda: func['controls_housing_cost_electricity'](controls_housing,ppl), unit, cat),

        Metric('dryer_hhx_cost_capex', lambda: func['dryer_hhx_cost_capex'](dryer_hhx,tea,ppl), unit, cat),
        Metric('dryer_hhx_cost_opex', lambda: func['dryer_hhx_cost_opex'](dryer_hhx,ppl), unit, cat),
        Metric('dryer_hhx_cost_electricity', lambda: func['dryer_hhx_cost_electricity'](dryer_hhx,ppl), unit, cat),

        Metric('ohx_cost_capex', lambda: func['ohx_cost_capex'](ohx,tea,ppl), unit, cat),
        Metric('ohx_cost_opex', lambda: func['ohx_cost_opex'](ohx,ppl), unit, cat),
        Metric('ohx_cost_electricity', lambda: func['ohx_cost_electricity'](ohx,ppl), unit, cat),

        Metric('carbonizer_base_cost_capex', lambda: func['carbonizer_base_cost_capex'](carbonizer_base,tea,ppl), unit, cat),
        Metric('carbonizer_base_cost_opex', lambda: func['carbonizer_base_cost_opex'](carbonizer_base,ppl), unit, cat),
        Metric('carbonizer_base_cost_electricity', lambda: func['carbonizer_base_cost_electricity'](carbonizer_base,ppl), unit, cat),

        Metric('pollution_control_cost_capex', lambda: func['pollution_control_cost_capex'](pollution_control,tea,ppl), unit, cat),
        Metric('pollution_control_cost_opex', lambda: func['pollution_control_cost_opex'](pollution_control,ppl), unit, cat),
        Metric('pollution_control_cost_electricity', lambda: func['pollution_control_cost_electricity'](pollution_control,ppl), unit, cat),
        ])

    unit = f'{GWP.unit}/cap/yr'
    cat = 'GHG breakdown'
    metrics.extend([
        Metric('GHG_ghg_front_end', lambda: func['front_end_ghg_capex'](front_end,lca,ppl), unit, cat),
        Metric('front_end_ghg_opex', lambda: func['front_end_ghg_opex'](front_end,lca,ppl), unit, cat),
        Metric('front_end_ghg_electricity', lambda: func['front_end_ghg_electricity'](front_end,ppl), unit, cat),
        Metric('front_end_ghg_direct', lambda: func['front_end_ghg_direct'](front_end,lca,ppl), unit, cat),

        Metric('transport_ghg_opex', lambda: func['transport_ghg_opex'](transport,lca,ppl), unit, cat),

        Metric('pretreatment_ghg_capex', lambda: func['pretreatment_ghg_capex'](pretreatment,lca,ppl), unit, cat),
        Metric('pretreatment_ghg_opex', lambda: func['pretreatment_ghg_opex'](pretreatment,lca,ppl), unit, cat),
        Metric('pretreatment_ghg_electricity', lambda: func['pretreatment_ghg_electricity'](pretreatment,ppl), unit, cat),
        Metric('pretreatment_ghg_direct', lambda: func['pretreatment_ghg_direct'](pretreatment,lca,ppl), unit, cat),

        Metric('liq_treatment_ghg_capex', lambda: func['liq_treatment_ghg_capex'](liq_treatment,lca,ppl), unit, cat),
        Metric('liq_treatment_ghg_opex', lambda: func['liq_treatment_ghg_opex'](liq_treatment,lca,ppl), unit, cat),
        Metric('liq_treatment_ghg_electricity', lambda: func['liq_treatment_ghg_electricity'](liq_treatment,ppl), unit, cat),
        Metric('liq_treatment_ghg_direct', lambda: func['liq_treatment_ghg_direct'](liq_treatment,lca,ppl), unit, cat),

        Metric('controls_housing_ghg_capex', lambda: func['controls_housing_ghg_capex'](controls_housing,lca,ppl), unit, cat),
        Metric('controls_housing_ghg_opex', lambda: func['controls_housing_ghg_opex'](controls_housing,lca,ppl), unit, cat),
        Metric('controls_housing_ghg_electricity', lambda: func['controls_housing_ghg_electricity'](controls_housing,ppl), unit, cat),
        Metric('controls_housing_ghg_direct', lambda: func['controls_housing_ghg_direct'](controls_housing,lca,ppl), unit, cat),

        Metric('dryer_hhx_ghg_capex', lambda: func['dryer_hhx_ghg_capex'](dryer_hhx,lca,ppl), unit, cat),
        Metric('dryer_hhx_ghg_opex', lambda: func['dryer_hhx_ghg_opex'](dryer_hhx,lca,ppl), unit, cat),
        Metric('dryer_hhx_ghg_electricity', lambda: func['dryer_hhx_ghg_electricity'](dryer_hhx,ppl), unit, cat),
        Metric('dryer_hhx_ghg_direct', lambda: func['dryer_hhx_ghg_direct'](dryer_hhx,lca,ppl), unit, cat),

        Metric('ohx_ghg_capex', lambda: func['ohx_ghg_capex'](ohx,lca,ppl), unit, cat),
        Metric('ohx_ghg_opex', lambda: func['ohx_ghg_opex'](ohx,lca,ppl), unit, cat),
        Metric('ohx_ghg_electricity', lambda: func['ohx_ghg_electricity'](ohx,ppl), unit, cat),
        Metric('ohx_ghg_direct', lambda: func['ohx_ghg_direct'](ohx,lca,ppl), unit, cat),

        Metric('carbonizer_base_ghg_capex', lambda: func['carbonizer_base_ghg_capex'](carbonizer_base,lca,ppl), unit, cat),
        Metric('carbonizer_base_ghg_opex', lambda: func['carbonizer_base_ghg_opex'](carbonizer_base,lca,ppl), unit, cat),
        Metric('carbonizer_base_ghg_electricity', lambda: func['carbonizer_base_ghg_electricity'](carbonizer_base,ppl), unit, cat),
        Metric('carbonizer_base_ghg_direct', lambda: func['carbonizer_base_ghg_direct'](carbonizer_base,lca,ppl), unit, cat),

        Metric('pollution_control_ghg_capex', lambda: func['pollution_control_ghg_capex'](pollution_control,lca,ppl), unit, cat),
        Metric('pollution_control_ghg_opex', lambda: func['pollution_control_ghg_opex'](pollution_control,lca,ppl), unit, cat),
        Metric('pollution_control_ghg_electricity', lambda: func['pollution_control_ghg_electricity'](pollution_control,ppl), unit, cat),
        Metric('pollution_control_ghg_direct', lambda: func['pollution_control_ghg_direct'](pollution_control,lca,ppl), unit, cat),
        ])

    return metrics


def batch_setting_unit_params(df, model, unit, exclude=()):
    for para in df.index:
        if para in exclude: continue
        b = getattr(unit, para)
        lower = float(df.loc[para]['low'])
        upper = float(df.loc[para]['high'])
        dist = df.loc[para]['distribution']
        if dist == 'uniform':
            D = shape.Uniform(lower=lower, upper=upper)
        elif dist == 'triangular':
            D = shape.Triangle(lower=lower, midpoint=b, upper=upper)
        elif dist == 'constant': continue
        else:
            raise ValueError(f'Distribution {dist} not recognized for unit {unit}.')
        model.parameter(setter=AttrSetter(unit, para),
                        name=para, element=unit, kind='coupled', units=df.loc[para]['unit'],
                        baseline=b, distribution=D)


# %%

# =============================================================================
# Pre-load data sheets as they will be used in multiple systems
# =============================================================================

join_path = lambda prefix, file_name: os.path.join(prefix, file_name)
su_data_path = join_path(data_path, 'sanunit_data')
excretion_data = load_data(join_path(su_data_path, '_excretion.tsv'))
control_data = load_data(join_path(su_data_path, '_control_box_op.tsv'))
housing_data = load_data(join_path(su_data_path, '_housing_biogenic_refinery.tsv'))
# Screw press for sysA/sysC, grinder for sysB
screw_data = load_data(join_path(su_data_path, '_screw_press.tsv'))
grinder_data = load_data(join_path(su_data_path, '_grinder.tsv'))
liquid_bed_data = load_data(join_path(su_data_path, '_liquid_treatment_bed.tsv'))
carbonizer_data = load_data(join_path(su_data_path, '_carbonizer_base.tsv'))
pcd_data = load_data(join_path(su_data_path, '_pollution_control_device.tsv'))
oil_hx_data = load_data(join_path(su_data_path, '_oil_heat_exchanger.tsv'))
hhx_data = load_data(join_path(su_data_path, '_hydronic_heat_exchanger.tsv'))
hhx_dryer_data = load_data(join_path(su_data_path, '_hhx_dryer.tsv'))


# %%

# =============================================================================
# Shared by systems A-C
# =============================================================================

def add_shared_parameters(sys, model, unit_dct, country_specific=False):
    param = model.parameter
    streams = sys_dct['stream_dct'][sys.ID]

    # Add these parameters if not running country-specific analysis,
    # in which they would be updated separately
    if not country_specific:
        b = systems.operator_daily_wage
        D = shape.Triangle(lower=(14.55), midpoint=b, upper=(43.68))
        @param(name='Operator daily wages', element='TEA', kind='cost', units='USD/d',
              baseline=b, distribution=D)
        def set_operator_daily_wage(i):
            sys._TEA.annual_labor = i*3*365

        # Discount factor for the excreta-derived fertilizers
        get_price_factor = systems.get_price_factor
        b = get_price_factor()
        D = shape.Uniform(lower=0.1, upper=0.4)
        @param(name='Price factor', element='TEA', kind='isolated', units='-',
                baseline=b, distribution=D)
        def set_price_factor(i):
            systems.price_factor = i

        D = shape.Uniform(lower=1.164, upper=2.296)
        @param(name='N fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
                baseline=1.507, distribution=D)
        def set_N_price(i):
            price_dct['N'] = streams['liq_N'] = streams['sol_N'] = i * get_price_factor()

        D = shape.Uniform(lower=2.619, upper=6.692)
        @param(name='P fertilizer price', element='TEA', kind='isolated', units='USD/kg P',
                baseline=3.983, distribution=D)
        def set_P_price(i):
            price_dct['P'] = streams['liq_P'] = streams['sol_P'] = i * get_price_factor()

        D = shape.Uniform(lower=1.214, upper=1.474)
        @param(name='K fertilizer price', element='TEA', kind='isolated', units='USD/kg K',
                baseline=1.333, distribution=D)
        def set_K_price(i):
            price_dct['K'] = streams['liq_K'] = streams['sol_K'] = i * get_price_factor()

        # Electricity price
        b = price_dct['Electricity']
        D = shape.Triangle(lower=0.08, midpoint=b, upper=0.14)
        @param(name='Electricity price', element='TEA', kind='isolated',
               units='$/kWh', baseline=b, distribution=D)
        def set_electricity_price(i):
            PowerUtility.price = i

        b = GWP_dct['Electricity']
        D = shape.Triangle(lower=0.6212, midpoint=b, upper=0.7592)
        @param(name='Electricity CF', element='LCA', kind='isolated',
                   units='kg CO2-eq/kWh', baseline=b, distribution=D)
        def set_electricity_CF(i):
            GWP_dct['Electricity'] = ImpactItem.get_item('e_item').CFs['GlobalWarming'] = i

        D = shape.Uniform(lower=(1.507*(14/17)*0.8), upper=(1.507*(14/17)*1.2))
        @param(name='NH3 fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
                baseline=(1.507*(14/17)*0.25), distribution=D)
        def set_con_NH3_price(i):
            price_dct['conc_NH3'] = streams['conc_NH3'] = streams['conc_NH3'] = i * get_price_factor()

        D = shape.Uniform(lower=(3.983*(31/245)*0.8), upper=(3.983*(31/245)*1.2))
        @param(name='Struvite fertilizer price', element='TEA', kind='isolated', units='USD/kg P',
                baseline=(3.983*(31/245)*0.25), distribution=D)
        def set_struvite_price(i):
            price_dct['struvite'] = streams['struvite'] = streams['struvite'] = i * get_price_factor()

        D = shape.Uniform(lower=(0.014*0.8), upper=(0.014*1.2))
        @param(name='Biochar  price', element='TEA', kind='isolated', units='USD/kg biochar',
                baseline=(0.014), distribution=D)
        def set_biochar_price(i):
            price_dct['biochar'] = streams['biochar'] = streams['biochar'] = i

        b = -GWP_dct['N']
        D = shape.Triangle(lower=1.8, midpoint=b, upper=8.9)
        @param(name='N fertilizer CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kg N', baseline=b, distribution=D)
        def set_N_fertilizer_CF(i):
            GWP_dct['N'] = systems.N_item.CFs['GlobalWarming'] = -i

        b = -GWP_dct['P']
        D = shape.Triangle(lower=4.3, midpoint=b, upper=5.4)
        @param(name='P fertilizer CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kg P', baseline=b, distribution=D)
        def set_P_fertilizer_CF(i):
            GWP_dct['P'] = systems.P_item.CFs['GlobalWarming'] = -i

        b = -GWP_dct['K']
        D = shape.Triangle(lower=1.1, midpoint=b, upper=2)
        @param(name='K fertilizer CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kg K', baseline=b, distribution=D)
        def set_K_fertilizer_CF(i):
            GWP_dct['K'] = systems.K_item.CFs['GlobalWarming'] = -i

        b = -GWP_dct['conc_NH3']
        D = shape.Triangle(lower=(1.8*(14/17)), midpoint=b, upper=(8.9*(14/17)))
        @param(name='conc_NH3 fertilizer CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kg conc_NH3', baseline=b, distribution=D)
        def set_conc_NH3_fertilizer_CF(i):
            GWP_dct['conc_NH3'] = systems.conc_NH3_item.CFs['GlobalWarming'] = -i

        b = -GWP_dct['struvite']
        D = shape.Triangle(lower=(4.3*(31/245)), midpoint=b, upper=(5.4*(31/245)))
        @param(name='struvite fertilizer CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kg struvite', baseline=b, distribution=D)
        def set_struvite_fertilizer_CF(i):
            GWP_dct['struvite'] = systems.struvite_item.CFs['GlobalWarming'] = -i

        b = -GWP_dct['biochar']
        D = shape.Triangle(lower=(0.2*0.9*(44/12)*.8), midpoint=b, upper=(0.2*0.9*(44/12)*1.2))
        @param(name='biochar CF', element='LCA', kind='isolated',
                units='kg CO2-eq/kg biochar', baseline=b, distribution=D)
        def set_biochar_CF(i):
            GWP_dct['biochar'] = systems.biochar_item.CFs['GlobalWarming'] = -i

    # ########## Specific units ##########
    # Diet and excretion
    excretion_unit = unit_dct['Excretion']
    exclude = ('e_cal','p_anim','p_veg') if country_specific else ()
    batch_setting_unit_params(excretion_data, model, excretion_unit, exclude)

    # Household size
    b = systems.household_size
    D = shape.Normal(mu=b, sigma=1.8)
    @param(name='Household size', element=excretion_unit, kind='coupled', units='cap/household',
           baseline=b, distribution=D)
    def set_household_size(i):
        systems.household_size = max(1, i)

    # Toilet
    toilet_unit = unit_dct['Toilet']
    b = systems.household_per_toilet
    D = shape.Uniform(lower=3, upper=5)
    @param(name='Toilet density', element=toilet_unit, kind='coupled', units='household/toilet',
            baseline=b, distribution=D)
    def set_toilet_density(i):
        systems.household_per_toilet = i

    # Control box (industrial control panel)
    control_unit = unit_dct['ControlBox']
    batch_setting_unit_params(control_data, model, control_unit)

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
    batch_setting_unit_params(carbonizer_data, model, carbonizer_unit)

    # Pollution control device
    pcd_unit = unit_dct['PCD']
    batch_setting_unit_params(pcd_data, model, pcd_unit)

    # Oil heat exchanger
    oil_hx_unit = unit_dct['OilHX']
    batch_setting_unit_params(oil_hx_data, model, oil_hx_unit)

    # Hydronic heat exchanger
    hhx_unit = unit_dct['HHX']
    batch_setting_unit_params(hhx_data, model, hhx_unit)

    # Dryer from HHX
    hhx_dryer_unit = unit_dct['HHXdryer']
    batch_setting_unit_params(hhx_dryer_data, model, hhx_dryer_unit)

    ##### Universal degradation parameters #####
    # Max methane emission
    toilet_unit = sys.path[1] # the first unit that involves degradation
    b = systems.max_CH4_emission
    D = shape.Triangle(lower=0.175, midpoint=b, upper=0.325)
    @param(name='Max CH4 emission', element=toilet_unit, kind='coupled', units='g CH4/g COD',
           baseline=b, distribution=D)
    def set_max_CH4_emission(i):
        systems.max_CH4_emission = i

    # Time to full degradation
    b = systems.tau_deg
    D = shape.Uniform(lower=1, upper=3)
    @param(name='Full degradation time', element=toilet_unit, kind='coupled', units='yr',
           baseline=b, distribution=D)
    def set_tau_deg(i):
        systems.tau_deg = i

    # Reduction at full degradation
    b = systems.log_deg
    D = shape.Uniform(lower=2, upper=4)
    @param(name='Log degradation', element=toilet_unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_log_deg(i):
        systems.log = i

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

    ######## General TEA settings ########
    # # Keeping discount rate constant
    # b = systems.discount_rate
    # D = shape.Uniform(lower=0.03, upper=0.06)
    # @param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
    #         baseline=b, distribution=D)
    # def set_discount_rate(i):
    #     systems.discount_rate = i

    b = price_dct['Polymer']
    D = shape.Uniform(lower=(b*0.95), upper=(b*1.05))
    @param(name='Polymer price', element='TEA', kind='isolated', units='USD/kg polymer',
            baseline=(b), distribution=D)
    def set_polymber_price(i):
        price_dct['Polymer'] = streams['polymer'].price = i

    b = price_dct['Resin']
    D = shape.Uniform(lower=(b*0.95), upper=(b*1.05))
    @param(name='Resin price', element='TEA', kind='isolated', units='USD/kg resin',
            baseline=(b), distribution=D)
    def set_resin_price(i):
        price_dct['Resin'] = streams['resin'].price = i

    b = price_dct['FilterBag']
    D = shape.Uniform(lower=(b*0.95), upper=(b*1.05))
    @param(name='Filter bag price', element='TEA', kind='isolated', units='USD/kg filter bag',
            baseline=(b), distribution=D)
    def set_filterbag_price(i):
        price_dct['FilterBag'] = streams['filter_bag'].price = i

    ######## General LCA settings ########
    b = GWP_dct['CH4']
    D = shape.Uniform(lower=28, upper=34)
    @param(name='CH4 CF', element='LCA', kind='isolated', units='kg CO2-eq/kg CH4',
           baseline=b, distribution=D)
    def set_CH4_CF(i):
        GWP_dct['CH4'] = systems.CH4_item.CFs['GlobalWarming'] = i

    b = GWP_dct['N2O']
    D = shape.Uniform(lower=265, upper=298)
    @param(name='N2O CF', element='LCA', kind='isolated', units='kg CO2-eq/kg N2O',
           baseline=b, distribution=D)
    def set_N2O_CF(i):
        GWP_dct['N2O'] = systems.N2O_item.CFs['GlobalWarming'] = i

    b = GWP_dct['Polymer']
    D = shape.Triangle(lower=b*0.95, midpoint=b, upper=b*1.05)
    @param(name='Polymer CF', element='LCA', kind='isolated',
            units='kg CO2-eq/kg N', baseline=b, distribution=D)
    def set_polymer_CF(i):
        GWP_dct['Polymer'] = systems.polymer_item.CFs['GlobalWarming'] = i

    b = GWP_dct['Resin']
    D = shape.Triangle(lower=b*0.95, midpoint=b, upper=b*1.05)
    @param(name='Resin CF', element='LCA', kind='isolated',
            units='kg CO2-eq/kg N', baseline=b, distribution=D)
    def set_resin_CF(i):
        GWP_dct['Resin'] = systems.resin_item.CFs['GlobalWarming'] = i

    b = GWP_dct['FilterBag']
    D = shape.Triangle(lower=b*0.95, midpoint=b, upper=b*1.05)
    @param(name='Filter bag CF', element='LCA', kind='isolated',
            units='kg CO2-eq/kg N', baseline=b, distribution=D)
    def set_filter_bag_CF(i):
        GWP_dct['FilterBag'] = systems.filter_bag_item.CFs['GlobalWarming'] = i

    item_path = join_path(br_data_path, 'impact_items.xlsx')
    data = load_data(item_path, sheet='GWP')
    for p in data.index:
        item = ImpactItem.get_item(p)
        b = item.CFs['GlobalWarming']
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
        model.parameter(name=p,
                        setter=DictAttrSetter(item, 'CFs', 'GlobalWarming'),
                        element='LCA', kind='isolated',
                        units=f'kg CO2-eq/{item.functional_unit}',
                        baseline=b, distribution=D)

    return model


# %%

# =============================================================================
# Shared by systems A, C, and D
# =============================================================================

path = join_path(su_data_path, '_toilet.tsv')
toilet_data = load_data(path)
path = join_path(su_data_path, '_pit_latrine.tsv')
pit_latrine_data = load_data(path)

MCF_lower_dct = dct_from_str(pit_latrine_data.loc['MCF_decay']['low'])
MCF_upper_dct = dct_from_str(pit_latrine_data.loc['MCF_decay']['high'])
N2O_EF_lower_dct = dct_from_str(pit_latrine_data.loc['N2O_EF_decay']['low'])
N2O_EF_upper_dct = dct_from_str(pit_latrine_data.loc['N2O_EF_decay']['high'])

def add_pit_latrine_parameters(sys, model, unit_dct):
    pit_unit = unit_dct['Toilet']
    param = model.parameter
    ######## Related to the toilet ########
    all_pit_data = pd.concat((toilet_data, pit_latrine_data))
    batch_setting_unit_params(all_pit_data, model, pit_unit, exclude=('MCF_decay', 'N2O_EF_decay'))

    kind = pit_unit._return_MCF_EF()
    if sys.ID != 'sysD': # Parameters not applicable to sysD
        ######## Related to conveyance ########
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

        b = systems.emptying_fee
        D = shape.Uniform(lower=0, upper=0.3)
        @param(name='Emptying fee', element=convey_unit, kind='coupled', units='USD',
                baseline=b, distribution=D)
        def set_emptying_fee(i):
            systems.emptying_fee = i

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

    return model


# %%

# =============================================================================
# Functions to create models
# =============================================================================

# System A
def create_modelA(country_specific=False):
    sysA = systems.sysA
    modelA = Model(sysA, add_metrics(sysA))

    # Shared parameters
    unit_dctA = {
        'Excretion': systems.A1,
        'Toilet': systems.A2,
        'Conveyance': systems.A3,
        'ControlBox': systems.A4, # industrial control panel
        'Housing': systems.A5,
        'Screw': systems.A6,
        'LiquidBed': systems.A7,
        'Carbonizer': systems.A8,
        'PCD': systems.A9, # pollution control device
        'OilHX': systems.A10,
        'HHX': systems.A11,
        'HHXdryer': systems.A12,
        }
    modelA = add_shared_parameters(sysA, modelA, unit_dctA, country_specific)

    # Pit latrine and conveyance
    modelA = add_pit_latrine_parameters(sysA, modelA, unit_dctA)

    return modelA


# System B
def create_modelB(country_specific=False):
    sysB = systems.sysB
    modelB = Model(sysB, add_metrics(sysB))
    paramB = modelB.parameter

    # Shared parameters
    unit_dctB = {
        'Excretion': systems.B1,
        'Toilet': systems.B2,
        'Conveyance': (systems.B3, systems.B4),
        'LiquidBed': systems.B7,
        'ControlBox': systems.B8, # industrial control panel
        'Housing': systems.B9,
        'Grinder': systems.B10,
        'Carbonizer': systems.B11,
        'PCD': systems.B12, # pollution control device
        'OilHX': systems.B13,
        'HHX': systems.B14,
        'HHXdryer': systems.B15,
        }
    modelB = add_shared_parameters(sysB, modelB, unit_dctB, country_specific)

    # UDDT
    B2 = systems.B2
    uddt_data = load_data(join_path(su_data_path, '_uddt.tsv'))
    all_uddt_data = pd.concat((toilet_data, uddt_data))
    batch_setting_unit_params(all_uddt_data, modelB, B2)

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
    B3 = systems.B3
    B4 = systems.B4
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

    b = systems.handcart_fee
    D = shape.Uniform(lower=0.004, upper=0.015)
    @paramB(name='Handcart fee', element=B3, kind='cost', units='USD',
           baseline=b, distribution=D)
    def set_handcart_fee(i):
        systems.handcart_fee = i

    b = systems.truck_fee
    D = shape.Uniform(lower=4.82, upper=8.5)
    @paramB(name='Truck fee', element=B3, kind='cost', units='USD',
           baseline=b, distribution=D)
    def set_truck_fee(i):
        systems.truck_fee = i

    # Struvite precipitation
    struvite_data = load_data(join_path(su_data_path, '_struvite_precipitation.tsv'))
    batch_setting_unit_params(struvite_data, modelB, systems.B5)

    # Ion exchange NH3
    ix_data = load_data(join_path(su_data_path, '_ion_exchange_NH3.tsv'))
    batch_setting_unit_params(ix_data, modelB, systems.B6)

    return modelB


# System C
def create_modelC(country_specific=False):
    sysC = systems.sysC
    modelC = Model(sysC, add_metrics(sysC))

    # Shared parameters
    unit_dctC = {
        'Excretion': systems.C1,
        'Toilet': systems.C2,
        'Conveyance': systems.C3,
        'ControlBox': systems.C4, # industrial control panel
        'Housing': systems.C5,
        'Screw': systems.C6,
        'LiquidBed': systems.C7,
        'Carbonizer': systems.C8,
        'PCD': systems.C9, # pollution control device
        'OilHX': systems.C10,
        'HHX': systems.C11,
        'HHXdryer': systems.C12,
        }
    modelC = add_shared_parameters(sysC, modelC, unit_dctC, country_specific)

    # Pit latrine and conveyance
    modelC = add_pit_latrine_parameters(sysC, modelC, unit_dctC)

    return modelC


# System D, the `country_specific` kwarg is just a placeholder to be consistent
# with other systems, isn't actually being used
def create_modelD(country_specific=False):
    sysD = systems.sysD
    modelD = Model(sysD, add_metrics(sysD))
    modelD = add_pit_latrine_parameters(sysD, modelD, unit_dct={'Toilet': systems.D2})
    anaerobic_lagoon_data = load_data(join_path(su_data_path, '_anaerobic_lagoon.tsv'))
    batch_setting_unit_params(anaerobic_lagoon_data, modelD, systems.D4)
    return modelD


# Wrapper function so that it'd work for all
def create_model(model_ID='A', country_specific=False):
    model_ID = model_ID.lstrip('model').lstrip('sys') # so that it'll work for "modelA"/"sysA"/"A"
    if model_ID == 'A': model = create_modelA(country_specific)
    elif model_ID == 'B': model = create_modelB(country_specific)
    elif model_ID == 'C': model = create_modelC(country_specific)
    else: model = create_modelD(country_specific)
    return model


# %%

# =============================================================================
# Functions to run simulation and generate plots
# =============================================================================

result_dct = {
        'sysA': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        'sysB': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        'sysC': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        'sysD': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman'))
        }

@time_printer
def run_uncertainty(model, seed=None, N=10000, rule='L',
                    percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                    print_time=False):
    global result_dct
    if seed:
        np.random.seed(seed)

    samples = model.sample(N, rule)
    model.load_samples(samples)
    model.evaluate()

    # Data organization
    dct = result_dct[model.system.ID]
    index_p = len(model.get_parameters())
    dct['parameters'] = model.table.iloc[:, :index_p].copy()
    dct['data'] = model.table.iloc[:, index_p:].copy()
    if percentiles:
        dct['percentiles'] = dct['data'].quantile(q=percentiles)
        dct['percentiles_parameters'] = dct['parameters'].quantile(q=percentiles)

    # Spearman's rank correlation
    spearman_metrics = model.metrics[:13]
    spearman_results = model.spearman(model.get_parameters(), spearman_metrics)
    spearman_results.columns = pd.Index([i.name_with_units for i in spearman_metrics])
    dct['spearman'] = spearman_results
    return dct


def save_uncertainty_results(model, dct={}, path=''):
    sys_ID = model.system.ID
    path = join_path(results_path, f'uncertainty{sys_ID[-1]}.xlsx') if path=='' else path
    dct = dct or result_dct[sys_ID]
    if dct['parameters'] is None:
        raise ValueError('No cached result, run model first.')
    with pd.ExcelWriter(path) as writer:
        dct['parameters'].to_excel(writer, sheet_name='Parameters')
        dct['data'].to_excel(writer, sheet_name='Uncertainty results')
        if 'percentiles' in dct.keys():
            dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
        dct['spearman'].to_excel(writer, sheet_name='Spearman')
        model.table.to_excel(writer, sheet_name='Raw data')