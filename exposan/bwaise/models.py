#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import os, pickle
import numpy as np
import pandas as pd
from chaospy import distributions as shape
from thermosteam.functional import V_to_rho, rho_to_V
from biosteam import PowerUtility
from biosteam.evaluation import Model, Metric
from qsdsan import currency, ImpactItem
from qsdsan.utils import (
    load_data, data_path, dct_from_str,
    AttrSetter, AttrFuncSetter, DictAttrSetter,
    FuncGetter,
    time_printer
    )
from exposan import bwaise as bw

c_path = bw._lca_data.c_path
lca_data_kind = bw.systems.lca_data_kind

__all__ = ('modelA', 'modelB', 'modelC', 'result_dct',
           'run_uncertainty', 'save_uncertainty_results')


# %%

# =============================================================================
# Functions for batch-making metrics and -setting parameters
# =============================================================================

systems = bw.systems
sys_dct = systems.sys_dct
price_dct = systems.price_dct
GWP_dct = systems.GWP_dct
get_summarizing_functions = systems.get_summarizing_functions

def add_LCA_metrics(system, metrics, kind):
    systems.update_lca_data(kind)
    lca = sys_dct['LCA'][system.ID]
    ppl = sys_dct['ppl'][system.ID]
    funcs = [
        lambda ID: lca.total_impacts[ID]/lca.lifetime/ppl,
        lambda ID: lca.total_construction_impacts[ID]/lca.lifetime/ppl,
        lambda ID: lca.total_transportation_impacts[ID]/lca.lifetime/ppl,
        lambda ID: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='direct_emission')[ID] \
            /lca.lifetime/ppl,
        lambda ID: lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='offset')[ID] \
            /lca.lifetime/ppl,
        lambda ID: lca.total_other_impacts[ID]/lca.lifetime/ppl
        ]

    for ind in lca.indicators:
        unit = f'{ind.unit}/cap/yr'
        cat = 'LCA results'
        metrics.extend([
            Metric(f'Net emission {ind.ID}', FuncGetter(funcs[0], (ind.ID,)), unit, cat),
            Metric(f'Construction {ind.ID}', FuncGetter(funcs[1], (ind.ID,)),unit, cat),
            Metric(f'Transportation {ind.ID}', FuncGetter(funcs[2], (ind.ID,)),unit, cat),
            Metric(f'Direct emission {ind.ID}', FuncGetter(funcs[3], (ind.ID,)),unit, cat),
            Metric(f'Offset {ind.ID}', FuncGetter(funcs[4], (ind.ID,)),unit, cat),
            Metric(f'Other {ind.ID}', FuncGetter(funcs[5], (ind.ID,)),unit, cat),
            ])

    return metrics

def add_metrics(system, kind):
    sys_ID = system.ID
    tea = sys_dct['TEA'][sys_ID]
    ppl = sys_dct['ppl'][sys_ID]
    func = get_summarizing_functions(system)

    metrics = []
    for i in ('COD', 'N', 'P', 'K'):
        cat = f'{i} recovery'
        metrics.extend([
            Metric(f'Liquid {i}', FuncGetter(func[f'get_liq_{i}_recovery'], (system, i)), '', cat),
            Metric(f'Solid {i}', FuncGetter(func[f'get_sol_{i}_recovery'], (system, i)), '', cat),
            Metric(f'Gas {i}', FuncGetter(func[f'get_gas_{i}_recovery'], (system, i)), '', cat),
            Metric(f'Total {i}', FuncGetter(func[f'get_tot_{i}_recovery'], (system, i)), '', cat)
            ])

    unit = f'{currency}/cap/yr'
    cat = 'TEA results'
    metrics.extend([
        Metric('Annual net cost', lambda: func['get_annual_net_cost'](tea, ppl), unit, cat),
        Metric('Annual cost', lambda: func['get_annual_cost'](tea, ppl), unit, cat),
        Metric('Annual CAPEX', lambda: func['get_annual_CAPEX'](tea, ppl), unit, cat),
        Metric('Annual OPEX', lambda: func['get_annual_OPEX'](tea, ppl), unit, cat),
        Metric('Annual sales', lambda: func['get_annual_sales'](tea, ppl), unit, cat)
        ])

    metrics = add_LCA_metrics(system, metrics, kind)

    return metrics


def update_metrics(model, kind):
    metrics = [i for i in model.metrics if i.element_name!='LCA results']
    model.metrics = add_LCA_metrics(model.system, metrics, kind)
    return model


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

        su_type = type(unit).__name__
        if su_type.lower() == 'lagoon':
            su_type = f'{unit.design_type.capitalize()} lagoon'
        name = f'{su_type} {para}'
        model.parameter(setter=AttrSetter(unit, para),
                        name=name, element=unit,
                        kind='coupled', units=df.loc[para]['unit'],
                        baseline=b, distribution=D)


# %%

# =============================================================================
# Shared by all three systems
# =============================================================================

su_data_path = os.path.join(data_path, 'sanunit_data/')
path = os.path.join(su_data_path, '_drying_bed.tsv')
drying_bed_data = load_data(path)
get_exchange_rate = systems.get_exchange_rate
get_decay_k = systems.get_decay_k
tau_deg = systems.tau_deg
log_deg = systems.log_deg

def add_shared_parameters(model, drying_bed_unit, crop_application_unit):
    ########## Related to multiple units ##########
    sys = model.system
    Excretion, Toilet = sys.path[0], sys.path[1]
    param = model.parameter
    streams = sys_dct['stream_dct'][sys.ID]
    tea = sys_dct['TEA'][sys.ID]

    # UGX-to-USD
    b = get_exchange_rate()
    D = shape.Triangle(lower=3600, midpoint=b, upper=3900)
    @param(name='Exchange rate', element=Excretion, kind='cost', units='UGX/USD',
           baseline=b, distribution=D)
    def set_exchange_rate(i):
        systems.exchange_rate = i

    ########## Related to human input ##########
    # Diet and excretion
    path = data_path + 'sanunit_data/_excretion.tsv'
    excretion_data = load_data(path)
    batch_setting_unit_params(excretion_data, model, Excretion)

    # Household size
    b = systems.household_size
    D = shape.Normal(mu=b, sigma=1.8)
    @param(name='Household size', element=Excretion, kind='coupled', units='cap/household',
           baseline=b, distribution=D,
            hook=lambda i: max(1, i)
           )
    def set_household_size(i):
        systems.household_size = i

    # Toilet
    b = systems.household_per_toilet
    D = shape.Uniform(lower=3, upper=5)
    @param(name='Toilet density', element=Toilet, kind='coupled', units='household/toilet',
           baseline=b, distribution=D)
    def set_toilet_density(i):
        systems.household_per_toilet = i

    path = data_path + 'sanunit_data/_toilet.tsv'
    toilet_data = load_data(path)
    batch_setting_unit_params(toilet_data, model, Toilet,
                              exclude=('desiccant_rho',)) # set separately

    toilet_type = type(Toilet).__name__
    WoodAsh = systems.cmps.WoodAsh
    b = V_to_rho(WoodAsh.V(298.15), WoodAsh.MW)
    D = shape.Triangle(lower=663, midpoint=b, upper=977)
    @param(name=f'{toilet_type} desiccant density', element=Toilet, kind='coupled',
           units='kg/m3', baseline=b, distribution=D)
    def set_desiccant_density(i):
        WoodAsh.V.local_methods['USER_METHOD'].value = rho_to_V(i, WoodAsh.MW)
        setattr(Toilet, 'desiccant_rho', i)

    b = WoodAsh.i_Mg
    D = shape.Triangle(lower=0.008, midpoint=b, upper=0.0562)
    @param(name=f'{toilet_type} desiccant Mg content', element=Toilet, kind='coupled',
           units='fraction', baseline=b, distribution=D)
    def set_desiccant_Mg(i):
        WoodAsh.i_Mg = i

    b = WoodAsh.i_Ca
    D = shape.Triangle(lower=0.0742, midpoint=b, upper=0.3716)
    @param(name=f'{toilet_type} desiccant Ca content', element=Toilet, kind='coupled',
           units='fraction', baseline=b, distribution=D)
    def set_desiccant_Ca(i):
        WoodAsh.i_Ca = i

    ##### Universal degradation parameters #####
    # Max methane emission
    unit = sys.path[1] # the first unit that involves degradation
    b = systems.max_CH4_emission
    D = shape.Triangle(lower=0.175, midpoint=b, upper=0.325)
    @param(name='Max CH4 emission', element=unit, kind='coupled', units='g CH4/g COD',
           baseline=b, distribution=D)
    def set_max_CH4_emission(i):
        systems.max_CH4_emission = i
        for unit in sys.units:
            if hasattr(unit, 'max_CH4_emission'):
                setattr(unit, 'max_CH4_emission', i)

    # Time to full degradation
    b = tau_deg
    D = shape.Uniform(lower=1, upper=3)
    @param(name='Full degradation time', element=unit, kind='coupled', units='yr',
           baseline=b, distribution=D)
    def set_tau_deg(i):
        tau_deg = i
        k = get_decay_k(tau_deg, log_deg)
        for unit in sys.units:
            if hasattr(unit, 'decay_k_COD'):
                setattr(unit, 'decay_k_COD', k)
            if hasattr(unit, 'decay_k_N'):
                setattr(unit, 'decay_k_N', k)

    # Reduction at full degradation
    b = systems.log_deg
    D = shape.Uniform(lower=2, upper=4)
    @param(name='Log degradation', element=unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_log_deg(i):
        systems.log = i
        k = get_decay_k(tau_deg, log_deg)
        for unit in sys.units:
            if hasattr(unit, 'decay_k_COD'):
                setattr(unit, 'decay_k_COD', k)
            if hasattr(unit, 'decay_k_N'):
                setattr(unit, 'decay_k_N', k)

    ##### Toilet material properties #####
    density = unit.density_dct
    b = density['Plastic']
    D = shape.Uniform(lower=0.31, upper=1.24)
    param(setter=DictAttrSetter(unit, 'density_dct', 'Plastic'),
          name='Plastic density', element=unit, kind='isolated', units='kg/m2',
          baseline=b, distribution=D)

    b = density['Brick']
    D = shape.Uniform(lower=1500, upper=2000)
    param(setter=DictAttrSetter(unit, 'density_dct', 'Brick'),
          name='Brick density', element=unit, kind='isolated', units='kg/m3',
          baseline=b, distribution=D)

    b = density['StainlessSteelSheet']
    D = shape.Uniform(lower=2.26, upper=3.58)
    param(setter=DictAttrSetter(unit, 'density_dct', 'StainlessSteelSheet'),
          name='SS sheet density', element=unit, kind='isolated', units='kg/m2',
          baseline=b, distribution=D)

    b = density['Gravel']
    D = shape.Uniform(lower=1520, upper=1680)
    param(setter=DictAttrSetter(unit, 'density_dct', 'Gravel'),
          name='Gravel density', element=unit, kind='isolated', units='kg/m3',
          baseline=b, distribution=D)

    b = density['Sand']
    D = shape.Uniform(lower=1281, upper=1602)
    param(setter=DictAttrSetter(unit, 'density_dct', 'Sand'),
          name='Sand density', element=unit, kind='isolated', units='kg/m3',
          baseline=b, distribution=D)

    b = density['Steel']
    D = shape.Uniform(lower=7750, upper=8050)
    param(setter=DictAttrSetter(unit, 'density_dct', 'Steel'),
          name='Steel density', element=unit, kind='isolated', units='kg/m3',
          baseline=b, distribution=D)

    ########## Drying bed ##########
    unit = drying_bed_unit
    D = shape.Uniform(lower=0, upper=0.1)
    batch_setting_unit_params(drying_bed_data, model, unit, exclude=('sol_frac', 'bed_H'))

    b = unit.sol_frac
    if unit.design_type == 'unplanted':
        D = shape.Uniform(lower=0.3, upper=0.4)
    elif unit.design_type == 'planted':
        D = shape.Uniform(lower=0.4, upper=0.7)
    param(setter=DictAttrSetter(unit, '_sol_frac', getattr(unit, 'design_type')),
          name='sol_frac', element=unit, kind='coupled', units='fraction',
          baseline=b, distribution=D)

    b = unit.bed_H['covered']
    D = shape.Uniform(lower=0.45, upper=0.75)
    param(setter=DictAttrSetter(unit, 'bed_H', ('covered', 'uncovered')),
          name='non_storage_bed_H', element=unit, kind='coupled', units='m',
          baseline=b, distribution=D)

    b = unit.bed_H['storage']
    D = shape.Uniform(lower=1.2, upper=1.8)
    param(DictAttrSetter(unit, 'bed_H', 'storage'),
          name='storage_bed_H', element=unit, kind='coupled', units='m',
          baseline=b, distribution=D)

    ########## Crop application ##########
    unit = crop_application_unit
    D = shape.Uniform(lower=0, upper=0.1)
    param(setter=DictAttrSetter(unit, 'loss_ratio', 'NH3'),
          name='NH3 application loss', element=unit, kind='coupled',
          units='fraction of applied', baseline=0.05, distribution=D)

    # Mg, Ca, C actually not affecting results
    D = shape.Uniform(lower=0, upper=0.05)
    param(setter=DictAttrSetter(unit, 'loss_ratio', ('NonNH3', 'P', 'K', 'Mg', 'Ca')),
          name='Other application losses', element=unit, kind='coupled',
          units='fraction of applied', baseline=0.02, distribution=D)

    # ######## Equipment lifetime ######## # added to test function, not included in Trimmer et al., 2020
    # for u in sys.units:
    #     if u.lifetime:
    #         if isinstance(u.lifetime, int): # add the lifetime of the unit
    #             b = u.lifetime
    #             D = shape.Uniform(lower=b*(1-0.25), upper=b*(1+0.25))
    #             param(setter=AttrSetter(u, 'lifetime'),
    #                   name=f'{u} lifetime', element=u, kind='isolated',
    #                   units='yr', baseline=b, distribution=D)
    #         else:
    #             for equip, lifetime in u._default_equipment_lifetime.items(): # add lifetime of all equipment
    #                 b = lifetime
    #                 D = shape.Uniform(lower=b*(1-0.25), upper=b*(1+0.25))
    #                 param(setter=DictAttrSetter(u, 'lifetime', equip),
    #                       name=f'{equip} lifetime', element=u, kind='isolated',
    #                       units='yr', baseline=b, distribution=D)


    ######## General TEA settings ########
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
        price_dct['N'] = streams['liq_N'].price = streams['sol_N'].price = i * get_price_factor()

    D = shape.Uniform(lower=2.619, upper=6.692)
    @param(name='P fertilizer price', element='TEA', kind='isolated', units='USD/kg P',
           baseline=3.983, distribution=D)
    def set_P_price(i):
        price_dct['P'] = streams['liq_P'].price = streams['sol_P'].price = i * get_price_factor()

    D = shape.Uniform(lower=1.214, upper=1.474)
    @param(name='K fertilizer price', element='TEA', kind='isolated', units='USD/kg K',
           baseline=1.333, distribution=D)
    def set_K_price(i):
        price_dct['K'] = streams['liq_K'].price = streams['sol_K'].price = i * get_price_factor()

    # Money discount rate
    b = systems.discount_rate
    D = shape.Uniform(lower=0.03, upper=0.06)
    @param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
           baseline=b, distribution=D)
    def set_discount_rate(i):
        systems.discount_rate = tea.discount_rate = i

    # Electricity price
    b = price_dct['Electricity']
    D = shape.Triangle(lower=0.08, midpoint=b, upper=0.21)
    @param(name='Electricity price', element='TEA', kind='isolated',
           units='$/kWh', baseline=b, distribution=D)
    def set_electricity_price(i):
        PowerUtility.price = i

    return model

get_biogas_factor = systems.get_biogas_factor
def add_LCA_CF_parameters(model, kind=bw._lca_data.lca_data_kind):
    # global lca_param_kind
    # lca_param_kind = kind
    param = model.parameter
    sys = model.system
    lca = sys_dct['LCA'][sys.ID]

    ######## LCA CF ########
    if kind == 'original':
        b = GWP_dct['CH4']
        D = shape.Uniform(lower=28, upper=34)
        @param(name='CH4 CF', element='LCA', kind='isolated', units='kg CO2-eq/kg CH4',
               baseline=b, distribution=D)
        def set_CH4_CF(i):
            GWP_dct['CH4'] = ImpactItem.get_item('CH4_item').CFs['GlobalWarming'] = i

        b = GWP_dct['N2O']
        D = shape.Uniform(lower=265, upper=298)
        @param(name='N2O CF', element='LCA', kind='isolated', units='kg CO2-eq/kg N2O',
               baseline=b, distribution=D)
        def set_N2O_CF(i):
            GWP_dct['N2O'] = ImpactItem.get_item('N2O_item').CFs['GlobalWarming'] = i

        b = GWP_dct['Electricity']
        D = shape.Uniform(lower=0.106, upper=0.121)
        @param(name='Electricity CF', element='LCA', kind='isolated',
               units='kg CO2-eq/kWh', baseline=b, distribution=D)
        def set_electricity_CF(i):
            GWP_dct['Electricity'] = ImpactItem.get_item('E_item').CFs['GlobalWarming'] = i

        b = -GWP_dct['N']
        D = shape.Triangle(lower=1.8, midpoint=b, upper=8.9)
        @param(name='N fertilizer CF', element='LCA', kind='isolated',
               units='kg CO2-eq/kg N', baseline=b, distribution=D)
        def set_N_fertilizer_CF(i):
            GWP_dct['N'] = ImpactItem.get_item('N_item').CFs['GlobalWarming'] = -i

        b = -GWP_dct['P']
        D = shape.Triangle(lower=4.3, midpoint=b, upper=5.4)
        @param(name='P fertilizer CF', element='LCA', kind='isolated',
               units='kg CO2-eq/kg P', baseline=b, distribution=D)
        def set_P_fertilizer_CF(i):
            GWP_dct['P'] = ImpactItem.get_item('P_item').CFs['GlobalWarming'] = -i

        b = -GWP_dct['K']
        D = shape.Triangle(lower=1.1, midpoint=b, upper=2)
        @param(name='K fertilizer CF', element='LCA', kind='isolated',
               units='kg CO2-eq/kg K', baseline=b, distribution=D)
        def set_K_fertilizer_CF(i):
            GWP_dct['K'] = ImpactItem.get_item('K_item').CFs['GlobalWarming'] = -i

        item_path = os.path.join(bw._lca_data.data_path, 'items_original.xlsx')
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
            model.parameter(name=p+'CF',
                            setter=DictAttrSetter(item, 'CFs', 'GlobalWarming'),
                            element='LCA', kind='isolated',
                            units=f'kg CO2-eq/{item.functional_unit}',
                            baseline=b, distribution=D)

        if sys.ID == 'sysB':
            D = shape.Uniform(lower=2.93, upper=3.05)
            @param(name='Liquid petroleum gas CF', element='LCA', kind='isolated', units='MJ/kg',
                   baseline=3, distribution=D)
            def set_LPG_CF(i):
                GWP_dct['Biogas'] = ImpactItem.get_item('Biogas_item').CFs['GlobalWarming'] = \
                    -i*get_biogas_factor()

    else:
        item_path = os.path.join(bw.data_path, 'cf_dct.pckl')
        f = open(item_path, 'rb')
        cf_dct = pickle.load(f)
        f.close()

        ind_new = load_data(os.path.join(bw._lca_data.data_path, 'indicators_new.tsv'))

        for p, df in cf_dct.items():
            item = ImpactItem.get_item(p)
            for ind in lca.indicators:
                full_name = ind_new[ind_new['indicator']==ind.ID]['full_name'].values.item().split("'")
                column = (full_name[1], full_name[3], full_name[5])
                ind_data = df[column]

                b = ind_data[df[df[('-', '-', 'activity name')]=='mean'].index].values.item()
                lower = ind_data[df[df[('-', '-', 'activity name')]=='min'].index].values.item()
                upper = ind_data[df[df[('-', '-', 'activity name')]=='max'].index].values.item()
                # All triangular distribution, set to baseline ±10% if only one data entry from ecoinvent
                name = f'{p} {ind.ID} CF'
                if lower==upper:
                    name += ' (±10%)'
                    lower = b * 0.9
                    upper = b * 1.1

                if item.ID in ('N_item', 'P_item', 'K_item'):
                    lower, b, upper = upper, b, lower
                elif item.ID == 'Biogas_item':
                    if not sys.ID == 'sysB':
                        continue
                    factor = get_biogas_factor()
                    lower, b, upper = upper*factor, b*factor, lower*factor
                D = shape.Triangle(lower=lower, midpoint=b, upper=upper)
                model.parameter(name=f'{p} {ind.ID} CF',
                                setter=DictAttrSetter(item, 'CFs', ind.ID),
                                element='LCA', kind='isolated',
                                units=f'{ind.unit}/{item.functional_unit}',
                                baseline=b, distribution=D)

    return model

def update_LCA_CF_parameters(model, kind):
    non_lca_params = [i for i in model.parameters if not 'CF' in i.name]
    model.set_parameters(non_lca_params)
    model = add_LCA_CF_parameters(model, kind)
    return model


# =============================================================================
# For the same processes in sysA and sysB
# =============================================================================

path = os.path.join(su_data_path, '_pit_latrine.tsv')
pit_latrine_data = load_data(path)

MCF_lower_dct = dct_from_str(pit_latrine_data.loc['MCF_decay']['low'])
MCF_upper_dct = dct_from_str(pit_latrine_data.loc['MCF_decay']['high'])
N2O_EF_lower_dct = dct_from_str(pit_latrine_data.loc['N2O_EF_decay']['low'])
N2O_EF_upper_dct = dct_from_str(pit_latrine_data.loc['N2O_EF_decay']['high'])

def add_pit_latrine_parameters(model):
    sys = model.system
    unit = sys.path[1]
    param = model.parameter
    ######## Related to the toilet ########
    batch_setting_unit_params(pit_latrine_data, model, unit,
                              exclude=('MCF_decay', 'N2O_EF_decay'))

    # MCF and N2O_EF decay parameters, specified based on the type of the pit latrine
    b = unit.MCF_decay
    kind = unit._return_MCF_EF()
    D = shape.Triangle(lower=MCF_lower_dct[kind], midpoint=b, upper=MCF_upper_dct[kind])
    param(setter=DictAttrSetter(unit, '_MCF_decay', kind),
          name='Pit latrine MCF decay', element=unit, kind='coupled',
          units='fraction of anaerobic conversion of degraded COD',
          baseline=b, distribution=D)

    b = unit.N2O_EF_decay
    D = shape.Triangle(lower=N2O_EF_lower_dct[kind], midpoint=b, upper=N2O_EF_upper_dct[kind])
    param(setter=DictAttrSetter(unit, '_N2O_EF_decay', kind),
          name='Pit latrine N2O EF decay', element=unit, kind='coupled',
          units='fraction of N emitted as N2O',
          baseline=b, distribution=D)

    # Costs
    b = unit.CAPEX
    D = shape.Uniform(lower=386, upper=511)
    param(setter=AttrSetter(unit, 'CAPEX'),
          name='Pit latrine capital cost', element=unit, kind='cost',
          units='USD/toilet', baseline=b, distribution=D)

    b = unit.OPEX_over_CAPEX
    D = shape.Uniform(lower=0.02, upper=0.08)
    param(setter=AttrSetter(unit, 'OPEX_over_CAPEX'),
          name='Pit latrine annual operating cost', element=unit, kind='cost',
          units='fraction of capital cost', baseline=b, distribution=D)

    ######## Related to conveyance ########
    unit = sys.path[2]
    b = unit.loss_ratio
    D = shape.Uniform(lower=0.02, upper=0.05)
    param(setter=AttrSetter(unit, 'loss_ratio'),
          name='Transportation loss', element=unit, kind='coupled', units='fraction',
          baseline=b, distribution=D)

    b = unit.single_truck.distance
    D = shape.Uniform(lower=2, upper=10)
    param(setter=AttrSetter(unit.single_truck, 'distance'),
          name='Transportation distance', element=unit, kind='coupled', units='km',
          baseline=b, distribution=D)

    b = systems.emptying_fee
    D = shape.Uniform(lower=0, upper=0.3)
    @param(name='Additional emptying fee', element=unit, kind='coupled', units='fraction of base cost',
           baseline=b, distribution=D)
    def set_emptying_fee(i):
        systems.emptying_fee = i

    return model

path = su_data_path + '_sludge_separator.tsv'
sludge_separator_data = load_data(path)
split_lower_dct = dct_from_str(sludge_separator_data.loc['split']['low'])
split_upper_dct = dct_from_str(sludge_separator_data.loc['split']['high'])
split_dist_dct = dct_from_str(sludge_separator_data.loc['split']['distribution'], dtype='str')

def add_sludge_separator_parameters(unit, model):
    param = model.parameter

    b = unit.settled_frac
    D = shape.Uniform(lower=0.1, upper=0.2)
    @param(name='Settled frac', element=unit, kind='coupled', units='fraction',
           baseline=b, distribution=D)
    def set_settled_frac(i):
        unit.settled_frac = i

    for key in split_lower_dct.keys():
        b = getattr(unit, 'split')[key]
        lower = split_lower_dct[key]
        upper = split_upper_dct[key]
        dist = split_dist_dct[key]
        if dist == 'uniform':
            D = shape.Uniform(lower=lower, upper=upper)
        elif dist == 'triangular':
            D = shape.Triangle(lower=lower, midpoint=b, upper=upper)
        param(setter=DictAttrSetter(unit, 'split', key),
              name='Frac of settled'+key, element=unit, kind='coupled',
              units='fraction',
              baseline=b, distribution=D)

    return model

def add_lagoon_parameters(unit, model):
    param = model.parameter
    b = systems.sewer_flow
    D = shape.Uniform(lower=2500, upper=3000)
    name = f'{unit.design_type.capitalize()} lagoon sewer flow'
    @param(name=name, element=unit, kind='coupled', units='m3/d',
           baseline=b, distribution=D)
    def set_sewer_flow(i):
        systems.sewer_flow = i
    return model

def add_existing_plant_parameters(toilet_unit, cost_unit, tea, model):
    param = model.parameter
    b = systems.ppl_exist_sewer
    D = shape.Uniform(lower=3e4, upper=5e4)
    @param(name='Sewer ppl', element=toilet_unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_sewer_ppl(i):
        systems.ppl_exist_sewer = i

    b = systems.ppl_exist_sludge
    D = shape.Triangle(lower=416667, midpoint=b, upper=458333)
    @param(name='Exist sludge ppl', element=toilet_unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_sludge_ppl(i):
        systems.ppl_exist_sludge = i

    b = cost_unit.lifetime
    D = shape.Triangle(lower=8, midpoint=b, upper=11)
    param(setter=AttrSetter(cost_unit, 'lifetime'),
          name='Plant lifetime', element='TEA/LCA', kind='isolated', units='yr',
          baseline=b, distribution=D)

    b = 3e6
    D = shape.Uniform(lower=1e6, upper=5e6)
    param(setter=AttrFuncSetter(tea, 'annual_labor',
                                lambda salary: salary*12*12/get_exchange_rate()),
          name='Staff salary', element='TEA', kind='isolated', units='MM UGX/cap/month',
          baseline=b, distribution=D)

    return model


# %%

# =============================================================================
# Scenario A (sysA)
# =============================================================================

sysA = systems.sysA
sysA.simulate()
modelA = Model(sysA, add_metrics(sysA, lca_data_kind))
paramA = modelA.parameter

# Shared parameters
modelA = add_shared_parameters(modelA, systems.A8, systems.A9)
modelA = add_LCA_CF_parameters(modelA)

# Pit latrine and conveyance
modelA = add_pit_latrine_parameters(modelA)

# WWTP costs
modelA = add_existing_plant_parameters(systems.A2, systems.A4, systems.teaA, modelA)

# Sedimentation tank
A5 = systems.A5
path = os.path.join(su_data_path, '_sedimentation_tank.tsv')
data = load_data(path)
batch_setting_unit_params(data, modelA, A5)
# The tank was based on a sludge separator
modelA = add_sludge_separator_parameters(A5, modelA)

# Anaerobic lagoon
A6 = systems.A6
path = os.path.join(su_data_path, '_anaerobic_lagoon.tsv')
anaerobic_lagoon_data = load_data(path)
batch_setting_unit_params(anaerobic_lagoon_data, modelA, A6)
modelA = add_lagoon_parameters(A6, modelA)

# Facultative lagoon
A7 = systems.A7
path = os.path.join(su_data_path, '_facultative_lagoon.tsv')
facultative_lagoon_data = load_data(path)
batch_setting_unit_params(facultative_lagoon_data, modelA, A7)
modelA = add_lagoon_parameters(A7, modelA)

all_paramsA = modelA.parameters


# %%

# =============================================================================
# Scenario B (sysB)
# =============================================================================

sysB, teaB = systems.sysB, systems.teaB
sysB.simulate()
modelB = Model(sysB, add_metrics(sysB, lca_data_kind))
paramB = modelB.parameter

# Shared parameters
modelB = add_shared_parameters(modelB, systems.B8, systems.B9)
modelB = add_LCA_CF_parameters(modelB)

# Pit latrine and conveyance
modelB = add_pit_latrine_parameters(modelB)

b = systems.ppl_alt
D = shape.Triangle(lower=45e3, midpoint=b, upper=55e3)
@paramB(name='Alt sludge ppl', element=systems.B2, kind='coupled', units='-',
        baseline=b, distribution=D)
def set_plant_ppl(i):
    systems.ppl_alt = i

# Anaerobic baffled reactor
B5 = systems.B5
path = os.path.join(su_data_path, '_anaerobic_baffled_reactor.tsv')
data = load_data(path)
batch_setting_unit_params(data, modelB, B5)

b = systems.biogas_energy
D = shape.Triangle(lower=802, midpoint=b, upper=870)
@paramB(name='Biogas energy', element=B5, kind='coupled', units='kJ/mol CH4',
        baseline=b, distribution=D)
def set_biogas_energy(i):
    systems.biogas_energy = i

# Cost of alternative plants
B4 = systems.B4
b = B4.CAPEX_dct['Lumped WWTP']
D = shape.Triangle(lower=303426, midpoint=b, upper=370854)
@paramB(name='Plant CAPEX', element=B4, kind='cost', units='USD',
        baseline=b, distribution=D)
def set_alt_plant_CAPEX(i):
    B4.CAPEX_dct['Lumped WWTP'] = i

b = B4.lifetime
D = shape.Triangle(lower=9, midpoint=b, upper=11)
@paramB(name='Plant lifetime', element='TEA/LCA', kind='isolated', units='yr',
        baseline=b, distribution=D)
def set_plant_lifetime(i):
    B4.lifetime = i

b = systems.get_unskilled_num()
D = shape.Uniform(lower=0, upper=10)
@paramB(name='Unskilled staff num', element='TEA', kind='isolated', units='-',
        baseline=b, distribution=D)
def set_unskilled_num(i):
    systems.unskilled_num = i
    teaB.annual_labor = systems.get_alt_salary()


b = 0.75e6
D = shape.Uniform(lower=0.5e6, upper=1e6)
@paramB(name='Unskilled staff salary', element='TEA', kind='isolated', units='MM UGX/cap/month',
        baseline=b, distribution=D)
def set_unskilled_salary(i):
    systems.unskilled_salary = i
    teaB.annual_labor = systems.get_alt_salary()

# Sludge separator
B6 = systems.B6
modelB = add_sludge_separator_parameters(B6, modelB)

# Liquid treatment bed
B7 = systems.B7
path = os.path.join(su_data_path, '_liquid_treatment_bed.tsv')
data = load_data(path)
batch_setting_unit_params(data, modelB, B7)

# Biogas combustion
B14 = systems.B14
b = B14.biogas_loss
D = shape.Uniform(lower=0, upper=0.2)
@paramB(name='Biogas loss ratio', element=B14, kind='coupled', units='fraction',
        baseline=b, distribution=D)
def set_biogas_loss(i):
    B14.biogas_loss = i

D = shape.Uniform(lower=6077, upper=6667)
@paramB(name='Liquid petroleum gas price', element='TEA', kind='isolated', units='UGX/kg',
        baseline=6500, distribution=D)
def set_LPG_price(i):
    price_dct['Biogas'] = sys_dct['stream_dct']['sysB']['biogas'].price = \
        i/get_exchange_rate()*get_biogas_factor()

b = systems.LPG_energy
D = shape.Uniform(lower=49.5, upper=50.4)
@paramB(name='Liquid petroleum gas energy', element='TEA/LCA', kind='isolated', units='MJ/kg',
        baseline=b, distribution=D)
def set_LPG_energy(i):
    old_LPG_energy = systems.LPG_energy
    systems.LPG_energy = i
    price_dct['Biogas'] = sys_dct['stream_dct']['sysB']['biogas'].price = \
        price_dct['Biogas'] / old_LPG_energy * i

all_paramsB = modelB.parameters


# %%

# =============================================================================
# Scenario C (sysC)
# =============================================================================

sysC = systems.sysC
sysC.simulate()
modelC = Model(sysC, add_metrics(sysC, lca_data_kind))
paramC = modelC.parameter

# Add shared parameters
modelC = add_shared_parameters(modelC, systems.C8, systems.C9)
modelC = add_LCA_CF_parameters(modelC)

# UDDT
C2 = systems.C2
path = os.path.join(su_data_path, '_uddt.tsv')
uddt_data = load_data(path)
batch_setting_unit_params(uddt_data, modelC, C2)

b = C2.CAPEX
D = shape.Uniform(lower=476, upper=630)
@paramC(name='UDDT capital cost', element=C2, kind='cost',
       units='USD/toilet', baseline=b, distribution=D)
def set_UDDT_CAPEX(i):
    C2.CAPEX = i

b = C2.OPEX_over_CAPEX
D = shape.Uniform(lower=0.05, upper=0.1)
@paramC(name='UDDT annual operating cost', element=C2, kind='cost',
       units='fraction of capital cost', baseline=b, distribution=D)
def set_UDDT_OPEX(i):
    C2.OPEX_over_CAPEX = i

# Conveyance
C3 = systems.C3
C4 = systems.C4
b = C3.loss_ratio
D = shape.Uniform(lower=0.02, upper=0.05)
@paramC(name='Transportation loss', element=C3, kind='coupled', units='fraction',
       baseline=b, distribution=D)
def set_trans_loss(i):
    C3.loss_ratio = C4.loss_ratio = i

b = C3.single_truck.distance
D = shape.Uniform(lower=2, upper=10)
@paramC(name='Transportation distance', element=C3, kind='coupled', units='km',
       baseline=b, distribution=D)
def set_trans_distance(i):
    C3.single_truck.distance = C4.single_truck.distance = i

b = systems.handcart_fee
D = shape.Uniform(lower=0.004, upper=0.015)
@paramC(name='Handcart fee', element=C3, kind='cost', units='USD/cap/d',
       baseline=b, distribution=D)
def set_handcart_fee(i):
    systems.handcart_fee = i

b = systems.truck_fee
D = shape.Uniform(lower=17e3, upper=30e3)
@paramC(name='Truck fee', element=C3, kind='cost', units='UGX/m3',
       baseline=b, distribution=D)
def set_truck_fee(i):
    systems.truck_fee = i

# WWTP costs
modelC = add_existing_plant_parameters(systems.C2, systems.C5, systems.teaC, modelC)

# Anaerobic lagoon
C6 = systems.C6
batch_setting_unit_params(anaerobic_lagoon_data, modelC, C6)
modelC = add_lagoon_parameters(C6, modelC)

# Facultative lagoon
C7 = systems.C7
batch_setting_unit_params(facultative_lagoon_data, modelC, C7)
modelC = add_lagoon_parameters(C7, modelC)

all_paramsC = modelC.parameters



# %%

# =============================================================================
# Functions to run simulation and generate plots
# =============================================================================

result_dct = {
        'sysA': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        'sysB': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        'sysC': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        }

@time_printer
def run_uncertainty(model, seed=None, N=1000, rule='L',
                    percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                    spearman_metrics='default'):
    global result_dct
    if seed:
        np.random.seed(seed)

    samples = model.sample(N, rule)

    model.load_samples(samples)
    model.evaluate()

    # Spearman's rank correlation,
    # metrics default to net cost, net emission, and total recoveries
    spearman_results = None
    if spearman_metrics:
        if spearman_metrics.lower() == 'default':
            spearman_metrics = [i for i in model.metrics
                                if 'net' in i.name.lower() or 'total' in i.name.lower()]

        spearman_results = model.spearman_r(model.parameters, spearman_metrics)[0]
        spearman_results.columns = pd.Index([i.name_with_units for i in spearman_metrics])

    dct = organize_uncertainty_results(model, percentiles, spearman_results)
    return dct


# Data organization
def organize_uncertainty_results(model, percentiles, spearman_results):
    dct = result_dct[model._system.ID]
    index_p = len(model.parameters)
    dct['parameters'] = model.table.iloc[:, :index_p].copy()
    dct['data'] = model.table.iloc[:, index_p:].copy()

    if percentiles:
        dct['percentiles'] = dct['data'].quantile(q=percentiles)

    if spearman_results:
        dct['spearman'] = spearman_results
    return dct


def save_uncertainty_results(model, path=''):
    if not path:
        path = os.path.join(c_path, 'results')

        if not os.path.isdir(path):
            os.mkdir(path)
        path = os.path.join(path, f'model{model._system.ID[-1]}.xlsx')

    elif not (path.endswith('xlsx') or path.endswith('xls')):
        extension = path.split('.')[-1]
        raise ValueError(f'Only "xlsx" and "xls" are supported, not {extension}.')

    dct = result_dct[model._system.ID]
    if dct['parameters'] is None:
        raise ValueError('No cached result, run model first.')
    with pd.ExcelWriter(path) as writer:
        dct['parameters'].to_excel(writer, sheet_name='Parameters')
        dct['data'].to_excel(writer, sheet_name='Uncertainty results')
        if 'percentiles' in dct.keys():
            dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
        dct['spearman'].to_excel(writer, sheet_name='Spearman')
        model.table.to_excel(writer, sheet_name='Raw data')