#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import os, pickle
from chaospy import distributions as shape
from thermosteam.functional import V_to_rho, rho_to_V
from qsdsan import currency, ImpactItem, PowerUtility, Model, Metric
from qsdsan.utils import (
    AttrFuncSetter,
    AttrSetter,
    DictAttrSetter,
    data_path,
    dct_from_str,
    load_data,
    ospath,
    )
from exposan.utils import batch_setting_unit_params, run_uncertainty as run
from exposan import bwaise as bw
from exposan.bwaise import (
    _load_components,
    create_system,
    get_alt_salary,
    get_biogas_factor,
    get_decay_k,
    get_LCA_metrics,
    get_TEA_metrics,
    get_recoveries,
    GWP_dct,
    price_dct,
    results_path,
    )

__all__ = ('create_model', 'run_uncertainty',)


# %%

# =============================================================================
# Functions for batch-making metrics and -setting parameters
# =============================================================================

def add_LCA_metrics(model, lca_kind):
    bw._load_lca_data(lca_kind)
    system = model.system
    lca = system.LCA
    metrics = []
    for ind in lca.indicators:
        unit = f'{ind.unit}/cap/yr'
        cat = 'LCA results'
        funcs = get_LCA_metrics(system, ind)
        metrics.extend([
            Metric(f'Net emission {ind.ID}', funcs[0], unit, cat),
            Metric(f'Construction {ind.ID}', funcs[1], unit, cat),
            Metric(f'Transportation {ind.ID}', funcs[2], unit, cat),
            Metric(f'Direct emission {ind.ID}', funcs[3], unit, cat),
            Metric(f'Offset {ind.ID}', funcs[4], unit, cat),
            Metric(f'Other {ind.ID}', funcs[5], unit, cat),
            ])
    model.metrics = [*model.metrics, *metrics]


def add_metrics(model, lca_kind):
    system = model.system
    metrics = []
    for i in ('COD', 'N', 'P', 'K'):
        funcs = get_recoveries(system, i)
        cat = f'{i} recovery'
        metrics.extend([
            Metric(f'Liquid {i}', funcs[0], '', cat),
            Metric(f'Solid {i}', funcs[1], '', cat),
            Metric(f'Gas {i}', funcs[2], '', cat),
            Metric(f'Total {i}', funcs[3], '', cat)
            ])

    unit = f'{currency}/cap/yr'
    cat = 'TEA results'
    funcs = get_TEA_metrics(system)
    metrics.extend([
        Metric('Annual net cost', funcs[0], unit, cat),
        Metric('Annual CAPEX', funcs[1], unit, cat),
        Metric('Annual OPEX', funcs[2], unit, cat),
        Metric('Annual sales', funcs[3], unit, cat)
        ])
    model.metrics = metrics
    add_LCA_metrics(model, lca_kind)


def update_metrics(model, lca_kind):
    model.metrics = [i for i in model.metrics if i.element_name!='LCA results']
    add_LCA_metrics(model, lca_kind)


# %%

# =============================================================================
# Data sheets
# =============================================================================

join_path = lambda prefix, file_name: os.path.join(prefix, file_name)

su_data_path = ospath.join(data_path, 'sanunit_data')
excretion_data = load_data(join_path(su_data_path, '_excretion.tsv'))
toilet_data = load_data(join_path(su_data_path, '_toilet.tsv'))
pit_latrine_data = load_data(join_path(su_data_path, '_pit_latrine.tsv'))
uddt_data = load_data(join_path(su_data_path, '_uddt.tsv'))
drying_bed_data = load_data(join_path(su_data_path, '_drying_bed.tsv'))
liquid_bed_data = load_data(join_path(su_data_path, '_liquid_treatment_bed.tsv'))
sedimentation_tank_data = load_data(join_path(su_data_path, '_sedimentation_tank.tsv'))
anaerobic_lagoon_data = load_data(join_path(su_data_path, '_anaerobic_lagoon.tsv'))
facultative_lagoon_data = load_data(join_path(su_data_path, '_facultative_lagoon.tsv'))
sludge_separator_data = load_data(join_path(su_data_path, '_sludge_separator.tsv'))
abr_data = load_data(join_path(su_data_path, '_anaerobic_baffled_reactor.tsv'))


# %%

# =============================================================================
# Shared by all three systems
# =============================================================================

MCF_lower_dct = dct_from_str(pit_latrine_data.loc['MCF_decay']['low'])
MCF_upper_dct = dct_from_str(pit_latrine_data.loc['MCF_decay']['high'])
N2O_EF_lower_dct = dct_from_str(pit_latrine_data.loc['N2O_EF_decay']['low'])
N2O_EF_upper_dct = dct_from_str(pit_latrine_data.loc['N2O_EF_decay']['high'])

def add_shared_parameters(model, drying_bed_unit, main_crop_application_unit):
    ##### Related to multiple units #####
    sys = model.system
    sys_stream = sys.flowsheet.stream
    tea = sys.TEA
    param = model.parameter

    Excretion, Toilet = sys.path[0], sys.path[1]

    # UGX-to-USD
    b = bw.exchange_rate
    D = shape.Triangle(lower=3600, midpoint=b, upper=3900)
    @param(name='Exchange rate', element=Excretion, kind='cost', units='UGX/USD',
           baseline=b, distribution=D)
    def set_exchange_rate(i):
        bw.exchange_rate = i

    ##### Related to human input #####
    # Diet and excretion
    batch_setting_unit_params(excretion_data, model, Excretion)

    # Household size
    b = bw.household_size
    D = shape.Trunc(shape.Normal(mu=b, sigma=1.8), lower=1)
    @param(name='Household size', element=Toilet, kind='coupled', units='cap/household',
           baseline=b, distribution=D)
    def set_household_size(i):
        bw.household_size = i

    # Toilet
    b = bw.household_per_toilet
    D = shape.Uniform(lower=3, upper=5)
    @param(name='Toilet density', element=Toilet, kind='coupled', units='household/toilet',
           baseline=b, distribution=D)
    def set_toilet_density(i):
        bw.household_per_toilet = i

    batch_setting_unit_params(toilet_data, model, Toilet,
                              exclude=('desiccant_rho',)) # set separately

    toilet_type = type(Toilet).__name__
    WoodAsh = bw.components.WoodAsh
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
    b = bw.max_CH4_emission
    D = shape.Triangle(lower=0.175, midpoint=b, upper=0.325)
    @param(name='Max CH4 emission', element=unit, kind='coupled', units='g CH4/g COD',
           baseline=b, distribution=D)
    def set_max_CH4_emission(i):
        bw.max_CH4_emission = i
        for unit in sys.units:
            if hasattr(unit, 'max_CH4_emission'):
                setattr(unit, 'max_CH4_emission', i)

    # Time to full degradation
    b = bw.tau_deg
    D = shape.Uniform(lower=1, upper=3)
    @param(name='Full degradation time', element=unit, kind='coupled', units='yr',
           baseline=b, distribution=D)
    def set_tau_deg(i):
        bw.tau_deg = i
        k = get_decay_k(i, bw.log_deg)
        for unit in sys.units:
            if hasattr(unit, 'decay_k_COD'):
                setattr(unit, 'decay_k_COD', k)
            if hasattr(unit, 'decay_k_N'):
                setattr(unit, 'decay_k_N', k)

    # Reduction at full degradation
    b = bw.log_deg
    D = shape.Uniform(lower=2, upper=4)
    @param(name='Log degradation', element=unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_log_deg(i):
        bw.log_deg = i
        k = get_decay_k(bw.tau_deg, i)
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

    ##### Drying bed #####
    unit = drying_bed_unit
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

    ##### Crop application #####
    unit = main_crop_application_unit
    D = shape.Uniform(lower=0, upper=0.1)
    param(setter=DictAttrSetter(unit, 'loss_ratio', 'NH3'),
          name='NH3 application loss', element=unit, kind='coupled',
          units='fraction of applied', baseline=0.05, distribution=D)

    # Mg, Ca, C actually not affecting results
    D = shape.Uniform(lower=0, upper=0.05)
    param(setter=DictAttrSetter(unit, 'loss_ratio', ('NonNH3', 'P', 'K', 'Mg', 'Ca')),
          name='Other application losses', element=unit, kind='coupled',
          units='fraction of applied', baseline=0.02, distribution=D)

    # ##### Equipment lifetime #####
    # # DO NOT DELETE
    # # Added to test function, not included in Trimmer et al., 2020
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


    ##### General TEA settings #####
    # Discount factor for the excreta-derived fertilizers
    b = bw.price_factor
    D = shape.Uniform(lower=0.1, upper=0.4)
    @param(name='Price factor', element='TEA', kind='isolated', units='-',
           baseline=b, distribution=D)
    def set_price_factor(i):
        bw.price_factor = i

    D = shape.Uniform(lower=1.164, upper=2.296)
    @param(name='N fertilizer price', element='TEA', kind='isolated', units='USD/kg N',
           baseline=1.507, distribution=D)
    def set_N_price(i):
        price_dct['N'] = sys_stream.liq_N.price = sys_stream.sol_N.price = i * bw.price_factor

    D = shape.Uniform(lower=2.619, upper=6.692)
    @param(name='P fertilizer price', element='TEA', kind='isolated', units='USD/kg P',
           baseline=3.983, distribution=D)
    def set_P_price(i):
        price_dct['P'] = sys_stream.liq_P.price = sys_stream.sol_P.price = i * bw.price_factor

    D = shape.Uniform(lower=1.214, upper=1.474)
    @param(name='K fertilizer price', element='TEA', kind='isolated', units='USD/kg K',
           baseline=1.333, distribution=D)
    def set_K_price(i):
        price_dct['K'] = sys_stream.liq_K.price = sys_stream.sol_K.price = i * bw.price_factor

    # Money discount rate
    b = bw.discount_rate
    D = shape.Uniform(lower=0.03, upper=0.06)
    @param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
           baseline=b, distribution=D)
    def set_discount_rate(i):
        bw.discount_rate = tea.discount_rate = i

    # Electricity price
    b = price_dct['Electricity']
    D = shape.Triangle(lower=0.08, midpoint=b, upper=0.21)
    @param(name='Electricity price', element='TEA', kind='isolated',
           units='$/kWh', baseline=b, distribution=D)
    def set_electricity_price(i):
        PowerUtility.price = i


def add_LCA_CF_parameters(model, lca_kind):
    param = model.parameter
    sys = model.system
    lca = sys.LCA

    ##### LCA CF #####
    if lca_kind == 'original':
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

        item_path = ospath.join(bw._lca_data.data_path, 'items_original.xlsx')
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

        if sys.ID == 'sysB' and ImpactItem.get_item('Biogas_item'):
            D = shape.Uniform(lower=2.93, upper=3.05)
            @param(name='Liquid petroleum gas CF', element='LCA', kind='isolated', units='MJ/kg',
                   baseline=3, distribution=D)
            def set_LPG_CF(i):
                GWP_dct['Biogas'] = ImpactItem.get_item('Biogas_item').CFs['GlobalWarming'] = \
                    -i*get_biogas_factor()

    else:
        item_path = ospath.join(bw.data_path, 'cf_dct.pckl')
        f = open(item_path, 'rb')
        cf_dct = pickle.load(f)
        f.close()

        ind_new = load_data(ospath.join(bw._lca_data.data_path, 'indicators_new.tsv'))

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


def update_LCA_CF_parameters(model, lca_kind):
    non_lca_params = [i for i in model.parameters if not ' CF' in i.name] # just 'CF' will exclude 'MCF' as well
    model.set_parameters(non_lca_params)
    add_LCA_CF_parameters(model, lca_kind)


# =============================================================================
# For the same processes in sysA and sysB
# =============================================================================

def add_pit_latrine_parameters(model):
    sys = model.system
    unit = sys.path[1]
    param = model.parameter
    ##### Related to the toilet #####
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

    ##### Related to conveyance #####
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

    b = bw.emptying_fee
    D = shape.Uniform(lower=0, upper=0.3)
    @param(name='Additional emptying fee', element=unit, kind='coupled', units='fraction of base cost',
           baseline=b, distribution=D)
    def set_emptying_fee(i):
        bw.emptying_fee = i


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


def add_lagoon_parameters(unit, model):
    param = model.parameter
    b = bw.sewer_flow
    D = shape.Uniform(lower=2500, upper=3000)
    name = f'{unit.design_type.capitalize()} lagoon sewer flow'
    @param(name=name, element=unit, kind='coupled', units='m3/d',
           baseline=b, distribution=D)
    def set_sewer_flow(i):
        bw.sewer_flow = i


def add_existing_plant_parameters(toilet_unit, cost_unit, model):
    param = model.parameter
    b = bw.ppl_exist_sewer
    D = shape.Uniform(lower=3e4, upper=5e4)
    @param(name='Sewer ppl', element=toilet_unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_sewer_ppl(i):
        bw.ppl_exist_sewer = i

    b = bw.ppl_exist_sludge
    D = shape.Triangle(lower=375000, midpoint=b, upper=458333)
    @param(name='Exist sludge ppl', element=toilet_unit, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_sludge_ppl(i):
        bw.ppl_exist_sludge = i

    b = cost_unit.lifetime
    D = shape.Triangle(lower=8, midpoint=b, upper=11)
    param(setter=AttrSetter(cost_unit, 'lifetime'),
          name='Plant lifetime', element='TEA/LCA', kind='isolated', units='yr',
          baseline=b, distribution=D)

    tea = model.system.TEA
    b = 3e6
    D = shape.Uniform(lower=1e6, upper=5e6)
    param(setter=AttrFuncSetter(tea, 'annual_labor',
                                lambda salary: salary*12*12/bw.exchange_rate),
          name='Staff salary', element='TEA', kind='isolated', units='MM UGX/cap/month',
          baseline=b, distribution=D)


# %%

# =============================================================================
# Scenario A (sysA)
# =============================================================================

def create_modelA(lca_kind='original', **model_kwargs):
    sysA = create_system('A', lca_kind=lca_kind)
    unitA = sysA.flowsheet.unit

    # Add shared metrics/parameters
    modelA = Model(sysA, **model_kwargs)
    add_metrics(modelA, lca_kind)
    add_shared_parameters(modelA, unitA.A8, unitA.A9)
    add_LCA_CF_parameters(modelA, lca_kind=lca_kind)

    # Pit latrine and conveyance
    add_pit_latrine_parameters(modelA)

    # WWTP costs
    add_existing_plant_parameters(unitA.A2, unitA.A4, modelA)

    # Sedimentation tank
    A5 = unitA.A5
    batch_setting_unit_params(sedimentation_tank_data, modelA, A5)
    # The tank was based on a sludge separator
    add_sludge_separator_parameters(A5, modelA)

    # Anaerobic lagoon
    A6 = unitA.A6
    batch_setting_unit_params(anaerobic_lagoon_data, modelA, A6)
    add_lagoon_parameters(A6, modelA)

    # Facultative lagoon
    A7 = unitA.A7
    batch_setting_unit_params(facultative_lagoon_data, modelA, A7)
    add_lagoon_parameters(A7, modelA)

    # # DO NOT DELETE
    # # Legacy codes to look at recoveries
    # A1 = unitA.A1
    # get_recovery = bw.get_recovery
    # metricsA = [m for m in modelA.metrics]
    # metricsA.extend([
    #     # Metric(f'Net emission {ind.ID}', FuncGetter(funcs[0], (ind.ID,)), unit, cat),
    #     Metric('A1', lambda: get_recovery(A1, systems.A2.ins, get_ppl('a'))['N'], '%', 'N'),
    #     Metric('A2', lambda: get_recovery(A1, systems.A3.ins, get_ppl('a'))['N'], '%', 'N'),
    #     Metric('A3', lambda: get_recovery(A1, systems.A4.ins, get_ppl('a'))['N'], '%', 'N'),
    #     Metric('A4', lambda: get_recovery(A1, systems.A5.ins, get_ppl('a'))['N'], '%', 'N'),
    #     Metric('A5', lambda: get_recovery(A1, systems.A6.ins, get_ppl('a'))['N'], '%', 'N'),
    #     Metric('A6', lambda: get_recovery(A1, systems.A7.ins, get_ppl('a'))['N'], '%', 'N'),
    #     Metric('A7', lambda: get_recovery(A1, systems.A8.ins, get_ppl('a'))['N'], '%', 'N'),
    #     Metric('A8', lambda: get_recovery(A1, systems.A12.ins, get_ppl('a'))['N'], '%', 'N'),
    #     Metric('A9', lambda: get_recovery(A1, systems.A13.ins, get_ppl('a'))['N'], '%', 'N'),
    #     ])
    # modelA.metrics = metricsA

    return modelA


# %%

# =============================================================================
# Scenario B (sysB)
# =============================================================================

def create_modelB(lca_kind='original', **model_kwargs):
    sysB = create_system('B', lca_kind=lca_kind)
    unitB = sysB.flowsheet.unit
    teaB = sysB.TEA

    # Add shared metrics/parameters
    modelB = Model(sysB, **model_kwargs)
    paramB = modelB.parameter
    add_metrics(modelB, lca_kind)
    add_shared_parameters(modelB, unitB.B8, unitB.B9)
    add_LCA_CF_parameters(modelB, lca_kind=lca_kind)

    # Pit latrine and conveyance
    add_pit_latrine_parameters(modelB)

    b = bw.ppl_alt
    D = shape.Triangle(lower=45e3, midpoint=b, upper=55e3)
    @paramB(name='Alt sludge ppl', element=unitB.B2, kind='coupled', units='-',
            baseline=b, distribution=D)
    def set_plant_ppl(i):
        bw.ppl_alt = i

    # Anaerobic baffled reactor
    B5 = unitB.B5
    batch_setting_unit_params(abr_data, modelB, B5)

    b = bw.biogas_energy
    D = shape.Triangle(lower=802, midpoint=b, upper=870)
    @paramB(name='Biogas energy', element=B5, kind='coupled', units='kJ/mol CH4',
            baseline=b, distribution=D)
    def set_biogas_energy(i):
        bw.biogas_energy = i

    # Cost of alternative plants
    B4 = unitB.B4
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

    b = bw.unskilled_num
    D = shape.Uniform(lower=0, upper=10)
    @paramB(name='Unskilled staff num', element='TEA', kind='isolated', units='-',
            baseline=b, distribution=D)
    def set_unskilled_num(i):
        bw.unskilled_num = i
        teaB.annual_labor = get_alt_salary()

    b = bw.unskilled_salary
    D = shape.Uniform(lower=0.5e6, upper=1e6)
    @paramB(name='Unskilled staff salary', element='TEA', kind='isolated', units='MM UGX/cap/month',
            baseline=b, distribution=D)
    def set_unskilled_salary(i):
        bw.unskilled_salary = i
        teaB.annual_labor = get_alt_salary()

    # Sludge separator
    add_sludge_separator_parameters(unitB.B6, modelB)

    # Liquid treatment bed
    B7 = unitB.B7
    batch_setting_unit_params(liquid_bed_data, modelB, B7)

    # Biogas combustion
    B14 = unitB.B14
    b = B14.biogas_loss
    D = shape.Uniform(lower=0, upper=0.2)
    @paramB(name='Biogas loss ratio', element=B14, kind='coupled', units='fraction',
            baseline=b, distribution=D)
    def set_biogas_loss(i):
        B14.biogas_loss = i

    biogas = sysB.flowsheet.stream.biogas
    D = shape.Uniform(lower=6077, upper=6667)
    @paramB(name='Liquid petroleum gas price', element='TEA', kind='isolated', units='UGX/kg',
            baseline=6500, distribution=D)
    def set_LPG_price(i):
        price_dct['Biogas'] = biogas.price = i/bw.exchange_rate*get_biogas_factor()

    b = bw.LPG_energy
    D = shape.Uniform(lower=49.5, upper=50.4)
    @paramB(name='Liquid petroleum gas energy', element='TEA/LCA', kind='isolated', units='MJ/kg',
            baseline=b, distribution=D)
    def set_LPG_energy(i):
        old_LPG_energy = bw.LPG_energy
        bw.LPG_energy = i
        price_dct['Biogas'] = biogas.price = price_dct['Biogas'] / old_LPG_energy * i

    return modelB


# %%

# =============================================================================
# Scenario C (sysC)
# =============================================================================

def create_modelC(lca_kind='original', **model_kwargs):
    sysC = create_system('C', lca_kind=lca_kind)
    unitC = sysC.flowsheet.unit

    # Add shared metrics/parameters
    modelC = Model(sysC, **model_kwargs)
    paramC = modelC.parameter
    add_metrics(modelC, lca_kind)
    add_shared_parameters(modelC, unitC.C8, unitC.C9)
    add_LCA_CF_parameters(modelC, lca_kind=lca_kind)

    # UDDT
    C2 = unitC.C2
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
    C3 = unitC.C3
    C4 = unitC.C4
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

    b = bw.handcart_fee
    D = shape.Uniform(lower=0.004, upper=0.015)
    @paramC(name='Handcart fee', element=C3, kind='cost', units='USD/cap/d',
           baseline=b, distribution=D)
    def set_handcart_fee(i):
        bw.handcart_fee = i

    b = bw.truck_fee
    D = shape.Uniform(lower=17e3, upper=30e3)
    @paramC(name='Truck fee', element=C3, kind='cost', units='UGX/m3',
           baseline=b, distribution=D)
    def set_truck_fee(i):
        bw.truck_fee = i

    # WWTP costs
    add_existing_plant_parameters(unitC.C2, unitC.C5, modelC)

    # Anaerobic lagoon
    C6 = unitC.C6
    batch_setting_unit_params(anaerobic_lagoon_data, modelC, C6)
    add_lagoon_parameters(C6, modelC)

    # Facultative lagoon
    C7 = unitC.C7
    batch_setting_unit_params(facultative_lagoon_data, modelC, C7)
    add_lagoon_parameters(C7, modelC)

    return modelC


# Wrapper functions
country_params = {
    'Caloric intake': 'Excretion e cal',
    'Vegetable protein intake': 'Excretion p veg',
    'Animal protein intake': 'Excretion p anim',
    'N fertilizer price': 'N fertilizer price',
    'P fertilizer price': 'P fertilizer price',
    'K fertilizer price': 'K fertilizer price',
    'Food waste ratio': 'Food waste ratio', # not in the original model
    'Price level ratio': 'Price level ratio', # not in the original model
    'Income tax': 'Income tax', # not in the original model
    }
def create_model(model_ID='A', country_specific=False, **model_kwargs):
    _load_components()
    model_ID = model_ID.lstrip('model').lstrip('sys').upper() # so that it'll work for "modelA"/"sysA"/"A"
    if model_ID == 'A': model = create_modelA(**model_kwargs)
    elif model_ID == 'B': model = create_modelB(**model_kwargs)
    elif model_ID == 'C': model = create_modelC(**model_kwargs)
    else: raise ValueError(f'`model_ID` can only be "A", "B", or "C", not "{model_ID}".')

    if country_specific: # add the remaining three more country-specific parameters
        param = model.parameter
        system = model.system

        unit = system.path[0]
        b = unit.waste_ratio
        D = shape.Uniform(lower=b*0.9, upper=b*1.1)
        @param(name='Food waste ratio', element=unit, kind='cost', units='fraction',
               baseline=b, distribution=D)
        def set_food_waste_ratio(i):
            unit.waste_ratio = i

        b = bw.systems.price_ratio
        D = shape.Uniform(lower=b*0.9, upper=b*1.1)
        @param(name='Price level ratio', element='TEA', kind='cost', units='',
               baseline=b, distribution=D)
        def set_price_ratio(i):
            bw.systems.price_ratio = i

        tea = system.TEA
        b = tea.income_tax
        D = shape.Uniform(lower=b*0.9, upper=b*1.1)
        @param(name='Income tax', element='TEA', kind='cost', units='fraction',
               baseline=b, distribution=D)
        def set_income_tax(i):
            tea.income_tax = i

    return model


def run_uncertainty(model, path='', **kwargs):
    kwargs['path'] = os.path.join(results_path, f'sys{model.system.ID[-1]}_model.xlsx') if path=='' else path
    run(model=model, **kwargs)
    return