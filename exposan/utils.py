#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License. Please refer to
https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np, pandas as pd
from math import log
from sklearn.linear_model import LinearRegression as LR
from chaospy import distributions as shape
<<<<<<< HEAD
=======
from thermosteam.functional import rho_to_V
>>>>>>> 39af7d1ccc5478a59bd0bf765ea48bcbcf694bb0
from qsdsan import ImpactItem, sanunits as su
from qsdsan.utils import time_printer, AttrSetter
from . import es_path



__all__ = (
    'add_fugitive_items',
<<<<<<< HEAD
=======
    'add_V_from_rho',
>>>>>>> 39af7d1ccc5478a59bd0bf765ea48bcbcf694bb0
    'batch_setting_unit_params',
    'clear_unit_costs',
    'get_decay_k',
    'get_generic_scaled_capital',
    'get_generic_tanker_truck_fee',
    'run_uncertainty',
    )


def add_fugitive_items(unit, item_ID):
    unit._run()
    for i in unit.ins:
        i.stream_impact_item = ImpactItem.get_item(item_ID).copy(set_as_source=True)


<<<<<<< HEAD
=======
def add_V_from_rho(cmp, rho):
    V_model = rho_to_V(rho, cmp.MW)
    try: cmp.V.add_model(V_model)
    except:
        handle = getattr(cmp.V, cmp.locked_state)
        handle.add_model(V_model)


>>>>>>> 39af7d1ccc5478a59bd0bf765ea48bcbcf694bb0
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


# Costs of WWTP units have been considered in the lumped unit
def clear_unit_costs(sys):
    for i in sys.units:
        if isinstance(i, su.LumpedCost): continue
        i.purchase_costs.clear()
        i.installed_costs.clear()


# Get reduction rate constant k for COD and N, use a function so that k can be
# changed during uncertainty analysis
def get_decay_k(tau_deg=2, log_deg=3):
    k = (-1/tau_deg)*np.log(10**-log_deg)
    return k


def get_generic_scaled_capital(tea, percent_CAPEX_to_scale, number_of_units,
                               percent_limit, learning_curve_percent):
    '''
    Scale capital cost based for the Nth system
    (would be lower than the cost for a single system due to scaling effect).

    Parameters
    ----------
    tea : obj
        TEA obj for the system of interest.
    percent_CAPEX_to_scale : float
        The fraction of the cost of specialty parts/cost of total parts.
    number_of_units : int
        Number of units to be constructed.
    percent_limit : float
        Percent of the lowest cost of the normal cost of a single system.
    learning_curve_percent : float
        The percent factor of the learning curve.
    '''
    CAPEX_to_scale = tea.annualized_CAPEX * percent_CAPEX_to_scale
    CAPEX_not_scaled = tea.annualized_CAPEX - CAPEX_to_scale
    scaled_limited = CAPEX_to_scale * percent_limit
    b = log(learning_curve_percent)/log(2)
    scaled_CAPEX_annualized  = (CAPEX_to_scale - scaled_limited)*number_of_units**b + scaled_limited
    new_CAPEX_annualized = scaled_CAPEX_annualized + CAPEX_not_scaled
    return new_CAPEX_annualized


fitting_dct = {
    3: 21.62,
    4.5: 32.43,
    8: 54.05,
    15: 67.57,
}
def get_generic_tanker_truck_fee(capacity,
                                 fitting_dct=fitting_dct,
                                 emptying_fee=0.15,
                                 exchange_rate=1):
    '''
    Exponential fitting to get the tanker truck fee based on capacity.

    cost = a*capacity**b -> ln(price) = ln(a) + bln(capacity)

    Parameters
    ----------
    capacity : float
        The capacity at which the tanker truck fee will be calculated.
    fitting_dct : dict(float, float)
        Capacity-based cost to develop the exponential fitting correlation,
        keys should be the capacities and values should be the corresponding costs.
        Capacities for fitting.
    emptying_fee : float
        Additional fraction of fee that will be added on top of the given prices.
    exchange_rate : float
        Exchange that will be multiplied to the prices.
    '''
    capacities = np.array(tuple(fitting_dct.keys()))
    costs = np.array(tuple(fitting_dct.values()))
    costs *= (1+emptying_fee)*exchange_rate
    ln_p = np.log(costs)
    ln_cap = np.log(np.array(capacities))
    model = LR().fit(ln_cap.reshape(-1,1), ln_p.reshape(-1,1))
    predicted = model.predict(np.array((np.log(capacity))).reshape(1, -1)).item()
    fee = np.exp(predicted)
    return fee


@time_printer
def run_uncertainty(model, seed=None, N=1000, rule='L',
                    percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                    path='', print_time=False):
    if seed: np.random.seed(seed)

    samples = model.sample(N, rule)
    model.load_samples(samples)
    model.evaluate()
    # Spearman's rank correlation
    spearman_results = model.spearman()
    spearman_results.columns = pd.Index([i.name_with_units for i in model.metrics])
    organize_and_save_results(model=model, percentiles=percentiles,
                              spearman_results=spearman_results, path=path)


def organize_and_save_results(
        model, percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
        spearman_results=None, path=''):
    dct = {}
    index_p = len(model.get_parameters())
    dct['parameters'] = model.table.iloc[:, :index_p].copy()
    dct['data'] = model.table.iloc[:, index_p:].copy()
    if percentiles:
        dct['percentiles_parameters'] = dct['parameters'].quantile(q=percentiles)
        dct['percentiles_results'] = dct['data'].quantile(q=percentiles)
    dct['spearman'] = spearman_results

    path = os.path.join(es_path, f'sys{model.system.ID[-1]}_model.xlsx') if path=='' else path
    with pd.ExcelWriter(path) as writer:
        dct['parameters'].to_excel(writer, sheet_name='Parameters')
        dct['data'].to_excel(writer, sheet_name='Uncertainty results')
        if percentiles:
            dct['percentiles_parameters'].to_excel(writer, sheet_name='Parameter percentiles')
            dct['percentiles_results'].to_excel(writer, sheet_name='Result percentiles')
        if spearman_results is not None: dct['spearman'].to_excel(writer, sheet_name='Spearman')
        model.table.to_excel(writer, sheet_name='Raw data')