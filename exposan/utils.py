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
from chaospy import distributions as shape
from qsdsan import ImpactItem, sanunits as su
from qsdsan.utils import time_printer, AttrSetter
from . import es_path



__all__ = (
    'add_fugitive_items',
    'batch_setting_unit_params',
    'clear_unit_costs',
    'run_uncertainty',
    )


def add_fugitive_items(unit, item_ID):
    unit._run()
    for i in unit.ins:
        i.stream_impact_item = ImpactItem.get_item(item_ID).copy(set_as_source=True)


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