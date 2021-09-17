#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

This script contains functions and codes used in the tutorial of
the `qsdsan.stats` module.

'''

import os
import pandas as pd
from qsdsan import stats as s
from exposan import bwaise as bw
from exposan.bwaise import results_path, figures_path

evaluate = bw.analyses.evaluate
get_key_metrics = bw.analyses.get_key_metrics

seed = 3221

def filter_parameters(model, df, threshold):
    new_df = pd.concat((df[df>=threshold], df[df<=-threshold]))
    filtered = new_df.dropna(how='all')
    param_dct = {p.name_with_units:p for p in model.get_parameters()}
    parameters = set(param_dct[i[1]] for i in filtered.index)
    return list(parameters)



#!!! Need to review and update with new qsdsan

# =============================================================================
# Pearson and Spearman
# =============================================================================

def run_plot_spearman(model, N, seed=seed, metrics=None, parameters=None,
                      plot_rho_cutoff=0.5, file_prefix='default'):
    suffix = model.system.ID[-1] if file_prefix=='default' else ''
    metrics = metrics if metrics else get_key_metrics(model)

    if file_prefix=='default':
        suffix = model.system.ID[-1]
        dct_file = os.path.join(results_path, f'Spearman{suffix}.xlsx')
        fig_file = os.path.join(figures_path, f'Spearman{suffix}.png')
    else:
        dct_file = f'{file_prefix}.xlsx' if file_prefix else ''
        fig_file = f'{file_prefix}.png' if file_prefix else ''

    m.run_uncertainty(model, N=N, seed=seed, rule='L',
                      percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

    spearman_rho, spearman_p = s.get_correlations(model, kind='Spearman',
                                                  nan_policy='raise',
                                                  file=dct_file)

    parameters = filter_parameters(model, spearman_rho, plot_rho_cutoff)

    fig, ax = s.plot_correlations(spearman_rho, parameters=parameters,
                                  metrics=metrics, file=fig_file)

    return spearman_rho, fig, ax


# =============================================================================
# Morris One-at-A-Time
# =============================================================================

def run_plot_morris(model, N, seed=seed, test_convergence=False,
                    metrics=None, plot_metric=None, parameters=None,
                    file_prefix='default'):
    inputs = s.define_inputs(model)
    metrics = metrics if metrics else get_key_metrics(model)
    plot_metric = plot_metric if plot_metric else metrics[0]
    if parameters:
        model.parameters = parameters

    if file_prefix=='default':
        suffix = model.system.ID[-1]
        conv = '_conv' if test_convergence else ''
        suffix += conv
    else:
        suffix = ''

    if file_prefix=='default':
        suffix = model.system.ID[-1]
        dct_file = os.path.join(results_path, f'Morris{suffix}.xlsx')
        fig_file = os.path.join(figures_path, f'Morris{suffix}.png')
    else:
        dct_file = f'{file_prefix}.xlsx' if file_prefix else ''
        fig_file = f'{file_prefix}.png' if file_prefix else ''

    if not test_convergence:
        morris_samples = s.generate_samples(inputs, kind='Morris', N=N, seed=seed)

        evaluate(model, morris_samples)

        dct = s.morris_analysis(model, inputs, metrics=metrics, seed=seed,
                                nan_policy='fill_mean', file=dct_file)
        fig, ax = s.plot_morris_results(dct, metric=plot_metric, file=fig_file)

    else:
        dct = s.morris_till_convergence(model, inputs, metrics=metrics, seed=seed,
                                        N_max=N, file=dct_file)
        fig, ax = s.plot_morris_convergence(dct, metric=plot_metric, plot_rank=True,
                                            file=fig_file)

    return dct, fig, ax


# =============================================================================
# (e)FAST and RBD-FAST
# =============================================================================

def run_plot_fast(model, kind, N, M, seed=seed, metrics=None, plot_metric=None,
                  parameters=None, file_prefix='default'):
    inputs = s.define_inputs(model)
    metrics = metrics if metrics else get_key_metrics(model)
    plot_metric = plot_metric if plot_metric else metrics[0]
    if parameters:
        model.parameters = parameters

    suffix = model.system.ID[-1] if file_prefix=='default' else ''
    if file_prefix=='default':
        suffix = model.system.ID[-1]
        dct_file = os.path.join(results_path, f'{kind}{suffix}.xlsx')
        fig_file = os.path.join(figures_path, f'{kind}{suffix}.png')
    else:
        dct_file = f'{file_prefix}.xlsx' if file_prefix else ''
        fig_file = f'{file_prefix}.png' if file_prefix else ''

    if kind.upper() in ('FAST', 'EFAST'):
        fast_samples = s.generate_samples(inputs, kind=kind, N=N, M=M, seed=seed)
    else:
        fast_samples = s.generate_samples(inputs, kind=kind, N=N, seed=seed)

    evaluate(model, fast_samples)

    dct = s.fast_analysis(model, inputs, kind=kind, metrics=metrics,
                          M=M, seed=seed, nan_policy='fill_mean', file=dct_file)

    fig, ax = s.plot_fast_results(dct, metric=plot_metric, file=fig_file)
    return dct, fig, ax


# =============================================================================
# Sobol
# =============================================================================

def run_plot_sobol(model, N, seed=seed, metrics=None, plot_metric=None,
                   parameters=None, file_prefix='default'):
    inputs = s.define_inputs(model)
    metrics = metrics if metrics else get_key_metrics(model)
    plot_metric = plot_metric if plot_metric else metrics[0]
    if parameters:
        model.parameters = parameters

    sobol_samples = s.generate_samples(inputs, kind='Sobol', N=N, seed=seed,
                                        calc_second_order=True)

    evaluate(model, sobol_samples)

    if file_prefix=='default':
        suffix = model.system.ID[-1]
        dct_file = os.path.join(results_path, f'Sobol{suffix}.xlsx')
        fig_file = os.path.join(figures_path, f'Sobol{suffix}.png')
    else:
        dct_file = f'{file_prefix}.xlsx' if file_prefix else ''
        fig_file = f'{file_prefix}.png' if file_prefix else ''

    sobol_dct = s.sobol_analysis(model, inputs,
                                 metrics=metrics, seed=seed,
                                 calc_second_order=True, conf_level=0.95,
                                 nan_policy='fill_mean', file=dct_file)

    fig, ax = s.plot_sobol_results(sobol_dct, metric=plot_metric,
                                    error_bar=True, annotate_heatmap=False,
                                    file=fig_file)

    return sobol_dct, fig, ax


#!!!! Update this scripts and tutorial of the `stats`
# and point the users to look at the documentation on the module
if __name__ == '__main__':
    model = bw.modelA
    m = bw.models

    spearman_rho, fig, ax = run_plot_spearman(model, N=100,
                                              file_prefix='/Users/yalinli_cabbi/downloads/a')

    ########## Uncommnet if want to test the convergence of Morris and plot ##########
    # morris_dct_conv, fig, ax = a.run_plot_morris(model, N=100, test_convergence=True)
    # fig, ax = s.plot_morris_convergence(morris_dct_conv, parameters=parameters,
    #                                     metric=key_metrics[0], plot_rank=True)

    ########## Uncomment if want to run FAST and plot results ##########
    # fast_dct, fig, ax = run_plot_fast(model, 'FAST', N=100, M=4)
    # fig, ax = s.plot_fast_results(fast_dct, key_metrics[0])

    ########## Uncomment if want to run RBD-FAST and plot results ##########
    # rbd_dct, fig, ax = run_plot_fast(model, 'RBD', N=100, M=10)
    # fig, ax = s.plot_fast_results(rbd_dct, key_metrics[0])

    ########## Uncomment if want to run Sobol ##########
    # sobol_dct, fig, ax = run_plot_sobol(model, N=10)
    # fig, ax = s.plot_sobol_results(sobol_dct, metric=key_metrics[0], kind='STS2',
    #                                plot_in_diagonal='ST')