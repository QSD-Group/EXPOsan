#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

TODO: add a function to pick model tables and reload model
'''


# %%

import os
import numpy as np, pandas as pd
from qsdsan import stats as s
from qsdsan.utils import time_printer, copy_samples, colors
from exposan import bwaise as bw
from exposan.bwaise import results_path, figures_path

# Comment these out if want to see all warnings
import warnings
warnings.filterwarnings(action='ignore')

modelA, modelB, modelC = bw.modelA, bw.modelB, bw.modelC
seed = 3221 # for numpy seeding and result consistency


# %%

# Net cost, net GWP, and total COD/N/P/K recovery
def get_key_metrics(model):
    key_metrics = [i for i in model.metrics if 'net' in i.name.lower()]
    key_metrics += [i for i in model.metrics if 'total' in i.name.lower()]
    return key_metrics


def filter_parameters(model, df, threshold):
    new_df = pd.concat((df[df>=threshold], df[df<=-threshold]))
    filtered = new_df.dropna(how='all')
    param_dct = {p.name_with_units:p for p in model.get_parameters()}
    parameters = set(param_dct[i[1]] for i in filtered.index)
    return list(parameters)


@time_printer
def evaluate(model, samples=None):
    if samples is not None:
        model.load_samples(samples)
    model.evaluate()


# %%

# =============================================================================
# Uncertainty analysis
# =============================================================================

# Save system, TEA, and LCA reports
# note that system and TEA reports will be in one Excel file,
# while LCA report will be in a standalone Excel file.
def save_reports(system):
    system.simulate()
    system.save_report(os.path.join(results_path, f'{system.ID}.xlsx'))
    system.LCA.save_report(os.path.join(results_path, f'{system.ID}_lca.xlsx'))


# Plot recoveries as 1D-density plots
def plot_box(model, color, metrics, kind='horizontal'):
    if kind == 'horizontal':
        fig, ax = s.plot_uncertainties(model,
                                       y_axis=metrics,
                                       kind='box',
                                       center_kws={'color': color,
                                                   'width': 0.6})
        ax.get_legend().remove()
        fig.set_figheight(2.5)
        ax.set_yticklabels([i.name.lstrip('Total ') for i in metrics], fontsize=20)
        ax.set_xlim(0, 1)
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.set_xticklabels([f'{i:.0%}' for i in np.linspace(0, 1, 6)], fontsize=20)

    else:
        fig, ax = s.plot_uncertainties(model,
                                       x_axis=metrics,
                                       kind='box',
                                       center_kws={'color': color,
                                                   'width': 0.6})
        ax.get_legend().remove()
        fig.set_figwidth(2.5)
        ax.set_xticklabels([i.name.lstrip('Total ') for i in metrics], fontsize=20)
        for label in ax.get_xticklabels():
            label.set_rotation(30)
        ax.set_ylim(0, 1)
        ax.set_yticks(np.linspace(0, 1, 6))
        ax.set_yticklabels([f'{i:.0%}' for i in np.linspace(0, 1, 6)], fontsize=20)

    fig.subplots_adjust(bottom=0.25)

    return fig, ax


# Plot net cost and emissions as 2D-density plots
def plot_kde(model, color, metrics):
    fig, ax = s.plot_uncertainties(model,
                                   x_axis=metrics[0],
                                   y_axis=metrics[1],
                                   kind='kde-box',
                                   center_kws={'fill': True, 'color': color},
                                   margin_kws={'color': color})
    for i in (ax[0], ax[1]):
        i.set_xlim(0, 100)
        i.set_xticks((np.linspace(0, 100, 6)))
    for i in (ax[0], ax[2]):
        i.set_ylim(0, 150)
        i.set_yticks((np.linspace(0, 150, 6)))
    for txt in ax[0].get_xticklabels()+ax[0].get_yticklabels():
        txt.set_fontsize(20)
    fig.subplots_adjust(left=0.15)
    return fig, ax



# %%

# =============================================================================
# Code used to generate figures and data used in the manuscript
# =============================================================================

if __name__ == '__main__':
    # # This is the new, recommended method for setting seed,
    # # but not seems to be widely used/supported
    # import numpy as np
    # rng = np.random.default_rng(3221)
    # rng.random(5)
    np.random.seed(seed)

    # This saves reports for the system, TEA, and LCA
    for sys in (bw.sysA, bw.sysB, bw.sysC):
        save_reports(sys)

    ########## Uncertainty analysis ##########
    figs, axs = {}, {}
    for model in (modelA, modelB, modelC):
        model.metrics = key_metrics = get_key_metrics(model)
        ID = model.system.ID[-1]
        figs[ID], axs[ID] = {}, {}

        samples = model.sample(N=1000, rule='L', seed=seed)
        model.load_samples(samples)
        if model is modelB:
            copy_samples(modelA, model)
            color = colors.Guest.green.RGBn
        elif model is modelC:
            copy_samples(modelA, model)
            copy_samples(modelB, model, exclude=modelA.parameters)
            color = colors.Guest.blue.RGBn
        else:
            color = colors.Guest.orange.RGBn

        evaluate(model)

        figs[ID]['box'], axs[ID]['box'] = \
            plot_box(model, color, key_metrics[2:], 'horizontal')
        figs[ID]['box'].savefig(os.path.join(figures_path, f'recoveries{ID}.png'), dpi=300)


        figs[ID]['kde'], axs[ID]['kde'] = plot_kde(model, color, key_metrics[:2])
        figs[ID]['kde'].savefig(os.path.join(figures_path, f'cost_emission{ID}.png'), dpi=300)

    ########## Morris One-at-A-Time ##########
    mu_star_origin_dct = {}
    mu_star_norm_dct = {}
    for model in (modelA, modelB, modelC):
        inputs = s.define_inputs(model)
        ID = model.system.ID[-1]
        morris_samples = s.generate_samples(inputs, kind='Morris', N=10, seed=seed)

        evaluate(model, morris_samples)

        morris_dct = s.morris_analysis(model, inputs, seed=seed,
                                       nan_policy='fill_mean',
                                       file=os.path.join(results_path, f'Morris{ID}.xlsx'))

        origin = []
        filtered = []
        for i in model.metrics:
            df = morris_dct[i.name]
            df.sort_values(by=['mu_star'], ascending=False, inplace=True)
            origin.append(df.mu_star)
            df_filtered = df.mu_star.iloc[0:5] # select the top five
            filtered.append(df_filtered)

        mu_star_origin = pd.concat(origin, axis=1)
        mu_star_filtered = pd.concat(filtered, axis=1)
        mu_star_origin.columns = mu_star_filtered.columns = [i.name for i in model.metrics]
        mu_star_norm = mu_star_filtered/mu_star_filtered.max() # normalize

        mu_star_origin_dct[f'{model.system.ID}'] = mu_star_origin
        mu_star_norm_dct[f'{model.system.ID}'] = mu_star_norm

    columns = []
    data = []
    for i in key_metrics:
        for ID, df in mu_star_norm_dct.items():
            columns.append(f'{i.name}-{ID}')
            data.append(df[i.name])

    combined = pd.concat(data, axis=1)
    combined.columns = columns

    writer = pd.ExcelWriter(os.path.join(results_path, 'Morris.xlsx'))
    combined.to_excel(writer, sheet_name='Combined')
    for ID, df in mu_star_origin_dct.items():
        df.to_excel(writer, sheet_name=ID)
    writer.save()


# %%

########## Below are functions used for tutorial purpose ##########

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
        suffix = model._system.ID[-1]
        conv = '_conv' if test_convergence else ''
        suffix += conv
    else:
        suffix = ''

    if file_prefix=='default':
        suffix = model._system.ID[-1]
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


# %%

# =============================================================================
# Pearson and Spearman
# =============================================================================

def run_plot_spearman(model, N, seed=seed, metrics=None, parameters=None,
                      file_prefix='default'):
    suffix = model._system.ID[-1] if file_prefix=='default' else ''
    metrics = metrics if metrics else get_key_metrics(model)

    if file_prefix=='default':
        suffix = model._system.ID[-1]
        dct_file = os.path.join(results_path, f'Spearman{suffix}.xlsx')
        fig_file = os.path.join(figures_path, f'Spearman{suffix}.png')
    else:
        dct_file = f'{file_prefix}.xlsx' if file_prefix else ''
        fig_file = f'{file_prefix}.png' if file_prefix else ''

    spearman_rho, spearman_p = s.get_correlations(model, kind='Spearman',
                                                  nan_policy='raise',
                                                  file=dct_file)

    fig, ax = s.plot_correlations(spearman_rho, parameters=parameters,
                                  metrics=metrics, file=fig_file)

    return spearman_rho, fig, ax


# %%

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

    suffix = model._system.ID[-1] if file_prefix=='default' else ''
    if file_prefix=='default':
        suffix = model._system.ID[-1]
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


# %%

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
        suffix = model._system.ID[-1]
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

        #!!! Need to figure out parameter filtering
        # if auto_filter_parameters:
        #     key_params = filter_parameters(model, spearman_rho, threshold)
        #     if len(key_params) > 5:
        #         model.set_parameters(key_params)

        # # Filter parameters based on Spearman's rho
        # parameters = filter_parameters(model, spearman_rho, 0.5)
        # fig, ax = s.plot_correlations(spearman_rho, parameters=parameters,
        #                               metrics=key_metrics[0])


        #!!!! Update the `stats` module based on the following,
        # and point the users to look at the documentation on the module
        ########## Uncomment if want to run Spearman and plot results ##########
        # spearman_rho, fig, ax, all_params = run_plot_spearman(model, N=100)

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