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

import os
import numpy as np, pandas as pd
from qsdsan import stats as s
from qsdsan.utils import time_printer, copy_samples, colors, load_pickle, save_pickle
from exposan import bwaise as bw
from exposan.bwaise import results_path, figures_path
from adjustText import adjust_text

# Comment these out if want to see all warnings
import warnings
warnings.filterwarnings(action='ignore')

modelA, modelB, modelC = bw.modelA, bw.modelB, bw.modelC
RGBs = {
    'A': colors.Guest.orange.RGBn,
    'B': colors.Guest.green.RGBn,
    'C': colors.Guest.blue.RGBn
    }    
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
def plot_box(model, ID, color, metrics, kind='horizontal'):
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
    fig.savefig(os.path.join(figures_path, f'recoveries{ID}.png'), dpi=300)


# Plot net cost and emissions as 2D-density plots
def plot_kde(model, ID, color, metrics):
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
    fig.savefig(os.path.join(figures_path, f'cost_emission{ID}.png'), dpi=300)


param_set = set(modelA.parameters).union(set(modelB.parameters), set(modelC.parameters))
param_names = tuple(p.name for p in param_set)

# Plot morris results
#!!! Make it possible to have 
def plot_morris(model, ID, morris_dct, color, annotate=True):
    ID = model.system.ID[-1]
    normalized = {}
    k3s = {}
    for metric in model.metrics:
        df = morris_dct[metric.name].copy()

        #!!! Remove this once results are refreshed
        try: int(df.index[0])
        except:
            df.reset_index(inplace=True)
            df.rename(columns={'index':'parameter'}, inplace=True)

        labels = [param_names.index(p) for p in df.parameter]
        df['label'] = labels
        
        k3s[metric.name] = round(df.mu_star.max()/df.sigma.max(), 1)
        df.mu_star_conf /= df.mu_star.max()
        df.mu_star /= df.mu_star.max()
        df.sigma /= df.sigma.max()
        normalized[metric.name] = df
    
        fig, ax = s.plot_morris_results(normalized, metric, label_kind=None,
                                        k1=k3s[metric.name]*0.1, k2=None,
                                        k3=k3s[metric.name],
                                        color=color)
        fig.set(figheight=3, figwidth=3)

        if annotate:
            legend = ax.get_legend()
            if legend:
                legend.texts[0].set_text('$\\sigma/\\mu^*$=1')
                legend.texts[1].set_text('$\\sigma/\\mu^*$=0.1')
            xlabel = r'$\mu^*$/$\mu^*_{max}$'
            ylabel = r'$\sigma$/$\sigma_{max}$'
        else:
            legend = ax.get_legend()
            if legend:
                legend.remove()
            xlabel = ylabel = ''

        ticks = [0, 0.25, 0.5, 0.75, 1]
        ax.set(xlim=(0, 1), ylim=(0, 1), xbound=(0, 1.1), ybound=(0, 1.1),
               xticks=ticks, yticks=ticks, xticklabels=ticks, yticklabels=ticks,
               xlabel=xlabel, ylabel=ylabel)
    
        # Only label the ones with mu_star/mu_star_max>0.1 or sigma/sigma_max>0.1
        top = df[(df.mu_star>=0.1)|(df.sigma>=0.1)]
        
        # Only markout the top five mu_star/mu_star_max parameters
        if top.shape[0] > 5:
            top.sort_values(by=['mu_star'], ascending=False, inplace=True)
            top = top[:5]
            
        labels = []
        for x, y, label in zip(top.mu_star, top.sigma, top.label):
            labels.append(ax.text(x, y, label))
        adjust_text(labels)
        
        # # Alternative codes that doesn't work as well as the adjustText package
        # for x, y, label in zip(top.mu_star, top.sigma, top.label):
        #     ax.annotate(label, (x, y), 
        #                 xytext=(2, 2), textcoords='offset points',
        #                 ha='center')
        
        
        fig.subplots_adjust(bottom=0.2, left=0.25)
        path = os.path.join(figures_path, f'Morris{ID}_{metric.name}.png')
        fig.savefig(path, dpi=300)
        

def run(from_record=True, annotate_morris=True):
    # # This is the new, recommended method for setting seed,
    # # but not seems to be widely used/supported
    # import numpy as np
    # rng = np.random.default_rng(3221)
    # rng.random(5)
    np.random.seed(seed)

    # A dict with all results for future rebuilding of the models
    global table_dct
    if not from_record:
        # This saves reports for the system, TEA, and LCA
        for sys in (bw.sysA, bw.sysB, bw.sysC):
            save_reports(sys)
        table_dct = dict(uncertainty={}, morris={})
    else:
        path = os.path.join(results_path, 'table_dct.pckl')
        table_dct = load_pickle(path)
        
    ########## Uncertainty analysis ##########
    for model in (modelA, modelB, modelC):
        model.metrics = key_metrics = get_key_metrics(model)
        ID = model.system.ID[-1]

        if not from_record:
            # Want to use a larger N (maybe 5000 or 10000)
            samples = model.sample(N=5000, rule='L', seed=seed)
            model.load_samples(samples)
            if model is modelB:
                copy_samples(modelA, model)
            elif model is modelC:
                copy_samples(modelA, model)
                copy_samples(modelB, model, exclude=modelA.parameters)
    
            evaluate(model)
            table_dct['uncertainty'][ID] = model.table.copy()
        else:
            model.table = table_dct['uncertainty'][ID]

        # Make the plots
        plot_box(model, ID, RGBs[ID], key_metrics[2:], 'horizontal')
        plot_kde(model, ID, RGBs[ID], key_metrics[:2])
        

    ########## Morris One-at-A-Time ##########
    mu_star_origin_dct = {}
    mu_star_norm_dct = {}
    for model in (modelA, modelB, modelC):
        ID = model.system.ID[-1]
        if not from_record:
            inputs = s.define_inputs(model)
            # Want to use a larger N (maybe 100)
            morris_samples = s.generate_samples(inputs, kind='Morris', N=100, seed=seed)
    
            evaluate(model, morris_samples)
    
            # These are the unprocessed data
            morris_dct = s.morris_analysis(model, inputs, seed=seed,
                                           nan_policy='fill_mean',
                                           file=os.path.join(results_path, f'Morris{ID}.xlsx'))
            
            table_dct['morris'][ID] = model.table.copy()
            table_dct['morris_dct'][ID] = morris_dct.copy()

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
        else:
            model.table = table_dct['morris'][ID]
            morris_dct = table_dct['morris_dct'][ID]

        # Make the plots
        plot_morris(model, ID, morris_dct, RGBs[ID], annotate_morris)
        
    # Data processing
    if not from_record:
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
        
        pickle_path = os.path.join(results_path, 'table_dct.pckl')
        save_pickle(table_dct, pickle_path)


# %%

# =============================================================================
# Code used to generate figures and data used in the manuscript
# =============================================================================

if __name__ == '__main__':
    run(from_record=True, annotate_morris=False)
    # run(from_record=True, annotate_morris=True)


# %%

########## Below are functions used for tutorial purpose ##########

#!!! Need to review and update with new qsdsan

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