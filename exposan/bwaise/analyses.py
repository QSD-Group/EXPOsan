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
import numpy as np, pandas as pd, seaborn as sns, matplotlib.colors as mcolors
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from qsdsan import stats as s
from qsdsan.utils import copy_samples, colors, load_pickle, save_pickle
from exposan import bwaise as bw
from exposan.bwaise import results_path, figures_path, evaluate, get_key_metrics
from exposan.bwaise.models import organize_uncertainty_results, save_uncertainty_results

# Comment these out if want to see all warnings
import warnings
warnings.filterwarnings(action='ignore')

models = modelA, modelB, modelC = bw.modelA, bw.modelB, bw.modelC

RGBs = {
    'A': colors.Guest.orange.RGBn,
    'B': colors.Guest.green.RGBn,
    'C': colors.Guest.blue.RGBn,
    }

# These label names will affect the order of the parameters in Morris plots,
# but not the value
alt_names = {
    # 'Total COD': 'COD recovery',
    # 'Total N': 'N recovery',
    # 'Total P': 'P recovery',
    # 'Total K': 'K recovery',
    'Annual net cost': 'Cost',
    'Net emission GlobalWarming': 'GWP',
    }
# Change both 'Total X' and 'X recovery' to 'X'
strip_label = lambda label: label.lstrip('Total ').rstrip(' recovery')

seed = 3221 # for numpy seeding and result consistency

def plot_minor_ticks_and_grid(ax, nx=None, ny=None,
                              color='k', linestyle='--', linewidth=0.3):
    if nx:
        ax.xaxis.set_minor_locator(AutoMinorLocator(nx))
    if ny:
        ax.yaxis.set_minor_locator(AutoMinorLocator(ny))
    ax.grid(color=color, linestyle=linestyle, linewidth=linewidth)
    return ax



# %%

# =============================================================================
# Uncertainty analysis
# =============================================================================

# Save system, TEA, and LCA reports
# note that system and TEA reports will be in one Excel file,
# while LCA report will be in a standalone Excel file.
def save_reports(system):
    system.simulate()
    system.save_report(os.path.join(results_path, f'{system.ID}_report.xlsx'))
    system.LCA.save_report(os.path.join(results_path, f'{system.ID}_lca.xlsx'))


# Plot recoveries as 1D-density plots
def plot_box(model, ID, color, metrics, ax_dct, kind='horizontal', adjust_hue=False,
             whis=(5, 95), sym=''):
    if kind == 'horizontal':
        fig, ax = s.plot_uncertainties(model, y_axis=metrics, kind='box',
                                       adjust_hue=adjust_hue,
                                       center_kws={'color': color, 'width': 0.5,
                                                   'whis': whis, 'sym': sym})
        if adjust_hue:
            ax.get_legend().remove()
        fig.set_figheight(2.5)
        ticks = np.arange(0, 1.2, 0.2)
        ax.set(xlabel='', xlim=(0, 1), xticks=ticks)
        ax.set_yticklabels([strip_label(i.name) for i in metrics], fontsize=20)
        ax.set_xticklabels([f'{i:.0%}' for i in ticks], fontsize=20)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(axis='both', which='both', direction='inout')

        ax2 = ax.secondary_xaxis('top')
        ax2.xaxis.set_major_formatter(plt.NullFormatter())
        ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax2.tick_params(axis='x', which='both', direction='in')
        fig.subplots_adjust(left=0.2, bottom=0.25)
    else:
        fig, ax = s.plot_uncertainties(model, x_axis=metrics, kind='box',
                                       center_kws={'color': color, 'width': 0.5,
                                                   'whis': whis, 'sym': sym})
        if adjust_hue:
            ax.get_legend().remove()
        fig.set_figwidth(2.5)
        ticks = np.arange(0, 1.2, 0.2)
        ax.set_xticklabels([strip_label(i.name) for i in metrics], fontsize=20)
        for label in ax.get_xticklabels():
            label.set_rotation(90)
        ax.set(ylabel='', ylim=(0, 1), yticks=ticks)
        ax.set_yticklabels([f'{i:.0%}' for i in ticks], fontsize=20)
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(axis='both', which='both', direction='inout')

        ax2 = ax.secondary_yaxis('right')
        ax2.yaxis.set_major_formatter(plt.NullFormatter())
        ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax2.tick_params(axis='y', which='both', direction='in')
        fig.subplots_adjust(left=0.35, bottom=0.3)

    fig.savefig(os.path.join(figures_path, f'sys{ID}_recoveries.png'), dpi=300)

    return fig, ax


# Plot net cost and emissions as 2D-density plots
def plot_kde(model, ID, color, metrics):
    fig, ax = s.plot_uncertainties(model,
                                   x_axis=metrics[0],
                                   y_axis=metrics[1],
                                   kind='kde-box',
                                   center_kws={'fill': True, 'color': color},
                                   margin_kws={'color': color,
                                               'whis': (5, 95),
                                               'width': 0.4,
                                               'sym': ''})

    ax0, ax1, ax2 = fig.axes # KDE, top box, right box
    plot_minor_ticks_and_grid(ax0, nx=2, ny=2)
    ax0x = ax0.secondary_xaxis('top')
    ax0x.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax0y = ax0.secondary_yaxis('right')
    ax0y.yaxis.set_minor_locator(AutoMinorLocator(2))
    for ax in (ax0, ax0x, ax0y):
        ax.tick_params(axis='both', which='both', direction='inout', width=0.5)

    for i in (ax0, ax1):
        i.set_xlim(0, 60)
        i.set_xticks((np.arange(0, 70, 10)))
    for i in (ax0, ax2):
        i.set_ylim(-20, 140)
        i.set_yticks((np.arange(-20, 160, 20)))
    for txt in ax0.get_xticklabels()+ax0.get_yticklabels():
        txt.set_fontsize(20)
    ax0x.set_xticklabels('')
    ax0y.set_yticklabels('')

    ax1.xaxis.set_visible(False)
    ax1.spines.clear()
    ax2.spines.clear()
    ax2.yaxis.set_visible(False)
    fig.subplots_adjust(left=0.15)
    # ax[0].yaxis.get_major_ticks()[-1].gridline.set_visible(False) # if want to remove a certain gridline
    fig.savefig(os.path.join(figures_path, f'sys{ID}_cost_emission.png'), dpi=300)

    return fig, ax


param_set = set(modelA.parameters).union(set(modelB.parameters), set(modelC.parameters))
param_names = tuple(p.name for p in param_set)

# Plot Morris results as scatter plot for each metric of each model
def plot_morris_scatter(model, morris_dct, combined_morris, color,
                        label_lines=True, label_points=True,
                        plot_combined=False, axs=None):
    ID = model.system.ID[-1]
    normalized = {}
    ax_dct = {}
    for n, metric in enumerate(model.metrics):
        df = morris_dct[metric.name].copy()
        df.index = df.parameter

        df.mu_star_conf /= df.mu_star.max()
        df.sigma /= df.mu_star.max()
        # df.sigma /= df.sigma.max() # if want to normalize by sigma
        df.mu_star /= df.mu_star.max()
        normalized[metric.name] = df

        axs_n = None if axs is None else axs[n]
        fig, ax = s.plot_morris_results(normalized, metric,
                                        ax=axs_n,
                                        label_kind=None,
                                        k1=0.1, k2=None, k3=1,
                                        color=color)
        if not plot_combined:
            fig.set(figheight=3, figwidth=3)
            xlabel = r'Normalized $\mu^*$ [$\mu^*$/$\mu^*_{max}$]'
            ylabel = r'Normalized $\sigma$ [$\sigma$/$\mu^*_{max}$]'
        else:
            label_lines = False
            xlabel = ylabel = ''

        if label_lines:
            legend = ax.get_legend()
            if legend:
                legend.texts[0].set_text('$\\sigma/\\mu^*$=1')
                legend.texts[1].set_text('$\\sigma/\\mu^*$=0.1')
        else:
            legend = ax.get_legend()
            if legend:
                legend.remove()

        ticks = np.arange(0, 1.4, 0.2)
        ticklabels = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2]
        ax.set(xlim=(0, 1.2), ylim=(0, 1.2), xbound=(0, 1.2), ybound=(0, 1.2),
               xticks=ticks, yticks=ticks, xticklabels=ticklabels, yticklabels=ticklabels,
               xlabel=xlabel, ylabel=ylabel)

        # Only label the ones with mu_star/mu_star_max>=0.1
        top = df[(df.mu_star>=0.1)]
        # # If want mu_star/mu_star_max or sigma/mu_star_max >=0.1
        # top = df[(df.mu_star>=0.1)|(df.sigma>=0.1)]

        # Only mark out the top five mu_star/mu_star_max parameters
        if top.shape[0] > 5:
            top.sort_values(by=['mu_star'], ascending=False, inplace=True)
            top = top[:5]

        if plot_combined: # `adjust_text` won't work here
            if label_points:
                for idx in top.index:
                    label = combined_morris[combined_morris.parameter==idx].index.values.item()
                    ax.annotate(label, (top.loc[idx].mu_star, top.loc[idx].sigma),
                                xytext=(2, 2), textcoords='offset points',
                                ha='center', fontsize=14)
            labelsize = 'large'
        else:
            if label_points:
                labels = []
                for idx in top.index:
                    label = combined_morris[combined_morris.parameter==idx].index.values.item()
                    labels.append(ax.text(top.loc[idx].mu_star, top.loc[idx].sigma, label))
                from adjustText import adjust_text
                adjust_text(labels)
            labelsize = None

        plot_minor_ticks_and_grid(ax, nx=1, ny=1)
        ax.tick_params(axis='both', which='both', direction='inout', width=0.5,
                       labelsize=labelsize)

        if not plot_combined:
            fig.subplots_adjust(bottom=0.2, left=0.25)
            path = os.path.join(figures_path, f'sys{ID}_morris{n+1}_{metric.name}.png')
            fig.savefig(path, dpi=300)
            ax_dct[metric.name] = ax

    return ax_dct


def plot_morris_blank(ax_dct):
    temp = ax_dct['morris_scatter']['A']['Cost']
    xrng = np.random.default_rng(seed-1)
    yrng = np.random.default_rng(seed+1)
    plt.scatter(xrng.random(5)*0.1, yrng.random(5)*0.1, color='k')
    plt.scatter(xrng.random(5)*0.5+0.5, yrng.random(5)*0.5+0.5, color='k')
    plt.scatter(xrng.random(5)*0.2+0.8, yrng.random(5)*0.2, color='k')
    xline = np.array([0, 2])
    plt.plot(xline, xline, color='k', linestyle='-')
    plt.plot(xline, 0.1*xline, color='k', linestyle='-.')
    ax, fig = plt.gca(), plt.gcf()
    def copy_prop(prop, origin, new):
        f = getattr(new, f'set_{prop}')
        f(getattr(origin, f'get_{prop}')())
    for suffix in ('lim', 'bound', 'ticks',
                   # 'ticklabels', 'label'
                    ):
        for xy in ('x', 'y'):
            copy_prop(f'{xy}{suffix}', temp, ax)
    copy_prop('size_inches', temp.figure, fig)
    par_dct = {k: getattr(temp.figure.subplotpars, k)
               for k in ('left', 'right', 'top', 'bottom')}
    fig.subplots_adjust(**par_dct)
    ax.set(xticks=(), yticks=(), xticklabels=(), yticklabels=())
    ax_dct['morris_blank'] = ax
    fig.savefig(os.path.join(figures_path, 'morris_blank.png'), dpi=300)
    return ax_dct


def plot_morris_scatter_all(models, morris_dct, combineds, RGBs):
    fig = plt.figure(figsize=(12, 18))
    tot_models = len(models)
    tot_metrics = len(models[0].metrics)
    fig.subplots(tot_metrics, tot_models, sharex=True, sharey=True)

    label_colors = (# A/B/C
        mcolors.to_rgba('#7F533E'),
        mcolors.to_rgba('#5A7F60'),
        mcolors.to_rgba('#497980'))
    for n_model, model in enumerate(models):
        ID = model.system.ID[-1]
        axs = [fig.axes[tot_models*i+n_model] for i in range(tot_metrics)]
        plot_morris_scatter(
            model, morris_dct[ID], combineds['mu_star'], RGBs[ID],
            label_lines=False, label_points=True,
            plot_combined=True, axs=axs)

        if n_model == 0:
            for n, ax in enumerate(axs):
                name = model.metrics[n].name
                if not name in ('Cost', 'GWP'):
                    ylabel = strip_label(name)+' recovery'
                else:
                    ylabel = name
                ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')
                ax.yaxis.labelpad = 15
        axs[0].set_xlabel(f'System {ID}', fontsize=14, fontweight='bold',
                          color=label_colors[n_model])
        axs[0].xaxis.set_label_position('top')
        axs[0].xaxis.labelpad = 15

    fig.tight_layout()
    fig.supxlabel(r'Normalized $\mu^*$ [$\mu^*$/$\mu^*_{max}$]', y=0.01, fontsize=20, fontweight='bold')
    fig.supylabel(r'Normalized $\sigma$ [$\sigma$/$\mu^*_{max}$]', x=0, fontsize=20, fontweight='bold')
    fig.subplots_adjust(left=0.12, bottom=0.05)
    fig.savefig(os.path.join(figures_path, 'morris_normalized_combined.png'), dpi=300)


def plot_morris_bubble(combineds):
    sns.set_theme(style='ticks')

    dct = {'mu_star': combineds['mu_star'].copy(),
           'sigma': combineds['sigma'].copy()}
    for k, df in dct.items():
        columns = [f'{i[0]}-{i[1]}' for i in df.columns]
        columns[0] = 'parameter'
        df.columns = columns
        df.index = df.parameter
        df = df.drop('parameter', axis=1)
        dct[k] = df.stack(dropna=False).to_frame().rename(columns={0: k})

    plot_df = pd.concat(dct.values(), axis=1).reset_index()
    plot_df.rename(columns={'level_1': 'metric'}, inplace=True)

    g = sns.relplot(
        data=plot_df,
        x='metric', y='parameter', hue='sigma', size='mu_star',
        # palette="vlag",
        hue_norm=(0, 1), # edgecolors='0.7',
        height=15, sizes=(50, 250), size_norm=(-.2, .8),
    )

    g.ax.grid(color='k', linestyle='--', linewidth=0.5)
    for spine in ('top', 'right'):
        g.ax.spines[spine].set_visible(True)

    g.ax.set_xlabel('Metric', fontsize=20, fontweight='bold')
    g.ax.set_ylabel('Parameter', fontsize=20, fontweight='bold')
    g.ax.tick_params(axis='both', which='both', direction='inout')

    for label in g.ax.get_xticklabels():
        label.set_rotation(60)
        label.set_fontsize(14)
    for label in g.ax.get_yticklabels():
        label.set_fontsize(14)
    for text in g.legend.get_texts():
        text.set_fontsize(14)

    g.fig.subplots_adjust(left=0.22, bottom=0.15)
    g.fig.savefig(os.path.join(figures_path, 'morris_normalized_bubble.png'), dpi=300)

    return g.ax


def run(N_uncertainty=5000, N_morris=50, from_record=True,
        label_morris_lines=True, label_morris_points=True):
    np.random.seed(seed)

    # A dict with all results for future rebuilding of the models
    global table_dct, ax_dct
    ax_dct = dict(box={}, kde={}, morris_scatter={}, morris_bubble=None)
    if not from_record:
        # This saves reports for the system, TEA, and LCA
        for sys in (bw.sysA, bw.sysB, bw.sysC):
            save_reports(sys)
        table_dct = dict(uncertainty={}, morris_table={}, morris_dct={})
    else:
        path = os.path.join(results_path, 'table_dct.pckl')
        table_dct = load_pickle(path)

    ########## Uncertainty analysis ##########
    for model in models:
        # Net cost, net GWP, and total COD/N/P/K recovery
        model.metrics = key_metrics = get_key_metrics(model, alt_names)
        ID = model.system.ID[-1]
        if not from_record:
            samples = model.sample(N=N_uncertainty, rule='L', seed=seed)
            model.load_samples(samples)
            if model is modelB:
                copy_samples(modelA, model)
            elif model is modelC:
                copy_samples(modelA, model)
                copy_samples(modelB, model, exclude=modelA.parameters)

            evaluate(model)

            table_dct['uncertainty'][ID] = model.table.copy()
            organized_dct = organize_uncertainty_results(
                model=model, spearman_results=model.spearman())
            save_uncertainty_results(model, organized_dct)

        else:
            model.table = table_dct['uncertainty'][ID]

        # Make the plots
        _, ax_dct['box'][ID] = plot_box(model, ID, RGBs[ID], key_metrics[:-2], 'horizontal')
        _, ax_dct['kde'][ID] = plot_kde(model, ID, RGBs[ID], key_metrics[-2:])

    ########## Morris One-at-A-Time ##########
    origin_dct = dict(mu_star={}, sigma={})
    norm_dct = dict(mu_star={}, sigma={})
    for model in models:
        ID = model.system.ID[-1]
        if not from_record:
            inputs = s.define_inputs(model)
            # Want to use a larger N (maybe 100)
            morris_samples = s.generate_samples(inputs, kind='Morris', N=N_morris, seed=seed)

            evaluate(model, morris_samples)

            # These are the unprocessed data
            morris_dct = s.morris_analysis(model, inputs, seed=seed,
                                           nan_policy='fill_mean',
                                           file=os.path.join(results_path, f'sys{ID}_morris.xlsx'))

            table_dct['morris_table'][ID] = model.table.copy()
            table_dct['morris_dct'][ID] = morris_dct.copy()

            origins = dict(mu_star=[], sigma=[])
            filtereds = dict(mu_star=[], sigma=[])

            for i in model.metrics:
                df = morris_dct[i.name]
                df.sort_values(by=['mu_star'], ascending=False, inplace=True)
                origins['mu_star'].append(df.mu_star)
                origins['sigma'].append(df.sigma)

                df_filtered = df.iloc[0:5] # select the top five
                filtereds['mu_star'].append(df_filtered.mu_star)
                filtereds['sigma'].append(df_filtered.sigma)

            mu_star_origin = pd.concat([df.parameter, *origins['mu_star']], axis=1)
            sigma_origin = pd.concat([df.parameter, *origins['sigma']], axis=1)
            mu_star_filtered = pd.concat(filtereds['mu_star'], axis=1)
            sigma_filtered = pd.concat(filtereds['sigma'], axis=1)

            # Won't be able to divde if having different column names
            sigma_filtered.columns = mu_star_filtered.columns = [i.name for i in model.metrics]

            # Normalize
            sigma_norm = sigma_filtered/mu_star_filtered.max()
            mu_star_norm = mu_star_filtered/mu_star_filtered.max()

            for i in (mu_star_filtered, mu_star_norm, sigma_filtered, sigma_norm):
                i.insert(0, 'parameter', mu_star_origin.parameter)

            columns = ['parameter'] + [i.name for i in model.metrics]
            for df in (mu_star_origin, mu_star_filtered, mu_star_norm,
                       sigma_origin, sigma_filtered, sigma_origin):
                df.columns = columns

            origin_dct['mu_star'][f'{model.system.ID}'] = mu_star_origin
            origin_dct['sigma'][f'{model.system.ID}'] = sigma_origin
            norm_dct['mu_star'][f'{model.system.ID}'] = mu_star_norm
            norm_dct['sigma'][f'{model.system.ID}'] = sigma_norm
        else:
            model.table = table_dct['morris_table'][ID]
            morris_dct = table_dct['morris_dct'][ID]

    # Process data
    if not from_record:
        combineds = {}
        writer = pd.ExcelWriter(os.path.join(results_path, 'morris_normalized_combined.xlsx'))

        for k in ('mu_star', 'sigma'):
            columns = []
            data = []
            for i in key_metrics:
                for ID, df in norm_dct[k].items():
                    df.index = df.parameter
                    columns.append(f'{i.name}-{ID}')
                    data.append(df[i.name])

            combined = pd.concat(data, axis=1)
            columns_mi = pd.MultiIndex.from_tuples([i.split('-') for i in columns])
            combined.columns = columns_mi
            combined.reset_index(inplace=True)
            combineds[k] = combined

            combined.to_excel(writer, sheet_name=f'{k}_combined')
            for ID, df in origin_dct[k].items():
                df.to_excel(writer, sheet_name=f'{k}_all_{ID}')

        writer.save()

        table_dct['morris_combined'] = combineds
        pickle_path = os.path.join(results_path, 'table_dct.pckl')
        save_pickle(table_dct, pickle_path)
    else:
        combineds = table_dct['morris_combined']

    # Make Morris plots
    for model in models:
        ID = model.system.ID[-1]
        morris_dct = table_dct['morris_dct'][ID]
        ax_dct['morris_scatter'][ID] = plot_morris_scatter(
            model, morris_dct, combineds['mu_star'], RGBs[ID],
            label_morris_lines, label_morris_points)
    plot_morris_blank(ax_dct) # plot a blank one for annotating in the manuscript
    plot_morris_scatter_all(models, table_dct['morris_dct'], combineds, RGBs)
    ax_dct['morris_bubble'] = plot_morris_bubble(combineds)


# %%

# =============================================================================
# Acutally run the functions
# =============================================================================

if __name__ == '__main__':
    # Run from scratch
    # run(N_uncertainty=5000, N_morris=50, from_record=False,
    #     label_morris_lines=False, label_morris_points=True)

    # Just make the plots
    run(from_record=True, label_morris_lines=False, label_morris_points=True)