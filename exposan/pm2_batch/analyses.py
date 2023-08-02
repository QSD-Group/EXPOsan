# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from time import time
from qsdsan.utils import ospath, load_data
from qsdsan.stats import get_correlations, plot_correlations
from exposan.pm2_batch import (
    results_path,
    figures_path,
    create_model,
    run_uncertainty
    )
from biosteam.evaluation._utils import var_indices
from math import ceil

import os, numpy as np, pandas as pd, matplotlib as mpl, seaborn as sns, \
    matplotlib.pyplot as plt
    # matplotlib.ticker as tk

mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams["figure.autolayout"] = True
mpl.rcParams['xtick.minor.visible'] = True

N = 5000

T = 0.25
# T = 7    # exclude
t_step = 0.01
# rmse_thresholds = [25, 25, 25]
nrmse_thresholds = [None, 0.1, 0.1]

kind='include'
# kind='exclude'

analysis='uasa'
# analysis='cali'

#%%
def seed_RGT():
    files = os.listdir(results_path)
    seeds = [int(file_name[-3:]) for file_name in files if file_name.startswith('time_series_data')]
    seed = int(str(time())[-3:])
    if len(set(seeds)) >= 1000:
        raise RuntimeError('The program has run out of 3-digit seeds to use. Consider'
                           'clean up the results folder.')
    while seed in seeds:
        seed = (seed+1) % 1000
    return seed

#%%
def run_UA_SA(seed=None, N=N, T=T, t_step=t_step, thresholds=[], kind=kind, analysis=analysis, plot=False):
    seed = seed or seed_RGT()
    mdl = create_model(kind=kind, analysis=analysis)
    mdl = run_uncertainty(mdl, N, T, t_step, seed=seed)
    thresholds = update_thresholds(mdl, thresholds)
    D, p = get_correlations(mdl, kind='KS', thresholds=thresholds,
                            file=ospath.join(results_path, f'KS_test_{seed}.xlsx')
                            )
    if plot:
        plot_cdf_by_group(mdl, thresholds=thresholds)
    fig, ax = plot_correlations(D, close_fig=False,
                                file=ospath.join(figures_path, 'KS_test_D.png'))
    return mdl

def update_thresholds(mdl, thresholds, metrics=None, quantile=0.25):
    metrics = metrics or mdl.metrics
    thresholds = thresholds or [None]*len(metrics)
    data = mdl.table[var_indices(metrics)]
    for i, col in enumerate(data):
        if thresholds[i] is None:
            thresholds[i] = data[col].quantile(quantile)
    return thresholds

def plot_cdf_by_group(mdl=None, seed=None, thresholds=None, parameters=None, metrics=None, kind=kind, analysis=analysis):
    if mdl is None:
        # global mdl
        mdl = create_model(kind=kind, analysis=analysis)
        mdl.table = load_data(ospath.join(results_path, f'table_{seed}.xlsx'), header=[0, 1])
    metrics = metrics or mdl.metrics
    parameters = parameters or mdl.parameters
    thresholds = thresholds or update_thresholds(mdl, thresholds, metrics)
    x_df = mdl.table[var_indices(parameters)]
    y_df = mdl.table[var_indices(metrics)]
    ncol = 4
    nrow = ceil(x_df.shape[1]/ncol)
    for m, t in zip(y_df.items(), thresholds):
        y, err = m
        group = err <= t
        fig, axes = plt.subplots(nrow, ncol, #sharey=True,
                                 figsize=(ncol*4, nrow*4),
                                 layout='constrained')
        for col, ax in zip(x_df, axes.ravel()):
            sns.kdeplot(data=x_df[col][group], ax=ax,
                        cumulative=False, common_norm=True,
                        # cumulative=True,
                        label=f'{y[-1]} <= {round(t,2)}')
            sns.kdeplot(data=x_df[col][1-group], ax=ax,
                        cumulative=False, common_norm=True,
                        # cumulative=True,
                        label=f'{y[-1]} > {round(t,2)}')
            ax.tick_params(axis='both', which='both', direction='inout')
            ax.legend()
            ax.set_xlabel(col[-1])
            ax.set_ylabel('density')
        # fig.subplots_adjust(hspace=0.4, wspace=0.05, bottom=0.2)
        fig.savefig(ospath.join(figures_path, f'pdf_{y[-1]}_v2.png'),
                    dpi=300, facecolor='white')
        del fig, axes

def KS_test_var_thresholds(mdl=None, seed=None, kind=kind, analysis=analysis):
    if mdl is None:
        mdl = create_model(kind=kind, analysis=analysis)
        mdl.table = load_data(ospath.join(results_path, f'table_{seed}.xlsx'),
                              header=[0,1])
    sig = []
    thresholds = []
    quantiles = np.linspace(0.05, 0.5, 10)
    for q in quantiles:
        thrs = update_thresholds(mdl, [], quantile=q)
        D, p = mdl.kolmogorov_smirnov_d(thresholds=thrs)
        thresholds.append(thrs)
        sig.append(p < 0.05)
    out = {m: pd.DataFrame() for m in var_indices(mdl.metrics)}
    thresholds = np.asarray(thresholds).T

    for s, q in zip(sig, quantiles):
        for m, df in out.items():
            df[q] = s[m]
    with pd.ExcelWriter(ospath.join(results_path, f'sig_params_{seed}.xlsx')) as writer:
        for m, thrs in zip(out, thresholds):
            df = out[m]
            df.index = df.index.droplevel()
            df.columns = pd.MultiIndex.from_tuples([(col, t) for col,t \
                                                    in zip(df.columns, thrs)],
                                                   names=['quantile', 'NRMSE'])
            df.to_excel(writer, sheet_name=m[-1])

#%%
if __name__ == '__main__':

    # seed = 100          # include, N=1000, exo_vars_batch_may_unit
    # seed = 101          # include, N=2000, exo_vars_batch_may_unit
    # seed = 102          # include, N=3000, exo_vars_batch_may_unit
    # seed = 103          # include, N=5000, exo_vars_batch_may_unit
    # seed = 107          # include, N=10000, exo_vars_batch_may_unit
    # seed = 105          # include, N=15000, exo_vars_batch_may_unit
    # seed = 106          # include, N=20000, exo_vars_batch_may_unit

    # seed = 200          # include, N=10000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 201          # include, N=1000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 202          # include, N=2000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 203          # include, N=3000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 205          # include, N=5000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 206          # include, N=15000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 210          # include, N=20000, exo_vars_batch_may_unit, modified_pm2_kwargs

    # seed = 300          # include, N=10000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 301          # include, N=1000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 302          # include, N=2000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 303          # include, N=3000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 305          # include, N=5000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 306          # include, N=15000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 310          # include, N=20000, exo_vars_batch_may_unit, modified_pm2_kwargs

    # seed = 400          # include, N=10000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 502          # include, N=2000, exo_vars_batch_may_unit, modified_pm2_kwargs
    # seed = 503          # include, N=3000, exo_vars_batch_may_unit, modified_pm2_kwargs
    seed = 505          # include, N=5000, exo_vars_batch_may_unit, modified_pm2_kwargs

    # seed = 160          # include, N=10000, exo_vars_batch_may_unit
    # seed = 150          # include, N=10000, exo_vars_may_kinetic (wrong)
    # seed = 125          # include, N=1000, exo_vars_may_kinetic (wrong)

    # seed = 170          # exclude, N=10000 (final)
    # seed = 130          # exclude, N=1000,
    # seed = 119          # DOE verification ver.
    # seed = 111          # exclude, N=10000

    mdl = run_UA_SA(seed=seed)

    # thrs = [0.343, 0.05, 0.08]
    plot_cdf_by_group(seed=seed)

    KS_test_var_thresholds(mdl=mdl,seed=seed)