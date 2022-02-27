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
the `qsdsan.stats` module (https://qsdsan.readthedocs.io/en/latest/stats.html).

'''

# =============================================================================
# Setup and uncertainty
# =============================================================================

import pandas as pd
from qsdsan import stats as s
from exposan import bwaise as bw

m = bw.models
modelA = bw.modelA

# Total COD/N/P/K recovery and net cost/GWP
modelA.metrics = key_metrics = bw.get_key_metrics(
    modelA, alt_names={'Annual net cost': 'Cost',
                       'Net emission GlobalWarming': 'GWP'})

seed = 3221 # set numpy seed for sample reproducibility

# Run Monte Carlo uncertainty analysis and get Spearman rank correlations,
# here we use a small sample size for demonstrative purpose
m.run_uncertainty(modelA, N=100, seed=seed, rule='L',
                  percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

# Pass a path to `file` or use `fig.savefig` if want to save the figure
fig, ax = s.plot_uncertainties(modelA,
                               x_axis=key_metrics[:-2], # only recoveries
                               kind='box', file='')

# Trim figure
fig.subplots_adjust(bottom=0.25)
for label in ax.get_xticklabels():
    label.set_rotation(45)

# Kernel density curve can be added to the histogram,
# with a log scale, we can have all metric results in the same plot
fig, ax = s.plot_uncertainties(modelA, y_axis=key_metrics, kind='hist',
                               center_kws={'kde':True, 'log_scale': 10})

# We can also have 2D histogram plot
fig, axes = s.plot_uncertainties(modelA,
                                 x_axis=key_metrics[-2], # cost
                                 y_axis=key_metrics[-1], # GWP
                                 kind='hist-box')

# Similar to histogram plots, kernel density plots can be 1D
fig, ax = s.plot_uncertainties(modelA, x_axis=key_metrics, kind='kde',
                               center_kws={'fill': True, 'log_scale': 2})

fig.subplots_adjust(bottom=0.25)

# Or 2D with different kinds of margins
fig, axes = s.plot_uncertainties(modelA, x_axis=key_metrics[-2],
                                 y_axis=key_metrics[-1], kind='kde-kde',
                                 center_kws={'fill': True})

fig, axes = s.plot_uncertainties(modelA, x_axis=key_metrics[-2],
                                 y_axis=key_metrics[-1], kind='kde-hist',
                                 center_kws={'fill': True},
                                 margin_kws={'kde': True, 'fill': False})


# %%

# =============================================================================
# Spearman
# =============================================================================

spearman_rho, spearman_p = s.get_correlations(
    modelA, kind='Spearman', nan_policy='raise',
    file='') # pass a path to `file` if you want to save the results as an Excel

# Filter out parameters that only meet a certain threshold
def filter_parameters(model, df, threshold):
    new_df = pd.concat((df[df>=threshold], df[df<=-threshold]))
    filtered = new_df.dropna(how='all')
    param_dct = {p.name_with_units:p for p in model.get_parameters()}
    parameters = set(param_dct[i[1]] for i in filtered.index)
    return list(parameters)

# Only want parameters with Spearman's rho >= 0.4 or <= -0.4
modelA.parameters = key_parameters = \
    filter_parameters(modelA, spearman_rho, threshold=0.4)

fig, ax = s.plot_correlations(spearman_rho, parameters=key_parameters,
	                          metrics=key_metrics[-2])

fig.subplots_adjust(left=0.25)


fig, ax = s.plot_correlations(
    spearman_rho, parameters=key_parameters, metrics=key_metrics)


#%%

# =============================================================================
# Morris
# =============================================================================

# Run Morris analysis without testing the convergence,
# here we use a small sample size for demonstrative purpose
inputs = s.define_inputs(modelA)
morris_samples = s.generate_samples(inputs, kind='Morris', N=10, seed=seed)

evaluate = bw.evaluate
evaluate(modelA, morris_samples)

dct = s.morris_analysis(modelA, inputs, metrics=key_metrics, seed=seed,
                        nan_policy='fill_mean')

# Unfortunately the auto-labeling is not good when you have close points,
# so you'll have to do some manual manipulation
fig, ax = s.plot_morris_results(dct, key_metrics[-2])

fig.subplots_adjust(bottom=0.3)


# Test if mu_star can converge within 10 trajectories
# (spoiler: it cannot because we already sort of selected the key parameters,
# and you will get a message prompt)
dct = s.morris_till_convergence(modelA, inputs, metrics=key_metrics, seed=seed,
                                N_max=100)

# Look at mu_star values for two parameters with regard to cost
fig, ax = s.plot_morris_convergence(dct,
                                    parameters=key_parameters[:2],
                                    metric=key_metrics[-2], plot_rank=False)


# Look at ranks of mu_star values for all parameters with regard to cost
fig, ax = s.plot_morris_convergence(dct, parameters=key_parameters,
                                    metric=key_metrics[-2], plot_rank=True)


# %%

# =============================================================================
# FAST
# =============================================================================

# Total and main effects from FAST analysis,
# here we use a small sample size for demonstrative purpose
fast_samples = s.generate_samples(inputs, kind='FAST', N=100, M=4, seed=seed)

evaluate(modelA, fast_samples)

dct = s.fast_analysis(modelA, inputs, kind='FAST', metrics=key_metrics,
                      M=4, seed=seed, nan_policy='fill_mean')

fig, ax = s.plot_fast_results(dct, metric=key_metrics[-2])

fig.subplots_adjust(left=0.4)


# Main effects from RBD-FAST analysis,
# here we use a small sample size for demonstrative purpose
fast_samples = s.generate_samples(inputs, kind='RBD', N=100, seed=seed)

evaluate(modelA, fast_samples)

dct = s.fast_analysis(modelA, inputs, kind='RBD', metrics=key_metrics,
                      seed=seed, nan_policy='fill_mean')

fig, ax = s.plot_fast_results(dct, metric=key_metrics[-2])

fig.subplots_adjust(left=0.4)


# %%

# =============================================================================
# Sobol
# =============================================================================

# Run Sobol analysis, here we use a small sample size for demonstrative purpose
sobol_samples = s.generate_samples(inputs, kind='Sobol', N=10,
                                   calc_second_order=True)

evaluate(modelA, sobol_samples)

dct = s.sobol_analysis(modelA, inputs, metrics=key_metrics, seed=seed,
                       calc_second_order=True, conf_level=0.95,
                       nan_policy='fill_mean')

fig, ax = s.plot_sobol_results(dct, metric=key_metrics[-1], kind='STS1')

fig.subplots_adjust(left=0.4, top=0.95)


fig, ax = s.plot_sobol_results(dct, metric=key_metrics[-1], kind='STS2',
                               plot_in_diagonal='ST')

for label in ax.get_xticklabels():
    label.set_rotation(45)

fig.subplots_adjust(left=0.4, bottom=0.4)


fig, ax = s.plot_sobol_results(dct, metric=key_metrics[-1], kind='all')