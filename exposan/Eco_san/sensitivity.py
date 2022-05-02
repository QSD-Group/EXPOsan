#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  8 12:03:05 2021

@author: lewisrowles
"""

# =============================================================================
# Sample codes to run the analyses
# (copy and put in another script, otherwise these scripts will run everytime
# you load the module)
# =============================================================================

from qsdsan import stats as s
from exposan import Eco_san as es

m = es.models
modelC= m.modelC

a = es.analyses
key_metrics = a.key_metrics


#spearman_rho, fig, ax, all_params = a.run_plot_spearman(modelA, N=10000)

uncertainty = m.run_uncertainty(modelC, N=10000)
m.save_uncertainty_results(modelC)

# morris_dct, fig, ax = a.run_plot_morris(modelB, 10, test_convergence=False)

# morris_dct_conv, fig, ax = a.run_plot_morris(modelB, 100, test_convergence=True)

# fast_dct, fig, ax = a.run_plot_fast(modelB, 'FAST', 100, M=4)

# rbd_dct, fig, ax = a.run_plot_fast(modelB, 'RBD', 100, M=10)

# sobol_dct, fig, ax = a.run_plot_sobol(modelB, 10, file_prefix='')


# fig, ax = s.plot_uncertainties(modelB, metrics=key_metrics)

# fig, ax = s.plot_correlations(spearman_rho, parameters=modelB.get_parameters(),
#                               metrics=key_metrics[0])

# fig, ax = s.plot_correlations(spearman_rho, parameters=modelB.get_parameters(),
#                               metrics=key_metrics)

# fig, ax = s.plot_morris_results(morris_dct, key_metrics[0], label_kind='name')

# fig, ax = s.plot_morris_convergence(morris_dct_conv,
#                                     parameters=modelB.get_parameters(),
#                                     metric=key_metrics[0], plot_rank=True)

# fig, ax = s.plot_fast_results(fast_dct, key_metrics[0])

# fig, ax = s.plot_fast_results(rbd_dct, key_metrics[0])

# fig, ax = s.plot_sobol_results(sobol_dct, metric=key_metrics[0], kind='STS2',
#                                 plot_in_diagonal='ST')