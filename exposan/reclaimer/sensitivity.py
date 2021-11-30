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
from exposan import reclaimer as R
from qsdsan import Model


m = R.models
modelA= m.modelA

#a = R.analyses
#key_metrics = a.key_metrics


#spearman_rho, fig, ax, all_params = a.run_plot_spearman(modelA, N=10000) 

uncertainty = m.run_uncertainty(modelA, N=10000)
m.save_uncertainty_results(modelA)

# # Filter out parameters that only meet a certain threshold
# def filter_parameters(model, df, threshold):
#     new_df = pd.concat((df[df>=threshold], df[df<=-threshold]))
#     filtered = new_df.dropna(how='all')
#     param_dct = {p.name_with_units:p for p in model.get_parameters()}
#     parameters = set(param_dct[i[1]] for i in filtered.index)
#     return list(parameters)

# # Only want parameters with Spearman's rho >= 0.4 or <= -0.4
# modelA.parameters = key_parameters = \
#     filter_parameters(modelA, spearman_rho, threshold=0.4)
