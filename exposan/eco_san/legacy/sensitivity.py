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
from exposan import eco_san as es

m = es.models
modelC= m.modelC

a = es.analyses
key_metrics = a.key_metrics


#spearman_rho, fig, ax, all_params = a.run_plot_spearman(modelA, N=10000)

uncertainty = m.run_uncertainty(modelC, N=10000)
m.save_uncertainty_results(modelC)