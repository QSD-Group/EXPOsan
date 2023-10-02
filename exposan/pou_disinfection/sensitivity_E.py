# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 17:17:30 2022

@author: be05055
"""

from qsdsan import stats as a
from exposan import POU_dis as pou 

m = pou.models
modelE = m.modelE


uncertainty = m.run_uncertainty(modelE, seed=5, N=10)
m.save_uncertainty_results(modelE)
