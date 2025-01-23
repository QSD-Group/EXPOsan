#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Created by Yuyao Huang, Siqi Tang, and Aaron Marszewski for uncertainty and sensitivity analysis of the EL system
'''

from exposan import enviroloo as el
from exposan.enviroloo import create_model, run_uncertainty

def run(model_ID, seed = None, N = 1000, country_specific =False, **model_kwargs):
          model = create_model(model_ID, country_specific = country_specific, **model_kwargs);
          run_uncertainty(model, seed = seed, N = N);

if _name_ == '_main_':
          el.INCLUDED_RESOURCE_RECOVERY = False;# does not include resource recovery at the current stage
          run('sysEL', seed = 5, N = 50);