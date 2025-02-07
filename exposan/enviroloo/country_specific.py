#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Created by Yuyao Huang and Siqi Tang for implication of Enviroloo Clear Reinvented Toilet for specific country of interest.
'''

from chaospy import distributions as shape
from qsdsan import ImpactItem, PowerUtility
from exposan.utils import general_country_specific_inputs, run_module_country_specific
from exposan import enviroloo as el
from exposan.enviroloo import (
    create_model,
    run_uncertainty,
    results_path,
    #update_resource_recovery_settings,
    )

__all__ = ('create_country_specific_model',)


