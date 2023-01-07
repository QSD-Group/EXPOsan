#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>
    Lane To <lane20@illinois.edu>
    Lewis Rowles <stetsonsc@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

from chaospy import distributions as shape
from exposan import biogenic_refinery as br
from exposan.biogenic_refinery import (
    create_model,
    run_uncertainty,
    results_path,
    update_resource_recovery_settings,
    )

__all__ = ('create_pyrolysis_model', 'get_dry_basis_yield', 'get_ash_free_yield', 'get_fixed_carbon', 'get_recalcitrance_pontential',)


# %%

# =============================================================================
# Create copy of model for biochar analysis
# =============================================================================


def get_default_uniform(b, ratio, lb=None, ub=None): # lb/ub for upper/lower bounds
    lower = max(b*(1-ratio), lb) if lb else b*(1-ratio)
    upper = min(b*(1+ratio), ub) if ub else b*(1+ratio)
    return shape.Uniform(lower=lower, upper=upper)

def create_pyrolysis_model(ID, model=None, ):
    