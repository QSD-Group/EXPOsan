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
import math
from exposan import biogenic_refinery as br
from exposan.biogenic_refinery import (
    create_model,
    run_uncertainty,
    results_path,
    )

__all__ = ('get_recalcitrance_pontential',)


# %%

# =============================================================================
# Create module for biochar analysis
# =============================================================================


def get_default_uniform(b, ratio, lb=None, ub=None): # lb/ub for upper/lower bounds
    lower = max(b*(1-ratio), lb) if lb else b*(1-ratio)
    upper = min(b*(1+ratio), ub) if ub else b*(1+ratio)
    return shape.Uniform(lower=lower, upper=upper)

def get_recalcitrance_pontential(ID, model=None):
    model = model.copy()
    sys = model.system
    pyrolysis_unit = sys.path[7] #carbonizer base in biogenic refinery system
    
    def set_pyrolysis_temp(i):
        pyrolysis_unit.pyrolysis_temp = i
    
    f_AC_dec = pyrolysis_unit.f_ash_content/10 #converts % ash content of feedstock to decimal
    
    # predictive equation for biochar yield based on feedstock ash content and pyrolysis temperature
    dry_basis_yield = 1.18 * f_AC_dec ** 0.843 + (1 - f_AC_dec) * 2.106 * math.exp(-0.0066 * pyrolysis_unit.pyrolysis_temp)
    
    # predictive equation for ash-free biochar yield (Neves et al. 2011)
    ash_free_yield = 100 * (0.106 + 2.43 * math.exp(-0.0066 * pyrolysis_unit.pyrolysis_temp))
    
    # predictive equation for biochar fixed carbon content 
    b_fixed_carbon = 87.786 * ash_free_yield ** -0.483
    
    # calculate biochar volatile matter and ash content 
    b_ash_content = (dry_basis_yield - ash_free_yield) * 100 / dry_basis_yield
    b_volatile_matter = 100 - b_ash_content - b_fixed_carbon
    
    # predictive equation for carbon recalcitrance potential (Klasson 2017)
    recalcitrance_potential = 0.17 * (0.474 * b_volatile_matter + 0.963 * b_fixed_carbon + 0.067 * b_ash_content) / (100 - b_ash_content) + 0.00479
    return recalcitrance_potential
    
        