#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import thermosteam as tmo, biosteam as bst
from exposan.htl import HTL_TEA

__all__ = ('create_tea',)

#!!! Need to see if we can follow all assumptions as in Jianan's paper

class TEA(HTL_TEA):
    '''
    With only minor modifications to the TEA class in the HTL module.
    '''
    
    @property
    def labor_cost(self):
        if callable(self._labor_cost): return self._labor_cost()
        return self._labor_cost
    @labor_cost.setter
    def labor_cost(self, i):
        self._labor_cost = i
    
def create_tea(sys, IRR_value=0.1, income_tax_value=0.21, finance_interest_value=0.08, labor_cost=1e6):
    OSBL_units = bst.get_OSBL(sys.cost_units)
    try:
        BT = tmo.utils.get_instance(OSBL_units, (bst.BoilerTurbogenerator, bst.Boiler))
    except:
        BT = None
    tea = TEA(
        system=sys, 
        IRR=IRR_value,
        duration=(2020, 2050),
        depreciation='MACRS7', # Jones et al. 2014
        income_tax=income_tax_value, # Davis et al. 2018
        operating_days=sys.operating_hours/24, # Jones et al. 2014
        lang_factor=None, # related to expansion, not needed here
        construction_schedule=(0.08, 0.60, 0.32), # Jones et al. 2014
        startup_months=6, # Jones et al. 2014
        startup_FOCfrac=1, # Davis et al. 2018
        startup_salesfrac=0.5, # Davis et al. 2018
        startup_VOCfrac=0.75, # Davis et al. 2018
        WC_over_FCI=0.05, # Jones et al. 2014
        finance_interest=finance_interest_value, # use 3% for waste management, use 8% for biofuel
        finance_years=10, # Jones et al. 2014
        finance_fraction=0.6, # debt: Jones et al. 2014
        OSBL_units=OSBL_units,
        warehouse=0.04, # Knorr et al. 2013
        site_development=0.09, # Knorr et al. 2013
        additional_piping=0.045, # Knorr et al. 2013
        proratable_costs=0.10, # Knorr et al. 2013
        field_expenses=0.10, # Knorr et al. 2013
        construction=0.20, # Knorr et al. 2013
        contingency=0.10, # Knorr et al. 2013
        other_indirect_costs=0.10, # Knorr et al. 2013
        labor_cost=labor_cost, # use default value
        labor_burden=0.90, # Jones et al. 2014 & Davis et al. 2018
        property_insurance=0.007, # Jones et al. 2014 & Knorr et al. 2013
        maintenance=0.03, # Jones et al. 2014 & Knorr et al. 2013
        steam_power_depreciation='MACRS20',
        boiler_turbogenerator=BT)
    return tea