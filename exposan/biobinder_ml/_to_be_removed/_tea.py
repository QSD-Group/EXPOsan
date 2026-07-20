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
import numpy as np, pandas as pd

__all__ = ('create_tea',)

#!!! Need to see if we can follow all assumptions as in Jianan's paper
#PNNL 32371 contains land costs for HTL & Upgrading 

class CAPEXTableBuilder:
    __slots__ = ('index', 'data')
    
    def __init__(self):
        self.index = []
        self.data =[]
    
    def entry(self, index: str, cost: list, notes: str = '-'):
        self.index.append(index)
        self.data.append([notes, *cost])

    @property
    def total_costs(self):
        N = len(self.data[0])
        return [sum([i[index] for i in self.data]) for index in range(1, N)]
    
    def table(self, names):
        return pd.DataFrame(self.data, 
                            index=self.index,
                            columns=('Notes', *[i + ' [MM$]' for i in names])
        )
class CostTableBuilder:
    __slots__ = ('index', 'data')
    
    def __init__(self):
        self.index = []
        self.data = []
    
    def entry(self, material_cost: list, utility_cost: list, notes: str = '-'):
        # Make sure to store both costs properly
        self.index.append(notes)
        self.data.append([*material_cost, *utility_cost])  # Flattening the costs into a single list

    def table(self, names):
        # Calculate number of utility costs based on data
        num_material_costs = len(self.data[0]) // 2  # Assuming equal split between material and utility costs
        num_utility_costs = len(self.data[0]) - num_material_costs
        
        columns = ['Notes'] + names + [f'Utility Cost {i + 1}' for i in range(num_utility_costs)]
        
        return pd.DataFrame(self.data, index=self.index, columns=columns)


class TEA(HTL_TEA):
    '''
    With only minor modifications to the TEA class in the HTL module.
    '''
    # __slots__ = (*HTL_TEA.__slots__, 'land')
    
    land = 0.
    
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule,
                 startup_months, startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac, WC_over_FCI,  finance_interest,
                 finance_years, finance_fraction, OSBL_units, warehouse,
                 site_development, additional_piping, proratable_costs,
                 field_expenses, construction, contingency,
                 other_indirect_costs, labor_cost, labor_burden,
                 property_insurance, maintenance, steam_power_depreciation,
                 boiler_turbogenerator, **kwargs):
        HTL_TEA.__init__(self, system, IRR, duration, depreciation, income_tax,
                     operating_days, lang_factor, construction_schedule,
                     startup_months, startup_FOCfrac, startup_VOCfrac,
                     startup_salesfrac, WC_over_FCI,  finance_interest,
                     finance_years, finance_fraction, OSBL_units, warehouse,
                     site_development, additional_piping, proratable_costs,
                     field_expenses, construction, contingency,
                     other_indirect_costs, labor_cost, labor_burden,
                     property_insurance, maintenance, steam_power_depreciation,
                     boiler_turbogenerator)
        for attr, val in kwargs.items():
            setattr(self, attr, val)
    
    @property
    def labor_cost(self):
        if callable(self._labor_cost): return self._labor_cost()
        return self._labor_cost
    @labor_cost.setter
    def labor_cost(self, i):
        self._labor_cost = i
    
def create_tea(sys, **kwargs):
    OSBL_units = bst.get_OSBL(sys.cost_units)
    try:
        BT = tmo.utils.get_instance(OSBL_units, (bst.BoilerTurbogenerator, bst.Boiler))
    except:
        BT = None
    default_kwargs = {
        'IRR': 0.1,
        'duration': (2020, 2050),
        'depreciation': 'MACRS7', # Jones et al. 2014
        'income_tax': 0.21, # Davis et al. 2018
        'operating_days': sys.operating_hours/24, # Jones et al. 2014
        'lang_factor': None, # related to expansion, not needed here
        'construction_schedule': (0.08, 0.60, 0.32), # Jones et al. 2014
        'startup_months': 6, # Jones et al. 2014
        'startup_FOCfrac': 1, # Davis et al. 2018
        'startup_salesfrac': 0.5, # Davis et al. 2018
        'startup_VOCfrac': 0.75, # Davis et al. 2018
        'WC_over_FCI': 0.05, # Jones et al. 2014
        'finance_interest': 0.08, # use 3% for waste management, use 8% for biofuel
        'finance_years': 10, # Jones et al. 2014
        'finance_fraction': 0.6, # debt: Jones et al. 2014
        'OSBL_units': OSBL_units,
        'warehouse': 0.04, # Knorr et al. 2013
        'site_development': 0.10, # Snowden-Swan et al. 2022
        'additional_piping': 0.045, # Knorr et al. 2013
        'proratable_costs': 0.10, # Knorr et al. 2013
        'field_expenses': 0.10, # Knorr et al. 2013
        'construction': 0.20, # Knorr et al. 2013
        'contingency': 0.10, # Knorr et al. 2013
        'other_indirect_costs': 0.10, # Knorr et al. 2013
        'labor_cost': 1e6, # use default value
        'labor_burden': 0.90, # Jones et al. 2014 & Davis et al. 2018
        'property_insurance': 0.007, # Jones et al. 2014 & Knorr et al. 2013
        'maintenance': 0.03, # Jones et al. 2014 & Knorr et al. 2013
        'steam_power_depreciation':'MACRS20',
        'boiler_turbogenerator': BT,
        'land':0
        }
    default_kwargs.update(kwargs)
    
    tea = TEA(
        system=sys, **default_kwargs)
        # IRR=IRR_value,
        # duration=(2020, 2050),
        # depreciation='MACRS7', # Jones et al. 2014
        # income_tax=income_tax_value, # Davis et al. 2018
        # operating_days=sys.operating_hours/24, # Jones et al. 2014
        # lang_factor=None, # related to expansion, not needed here
        # construction_schedule=(0.08, 0.60, 0.32), # Jones et al. 2014
        # startup_months=6, # Jones et al. 2014
        # startup_FOCfrac=1, # Davis et al. 2018
        # startup_salesfrac=0.5, # Davis et al. 2018
        # startup_VOCfrac=0.75, # Davis et al. 2018
        # WC_over_FCI=0.05, # Jones et al. 2014
        # finance_interest=finance_interest_value, # use 3% for waste management, use 8% for biofuel
        # finance_years=10, # Jones et al. 2014
        # finance_fraction=0.6, # debt: Jones et al. 2014
        # OSBL_units=OSBL_units,
        # warehouse=0.04, # Knorr et al. 2013
        # site_development=0.10, # Snowden-Swan et al. 2022
        # additional_piping=0.045, # Knorr et al. 2013
        # proratable_costs=0.10, # Knorr et al. 2013
        # field_expenses=0.10, # Knorr et al. 2013
        # construction=0.20, # Knorr et al. 2013
        # contingency=0.10, # Knorr et al. 2013
        # other_indirect_costs=0.10, # Knorr et al. 2013
        # labor_cost=labor_cost, # use default value
        # labor_burden=0.90, # Jones et al. 2014 & Davis et al. 2018
        # property_insurance=0.007, # Jones et al. 2014 & Knorr et al. 2013
        # maintenance=0.03, # Jones et al. 2014 & Knorr et al. 2013
        # steam_power_depreciation='MACRS20',
        # boiler_turbogenerator=BT,
        # land=land)
    return tea

def capex_table(teas, names=None):
    if isinstance(teas, bst.TEA): teas = [teas]
    capex = CAPEXTableBuilder()
    tea, *_ = teas
    ISBL_installed_equipment_costs = np.array([i.ISBL_installed_equipment_cost / 1e6 for i in teas])
    OSBL_installed_equipment_costs = np.array([i.OSBL_installed_equipment_cost / 1e6 for i in teas])
    capex.entry('ISBL installed equipment cost', ISBL_installed_equipment_costs)
    capex.entry('OSBL installed equipment cost', OSBL_installed_equipment_costs)
    ISBL_factor_entry = lambda name, value: capex.entry(name, ISBL_installed_equipment_costs * value, f"{value:.1%} of ISBL")
    ISBL_factor_entry('Warehouse', tea.warehouse)
    ISBL_factor_entry('Site development', tea.site_development)
    ISBL_factor_entry('Additional piping', tea.additional_piping)
    TDC = np.array(capex.total_costs)
    capex.entry('Total direct cost (TDC)', TDC)
    TDC_factor_entry = lambda name, value: capex.entry(name, TDC * value, f"{value:.1%} of TDC")
    TDC_factor_entry('Proratable costs', tea.proratable_costs)
    TDC_factor_entry('Field expenses', tea.field_expenses)
    TDC_factor_entry('Construction', tea.construction)
    TDC_factor_entry('Contingency', tea.contingency)
    TDC_factor_entry('Other indirect costs (start-up, permits, etc.)', tea.other_indirect_costs)
    TIC = np.array(capex.total_costs) - 2 * TDC
    capex.entry('Total indirect cost', TIC)
    FCI = TDC + TIC
    capex.entry('Fixed capital investment (FCI)', FCI)
    working_capital = FCI * tea.WC_over_FCI
    capex.entry('Working capital', working_capital, f"{tea.WC_over_FCI:.1%} of FCI")
    TCI = FCI + working_capital
    capex.entry('Total capital investment (TCI)', TCI)
    if names is None: names = [i.system.ID for i in teas]
    names = [i for i in names]
    return capex.table(names)
voc_table = bst.report.voc_table

def foc_table(teas, names=None):
    if isinstance(teas, bst.TEA): teas = [teas]
    tea, *_ = teas
    foc = bst.report.FOCTableBuilder()
    ISBL = np.array([i.ISBL_installed_equipment_cost / 1e6 for i in teas])
    labor_cost = np.array([i.labor_cost / 1e6 for i in teas])
    foc.entry('Labor salary', labor_cost)
    foc.entry('Labor burden', tea.labor_burden * labor_cost, '90% of labor salary')
    foc.entry('Maintenance', tea.maintenance * ISBL, f'{tea.maintenance:.1%} of ISBL')
    foc.entry('Property insurance', tea.property_insurance * ISBL, f'{tea.property_insurance:.1%} of ISBL')
    if names is None: names = [i.system.ID for i in teas]
    names = [i + ' MM$/yr' for i in names]
    return foc.table(names)
def cost_table(teas, names=None):
    if isinstance(teas, bst.TEA):
        teas = [teas]
    
    cost_builder = CostTableBuilder()
    tea, *_ = teas  # Get the first TEA object for shared attributes

    material_costs = np.array([i.material_cost / 1e6 for i in teas])  # Convert to MM$
    utility_costs = np.array([i.utility_cost / 1e6 for i in teas])    # Convert to MM$

    note = tea.name if hasattr(tea, 'name') else 'Cost Summary'
    
    cost_builder.entry(material_costs, utility_costs, note)

    if names is None:
        names = [f'TEA {i+1}' for i in range(len(teas))]  # Generate names if not provided

    return cost_builder.table(names)






# # Example usage in the main block
# if __name__ == '__main__':
#     your_tea_instance = ...  # Replace with actual TEA instance(s)
#     cost_df = cost_table(your_tea_instance)  # Call the function
#     print(cost_df)  # Print the generated DataFrame
