# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from qsdsan.sanunits import AnaerobicCSTR
from biosteam._tea import add_replacement_cost_to_cashflow_array

__all__ = ('categorize_construction_cashflow',
           'categorize_replacement_cashflows',
           'categorize_cashflow')

#%%

def categorize_replacement_cashflows(vessel, beads, dm, others, units, years, start):
    for u in units:
        for name, installed_cost in u.installed_costs.items():
            lifetime = u.equipment_lifetime.get(name)
            if lifetime:
                if '_beads' in name: arr = beads
                elif ' - ' in name: arr = others
                else: arr = dm if u.ID.startswith('DM') else vessel
                add_replacement_cost_to_cashflow_array(installed_cost, lifetime, arr, years, start)

def categorize_construction_cashflow(vessel, beads, dm, others, units, start, construction_schedule):
    for u in units:
        for name, installed_cost in u.installed_costs.items():
            if '_beads' in name: arr = beads
            elif ' - ' in name: arr = others
            else: arr = dm if u.ID.startswith('DM') else vessel
            arr[:start] += installed_cost * construction_schedule

def categorize_cashflow(tea):
    C_D, C_FC, C_WC, D, L, LI, LP, LPl, C, S, T, I, NE, CF, DF, NPV, CNPV = tea.get_cashflow_table().to_numpy().transpose() * 1e6
    DF /= 1e6
    start = tea._start
    years = tea._years
    length = start + years
    vessel, beads, dm, others, electricity, heat_onsite, chemicals, biogas_offset = data = np.zeros((8, length))
    system = tea.system
    
    # C_FC categorized
    categorize_construction_cashflow(vessel, beads, dm, others, system.cost_units, start, tea._construction_schedule)
    categorize_replacement_cashflows(vessel, beads, dm, others, system.cost_units, years, start)
    
    # sales
    biogas_offset[:] = -S
    
    # VOC
    electricity[start:] = sum(u.power_utility.cost for u in system.cost_units) * system.operating_hours
    heat_onsite[start:] = sum(sum(hu.cost for hu in u.heat_utilities) for u in system.cost_units) * system.operating_hours
    
    # FOC
    chemicals[start:] = tea.FOC 
    
    out = np.dot(data, DF)
    lcc = sum(out)
    out = dict(zip(['vessel', 'beads', 'dm', 'others', 'electricity', 
                    'heat_onsite', 'chemicals', 'biogas_offset'], out))
    out['total'] = lcc
    out['cnpv'] = CNPV[-1]
    return out
