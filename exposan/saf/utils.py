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

import numpy as np, pandas as pd

__all__ = (
    'find_Lr_Hr',
    'get_mass_energy_balance',
    )

# To find Lr/Hr of a distillation column
Lr_trial_range = Hr_trial_range = np.linspace(0.05, 0.95, 19)
def find_Lr_Hr(unit, target_light_frac=None, Lr_trial_range=Lr_trial_range, Hr_trial_range=Hr_trial_range):
    results = {}
    outs0, outs1 = unit.outs
    F_mass_in = unit.F_mass_in
    _Lr, _Hr = unit.Lr, unit.Hr
    for Lr in Lr_trial_range:
        unit.Lr = round(Lr,2)
        Hr_results = {}
        for Hr in Hr_trial_range:
            unit.Hr = round(Hr,2)
            try: 
                unit.simulate()
                Hr_results[Hr] = outs0.F_mass/F_mass_in
            except:
                Hr_results[Hr] = None
        results[Lr] = Hr_results
    results_df = pd.DataFrame.from_dict(results) # columns are Lr, rows are Hr
    unit.Lr, unit.Hr = _Lr, _Hr
    try: unit.simulate()
    except: pass
    if not target_light_frac:
        return results_df
    try:
        diff_df = (results_df-target_light_frac).abs()
        where = np.where(diff_df==diff_df.min(None))
        Lr = results_df.columns[where[1]].to_list()[0]
        Hr = results_df.index[where[0]].to_list()[0]
    except: Lr = Hr = None
    return results_df, Lr, Hr


def get_mass_energy_balance(sys):
    IDs = []
    F_mass_ins = []
    F_mass_outs = []
    mass_ratios = []
    Hnets = []
    duties = []
    energy_ratios = []
    for i in sys.path:
        IDs.append(i.ID)
        F_mass_ins.append(i.F_mass_in)
        F_mass_outs.append(i.F_mass_out)
        mass_ratios.append(i.F_mass_out/i.F_mass_in)
        Hnets.append(i.Hnet)
        duty = sum(hu.duty for hu in i.heat_utilities)
        duties.append(duty)
        if duty == 0:
            energy_ratios.append(0)
        else:
            energy_ratios.append(i.Hnet/duty)

    df = pd.DataFrame({
        'ID': IDs,
        'F_mass_in': F_mass_ins,
        'F_mass_out': F_mass_outs,
        'mass_ratio': mass_ratios,
        'Hnet': Hnets,
        'duty': duties,
        'energy_ratio': energy_ratios,
        })
    
    return df