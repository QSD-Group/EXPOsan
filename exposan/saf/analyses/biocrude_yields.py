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

# !!! Temporarily ignoring warnings
import warnings
warnings.filterwarnings('ignore')

import os, numpy as np, pandas as pd, qsdsan as qs
from exposan.saf import (
    data_path,
    results_path,
    create_system,
    get_MFSP,
    )

data_path = os.path.join(data_path, 'biocrude_yields.csv')
df = pd.read_csv(data_path)

def MFSP_across_biocrude_yields(yields=[], **config_kwargs):
    sys = create_system(**config_kwargs)
    unit = sys.flowsheet.unit
    stream = sys.flowsheet.stream
    
    HTL = unit.HTL
    default_yield = HTL.dw_yields.copy()
    
    CrudeSplitter = unit.CrudeSplitter
    default_fracs = CrudeSplitter.cutoff_fracs.copy()
    crude_and_char0 = sum(default_fracs[1:])
    
    gas0, aq0 = default_yield['gas'], default_yield['aqueous']
    char0 = default_yield['biocrude'] * default_fracs[-1]/crude_and_char0
    non_crudes = [gas0, aq0, char0]
    non_crude0 = sum(non_crudes)
    non_crudes = [i/non_crude0 for i in non_crudes]
    
    def adjust_yield(y_crude):
        non_crude = 1 - y_crude
        return [i*non_crude for i in non_crudes]

    feedstock = stream.feedstock
    mixed_fuel = stream.mixed_fuel
    dry_feedstock = feedstock.F_mass - feedstock.imass['Water']
    
    MFSPs = []
    fuel_yields = []
    for y in yields:
        print(f'yield: {y}')
        sys.reset_cache()
        crude = y/100
        gas, aq, char = adjust_yield(crude)

        dw_yields = default_yield.copy()
        dw_yields['gas'] = gas
        dw_yields['aqueous'] = aq
        dw_yields['biocrude'] = crude+char
        HTL.dw_yields = dw_yields

        CrudeSplitter.cutoff_fracs = [
            1-crude_and_char0,
            crude_and_char0*crude/(crude+char),
            crude_and_char0*char/(crude+char),
            ]

        try: 
            sys.simulate()
            MFSP = get_MFSP(sys, print_msg=False)
            fuel_yield = mixed_fuel.F_mass/dry_feedstock
            print(f'MFSP: ${MFSP:.2f}/GGE; fuel yields {fuel_yield:.2%}.\n')
        except:
            MFSP = fuel_yield = None
            print('Simulation failed.\n')

        MFSPs.append(MFSP)
        fuel_yields.append(fuel_yield)

    return MFSPs, fuel_yields


if __name__ == '__main__':
    # config = {'include_PSA': False, 'include_EC': False,}
    config = {'include_PSA': True, 'include_EC': False,}
    # config = {'include_PSA': True, 'include_EC': True,}
    flowsheet = qs.main_flowsheet
    dct = globals()
    dct.update(flowsheet.to_dict())
    
    # single = [67.3] # normalized from the 80.2 biocrude+char
    # single = [20]
    # results = MFSP_across_biocrude_yields(yields=single, **config)
    
    yields_results = df.copy()
    tested = MFSP_across_biocrude_yields(yields=df.y_test, **config)
    yields_results['y_test_MFSP'] = tested[0]
    yields_results['y_test_yields'] = tested[1]

    predicted = MFSP_across_biocrude_yields(yields=df.y_pred, **config)
    yields_results['y_pred_MFSP'] = predicted[0]
    yields_results['y_pred_yields'] = predicted[1]
    
    outputs_path = os.path.join(results_path, f'biocrude_yields_{flowsheet}.csv')
    yields_results.to_csv(outputs_path)
