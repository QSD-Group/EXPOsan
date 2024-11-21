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
from qsdsan.utils import time_printer
from exposan.saf import (
    create_system,
    data_path,
    HTL_yields,
    get_GWP,
    get_MFSP,
    results_path,
    )

data_path = os.path.join(data_path, 'biocrude_yields.csv')
df = pd.read_csv(data_path)

@time_printer
def evaluation_across_biocrude_yields(yields=[], **config_kwargs):
    sys = create_system(**config_kwargs)
    unit = sys.flowsheet.unit
    stream = sys.flowsheet.stream
    
    HTL = unit.HTL
    CrudeSplitter = unit.CrudeSplitter
    default_fracs = CrudeSplitter.cutoff_fracs.copy()
    crude_and_char0 = sum(default_fracs[1:])
    
    gas0, aq0 = HTL_yields['gas'], HTL_yields['aqueous']
    # HTL_yields['biocrude'] includes the residuals/char after the first distillation
    char0 = HTL_yields['char'] + HTL_yields['biocrude']/crude_and_char0*default_fracs[-1]
    non_crudes = [gas0, aq0, char0]
    non_crude0 = sum(non_crudes)
    non_crudes = [i/non_crude0 for i in non_crudes]
    
    def adjust_yield(y_crude):
        non_crude = 1 - y_crude
        return [i*non_crude for i in non_crudes]

    feedstock = stream.feedstock
    mixed_fuel = stream.mixed_fuel
    dry_feedstock = feedstock.F_mass - feedstock.imass['Water']
    
    crudes = []
    fuel_yields = []
    MFSPs = []
    GWPs = []
    for y in yields:
        sys.reset_cache()
        crude = y/100 if y>1 else y
        print(f'yield: {crude:.2%}') 
        gas, aq, char = adjust_yield(crude)

        dw_yields = HTL_yields.copy()
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
            fuel_yield = mixed_fuel.F_mass/dry_feedstock
            MFSP = get_MFSP(sys, print_msg=False)
            GWP = get_GWP(sys, print_msg=False)
            print(f'Fuel yield: {fuel_yield:.2%}; MFSP: ${MFSP:.2f}/GGE; GWP: {GWP:.2f} kg CO2e/GGE.\n')
        except:
            fuel_yield = MFSP = GWP = None
            print('Simulation failed.\n')
        
        crudes.append(crude)
        fuel_yields.append(fuel_yield)
        MFSPs.append(MFSP)
        GWPs.append(GWP)
    
    df = pd.DataFrame({
        'biocrude_yield': crudes,
        'fuel_yield': fuel_yields,
        'MFSP': MFSPs,
        'GWP': GWPs,
        })

    return df


if __name__ == '__main__':
    # config_kwargs = {'include_PSA': False, 'include_EC': False,}
    # config_kwargs = {'include_PSA': True, 'include_EC': False,}
    config_kwargs = {'include_PSA': True, 'include_EC': True,}
    flowsheet = qs.main_flowsheet
    dct = globals()
    dct.update(flowsheet.to_dict())
    
    # normalized from the 80.2% biocrude+char with 3.39+81.04% light/medium distillate
    single = [0.802*(0.0339+0.8104)]
    df = evaluation_across_biocrude_yields(yields=single, **config_kwargs)
    
    # yields = np.arange(1, 100, 1)
    # df = evaluation_across_biocrude_yields(yields=yields, **config_kwargs)
    # outputs_path = os.path.join(results_path, f'biocrude_yields_{flowsheet}.csv')
    # df.to_csv(outputs_path)

    # test_df = evaluation_across_biocrude_yields(yields=df.y_test, **config_kwargs)
    # outputs_path = os.path.join(results_path, f'biocrude_yields_{flowsheet}_test.csv')
    # test_df.to_csv(outputs_path)

    # pred_df = evaluation_across_biocrude_yields(yields=df.y_pred, **config_kwargs)
    # outputs_path = os.path.join(results_path, f'biocrude_yields_{flowsheet}_pred.csv')
    # pred_df.to_csv(outputs_path)