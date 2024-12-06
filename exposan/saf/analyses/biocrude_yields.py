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
    config_baseline,
    config_EC,
    config_EC_future,
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
    non_crudes0 = [HTL_yields['gas'], HTL_yields['aqueous'], HTL_yields['char']]
    non_crude0 = sum(non_crudes0)
    non_crudes0 = [i/non_crude0 for i in non_crudes0]
    
    def adjust_yield(y_crude):
        non_crude = 1 - y_crude
        return [i*non_crude for i in non_crudes0]

    feedstock = stream.feedstock
    mixed_fuel = stream.mixed_fuel
    dry_feedstock = feedstock.F_mass - feedstock.imass['Water']
    
    crudes = []
    fuel_yields = []
    MFSPs = []
    GWPs = []
    for y in yields:
        sys.reset_cache()
        crude = y/100 if y>=1 else y
        print(f'yield: {crude:.2%}') 
        gas, aq, char = adjust_yield(crude)

        HTL.dw_yields = {
            'gas': gas,
            'aqueous': aq,
            'biocrude': crude,
            'char': char,
            }

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
    # config_kwargs = config_baseline
    # config_kwargs = config_EC
    config_kwargs = config_EC_future
    
    flowsheet = qs.main_flowsheet
    dct = globals()
    dct.update(flowsheet.to_dict())
    
    # Original setting, char subtracted
    # single = [0.802-0.0614]
    # df = evaluation_across_biocrude_yields(yields=single, **config_kwargs)
    
    yields = np.arange(1, 100, 1)
    df = evaluation_across_biocrude_yields(yields=yields, **config_kwargs)
    outputs_path = os.path.join(results_path, f'biocrude_yields_{flowsheet}.csv')
    df.to_csv(outputs_path)

    ### Below are for the ML paper ###
    # test_df = evaluation_across_biocrude_yields(yields=df.y_test, **config_kwargs)
    # outputs_path = os.path.join(results_path, f'biocrude_yields_{flowsheet}_test.csv')
    # test_df.to_csv(outputs_path)

    # pred_df = evaluation_across_biocrude_yields(yields=df.y_pred, **config_kwargs)
    # outputs_path = os.path.join(results_path, f'biocrude_yields_{flowsheet}_pred.csv')
    # pred_df.to_csv(outputs_path)