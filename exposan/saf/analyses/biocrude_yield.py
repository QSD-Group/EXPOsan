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

import os, pandas as pd, qsdsan as qs
from exposan.saf import (
    data_path,
    results_path,
    create_system,
    get_MFSP,
    )

data_path = os.path.join(data_path, 'biocrude_yields.csv')
df = pd.read_csv(data_path)

def MFSP_across_biocrude_yields(yields=[]):
    sys = create_system()
    unit = sys.flowsheet.unit
    
    HTL = unit.HTL
    default_yield = HTL.dw_yields.copy()
    
    GGEs = []
    for y in yields:
        print(f'yield: {y}')
        sys.reset_cache()
        dct = default_yield.copy()
        dct['biocrude'] = y/100
        HTL.dw_yields = dct
        try: 
            sys.simulate()
            MFSP = get_MFSP(sys, print_msg=False)
            print(f'MFSP: ${MFSP:.2f}/GGE\n')
        except:
            MFSP = None
            print('simulation failed.\n')
        GGEs.append(MFSP)
    
    return GGEs


if __name__ == '__main__':
    flowsheet = qs.main_flowsheet
    dct = globals()
    dct.update(flowsheet.to_dict())
    
    # GGEs = MFSP_across_biocrude_yields(yields=[80.2])
    GGEs = MFSP_across_biocrude_yields(yields=df.y_pred[:10])
    # df['$/GGE'] = GGEs
    # result_path = os.path.join(results_path, 'biocrude_yields_results.csv')
    # df.to_csv(result_path)
