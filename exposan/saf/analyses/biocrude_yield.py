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
    gasoline = stream.gasoline
    jet = stream.jet
    diesel = stream.diesel
    dry_feedstock = feedstock.F_mass - feedstock.imass['Water']
    
    GGEs = []
    y_gasolines = []
    y_jets = []
    y_diesels = []
    
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
            y_gasoline = gasoline.F_mass/dry_feedstock
            y_jet = jet.F_mass/dry_feedstock
            y_diesel = diesel.F_mass/dry_feedstock
            print(f'MFSP: ${MFSP:.2f}/GGE\n')
        except:
            MFSP = y_gasoline = y_jet = y_diesel = None
            print('simulation failed.\n')

        GGEs.append(MFSP)
        y_gasolines.append(y_gasoline)
        y_jets.append(y_jet)
        y_diesels.append(y_diesel)
    
    return GGEs, y_gasolines, y_jets, y_diesels


if __name__ == '__main__':
    flowsheet = qs.main_flowsheet
    dct = globals()
    dct.update(flowsheet.to_dict())
    
    # single=[67.3] # normalized from the 80.2 biocrude+char, $2.67/GGE
    # single=[22.304035] # $3.91/GGE
    # results = MFSP_across_biocrude_yields(yields=single)
    
    tested_results = MFSP_across_biocrude_yields(yields=df.y_test)
    tested = df.copy()
    tested['$/GGE'] = tested_results[0]
    tested['y_gasoline'] = tested_results[1]
    tested['y_jet'] = tested_results[2]
    tested['y_diesel'] = tested_results[3]
    tested_path = os.path.join(results_path, 'tested_results.csv')
    tested.to_csv(tested_path)
    
    predicted_results = MFSP_across_biocrude_yields(yields=df.y_pred)
    predicted = df.copy()
    predicted['$/GGE'] = predicted_results[0]
    predicted['y_gasoline'] = predicted_results[1]
    predicted['y_jet'] = predicted_results[2]
    predicted['y_diesel'] = predicted_results[3]
    predicted_path = os.path.join(results_path, 'predicted_results.csv')
    predicted.to_csv(predicted_path)
