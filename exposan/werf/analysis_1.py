# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 08:41:06 2025

@author: joy_c
"""
import time as tm, pandas as pd, os
from exposan.werf import create_system, add_performance_metrics, results_path
from exposan.werf.utils import plantwide_N_mass_flows
from qsdsan import Model, sanunits as su
from biosteam.evaluation._utils import var_columns

import warnings
warnings.filterwarnings("ignore")

N_mass_dfs = {}
metrics = {}
pe = 108918     # person equivalent for 10 MGD WRRF

for ID in (
        'B1', 'B2', 'B3', 
        'C1', 'C2', 'C3', 
        'E2', 'E2P', 
        'F1', 
        'G1', 'G2', 'G3', 
        'H1', 
        'I1', 'I2', 'I3', 
        'N1', 'N2'
        ):
    txt = f"System {ID}"
    print(txt)
    print("="*len(txt))
    sys = create_system(ID)
    fs = sys.flowsheet.stream
    mdl = Model(sys)
    add_performance_metrics(mdl)

    start = tm.time()
    print("Baseline start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
    sys.simulate(t_span=(0,300), method='BDF')
    end = tm.time()
    print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end-start)), '\n')
    N_mass_dfs[ID] = plantwide_N_mass_flows(sys)
    metrics[ID] = [m() for m in mdl.metrics]

    if 'urine' not in fs:
        EX = su.ExcretionmASM2d('EX', outs=('urine', 'feces'))
        EX.simulate()
        fs.urine.scale(pe)
    fs.RWW.separate_out(fs.urine)
    
    start = tm.time()
    print("With UD start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
    sys.simulate(state_reset_hook='reset_cache', t_span=(0,400), method='BDF')
    end = tm.time()
    print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end-start)), '\n')
    N_mass_dfs[ID+'_UD'] = plantwide_N_mass_flows(sys)
    metrics[ID+'_UD'] = [m() for m in mdl.metrics]

with pd.ExcelWriter(os.path.join(results_path, 'N_mass.xlsx')) as writer:
    for k, v in N_mass_dfs.items():
        v.to_excel(writer, sheet_name=k)

metrics = pd.DataFrame.from_dict(metrics, orient='index', columns=var_columns(mdl.metrics))
metrics.to_excel(os.path.join(results_path, 'ana1_performance.xlsx'))