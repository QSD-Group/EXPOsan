# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""
import time as tm, pandas as pd, os
from exposan.werf import create_system, add_performance_metrics, results_path
from exposan.werf.utils import plantwide_N_mass_flows, plantwide_P_mass_flows
from qsdsan import Model, WasteStream, sanunits as su, processes as pc
from biosteam.evaluation._utils import var_columns

import warnings
warnings.filterwarnings("ignore")

metrics = {'0':{}, '10':{}, '30':{}, '100':{}}
pe = 108918     # person equivalent for 10 MGD WRRF

pc.create_masm2d_cmps()
EX = su.ExcretionmASM2d('EX')
EX.simulate()
urine, feces = EX.outs
urine.scale(pe*0.1)
_urine = WasteStream()

mode = 'w'
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
    print(f"System {ID}")
    print("="*30)
    sys = create_system(ID)
    fs = sys.flowsheet.stream
    _urine.copy_like(urine)
    mdl = Model(sys)
    add_performance_metrics(mdl)

    for pud, multiple in zip((0, 10, 30, 100), (0,1,2,3.5)):
        try:
            if pud > 0:
                _urine.scale(multiple)
                fs.RWW.separate_out(_urine)
            start = tm.time()
            print(f"{pud}% UD start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
            sys.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
            end = tm.time()
            print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end-start)), '\n')

            ndf = plantwide_N_mass_flows(sys)
            with pd.ExcelWriter(os.path.join(results_path, f'N_mass_{pud}UD.xlsx'), mode=mode) as writer:
                ndf.to_excel(writer, sheet_name=ID)
            pdf = plantwide_P_mass_flows(sys)
            with pd.ExcelWriter(os.path.join(results_path, f'P_mass_{pud}UD.xlsx'), mode=mode) as writer:
                pdf.to_excel(writer, sheet_name=ID)
            metrics[f'{pud}'][ID] = [m() for m in mdl.metrics]
        except Exception as exc:
            print(exc, '\n')
    
    mode = 'a'
    
    del sys, fs

#%%


with pd.ExcelWriter(os.path.join(results_path, 'UD_performance.xlsx')) as writer:
    for k, v in metrics.items():
        df = pd.DataFrame.from_dict(v, orient='index', columns=var_columns(mdl.metrics))
        df.to_excel(writer, sheet_name=k)
