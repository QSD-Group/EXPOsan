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
from exposan.werf import create_system, add_performance_metrics, add_OPEX_metrics, results_path, baseline_underflows
from exposan.werf.utils import plantwide_N_mass_flows, plantwide_P_mass_flows, load_state
from qsdsan import Model, WasteStream, sanunits as su, processes as pc
from biosteam.evaluation._utils import var_columns

import warnings
warnings.filterwarnings("ignore")

metrics = {'10':{}, '30':{}, '100':{}}
pe = 108918     # person equivalent for 10 MGD WRRF

pc.create_masm2d_cmps()
EX = su.ExcretionmASM2d('EX')
EX.simulate()
urine, feces = EX.outs
urine.scale(pe*0.1)
_urine = WasteStream()

kwargs = dict(
    mode = 'w',
    if_sheet_exists=None
    )
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
    load_state(sys, folder='steady_states/baseline_unopt')
    s = sys.flowsheet.stream
    u = sys.flowsheet.unit
    _urine.copy_like(urine)
    mdl = Model(sys)
    add_performance_metrics(mdl)
    add_OPEX_metrics(mdl)
    

    if "PS" in s: 
        if 'AD' in u: cake_tss = 18e4
        elif 'AED' in u: cake_tss = 17e4
        else: cake_tss = 20e4
    else: 
        cake_tss = 17e4
    
    if "thickened_WAS" in s: 
        thickened = s.thickened_WAS
        thickener = u.MT
    else: 
        thickened = s.thickened_sludge
        thickener = u.GT
    
    thickener.sludge_flow_rate, u.DW.sludge_flow_rate = baseline_underflows[ID]

    for pud, multiple in zip((10, 30, 100), (1,2,3.5)):
        try:
            _urine.scale(multiple)
            s.RWW.separate_out(_urine)
            s.RWW._init_state()
            start = tm.time()
            print(f"{pud}% UD start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
            try: sys.simulate(t_span=(0,300), method='BDF')
            except: sys.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
            end = tm.time()
            print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end-start)))
            print('Adjusting TS% ...')
            r_thick = thickened.get_TSS()/5e4
            r_cake = s.cake.get_TSS()/cake_tss
            while abs(r_thick - 1) > 0.01 or abs(r_cake - 1) > 0.01:
                print(f"{r_thick:.3f}  {r_cake:.3f}")
                thickener.sludge_flow_rate *= r_thick
                u.DW.sludge_flow_rate *= r_cake
                try: sys.simulate(t_span=(0,300), method='BDF')
                except: sys.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
                r_thick = thickened.get_TSS()/5e4
                r_cake = s.cake.get_TSS()/cake_tss
            end2 = tm.time()
            print("Final underflows: ", f"{thickener.sludge_flow_rate:.2f}  {u.DW.sludge_flow_rate:.2f}")
            print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end2-end)), '\n')
            
            ndf = plantwide_N_mass_flows(sys)
            with pd.ExcelWriter(os.path.join(results_path, f'N_mass_{pud}UD.xlsx'), **kwargs) as writer:
                ndf.to_excel(writer, sheet_name=ID)
            pdf = plantwide_P_mass_flows(sys)
            with pd.ExcelWriter(os.path.join(results_path, f'P_mass_{pud}UD.xlsx'), **kwargs) as writer:
                pdf.to_excel(writer, sheet_name=ID)
            metrics[f'{pud}'][ID] = [m() for m in mdl.metrics]
        except Exception as exc:
            print(exc, '\n')
    
    kwargs['mode'] = 'a'
    kwargs['if_sheet_exists'] = 'replace'
    sys.flowsheet.clear()
    del sys
    
#%%


with pd.ExcelWriter(os.path.join(results_path, 'UD_performance.xlsx')) as writer:
    for k, v in metrics.items():
        df = pd.DataFrame.from_dict(v, orient='index', columns=var_columns(mdl.metrics))
        df.to_excel(writer, sheet_name=k)
