# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

import time as tm, pandas as pd, numpy as np
from exposan.werf import create_system, add_performance_metrics, add_OPEX_metrics, baseline_underflows, opt_underflows
from exposan.werf.utils import load_state
from qsdsan import Model
from qsdsan.utils import get_SRT
from biosteam.evaluation._utils import var_columns

import warnings
warnings.filterwarnings("ignore")
#%%
def display_metrics(model):
    vals = [m() for m in model.metrics]
    print(f'{vals[-1]:.2f}')
    idx = var_columns(model.metrics)
    df = pd.DataFrame(vals, index=idx, columns=[model._system.ID])
    return df

MGD2cmd = 3785.412
#%%

ID = 'G1'

sys = create_system(ID)
s = sys.flowsheet.stream
u = sys.flowsheet.unit
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
# load_state(sys, folder='steady_states/baseline_unopt')
    
#%%
def set_do(c):
    for unit in (u.O5, u.O6):
        unit.aeration = c
        unit._ODE = None
    sys._DAE = None

def set_carb(c):
    s.carbon.imass['S_A'] = c
    s.carbon._init_state()

#%%
print(f"System {ID} Operation Adjusted")
c0 = s.carbon.imass['S_A']
cmps = s.SE.components
carbs = np.linspace(c0, c0-40, 5)
dos = np.linspace(2.0, 0.5, 6)
outs_all = []
for do in dos[1:]:
    set_do(do)
    print("="*30)
    print(f"{do} mg/L")
    outs = []
    for c in carbs[1:]:
        set_carb(c)
        start = tm.time()
        print(f"{c} kg/hr Start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
        try:
            try: sys.simulate(t_span=(0,300), method='BDF')
            except: sys.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
            end = tm.time()
            print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end-start)))
            print('Adjusting TS% ...')
            r_thick = thickened.get_TSS()/5e4
            r_cake = s.cake.get_TSS()/cake_tss
            while 1 > abs(r_thick - 1) > 0.01 or 1 > abs(r_cake - 1) > 0.01:
                print(f"{r_thick:.3f}  {r_cake:.3f}")
                thickener.sludge_flow_rate *= r_thick
                u.DW.sludge_flow_rate *= r_cake
                sys.simulate(t_span=(0,300), method='BDF')
                try: sys.simulate(t_span=(0,300), method='BDF')
                except: sys.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
                r_thick = thickened.get_TSS()/5e4
                r_cake = s.cake.get_TSS()/cake_tss
            if abs(r_thick - 1) > 1 or abs(r_cake - 1) > 1: 
                outs.append([])
                print(f"{r_thick:.3f}  {r_cake:.3f}")
                continue
                
            end2 = tm.time()
            print("Final underflows: ", f"{thickener.sludge_flow_rate:.2f}  {u.DW.sludge_flow_rate:.2f}")
            print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end2-end)), '\n')
            
            outs.append([m() for m in mdl.metrics])
            srt = get_SRT(sys, ('X_H', 'X_PAO', 'X_AUT'), wastage=[s.WAS], 
                          active_unit_IDs=('A1', 'A2', 'A3', 'A4', 'O5', 'O6'))
            print(f'SRT = {srt:.2f} d')
            if 'ASR' in u: arr = u.ASR.state.iloc[:,:-1].to_numpy()
            else: arr = [unit._state[:-1] for unit in (u.A1, u.A2, u.A3, u.A4, u.O5, u.O6)]
            mlss = np.sum(cmps.i_mass * cmps.x * arr, axis=1)
            print(f'MLSS ~ {np.mean(mlss):.0f} mg/L\n')
        
        except Exception as exc:
            print(exc, '\n')
            outs.append([])
            
    df = pd.DataFrame(outs, index=carbs[1:], columns=var_columns(mdl.metrics))
    outs_all.append(df)
#%%
dfs = pd.concat(outs_all, keys=dos)
dfs.to_clipboard()
