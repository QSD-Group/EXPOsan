# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

import time as tm, pandas as pd
from exposan.werf import create_system, add_performance_metrics, add_OPEX_metrics, results_path, baseline_underflows
from exposan.werf.utils import load_state
from qsdsan import Model
from biosteam.evaluation._utils import var_columns

#%%
def display_metrics(model):
    vals = [m() for m in model.metrics]
    print(f'{vals[-1]:.2f}')
    idx = var_columns(model.metrics)
    df = pd.DataFrame(vals, index=idx, columns=[model._system.ID])
    return df

# underflows = {
#     'C1': (136.64, 19.67),
#     'C2': (131.76, 6.69),
#     'C3': (138.19, 37.70),
#     'E2': (107.95, 5.26),
#     'E2P': (53.92, 19.93),
#     'F1': (54.64, 19.42),
#     'G2': (89.71, 15.37),
#     'I3': (116.21, 31.65)
#     }

MGD2cmd = 3785.412
#%%

# ID = 'C1'
# ID = 'C2'
# ID = 'C3'
# ID = 'E2'
# ID = 'E2P'
# ID = 'F1'
# ID = 'G2'
# ID = 'I3'

# ID = 'B1'
# ID = 'B2'
# ID = 'B3'
# ID = 'G1'
# ID = 'G3'
# ID = 'H1'
# ID = 'I1'
# ID = 'I2'
# ID = 'N1'
ID = 'N2'

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

#%%
try:
    print(f"System {ID}")
    print("="*30)
    start = tm.time()
    print("Start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
    load_state(sys, folder='steady_states/baseline_unopt')
    thickener.sludge_flow_rate, u.DW.sludge_flow_rate = baseline_underflows[ID]
    try: sys.simulate(t_span=(0,300), method='BDF')
    except: sys.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
    end = tm.time()
    print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end-start)))
    display_metrics(mdl)
    
except Exception as exc:
    print(exc, '\n')
    
#%%
# u.ASR.DO_setpoints *= 0
# u.ASR.DO_setpoints += 1
# u.ASR.DO_setpoints[:] = [0,0,1,1,0,1]
u.ASR.DO_setpoints[:] = [0,0,1,1,0]
# Vs = [0.63, 1.5, 2.0, 2.0, 2.3, 0.21] # MG
# Vs = [0.84, 1.5, 2.0, 2.0, 2.1, 0.2]
# Vs = [0.99, 1.35, 2.0, 2.0, 2.1, 0.2]
# Vs = [0.3, 1.5, 2.03, 2.1, 2.3, 0.41]
# u.ASR.V_tanks[:] = [v * MGD2cmd for v in Vs]
V_tot = 2.61 * MGD2cmd
# fr_V = [0.12, 0.18, 0.24, 0.24, 0.18, 0.04]
fr_V = [0.16, 0.16, 0.24, 0.24, 0.17, 0.03]     # larger anaerobic zone seems better for EBPR
u.ASR.V_tanks[:] = [v * V_tot for v in fr_V[:-1]]
# u.ASR.internal_recycles[0] = (3,1,20*MGD2cmd)
u.ASR._ODE = None

u.MBR.V_max = fr_V[-1] * V_tot
u.MBR.aeration = 1.0
u.MBR._ODE = None


# for unit in (u.O5, u.O6):
#     unit.aeration = 1.0
#     unit._ODE = None

s.carbon.imass['S_A'] = 70
s.carbon._init_state()
# u.MD.metal_dosage = 12
# u.MD._AE = None
# u.FC.underflow = 0.5 * 10 * MGD2cmd
# u.FC.wastage = 0.3 * MGD2cmd
# u.FC._ODE = None

# u.AED.V_max = 0.5 * MGD2cmd
# u.AED._ODE = None

sys._DAE = None

#%%
print(f"System {ID} Operation Adjusted")
print("="*30)
start = tm.time()
print("Start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
sys.simulate(t_span=(0,300), method='BDF')
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

display_metrics(mdl)

#%%
df = display_metrics(mdl)
df.T.to_clipboard()
