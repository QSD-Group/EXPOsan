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
from exposan.werf import create_system, add_performance_metrics, add_OPEX_metrics, opt_underflows
from exposan.werf.utils import load_state, cache_state#, results_path
from qsdsan import Model, sanunits as su, processes as pc
from qsdsan.utils import get_SRT#, ospath, load_data
from biosteam.evaluation._utils import var_columns

#%%
def display_metrics(model):
    vals = [m() for m in model.metrics]
    print(f'OPEX = {vals[-1]:.2f} USD/d')
    idx = var_columns(model.metrics)
    df = pd.DataFrame(vals, index=idx, columns=[model._system.ID])
    return df

MGD2cmd = 3785.412

pct = 0.1
# pct = 0.3
# pct = 1.0
pe = 108918     # person equivalent for 10 MGD WRRF
try: EX = su.ExcretionmASM2d('EX')
except: 
    pc.create_masm2d_cmps()
    EX = su.ExcretionmASM2d('EX')
EX.simulate()
urine, feces = EX.outs
urine.scale(pe*pct)

# ufs = load_data(ospath.join(results_path, 'UD_opt_performance.xlsx'), 
#                 sheet='opt_commands', skiprows=[2], header=[0,1], index_col=0)
#%%

ID = 'B1'
# ID = 'B2'
# ID = 'B3'
# ID = 'C1'
# ID = 'C2'
# ID = 'C3'
# ID = 'E2'
# ID = 'E2P'
# ID = 'F1'
# ID = 'G1'
# ID = 'G2'
# ID = 'G3'
# ID = 'H1'
# ID = 'I1'
# ID = 'I2'
# ID = 'I3'
# ID = 'N1'
# ID = 'N2'

sys = create_system(ID)
s = sys.flowsheet.stream
u = sys.flowsheet.unit
mdl = Model(sys)
add_performance_metrics(mdl)
add_OPEX_metrics(mdl)

s.RWW.separate_out(urine)

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

# thickener.sludge_flow_rate = ufs.at[ID, (int(pct*100), 'Thickener')]
# u.DW.sludge_flow_rate = ufs.at[ID, (int(pct*100), 'Dewatering')]

thickener.sludge_flow_rate, u.DW.sludge_flow_rate = opt_underflows[ID]
# thickener.sludge_flow_rate, u.DW.sludge_flow_rate = (62.37, 11.34)
load_state(sys, folder='steady_states/baseline_unopt')
# load_state(sys, folder='steady_states/UD100_opt')
    
#%%
u.ASR.DO_setpoints *= 0
u.ASR.DO_setpoints += 1
# u.ASR.DO_setpoints[:] = [0,0,0,0,1,1]
# u.ASR.DO_setpoints[:] = [0.5,0,2,2,0,1]
# u.ASR.DO_setpoints[:] = [0,0,1,1,0,1]
# u.ASR.DO_setpoints[:] = [0,0,1,1,0]
# Vs = [0.63, 1.5, 2.0, 2.0, 2.3, 0.21] # MG
# Vs = [0.84, 1.5, 2.0, 2.0, 2.1, 0.2]
# Vs = [0.99, 1.35, 2.0, 2.0, 2.1, 0.2]
# Vs = [0.33, 1.5, 2.0, 2.1, 2.3, 0.41]
# u.ASR.V_tanks[:] = [v * MGD2cmd for v in Vs]
# V_tot = 2.61 * MGD2cmd
# fr_V = [0.12, 0.18, 0.24, 0.24, 0.18, 0.04]
# fr_V = [0.18, 0.14, 0.24, 0.24, 0.16, 0.04]
# fr_V = [0.16, 0.16, 0.24, 0.24, 0.17, 0.03]     # larger anaerobic zone seems better for EBPR
# u.ASR.V_tanks[:] = [v * V_tot for v in fr_V[:-1]]
# u.ASR.internal_recycles[0] = (3,1,20*MGD2cmd)
u.ASR._ODE = None

# Q_ras = 2 * 10 * MGD2cmd
# Q_was = 0.2 * MGD2cmd
# u.MBR.pumped_flow = Q_ras + Q_was
# u.S1.split = Q_ras / (Q_ras + Q_was)
# !!! must reset cache

# u.MBR.V_max = fr_V[-1] * V_tot
# u.MBR.aeration = 1.0
# u.MBR._ODE = None

# for unit in (u.O5, u.O6):
#     unit.aeration = 1.0
#     unit._ODE = None

# s.carbon.imass['S_A'] = 50
# s.carbon._init_state()
# u.MD.metal_dosage = 6
# u.MD._AE = None
u.FC.underflow = 0.4 * 10 * MGD2cmd
u.FC.wastage = 0.2 * MGD2cmd
u.FC._ODE = None

# u.AED.V_max = 0.5 * MGD2cmd
# u.AED._ODE = None

sys._DAE = None

#%%
print(f"System {ID} Operation Adjusted")
print("="*30)
start = tm.time()
print("Start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
# sys.simulate(t_span=(0,300), method='BDF')
sys.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF') # N1, N2
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

df = display_metrics(mdl)
srt = get_SRT(sys, ('X_H', 'X_PAO', 'X_AUT'), wastage=[s.WAS], 
              active_unit_IDs=('ASR', 'MBR', 'A1', 'A2', 'A3', 'A4', 'O5', 'O6'))
print(f'SRT = {srt:.2f} d')
cmps = s.SE.components
if 'ASR' in u: arr = u.ASR.state.iloc[:,:-1].to_numpy()
else: arr = [unit._state[:-1] for unit in (u.A1, u.A2, u.A3, u.A4, u.O5, u.O6)]
mlss = np.sum(cmps.i_mass * cmps.x * arr, axis=1)
print(f'MLSS ~ {np.mean(mlss):.0f} mg/L')

#%%
df = display_metrics(mdl)
df.T.to_clipboard()
cache_state(sys, 'steady_states/UD10_opt')
