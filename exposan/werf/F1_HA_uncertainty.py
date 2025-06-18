# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import time as tm, pandas as pd, numpy as np
from exposan.werf import (
    create_system, 
    SelectiveRecovery,
    add_performance_metrics, 
    add_OPEX_metrics, 
    add_NH4_recovery_metric,
    opt_underflows,
    results_path
    )
from exposan.werf.utils import cache_state, load_state
from qsdsan import Model, System, processes as pc
from qsdsan.utils import get_SRT, ospath
from biosteam.evaluation._utils import var_columns

import warnings
warnings.filterwarnings("ignore")

MGD2cmd = 3785.412
f_rmv = 0.7 # (0.5, 0.8)

# influent strength: low, medium, high
# tech performance: % NH4-N recovery (50-80% reasonable range for ionic strength typical in domestic ww)
# maintain operation as adjusted unless violation?

# contextual uncertainties
# %%
ID = 'F1'
sys = create_system(ID)
s = sys.flowsheet.stream
u = sys.flowsheet.unit
cmps = s.RWW.components
for c in cmps:
    if c.organic:
        if c.i_N: c.i_N *= 0.7
        if c.i_P: c.i_P *= 1.1
cmps.S_I.f_Vmass_Totmass = cmps.X_I.f_Vmass_Totmass = 0.6
cmps.X_S.i_mass = cmps.X_I.i_mass = 0.63
cmps.refresh_constants()
fr_orgs = dict(fr_SI=0.075, fr_SF=0.24, fr_SA=0.05, fr_XI=0.3)
wws = {}
wws['low'] = pc.create_masm2d_inf(
    'low', 10, 'MGD', 
    COD=339, NH4_N=14, PO4_P=1.6, S_K=11, S_Cl=39, 
    **fr_orgs
    )
wws['mid'] = pc.create_masm2d_inf(
    'mid', 10, 'MGD', 
    COD=508, NH4_N=20, PO4_P=2.4, S_K=16, S_Cl=59, 
    **fr_orgs
    )
wws['high'] = pc.create_masm2d_inf(
    'high', 10, 'MGD', 
    COD=1016, NH4_N=41, PO4_P=4.7, S_K=32, S_Cl=118, 
    **fr_orgs
    )

#%%
location_insert = s.PE
downstream_unit = location_insert.sink
i_inlet = downstream_unit.ins.index(location_insert)

i_path = sys.path.index(downstream_unit)
location_insert.disconnect_sink()
HA = SelectiveRecovery(
    'HydrogelAbsorbent', ins=location_insert, 
    outs=('Recovered_NH4', 'HA_eff'), 
    split={'S_NH4': f_rmv}, init_with='WasteStream'
    )
HA-1-i_inlet-downstream_unit
sys_ha = System(ID+'ha', 
                path=(*sys.path[:i_path], HA, *sys.path[i_path:]),
                recycle=sys.recycle)
sys_ha.set_dynamic_tracker(*sys.scope.subjects)
mdl = Model(sys_ha)
add_performance_metrics(mdl)
add_OPEX_metrics(mdl)
add_NH4_recovery_metric(mdl)

cake_tss = 18e4
thickened = s.thickened_WAS
thickener = u.MT

thickener.sludge_flow_rate, u.DW.sludge_flow_rate = opt_underflows[ID]
u.ASR.DO_setpoints *= 0
u.ASR.DO_setpoints += 1

load_state(sys_ha, folder='steady_states/HA_opt')

#%%
removal_efficiencies = np.arange(0.5, 0.8, 0.05)
metrics = {}
for strength, qwas in zip(('low', 'mid', 'high'), (0.15, 0.2, 0.3)):
    s.RWW.copy_flow(wws[strength])
    s.RWW._init_state()
    u.FC.wastage = qwas * MGD2cmd
    u.FC._ODE = None
    sys_ha._DAE = None
    out = []
    print(f"System {ID} + HA, {strength}-strength ww")
    print("="*35)
    for f_rmv in removal_efficiencies:
        HA.split[cmps.index('S_NH4')] = f_rmv
        start = tm.time()
        print(f"{f_rmv:.0%} Start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
        try: sys_ha.simulate(t_span=(0,300), method='BDF')
        except: sys_ha.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
        end = tm.time()
        print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end-start)))
        print('Adjusting TS% ...')
        r_thick = thickened.get_TSS()/5e4
        r_cake = s.cake.get_TSS()/cake_tss
        while abs(r_thick - 1) > 0.01 or abs(r_cake - 1) > 0.01:
            print(f"{r_thick:.3f}  {r_cake:.3f}")
            thickener.sludge_flow_rate *= r_thick
            u.DW.sludge_flow_rate *= r_cake
            try: sys_ha.simulate(t_span=(0,300), method='BDF')
            except: sys_ha.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
            r_thick = thickened.get_TSS()/5e4
            r_cake = s.cake.get_TSS()/cake_tss
        end2 = tm.time()
        print("Final underflows: ", f"{thickener.sludge_flow_rate:.2f}  {u.DW.sludge_flow_rate:.2f}")
        print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end2-end)))
        srt = get_SRT(sys_ha, ('X_H', 'X_PAO', 'X_AUT'), wastage=[s.WAS], 
                      active_unit_IDs=('ASR', ))
        print(f'SRT = {srt:.2f} d')
        arr = u.ASR.state.iloc[:,:-1].to_numpy()
        mlss = np.sum(cmps.i_mass * cmps.x * arr, axis=1)
        print(f'MLSS = {np.mean(mlss):.0f} mg/L', '\n')
        cache_state(sys_ha, f'steady_states/HA_F1/{strength}', f'{f_rmv*100:.0f}')
        out.append([m() for m in mdl.metrics])
    metrics[strength] = pd.DataFrame(out, index=removal_efficiencies, 
                                     columns=var_columns(mdl.metrics))

#%%
df = pd.concat(metrics.values(), keys=metrics.keys(), 
               names=['WW strength', 'Recovery efficiency'])
df.to_excel(ospath.join(results_path, 'HA_F1_upstream_uncertainty.xlsx'))