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
from exposan.werf import (
    create_system, 
    SelectiveRecovery,
    add_performance_metrics, 
    add_NH4_recovery_metric, 
    results_path
    )
from exposan.werf.utils import plantwide_N_mass_flows, plantwide_P_mass_flows
from qsdsan import System, Model
from biosteam.evaluation._utils import var_columns

import warnings
warnings.filterwarnings("ignore")

f_rmv = 0.9

mode = 'w'
metrics = {}
for ID in (
        'B1', 'C1', 'F1', 
        'G1', 'H1', 'I1', 'N1', 
        ):
    print(f"System {ID} w Electrochemical Stripping")
    print("="*40)
    sys = create_system(ID)
    u = sys.flowsheet.unit
    s = sys.flowsheet.stream
    
    location_insert = u.DW.outs[0] 
    downstream_unit = location_insert.sink
    i_inlet = downstream_unit.ins.index(location_insert)
    i_path = sys.path.index(downstream_unit)
    
    location_insert.disconnect_sink()
    ECS = SelectiveRecovery(
        'ElectrochemStripping', ins=location_insert, 
        outs=('Recovered_NH4', 'ECS_eff'), 
        split={'S_NH4': f_rmv}, init_with='WasteStream'
        )
    ECS-1-i_inlet-downstream_unit
    sys_ecs = System(ID+'ecs', 
                    path=(*sys.path[:i_path], ECS, *sys.path[i_path:]),
                    recycle=sys.recycle)
    sys_ecs.set_dynamic_tracker(*sys.scope.subjects)
    mdl = Model(sys_ecs)
    add_performance_metrics(mdl)
    add_NH4_recovery_metric(mdl)
    if "PC" in u: cake_tss = 18e4
    else: cake_tss = 17e4
    
    if "thickened_WAS" in s: 
        thickened = s.thickened_WAS
        thickener = u.MT
    else: 
        thickened = s.thickened_sludge
        thickener = u.GT
    
    try:
        start = tm.time()
        print("Start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
        try: sys_ecs.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
        except: sys_ecs.simulate(t_span=(0,300), method='BDF')
        end = tm.time()
        print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end-start)))
        print('Adjusting TS% ...')
        r_thick = thickened.get_TSS()/5e4
        r_cake = s.cake.get_TSS()/cake_tss
        while abs(r_thick - 1) > 0.01 or abs(r_cake - 1) > 0.01:
            print(f"{r_thick:.3f}  {r_cake:.3f}")
            thickener.sludge_flow_rate *= r_thick
            u.DW.sludge_flow_rate *= r_cake
            sys_ecs.simulate(t_span=(0,300), method='BDF')
            r_thick = thickened.get_TSS()/5e4
            r_cake = s.cake.get_TSS()/cake_tss
        end2 = tm.time()
        print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end2-end)), '\n')

        # sys_ecs.diagram()
        ndf = plantwide_N_mass_flows(sys_ecs)
        with pd.ExcelWriter(os.path.join(results_path, 'N_mass_ECS.xlsx'), mode=mode) as writer:
            ndf.to_excel(writer, sheet_name=ID)
        pdf = plantwide_P_mass_flows(sys_ecs)
        with pd.ExcelWriter(os.path.join(results_path, 'P_mass_ECS.xlsx'), mode=mode) as writer:
            pdf.to_excel(writer, sheet_name=ID)
        metrics[ID] = [m() for m in mdl.metrics]
    except Exception as exc:
        print(exc, '\n')
    
    mode = 'a'
    sys.flowsheet.clear()
    sys_ecs.flowsheet.clear()
    del sys, sys_ecs

metrics = pd.DataFrame.from_dict(metrics, orient='index', columns=var_columns(mdl.metrics))
metrics.to_excel(os.path.join(results_path, 'ECS_performance.xlsx'))