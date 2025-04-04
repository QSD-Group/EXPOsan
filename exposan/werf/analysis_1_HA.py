# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

#%%
import time as tm, pandas as pd, os
from exposan.werf import (
    create_system, 
    SelectiveRecovery,
    add_performance_metrics, 
    add_NH4_recovery_metric, 
    results_path
    )
from exposan.werf.utils import plantwide_N_mass_flows, plantwide_P_mass_flows
from qsdsan import System, Model, sanunits as su
from biosteam.evaluation._utils import var_columns

import warnings
warnings.filterwarnings("ignore")

f_rmv = 0.7

mode = 'w'
metrics = {}
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
    print(f"System {ID} w Hydrogel Absorbent")
    print("="*30)
    sys = create_system(ID)
    u = sys.flowsheet.unit
    s = sys.flowsheet.stream
    
    if "PE" in s: 
        location_insert = s.PE
        downstream_unit = location_insert.sink
        i_inlet = downstream_unit.ins.index(location_insert)
    else: 
        downstream_unit = s.reject.source
        assert isinstance(downstream_unit, su.HydraulicDelay)
        i_inlet = 0
        location_insert = downstream_unit.ins[i_inlet]
        s.RWW.imass['S_NH4'] *= (1-f_rmv)

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
    add_NH4_recovery_metric(mdl)
    
    try:
        start = tm.time()
        print("Start time: ", tm.strftime('%H:%M:%S', tm.localtime()))
        try: sys_ha.simulate(state_reset_hook='reset_cache', t_span=(0,300), method='BDF')
        except: sys_ha.simulate(t_span=(0,300), method='BDF')
        end = tm.time()
        print('Duration: ', tm.strftime('%H:%M:%S', tm.gmtime(end-start)), '\n')
        
        sys_ha.diagram()
        ndf = plantwide_N_mass_flows(sys_ha)
        with pd.ExcelWriter(os.path.join(results_path, 'N_mass_HA.xlsx'), mode=mode) as writer:
            ndf.to_excel(writer, sheet_name=ID)
        pdf = plantwide_P_mass_flows(sys_ha)
        with pd.ExcelWriter(os.path.join(results_path, 'P_mass_HA.xlsx'), mode=mode) as writer:
            pdf.to_excel(writer, sheet_name=ID)
        metrics[ID] = [m() for m in mdl.metrics]
    except Exception as exc:
        print(exc, '\n')
    
    mode = 'a'
    sys.flowsheet.clear()
    sys_ha.flowsheet.clear()
    del sys, sys_ha

metrics = pd.DataFrame.from_dict(metrics, orient='index', columns=var_columns(mdl.metrics))
metrics.to_excel(os.path.join(results_path, 'HA_performance.xlsx'))