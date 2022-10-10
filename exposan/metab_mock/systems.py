# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
import qsdsan as qs
from qsdsan import sanunits as su, processes as pc, WasteStream, System
from qsdsan.utils import time_printer, ospath
from chemicals.elements import molecular_weight as get_mw
from exposan.metab_mock import DegassingMembrane as DM, METAB_AnCSTR as AB

folder = ospath.dirname(__file__)

__all__ = (
    'create_systems', 
    'default_inf_concs',
    'default_R1_init_conds',
    'default_R2_init_conds',
    'yields_bl', 'mus_bl', 'Ks_bl'
    )

#%% default values
Q = 5           # influent flowrate [m3/d]
T1 = 273.15+35  # temperature [K]
Vl1 = 5         # liquid volume [m^3]
Vg1 = 0.556     # headspace volume [m^3]
split_1 = 0.75  # split ratio to side-stream
tau_1 = 0.021   # degassing membrane retention time [d]

T2 = 273.15+25    
Vl2 = 75
Vg2 = 5
split_2 = 0.75
tau_2 = 0.021

C_mw = get_mw({'C':1})
N_mw = get_mw({'N':1})

default_inf_concs = {
    'S_su':3.0,
    'S_aa':0.6,
    'S_fa':0.4,
    'S_va':0.4,
    'S_bu':0.4,
    'S_pro':0.4,
    'S_ac':0.4,
    'S_h2':5e-9,
    'S_ch4':5e-6,
    'S_IC':0.04*C_mw,
    'S_IN':0.01*N_mw,
    'S_I':0.02,
    'X_c':0.1,
    'X_ch':0.3,
    'X_pr':0.5,
    'X_li':0.25,
    'X_aa':1e-3,
    'X_fa':1e-3,
    'X_c4':1e-3,
    'X_pro':1e-3, 
    'X_ac':1e-3, 
    'X_h2':1e-3, 
    'X_I':0.025, 
    'S_cat':0.04, 
    'S_an':0.02
    }

yields_bl = {
         'Y_su': 0.1,
         'Y_aa': 0.08,
         'Y_fa': 0.06,
         'Y_c4': 0.06,
         'Y_pro': 0.04,
         'Y_ac': 0.05,
         'Y_h2': 0.06
         }

mus_bl = np.array([5.0e-01, 1.0e+01, 1.0e+01, 1.0e+01, 3.0e+01, 5.0e+01, 6.0e+00,
                   2.0e+01, 2.0e+01, 1.3e+01, 8.0e+00, 3.5e+01, 2.0e-02, 2.0e-02,
                   2.0e-02, 2.0e-02, 2.0e-02, 2.0e-02, 2.0e-02])

Ks_bl = np.array([5.0e-01, 3.0e-01, 4.0e-01, 2.0e-01, 
                  2.0e-01, 1.0e-01, 1.5e-01, 7.0e-06])

default_R1_init_conds = {
    'S_su': 0.0124*1e3,
    'S_aa': 0.0055*1e3,
    'S_fa': 0.1074*1e3,
    'S_va': 0.0123*1e3,
    'S_bu': 0.0140*1e3,
    'S_pro': 0.0176*1e3,
    'S_ac': 0.0893*1e3,
    'S_h2': 2.5055e-7*1e3,
    'S_ch4': 0.0555*1e3,
    'S_IC': 0.0951*C_mw*1e3,
    'S_IN': 0.0945*N_mw*1e3,
    'S_I': 0.1309*1e3,
    'X_ch': 0.0205*1e3,
    'X_pr': 0.0842*1e3,
    'X_li': 0.0436*1e3,
    'X_su': 1.87*1e3,
    'X_aa': 5.58*1e3,
    'X_fa': 2.03*1e3,
    'X_c4': 2.15*1e3,
    'X_pro': 1.00*1e3,
    }

default_R2_init_conds = {
    'S_su': 0.0124*1e3,
    'S_aa': 0.0055*1e3,
    'S_fa': 0.1074*1e3,
    'S_va': 0.0123*1e3,
    'S_bu': 0.0140*1e3,
    'S_pro': 0.0176*1e3,
    'S_ac': 0.0893*1e3,
    'S_h2': 2.5055e-7*1e3,
    'S_ch4': 0.0555*1e3,
    'S_IC': 0.0951*C_mw*1e3,
    'S_IN': 0.0945*N_mw*1e3,
    'S_I': 0.1309*1e3,
    'X_ac': 8.80*1e3,
    'X_h2': 3.70*1e3,
    }

# H2E = su.AnaerobicCSTR('H2E', ins=[brewery_ww, 'return_1'], outs=('headspace_1', ''), 
#                       V_liq=5, V_gas=0.556, T=T1, model=adm1, 
#                       retain_cmps=('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro'))
# S1 = su.Splitter('S1', ins=H2E-1, outs=('sidestream_1', ''),
#                   split=split_1, isdynamic=True)
# DM1 = DM('DM1', ins=S1-0, outs=(bg1, 1-H2E), 
#           H2_degas_efficiency=0.5,
#           CH4_degas_efficiency=0, 
#           CO2_degas_efficiency=0.2)

# CH4E = su.AnaerobicCSTR('CH4E', ins=[S1-1, 'return_2'], outs=('headspace_2', ''), 
#                         V_liq=75, V_gas=5, T=T2, model=adm1,
#                         retain_cmps=('X_ac', 'X_h2'))
# S2 = su.Splitter('S2', ins=CH4E-1, outs=('sidestream_2', eff),
#                   split=split_2, isdynamic=True)
# DM2 = DM('DM2', ins=S2-0, outs=(bg2, 1-CH4E),
#           H2_degas_efficiency=0,
#           CH4_degas_efficiency=0.5, 
#           CO2_degas_efficiency=0.2)

# sys = System('mock_METAB', path=(H2E, S1, DM1, CH4E, S2, DM2),
#               recycle=(DM1-1, DM2-1))

# %%
# =============================================================================
# Preliminary analyses with mock METAB configuration
# =============================================================================

def create_systems(flowsheet_A=None, flowsheet_B=None, 
                   inf_concs={}, R1_init_conds={}, R2_init_conds={}):
    flowsheet_A = flowsheet_A or qs.Flowsheet('METAB_sysA')
    qs.main_flowsheet.set_flowsheet(flowsheet_A)
    
    ############# load components and set thermo #############
    pc.create_adm1_cmps()

    ############# create WasteStream objects #################
    inf_concs = inf_concs or default_inf_concs
    brewery_ww = WasteStream('BreweryWW_A', T=T1)
    brewery_ww.set_flow_by_concentration(Q, concentrations=inf_concs, units=('m3/d', 'kg/m3'))
    eff_A = WasteStream('Effluent_A', T=T2)
    bg1_A = WasteStream('biogas_1A', phase='g')
    bg2_A = WasteStream('biogas_2A', phase='g')
    
    ############# load process model ###########################
    adm1 = pc.ADM1()
    
    ############# sysA unit operation ########################
    R1_init_conds = R1_init_conds or default_R1_init_conds
    R2_init_conds = R2_init_conds or default_R2_init_conds
    
    H2E = AB('H2E', ins=[brewery_ww, 'return_1'], outs=('sidestream_1', ''), 
            split=(split_1, 1-split_1), V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
            retain_cmps=('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro'))
    DM1 = DM('DM1', ins=H2E-0, outs=(bg1_A, 1-H2E), tau=tau_1)
    CH4E = AB('CH4E', ins=[H2E-1, 'return_2'], outs=('sidestream_2', eff_A), 
            split=(split_2, 1-split_2), V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
            retain_cmps=('X_ac', 'X_h2'))
    DM2 = DM('DM2', ins=CH4E-0, outs=(bg2_A, 1-CH4E), tau=tau_2)
    H2E.set_init_conc(**R1_init_conds)
    CH4E.set_init_conc(**R2_init_conds)
    sysA = System('mock_METAB', 
                  path=(H2E, DM1, CH4E, DM2),
                  recycle=(DM1-1, DM2-1))
    sysA.set_dynamic_tracker(H2E, CH4E, bg1_A, bg2_A)
    
    #***************************************************
    flowsheet_B = flowsheet_A or qs.Flowsheet('METAB_sysB')
    qs.main_flowsheet.set_flowsheet(flowsheet_B)

    ############# sysB streams ########################
    inf_b = brewery_ww.copy('BreweryWW_B')
    eff_B = WasteStream('Effluent_B', T=T2)
    bg1_B = WasteStream('biogas_1B', phase='g')
    bg2_B = WasteStream('biogas_2B', phase='g')
    
    ############# sysB unit operation #################
    AnR1 = su.AnaerobicCSTR('AnR1', ins=inf_b, outs=(bg1_B, ''), 
                            V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1, 
                            retain_cmps=('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro'))
    AnR2 = su.AnaerobicCSTR('AnR2', ins=AnR1-1, outs=(bg2_B, eff_B), 
                            V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1,
                            retain_cmps=('X_ac', 'X_h2'))
    AnR1.set_init_conc(**R1_init_conds)
    AnR2.set_init_conc(**R2_init_conds)
    sysB = System('baseline', path=(AnR1, AnR2))
    sysB.set_dynamic_tracker(AnR1, AnR2, bg1_B, bg2_B)
    
    return sysA, sysB

#%%
@time_printer
def run(t, t_step, method=None, **kwargs):
    sysA, sysB = create_systems()
    print(f'Simulating {sysA.ID}...')
    sysA.simulate(state_reset_hook='reset_cache',
                  t_span=(0,t),
                  t_eval=np.arange(0, t+t_step, t_step),
                  method=method,
                  # export_state_to=ospath.join(folder, f'results/{method}_{t}d_sysA.xlsx'),
                  **kwargs)
    print(f'Simulating {sysB.ID}...')
    sysB.simulate(state_reset_hook='reset_cache',
                  t_span=(0,t),
                  t_eval=np.arange(0, t+t_step, t_step),
                  method=method,
                  # export_state_to=ospath.join(folder, f'results/{method}_{t}d_sysB.xlsx'),
                  **kwargs)

if __name__ == '__main__':
    t = 120
    t_step = 3
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}') # long live OCD!
    print(f'Time span 0-{t}d \n')
    run(t, t_step, method=method)
