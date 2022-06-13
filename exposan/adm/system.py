# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np
from qsdsan import sanunits as su, processes as pc, set_thermo, WasteStream, System
from qsdsan.utils import time_printer

adm_path = os.path.dirname(__file__)

# =============================================================================
# Validation & Verification of ADM1
# =============================================================================

############# load components and set thermo #############
cmps = pc.load_adm1_cmps()
set_thermo(cmps)

############# create WasteStream objects #################
Q = 170           # influent flowrate [m3/d]
Temp = 273.15+35    # temperature [K]

inf = WasteStream('Influent', T=Temp)
inf_kwargs = {
    'concentrations': {
        'S_su':0.01,
        'S_aa':1e-3,
        'S_fa':1e-3,
        'S_va':1e-3,
        'S_bu':1e-3,
        'S_pro':1e-3,
        'S_ac':1e-3,
        'S_h2':1e-8,
        'S_ch4':1e-5,
        'S_IC':0.04*12,
        'S_IN':0.01*14,
        'S_I':0.02,
        'X_c':2.0,
        'X_ch':5.0,
        'X_pr':20.0,
        'X_li':5.0,
        'X_aa':1e-2,
        'X_fa':1e-2,
        'X_c4':1e-2,
        'X_pro':1e-2, 
        'X_ac':1e-2, 
        'X_h2':1e-2, 
        'X_I':25, 
        'S_cat':0.04, 
        'S_an':0.02
        },
    'units': ('m3/d', 'kg/m3'),
    }
inf.set_flow_by_concentration(Q, **inf_kwargs)

eff = WasteStream('Effluent', T=Temp)
bg = WasteStream('Biogas')

############# load process model ###########################
adm1 = pc.ADM1()

############# create unit operation ########################
AD = su.AnaerobicCSTR('AD', ins=inf, outs=(bg, eff), model=adm1)
AD.set_init_conc(
    S_su =  0.0124*1e3,
    S_aa =  0.0055*1e3,
    S_fa =  0.1074*1e3,
    S_va =  0.0123*1e3,
    S_bu = 0.0140*1e3,
    S_pro = 0.0176*1e3,
    S_ac =  0.0893*1e3,
    S_h2 =  2.5055e-7*1e3,
    S_ch4 = 0.0555*1e3,
    S_IC = 0.0951*1e3,
    S_IN = 0.0945*1e3,
    S_I = 0.1309*1e3,
    X_c = 0.0*1e3,
    X_ch = 0.0205*1e3,
    X_pr = 0.0842*1e3,
    X_li = 0.0436*1e3,
    X_su = 0.3122*1e3,
    X_aa = 0.9317*1e3,
    X_fa = 0.3384*1e3,
    X_c4 = 0.3258*1e3,
    X_pro = 0.1011*1e3,
    X_ac = 0.6772*1e3,
    X_h2 =  0.2848*1e3,
    X_I =  17.2162*1e3
    )

sys = System('ADM1_test', path=(AD,))
sys.set_dynamic_tracker(AD, eff, bg)

#%%
@time_printer
def run(t, t_step, method=None, **kwargs):
    sys.simulate(state_reset_hook='reset_cache',
                 t_span=(0,t),
                 t_eval=np.arange(0, t+t_step, t_step),
                 method=method,
                 export_state_to=f'results/sol2_{t}d_{method}.xlsx',
                 **kwargs)

if __name__ == '__main__':
    t = 50
    t_step = 1
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
