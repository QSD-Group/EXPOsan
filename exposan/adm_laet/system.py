# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from chemicals.elements import molecular_weight as get_mw
from qsdsan import sanunits as su, processes as pc, WasteStream, System
from qsdsan.utils import time_printer
# from exposan.adm import data_path


__all__ = (
    'create_system',
    'brewery_inf_kwargs',
    'default_inf_kwargs',
    'default_init_conds',
    )


# =============================================================================
# Parameters
# =============================================================================

Q = 170           # influent flowrate [m3/d]
Temp = 273.15+35    # temperature [K]
C_mw = get_mw({'C':1})
N_mw = get_mw({'N':1})

default_inf_kwargs = {
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
        'S_IC':0.04*C_mw,
        'S_IN':0.01*N_mw,
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
        'S_an':0.02,
        },
    'units': ('m3/d', 'kg/m3'),
    }
    
brewery_inf_kwargs = {
    'concentrations': {
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
        'S_an':0.02,
        },
    'units': ('m3/d', 'kg/m3'),
    }

default_init_conds = {
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
    'X_su': 0.3122*1e3,
    'X_aa': 0.9317*1e3,
    'X_fa': 0.3384*1e3,
    'X_c4': 0.3258*1e3,
    'X_pro': 0.1011*1e3,
    'X_ac': 0.6772*1e3,
    'X_h2': 0.2848*1e3,
    'X_I': 17.2162*1e3
    }


# %%

# =============================================================================
# Validation & Verification of ADM1
# =============================================================================

def create_system(flowsheet=None, inf_kwargs={}, adm_kwargs={}, init_conds={}):
    flowsheet = flowsheet or qs.Flowsheet('adm')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and streams
    pc.create_adm1_cmps()   
    inf = WasteStream('Influent', T=Temp)
    inf_kwargs = inf_kwargs or default_inf_kwargs
    inf.set_flow_by_concentration(Q, **inf_kwargs)
    
    eff = WasteStream('Effluent', T=Temp)
    bg = WasteStream('Biogas')
    
    # Process model
    adm1 = pc.ADM1(**adm_kwargs)
    
    # System setup
    AD = su.AnaerobicCSTR('AD', ins=inf, outs=(bg, eff), model=adm1, T=Temp)
    
    init_conds = init_conds or default_init_conds
    AD.set_init_conc(**init_conds)
    
    sys = System('ADM1_test', path=(AD,))
    sys.set_dynamic_tracker(AD, bg)
    
    return sys


# %%

@time_printer
def run(t, t_step, method=None, **kwargs):
    sys = create_system()
    # adm_kwargs = dict(path=os.path.join(data_path, '_adm1.tsv'))
    # sys = create_system(inf_kwargs=brewry_inf_kwargs, adm_kwargs=adm_kwargs)
    sys.simulate(state_reset_hook='reset_cache',
                 t_span=(0,t),
                 t_eval=np.arange(0, t+t_step, t_step),
                 method=method,
                 # export_state_to=f'results/sol2_{t}d_{method}_Phead.xlsx',
                 **kwargs)

if __name__ == '__main__':
    t = 200
    t_step = 5
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