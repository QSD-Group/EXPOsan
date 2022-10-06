# -*- coding: utf-8 -*-
'''
Created on Thu Aug  4 13:45:29 2022

@author: joy_c
'''

import numpy as np
from qsdsan import sanunits as su, processes as pc, WasteStream, System
from qsdsan.utils import time_printer, ospath
from chemicals.elements import molecular_weight as get_mw
from exposan.metab_mock import DegassingMembrane as DM, METAB_AnCSTR as AB

folder = ospath.dirname(__file__)

#%%
# =============================================================================
# Test ADM1 with mock METAB configuration
# =============================================================================

############# load components and set thermo #############
cmps = pc.create_adm1_cmps()

############# create WasteStream objects #################
Q = 5           # influent flowrate [m3/d]
T1 = 273.15+35
T2 = 273.15+25    # temperature [K]
C_mw = get_mw({'C':1})
N_mw = get_mw({'N':1})

brewery_ww = WasteStream('BreweryWW', T=T1)

conc = {
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

brewery_ww.set_flow_by_concentration(Q, concentrations=conc, units=('m3/d', 'kg/m3'))

eff = WasteStream('Effluent', T=T2)
# hs1 = WasteStream('headspace_1')
# hs2 = WasteStream('headspace_2')
# ss1 = WasteStream('sidestream_1')
# ss2 = WasteStream('sidestream_2')
# rs1 = WasteStream('return_1')
# rs2 = WasteStream('return_2')
bg1 = WasteStream('biogas_1', phase='g')
bg2 = WasteStream('biogas_2', phase='g')

############# load process model ###########################
# adm1 = pc.ADM1(path=os.path.join(adm_path, '_adm1.tsv'))
adm1 = pc.ADM1()

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

############# create unit operation ########################
# H2E = su.AnaerobicCSTR('H2E', ins=[brewery_ww, 'return_1'], outs=('headspace_1', ''), 
#                       V_liq=5, V_gas=0.556, T=T1, model=adm1, 
#                       retain_cmps=('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro'),
#                       pipe_resistance=5.0e4, fixed_headspace_P=True)
# S1 = su.Splitter('S1', ins=H2E-1, outs=('sidestream_1', ''),
#                  split=0.1, isdynamic=True)
# DM1 = DM('DM1', ins=S1-0, outs=(bg1, 1-H2E))

# CH4E = su.AnaerobicCSTR('CH4E', ins=[S1-1, 'return_2'], outs=('headspace_2', ''), 
#                        V_liq=75, V_gas=5, T=T2, model=adm1,
#                        retain_cmps=('X_ac', 'X_h2'),
#                        pipe_resistance=5.0e4, fixed_headspace_P=True)
# S2 = su.Splitter('S2', ins=CH4E-1, outs=('sidestream_2', eff),
#                  split=0.1, isdynamic=True)
# DM2 = DM('DM2', ins=S2-0, outs=(bg2, 1-CH4E))

H2E = AB('H2E', ins=[brewery_ww, 'return_1'], outs=('sidestream_1', ''), 
        split=(0.1, 0.9), V_liq=5, V_gas=0.556, T=T1, model=adm1, 
        retain_cmps=('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro'))
DM1 = DM('DM1', ins=H2E-0, outs=(bg1, 1-H2E))
CH4E = AB('CH4E', ins=[H2E-1, 'return_2'], outs=('sidestream_2', ''), 
        split=(0.1, 0.9), V_liq=75, V_gas=5, T=T2, model=adm1,
        retain_cmps=('X_ac', 'X_h2'))
DM2 = DM('DM2', ins=CH4E-0, outs=(bg2, 1-CH4E))

_ic1 = {
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

_ic2 = {
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

H2E.set_init_conc(**_ic1)
CH4E.set_init_conc(**_ic2)

# sys = System('mock_METAB', path=(H2E, S1, DM1, CH4E, S2, DM2),
#              recycle=(DM1-1, DM2-1))
sys = System('mock_METAB', path=(H2E, DM1, CH4E, DM2),
             recycle=(DM1-1, DM2-1))
sys.set_dynamic_tracker(H2E, CH4E)

__all__ = (
    'cmps', 'Q', 'T1', 'T2', 'conc',
    'brewery_ww', 'eff', 'bg1', 'bg2',
    'adm1', 'yields_bl', 'mus_bl', 'Ks_bl', 
    'H2E', 'CH4E', 
    '_ic1', '_ic2', 'sys', 
    )

#%%
@time_printer
def run(t, t_step, method=None, **kwargs):
    sys.simulate(state_reset_hook='reset_cache',
                 t_span=(0,t),
                 t_eval=np.arange(0, t+t_step, t_step),
                 method=method,
                 export_state_to=ospath.join(folder, f'results/sol_{t}d.xlsx'),
                 **kwargs)

if __name__ == '__main__':
    t = 12
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
