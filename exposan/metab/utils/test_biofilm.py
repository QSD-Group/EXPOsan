# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 10:14:48 2023

@author: joy_c
"""
from qsdsan.utils import SanUnitScope
from exposan.metab import create_system, fermenters, methanogens
import numpy as np, pandas as pd

#%%
ecp = create_system(reactor_type='FB')
ue = ecp.flowsheet.unit
se = ecp.flowsheet.stream
Q = se.inf.F_vol * 24

#%%
ue.R1.voidage = 0.39
ue.R1.n_dz = 10
ue.R1._concs = np.tile(ue.R1._concs[:27], 11)
rcod = []
for db in (10, 2, 1, 0.2, 0.02):
    print(f'\n{db} mm\n{"="*10}')
    ue.R1.bead_diameter = db
    arr = []
    for tau in (1, 1/2, 10/24, 8/24, 4/24, 2/24, 1/24):
        print(f'tau = {tau:.3f} d')
        ue.R1.V_liq = tau * Q
        ue.R1.V_gas = tau * Q * 0.1
        ue.R1._prep_model()
        ue.R1._compile_ODE()
        ue.R1.scope = SanUnitScope(ue.R1)
        ecp.simulate(state_reset_hook='reset_cache', t_span=(0,400), method='BDF')
        arr.append(1-se.eff.COD/se.inf.COD)
    rcod.append(arr)

#%%
cmps = se.inf.components
ue.R1.n_dz = 5
ue.R1._concs = np.tile(ue.R1._concs[:27], 6)
ue.R1.bead_diameter = 0.02
ecp_rcod = []
srts = []
for void in np.linspace(0.39, 0.9, 10):
    ue.R1.voidage = void
    print(f'\nvoid {void}\n{"="*10}')
    for tau in (1, 1/2, 10/24, 8/24, 4/24, 2/24, 1/24):
        print(f'tau = {tau:.3f} d')
        ue.R1.V_liq = tau * Q
        ue.R1.V_gas = tau * Q * 0.1
        ue.R1._prep_model()
        ue.R1._compile_ODE()
        ue.R1.scope = SanUnitScope(ue.R1)
        ecp.simulate(state_reset_hook='reset_cache', t_span=(0,400), method='BDF')
        ecp_rcod.append(1-se.eff.COD/se.inf.COD)

        bk, en = ue.R1.biomass_tss([*fermenters, *methanogens])
        srt = (1 + (1-void)/void*en/bk) * tau
        srts.append(srt)

#%%
sys = create_system()
us = sys.flowsheet.unit
ss = sys.flowsheet.stream
allx_rcod = []
bm_rcod = []
hrts = np.tile((1, 1/2, 10/24, 8/24, 4/24, 2/24, 1/24), 10)
for srt, tau in zip(srts, hrts):
    print(f'{srt}, {tau}\n')
    us.R1.V_liq = tau * Q
    us.R1.V_gas = tau * Q * 0.1
    us.R1._f_retain[:] = (1-tau/srt) * cmps.x
    us.R1._ODE = None
    sys.simulate(state_reset_hook='reset_cache', t_span=(0,400), method='BDF')
    allx_rcod.append(1-ss.eff.COD/ss.inf.COD)
    
    us.R1._f_retain[cmps.indices([*fermenters, *methanogens])] = 1-tau/srt
    us.R1._ODE = None
    sys.simulate(state_reset_hook='reset_cache', t_span=(0,400), method='BDF')
    bm_rcod.append(1-ss.eff.COD/ss.inf.COD)       
