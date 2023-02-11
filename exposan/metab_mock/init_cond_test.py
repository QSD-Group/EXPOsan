# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 16:55:35 2023

@author: joy_c
"""

from exposan.metab_mock import create_systems, biomass_IDs
import numpy as np

sys, = create_systems(which='C')
u = sys.units[0]
bg, eff = u.outs
inf, = u.ins
idx = inf.components.indices(biomass_IDs)
_ic = u._concs[idx].copy()
init_bm = np.linspace(1, 10, 10)

u.V_liq = 2.5
u.V_gas = 0.25

rcod = []
bm = []
bm_dist = []
fug_ch4 = []

for i in init_bm:
    u._concs[idx] = _ic*i
    sys.simulate(state_reset_hook='reset_cache', t_span=(0, 400), method='BDF')
    rcod.append(1-eff.COD/inf.COD)
    bm.append(sum(u._state[idx]))
    bm_dist.append(u._state[idx]/sum(u._state[idx])*100)
    fug_ch4.append(eff.iconc['S_ch4'])
