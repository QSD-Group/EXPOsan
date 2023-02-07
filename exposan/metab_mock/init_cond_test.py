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

u.V_liq = 5
u.V_gas = 0.556
idx = inf.components.indices(biomass_IDs)

_ic = u._concs[idx].copy()
rcod = []
bm = []
fug_ch4 = []

init_bm = np.linspace(11, 20, 10)
for i in init_bm:
    u._concs[idx] = _ic*i
    sys.simulate(state_reset_hook='reset_cache', t_span=(0, 400), method='BDF')
    rcod.append(1-eff.COD/inf.COD)
    bm.append(sum(u._state[idx]))
    fug_ch4.append(eff.iconc['S_ch4'])
