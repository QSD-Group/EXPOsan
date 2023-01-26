# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 12:48:54 2023

@author: joy_c
"""

from qsdsan import processes as pc, sanunits as su, WasteStream, System
from qsdsan.utils import ExogenousDynamicVariable as EDV
from exposan.metab_mock import (
    rhos_adm1_ph_ctrl, 
    default_inf_concs,
    default_R1_init_conds,
    R1_ss_conds,
    biomass_IDs,
    )
import numpy as np

pc.create_adm1_cmps()
adm1 = pc.ADM1()
# dyn_params = adm1.rate_function.params.copy()
# adm1.set_rate_function(rhos_adm1_ph_ctrl)
# adm1.rate_function._params = dyn_params

Q = 5
T1 = 298.15
T2 = 298.15
ph1 = 5.8
ph2 = 7.2

Temp = EDV('T1', function=lambda t: T1)
Temp2 = EDV('T2', function=lambda t: T2)
pH1 = EDV('pH1', function=lambda t: ph1)
pH2 = EDV('pH2', function=lambda t: ph2)

ics = {
 'S_su': 0.32679646003805314,
 'S_aa': 0.3801610819236484,
 'S_fa': 12.437222319633748,
 'S_va': 0.3719673543175543,
 'S_bu': 0.47283246583627593,
 'S_pro': 0.3946420365926535,
 'S_ac': 10.182894473261367,
 'S_h2': 1.1655700001506622e-05,
 'S_ch4': 67.17348627201263,
 'S_IC': 846.4879614661522,
 'S_IN': 222.13725282096297,
 'S_I': 110.71467278942289,
 'X_c': 107.43132381172228,
 'X_ch': 1.2600235711799973,
 'X_pr': 1.3804329631122664,
 'X_li': 1.7696259648387357,
 'X_su': 732.9760678333023,
 'X_aa': 224.81751931525334,
 'X_fa': 126.7301174776879,
 'X_c4': 227.8726398428066,
 'X_pro': 140.2738127019708,
 'X_ac': 669.4626559278454,
 'X_h2': 245.67774602566578,
 'X_I': 206.42934561053158,
 'S_cat': 40.0,
 'S_an': 20.0,
 }

inf = WasteStream('inf')
eff = WasteStream('eff')
bg = WasteStream('bg', phase='g')
inf.set_flow_by_concentration(Q, concentrations=default_inf_concs, 
                                     units=('m3/d', 'kg/m3'))

R1 = su.AnaerobicCSTR('R1', ins=inf, outs=(bg, eff), 
                      # V_liq=50, V_gas=5, 
                      V_liq=20, V_gas=2, 
                      model=adm1, T=T1, 
                      retain_cmps=biomass_IDs, 
                      # exogenous_vars=(Temp1, pH1),
                      )
# R1.set_init_conc(**R1_ss_conds)
R1.set_init_conc(**ics)

sys = System('sys', path=(R1,))
sys.set_dynamic_tracker(R1, eff, bg)

#%%
HRTs = 2**np.linspace(-1, 4)
rCOD = {}
for i, tau in enumerate(HRTs):
    R1.V_liq = Q*tau
    R1.V_gas = R1.V_liq/(i+8)
    sys.simulate(state_reset_hook='reset_cache',
                 t_span=(0, 200),
                 method='BDF'
                 )
    assert np.all(R1._state >= 0)
    R1.scope.plot_time_series(('S_h2', 'S_ch4'))
    rCOD[tau] = 1-eff.COD/inf.COD