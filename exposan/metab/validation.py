# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 09:19:24 2023

@author: Joy Zhang <joycheung1994@gmail.com>
"""

import qsdsan as qs, matplotlib.pyplot as plt#, numpy as np
from qsdsan.utils import ospath
from qsdsan import processes as pc
from exposan.metab import flex_rhos_adm1, METAB_BatchExp, figures_path

#%%
cmps = pc.create_adm1_cmps()
adm1 = pc.ADM1(flex_rate_function=flex_rhos_adm1)
adm1.rate_function.params['rate_constants'][1:4] /= 16
# adm1.rate_function.params['rate_constants'][4:12] *= 0.8
adm1.rate_function.params['rate_constants'][12:19] *= 2

fy = 0.3
adm1.set_parameters(
    # f_bu_su=0.18, f_pro_su=0.35, f_ac_su=0.33, f_h2_su=0.14,
    # Y_su=0.1*fy,
    # Y_aa=0.08*fy,
    # Y_fa=0.06*fy,
    # Y_c4=0.06*fy,
    # Y_pro=0.04*fy,
    Y_ac=0.05*fy,
    # Y_h2=0.06*fy
    )

syn_ww = dict(
    # X_ch=2.064,
    # X_pr=0.118,
    # X_li=2.1e-3,
    # S_su=0.541,
    X_ch=2.239,
    X_pr=0.586,
    X_li=3.5e-3,
    S_su=3.461,
    S_aa=0.159,
    S_fa=2.755,
    S_IC=0.014,
    S_IN=0.014,
    S_cat=0.93,
    S_an=0.83
    )

#%%
# yields ~ 0.1g/L VSS
h2r_seed = dict(
    X_su=4.69/59,
    X_aa=1.23/59,
    X_fa=0.484/59,
    X_c4=1.42/59,
    X_pro=0.898/59
    )

bg1 = qs.WasteStream('bg1', phase='g')
BE1 = METAB_BatchExp('BE1', outs=(bg1,), V_liq=25, V_gas=35, V_beads=10, 
                    T=273.15+37, pH_ctrl=False, 
                    model=adm1, fixed_headspace_P=True)
# BE1.f_diff = 0.8
BE1.set_bulk_init(**syn_ww)
BE1.set_encap_init(**h2r_seed)

sys1 = qs.System('batch_1', path=(BE1,))
sys1.set_dynamic_tracker(BE1, bg1)

sys1.simulate(state_reset_hook='reset_cache', t_span=(0,8), method='BDF')
rh2 = bg1.scope.record[:, cmps.index('S_h2')] * bg1.scope.record[:, -1] # mg/L * m3/d = g(COD)/d
rh2 *= 1e3 # g(COD)/d -> mg(COD)/d
rh2 /= 1e6 # m3 reactor -> mL reactor
rh2 = rh2 * cmps.S_h2.i_mass / cmps.S_h2.chem_MW # mg/d -> mmol/d
t1 = bg1.scope.time_series
# bg1.scope.plot_time_series(('S_h2', 'S_IC', 'S_ch4', 'Q'))

#%%
tex1 = [0,1,2,3,4,5,7]
yex1 = [0, 0.18, 0.225, 0.145, 0.125, 0.13, 0.035]
yerr1 = [0, 0.012, 0.005, 0.08, 0.07, 0.09, 0.03]
fig, ax = plt.subplots()
ax.plot(t1, rh2, color='black')
ax.errorbar(tex1, yex1, yerr=yerr1, 
            marker='o', color='red', linestyle='dashed',
            capsize=2)
ax.set_xlabel('Days')
ax.set_ylabel('H2 production rate [mmol/d]')
ax.tick_params(axis='both', which='major', direction='inout')
fig.savefig(ospath.join(figures_path, 'validation_H2E_batch.png'),
            dpi=300, facecolor='white')

#%%
# yields ~ 1.0g/L VSS
ch4r_seed = dict(
    X_su=4.708e-1, 
    # X_su=0.45,
    X_aa=0.13, 
    X_fa=0.08, 
    # X_aa=4.838e-2,
    # X_fa=1.239e-1,
    # X_c4=1.423e-1, 
    X_c4=0.2,
    # X_pro=8.978e-2, 
    X_pro=0.1,
    # X_ac=2.959e-1, 
    # X_h2=1.467e-1, 
    X_ac=0.25,
    X_h2=0.15
    )

# ch4r_seed = {k:v*0.9 for k,v in ch4r_seed.items()}

bg2 = qs.WasteStream('bg2', phase='g')
BE2 = METAB_BatchExp('BE2', outs=(bg2,), V_liq=25, V_gas=35, V_beads=10, 
                    T=273.15+37, pH_ctrl=False, 
                    # max_encapsulation_tss=6,
                    model=adm1, fixed_headspace_P=True)
# BE2.f_diff = 0.8
BE2.set_bulk_init(**syn_ww)
BE2.set_encap_init(**ch4r_seed)

sys2 = qs.System('batch_2', path=(BE2,))
sys2.set_dynamic_tracker(BE2, bg2)

sys2.simulate(state_reset_hook='reset_cache', t_span=(0,8), method='BDF')
# bg2.scope.plot_time_series(('S_h2', 'S_IC', 'S_ch4', 'Q'))
rch4 = bg2.scope.record[:, cmps.index('S_ch4')] * bg2.scope.record[:, -1] # mg/L * m3/d = g(COD)/d
rch4 *= 1e3 # g(COD)/d -> mg(COD)/d
rch4 /= 1e6 # m3 reactor -> mL reactor
rch4 = rch4 * cmps.S_ch4.i_mass / cmps.S_ch4.chem_MW # mg/d -> mmol/d
t2 = bg2.scope.time_series

#%%
tex2 = [0,1,2,3,5,7,8]
yex2 = [0, 0.025, 0.042, 0.052, 0.065, 0.055, 0.062]
yerr2 = [0, 5e-3, 7.5e-3, 0.015, 0.02, 0.012, 0.02]
tex3 = [0,1,2,3,5,6,7]
yex3 = [0,0,0,0.019, 0.07, 0.075, 0.108]
yerr3 = [0,0,0, 6e-3, 6e-3, 0.01, 0.018]
fig, ax = plt.subplots()
ax.plot(t2, rch4, color='black')
ax.errorbar(tex2, yex2, yerr=yerr2, 
            marker='o', color='red', linestyle='dashed',
            capsize=2)
ax.errorbar(tex3, yex3, yerr=yerr3, 
            marker='^', color='blue', linestyle='dashed',
            capsize=2)
ax.set_xlabel('Days')
ax.set_ylabel('CH4 production rate [mmol/d]')
ax.tick_params(axis='both', which='major', direction='inout')
fig.savefig(ospath.join(figures_path, 'validation_CH4E_batch.png'),
            dpi=300, facecolor='white')