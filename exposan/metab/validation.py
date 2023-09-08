# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 09:19:24 2023

@author: Joy Zhang <joycheung1994@gmail.com>
"""

import qsdsan as qs, matplotlib.pyplot as plt, numpy as np
from qsdsan.utils import ospath
from qsdsan import processes as pc
from exposan.metab import flex_rhos_adm1, METAB_BatchExp, figures_path

#%%
cmps = pc.create_adm1_cmps()
adm1 = pc.ADM1(flex_rate_function=flex_rhos_adm1)
# adm1.rate_function.params['rate_constants'][1:4] /= 16
# adm1.rate_function.params['rate_constants'][4:12] *= 1.5
# adm1.rate_function.params['half_sat_coeffs'][:-2] *= 0.5
# adm1.rate_function.params['rate_constants'][10:12] *= 1.5
# adm1.rate_function.params['half_sat_coeffs'][-2:] *= 0.5
# adm1.rate_function.params['rate_constants'][12:19] *= 2

# adm1.set_parameters(
#     # f_bu_aa=0.18, f_pro_aa=0.35, f_ac_aa=0.33, f_va_aa=0.14
#     # 'f_va_aa': 0.23,
#     # 'f_bu_aa': 0.26,
#     # 'f_pro_aa': 0.05,
#     # 'f_ac_aa': 0.4,
#     # 'f_h2_aa': 0.06,
#     # Y_aa=0.04,
#     # Y_ac=0.025
#     )

syn_ww = dict(
    X_ch=2.065,
    X_pr=5.832,
    S_aa=0.159,
    S_fa=2.755,
    S_IC=0.0143,
    S_IN=0.014,
    S_cat=0.0075,
    S_an=0.0058
    )

#%%
from exposan.metab import create_system
sys = create_system(1, inf_concs=syn_ww, tot_HRT=8)
# mdl = sys.units[0].model
# mdl.rate_function.params['rate_constants'][1:4] /= 16
# mdl.rate_function.params['rate_constants'][4:12] *= 1.5
# mdl.rate_function.params['half_sat_coeffs'][:-2] *= 0.5
sys.simulate(state_reset_hook='reset_cache', t_span=(0,200), method='BDF')

#%%
# yields ~ 0.1g/L VSS
# h2r_seed = dict(
#     X_su=4.69/59,
#     X_aa=1.23/59,
#     X_fa=0.484/59,
#     X_c4=1.42/59,
#     X_pro=0.898/59
#     )

h2r_seed = {
    'X_su': 1.956e-2,
    'X_aa': 4.229e-2,
    'X_fa': 1.546e-2,
    'X_c4': 1.561e-2,
    'X_pro': 0.511e-2,
    # 'X_su': 1.1795398908390478/100,
    # 'X_aa': 2.5184674870545223/100,
    # 'X_fa': 0.9329873586793094/100,
    # 'X_c4': 0.9313043558027272/100,
    # 'X_pro': 0.3060024744735877/100,
    }

bg1 = qs.WasteStream('bg1', phase='g')
BE1 = METAB_BatchExp('BE1', outs=(bg1,), V_liq=25, V_gas=70, V_beads=10, 
                    T=273.15+37, pH_ctrl=False, 
                    model=adm1, fixed_headspace_P=True)
BE1.f_diff = 0.8
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
dt1 = t1[1:] - t1[:-1]
ch2 = np.cumsum(rh2[1:]*dt1)
ch2 = np.insert(ch2, 0, 0)

#%%
tex1 = [0,1,2,3,4,5,7,9]
# yex1 = [0, 1.5541e-1, 3.2176e-2, 0, 0, 2.6512e-3, 0]   # mmol/mg protein/d
yex1 = [0, 0.1554, 0.1876, 0.1226, 0.1081, 0.1108, 0.0305, 0.0738]
# yerr1 = [0, 1.2e-1, 1.33e-1, 0,0,0,0]
yerr1 = [0, 0.01231, 0.00414, 0.06947, 0.06282, 0.07872, 0.02465, 0.07875]
yex1 = [x*0.884 for x in yex1]    # 88.4 mg protein/L * 10 mL = 0.884 mg
yerr1 = [x*0.884 for x in yerr1]
fig, ax = plt.subplots()
# ax.plot(t1, rh2, color='black', label='simulated')
ax.plot(t1, ch2, color='black', label='simulated')
ax.errorbar(tex1, yex1, yerr=yerr1, label='PEG',
            marker='o', color='red', linestyle='dashed',
            capsize=2)
ax.legend()
ax.set_xlabel('Days')
ax.set_ylabel('Cumulative H2 production [mmol]')
ax.tick_params(axis='both', which='major', direction='inout')
fig.savefig(ospath.join(figures_path, 'validation_H2E_batch.png'),
            dpi=300, facecolor='white')

#%%
# yields ~ 1.0g/L VSS
# ch4r_seed = dict(
#     X_su=4.708e-1, 
#     # X_aa=4.838e-2,
#     # X_fa=1.239e-1,
#     X_c4=1.423e-1, 
#     X_pro=8.978e-2, 
#     # X_ac=2.959e-1, 
#     # X_h2=1.467e-1, 
#     # X_su=0.5,
#     X_aa=0.13, 
#     X_fa=0.08, 
#     # X_c4=0.2,
#     # X_pro=0.1,
#     X_ac=0.25,
#     X_h2=0.15
#     )

ch4r_seed = {
    'X_su': 1.956004638872764,
    'X_aa': 4.229473465645631,
    'X_fa': 1.5466419486676024,
    'X_c4': 1.5615615740437239,
    'X_pro': 0.511117652665559,
    'X_ac': 3.2327967588577007,
    'X_h2': 1.390562414197235,
    # 'X_su': 1.1795398908390478,
    # 'X_aa': 2.5184674870545223,
    # 'X_fa': 0.9329873586793094,
    # 'X_c4': 0.9313043558027272,
    # 'X_pro': 0.3060024744735877,
    # 'X_ac': 1.9359643800709951,
    # 'X_h2': 0.8342225195583577
    }



ch4r_seed = {k: v/10 for k,v in ch4r_seed.items()}

bg2 = qs.WasteStream('bg2', phase='g')
BE2 = METAB_BatchExp('BE2', outs=(bg2,), V_liq=25, V_gas=35, V_beads=10, 
                    T=273.15+37, pH_ctrl=False, 
                    # max_encapsulation_tss=6,
                    model=adm1, fixed_headspace_P=True)
BE2.f_diff = 0.8
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
dt2 = t2[1:] - t2[:-1]
cch4 = np.cumsum(rch4[1:]*dt2)
cch4 = np.insert(cch4, 0, 0)


#%%
# PEG, suspended AS
tex2 = [0,1,2,3,5,7,8]
yex2 = [0, 0.026, 0.043, 0.052, 0.067, 0.056, 0.065]
yerr2 = [0, 3.6e-3, 7.2e-3, 0.0156, 0.02, 0.013, 0.02]

# PEG, granular AS
tex4 = [0,1,2,3,4,7,8]
yex4 = [0, 5.75e-3, 9.41e-3, 14.9e-3, 13.8e-3, 16.3e-3, 18.7e-3]
yex4 = [x*13 for x in yex4]
yerr4 = [0, 1.44e-3, 1.58e-3, 1.14e-3, 1.21e-3, 4.71e-4, 2.59e-3]
yerr4 = [x*13 for x in yerr4]

# PEG, PAC-supported biofilm
tex5 = [0,1,2,4,6]
yex5 = [0, 1.03e-2, 2.53e-2, 8.4e-2, 9.67e-2]
yex5 = [x*17 for x in yex5]
yerr5 = [0, 3.1e-4, 3.46e-3, 6.47e-3, 7.76e-3]
yerr5 = [x*17 for x in yerr5]

# PEG, PAC-supported biofilm, replicate
tex6 = [0,1,2,3,4,7,8]
yex6 = [0, 0.022, 0.045, 0.095, 0.121, 0.140, 0.153]
yex6 = [x*11 for x in yex6]
yerr6 = [0, 3.5e-4, 7.38e-3, 4.91e-3, 3.91e-3, 8.56e-3, 2.1e-3]
yerr6 = [x*11 for x in yerr6]

fig, ax = plt.subplots()
# ax.plot(t2, rch4, color='black', label='simulated')
ax.plot(t2, cch4, color='black', label='simulated')
ax.errorbar(tex2, yex2, yerr=yerr2, label='PEG, suspended AS',
            marker='s', linestyle='dashed', capsize=2)
ax.errorbar(tex4, yex4, yerr=yerr4, label='PEG, granular AS',
            marker='D', linestyle='dashed', capsize=2)
ax.errorbar(tex5, yex5, yerr=yerr5, label='PEG, PAC-supported biofilm',
            marker='^', linestyle='dashed', capsize=2)
ax.errorbar(tex6, yex6, yerr=yerr6, label='PEG, biofilm, replicate',
            marker='o', linestyle='dashed', capsize=2)

ax.legend(loc='upper left')
# ax.set_ylim(-0.1, 2.6)
ax.set_xlabel('Days')
ax.set_ylabel('Cumulative CH4 production [mmol]')
ax.tick_params(axis='both', which='major', direction='inout')
fig.savefig(ospath.join(figures_path, 'validation_CH4E_batch.png'),
            dpi=300, facecolor='white')