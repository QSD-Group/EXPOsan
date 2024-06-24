# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from chemicals.elements import molecular_weight as get_mw
from qsdsan import sanunits as su, processes as pc, WasteStream, System, get_thermo
from qsdsan.utils import load_data, ospath, time_printer
from exposan.bsm2 import data_path


# =============================================================================
# Parameters
# =============================================================================

Q = 190           # influent flowrate [m3/d]
HRT = 20
V_liq = Q*HRT
V_gas = 0.088*V_liq
Temp = 273.15+35    # temperature [K]
C_mw = get_mw({'C':1})
N_mw = get_mw({'N':1})
P_mw = get_mw({'P':1})
adm1init = load_data(ospath.join(data_path, 'adm1init.csv'), index_col=0).to_dict('index')

# Table 1.1 [mg/L]
inf_asm2d = dict(
    S_O2=0,
    S_F=26.44,
    S_A=17.66,
    S_I=27.23,
    S_NH4=18.58,
    S_N2=5.07,
    S_NO3=0.02,
    S_PO4=4.69,
    S_IC=78.99,
    X_I=10964.41,
    X_S=19084.76,
    X_H=9479.39,
    X_PAO=3862.20,
    X_PP=450.87,
    X_PHA=24.64,
    X_AUT=333.79,
    S_K=19.79,
    S_Mg=189.87,
    S_Na=70,
    S_Cl=1035,
    S_Ca=300,
    )
    
# Table 1.3 [kg/m3]
inf_adm1p = dict(
    S_su=0.018,
    S_aa=0.008,
    S_ac=0.018,
    S_IC=0.021*C_mw,
    S_IN=0.036*N_mw,
    S_IP=0.006*P_mw,
    S_I=0.027,
    X_ch=8.020,
    X_pr=8.481,
    X_li=11.416,
    X_I=11.946,
    X_PHA=0.025,
    X_PP=0.015*P_mw,
    X_PAO=3.862,
    S_K=0.001*39,
    S_Mg=0.008*24.3,
    S_Ca=0.007*40,
    S_Na=0.003*23,
    S_Cl=0.029*35.5,
    # S_N2=0.0004*14
    )

# Table 1.4 [kg/m3]
out_adm1p = dict(
    S_su=0.013,
    S_aa=0.006,
    S_fa=0.116,
    S_va=0.012,
    S_bu=0.016,
    S_pro=0.019,
    S_ac=0.055,
    S_h2=2.65e-7,
    S_ch4=0.052,
    S_IC=0.059*C_mw,
    S_IN=0.080*N_mw,
    S_IP=0.007*P_mw,
    S_I=0.027,
    X_ch=1.441,
    X_pr=1.513,
    X_li=2.025,
    X_I=12.345,
    X_PHA=0.252,
    X_PP=8.05e-6*P_mw,
    # X_biomass=3.600,
    X_su=3.600,
    S_K=0.005*39,
    S_Mg=0.001*24.3,
    S_Ca=0.001*40,
    X_ACP=0.002*310.176722,
    X_struv=0.011*245.406502,
    S_Na=0.003*23,
    S_Cl=0.029*35.5,
    # S_N2=0.0004*14
    )

# Table 1.5 [mg/L]
out_asm2d = dict(
    S_NH4=1291.68,
    S_PO4=298.09,
    S_F=134.43,
    S_A=353.82,
    S_I=27.23,
    S_IC=885.27,
    S_K=208.84,
    S_Mg=28.29,
    X_I=12704.93,
    X_S=8218.94,
    S_Na=70,
    S_Cl=1035,
    S_Ca=20.45,
    X_ACP=722.17,
    X_struv=1578.52
    )

# default_init_conds = {
#     'S_su': 0.0124*1e3,
#     'S_aa': 0.0055*1e3,
#     'S_fa': 0.1074*1e3,
#     'S_va': 0.0123*1e3,
#     'S_bu': 0.0140*1e3,
#     'S_pro': 0.0176*1e3,
#     'S_ac': 0.0893*1e3,
#     'S_h2': 2.5055e-7*1e3,
#     'S_ch4': 0.0555*1e3,
#     'S_IC': 0.0951*C_mw*1e3,
#     'S_IN': 0.0945*N_mw*1e3,
#     'S_I': 0.1309*1e3,
#     'X_ch': 0.0205*1e3,
#     'X_pr': 0.0842*1e3,
#     'X_li': 0.0436*1e3,
#     'X_su': 0.3122*1e3,
#     'X_aa': 0.9317*1e3,
#     'X_fa': 0.3384*1e3,
#     'X_c4': 0.3258*1e3,
#     'X_pro': 0.1011*1e3,
#     'X_ac': 0.6772*1e3,
#     'X_h2': 0.2848*1e3,
#     'X_I': 17.2162*1e3
#     }

default_init_conds = {
    'S_su': 0.014*1e3,
    'S_aa': 0.0062*1e3,
    'S_fa': 0.126*1e3,
    'S_va': 0.0129*1e3,
    'S_bu': 0.0168*1e3,
    'S_pro': 0.0204*1e3,
    'S_ac': 0.0588*1e3,
    'S_h2': 2.8309e-7*1e3,
    'S_ch4': 0.0544*1e3,
    'S_IC': 0.089*12*1e3,
    'S_IN': 0.0663*14*1e3,
    'S_IP': 0.028*31*1e3,
    'S_I': 0.1309*1e3,
    'X_ch': 1.302*1e3,
    'X_pr': 1.3613*1e3,
    'X_li': 1.8127*1e3,
    'X_su': 0.5146*1e3,
    'X_aa': 0.4017*1e3,
    'X_fa': 0.3749*1e3,
    'X_c4': 0.1596*1e3,
    'X_pro': 0.0896*1e3,
    'X_ac': 0.5006*1e3,
    'X_h2': 0.258*1e3,
    'X_I': 12.9232*1e3,
    'X_PHA': 0.6697*1e3,
    'X_PAO': 0.9154*1e3,
    'S_K': 0.0129*1e3,
    'S_Mg': 0.0001*1e3,
    'S_Ca': 2e-4*1e3,
    'X_struv':0.0161*1e3,
    'X_ACP': 9e-4*1e3,
    'X_FePO4': 0.001*1e3,
    'S_Na': 0.061*1e3,
    'S_Cl': 0.0126*1e3
    }

#%%

cmps_asm = pc.create_masm2d_cmps()
inf_asm = WasteStream('inf_asm', T=Temp)
inf_asm.set_flow_by_concentration(
    flow_tot=Q, 
    concentrations=inf_asm2d,
    units=('m3/d', 'mg/L')
    )
asm = pc.mASM2d()
thermo_asm = get_thermo()
cmps_adm = pc.create_adm1p_cmps()
alt_inf_adm = WasteStream('alt_inf_adm', T=Temp)
alt_inf_adm.set_flow_by_concentration(
    flow_tot=Q,
    concentrations=inf_adm1p,
    units=('m3/d', 'kg/m3')
    )
alt_eff_adm = WasteStream('alt_eff_adm', T=Temp)
alt_eff_adm.set_flow_by_concentration(
    flow_tot=Q, 
    concentrations=out_adm1p,
    units=('m3/d', 'kg/m3')
    )
adm = pc.ADM1p(
    f_bu_su=0.1328, f_pro_su=0.2691, f_ac_su= 0.4076,
    q_ch_hyd=0.3, q_pr_hyd=0.3, q_li_hyd=0.3, 
    )
thermo_adm = get_thermo()

J1 = su.mASM2dtoADM1p('J1', upstream=inf_asm, downstream='inf_adm', 
                      thermo=thermo_adm, isdynamic=True, 
                      adm1_model=adm, asm2d_model=asm)
J1.xs_to_li = 0.6
AD = su.AnaerobicCSTR('AD', ins=alt_inf_adm, outs=('biogas', 'eff_adm'), isdynamic=True, 
# AD = su.AnaerobicCSTR('AD', ins=J1-0, outs=('biogas', 'eff_adm'), isdynamic=True, 
                      V_liq=V_liq, V_gas=V_gas, T=Temp, model=adm)
# AD.algebraic_h2 = True
AD.algebraic_h2 = False
# AD.set_init_conc(**adm1init['AD1'])
AD.set_init_conc(**default_init_conds)
J2 = su.ADM1ptomASM2d('J2', upstream=AD-1, downstream='eff_asm', 
# J2 = su.ADM1ptomASM2d('J2', upstream=alt_eff_adm, downstream='eff_asm', 
                      thermo=thermo_asm, isdynamic=True, 
                      adm1_model=adm, asm2d_model=asm)

sys = System(path=(J1, AD, J2))
fs = sys.flowsheet.stream
sys.set_dynamic_tracker(AD, fs.biogas, fs.eff_adm)

#%%
@time_printer
def run():
    sys.simulate(state_reset_hook='reset_cache', t_span=(0, 400), method='BDF')

run()

# fs.inf_adm.conc vs. inf_adm1p
# fs.eff_asm.conc vs. out_asm2d