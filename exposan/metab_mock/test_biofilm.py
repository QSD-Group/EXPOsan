# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 15:53:18 2023

@author: joy_c
"""
import numpy as np
import qsdsan.processes as pc
from qsdsan.utils import auom, ospath

Q = 5 # m3/d

r_beads = 5e-3  # m; 5 mm
l_bl = 1e-5     # m, boundary layer thickness
n_dz = 10
dz = r_beads / n_dz
zs = np.linspace(dz, r_beads, n_dz)

cmps = pc.create_adm1_cmps()
S_idx = cmps.indices([i for i in cmps.IDs if i.startswith('S_')])
adm1 = pc.ADM1()
stoi_bk = adm1.stoichio_eval()
stoi_en = stoi_bk[:-3]  # no liquid-gas transfer
Rho_bk = adm1.rate_function

#%%
from qsdsan.processes._adm1 import (
    mass2mol_conversion,
    T_correction_factor,
    acid_base_rxn,
    substr_inhibit,
    Hill_inhibit,
    non_compet_inhibit
    )
from scipy.optimize import brenth

rhos = np.zeros(19) # 19 kinetic processes, no gas-liquid transfer
Cs = np.empty(19)
def Rho_en(state_arr, params=Rho_bk._params, T_op=298.15):
    ks = params['rate_constants']
    Ks = params['half_sat_coeffs']
    cmps = params['components']
    pH_ULs = params['pH_ULs']
    pH_LLs = params['pH_LLs']
    KS_IN = params['KS_IN']
    KI_nh3 = params['KI_nh3']
    KIs_h2 = params['KIs_h2']
    Kab = params['Ka_base']
    Ka_dH = params['Ka_dH']
    T_base = params['T_base']

    Cs[:8] = state_arr[12:20]
    Cs[8:12] = state_arr[19:23]
    Cs[12:] = state_arr[16:23]

    substrates = state_arr[:8]

    S_va, S_bu, S_h2, S_IN = state_arr[[3,4,7,10]]
    unit_conversion = mass2mol_conversion(cmps)
    cmps_in_M = state_arr[:27] * unit_conversion
    weak_acids = cmps_in_M[[24, 25, 10, 9, 6, 5, 4, 3]]

    Kas = Kab * T_correction_factor(T_base, T_op, Ka_dH)

    rhos[:] = ks * Cs
    Monod = substr_inhibit(substrates, Ks)
    rhos[4:12] *= Monod
    if S_va > 0: rhos[7] *= 1/(1+S_bu/S_va)
    if S_bu > 0: rhos[8] *= 1/(1+S_va/S_bu)

    try:
        h = brenth(acid_base_rxn, 1e-14, 1.0,
                args=(weak_acids, Kas),
                xtol=1e-12, maxiter=100)
    except:
        h = 10**(-6.5)
        
    nh3 = Kas[1] * weak_acids[2] / (Kas[1] + h)
    
    Iph = Hill_inhibit(h, pH_ULs, pH_LLs)
    Iin = substr_inhibit(S_IN, KS_IN)
    Ih2 = non_compet_inhibit(S_h2, KIs_h2)
    Inh3 = non_compet_inhibit(nh3, KI_nh3)
    rhos[4:12] *= Iph * Iin
    rhos[6:10] *= Ih2
    rhos[10] *= Inh3
    return rhos

n_cmps = len(cmps)
Diff = np.zeros(n_cmps)
Diff[:12] = [4.56e-6, 8.62e-6, 5.33e-6, 5.00e-6, 5.04e-6, 6.00e-6, 
             6.48e-6, 4.02e-4, 1.36e-4, 1.71e-4, 1.52e-4, 6.00e-6]
Diff[-3:-1] = 1.17e-4

f_diff = 0.75       # bead-to-water diffusivity fraction
D = Diff * f_diff   # Diffusivity in beads
k = Diff/l_bl       # mass transfer coeffient through liquid boundary layer

Cs_in = np.array([
    3.000e+03, 6.000e+02, 4.000e+02, 4.000e+02, 4.000e+02, 4.000e+02,
    4.000e+02, 5.000e-06, 5.000e-03, 4.804e+02, 1.401e+02, 2.000e+01,
    1.000e+02, 3.000e+02, 5.000e+02, 2.500e+02, 0.000e+00, 1.000e+00,
    1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 2.500e+01,
    4.000e+01, 2.000e+01, 9.900e+05
    ])*1e-3 # mg/L to kg/m3

_R = 8.3145e-2 # Universal gas constant, [bar/M/K]
_k_p = 5.0e4
_P_atm = 1.013
Pa_to_bar = auom('Pa').conversion_factor('bar')
p_vapor = lambda T: cmps.H2O.Psat(T)*Pa_to_bar
gas_mass2mol_conversion = (cmps.i_mass / cmps.chem_MW)[cmps.indices(('S_h2', 'S_ch4', 'S_IC'))]

def f_qgas(S_gas, T):
    p_gas = S_gas * _R * T
    P = sum(p_gas) + p_vapor(T) 
    q_gas = max(0, _k_p * (P - _P_atm))
    return q_gas

def dydt(t, y, HRT):
    Q, T = y[-2:]
    V_liq = HRT * Q
    V_bed = V_liq / 0.1
    V_beads = V_bed - V_liq
    V_gas = V_bed * 0.1
    A_beads = 3 * V_beads / r_beads # m2, total bead surface area

    S_gas = y[-5:-2]
    Cs_bk = y[-(n_cmps+5):-5]        # bulk liquid concentrations
    Cs_en = y[:n_dz*n_cmps].reshape((n_dz, n_cmps))  # each row is one control volume
    
    # Transformation
    rhos_en = np.apply_along_axis(Rho_en, 1, Cs_en, T_op=T)
    Rs_en = rhos_en @ stoi_en       # n_dz * n_cmps
    rhos_bk = Rho_bk(y[-(n_cmps+5):])
    Rs_bk = np.dot(stoi_bk.T, rhos_bk) # n_cmps+5
    q_gas = f_qgas(S_gas, T)
    gas_transfer = - q_gas*S_gas/V_gas + rhos_bk[-3:] * V_liq/V_gas * gas_mass2mol_conversion
    
    # Mass transfer (centered differences)
    dCdz = Cs_en * 0.
    d2Cdz2 = Cs_en * 0.
    dCdz[1:-1] = (Cs_en[2:] - Cs_en[:-2])/(2*dz)
    d2Cdz2[1:-1] = (Cs_en[2:] + Cs_en[:-2] - 2*Cs_en[1:-1])/(dz**2)
    d2Cdz2[0,:] = (Cs_en[1] - Cs_en[0])*(1+dz/(2*zs[0]))**2/dz  # inner boundary condition, forward difference
    C_lf = Cs_en[-1].copy()
    C_lf[S_idx] = (D*Cs_en[-2] + k*Cs_bk*dz)[S_idx]/(D+k*dz)[S_idx]  # Only solubles based on outer b.c.
    d2Cdz2[-1,:] = ((1+dz/zs[-1])**2*(C_lf-Cs_en[-2]) - (1-dz/zs[-1])**2*(Cs_en[-2]-Cs_en[-3]))/(dz**2)  # backward difference
    
    # Mass balance
    dCdt_en = D*(d2Cdz2 + np.diag(2/zs) @ dCdz) + Rs_en
    dCdt_bk = (Cs_in-Cs_bk)/HRT - A_beads/V_liq*k*(Cs_bk-C_lf) + Rs_bk
    
    dy = y*0.
    dy[:n_dz*n_cmps] = dCdt_en.flatten()
    dy[-(n_cmps+5):-5] = dCdt_bk
    dy[-5:-2] = gas_transfer
    
    return dy

#%%
from scipy.integrate import solve_ivp
import pandas as pd
from exposan.metab_mock import results_path, biomass_IDs, create_systems
# HRT = 1
# C0 = np.array([
#     1.204e-02, 5.323e-03, 9.959e-02, 1.084e-02, 1.411e-02, 1.664e-02,
#     4.592e-02, 2.409e-07, 7.665e-02, 5.693e-01, 1.830e-01, 3.212e-02,
#     2.424e-01, 2.948e-02, 4.766e-02, 2.603e-02, 4.708e+00, 1.239e+00,
#     4.838e-01, 1.423e+00, 8.978e-01, 2.959e+00, 1.467e+00, 4.924e-02,
#     4.000e-02, 2.000e-02, 9.900e+02
#     ])
# y0_bulk = np.array([
#     1.204e-02, 5.323e-03, 9.959e-02, 1.084e-02, 1.411e-02, 1.664e-02,
#     4.592e-02, 2.409e-07, 7.665e-02, 5.693e-01, 1.830e-01, 3.212e-02,
#     2.424e-01, 2.948e-02, 4.766e-02, 2.603e-02, 0, 0,
#     0, 0, 0, 0, 0, 0,   # no biomass in bulk liquid initially
#     4.000e-02, 2.000e-02, 9.900e+02, 3.922e-07, 2.228e-02, 1.732e-02,
#     Q, 298.15 
#     ])
# y0 = np.tile(C0, n_dz)
# y0 = np.append(y0, y0_bulk)

# sol = solve_ivp(dydt, t_span=(0, 100), y0=y0, method='BDF',
#                 t_eval=np.arange(0, 101, 5), args=(HRT,))

# loc = [f'bead_z{i}' for i in np.arange(n_dz)] + ['bulk']
# hdr = pd.MultiIndex.from_product([loc, cmps.IDs], names=['location', 'variable'])
# hdr = hdr.append(pd.MultiIndex.from_product([['bulk'], ['S_h2_gas', 'S_ch4_gas', 'S_IC_gas', 'Q', 'T']], 
#                                             names=['location', 'variable']))
# out = pd.DataFrame(sol.y.T, index=sol.t, columns=hdr)
# out.to_excel(ospath.join(results_path, f'sol_bfm_{HRT}d.xlsx'))

#%%
bm_idx = cmps.indices(biomass_IDs)

sys, = create_systems(which='C')
u, = sys.units
u.encapsulate_concentration = 1e5
inf, = sys.feeds
bg, eff = sys.products
_ic = u._concs[bm_idx].copy()
init_bm = np.linspace(1, 10, 10)

def suspended_vs_attached(HRT):
    print('\n')
    print('='*15)
    print(f'HRT: {HRT:.4f} d')
    u.V_liq = HRT * 5
    u.V_gas = HRT * 0.5
    sys.simulate(state_reset_hook='reset_cache', t_span=(0, 400), method='BDF')
    sus_rcod = 1-eff.COD/inf.COD
    i = 0
    while sus_rcod < 0.7:
        i += 1
        print('add VSS...')
        u._concs[bm_idx] *= 2
        sys.simulate(state_reset_hook='reset_cache', t_span=(0, 400), method='BDF')
        sus_rcod = 1-eff.COD/inf.COD
        if i > 10: return
    sus_bm = sum(u._state[bm_idx] * cmps.i_mass[bm_idx])
    y0_bulk = np.append(u._state, 298.15)
    y0_bulk[bm_idx] = 0.
    y0 = np.append(np.tile(u._state[:n_cmps], n_dz), y0_bulk)
    try: sol = solve_ivp(dydt, t_span=(0, 100), y0=y0, method='BDF', args=(HRT,))
    except: return sus_rcod, np.nan, sus_bm, np.nan
    bk_Cs = sol.y.T[-1][-(n_cmps+5):-5]
    en_Cs = sol.y.T[-2][-(n_cmps+5):-5]
    att_COD = sum(bk_Cs * cmps.i_COD * (1-cmps.g)) * 1e3
    att_rcod = 1-att_COD/inf.COD
    att_bm = sum(en_Cs[bm_idx] * cmps.i_mass[bm_idx]) * 1e3
    return sus_rcod, att_rcod, sus_bm, att_bm

#%%
record = []
rows = []
HRTs = [1, 0.5, 10/24, 8/24, 4/24, 2/24, 1/24]
for tau in HRTs:
    out = suspended_vs_attached(tau)
    if out: 
        record.append(out)
        rows.append(tau)

df = pd.DataFrame(record, index=rows, 
                  columns=pd.MultiIndex.from_product(
                      [['rCOD', 'biomass [mg TSS/L]'], ['Suspended', 'Encapsulated']]
                      ))
df.index.name = 'HRT [d]'
