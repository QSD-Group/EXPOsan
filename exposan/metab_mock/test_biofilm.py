# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 15:53:18 2023

@author: joy_c
"""
import numpy as np
import qsdsan.processes as pc
from qsdsan.utils import auom

#%%
HRT = 1/24 # d
Q = 5 # m3/d
V_liq = HRT * Q
V_bed = V_liq / 0.1
V_beads = V_bed - V_liq
V_gas = V_bed * 0.1

r_beads = 5e-3  # m; 5 mm
A_beads = 3 * V_beads / r_beads # m2, total bead surface area
l_bl = 1e-5     # m, boundary layer thickness

cmps = pc.create_adm1_cmps()
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

    rhos[:-3] = ks * Cs
    Monod = substr_inhibit(substrates, Ks)
    rhos[4:12] *= Monod
    if S_va > 0: rhos[7] *= 1/(1+S_bu/S_va)
    if S_bu > 0: rhos[8] *= 1/(1+S_va/S_bu)

    h = brenth(acid_base_rxn, 1e-14, 1.0,
            args=(weak_acids, Kas),
            xtol=1e-12, maxiter=100)

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

y0_bulk = np.ones(31)
n_dz = 10
dz = r_beads / n_dz
zs = np.linspace(dz, r_beads, n_dz)
y0 = np.tile(y0_bulk, n_dz+1)

_R = 8.3145e-2 # Universal gas constant, [bar/M/K]
_k_p = 5.0e4
_P_atm = 1.013
Pa_to_bar = auom('Pa').conversion_factor('bar')
p_vapor = lambda T: cmps.H2O.Psat(T)*Pa_to_bar
gas_mass2mol_conversion = (cmps.i_mass / cmps.chem_MW)[cmps.indices(('S_h2', 'S_ch4', 'S_IC'))]

def f_qgas(rhoTs, S_gas, T):
    p_gas = S_gas * _R * T
    P = sum(p_gas) + p_vapor(T) 
    q_gas = max(0, _k_p * (P - _P_atm))
    return q_gas

def dydt(t, y):
    Q, T = y[-2:]
    S_gas = y[-5:-2]
    Cs = y[:-5].reshape(n_dz, n_cmps)  # each row is one control volume, last row is bulk liquid
    
    # Transformation
    rhos_en = np.apply_along_axis(Rho_en, 1, Cs[:-1], T=T)
    Rs_en = rhos_en @ stoi_en       # n_dz * n_cmps
    rhos_bk = Rho_bk(y[-(n_cmps+5):])
    Rs_bk = np.dot(stoi_bk.T, rhos_bk) # n_cmps+5
    q_gas = f_qgas(rhos_bk[-3:], S_gas, T)
    gas_transfer = - q_gas*S_gas/V_gas + rhos_bk[-3:] * V_liq/V_gas * gas_mass2mol_conversion
    
    # Mass transfer
    
