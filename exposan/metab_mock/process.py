# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from qsdsan.processes._adm1 import (
    R,
    mass2mol_conversion,
    T_correction_factor,
    acid_base_rxn,
    substr_inhibit,
    Hill_inhibit,
    non_compet_inhibit
    )

__all__ = ('rhos_adm1_ph_ctrl',)


#%%
rhos = np.zeros(22) # 22 kinetic processes
Cs = np.empty(19)

def rhos_adm1_ph_ctrl(state_arr, params):
    ks = params['rate_constants']
    Ks = params['half_sat_coeffs']
    cmps = params['components']
    pH_ULs = params['pH_ULs']
    pH_LLs = params['pH_LLs']
    KS_IN = params['KS_IN']
    KI_nh3 = params['KI_nh3']
    KIs_h2 = params['KIs_h2']
    KHb = params['K_H_base']
    Kab = params['Ka_base']
    KH_dH = params['K_H_dH']
    Ka_dH = params['Ka_dH']
    kLa = params['kLa']
    T_base = params['T_base']
    root = params['root']
    

    Cs[:8] = state_arr[12:20]
    Cs[8:12] = state_arr[19:23]
    Cs[12:] = state_arr[16:23]
    substrates = state_arr[:8]
    S_va, S_bu, S_h2, S_IN = state_arr[[3,4,7,10]]
    unit_conversion = mass2mol_conversion(cmps)
    cmps_in_M = state_arr[:27] * unit_conversion
    weak_acids = cmps_in_M[[24, 25, 10, 9, 6, 5, 4, 3]]

    T_op, pH = state_arr[-2:]   #!!! change in state_arr
    biogas_S = state_arr[7:10].copy()
    biogas_p = R * T_op * state_arr[27:30]
    Kas = Kab * T_correction_factor(T_base, T_op, Ka_dH)
    KH = KHb * T_correction_factor(T_base, T_op, KH_dH) / unit_conversion[7:10]

    rhos[:-3] = ks * Cs
    Monod = substr_inhibit(substrates, Ks)
    rhos[4:12] *= Monod
    if S_va > 0: rhos[7] *= 1/(1+S_bu/S_va)
    if S_bu > 0: rhos[8] *= 1/(1+S_va/S_bu)

    h = 10**(-pH)
    delta = acid_base_rxn(h, weak_acids, Kas)
    S_cat = weak_acids[0] - delta

    nh3 = Kas[1] * weak_acids[2] / (Kas[1] + h)
    co2 = weak_acids[3] - Kas[2] * weak_acids[3] / (Kas[2] + h)
    biogas_S[-1] = co2 / unit_conversion[9]
    
    Iph = Hill_inhibit(h, pH_ULs, pH_LLs)
    Iin = substr_inhibit(S_IN, KS_IN)
    Ih2 = non_compet_inhibit(S_h2, KIs_h2)
    Inh3 = non_compet_inhibit(nh3, KI_nh3)
    rhos[4:12] *= Iph * Iin
    rhos[6:10] *= Ih2
    rhos[10] *= Inh3
    rhos[-3:] = kLa * (biogas_S - KH * biogas_p)
    root.data = {
        'pH':pH, 
        'Iph':Iph, 
        'Ih2':Ih2, 
        'Iin':Iin, 
        'Inh3':Inh3,
        'Monod':Monod,
        'rhos':rhos[4:12].copy(),
        'gas_transfer':rhos[-3:].copy(),
        'S_cat':S_cat
        }
    return rhos