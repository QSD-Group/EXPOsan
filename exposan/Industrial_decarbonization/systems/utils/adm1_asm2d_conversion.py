# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 13:52:16 2023

@author: joy_c
"""
import numpy as np
from warnings import warn
from math import isclose
import qsdsan.processes as pc

__all__ = (
    'industrial_wws_adm1',
    'adm1_2_asm2d',
       )

#%%
# adm1 = [S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN,
#         S_I, X_c, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, 
#         S_cat, S_an]

# Without pretreatment 

bm2_ind = np.array([180, 180, 180, 150, 150, 150, 150,   0.,  0.,  1, 0,
                    24, 120, 360, 600, 300, 0, 0, 0, 0, 0, 0, 0, 30, 0, 0])

bm1_ind = np.array([520, 520, 500, 500, 500, 500, 500,   0.,  0.,  1, 0,
                    24, 120, 360, 600, 300, 0, 0, 0, 0, 0, 0, 0, 30, 0, 0])


b_ind  = np.array([750, 750, 750, 750, 750, 720, 720,   0.,  0.,  1, 0,
                   24, 120, 720, 900, 600, 0, 0, 0, 0, 0, 0, 0, 30, 0, 0])

bp1_ind = np.array([1100, 1100, 1100, 1100, 1100, 1100, 1000,   0.,  0.,  1, 0,
                    24, 120, 720, 900, 600, 0, 0, 0, 0, 0, 0, 0, 30, 0, 0])

bp2_ind = np.array([1450, 1450, 1450, 1450, 1450, 1450, 1400,   0.,  0.,  1, 0,
                    24, 120, 720, 900, 600, 0, 0, 0, 0, 0, 0, 0, 30, 0, 0])

# With pretreatment 
# Ignoring CH4 in effluents 
bm2_metab = np.array([
    4.283e+00, 2.075e+00, 2.672e+01, 3.471e+00, 4.172e+00, 4.562e+00,
    1.882e+01, 0, 0, 1.920e+02, 6.133e+01, 4.862e+01,
    8.427e+01, 1.085e+01, 1.751e+01, 9.605e+00, 4.936e+01, 5.795e+01,
    4.773e+00, 3.127e+01, 1.124e+01, 2.019e+01, 2.099e+01, 7.816e+01,
    4.000e+01, 2.000e+01])


bm1_metab = np.array([
    4.376e+00, 2.096e+00, 2.765e+01, 3.615e+00, 4.219e+00, 4.691e+00,
    1.979e+01, 0, 0, 2.872e+02, 7.562e+01, 5.392e+01,
    9.475e+01, 1.095e+01, 1.762e+01, 9.761e+00, 1.028e+02, 1.061e+02,
    1.365e+01, 1.020e+02, 4.199e+01, 6.123e+01, 6.176e+01, 8.847e+01,
    4.000e+01, 2.000e+01])

b_metab = np.array([
    5.2, 2.4, 31, 3.9, 4.7, 5.2, 
    22.3,   0.,  0,  357.9, 103.3,  58.4, 
    106.4, 21.07, 26.1, 18.27, 188, 169.3, 
    29.9, 169, 68.3, 108, 105.2, 98.3, 
    0, 0])

bp1_metab = np.array([
    5.265e+00, 2.434e+00, 3.225e+01, 4.074e+00, 4.795e+00, 5.442e+00,
    2.340e+01, 0, 0, 4.073e+02, 1.161e+02, 6.003e+01,
    1.163e+02, 2.117e+01, 2.617e+01, 1.841e+01, 2.417e+02, 2.134e+02,
    4.257e+01, 2.457e+02, 1.016e+02, 1.595e+02, 1.495e+02, 1.015e+02,
    4.000e+01, 2.000e+01])

bp2_metab = np.array([
    5.368e+00, 2.462e+00, 3.364e+01, 4.206e+00, 4.899e+00, 5.671e+00,
    2.463e+01, 0, 0, 4.527e+02, 1.285e+02, 6.130e+01,
    1.258e+02, 2.126e+01, 2.626e+01, 1.855e+01, 2.945e+02, 2.563e+02,
    5.556e+01, 3.228e+02, 1.319e+02, 2.158e+02, 1.918e+02, 1.041e+02,
    4.000e+01, 2.000e+01])

industrial_wws_adm1 = {
    'bm2_ind':bm2_ind,
    'bm1_ind':bm1_ind,
    'b_ind':b_ind,
    'bp1_ind':bp1_ind,
    'bp2_ind':bp2_ind,
    'bm2_metab':bm2_metab,
    'bm1_metab':bm1_metab,
    'b_metab':b_metab,
    'bp1_metab':bp1_metab,
    'bp2_metab':bp2_metab
    }

# # Just for testing

# random_adm1_vals = np.array([110, 30, 10.,  10.,  10., 10.,  10.,   0.,  0.,  10., 24,
#                              112, 124, 36, 87, 84, 10.,  10.,  10.,  10., 10.,  10.,  10., 63,  10.,  10.])

# werf_adm1_vals = np.array([31.51, 0., 0., 0., 0., 0.,  0.,  0.,  0., 0., 23.5, 30.,  
#                             42.08917959, 96.75673469, 98.2080857,  0.,  0., 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0 ])

# metro_adm1_vals = np.array([110, 30, 0.,  0.,  0., 0.,  0.,   0.,  0.,  0., 24,
#                              112, 124, 36, 87, 84, 0.,  0.,  0.,  0., 0.,  0.,  0., 63,  0.,  0.])

#%%
casm = pc.create_asm2d_cmps(False)
cadm = pc.create_adm1_cmps(False)

asm2d_i_COD = casm.i_COD[:-1]
asm2d_f_B2COD = casm.f_BOD5_COD[:-1]
asm2d_i_N = casm.i_N[:-1]

X_S_i_N = casm.X_S.i_N
S_F_i_N = casm.S_F.i_N
asm_X_I_i_N = casm.X_I.i_N
asm_S_I_i_N = casm.S_I.i_N
    
adm1_i_COD = cadm.i_COD[:-1]
adm1_f_B2COD = cadm.f_BOD5_COD[:-1]
adm1_i_N = cadm.i_N[:-1]

X_c_i_N = cadm.X_c.i_N
X_pr_i_N = cadm.X_pr.i_N
S_aa_i_N = cadm.S_aa.i_N
adm_X_I_i_N = cadm.X_I.i_N
adm_S_I_i_N = cadm.S_I.i_N

def adm1_2_asm2d(adm_vals, return_dict=True):

    S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I, X_c, \
    X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, S_cat, S_an = adm_vals
    
    bio_cod = X_su + X_aa + X_fa + X_c4 + X_pro + X_ac + X_h2
    bio_n = np.sum((adm_vals*adm1_i_N)[16:23])

    #!!! In default ASM2d stoichiometry, biomass decay (cell lysis)
    #!!! yields 90% particulate substrate + 10% X_I
    #!!! so: convert both biomass and X_I in adm to X_S and X_I in asm
    xi_n = X_I*adm_X_I_i_N
    
    
    xs_cod = bio_cod * 0.9
    xs_ndm = xs_cod * X_S_i_N
    xi_cod = bio_cod * 0.1 + X_I
    xi_ndm = xi_cod * asm_X_I_i_N
    
    if xs_ndm > bio_n:
        warn('Not enough biomass N to map the specified proportion of '
             'biomass COD into X_S. Rest of the biomass COD goes to S_A')
        X_S = bio_n / X_S_i_N
        xs_cod -= X_S
        bio_n = 0
    else:
        X_S = xs_cod
        xs_cod = 0
        bio_n -= xs_ndm
    if xi_ndm > bio_n + xi_n + S_IN:
        warn('Not enough N in biomass and X_I to map the specified proportion of '
             'biomass COD into X_I. Rest of the biomass COD goes to S_A')
        X_I = (bio_n + xi_n + S_IN) / asm_X_I_i_N
        xi_cod -= X_I
        bio_n = xi_n = S_IN = 0
    else:
        X_I = xi_cod
        xi_cod = 0
        xi_n -= xi_ndm
        if xi_n < 0:
            bio_n += xi_n
            xi_n = 0
            if bio_n < 0:
                S_IN += bio_n
                bio_n = 0
    
    xsub_cod = X_c + X_ch + X_pr + X_li 
    xsub_n = X_c*X_c_i_N + X_pr*X_pr_i_N
    
    xs_ndm = xsub_cod * X_S_i_N
    if xs_ndm > xsub_n + bio_n:
        X_S_temp = (xsub_n + bio_n)/X_S_i_N
        X_S += X_S_temp
        xsub_cod -= X_S_temp
        xsub_n = bio_n = 0
    else:
        X_S += xsub_cod
        xsub_cod = 0
        xsub_n -= xs_ndm
        if xsub_n < 0: 
            bio_n += xsub_n
            xsub_n = 0
        
    ssub_cod = S_su + S_aa + S_fa
    ssub_n = S_aa * S_aa_i_N
    sf_ndm = ssub_cod * S_F_i_N
    if sf_ndm > ssub_n + xsub_n + bio_n:
        S_F = (ssub_n + xsub_n + bio_n) / S_F_i_N
        ssub_cod -= S_F
        ssub_n = xsub_n = bio_n = 0
    else:
        S_F = ssub_cod
        ssub_cod = 0
        ssub_n -= sf_ndm
        if ssub_n < 0:
            xsub_n += ssub_n
            ssub_n = 0
            if xsub_n < 0:
                bio_n += xsub_n
                xsub_n = 0
    
    S_A = S_ac + S_pro + S_bu + S_va

    si_cod = S_I    
    si_n = S_I * adm_S_I_i_N
    si_ndm = si_cod * asm_S_I_i_N
    if si_ndm > si_n + xi_n + S_IN:
        warn('Not enough N in S_I and X_I to map all S_I from ADM1 to ASM2d. '
             'Rest of the S_I COD goes to S_A')
        S_I = (si_n + xi_n + S_IN) / asm_S_I_i_N
        si_cod -= S_I
        si_n = xi_n = S_IN = 0
    else:
        S_I = si_cod
        si_cod = 0
        si_n -= si_ndm
        if si_n < 0:
            xi_n += si_n
            si_n = 0
            if xi_n < 0:
                S_IN += xi_n
                xi_n = 0
    
    S_NH4 = S_IN + si_n + ssub_n + xsub_n + xi_n + bio_n
    S_A += si_cod + ssub_cod + xsub_cod + xi_cod + xs_cod
    S_ALK = S_IC
    
    asm_vals = np.array(([
        0, 0, # S_O2, S_N2,
        S_NH4, 
        0, 
        0, 
        S_F, S_A, S_I, S_ALK,
        X_I, X_S, 
        0,  # X_H,
        0, 0, 0, 
        0, # X_AUT,
        0, 0]))
    
    if S_h2 > 0 or S_ch4 > 0:
        warn('Ignored dissolved H2 or CH4.')
        
    asm2d_COD =  np.sum(asm_vals*asm2d_i_COD)
    asm2d_BOD =  np.sum(asm_vals*asm2d_i_COD*asm2d_f_B2COD)
    asm2d_TN = np.sum(asm_vals*asm2d_i_N)
    
    adm1_COD =  np.sum(adm_vals*adm1_i_COD)
    adm1_BOD = np.sum(adm_vals*adm1_i_COD*adm1_f_B2COD)
    adm1_TN = np.sum(adm_vals*adm1_i_N)
    
    print(f'BOD\n'
          f'ASM2d: {asm2d_BOD}\n'
          f'ADM1:  {adm1_BOD}\n')

    if not isclose(asm2d_COD, adm1_COD, rel_tol=1e-2, abs_tol=1e-6):        
        print('COD is not balanced,\n'
              f'ASM2D COD is {asm2d_COD}.\n'
              f'while ADM1 COD is {adm1_COD}.\n'
              f'difference is {asm2d_COD - adm1_COD}.\n')
        
    if not isclose(asm2d_TN, adm1_TN, rel_tol=1e-2, abs_tol=1e-6):
        print('TN is not balanced,\n'
              f'ASM2D TKN is {asm2d_TN}.\n'
              f'while ADM1 TKN is {adm1_TN}.\n'
              f'difference is {asm2d_TN - adm1_TN}.\n'
              f'ASM2d TN array {asm_vals*asm2d_i_N}.\n'
              f'ADM1 TN array {adm_vals*adm1_i_N}.\n'
              )
    
    if return_dict:
        return {
            'S_NH4':S_NH4,
            'S_F':S_F,
            'S_A':S_A,
            'S_I':S_I,
            'S_ALK':S_ALK,
            'X_I':X_I,
            'X_S':X_S
            }
    else: return asm_vals
 
asm2d_cmps = adm1_2_asm2d(bp2_metab, return_dict=False)

#%%
if __name__ == '__main__':
    concs_asm = {}
    for k, vals in industrial_wws_adm1.items():
        concs_asm[k] = adm1_2_asm2d(vals, return_dict=True)