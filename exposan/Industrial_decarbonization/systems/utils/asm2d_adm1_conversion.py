#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:40:23 2024

@author: saumitrarai
"""

import numpy as np
from warnings import warn
from math import isclose
import qsdsan.processes as pc

#%%

cmps_asm = pc.create_asm2d_cmps(False)
cmps_adm = pc.create_adm1_cmps(False)

asm2d_i_COD = cmps_asm.i_COD[:-1]
asm2d_f_B2COD = cmps_asm.f_BOD5_COD[:-1]

# NEW
# asm2d_i_N = cmps_asm.i_N[:-1]

# OLD
### Change array while checking N balance ###

asm2d_i_TN = np.zeros(18)
asm2d_i_TN[0:9] = cmps_asm.i_N[:9]
asm2d_i_TN[9] = 0.06
asm2d_i_TN[10:] = cmps_asm.i_N[10:-1]
asm2d_i_N = asm2d_i_TN

### Change array while checking N balance ###

adm1_i_COD = cmps_adm.i_COD[:-1]
adm1_f_B2COD = cmps_adm.f_BOD5_COD[:-1]
adm1_i_N = cmps_adm.i_N[:-1]

# User defined values
xs_to_li = 0.7
bio_to_li = 0.4
frac_deg = 0.68
rtol = 1e-2
atol = 1e-6

S_NO3_i_COD = cmps_asm.S_NO3.i_COD
X_H_i_N = cmps_asm.X_H.i_N
X_AUT_i_N = cmps_asm.X_AUT.i_N
X_PAO_i_N = cmps_asm.X_PAO.i_N
S_F_i_N = cmps_asm.S_F.i_N
X_S_i_N = cmps_asm.X_S.i_N

# While using this interface X_I.i_N in ASM2d should be 0.06, instead of 0.02. 

# NEW
# asm_X_I_i_N = cmps_asm.X_I.i_N

# OLD
asm_X_I_i_N = cmps_adm.X_I.i_N # which is 0.06

asm_S_I_i_N = cmps_asm.S_I.i_N

S_aa_i_N = cmps_adm.S_aa.i_N
X_pr_i_N = cmps_adm.X_pr.i_N
adm_X_I_i_N = cmps_adm.X_I.i_N
adm_S_I_i_N = cmps_adm.S_I.i_N

def asm2adm(asm_vals):
    # S_I, S_S, X_I, X_S, X_BH, X_BA, X_P, S_O, S_NO, S_NH, S_ND, X_ND, S_ALK, S_N2, H2O = asm_vals
    
    S_O2, S_N2, S_NH4, S_NO3, S_PO4, S_F, S_A, S_I, S_ALK, X_I, X_S, X_H, \
        X_PAO, X_PP, X_PHA, X_AUT, X_MeOH, X_MeP = asm_vals
      
    # Step 1: remove any remaining COD demand
    O2_coddm = S_O2
    NO3_coddm = -S_NO3*S_NO3_i_COD
    
    # cod_spl = S_S + X_S + X_BH + X_BA
    # Replacing S_S with S_F + S_A (IWA ASM textbook)
    
    cod_spl = (S_A + S_F) + (X_S + X_PHA) + (X_H + X_AUT + X_PAO)
    
    bioN = X_H*X_H_i_N + X_AUT*X_AUT_i_N + X_PAO*X_PAO_i_N
    
    # To be used in Step 2
    S_F_N = S_F*S_F_i_N   #S_ND (in asm1) equals the N content in S_F
    # To be used in Step 3
    X_S_N = X_S*X_S_i_N   #X_ND (in asm1) equals the N content in X_S
    
    if cod_spl <= O2_coddm:
        S_O2 = O2_coddm - cod_spl
        S_F = S_A =  X_S = X_H = X_AUT = 0
    elif cod_spl <= O2_coddm + NO3_coddm:
        S_O2 = 0
        S_NO3 = -(O2_coddm + NO3_coddm - cod_spl)/S_NO3_i_COD
        S_A = S_F = X_S = X_H = X_AUT = 0
    else:
        S_A -= O2_coddm + NO3_coddm
        if S_A < 0:
            S_F += S_A
            S_A = 0
            if S_F < 0:
                X_S += S_F
                S_F = 0
                if X_S < 0:
                    X_PHA += X_S
                    X_S = 0
                    if X_PHA < 0:
                        X_H += X_PHA
                        X_PHA = 0
                        if X_H < 0:
                            X_AUT += X_H
                            X_H = 0
                            if X_AUT < 0:
                                X_PAO += X_AUT
                                X_AUT = 0
                    
        S_O2 = S_NO3 = 0
    
    # Step 2: convert any readily biodegradable 
    # COD and TKN into amino acids and sugars
    
    # S_S (in asm1) equals to the sum of S_F and S_A (pg. 82 IWA ASM models handbook)
    S_S_asm1 = S_F + S_A 
    
    # First we calculate the amount of amino acid required in ADM1
    # if all available soluble organic N can be mapped to amino acid
    req_scod = S_F_N / S_aa_i_N
    
    # if available S_S is not enough to fulfill that amino acid requirement 
    if S_S_asm1 < req_scod: 
        # then all available S_S is mapped to amino acids 
        S_aa = S_S_asm1
        # and no S_S would be available for conversion to sugars
        S_su = 0
        # This needs to be followed by a corresponding loss in soluble organic N 
        S_F_N -= S_aa * S_aa_i_N
    # if available S_S is more than enough to fulfill that amino acid requirement 
    else:
        # All soluble organic N will be mapped to amino acid
        S_aa = req_scod
        # The line above implies that a certain portion of S_S would also be consumed to form amino acid
        # The S_S which is left would form sugar 
        # In simpler terms; S_S = S_S - S_aa; S_su = S_S 
        S_su = S_S_asm1 - S_aa
        # All soluble organic N would thus be consumed in amino acid formation
        S_F_N = 0

    # Step 3: convert slowly biodegradable COD and TKN
    # into proteins, lipids, and carbohydrates
    
    # First we calculate the amount of protein required in ADM1
    # if all available particulate organic N can be mapped to protein
    req_xcod = X_S_N / X_pr_i_N
    
    # if available X_S is not enough to fulfill that protein requirement
    if X_S < req_xcod:
        # then all available X_S is mapped to amino acids
        X_pr = X_S
        # and no X_S would be available for conversion to lipid or carbohydrates 
        
        X_li = xs_to_li * X_PHA
        X_ch = (1 - xs_to_li)*X_PHA 
       
        # This needs to be followed by a corresponding loss in particulate organic N 
        X_S_N -= X_pr * X_pr_i_N
        
    # if available X_S is more than enough to fulfill that protein requirement
    else:
        # All particulate organic N will be mapped to amino acid
        X_pr = req_xcod
        # The line above implies that a certain portion of X_S would also be consumed to form protein
        # The X_S which is left would form lipid and carbohydrates in a percentage define by the user  
        X_li = xs_to_li * (X_S + X_PHA - X_pr)
        X_ch = (X_S + X_PHA - X_pr) - X_li
        # All particulate organic N would thus be consumed in amino acid formation
        X_S_N = 0
    
    # Step 4: convert active biomass into protein, lipids, 
    # carbohydrates and potentially particulate TKN
    
    # First the amount of biomass N available for protein, lipid etc is determined
    # For this calculation, from total biomass N available the amount 
    # of particulate inert N expected in ADM1 is subtracted 
    
    available_bioN = bioN - (X_H + X_AUT + X_PAO) * (1-frac_deg) * adm_X_I_i_N
    
    if available_bioN < 0:
        raise RuntimeError('Not enough N in X_H, X_AUT and X_PAO to fully convert '
                           'the non-biodegradable portion into X_I in ADM1.')
        
    # Then the amount of biomass N required for biomass conversion to protein is determined
    req_bioN = (X_H + X_AUT + X_PAO) * frac_deg * X_pr_i_N
    # req_bioP = (X_H + X_AUT) * frac_deg * X_pr_i_P
    
    # If available biomass N and particulate organic N is greater than 
    # required biomass N for conversion to protein
    if available_bioN + X_S_N >= req_bioN:
        # then all biodegradable biomass N (corrsponding to protein demand) is converted to protein
        X_pr += (X_H + X_AUT + X_PAO) * frac_deg
        # the remaining biomass N is transfered as organic N
        X_S_N += available_bioN - req_bioN 
    else:
        # all available N and particulate organic N is converted to protein
        bio2pr = (available_bioN + X_S_N)/X_pr_i_N
        X_pr += bio2pr
        # Biodegradable biomass available after conversion to protein is calculated 
        bio_to_split = (X_H + X_AUT + X_PAO) * frac_deg - bio2pr
        # Part of the remaining biomass is mapped to lipid based on user defined value 
        bio_split_to_li = bio_to_split * bio_to_li
        X_li += bio_split_to_li
        # The other portion of the remanining biomass is mapped to carbohydrates 
        X_ch += (bio_to_split - bio_split_to_li)
        # Since all organic N has been mapped to protein, none is left
        X_S_N = 0
    
    # Step 5: map particulate inerts
    
    # 5 (a)
    # First determine the amount of particulate inert N available from ASM2d
    xi_nsp_asm2d = X_I * asm_X_I_i_N
    
    # Then determine the amount of particulate inert N that could be produced 
    # in ADM1 given the ASM1 X_I
    xi_ndm = X_I * adm_X_I_i_N

    # if particulate inert N available in ASM1 is greater than ADM1 demand
    if xi_nsp_asm2d + X_S_N >= xi_ndm:
        deficit = xi_ndm - xi_nsp_asm2d
        # COD balance 
        X_I += (X_H + X_AUT + X_PAO) * (1-frac_deg)
        # N balance 
        X_S_N -= deficit
    elif isclose(xi_nsp_asm2d+X_S_N, xi_ndm, rel_tol=rtol, abs_tol=atol):
        # COD balance 
        X_I += (X_H + X_AUT + X_PAO) * (1-frac_deg)
        # N balance 
        X_S_N = 0
        
    ## NEW ###
    # elif xi_nsp_asm2d + X_ND_asm1 <= xi_ndm:
        
    #     excess = xi_ndm - xi_nsp_asm2d
        
    #     # COD balance 
    #     X_I += (X_H + X_AUT + X_PAO) * (1-frac_deg)
        
    #     # N balance 
    #     X_ND_asm1 -= excess  # Negative X_ND_asm1
        
    ## ### ###
    
    # OLD #
    else:
        raise RuntimeError('Not enough N in X_I, X_S to fully '
                            'convert X_I in ASM2d into X_I in ADM1.')
        
    # 5(b)
    
    # Then determine the amount of soluble inert N that could be produced 
    # in ADM1 given the ASM1 X_I
    req_sn = S_I * adm_S_I_i_N
    supply_inert_n_asm2d = S_I * asm_S_I_i_N
    
    # N balance 
    if req_sn <= S_F_N + supply_inert_n_asm2d:
        S_F_N -= (req_sn - supply_inert_n_asm2d)
        supply_inert_n_asm2d = 0 
    # N balance
    elif req_sn <= S_F_N + X_S_N + supply_inert_n_asm2d:
        X_S_N -= (req_sn - S_F_N - supply_inert_n_asm2d)
        S_F_N = supply_inert_n_asm2d = 0
    # N balance
    elif req_sn <= S_F_N + X_S_N + S_NH4 + supply_inert_n_asm2d:
        S_NH4 -= (req_sn - S_F_N - X_S_N - supply_inert_n_asm2d)
        S_F_N = X_S_N = supply_inert_n_asm2d = 0
    else:
        warn('Additional soluble inert COD is mapped to S_su.')
        SI_cod = (S_F_N + X_S_N + S_NH4 + supply_inert_n_asm2d)/adm_S_I_i_N
        S_su += S_I - SI_cod
        S_I = SI_cod
        S_F_N = X_S_N = S_NH4 = supply_inert_n_asm2d = 0
        
    # Step 6: Step map any remaining TKN/P
    S_IN = S_F_N + X_S_N + S_NH4 + supply_inert_n_asm2d
    
    # Step 8: check COD and TKN balance
    # has TKN: S_aa, S_IN, S_I, X_pr, X_I
    
    # these will not be zero if we do charge balance
    S_IC = S_cat = S_an = 0
    
    adm_vals = np.array([
        S_su, S_aa, 
        0, 0, 0, 0, 0, # S_fa, S_va, S_bu, S_pro, S_ac, 
        0, 0, # S_h2, S_ch4,
        S_IC, S_IN, S_I, 
        0, # X_c, 
        X_ch, X_pr, X_li, 
        0, 0, 0, 0, 0, 0, 0, # X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2,
        X_I, S_cat, S_an])
    
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
    
    return adm_vals

#%%

# asm2d_cmps = np.array([ 0.,  0.,  196.04,    0.,    0.,   41.47, 39.406,   61.3,  452.7,  250.966, 1513.664,  10.,  10.,
#                        0.,    0.,    10.,    0.,    0.])

# ([S_O2, S_N2, S_NH4, S_NO3, S_PO4, S_F, S_A, S_I, S_ALK, X_I, X_S, X_H, X_PAO, 
#   X_PP, X_PHA, X_AUT, X_MeOH, X_MeP])

asm2d_cmps = np.array([ 10.,  0.,  196.04,    0,    10,   41.47, 39.406,   61.3,  452.7,  250.966, 1513.664,  10.,  10.,
                       10.,    100.,    10.,    10.,    10.])

adm1_cmps = asm2adm(asm2d_cmps)
print(f'adm1 components = {adm1_cmps}')