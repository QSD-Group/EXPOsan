# -*- coding: utf-8 -*-
"""
Created on Fri May  5 08:30:48 2023

@author: raisa
"""
import numpy as np
from warnings import warn
from math import isclose

# metro_asm2d_inf_kwargs = {
#     'concentrations': {
#         'S_I': 112.51, / 
#         'X_I': 50, /
#         'S_F': 100, / 
#         'S_A': 41.21, / 
#         'X_S': 200, / 
#         'S_NH4': 29.3, / 
#         'S_N2': 0, / 
#         'S_NO3': 0.28, / 
#         'S_PO4': 2.8, / 
#         'X_PP': 0.4,/
#         'X_PHA': 174.56, /
#         'X_H': 5, /
#         'X_AUT': 5,/
#         'X_PAO': 5, /
#         'S_ALK':7*12, / 
#         },
#     'units': ('m3/d', 'mg/L'),
#     }

# werf_asm2d_inf_kwargs = {
#     'concentrations': {
#         'S_I': 30, /
#         'X_I': 0, / 
#         'S_F': 0, /
#         'S_A': 31.504, / 
#         'X_S': 237.054, /
#         'S_NH4': 25, /
#         'S_N2': 0, /
#         'S_NO3': 0, /
#         'S_PO4': 5, /
#         'X_PP': 0,
#         'X_PHA': 54.754,
#         'X_H': 0, /
#         'X_AUT': 0, 
#         'X_PAO': 0, /
#         'S_ALK':7*12, /
#         },
#     'units': ('m3/d', 'mg/L'),
#     }

random_asm2d_vals = np.array([0, 0, 40, 0, 10, 10, 31.51, 100, 7*12, 0, 237.054, 0, 
            0, 0, 54.754, 5, 0, 0])

# Typical influent wastewater. Such wastewater wouldn't actually be influent to 
# an anaerobic digestor. Typical influent to digestor would have nitrate, oxygen, etc. 
werf_asm2d_vals = np.array([0, 0, 25, 0, 5, 0, 31.51, 30, 7*12, 0, 237.054, 0, 
            0, 0, 54.754, 0, 0, 0])

metro_asm2d_vals = np.array([0, 0, 29.3, 0.28, 2.8, 100, 41.21, 112.51, 7*12, 50, 200, 5, 5, 0.4, 174.56, 
                             5, 0, 0])

metro2_asm2d_vals = np.array([100, 0, 29.3, 0.28, 2.8, 100, 41.21, 50, 7*12, 50, 200, 5, 5, 40, 174.56, 
                             5, 10, 10])

# metro_asm2d_vals = np.array([S_O2 = 0, S_N2 = 0, S_NH4 = 29.3, S_NO3 = 0.28, S_PO4 = 2.8,  S_F = 100, 
#                         S_A = 41.21, S_I = 112.51, S_ALK = 7*12, X_I = 50, X_S = 2OO, X_H = 5, 
#             X_PAO = 5, X_PP = 0.4, X_PHA = 174.56, X_AUT = 5, X_MeOH = 0, X_MeP = 0])

def test_asm2adm(asm_vals):
    
    # asm2d = [S_O2, S_N2, S_NH4, S_NO3, S_PO4, S_F, S_A, S_I, S_ALK, X_I, X_S, X_H, 
    #          X_PAO, X_PP, X_PHA, X_AUT, X_MeOH, X_MeP]
    
    # Have removed water from all arrays
    
    asm2d_i_COD = np.array([-1.   ,  0.   ,  0.   , -4.569,  0.   ,  1.   ,  1.   ,  1.   ,
            0.   ,  1.   ,  1.   ,  1.   ,  1.   ,  0.   ,  1.   ,  1.   ,
            0.   ,  0. ])

    # X_I.i_N value was changed from 0.02 to 0.06, because of mass balance issues 05/17 (Jeremy)
    # Made S_NO3.i_N = 0, since TKN is being balanced not TN 05/18
    asm2d_i_N = np.array([0.  , 1.  , 1.  ,  0, 0.  , 0.03, 0.  , 0.01, 0.  , 0.06, 0.04,
           0.07, 0.07, 0.  , 0.  , 0.07, 0.  , 0.])
    
    asm2d_i_P = np.array([0.   , 0.   , 0.   , 0.   , 1.   , 0.01 , 0.   , 0.   , 0.   ,
           0.01 , 0.01 , 0.02 , 0.02 , 1.   , 0.   , 0.02 , 0.   , 0.205])

    # adm1 = [S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_IP, 
    #         S_I, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, 
    #         X_PHA, X_PP, X_PAO, S_K, S_Mg, X_MeOH, X_MeP, S_cat, S_an]

    adm1_i_COD = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 1., 1., 1., 1., 1.,
           1., 1., 1., 1., 1., 1., 1., 1., 0., 1., 0., 0., 0., 0., 0., 0.])

    # X_PAO.i_N value was changed from 0.068 to 0.07, because X_PAO is directly mapped 
    adm1_i_N = np.array([0.   , 0.098, 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
           0.   , 1.   , 0.   , 0.06 , 0.   , 0.098, 0.   , 0.08 , 0.08 ,
           0.08 , 0.08 , 0.08 , 0.08 , 0.08 , 0.06 , 0.   , 0.   , 0.07,
           0.   , 0.   , 0.   , 0.   , 0.   , 0.])
    
    # X_PAO.i_P value was changed from 0.019 to 0.02, because X_PAO is directly mapped 
    adm1_i_P = np.array([0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
           0.   , 0.   , 1.   , 0.01 , 0.   , 0.01 , 0.   , 0.019, 0.019,
           0.019, 0.019, 0.019, 0.019, 0.019, 0.01 , 0.   , 1.   , 0.02,
           0.   , 0.   , 0.   , 0.205, 0.   , 0.])
    
    S_O2, S_N2, S_NH4, S_NO3, S_PO4, S_F, S_A, S_I, S_ALK, X_I, X_S, X_H, \
        X_PAO, X_PP, X_PHA, X_AUT, X_MeOH, X_MeP = asm_vals

    # adm_ions_idx = cmps_adm.indices(['S_IN', 'S_IC', 'S_cat', 'S_an'])
    # adm_ions_idx = np.array([10, 9, 31, 32])
    
    # from ASM2d
    S_NO3_i_COD = -4.569
    X_H_i_N = 0.07
    X_AUT_i_N = 0.07
    S_F_i_N = 0.03
    X_S_i_N = 0.04
    # asm_X_I_i_N = 0.02 (original value)
    asm_X_I_i_N = 0.06 # to fix mass balance 05/17 (Jeremy)
    
    X_H_i_P = 0.02
    X_AUT_i_P = 0.02
    S_F_i_P = 0.01
    X_S_i_P = 0.01
    asm_X_I_i_P = 0.01
    
    # from ADM1
    S_aa_i_N = 0.098
    X_pr_i_N = 0.098
    adm_X_I_i_N = 0.06
    asm_S_I_i_N = 0.01
    S_I_i_N = 0.06
    
    X_pr_i_P = 0.01
    adm_X_I_i_P = 0.01
    adm_S_I_i_P = 0.01
    
    
    # For process splits
    xs_to_li = 0.7
    frac_deg = 0.68
    bio_to_li = 0.4
    
    # # Consideration of pH
    # pKw = 14
    # pH = 7
    # pKa_IN = 9.3
    # pKa_IC = 2.12
    # proton_charge = 10**(-pKw+pH) - 10**(-pH)
    # alpha_IN = 10**(pKa_IN-pH)/(1+10**(pKa_IN-pH))/14
    # alpha_IC = -1/(1+10**(pKa_IC-pH))/12
    
    # # Step 0: charged component snapshot (# pg. 84 of IWA ASM textbook)
    # _sno3 = S_NO3
    # _snh4 = S_NH4
    # _salk = S_ALK  
    # _spo4 = S_PO4
    # _sa = S_A 
    # _xpp = X_PP 
      
    # Step 1: remove any remaining COD demand
    O2_coddm = S_O2
    NO3_coddm = -S_NO3*S_NO3_i_COD
    
    # cod_spl = S_S + X_S + X_BH + X_BA
    # Replacing S_S with S_F + S_A (IWA ASM textbook)
    
    cod_spl = (S_A + S_F) + X_S + (X_H + X_AUT)
    
    # bioN = X_BH*X_BH_i_N + X_BA*X_BA_i_N
    
    bioN = X_H*X_H_i_N + X_AUT*X_AUT_i_N
    bioP = X_H*X_H_i_P + X_AUT*X_AUT_i_P
    
    # To be used in Step 2
    S_ND_asm1 = S_F*S_F_i_N   #S_ND (in asm1) equals the N content in S_F (Joy)
    
    # To be used in Step 3
    X_ND_asm1 = X_S*X_S_i_N   #X_ND (in asm1) equals the N content in X_S (Joy)
    
    # To be used in Step 5 (a)
    X_S_P = X_S*X_S_i_P
    # To be used in Step 5 (b)
    S_F_P = S_F*S_F_i_P 
    
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
                    X_H += X_S
                    X_S = 0
                    if X_H < 0:
                        X_AUT += X_H
                        X_H = 0
        S_O2 = S_NO3 = 0
    
    # Step 2: convert any readily biodegradable 
    # COD and TKN into amino acids and sugars
    
    # S_S (in asm1) equals to the sum of S_F and S_A (pg. 82 IWA ASM models handbook)
    S_S_asm1 = S_F + S_A 
    
    # First we calculate the amount of amino acid required in ADM1
    # if all available soluble organic N can be mapped to amino acid
    req_scod = S_ND_asm1 / S_aa_i_N
    
    # if available S_S is not enough to fulfill that amino acid requirement 
    if S_S_asm1 < req_scod: 
        # then all available S_S is mapped to amino acids 
        S_aa = S_S_asm1
        # and no S_S would be available for conversion to sugars
        S_su = 0
        # This needs to be followed by a corresponding loss in soluble organic N 
        S_ND_asm1 -= S_aa * S_aa_i_N
    # if available S_S is more than enough to fulfill that amino acid requirement 
    else:
        # All soluble organic N will be mapped to amino acid
        S_aa = req_scod
        # The line above implies that a certain portion of S_S would also be consumed to form amino acid
        # The S_S which is left would form sugar 
        # In simpler terms; S_S = S_S - S_aa; S_su = S_S 
        S_su = S_S_asm1 - S_aa
        # All soluble organic N would thus be consumed in amino acid formation
        S_ND_asm1 = 0

    # Step 3: convert slowly biodegradable COD and TKN
    # into proteins, lipids, and carbohydrates
    
    # First we calculate the amount of protein required in ADM1
    # if all available particulate organic N can be mapped to amino acid
    req_xcod = X_ND_asm1 / X_pr_i_N
    # Since X_pr_i_N >> X_pr_i_P there's no need to check req_xcod for N and P separately (CONFIRM LATER 05/16)
    
    # if available X_S is not enough to fulfill that protein requirement
    if X_S < req_xcod:
        # then all available X_S is mapped to amino acids
        X_pr = X_S
        # and no X_S would be available for conversion to lipid or carbohydrates 
        X_li = X_ch = 0
        # This needs to be followed by a corresponding loss in particulate organic N 
        X_ND_asm1 -= X_pr * X_pr_i_N
        
        # For P balance (CONFIRM LATER 05/16)
        # This needs to be followed by a corresponding loss in particulate organic P
        X_S_P -= X_pr * X_pr_i_P
        
    # if available X_S is more than enough to fulfill that protein requirement
    else:
        # All particulate organic N will be mapped to amino acid
        X_pr = req_xcod
        # The line above implies that a certain portion of X_S would also be consumed to form protein
        # The X_S which is left would form lipid and carbohydrates in a percentage define by the user  
        X_li = xs_to_li * (X_S - X_pr)
        X_ch = (X_S - X_pr) - X_li
        # All particulate organic N would thus be consumed in amino acid formation
        X_ND_asm1 = 0
        
        # For P balance (CONFIRM LATER 05/16)
        # This needs to be followed by a corresponding loss in particulate organic N 
        X_S_P -= X_pr * X_pr_i_P
    
    # Step 4: convert active biomass into protein, lipids, 
    # carbohydrates and potentially particulate TKN
    
    # First the amount of biomass N/P available for protein, lipid etc is determined
    # For this calculation, from total biomass N available the amount 
    # of particulate inert N/P expected in ADM1 is subtracted 
    
    available_bioN = bioN - (X_H + X_AUT) * (1-frac_deg) * adm_X_I_i_N
    if available_bioN < 0:
        raise RuntimeError('Not enough N in X_H and X_AUT to fully convert '
                           'the non-biodegradable portion into X_I in ADM1.')
        
    available_bioP = bioP - (X_H + X_AUT) * (1-frac_deg) * adm_X_I_i_P
    if available_bioP < 0:
        raise RuntimeError('Not enough P in X_H and X_AUT to fully convert '
                           'the non-biodegradable portion into X_I in ADM1.')
        
    # Then the amount of biomass N/P required for biomass conversion to protein is determined
    req_bioN = (X_H + X_AUT) * frac_deg * X_pr_i_N
    req_bioP = (X_H + X_AUT) * frac_deg * X_pr_i_P
    
    # Case I: if both available biomass N/P and particulate organic N/P is greater than 
    # required biomass N/P for conversion to protein
    if available_bioN + X_ND_asm1 >= req_bioN and available_bioP + X_S_P >= req_bioP:
        print('radhe shyam1')
        # then all biodegradable biomass N/P (corrsponding to protein demand) is converted to protein
        X_pr += (X_H + X_AUT) * frac_deg
        # the remaining biomass N/P is transfered as organic N/P
        X_ND_asm1 += available_bioN - req_bioN 
        X_S_P += available_bioP - req_bioP   
    
    # Case II: if available biomass N and particulate organic N is less than 
    # required biomass N for conversion to protein, but available biomass P and  
    # particulate organic P is greater than required biomass P for conversion to protein
    
    # Case III: if available biomass P and particulate organic P is less than 
    # required biomass P for conversion to protein, but available biomass N and  
    # particulate organic N is greater than required biomass N for conversion to protein
    
    # Case IV: if both available biomass N/P and particulate organic N/P is less than 
    # required biomass N/P for conversion to protein
    else:
        
        if (available_bioP + X_S_P)/X_pr_i_P < (available_bioN + X_ND_asm1)/X_pr_i_N:
            # all available P and particulate organic P is converted to protein
            bio2pr = (available_bioP + X_S_P)/X_pr_i_P
            X_pr += bio2pr
            # Biodegradable biomass available after conversion to protein is calculated 
            bio_to_split = (X_H + X_AUT) * frac_deg - bio2pr
            # Part of the remaining biomass is mapped to lipid based on user defined value 
            bio_split_to_li = bio_to_split * bio_to_li
            X_li += bio_split_to_li
            # The other portion of the remanining biomass is mapped to carbohydrates 
            X_ch += (bio_to_split - bio_split_to_li)
            # Since all organic P has been mapped to protein, none is left
            X_S_P = 0
            
            # the remaining biomass N is transfered as organic N
            X_ND_asm1 += available_bioN - (bio2pr*X_pr_i_N)
        
        else:
            # all available N and particulate organic N is converted to protein
            bio2pr = (available_bioN + X_ND_asm1)/X_pr_i_N
            X_pr += bio2pr
            # Biodegradable biomass available after conversion to protein is calculated 
            bio_to_split = (X_H + X_AUT) * frac_deg - bio2pr
            # Part of the remaining biomass is mapped to lipid based on user defined value 
            bio_split_to_li = bio_to_split * bio_to_li
            X_li += bio_split_to_li
            # The other portion of the remanining biomass is mapped to carbohydrates 
            X_ch += (bio_to_split - bio_split_to_li)
            # Since all organic N has been mapped to protein, none is left
            X_ND_asm1 = 0
            
            # the remaining biomass P is transfered as organic P
            X_S_P += available_bioP - (bio2pr*X_pr_i_P)
    
    # Step 5: map particulate inerts
    
    # 5 (a)
    # First determine the amount of particulate inert N/P available from ASM2d
    xi_nsp_asm2d = X_I * asm_X_I_i_N
    xi_psp_asm2d = X_I * asm_X_I_i_P
    
    # Then determine the amount of particulate inert N/P that could be produced 
    # in ADM1 given the ASM1 X_I
    xi_ndm = X_I * adm_X_I_i_N
    xi_pdm = X_I * adm_X_I_i_P

    # if particulate inert N available in ASM2d is greater than ADM1 demand
    if xi_nsp_asm2d + X_ND_asm1 >= xi_ndm:
        print('durga ma1')
        deficit = xi_ndm - xi_nsp_asm2d
        # COD balance 
        X_I += (X_H+X_AUT) * (1-frac_deg)
        # N balance 
        X_ND_asm1 -= deficit
        # P balance 
        if xi_psp_asm2d + X_S_P >= xi_pdm:
            deficit = xi_pdm - xi_psp_asm2d
            X_S_P -= deficit
        elif isclose(xi_psp_asm2d+X_S_P, xi_pdm, rel_tol=1e-2, abs_tol=1e-6):
            X_S_P  = 0
        else:
            raise RuntimeError('Not enough P in X_I, X_S to fully '
                               'convert X_I in ASM2d into X_I in ADM1.')
    elif isclose(xi_nsp_asm2d+X_ND_asm1, xi_ndm, rel_tol=1e-2, abs_tol=1e-6):
        print('durga ma2')
        # COD balance 
        X_I += (X_H+X_AUT) * (1-frac_deg)
        # N balance 
        X_ND_asm1 = 0
        # P balance
        if xi_psp_asm2d + X_S_P >= xi_pdm:
            deficit = xi_pdm - xi_psp_asm2d
            X_S_P -= deficit
        elif isclose(xi_psp_asm2d+X_S_P, xi_pdm, rel_tol=1e-2, abs_tol=1e-6):
            X_S_P  = 0
        else:
            raise RuntimeError('Not enough P in X_I, X_S to fully '
                               'convert X_I in ASM2d into X_I in ADM1.')
    else:
        raise RuntimeError('Not enough N in X_I, X_ND_asm1, X_S_P to fully '
                           'convert X_I in ASM2d into X_I in ADM1.')
        
    # 5(b)
    
    # Then determine the amount of soluble inert N/P that could be produced 
    # in ADM1 given the ASM2d X_I
    # S_I_i_N is for ADM1 
    
    req_sn = S_I * S_I_i_N
    req_sp = S_I * adm_S_I_i_P
    
    supply_inert_n_asm2d = S_I * asm_S_I_i_N
    
    # N balance 
    if req_sn <= S_ND_asm1 + supply_inert_n_asm2d:
        print('sita ram1')
        S_ND_asm1 -= (req_sn - supply_inert_n_asm2d)
        supply_inert_n_asm2d = 0 
        # P balance 
        if req_sp <= S_F_P:
            S_F_P -= req_sp
        elif req_sp <= S_F_P + X_S_P:
            X_S_P -= (req_sp - S_F_P)
            S_F_P = 0
        elif req_sp <= S_F_P + X_S_P + S_PO4:
            S_PO4 -= (req_sp - S_F_P - X_S_P)
            S_F_P = X_S_P = 0
        else:
            S_PO4 -= (req_sp - S_F_P - X_S_P)
            S_F_P =  X_S_P = 0
    # N balance
    elif req_sn <= S_ND_asm1 + X_ND_asm1 + supply_inert_n_asm2d:
        print('sita ram2')
        X_ND_asm1 -= (req_sn - S_ND_asm1 - supply_inert_n_asm2d)
        S_ND_asm1 = supply_inert_n_asm2d = 0
        # P balance
        if req_sp <= S_F_P:
            S_F_P -= req_sp
        elif req_sp <= S_F_P + X_S_P:
            X_S_P -= (req_sp - S_F_P)
            S_F_P = 0
        elif req_sp <= S_F_P + X_S_P + S_PO4:
            S_PO4 -= (req_sp - S_F_P - X_S_P)
            S_F_P = X_S_P = 0
        else:
            S_PO4 -= (req_sp - S_F_P - X_S_P)
            S_F_P =  X_S_P = 0
    # N balance
    elif req_sn <= S_ND_asm1 + X_ND_asm1 + S_NH4:
        print('sita ram3')
        S_NH4 -= (req_sn - S_ND_asm1 - X_ND_asm1 - supply_inert_n_asm2d)
        S_ND_asm1 = X_ND_asm1 = supply_inert_n_asm2d = 0
        # P balance 
        if req_sp <= S_F_P:
            print('saibaba1')
            S_F_P -= req_sp
        elif req_sp <= S_F_P + X_S_P:
            print('saibaba2')
            X_S_P -= (req_sp - S_F_P)
            S_F_P = 0
        elif req_sp <= S_F_P + X_S_P + S_PO4:
            print('saibaba3')
            S_PO4 -= (req_sp - S_F_P - X_S_P)
            S_F_P = X_S_P = 0
        else:
            S_PO4 -= (req_sp - S_F_P - X_S_P)
            S_F_P =  X_S_P = 0
            
    elif req_sp <= S_F_P or req_sp <= S_F_P + X_S_P or req_sp <= S_F_P + X_S_P + S_PO4:
        
        S_NH4 -= (req_sn - S_ND_asm1 - X_ND_asm1 - supply_inert_n_asm2d)
        S_ND_asm1 = X_ND_asm1 = supply_inert_n_asm2d = 0
        
        if req_sp <= S_F_P:
            S_F_P -= req_sp
        elif req_sp <= S_F_P + X_S_P:
            X_S_P -= (req_sp - S_F_P)
            S_F_P = 0
        elif req_sp <= S_F_P + X_S_P + S_PO4:
            S_PO4 -= (req_sp - S_F_P - X_S_P)
            S_F_P = X_S_P = 0
    else:
        print('sita ram5')
        if (S_ND_asm1 + X_ND_asm1 + S_NH4 + supply_inert_n_asm2d)/S_I_i_N < (S_F_P + X_S_P + S_PO4)/adm_S_I_i_P:
            warn('Additional soluble inert COD is mapped to S_su.')
            SI_cod = (S_ND_asm1 + X_ND_asm1 + S_NH4 + supply_inert_n_asm2d)/S_I_i_N
            S_su += S_I - SI_cod
            S_I = SI_cod
            S_ND_asm1 = X_ND_asm1 = S_NH4 = supply_inert_n_asm2d = 0
            
            req_sp = S_I * adm_S_I_i_P
            S_PO4 -= (req_sp - S_F_P - X_S_P)
            S_F_P = X_S_P = 0
        else:
            warn('Additional soluble inert COD is mapped to S_su.')
            SI_cod = (S_F_P + X_S_P + S_PO4)/adm_S_I_i_P
            S_su += S_I - SI_cod
            S_I = SI_cod
            S_F_P = X_S_P = S_PO4 = 0
            
            req_sn = S_I * S_I_i_N
            S_NH4 -= (req_sn - S_ND_asm1 - X_ND_asm1 - supply_inert_n_asm2d)
            S_ND_asm1 = X_ND_asm1 = supply_inert_n_asm2d = 0
        
    # Step 6: Step map any remaining TKN/P
    S_IN = S_ND_asm1 + X_ND_asm1 + S_NH4
    S_IP = S_F_P + X_S_P + S_PO4 

    # Step 8: check COD and TKN balance
    # has TKN: S_aa, S_IN, S_I, X_pr, X_I
    S_IC = S_cat = S_an = 0
    
    # Step 9: Mapping common state variables directly    
    # The next three commented lines are executed when outputting
    # array of ADM1 components 
    # X_PAO (ADM1) = X_PAO (ASM2d)
    # X_PP (ADM1) = X_PP (ASM2d)
    # X_PHA (ADM1) = X_PHA (ASM2d)
    # X_MeOH (ADM1) = X_MeOH (ASM2d)
    # X_MeP (ADM1) = X_MeP (ASM2d)
    # S_IP = S_PO4 # correct, but not using to balance P 
        
    # Step 6: map any remaining TKN
    S_IN = S_ND_asm1 + X_ND_asm1 + S_NH4
    S_IP = S_F_P + X_S_P + S_PO4
    #X_PP = X_PP + X_H*frac_deg*X_H_i_P + X_AUT*frac_deg*X_AUT_i_P     
    
    # Step 8: check COD and TKN balance
    # has TKN: S_aa, S_IN, S_I, X_pr, X_I
    S_IC = S_cat = S_an = 0
    
    # adm_vals = np.array([
    #     S_su, S_aa, 
    #     0, 0, 0, 0, 0, # S_fa, S_va, S_bu, S_pro, S_ac, 
    #     0, 0, # S_h2, S_ch4,
    #     S_IC, S_IN, S_I, 
    #     0, # X_c, 
    #     X_ch, X_pr, X_li, 
    #     0, 0, 0, 0, 0, 0, 0, # X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2,
    #     X_I, S_cat, S_an, H2O])
    
    adm_vals = np.array([
        S_su, S_aa, 
        0, 0, 0, 0, 0, # S_fa, S_va, S_bu, S_pro, S_ac, 
        0, 0, # S_h2, S_ch4,
        S_IC, S_IN, S_IP, S_I, 
        X_ch, X_pr, X_li, 
        0, 0, 0, 0, 0, 0, 0, # X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2,
        X_I, X_PHA, X_PP, X_PAO, 
        0, 0,  # S_K, S_Mg,
        X_MeOH, X_MeP,
        S_cat, S_an])
    
    # adm_vals = f_corr(asm_vals, adm_vals)
    
    # # Step 7: charge balance
    # asm_charge_tot = - _sa/64 + _snh4/14 - _sno3/14 - 1.5*_spo4/31 - _salk - _xpp/31 #Based on page 84 of IWA ASM handbook
    # #!!! charge balance should technically include VFAs, 
    # # but VFAs concentrations are assumed zero per previous steps??
    # S_IN = adm_vals[adm_ions_idx[0]]
    # S_IC = (asm_charge_tot -  S_IN*alpha_IN)/alpha_IC
    # net_Scat = asm_charge_tot + proton_charge
    # if net_Scat > 0:  
    #     S_cat = net_Scat
    #     S_an = 0
    # else:
    #     S_cat = 0
    #     S_an = -net_Scat
    
    # adm_vals[adm_ions_idx[1:]] = [S_IC, S_cat, S_an]
    
    asm2d_COD =  np.sum(asm_vals*asm2d_i_COD)
    asm2d_TN = np.sum(asm_vals*asm2d_i_N)
    asm2d_TP = np.sum(asm_vals*asm2d_i_P)
    
    adm1_COD =  np.sum(adm_vals*adm1_i_COD)
    adm1_TN = np.sum(adm_vals*adm1_i_N)
    adm1_TP = np.sum(adm_vals*adm1_i_P)
        
    if isclose(asm2d_COD, adm1_COD, rel_tol=1e-2, abs_tol=1e-6):
        print('COD is balanced\n'
              f'ASM2D COD is {asm2d_COD}.\n'
              f'while ADM1 COD is {adm1_COD}.\n'
              f'difference is {asm2d_COD - adm1_COD}.\n')
    else:
        print('COD is not balanced,\n'
              f'ASM2D COD is {asm2d_COD}.\n'
              f'while ADM1 COD is {adm1_COD}.\n'
              f'difference is {asm2d_COD - adm1_COD}.\n')
        
    if isclose(asm2d_TN, adm1_TN, rel_tol=1e-2, abs_tol=1e-6):
        print('TN is balanced\n'
              f'ASM2D TN is {asm2d_TN}.\n'
              f'while ADM1 TN is {adm1_TN}.\n'
              f'difference is {asm2d_TN - adm1_TN}.\n'
              )
    else:
        print('TN is not balanced,\n'
              f'ASM2D TKN is {asm2d_TN}.\n'
              f'while ADM1 TKN is {adm1_TN}.\n'
              f'difference is {asm2d_TN - adm1_TN}.\n'
              f'ASM2d TN array {asm_vals*asm2d_i_N}.\n'
              f'ADM1 TN array {adm_vals*adm1_i_N}.\n'
              )
        
    if isclose(asm2d_TP, adm1_TP, rel_tol=1e-2, abs_tol=1e-6):
        print('TP is balanced\n'
              f'ASM2D TP is {asm2d_TP}.\n'
              f'while ADM1 TP is {adm1_TP}.\n'
              f'difference is {asm2d_TP - adm1_TP}.\n')
    else:
        print('TP is not balanced,\n'
              f'ASM2D TP is {asm2d_TP}.\n'
              f'while ADM1 TP is {adm1_TP}.\n'
              f'difference is {asm2d_TP - adm1_TP}.\n'
              f'ASM2d TP array {asm_vals*asm2d_i_P}.\n'
              f'ADM1 TP array {adm_vals*adm1_i_P}.\n'
              )
    
    return

mymy = test_asm2adm(metro2_asm2d_vals)