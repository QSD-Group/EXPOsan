# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 13:10:02 2022

@author: joy_c

Nopens, I.; Batstone, D. J.; Copp, J. B.; Jeppsson, U.; Volcke, E.; 
Alex, J.; Vanrolleghem, P. A. An ASM/ADM Model Interface for Dynamic 
Plant-Wide Simulation. Water Res. 2009, 43, 1913–1923.
"""
from qsdsan import processes as pc, Components, Process, set_thermo
from warnings import warn
import numpy as np


cmps_asm = pc.create_asm1_cmps(False)
cmps_adm = pc.create_adm1_cmps(False)
# cmps = Components([*cmps_asm, *cmps_adm])
# cmps.compile()
# set_thermo(cmps)

# p1 = Process('biomass_convert', 
#              reaction='X_BH + X_ND -> [?]X_pr + [?]X_li + [?]X_ch + [?]X_I',
#              ref_component='X_BH',
#              conserved_for=('COD', 'N', 'P', 'mass'))

p1 = Process('biomass_convert', 
              reaction='X_BH + [?]X_ND -> X_pr + [0.32]X_I',
              ref_component='X_BH',
              conserved_for=('N',))

li_ch_split_XS = [0.7, 0.3]
li_ch_split_bio = [0.4, 0.6]
frac_deg = 0.68

# _acid_base_pairs = (('H+', 'OH-'), ('NH4+', 'NH3'), ('CO2', 'HCO3-'),
#                     ('HAc', 'Ac-'), ('HPr', 'Pr-'),
#                     ('HBu', 'Bu-'), ('HVa', 'Va-'))

#!!! these values should be retrived from the ADM1 object's rate function parameter set
T_base = 273.15 + 25
pKa_base = np.array([14, 9.25, 6.35, 4.76, 4.88, 4.82, 4.86])
Ka_dH = np.array([55900, 51965, 7646, 0, 0, 0, 0])
vfa_cod = np.array([64, 112, 160, 208])

def calc_pKa(T):
    return pKa_base - np.log10(np.exp(pc.T_correction_factor(T_base, T, Ka_dH)))

def calc_vfa_alpha(pKa_vfa, pH):
    return 1.0/vfa_cod*(-1.0/(1.0 + 10**(pKa_vfa - pH)))

def asm2adm(asm_vals, T, pH):
    
    S_I, S_S, X_I, X_S, X_BH, X_BA, X_P, S_O, S_NO, S_NH, S_ND, X_ND, S_ALK, S_N2, H2O = asm_vals
    
    S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, \
        X_c, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2 = 0
    
    # Step 0: charged component snapshot
    _sno = S_NO
    _snh = S_NH
    _salk = S_ALK
    
    # Step 1: remove any remaining COD demand
    O_coddm = S_O
    NO_coddm = -S_NO*cmps_asm.S_NO.i_COD
    cod_spl = S_S + X_S + X_BH + X_BA
    
    if cod_spl <= O_coddm:
        S_O = O_coddm - cod_spl
        S_S, X_S, X_BH, X_BA = 0
    elif cod_spl <= O_coddm + NO_coddm:
        S_O = 0
        S_NO = -(O_coddm + NO_coddm - cod_spl)/cmps_asm.S_NO.i_COD
        S_S, X_S, X_BH, X_BA = 0
    else:
        S_S -= O_coddm + NO_coddm
        if S_S < 0:
            X_S += S_S
            S_S = 0
        if X_S < 0:
            X_BH += X_S
            X_S = 0
        if X_BH < 0:
            X_BA += X_BH
            X_BH = 0
        S_O, S_NO = 0
    
    # Step 2: convert any readily biodegradable 
    # COD and TKN into amino acids and sugars
    req_scod = S_ND / cmps_adm.S_aa.i_N
    if S_S < req_scod:
        S_aa = S_S
        S_su = 0
        S_ND -= S_aa * cmps_adm.S_aa.i_N
    else:
        S_aa = req_scod
        S_su = S_S - S_aa
        S_ND = 0
    S_S = 0
    
    # Step 3: convert slowly biodegradable COD and TKN
    # into proteins, lipids, and carbohydrates
    req_xcod = X_ND / cmps_adm.X_pr.i_N
    if X_S < req_xcod:
        X_pr = X_S
        X_li, X_ch = 0
        X_ND -= X_pr * cmps_adm.X_pr.i_N
    else:
        X_pr = req_xcod
        X_li, X_ch = [frac*(X_S - X_pr) for frac in li_ch_split_XS]
        X_ND = 0
    X_S = 0
    
    # Step 4: convert active biomass into protein, lipids, 
    # carbohydrates and potentially particulate TKN
    available_bioN = X_BH * cmps_asm.X_BH.i_N \
        + X_BA * cmps_asm.X_BA.i_N \
        - (X_BH+X_BA) * (1-frac_deg) * cmps_adm.X_I.i_N
    if available_bioN < 0:
        raise RuntimeError('Not enough N in X_BA and X_BH to fully convert the non-biodegrable'
                           'portion into X_I in ADM1.')
    req_bioN = (X_BH+X_BA) * frac_deg * cmps_adm.X_pr.i_N
    if available_bioN + X_ND >= req_bioN:
        X_pr += (X_BH+X_BA) * frac_deg
        X_ND += available_bioN - req_bioN
    else:
        bio2pr = (available_bioN + X_ND)/cmps_adm.X_pr.i_N
        X_pr += bio2pr
        X_li += ((X_BH+X_BA) * frac_deg - bio2pr) * li_ch_split_bio[0]
        X_ch += ((X_BH+X_BA) * frac_deg - bio2pr) * li_ch_split_bio[1]
        X_ND = 0
    X_BH, X_BA = 0
    
    # Step 5: map particulate inerts
    if cmps_asm.X_P.i_N * X_P + cmps_asm.X_I.i_N * X_I + X_ND < (X_P+X_I) * cmps_adm.X_I.i_N:
        raise RuntimeError('Not enough N in X_I, X_P, X_ND to fully convert X_I and X_P'
                           'into X_I in ADM1.')
    deficit = (X_P+X_I) * cmps_adm.X_I.i_N - cmps_asm.X_P.i_N * X_P + cmps_asm.X_I.i_N * X_I
    X_I = X_I + X_P + (X_BH+X_BA) * (1-frac_deg)
    X_ND -= deficit
    
    req_sn = S_I * cmps_adm.S_I.i_N
    if req_sn <= S_ND:
        S_ND -= req_sn
    elif req_sn <= S_ND + X_ND:
        X_ND -= (req_sn - S_ND)
        S_ND = 0
    elif req_sn <= S_ND + X_ND + S_NH:
        S_NH -= (req_sn - S_ND - X_ND)
        S_ND, X_ND = 0
    else:
        warn('Additional soluble inert COD is mapped to S_su.')
        SI_cod = (S_ND + X_ND + S_NH)/cmps_adm.S_I.i_N
        S_su += S_I - SI_cod
        S_I = SI_cod
        S_ND, X_ND, S_NH = 0
        
    # Step 6: maps any remaining nitrogen
    S_IN = S_ND + X_ND + S_NH
    
    # Step 7: charge balance
    asm_charge_tot = _snh/14 - _sno/14 - _salk/12
    pKw, pKa_IN, pKa_IC = calc_pKa(T)[:2]
    alpha_IN = 10**(pKa_IN-pH)/(1+10**(pKa_IN-pH))/14 # charge per g N
    alpha_IC = -1/(1+10**(pKa_IC-pH))/12 # charge per g C
    #!!! charge balance should technically include VFAs, 
    # but VFAs concentrations are assumed zero per previous steps??
    S_IC = (asm_charge_tot -  S_IN*alpha_IN)/alpha_IC
    net_Scat = asm_charge_tot + 10**(-pKw + pH) - 10**(-pH)   
    if net_Scat > 0:  
        S_cat = net_Scat
        S_an = 0
    else:
        S_cat = 0
        S_an = -net_Scat
    
    # Step 8: check COD and TKN balance
    adm_vals = np.array([S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, 
                        S_h2, S_ch4, S_IC, S_IN, S_I, X_c, X_ch, 
                        X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, 
                        X_ac, X_h2, X_I, S_cat, S_an, H2O])
    
    assert sum(asm_vals*cmps_asm.i_COD) == sum(adm_vals*cmps_adm.i_COD)
    assert sum(asm_vals*cmps_asm.i_N) - sum(asm_vals[cmps_asm.indices(('S_NO', 'S_N2'))]) \
        == sum(adm_vals*cmps_adm.i_N)
    
    # unit conversion from mg/L (ASM) to kg/m3 (ADM)
    return adm_vals/1000

xs_xp_split_bio = [0.7, 0.3]

def adm2asm(adm_vals, T, pH):    
    S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I, \
        X_c, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, \
        S_cat, S_an, H2O = adm_vals
    
    X_BH, X_BA, S_O, S_NO, S_N2 = 0
    
    # Step 0: snapshot of charged components
    _ions = np.array([S_IN, S_IC, S_ac, S_pro, S_bu, S_va])
    
    # Step 1a: convert biomass into X_S+X_ND and X_P
    bio_cod = X_su + X_aa + X_fa + X_c4 + X_pro + X_ac + X_h2
    bio_n = sum((adm_vals*cmps_adm.i_N)[cmps_adm.indices(('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2'))])
    xs_cod, xp_cod = [frac*bio_cod for frac in xs_xp_split_bio]
    xp_ndm = xp_cod*cmps_asm.X_P.i_N
    if xp_ndm > bio_n:
        X_P = bio_n/cmps_asm.i_N
        bio_n = 0
    else:
        X_P = xp_cod
        bio_n -= xp_ndm
    X_S = bio_cod - X_P
    xs_ndm = X_S*cmps_asm.X_S.i_N
    if xs_ndm <= bio_n:
        X_ND = bio_n - xs_ndm
        bio_n = 0
    elif xs_ndm <= bio_n + S_IN:
        X_ND = 0
        S_IN -= (xs_ndm - bio_n)
        bio_n = 0
    else:
    #     X_S = (bio_n + S_IN)/cmps_asm.X_S.i_N
    #     X_ND = 0
    #     S_IN, bio_n = 0
    # bio_cod -= (X_S + X_P)
    # if bio_cod > 0:
        raise RuntimeError('Not enough nitrogen (S_IN + biomass) to map '
                            'all biomass COD into X_P and X_S')
    
    # Step 1b: convert particulate substrates into X_S + X_ND
    xsub_cod = X_c + X_ch + X_pr + X_li
    xsub_n = X_c*cmps_adm.X_c.i_N + X_pr*cmps_adm.X_pr.i_N    
    X_S += xsub_cod
    X_ND += xsub_n - xsub_cod*cmps_asm.X_S.i_N  # X_S.i_N should technically be zero
    if X_ND < 0:
        raise RuntimeError('Not enough nitrogen (substrate + excess X_ND) '
                           'to map all particulate substrate COD into X_S')
    
    # Step 2: map all X_I from ADM to ASM
    excess_XIn = X_I * (cmps_adm.X_I.i_N - cmps_asm.X_I.i_N)
    S_IN += excess_XIn
    if S_IN < 0:
        raise RuntimeError('Not enough nitrogen (X_I + S_IN) to map '
                           'all ADM X_I into ASM X_I')
    
    # Step 3: map ADM S_I into ASM S_I and S_NH
    excess_SIn = S_I * (cmps_adm.S_I.i_N - cmps_asm.S_I.i_N)
    if excess_SIn > 0:
        S_NH = excess_SIn
    else:
        S_NH = 0
        S_IN += excess_SIn
        if S_IN < 0:
            raise RuntimeError('Not enough nitrogen (S_I + S_IN) to map '
                               'all ADM S_I into ASM S_I')
        
    # Step 4: map all soluble substrates into S_S and S_ND
    ssub_cod = S_su + S_aa + S_fa + S_va + S_bu + S_pro + S_ac
    ssub_n = S_aa * cmps_adm.S_aa.i_N
    if ssub_cod*cmps_asm.S_S.i_N <= ssub_n:
        S_S = ssub_cod
        S_ND = ssub_n - S_S/cmps_asm.S_S.i_N # S_S.i_N should technically be zero
    else:
        raise RuntimeError('Not enough nitrogen to map all soluble substrates into ASM S_S')
    
    # Step 5: charge balance for alkalinity
    pKa = calc_pKa(T)
    pKw, pKa_IN, pKa_IC = pKa[:3]
    pKa_vfa = pKa[3:]
    alpha_IN = 10**(pKa_IN-pH)/(1+10**(pKa_IN-pH))/14 # charge per g N
    alpha_IC = -1/(1+10**(pKa_IC-pH))/12
    alpha_vfa = calc_vfa_alpha(pKa_vfa, pH)
    S_ALK = (sum(_ions * np.append([alpha_IN, alpha_IC], alpha_vfa)) - S_NH/14)*(-12)
    
    # Step 6: check COD and TKN balance
    asm_vals = np.array(([S_I, S_S, X_I, X_S, X_BH, X_BA, X_P, S_O, S_NO, S_NH, S_ND, X_ND, S_ALK, S_N2, H2O]))
    assert sum(asm_vals*cmps_asm.i_COD) == sum(adm_vals*cmps_adm.i_COD)
    assert sum(asm_vals*cmps_asm.i_N) == sum(adm_vals*cmps_adm.i_N)
    
    # convert from kg/m3 (ADM) to mg/L(ASM)
    return asm_vals*1000