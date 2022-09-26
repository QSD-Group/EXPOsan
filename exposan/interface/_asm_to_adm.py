# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
     
    Joy Zhang <joycheung1994@gmail.com>
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from warnings import warn
from qsdsan import processes as pc
from . import Junction

__all__ = ('ASMtoADM',)

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

def calc_pKa(T):
    return pKa_base - np.log10(np.exp(pc.T_correction_factor(T_base, T, Ka_dH)))


class ASMtoADM(Junction):
    '''
    Interface unit to convert activated sludge model (ASM) components
    to anaerobic digestion model (ADM) components.
    
    Parameters
    ----------
    upstream : stream or str
        Influent stream with ASM components.
    downstream : stream or str
        Effluent stream with ADM components.
    
    References
    ----------
    [1] Nopens, I.; Batstone, D. J.; Copp, J. B.; Jeppsson, U.; Volcke, E.; 
    Alex, J.; Vanrolleghem, P. A. An ASM/ADM Model Interface for Dynamic 
    Plant-Wide Simulation. Water Res. 2009, 43, 1913–1923.
    
    See Also
    --------
    :class:`qsdsan.sanunits.Junction`_
    
    :class:`qsdsan.sanunits.ADMtoASM`_    
    '''   
    def _parse_reactions(self, rxns):
        raise RuntimeError('Reactions are automatically compiled.')
    
    def _compile_reactions(self):
        # Retrieve constants
        ins = self.ins[0]
        outs = self.outs[0]
        pH = ins.pH
        
        cmps_asm = ins.components
        S_NO_i_COD = cmps_asm.S_NO.i_COD
        X_BH_i_N = cmps_asm.X_BH.i_N
        X_BA_i_N = cmps_asm.X_BA.i_N
        asm_X_I_i_N = cmps_asm.X_I.i_N
        X_P_i_N = cmps_asm.X_P.i_N
        asm_i_COD = cmps_asm.i_COD
        asm_i_N = cmps_asm.i_N
        asm_N_gas_indices = cmps_asm.indices(('S_NO', 'S_N2'))
        
        cmps_adm = outs.components
        S_aa_i_N = cmps_adm.S_aa.i_N
        X_pr_i_N = cmps_adm.X_pr.i_N
        S_I_i_N = cmps_adm.S_I.i_N
        adm_X_I_i_N = cmps_adm.X_I.i_N
        adm_i_COD = cmps_adm.i_COD
        adm_i_N = cmps_adm.i_N
        
        def asm2adm(asm_vals):
            S_I, S_S, X_I, X_S, X_BH, X_BA, X_P, S_O, S_NO, S_NH, S_ND, X_ND, S_ALK, S_N2, H2O = asm_vals
            
            # Step 0: charged component snapshot
            _sno = S_NO
            _snh = S_NH
            _salk = S_ALK
            
            # Step 1: remove any remaining COD demand
            O_coddm = S_O
            NO_coddm = -S_NO*S_NO_i_COD
            cod_spl = S_S + X_S + X_BH + X_BA
            
            if cod_spl <= O_coddm:
                S_O = O_coddm - cod_spl
                S_S, X_S, X_BH, X_BA = 0
            elif cod_spl <= O_coddm + NO_coddm:
                S_O = 0
                S_NO = -(O_coddm + NO_coddm - cod_spl)/S_NO_i_COD
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
            req_scod = S_ND / S_aa_i_N
            if S_S < req_scod:
                S_aa = S_S
                S_su = 0
                S_ND -= S_aa * S_aa_i_N
            else:
                S_aa = req_scod
                S_su = S_S - S_aa
                S_ND = 0
            S_S = 0
            
            # Step 3: convert slowly biodegradable COD and TKN
            # into proteins, lipids, and carbohydrates
            req_xcod = X_ND / X_pr_i_N
            if X_S < req_xcod:
                X_pr = X_S
                X_li, X_ch = 0
                X_ND -= X_pr * X_pr_i_N
            else:
                X_pr = req_xcod
                X_li, X_ch = [frac*(X_S - X_pr) for frac in li_ch_split_XS]
                X_ND = 0
            X_S = 0
            
            # Step 4: convert active biomass into protein, lipids, 
            # carbohydrates and potentially particulate TKN
            available_bioN = X_BH * X_BH_i_N \
                + X_BA * X_BA_i_N \
                - (X_BH+X_BA) * (1-frac_deg) * adm_X_I_i_N
            if available_bioN < 0:
                raise RuntimeError('Not enough N in X_BA and X_BH to fully convert the non-biodegrable'
                                   'portion into X_I in ADM1.')
            req_bioN = (X_BH+X_BA) * frac_deg * X_pr_i_N
            if available_bioN + X_ND >= req_bioN:
                X_pr += (X_BH+X_BA) * frac_deg
                X_ND += available_bioN - req_bioN
            else:
                bio2pr = (available_bioN + X_ND)/X_pr_i_N
                X_pr += bio2pr
                X_li += ((X_BH+X_BA) * frac_deg - bio2pr) * li_ch_split_bio[0]
                X_ch += ((X_BH+X_BA) * frac_deg - bio2pr) * li_ch_split_bio[1]
                X_ND = 0
            X_BH, X_BA = 0
            
            # Step 5: map particulate inerts
            if X_P_i_N * X_P + asm_X_I_i_N * X_I + X_ND < (X_P+X_I) * adm_X_I_i_N:
                raise RuntimeError('Not enough N in X_I, X_P, X_ND to fully convert X_I and X_P'
                                   'into X_I in ADM1.')
            deficit = (X_P+X_I) * adm_X_I_i_N - X_P_i_N * X_P + asm_X_I_i_N * X_I
            X_I = X_I + X_P + (X_BH+X_BA) * (1-frac_deg)
            X_ND -= deficit
            

            req_sn = S_I * S_I_i_N
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
                SI_cod = (S_ND + X_ND + S_NH)/S_I_i_N
                S_su += S_I - SI_cod
                S_I = SI_cod
                S_ND, X_ND, S_NH = 0
                
            # Step 6: maps any remaining nitrogen
            S_IN = S_ND + X_ND + S_NH
            
            # Step 7: charge balance
            asm_charge_tot = _snh/14 - _sno/14 - _salk/12
            pKw, pKa_IN, pKa_IC = calc_pKa()[:2]
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
            adm_vals = np.array([
                S_su, S_aa,
                0, 0, 0, 0, 0, # S_fa, S_va, S_bu, S_pro, S_ac, 
                0, 0, # S_h2, S_ch4,
                S_IC, S_IN, S_I, 
                0, # X_c, 
                X_ch, X_pr, X_li, 
                0, 0, 0, 0, 0, 0, 0, # X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2,
                X_I, S_cat, S_an, H2O])
                      
            assert sum(asm_vals*asm_i_COD) == sum(adm_vals*adm_i_COD), 'COD not balanced.'
            assert sum(asm_vals*asm_i_N) - sum(asm_vals[asm_N_gas_indices]) \
                == sum(adm_vals*adm_i_N), 'N not balanced.'
            
            # # Aren't all conc in mg/L?
            # # unit conversion from mg/L (ASM) to kg/m3 (ADM)
            # return adm_vals/1000
            return adm_vals
        
        self._reactions = asm2adm
        

    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        asm2adm = self.reactions
        
        def yt(t, QC_ins, dQC_ins):
            # S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, \
            #     X_c, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2 = 0
            
            # adm_vals = np.array([
            #     S_su, S_aa,
            #     0, 0, 0, 0, 0, # S_fa, S_va, S_bu, S_pro, S_ac, 
            #     0, 0, # S_h2, S_ch4,
            #     S_IC, S_IN, S_I, 
            #     0, # X_c, 
            #     X_ch, X_pr, X_li, 
            #     0, 0, 0, 0, 0, 0, 0, # X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2,
            #     X_I, S_cat, S_an, H2O])
            
            for i, j in zip((QC_ins, dQC_ins), (_state, _dstate)):                             
                asm_vals = i[0][:-1] # shape = (1, num_upcmps)
                adm_vals = asm2adm(asm_vals)
                Q = adm_vals.sum() #!!! what's the unit of Q?
                j[:-1] = adm_vals
                j[-1] = Q

            _update_state()
            _update_dstate()
        
        self._AE = yt