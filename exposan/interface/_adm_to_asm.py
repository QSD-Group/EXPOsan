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
from . import Junction
from ._asm_to_adm import calc_pKa #!!! will it be different for asm/adm modules?

__all__ = ('ADMtoASM',)

#!!! these values should be retrieved from the ADM1 object's rate function parameter set
vfa_cod = np.array([64, 112, 160, 208])
def calc_vfa_alpha(pKa_vfa, pH):
    return 1.0/vfa_cod*(-1.0/(1.0 + 10**(pKa_vfa - pH)))

xs_xp_split_bio = [0.7, 0.3]

class ADMtoASM(Junction):
    '''
    Interface unit to convert anaerobic digestion model (ADM) components
    to activated sludge model (ASM) components.
    
    Parameters
    ----------
    upstream : stream or str
        Influent stream with ADM components.
    downstream : stream or str
        Effluent stream with ASM components.
    
    References
    ----------
    [1] Nopens, I.; Batstone, D. J.; Copp, J. B.; Jeppsson, U.; Volcke, E.; 
    Alex, J.; Vanrolleghem, P. A. An ASM/ADM Model Interface for Dynamic 
    Plant-Wide Simulation. Water Res. 2009, 43, 1913–1923.
    
    See Also
    --------
    :class:`qsdsan.sanunits.Junction`_
    
    :class:`qsdsan.sanunits.ASMtoADM`_    
    '''
    def __init__(self, ID='', upstream=None, downstream=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False):
        super().__init__(self, ID=ID, upstream=upstream, downstream=downstream,
                         init_with=init_with, F_BM_default=F_BM_default,
                         isdynamic=isdynamic)

    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        
        # Retrieve constants
        ins = self.ins[0]
        outs = self.outs[0]
        T = ins.T
        pH = ins.pH
        
        cmps_adm = ins.components
        X_c_i_N = cmps_adm.X_c.i_N
        X_pr_i_N = cmps_adm.X_pr.i_N
        S_aa_i_N = cmps_adm.S_aa.i_N
        adm_X_I_i_N = cmps_adm.X_I.i_N
        adm_S_I_i_N = cmps_adm.S_I.i_N
        adm_i_COD = cmps_adm.i_COD
        adm_i_N = cmps_adm.i_N
        adm_bio_N_indices = cmps_adm.indices(('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2'))

        cmps_asm = outs.components
        X_P_i_N = cmps_asm.X_P.i_N
        X_S_i_N = cmps_asm.X_S.i_N
        S_S_i_N = cmps_asm.S_S.i_N
        asm_X_I_i_N = cmps_asm.X_I.i_N
        asm_S_I_i_N = cmps_asm.S_I.i_N
        asm_i_COD = cmps_asm.i_COD
        asm_i_N = cmps_asm.i_N

        def adm2asm(adm_vals):    
            S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I, \
                X_c, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, \
                S_cat, S_an, H2O = adm_vals
                       
            # Step 0: snapshot of charged components
            _ions = np.array([S_IN, S_IC, S_ac, S_pro, S_bu, S_va])
            
            # Step 1a: convert biomass into X_S+X_ND and X_P
            bio_cod = X_su + X_aa + X_fa + X_c4 + X_pro + X_ac + X_h2
            bio_n = sum((adm_vals*adm_i_N)[adm_bio_N_indices])
            xs_cod, xp_cod = [frac*bio_cod for frac in xs_xp_split_bio]
            xp_ndm = xp_cod*X_P_i_N
            if xp_ndm > bio_n:
                X_P = bio_n/asm_i_N
                bio_n = 0
            else:
                X_P = xp_cod
                bio_n -= xp_ndm
            X_S = bio_cod - X_P
            xs_ndm = X_S*X_S_i_N
            if xs_ndm <= bio_n:
                X_ND = bio_n - xs_ndm
                bio_n = 0
            elif xs_ndm <= bio_n + S_IN:
                X_ND = 0
                S_IN -= (xs_ndm - bio_n)
                bio_n = 0
            else:
            #     X_S = (bio_n + S_IN)/X_S_i_N
            #     X_ND = 0
            #     S_IN, bio_n = 0
            # bio_cod -= (X_S + X_P)
            # if bio_cod > 0:
                raise RuntimeError('Not enough nitrogen (S_IN + biomass) to map '
                                    'all biomass COD into X_P and X_S')
            
            # Step 1b: convert particulate substrates into X_S + X_ND
            xsub_cod = X_c + X_ch + X_pr + X_li
            xsub_n = X_c*X_c_i_N + X_pr*X_pr_i_N
            X_S += xsub_cod
            X_ND += xsub_n - xsub_cod*X_S_i_N  # X_S.i_N should technically be zero
            if X_ND < 0:
                raise RuntimeError('Not enough nitrogen (substrate + excess X_ND) '
                                   'to map all particulate substrate COD into X_S')
            
            # Step 2: map all X_I from ADM to ASM           
            excess_XIn = X_I * (adm_X_I_i_N - asm_X_I_i_N)
            S_IN += excess_XIn
            if S_IN < 0:
                raise RuntimeError('Not enough nitrogen (X_I + S_IN) to map '
                                   'all ADM X_I into ASM X_I')
            
            # Step 3: map ADM S_I into ASM S_I and S_NH
            excess_SIn = S_I * (adm_S_I_i_N - asm_S_I_i_N)
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
            ssub_n = S_aa * S_aa_i_N
            if ssub_cod*S_S_i_N <= ssub_n:
                S_S = ssub_cod
                S_ND = ssub_n - S_S/S_S_i_N # S_S.i_N should technically be zero
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
            asm_vals = np.array(([
                S_I, S_S, X_I, X_S, 
                0, 0, # X_BH, X_BA, 
                X_P, 
                0, 0, # S_O, S_NO, 
                S_NH, S_ND, X_ND, S_ALK, 
                0, # S_N2, 
                H2O]))
            assert sum(asm_vals*asm_i_COD) == sum(adm_vals*adm_i_COD)
            assert sum(asm_vals*asm_i_N) == sum(adm_vals*adm_i_N)
            
            # # Aren't all conc in mg/L?
            # # convert from kg/m3 (ADM) to mg/L(ASM)
            # return asm_vals*1000
            return asm_vals
        
        def yt(t, QC_ins, dQC_ins):
            # X_BH, X_BA, S_O, S_NO, S_N2 = 0
            
            # asm_vals = np.array(([
            #     S_I, S_S, X_I, X_S, 
            #     0, 0, # X_BH, X_BA, 
            #     X_P, 
            #     0, 0, # S_O, S_NO, 
            #     S_NH, S_ND, X_ND, S_ALK, 
            #     0, # S_N2, 
            #     H2O]))
            
            for i, j in zip((QC_ins, dQC_ins), (_state, _dstate)):                             
                adm_vals = i[0][:-1] # shape = (1, num_upcmps)
                asm_vals = adm2asm(adm_vals)
                Q = asm_vals.sum() #!!! what's the unit of Q?
                j[:-1] = asm_vals
                j[-1] = Q

            _update_state()
            _update_dstate()
        self._AE = yt