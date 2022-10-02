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

import numpy as np, qsdsan as qs
from warnings import warn
from math import isclose
from biosteam.units import Junction as BSTjunction
from qsdsan import SanUnit, processes as pc

__all__ = (
    'Junction',
    'ADMjunction', 'ADMtoASM', 'ASMtoADM',
    )

class Junction(SanUnit):
    '''
    A non-reactive class that serves to convert a stream with one set of components
    and into another.
    
    Thermal conditions of the downstream (T, P) will be copied from those of the upstream.


    Parameters
    ----------
    upstream : stream or str
        Influent stream.
    downstream : stream or str
        Effluent stream.
    reactions : iterable(dict) | callable
        Iterable of dict that has the conversion of upstream components to
        downstream components,
        or a function that will return the concentration of the effluent
        with influent concentration as the input.
        
        If given as an iterable of dict, keys of each of the dict should be 
        the ID or alias of components,
        values should be the conversion/yield,
        which should be negative for reactants and positive for products.
    '''
    _graphics = BSTjunction._graphics

    def __init__(self, ID='', upstream=None, downstream=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 reactions=None, **kwargs):
        thermo = downstream.thermo if downstream else thermo
        SanUnit.__init__(self, ID, ins=upstream, outs=downstream, thermo=thermo,
                         init_with=init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic,
                         skip_property_package_check=True)
        if reactions: self.reactions = reactions
        else: self.reactions = None
        
        for key, val in kwargs.items():
            setattr(self, key, val)

    def _no_parse_reactions(self, rxns):
        if rxns is None: return
        raise RuntimeError('Reactions are automatically compiled.')


    def _parse_reactions(self, rxns):
        cmps_in = self.ins[0].components
        cmps_outs = self.outs[0].components

        num_rxns = len(rxns)
        num_upcmps = len(cmps_in)
        num_downcmps = len(cmps_outs)
        RX = np.zeros(shape=(num_rxns, num_upcmps)) # each row is a reaction
        RY = np.zeros(shape=(num_rxns, num_downcmps))

        for n, dct in enumerate(rxns):
            RXn, RYn = RX[n], RY[n]
            for cmp, val in dct.items():
                if val < 0: RXn[cmps_in.index(cmp)] = val
                elif val > 0: RYn[cmps_outs.index(cmp)] = val

        # Transfer overlapped components
        overlapped = set(cmps_in.IDs).intersection(set(cmps_outs.IDs))
        RXsum = RX.sum(axis=0)
        RXnew = []
        RYnew = []
        for ID in overlapped:
            idx_up = cmps_in.index(ID)
            idx_down= cmps_outs.index(ID)
            if RXsum[idx_up] != -1:
                newRXn = np.zeros(num_upcmps)
                newRYn = np.zeros(num_downcmps)
                newRXn[idx_up] = -1 - RXsum[idx_up]
                newRYn[idx_down] = -newRXn[idx_up]
                RXnew.append(newRXn)
                RYnew.append(newRYn)
        RX = np.concatenate((RX, np.array(RXnew)))
        RY = np.concatenate((RY, np.array(RYnew)))

        # Check if all upstream components are converted
        RXsum = RX.sum(axis=0)
        if np.where(RXsum!=-1)[0].any():
            index = np.where(RXsum!=-1)[0][0]
            raise ValueError(f'Conversion for Component "{cmps_in.IDs[index]}" '
                             f'is {abs(RXsum[index])}, not 100%.')
        self._RX = RX
        self._RY = RY
        
        
    def _compile_reactions(self):
        def reactions(X):
            X.reshape(1, len(X))
            Yarr = -(self._RX*X).T @ self._RY # _RX: (num_rxns, num_upcmps); _RY: (num_rxns, num_downcmps)
            Y = Yarr.sum(axis=0) # Yarr: (num_upcmps, num_downcmps)
            return Y.reshape(Y.shape[1],)
        self._reactions = reactions
        

    def _run(self):
        ins = self.ins[0]
        rxns = self.reactions
        X = ins.conc.value
        Y = rxns(X)
        outs = self.outs[0]
        outs.thermal_condition.copy_like(ins.thermal_condition)
        outs.set_flow_by_concentration(
            flow_tot=ins.F_vol,
            concentrations=dict(zip(outs.components.IDs, Y)),
            units=('m3/hr', 'mg/L'))


    # Below are dynamic simulation-related properties
    @property
    def state(self):
        '''The state of the Junction, including component concentrations [mg/L] and flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    def _init_dynamic(self):
        super()._init_dynamic()
        # Need to use ins' components, otherwise _ins_QC will follow the shape of
        # the unit's (i.e., downstream) components
        self._ins_QC = np.zeros((len(self._ins), len(self.ins[0].components)+1))
        self._ins_dQC = self._ins_QC.copy()

    def _init_state(self):
        '''
        Initialize state by specifying or calculating component concentrations
        based on influents. Total flow rate is always initialized as the sum of
        influent wastestream flows.
        '''
        self._state = np.append(self.outs[0].conc, self.outs[0].F_vol*24)
        self._dstate = self._state * 0.

    def _update_state(self):
        '''
        Updates conditions of output stream based on conditions of the Junction.
        '''
        self._outs[0].state = self._state

    def _update_dstate(self):
        '''
        Updates rates of change of output stream from rates of change of the Junction.
        '''
        self._outs[0].dstate = self._dstate

    # The unit's state should be the same as the effluent state
    # react the state arr and dstate arr
    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        rxns = self.reactions
        def yt(t, QC_ins, dQC_ins):
            for i, j in zip((QC_ins, dQC_ins), (_state, _dstate)):
                X = i[0][:-1] # shape = (1, num_upcmps)
                Y = rxns(X)
                j[:-1] = Y
                j[-1] = i[0][-1] # volumetric flow of outs should equal that of ins
            _update_state()
            _update_dstate()
        self._AE = yt

    
    @property
    def upstream(self):
        '''[qsdsan.WasteStream] Influent.'''
        return self.ins[0]
    @upstream.setter
    def upstream(self, upstream):
        self.ins[0] = upstream

    @property
    def downstream(self):
        '''[qsdsan.WasteStream] Effluent.'''
        return self.outs[0]
    @downstream.setter
    def downstream(self, downstream):
        self.outs[0] = downstream
    
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE
    
    @property
    def reactions(self):
        '''
        [callable] Function that takes the concentration array of the influent
        and convert to the concentration array of the effluent.
        '''
        if not self._reactions: self._compile_reactions()
        return self._reactions
    @reactions.setter
    def reactions(self, i):
        if callable(i): self._reactions = i
        else:
            self._parse_reactions(i)
            self._compile_reactions()
            
            
# %%

#TODO: add a `rtol` kwargs for error checking
class ADMjunction(Junction):
    '''
    An abstract superclass holding common properties of ADM interface classes.
    Users should use its subclasses (e.g., ``ASMtoADM``, ``ADMtoASM``) instead.
    
    See Also
    --------
    :class:`qsdsan.sanunits.Junction`
    
    :class:`qsdsan.sanunits.ADMtoASM`
    
    :class:`qsdsan.sanunits.ASMtoADM`
    '''
    _parse_reactions = Junction._no_parse_reactions
    tolerance = 1e-3
    
    def __init__(self, ID='', upstream=None, downstream=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 adm1_model=None):
        self.adm1_model = adm1_model # otherwise there won't be adm1_model when `_compile_reactions` is called
        if thermo is None:
            warn('No `thermo` object is provided and is prone to raise error. '
                 'If you are not sure how to get the `thermo` object, '
                 'use `thermo = qsdsan.set_thermo` after setting thermo with the `Components` object.')
        super().__init__(ID=ID, upstream=upstream, downstream=downstream,
                         thermo=thermo, init_with=init_with, 
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        
   
    @property
    def T(self):
        '''[float] Temperature of the upstream/downstream [K].'''
        return self.ins[0].T
    @T.setter
    def T(self, T):
        self.ins[0].T = self.outs[0].T = T
    
    @property
    def pH(self):
        '''[float] pH of the upstream/downstream.'''
        return self.ins[0].pH
    
    @property
    def adm1_model(self):
        '''[qsdsan.Process] ADM process model.'''
        return self._adm1_model
    @adm1_model.setter
    def adm1_model(self, model):
        if not isinstance(model, pc.ADM1):
            raise ValueError('`adm1_model` must be an `AMD1` object, '
                             f'the given object is {type(model).__name__}.')
        self._adm1_model = model
        
    @property
    def T_base(self):
        '''[float] Base temperature in the ADM1 model.'''
        return self.adm1_model.rate_function.params['T_base']
    
    @property
    def pKa_base(self):
        '''[float] pKa of the acid-base pairs at the base temperature in the ADM1 model.'''
        Ka_base = self.adm1_model.rate_function.params['Ka_base']
        return -np.log10(Ka_base)
    
    @property
    def Ka_dH(self):
        '''[float] Heat of reaction for Ka.'''
        return self.adm1_model.rate_function.params['Ka_dH']
    
    @property
    def pKa(self):
        '''
        [numpy.array] pKa array of the following acid-base pairs:
        ('H+', 'OH-'), ('NH4+', 'NH3'), ('CO2', 'HCO3-'),
        ('HAc', 'Ac-'), ('HPr', 'Pr-'), ('HBu', 'Bu-'), ('HVa', 'Va-')
        '''
        return self.pKa_base-np.log10(np.exp(pc.T_correction_factor(self.T_base, self.T, self.Ka_dH)))
    
    @property
    def alpha_IC(self):
        '''[float] Charge per g of C.'''
        pH = self.pH
        pKa_IC = self.pKa[2]
        return -1/(1+10**(pKa_IC-pH))/12

    @property
    def alpha_IN(self):
        '''[float] Charge per g of N.'''
        pH = self.pH
        pKa_IN = self.pKa[1]
        return 10**(pKa_IN-pH)/(1+10**(pKa_IN-pH))/14
    
    
# %%

class ADMtoASM(ADMjunction):
    '''
    Interface unit to convert anaerobic digestion model (ADM) components
    to activated sludge model (ASM) components.
    
    Parameters
    ----------
    upstream : stream or str
        Influent stream with ADM components.
    downstream : stream or str
        Effluent stream with ASM components.
    adm1_model : obj
        The anaerobic digestion process model (:class:`qsdsan.processes.ADM1`).
    bio_to_xs : float
        Split of the total biomass COD to slowly biodegradable substrate (X_S),
        the rest is assumed to be mapped into X_P.
    tolerance : float
        Error tolerance,
        relative when checking mass balance,
        absolute when checking material abundance.
    
    References
    ----------
    [1] Nopens, I.; Batstone, D. J.; Copp, J. B.; Jeppsson, U.; Volcke, E.; 
    Alex, J.; Vanrolleghem, P. A. An ASM/ADM Model Interface for Dynamic 
    Plant-Wide Simulation. Water Res. 2009, 43, 1913–1923.
    
    See Also
    --------
    :class:`qsdsan.sanunits.ADMjunction`
    
    :class:`qsdsan.sanunits.ASMtoADM`  
    '''
    # User defined values
    bio_to_xs = 0.7
    
    # Should be constants
    cod_vfa = np.array([64, 112, 160, 208])
        
    def _compile_reactions(self):
        # Retrieve constants
        ins = self.ins[0]
        outs = self.outs[0]
        tol = self.tolerance
        
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
        
        alpha_IN = self.alpha_IN
        alpha_IC = self.alpha_IC
        alpha_vfa = self.alpha_vfa

        def adm2asm(adm_vals):    
            S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I, \
                X_c, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, \
                S_cat, S_an, H2O = adm_vals
                       
            # Step 0: snapshot of charged components
            _ions = np.array([S_IN, S_IC, S_ac, S_pro, S_bu, S_va])
            
            # Step 1a: convert biomass into X_S+X_ND and X_P
            bio_cod = X_su + X_aa + X_fa + X_c4 + X_pro + X_ac + X_h2
            bio_n = sum((adm_vals*adm_i_N)[adm_bio_N_indices])
            xp_cod = bio_cod * (1-self.bio_to_xs)
            xp_ndm = xp_cod*X_P_i_N
            if xp_ndm > bio_n:
                warn('Not enough biomass N to map the specified proportion of '
                     'biomass COD into X_P. Mapped as much COD as possible, the rest '
                     'goes to X_S.')
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
            elif xs_ndm*(1+tol) <= bio_n + S_IN:
                X_ND = 0
                S_IN -= (xs_ndm - bio_n)
                bio_n = 0
            else:
                raise RuntimeError('Not enough nitrogen (S_IN + biomass) to map '
                                    'all biomass COD into X_P and X_S')
            
            # Step 1b: convert particulate substrates into X_S + X_ND
            xsub_cod = X_c + X_ch + X_pr + X_li
            xsub_n = X_c*X_c_i_N + X_pr*X_pr_i_N
            X_S += xsub_cod
            X_ND += xsub_n - xsub_cod*X_S_i_N  # X_S.i_N should technically be zero
            if X_ND+tol < 0:
                raise RuntimeError('Not enough nitrogen (substrate + excess X_ND) '
                                   'to map all particulate substrate COD into X_S')
            else: X_ND = max(0, X_ND) # reset small error
            
            # Step 2: map all X_I from ADM to ASM           
            excess_XIn = X_I * (adm_X_I_i_N - asm_X_I_i_N)
            S_IN += excess_XIn
            if S_IN+tol < 0:
                breakpoint()
                raise RuntimeError('Not enough nitrogen (X_I + S_IN) to map '
                                   'all ADM X_I into ASM X_I')
            else: S_IN = max(0, S_IN) # reset small error
            
            # Step 3: map ADM S_I into ASM S_I and S_NH
            excess_SIn = S_I * (adm_S_I_i_N - asm_S_I_i_N)
            if excess_SIn > 0:
                S_NH = excess_SIn
            else:
                S_NH = 0
                S_IN += excess_SIn
                if S_IN+tol < 0:
                    raise RuntimeError('Not enough nitrogen (S_I + S_IN) to map '
                                       'all ADM S_I into ASM S_I')
                else: S_IN = max(0, S_IN) # reset small error
            S_NH += S_IN
                
            # Step 4: map all soluble substrates into S_S and S_ND        
            ssub_cod = S_su + S_aa + S_fa + S_va + S_bu + S_pro + S_ac
            ssub_n = S_aa * S_aa_i_N
            if ssub_n<0 and abs(ssub_n) < tol: ssub_n = 0
            if (ssub_cod*S_S_i_N)*(1+tol) <= ssub_n:
                S_S = ssub_cod
                S_ND = ssub_n
                if S_S_i_N != 0:
                    S_ND -= S_S/S_S_i_N # S_S.i_N should technically be zero
            else:
                raise RuntimeError('Not enough nitrogen to map all soluble substrates into ASM S_S')
            
            # Step 5: charge balance for alkalinity
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
            
            if S_h2 > tol or S_ch4 > tol:
                warn('Ignored dissolved H2 or CH4.')
            
            lhs = sum(adm_vals*adm_i_COD) - S_h2 - S_ch4
            rhs = sum(asm_vals*asm_i_COD)
            if not isclose(lhs, rhs, rel_tol=tol):
                raise RuntimeError('COD not balanced, '
                                   f'influent (ADM) COD is {lhs}, '
                                   f'effluent (ASM) COD is {rhs}.')

            lhs = sum(adm_vals*adm_i_N)
            rhs = sum(asm_vals*asm_i_N)
            if not isclose(lhs, rhs, rel_tol=tol):
                raise RuntimeError('TKN not balanced, '
                                   f'influent (ASM) TKN is {lhs}, '
                                   f'effluent (ADM) TKN is {rhs}.')
            
            return asm_vals
        
        self._reactions = adm2asm
    
    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        adm2asm = self.reactions
               
        def yt(t, QC_ins, dQC_ins):          
            for i, j in zip((QC_ins, dQC_ins), (_state, _dstate)):                             
                adm_vals = i[0][:-1] # shape = (1, num_upcmps)
                asm_vals = adm2asm(adm_vals)
                j[:-1] = asm_vals
                j[-1] = i[0][-1] # volumetric flow of outs should equal that of ins

            _update_state()
            _update_dstate()
        
        self._AE = yt

    @property
    def alpha_vfa(self):
        return 1.0/self.cod_vfa*(-1.0/(1.0 + 10**(self.pKa[3:]-self.pH)))
        
        
# %%

class ASMtoADM(ADMjunction):
    '''
    Interface unit to convert activated sludge model (ASM) components
    to anaerobic digestion model (ADM) components.
    
    Parameters
    ----------
    upstream : stream or str
        Influent stream with ASM components.
    downstream : stream or str
        Effluent stream with ADM components.
    adm1_model : obj
        The anaerobic digestion process model (:class:`qsdsan.processes.ADM1`).
    xs_to_li : float
        Split of slowly biodegradable substrate COD to lipid, 
        after all N is mapped into protein.
    bio_to_li : float
        Split of biomass COD to lipid, after all biomass N is
        mapped into protein.
    frac_deg : float
        Biodegradable fraction of biomass COD.
    tolerance : float
        Error tolerance,
        relative when checking mass balance,
        absolute when checking material abundance.
    
    References
    ----------
    [1] Nopens, I.; Batstone, D. J.; Copp, J. B.; Jeppsson, U.; Volcke, E.; 
    Alex, J.; Vanrolleghem, P. A. An ASM/ADM Model Interface for Dynamic 
    Plant-Wide Simulation. Water Res. 2009, 43, 1913–1923.
    
    See Also
    --------
    :class:`qsdsan.sanunits.ADMjunction`
    
    :class:`qsdsan.sanunits.ADMtoASM`  
    '''    
    # User defined values
    xs_to_li = 0.7
    bio_to_li = 0.4
    frac_deg = 0.68
    
    def _compile_reactions(self):
        # Retrieve constants
        ins = self.ins[0]
        outs = self.outs[0]
        tol = self.tolerance

        cmps_asm = ins.components
        S_NO_i_COD = cmps_asm.S_NO.i_COD
        X_BH_i_N = cmps_asm.X_BH.i_N
        X_BA_i_N = cmps_asm.X_BA.i_N
        asm_X_I_i_N = cmps_asm.X_I.i_N
        X_P_i_N = cmps_asm.X_P.i_N
        asm_i_COD = cmps_asm.i_COD
        asm_i_N = cmps_asm.i_N
        asm_nonTKN_indices = cmps_asm.indices(('S_NO', 'S_N2'))
        
        cmps_adm = outs.components
        S_aa_i_N = cmps_adm.S_aa.i_N
        X_pr_i_N = cmps_adm.X_pr.i_N
        S_I_i_N = cmps_adm.S_I.i_N
        adm_X_I_i_N = cmps_adm.X_I.i_N
        adm_i_COD = cmps_adm.i_COD
        adm_i_N = cmps_adm.i_N
        
        frac_deg = self.frac_deg
        alpha_IN = self.alpha_IN
        alpha_IC = self.alpha_IC
        proton_charge = 10**(-self.pKa[0]+self.pH) - 10**(-self.pH) # self.pKa[0] is pKw

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
                S_S = X_S = X_BH = X_BA = 0
            elif cod_spl <= O_coddm + NO_coddm:
                S_O = 0
                S_NO = -(O_coddm + NO_coddm - cod_spl)/S_NO_i_COD
                S_S = X_S = X_BH = X_BA = 0
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
                S_O = S_NO = 0
            
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

            # Step 3: convert slowly biodegradable COD and TKN
            # into proteins, lipids, and carbohydrates
            req_xcod = X_ND / X_pr_i_N
            if X_S < req_xcod:
                X_pr = X_S
                X_li = X_ch = 0
                X_ND -= X_pr * X_pr_i_N
            else:
                X_pr = req_xcod
                X_li = self.xs_to_li * (X_S - X_pr)
                X_ch = (X_S - X_pr) - X_li
                X_ND = 0
            
            # Step 4: convert active biomass into protein, lipids, 
            # carbohydrates and potentially particulate TKN
            available_bioN = (
                X_BH * X_BH_i_N
                + X_BA * X_BA_i_N
                - (X_BH+X_BA) * (1-frac_deg) * adm_X_I_i_N
                )
            if available_bioN+tol < 0:
                raise RuntimeError('Not enough N in X_BA and X_BH to fully convert the non-biodegrable '
                                   'portion into X_I in ADM1.')
            req_bioN = (X_BH+X_BA) * frac_deg * X_pr_i_N
            if available_bioN + X_ND >= req_bioN:
                X_pr += (X_BH+X_BA) * frac_deg
                X_ND += available_bioN - req_bioN
            else:
                bio2pr = (available_bioN + X_ND)/X_pr_i_N
                X_pr += bio2pr
                bio_to_split = (X_BH+X_BA) * frac_deg - bio2pr
                bio_split_to_li = bio_to_split * self.bio_to_li
                X_li += bio_split_to_li
                X_ch += (bio_to_split - bio_split_to_li)
                X_ND = 0
            
            # Step 5: map particulate inerts
            xi_nsp = X_P_i_N * X_P + asm_X_I_i_N * X_I
            xi_ndm = (X_P+X_I) * adm_X_I_i_N
            if (xi_nsp + X_ND)*(1+tol) < xi_ndm: # allow for minor rounding error
                raise RuntimeError('Not enough N in X_I, X_P, X_ND to fully convert X_I and X_P '
                                   'into X_I in ADM1.')
            deficit = xi_ndm - xi_nsp
            X_I += X_P + (X_BH+X_BA) * (1-frac_deg)
            X_ND -= deficit            

            req_sn = S_I * S_I_i_N
            if req_sn <= S_ND:
                S_ND -= req_sn
            elif req_sn <= S_ND + X_ND:
                X_ND -= (req_sn - S_ND)
                S_ND = 0
            elif req_sn*(1+tol) <= S_ND + X_ND + S_NH:
                S_NH -= (req_sn - S_ND - X_ND)
                S_ND = X_ND = 0
            else:
                warn('Additional soluble inert COD is mapped to S_su.')
                SI_cod = (S_ND + X_ND + S_NH)/S_I_i_N
                S_su += S_I - SI_cod
                S_I = SI_cod
                S_ND = X_ND = S_NH = 0
                
            # Step 6: map any remaining nitrogen
            S_IN = S_ND + X_ND + S_NH
            
            # Step 7: charge balance
            asm_charge_tot = _snh/14 - _sno/14 - _salk/12
            #!!! charge balance should technically include VFAs, 
            # but VFAs concentrations are assumed zero per previous steps??
            S_IC = (asm_charge_tot -  S_IN*alpha_IN)/alpha_IC
            net_Scat = asm_charge_tot + proton_charge
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
            
            lhs = sum(asm_vals*asm_i_COD)
            rhs = sum(adm_vals*adm_i_COD)
            if not isclose(lhs, rhs, rel_tol=tol):
                raise RuntimeError('COD not balanced, '
                                   f'influent (ASM) COD is {lhs}, '
                                   f'effluent (ADM) COD is {rhs}.')
                
            lhs = sum(asm_vals*asm_i_N) - sum(asm_vals[asm_nonTKN_indices])
            rhs = sum(adm_vals*adm_i_N)
            if not isclose(lhs, rhs, rel_tol=tol):
                raise RuntimeError('TKN not balanced, '
                                   f'influent (ASM) TKN is {lhs}, '
                                   f'effluent (ADM) TKN is {rhs}.')
            return adm_vals
        
        self._reactions = asm2adm
        

    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        asm2adm = self.reactions
        
        def yt(t, QC_ins, dQC_ins):
            for i, j in zip((QC_ins, dQC_ins), (_state, _dstate)):   
                asm_vals = i[0][:-1] # shape = (1, num_upcmps)
                adm_vals = asm2adm(asm_vals)
                j[:-1] = adm_vals                
                j[-1] = i[0][-1] # volumetric flow of outs should equal that of ins

            _update_state()
            _update_dstate()
        
        self._AE = yt        
    