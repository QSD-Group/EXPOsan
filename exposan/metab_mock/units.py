# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 12:17:42 2022

@author: joy_c
"""
from qsdsan import SanUnit
from qsdsan.sanunits import AnaerobicCSTR, CSTR
import numpy as np
from warnings import warn

__all__ = ('DegassingMembrane',
           'METAB_AnCSTR')

#%%
class DegassingMembrane(SanUnit):
    
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=True,
                 tau=0.01,
                 H2_degas_efficiency=0.8, CH4_degas_efficiency=0.8, 
                 CO2_degas_efficiency=0.8, gas_IDs=('S_h2', 'S_ch4', 'S_IC')):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo,
                         init_with=init_with, F_BM_default=F_BM_default,
                         isdynamic=isdynamic)
        self.tau = tau
        self.H2_degas_efficiency = H2_degas_efficiency
        self.CH4_degas_efficiency = CH4_degas_efficiency
        self.CO2_degas_efficiency = CO2_degas_efficiency
        self.gas_IDs = gas_IDs
        self._split = np.zeros(len(self.components))
        self._gas_idx = self.components.indices(gas_IDs)
    
    @property
    def H2_degas_efficiency(self):
        return self._h2_ermv
    
    @H2_degas_efficiency.setter
    def H2_degas_efficiency(self, e):
        if e > 1 or e < 0:
            raise ValueError(f'degassing efficiency must be within [0, 1], not {e}')
        self._h2_ermv = e
        
    @property
    def CH4_degas_efficiency(self):
        return self._ch4_ermv
    
    @CH4_degas_efficiency.setter
    def CH4_degas_efficiency(self, e):
        if e > 1 or e < 0:
            raise ValueError(f'degassing efficiency must be within [0, 1], not {e}')
        self._ch4_ermv = e
    
    @property
    def CO2_degas_efficiency(self):
        return self._co2_ermv
    
    @CO2_degas_efficiency.setter
    def CO2_degas_efficiency(self, e):
        if e > 1 or e < 0:
            raise ValueError(f'degassing efficiency must be within [0, 1], not {e}')
        self._co2_ermv = e
    
    @property
    def split(self):
        s = self._split * 0
        s[self._gas_idx] = [self._h2_ermv, self._ch4_ermv, self._co2_ermv]
        return s   
    
    def _run(self):
        inf, = self.ins
        gas, liquid = self.outs
        s = self.split
        inf.split_to(gas, liquid, s)
        gas.phase = 'g'

    def _init_state(self):
        self._state = self._ins_QC[0]
        self._dstate = self._state * 0.
        
    def _update_state(self):
        arr = self._state
        gas, liquid = self.outs
        s = self.split
        Q_liq = arr[-1]
        if gas.state is None: gas.state = arr*0.0
        gas.state[:-1] = s * arr[:-1] * Q_liq
        gas.state[-1] = 1
        if liquid.state is None: liquid.state = arr*0.0    
        liquid.state[:-1] = (1-s) * arr[:-1]
        liquid.state[-1] = Q_liq

    def _update_dstate(self):
        arr = self._dstate
        gas, liquid = self.outs
        s = self.split
        Q_liq = self._state[-1]
        C_liq = self._state[:-1]
        if gas.dstate is None: gas.dstate = arr*0.0
        gas.dstate[:-1] = s * (arr[:-1] * Q_liq + C_liq * arr[-1])
        gas.dstate[-1] = 0
        if liquid.dstate is None: liquid.dstate = arr*0.0    
        liquid.dstate[:-1] = (1-s) * arr[:-1]
        liquid.dstate[-1] = arr[-1]
     
    # @property
    # def AE(self):
    #     if self._AE is None:
    #         self._compile_AE()
    #     return self._AE

    # def _compile_AE(self):
    #     _state = self._state
    #     _dstate = self._dstate
    #     _update_state = self._update_state
    #     _update_dstate = self._update_dstate
    #     def yt(t, QC_ins, dQC_ins):
    #         _state[:] = QC_ins[0]
    #         _dstate[:] = dQC_ins[0]
    #         _update_state()
    #         _update_dstate()
    #     self._AE = yt

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE
    
    def _compile_ODE(self):
        _dstate = self._dstate
        _update_dstate = self._update_dstate
        tau = self.tau
        def dy_dt(t, QC_ins, QC, dQC_ins):
            _dstate[:] = (QC_ins[0] - QC)/tau
            _dstate[-1] = dQC_ins[0,-1]
            _update_dstate()
        self._ODE = dy_dt

#%%

class METAB_AnCSTR(AnaerobicCSTR, CSTR):
    
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    _R = 8.3145e-2 # Universal gas constant, [bar/M/K]
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', V_liq=3400, V_gas=300, split=None,
                 model=None, T=308.15, max_headspace_P=30, 
                 retain_cmps=(), fraction_retain=0.95,
                 isdynamic=True, exogenous_vars=(), **kwargs):
        AnaerobicCSTR.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, 
                        init_with=init_with, V_liq=V_liq, V_gas=V_gas, model=model,  
                        T=T, retain_cmps=retain_cmps, fraction_retain=fraction_retain,
                        isdynamic=isdynamic, exogenous_vars=exogenous_vars, **kwargs)
        self.max_headspace_P = max_headspace_P
        self.split = split
        self._fixed_P_gas = None
        self.f_q_gas_fixed_P_headspace = self.f_q_gas_var_P_headspace = lambda *args: 0.
        self._biogas = None
    
    # split = property(CSTR.split.fget, CSTR.split.fset)
    split = CSTR.split
        
    headspace_P = property(AnaerobicCSTR.headspace_P.fget)
    @headspace_P.setter
    def headspace_P(self, P):
        pass
        
    @property
    def external_P(self):
        return None
    @external_P.setter
    def external_P(self, P):
        pass
    
    @property
    def pipe_resistance(self):
        return None
    @pipe_resistance.setter
    def pipe_resistance(self, k):
        pass

    @property
    def fixed_headspace_P(self):
        return None
    @fixed_headspace_P.setter
    def fixed_headspace_P(self, b):
        pass
    
    @property
    def max_headspace_P(self):
        '''Maximum headspace pressure [bar]'''
        return self._Pmax_hs
    @max_headspace_P.setter
    def max_headspace_P(self, p):
        p_vap = self.p_vapor()
        if p < p_vap:
            raise ValueError(f'max_headspace_P must at least be the saturated '
                             f'vapor pressure {p_vap} bar, not {p}')
        self._Pmax_hs = p
    
    _run = CSTR._run

    def _update_state(self):
        arr = self._state
        f_rtn = self._f_retain
        n_cmps = len(self.components)
        Cs = arr[:n_cmps]*(1-f_rtn)*1e3  # kg/m3 to mg/L
        S_gas = arr[n_cmps: (n_cmps+self._n_gas)]
        p_gas = sum(self._R * self.T * S_gas) + self.p_vapor()
        if p_gas> self._Pmax_hs:
            warn(f'headspace pressure is {p_gas} bar, exceeding design maximum {self._Pmax_hs}')
        if self.split is None:
            out, = self.outs
            if out.state is None:
                out.state = np.append(Cs, arr[-1])
            else: 
                out.state[:n_cmps] = Cs
                out.state[-1] = arr[-1]
        else:
            for ws, spl in zip(self._outs, self.split):
                if ws.state is None: ws.state = np.empty(n_cmps+1)
                ws.state[:n_cmps] = Cs
                ws.state[-1] = arr[-1]*spl

    def _update_dstate(self):
        arr = self._dstate
        f_rtn = self._f_retain
        n_cmps = len(self.components)
        dCs = arr[:n_cmps]*(1-f_rtn)*1e3  # kg/m3 to mg/L
        if self.split is None:
            out, = self.outs
            if out.dstate is None:
                out.dstate = np.append(dCs, arr[-1])
            else: 
                out.dstate[:n_cmps] = dCs
                out.dstate[-1] = arr[-1]
        else:
            for ws, spl in zip(self._outs, self.split):
                if ws.dstate is None: ws.dstate = np.empty(n_cmps+1)
                ws.dstate[:n_cmps] = dCs
                ws.dstate[-1] = arr[-1]*spl
    