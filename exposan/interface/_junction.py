# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from biosteam.units import Junction as BSTjunction
from qsdsan import SanUnit

__all__ = ('Junction',)

class Junction(SanUnit):
    '''
    A non-reactive class that serves to convert a stream with one set of components
    and into another.

    (i.e., all the outs at the same as the ins) unit that
    is used in dynamic simulation to record the unit/stream states.


    Parameters
    ----------
    upstream : stream or str
        Influent stream.
    downstream : stream or str
        Effluent stream.
    conversions : iterable(dict) | callable
        Iterable of dict that has the conversion of upstream components to
        downstream components,
        or a function that will return the concentration of the effluent
        with influent concentration as the input.
        
        If given as an iterable of dict, keys of each of the dict should be 
        the ID or alias of components,
        values should be the conversion/yield,
        which should be negative for reactants and positive for products.

    .. note::
        `ins` and `outs` have all the components,
        `upstream` and `downstream` only have their individual components;
        but the mass flow of `ins` equals that of `upstream`,
        and the mass flow of `outs` equals that of `downstream`.
    '''
    _graphics = BSTjunction._graphics

    def __init__(self, ID='', upstream=None, downstream=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 reactions=None):
        thermo = downstream.thermo if downstream else thermo
        SanUnit.__init__(self, ID, ins=upstream, outs=downstream, thermo=thermo,
                         init_with=init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic,
                         skip_property_package_check=True)
        if reactions: self.reactions = reactions
        else: self.reactions = None

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
        self._state = np.zeros(len(self.components)+1)
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

# =============================================================================
# Legacy codes
# =============================================================================

    # def __init__(self, ID='', upstream=None, outs=None, downstream=None, isdynamic=False,
    #              reactions=()):
    #     self._register(ID)
    #     self._init_specification()
    #     self._upstream = upstream
    #     self._downstream = downstream
    #     self.isdynamic = isdynamic
    #     self.reactions = reactions

    #     cmps = self._compile_components(upstream.components, downstream.components)
    #     thermo = self._load_thermo(cmps)
    #     self._ins = Inlets(self, 1, upstream, thermo, True, self._stacklevel)
    #     self._outs = Outlets(self, 1, downstream, thermo, True, self._stacklevel)
    #     self._update_streams(upstream=upstream, downstream=downstream)

    # def _compile_components(self, cmps1=None, cmps2=None):
    #     '''Compile components from the influent and effluent streams into one.'''
    #     cmps1 = cmps1 or self.upstream.components
    #     cmps2 = cmps2 or self.downstream.components
    #     cmps = qs.Components((*cmps1, *cmps2))
    #     cmps.compile()
    #     self._components = cmps
    #     return cmps


    # @ignore_docking_warnings
    # def _update_streams(self, upstream=None, downstream=None):
    #     '''
    #     Update `upstream` and `downstream` to have the same flows as `ins` and `outs`,
    #     but with different components.
    #     '''
    #     upstream = upstream or self.upstream
    #     downstream = downstream or self.downstream

    #     cmps = self._compile_components(cmps1=upstream.components)
    #     qs.set_thermo(cmps)
    #     self._reset_thermo(cmps)

    #     ins = qs.WasteStream(ID=self._ins[0].ID)
    #     ins.copy_like(upstream)
    #     self._ins[0] = ins

    #     outs = qs.WasteStream(ID=self._outs[0].ID)
    #     outs.copy_like(ins)
    #     self.reactions(outs.mol)
    #     self._outs[0] = outs

    #     downstream = self.downstream
    #     downstream._thermal_condition = outs._thermal_condition
    #     for ID in outs.components.IDs:
    #         downstream.imass[ID] = outs.imass[ID]


    # @property
    # def upstream(self):
    #     return self._upstream
    # @upstream.setter
    # def upstream(self, upstream):
    #     if self._upstream.thermo is not upstream.thermo:
    #         self._update_stream(upstream=upstream)
    #     self._upstream = upstream

    # @property
    # def downstream(self):
    #     return self._downstream
    # @downstream.setter
    # def downstream(self, downstream):
    #     if self._downstream.thermo is not downstream.thermo:
    #         self._update_stream(downstream=downstream)
    #     self._downstream = downstream