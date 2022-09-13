# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from biosteam.units import Junction as BSTjunction
from qsdsan import SanUnit, System
from qsdsan.sanunits import DynamicInfluent

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
    reaction : Iterable(dict)
        Iterable of dict that has the conversion of upstream components to
        downstream components.
        Keys of the dict should be the ID or alias of components,
        values should be the conversion/yield,
        negative values for reactants and positive values for products.

    .. note::
        `ins` and `outs` have all the components,
        `upstream` and `downstream` only have their individual components;
        but the mass flow of `ins` equals that of `upstream`,
        and the mass flow of `outs` equals that of `downstream`.
    '''

    # _stacklevel = Unit._stacklevel
    _graphics = BSTjunction._graphics
    # heat_utilities = ()
    # power_utility = PowerUtility()
    # baseline_purchase_cost = 0.
    # baseline_purchase_costs = {}
    # purchase_cost = 0.
    # purchase_costs = {}
    # installed_cost = 0.
    # installed_costs = {}
    # utility_cost = 0.


    def __init__(self, ID='', upstream=None, downstream=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 *, reactions):
        thermo = downstream.thermo
        SanUnit.__init__(self, ID, upstream, downstream, thermo, init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic,
                         skip_property_package_check=True)
        self.reactions = reactions
        self._parse()



    def _parse(self):
        rxns = self.reactions
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

        self.RX = RX
        self.RY = RY


    def _run(self):
        ins = self.ins[0]
        MXsum = ins.mass.reshape(1, len(ins.components))
        MY = (self.RX*MXsum).T @ self.RY # _RX: (num_rxns, num_upcmps); _RY: (num_rxns, num_downcmps)
        MYsum = MY.sum(axis=0) # MY: (num_upcmps, num_downcmps)
        self.outs[0].mass = MYsum


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
        _RX = self.RX
        _RY = self.RY
        def yt(t, QC_ins, dQC_ins):
            for i, j in zip((QC_ins, dQC_ins), (_state, _dstate)):
                MXsum = i[0][:-1] # shape = (1, num_upcmps)
                MY = -(_RX*MXsum).T @ _RY # _RX: (num_rxns, num_upcmps); _RY: (num_rxns, num_downcmps)
                MYsum = MY.sum(axis=0) # MY: (num_upcmps, num_downcmps)
                Q = MY.sum()
                j[:-1] = MYsum
                j[-1] = Q

            # if _n_ins > 1:
            #     Q_ins = QC_ins[:, -1]
            #     C_ins = QC_ins[:, :-1]
            #     dQ_ins = dQC_ins[:, -1]
            #     dC_ins = dQC_ins[:, :-1]
            #     Q = Q_ins.sum()
            #     C = Q_ins @ C_ins / Q
            #     _state[-1] = Q
            #     _state[:-1] = C
            #     Q_dot = dQ_ins.sum()
            #     C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
            #     _dstate[-1] = Q_dot
            #     _dstate[:-1] = C_dot

            # else:
            # _state[:] = QC_ins[0]
            # _dstate[:] = dQC_ins[0]

            _update_state()
            _update_dstate()
        self._AE = yt

        
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE


# %%

import qsdsan as qs


cmps_asm1 = qs.processes.create_asm1_cmps()
thermo_asm1 = qs.get_thermo()
s1 = qs.WasteStream('s1')
for ID in cmps_asm1.IDs: s1.imass[ID] = 1


cmps_adm1 = qs.processes.create_adm1_cmps()
s2 = qs.WasteStream('s2')
for ID in cmps_adm1.IDs: s2.imass[ID] = 2


reactions = []
for n, ID in enumerate(cmps_asm1.IDs):
    if n < 13:
        rxn = {ID: -1,
               cmps_adm1.IDs[2*n]: 0.5,
               cmps_adm1.IDs[2*n+1]: 0.5,
               }
        reactions.append(rxn)
reactions.append({'S_N2': -1})



DI = DynamicInfluent('DI', thermo=thermo_asm1)

J1 = Junction('J1', upstream=DI.outs[0], downstream=s2, reactions=reactions, isdynamic=True)

sys = System('sys', path=(DI, J1,))
sys.set_dynamic_tracker(DI, s2, J1)
sys.simulate(
    state_reset_hook='reset_cache',
    t_span=(0, 10),
    t_eval=np.arange(0, 10.5, 0.5),
    )



# %%

# # Probably shouldn't use reactions at all,
# # should directly use the array

# from thermosteam.reaction import Reaction as Rxn, ParallelReaction as PRxn
# cmps_compiled = qs.Components((*cmps_asm1, *cmps_adm1))
# cmps_compiled.compile()
# qs.set_thermo(cmps_compiled)

# lst = []
# for n, ID in enumerate(cmps_asm.IDs):
#     if n < 13:
#         lst.append(Rxn(
#             {
#                 ID: -1,
#                 cmps_adm1.IDs[2*n]: 0.5,
#                 cmps_adm1.IDs[2*n+1]: 0.5,
#             },
#             reactant=ID,
#             X=1,
#             ))
# asm2adm1 = PRxn(lst)


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