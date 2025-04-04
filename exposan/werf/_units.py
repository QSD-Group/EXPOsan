# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 14:43:51 2025

@author: joy_c
"""

import qsdsan as qs
from qsdsan.sanunits import Splitter

__all__ = ('SelectiveRecovery',)

class SelectiveRecovery(Splitter):
    '''
    Selective removal and recovery of certain wastewater constituents (e.g., Ammonia). 
    Specify recovery efficiency with `split` between two outlet streams: 
    [0] recovered constituents and [1] wastewater.
    
    '''
    def _init_state(self):
        self._state = self._ins_QC[0]
        self._dstate = self._state * 0.

    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Splitter'''
        arr = self._state
        Cs = arr[:-1]
        Q = arr[-1]
        sp = self.split
        recv, eff = self.outs
        if recv.state is None: recv.state = 0 * arr
        if eff.state is None: eff.state = 1 * arr
        recv.state[-1] = 1 # assume recovered in a liquid solution of 1 m3/d
        recv.state[:-1] = sp * Cs * Q  # mass flows in g/d, equivalent to g/m3 concentrations in 1m3/d solution
        eff.state[:-1] = (1-sp) * Cs
        eff.state[-1] = Q

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Splitter'''
        pass
