#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import SanUnit

__all__ = ('PhaseChanger',)

# =============================================================================
# PhaseChanger
# =============================================================================

class PhaseChanger(SanUnit):
    '''
    Change the effluent phase to the desired one, also allow the switch between stream types.
    
    Parameters
    ----------
    ins : Iterable(stream)
        influent
    outs : Iterable(stream)
        effluent
    phase : str
        Desired phase, can only be one of ("g", "l", or "s").
    '''
    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', phase='l'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.phase = phase

        
    def _run(self):
        
        influent = self.ins[0]
        effluent = self.outs[0]
        effluent.copy_like(influent)
        effluent.phase = self.phase
        
    @property
    def phase(self):
        return self._phase
    @phase.setter
    def phase(self, i):
        if not i in ('g', 'l', 's'):
            raise ValueError('`phase` must be one of ("g", "l", or "s"), '
                             f'not "{i}".')
        self._phase = i