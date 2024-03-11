#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# import biosteam as bst
from math import ceil, log
from qsdsan import SanUnit
from qsdsan.sanunits import Reactor
from qsdsan.utils import auom

__all__ = (
    'XXX',
    'XXX',
    )

CO2 absorber/stripper

electrochemical cells

# =============================================================================
# XXX
# =============================================================================

class XXX():
    '''
    XXX.
    
    Parameters
    ----------
    ins : Iterable(stream)
        XXX, XXX.
    outs : Iterable(stream)
        XXX, XXX.
        
    References
    ----------

    '''
    _N_ins = 
    _N_outs = 
        
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.XXX = XXX

    def _run(self):
        
        XXX, XXX = self.ins
        XXX, XXX = self.outs

    def _design(self):
        
    def _cost(self):
        
    @property
    def XXX(self):
        return