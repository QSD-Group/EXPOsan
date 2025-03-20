#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import SanUnit

__all__ = (
    # 'AerobicDigester',
    # 'HeatDrying',
    # 'LandApplication',
    'Landfilling',
    # 'LimeStabilization',
    )

# =============================================================================
# Landfilling
# =============================================================================

class Landfilling(SanUnit):
    '''
    Sludge / biosolids disposal through landfilling, including hauling.
    Costs based on [1], and fugutive emissions from landfilling based on [2].
    
    Parameters
    ----------
    ins : iterable
        XXX, XXX.
    outs : iterable
        XXX, XXX, XXX.
    hauling_distance : float
        The hauling distance from wastewater treatment plants to landfills, [km].
    timescale : float
        XXX.
    '''
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', hauling_distance=13,
                 timescale=1):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.hauling_distance = hauling_distance
        self.timescale = timescale
    
    def _run(self):
        pass
    
    def _design(self):
        pass
    
    def _cost(self):
        pass