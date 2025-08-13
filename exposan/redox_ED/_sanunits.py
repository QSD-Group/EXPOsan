#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:27:03 2025

@author: blues
"""

from qsdsan import SanUnit,Construction, WasteStream

__all__ = ('redox_ED',)

#%%
class redox_ED(SanUnit):
    _N_ins = 2
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
    def _run(self):
        fc_in, ac_in = self.ins
        fc_out, ac_out = self.outs
        
        
from exposan import g2rt
g2rt._components_loaded
g2rt.load()
g2rt._components_loaded
g2rt.components