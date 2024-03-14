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
import biosteam as bst
from math import ceil, log
from qsdsan import SanUnit
from qsdsan.sanunits import Reactor
from qsdsan.utils import auom

__all__ = (
    'ALFProduction',
    'ALFCrystallizer',
    )

# =============================================================================
# CO2 absorber/stripper
# 
# electrochemical cells
# =============================================================================

# =============================================================================
# ALFProduction
# =============================================================================

class ALFProduction(bst.CSTR):
    
    _N_ins = 1
    _N_outs = 1
    # TODO: confirm T and tau
    T_default = 60 + 273.15
    P_default = 101325
    tau_default = 8 # hr
    
    def _setup(self):
        super()._setup()
        self.ALF_production = bst.Reaction('AlH3O3,s + 3HCOOH,l -> C3H3AlO6,s + 3H2O,l', 'AlH3O3', 1)
    
    def _run(self):
        effluent = self.outs[0]
        effluent.copy_like(self.ins[0])
        self.ALF_production(effluent)
        effluent.T = self.T
        effluent.P = self.P

# =============================================================================
# ALFCrystallizer
# =============================================================================
class ALFCrystallizer(bst.BatchCrystallizer):
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 T, crystal_ALF_yield=1, order=None):
        bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                                       tau=5, V=1e6, T=T)
        self.crystal_ALF_yield = crystal_ALF_yield

    @property
    def Hnet(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        if 's' in feed.phases:
            H_in = - sum([i.Hfus * j for i,j in zip(self.chemicals, feed['s'].mol) if i.Hfus])
        else:
            H_in = 0.
        solids = effluent['s']
        H_out = - sum([i.Hfus * j for i,j in zip(self.chemicals, solids.mol) if i.Hfus])
        return H_out - H_in
        
    def _run(self):
        outlet = self.outs[0]
        outlet.phases = ('s', 'l')
        crystal_ALF_yield = self.crystal_ALF_yield
        feed = self.ins[0]
        ALF = feed.imass['C3H3AlO6']
        outlet.empty()
        
        outlet.imass['s', 'C3H3AlO6'] = ALF*crystal_ALF_yield
        outlet.imass['l', ('C3H3AlO6','HCOOH','H2O')] = [ALF*(1-crystal_ALF_yield), feed.imass['HCOOH'], feed.imass['H2O']]
        
        outlet.T = self.T

# =============================================================================
# class XXX():
#     '''
#     XXX.
#     
#     Parameters
#     ----------
#     ins : Iterable(stream)
#         XXX, XXX.
#     outs : Iterable(stream)
#         XXX, XXX.
#         
#     References
#     ----------
# 
#     '''
#     _N_ins = 
#     _N_outs = 
#         
#     def __init__(self, ID='', ins=None, outs=(), thermo=None,
#                  init_with='WasteStream'):
#         
#         SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
#         self.XXX = XXX
# 
#     def _run(self):
#         
#         XXX, XXX = self.ins
#         XXX, XXX = self.outs
# 
#     def _design(self):
#         
#     def _cost(self):
#         
#     @property
#     def XXX(self):
#         return
# =============================================================================