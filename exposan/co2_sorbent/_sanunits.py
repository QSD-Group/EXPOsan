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
    # 'XXX',
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
# class ALFCrystallizer(bst.BatchCrystallizer):
#     
#     def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
#                  T, crystal_ALF_purity=0.95, melt_AcTAG_purity=0.90,
#                  order=None):
#         bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
#                                        tau=5, V=1e6, T=T)
#         self.melt_AcTAG_purity = melt_AcTAG_purity
#         self.crystal_TAG_purity = crystal_TAG_purity
# 
#     @property
#     def Hnet(self):
#         feed = self.ins[0]
#         effluent = self.outs[0]
#         if 's' in feed.phases:
#             H_in = - sum([i.Hfus * j for i,j in zip(self.chemicals, feed['s'].mol) if i.Hfus])
#         else:
#             H_in = 0.
#         solids = effluent['s']
#         H_out = - sum([i.Hfus * j for i,j in zip(self.chemicals, solids.mol) if i.Hfus])
#         return H_out - H_in
#         
#     def _run(self):
#         outlet = self.outs[0]
#         outlet.phases = ('s', 'l')
#         crystal_TAG_purity = self.crystal_TAG_purity
#         melt_AcTAG_purity = self.melt_AcTAG_purity
#         feed = self.ins[0]
#         TAG, AcTAG = feed.imass['TAG', 'AcTAG']
#         total = TAG + AcTAG
#         minimum_melt_purity = AcTAG / total
#         minimum_crystal_purity = TAG / total
#         outlet.empty()
#         if crystal_TAG_purity < minimum_crystal_purity:
#             outlet.imol['s'] = feed.mol
#         elif melt_AcTAG_purity < minimum_melt_purity:
#             outlet.imol['l'] = feed.mol
#         else: # Lever rule
#             crystal_AcTAG_purity = (1. - crystal_TAG_purity)
#             melt_fraction = (minimum_melt_purity - crystal_AcTAG_purity) / (melt_AcTAG_purity - crystal_AcTAG_purity)
#             melt = melt_fraction * total
#             AcTAG_melt = melt * melt_AcTAG_purity
#             TAG_melt = melt - AcTAG
#             outlet.imass['l', ('AcTAG', 'TAG')] = [AcTAG_melt, TAG_melt]
#             outlet.imol['s'] = feed.mol - outlet.imol['l']
#         outlet.T = self.T
# =============================================================================

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