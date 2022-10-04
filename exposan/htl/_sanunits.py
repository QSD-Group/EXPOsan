#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Jianan Feng <jiananf2@illinois.edu>
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''

from warnings import warn
from qsdsan import SanUnit

__all__ = (
    'HTL',
    'HT',
    'CHG',
    'StruvitePrecipitation',
    'MembraneDistillation')

class HTL(SanUnit):
    '''
    HTL converts dewatered sludge to biocrude, aqueous, off-gas, and biochar under
    elevated temperature and pressure. The products percentage (wt%) can be evaluatad
    using revised MCA model (Li et al., 2017, Leow et al., 2018) with known sludge
    composition (protein%, lipid%, and carbohydrate%)
    
    
    
    
    Parameters
    ----------
    

    
    References
    ----------
    
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)

        
    _N_ins=1 #One wastestream with four components (H2O, Sludge_lipid, Sludge_protein, and Sludge_carbon)
    _N_outs=4
        
        
    def _run(self):
        pass
        # Design a function to check if the sum of sludge compositions is one. If not, warn('').
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
        
        
        
        
    