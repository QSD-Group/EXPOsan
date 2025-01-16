#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import biosteam as bst
from warnings import warn
from math import ceil
import numpy as np
from qsdsan.sanunits._abstract import Mixer
from qsdsan.sanunits import IsothermalCompressor
from qsdsan.processes import Decay
from qsdsan import SanUnit,Construction, WasteStream
from qsdsan.sanunits import SludgeThickening, Copier
from biosteam.units.design_tools import flash_vessel_design
from qsdsan.utils import ospath, data_path, load_data, price_ratio
import CoolProp.CoolProp as CP
from exposan.g2rt._sanunits import VolumeReductionCombustor
g2rt_su_data_path = ospath.join(data_path, 'sanunit_data/g2rt')


__all__ = ('BiomassCombustion',
           )

#%%
@price_ratio()
class BiomassCombustion(VolumeReductionCombustor):
    '''
    Combust dry biomass (feces and wood pellets) to generate heat (in the form of power) 
    to offset the electrical heating demand for single unit reinvented toilet. 
    
    Parameters
    ----------
    ins : Iterable(WasteStream)
        Dewatered feces solid cakes, wood pellets, air
    outs : Iterable(WasteStream)
        Wood Ash, hot gas, fugitive N2O, fugitive CH4, fugitive NO, fugitive SO2
    energy_recovery_efficiency: float from 0 to 1
        The efficiency of the bioler to reuse the combustion heat to offset system
        heating demand
    biofuel_to_solids: float
        g-wood pellets/g-feces solids
    '''
    _N_ins = 3
    _N_outs = 7
    
    

