#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Bright Elijah <be05055@georgiasouthern.edu & brightcarlelijah@gmail.com>
    
Department of Civil Engineering and Construction
Georgia Southern University

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# %%

from .. import SanUnit
from ..utils import ospath, load_data, data_path

__all__ = ('RawWater',)

raw_water_path = ospath.join(data_path, 'sanunit_data/_raw_water1.tsv')
#raw_water_path = ospath.join(data_path, 'sanunit_data/_raw_water2.tsv')


# %%

class RawWater(SanUnit):
    '''
    References
    ----------
Iii, R.; Stetson, L. Socially Embedded and Sustained Point-of-Use Disinfection : Enhancing Silver Nanoparticle Enabled Ceramic Water Filters with a Navajo Pottery Technique. Thesis, 2020. https://doi.org/10.26153/tsw/13860.   

Stumm, W.; Morgan, J. J. Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters; Wiley, 1996.

Maciel, P. M. F.; Fava, N. de M. N.; Lamon, A. W.; Fernandez-Ibañez, P.; Byrne, J. A.; Sabogal-Paz, L. P. Household Water Purification System Comprising Cartridge Filtration, UVC Disinfection and Chlorination to Treat Turbid Raw Water. Journal of Water Process Engineering 2021, 43, 102203. https://doi.org/10.1016/j.jwpe.2021.102203.
    
Wilhelm, N.; Kaufmann, A.; Blanton, E.; Lantagne, D. Sodium Hypochlorite Dosage for Household and Emergency Water Treatment: Updated Recommendations. Journal of Water and Health 2017, 16 (1), 112–125. https://doi.org/10.2166/wh.2017.012.
    Parameters
    ----------
       ins :none
    outs : Raw water

   

    
    '''

    _N_ins = 0
    _N_outs = 1

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 household_size=6, number_of_households=1, water_demand=3.7, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.household_size = household_size
        self.number_of_households = number_of_households
        self.water_demand = water_demand


        data = load_data(path=raw_water_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)                                                                                  

    def _run(self):
        water = self.outs[0]
        water.empty() 
    
        factor = self.household_size * self.water_demand / 24 # from L per day of water to kg per hour
        water.imass['Ecoli'] = self.E_coli/1000*factor 
        ###! units are MPN/hr need to confirm this is done correctly

        water._TOC = self.TOC # mg/L
        water.imass['Ca'] = self.Ca/1000*factor # kg/hr
        water.imass['Mg'] = self.Ca/1000*factor # kg/hr
        water._turbidity = self.Turbidity # NTU
        water.F_vol = factor
        # water._UVT = self.UVT # %
        
       
        
      