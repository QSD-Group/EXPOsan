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

import qsdsan as qs
from qsdsan import Construction

__all__ = ('CombinedHeatPower',)

# =============================================================================
# HTLCHP
# =============================================================================
    
class CombinedHeatPower(qs.sanunits.CHP):
    '''
    Similar to qsdsan.sanunits.CHP, but can calculate material usage.
    References
    ----------
    .. [1] Havukainen, J.; Nguyen, M. T.; Väisänen, S.; Horttanainen, M.
           Life Cycle Assessment of Small-Scale Combined Heat and Power Plant:
           Environmental Impacts of Different Forest Biofuels and Replacing
           District Heat Produced from Natural Gas. Journal of Cleaner
           Production 2018, 172, 837–846.
           https://doi.org/10.1016/j.jclepro.2017.10.241.
    '''
    
    _units = {'Steel': 'kg',
              'Furnace': 'kg',
              'Concrete': 'kg',
              'Reinforcing steel': 'kg'}
    
    def _init_lca(self):
        self.construction = [
            Construction('carbon_steel', linked_unit=self, item='Carbon_steel', quantity_unit='kg'),
            Construction('furnace', linked_unit=self, item='Furnace', quantity_unit='kg'),
            Construction('concrete', linked_unit=self, item='Concrete', quantity_unit='kg'),
            Construction('reinforcing_steel', linked_unit=self, item='Reinforcing_steel', quantity_unit='kg'),
            ]
    
    def _design(self):
        super()._design()
        D = self.design_results
        constr = self.construction
        
        # material calculation based on [1], linearly scaled on power (kW)
        # in [1], a 580 kW CHP:
        # steel: 20098 kg
        # furnace: 12490 kg
        # reinforced concrete: 15000 kg (concrete + reinforcing steel)
        # 1 m3 reinforced concrete: 98 v/v% concrete with a density of 2500 kg/m3 (2450 kg)
        #                            2 v/v% reinforcing steel with a density of 7850 kg/m3 (157kg)
        factor = self.H_net_feeds/3600/580
        constr[0].quantity = D['Steel'] = factor*20098
        constr[1].quantity = D['Furnace'] = factor*12490
        constr[2].quantity = D['Concrete'] = factor*15000*2450/(2450 + 157)
        constr[3].quantity = D['Reinforcing steel'] = factor*15000*157/(2450 + 157)