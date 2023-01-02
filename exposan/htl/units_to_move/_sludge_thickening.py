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

import biosteam as bst, qsdsan as qs
from math import floor
from qsdsan import Construction
from qsdsan.utils import auom
from . import Pump

__all__ = ('SludgeCentrifuge',)

_lb_to_kg = auom('lb').conversion_factor('kg')
_m3_to_gal = auom('m3').conversion_factor('gallon')

# =============================================================================
# HTL_sludge_centrifuge
# =============================================================================
    
class SludgeCentrifuge(qs.sanunits.SludgeThickening, bst.units.SolidsCentrifuge):
    '''
    Similar to qsdsan.sanunits.SludgeCentrifuge, but can calculate material usage.
    
    References
    ----------
    .. [1] https://dolphincentrifuge.com/wastewater-centrifuge/ (accessed 12-4-2022).
    '''
    
    _units = {'Total pump stainless steel': 'kg',
              'Total pipe stainless steel': 'kg',
              'Centrifige stainless steel': 'kg',
              'Total stainless steel': 'kg'}

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 sludge_moisture=0.8, solids=(),
                 centrifuge_type='scroll_solid_bowl'):
        qs.sanunits.SludgeThickening.__init__(self, ID, ins, outs, thermo, init_with,
                                sludge_moisture=sludge_moisture,
                                solids=solids)
        self.centrifuge_type = centrifuge_type
        ID = self.ID
        eff = self.outs[0].proxy(f'{ID}_eff')
        sludge = self.outs[1].proxy(f'{ID}_sludge')
        self.effluent_pump = Pump(f'.{ID}_eff_pump', ins=eff, init_with=init_with)
        self.sludge_pump = Pump(f'.{ID}_sludge_pump', ins=sludge, init_with=init_with)

    _run = qs.sanunits.SludgeThickening._run

    def _design(self):
        bst.units.SolidsCentrifuge._design(self)
        D = self.design_results
        self.effluent_pump.simulate()
        self.sludge_pump.simulate()
        D['Total pump stainless steel'] = self.effluent_pump.design_results['Pump stainless steel'] +\
                                          self.sludge_pump.design_results['Pump stainless steel']
        D['Total pipe stainless steel'] = self.effluent_pump.design_results['Pump pipe stainless steel'] +\
                                          self.sludge_pump.design_results['Pump pipe stainless steel']
        # based on [1]:
        # when rated capacity <= 80 GPM: weight = 2500 lb
        # when rated capacity (80, 170]: weight = 4000 lb
        # when rated capacity > 170 GPM, use a combination of large and small centrifuges
        
        D['Number of large centrifuge'] = floor(self.F_vol_in*_m3_to_gal/60/170)
        D['Number of small centrifuge'] = 0
        if self.F_vol_in*_m3_to_gal/60 - D['Number of large centrifuge']*170 <= 80:
            D['Number of small centrifuge'] = 1
        else:
            D['Number of large centrifuge'] += 1
        
        D['Centrifige stainless steel'] = (4000*D['Number of large centrifuge'] + 2500*D['Number of small centrifuge'])*_lb_to_kg
        total_steel = D['Total stainless steel'] = D['Total pump stainless steel'] + D['Total pipe stainless steel'] + D['Centrifige stainless steel']
        
        construction = getattr(self, 'construction', ())
        if construction: construction[0].quantity = total_steel
        else:
            self.construction = (
                Construction('stainless_steel', linked_unit=self, item='Stainless_steel', 
                             quantity=total_steel, quantity_unit='kg'),
                )
        
    def _cost(self):
        qs.sanunits.SludgeThickening._cost(self)