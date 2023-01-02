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
from warnings import warn
from math import pi
from biosteam.exceptions import bounds_warning, DesignWarning
from biosteam.units.design_tools import flash_vessel_design
from biosteam.units.design_tools.specification_factors import material_densities_lb_per_ft3
from qsdsan import Construction
from qsdsan.utils import auom

__all__ = ('StorageTank',)

_lb_to_kg = auom('lb').conversion_factor('kg')
_m_to_ft = auom('m').conversion_factor('ft')
_Pa_to_psi = auom('Pa').conversion_factor('psi')

# =============================================================================
# HTL_storage_tank
# =============================================================================
        
class StorageTank(qs.sanunits.StorageTank):
    '''
    Similar to qsdsan.sanunits.StorageTank, but can calculate material usage.
    '''
    
    _units = {'Diameter': 'ft',
              'Length': 'ft',
              'Wall thickness': 'in',
              'Weight': 'kg'}
    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Horizontal vessel diameter': (3, 21),
               'Vertical vessel length': (12, 40)}
    _vessel_material = 'Stainless steel'
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  vessel_type=None, tau=None, V_wf=None,
                  vessel_material=None, kW_per_m3=0.,
                  init_with='WasteStream', F_BM_default=None, length_to_diameter=2):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo,
                      vessel_type=vessel_type, tau=tau, V_wf=V_wf,
                      vessel_material=vessel_material, kW_per_m3=kW_per_m3,
                      init_with=init_with, F_BM_default=F_BM_default)
        self.length_to_diameter = length_to_diameter
        
    def _init_lca(self):
        item_name = self.vessel_material.replace(' ', '_')
        self.construction = [
            Construction(item_name.lower(), linked_unit=self, item=item_name, quantity_unit='kg'),
            ]
        
    
    def _design(self):
        super()._design()
        D = self.design_results
        
        Diameter = (4*D['Total volume']/pi/self.length_to_diameter)**(1/3)
        Diameter *= _m_to_ft # convert from m to ft
        L = Diameter * self.length_to_diameter # ft

        Tank_design = self._horizontal_vessel_design(self.ins[0].P*_Pa_to_psi, Diameter, L)
        
        D['Diameter'] = Diameter
        D['Length'] = L
        D['Wall thickness'] = Tank_design['Wall thickness']
        D['Material'] = self.vessel_material
        self.construction[0].quantity = D['Weight'] = Tank_design['Weight']*_lb_to_kg


    def _horizontal_vessel_design(self, pressure, diameter, length) -> dict:
        pressure = pressure
        diameter = diameter
        length = length
        # Calculate vessel weight and wall thickness
        if self.vessel_material == 'Carbon steel':
            rho_M = material_densities_lb_per_ft3[self.vessel_material]
        else:
            rho_M = material_densities_lb_per_ft3['Stainless steel 304']
        if pressure < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=DesignWarning)
        VW, VWT = flash_vessel_design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M)
        bounds_warning(self, 'Horizontal vessel weight', VW, 'lb',
                       self._bounds['Horizontal vessel weight'], 'cost')
        bounds_warning(self, 'Horizontal vessel diameter', diameter, 'ft',
                       self._bounds['Horizontal vessel diameter'], 'cost')
        Design = {}
        Design['Vessel type'] = 'Horizontal'
        Design['Length'] = length # ft
        Design['Diameter'] = diameter # ft
        Design['Weight'] = VW # lb
        Design['Wall thickness'] = VWT # in
        return Design

    @property
    def vessel_material(self):
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, i):
        exist_material = getattr(self, '_vessel_material', None)
        qs.sanunits.StorageTank.vessel_material.fset(self, i)
        if i and exist_material == i: return # type doesn't change, no need to reload construction items
        self._init_lca()