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
from warnings import warn
from math import ceil, pi
from biosteam.exceptions import bounds_warning, DesignWarning
from biosteam.units.design_tools import flash_vessel_design
from biosteam.units.design_tools.specification_factors import material_densities_lb_per_ft3
from qsdsan import Construction
from qsdsan.utils import auom

__all__ = ('HXutility',)

_in_to_ft = auom('in').conversion_factor('ft')
_lb_to_kg = auom('lb').conversion_factor('kg')
_Pa_to_psi = auom('Pa').conversion_factor('psi')

# =============================================================================
# HTLHX
# =============================================================================

class HXutility(qs.sanunits.HXutility):
    '''
    Similar to qsdsan.sanunits.HXutility, but can calculate material usage.
    
    References
    ----------
    .. [1] Seider, W. D., Lewin, D. R., Seader, J. D., Widagdo, S., Gani, R., &
           Ng, M. K. (2017). Product and Process Design Principles. Wiley.
           Chapter 12: Heat Exchanger Design.
    '''
    
    line = qs.unit.HXutility.line
    _graphics = qs.unit.HXutility._graphics
    _units = {'Area': 'ft^2',
              'Total tube length': 'ft',
              'Inner pipe weight': 'kg',
              'Outer pipe weight': 'kg',
              'Total steel weight': 'kg',
              'Shell length': 'ft',
              'Shell diameter': 'ft',
              'Shell steel weight': 'kg',
              'Tube weight': 'kg'}
    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Horizontal vessel diameter': (3, 21),
               'Vertical vessel length': (12, 40)}
    
    def _design(self, duty=None):
        # Set duty and run heat utility
        if duty is None: duty = self.Hnet # Includes heat of formation
        inlet = self.ins[0]
        outlet = self.outs[0] 
        T_in = inlet.T
        T_out = outlet.T
        iscooling = duty < 0.
        if iscooling: # Assume there is a pressure drop before the heat exchanger
            if T_out > T_in: T_in = T_out
        else:
            if T_out < T_in: T_out = T_in
        self.add_heat_utility(duty, T_in, T_out, 
                              heat_transfer_efficiency=self.heat_transfer_efficiency,
                              hxn_ok=True)
        bst.units.HX._design(self)

        D = self.design_results
        if D['Area'] < 150: # double pipe
            # Assume use 1 1/4 nominal size of inner tube, based on [1] page 365
            # Table 12.3, when use Schedule 40, surface area per foot is 0.435 ft2
            # and weight is 2.28 lb steel per foot
            D['Total tube length'] = D['Area']/0.435
            D['Inner pipe weight'] = D['Total tube length']*2.28*_lb_to_kg
            # Assume use 2 nominal size of outer tube, same length as inner tube
            # the weight is 3.66 lb steel per foot
            D['Outer pipe weight'] = D['Total tube length']*3.66*_lb_to_kg
            total_steel = D['Total steel weight'] = D['Inner pipe weight'] + D['Outer pipe weight']
        else: # shell and tube
            # assume all tubes are 16 ft long and 3/4 in O.D., 16 BWG
            D['Shell length'] = 16
            # D['Total tube length'] = 'N/A'
            single_shell_area = D['Area']/self.N_shells
            if single_shell_area <= 100:
                D['Shell diameter'] = 1
            elif single_shell_area <= 400:
                D['Shell diameter'] = 2
            elif single_shell_area < 1100:
                D['Shell diameter'] = 3
            else:
                D['Shell diameter'] = 3
                self.N_shells = ceil(D['Area']/1100) # ceil shell number
                # if increase the number of N_shells, ft (correction factor) will increase (max = 1),
                # then the required area will decrease, so the calculation here is conservative.
                single_shell_area = D['Area']/self.N_shells

            Shell_design = self._horizontal_vessel_design(self.ins[0].P*_Pa_to_psi, D['Shell diameter'], D['Shell length'])
            D['Shell steel weight'] = Shell_design['Weight']*self.N_shells*_lb_to_kg
            
            single_tube_area = pi*(3/4)*_in_to_ft*D['Shell length']
            
            D['Total tube length'] = D['Shell length']*single_shell_area/single_tube_area*self.N_shells
            
            # according to [1] page 367, Table 12.4, 3/4 in O.D., 16 BWG tube: 0.520 lb/ft
            D['Tube weight'] = D['Total tube length']*0.520*_lb_to_kg
            
            total_steel = D['Total steel weight'] = D['Shell steel weight'] + D['Tube weight']
            
        construction = getattr(self, 'construction', ())
        if construction: construction[0].quantity = total_steel
        else:
            self.construction = (
                Construction('stainless_steel', linked_unit=self, item='Stainless_steel', 
                             quantity=total_steel, quantity_unit='kg'),
                )
            
    def _horizontal_vessel_design(self, pressure, diameter, length) -> dict:
        pressure = pressure
        diameter = diameter
        length = length
        # Calculate vessel weight and wall thickness
        rho_M = material_densities_lb_per_ft3['Carbon steel']
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