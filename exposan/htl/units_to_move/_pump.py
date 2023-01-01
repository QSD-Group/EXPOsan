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
from math import pi
from qsdsan import Construction
from qsdsan.utils import auom, select_pipe

__all__ = ('Pump',)

_ft3_to_gal = auom('ft3').conversion_factor('gallon')
_lb_to_kg = auom('lb').conversion_factor('kg')
_m3_to_gal = auom('m3').conversion_factor('gallon')


# =============================================================================
# HTLpump            
# =============================================================================

class Pump(qs.sanunits.Pump):
    '''
    Pumps used in HTL system
    See qsdsan.sanunits.WWTpump for pipe and pump weight calculation
    See bst.units.Pump for other functions
    All pumps are assumed to be made of stainless steel and specific for sludge.
    
    Parameters
    ----------
    P : float
        pump pressure
        
    References
    ----------
    .. [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
        Valorization of Dilute Organic Carbon Waste Streams.
        Energy Environ. Sci. 2016, 9 (3), 1102–1112.
        https://doi.org/10.1039/C5EE03715H.
        
    See Also
    --------
    :class:`qsdsan.sanunits.WWTpump`
    
    :class:`biosteam.units.Pump`
    '''
    
    _N_pump = 1
    _H_ts = 0. # total static head
    _H_p = 0. # total pressure head
    _v = 3 # fluid velocity, [ft/s]
    _C = 110 # Hazen-Williams coefficient for stainless steel (SS)
    _SS_per_pump = 725 * 0.5
    _units = {'Pump pipe stainless steel': 'kg',
              'Pump stainless steel': 'kg'}

    def _design(self):
        super()._design()
        
        pipe, pumps, hdpe = self.design_sludge()

        D = self.design_results
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pumps
        construction = getattr(self, 'construction', ()) # would work for both biosteam/qsdsan units
        if construction: construction[0].quantity = pipe + pumps
        else:
            self.construction = [
                Construction('stainless_steel', linked_unit=self, item='Stainless_steel', 
                             quantity=pipe + pumps, quantity_unit='kg'),
                ]

    def design_sludge(self, Q_mgd=None, N_pump=None, **kwargs):
        '''
        Design pump for handling waste sludge.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_pump : int
            Number of the pumps.
        kwargs : dict
            Additional attribute values to set (e.g., `L_s`, `H_ts`),
            this will overwrite the default values.
        '''
        Q_mgd = Q_mgd or self.Q_mgd
        N_pump = N_pump or 1

        val_dct = dict(
            L_s=50, # length of suction pipe, [ft]
            L_d=50, # length of discharge pipe, [ft]
            H_ts=0., # H_ds_LIFT (D) - H_ss_LIFT (0)
            H_p=0. # no pressure
            )
        val_dct.update(kwargs)

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd, N_pump=N_pump, **val_dct)

        return M_SS_IR_pipe, M_SS_IR_pump, 0
    
    def _design_generic(self, Q_mgd, N_pump=None, L_s=0., L_d=0., H_ts=0., H_p=0.):
        self.Q_mgd = Q_mgd
        self._H_ts = H_ts or self.H_ts
        self._H_p = H_p or self.H_p
        N_pump = N_pump or self.N_pump

        v, C, Q_cfs = self.v, self.C, self.Q_cfs # [ft/s], -, [ft3/s]

        ### Suction side ###
        # Suction pipe (permeate header) dimensions
        OD_s, t_s, ID_s = select_pipe(Q_cfs/N_pump, v) # [in]

        # Suction friction head, [ft]
        self._H_sf = 3.02 * L_s * (v**1.85) * (C**(-1.85)) * ((ID_s/12)**(-1.17))

        ### Discharge side ###
        # Discharge pipe (permeate collector) dimensions
        OD_d, t_d, ID_d = select_pipe(Q_cfs, v)

        # Discharge friction head, [ft]
        self._H_df = 3.02 * L_d * (v**1.85) * (C**(-1.85)) * ((ID_d/12)**(-1.17))

        ### Material usage ###
        # Pipe SS, assume stainless steel, density = 0.29 lbs/in3
        # SS volume for suction, [in3]
        self._N_pump = N_pump
        V_s = N_pump * pi/4*((OD_s)**2-(ID_s)**2) * (L_s*12)
        # SS volume for discharge, [in3]
        V_d = pi/4*((OD_d)**2-(ID_d)**2) * (L_d*12)

        # Total SS mass, [kg]
        M_SS_pipe = 0.29 * (V_s+V_d) * _lb_to_kg
        M_SS_pump = N_pump * self.SS_per_pump
        return M_SS_pipe, M_SS_pump
    
    @property
    def Q_mgd(self):
        '''
        [float] Volumetric flow rate in million gallon per day, [mgd].
        Will use total volumetric flow through the unit if not provided.
        '''
        return self.F_vol_in*_m3_to_gal*24/1e6
    @Q_mgd.setter
    def Q_mgd(self, i):
        self._Q_mgd = i
        
    @property
    def Q_cfs(self):
        '''[float] Volumetric flow rate in cubic feet per second, [cfs].'''
        return self.Q_mgd*1e6/24/60/60/_ft3_to_gal
        
    @property
    def H_ts(self):
        '''[float] Total static head, [ft].'''
        return self._H_ts
    
    @property
    def H_p(self):
        '''[float] Pressure head, [ft].'''
        return self._H_p
    
    @property
    def v(self):
        '''[float] Fluid velocity, [ft/s].'''
        return self._v
    @v.setter
    def v(self, i):
        self._v = i
        
    @property
    def C(self):
        '''[float] Hazen-Williams coefficient to calculate fluid friction.'''
        return self._C
    @C.setter
    def C(self, i):
        self._C = i
        
    @property
    def SS_per_pump(self):
        '''[float] Quantity of stainless steel per pump, [kg/ea].'''
        return self._SS_per_pump
    @SS_per_pump.setter
    def SS_per_pump(self, i):
        self._SS_per_pump = i