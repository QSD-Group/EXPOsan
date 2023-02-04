# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan.utils import auom
from biosteam.units.design_tools.flash_vessel_design import compute_vessel_weight_and_wall_thickness

__all__ = ('heat_transfer_U', 'stainless_steel_wall_thickness')

#%%

_h_air = 37.5       # W/m2/K, convective heat transfer coefficient, assume at free air relative velocity = 10 m/s, https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html
_h_water = 3000     # W/m2/K, assume moderate forced flow, https://www.engineersedge.com/heat_transfer/convective_heat_transfer_coefficients__13378.htm
_k_rockwool = 0.038 # 0.035 – 0.040 W/m/K, thermal conductivity
_k_concrete = 2.4   # 1.6-3.2 W/m/K  # http://www.jett.dormaj.com/docs/Volume8/Issue%203/Investigation%20of%20Thermal%20Properties%20of%20Normal%20Weight%20Concrete%20for%20Different%20Strength%20Classes.pdf
_k_ssteel = 20      # 16-24 W/m/K
_k_alum = 230       # 205-250
_k_csteel = 45

_rho_ssteel = 7930  # kg/m3

def heat_transfer_U(twall, tinsl, tface, tbase, concrete=True):
    '''
    Vessel body heat transfer coefficients.

    Parameters
    ----------
    twall : float
        Vessel wall thickness [m].
    tinsl : float
        Vessel insulation layer thickness [m].
    tface : float
        Thickness of facing over the insulation laye [m].
    tbase : float
        Base thickness [m].
    concrete : bool, optional
        Made of concrete of stainless steel. The default is concrete.

    '''
    if concrete: kwall = kbase = _k_concrete
    else: kwall = kbase = _k_ssteel
    Uwall = 1/(1/_h_water + 1/_h_air \
                + twall/kwall \
                + tinsl/_k_rockwool \
                + tface/_k_csteel)
    Ucover = 1/(2/_h_air + twall/kwall)
    if concrete: 
        Ubase = 1/(1/_h_water + tbase/kbase)
    else:
        Ubase = 1/(1/_h_water \
                   + tbase/kbase \
                   + tinsl/_k_rockwool \
                   + tface/_k_csteel)
    return Uwall, Ucover, Ubase

_in2m = auom('inch').conversion_factor('m')
_m2ft = auom('m').conversion_factor('ft')
_Pa2psi = auom('Pa').conversion_factor('psi')

def stainless_steel_wall_thickness(P, D, h):
    '''
    Vessel wall thickness [m] if made of stainless steel.    

    Parameters
    ----------
    P : float
        Vessel absolute internal pressure [Pa].
    D : float
        Vessel inner diameter [m].
    h : float
        Liquid depth [m].
    '''
    P *= _Pa2psi
    D *= _m2ft
    h *= _m2ft
    return compute_vessel_weight_and_wall_thickness(P, D, h, _rho_ssteel)[1] * _in2m
    