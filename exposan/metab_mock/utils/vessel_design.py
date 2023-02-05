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
from math import pi

__all__ = ('heat_transfer_U', 'stainless_steel_wall_thickness', 'UASB_sizing')

#%%

_h_air = 10         # W/m2/K, convective heat transfer coefficient, assume free convection, https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html
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

def UASB_sizing(Q, Vliq, Vgas, max_h2d, upflow_velocity, pipe_velocity):
    '''
    Determine aspect ratio of the reactor based on design upflow velocity.

    Parameters
    ----------
    Q : float
        Total influent flowrate [m3/d].
    Vliq : float
        Liquid phase volume [m3].
    Vgas : float
        Gas phase volume [m3].
    max_h2d : float
        Maximum liquid depth to diameter ratio.
    upflow_velocity : float
        Design upflow velocity [m/h].
    pipe_velocity : float
        Flow velocity through the influent pipe [m/h]

    Returns
    -------
    Vtot : float
        Reactor total volume [m3].
    H : float
        Reactor height [m].
    dia : float
        Reactor diameter [m].
    velocity_head : float
        Additional velocity head for pumping.
    '''
    hrt = Vliq/Q  # d
    vlim_h2d = (4*Q*max_h2d**2/pi/hrt**2)**(1/3) # m/d
    vpipe = pipe_velocity/3600 # m/s
    if vlim_h2d < upflow_velocity*24:
        h = hrt * vlim_h2d
        dia = h/max_h2d
        velocity_head = ((vpipe * upflow_velocity*24/vlim_h2d)**2 - vpipe**2)/2/9.80665 # m
    else:
        h = hrt * upflow_velocity*24
        dia = (4*Q/pi/(upflow_velocity*24))**(1/2)
        velocity_head = 0.
    H = h*(1+Vgas/Vliq)
    Vtot = Vgas + Vliq
    return Vtot, H, dia, velocity_head