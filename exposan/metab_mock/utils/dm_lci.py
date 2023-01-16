# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 12:33:21 2022

@author: joy_c
"""

# from qsdsan.utils import auom
from math import pi

__all__ = ('DuPont_items', 'DuPont_input')
#%%

DuPont_items = {
    'PP': ['polypropylene production', 'kg'],
    'PVC': ['polyvinylchloride production, emulsion polymerisation', 'kg'],
    'PS': ['polysulfone production, for membrane filtration production', 'kg'],
    'Epoxy': ['epoxy resin production, liquid', 'kg'],
    'electric': ['electricity, high voltage, production mix', 'kWh'],
    'molding': ['injection moulding', 'kg'],
    'extrusion': ['extrusion, plastic pipes', 'kg']
    }

density = {
    'PMP': 833,    # kg/m3
    'PVC': 1380,
    'PS': 1240,
    'Epoxy': 1100, 
    }

# https://www.dupont.com/content/dam/dupont/amer/us/en/water-solutions/public/documents/en/MDG-Ligasep-LDM-040-PDS-45-D00501-en.pdf
# od_fiber = 210e-6           # 180-240 um
# dw_fiber = 35e-6            # 30-40 um
# l_fiber = l_shell = 536e-3  # 536 mm
# od_shell = 165e-3           # 165 mm
# dw_shell = 2.5e-3           # assume 2.5 mm thick housing plastic
# V_liq = 6.5e-3              # 6.5 L
# od_potting = 180e-3         # 180 mm
# l_potting = 36e-3 * 2       # assume 36 mm of potting length on each end
# od_pipe = 108e-3            # 108 mm
# h_pipe = (159-165/2)*1e-3   # mm
# m_unit = 10                 # kg

def V_cylinder(od, h):
    return pi*(od/2)**2*h

def V_cylindrical_shell(od, dw, h):
    return pi * ((od/2)**2 - (od/2-dw)**2) * h

def S_cylinder_wall(od, dw, h):
    return pi * (od-dw/2) * h

def specific_area(od, dw):
    '''membrane area per unit volume of cylindrical fiber [m2/m3].'''
    return (od-dw/2) / (od/2)**2

def DuPont_input(od_fiber=2.1e-4, dw_fiber=3.5e-5, l_fiber=0.536,
                 od_shell=0.165, dw_shell=2.5e-3, l_shell=0.536,
                 od_potting=0.180, l_potting=0.072, od_pipe=0.108, 
                 h_pipe=0.0765, V_liq=6.5e-3, total_mass=10):
    '''
    Estimate the technosphere input to manufacture a 
    DuPont Ligasep LDM-040 membrane module.

    Parameters
    ----------
    od_fiber : float, optional
        Outer diameter of a single fiber [m]. The default is 2.1e-4.
    dw_fiber : float, optional
        Membrane thickness [m]. The default is 3.5e-5.
    l_fiber : float, optional
        Fiber length [m]. The default is 0.536.
    od_shell : float, optional
        Outer diameter of the housing shell [m]. The default is 0.165.
    dw_shell : float, optional
        Shell thickness [m]. The default is 2.5e-3.
    l_shell : float, optional
        Shell length [m]. The default is 0.536.
    od_potting : float, optional
        Outer diameter of the potting section [m]. The default is 0.180.
    l_potting : float, optional
        Potting depth [m]. The default is 0.072.
    od_pipe : float, optional
        Approximate pipe outer diameter [m]. The default is 0.108.
    h_pipe : float, optional
        Approximate pipe length [m]. The default is 0.0765.
    V_liq : float, optional
        Liquid phase volume [m3]. The default is 6.5e-3.
    total_mass : float, optional
        Total mass of the module [kg]. The default is 10.
    '''
    V_fibers = V_cylinder(od_shell-dw_shell*2, l_shell) - V_liq
    S_membrane = V_fibers * specific_area(od_fiber, dw_fiber)
    electricity = 0.0408 * S_membrane   # kWh
    
    V_pvc = V_cylindrical_shell(od_shell, dw_shell, l_shell)
    m_pvc = V_pvc * density['PVC']
    V_cap = pi * (od_shell/2)**2 * dw_shell * 2
    V_pipe = pi * (od_pipe/2)**2 * h_pipe * 4
    m_ps = (V_cap + V_pipe) * density['PS']
    V_epoxy = V_cylinder(od_potting-dw_shell*2, l_potting) - V_fibers * l_potting/l_fiber
    m_epoxy = V_epoxy * density['Epoxy']
    m_pp = total_mass - (m_pvc + m_ps + m_ps)
    m_mold = m_pvc + V_cap * density['PS']
    m_extru = m_pp + V_pipe * density['PS']
    
    return dict(zip(DuPont_items.keys(), [m_pp, m_pvc, m_ps, m_epoxy, electricity, m_mold, m_extru]))