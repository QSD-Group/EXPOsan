# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 14:18:10 2022

@author: joy_c
"""

from math import pi

__all__ = ('gas_holder_items',
           'iron_sponge_items',
           'gas_holder_input',
           'iron_sponge_input')
#%%

gas_holder_items = {
    'PE': ['textile production, nonwoven polyester, needle-punched', 'kg'],
    'PVC': ['polyvinylchloride production, emulsion polymerisation', 'kg'],
    'Varnish': ['acrylic varnish production, product in 87.5% solution state', 'kg'],
    'Slab concrete': ['concrete slab production', 'm3']
    }

# http://mvseer.com/wp-content/uploads/MV_IRON_SPONGE-MSDS.pdf
iron_sponge_items = {
    'Wood chips': ['wood chips production, softwood, at sawmill', 'kg'],
    'Iron oxide': ['portafer production', 'kg'],
    'Soda ash': ['soda ash, dense, to generic market for neutralising agent', 'kg'],
    'CaCO3': ['calcium carbonate production, precipitated', 'kg'],
    # 'Carbon steel': ['reinforcing steel production', 'kg']
    }

density = {
    'PE': 1300,  # 1230-1380 kg/m3
    'PVC': 1380, # kg/m3
    'varnish': 900, # https://www.industrialcoatingsltd.com/app/uploads/2020/10/Selett-Clear-Varnish-Data-Sheet.pdf
    'iron_sponge': 800, 
    'steel': 7850  # https://intotheboxes.com/unit-weight-for-steel/517/#:~:text=Density%20of%20steel%20rebar%20in,is%207850%20kg%2Fm3.
    }

def cap_radius_from_V(V, h2r=0.7):
    c = (4 - h2r**2 * (3-h2r))/3
    return (V/(c*pi))**(1/3)

def A_cap(r, h2r=0.7):
    return 2*pi*(2 - h2r)*r**2

# V_tank = 10 # m3, biogas production is close to 20 m3/d at baseline

# d_pe = 1e-3     # assume thickness 1 mm
# d_pvc = 75e-6   # assume surface treatment thickness = 75 um
# d_slab = 0.2    # assume 20 cm

# # https://www.industrialcoatingsltd.com/app/uploads/2020/10/Selett-Clear-Varnish-Data-Sheet.pdf
# varnish_spreading_rate = 10.65e3 # 11.6 - 9.7 m2/L 

def gas_holder_input(V_tank, d_pe=1e-3, d_pvc=75e-6, d_slab=0.2, 
                     varnish_spreading_rate=10.65e3):
    R_tank = cap_radius_from_V(V_tank)
    A_membrane = A_cap(R_tank)
    
    m_pe = 2 * A_membrane * d_pe * density['PE']
    m_pvc = 4 * A_membrane * d_pe * density['PVC']
    m_varnish = 4 * A_membrane / varnish_spreading_rate * density['varnish']
    V_slab = pi * R_tank**2 * d_slab
    
    return dict(zip(gas_holder_items.keys(), [m_pe, m_pvc, m_varnish, V_slab]))

# V_sponge = 1.5 # m3

def iron_sponge_input(density=800, Fe2O3_content=15):
    '''mass fractions of media materials, in kg/kg iron sponge'''
    m_iron_oxide = Fe2O3_content * 12.872 / density
    m_wood = 0.25
    m_soda_ash = 0.03 / 1.325     # Convert from Na2CO3 to 2*(NaOH)
    m_caco3 = 0.015

    return dict(zip(iron_sponge_items.keys(), 
                    [m_wood, m_iron_oxide, m_caco3, m_soda_ash]))