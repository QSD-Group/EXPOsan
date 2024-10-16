#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho
from exposan import htl

__all__ = ('create_components',)


def create_components(set_thermo=True):
    htl_cmps = htl.create_components()
    
    # Components in the feedstock
    Lipids = htl_cmps.Sludge_lipid.copy('Lipids')
    Proteins = htl_cmps.Sludge_protein.copy('Proteins')
    Carbohydrates = htl_cmps.Sludge_carbo.copy('Carbohydrates')
    Ash = htl_cmps.Sludge_ash.copy('Ash')
    saf_cmps = Components([
        Lipids, Proteins, Carbohydrates, Ash,
        ])
    
    # Generic components for HTL products
    Biocrude = htl_cmps.Biocrude
    HTLaqueous = htl_cmps.HTLaqueous
    Hydrochar = htl_cmps.Hydrochar
    saf_cmps.extend([Biocrude, HTLaqueous, Hydrochar])
    
    # Components in the biocrude/bio-oil
    biooil_IDs = [
        'C4H10', 'TWOMBUTAN', 'NPENTAN', 'TWOMPENTA', 'CYCHEX',
        'HEXANE', 'TWOMHEXAN', 'HEPTANE', 'CC6METH', 'PIPERDIN',
        'TOLUENE', 'THREEMHEPTA', 'OCTANE', 'ETHCYC6', 'ETHYLBEN',
        'OXYLENE', 'C9H20', 'PROCYC6', 'C3BENZ', 'FOURMONAN', 'C10H22',
        'C4BENZ', 'C11H24', 'C10H12', 'C12H26', 'C13H28', 'C14H30',
        'OTTFNA', 'C6BENZ', 'OTTFSN', 'C7BENZ', 'C8BENZ', 'C10H16O4',
        'C15H32', 'C16H34', 'C17H36', 'C18H38', 'C19H40', 'C20H42', 'C21H44',
        'TRICOSANE', 'C24H38O4', 'C26H42O4', 'C30H62',
        ]
    saf_cmps.extend([i for i in htl_cmps if i.ID in biooil_IDs])
    
    # Components in the aqueous product
    aq_kwargs = {
        'phase': 'l',
        'particle_size': 'Soluble',
        'degradability': 'Undegradable',
        'organic': False,
        }
    H2O = htl_cmps.H2O
    C = Component('C', search_ID='Carbon', **aq_kwargs)
    N = Component('N', search_ID='Nitrogen', **aq_kwargs)
    P = Component('P', search_ID='Phosphorus', **aq_kwargs)
    K = Component('K', search_ID='Potassium', **aq_kwargs)
    KH2PO4= Component('KH2PO4', **aq_kwargs)
    saf_cmps.extend([H2O, C, N, P, K, KH2PO4,])
    
    # Components in the gas product
    CO2 = htl_cmps.CO2
    CH4 = htl_cmps.CH4
    C2H6 = htl_cmps.C2H6
    O2 = htl_cmps.O2
    N2 = htl_cmps.N2
    CO = htl_cmps.CO
    H2 = htl_cmps.H2
    NH3 = htl_cmps.NH3
    saf_cmps.extend([CO2, CH4, C2H6, O2, N2, CO, H2, NH3])
    
    # Surrogate compounds based on the carbon range
    org_kwargs = {
        'particle_size': 'Soluble',
        'degradability': 'Slowly',
        'organic': True,
        }
    # Tb = 391.35 K (118.2°C)
    Gasoline = Component('Gasoline', search_ID='C8H18', phase='l', **org_kwargs)
    # Tb = 526.65 K (253.5°C)
    SAF = Component('SAF', search_ID='C14H30', phase='l', **org_kwargs)
    # Tb = 632.15 K (359°C)
    Diesel = Component('Diesel', search_ID='C21H44', phase='l', **org_kwargs)
    saf_cmps.extend([Gasoline, SAF, Diesel])

    # Consumables only for cost purposes, thermo data for these components are made up
    sol_kwargs = {
        'phase': 's',
        'particle_size': 'Particulate',
        'degradability': 'Undegradable',
        'organic': False,
        }
    HC_catalyst = Component('HC_catalyst', **sol_kwargs) # Fe-ZSM5
    add_V_from_rho(HC_catalyst, 1500)
    HC_catalyst.copy_models_from(Component('CaCO3'),('Cn',))
    
    HT_catalyst = HC_catalyst.copy('HT_catalyst') # Pd-Al2O3
    
    EOmembrane = HC_catalyst.copy('EOmembrane')
    EOanode = HC_catalyst.copy('EOanode')
    EOcathode = HC_catalyst.copy('EOcathode')
    
    ECmembrane = HC_catalyst.copy('ECmembrane')
    ECanode = HC_catalyst.copy('ECanode')
    ECcathode = HC_catalyst.copy('ECcathode')
    
    saf_cmps.extend([
        HC_catalyst, HT_catalyst,
        EOmembrane, EOanode, EOcathode,
        ECmembrane, ECanode, ECcathode,
        ])
    
    for i in saf_cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)
        i.default() # default properties to those of water

    saf_cmps.compile()
    saf_cmps.set_alias('H2O', 'Water')
    saf_cmps.set_alias('H2O', '7732-18-5')
    
    # So that HTL units can run   
    saf_cmps.set_alias('Lipids', 'Sludge_lipid')
    saf_cmps.set_alias('Proteins', 'Sludge_protein')
    saf_cmps.set_alias('Carbohydrates', 'Sludge_carbo')
    saf_cmps.set_alias('Carbohydrates', 'Carbs')
    saf_cmps.set_alias('Ash', 'Sludge_ash')
    
    saf_cmps.set_alias('C', 'Carbon')
    saf_cmps.set_alias('N', 'Nitrogen')
    saf_cmps.set_alias('P', 'Phosphorus')
    saf_cmps.set_alias('P', 'Potassium')

    if set_thermo: qs_set_thermo(saf_cmps)

    return saf_cmps