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
# from exposan.utils import add_V_from_rho
from exposan.saf import create_components as create_saf_components

__all__ = ('create_components',)


# def estimate_heating_values(component):
#     '''
#     Estimate the HHV of a component based on the Dulong's equation (MJ/kg):
        
#         HHV [kJ/g] = 33.87*C + 122.3*(H-O/8) + 9.4*S
        
#     where C, H, O, and S are the wt% of these elements.
        
#     Estimate the LHV based on the HHV as:
        
#         LHV [kJ/g] = HHV [kJ/g] – 2.51*(W + 9H)/100
        
#     where W and H are the wt% of moisture and H in the fuel
    
#     References
#     ----------
#     [1] https://en.wikipedia.org/wiki/Heat_of_combustion
#     [2] https://www.sciencedirect.com/science/article/abs/pii/B9780128203606000072
        
#     '''
#     atoms = component.atoms
#     MW = component.MW
#     HHV = (33.87*atoms.get('C', 0)*12 +
#            122.3*(atoms.get('H', 0)-atoms.get('O', 0)/8) +
#            9.4*atoms.get('S', 0)*32
#            )/MW
#     LHV = HHV - 2.51*(9*atoms.get('H', 0)/MW)
    
#     return HHV*MW*1000, LHV*MW*1000

def create_components(set_thermo=True):
    saf_cmps = create_saf_components(set_thermo=False)
    biobinder_cmps = Components([i for i in saf_cmps])
    
    # Other needed components
    Biofuel = saf_cmps.C16H34.copy('Biofuel') # Tb = 559 K
    Biobinder = saf_cmps.TRICOSANE.copy('Biobinder') # Tb = 654 K
    
    biobinder_cmps.extend([Biofuel, Biobinder])
    
    # htl_cmps = htl.create_components()
    
    # # Components in the feedstock
    # Lipids = htl_cmps.Sludge_lipid.copy('Lipids')
    # Proteins = htl_cmps.Sludge_protein.copy('Proteins')
    # Carbohydrates = htl_cmps.Sludge_carbo.copy('Carbohydrates')
    # Ash = htl_cmps.Sludge_ash.copy('Ash')
    
    # # Generic components for HTL products
    # Biocrude = htl_cmps.Biocrude
    # HTLaqueous = htl_cmps.HTLaqueous
    # Hydrochar = htl_cmps.Hydrochar
    
    # # Components in the biocrude
    # org_kwargs = {
    #     'particle_size': 'Soluble',
    #     'degradability': 'Slowly',
    #     'organic': True,
    #     }
    # biocrude_dct = { # ID, search_ID (CAS#)
    #     '1E2PYDIN':     '2687-91-4',
    #     # 'C5H9NS':       '10441-57-3',
    #     'ETHYLBEN':     '100-41-4',
    #     '4M-PHYNO':     '106-44-5',
    #     '4EPHYNOL':     '123-07-9',
    #     'INDOLE':       '120-72-9',
    #     '7MINDOLE':     '933-67-5',
    #     'C14AMIDE':     '638-58-4',
    #     'C16AMIDE':     '629-54-9',
    #     'C18AMIDE':     '124-26-5',
    #     'C16:1FA':      '373-49-9',
    #     'C16:0FA':      '57-10-3',
    #     'C18FACID':     '112-80-1',
    #     'NAPHATH':      '91-20-3',
    #     'CHOLESOL':     '57-88-5',
    #     'AROAMINE':     '74-31-7',
    #     'C30DICAD':     '3648-20-2',
    #     }
    # biocrude_cmps = {}
    # for ID, search_ID in biocrude_dct.items():
    #     cmp = Component(ID, search_ID=search_ID, **org_kwargs)
    #     if not cmp.HHV or not cmp.LHV: 
    #         HHV, LHV = estimate_heating_values(cmp)
    #         cmp.HHV = cmp.HHV or HHV
    #         cmp.LHV = cmp.LHV or LHV
    #     biocrude_cmps[ID] = cmp
        
    # # # Add missing properties
    # # # http://www.chemspider.com/Chemical-Structure.500313.html?rid=d566de1c-676d-4064-a8c8-2fb172b244c9
    # # C5H9NS = biocrude_cmps['C5H9NS']
    # # C5H9NO = Component('C5H9NO')
    # # C5H9NS.V.l.add_method(C5H9NO.V.l)
    # # C5H9NS.copy_models_from(C5H9NO) #!!! add V.l.
    # # C5H9NS.Tb = 273.15+(151.6+227.18)/2 # avg of ACD and EPIsuite
    # # C5H9NS.Hvap.add_method(38.8e3) # Enthalpy of Vaporization, 38.8±3.0 kJ/mol
    # # C5H9NS.Psat.add_method((3.6+0.0759)/2*133.322) # Vapour Pressure, 3.6±0.3/0.0756 mmHg at 25°C, ACD/EPIsuite
    # # C5H9NS.Hf = -265.73e3 # C5H9NO, https://webbook.nist.gov/cgi/cbook.cgi?ID=C872504&Mask=2

    # # Rough assumption based on the formula
    # biocrude_cmps['7MINDOLE'].Hf = biocrude_cmps['INDOLE'].Hf
    # biocrude_cmps['C30DICAD'].Hf = biocrude_cmps['CHOLESOL'].Hf
    
    # # Components in the aqueous product
    # H2O = htl_cmps.H2O
    # C = Component('C', search_ID='Carbon', particle_size='Soluble',
    #               degradability='Undegradable', organic=False)
    # N = Component('N', search_ID='Nitrogen', particle_size='Soluble',
    #               degradability='Undegradable', organic=False)
    # NH3 = htl_cmps.NH3
    # P = Component('P', search_ID='Phosphorus', particle_size='Soluble',
    #               degradability='Undegradable', organic=False)
    # for i in (C, N, P): i.at_state('l')
    
    # # Components in the gas product
    # CO2 = htl_cmps.CO2
    # CH4 = htl_cmps.CH4
    # C2H6 = htl_cmps.C2H6
    # O2 = htl_cmps.O2
    # N2 = htl_cmps.N2

    # # Other needed components
    # Biofuel = htl_cmps.C16H34.copy('Biofuel') # Tb = 559 K
    # Biobinder = htl_cmps.TRICOSANE.copy('Biobinder') # Tb = 654 K

    # # Compile components
    # biobinder_cmps = Components([
    #     Lipids, Proteins, Carbohydrates, Ash,
    #     Biocrude, HTLaqueous, Hydrochar,
    #     *biocrude_cmps.values(),
    #     H2O, C, N, NH3, P,
    #     CO2, CH4, C2H6, O2, N2,
    #     Biofuel, Biobinder,
    #     ])
    
    # for i in biobinder_cmps:
    #     for attr in ('HHV', 'LHV', 'Hf'):
    #         if getattr(i, attr) is None: setattr(i, attr, 0)
    #     i.default() # default properties to those of water

    biobinder_cmps.compile()
    biobinder_cmps.set_alias('H2O', 'Water')
    biobinder_cmps.set_alias('H2O', '7732-18-5')
    biobinder_cmps.set_alias('Carbohydrates', 'Carbs')
    biobinder_cmps.set_alias('C', 'Carbon')
    biobinder_cmps.set_alias('N', 'Nitrogen')
    biobinder_cmps.set_alias('P', 'Phosphorus')

    if set_thermo: qs_set_thermo(biobinder_cmps)

    return biobinder_cmps