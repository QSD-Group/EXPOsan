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
from exposan.saf import feedstock_composition, HTL_yields

__all__ = ('create_components',)

def estimate_heating_values(component):
    '''
    Estimate the HHV of a component based on the Dulong's equation (MJ/kg):
        
        HHV [kJ/g] = 33.87*C + 122.3*(H-O/8) + 9.4*S
        
    where C, H, O, and S are the wt% of these elements.
        
    Estimate the LHV based on the HHV as:
        
        LHV [kJ/g] = HHV [kJ/g] – 2.51*(W + 9H)/100
        
    where W and H are the wt% of moisture and H in the fuel
    
    References
    ----------
    [1] https://en.wikipedia.org/wiki/Heat_of_combustion
    [2] https://www.sciencedirect.com/science/article/abs/pii/B9780128203606000072
        
    '''
    atoms = component.atoms
    MW = component.MW
    HHV = (33.87*atoms.get('C', 0)*12 +
           122.3*(atoms.get('H', 0)-atoms.get('O', 0)/8) +
           9.4*atoms.get('S', 0)*32
           )/MW
    LHV = HHV - 2.51*(9*atoms.get('H', 0)/MW)
    
    return HHV*MW*1000, LHV*MW*1000


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
    HTLbiocrude = htl_cmps.Biocrude
    HTLaqueous = htl_cmps.HTLaqueous
    # 43040 mg/L COD
    moisture = feedstock_composition['Water']
    HTLaqueous.i_COD = 43040*moisture/1e6/((1-moisture)*HTL_yields['aqueous'])
    
    HTLchar = htl_cmps.Hydrochar.copy('HTLchar')
    saf_cmps.extend([HTLbiocrude, HTLaqueous, HTLchar])
    
    # Components in the biocrude
    org_kwargs = {
        'particle_size': 'Soluble',
        'degradability': 'Slowly',
        'organic': True,
        }
    biocrude_dct = { # ID, search_ID (CAS#)
        '1E2PYDIN':     '2687-91-4',
        # 'C5H9NS':       '10441-57-3',
        'ETHYLBEN':     '100-41-4',
        '4M-PHYNO':     '106-44-5',
        '4EPHYNOL':     '123-07-9',
        'INDOLE':       '120-72-9',
        '7MINDOLE':     '933-67-5',
        'C14AMIDE':     '638-58-4',
        'C16AMIDE':     '629-54-9',
        'C18AMIDE':     '124-26-5',
        'C16:1FA':      '373-49-9',
        'C16:0FA':      '57-10-3',
        'C18FACID':     '112-80-1',
        'NAPHATH':      '91-20-3',
        'CHOLESOL':     '57-88-5',
        'AROAMINE':     '74-31-7',
        'C30DICAD':     '3648-20-2',
        }
    for ID, search_ID in biocrude_dct.items():
        cmp = Component(ID, search_ID=search_ID, **org_kwargs)
        if not cmp.HHV or not cmp.LHV: 
            HHV, LHV = estimate_heating_values(cmp)
            cmp.HHV = cmp.HHV or HHV
            cmp.LHV = cmp.LHV or LHV
        saf_cmps.append(cmp)
        
    # # Add missing properties
    # # http://www.chemspider.com/Chemical-Structure.500313.html?rid=d566de1c-676d-4064-a8c8-2fb172b244c9
    # C5H9NS = biocrude_cmps['C5H9NS']
    # C5H9NO = Component('C5H9NO')
    # C5H9NS.V.l.add_method(C5H9NO.V.l)
    # C5H9NS.copy_models_from(C5H9NO) #!!! add V.l.
    # C5H9NS.Tb = 273.15+(151.6+227.18)/2 # avg of ACD and EPIsuite
    # C5H9NS.Hvap.add_method(38.8e3) # Enthalpy of Vaporization, 38.8±3.0 kJ/mol
    # C5H9NS.Psat.add_method((3.6+0.0759)/2*133.322) # Vapour Pressure, 3.6±0.3/0.0756 mmHg at 25°C, ACD/EPIsuite
    # C5H9NS.Hf = -265.73e3 # C5H9NO, https://webbook.nist.gov/cgi/cbook.cgi?ID=C872504&Mask=2

    # Components in the biooil
    biooil_IDs = {
        'C4H10', 'TWOMBUTAN', 'NPENTAN', 'TWOMPENTA', 'CYCHEX',
        'HEXANE', 'TWOMHEXAN', 'HEPTANE', 'CC6METH', 'PIPERDIN',
        'TOLUENE', 'THREEMHEPTA', 'OCTANE', 'ETHCYC6', 'ETHYLBEN',
        'OXYLENE', 'C9H20', 'PROCYC6', 'C3BENZ', 'FOURMONAN', 'C10H22',
        'C4BENZ', 'C11H24', 'C10H12', 'C12H26', 'C13H28', 'C14H30',
        'OTTFNA', 'C6BENZ', 'OTTFSN', 'C7BENZ', 'C8BENZ', 'C10H16O4',
        'C15H32', 'C16H34', 'C17H36', 'C18H38', 'C19H40', 'C20H42', 'C21H44',
        'TRICOSANE', 'C24H38O4', 'C26H42O4', 'C30H62',
        }.difference(biocrude_dct.keys())
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
    C3H8 = htl_cmps.C3H8
    O2 = htl_cmps.O2
    N2 = htl_cmps.N2
    CO = htl_cmps.CO
    H2 = htl_cmps.H2
    NH3 = htl_cmps.NH3
    saf_cmps.extend([CO2, CH4, C2H6, C3H8, O2, N2, CO, H2, NH3])
    
    # Surrogate compounds based on the carbon range
    # Tb = 391.35 K (118.2°C)
    Gasoline = Component('Gasoline', search_ID='C8H18', **org_kwargs)
    # Tb = 526.65 K (253.5°C)
    Jet = Component('Jet', search_ID='C14H30', **org_kwargs)
    # Tb = 632.15 K (359°C)
    Diesel = Component('Diesel', search_ID='C21H44', **org_kwargs)
    saf_cmps.extend([Gasoline, Jet, Diesel])

    # Consumables only for cost purposes, thermo data for these components are made up
    sol_kwargs = {
        'phase': 's',
        'particle_size': 'Particulate',
        'degradability': 'Undegradable',
        'organic': False,
        }
    HCcatalyst = Component('HCcatalyst', **sol_kwargs) # Fe-ZSM5
    add_V_from_rho(HCcatalyst, 1500)
    HCcatalyst.copy_models_from(Component('CaCO3'),('Cn',))
    
    HTcatalyst = HCcatalyst.copy('HTcatalyst') # Pd-Al2O3
    
    EOmembrane = HCcatalyst.copy('EOmembrane')
    EOanode = HCcatalyst.copy('EOanode')
    EOcathode = HCcatalyst.copy('EOcathode')
    
    ECmembrane = HCcatalyst.copy('ECmembrane')
    ECanode = HCcatalyst.copy('ECanode')
    ECcathode = HCcatalyst.copy('ECcathode')
    
    saf_cmps.extend([
        HCcatalyst, HTcatalyst,
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
    saf_cmps.set_alias('C', 'Carbon')
    saf_cmps.set_alias('N', 'Nitrogen')
    saf_cmps.set_alias('P', 'Phosphorus')
    saf_cmps.set_alias('K', 'Potassium')
    saf_cmps.set_alias('Biocrude', 'HTLbiocrude')
    saf_cmps.set_alias('HTLchar', 'Hydrochar')

    if set_thermo: qs_set_thermo(saf_cmps)

    return saf_cmps