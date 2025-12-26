#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Xuan Wang <easonbiubiu99@gmail.com>
    
    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_components',)

def create_components(set_thermo=True):

    Fe_coagulant = Component(
        ID='Fe_coagulant',
        search_ID='7705-08-0',
        particle_size='particulate',
        degradability='Undegradable',
        organic=False
    )
    
    Food_waste = Component(
        ID='Food_waste',
        search_ID='50-99-7',
        particle_size='Particulate',
        degradability='Biodegardable',
        organic=True
    )
    
# =============================================================================
# gases emissions during the fermentation of sludge and food waste
# =============================================================================
    
    CO2 = Component(
        ID='CO2',
        phase='g',
        particle_size='Dissolved gas',
        degradability='Undegradable',
        organic=False
    )
    
    H2 = Component(
        ID='H2',
        phase='g',
        particle_size='Dissolved gas',
        degradability='Undegradable',
        organic=False
    )
    
    H2S = Component(
        ID='H2S',
        phase='g',
        particle_size='Dissolved gas',
        degradability= 'Undegradable',
        organic=False
    )
    
    CH4 = Component(
        ID='CH4',
        phase='g',
        particle_size='Dissolved gas',
        degradability='Undegradable',
        organic=False
    )
    
# =============================================================================
# water
# =============================================================================
    
    H2O = Component(
        ID='Water',
        particle_size='Soluble',
        degradability='Undegradable',
        organic=False
    )
    
# =============================================================================
# inorganic ions dissolved during the fermentation
# =============================================================================
    
    Fe2 = Component(
        ID='Fe2',
        search_ID='Water',
        phase='l',
        particle_size='Soluble',
        degradability='Undegradable',
        organic=False
    )
    
    Fe3 = Component(
        ID='Fe3',
        search_ID='Water',
        phase='l',
        particle_size='Soluble',
        degradability='Undegradable',
        organic=False
    )
    
    PO4 = Component(
        ID='PO4',
        search_ID='Water',
        phase='l',
        particle_size='Soluble',
        degradability='Undegradable',
        organic=False
    )
    
    Ca2 = Component(
        ID='Ca2',
        search_ID='Water',
        phase='l',
        particle_size='Soluble',
        degradability='Undegradable',
        organic=False
    )
    
    Mg2 = Component(
        ID='Mg2',
        search_ID='Water',
        phase='l',
        particle_size='Soluble',
        degradability='Undegradable',
        organic=False
    )
    
# =============================================================================
# organics released during the fermentation
# =============================================================================
    
    Org = Component(
        ID='Org',
        search_ID='50-99-7',
        particle_size='Soluble',
        degradability='Biodegradable',
        organic=True
    )
    
# =============================================================================
# VFAs, the key acidification products
# =============================================================================
    
    Ac = Component(
        ID='Acetic_acid',
        search_ID='64-19-7',
        particle_size='Soluble',
        degradability='Readily',
        organic=True
    )
    
    Pr = Component(
        ID='Propionic_acid',
        search_ID='79-09-4',
        particle_size='Soluble',
        degradability='Readily',
        organic=True
    )
    
    Bu = Component(
        ID='Butyric_acid',
        search_ID='107-92-6',
        particle_size='Soluble',
        degradability='Readily',
        organic=True
    )
    
    Va = Component(
        ID='Valeric_acid',
        search_ID='109-52-4',
        particle_size='Soluble',
        degradability='Readily',
        organic=True
    )
    
    Lac = Component(
        ID='Lactic_acid',
        search_ID='50-21-5',
        particle_size='Soluble',
        degradability='Readily',
        organic=True
    )
    
    Etoh = Component(
        ID='Ethanol',
        search_ID='64-17-5',
        particle_size='Soluble',
        degradability='Readily',
        organic=True
    )
    
# =============================================================================
# sludge residue to landfill after fermentation
# =============================================================================
    
    Inert = Component(
        ID='Inert',
        phase='s',
        particle_size='Particulate',
        degradability='Undegradable',
        organic=False
    )
    
    # assume 1500 kg/m3
    add_V_from_rho(Inert, 1500)
    Inert.copy_models_from(Chemical('CaCO3'),('Cn',))
    # the stream is quite diluted; assume similar to water
    Inert.mu.add_model(1e-3)
    
    Residue = Component(
        ID='Residue',
        phase='s',
        particle_size='Particulate',
        degradability='Undegradable',
        organic=False
    )
    
    # assume 1500 kg/m3
    add_V_from_rho(Residue, 1500)
    Residue.copy_models_from(Chemical('CaCO3'),('Cn',))
    # the stream is quite diluted; assume similar to water
    Residue.mu.add_model(1e-3)

# =============================================================================
# gases emission dring the calcination, including input (air-O2+N2) and output, CO2 was defined before
# =============================================================================
    
    O2 = Component(
        ID='O2',
        phase='g',
        particle_size='Dissolved gas',
        degradability='Undegradable',
        organic=False
    )
    
    N2 = Component(
        ID='N2',
        phase='g',
        particle_size='Dissolved gas',
        degradability='Undegradable',
        organic=False
    )
    
    SO2 = Component(
        ID='SO2',
        phase='g',
        particle_size='Dissolved gas',
        degradability='Undegradable',
        organic=False        
    )
    
# =============================================================================
# chemicals used during the FePO4 precipitation
# =============================================================================
    
    H2SO4 = Component(
        ID='H2SO4',
        particle_size='Soluble',
        degradability='Undegradable',
        organic=False
    )
    
    H2O2 = Component(
        ID='H2O2',
        particle_size='Soluble',
        degradability='Undegradable',
        organic=False
    )
    
# =============================================================================
# two kinds of chemicals used as the P source
# =============================================================================
    
    H3PO4 = Component(
        ID='H3PO4',
        particle_size='Soluble',
        degradability='Undegradable',
        organic=False
    )
    
    NH4H2PO4 = Component(
        ID='NH4H2PO4',
        particle_size='Soluble',
        formula='NH4H2PO4',
        degradability='Undegradable',
        organic=False
    )
    
# =============================================================================
# two kinds of chemicals used as the Fe source
# =============================================================================
    
    FeCl3 = Component(
        ID='FeCl3',
        particle_size='Soluble',
        formula='FeCl3',
        degradability='Undegradable',
        organic=False
    )
    
    FeSO4_7H2O = Component(
        ID='FeSO4_7H2O',
        search_ID='7782-63-0',
        particle_size='Soluble',
        degradability='Undegradable',
        organic=False
    )
    
# =============================================================================
# output before/after the purification
# =============================================================================
    
    FePO4 = Component(
        ID='FePO4',
        search_ID='10045-86-0',
        phase='s',
        particle_size='Particulate',
        degradability='Undegradable',
        organic=False
    )
    
    FePO4_2H2O = Component(
        ID='FePO4_2H2O',
        formula='FePO4H4O2',    
        phase='s',
        particle_size='Particulate',
        degradability='Undegradable',
        organic=False
    )
    
    # 2770 kg/m3; https://materials.springer.com/isp/crystallographic/docs/sd_0381967 (accessed 2025-12-22)
    add_V_from_rho(FePO4_2H2O, 2770)
    FePO4_2H2O.copy_models_from(Chemical('FePO4'),('Cn',))
    # the stream is quite diluted; assume similar to water
    FePO4_2H2O.mu.add_model(1e-3)

# =============================================================================
# compilation
# =============================================================================
    
    cmps = Components([
        Fe_coagulant, Food_waste, 
        CO2, H2, H2S, CH4,
        H2O, 
        Fe2, Fe3, PO4, Ca2, Mg2,
        Org, Ac, Pr, Bu, Va, Lac, Etoh, Inert, Residue,
        O2, N2, SO2, 
        H2SO4, H2O2, H3PO4, NH4H2PO4, FeCl3, FeSO4_7H2O,
        FePO4, FePO4_2H2O
    ])
    
    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)
    
    cmps.compile()
    cmps.set_alias('Water','H2O')
    if set_thermo: qs_set_thermo(cmps)
    
    return cmps