#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    BlueS928

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = (
    'ammonium_sulfate_solution_molarity',
    'NH3_solution_molarity',
    'create_components',
    )

ammonium_sulfate_solution_molarity = 0.835 # mol/L
# this is the average of the 0.4 M to 1.3 M TAN in ammonium salt concetration from Kogler et al.
NH3_solution_molarity = 1.76 # mol/L
# NH3 product solution concentration from Kogler et al.

def create_components(set_thermo=True):
    H2O = Component('H2O', particle_size='Soluble',
                    degradability='Undegradable', organic=False)

    NH3 = Component('NH3', particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)

    O2 = Component('O2', phase='g', particle_size='Dissolved gas',
                   degradability='Undegradable', organic=False)

    N2 = Component('N2', phase='g', particle_size='Dissolved gas',
                   degradability='Undegradable', organic=False)

    CH4 = Component('CH4', phase='g', particle_size='Dissolved gas',
                    degradability='Slowly', organic=True)

    CO2 = Component('CO2', phase='g', particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)

    NH42SO4 = Component('NH42SO4', search_ID='AmmoniumSulfate',
                        phase='l', particle_size='Soluble',
                        degradability='Undegradable', organic=False)
    add_V_from_rho(NH42SO4, 1770) # density in kg/m3 from PubChem

    H2SO4 = Component('H2SO4', phase='l', particle_size='Soluble',
                      degradability='Undegradable', organic=False)

    CitricAcid = Component('CitricAcid', search_ID='CitricAcid',
                           phase='l', particle_size='Soluble',
                           degradability='Undegradable', organic=True)

    LacticAcid = Component('LacticAcid', search_ID='LacticAcid',
                           phase='l', particle_size='Soluble',
                           degradability='Undegradable', organic=True)

    MaleicAcid = Component('MaleicAcid', search_ID='MaleicAcid',
                           phase='l', particle_size='Soluble',
                           degradability='Undegradable', organic=True)

    NH4Citrate = Component('NH4Citrate', phase='l',
                           formula='C6H14N2O7', particle_size='Soluble',
                           degradability='Undegradable', organic=True)
    # Diammonium citrate density: 1.48 g/cm3, MediaDive ingredient data.
    add_V_from_rho(NH4Citrate, 1480)
    # Provisional S from parent citric acid
    NH4Citrate.copy_models_from(Chemical('CitricAcid'), ('mu', 'Cn'))
    
    NH4Lactate = Component('NH4Lactate', phase='l',
                           formula='C3H9NO3', particle_size='Soluble',
                           degradability='Undegradable', organic=True)
    # Ammonium lactate specific gravity: 1.2 at 59 F, NOAA CAMEO/PubChem.
    add_V_from_rho(NH4Lactate, 1200)
    # Provisional S from parent lactic acid
    NH4Lactate.copy_models_from(Chemical('LacticAcid'), ('mu', 'Cn'))
    
    NH4Maleate = Component('NH4Maleate', phase='l',
                           formula='C4H10N2O4', particle_size='Soluble',
                           degradability='Undegradable', organic=True)
    # Diammonium maleate density from ALFA Chemistry product specs
    add_V_from_rho(NH4Maleate, 1641)
    # Provisional S from parent maleic acid
    NH4Maleate.copy_models_from(Chemical('MaleicAcid'), ('mu', 'Cn'))
    
    NaOH = Component('NaOH', phase='s', particle_size='Particulate',
                     degradability='Undegradable', organic=False)

    cmps = Components([H2O, NH3, O2, N2, CH4, CO2, NH42SO4, H2SO4,
                       CitricAcid, LacticAcid, MaleicAcid, NH4Citrate,
                       NH4Lactate, NH4Maleate, NaOH])

    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_alias('H2O', 'Water')
    cmps.set_alias('NH42SO4', 'AmmoniumSulfate')
    cmps.set_alias('NH42SO4', '(NH4)2SO4')
    cmps.set_alias('NH4Citrate', 'AmmoniumCitrate')
    cmps.set_alias('NH4Citrate', 'DiammoniumCitrate')
    cmps.set_alias('NH4Citrate', '(NH4)2HC6H5O7')
    cmps.set_alias('NH4Lactate', 'AmmoniumLactate')
    cmps.set_alias('NH4Lactate', 'NH4C3H5O3')
    cmps.set_alias('NH4Maleate', 'AmmoniumMaleate')
    cmps.set_alias('NH4Maleate', '(NH4)2C4H2O4')
    cmps.set_alias('NaOH', 'SodiumHydroxide')

    if set_thermo: qs_set_thermo(cmps)

    return cmps
