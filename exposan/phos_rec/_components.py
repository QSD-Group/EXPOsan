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
from qsdsan import Component, Components

Fe_coagulant = Component(
    ID='Fe_coagulant',
    phase='s',
    particle_size='particulate',
    formula='FeCl3',
    degradability='Undegradable',
    organic=False
)

    #TODO: define the Fe_sludge latter Fe_sludge = Component(ID='Fe_slduge',phase='s',particle_size='particulate',degradability='Bidegradable',organic=True)
    #TODO: CEPS_supernatant need to be defined
    #TODO: Residues after fermentaion need to be defined
    
Food_waste = Component(
    ID='Food_waste',
    phase='s',
    particle_size='Particulate',
    # TODO: change the food waste formula
    formula= 'C6H12O6',
    degradability='Biodegardable',
    organic=True
)

    # TODO:Gases emissions during the fermentation of sludge and food_waste
    
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
    particle_size='Dissloved gas',
    degradability='Undegradable',
    organic=False
)

H2S = Component(
    ID='H2S',
    phase='g',
    particle_size='Dissloved gas',
    degradability= 'Undegradable',
    organic=False
)

CH4 = Component(
    ID='CH4',
    phase='g',
    particle_size='Dissloved gas',
    degradability='Undegradable',
    organic=False
)

    # TODO:Water
    
H2O = Component(
    ID='Water',
    phase='l',
    particle_size='Soluble',
    degradability='Undegradable',
    organic=False
)
    #TODO:Inorganic ions dissloved during the fermentation
    
Fe2 = Component(
    ID='Fe2',
    phase='l',
    particle_size='Solbule',
    formula='Fe',
    degradability='Undegradable',
    organic=False
)

Fe3 = Component(
    ID='Fe3',
    phase='l',
    particle_size='Solbule',
    formula='Fe',
    degradability='Undegradable',
    organic=False
)

PO4 = Component(
    ID='PO4',
    phase='l',
    particle_size='Soluble',
    formula='PO4',
    degradability='Undegradable',
    organic=False
)

Ca2 = Component(
    ID='Ca2',
    phase='l',
    particle_size='Soluble',
    formula='Ca',
    degradability='Undegradable',
    organic=False
)

Mg1 = Component(
    ID='Mg1',
    phase='l',
    particle_size='Soluble',
    formula='Mg',
    degradability='Undegradable',
    organic=False
)

    # TODO:Organics released during the fermentation
    
Org = Component(
    ID='Org',
    phase='l',
    particle_size='Soluble',
    formula='CH2O',
    degradability='Biodegardable',
    organic=True
)
   #TODO:VFAs, the key acidification products
   
Ac = Component(
    ID='Ac',
    phase='l',
    particle_size='Soluble',
    formula='C2H4O2',
    name='Acetate',
    degradability='Readily',
    organic=True
)

Pr = Component(
    ID='Pr',
    phase='l',
    particle_size='Soluble',
    formula='C3H6O2',
    name='Propionate',
    degradability='Readily',
    organic=True
)

Bu = Component(
    ID='Bu',
    phase='l',
    particle_size='Soluble',
    formula='C4H8O2',
    degradability='Readily',
    organic=True
)

Va = Component(
    ID='Va',
    phase='l',
    particle_size='Soluble',
    formula='C5H10O2',
    name='Valerate',
    degradability='Readily',
    organic=True
)

Lac = Component(
    ID='Lac',
    phase='l',
    particle_size='Soluble',
    formula='C3H6O3',
    name='Lactate',
    degradability='Readily',
    organic=True
)

Etoh = Component(
    ID='Etoh',
    phase='l',
    particle_size='Soluble',
    formula='C2H6O',
    name='Ethanol',
    degradability='Readily',
    organic=True
)
    # TODO:sludge residue to landfill after fermentation
    
residue = Component(
    ID='Sludge',
    phase='s',
    particle_size='Particulate',
    degradability='slowly',
    organic=True
)

    # TODO:Gases emission dring the calcination,including input (air-O2) and output
    
O2 = Component(
    ID='O2',
    phase='g',
    particle_size='Dissloved gas',
    degradability='Undegradable',
    organic=False
)

SO2 = Component(
    ID='SO2',
    phase='g',
    particle_size='Dissloved gas',
    degradability='Undegradable',
    organic=False        
)

Gas_H2O = Component(
    ID='Water',
    phase='g',
    particle_size='Dissloved gas',
    degradability='Undegradable',
    organic=False
)

    # TODO:chemicals used during the FePO4 precipitaion
    
H2SO4 = Component(
    ID='H2SO4',
    phase='l',
    particle_size='Soluble',
    degradability='Undegradable',
    organic=False
)

H2O2 = Component(
    ID='H2O2',
    phase='l',
    particle_size='Soluble',
    degradability='Undegradable',
    organic=False
)

    # TODO:two kinds of chemicals used as the P source

H3PO4 = Component(
    ID='H3PO4',
    phase='l',
    particle_size='Soluble',
    degradability='Undegradable',
    organic=False
)

NH4H2PO4 = Component(
    ID='NH4H2PO4',
    phase='Partculate',
    particle_size='Soluble',
    formula='NH4H2PO4',
    degradability='Undegradable',
    organic=False
)

    # TODO:two kinds of chemicals used as the Fe source

FeCl3 = Component(
    ID='FeCl3',
    phase='Particulate',
    particle_size='Soluble',
    formula='FeCl3',
    degradability='Undegradable',
    organic=False
)

FeSO4_7H2O = Component(
    ID='FeSO4_7H2O',
    phase='Particulate',
    particle_size='Soluble',
    formula='FeSO4_7H2O',
    degradability='Undegradable',
    organic=False
)

    # TODO:output before/after the purification
    
Fe_P_precipitate = Component(
    ID='Fe_P_precipitate',
    phase='s',
    particle_size='Particulate',
    degradability='Undegradable',
    organic=True
)

FePO4 = Component(
    ID='FePO4',
    phase='s',
    particle_size='Particulate',
    degradability='Undegradable',
    organic=False
)

cmps = Components([Fe_coagulant, Food_waste, 
                   CO2, H2, H2S, CH4,
                   H2O, 
                   Fe2, Fe3, PO4, Ca2, Mg1,
                   Org, Ac, Pr, Bu, Va, Lac, Etoh, residue,
                   O2, SO2, Gas_H2O, 
                   H2SO4, H2O2, H3PO4, NH4H2PO4, FeCl3, FeSO4_7H2O,
                   Fe_P_precipitate, FePO4])
cmps.compile()
cmps.set_alias('H2O','Water')
