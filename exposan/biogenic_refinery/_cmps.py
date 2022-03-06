#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# %%

from thermosteam.functional import rho_to_V
from qsdsan import Chemical, Component, Components
from exposan.bwaise import cmps as bw_cmps

__all__ = ('cmps', )

C = Component('C', phase='l', particle_size='Soluble',
              degradability='Undegradable', organic=False)


# NH3 = Component('NH3', measured_as='N',
#                 phase='l', particle_size='Soluble',
#                 degradability='Undegradable', organic=False)

# NonNH3 = Component('NonNH3', search_ID='N', measured_as='N',
#                    phase='l', particle_size='Soluble',
#                    degradability='Undegradable', organic=False,
#                    description='Non-NH3 nitrogen')

# P = Component('P', phase='l', particle_size='Soluble',
#               degradability='Undegradable', organic=False)

# K = Component('K', phase='l', particle_size='Soluble',
#               degradability='Undegradable', organic=False)

# Mg = Component('Mg', phase='l', particle_size='Soluble',
#                degradability='Undegradable', organic=False)

# Ca = Component('Ca', phase='l', particle_size='Soluble',
#                degradability='Undegradable', organic=False)

# H2O = Component('H2O', phase='l', particle_size='Soluble',
#                 degradability='Undegradable', organic=False)

# OtherSS = Component('OtherSS', phase='l', particle_size='Soluble',
#                     degradability='Undegradable', organic=False,
#                     description='Unspecified soluble solids')

# N2O = Component('N2O', phase='g', particle_size='Dissolved gas',
#                 degradability='Undegradable', organic=False)

# CH4 = Component('CH4', phase='g', particle_size='Dissolved gas',
#                 degradability='Slowly', organic=True)

# # Below three are for combustion reactions
# O2 = Component('O2', phase='g', particle_size='Dissolved gas',
#                degradability='Undegradable', organic=False)

# N2 = Component('N2', phase='g', particle_size='Dissolved gas',
#                degradability='Undegradable', organic=False)

# CO2 = Component('CO2', phase='g', particle_size='Dissolved gas',
#                 degradability='Undegradable', organic=False)

# P4O10 = Component.from_chemical('P4O10', Chemical('P4O10'),
#                                 phase='s', particle_size='Particulate',
#                                 degradability='Undegradable', organic=False)
# # The following will lead to an error as it won't be able to copy the V model
# # P4O10 = Component('P4O10', phase='s', particle_size='Particulate',
# #                   degradability='Undegradable', organic=False)

def add_V_from_rho(cmp, rho):
    V_model = rho_to_V(rho, cmp.MW)
    try: cmp.V.add_model(V_model)
    except:
        handle = getattr(cmp.V, cmp.locked_state)
        handle.add_model(V_model)

FilterBag = Component('FilterBag', search_ID='Poly(hexamethylene adipamide)',
                      formula='C12H20N2O2', phase='s', particle_size='Particulate',
                      degradability='Undegradable', organic=False)

# Tissue = Component('Tissue', MW=1, phase='s', particle_size='Particulate',
#                     degradability='Undegradable', organic=False,
#                     description='Tissue for toilet paper')
# # 375 kg/m3 is the average of 250-500 for tissue from
# # https://paperonweb.com/density.htm (accessed 2020-11-12)
# add_V_from_rho(Tissue, 375)

# WoodAsh = Component('WoodAsh', MW=1, phase='s', i_Mg=0.0224, i_Ca=0.3034,
#                     particle_size='Particulate', degradability='Undegradable',
#                     organic=False, description='Wood ash for desiccant')
# add_V_from_rho(WoodAsh, 760)

# for i in (Tissue, WoodAsh):
#     i.copy_models_from(Chemical('Glucose'), ('Cn', 'mu'))

# Struvite = Component('Struvite', search_ID='MagnesiumAmmoniumPhosphate',
#                      formula='NH4MgPO4·H12O6',
#                      phase='s', particle_size='Particulate',
#                      degradability='Undegradable', organic=False)
# # http://www.chemspider.com/Chemical-Structure.8396003.html (accessed 2020-11-19)
# add_V_from_rho(Struvite, 1711)

# HAP = Component('HAP', search_ID='Hydroxyapatite',
#                 phase='s', particle_size='Particulate',
#                 degradability='Undegradable', organic=False)
# # Taking the average of 3.1-3.2 g/cm3 from
# # https://pubchem.ncbi.nlm.nih.gov/compound/Hydroxyapatite (accessed 2020-11-19)
# add_V_from_rho(HAP, 3150)

# Polymer for dewatering (polyacrylamide)
Polymer = Component('Polyacrylamide', formula='C3H5NO',
                    phase='s', particle_size='Particulate',
                    degradability='Slowly', organic=False)

Resin = Component('Polystyrene', search_ID='Polystyrene',
                  formula='C8H8', phase='s', particle_size='Particulate',
                  degradability='Undegradable', organic=False)
# 1200 kg/m3
#
# add_V_from_rho(Resin, 1200)

H2SO4 = Component('H2SO4', phase='l', particle_size='Soluble',
                  degradability='Undegradable', organic=False)

MgOH2 = Component('MagnesiumHydroxide', formula='Mg(OH)2',
                  phase='s', particle_size='Particulate',
                  degradability='Slowly', organic=False)

MgCO3 = Component('MagnesiumCarbonate', formula='MgCO3',
                  phase='s', particle_size='Particulate',
                  degradability='Slowly', organic=False)
# Agricultural Residues
# data collect from source below unless noted
# https://www.sciencedirect.com/science/article/abs/pii/S0016236101001314?via%3Dihub
# Bamboo wood
# moisture content = 20.55 %, caloric value (HHV) = 11.5 MJ/kg
BambooWood = Component('BambooWood', MW=1, phase='s', i_C = 0.4876, i_N = 0.002,
                       particle_size='Particulate',
                       degradability='Undegradable', organic=False)
# 406 kg/m3 is average from:
# https://ojs.cnr.ncsu.edu/index.php/BioRes/article/view/10244 (accessed 2021-04-17)
add_V_from_rho(BambooWood, 406)
BambooWood.HHV = 11.5

# Coconut shell
# moisture content = 20.5 %, caloric value (HHV) = 8.27 MJ/kg
CoconutShell = Component('CoconutShell', phase='s', i_C = 0.5022, i_N = 0.0001,
                         particle_size='Particulate',
                         degradability='Undegradable', organic=False)
# 700 kg/m3 is average from:
# https://aip.scitation.org/doi/pdf/10.1063/1.5127145#:~:text=Actual%20density%20of%20coconut%20shell,600%2D800%20kg%2Fm3. (accessed 2021-04-17)
add_V_from_rho(CoconutShell, 700)
CoconutShell.HHV = 8.27

# Coconut husk
# moisture content = 18.067 %, caloric value (HHV) = 19.77 MJ/kg
CoconutHusk = Component('CoconutHusk', phase='s', i_C = 0.4876, i_N = 0.002,
                        particle_size='Particulate',
                        degradability='Undegradable', organic=False)
# 69 kg/m3 is average for crushed husk from:
# https://demelenterprises.com/Cocopeat/Cocopeat.pdf (accessed 2021-04-17)
add_V_from_rho(CoconutHusk, 69)
CoconutHusk.HHV = 19.77

# Rice husk
# moisture content = 14.693 %, caloric value (HHV) = 8.47 MJ/kg
RiceHusk = Component('RiceHusk', phase='s', i_C = 0.385, i_N = 0.0045,
                      particle_size='Particulate',
                      degradability='Undegradable', organic=False)
# 358 kg/m3 is average from:
# https://thescipub.com/pdf/ajassp.2012.1757.1768.pdf (accessed 2021-04-17)
add_V_from_rho(RiceHusk, 358)
RiceHusk.HHV = 8.47

# Corn stover
# moisture content = 18.5 %, caloric value (HHV) = 46.5 MJ/kg
CornStover = Component('CornStoverk', phase='s', i_C = 0.465, i_N = 0.0056,
                       particle_size='Particulate',
                       degradability='Undegradable', organic=False)
# 90 kg/m3 is average from:
# https://core.ac.uk/download/pdf/38931685.pdf (accessed 2021-04-17)
add_V_from_rho(CornStover, 90)
CornStover.HHV = 46.5

for i in  (BambooWood, CoconutShell, CoconutHusk, RiceHusk, CornStover):
    i.copy_models_from(Chemical('Glucose'), ('Cn', 'mu'))


# for cmp in (NonNH3, C, P, K, Mg, Ca, OtherSS):
#     cmp.default()
#     cmp.copy_models_from(H2O, ('sigma', 'epsilon', 'kappa', 'Cn', 'mu'))
#     add_V_from_rho(cmp, 1e3) # assume the same density as water

cmps = Components((*bw_cmps, C, FilterBag, Polymer, Resin, H2SO4, MgOH2, MgCO3,
                   BambooWood, CoconutShell, CoconutHusk, RiceHusk, CornStover))

# cmps = Components((NH3, NonNH3, C, P, K, Mg, Ca, H2O, OtherSS, N2O, CH4, O2, N2,
#                    CO2, P4O10, FilterBag, Tissue, WoodAsh, Struvite, HAP,
#                    Polymer, Resin, H2SO4, MgOH2, MgCO3,
#                    BambooWood, CoconutShell, CoconutHusk, RiceHusk, CornStover))

for i in cmps:
    if i.HHV is None: i.HHV = 0
    if i.LHV is None: i.LHV = 0
    if i.Hf is None: i.Hf = 0
cmps.compile()

cmps.set_synonym('H2O', 'Water')