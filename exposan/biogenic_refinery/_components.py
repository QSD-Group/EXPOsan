#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>
    Lane To <lane20@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''


# %%

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_components', )

def create_components(set_thermo=True):

    # NH3 through HAP are copied from Bwaise components    

    NH3 = Component('NH3', measured_as='N',
                    phase='l', particle_size='Soluble',
                    degradability='Undegradable', organic=False)

    NonNH3 = Component('NonNH3', formula='N', measured_as='N',
                       phase='l', particle_size='Soluble',
                       degradability='Undegradable', organic=False,
                       description='Non-NH3 nitrogen')

    P = Component('P', phase='l', particle_size='Soluble',
                  degradability='Undegradable', organic=False)

    K = Component('K', phase='l', particle_size='Soluble',
                  degradability='Undegradable', organic=False)

    Mg = Component('Mg', phase='l', particle_size='Soluble',
                   degradability='Undegradable', organic=False)

    Ca = Component('Ca', phase='l', particle_size='Soluble',
                   degradability='Undegradable', organic=False)

    H2O = Component('H2O', phase='l', particle_size='Soluble',
                    degradability='Undegradable', organic=False)

    OtherSS = Component('OtherSS', phase='l', particle_size='Soluble',
                        degradability='Undegradable', organic=False,
                        description='Unspecified soluble solids')

    N2O = Component('N2O', phase='g', particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)

    CH4 = Component('CH4', phase='g', particle_size='Dissolved gas',
                    degradability='Slowly', organic=True)

    # Below three are for combustion reactions
    O2 = Component('O2', phase='g', particle_size='Dissolved gas',
                   degradability='Undegradable', organic=False)

    N2 = Component('N2', phase='g', particle_size='Dissolved gas',
                   degradability='Undegradable', organic=False)

    CO2 = Component('CO2', phase='g', particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)

    P4O10 = Component.from_chemical('P4O10', Chemical('P4O10'),
                                    phase='s', particle_size='Particulate',
                                    degradability='Undegradable', organic=False)
    # https://en.wikipedia.org/wiki/Phosphorus_pentoxide (accessed 2022-09-08)
    add_V_from_rho(P4O10, 2390)
    # The following will lead to an error as it won't be able to copy the V model
    # P4O10 = Component('P4O10', phase='s', particle_size='Particulate',
    #                   degradability='Undegradable', organic=False)

    Tissue = Component('Tissue', MW=1, phase='s', particle_size='Particulate',
                        degradability='Undegradable', organic=False,
                        description='Tissue for toilet paper')
    # 375 kg/m3 is the average of 250-500 for tissue from
    # https://paperonweb.com/density.htm (accessed 2020-11-12)
    add_V_from_rho(Tissue, 375)

    WoodAsh = Component('WoodAsh', MW=1, phase='s', i_Mg=0.0224, i_Ca=0.3034,
                        particle_size='Particulate', degradability='Undegradable',
                        organic=False, description='Wood ash for desiccant')
    add_V_from_rho(WoodAsh, 760)

    for i in (Tissue, WoodAsh):
        i.copy_models_from(Chemical('Glucose'), ('Cn', 'mu'))

    Struvite = Component('Struvite', search_ID='MagnesiumAmmoniumPhosphate',
                         formula='NH4MgPO4·H12O6',
                         phase='s', particle_size='Particulate',
                         degradability='Undegradable', organic=False)
    # http://www.chemspider.com/Chemical-Structure.8396003.html (accessed 2020-11-19)
    add_V_from_rho(Struvite, 1711)

    HAP = Component('HAP', search_ID='Hydroxyapatite',
                    phase='s', particle_size='Particulate',
                    degradability='Undegradable', organic=False)
    # Taking the average of 3.1-3.2 g/cm3 from
    # https://pubchem.ncbi.nlm.nih.gov/compound/Hydroxyapatite (accessed 2020-11-19)
    add_V_from_rho(HAP, 3150)


    C = Component('C', phase='l', particle_size='Soluble',
                  degradability='Undegradable', organic=False)
    
    # Ash, VolatileMatter, and FixedCarbon are components exclusive to the biochar stream in the biogenic refinery carbonizer base
    
    Ash = Component('Ash', phase='s', particle_size='Particulate', 
                    degradability='Undegradable', organic=False)
    
    VolatileMatter = Component('VolatileMatter', phase='s', particle_size='Particulate', 
                    degradability='Undegradable', organic=False)
    
    FixedCarbon = Component('FixedCarbon', phase='s', particle_size='Particulate', 
                    degradability='Undegradable', organic=False)

    FilterBag = Component('FilterBag', search_ID='Poly(hexamethylene adipamide)',
                          formula='C12H20N2O2', phase='s', particle_size='Particulate',
                          degradability='Undegradable', organic=False)
    # 969 kg/m3 is average from:
    # https://www.chembk.com/en/chem/Poly(hexamethylene%20adipamide) (accessed 2022-07-25)
    add_V_from_rho(FilterBag, 969)

    # Polymer for dewatering (polyacrylamide)
    Polymer = Component('Polyacrylamide', formula='C3H5NO',
                        phase='s', particle_size='Particulate',
                        degradability='Slowly', organic=False)

    Resin = Component('Polystyrene', search_ID='Polystyrene',
                      formula='C8H8', phase='s', particle_size='Particulate',
                      degradability='Undegradable', organic=False)
    # See this issue for density: https://github.com/QSD-Group/EXPOsan/issues/15
    add_V_from_rho(Resin, 750)

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
    # moisture content = 11.50 %, caloric value (HHV) = 20.547 MJ/kg
    BambooWood = Component('BambooWood', MW=1, phase='s', i_C = 0.4876, i_N = 0.002,
                           f_Amass_Totmass = 0.0195, particle_size='Particulate',
                           degradability='Undegradable', organic=False)
    # 406 kg/m3 is average from:
    # https://ojs.cnr.ncsu.edu/index.php/BioRes/article/view/10244 (accessed 2021-04-17)
    add_V_from_rho(BambooWood, 406)
    BambooWood.HHV = 20.547

    # Coconut shell
    # moisture content = 8.27 %, caloric value (HHV) = 20.498 MJ/kg
    CoconutShell = Component('CoconutShell', phase='s', i_C = 0.5022, i_N = 0.0001,
                             f_Amass_Totmass = 0.0071, particle_size='Particulate',
                             degradability='Undegradable', organic=False)
    # 700 kg/m3 is average from:
    # https://aip.scitation.org/doi/pdf/10.1063/1.5127145#:~:text=Actual%20density%20of%20coconut%20shell,600%2D800%20kg%2Fm3. (accessed 2021-04-17)
    add_V_from_rho(CoconutShell, 700)
    CoconutShell.HHV = 20.498

    # Coconut husk
    # moisture content = 19.777 %, caloric value (HHV) = 18.067 MJ/kg
    CoconutHusk = Component('CoconutHusk', phase='s', i_C = 0.4876, i_N = 0.002,
                            f_Amass_Totmass = 0.0516, particle_size='Particulate',
                            degradability='Undegradable', organic=False)
    # 69 kg/m3 is average for crushed husk from:
    # https://demelenterprises.com/Cocopeat/Cocopeat.pdf (accessed 2021-04-17)
    add_V_from_rho(CoconutHusk, 69)
    CoconutHusk.HHV = 18.067

    # Rice husk
    # moisture content = 8.47 %, caloric value (HHV) = 14.693 MJ/kg
    RiceHusk = Component('RiceHusk', phase='s', i_C = 0.385, i_N = 0.0045,
                          f_Amass_Totmass = 0.2124, particle_size='Particulate',
                          degradability='Undegradable', organic=False)
    # 358 kg/m3 is average from:
    # https://thescipub.com/pdf/ajassp.2012.1757.1768.pdf (accessed 2021-04-17)
    add_V_from_rho(RiceHusk, 358)
    RiceHusk.HHV = 14.693

    # Corn stover
    # moisture content = 18.5 %, caloric value (HHV) = 46.5 MJ/kg
    CornStover = Component('CornStover', phase='s', i_C = 0.465, i_N = 0.0056,
                           f_Amass_Totmass = 0.109, particle_size='Particulate',
                           degradability='Undegradable', organic=False)
    # Ash content is from: https://doi.org/10.1111/gcbb.12183
    # 90 kg/m3 is average from:
    # https://core.ac.uk/download/pdf/38931685.pdf (accessed 2021-04-17)
    add_V_from_rho(CornStover, 90)
    CornStover.HHV = 46.5

    for i in (BambooWood, CoconutShell, CoconutHusk, RiceHusk, CornStover):
        i.copy_models_from(Chemical('Glucose'), ('Cn', 'mu'))

    # Reuse components in the bwaise module for consistency
    cmps = Components((NH3, NonNH3, P, K, Mg, Ca, H2O, OtherSS, N2O, CH4, O2, N2,
                       CO2, P4O10, Tissue, WoodAsh, Struvite, HAP, C, Ash, 
                       VolatileMatter, FixedCarbon, FilterBag, Polymer, Resin, 
                       H2SO4, MgOH2, MgCO3, BambooWood, CoconutShell, CoconutHusk, 
                       RiceHusk, CornStover))

    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_alias('H2O', 'Water')
    if set_thermo: qs_set_thermo(cmps)

    return cmps