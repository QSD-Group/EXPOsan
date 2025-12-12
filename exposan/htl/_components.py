#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References:

(1) Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
    Quantitative Evaluation of an Integrated System for Valorization of
    Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
    Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
    https://doi.org/10.1021/acs.est.8b04035.
    
(2) Leow, S.; Witter, J. R.; Vardon, D. R.; Sharma, B. K.;
    Guest, J. S.; Strathmann, T. J. Prediction of Microalgae Hydrothermal
    Liquefaction Products from Feedstock Biochemical Composition.
    Green Chem. 2015, 17 (6), 3584–3599. https://doi.org/10.1039/C5GC00574D.
    
(3) Snowden-Swan, L. J.; Zhu, Y.; Jones, S. B.; Elliott, D. C.; Schmidt, A. J.; 
    Hallen, R. T.; Billing, J. M.; Hart, T. R.; Fox, S. P.; Maupin, G. D. 
    Hydrothermal Liquefaction and Upgrading of Municipal Wastewater Treatment 
    Plant Sludge: A Preliminary Techno-Economic Analysis; 
    PNNL--25464, 1258731; 2016; https://doi.org/10.2172/1258731.

(4) Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
    Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
    Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
    Process Design and Economics for the Conversion of Algal Biomass to
    Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
    PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
'''

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_components',)

def create_components(set_thermo=True):
    
    Sludge_lipid = Component('Sludge_lipid', phase='s',
                              particle_size='Particulate',
                              formula='C56H95O24N9P',
                              degradability='Undegradable',
                              organic=True)
    add_V_from_rho(Sludge_lipid, 1400)
    # https://www.climate-policy-watcher.org/wastewater-sludge/physical-
    # and-biological-properties.html (accessed 2022-10-23)
    # https://web.deu.edu.tr/atiksu/ana52/wdesign06.html (accessed 2022-10-23)
    Sludge_lipid.HHV = 22.0*10**6*Sludge_lipid.MW/1000  #Li et al., 2018
    Sludge_lipid.Cn.add_model(1.25*10**3*Sludge_lipid.MW/1000) 
    # Leow et al., 2015
    Sludge_lipid.mu.add_model(6000)
    # made up value, so that HTL.ins[0].nu = 0.03 m2/s ~30000 cSt
    # (NREL 2013 appendix B)
    
    Sludge_protein = Component('Sludge_protein', phase='s',
                                particle_size='Particulate',
                                formula='C56H95O24N9P',
                                degradability='Undegradable',
                                organic=True)
    add_V_from_rho(Sludge_protein, 1400)
    Sludge_protein.HHV = 22.0*10**6*Sludge_protein.MW/1000
    Sludge_protein.Cn.add_model(1.25*10**3*Sludge_protein.MW/1000)
    Sludge_protein.mu.add_model(6000)
    # made up value, so that HTL.ins[0].nu = 0.03 m2/s ~30000 cSt
    # (NREL 2013 appendix B)
    
    Sludge_carbo = Component('Sludge_carbo', phase='s',
                              particle_size='Particulate',
                              formula='C56H95O24N9P',
                              degradability='Undegradable',
                              organic=True)
    add_V_from_rho(Sludge_carbo, 1400)
    Sludge_carbo.HHV = 22.0*10**6*Sludge_carbo.MW/1000
    Sludge_carbo.Cn.add_model(1.25*10**3*Sludge_carbo.MW/1000)
    Sludge_carbo.mu.add_model(6000)
    # made up value, so that HTL.ins[0].nu = 0.03 m2/s ~30000 cSt
    # (NREL 2013 appendix B)
    
    # ash in sludge (not separated), therefore set organic to True
    Sludge_ash = Component('Sludge_ash', phase='s',
                            particle_size='Particulate',
                            formula='C56H95O24N9P',
                            degradability='Undegradable',
                            organic=True)
    add_V_from_rho(Sludge_ash, 1400)
    Sludge_ash.HHV = 22.0*10**6*Sludge_ash.MW/1000
    Sludge_ash.Cn.add_model(1.25*10**3*Sludge_ash.MW/1000)
    Sludge_ash.mu.add_model(6000)
    # made up value, so that HTL.ins[0].nu = 0.03 m2/s ~30000 cSt
    # (NREL 2013 appendix B)

    Hydrochar = Component('Hydrochar', phase='s', particle_size='Particulate',
                        degradability='Undegradable', organic=False)
    # assume a bulk density of 700 kg/m3; the density does not affect transportation (if any) which is based on weight
    add_V_from_rho(Hydrochar, 700)
    Hydrochar.copy_models_from(Chemical('CaCO3'),('Cn',))
    
    # when calculating price using HHV, use a separate HHV value but not directly from this Componenet
    # use palmitamide to represent biocrude
    Biocrude = Component('Biocrude', search_ID='629-54-9',
                         particle_size='Soluble', degradability='Slowly',
                         organic=True)
    
    HTLaqueous = Component('HTLaqueous', search_ID='water', particle_size='Soluble',
                           degradability='Undegradable', organic=False)
    
    Struvite = Component('Struvite', search_ID='MagnesiumAmmoniumPhosphate',
                         formula='NH4MgPO4·H12O6', phase='s',
                         particle_size='Particulate', 
                         degradability='Undegradable',
                         organic=False)
    add_V_from_rho(Struvite, 1710)
    # http://webmineral.com/data/Struvite.shtml#.YzYvqOzMIiM
    # (accessed 2022-9-30)
    
    Residue = Component('Residue', phase='s', particle_size='Particulate',
                        degradability='Undegradable', organic=False)
    add_V_from_rho(Residue, 1500)  # assume 1500kg/m3
    Residue.copy_models_from(Chemical('CaCO3'),('Cn',)) #CaCO3?
    
    H2O = Component('H2O', particle_size='Soluble',
                    degradability='Undegradable', organic=False)
    
    # intentionally using the molecular formula and MW of water for C/N/P
    C = Component('C', search_ID='water', particle_size='Soluble',
                  degradability='Undegradable', organic=False)
    C._CAS = 'C'
    
    N = Component('N', search_ID='water', particle_size='Soluble',
                  degradability='Undegradable', organic=False)
    N._CAS = 'N'
    
    P = Component('P', search_ID='water', particle_size='Soluble',
                  degradability='Undegradable', organic=False)
    P._CAS = 'P'
    
    O2 = Component('O2', phase='g', particle_size='Dissolved gas',
                   degradability='Undegradable', organic=False)
    
    N2 = Component('N2', phase='g', particle_size='Dissolved gas',
                   degradability='Undegradable', organic=False)
    
    N2O = Component('N2O', phase='g', particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)
    
    CH4 = Component('CH4', phase='g', particle_size='Dissolved gas',
                    degradability='Slowly', organic=True)
    
    C2H6 = Component('C2H6', phase='g', search_ID='ethane',
                     particle_size='Dissolved gas', degradability='Slowly',
                     organic=True)
    
    C3H8 = Component('C3H8', phase='g', particle_size='Dissolved gas',
                     degradability='Slowly', organic=True)
    
    # CH4, C2H6, and C3H8 are gas phase only
    
    CO2 = Component('CO2', phase='g', particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)
    
    CO = Component('CO', phase='g', particle_size='Dissolved gas',
                   degradability='Undegradable', organic=False)
    
    H2 = Component('H2', phase='g', particle_size='Dissolved gas',
                   degradability='Undegradable', organic=False)
    
    NH3 = Component('NH3', particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)
    
    H2SO4 = Component('H2SO4', phase='l', particle_size='Soluble',
                      degradability='Undegradable', organic=False)
    
    H3PO4 = Component('H3PO4', phase='l', particle_size='Soluble',
                      degradability='Undegradable', organic=False)
    
    CaO = Component('CaO', phase='s', particle_size='Particulate',
                    degradability='Undegradable', organic=False)
    
    CaOH2 = Component('CaOH2', search_ID='1305-62-0',
                      phase='s', particle_size='Particulate',
                      degradability='Undegradable', organic=False)
    
    CaCO3 = Component('CaCO3', phase='s', particle_size='Particulate',
                      degradability='Undegradable', organic=False)
    
    MgCl2 = Component('MgCl2', phase='l', particle_size='Soluble',
                      degradability='Undegradable', organic=False)
    
    MgO = Component('MgO', phase='l', particle_size='Particulate',
                      degradability='Undegradable', organic=False)
    
    NaOH = Component('NaOH', phase='l', particle_size='Soluble',
                     degradability='Undegradable', organic=False)
    
    NH42SO4 = Component('NH42SO4', phase='l', particle_size='Soluble',
                        degradability='Undegradable', organic=False)
    add_V_from_rho(NH42SO4, 1770)
    # https://en.wikipedia.org/wiki/Ammonium_sulfate (accessed 2022-9-30)
    
    NH4Cl = Component('NH4Cl', phase='l', particle_size='Soluble',
                        degradability='Undegradable', organic=False)

    # below are biooil-related components
    # CH4, C2H6, and C3H8 are gas phase only
    C4H10 = Component('C4H10', particle_size='Dissolved gas',
                   degradability='Slowly', organic=True)
    
    TWOMBUTAN = Component('TWOMBUTAN', search_ID='78-78-4',
                          particle_size='Soluble', degradability='Slowly',
                          organic=True)
    
    NPENTAN = Component('NPENTAN', search_ID='109-66-0',
                        particle_size='Soluble', degradability='Slowly',
                        organic=True)

    TWOMPENTA = Component('TWOMPENTA', search_ID='107-83-5',
                          particle_size='Soluble', degradability='Slowly',
                          organic=True)
    
    CYCHEX = Component('CYCHEX', search_ID='110-82-7', particle_size='Soluble',
                        degradability='Slowly', organic=True)
    
    
    HEXANE = Component('HEXANE', search_ID='110-54-3', particle_size='Soluble',
                        degradability='Slowly', organic=True)

    TWOMHEXAN = Component('TWOMHEXAN', search_ID='591-76-4',
                          particle_size='Soluble', degradability='Slowly',
                          organic=True)

    HEPTANE = Component('HEPTANE', search_ID='142-82-5',
                        particle_size='Soluble', degradability='Slowly',
                        organic=True)

    CC6METH = Component('CC6METH', search_ID='108-87-2',
                        particle_size='Soluble', degradability='Slowly',
                        organic=True)
    
    PIPERDIN = Component('PIPERDIN', search_ID='110-89-4',
                         particle_size='Soluble', degradability='Slowly',
                         organic=True)
    
    TOLUENE = Component('TOLUENE', search_ID='108-88-3',
                        particle_size='Soluble', degradability='Slowly',
                        organic=True)
    
    THREEMHEPTA = Component('THREEMHEPTA', search_ID='589-81-1',
                            particle_size='Soluble', degradability='Slowly',
                            organic=True)

    OCTANE = Component('OCTANE', search_ID='111-65-9', particle_size='Soluble',
                        degradability='Slowly', organic=True)

    ETHCYC6 = Component('ETHCYC6', search_ID='1678-91-7',
                        particle_size='Soluble', degradability='Slowly',
                        organic=True)

    ETHYLBEN = Component('ETHYLBEN', search_ID='100-41-4',
                         particle_size='Soluble', degradability='Slowly',
                         organic=True)

    OXYLENE = Component('OXYLENE', search_ID='95-47-6',
                        particle_size='Soluble', degradability='Slowly',
                        organic=True)

    C9H20 = Component('C9H20', search_ID='111-84-2', particle_size='Soluble',
                      degradability='Slowly', organic=True)

    PROCYC6 = Component('PROCYC6', search_ID='1678-92-8',
                        particle_size='Soluble', degradability='Slowly',
                        organic=True)

    C3BENZ = Component('C3BENZ', search_ID='103-65-1', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    FOURMONAN = Component('FOURMONAN', search_ID='17301-94-9',
                          particle_size='Soluble', degradability='Slowly',
                          organic=True)

    C10H22 = Component('C10H22', search_ID='124-18-5', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    C4BENZ = Component('C4BENZ', search_ID='104-51-8', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    C11H24 = Component('C11H24', search_ID='1120-21-4',
                       particle_size='Soluble', degradability='Slowly',
                       organic=True)

    C10H12 = Component('C10H12', search_ID='119-64-2', particle_size='Soluble',
                       degradability='Slowly', organic=True)
    # use 1234NA (OTTFNA)'s search_ID to replace since they are both C10H12

    C12H26 = Component('C12H26', search_ID='112-40-3', particle_size='Soluble',
                       degradability='Slowly', organic=True)
    
    C13H28 = Component('C13H28', search_ID='629-50-5', particle_size='Soluble',
                       degradability='Slowly', organic=True)
    
    C14H30 = Component('C14H30', search_ID='629-59-4', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    OTTFNA = Component('OTTFNA', search_ID='119-64-2', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    C6BENZ = Component('C6BENZ', search_ID='1077-16-3',
                       particle_size='Soluble', degradability='Slowly',
                       organic=True)

    OTTFSN = Component('OTTFSN', search_ID='1680-51-9',
                       particle_size='Soluble', degradability='Slowly',
                       organic=True)

    C7BENZ = Component('C7BENZ', search_ID='1078-71-3',
                       particle_size='Soluble', degradability='Slowly',
                       organic=True)
    C6BENZ.copy_models_from(C7BENZ, ['mu'])

    C8BENZ = Component('C8BENZ', search_ID='2189-60-8',
                       particle_size='Soluble', degradability='Slowly',
                       organic=True)

    C10H16O4 = Component('C10H16O4', search_ID='94-60-0',
                         particle_size='Soluble', degradability='Slowly',
                         organic=True)

    C15H32 = Component('C15H32', search_ID='629-62-9', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    C16H34 = Component('C16H34', search_ID='544-76-3', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    C17H36 = Component('C17H36', search_ID='629-78-7', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    C18H38 = Component('C18H38', search_ID='593-45-3', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    C19H40 = Component('C19H40', search_ID='629-92-5', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    C20H42 = Component('C20H42', search_ID='638-36-8', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    C21H44 = Component('C21H44', search_ID='629-94-7', particle_size='Soluble',
                       degradability='Slowly', organic=True)

    TRICOSANE = Component('TRICOSANE', search_ID='638-67-5',
                          particle_size='Soluble', degradability='Slowly',
                          organic=True)

    C24H38O4 = Component('C24H38O4', search_ID='27554-26-3',
                         particle_size='Soluble', degradability='Slowly',
                         organic=True)

    C26H42O4 = Component('C26H42O4', search_ID='27554-26-3',
                         particle_size='Soluble', degradability='Slowly',
                         organic=True)
    # use C24H38O4 to replace since they are similar

    C30H62 = Component('C30H62', search_ID='638-68-6', particle_size='Soluble',
                      degradability='Slowly', organic=True)
    
    Gasoline = Component('Gasoline', search_ID='544-76-3', phase='l', particle_size='Soluble',
                         degradability='Slowly', organic=True)
    # Gasoline copies C16H34, do not need to be precise, since Gasoline is just used to calculte fuel production amount

    Diesel = Component('Diesel', search_ID='629-94-7', phase='l', particle_size='Soluble',
                       degradability='Slowly', organic=True)
    # Diesel copies C21H44, do not need to be precise, since Diesel is just used to calculte fuel production amount
    
    # values for all catalysts are made up since we just need to calculate cost for catalysts
    CHG_catalyst = Component('CHG_catalyst', phase='s', particle_size='Particulate',
                             degradability='Undegradable', organic=False)
    add_V_from_rho(CHG_catalyst, 1500)
    CHG_catalyst.copy_models_from(Chemical('CaCO3'),('Cn',))
    
    HT_catalyst = Component('HT_catalyst', phase='s', particle_size='Particulate',
                             degradability='Undegradable', organic=False)
    add_V_from_rho(HT_catalyst, 1500)
    HT_catalyst.copy_models_from(Chemical('CaCO3'),('Cn',))
    
    HC_catalyst = Component('HC_catalyst', phase='s', particle_size='Particulate',
                             degradability='Undegradable', organic=False)
    add_V_from_rho(HC_catalyst, 1500)
    HC_catalyst.copy_models_from(Chemical('CaCO3'),('Cn',))
    
    # values for Membrane are made up since we just need to calculate cost for membrane
    Membrane = Component('Membrane', phase='s', particle_size='Particulate',
                             degradability='Undegradable', organic=False)
    add_V_from_rho(Membrane, 1500)
    Membrane.copy_models_from(Chemical('CaCO3'),('Cn',))
    
    DAP = Component('DAP', search_ID='7783-28-0', phase='s', particle_size='Soluble',
                    degradability='Undegradable', organic=False)
    
    MEA = Component('MEA', search_ID='141-43-5', phase='l', particle_size='Soluble',
                    degradability='Slowly', organic=True)
    
    Urea = Component('Urea', search_ID='57-13-6', phase='s', particle_size='Soluble',
                     degradability='Slowly', organic=True)
    
    HNO3 = Component('HNO3', search_ID='7697-37-2', phase='l', particle_size='Soluble',
                     degradability='Undegradable', organic=False)
    
    UAN = Component('UAN', formula='CH6N4O4', phase='s', particle_size='Soluble',
                    degradability='Slowly', organic=True)
    UAN.copy_models_from(Chemical('Urea'),('Cn',))
    # https://www.cfindustries.com/globalassets/cf-industries/media/documents/\
    # product-specification-sheets/uan---north-america/urea-ammonium-nitrate-solution-\
    # 28-30-32.pdf (accessed 2025-02-09)
    add_V_from_rho(UAN, 1300)
    
    # use the CAS number of acrylamide instead of polyacrylamide (not in database)
    PAM = Component('PAM', search_ID='79-06-1', phase='s', particle_size='Soluble',
                    degradability='Slowly', organic=True)
    
    # assume Sawdust have the same composition as Sludge_carbo (carbohydrate)
    Sawdust = Component('Sawdust', phase='s',
                        particle_size='Particulate',
                        formula='C56H95O24N9P',
                        degradability='Undegradable',
                        organic=True)
    add_V_from_rho(Sawdust, 1400)
    Sawdust.HHV = 22.0*10**6*Sawdust.MW/1000
    Sawdust.Cn.add_model(1.25*10**3*Sawdust.MW/1000)
    Sawdust.mu.add_model(6000)
    # made up value, so that HTL.ins[0].nu = 0.03 m2/s ~30000 cSt
    # (NREL 2013 appendix B)
    
    # assume Compost have the same composition as Sludge_carbo (carbohydrate)
    Compost = Component('Compost', phase='s',
                        particle_size='Particulate',
                        formula='C56H95O24N9P',
                        degradability='Undegradable',
                        organic=True)
    add_V_from_rho(Compost, 1400)
    Compost.HHV = 22.0*10**6*Compost.MW/1000
    Compost.Cn.add_model(1.25*10**3*Compost.MW/1000)
    Compost.mu.add_model(6000)
    # made up value, so that HTL.ins[0].nu = 0.03 m2/s ~30000 cSt
    # (NREL 2013 appendix B)
    
    HCl = Component('HCl', phase='l', particle_size='Soluble',
                    degradability='Undegradable', organic=False)
    
    # when calculating price using HHV, use a separate HHV value but not directly from this Componenet
    # assume Biooil is the same as Biocrude, use palmitamide to represent both
    Biooil = Component('Biooil', search_ID='629-54-9',
                       particle_size='Soluble', degradability='Slowly',
                       organic=True)
    
    # use anthracene to represent it
    Tar = Component('Tar', search_ID='120-12-7',
                    particle_size='Soluble', degradability='Slowly',
                    organic=True)
    
    Biochar = Component('Biochar', phase='s', particle_size='Particulate',
                        degradability='Undegradable', organic=False)
    # assume a bulk density of 400 kg/m3; the density does not affect transportation (if any) which is based on weight
    add_V_from_rho(Biochar, 400)
    Biochar.copy_models_from(Chemical('CaCO3'),('Cn',))
    
    cmps = Components([Sludge_lipid, Sludge_protein, Sludge_carbo, Sludge_ash,
                       Struvite, Hydrochar, Residue, Biocrude, HTLaqueous, H2O,
                       C, N, P, O2, N2, N2O, CH4, C2H6, C3H8, CO2, CO, H2, NH3,
                       H2SO4, H3PO4, CaO, CaOH2, CaCO3, MgCl2, MgO, NaOH,
                       NH42SO4, NH4Cl, C4H10, TWOMBUTAN, NPENTAN, TWOMPENTA,
                       CYCHEX, HEXANE, TWOMHEXAN, HEPTANE, CC6METH, PIPERDIN,
                       TOLUENE, THREEMHEPTA, OCTANE, ETHCYC6, ETHYLBEN, OXYLENE,
                       C9H20, PROCYC6, C3BENZ, FOURMONAN, C10H22, C4BENZ,
                       C11H24, C10H12, C12H26, C13H28, C14H30, OTTFNA, C6BENZ,
                       OTTFSN, C7BENZ, C8BENZ, C10H16O4, C15H32, C16H34, C17H36,
                       C18H38, C19H40, C20H42, C21H44, TRICOSANE, C24H38O4,
                       C26H42O4, C30H62, Gasoline, Diesel, CHG_catalyst,
                       HT_catalyst, HC_catalyst, Membrane, DAP, MEA, Urea, HNO3,
                       UAN, PAM, Sawdust, Compost, HCl, Biooil, Tar, Biochar])
    
    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_alias('H2O','Water')
    if set_thermo: qs_set_thermo(cmps)

    return cmps