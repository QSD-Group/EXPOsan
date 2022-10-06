#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Jianan Feng <jiananf2@illinois.edu>
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References:

(1) Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
    Quantitative Evaluation of an Integrated System for Valorization of Wastewater 
    Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
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
    Process Design and Economics for the Conversion of Algal Biomass to Hydrocarbons: 
    Whole Algae Hydrothermal Liquefaction and Upgrading; PNNL--23227, 1126336; 2014; 
    https://doi.org/10.2172/1126336.
'''

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_components',)

def create_components(set_thermo=True):
    
    #Solids components
    
    #Sludge is separated into lipid, protein, and carbohydrate for the convenience
    #during sanunits setup.
    C_s = Component('C_s',phase='s',particle_size='Particulate',degradability='Undegradable',
                  organic=False)
    
    N_s = Component('N_s',phase='s',particle_size='Particulate',degradability='Undegradable',
                  organic=False)
    
    P_s = Component('P_s',phase='s',particle_size='Particulate',degradability='Undegradable',
                  organic=False)
    
    Sludge_lipid = Component('Sludge_lipid', phase='s', formula='C56H95O27N6P',
                       particle_size='Particulate',
                       degradability='Undegradable',organic=False)
    add_V_from_rho(Sludge_lipid,721)
    #https://www.aqua-calc.com/page/density-table/substance/sewage-coma-and-blank-sludge (accessed 2022-9-30)
    Sludge_lipid.HHV = 22.0*10**6*Sludge_lipid.MW/1000  #Li et al., 2018
    Sludge_lipid.Cn.add_model(1.25*10**3*Sludge_lipid.MW/1000) # Leow et al., 2015
    
    Sludge_protein = Component('Sludge_protein', phase='s', formula='C56H95O27N6P',
                       particle_size='Particulate',
                       degradability='Undegradable',organic=False)
    add_V_from_rho(Sludge_protein,721)
    Sludge_protein.HHV = 22.0*10**6*Sludge_protein.MW/1000
    Sludge_protein.Cn.add_model(1.25*10**3*Sludge_protein.MW/1000)
    
    Sludge_carbo = Component('Sludge_carbo', phase='s', formula='C56H95O27N6P',
                       particle_size='Particulate',
                       degradability='Undegradable',organic=False)
    add_V_from_rho(Sludge_carbo,721)
    Sludge_carbo.HHV = 22.0*10**6*Sludge_carbo.MW/1000
    Sludge_carbo.Cn.add_model(1.25*10**3*Sludge_carbo.MW/1000)
    
    Struvite = Component('Struvite',search_ID='MagnesiumAmmoniumPhosphate',
                         formula='NH4MgPO4·H12O6',phase='s',
                         particle_size='Particulate', 
                         degradability='Undegradable',
                         organic=False)
    add_V_from_rho(Struvite, 1710)
    #http://webmineral.com/data/Struvite.shtml#.YzYvqOzMIiM (accessed 2022-9-30)
    
    # Biochar = Component('Biochar',phase='s',particle_size='Particulate',
    #                     degradability='Undegradable',organic=False,i_C=0.399,
    #                     i_N=0.026,i_P=0.184)
    # add_V_from_rho(Biochar, 1500)  #Assume 1500kg/m3
    # Biochar.copy_models_from(Chemical('CaCO3'),('Cn',))
    
    Residual = Component('Residual',phase='s',particle_size='Particulate',
                        degradability='Undegradable',organic=False,i_C=0.484,
                        i_N=0.031,i_P=0.011)
    add_V_from_rho(Residual, 1500)  #Assume 1500kg/m3
    Residual.copy_models_from(Chemical('CaCO3'),('Cn',))  #CaCO3?
    
    Char = Component('Char',phase='s',particle_size='Particulate',
                        degradability='Undegradable',organic=False,i_C=0.515,
                        i_N=0.150)
    add_V_from_rho(Char, 1500)  #Assume 1500kg/m3
    Char.copy_models_from(Chemical('CaCO3'),('Cn',))
    
    #Oil components
    
    Biocrude = Component('Biocrude',phase='l',formula = 'C14H21O1.8N',
                         particle_size='Soluble',
                         degradability='Slowly',organic=True)
    Biocrude.HHV = 34.9*10**6*Biocrude.MW/1000  #Li et al., 2018
    add_V_from_rho(Biocrude, 980)  #SS et al., PNNL 2021
    Biocrude.copy_models_from(Chemical('palmitamide'),('Cn',))  #Jones et al., 2014
    
    Biooil = Component('Biooil',phase='l',formula='C100H165O1.5N',
                       particle_size='Soluble',
                       degradability='Slowly',organic=True)
    Biooil.HHV = 45.4*10**6*Biooil.MW/1000  #Li et al., 2018
    add_V_from_rho(Biooil, 794)  #SS et al., PNNL 2016
    Biooil.copy_models_from(Chemical('hexadecane'),('Cn',))  #Jones et al., 2014
    
    #Aqueous components
    C_l = Component('C_l',phase='l',particle_size='Soluble',degradability='Undegradable',
                  organic=False)
    
    N_l = Component('N_l',phase='l',particle_size='Soluble',degradability='Undegradable',
                  organic=False)
    
    P_l = Component('P_l',phase='l',particle_size='Soluble',degradability='Undegradable',
                  organic=False)
    
    HTLaqueous = Component('HTLaqueous',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.250,
                        i_N=0.103,i_P=0.028,i_COD=0.602)
    add_V_from_rho(HTLaqueous, 1000)
    HTLaqueous.copy_models_from(Chemical('H2O'),('Cn',))
    
    Mixture = Component('Mixture',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.213,
                        i_N=0.088,i_P=0.071,i_COD=0.512)
    add_V_from_rho(Mixture, 1000)
    Mixture.copy_models_from(Chemical('H2O'),('Cn',))
    
    CHGfeed = Component('CHGfeed',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.316,
                        i_N=0.085,i_P=0.005,i_COD=0.759)
    add_V_from_rho(CHGfeed, 1000)
    CHGfeed.copy_models_from(Chemical('H2O'),('Cn',))
    
    CHGeffluent = Component('CHGeffluent',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.320,
                        i_N=0.107,i_P=0.007,i_COD=0.015)
    add_V_from_rho(CHGeffluent, 1000)
    CHGeffluent.copy_models_from(Chemical('H2O'),('Cn',))
    
    WW = Component('WW',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.418,
                        i_N=0.014,i_P=0.009,i_COD=0.020)
    add_V_from_rho(WW, 1000)
    WW.copy_models_from(Chemical('H2O'),('Cn',))
    
    #Gas components
    
    CH4 = Component('CH4',phase='g',particle_size='Dissolved gas',
                    degradability='Slowly',organic=True)
    
    C2H6 = Component('C2H6',phase='g',particle_size='Dissolved gas',
                   degradability='Slowly',organic=True)
    
    C3H8 = Component('C3H8',phase='g',particle_size='Dissolved gas',
                   degradability='Slowly',organic=True)
    
    C5H12 = Component('C5H12',phase='g',particle_size='Dissolved gas',
                   degradability='Slowly',organic=True)
    
    CO2 = Component('CO2',phase='g',particle_size='Dissolved gas',
                    degradability='Undegradable',organic=False)
    
    CO = Component('CO',phase='g',particle_size='Dissolved gas',
                   degradability='Undegradable',organic=False)
    
    H2 = Component('H2',phase='g',particle_size='Dissolved gas',
                   degradability='Undegradable',organic=False)
    
    #Other components
    
    H2SO4 = Component('H2SO4',phase='l',particle_size='Soluble',
                      degradability='Undegradable',organic=False)
    
    H3PO4 = Component('H3PO4',phase='l',particle_size='Soluble',
                      degradability='Undegradable',organic=False)
    
    MgCl2 = Component('MgCl2',phase='l',particle_size='Soluble',
                      degradability='Undegradable',organic=False)
    
    NaOH = Component('NaOH',phase='l',particle_size='Soluble',
                     degradability='Undegradable',organic=False)
    
    NH42SO4 = Component('NH42SO4',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False)
    add_V_from_rho(NH42SO4, 1770)
    #https://en.wikipedia.org/wiki/Ammonium_sulfate (accessed 2022-9-30)
    
    H2O = Component('H2O',phase='l',particle_size='Soluble',
                    degradability='Undegradable',organic=False)
         
    NH3 = Component('NH3',phase='g',particle_size='Dissolved gas',
                    degradability='Undegradable',organic=False)

    for i in (C_s, N_s, P_s):
        i.default()
        i.copy_models_from(Chemical('CaCO3'),('sigma','epsilon','kappa','Cn','mu'))
        add_V_from_rho(i, 1500)

    for i in (C_l, N_l, P_l):
        i.default()
        i.copy_models_from(H2O,('sigma','epsilon','kappa','Cn','mu'))
        add_V_from_rho(i, 1000)

    cmps = Components([C_s,N_s,P_s,Sludge_lipid,Sludge_protein,Sludge_carbo,Struvite,Residual,Char,#Biochar,
                       Biocrude,Biooil,
                       C_l,N_l,P_l,HTLaqueous,Mixture,CHGfeed,CHGeffluent,WW,
                       CH4,C2H6,C3H8,C5H12,CO2,CO,H2,
                       H2SO4,H3PO4,MgCl2,NaOH,NH42SO4,H2O,NH3])
    
    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_alias('H2O', 'Water')
    if set_thermo: qs_set_thermo(cmps)

    return cmps