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
    
(2) Snowden-Swan, L. J.; Zhu, Y.; Jones, S. B.; Elliott, D. C.; Schmidt, A. J.; 
    Hallen, R. T.; Billing, J. M.; Hart, T. R.; Fox, S. P.; Maupin, G. D. 
    Hydrothermal Liquefaction and Upgrading of Municipal Wastewater Treatment 
    Plant Sludge: A Preliminary Techno-Economic Analysis; 
    PNNL--25464, 1258731; 2016; https://doi.org/10.2172/1258731.


(3) Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
    Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
    Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
    Process Design and Economics for the Conversion of Algal Biomass to Hydrocarbons: 
    Whole Algae Hydrothermal Liquefaction and Upgrading; PNNL--23227, 1126336; 2014; 
    https://doi.org/10.2172/1126336.

#!!! Add full citation here (you can't find the report with just Jones et al., 2014)
'''

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_components',)

def create_components(set_thermo=True):
    
    #Solids components
    #!!! we usually use upper camel case when naming components/units,
    # e.g., sludge would be Sludge, biooil would be BioOil
    # but we use lower case with underscores for streams
    # e.g., upgraded_biocrude
    Sludge = Component('Sludge', phase='s', formula='C56H95O27N6P',
                       particle_size='Particulate',
                       degradability='Undegradable',organic=False)
    add_V_from_rho(Sludge,721)   
    #https://www.aqua-calc.com/page/density-table/substance/sewage-coma-and-blank-sludge (accessed 2022-9-30)
    #!!! for references, add "accessed"
    Sludge.HHV = 22.0  #Li et al 2018 #!!! double-check unit, 22 is probably MJ/kg, here is J/mol #kg->mol pending
    # Sludge.copy_models_from(Chemical('glucose'),('Cn',)) #glucose will be replaced
    Sludge.add_model(constant) # convert the 1.25 kJ/kg/K from Leow et al., 2015 (https://pubs.rsc.org/en/content/articlelanding/2015/gc/c5gc00574d) to J/mol/K
    
    Struvite = Component('struvite',search_ID='MagnesiumAmmoniumPhosphate',
                         formula='NH4MgPO4·H12O6',phase='s',
                         particle_size='Particulate', 
                         degradability='Undegradable',
                         organic=False)
    add_V_from_rho(Struvite, 1710)
    #http://webmineral.com/data/Struvite.shtml#.YzYvqOzMIiM (accessed 2022-9-30)
    
    Biochar = Component('biochar',phase='s',particle_size='Particulate',
                        degradability='Undegradable',organic=False,i_C=0.399,
                        i_N=0.026,i_P=0.184)
    add_V_from_rho(Biochar, 1500)  #Assume 1500kg/m3
    Biochar.copy_models_from(Chemical('CaCO3'),('Cn',))  #CaCO3?
    
    Residual = Component('residual',phase='s',particle_size='Particulate',
                        degradability='Undegradable',organic=False,i_C=0.484,
                        i_N=0.031,i_P=0.011)
    add_V_from_rho(Residual, 1500)  #Assume 1500kg/m3
    Residual.copy_models_from(Chemical('CaCO3'),('Cn',))  #CaCO3?
    
    HTchar = Component('HTchar',phase='s',particle_size='Particulate',
                        degradability='Undegradable',organic=False,i_C=0.515,
                        i_N=0.150)
    add_V_from_rho(HTchar, 1500)  #Assume 1500kg/m3
    HTchar.copy_models_from(Chemical('CaCO3'),('Cn',))  #CaCO3?
    
    #Oil components
    
    Biocrude = Component('biocrude',phase='l',formula = 'C14H21O1.8N',
                         particle_size='Soluble',
                         degradability='Slowly',organic=True)
    Biocrude.HHV = 34.9  #Li et al 2018 #kg->mol pending
    add_V_from_rho(Biocrude, 980)  #SS et al PNNL 2021
    Biocrude.copy_models_from(Chemical('palmitamide'),('Cn',))  #Jones et al 2014
    
    Biooil = Component('biooil',phase='l',formula='C100H165O1.5N',
                       particle_size='Soluble',
                       degradability='Slowly',organic=True)
    Biooil.HHV = 45.4  #Li et al 2018 #kg->mol pending
    add_V_from_rho(Biooil, 794)  #SS et al PNNL 2016
    Biooil.copy_models_from(Chemical('hexadecane'),('Cn',))  #Jones et al 2014
    
    #Aqueous components
    
    HTLaqueous = Component('HTLaqueous',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.250,
                        i_N=0.103,i_P=0.028,i_COD=0.602)
    add_V_from_rho(HTLaqueous, 1000)
    HTLaqueous.copy_models_from(Chemical('H2O'),('Cn',))
    
    Mixture = Component('mixture',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.213,
                        i_N=0.088,i_P=0.071,i_COD=0.512)
    add_V_from_rho(Mixture, 1000)
    Mixture.copy_models_from(Chemical('H2O'),('Cn',))
    
    CHGfeed = Component('CHGfeed',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.316,
                        i_N=0.130,i_P=0.005,i_COD=0.759)
    add_V_from_rho(CHGfeed, 1000)
    CHGfeed.copy_models_from(Chemical('H2O'),('Cn',))
    
    CHGeffluent = Component('CHGeffluent',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.320,
                        i_N=0.107,i_P=0.007,i_COD=0.015)
    add_V_from_rho(CHGeffluent, 1000)
    CHGeffluent.copy_models_from(Chemical('H2O'),('Cn',))
    
    WW = Component('WW',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.365,
                        i_N=0.012,i_P=0.008,i_COD=0.017)
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


    cmps = Components([Sludge,Struvite,Biochar,HTchar,
                       Biocrude,Biooil,
                       HTLaqueous,CHGfeed,CHGeffluent,WW,
                       CH4,C2H6,C3H8,C5H12,CO2,CO,H2,
                       H2SO4,H3PO4,MgCl2,NaOH,NH42SO4,H2O,NH3])
    
    #!!! This means that if there's no HHV, then you HHV, LHV, and Hf to 0
    # this may or may not be what we wanted, e.g., now cmps.biochar.HHV == 0 #pending
    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_alias('H2O', 'Water')
    if set_thermo: qs_set_thermo(cmps)

    return cmps