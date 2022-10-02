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
    sludge = Component('sludge', phase='s', formula='C56H95O27N6P',
                       particle_size='Particulate',
                       degradability='Undegradable',organic=False)
    add_V_from_rho(sludge,721)   
    #https://www.aqua-calc.com/page/density-table/substance/sewage-coma-and-blank-sludge
    #!!! for references, add "accessed"
    sludge.HHV = 22.0  #Li et al 2018 #!!! double-check unit, 22 is probably MJ/kg, here is J/mol
    sludge.copy_models_from(Chemical('glucose'),('Cn',)) #glucose will be replaced
    
    struvite = Component('struvite',search_ID='MagnesiumAmmoniumPhosphate',
                         formula='NH4MgPO4·H12O6',phase='s',
                         particle_size='Particulate', 
                         degradability='Undegradable',
                         organic=False)
    add_V_from_rho(struvite, 1710)
    #http://webmineral.com/data/Struvite.shtml#.YzYvqOzMIiM
    
    biochar = Component('biochar',phase='s',particle_size='Particulate',
                        degradability='Undegradable',organic=False,i_C=0.399,
                        i_N=0.026,i_P=0.184)
    add_V_from_rho(biochar, 1500)  #Assume 1500kg/m3
    biochar.copy_models_from(Chemical('CaCO3'),('Cn',))  #CaCO3?
    
    residual = Component('residual',phase='s',particle_size='Particulate',
                        degradability='Undegradable',organic=False,i_C=0.484,
                        i_N=0.031,i_P=0.011)
    add_V_from_rho(residual, 1500)  #Assume 1500kg/m3
    residual.copy_models_from(Chemical('CaCO3'),('Cn',))  #CaCO3?
    
    HTchar = Component('HTchar',phase='s',particle_size='Particulate',
                        degradability='Undegradable',organic=False,i_C=0.515,
                        i_N=0.150)
    add_V_from_rho(HTchar, 1500)  #Assume 1500kg/m3
    HTchar.copy_models_from(Chemical('CaCO3'),('Cn',))  #CaCO3?
    
    #Oil components
    
    biocrude = Component('biocrude',phase='l',formula = 'C14H21O1.8N',
                         particle_size='Soluble',
                         degradability='Slowly',organic=True)
    biocrude.HHV = 34.9  #Li et al 2018
    add_V_from_rho(biocrude, 980)  #SS et al PNNL 2021
    biocrude.copy_models_from(Chemical('palmitamide'),('Cn',))  #Jones et al 2014
    
    biooil = Component('biooil',phase='l',formula='C100H165O1.5N',
                       particle_size='Soluble',
                       degradability='Slowly',organic=True)
    biooil.HHV = 45.4  #Li et al 2018
    add_V_from_rho(biooil, 794)  #SS et al PNNL 2016
    biooil.copy_models_from(Chemical('hexadecane'),('Cn',))  #Jones et al 2014
    
    #Aqueous components
    
    HTLaqueous = Component('HTLaqueous',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.250,
                        i_N=0.103,i_P=0.028,i_COD=0.602)
    add_V_from_rho(HTLaqueous, 1000)
    HTLaqueous.copy_models_from(Chemical('H2O'),('Cn',))
    
    mixture = Component('mixture',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False,i_C=0.213,
                        i_N=0.088,i_P=0.071,i_COD=0.512)
    add_V_from_rho(mixture, 1000)
    mixture.copy_models_from(Chemical('H2O'),('Cn',))
    
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
    #https://en.wikipedia.org/wiki/Ammonium_sulfate
    
    H2O = Component('H2O',phase='l',particle_size='Soluble',
                    degradability='Undegradable',organic=False)
         
    NH3 = Component('NH3',phase='g',particle_size='Dissolved gas',
                    degradability='Undegradable',organic=False)


    cmps = Components([sludge,struvite,biochar,HTchar,
                       biocrude,biooil,
                       HTLaqueous,CHGfeed,CHGeffluent,WW,
                       CH4,C2H6,C3H8,C5H12,CO2,CO,H2,
                       H2SO4,H3PO4,MgCl2,NaOH,NH42SO4,H2O,NH3])
    
    #!!! This means that if there's no HHV, then you HHV, LHV, and Hf to 0
    # this may or may not be what we wanted, e.g., now cmps.biochar.HHV == 0
    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_alias('H2O', 'Water')
    if set_thermo: qs_set_thermo(cmps)

    return cmps