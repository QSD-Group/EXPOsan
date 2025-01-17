#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.reclaimer import create_components as create_re_components
from exposan.utils import add_V_from_rho

__all__ = ('create_components',)

def create_components(set_thermo=True):
    #Reuse components in the bwaise and reclaimer for consistency
    #Modify Tissue so that it becomes combustaible
    #HHV 25MJ/lg, https://doi.org/10.1016/j.jenvman.2019.01.005
    re_cmps = create_re_components(set_thermo = True)
    re_cmps = Components(re_cmps)
    re_cmps.Tissue.degradability= 'Slowly' 
    re_cmps.Tissue.HHV= 25*1e3 #kJ/kg
    re_cmps.Tissue.organic= True
    re_cmps.Tissue.f_Vmass_Totmass = 0.9
    
    #create components for sCOD(soluble COD) and xCOD (particulate COD)
    sCOD = Component('sCOD', phase='l', particle_size='Soluble',measured_as='COD',
                     i_C=0.35,
                    degradability='Readily', organic=True, 
                    f_Vmass_Totmass = 1
                    )
    
    xCOD = Component('xCOD', phase='s', particle_size='Particulate',measured_as='COD',
                     i_C=0.35,
                    degradability='Readily', organic=True, 
                    f_Vmass_Totmass = 1
                    )
    
    #https://doi.org/10.1016/j.fuel.2016.07.077
    # HHV is from the above literature
    add_V_from_rho(sCOD, 1560)
    sCOD.HHV = 24.7 * 1e3 *2 #kJ/kg, enlarge 2 times to compensate for OtherSS so that final feces HHV matches 24.7 MJ/kg
    sCOD.copy_models_from(Chemical('Glucose'), ('Cn', 'mu'))
    
    add_V_from_rho(xCOD, 1560)
    xCOD.HHV = 24.7 * 1e3 *2
    xCOD.copy_models_from(Chemical('Glucose'), ('Cn', 'mu'))
    
    NO = Component('NO', phase='g', particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)
    SO2 = Component('SO2', phase='g', particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)
    
    # Wood pellets
    # moisture content = 18.5 %, caloric value (HHV) = 19.16 MJ/kg
    WoodPellet = Component('WoodPellet', phase='s', i_C = 0.474, i_N = 0.0031,
                           particle_size='Particulate',
                           degradability='Undegradable', organic=False, description='Wood pellets as biofuel for combustion')
    # 700 kg/m3 is average from:
    # https://www.biofuelmachines.com/wood-pellets-quality-standards-study.html (accessed 2024-11-20)
    # https://doi.org/10.1016/j.biombioe.2023.106951
    add_V_from_rho(WoodPellet, 700)
    WoodPellet.HHV = 19.16
    WoodPellet.copy_models_from(Chemical('Glucose'), ('Cn', 'mu'))
    
    #remove components that are not used
    not_used_cmps = ('KCl','GAC','Zeolite','LPG','HCl','MgOH2','Struvite')
    
    cmps = Components([cmp for cmp in re_cmps if cmp.ID not in not_used_cmps]+ 
                      [sCOD, xCOD,WoodPellet,NO,SO2]) #new components
    
    #assign value 0 to heat value that does not exist
    for cmp in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(cmp, attr) is None: 
                setattr(cmp, attr, 0)

    cmps.default_compile()
    
    cmps.set_alias('H2O','Water')
    
    if set_thermo: qs_set_thermo(cmps)
    
    return cmps
    
    
    

