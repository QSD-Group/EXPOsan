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
    
    
    Sludge_lipid = Component('Sludge_lipid', phase='s',particle_size='Particulate',
                       formula='C56H95O24N9P',degradability='Undegradable',organic=False)
    add_V_from_rho(Sludge_lipid,721)
    #https://www.aqua-calc.com/page/density-table/substance/sewage-coma-and-blank-sludge (accessed 2022-9-30)
    Sludge_lipid.HHV = 22.0*10**6*Sludge_lipid.MW/1000  #Li et al., 2018
    Sludge_lipid.Cn.add_model(1.25*10**3*Sludge_lipid.MW/1000) # Leow et al., 2015
    Sludge_lipid.mu.add_model(0.006) #https://www.researchgate.net/figure/Apparent
    # -viscosity-of-activated-sludge-in-various-MLSS_fig2_254226541 (accessed 10-7-2022, MLSS = 10 g/L)
    
    
    Sludge_protein = Component('Sludge_protein', phase='s',particle_size='Particulate',
                       formula='C56H95O24N9P',degradability='Undegradable',organic=False)
    add_V_from_rho(Sludge_protein,721)
    Sludge_protein.HHV = 22.0*10**6*Sludge_protein.MW/1000
    Sludge_protein.Cn.add_model(1.25*10**3*Sludge_protein.MW/1000)
    Sludge_protein.mu.add_model(0.006)
    
    
    Sludge_carbo = Component('Sludge_carbo', phase='s',particle_size='Particulate',
                       formula='C56H95O24N9P',degradability='Undegradable',organic=False)
    add_V_from_rho(Sludge_carbo,721)
    Sludge_carbo.HHV = 22.0*10**6*Sludge_carbo.MW/1000
    Sludge_carbo.Cn.add_model(1.25*10**3*Sludge_carbo.MW/1000)
    Sludge_carbo.mu.add_model(0.006)
    
    Sludge_ash = Component('Sludge_ash', phase='s',particle_size='Particulate',
                       formula='C56H95O24N9P',degradability='Undegradable',organic=False)
    add_V_from_rho(Sludge_ash,721)
    Sludge_ash.HHV = 22.0*10**6*Sludge_ash.MW/1000
    Sludge_ash.Cn.add_model(1.25*10**3*Sludge_ash.MW/1000)
    Sludge_ash.mu.add_model(0.006)


    H2O = Component('H2O',phase='l',particle_size='Soluble',
                    degradability='Undegradable',organic=False)
    
    Biochar = Component('Biochar',phase='s',particle_size='Particulate',
                        degradability='Undegradable',organic=False)
    add_V_from_rho(Biochar, 1500)  #Assume 1500kg/m3
    Biochar.copy_models_from(Chemical('CaCO3'),('Cn',))
    
    Biocrude = Component('Biocrude',phase='l',particle_size='Soluble',
                          formula='C19H28O2.7N',degradability='Undegradable',organic=False)
    Biocrude.HHV = 34.9*10**6*Biocrude.MW/1000  #Li et al., 2018
    add_V_from_rho(Biocrude, 980)  #SS et al., PNNL 2021
    Biocrude.copy_models_from(Chemical('palmitamide'),('Cn',))  #Jones et al., 2014
    
    HTLaqueous = Component('HTLaqueous',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False)
    add_V_from_rho(HTLaqueous, 1000)
    HTLaqueous.copy_models_from(Chemical('H2O'),('Cn',))  #HTLaqueous referd to TDS in HTL aqueous phase
    
    HTaqueous = Component('HTaqueous',phase='l',particle_size='Soluble',
                        degradability='Undegradable',organic=False)
    add_V_from_rho(HTaqueous, 1000)
    HTaqueous.copy_models_from(Chemical('H2O'),('Cn',))  #HTaqueous referd to HT aqueous waste
    
    Biooil = Component('Biooil',phase='l',particle_size='Soluble',
                        formula='C100H165O1.5N',degradability='Undegradable',organic=False)
    Biooil.HHV = 45.4*10**6*Biooil.MW/1000  #Li et al., 2018
    add_V_from_rho(Biooil, 794)  #SS et al., PNNL 2016
    Biooil.copy_models_from(Chemical('hexadecane'),('Cn',))  #Jones et al., 2014

    Struvite = Component('Struvite',search_ID='MagnesiumAmmoniumPhosphate',
                         formula='NH4MgPO4·H12O6',phase='s',
                         particle_size='Particulate', 
                         degradability='Undegradable',
                         organic=False)
    add_V_from_rho(Struvite, 1710)
    #http://webmineral.com/data/Struvite.shtml#.YzYvqOzMIiM (accessed 2022-9-30)
    
    Residual = Component('Residual',phase='s',particle_size='Particulate',
                        degradability='Undegradable',organic=False)
    add_V_from_rho(Residual, 1500)  #Assume 1500kg/m3
    Residual.copy_models_from(Chemical('CaCO3'),('Cn',))  #CaCO3?
    
    C = Component('C',phase='l',particle_size='Soluble',degradability='Undegradable',
                  organic=False)
    
    N = Component('N',phase='l',particle_size='Soluble',degradability='Undegradable',
                  organic=False)
    
    P = Component('P',phase='l',particle_size='Soluble',degradability='Undegradable',
                  organic=False)
    
    for i in(C,N,P):
        i.default()
        i.copy_models_from(H2O,('sigma','epsilon','kappa','Cn','mu'))
        add_V_from_rho(i, 1000)
    
    
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
         
    NH3 = Component('NH3',phase='g',particle_size='Dissolved gas',
                    degradability='Undegradable',organic=False)
    
    #heating agent
    p-Terphenyl = Component('Terphenyl',CAS='92-94-4')

    cmps = Components([Sludge_lipid,Sludge_protein,Sludge_carbo,Sludge_ash,
                       Struvite,Biochar,Residual,
                       Biocrude,Biooil,
                       HTLaqueous,HTaqueous,C,N,P,
                       CH4,C2H6,C3H8,C5H12,CO2,CO,H2,
                       H2SO4,H3PO4,MgCl2,NaOH,NH42SO4,H2O,NH3])
    
    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()
    cmps.set_alias('H2O', 'Water')
    if set_thermo: qs_set_thermo(cmps)

    return cmps