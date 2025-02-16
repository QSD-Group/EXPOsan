#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created by Yuyao Huang and Siqi Tang for Enviroloo Clear Toilet system
'''
from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho
from exposan.bwaise import create_components as create_bw_components
import thermosteam as tmo

__all__ = ('create_components', )

## Pre-define the involving components in the toilet system of concern such as Enviroloo Clear
def create_components(set_thermo = True, 
                      #adjust_MW_to_measured_as=False
                      ):
    bw_cmps = create_bw_components(set_thermo=False)

    C = Component('C', phase='l', particle_size='Soluble', degradability='Undegradable', organic=False)
    
    SolubleCH4 = Component('SolubleCH4', search_ID='CH4', phase='l', particle_size='Soluble', degradability='Slowly', organic=True)

    PAC = Component('PAC', search_ID='10124-27-3', phase='s', particle_size='Particulate', degradability='Slowly', organic=False)
    add_V_from_rho(PAC, rho=2800)
                    
    Glucose = Component('Glucose', search_ID='50-99-7', phase='s', particle_size='Particulate', degradability='Readily', organic=False)
    add_V_from_rho(Glucose, rho=1560)
          
    O3 = Component('O3', search_ID='10028-15-6', phase='g', particle_size='Dissolved gas', degradability='Readily', organic=False)
          
          
    NaOH = Component('NaOH', search_ID='1310-73-2', phase='s', particle_size='Particulate', degradability='Readily', organic=False)
    add_V_from_rho(NaOH, rho=2130)
          
    NaClO = Component('NaClO', search_ID='7681-52-9', phase='s', particle_size='Particulate', degradability='Readily', organic=False)
    add_V_from_rho(NaClO, rho=1250)
    
    NO3 = Component('NO3', measured_as = 'N', phase='l', particle_size='Soluble', degradability='Undegradable', organic=False)
    add_V_from_rho(NO3, rho=1.523) # need check
          
    #NH3_l = Component('NH3_l', measured_as = 'N', phase='l', particle_size='Soluble', degradability='Undegradable', organic=False)
          
    #NonNH3 = Component('NonNH3', formula = 'N', measured_as = 'N', phase='l', particle_size='Soluble', degradability='Undegradable', organic=False, description='Non-NH3');
          
    air = Component('air', search_ID='17778-88-0', phase='g', particle_size='Dissolved gas',
                    degradability='Readily', organic=False)
                    # 1.204 kg/m3, cited from https://en.wikipedia.org/wiki/Density_of_air#:~:text=Air%20density%2C%20like%20air%20pressure,International%20Standard%20Atmosphere%20(ISA).
    add_V_from_rho(air, rho=1.204)



          #allowed_values = {
          #'particle_size': ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'),
          #'degradability': ('Readily', 'Slowly', 'Undegradable'),
          #'organic': (True, False)}
          
    cmps = Components((*bw_cmps, C, SolubleCH4, NO3,
                       #H2O, CO2, CH4, N2O, NH3
                       Glucose, O3, air, PAC, NaOH, NaClO))
    
    for i in cmps:
        for attr in ('HHV', 'LHV', 'Hf'):
            if getattr(i, attr) is None: setattr(i, attr, 0)

    cmps.compile()

    cmps.set_alias('H2O', 'Water')
    #cmps.set_alias('CO2', 'Carbon Dioxide')
    cmps.set_alias('CH4', 'Methane')
    if set_thermo: qs_set_thermo(cmps)

    return cmps
