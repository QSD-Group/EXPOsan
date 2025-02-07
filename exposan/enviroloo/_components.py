#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created by Yuyao Huang and Siqi Tang for Enviroloo Clear Toilet system
'''
from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_components', )

## Pre-define the involving components in the toilet system of concern such as Enviroloo Clear
def create_components(set_thermo = True):
          H2O = Component('H2O', search_ID='H2O', particle_size='Soluble',
                     degradability='Undegradable', organic=False)

          CO2 = Component('CO2', search_ID='CO2', phase='g',
                    particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)

          N2O = Component('N2O', search_ID='N2O', phase='g',
                    particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)

          CH4 = Component('CH4', search_ID='CH4', phase='g', particle_size='Dissolved gas',
                     degradability='Readily', organic=True)
          
          PAC = Component('PAC', search_ID='10124-27-3', phase='s', particle_size='Particulate',
                    degradability='Slowly', organic=False)
                    
          Glucose = Component('Glucose', search_ID='50-99-7', phase='l', particle_size='Soluble',
                         degradability='Readily', organic=True)
          
          O3 = Component('O3', search_ID='10028-15-6', phase='g', particle_size='Dissolved gas',
                         degradability='Readily', organic=False)
          
          NH3 = Component('NH3', search_ID='7664-41-7', phase='g', particle_size='Dissolved gas',
                         degradability='Readily', organic=False)
          
          NaOH = Component('NaOH', search_ID='1310-73-2', phase='s', particle_size='Particulate',
                           degradability='Readily', organic=False)
          
          NaClO = Component('NaClO', search_ID='7681-52-9', phase='s', particle_size='Particulate',
                            degradability='Readily', organic=False)
          
          #NH3_l = Component('NH3_l', measured_as = 'N', phase='l', particle_size='Soluble',
                         #degradability='Undegradable', organic=False)
          
          
          #NonNH3 = Component('NonNH3', formula = 'N', measured_as = 'N',
                             #phase='l', particle_size='Soluble',
                             #degradability='Undegradable', organic=False,
                             #description='Non-NH3');
          
          air = Component('air', MW=29, phase='g', particle_size='Dissolved gas', 
                          degradability='Readily', organic=False)
          # 1.204 kg/m3, cited from https://en.wikipedia.org/wiki/Density_of_air#:~:text=Air%20density%2C%20like%20air%20pressure,International%20Standard%20Atmosphere%20(ISA).
          add_V_from_rho(air, 1.204)


          #allowed_values = {
          #'particle_size': ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'),
          #'degradability': ('Readily', 'Slowly', 'Undegradable'),
          #'organic': (True, False)}
          
          cmps = Components((H2O, CO2, CH4, N2O, Glucose, O3, NH3, air, PAC, NaOH, NaClO))
          
          for cmp in cmps:
                    cmp.default()

          cmps.compile()

          cmps.set_alias('H2O', 'Water')
          cmps.set_alias('CO2', 'Carbon Dioxide')
          cmps.set_alias('CH4', 'Methane')

          if set_thermo: qs_set_thermo(cmps)

          return cmps
