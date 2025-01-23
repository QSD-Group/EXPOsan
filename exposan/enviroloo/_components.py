#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created by Yuyao Huang and Siqi Tang for Enviroloo Clear Toilet system
'''
from qsdsan import Chemical, Component, Components, set_thermo as qs_set_thermo

__all__ = ('create_components',)

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
          
          ''' valid CAS but not included in database ????'''
          #PAC = Component('PAC', search_ID='1327-41-9', phase='s', particle_size='Particulate',
                    #degradability='Slowly', organic=False); 
                    
          
          Glucose = Component('Glucose', search_ID='50-99-7', phase='l', particle_size='Soluble',
                         degradability='Readily', organic=True)
          
          O3 = Component('O3', search_ID='10028-15-6', phase='g', particle_size='Dissolved gas',
                         degradability='Readily', organic=False)
          
          NH3 = Component('NH3', search_ID='7664-41-7', phase='g', particle_size='Dissolved gas',
                         degradability='Readily', organic=False)
          
          #NH3_l = Component('NH3_l', measured_as = 'N', phase='l', particle_size='Soluble',
                         #degradability='Undegradable', organic=False)
          
          
          #NonNH3 = Component('NonNH3', formula = 'N', measured_as = 'N',
                             #phase='l', particle_size='Soluble',
                             #degradability='Undegradable', organic=False,
                             #description='Non-NH3');
          
         # air = Component('air', search_ID='air', phase='g', particle_size='Dissolved gas',
                         #degradability='Readily', organic=False)

          #allowed_values = {
          #'particle_size': ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'),
          #'degradability': ('Readily', 'Slowly', 'Undegradable'),
          #'organic': (True, False)}
          
          cmps = Components((H2O, CO2, CH4,  N2O, Glucose, O3, NH3))
          
          for cmp in cmps:
                    cmp.default()

          cmps.compile()

          cmps.set_alias('H2O', 'Water')
          cmps.set_alias('CO2', 'Carbon Dioxide')
          cmps.set_alias('CH4', 'Methane')

          if set_thermo: qs_set_thermo(cmps)

          return cmps
