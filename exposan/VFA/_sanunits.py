#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:27:03 2025

@author: blues
"""

from qsdsan import SanUnit,Construction, WasteStream, Chemical, Component, \
Components, set_thermo as qs_set_thermo
from qsdsan.utils import ospath, data_path, load_data, price_ratio
from numpy import exp

__all__ = ('redox_ED',)

#%%
# creating a fake wastestream to verify the run function of redox-ED, to be 
# replaced by upstream units

Na = Component('Na', phase='l', particle_size='Soluble',
                  degradability='Undegradable', organic=False)
K = Component('K', phase='l', particle_size='Soluble',
                  degradability='Undegradable', organic=False)
H2O = Component('H2O', phase='l', particle_size='Soluble',
                    degradability='Undegradable', organic=False)
Cl = Component('Cl', phase='l', particle_size='Soluble',
                    degradability='Undegrdable', organic=False)
Propionate = Component('Propionate', formula='C3H6O2', phase='l', particle_size='Soluble',
                           degradability='Readily', organic=True)
Butyrate = Component('Butyrate', phase='l', particle_size='Soluble',
                         degradability='Readily', organic=True)
Hexanoate = Component('Hexanoate', phase='l', particle_size='Soluble',
                         degradability='Readily', organic=True)

cmps = Components((Na, K, H2O,Cl,Propionate,Butyrate, Hexanoate))

qs_set_thermo(cmps)

ws1 = WasteStream('ws1', Na=0.0375, K=0, H2O=1662.28, Cl=0, 
                     Propionate=0.0125, Butyrate=0.0125, Hexanoate=0.0125, units='mmol/hr')
# !!! Double check concentrations taken from raw data: C3C4C6 = 25/25/25 mM,
# volumetric flowrate is 0.5 ml/min. Assume water density at 25 C is 0.9982 g/ml,
# and water molar weight is 18.015 g/mol.

ws2 = WasteStream('ws2', Na=0.0375, K=0, H2O=1662.28, Cl=0.0375, 
                     Propionate=0, Butyrate=0, Hexanoate=0, units='mmol/hr')
# concentrations taken from raw data: NaCl is 75 mM/min, volumetric flowrate
# is 0.5 ml/min. Assume water density at 25 C. Assume water density at 25 C is 0.9982 g/ml,
# and water molar weight is 18.015 g/mol.

    
#%%
redox_ED_path = ospath.join(data_path, 'sanunit_data/VFA/_redox_ED.csv')
class redox_ED(SanUnit):
    '''
    Redox-electrodialysis unit to upconcentration volatile fatty acids from waste.
    
    Parameters
    ----------
    ins: Iterable(stream)
        fc_in is feeding channel inflow
        ac_in is accumulating channel inflow
    outs: Iterable(stream)
        fc_out is feeding channel effluent
        ac_out is accumulating channel effluent
    voltage: float
        cell voltage (V)
    op_time: float
        operational time (hr)
    fc_c3: float
        feeding channel propionate inflow concentration (M)
    fc_c4: float
        feeding channel butyrate inflow concentration (M)
    fc_c6: float
        feeding channel hexanoate inflow concentration (M)
    
    Unit conventions:
    -----------------
    concentration: mol/L or M
    volume: L
    voltage: V
    current: A
    time: hr
    
    Concentration calculations:
    ---------------------------
    Assume a balck box model for the ionic flux, continuously-stirred tank reactors
    for the feeding and accumulating channels. Model membrane transfer as 1st order
    reaction and solve a steady state mass balance for concentration.
    
    '''
    _N_ins = 2
    _N_outs = 2
    
    # !!!
    # Linear regression coefficients to extrapolate flux, calculated in excel
    # and imported as instance attributes. Linear regression coefficients
    # should be updated when new experimental data is used for future instances
    # !!! Check units for the mass balance and concentrations!
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 voltage=None, op_time=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        # Defining instance attributes, operation related for redox-ED
        # !!! When do the users input values? because this script is only run 
        # implicitly when a system is created
        
        data = load_data(path=redox_ED_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        
    def _run(self):
        fc_in, ac_in = self.ins
        # !!! need to edit the inputs, fc_in is from the wastestream, ac_in is
        # a supporting electrolyte solution that only inlcudes sodium, potassium,
        # and chloride
        fc_out, ac_out = self.outs
        # Calculate the average flux from cell voltage
        # !!! how to reduce repetitiveness? when flux extrpolation constants
        # are imported in other ways instead of as class attributes in future,
        # need to change ways to retrieve this data as well
        c3_flux = self.c3_slope*self.voltage + self.c3_const
        c4_flux = self.c4_slope*self.voltage + self.c4_const
        c6_flux = self.c6_slope*self.voltage + self.c6_const
        
        area = self.dimension**2
        
        # Calculate accumulating channel effluent concentration
        ac_out.imass['propionate'] = ac_in.imass['propionate'] + \
        c3_flux * area / self.flowrate 
        ac_out.imass['butyrate'] = ac_in.imass['butyrate'] + \
        c4_flux * area / self.flowrate
        ac_out.imass['hexanoate'] = ac_in.imass['hexanoate'] + \
        c6_flux * area / self.flowrate
        
        # Calcualte feeding channel effluent concentration        
        fc_out.imass['propionate'] = fc_in.imass['propionate'] - \
        c3_flux * area / self.flowrate
        fc_out.imass['butyrate'] = fc_in.imass['butyrate'] - \
        c4_flux * area / self.flowrate
        fc_out.imass['hexanoate'] = fc_in.imass['hexanoate'] - \
        c6_flux * area / self.flowrate
        
    def _init_lca(self):
        pass
    
    def _design(self):
        pass
    
    def _cost(self):
        pass

# =============================================================================
#     def _init_lca(self): 
#         self.construction = [Construction('carbon_cloth', linked_unit=self, 
#                                           item='CarbonCloth', 
#                                           quantity_unit='kg'),
#                              Construction('current_collector', linked_unit=self,
#                                           item='CurrentCollector', 
#                                           quantity_unit='kg'),
#                              Construction('silicon_spacer', linked_unit=self,
#                                           item='SiliconSpacer',
#                                           quantity_unit='kg'),
#                              Construction('electrolyte_membrane', linked_unit=self,
#                                           item='ElectrolyteMembrane',
#                                           quantity_unit='m2'),    #is there a common format for units?
#                              Construction('housing', linked_unit=self,
#                                           item='Housing',
#                                           quantity_unit='kg'),
#                              Construction('piping', linked_unit=self,
#                                           item='Piping',
#                                           quantity_unit='m')] 
#             
#     def _design(self):
#         design = self.design_results
#         constr = self.construction
# =============================================================================
        
        
        
        
        