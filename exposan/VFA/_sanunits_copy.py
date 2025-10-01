#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:27:03 2025

@author: blues
"""

from qsdsan import SanUnit,Construction, WasteStream, Chemical, Component, \
Components, set_thermo as qs_set_thermo
from qsdsan.utils import ospath, data_path, load_data, price_ratio
# from exposan.g2rt._sanunits import G2RTSolidsSeparation
from numpy import exp

__all__ = ('redox_ED',)


#%% Pretreatment: solids separation using centrifugation
centrifuge_path = ospath.join(data_path, 'sanunit_data/VFA/_centrifuge.csv')
class centrifuge(SanUnit):
    '''
    Reference units are G2RTSolidsSeparation developed by Zixuan Wang and the SolidCetrifuge
    developed by Yoel.
    
    '''
# !!! This is the first unit of the system, the in flow would be AD effluent   
    _N_ins = 1
    _N_outs = 2
    # _units = {}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
# !!! Check what F_BM_default = 1 means, deleted scale up, ppl, and estreme arguments from parent class
        data = load_data(path=centrifuge_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
# !!! in the line below, what components is it referring to? Should import from _components.py?        
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])

    def _init_lca(self):
        self.construction = [Construction("stainless_steel", linked_unit=self,
                                          item = "StainlessSteel", 
                                          quantity_unit= "kg"),
                             Construction("polyethylene", linked_unit=self,
                                          item = "Polyethylene",
                                          quantity_unit= "kg"),
                             Construction("electric_motor", linked_unit=self,
                                          item = "ElectricMotor",
                                          quantity_unit= "kg"),
                             Construction("pump", linked_unit=self,
                                          item = "Pump",
                                          quantity_unit= "ea"),
                             ]        
        
    def _run(self):
        AD_eff, = self.ins # index or comma needed when there is only one stream in inlets
        liquid_stream, solid_stream = self.outs
# This following is defining the subset of components that later could be accessed through []
        solubles, solids = self.solubles, self.solids
        
        TL_in = AD_eff.F_mass - AD_eff.imass[solids].sum()
        mc_in = TL_in / AD_eff.F_mass # including both water and solubles in the moisture content
        mc_out = self.moisture_content_out # this moisture data should assume the total mass of water and solubles
# !!! Below the logic check is for moisture content in the solids but the focus here is the liquid stream
        if mc_in < mc_out*0.999:
            mc_out = mc_in
            
        liquid_stream.imass[solids] = AD_eff.imass[solids] * (1- self.TSS_removal)
        TS_out = liquid_stream.imass[solids].sum() # total soilds in the liquid stream
        TL_out = TS_out / (1 - mc_out) * mc_out # total solubles and water mass flowrate in the liquid stream
        liquid_stream.imass[solubles] = AD_eff.imass[solubles] * TL_out / TL_in 
        breakpoint()
        liquid_stream.imass['H2O'] = TL_out - liquid_stream.imass[solubles].sum()
        # the above line assumes that the solubles (including water and soluble chemicals) partition together and by the same ratio
        #liquid_stream.imass['H2O'] = TS_out*(1-mc_out) * mc_out - liquid_stream.imass[solubles].sum()
        
        solid_stream.imass[solids] = AD_eff.imass[solids] * self.TSS_removal
        solid_stream.imass[solubles] = AD_eff.imass[solubles] - liquid_stream.imass[solubles]
        solid_stream.imass['H2O'] = AD_eff.imass['H2O'] - liquid_stream.imass['H2O']
    
    def _design(self):
        pass
    
    def _cost(self):
        pass
        
        
        
    
#%% Separation: redox-ED
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
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
        self.voltage = voltage
        self.area = self.compartment_length * self.compartment_width
        
    def _run(self):
        fc_in, ac_in = self.ins
        # !!! need to edit the inputs, fc_in is from the wastestream, ac_in is
        # a supporting electrolyte solution that only inlcudes sodium, potassium,
        # and chloride
        fc_out, ac_out = self.outs
        # Calculate the average flux from cell voltage
        # !!! how to reduce repetitiveness? when flux extrapolation constants
        # are imported in other ways instead of as class attributes in future,
        # need to change ways to retrieve this data as well
        c3_flux = self.c3_slope*self.voltage + self.c3_const
        c4_flux = self.c4_slope*self.voltage + self.c4_const
        c6_flux = self.c6_slope*self.voltage + self.c6_const
        
        # Calculate accumulating channel effluent concentration
        ac_out.imass['Propionate'] = ac_in.imass['Propionate'] + \
        c3_flux * self.area / self.flowrate 
        ac_out.imass['Butyrate'] = ac_in.imass['Butyrate'] + \
        c4_flux * self.area / self.flowrate
        ac_out.imass['Hexanoate'] = ac_in.imass['Hexanoate'] + \
        c6_flux * self.area / self.flowrate
        
        # Calcualte feeding channel effluent concentration        
        fc_out.imass['Propionate'] = fc_in.imass['Propionate'] - \
        c3_flux * self.area / self.flowrate
        fc_out.imass['Butyrate'] = fc_in.imass['Butyrate'] - \
        c4_flux * self.area / self.flowrate
        fc_out.imass['Hexanoate'] = fc_in.imass['Hexanoate'] - \
        c6_flux * self.area / self.flowrate
        
# !!! need to consider the mass transport of water?
        
# =============================================================================
#     def _init_lca(self):
#         self.construction = [Construction('carbon_cloth', linked_unit=self, 
#                                            item='CarbonCloth', 
#                                            quantity_unit='kg'),
#                              Construction('current_collector', linked_unit=self,
#                                            item='CurrentCollector', 
#                                            quantity_unit='kg'),
#                              Construction('silicone_spacer', linked_unit=self,
#                                            item='SiliconeSpacer',
#                                            quantity_unit='kg'),
#                              Construction('IEM', linked_unit=self,
#                                            item='IEM',
#                                            quantity_unit='m2'),    #is there a common format for units?
#                              Construction('housing', linked_unit=self,
#                                            item='Housing',
#                                            quantity_unit='kg'),
#                              Construction('piping', linked_unit=self,
#                                            item='Piping',
#                                            quantity_unit='m')]
# 
# # !!! change construction items when scaling up system. Now the materials are
# # based on bench top experiment.
#     
#     def _design(self):
#         design = self.design_results
#         constr = self.construction
#         design['CarbonCloth'] = constr[0].quantity = self.area
#         design['CurrentCollector'] = constr[1].quantity = self.area*self.num_current_collector
#         design['SiliconSpacer'] = constr[2].quantity = self.area*self.num_spacer
#         design['IEM'] = constr[3].quantity = self.area*self.num_IEM
#         design['Housing'] = constr[4].quantity # where would I input quantity value?
#         design['Piping'] = constr[5].quantity
#     
#     def _cost(self):
#         pass
# 
# =============================================================================
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
        
        
        
        
        