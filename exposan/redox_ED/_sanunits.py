#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:27:03 2025

@author: blues
"""

from qsdsan import SanUnit,Construction, WasteStream
from qsdsan.utils import ospath, data_path, load_data, price_ratio
from numpy import exp

__all__ = ('redox_ED',)

#%%
redox_ED_path = ospath.join(data_path, 'sanunit_data/_redox_ED.csv')
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
    reaction and solve a non steady state mass balance for concentration.
    
    '''
    _N_ins = 2
    _N_outs = 2
    
    # !!!
    # Linear regression coefficients to extrpolate flux, calculated in excel
    # and imported as instance attributes. Linear regression coefficients
    # should be updated when new experimental data is used for future instances
    
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 voltage=None, op_time=None, fc_c3=None, fc_c4=None, fc_c6=None, 
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        # Defining instance attributes, operation related for redox-ED
        # When do the users input values? because this script is only run 
        # implicitly when a system is created
        
        self.voltage = voltage
        self.operation_time = op_time
        self.fc_c3 = fc_c3
        self.fc_c4 = fc_c4
        self.fc_c6 = fc_c6
        
        data = load_data(path=redox_ED_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        
    def _run(self):
        fc_in, ac_in = self.ins
        fc_out, ac_out = self.outs
        fc_in.imass['propionate'] = self.fc_c3
        fc_in.imass['butyrate'] = self.fc_c4
        fc_in.imass['hexanoate'] = self.fc_c6
        ac_in.imass['propionate'] = self.ac_c3
        ac_in.imass['butyrate'] = self.ac_c4
        ac_in.imass['hexanoate'] = self.ac_c6
        # Calculate the average flux from cell voltage
        # !!! how to reduce repetitiveness? when flux extrpolation constants
        # are imported in other ways instead of as class attributes in future,
        # need to change ways to retrieve this data as well
        c3_flux = self.c3_slope*self.voltage + self.c3_const
        c4_flux = self.c4_slope*self.voltage + self.c4_const
        c6_flux = self.c6_slope*self.voltage + self.c6_const
        
        area = self.dimension**2
        volume = self.dimension ** 2 * self.thickness
        
        # Calculate accumulating channel effluent concentration
        ac_out.imass['propionate'] = ac_in.imass['propionate'] + \
        c3_flux * area / self.flowrate - \
        exp(-self.flowrate / volume 
            *self.operation_time + c3_flux * self.area - 1)
        ac_out.imass['butyrate'] = ac_in.imass['butyrate'] + \
        c4_flux * area / self.flowrate - \
        exp(-self.flowrate / volume 
            *self.operation_time + c4_flux * self.area - 1)
        ac_out.imass['hexanoate'] = ac_in.imass['hexanoate'] + \
        c6_flux * area / self.flowrate - \
        exp(-self.flowrate / volume 
            *self.operation_time + c6_flux * self.area - 1)
        
        # Calcualte feeding channel effluent concentration        
        fc_out.imass['propionate'] = fc_in.imass['propionate'] - \
        c3_flux * area / self.flowrate - \
        exp(-self.flowrate / volume 
            *self.operation_time - c3_flux * self.area - 1)
        fc_out.imass['butyrate'] = fc_in.imass['butyrate'] - \
        c4_flux * area / self.flowrate - \
        exp(-self.flowrate / volume 
            *self.operation_time - c4_flux * self.area - 1)
        fc_out.imass['hexanoate'] = fc_in.imass['hexanoate'] - \
        c6_flux * area / self.flowrate - \
        exp(-self.flowrate / volume 
            *self.operation_time - c6_flux * self.area - 1)

    def _init_lca(self): 
        self.construction = [Construction()]
    def _design(self):
        design = self.design_results
        constr = self.construction
        
        
        
        
        