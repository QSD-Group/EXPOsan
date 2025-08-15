#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:27:03 2025

@author: blues
"""

from qsdsan import SanUnit,Construction, WasteStream
from numpy import exp

__all__ = ('redox_ED',)

#%%
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
    # and imported as class variables here. Change way to calculate and import 
    # in future
    c3_slope = 0.149374138
    c3_const = 0.307055862
    c4_slope = 0.267336207
    c4_const = 0.155695793
    c6_slope = 0.236910345
    c6_const = 0.088617655
    dimension = 4 # Unit is cm
    thickness = 3 # Unit is cm
    flowrate = 0.5 # ml/min
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 voltage=None, op_time=None, fc_c3=None, fc_c4=None, fc_c6=None, 
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        # Defining instance attributes, operation related for redox-ED
        # When are they assigned, because this script is only run implicitly
        # when a system is created
        self.voltage = voltage
        self.operation_time = op_time
        self.fc_c3 = fc_c3
        self.fc_c4 = fc_c4
        self.fc_c6 = fc_c6
        
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
        c3_flux = redox_ED.c3_slope*self.voltage + redox_ED.c3_const
        c4_flux = redox_ED.c4_slope*self.voltage + redox_ED.c4_const
        c6_flux = redox_ED.c6_slope*self.voltage + redox_ED.c6_const
        
        # Calculate concentration
        fc_out.imass['propionate'] = fc_in.imass['propionate'] + \
        c3_flux * redox_ED.dimension**2 / redox_ED.flowrate - \
        exp(-redox_ED.flowrate / (redox_ED.dimension ** 2 * redox_ED.thickness) 
            *self.operation_time + c3_flux * redox_ED.dimension ** 2 - 1)
        
        fc_out.imass['butyrate'] = fc_in.imass['butyrate'] + \
        c4_flux * redox_ED.dimension**2 / redox_ED.flowrate - \
        exp(-redox_ED.flowrate / (redox_ED.dimension ** 2 * redox_ED.thickness) 
            *self.operation_time + c4_flux * redox_ED.dimension ** 2 - 1)
        
        fc_out.imass['hexanoate'] = fc_in.imass['hexanoate'] + \
        c6_flux * redox_ED.dimension**2 / redox_ED.flowrate - \
        exp(-redox_ED.flowrate / (redox_ED.dimension ** 2 * redox_ED.thickness) 
            *self.operation_time + c6_flux * redox_ED.dimension ** 2 - 1)
        
        
        
        
        