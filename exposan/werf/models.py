# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from qsdsan import WasteStream, sanunits as su
from qsdsan.utils import AttrGetter, auom
from ._units import SelectiveRecovery
from exposan.werf.utils import plantwide_aeration_demand, plantwide_aeration_energy


__all__ = (
    'add_performance_metrics',
    'add_NH4_recovery_metric',
    )

#%%
def add_performance_metrics(model):
    
    metric = model.metric
    sys = model.system
    s = sys.flowsheet.stream
    cmps = s.RWW.components
    
    kwargs = dict(units='mg/L', element='Effluent')
    metric(getter=AttrGetter(s.SE, 'COD'), name='COD', **kwargs)
    metric(getter=AttrGetter(s.SE, 'BOD'), name='BOD', **kwargs)
    
    @metric(name='TSS', **kwargs)
    def get_TSS():
        return s.SE.get_TSS()
    
    metric(getter=AttrGetter(s.SE, 'TN'), name='TN', **kwargs)
   
    @metric(name='NH4_N', **kwargs)
    def get_NH4_N():
        return s.SE.iconc['S_NH4']
    
    metric(getter=AttrGetter(s.SE, 'TP'), name='TP', **kwargs)
    
    @metric(name='ortho_P', **kwargs)
    def get_orthoP():
        return s.SE.iconc['S_PO4']
    
    metric(getter=AttrGetter(s.SE, 'TOC'), name='TOC', **kwargs)
    
    @metric(name='CH4 production', units='kg/hr', element='Biogas')
    def get_CH4_production():
        if 'biogas' in s: return s.biogas.imass['S_ch4'] * s.biogas.components.S_ch4.i_mass
        else: return np.nan
    
    @metric(name='CH4 content', units='%', element='Biogas')
    def get_CH4_content():
        if 'biogas' in s: return s.biogas.imol['S_ch4']/sum(s.biogas.mol) * 100
        else: return np.nan
    
    @metric(name='sludge production', units='tonne/d', element='Sludge')
    def get_sludge_production():
        return sum(s.cake.mass * cmps.i_mass) * 24e-3

    add_aeration_metrics(model, False)

#%%

def add_aeration_metrics(model, energy=True):
    
    metric = model.metric
    sys = model.system
    
    _cached_aer = {}
    @metric(name='liquid aeration flowrate', units='m3/d', element='Aeration')
    def get_liquid_qair():
        _cached_aer.update(plantwide_aeration_demand(sys))
        qair = 0.
        for k,v in _cached_aer.items():
            if k == 'AED': continue
            else: qair += v
        return qair
    
    @metric(name='sludge aeration flowrate', units='m3/d', element='Aeration')
    def get_aed_qair():
        if 'AED' in _cached_aer: qair = _cached_aer
        else: qair = np.nan
        if not energy: _cached_aer.clear()
        return qair

    if energy:
        blower_energy = {}
        @metric(name='aeration energy', units='kW', element='Aeration')
        def get_aer_energy():
            blower_energy.clear()
            blower_energy.update(plantwide_aeration_energy(sys, _cached_aer))
            _cached_aer.clear()
            return sum(blower_energy.values())

        @metric(name='aeration energy cost', units='USD/d', element='OPEX')
        def get_aer_cost():
            return sum(blower_energy.values()) * 24 * sys.power_utility.price 

#%%    
def add_NH4_recovery_metric(model):
    
    metric = model.metric
    sys = model.system
    s = sys.flowsheet.stream
    u = sys.flowsheet.unit
    cmps = s.RWW.components
    nh4_idx = cmps.index('S_NH4')
    for unit in u:
        if isinstance(unit, SelectiveRecovery): 
            SR = unit
            break
    
    @metric(name='NH4_recovery', units='kg-N/d', element='System')
    def get_NH4_recovery():
        recovery = s.Recovered_NH4.imass['S_NH4']
        if "HA_eff" in s and "PC" not in u:
            r = SR.split[nh4_idx]
            recovery += s.RWW.imass['S_NH4']/(1-r)*r
        return recovery * 24        
    
def add_OPEX_metrics(model):
    
    metric = model.metric
    sys = model.system
    s = sys.flowsheet.stream
    u = sys.flowsheet.unit
    
    pumpin = WasteStream('pumpin', T=s.RWW.T, P=s.RWW.P)
    pump = su.Pump('pump', ins=pumpin, ignore_NPSH=False, 
                   init_with='WasteStream', isdynamic=False)
    pump_energy = {}
    
    def get_dP(hydraulic_head, headloss, unit='ft'):
        H = (hydraulic_head + headloss) * auom(unit).conversion_factor('m')
        return 1000 * 9.81 * H      # in Pa
    
    
    @metric(name='influent pumping energy', units='kW', element='Pumping')
    def get_inf_pump_power():
        pump.dP_design = get_dP(50, 1.0)    #!!! update input
        pumpin.copy_like(s.RWW)
        pump.simulate()
        pump_energy['inf'] = power = pump.power_utility.consumption
        return power
        
    @metric(name='RAS pumping energy', units='kW', element='Pumping')
    def get_ras_pump_power():
        pump.dP_design = get_dP(50, 1.0)    #!!! update input
        pumpin.copy_like(s.RAS)
        pump.simulate()
        pump_energy['RAS'] = power = pump.power_utility.consumption
        return power
    
    # internal recirculations
    # reject waters 
    @metric(name='pumping energy cost', units='USD/d', element='OPEX')
    def get_pump_cost():
        return sum(pump_energy.values()) * 24 * sys.power_utility.price
    
    @metric(name='sludge disposal cost', units='USD/d', element='OPEX')
    def get_sludge_disposal_cost():
        return 