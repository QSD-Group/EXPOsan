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
    'add_OPEX_metrics'
    )

#%%
def add_performance_metrics(model, aeration_energy=False):
    
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

    add_aeration_metrics(model, aeration_energy)

#%%

def add_aeration_metrics(model, energy=True):
    
    metric = model.metric
    sys = model.system
    
    _cached_aer = {}
    @metric(name='liquid aeration flowrate', units='m3/d', element='Aeration')
    def get_liquid_qair():
        _cached_aer.clear()
        _cached_aer.update(plantwide_aeration_demand(sys))
        qair = 0.
        for k,v in _cached_aer.items():
            if k == 'AED': continue
            else: qair += v
        return qair
    
    @metric(name='sludge aeration flowrate', units='m3/d', element='Aeration')
    def get_aed_qair():
        if 'AED' in _cached_aer: qair = _cached_aer['AED']
        else: qair = np.nan
        return qair

    if energy:
        blower_energy = {}
        @metric(name='aeration energy', units='kW', element='Aeration')
        def get_aer_energy():
            blower_energy.clear()
            blower_energy.update(plantwide_aeration_energy(sys, _cached_aer))
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
    has_sr = False
    for unit in u:
        if isinstance(unit, SelectiveRecovery): 
            SR = unit
            has_sr = True
            break
    
    @metric(name='NH4_recovery', units='kg-N/d', element='System')
    def get_NH4_recovery():
        if has_sr:
            recovery = s.Recovered_NH4.imass['S_NH4']
            if "HA_eff" in s and "PC" not in u:
                r = SR.split[nh4_idx]
                recovery += s.RWW.imass['S_NH4']/(1-r)*r
            return recovery * 24
        return np.nan
    
def add_OPEX_metrics(model):
    
    metric = model.metric
    sys = model.system
    s = sys.flowsheet.stream
    u = sys.flowsheet.unit
    
    pumpin = WasteStream('pumpin', T=s.RWW.T, P=s.RWW.P)
    pump = su.Pump('pump', ins=pumpin, ignore_NPSH=False, 
                   init_with='WasteStream', isdynamic=False)
    pump_energy = {}
    
    ft2m = auom('ft').conversion_factor('m')
    def get_dP(hydraulic_head, headloss, unit='ft'):
        H = (hydraulic_head + headloss) * ft2m
        if unit != 'ft': H *= auom(unit).conversion_factor('m')
        return 1000 * 9.81 * H      # in Pa
    
    @metric(name='influent pumping energy', units='kW', element='Pumping')
    def get_inf_pump_power():
        pump.dP_design = get_dP(40, 0)
        pumpin.copy_like(s.RWW)
        pump.simulate()
        pump_energy['inf'] = power = pump.power_utility.consumption
        return power
    
    @metric(name='RAS pumping energy', units='kW', element='Pumping')
    def get_ras_pump_power():
        if 'RAS' in s:
            pump.dP_design = get_dP(25, 0)
            pumpin.copy_like(s.RAS)
            pump.simulate()
            power = pump.power_utility.consumption
        else: power = 0
        pump_energy['RAS'] = power
        return power
    
    @metric(name='WAS pumping energy', units='kW', element='Pumping')
    def get_was_pump_power():
        if 'WAS' in s:
            pump.dP_design = get_dP(45, 0)
            pumpin.copy_like(s.WAS)
            pump.simulate()
            power = pump.power_utility.consumption
        else: power = 0
        pump_energy['WAS'] = power
        return power
    
    @metric(name='internal recycle pumping energy', units='kW', element='Pumping')
    def get_intr_pump_power():
        power = 0
        if 'ASR' in u:
            for i, j, q in u.ASR.internal_recycles:
                y =  u.ASR.state.iloc[i,:].to_numpy()
                pumpin.mass = y[:-1] * q / 24e3
                pumpin.F_vol = q / 24
                pump.dP_design = get_dP(10, 0)
                pump.simulate()
                power += pump.power_utility.consumption
        if 'intr' in s:
            pumpin.copy_like(s.intr)
            pump.dP_design = get_dP(10, 0)
            pump.simulate()
            power += pump.power_utility.consumption
        pump_energy['intr'] = power
        return power
    
    @metric(name='permeate pumping energy', units='kW', element='Pumping')
    def get_permeate_pump_power():
        if 'MBR' in u:
            pumpin.copy_like(u.MBR.outs[0])
            pump.dP_design = get_dP(10, 0)
            pump.simulate()
            power = pump.power_utility.consumption
        else: power = 0
        pump_energy['Permeate'] = power
        return power
    
    @metric(name='GT underflow pumping energy', units='kW', element='Pumping')
    def get_gtuf_pump_power():
        if 'GT' in u:
            pumpin.copy_like(u.GT.outs[1])
            pump.dP_design = get_dP(40, 0)
            pump.simulate()
            power = pump.power_utility.consumption
        else: power = 0
        pump_energy['GT_uf'] = power
        return power

    @metric(name='MT underflow pumping energy', units='kW', element='Pumping')
    def get_mtuf_pump_power():
        if 'MT' in u:
            pumpin.copy_like(u.MT.outs[1])
            pump.dP_design = get_dP(50, 0)
            pump.simulate()
            power = pump.power_utility.consumption
        else: power = 0
        pump_energy['MT_uf'] = power
        return power
        
    @metric(name='AD pumping energy', units='kW', element='Pumping')
    def get_ad_pump_power():
        if 'AD' in u:
            pumpin.copy_like(u.J1.ins[0])
            pump.dP_design = get_dP(30, 0)
            pump.simulate()
            power = pump.power_utility.consumption
        else: power = 0
        pump_energy['AD'] = power
        return power
    
    @metric(name='pumping energy cost', units='USD/d', element='OPEX')
    def get_pump_cost():
        return sum(pump_energy.values()) * 24 * sys.power_utility.price
    
    mixing_energy = {}
    @metric(name='ASR mixing energy', units='kW', element='Mixing')
    def get_asr_mixing_power():
        power = 0
        for unit in u:
            if unit.ID[0] == 'A' and unit.ID[1].isdigit():
                if unit.aeration is None:
                    power += unit.V_max * 3e-3
        if 'ASR' in u:
            V = sum(u.ASR.V_tanks[i] for i in range(u.ASR.N_tanks_in_series) if u.ASR.DO_setpoints[i] == 0)
            power += V * 3e-3 # 3W/m3 mixing power requirement
        mixing_energy['ASR'] = power
        return power
    
    @metric(name='AD mixing energy', units='kW', element='Mixing')
    def get_ad_mixing_power():
        if 'AD' in u: power = u.AD.V_liq * 3e-3
        else: power = 0
        mixing_energy['AD'] = power
        return power
    
    @metric(name='mixing energy cost', units='USD/d', element='OPEX')
    def get_mixing_cost():
        return sum(mixing_energy.values()) * 24 * sys.power_utility.price
    
    @metric(name='chemical cost', units='USD/d', element='OPEX')
    def get_chemical_cost():
        cost = 0
        if 'carbon' in s:
            s.carbon.price = 1.8 * s.carbon.components.S_A.i_mass   # 1.8 $/kg 100% acetic acid, GPS-X default
            cost += s.carbon.cost * 24
        if 'MD' in u:
            cost += u.MD.add_OPEX['Coagulant'] * 24
        if ('AD' not in u) and ('AED' not in u):    # class B lime stabilization in solid trains 3
            cmps = s.cake.components
            tss = s.cake.get_TSS() * 1e-4 # in TS%
            dose = 50 + 4.0*(tss-10)    # in lb CaO per wet ton, linearly correlated w TS%, MOP8 Fig 23.79
            dose *= 0.5 # convert from lb/ton to kg/tonne
            cost += sum(s.cake.mass * cmps.i_mass) * 24e-3 * dose / 0.9 * 0.124 # assume 90% purity, 124 USD/tonne https://www.imarcgroup.com/quicklime-pricing-report 
        return cost
    
    @metric(name='lime stablization energy cost', units='USD/d', element='Misc')
    def get_stabilization_power():
        if ('AD' not in u) and ('AED' not in u):    # class B lime stabilization in solid trains 3
            cmps = s.cake.components
            return sum(s.cake.mass * cmps.i_mass) * 24e-3 * 4.85 * sys.power_utility.price  # 4.4 kW/wet ton, Tarallo et al. 2015 --> makes more sense to be kW/(tonne/d)
        else:
            return 0
    
    @metric(name='sludge disposal cost', units='USD/d', element='OPEX')
    def get_sludge_disposal_cost():
        cmps = s.cake.components
        return sum(s.cake.mass * cmps.i_mass) * 24e-3 * 80 # 80 USD/tonne, GPS-X default