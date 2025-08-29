# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from qsdsan import WasteStream, sanunits as su
from qsdsan.utils import AttrGetter, auom, load_data, ospath, get_SRT
from . import data_path
from ._units import SelectiveRecovery
from exposan.werf.utils import plantwide_aeration_demand, plantwide_aeration_energy
from chaospy import distributions as shape


__all__ = (
    'add_performance_metrics',
    'add_NH4_recovery_metric',
    'add_OPEX_metrics',
    'add_downstream_uncertainty'
    )

#%%
def add_performance_metrics(model, effluent_quality=True, SRT=True, biogas=True,
                            sludge=True, aeration=False, aeration_energy=False):
    
    metric = model.metric
    sys = model.system
    s = sys.flowsheet.stream
    cmps = s.RWW.components
    
    if effluent_quality:
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
    
    if SRT:
        @metric(name='SRT', units='d', element='SRT')
        def get_srt():
            return get_SRT(
                sys, ('X_H', 'X_PAO', 'X_AUT'), wastage=[s.WAS],
                active_unit_IDs=('A1', 'A2', 'A3', 'A4', 'O5', 'O6', 'ASR', 'MBR')
                )
    
    if biogas:
        @metric(name='CH4 production', units='kg/hr', element='Biogas')
        def get_CH4_production():
            if 'biogas' in s: return s.biogas.imass['S_ch4'] * s.biogas.components.S_ch4.i_mass
            else: return np.nan
        
        @metric(name='CH4 content', units='%', element='Biogas')
        def get_CH4_content():
            if 'biogas' in s: return s.biogas.imol['S_ch4']/sum(s.biogas.mol) * 100
            else: return np.nan
    
    if sludge:
        @metric(name='sludge production', units='tonne/d', element='Sludge')
        def get_sludge_production():
            return sum(s.cake.mass * cmps.i_mass) * 24e-3

    if aeration: add_aeration_metrics(model, aeration_energy)

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

    
x_aer = 1
x_pump = 1

mixing_power = 3e-3     # 3W/m3 = 3e-3 kW/m3 mixing power requirement
stablize_power = 4.85   # 4.4 kW/wet ton, Tarallo et al. 2015 --> makes more sense to be kWh/wet ton --> 4.85 kWh/wet tonne

acetic_acid_price = 1.8 # 1.8 $/kg 100% acetic acid, GPS-X default
lime_price = 0.124      # assume 90% purity, 124 USD/tonne https://www.imarcgroup.com/quicklime-pricing-report
disposal_price = 68.5   # "land application price per wet ton at WRRF gate", based on National Biosolids Data Project report; GPS-X default is 80 USD/tonne

def add_OPEX_metrics(model):
    
    metric = model.metric
    sys = model.system
    s = sys.flowsheet.stream
    u = sys.flowsheet.unit
    
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

    blower_energy = {}
    @metric(name='liquid aeration energy', units='kW', element='Aeration')
    def get_liq_aer_energy():
        blower_energy.clear()
        blower_energy.update(plantwide_aeration_energy(sys, _cached_aer))
        if 'AED' in blower_energy: not_liq = blower_energy['AED']
        else: not_liq = 0.
        return sum(blower_energy.values()) - not_liq

    @metric(name='total aeration energy', units='kW', element='Aeration')
    def get_aer_energy():
        return sum(blower_energy.values()) * x_aer
    
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
            if 'MBR' in u: h = 10
            else: h = 25
            pump.dP_design = get_dP(h, 0)
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
    
    mixing_energy = {}
    @metric(name='ASR mixing energy', units='kW', element='Mixing')
    def get_asr_mixing_power():
        power = 0
        for unit in u:
            if unit.ID[0] == 'A' and unit.ID[1].isdigit():
                if unit.aeration is None:
                    power += unit.V_max * mixing_power
        if 'ASR' in u:
            V = sum(u.ASR.V_tanks[i] for i in range(u.ASR.N_tanks_in_series) if u.ASR.DO_setpoints[i] == 0)
            power += V * mixing_power
        mixing_energy['ASR'] = power
        return power
    
    @metric(name='AD mixing energy', units='kW', element='Mixing')
    def get_ad_mixing_power():
        if 'AD' in u: power = u.AD.V_liq * mixing_power
        else: power = 0
        mixing_energy['AD'] = power
        return power
    
    opex = {}

    @metric(name='aeration energy cost', units='USD/d', element='OPEX')
    def get_aer_cost():
        opex['aeration'] = c = sum(blower_energy.values()) * x_aer * 24 * sys.power_utility.price
        return c
   
    @metric(name='pumping energy cost', units='USD/d', element='OPEX')
    def get_pump_cost():
        opex['pumping'] = c =  sum(pump_energy.values()) * x_pump * 24 * sys.power_utility.price
        return c
    
    @metric(name='mixing energy cost', units='USD/d', element='OPEX')
    def get_mixing_cost():
        opex['mixing'] = c =  sum(mixing_energy.values()) * 24 * sys.power_utility.price
        return c
    
    @metric(name='external carbon cost', units='USD/d', element='OPEX')
    def get_carbon_cost():
        if 'carbon' in s:
            s.carbon.price = acetic_acid_price * s.carbon.components.S_A.i_mass
            opex['carbon'] = c = s.carbon.cost * 24
        else: opex['carbon'] = c = 0
        return c
    
    @metric(name='coagulant cost', units='USD/d', element='OPEX')
    def get_coagulant_cost():
        if 'MD' in u:
            opex['coagulant'] = c = u.MD.add_OPEX['Coagulant'] * 24
        else: opex['coagulant'] = c = 0
        return c    

    @metric(name='lime cost', units='USD/d', element='OPEX')
    def get_lime_cost():
        if ('AD' not in u) and ('AED' not in u):    # class B lime stabilization in solid trains 3
            cmps = s.cake.components
            tss = s.cake.get_TSS() * 1e-4 # in TS%
            dose = 50 + 4.0*(tss-10)    # in lb CaO per wet ton, linearly correlated w TS%, MOP8 Fig 23.79
            dose *= 0.5 / 0.9 # convert from lb/ton to kg/tonne, 90% purity
            opex['lime'] = c =  sum(s.cake.mass * cmps.i_mass) * 24e-3 * dose * lime_price  
        else: opex['lime'] = c = 0
        return c

    @metric(name='lime stablization energy cost', units='USD/d', element='OPEX')
    def get_stabilization_power_cost():
        if ('AD' not in u) and ('AED' not in u):    # class B lime stabilization in solid trains 3
            cmps = s.cake.components
            opex['stablization_energy'] = c = sum(s.cake.mass * cmps.i_mass) * 24e-3 * stablize_power * sys.power_utility.price
        else: opex['stablization_energy'] = c = 0
        return c

    @metric(name='sludge disposal cost', units='USD/d', element='OPEX')
    def get_sludge_disposal_cost():
        cmps = s.cake.components
        opex['sludge disposal'] = c = sum(s.cake.mass * cmps.i_mass) * 24e-3 * disposal_price
        return c
    
    @metric(name='total OPEX', units='USD/d', element='OPEX')
    def get_opex():
        return sum(opex.values())

# %%

def add_downstream_uncertainty(model):
    param = model.parameter
        
    b = 1
    D = shape.Uniform(0.85, 1.15)   # based on % error between QSDsan and GPS-X results at baseline (unoptimized)
    @param(name='Aeration energy uncertainty', units='-', element='System', 
           baseline=b, distribution=D)
    def set_aer_energy_scalar(x):
        global x_aer
        x_aer = x
    
    b = 1
    D = shape.Uniform(0.85, 1.3)    # based on % error between QSDsan and GPS-X results at baseline (unoptimized)
    @param(name='Pumping energy uncertainty', units='-', element='System', 
           baseline=b, distribution=D)
    def set_pump_energy_scalar(x):
        global x_pump
        x_pump = x

    b = 3e-3
    D = shape.Uniform(b*0.9, b*1.1)
    @param(name='Mixing power', units='kW/m3', element='System', 
           baseline=b, distribution=D)
    def set_mixing_power(p):
        global mixing_power
        mixing_power = p
        
    eia_df = load_data(
        ospath.join(data_path, 'eia_electricity_price.xlsx'), # https://www.eia.gov/electricity/monthly/epm_table_grapher.php?t=table_5_06_b
        header=[0,1], index_col=0
        )
    data = eia_df.loc[:,('All Sectors', 'March 2025 YTD')].to_numpy()

    b = eia_df.loc['U.S. Total', ('All Sectors', 'March 2025 YTD')] 
    D = shape.TruncNormal(data.min()-0.5, data.max()+0.5, mu=b, sigma=data.std())
    @param(name='Electricity price', units='cents/kWh', element='System', 
           baseline=b, distribution=D)
    def set_electricity_price(p):
        qs.PowerUtility.price = p/100
    
    b = 68.5 # land application price per wet ton at WRRF gate, based on National Biosolids Data Project report
    D = shape.Triangle(b*0.5, b, b*1.5)
    @param(name='Sludge disposal price', units='USD/wet tonne', element='System', 
           baseline=b, distribution=D)
    def set_disposal_price(p):
        global disposal_price
        disposal_price = p
