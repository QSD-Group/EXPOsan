# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import math
import numpy as np, qsdsan as qs
from qsdsan import WasteStream, unit_operations as su
from qsdsan.utils import AttrGetter, auom, load_data, ospath, get_SRT
from exposan.werf import data_path
from exposan.werf._units import SelectiveRecovery
from exposan.werf.utils import plantwide_aeration_demand, plantwide_aeration_energy
from chaospy import distributions as shape


__all__ = (
    'add_performance_metrics',
    'add_NH4_recovery_metric',
    'add_OPEX_metrics',
    'get_AD_heat_demand',
    'get_CHP_outputs',
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
# Heat price as $/kWh_th (useful heat) derived from EIA industrial natural gas price ($/Mcf)
# Sources / conversions (documented here; multiplication happens outside the functions):
# - EIA industrial natural gas price: https://www.eia.gov/dnav/ng/ng_pri_sum_dcu_nus_a.htm  [3.93 $ per Mcf]
# - Typical heat content: 1 Mcf ≈ 1.038 MMBtu (EIA heat content convention)
# - Unit conversion: 1 MMBtu = 293.071 kWh (since 1 kWh = 3412 Btu)
# - Useful heat accounts for boiler efficiency (e.g., 0.80), so:
#   $/kWh_th_useful = ($/Mcf) / (1.038 MMBtu/Mcf) / (293.071 kWh/MMBtu) / (boiler_eff)
heat_price = 0.0160     # USD/kWh_th useful heat, industrial natural-gas-based proxy

CHP_efficiency_e = 0.40
CHP_efficiency_th = 0.45
CHP_fraction_flared = 0.10
CHP_LHV_kWh_per_kg = 13.9

def get_AD_heat_demand(
        m_in_kg_day, Q_wastewater_m3_day, T_in_C,
        V_digestor_total_m3, V_digestor_headspace_m3, *,
        T_d_C=35, cp_J_kgC=4200.0, U_wall=0.68, U_floor=2.85,
        U_roof=1.5, T_earth_wall=0.0, T_earth_floor=5.0,
        T_air=-5.0):
    if V_digestor_headspace_m3 <= 0:
        raise ValueError('V_digestor_headspace_m3 must be > 0.')
    if V_digestor_total_m3 <= V_digestor_headspace_m3:
        raise ValueError('V_digestor_total_m3 must be > V_digestor_headspace_m3.')
    if Q_wastewater_m3_day <= 0:
        raise ValueError('Q_wastewater_m3_day must be > 0.')
    if m_in_kg_day < 0:
        raise ValueError('m_in_kg_day must be >= 0.')

    V_head = V_digestor_headspace_m3
    V_work = V_digestor_total_m3 - V_head
    r = (3.0 * V_head / (2.0 * math.pi)) ** (1.0 / 3.0)
    V_cone = 0.15 * V_work
    h_cone = (3.0 * V_cone) / (math.pi * r**2)
    V_cyl = V_work - V_cone
    h_cyl = V_cyl / (math.pi * r**2)

    Q_sensible = m_in_kg_day * cp_J_kgC * (T_d_C - T_in_C)
    A_wall = 2.0 * math.pi * r * h_cyl
    A_floor = math.pi * r * math.sqrt(r**2 + h_cone**2)
    A_roof = 2.0 * math.pi * r**2

    sec_per_day = 86400
    Q_loss = (
        U_wall * A_wall * (T_d_C - T_earth_wall)
        + U_floor * A_floor * (T_d_C - T_earth_floor)
        + U_roof * A_roof * (T_d_C - T_air)
        ) * sec_per_day
    return (Q_sensible + Q_loss) / 3.6e6 / Q_wastewater_m3_day

def get_CHP_outputs(
        ch4_kg_hr, Q_wastewater_m3_day, *,
        efficiency_e=CHP_efficiency_e, efficiency_th=CHP_efficiency_th,
        fraction_flared=CHP_fraction_flared,
        LHV_kWh_per_kg=CHP_LHV_kWh_per_kg):
    if ch4_kg_hr < 0:
        raise ValueError('ch4_kg_hr must be >= 0.')
    if Q_wastewater_m3_day <= 0:
        raise ValueError('Q_wastewater_m3_day must be > 0.')
    if not 0 <= fraction_flared <= 1:
        raise ValueError('fraction_flared must be between 0 and 1.')
    if not 0 <= efficiency_e <= 1 or not 0 <= efficiency_th <= 1:
        raise ValueError('efficiency_e and efficiency_th must be between 0 and 1.')

    ch4_to_chp_kg_day = (1.0 - fraction_flared) * ch4_kg_hr * 24
    energy_input_kWh_day = ch4_to_chp_kg_day * LHV_kWh_per_kg
    electricity_kWh_day = energy_input_kWh_day * efficiency_e
    heat_kWh_day = energy_input_kWh_day * efficiency_th
    return electricity_kWh_day / Q_wastewater_m3_day, heat_kWh_day / Q_wastewater_m3_day

def _has_AD(sys):
    return 'AD' in sys.flowsheet.unit

def _has_CHP(sys):
    return sys.ID.upper().endswith('E') and 'AD' in sys.flowsheet.unit

def _get_CH4_kg_hr(s):
    return s.biogas.imass['S_ch4'] * s.biogas.components.S_ch4.i_mass

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
    chp = {}

    if _has_AD(sys):
        @metric(name='AD heating demand', units='kW', element='AD')
        def get_ad_heating_demand():
            Q = s.RWW.F_vol * 24
            uAD = u.AD
            heat_kWh_m3 = get_AD_heat_demand(
                m_in_kg_day=uAD.ins[0].F_mass * 24,
                Q_wastewater_m3_day=Q,
                T_in_C=uAD.ins[0].T - 273.15,
                V_digestor_total_m3=uAD.V_liq + uAD.V_gas,
                V_digestor_headspace_m3=uAD.V_gas,
                )
            chp['AD_heat_kWh_m3'] = heat_kWh_m3
            return heat_kWh_m3 * Q / 24

    if _has_CHP(sys):
        @metric(name='CHP electricity recovery', units='kW', element='CHP')
        def get_chp_electricity_recovery():
            Q = s.RWW.F_vol * 24
            electricity_kWh_m3, heat_kWh_m3 = get_CHP_outputs(_get_CH4_kg_hr(s), Q)
            chp['electricity_kWh_m3'] = electricity_kWh_m3
            chp['heat_kWh_m3'] = heat_kWh_m3
            return -electricity_kWh_m3 * Q / 24

        @metric(name='CHP heat recovery', units='kW', element='CHP')
        def get_chp_heat_recovery():
            if 'heat_kWh_m3' in chp:
                heat_kWh_m3 = chp['heat_kWh_m3']
            else:
                Q = s.RWW.F_vol * 24
                _, heat_kWh_m3 = get_CHP_outputs(_get_CH4_kg_hr(s), Q)
                chp['heat_kWh_m3'] = heat_kWh_m3
            return -heat_kWh_m3 * s.RWW.F_vol

        @metric(name='CHP total recovery', units='kW', element='CHP')
        def get_chp_total_recovery():
            if 'electricity_kWh_m3' not in chp or 'heat_kWh_m3' not in chp:
                Q = s.RWW.F_vol * 24
                electricity_kWh_m3, heat_kWh_m3 = get_CHP_outputs(_get_CH4_kg_hr(s), Q)
                chp['electricity_kWh_m3'] = electricity_kWh_m3
                chp['heat_kWh_m3'] = heat_kWh_m3
            return -(chp['electricity_kWh_m3'] + chp['heat_kWh_m3']) * s.RWW.F_vol

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

    if _has_AD(sys):
        @metric(name='AD heating cost', units='USD/d', element='OPEX')
        def get_ad_heating_cost():
            if 'AD_heat_kWh_m3' in chp:
                heat_kWh_m3 = chp['AD_heat_kWh_m3']
            else:
                uAD = u.AD
                Q = s.RWW.F_vol * 24
                heat_kWh_m3 = get_AD_heat_demand(
                    m_in_kg_day=uAD.ins[0].F_mass * 24,
                    Q_wastewater_m3_day=Q,
                    T_in_C=uAD.ins[0].T - 273.15,
                    V_digestor_total_m3=uAD.V_liq + uAD.V_gas,
                    V_digestor_headspace_m3=uAD.V_gas,
                    )
                chp['AD_heat_kWh_m3'] = heat_kWh_m3
            opex['AD heating'] = c = heat_kWh_m3 * s.RWW.F_vol * 24 * heat_price
            return c

    if _has_CHP(sys):
        @metric(name='CHP electricity credit', units='USD/d', element='OPEX')
        def get_chp_electricity_credit():
            if 'electricity_kWh_m3' in chp:
                electricity_kWh_m3 = chp['electricity_kWh_m3']
            else:
                Q = s.RWW.F_vol * 24
                electricity_kWh_m3, heat_kWh_m3 = get_CHP_outputs(_get_CH4_kg_hr(s), Q)
                chp['electricity_kWh_m3'] = electricity_kWh_m3
                chp['heat_kWh_m3'] = heat_kWh_m3
            opex['CHP electricity'] = c = -electricity_kWh_m3 * s.RWW.F_vol * 24 * sys.power_utility.price
            return c

        @metric(name='CHP heat credit', units='USD/d', element='OPEX')
        def get_chp_heat_credit():
            if 'heat_kWh_m3' in chp:
                heat_kWh_m3 = chp['heat_kWh_m3']
            else:
                Q = s.RWW.F_vol * 24
                _, heat_kWh_m3 = get_CHP_outputs(_get_CH4_kg_hr(s), Q)
                chp['heat_kWh_m3'] = heat_kWh_m3
            opex['CHP heat'] = c = -heat_kWh_m3 * s.RWW.F_vol * 24 * heat_price
            return c

        @metric(name='CHP total credit', units='USD/d', element='OPEX')
        def get_chp_total_credit():
            return opex.get('CHP electricity', 0) + opex.get('CHP heat', 0)

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

    b = eia_df.loc['U.S. Total', ('Industrial', 'March 2025 YTD')]
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

    if _has_CHP(model.system):
        b = 0.0160
        D = shape.Triangle(b*0.5, b, b*1.5)
        @param(name='Heat price', units='USD/kWh_th', element='System',
               baseline=b, distribution=D)
        def set_heat_price(p):
            global heat_price
            heat_price = p
