#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Xuan Wang <easonbiubiu99@gmail.com>
    
    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, qsdsan as qs, biosteam as bst, pandas as pd
from chaospy import distributions as shape
from qsdsan.utils import auom, DictAttrSetter
from exposan.phos_rec import get_precipitation_supernatant_price

__all__ = ('create_model',)

_MMBTU_to_MJ = auom('BTU').conversion_factor('MJ')*10**6
_ton_to_kg = auom('ton').conversion_factor('kg')
_tonne_to_kg = auom('tonne').conversion_factor('kg')

MW_CH4 = 16

folder = os.path.dirname(__file__)

def create_model(system=None,
                 perspective='FePO4',
                 exclude_context=False,
                 exclude_IRR=False):
    
    if perspective not in ['FePO4','sludge']:
        raise KeyError('perspective can only be FePO4 or sludge')
    
    flowsheet = system.flowsheet
    unit = flowsheet.unit
    stream = flowsheet.stream
    model = qs.Model(system)
    param = model.parameter
    
    # =============================================================================
    # AcidogenicFermenter
    # =============================================================================
    AF = unit.AF 
    
    df = pd.read_excel(os.path.join(folder, 'data/ratio_reaction_time.xlsx'))
    for _, row in df.iterrows():        
        dist = shape.Uniform(row['baseline']*0.8, row['baseline']*1.2)
        @param(name=row['name'],
               element='AF',
               kind='coupled',
               units='-',
               baseline=row['baseline'],
               distribution=dist)
        def set_param(i, name=row['name']):
            setattr(AF, row['name'], i)     
    
    dist = shape.Uniform(0.68,0.8)
    @param(name='food_waste_moisture',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.74,
           distribution=dist)
    def set_sludge_moisture(i):
        AF.food_waste_moisture=i
    
    dist = shape.Uniform(0.018,0.022)
    @param(name='org_to_gas',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.02,
           distribution=dist)
    def set_org_to_gas(i):
        AF.org_to_gas=i
        
    dist = shape.Uniform(0.54,0.66)
    @param(name='org_to_vfa',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.60,
           distribution=dist)
    def set_org_to_vfa(i):
        AF.org_to_vfa=i
    
    dist = shape.Uniform(0.018,0.022)
    @param(name='org_to_ethanol',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.02,
           distribution=dist)
    def set_org_to_ethanol(i):
        AF.org_to_ethanol=i
    
    dist = shape.Uniform(0.36,0.44)
    @param(name='org_to_residue',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.40,
           distribution=dist)
    def set_org_to_residue(i):
        if AF.org_to_gas + AF.org_to_vfa + AF.org_to_ethanol + i > 1:
            AF.org_to_residue = max(0, 1 - AF.org_to_gas - AF.org_to_vfa - AF.org_to_ethanol)
        else:
            AF.org_to_residue=i
    
    dist = shape.Triangle(0.95,0.98,1)
    @param(name='Fe_reduction',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.98,
           distribution=dist)
    def set_Fe_reduction(i):
        AF.Fe_reduction=i
    
    # =============================================================================
    # SelectivePrecipitation    
    # =============================================================================
    SP = unit.SP
    dist = shape.Uniform(0.0027,0.0033)
    @param(name='acid_dose',
           element='SP',
           kind='coupled',
           units='-',
           baseline=0.003,
           distribution=dist)
    def set_acid_dose(i):
        SP.acid_dose=i
    
    dist = shape.Triangle(1,2,2.2)
    @param(name='oxidant_excess',
           element='SP',
           kind='coupled',
           units='-',
           baseline=2,
           distribution=dist)
    def set_oxidant_excess(i):
        SP.oxidant_excess=i
 
    # =============================================================================
    # HeatDrying
    # =============================================================================
    HD=unit.HD
    
    dist = shape.Uniform(4.5*0.8,4.5*1.2)
    @param(name='unit_heat',
           element='HD',
           kind='coupled',
           units='-',
           baseline=4.5,
           distribution=dist)
    def set_unit_heat(i):
        HD.unit_heat=i
    
    dist = shape.Uniform(214*0.8,214*1.2)
    @param(name='unit_HD_electricity',
           element='HD',
           kind='coupled',
           units='kWh/t',  
           baseline=214,
           distribution=dist)
    def set_HD_unit_electricity(i):
        HD.unit_electricity=i

    # =============================================================================
    # Sintering
    # =============================================================================
    SI=unit.SI
    dist = shape.Uniform(0.8*0.8,0.8*1.2)
    @param(name='combustion_eff',
           element='SI',
           kind='coupled',
           units='-',
           baseline=0.8,
           distribution=dist)
    def set_combustion_eff(i):
        SI.combustion_eff=i
    
    dist = shape.Uniform(30*0.8,30*1.2)
    @param(name='unit_SI_electricity',
           element='SI',
           kind='coupled',
           units='kWh/t',
           baseline=30,
           distribution=dist)
    def set_SI_unit_electricity(i):
        SI.unit_electricity=i       
    
    # =============================================================================
    # TEA
    # =============================================================================
    food_waste = stream.food_waste
    dist = shape.Uniform(75.8*0.8,75.8*1.2)
    @param(name='food_waste_credit',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=75.8,
           distribution=dist)
    def set_food_waste_credit(i):
        food_waste_dry_kg_hr = food_waste.F_mass - food_waste.imass['Water']
        food_waste.price = -i*(food_waste_dry_kg_hr/1000)/food_waste.F_mass
    
    H2SO4 = stream.acid
    dist = shape.Uniform(0.08*0.8,0.08*1.2)
    @param(name='H2SO4_price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=0.08,
           distribution=dist)
    def set_H2SO4_price(i):
        H2SO4.price=i
    
    H2O2 = stream.oxidant
    dist = shape.Uniform(1.46*0.8,1.46*1.2)
    @param(name='H2O2_price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=1.46,
           distribution=dist)
    def set_H2O2_price(i):
        H2O2.price=i
    
    FePO4 = stream.product
    dist = shape.Uniform(2.14*0.8,2.14*1.2)
    @param(name='FePO4_price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=2.14,
           distribution=dist)
    def set_FePO4_price(i):
        FePO4.price=i
    
    HD_NG = stream.heat_drying_natural_gas
    SI_NG = stream.sintering_natural_gas
    dist = shape.Uniform(3.49672/MW_CH4*0.9,3.49672/MW_CH4*1.1)
    @param(name='NG_price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=3.49672/MW_CH4,
           distribution=dist)
    def set_NG_price(i):
        HD_NG.price=SI_NG.price=i
    
    dist = shape.Uniform(0.0782*0.8,0.0782*1.2)
    @param(name='electricity_price',
           element='TEA',
           kind='isolated',
           units='$/kWh',
           baseline=0.0782,
           distribution=dist)
    def set_electricity_price(i):
        qs.PowerUtility.price = i
    
    for heating_agent in bst.HeatUtility.heating_agents:
        if heating_agent.ID == 'low_pressure_steam':
            dist = shape.Uniform(0.2378*0.8,0.2378*1.2)
            @param(name='low_pressure_steam_price',
                   element='TEA',
                   kind='isolated',
                   units='$/kmol',
                   baseline=0.2378,
                   distribution=dist)
            def set_low_pressure_steam_price(i):
                heating_agent.regeneration_price = i
        
        if heating_agent.ID == 'medium_pressure_steam':
            dist = shape.Uniform(0.2756*0.8,0.2756*1.2)
            @param(name='medium_pressure_steam_price',
                   element='TEA',
                   kind='isolated',
                   units='$/kmol',
                   baseline=0.2756,
                   distribution=dist)
            def set_medium_pressure_steam_price(i):
                heating_agent.regeneration_price = i
        
        if heating_agent.ID == 'high_pressure_steam':
            dist = shape.Uniform(0.3171*0.8,0.3171*1.2)
            @param(name='high_pressure_steam_price',
                   element='TEA',
                   kind='isolated',
                   units='$/kmol',
                   baseline=0.3171,
                   distribution=dist)
            def set_high_pressure_steam_price(i):
                heating_agent.regeneration_price = i
        
        if heating_agent.ID == 'natural_gas':
            dist = shape.Uniform(3.49672*0.8,3.49672*1.2)
            @param(name='natural_gas_utility_price',
                   element='TEA',
                   kind='isolated',
                   units='$/kmol',
                   baseline=3.49672,
                   distribution=dist)
            def set_natural_gas_utility_price_price(i):
                heating_agent.regeneration_price = i
    
    if not exclude_context:
        precipitation_supernatant = stream.precipitation_supernatant
        precipitation_supernatant_baseline_price = get_precipitation_supernatant_price(AF.food_sludge_ratio, AF.reaction_time)
        dist = shape.Uniform(precipitation_supernatant_baseline_price*0.8,precipitation_supernatant_baseline_price*1.2)
        @param(name='precipitation_supernatant_price',
               element='TEA',
               kind='isolated',
               units='$/kg',
               baseline=precipitation_supernatant_baseline_price,
               distribution=dist)
        def set_precipitation_supernatant_price(i):
            precipitation_supernatant.price = i
        
        residue = stream.residue
        dist = shape.Uniform(-62.28/_ton_to_kg*1.2,-62.28/_ton_to_kg*0.8)
        @param(name='residue_price',
               element='TEA',
               kind='isolated',
               units='$/kg',
               baseline=-62.28/_ton_to_kg,
               distribution=dist)
        def set_landfill_price(i):
            residue.price=i
    
    if not exclude_IRR:
        if perspective == 'FePO4':
            dist = shape.Triangle(0.05,0.1,0.15)
            @param(name='IRR',
                   element='TEA',
                   kind='isolated',
                   units='-',
                   baseline=0.1,
                   distribution=dist)
            def set_IRR(i):
                tea.IRR=i
        else:
            dist = shape.Triangle(0,0.03,0.05)
            @param(name='IRR',
                   element='TEA',
                   kind='isolated',
                   units='-',
                   baseline=0.03,
                   distribution=dist)
            def set_IRR(i):
                tea.IRR=i
    
    # =========================================================================
    # LCA (unifrom ± 10%)
    # =========================================================================
    for item in qs.ImpactItem.get_all_items().keys():
        if qs.ImpactItem.get_item(item).CFs and qs.ImpactItem.get_item(item).CFs['GlobalWarming'] != 0:
            abs_small = 0.9*qs.ImpactItem.get_item(item).CFs['GlobalWarming']
            abs_large = 1.1*qs.ImpactItem.get_item(item).CFs['GlobalWarming']
            dist = shape.Uniform(min(abs_small,abs_large),max(abs_small,abs_large))
            @param(name=f'{item}_GlobalWarming',
                   setter=DictAttrSetter(qs.ImpactItem.get_item(item), 'CFs', 'GlobalWarming'),
                   element='LCA',
                   kind='isolated',
                   units=qs.ImpactIndicator.get_indicator('GlobalWarming').unit,
                   baseline=qs.ImpactItem.get_item(item).CFs['GlobalWarming'],
                   distribution=dist)
            def set_LCA(i):
                qs.ImpactItem.get_item(item).CFs['GlobalWarming']=i
    
    # =========================================================================
    # metrics
    # =========================================================================
    metric = model.metric
    tea = system.TEA
    lca = system.LCA
    fe_sludge = stream.fe_sludge
    product = stream.product
    SI = unit.SI
    
    @metric(name='Fe_recovery', units='-', element='Performance')
    def get_Fe_recovery():
        return SI.Fe_recovery

    @metric(name='P_recovery', units='-', element='Performance')
    def get_P_recovery():
        return SI.P_recovery
    
    if perspective == 'FePO4':
        @metric(name='FePO4_MSP',units='$/kg',element='TEA')
        def get_FePO4_MSP():
            return tea.solve_price(product)
        
        @metric(name='FePO4_GWP',units='kg_CO2_eq/kg',element='LCA')
        def get_FePO4_GWP():
            return lca.get_total_impacts(operation_only=True,
                                         exclude=(product,),
                                         annual=True)['GlobalWarming']/product.F_mass/system.operating_hours
    else:
        @metric(name='sludge_management_cost',units='$/tonne',element='TEA')
        def get_sludge_management_cost():
            return -tea.solve_price(fe_sludge)*_tonne_to_kg*fe_sludge.F_mass/(fe_sludge.F_mass - fe_sludge.imass['H2O'])
        
        @metric(name='sludge_management_GWP',units='kg_CO2_eq/tonne',element='LCA')
        def get_sludge_management_GWP():
            return lca.get_total_impacts(operation_only=True,
                                         exclude=(fe_sludge,),
                                         annual=True)['GlobalWarming']/(fe_sludge.F_mass - fe_sludge.imass['H2O'])/system.operating_hours*_tonne_to_kg
    
    return model