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

__all__ = ('create_model',)

_MMBTU_to_MJ = auom('BTU').conversion_factor('MJ')*10**6
_ton_to_kg = auom('ton').conversion_factor('kg')
_tonne_to_kg = auom('tonne').conversion_factor('kg')

MW_CH4 = 16

folder = os.path.dirname(__file__)

def create_model(system=None,
                 perspective='FePO4',
                 exclude_disposal=False,
                 exclude_IRR=False,
                 exclude_breakdown=False):
    
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
        name = row['name']
        baseline = row['baseline']
        
        dist = shape.Uniform(baseline*0.8, min(1, baseline*1.2))
        @param(
            name=name,
            element='AF',
            kind='coupled',
            units='-',
            baseline=baseline,
            distribution=dist
        )
        def set_param(i, name=name):
            setattr(AF, name, i)
    
    dist = shape.Uniform(0.74*0.8,0.74*1.2)
    @param(name='food_waste_moisture',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.74,
           distribution=dist)
    def set_sludge_moisture(i):
        AF.food_waste_moisture=i
    
    dist = shape.Uniform(0.03*0.8,0.03*1.2)
    @param(name='org_unconverted',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.03,
           distribution=dist)
    def set_org_unconverted(i):
        AF.org_unconverted=i
    
    dist = shape.Uniform(0.05*0.8,0.05*1.2)
    @param(name='org_to_gas',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.05,
           distribution=dist)
    def set_org_to_gas(i):
        AF.org_to_gas=i
        
    dist = shape.Uniform(0.65*0.8,0.65*1.2)
    @param(name='org_to_vfa',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.65,
           distribution=dist)
    def set_org_to_vfa(i):
        AF.org_to_vfa=i
    
    dist = shape.Uniform(0.02*0.8,0.02*1.2)
    @param(name='org_to_ethanol',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.02,
           distribution=dist)
    def set_org_to_ethanol(i):
        AF.org_to_ethanol=i
    
    dist = shape.Uniform(0.25*0.8,0.25*1.2)
    @param(name='org_to_residue',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.25,
           distribution=dist)
    def set_org_to_residue(i):
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
    dist = shape.Uniform(0.003*0.8,0.003*1.2)
    @param(name='acid_dose',
           element='SP',
           kind='coupled',
           units='-',
           baseline=0.003,
           distribution=dist)
    def set_acid_dose(i):
        SP.acid_dose=i
    
    dist = shape.Uniform(2*0.8,2*1.2)
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
    HD = unit.HD
    
    dist = shape.Uniform(4.5*0.8,4.5*1.2)
    @param(name='unit_heat',
           element='HD',
           kind='coupled',
           units='-',
           baseline=4.5,
           distribution=dist)
    def set_unit_heat(i):
        HD.unit_heat=i
    
    dist = shape.Uniform(0.8*0.8,0.8*1.2)
    @param(name='combustion_eff_HD',
           element='SI',
           kind='coupled',
           units='-',
           baseline=0.8,
           distribution=dist)
    def set_combustion_eff_HD(i):
        SI.combustion_eff_HD=i
    
    dist = shape.Uniform(214*0.8,214*1.2)
    @param(name='unit_electricity_HD',
           element='HD',
           kind='coupled',
           units='kWh/t',  
           baseline=214,
           distribution=dist)
    def set_unit_electricity_HD(i):
        HD.unit_electricity_HD=i

    # =============================================================================
    # Sintering
    # =============================================================================
    SI = unit.SI
    dist = shape.Uniform(0.8*0.8,0.8*1.2)
    @param(name='combustion_eff_SI',
           element='SI',
           kind='coupled',
           units='-',
           baseline=0.8,
           distribution=dist)
    def set_combustion_eff_SI(i):
        SI.combustion_eff_SI=i
    
    dist = shape.Uniform(30*0.8,30*1.2)
    @param(name='unit_electricity_SI',
           element='SI',
           kind='coupled',
           units='kWh/t',
           baseline=30,
           distribution=dist)
    def set_unit_electricity_SI(i):
        SI.unit_electricity_SI=i       
    
    # =============================================================================
    # TEA
    # =============================================================================
    
    residue = stream.residue
    if not exclude_disposal:
        dist = shape.Uniform(-62.28/_ton_to_kg*1.2,-62.28/_ton_to_kg*0.8)
        @param(name='residue_price',
               element='TEA',
               kind='isolated',
               units='$/kg',
               baseline=-62.28/_ton_to_kg,
               distribution=dist)
        def set_landfill_price(i):
            residue.price=i
    
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
    
    HD_NG = stream.heat_drying_natural_gas
    SI_NG = stream.sintering_natural_gas
    dist = shape.Uniform(3.49672/MW_CH4*0.8,3.49672/MW_CH4*1.2)
    @param(name='NG_price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=3.49672/MW_CH4,
           distribution=dist)
    def set_NG_price(i):
        HD_NG.price=SI_NG.price=i
    
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
        if item != 'VFA_credit':
            if qs.ImpactItem.get_item(item).CFs and qs.ImpactItem.get_item(item).CFs['GlobalWarming'] != 0:
                abs_small = qs.ImpactItem.get_item(item).CFs['GlobalWarming']*0.9
                abs_large = qs.ImpactItem.get_item(item).CFs['GlobalWarming']*1.1
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
    
    if not exclude_breakdown:
        # CAPEX metrics (as total installed cost, TIC)
        @metric(name='TIC',units='$',element='TEA')
        def get_TIC():
            return system.installed_equipment_cost
        
        P1 = unit.P1
        FC = unit.FC
        PC = unit.PC
        
        @metric(name='AF_TIC',units='$',element='TEA')
        def get_AF_TIC():
            return AF.installed_cost
        
        @metric(name='SP_TIC',units='$',element='TEA')
        def get_SP_TIC():
            return SP.installed_cost
        
        @metric(name='thermal_treatment_TIC',units='$',element='TEA')
        def get_thermal_treatment_TIC():
            return HD.installed_cost+SI.installed_cost
        
        @metric(name='others_TIC',units='$',element='TEA')
        def get_others_TIC():
            return sum(i.installed_cost for i in [P1, FC, PC])
        
        @metric(name='AOC',units='$/yr',element='TEA')
        def get_AOC():
            return tea.AOC
        
        @metric(name='FOC',units='$/yr',element='TEA')
        def get_FOC():
            return tea.FOC
        
        @metric(name='VOC',units='$/yr',element='TEA')
        def get_VOC():
            return tea.VOC
        
        @metric(name='material_VOC',units='$/yr',element='TEA')
        def get_material_VOC():
            return system.material_cost-(HD_NG.cost+SI_NG.cost)*system.operating_hours
        
        @metric(name='fe_sludge_VOC',units='$/yr',element='TEA')
        def get_fe_sludge_VOC():
            return fe_sludge.cost*system.operating_hours
        
        @metric(name='residue_VOC',units='$/yr',element='TEA')
        def get_residue_VOC():
            return -residue.cost*system.operating_hours
        
        @metric(name='H2SO4_VOC',units='$/yr',element='TEA')
        def get_H2SO4_VOC():
            return H2SO4.cost*system.operating_hours
        
        @metric(name='H2O2_VOC',units='$/yr',element='TEA')
        def get_H2O2_VOC():
            return H2O2.cost*system.operating_hours
        
        @metric(name='product_VOC',units='$/yr',element='TEA')
        def get_product_VOC():
            return -product.cost*system.operating_hours
        
        @metric(name='utility_VOC',units='$/yr',element='TEA')
        def get_utility_VOC():
            return system.utility_cost+(HD_NG.cost+SI_NG.cost)*system.operating_hours
        
        @metric(name='AF_utility_VOC',units='$/yr',element='TEA')
        def get_AF_utility_VOC():
            return AF.utility_cost*system.operating_hours
        
        @metric(name='SP_utility_VOC',units='$/yr',element='TEA')
        def get_SP_utility_VOC():
            return SP.utility_cost*system.operating_hours
        
        @metric(name='thermal_treatment_utility_VOC',units='$/yr',element='TEA')
        def get_thermal_treatment_utility_VOC():
            return HD.utility_cost*system.operating_hours+\
                   HD_NG.cost*system.operating_hours+\
                   SI.utility_cost*system.operating_hours+\
                   SI_NG.cost*system.operating_hours
        
        @metric(name='others_utility_VOC',units='$/yr',element='TEA')
        def get_others_utility_VOC():
            return sum(i.utility_cost for i in [P1, FC, PC])*system.operating_hours

        @metric(name='stream_GWP',units='kg CO2 eq',element='LCA')
        def get_stream_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return lca.get_stream_impacts()['GlobalWarming']-\
                   table_stream['heat_drying_natural_gas']-\
                   table_stream['sintering_natural_gas']
        
        @metric(name='fe_sludge_GWP',units='kg CO2 eq',element='LCA')
        def get_fe_sludge_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return table_stream['fe_sludge']
        
        @metric(name='residue_GWP',units='kg CO2 eq',element='LCA')
        def get_residue_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return table_stream['residue']
        
        @metric(name='H2SO4_GWP',units='kg CO2 eq',element='LCA')
        def get_H2SO4_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return table_stream['acid']
        
        @metric(name='H2O2_GWP',units='kg CO2 eq',element='LCA')
        def get_H2O2_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return table_stream['oxidant']
        
        @metric(name='product_GWP',units='kg CO2 eq',element='LCA')
        def get_product_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return table_stream['product']
        
        @metric(name='utility_GWP',units='kg CO2 eq',element='LCA')
        def get_utility_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return lca.get_other_impacts()['GlobalWarming']+\
                   table_stream['heat_drying_natural_gas']+\
                   table_stream['sintering_natural_gas']
    
        @metric(name='electricity_GWP',units='kg CO2 eq',element='LCA')
        def get_electricity_GWP():
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            return table_other['Electricity']
        
        @metric(name='steam_GWP',units='kg CO2 eq',element='LCA')
        def get_steam_GWP():
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            return table_other['Steam']
        
        @metric(name='natural_gas_GWP',units='kg CO2 eq',element='LCA')
        def get_natural_gas_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            return table_stream['heat_drying_natural_gas'] + table_stream['sintering_natural_gas']
    
        @metric(name='AF_utility_GWP',units='kg CO2 eq',element='LCA')
        def get_AF_utility_GWP():
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            duty_per_year = sum(i.duty for i in AF.heat_utilities if i.duty>0)*system.operating_hours
            return table_other['Steam']/system.get_heating_duty()*duty_per_year+\
                   table_other['Electricity']/(system.get_electricity_consumption()-system.get_electricity_production())*\
                   AF.power_utility.consumption*system.operating_hours
        
        @metric(name='SP_utility_GWP',units='kg CO2 eq',element='LCA')
        def get_SP_utility_GWP():
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            duty_per_year = sum(i.duty for i in SP.heat_utilities if i.duty>0)*system.operating_hours
            return table_other['Steam']/system.get_heating_duty()*duty_per_year+\
                   table_other['Electricity']/(system.get_electricity_consumption()-system.get_electricity_production())*\
                   SP.power_utility.consumption*system.operating_hours
        
        @metric(name='thermal_treatment_utility_GWP',units='kg CO2 eq',element='LCA')
        def get_thermal_treatment_utility_GWP():
            table_stream = lca.get_impact_table('Stream')['GlobalWarming [kg CO2-eq]']
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            duty_per_year = sum(i.duty for j in [HD, SI] for i in j.heat_utilities if i.duty>0)*system.operating_hours
            return table_other['Steam']/system.get_heating_duty()*duty_per_year+\
                   table_other['Electricity']/(system.get_electricity_consumption()-system.get_electricity_production())*\
                   HD.power_utility.consumption*system.operating_hours+\
                   table_stream['heat_drying_natural_gas']+table_stream['sintering_natural_gas']
        
        @metric(name='others_utility_GWP',units='kg CO2 eq',element='LCA')
        def get_others_utility_GWP():
            table_other = lca.get_impact_table('Other')['GlobalWarming [kg CO2-eq]']
            duty_per_year = sum(i.duty for j in [P1, FC, PC] for i in j.heat_utilities if i.duty>0)*system.operating_hours
            return table_other['Steam']/system.get_heating_duty()*duty_per_year+\
                   table_other['Electricity']/(system.get_electricity_consumption()-system.get_electricity_production())*\
                   sum(i.power_utility.consumption for i in [P1, FC, PC])*system.operating_hours
    
    return model