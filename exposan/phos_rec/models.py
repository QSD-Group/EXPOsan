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

import qsdsan as qs, biosteam as bst
from chaospy import distributions as shape
from qsdsan.utils import auom, DictAttrSetter

__all__ = ('create_model',)

_MMBTU_to_MJ = auom('BTU').conversion_factor('MJ') * 10**6
_ton_to_kg = auom('ton').conversion_factor('kg')

MW_CH4 = 16

def create_model(system=None, perspective='FePO4'): #revised version：def create_model(system): flowsheet = system.flowsheet  for the creat_model() to a AttributeError
    
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
    
    # # TODO: remove this; add this add a scenario analysis instead
    # #discrete scenarcios for food_slduge_ratio
    # _choices = [0, 1/3, 2/3, 1, 4/3] #5 discrete points
    # dist = shape.DiscreteUniform(0,len(_choices)-1)
    # @param(name='food_sludge_ratio',
    #        element='AF',
    #        kind='coupled',
    #        units='-',
    #        baseline=3, # TODO: the value of i
    #        distribution=dist)
    # def set_food_sludge_ratio(i):
    #     AF.food_sludge_ratio=_choices[int(i)]  
    
    dist = shape.Uniform(0.02124263*0.8,0.02124263*1.2)
    @param(name='metal_release_ratio_0',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.02124263,
           distribution=dist)
    def set_metal_release_ratio_0(i):
        AF.metal_release_ratio_0=i
    
    dist = shape.Uniform(0.2113664*0.8,0.2113664*1.2)
    @param(name='metal_release_ratio_1_3',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.2113664,
           distribution=dist)
    def set_metal_release_ratio_1_3(i):
        AF.metal_release_ratio_1_3=i
    
    dist = shape.Uniform(0.6067351*0.8,0.6067351*1.2)
    @param(name='metal_release_ratio_2_3',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.6067351,
           distribution=dist)
    def set_metal_release_ratio_2_3(i):
        AF.metal_release_ratio_2_3=i
    
    dist = shape.Uniform(0.8307559*0.8,0.8307559*1.2)
    @param(name='metal_release_ratio_1',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.8307559,
           distribution=dist)
    def set_metal_release_ratio_1(i):
        AF.metal_release_ratio_1=i
    
    dist = shape.Uniform(0.85*0.8,0.85*1.2)
    @param(name='metal_release_ratio_4_3',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.85,
           distribution=dist)
    def set_metal_release_ratio_4_3(i):
        AF.metal_release_ratio_4_3=i
    
    dist = shape.Uniform(0.110693*0.8,0.110693*1.2)
    @param(name='P_release_ratio_0',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.110693,
           distribution=dist)
    def set_P_release_ratio_0(i):
        AF.P_release_ratio_0=i
    
    dist = shape.Uniform(0.4431886*0.8,0.4431886*1.2)
    @param(name='P_release_ratio_1_3',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.4431886,
           distribution=dist)
    def set_P_release_ratio_1_3(i):
        AF.P_release_ratio_1_3=i
    
    dist = shape.Uniform(0.7122673*0.8,0.7122673*1.2)
    @param(name='P_release_ratio_2_3',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.7122673,
           distribution=dist)
    def set_P_release_ratio_2_3(i):
        AF.P_release_ratio_2_3=i
    
    dist = shape.Uniform(0.8230645*0.8,0.8230645*1.2)
    @param(name='P_release_ratio_1',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.8230645,
           distribution=dist)
    def set_P_release_ratio_1(i):
        AF.P_release_ratio_1=i
    
    dist = shape.Uniform(0.83*0.8,0.83*1.2)
    @param(name='P_release_ratio_4_3',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.83,
           distribution=dist)
    def set_P_release_ratio_4_3(i):
        AF.P_release_ratio_4_3=i
    
    dist = shape.Uniform(0.68,0.8)
    @param(name='food_waste_moisture',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.74,
           distribution=dist)
    def set_sludge_moisture(i):
        AF.food_waste_moisture=i  
    
    # TODO: leave HRT in uncertainty analysis for now; but this can be removed if that makes sense
    dist = shape.Uniform(80,120)
    @param(name='AF_HRT',
           element='AF',
           kind='coupled',
           units='hr',
           baseline=100,
           distribution=dist)
    def set_AF_HRT(i):
        # TODO: check if the value actually changes: the HRT value changes in the result spreadsheet, however, it may be safer to write HRT as a parameter in the sanunit
        AF.HRT=i
    
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
    
    # TODO: leave HRT in uncertainty analysis for now; but this can be removed if that makes sense
    dist = shape.Uniform(8,12)
    @param(name='SP_HRT',
           element='SP',
           kind='coupled',
           units='hr',
           baseline=10,
           distribution=dist)
    def set_SP_HRT(i):
        # TODO: check if the value actually changes: the HRT value changes in the result spreadsheet, however, it may be safer to write HRT as a parameter in the sanunit
        SP.HRT=i
    
    # =============================================================================
    # TEA
    # =============================================================================
    H2SO4 = stream.acid
    # TODO: update if needed
    dist = shape.Uniform(0.08*0.9,0.08*1.1)
    @param(name='H2SO4_price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=0.08,
           distribution=dist)
    def set_H2SO4_price(i):
        H2SO4.price=i
    
    H2O2 = stream.oxidant
    # TODO: update if needed
    dist = shape.Uniform(1.46*0.9,1.46*1.1)
    @param(name='H2O2_price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=1.46,
           distribution=dist)
    def set_H2O2_price(i):
        H2O2.price=i
    
    FePO4 = stream.product
    # TODO: update if needed
    dist = shape.Uniform(0.2*0.9,0.2*1.1)
    @param(name='FePO4_price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=0.2,
           distribution=dist)
    def set_FePO4_price(i):
        FePO4.price=i
    
    # heat drying natural gas
    HD_NG = stream.heat_drying_natural_gas
    # sintering natural gas
    SI_NG = stream.sintering_natural_gas
    dist = shape.Uniform(bst.HeatUtility.heating_agents[-1].regeneration_price/MW_CH4*0.9,bst.HeatUtility.heating_agents[-1].regeneration_price/MW_CH4*1.1)
    @param(name='NG_price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=bst.HeatUtility.heating_agents[-1].regeneration_price/MW_CH4,
           distribution=dist)
    def set_NG_price(i):
        HD_NG.price=SI_NG.price=i
    
    landfill = stream.residue
    # TODO: update if needed
    dist = shape.Uniform(-62.28/_ton_to_kg*1.1,-62.28/_ton_to_kg*0.9)
    @param(name='landfill_price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=-62.28/_ton_to_kg, # landfill price: 62.28 $/ton
           distribution=dist)
    def set_landfill_price(i):
        landfill.price=i
    
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
        # TODO: may need update
        @metric(name='FePO4_MSP',units='$/kg',element='TEA')
        def get_FePO4_MSP():
            return tea.solve_price(product)
        
        # TODO: may need update
        @metric(name='FePO4_GWP',units='kg_CO2_eq/kg',element='LCA')
        def get_FePO4_GWP():
            return lca.get_total_impacts(operation_only=True,
                                         exclude=(product,),
                                         annual=True)['GlobalWarming']/product.F_mass/system.operating_hours
    else:
        # TODO: may need update (*1000 is to convert kg to tonne, need to add this as a constant, *100 is just a rough estimate assuming the moisture content of fe_sludge is 99%)
        @metric(name='sludge_management_cost',units='$/tonne',element='TEA')
        def get_sludge_management_cost():
            return -tea.solve_price(fe_sludge)*1000*100
        
        # TODO: may need update (*1000 is to convert kg to tonne, need to add this as a constant, *100 is just a rough estimate assuming the moisture content of fe_sludge is 99%)
        @metric(name='sludge_management_GWP',units='kg_CO2_eq/tonne',element='LCA')
        def get_sludge_management_GWP():
            return lca.get_total_impacts(operation_only=True,
                                         exclude=(fe_sludge,),
                                         annual=True)['GlobalWarming']/fe_sludge.F_mass/system.operating_hours*1000*100
    
    return model