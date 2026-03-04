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

import qsdsan as qs
from chaospy import distributions as shape
from qsdsan.utils import auom, DictAttrSetter

__all__ = ('create_model',)

_MMBTU_to_MJ = auom('BTU').conversion_factor('MJ') * 10**6
_ton_to_kg = auom('ton').conversion_factor('kg')

def create_model(system=None): #revised version：def create_model(system): flowsheet = system.flowsheet  for the creat_model() to a AttributeError
    flowsheet = system.flowsheet
    unit = flowsheet.unit
    stream = flowsheet.stream
    model = qs.Model(system)
    param = model.parameter
    
    # =============================================================================
    # AcidogenicFermenter
    # =============================================================================
    AF = unit.AF
    #discrete scenarcios for food_slduge_ratio
    _choices = [0,1/3,2/3,1,4/3] #5 discrete points
    dist = shape.DiscreteUniform(0,len(_choices)-1)
    @param(name='food_sludge_ratio',
           element='AF',
           kind='coupled',
           units='-',
           baseline=1,
           distribution=dist)
    def set_food_sludge_ratio(i):
        AF.food_sludge_ratio=_choices[int(i)]  
        
    dist = shape.Uniform(0.68,0.8)
    @param(name='food_waste_moisture',
           element='AF',
           kind='coupled',
           units='-',
           baseline=0.74,
           distribution=dist)
    def set_sludge_moisture(i):
        AF.food_waste_moisture=i  
    
    dist = shape.Uniform(80,120)
    @param(name='HRT', # whether to be the name='AF_HRT'
           element='AF',
           kind='coupled',
           units='hr',
           baseline=100,
           distribution=dist)
    def set_AF_HRT(i):
        # TODO: check if the value actually changes
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
    
    dist = shape.Uniform(8,12)
    @param(name='HRT', # whether to be the name='SP_HRT'
           element='SP',
           kind='coupled',
           units='hr',
           baseline=10,
           distribution=dist)
    def set_SP_HRT(i):
        # TODO: check if the value actually changes
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
    
    # TODO: biosteam has a default price
    # price of natural gas is 7.7$/MMBTU, from Everbatt2023, LHV_ng = 50.1 MJ/kg
    # natural gas consumption was calculated using the lower heating value (LHV = 50.1 MJ/kg), assuming that the latent heat of water vapor in the flue gas was not recovered
    # heat drying natural gas
    HD_NG = stream.heat_drying_natural_gas
    # sintering natural gas
    SI_NG = stream.sintering_natural_gas
    # TODO: update if needed
    dist = shape.Uniform(7.7 / _MMBTU_to_MJ * 50.1*0.9,7.7 / _MMBTU_to_MJ * 50.1*1.1)
    @param(name='NG_price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=7.7 / _MMBTU_to_MJ * 50.1,
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
    product = stream.product
    SI = unit.SI
    
    @metric(name='Fe_recovery', units='-', element='Performance')
    def get_Fe_recovery():
        return SI.Fe_recovery

    @metric(name='P_recovery', units='-', element='Performance')
    def get_P_recovery():
        return SI.P_recovery

    @metric(name='C_recovery', units='-', element='Performance')
    def get_C_recovery():
        return SI.C_recovery
    
    # TODO: may need update
    @metric(name='FePO4_minimum_price',units='$/kg',element='TEA')
    def get_FePO4_minimum_price():
        return tea.solve_price(product)
    
    # TODO: may need update
    @metric(name='FePO4_GWP',units='kg_CO2_eq/kg',element='LCA')
    def get_FePO4_GWP():
        return lca.get_total_impacts(operation_only=True,
                                     exclude=(product,),
                                     annual=True)['GlobalWarming']/product.F_mass/system.operating_hours
    
    return model