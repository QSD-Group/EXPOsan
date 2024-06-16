#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Jianan Feng <jiananf2@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs
from exposan.co2_sorbent import create_system_A, create_system_B, create_system_C
from chaospy import distributions as shape
from qsdsan.utils import DictAttrSetter

__all__ = (
    'create_model_A', # ALF production using Al(OH)3
    'create_model_B', # ALF production using Bauxite
    'create_model_C' # CO2 capture and utilization
    )

# =============================================================================
# ALF production: Al(OH)3 + HCOOH
# =============================================================================
def create_model_A(system=None):pass
    
    # sys = create_system_A() if not system else system
    # flowsheet = sys.flowsheet
    # unit = flowsheet.unit
    # stream = flowsheet.stream
    # model = qs.Model(sys)
    # param = model.parameter
    
    # # ALF production
    # ALF_production = unit.ALF_production
    
    # dist = shape.
    # @param(name='ALF_production_T',
    #        element=ALF_production,
    #        kind='coupled',
    #        units='K',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_production_T(i):
    #     ALF_production.T=i
    
    # dist = shape.
    # @param(name='ALF_production_tau',
    #        element=ALF_production,
    #        kind='coupled',
    #        units='h',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_production_tau(i):
    #     ALF_production.tau=i
    
    # # ALF crystallization
    # ALF_crystallizer = unit.ALF_crystallizer
    
    # dist = shape.
    # @param(name='ALF_crystallizer_tau',
    #        element=ALF_crystallizer,
    #        kind='coupled',
    #        units='h',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_crystallizer_tau(i):
    #     ALF_crystallizer.tau=i
    
    # dist = shape.
    # @param(name='ALF_crystallizer_yield',
    #        element=ALF_crystallizer,
    #        kind='coupled',
    #        units='-',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_crystallizer_yield(i):
    #     ALF_crystallizer.crystal_ALF_yield=i
    
    # # ALF pressure filter
    # ALF_filter = unit.ALF_filter
    
    # dist = shape.
    # @param(name='ALF_filter_moisture_content',
    #        element=ALF_filter,
    #        kind='coupled',
    #        units='-',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_filter_moisture_content(i):
    #     ALF_filter.moisture_content=i
    
    # # reverse osmosis
    # RO = unit.Reverse_osmosis
    
    # dist = shape.
    # @param(name='RO_water_recovery',
    #        element=RO,
    #        kind='coupled',
    #        units='-',
    #        baseline=,
    #        distribution=dist)
    # def set_RO_water_recovery(i):
    #     RO.water_recovery=i
    
    # # storage tank
    # ALF_storage_tank = unit.ALF_storage_tank
    
    # dist = shape.
    # @param(name='ALF_storage_tank_tau',
    #        element=ALF_storage_tank,
    #        kind='coupled',
    #        units='h',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_storage_tank_tau(i):
    #     ALF_storage_tank.tau=i
    
    # # TEA
    # tea = sys.TEA
    
    # dist = shape.
    # @param(name='IRR',
    #        element='TEA',
    #        kind='isolated',
    #        units='-',
    #        baseline=,
    #        distribution=dist)
    # def set_IRR(i):
    #     tea.IRR=i
    
    # # TODO: add a parameter to set the baseline electricity price?
    # dist = shape.
    # @param(name='electricity_price',
    #        element='TEA',
    #        kind='isolated',
    #        units='$/kg',
    #        baseline=,
    #        distribution=dist)
    # def set_electricity_price(i):
    #     bst.PowerUtility.price=i
    
    # aluminum_hydroxide = stream.aluminum_hydroxide
    
    # dist = shape.
    # @param(name='aluminum_hydroxide_price',
    #        element='TEA',
    #        kind='isolated',
    #        units='$/kg',
    #        baseline=,
    #        distribution=dist)
    # def set_aluminum_hydroxide_price(i):
    #     aluminum_hydroxide.price=i
    
    # dist = shape.
    # @param(name='water_price',
    #        element='TEA',
    #        kind='isolated',
    #        units='$/kg',
    #        baseline=,
    #        distribution=dist)
    # def set_water_price(i):
    #     bst.stream_prices['Reverse osmosis water']=i
    
    # formic_acid = stream.formic_acid
    
    # dist = shape.
    # @param(name='formic_acid_price',
    #        element='TEA',
    #        kind='isolated',
    #        units='$/kg',
    #        baseline=,
    #        distribution=dist)
    # def set_formic_acid_price(i):
    #     formic_acid.price=i
    
    # dist = shape.
    # @param(name='natural_gas_price',
    #        element='TEA',
    #        kind='isolated',
    #        units='$/kg',
    #        baseline=,
    #        distribution=dist)
    # def set_natural_gas_price(i):
    #     bst.stream_prices['Natural gas']=i
    
    # # LCA
    # qs.ImpactItem.get_all_items()
    # for item in qs.ImpactItem.get_all_items().keys():
    #     for CF in qs.ImpactIndicator.get_all_indicators().keys():
    #         abs_small = 0.9*qs.ImpactItem.get_item(item).CFs[CF]
    #         abs_large = 1.1*qs.ImpactItem.get_item(item).CFs[CF]
    #         dist = shape.Uniform(min(abs_small,abs_large),max(abs_small,abs_large))
    #         @param(name=f'{item}_{CF}',
    #                setter=DictAttrSetter(qs.ImpactItem.get_item(item), 'CFs', CF),
    #                element='LCA',
    #                kind='isolated',
    #                units=qs.ImpactIndicator.get_indicator(CF).unit,
    #                baseline=qs.ImpactItem.get_item(item).CFs[CF],
    #                distribution=dist)
    #         def set_LCA(i):
    #             qs.ImpactItem.get_item(item).CFs[CF]=i

# =============================================================================
# ALF production: Bauxite + HCOOH
# =============================================================================
def create_model_B(system=None): pass

    # sys = create_system_B() if not system else system
    # flowsheet = sys.flowsheet
    # unit = flowsheet.unit
    # stream = flowsheet.stream
    # model = qs.Model(sys)
    # param = model.parameter
    
    # # TODO: if it is necessary to add bauxite_Al2O3 and bauxite_SiO2 as uncertainty parameters,
    # # we can create a fake unit and and these as parameters
    
    # # ALF production
    # ALF_production = unit.ALF_production
    
    # dist = shape.
    # @param(name='ALF_production_T',
    #        element=ALF_production,
    #        kind='coupled',
    #        units='K',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_production_T(i):
    #     ALF_production.T=i
    
    # dist = shape.
    # @param(name='ALF_production_tau',
    #        element=ALF_production,
    #        kind='coupled',
    #        units='h',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_production_tau(i):
    #     ALF_production.tau=i
    
    # # solid_filter
    # solid_filter = unit.solid_filter
    
    # dist = shape.
    # @param(name='solid_filter_moisture_content',
    #        element=solid_filter,
    #        kind='coupled',
    #        units='-',
    #        baseline=,
    #        distribution=dist)
    # def set_solid_filter_moisture_content(i):
    #     solid_filter.moisture_content=i
    
    # # ALF crystallization
    # ALF_crystallizer = unit.ALF_crystallizer
    
    # dist = shape.
    # @param(name='ALF_crystallizer_tau',
    #        element=ALF_crystallizer,
    #        kind='coupled',
    #        units='h',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_crystallizer_tau(i):
    #     ALF_crystallizer.tau=i
    
    # dist = shape.
    # @param(name='ALF_crystallizer_yield',
    #        element=ALF_crystallizer,
    #        kind='coupled',
    #        units='-',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_crystallizer_yield(i):
    #     ALF_crystallizer.crystal_ALF_yield=i
    
    # # ALF pressure filter
    # ALF_filter = unit.ALF_filter
    
    # dist = shape.
    # @param(name='ALF_filter_moisture_content',
    #        element=ALF_filter,
    #        kind='coupled',
    #        units='-',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_filter_moisture_content(i):
    #     ALF_filter.moisture_content=i
    
    # # reverse osmosis
    # RO = unit.Reverse_osmosis
    
    # dist = shape.
    # @param(name='RO_water_recovery',
    #        element=RO,
    #        kind='coupled',
    #        units='-',
    #        baseline=,
    #        distribution=dist)
    # def set_RO_water_recovery(i):
    #     RO.water_recovery=i
    
    # # storage tank
    # ALF_storage_tank = unit.ALF_storage_tank
    
    # dist = shape.
    # @param(name='ALF_storage_tank_tau',
    #        element=ALF_storage_tank,
    #        kind='coupled',
    #        units='h',
    #        baseline=,
    #        distribution=dist)
    # def set_ALF_storage_tank_tau(i):
    #     ALF_storage_tank.tau=i
    
    # # TEA
    # tea = sys.TEA
    
    # dist = shape.
    # @param(name='IRR',
    #        element='TEA',
    #        kind='isolated',
    #        units='-',
    #        baseline=,
    #        distribution=dist)
    # def set_IRR(i):
    #     tea.IRR=i
    
    # # TODO: add a parameter to set the baseline electricity price?
    # dist = shape.
    # @param(name='electricity_price',
    #        element='TEA',
    #        kind='isolated',
    #        units='$/kg',
    #        baseline=,
    #        distribution=dist)
    # def set_electricity_price(i):
    #     bst.PowerUtility.price=i
    
    # bauxite_ore = stream.bauxite_ore
    
    # dist = shape.
    # @param(name='bauxite_ore_price',
    #        element='TEA',
    #        kind='isolated',
    #        units='$/kg',
    #        baseline=,
    #        distribution=dist)
    # def set_bauxite_ore_price(i):
    #     bauxite_ore.price=i
    
    # dist = shape.
    # @param(name='water_price',
    #        element='TEA',
    #        kind='isolated',
    #        units='$/kg',
    #        baseline=,
    #        distribution=dist)
    # def set_water_price(i):
    #     bst.stream_prices['Reverse osmosis water']=i
    
    # formic_acid = stream.formic_acid
    
    # dist = shape.
    # @param(name='formic_acid_price',
    #        element='TEA',
    #        kind='isolated',
    #        units='$/kg',
    #        baseline=,
    #        distribution=dist)
    # def set_formic_acid_price(i):
    #     formic_acid.price=i
    
    # # TODO: the price of solid_waste should be negative
    # solid_waste = stream.solid_waste
    
    # dist = shape.
    # @param(name='solid_waste_price',
    #        element='TEA',
    #        kind='isolated',
    #        units='$/kg',
    #        baseline=,
    #        distribution=dist)
    # def set_solid_waste_price(i):
    #     solid_waste.price=i
    
    # dist = shape.
    # @param(name='natural_gas_price',
    #        element='TEA',
    #        kind='isolated',
    #        units='$/kg',
    #        baseline=,
    #        distribution=dist)
    # def set_natural_gas_price(i):
    #     bst.stream_prices['Natural gas']=i
    
    # # LCA
    # qs.ImpactItem.get_all_items()
    # for item in qs.ImpactItem.get_all_items().keys():
    #     for CF in qs.ImpactIndicator.get_all_indicators().keys():
    #         abs_small = 0.9*qs.ImpactItem.get_item(item).CFs[CF]
    #         abs_large = 1.1*qs.ImpactItem.get_item(item).CFs[CF]
    #         dist = shape.Uniform(min(abs_small,abs_large),max(abs_small,abs_large))
    #         @param(name=f'{item}_{CF}',
    #                setter=DictAttrSetter(qs.ImpactItem.get_item(item), 'CFs', CF),
    #                element='LCA',
    #                kind='isolated',
    #                units=qs.ImpactIndicator.get_indicator(CF).unit,
    #                baseline=qs.ImpactItem.get_item(item).CFs[CF],
    #                distribution=dist)
    #         def set_LCA(i):
    #             qs.ImpactItem.get_item(item).CFs[CF]=i

# =============================================================================
# CO2 capture and utilization
# =============================================================================
def create_model_C(system=None, upgrade=False): pass
    
    # # note if a system is given, the parameter 'upgrade' will be ignored
    # sys = create_system_C(upgrade=upgrade) if not system else system
    # flowsheet = sys.flowsheet
    # unit = flowsheet.unit
    # stream = flowsheet.stream
    # model = qs.Model(sys)
    # param = model.parameter
    
    # # TODO: is it necessary to add CO2 purity and CO2 recovery as uncertainty parameters?
    # # if so, decide how to do this
    # # TODO: how to add electrolyzer price as an uncertainty parameter?
    
    # # TSA
    # TSA = unit.ALF_TSA
    
    # dist = shape.
    # @param(name='TSA_adsorbent_lifetime',
    #        element=TSA,
    #        kind='coupled',
    #        units='year',
    #        baseline=,
    #        distribution=dist)
    # def set_TSA_adsorbent_lifetime(i):
    #     TSA.adsorbent_lifetime=i
    
    # dist = shape.
    # @param(name='TSA_cycle_time',
    #        element=TSA,
    #        kind='coupled',
    #        units='h',
    #        baseline=,
    #        distribution=dist)
    # def set_TSA_cycle_time(i):
    #     TSA.cycle_time=i
    
    # dist = shape.
    # @param(name='TSA_adsorbent_capacity',
    #        element=TSA,
    #        kind='coupled',
    #        units='g CO2/g ALF',
    #        baseline=,
    #        distribution=dist)
    # def set_TSA_adsorbent_capacity(i):
    #     TSA.adsorbent_capacity=i
    
    # # CO2 electrolyzer
    # try: electrolyzer = unit.CO2_electrolyzer
    # except AttributeError: pass
    # else:
    #     dist = shape.
    #     @param(name='electrolyzer_current_density',
    #            element=electrolyzer,
    #            kind='coupled',
    #            units='A/cm2',
    #            baseline=,
    #            distribution=dist)
    #     def set_electrolyzer_current_density(i):
    #         electrolyzer.current_density=i
        
    #     dist = shape.
    #     @param(name='electrolyzer_cell_voltage',
    #            element=electrolyzer,
    #            kind='coupled',
    #            units='V',
    #            baseline=,
    #            distribution=dist)
    #     def set_electrolyzer_cell_voltage(i):
    #         electrolyzer.cell_voltage=i
        
    #     dist = shape.
    #     @param(name='electrolyzer_product_selectivity',
    #            element=electrolyzer,
    #            kind='coupled',
    #            units='-',
    #            baseline=,
    #            distribution=dist)
    #     def set_electrolyzer_product_selectivity(i):
    #         electrolyzer.product_selectivity=i
        
    #     dist = shape.
    #     @param(name='electrolyzer_conversion',
    #            element=electrolyzer,
    #            kind='coupled',
    #            units='-',
    #            baseline=,
    #            distribution=dist)
    #     def set_electrolyzer_conversion(i):
    #         electrolyzer.conversion=i
        
    #     # TEA
    #     dist = shape.
    #     @param(name='water_price',
    #            element='TEA',
    #            kind='isolated',
    #            units='$/kg',
    #            baseline=,
    #            distribution=dist)
    #     def set_water_price(i):
    #         bst.stream_prices['Reverse osmosis water']=i
            
    #     hydrogen = stream.hydrogen
        
    #     dist = shape.
    #     @param(name='hydrogen_price',
    #            element='TEA',
    #            kind='isolated',
    #            units='$/kg',
    #            baseline=,
    #            distribution=dist)
    #     def set_hydrogen_price(i):
    #         hydrogen.price=i
            
    #     # TODO: add adsorbent cost here
        
    #     # LCA
    #     # TODO: add adsorbent CI here
    
    # # TEA
    # tea = sys.TEA
    
    # dist = shape.
    # @param(name='IRR',
    #        element='TEA',
    #        kind='isolated',
    #        units='-',
    #        baseline=,
    #        distribution=dist)
    # def set_IRR(i):
    #     tea.IRR=i
    
    # # LCA
    # qs.ImpactItem.get_all_items()
    # for item in qs.ImpactItem.get_all_items().keys():
    #     for CF in qs.ImpactIndicator.get_all_indicators().keys():
    #         abs_small = 0.9*qs.ImpactItem.get_item(item).CFs[CF]
    #         abs_large = 1.1*qs.ImpactItem.get_item(item).CFs[CF]
    #         dist = shape.Uniform(min(abs_small,abs_large),max(abs_small,abs_large))
    #         @param(name=f'{item}_{CF}',
    #                setter=DictAttrSetter(qs.ImpactItem.get_item(item), 'CFs', CF),
    #                element='LCA',
    #                kind='isolated',
    #                units=qs.ImpactIndicator.get_indicator(CF).unit,
    #                baseline=qs.ImpactItem.get_item(item).CFs[CF],
    #                distribution=dist)
    #         def set_LCA(i):
    #             qs.ImpactItem.get_item(item).CFs[CF]=i