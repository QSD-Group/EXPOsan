#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, pandas as pd, qsdsan as qs
from chaospy import distributions as shape
from exposan.saf import (
    results_path,
    create_system,
    _HHV_per_GGE,
    )

__all__ = (
    'create_model',
    )

def create_model(**sys_kwargs):
    sys = create_system(**sys_kwargs)
    flowsheet = sys.flowsheet
    unit = flowsheet.unit
    stream = flowsheet.stream
    model = qs.Model(sys)
    param = model.parameter
    
    feedstock = stream.feedstock
       
    dist = shape.Uniform(0.846,1.034)
    @param(name='ww_2_dry_sludge',
           element=WWTP,
           kind='coupled',
           units='ton/d/MGD',
           baseline=0.94,
           distribution=dist)
    def set_ww_2_dry_sludge(i):
        WWTP.ww_2_dry_sludge=i
    
    
    # TEA
    tea = sys.TEA
    
    - feedstock tipping fee
    - electricity (net metering)
    - EC voltage
    - HTL capital
    - HC capital
    - HT capital
    - CHG capital
    - capital cost
    
    
    #!!! Potential codes for LCA
    # if include_CFs_as_metrics:
    #     qs.ImpactItem.get_all_items().pop('feedstock_item')
    #     for item in qs.ImpactItem.get_all_items().keys():
    #         for CF in qs.ImpactIndicator.get_all_indicators().keys():
    #             abs_small = 0.9*qs.ImpactItem.get_item(item).CFs[CF]
    #             abs_large = 1.1*qs.ImpactItem.get_item(item).CFs[CF]
    #             dist = shape.Uniform(min(abs_small,abs_large),max(abs_small,abs_large))
    #             @param(name=f'{item}_{CF}',
    #                    setter=DictAttrSetter(qs.ImpactItem.get_item(item), 'CFs', CF),
    #                    element='LCA',
    #                    kind='isolated',
    #                    units=qs.ImpactIndicator.get_indicator(CF).unit,
    #                    baseline=qs.ImpactItem.get_item(item).CFs[CF],
    #                    distribution=dist)
    #             def set_LCA(i):
    #                 qs.ImpactItem.get_item(item).CFs[CF]=i

    # =========================================================================
    # Metrics
    # =========================================================================   
    metric = model.metric
    
    # Mass balance
    gasoline = stream.gasoline
    @metric(name='Gasoline yield',units='dw',element='TEA')
    def get_gasoline_yield():
        return gasoline.F_mass/(feedstock.F_mass-feedstock.imass['Water'])

    jet = stream.jet
    @metric(name='Jet yield',units='dw',element='TEA')
    def get_jet_yield():
        return jet.F_mass/(feedstock.F_mass-feedstock.imass['Water'])
    
    diesel = stream.diesel
    @metric(name='Diesel yield',units='dw',element='TEA')
    def get_diesel():
        return diesel.F_mass/(feedstock.F_mass-feedstock.imass['Water'])
    
    #!!!

    mixed_fuel = stream.mixed_fuel
    @metric(name='Annual GGE',units='GGE/yr',element='TEA')
    def get_annual_GGE():
        return mixed_fuel.HHV/1e3/_HHV_per_GGE*sys.operating_hours

    # TEA
    @metric(name='MFSP',units='$/GGE',element='TEA')
    def get_MFSP():
        mixed_fuel.price = sys.TEA.solve_price(mixed_fuel)
        GGE = mixed_fuel.HHV/1e3/_HHV_per_GGE
        return mixed_fuel.cost/GGE
    
    
    #!!! LCA ones
    # @metric(name='GWP_sludge',units='kg CO2/tonne dry sludge',element='LCA')
    # def get_GWP_sludge():
    #     return lca.get_total_impacts(exclude=(raw_wastewater,))['GlobalWarming']/raw_wastewater.F_vol/_m3perh_to_MGD/WWTP.ww_2_dry_sludge/(sys.operating_hours/24)/lca.lifetime