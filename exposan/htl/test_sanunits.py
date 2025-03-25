#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, math, pandas as pd, qsdsan as qs
from qsdsan import SanUnit
from qsdsan.utils import auom

__all__ = ('create_geospatial_system',)

_mile_to_km = auom('mile').conversion_factor('km')

# GDPCTPI (Gross Domestic Product: Chain-type Price Index)
# https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20)
GDPCTPI = {2007: 86.352,
           2008: 87.977,
           2009: 88.557,
           2010: 89.619,
           2011: 91.466,
           2012: 93.176,
           2013: 94.786,
           2014: 96.436,
           2015: 97.277,
           2016: 98.208,
           2017: 100.000,
           2018: 102.290,
           2019: 104.008,
           2020: 105.407,
           2021: 110.220,
           2022: 117.995,
           2023: 122.284}

folder = os.path.dirname(__file__)

__all__ = (
    # 'AerobicDigester',
    # 'HeatDrying',
    # 'LandApplication',
    'Landfilling',
    # 'LimeStabilization',
    )

# =============================================================================
# Landfilling
# =============================================================================

state_tipping_fee = pd.read_excel(folder + '/data/landfilling_tipping_fee.xlsx', 'data')

# TODO: check
class Landfilling(SanUnit):
    '''
    Municipal solid waste (MWS, e.g., sludge, biosolids) disposal through landfilling,
    including hauling. Tipping fee based on [1]. Hauling fee based on [2].
    Fugutive emission based on [3]. Hauling-related emission based on [4].
    
    Parameters
    ----------
    ins : iterable
        MSW (municipal solid waste).
    outs : iterable
        Landfilled_MSW, emission.
    MSW_wet_density : float
        MSW wet density (assuming sludge with 80% moisture content), [kg·m-3].
    hauling_distance : float
        The hauling distance from wastewater treatment plants to landfills, [km].
    k_methane : float
        Methane generation rate, [year-1].
    L_methane : float
        Potential methane generation capacity, [m3 methane·tonne-1 wastewater solids]
    methane_ratio : float
        Methane ratio in landfill gas, [-].
    methane_density : float
        Methane density at 15 C and 1 atm, [kg·m-3].
        Default to 0.68 kg·m-3 based on [5].
    carbon_dioxide_density : float
        Carbon dioxide density, [kg·m-3].
    lifetime : float
        Lifetime of WRRF sending wastewater solids, [year].
    emission_time : float
        Time duration for emission accounting, can be float(inf), [year].
    other_pollutants : list
        Pollutant and their concentration and molecular weight should be in the format
        of ['pollutant','concentration','molecular_weight'].
        Concentration in [ppmv] and molecular_weight in [g·mol-1].
    US_state : str
        U.S. state where the landfill is located, use 'US' for the U.S. averge data,
        will be ignored if tipping_fee is given.
    tipping_fee : float
        Tipping fee, optional, [USD·tonne-1].
    operation_hours : float
        Yearly operation hour, [h·year-1].
    
    References
    ----------
    .. [1] Analyzing Municipal Solid Waste Landfill Tipping Fees — 2023;
        The Environmental Research & Education Foundation, 2024.
    .. [2] Marufuzzaman, M.; Ekşioğlu, S. D.; Hernandez, R.
        Truck versus Pipeline Transportation Cost Analysis of Wastewater Sludge.
        Transportation Research Part A: Policy and Practice 2015, 74, 14–30.
        https://doi.org/10.1016/j.tra.2015.02.001.
    .. [3] Landfill Gas Emissions Model (LandGEM), 2023.
        https://www.epa.gov/land-research/landfill-gas-emissions-model-landgem.
    .. [4] Ecoinvent 3.8 Database. Swiss Centre for Life Cycle Inventories
        (accessed 2023-01-06).
    .. [5] El Abbadi, S. H.; Feng, J.; Hodson, A. R.; Amouamouha, M.;
        Busse, M. M.; Polcuch, C.; Zhou, P.; Macknick, J.; Guest, J. S.;
        Stokes-Draut, J. R.; Dunn, J. B. Benchmarking Greenhouse Gas Emissions
        from U.S. Wastewater Treatment for Targeted Reduction. 2024.
        https://doi.org/10.31223/X5VQ59.
    '''
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 MSW_wet_density=1040, hauling_distance=13,
                 k_methane=0.05, L_methane=170,
                 methane_ratio=0.5, methane_density=0.68,
                 carbon_dioxide_density=0.68/16.04*44.009,
                 lifetime=30, emission_time=30,
                 other_pollutants=[],
                 US_state=None, tipping_fee=None,
                 operation_hours=8760):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.MSW_wet_density=MSW_wet_density
        self.hauling_distance = hauling_distance
        self.k_methane = k_methane
        self.L_methane = L_methane
        self.methane_ratio = methane_ratio
        self.methane_density = methane_density
        self.carbon_dioxide_density = carbon_dioxide_density
        self.lifetime = lifetime
        self.emission_time = emission_time
        self.other_pollutants = other_pollutants
        self.US_state = US_state
        self.tipping_fee = tipping_fee
        self.operation_hours = operation_hours
    
    def _run(self):
        MSW = self.ins[0]
        landfilled_MSW, emission = self.outs
        
        emission.phase = 'g'
        
        if self.lifetime > self.emission_time:
            self.lifetime = self.emission_time
        
        # methane emission is based on the total MSW weight (including moisture and ash)
        methane_m3_per_h = self.L_methane*MSW.F_mass/1000*(1 - (math.e)**(-self.k_methane*self.emission_time))
        
        emission.imass['CH4'] = methane_m3_per_h*self.methane_density
        
        carbon_dioxide_m3_per_h = methane_m3_per_h/self.methane_ratio*(1 - self.methane_ratio)
        
        emission.imass['CO2'] = carbon_dioxide_m3_per_h*self.carbon_dioxide_density
        
        for pollutant in self.other_pollutants:
            emission.imass[pollutant[0]] = (methane_m3_per_h + carbon_dioxide_m3_per_h)/\
                1000000*pollutant[1]*0.68/16.04*pollutant[2]

        for component in qs.get_components():
            landfilled_MSW.imass[component.ID] = MSW.imass[component.ID]*(MSW.F_mass - emission.F_mass)/MSW.F_mass
        
        if self.tipping_fee:
            if self.US_state:
                Warning('US_state will be ignored since tipping_fee is given.')
                # set price for landfilled_MSW, since MSW is usually from another unit
            landfilled_MSW.price = -self.tipping_fee/1000*MSW.F_mass/landfilled_MSW.F_mass
        elif self.US_state:
            try:
                # set price for landfilled_MSW, since MSW is usually from another unit
                landfilled_MSW.price = -state_tipping_fee[state_tipping_fee['state']==self.US_state]\
                    ['average_tipping_fee'].iloc[0]/1000*MSW.F_mass/landfilled_MSW.F_mass
            except IndexError:
                raise KeyError(f'{self.US_state} is not a valid U.S. stata name, use full name, ' +
                               'with first letter capitalized, use "US" for the U.S. average data')
        else:
            raise ValueError('Either "tipping_fee" or "US_state" must be given.')
    
    def _design(self):
        MSW = self.ins[0]
        landfilled_MSW, emission = self.outs
        
        GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                           method='TRACI',
                                           category='environmental impact',
                                           unit='kg CO2-eq',
                                           description='Global Warming Potential')
        
        try:
            qs.ImpactItem.get_all_items()['MSW_hauling']
            
        except KeyError:
            MSW_hauling = qs.ImpactItem('MSW_hauling', functional_unit='tonne*km')
            MSW_hauling.add_indicator(GlobalWarming, 0.13004958)
            # price for per functional unit
            # 4.56 $/m3, 0.072 $/m3/mile (likely 2015$, [2])
            MSW_hauling.price = (4.56/self.MSW_wet_density+0.072/_mile_to_km/\
                                 self.MSW_wet_density*self.hauling_distance)/\
                GDPCTPI[2015]*GDPCTPI[2022]/self.hauling_distance*1000
            
            MSW_transportation = qs.Transportation('MSW_transportation',
                                                   linked_unit=self,
                                                   item=MSW_hauling,
                                                   load_type='mass',
                                                   load=MSW.F_mass/1000,
                                                   load_unit='tonne',
                                                   distance=self.hauling_distance,
                                                   distance_unit='km',
                                                   # set to 1 h since load = tonne/h
                                                   interval='1',
                                                   interval_unit='h')
            self.transportation = MSW_transportation
        
        try:
            qs.ImpactItem.get_all_items()['landfill_methane']
        except KeyError:
            qs.StreamImpactItem(ID='landfill_methane',
                                linked_stream=emission,
                                GlobalWarming=29.8*emission.imass['CH4']/emission.F_mass)
    
    @property
    def cumulative_emission(self):
        # methane emission is based on the total MSW weight (including moisture and ash)
        return self.L_methane*self.ins[0].F_mass/1000*self.operation_hours*\
            (self.lifetime -\
             1/self.k_methane*(math.e)**(-self.k_methane*(self.emission_time - self.lifetime)) +\
             1/self.k_methane*(math.e)**(-self.k_methane*self.emission_time))