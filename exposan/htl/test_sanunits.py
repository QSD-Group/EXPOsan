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

folder = os.path.dirname(__file__)

__all__ = (
    # 'AerobicDigester',
    # 'HeatDrying',
    'LandApplication',
    'Landfilling',
    # 'LimeStabilization',
    )

# =============================================================================
# LandApplication
# =============================================================================

# TODO: check
class LandApplication(SanUnit):
    '''
    Biosolids application to e.g., agricultural fields, not including hauling.
    Selling prices for Class A-EQ and B biosolids are from [1]. The selling price
    for Class A biosolids is from [2]. Fugutive emission is based on [3]. 
    
    Parameters
    ----------
    ins : iterable
        biosolids.
    outs : iterable
        applied_biosolids, emission.
    N_content : float
        Nitrogen content in dry biosolids.
    N_emission_ratio : float
        Ratio of N emitted as N2O from N.
    other_pollutants : list
        Pollutant and their concentration and molecular weight should be in the format
        of ['pollutant','concentration']. Concentration in [ppmw].
    biosolids_grade : str
        The grade of biosolids, can be 'A-EQ', 'A', or 'B'.
    biosolids_price : float
        Biosolids selling price, [USD·wet tonne-1].
    
    References
    ----------
    # TODO: update citations
    
    .. [1] https://efotg.sc.egov.usda.gov/references/Public/FL/ECN_FL_9.pdf
        (accessed 2025-03-26).
    .. [2] https://www.unionleader.com/news/environment/with-fertilizer-costs-\
        skyrocketing-farmers-find-a-friend-in-new-hampshire-sludge/article_\
        4604d78a-1cbb-556c-a739-a5f17e4aaa51.html (accessed 2025-03-26).
    .. [3] El Abbadi, S. H.; Feng, J.; Hodson, A. R.; Amouamouha, M.;
        Busse, M. M.; Polcuch, C.; Zhou, P.; Macknick, J.; Guest, J. S.;
        Stokes-Draut, J. R.; Dunn, J. B. Benchmarking Greenhouse Gas Emissions
        from U.S. Wastewater Treatment for Targeted Reduction. 2024.
        https://doi.org/10.31223/X5VQ59.
    '''
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 N_content=0.05, N_emission_ratio=0.01,
                 other_pollutants=[], biosolids_grade='B',
                 biosolids_price=0):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.N_content = N_content
        self.N_emission_ratio = N_emission_ratio
        self.other_pollutants = other_pollutants
        self.biosolids_grade = biosolids_grade
        self.biosolids_price = biosolids_price
    
    def _run(self):
        biosolids = self.ins[0]
        applied_biosolids, emission = self.outs
        
        emission.phase = 'g'
        
        emission.imass['N2O'] = (biosolids.F_mass - biosolids.imass['H2O'])*\
            self.N_content*self.N_emission_ratio/28.0134*44.013
        
        for pollutant in self.other_pollutants:
            emission.imass[pollutant[0]] = emission.F_mass/1000000*pollutant[1]
        
        for component in qs.get_components():
            if component.ID in ['water','H2O']:
                applied_biosolids.imass[component.ID] = biosolids.imass[component.ID]
            else:
                applied_biosolids.imass[component.ID] = biosolids.imass[component.ID]*\
                    (biosolids.F_mass - biosolids.imass['H2O'] - emission.F_mass)/\
                        (biosolids.F_mass - biosolids.imass['H2O'])
        
        if self.biosolids_grade not in ['A-EQ','A','B']:
            raise ValueError('biosolids_grade must be one of the following: "A-EQ", "A", or "B".')
        
        if self.biosolids_price:
            if self.biosolids_grade:
                Warning('biosolids_grade will be ignored since biosolids_price is given.')
            # set price for applied_biosolids, since biosolids is usually from another unit
            applied_biosolids.price = self.biosolids_price/1000*biosolids.F_mass/applied_biosolids.F_mass
        elif self.biosolids_grade:
            if self.biosolids_grade == 'A-EQ':
                # TODO: is this price for dry solids or wet solids? probably for dry solids based on the table title
                # from [1]: 60 $/ton
                applied_biosolids.price = 66/1000*biosolids.F_mass/applied_biosolids.F_mass
            elif self.biosolids_grade == 'A':
                # TODO: is this price for dry solids or wet solids?
                # from [2]: 5-15 $/ton
                applied_biosolids.price = 11/1000*biosolids.F_mass/applied_biosolids.F_mass
            else:
                applied_biosolids.price = 0
        else:
            raise ValueError('Either "biosolids_grade" or "biosolids_price" must be given.')

# =============================================================================
# Landfilling
# =============================================================================

state_tipping_fee = pd.read_excel(folder + '/data/landfilling_tipping_fee.xlsx', 'data')

# TODO: check
class Landfilling(SanUnit):
    '''
    Wastewater solids (i.e., sludge, biosolids) disposal through landfilling,
    not including hauling. Methane emission is based on dry wastewater solids
    weight, since wastewater solids have high moisture content (e.g., ~80%),
    while the landfill-wise moisture content is around 20-30%, which is similar
    to ash content of wastewater solids, therefore, it is reasonable to not
    consider water but consider ash for methane emission calculation (as if
    water is considered but ash is ignored). Tipping fee is based on [1].
    Fugutive emission is based on [2].
    
    Parameters
    ----------
    ins : iterable
        wastewater_solids.
    outs : iterable
        landfilled_solids, emission.
    k_methane : float
        Methane generation rate, [year-1].
    L_methane : float
        Potential methane generation capacity, [m3 methane·tonne-1 wastewater solids]
    methane_ratio : float
        Methane ratio in landfill gas, [-].
    methane_density : float
        Methane density at 15 C and 1 atm, [kg·m-3].
        Default to 0.68 kg·m-3 based on [3].
    carbon_dioxide_density : float
        Carbon dioxide density, [kg·m-3].
    lifetime : float
        Lifetime of wastewater treatment plants producing wastewater solids, [year].
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
        Tipping fee, optional, [USD·wet tonne-1].
    operation_hours : float
        Yearly operation hour, [h·year-1].
    
    References
    ----------
    .. [1] Analyzing Municipal Solid Waste Landfill Tipping Fees — 2023;
        The Environmental Research & Education Foundation, 2024.
    .. [2] Landfill Gas Emissions Model (LandGEM), 2023.
        https://www.epa.gov/land-research/landfill-gas-emissions-model-landgem.
    .. [3] El Abbadi, S. H.; Feng, J.; Hodson, A. R.; Amouamouha, M.;
        Busse, M. M.; Polcuch, C.; Zhou, P.; Macknick, J.; Guest, J. S.;
        Stokes-Draut, J. R.; Dunn, J. B. Benchmarking Greenhouse Gas Emissions
        from U.S. Wastewater Treatment for Targeted Reduction. 2024.
        https://doi.org/10.31223/X5VQ59.
    '''
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 k_methane=0.05, L_methane=170,
                 methane_ratio=0.5, methane_density=0.68,
                 carbon_dioxide_density=0.68/16.04*44.009,
                 lifetime=30, emission_time=30,
                 other_pollutants=[],
                 US_state=None, tipping_fee=None,
                 operation_hours=8760):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
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
        wastewater_solids = self.ins[0]
        landfilled_solids, emission = self.outs
        
        emission.phase = 'g'
        
        methane_m3_per_h = self.L_methane*(wastewater_solids.F_mass - wastewater_solids.imass['H2O'])/1000*\
            (1 - (math.e)**(-self.k_methane*self.emission_time))
        
        emission.imass['CH4'] = methane_m3_per_h*self.methane_density
        
        carbon_dioxide_m3_per_h = methane_m3_per_h/self.methane_ratio*(1 - self.methane_ratio)
        
        emission.imass['CO2'] = carbon_dioxide_m3_per_h*self.carbon_dioxide_density
        
        for pollutant in self.other_pollutants:
            emission.imass[pollutant[0]] = (methane_m3_per_h + carbon_dioxide_m3_per_h)/\
                1000000*pollutant[1]*0.68/16.04*pollutant[2]
        
        for component in qs.get_components():
            if component.ID in ['water','H2O']:
                landfilled_solids.imass[component.ID] = wastewater_solids.imass[component.ID]
            else:
                landfilled_solids.imass[component.ID] = wastewater_solids.imass[component.ID]*\
                    (wastewater_solids.F_mass - wastewater_solids.imass['H2O'] - emission.F_mass)/(wastewater_solids.F_mass - wastewater_solids.imass['H2O'])
        
        if self.tipping_fee:
            if self.US_state:
                Warning('US_state will be ignored since tipping_fee is given.')
            # set price for landfilled_solids, since wastewater_solids is usually from another unit
            landfilled_solids.price = -self.tipping_fee/1000*wastewater_solids.F_mass/landfilled_solids.F_mass
        elif self.US_state:
            try:
                # set price for landfilled_solids, since wastewater_solids is usually from another unit
                landfilled_solids.price = -state_tipping_fee[state_tipping_fee['state']==self.US_state]\
                    ['average_tipping_fee'].iloc[0]/1000*wastewater_solids.F_mass/landfilled_solids.F_mass
            except IndexError:
                raise KeyError(f'{self.US_state} is not a valid U.S. stata name, use full name, ' +
                               'with first letter capitalized, use "US" for the U.S. average data')
        else:
            raise ValueError('Either "US_state" or "tipping_fee" must be given.')
    
    @property
    def cumulative_emission(self):
        if self.lifetime > self.emission_time:
            self.lifetime = self.emission_time
        return self.L_methane*(self.ins[0].F_mass - self.ins[0].imass['H2O'])/1000*self.operation_hours*\
            (self.lifetime -\
             1/self.k_methane*(math.e)**(-self.k_methane*(self.emission_time - self.lifetime)) +\
             1/self.k_methane*(math.e)**(-self.k_methane*self.emission_time))