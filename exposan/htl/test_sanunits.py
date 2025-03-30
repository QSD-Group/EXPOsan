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
from biosteam.units.decorators import cost
from biosteam.units.design_tools import CEPCI_by_year

_ton_to_tonne = auom('ton').conversion_factor('tonne')

folder = os.path.dirname(__file__)

__all__ = (
    # 'AerobicDigester',
    # 'HeatDrying',
    'LandApplication',
    'Landfilling',
    'LimeStabilization',
    )

# =============================================================================
# LandApplication
# =============================================================================

class LandApplication(SanUnit):
    '''
    Biosolids application to e.g., agricultural fields, not including hauling.
    Selling prices are from [1-3]. Fugutive emission is based on [4]. 
    [1]: Class A-EQ (AA) biosolids: $60/dry ton, Class B: 0
    [2]: Class A biosolids: $5-15/wet ton, Class B: 0
    [3]: Alkaline stabilized biosolids: $3-5/wet ton.
    Overall, assume Class A-EQ and Class A biosolids have the same selling price of
    $10/wet tonne and Class B biosolids is free.
    
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
        Pollutants and their concentrations, assume trace amounts,thus
        emission.F_mass = emission.imass['N2O'], should be in the format
        of [['pollutant','concentration'], ...]. Concentration in [ppmw].
    biosolids_grade : str
        The grade of biosolids, can be 'A-EQ', 'AA', 'A', or 'B', will be ignored
        if biosolids_price is given.
    biosolids_price : float
        Biosolids selling price, [USD·wet tonne-1].
    
    References
    ----------
    .. [1] https://efotg.sc.egov.usda.gov/references/Public/FL/ECN_FL_9.pdf
        (accessed 2025-03-26).
    .. [2] https://www.unionleader.com/news/environment/with-fertilizer-costs-\
        skyrocketing-farmers-find-a-friend-in-new-hampshire-sludge/article_\
        4604d78a-1cbb-556c-a739-a5f17e4aaa51.html (accessed 2025-03-26).
    .. [3] https://www.epa.gov/sites/default/files/2018-11/documents/alkaline-\
        stabilization-biosolids-factsheet.pdf (accessed 2025-03-27).
    .. [4] El Abbadi, S. H.; Feng, J.; Hodson, A. R.; Amouamouha, M.;
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
                 other_pollutants=[], biosolids_grade=None,
                 biosolids_price=None):
        
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
            emission.imass[pollutant[0]] = emission.imass['N2O']/1e6*pollutant[1]
        
        for component in qs.get_components():
            if component.ID in ['water','H2O']:
                applied_biosolids.imass[component.ID] = biosolids.imass[component.ID]
            else:
                applied_biosolids.imass[component.ID] = biosolids.imass[component.ID]/\
                    (biosolids.F_mass - biosolids.imass['H2O'])*\
                        (biosolids.F_mass - biosolids.imass['H2O'] - emission.F_mass)
        
        if self.biosolids_grade not in ['A-EQ','AA','A','B']:
            raise ValueError('biosolids_grade must be one of the following: "A-EQ", "AA", "A", or "B".')
        
        if self.biosolids_price:
            if self.biosolids_grade:
                Warning('biosolids_grade will be ignored since biosolids_price is given.')
            # set price for applied_biosolids, since biosolids is usually from another unit
            applied_biosolids.price = self.biosolids_price/1000*biosolids.F_mass/applied_biosolids.F_mass
        elif self.biosolids_grade:
            if self.biosolids_grade in ['A-EQ','AA','A']:
                applied_biosolids.price = 10/1000*biosolids.F_mass/applied_biosolids.F_mass
            else:
                applied_biosolids.price = 0
        else:
            raise ValueError('Either "biosolids_grade" or "biosolids_price" must be given.')

# =============================================================================
# Landfilling
# =============================================================================

# tipping fee in $·wet ton-1
state_tipping_fee = pd.read_excel(folder + '/data/landfilling_tipping_fee.xlsx', 'data')

class Landfilling(SanUnit):
    '''
    Wastewater solids (i.e., sludge, biosolids) co-disposal through landfilling
    (with other municipal solid wastes), not including hauling. Methane emission
    is based on dry wastewater solids weight, since wastewater solids have high
    moisture content (e.g., ~80%), while the landfill-wise moisture content is
    around 20-30%, which is similar to ash content of wastewater solids,
    therefore, the current model does not consider water but consider ash for
    methane emission calculation (as if water is considered but ash is ignored).
    
    It is assumed for a given amount of wastewater solids feed rate (M), the
    cumulative emission at emission_time (T) will be emitted at one time as
    the effluent, based on the following equation:
        ∫0->T k_methane*L_methane*M*dt = L_methane*M*(1 - e^(-k_methane*T))
    
    The property cumulative_emission represents the all emissions from the time
    landfills start to accept wastewater solids to the time stopping accounting
    for emissions (tau). If tau is less than the lifetime of wastewater treatment
    plants or landfills, depending on which is the studied object (both denoted
    as t), then t=tau. The cumulative_emission is calculated based on the
    following equations:
        ∫0->t∫y->tau k_methane*L_methane*M*e^(-k_methane*(x - y))*dxdy = 
        L_methane*M*(t - 1/k_methane*e^(-k_methane*(tau - t)) +
        1/k_methane*e^(-k_methane*tau))
    
    Tipping fee is based on [1]. Fugutive emission is based on [2],
    but using integrals.
    
    Parameters
    ----------
    ins : iterable
        wastewater_solids.
    outs : iterable
        landfilled_solids, emission.
    k_methane : float
        Methane generation rate, [year-1].
    L_methane : float
        Potential methane generation capacity, based on pseudo-dry wastewater
        solids, [m3 methane·tonne-1 wastewater solids].
    methane_ratio : float
        Methane ratio in landfill gas, [-].
    methane_density : float
        Methane density at 15 C and 1 atm, [kg·m-3].
        Default to 0.68 kg·m-3 based on [3].
    lifetime : float
        Lifetime of wastewater treatment plants producing wastewater solids or
        landfills, depending on which is the studied object, [year].
    emission_time : float
        Time duration for effluent one-time emission calculation, 
        an be float('inf'), [year].
    stop_time : float
        Stop time for cumulative emission calculation, can be float('inf'), [year].
    other_pollutants : list
        Pollutants and their concentrations and molecular weights, assume trace
        amounts,thus emission.F_vol = emission.ivol['CH4'] + emission.ivol['CO2'],
        should be in the format of [['pollutant','concentration','molecular_weight'],
        ...]. Concentration in [ppmv] and molecular_weight in [g·mol-1].
    US_state : str
        U.S. state where the landfill is located, use 'US' for the U.S. averge data,
        will be ignored if tipping_fee is given.
    tipping_fee : float
        Tipping fee, [USD·wet tonne-1].
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
                 lifetime=30, emission_time=1,
                 stop_time=float('inf'),
                 other_pollutants=[],
                 US_state=None, tipping_fee=None,
                 operation_hours=8760):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.k_methane = k_methane
        self.L_methane = L_methane
        self.methane_ratio = methane_ratio
        self.methane_density = methane_density
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
        
        # methane MW: 16.04 g/mol, CO2 MW: 44.009 g/mol
        emission.imass['CO2'] = carbon_dioxide_m3_per_h*self.methane_density/16.04*44.009
        
        for pollutant in self.other_pollutants:
            emission.imass[pollutant[0]] = (methane_m3_per_h + carbon_dioxide_m3_per_h)/\
                1e6*pollutant[1]*self.methane_density/16.04*pollutant[2]
        
        for component in qs.get_components():
            if component.ID in ['water','H2O']:
                landfilled_solids.imass[component.ID] = wastewater_solids.imass[component.ID]
            else:
                landfilled_solids.imass[component.ID] = wastewater_solids.imass[component.ID]/\
                    (wastewater_solids.F_mass - wastewater_solids.imass['H2O'])*\
                        (wastewater_solids.F_mass - wastewater_solids.imass['H2O'] - emission.F_mass)
        
        if self.tipping_fee:
            if self.US_state:
                Warning('US_state will be ignored since tipping_fee is given.')
            # set price for landfilled_solids, since wastewater_solids is usually from another unit
            landfilled_solids.price = -self.tipping_fee/1000*wastewater_solids.F_mass/landfilled_solids.F_mass
        elif self.US_state:
            try:
                # set price for landfilled_solids, since wastewater_solids is usually from another unit
                landfilled_solids.price = -state_tipping_fee[state_tipping_fee['state']==self.US_state]\
                    ['average_tipping_fee'].iloc[0]/_ton_to_tonne/1000*wastewater_solids.F_mass/landfilled_solids.F_mass
            except IndexError:
                raise KeyError(f'{self.US_state} is not a valid U.S. stata name, use full name, ' +
                               'with first letter capitalized, use "US" for the U.S. average data')
        else:
            raise ValueError('Either "US_state" or "tipping_fee" must be given.')
    
    @property
    def cumulative_emission_m3(self):
        if self.lifetime > self.stop_time:
            self.lifetime = self.stop_time
        return self.L_methane*(self.ins[0].F_mass - self.ins[0].imass['H2O'])/1000*self.operation_hours*\
            (self.lifetime -\
             1/self.k_methane*(math.e)**(-self.k_methane*(self.stop_time - self.lifetime)) +\
             1/self.k_methane*(math.e)**(-self.k_methane*self.stop_time))
    
    @property
    def cumulative_emission_tonne(self):
        if self.lifetime > self.stop_time:
            self.lifetime = self.stop_time
        cumulative_emission_m3 = self.L_methane*(self.ins[0].F_mass - self.ins[0].imass['H2O'])/1000*\
            self.operation_hours*(self.lifetime -\
             1/self.k_methane*(math.e)**(-self.k_methane*(self.stop_time - self.lifetime)) +\
             1/self.k_methane*(math.e)**(-self.k_methane*self.stop_time))
        return cumulative_emission_m3*self.methane_density/1000

# =============================================================================
# LimeStabilization
# =============================================================================
# TODO: check
@cost(ID='Quick lime all parts', basis='Quick lime dry solids flow',
      units='tonne/day', cost=311386, S=1, CE=CEPCI_by_year[2004], n=0.5623, BM=1)
@cost(ID='Hydrated lime all parts', basis='Hydrated lime dry solids flow',
      units='tonne/day', cost=311386*0.8, S=1, CE=CEPCI_by_year[2004], n=0.5623, BM=1)
class LimeStabilization(SanUnit):
    '''
    Lime stabilization of wastewater solids. Assume impurities in lime are CaCO3
    and all lime become CaCO3.
    
    Annualized capital cost (which likely means installed cost) in 2004$
    (30 years, 7% discount rate [r]) from [1-2]:
        1  MGD (~1 tonne dry solids per day)  $  22,600
        4  MGD (~4 tonne dry solids per day)  $  64,700
        40 MGD (~40 tonne dry solids per day) $ 187,500
    
   The total installed costs can be calculated using equations below:
        annualized installed cost = installed cost*r/(1 - (1 + r)^(-lifetime))
    
    The estimated total installed costs are:
        1  MGD (~1 tonne dry solids per day)  $   280,444
        4  MGD (~4 tonne dry solids per day)  $   802,865
        40 MGD (~40 tonne dry solids per day) $ 2,326,695
    
    The total installed cost follows a power function:
        capital cost = 311386*X^0.5623
    
    Assume all costs above are for cases when quick lime is used. For hydrated
    lime, the installed cost reduces 20% since slaking is not needed.
    
    Electricity requirement is based on [3-4].
    
    Parameters
    ----------
    ins : iterable
        wastewater_solids, lime, water, carbon_dioxide.
    outs : iterable
        stabilized_solids, vapor.
    lime_type : str
        Type of lime used, can only be 'quick_lime' or 'hydrated_lime'.
    lime_purity : float
        The purity of lime (i.e., CaO content in quick lime or Ca(OH)2 content
        in hydrated lime), [-].
    water_ratio : float
        Ratio between the weight of added water and the weight of lime, [-].
    vapor_ratio : float
        Ratio of water loss as vapor during quick lime slaking.
    CaOH2_ratio : float
        Lime amount as a ratio of dry wastewater solids, [-].
        From [1]:
            Primary sludge                  -> 0.10-0.15
            Activated sludge                -> 0.30-0.50
            Septage                         -> 0.10-0.30
            Alum-sludge                     -> 0.40-0.60
            Alum-sludge plus primary sludge -> 0.25-0.40
            Iron-sludge                     -> 0.35-0.60
        From [2]:
            Wastewater treatment plants     -> 0.20
    unit_electricity : float
        Electricity for wastewater solids and lime mixing,
        [kWh·dry tonne-1 wastewater solids].
    
    References
    ----------
    .. [1] Williford, C.; Chen, W.-Y.; Shamas, N. K.; Wang, L. K. Lime
        Stabilization. In Biosolids Treatment Processes; Wang, L. K.,
        Shammas, N. K., Hung, Y.-T., Eds.; Humana Press: Totowa, NJ,
        2007; pp 207–241. https://doi.org/10.1007/978-1-59259-996-7_7.
    .. [2] Process Design Manual for Sludge Treatment and Disposal;
        United States Environmental Protection Agency, 1979.
    .. [3] Murray, A.; Horvath, A.; Nelson, K. L. Hybrid Life-Cycle Environmental
        and Cost Inventory of Sewage Sludge Treatment and End-Use Scenarios: A
        Case Study from China. Environ. Sci. Technol. 2008, 42 (9), 3163–3169.
        https://doi.org/10.1021/es702256w.
    .. [4] Liu, H. Novel Approach on Reduction in GHG Emissions from Sludge Lime
        Stabilization as an Emergent and Regional Treatment in China. Sci Rep
        2018, 8 (1), 16564. https://doi.org/10.1038/s41598-018-35052-9.
    '''
    _N_ins = 4
    _N_outs = 2
    _units = {'Quick lime dry solids flow':'tonne/day',
              'Hydrated lime dry solids flow':'tonne/day'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 lime_type='hydrated_lime', lime_purity=None,
                 water_ratio=None, vapor_ratio=0.2,
                 CaOH2_ratio=0.2, unit_electricity=5):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.lime_type = lime_type
        self.lime_purity = lime_purity
        self.water_ratio = water_ratio
        self.vapor_ratio = vapor_ratio
        self.CaOH2_ratio = CaOH2_ratio
        self.unit_electricity = unit_electricity
    
    def _run(self):
        wastewater_solids, lime, water, carbon_dioxide = self.ins
        stabilized_solids, vapor = self.outs
        
        vapor.phase = 'g'
        
        if self.lime_type not in ['quick_lime','hydrated_lime']:
            raise ValueError('lime_type must be one of the following: "quick_lime", "hydrated_lime".')
        
        if self.lime_type == 'quick_lime':
            required_CaOH2 = (wastewater_solids.F_mass - wastewater_solids.imass['H2O'])*self.CaOH2_ratio
                     
            lime.imass['CaO'] = required_CaOH2/74.093*56.0774
            
            # quicklime contains 85% of CaO, [2]
            if self.lime_purity is None: self.lime_purity = 0.85
            
            lime.imass['CaCO3'] = lime.imass['CaO']/self.lime_purity*(1 - self.lime_purity)
            
            # quicklime can be slaked by mixing one part quicklime with two to three parts water, [1]
            if self.water_ratio is None: self.water_ratio = 3
            
            water.imass['H2O'] = lime.F_mass*self.water_ratio
            
            carbon_dioxide.imass['CO2'] = lime.imass['CaO']/56.0774*44.009
            
            stabilized_solids.copy_like(wastewater_solids)
            stabilized_solids.imass['CaCO3'] += lime.imass['CaO']/56.0774*100.0869 + lime.imass['CaCO3']
            
            vapor.imass['H2O'] = water.F_mass*self.vapor_ratio
            
            # assume no water is evaporized during the mixing stage and the
            # temperature does not increase
            # heat generated by slaking does not raise temperature significantly
            # unless the biosolids are dewatered and the lime dose is
            # high on the order of 0.2-0.4, [1]
            stabilized_solids.imass['H2O'] += water.F_mass - vapor.F_mass
        
        else:
            lime.imass['CaOH2'] = (wastewater_solids.F_mass - wastewater_solids.imass['H2O'])*self.CaOH2_ratio
            
            # hydrated lime contains 47% CaO, [2]
            if self.lime_purity is None: self.lime_purity = 0.47/56.0774*74.093
            
            lime.imass['CaCO3'] = lime.imass['CaOH2']/self.lime_purity*(1 - self.lime_purity)
            
            # hydrated lime is fed as a 6–18% Ca(OH)2 slurry by weight, [1]
            if self.water_ratio is None: self.water_ratio = self.lime_purity/0.18 - 1
            
            water.imass['H2O'] = lime.F_mass*self.water_ratio
            
            carbon_dioxide.imass['CO2'] = lime.imass['CaOH2']/74.093*44.009
            
            stabilized_solids.copy_like(wastewater_solids)
            stabilized_solids.imass['CaCO3'] += lime.imass['CaOH2']/74.093*100.0869 + lime.imass['CaCO3']
            stabilized_solids.imass['H2O'] += lime.imass['CaOH2']/74.093*18.01528 + water.F_mass
        
    def _design(self):
        if self.lime_type == 'quick_lime':
            self.design_results['Quick lime dry solids flow'] = (self.ins[0].F_mass - self.ins[0].imass['H2O'])/1000*24
            self.design_results['Hydrated lime dry solids flow'] = 0
        else:
            self.design_results['Quick lime dry solids flow'] = 0
            self.design_results['Hydrated lime dry solids flow'] = (self.ins[0].F_mass - self.ins[0].imass['H2O'])/1000*24
        
        self.add_power_utility(self.unit_electricity*(self.ins[0].F_mass - self.ins[0].imass['H2O'])/1000)
    
    @property
    def moisture_before(self):
        return self.ins[0].imass['H2O']/self.ins[0].F_mass
    
    @property
    def moisture_after(self):
        return self.outs[0].imass['H2O']/self.outs[0].F_mass