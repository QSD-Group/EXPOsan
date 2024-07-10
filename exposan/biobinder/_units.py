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

import math, biosteam as bst, qsdsan as qs
from biosteam.units.decorators import cost
from qsdsan import SanUnit, Stream, sanunits as qsu
from exposan.biobinder import CEPCI_by_year

__all__ = (
    'AqueousFiltration',
    'BiocrudeDeashing',
    'BiocrudeDewatering',
    'BiocrudeSplitter',
    'Disposal',
    # 'GasScrubber',
    'PilotHTL',
    'PreProcessing',
    'ShortcutColumn',
    'Transportation'
    )

salad_dressing_composition = {
    'Water': 0.7566,
    'Lipids': 0.2434*0.6245,
    'Proteins': 0.2434*0.0238,
    'Carbohydrates': 0.2434*0.2946,
    'Ash': 0.2434*0.0571,
    }

class PreProcessing(qsu.MixTank):
    '''
    Adjust the composition and moisture content of the feedstock.
    '''
    _N_ins = 2
    _centralized_dry_flowrate = _decentralized_dry_flowrate = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', F_BM_default=1,
                  feedstock_composition=salad_dressing_composition,
                  decentralized_dry_flowrate=1, # dry kg/hr
                  centralized_dry_flowrate=1, # dry kg/hr
                  target_HTL_solid_loading=0.2,
                  tau=1, **add_mixtank_kwargs,
                  ):
        mixtank_kwargs = add_mixtank_kwargs.copy()
        mixtank_kwargs['tau'] = tau
        qsu.MixTank.__init__(self, ID, ins, outs, thermo, 
                             init_with=init_with, F_BM_default=F_BM_default, **mixtank_kwargs)
        self.feedstock_composition = feedstock_composition
        self.decentralized_dry_flowrate = decentralized_dry_flowrate
        self.centralized_dry_flowrate = centralized_dry_flowrate
        self.target_HTL_solid_loading = target_HTL_solid_loading

    
    def _run(self):
        feedstock, htl_process_water = self.ins
        feedstock.empty()
        htl_process_water.empty()
        
        feedstock_composition = self.feedstock_composition
        for i, j in self.feedstock_composition.items():
            feedstock.imass[i] = j
        
        decentralized_dry_flowrate = self.decentralized_dry_flowrate
        feedstock.F_mass = decentralized_dry_flowrate/(1-feedstock_composition['Water']) # scale flowrate
        htl_wet_mass = decentralized_dry_flowrate/self.target_HTL_solid_loading
        required_water = htl_wet_mass - feedstock.imass['Water']
        htl_process_water.imass['Water'] = max(0, required_water)
        
        qsu.MixTank._run(self)
        
    def _cost(self):
        qsu.MixTank._cost(self)
        N = math.ceil(self.centralized_dry_flowrate / self.decentralized_dry_flowrate)
        self.parallel['self'] *= N
        # baseline_purchase_costs = self.baseline_purchase_costs
        # for i, j in baseline_purchase_costs.items():
        #     baseline_purchase_costs[i] *= N
        # self.power_utility.consumption *= N
        # self.power_utility.production *= N
        

#!!! TO BE UPDATED THROUGHOUT
base_feedstock_flowrate = 11.46 # kg/h
@cost(basis='Feedstock dry flowrate', ID='Feedstock Tank', units='kg/hr',
      cost=4330, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023], n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Feedstock Pump', units='kg/hr',
      cost=6180, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2.3)
@cost(basis='Feedstock dry flowrate', ID= 'Inverter', units='kg/hr',
      cost=240, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'High Pressure Pump', units='kg/hr',
      cost=1634, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2.3)
@cost(basis='Feedstock dry flowrate', ID= 'Reactor Core', units='kg/hr',
      cost=30740, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2)
@cost(basis='Feedstock dry flowrate', ID= 'Reactor Vessel', units='kg/hr',
      cost=4330, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Heat Transfer Putty', units='kg/hr',
      cost=2723, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Electric Heaters', units='kg/hr',
      cost=8400, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'J Type Thermocouples', units='kg/hr',
      cost=497, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Ceramic Fiber', units='kg/hr',
      cost=5154, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Steel Jacket', units='kg/hr',
      cost=22515, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Counterflow Heat Exchanger', units='kg/hr',
      cost=14355, S=base_feedstock_flowrate, CE=CEPCI_by_year[2013],n=0.77, BM=2.2)
@cost(basis='Feedstock dry flowrate', ID= 'Temperature Control and Data Logging Unit', units='kg/hr',
      cost=905, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'Pulsation Dampener', units='kg/hr',
      cost=3000, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'Fluid Accumulator', units='kg/hr',
      cost=995, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'Burst Rupture Discs', units='kg/hr',
      cost=1100, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023], n=0.77, BM=1.6)
@cost(basis='Feedstock dry flowrate', ID= 'Pressure Relief Vessel', units='kg/hr',
      cost=4363, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2)
@cost(basis='Feedstock dry flowrate', ID= 'Gas Scrubber', units='kg/hr',
      cost=1100, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'BPR', units='kg/hr',
      cost=4900, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.6)
@cost(basis='Feedstock dry flowrate', ID= 'Primary Collection Vessel', units='kg/hr',
      cost=7549, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Belt Oil Skimmer', units='kg/hr',
      cost=2632, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Bag Filter', units='kg/hr',
      cost=8800, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.7)
@cost(basis='Feedstock dry flowrate', ID= 'Oil Vessel', units='kg/hr',
      cost=4330, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Mobile HTL system', units='kg/hr',
      cost=23718, S=base_feedstock_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Non-scaling factor', ID='Magnotrol Valves Set', units='ea',
      cost=343, S=1, CE=CEPCI_by_year[2023], n=1, BM=1)
class PilotHTL(SanUnit):
    '''
    
    References
    ----------
    [1] Leow, S.; Witter, J. R.; Vardon, D. R.; Sharma, B. K.;
        Guest, J. S.; Strathmann, T. J. Prediction of Microalgae Hydrothermal
        Liquefaction Products from Feedstock Biochemical Composition.
        Green Chem. 2015, 17 (6), 3584–3599. https://doi.org/10.1039/C5GC00574D.
    [2] Li, Y.; Leow, S.; Fedders, A. C.; Sharma, B. K.; Guest, J. S.;
        Strathmann, T. J. Quantitative Multiphase Model for Hydrothermal
        Liquefaction of Algal Biomass. Green Chem. 2017, 19 (4), 1163–1174.
        https://doi.org/10.1039/C6GC03294J.
    [3] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J.
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products.
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727.
        https://doi.org/10.1021/acs.est.8b04035.
    [4] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    [5] Matayeva, A.; Rasmussen, S. R.; Biller, P. Distribution of Nutrients and
        Phosphorus Recovery in Hydrothermal Liquefaction of Waste Streams.
        BiomassBioenergy 2022, 156, 106323.
        https://doi.org/10.1016/j.biombioe.2021.106323.
    [6] Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels
        via Liquefaction - Hydrothermal Liquefaction Reactor Design:
        April 5, 2013; NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462,
        1111191. https://doi.org/10.2172/1111191.
    '''
    
    _N_ins = 1
    _N_outs = 4
    
    _units= {
        'Feedstock dry flowrate': 'kg/hr',
        'Non-scaling factor': 'ea',
        }
    
    _centralized_dry_flowrate = _decentralized_dry_flowrate = 1
    
    # ID of the components that will be used in mass flowrate calculations
    ash_ID = 'Ash'
    water_ID = 'Water'
    
    # Product condition adjustment based on ref [4]
    gas_composition = {
        'CH4':0.050,
        'C2H6':0.032,
        'CO2':0.918
        }
    
    # Product conditions per [4], pressure converted from psi to Pa
    biocrude_moisture_content = 0.063
    hydrochar_P = 3029.7*6894.76
    HTLaqueous_P = 30*6894.76
    biocrude_P = 30*6894.76
    offgas_P = 30*6894.76
    eff_T = 60+273.15
    

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream',
                  P=None, tau=15/60, V_wf=0.45,
                  decentralized_dry_flowrate=1, # dry kg/hr
                  centralized_dry_flowrate=1, # dry kg/hr
                  afdw_biocrude_yield=0.5219,
                  afdw_aqueous_yield=0.2925,
                  afdw_gas_yield=0.1756,
                  piping_cost_ratio=0.15,
                  accessory_cost_ratio=0.08,
                  **kwargs,
                  ):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        #!!! Need to compare the externally sourced HX cost and BioSTEAM default
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = qsu.HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=self.eff_T, rigorous=True)
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.decentralized_dry_flowrate = decentralized_dry_flowrate
        self.centralized_dry_flowrate = centralized_dry_flowrate
        self.afdw_biocrude_yield = afdw_biocrude_yield
        self.afdw_aqueous_yield = afdw_aqueous_yield
        self.afdw_gas_yield = afdw_gas_yield
        self.piping_cost_ratio = piping_cost_ratio
        self.accessory_cost_ratio = accessory_cost_ratio
        for attr, val in kwargs.items(): setattr(self, attr, val)
    

    def _run(self):        
        feedstock = self.ins[0]
        hydrochar, HTLaqueous, biocrude, offgas = self.outs
        for i in self.outs: i.empty()
        
        afdw_in = self.afdw_mass_in
        hydrochar.imass['Hydrochar'] = afdw_in * self.afdw_hydrochar_yield
        # HTLaqueous is TDS in aqueous phase
        HTLaqueous.imass['HTLaqueous'] = afdw_in * self.afdw_aqueous_yield
        
        gas_mass = afdw_in * self.afdw_gas_yield
        for name, ratio in self.gas_composition.items():
            offgas.imass[name] = gas_mass*ratio
            
        biocrude.imass['Biocrude'] = afdw_in * self.afdw_biocrude_yield
        biocrude.imass['H2O'] = biocrude.imass['Biocrude']/(1 -\
                                self.biocrude_moisture_content) -\
                                biocrude.imass['Biocrude']
                                
        HTLaqueous.imass['H2O'] = feedstock.F_mass - hydrochar.F_mass -\
                                  biocrude.F_mass - gas_mass - HTLaqueous.imass['HTLaqueous']
        # assume ash (all soluble based on Jones) goes to water
        
        hydrochar.phase = 's'
        offgas.phase = 'g'
        HTLaqueous.phase = biocrude.phase = 'l'
        
        hydrochar.P = self.hydrochar_P
        HTLaqueous.P = self.HTLaqueous_P
        biocrude.P = self.biocrude_P
        offgas.P = self.offgas_P
        
        self._refresh_parallel()
        for stream in self.outs:
            stream.T = self.heat_exchanger.T
            stream.F_mass *= self.parallel['self']
        

    def _design(self):
        Design = self.design_results
        Design['Feedstock dry flowrate'] = self.dry_mass_in
        Design['Non-scaling factor'] = 1
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from((self.outs[1], self.outs[2], self.outs[3]))
        hx_outs0.copy_like(hx_ins0)
        hx_ins0.T = self.ins[0].T # temperature before/after HTL are similar
        hx_outs0.T = hx.T
        hx_ins0.P = hx_outs0.P = self.outs[0].P # cooling before depressurized, heating after pressurized
        # in other words, both heating and cooling are performed under relatively high pressure
        hx_ins0.vle(T=hx_ins0.T, P=hx_ins0.P)
        hx_outs0.vle(T=hx_outs0.T, P=hx_outs0.P)
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)

        self.P = self.ins[0].P

        
    def _cost(self):
        self._refresh_parallel()
        self._decorated_cost()
        baseline_purchase_cost = self.baseline_purchase_cost
        self.baseline_purchase_costs['Piping'] = baseline_purchase_cost*self.piping_cost_ratio
        self.baseline_purchase_costs['Accessories'] = baseline_purchase_cost*self.accessory_cost_ratio
        
        
        # # If need to consider additional cost factors
        # purchase_costs = self.baseline_purchase_costs
        # for item in purchase_costs.keys():
        #     purchase_costs[item] *= self.CAPEX_factor
            
        # for aux_unit in self.auxiliary_units:
        #     purchase_costs = aux_unit.baseline_purchase_costs
        #     installed_costs = aux_unit.installed_costs
        #     for item in purchase_costs.keys():
        #         purchase_costs[item] *= self.CAPEX_factor
        #         installed_costs[item] *= self.CAPEX_factor
        

    def _refresh_parallel(self):
        self.parallel['self'] = math.ceil(self.centralized_dry_flowrate/self.decentralized_dry_flowrate)

    @property
    def dry_mass_in(self):
        '''Total dry mass of the feedstock.'''
        feedstock = self.ins[0]
        return feedstock.F_mass-feedstock.imass[self.water_ID]

    @property
    def afdw_mass_in(self):
        '''Total ash-free dry mass of the feedstock.'''
        feedstock = self.ins[0]
        return feedstock.F_mass-feedstock.imass[self.ash_ID]-feedstock.imass[self.water_ID]
    
    @property
    def afdw_hydrochar_yield(self):
        '''Hydrochar product yield on the ash-free dry weight basis of the feedstock.'''
        char_yield = 1-self.afdw_biocrude_yield-self.afdw_aqueous_yield-self.afdw_gas_yield
        if char_yield < 0:
            raise ValueError('Sum of biocrude, aqueous, and gas product exceeds 100%.')
        return char_yield
    
    @property
    def decentralized_dry_flowrate(self):
        '''Dry mass flowrate for the decentralized configuration.'''
        return self._decentralized_dry_flowrate
    @decentralized_dry_flowrate.setter
    def decentralized_dry_flowrate(self, i):
        self._decentralized_dry_flowrate = i
        self._refresh_parallel()

    @property
    def centralized_dry_flowrate(self):
        '''Dry mass flowrate for the centralzied configuration.'''
        return self._centralized_dry_flowrate
    @centralized_dry_flowrate.setter
    def centralized_dry_flowrate(self, i):
        self._centralized_dry_flowrate = i
        self._refresh_parallel()
        

# Jone et al., Table C-1
#!!! Might want to redo this part by adjusting the components.
default_biocrude_ratios = {
    '1E2PYDIN':     0.067912,
    # 'C5H9NS':       0.010257,
    'ETHYLBEN':     0.025467,
    '4M-PHYNO':     0.050934,
    '4EPHYNOL':     0.050934,
    'INDOLE':       0.050934,
    '7MINDOLE':     0.033956,
    'C14AMIDE':     0.033956,
    'C16AMIDE':     0.152801,
    'C18AMIDE':     0.067912,
    'C16:1FA':      0.135823,
    'C16:0FA':      0.101868,
    'C18FACID':     0.016978,
    'NAPHATH':      0.050934,
    'CHOLESOL':     0.016978,
    'AROAMINE':     0.081424,
    'C30DICAD':     0.050934,
    }


class BiocrudeSplitter(SanUnit):
    '''
    Split biocrude into the respective components.
    '''
    _N_ins = _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', F_BM_default=1,
                  cutoff_Tb=273.15+343, light_frac=0.5316,
                  biocrude_IDs=('Biocrude',),
                  biocrude_ratios=default_biocrude_ratios,
                  **kwargs,
                  ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
        self.cutoff_Tb = cutoff_Tb
        self.light_frac = light_frac
        self.biocrude_IDs = biocrude_IDs
        self.biocrude_ratios = biocrude_ratios
        for kw, arg in kwargs.items(): setattr(self, kw, arg)
        
    def _run(self):
        biocrude_in = self.ins[0]
        biocrude_out = self.outs[0]
        
        biocrude_IDs = self.biocrude_IDs
        total_crude = biocrude_in.imass[self.biocrude_IDs].sum()
        total_light = total_crude * self.light_frac
        total_heavy = total_crude - total_light
        
        # Firstly copy the non-biocrude components
        light_ratios, heavy_ratios = self.light_component_ratios, self.heavy_component_ratios
        biocrude_out.copy_like(biocrude_in)
        biocrude_out.imass[[*light_ratios, *heavy_ratios]] = 0
        
        # Set the mass for the biocrude components
        biocrude_out.imass[light_ratios] = [total_light*i for i in light_ratios.values()]
        biocrude_out.imass[heavy_ratios] = [total_heavy*i for i in heavy_ratios.values()]
        biocrude_out.imass[biocrude_IDs] = 0 # clear out biocrude

    def _update_component_ratios(self):
        '''Update the light and heavy ratios of the biocrude components.'''
        if not hasattr(self, 'cutoff_Tb'): return
        if not hasattr(self, 'biocrude_ratios'): return

        cmps = self.components
        Tb = self.cutoff_Tb
        ratios = self.biocrude_ratios
        
        light_ratios = {}
        for ID, ratio in ratios.items():
            if cmps[ID].Tb <= Tb:
                light_ratios[ID] = ratio
                light_key = ID
            else:
                heavy_key = ID
                break
        self._light_key = light_key
        self._heavy_key = heavy_key
        
        heavy_cmps = set(ratios).difference(set(light_ratios))
        heavy_ratios = {ID: ratios[ID] for ratio in heavy_cmps}
        
        # Normalize the ratios
        ratio_light_sum = sum(light_ratios.values())
        ratio_heavy_sum = sum(heavy_ratios.values())
        for ID, ratio in light_ratios.items():
            light_ratios[ID] = ratio/ratio_light_sum
        for ID, ratio in heavy_ratios.items():
            heavy_ratios[ID] = ratio/ratio_heavy_sum
            
        self._light_component_ratios = light_ratios    
        self._heavy_component_ratios = heavy_ratios

    @property
    def cutoff_Tb(self):
        '''[float] Cutoff of the boiling point of light and heavy fractions.'''
        return self._cutoff_Tb
    @cutoff_Tb.setter
    def cutoff_Tb(self, Tb):
        if hasattr(self, '_cutoff_Tb'):
            if Tb == self._cutoff_Tb: return # no need to do anything if Tb unchanged
        self._cutoff_Tb = Tb
        self._update_component_ratios()

    @property
    def light_component_ratios(self):
        '''Mass ratios of the components in the light fraction of the biocrude.'''
        return self._light_component_ratios

    @property
    def heavy_component_ratios(self):
        '''Mass ratios of the components in the heavy fraction of the biocrude.'''
        return self._heavy_component_ratios
    
    @property
    def light_key(self):
        '''ID of the component that has the highest boiling point in the light fraction of the biocrude.'''
        return self._light_key
    
    @property
    def heavy_key(self):
        '''ID of the component that has the lowest boiling point in the heavy fraction of the biocrude.'''
        return self._heavy_key
    
    @property
    def biocrude_ratios(self):
        '''[dict] Mass ratios of the components used to model the biocrude.'''
        return self._biocrude_ratios
    @biocrude_ratios.setter
    def biocrude_ratios(self, ratios):
        cmps = self.components
        # Sort the biocrude ratios by the boiling point
        ratios = {ID: ratio for ID, ratio in 
                  sorted(ratios.items(), key=lambda item: cmps[item[0]].Tb)}
        self._biocrude_ratios = ratios
        self._update_component_ratios()
        

# # Included in the HTL reactor
# class GasScrubber(qsu.Copier):
#     '''
#     Placeholder for the gas scrubber. All outs are copied from ins.
#     '''


base_ap_flowrate = 49.65 #kg/hr
@cost(basis='Aqueous flowrate', ID= 'Sand Filtration Unit', units='kg/hr',
      cost=318, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.7)
@cost(basis='Aqueous flowrate', ID= 'EC Oxidation Tank', units='kg/hr',
      cost=1850, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
@cost(basis='Aqueous flowrate', ID= 'Biological Treatment Tank', units='kg/hr',
      cost=4330, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
@cost(basis='Aqueous flowrate', ID= 'Liquid Fertilizer Storage', units='kg/hr',
      cost=7549, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
class AqueousFiltration(SanUnit):
    '''
    Placeholder for the aqueous filtration unit. All outs are copied from ins.
    '''
     
    _N_outs = 1
    _units= {
        'Aqueous flowrate': 'kg/hr',
        }
    def _run(self):
        HTL_aqueous = self.ins[0]
        treated_aq = self.outs
        
        #treated_aq.copy_like(HTL_aqueous)
        
    def _design(self):
        aqueous = self.ins[0]
        self.design_results['Aqueous flowrate'] = aqueous.F_mass
                
base_biocrude_flowrate = 5.64 # kg/hr
@cost(basis='Biocrude flowrate', ID= 'Deashing Tank', units='kg/hr',
      cost=4330, S=base_biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
class BiocrudeDeashing(SanUnit):
    '''
    Placeholder for the deashing unit.
    '''
    
    _N_outs = 2
    _units= {
        'Biocrude flowrate': 'kg/hr',
        }
    target_ash = 0.01 # dry weight basis
    
    def _run(self):
        biocrude = self.ins[0]
        deashed, ash = self.outs
        
        deashed.copy_like(biocrude)
        ash.empty()
        dw = deashed.F_mass - deashed.imass['Water']
        excess_ash = deashed.imass['Ash'] - dw * self.target_ash
        # Remove excess ash
        if excess_ash >= 0:
            deashed.imass['Ash'] -= excess_ash
            ash.imass['Ash'] = excess_ash
            
    def _design(self):
        biocrude = self.ins[0]
        self.design_results['Biocrude flowrate'] = biocrude.F_mass
            

@cost(basis='Biocrude flowrate', ID= 'Dewatering Tank', units='kg/hr',
      cost=4330, S=base_biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
class BiocrudeDewatering(SanUnit):
    '''
    Placeholder for the dewatering unit.
    '''
    
    _N_outs = 2
    _units= {
        'Biocrude flowrate': 'kg/hr',
        }
    target_moisture = 0.01 # weight basis
    
    def _run(self):
        biocrude = self.ins[0]
        dewatered, water = self.outs
        
        dewatered.copy_like(biocrude)
        water.empty()
        dw = dewatered.F_mass - dewatered.imass['Water']
        excess_water = dw/(1-self.target_moisture) - dw
        # Remove excess water
        if excess_water >= 0:
            dewatered.imass['Water'] -= excess_water
            water.imass['Water'] = excess_water
            
    def _design(self):
        biocrude = self.ins[0]
        self.design_results['Biocrude flowrate'] = biocrude.F_mass

            
class ShortcutColumn(bst.units.ShortcutColumn, qs.SanUnit):
    '''
    Similar to biosteam.units.ShortcutColumn.
    
    See Also
    --------
    `biosteam.units.ShortcutColumn <https://biosteam.readthedocs.io/en/latest/API/units/distillation.html>`_
    '''
            

class Transportation(qsu.Copier):    
    '''
    To account for transportation cost. All outs are copied from ins.
    
    Parameters
    ----------
    transportation_distance : float
        Transportation distance in km.
    transportation_price : float
        Transportation price in $/kg/km.
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', F_BM_default=1,
                  transportation_distance=0, # km
                  transportation_price=0, # $/kg/km
                  **kwargs,
                  ):
        qsu.Copier.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
        self.transportation_distance = transportation_distance
        self.transportation_price = transportation_price
        for kw, arg in kwargs.items(): setattr(self, kw, arg)
    
    def _cost(self):
        self.baseline_purchase_costs['Transportation'] = (
            self.F_mass_in *
            self.transportation_distance *
            self.transportation_price
            )

    
class Disposal(SanUnit):
    '''
    Update the cost for disposal, where the price is given for dry weights.
    '''
    
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', F_BM_default=1,
                  disposal_price=0,
                  exclude_components=('Water',),
                  **kwargs,
                  ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
        self.disposal_price = disposal_price
        self.exclude_components = exclude_components
        for kw, arg in kwargs.items(): setattr(self, kw, arg)
    
    def _run(self):
        inf = self.ins[0]
        waste, others = self.outs        
        
        waste.copy_like(inf)
        waste.imass[self.exclude_components] = 0
        
        others.copy_like(inf)
        others.imass[self.components.IDs] -= waste.imass[self.components.IDs]
        
    def _cost(self):
        self.outs[0].price = self.disposal_price
