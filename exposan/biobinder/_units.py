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

import biosteam as bst, qsdsan as qs
from math import ceil, log
from biosteam.units.decorators import cost
from qsdsan import SanUnit, Stream, sanunits as qsu
# from qsdsan.sanunits._hydrothermal import KnockOutDrum
# from qsdsan.utils import auom
from exposan.biobinder import CEPCI_by_year

__all__ = (
    'BiocrudeDeashing',
    'BiocrudeDewatering',
    'BiocrudeSplitter',
    # 'GasScrubber',
    'PilotHTL',
    'AqueousFiltration',
    'ShortcutColumn',
    'Transportation'
    )

#!!! TO BE UPDATED THROUGHOUT
pilot_flowrate = 11.46 # kg/h
@cost(basis='Feedstock dry flowrate', ID='Feedstock Tank', units='kg/h',
      cost=4330, S=pilot_flowrate, CE=CEPCI_by_year[2023], n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Feedstock Pump', units='kg/h',
      cost=6180, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2.3)
@cost(basis='Feedstock dry flowrate', ID= 'Inverter', units='kg/h',
      cost=240, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'High Pressure Pump', units='kg/h',
      cost=1634, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2.3)
@cost(basis='Feedstock dry flowrate', ID= 'Reactor Core', units='kg/h',
      cost=30740, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2)
@cost(basis='Feedstock dry flowrate', ID= 'Reactor Vessel', units='kg/h',
      cost=4330, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Heat Transfer Putty', units='kg/h',
      cost=2723, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Electric Heaters', units='kg/h',
      cost=8400, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'J Type Thermocouples', units='kg/h',
      cost=497, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Ceramic Fiber', units='kg/h',
      cost=5154, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Steel Jacket', units='kg/h',
      cost=22515, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Counterflow Heat Exchanger', units='kg/h',
      cost=14355, S=pilot_flowrate, CE=CEPCI_by_year[2013],n=0.77, BM=2.2)
@cost(basis='Feedstock dry flowrate', ID= 'Temperature Control and Data Logging Unit', units='kg/h',
      cost=905, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'Pulsation Dampener', units='kg/h',
      cost=3000, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'Fluid Accumulator', units='kg/h',
      cost=995, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'Burst Rupture Discs', units='kg/h',
      cost=1100, S=pilot_flowrate, CE=CEPCI_by_year[2023], n=0.77, BM=1.6)
@cost(basis='Feedstock dry flowrate', ID= 'Pressure Relief Vessel', units='kg/h',
      cost=4363, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2)
@cost(basis='Feedstock dry flowrate', ID= 'Gas Scrubber', units='kg/h',
      cost=1100, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'BPR', units='kg/h',
      cost=4900, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.6)
@cost(basis='Feedstock dry flowrate', ID= 'Primary Collection Vessel', units='kg/h',
      cost=7549, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Belt Oil Skimmer', units='kg/h',
      cost=2632, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Bag Filter', units='kg/h',
      cost=8800, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.7)
@cost(basis='Feedstock dry flowrate', ID= 'Oil Vessel', units='kg/h',
      cost=4330, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Mobile HTL system', units='kg/h',
      cost=23718, S=pilot_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
class PilotHTL(qsu.HydrothermalLiquefaction):   
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
        'Feedstock dry flowrate': 'kg/h',
        }
    
    # auxiliary_unit_names=('heat_exchanger','kodrum')

    _F_BM_default = {
        **qsu.HydrothermalLiquefaction._F_BM_default,
        # 'Feedstock Tank': 1.5,
        }
    
    # ID of the components that will be used in mass flowrate calculations
    ash_ID = 'Ash'
    water_ID = 'Water'
    
    # Product condition adjustment based on 
    gas_composition = {
        'CH4':0.050,
        'C2H6':0.032,
        'CO2':0.918
        } # [4]
    
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
                  N=4, V=None, auxiliary=False,
                  mixing_intensity=None, kW_per_m3=0,
                  wall_thickness_factor=1,
                  vessel_material='Stainless steel 316',
                  vessel_type='Horizontal',
                  CAPEX_factor=1,
                  HTL_steel_cost_factor=2.7, # so the cost matches [6]
                  **kwargs,
                  ):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        #!!! Need to compare the externally sourced HX cost and BioSTEAM default
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = qsu.HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=self.eff_T, rigorous=True)
        #!!! Probably need to add the knockout drum
        # self.kodrum = KnockOutDrum(ID=f'.{ID}_KOdrum')
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.CAPEX_factor = CAPEX_factor
        for attr, val in kwargs.items(): setattr(self, attr, val)
    

    def _run(self):        
        feedstock = self.ins[0]
        hydrochar, HTLaqueous, biocrude, offgas = self.outs
        
        #!!! Should allow users to update the yields and properties
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
        
        # hydrochar.P = self.hydrochar_P
        # HTLaqueous.P = self.HTLaqueous_P
        # biocrude.P = self.biocrude_P
        # offgas.P = self.offgas_P
        
        for stream in self.outs : stream.T = self.heat_exchanger.T
        
    @property
    def afdw_mass_in(self):
        '''Total ash-free dry mass of the feedstock.'''
        feedstock = self.ins[0]
        return feedstock.F_mass-feedstock.imass[self.ash_ID]-feedstock.imass[self.water_ID]

    @property
    def afdw_biocrude_yield(self):
        '''Biocrude product yield on the ash-free dry weight basis of the feedstock.'''
        return 0.5219

    @property
    def afdw_aqueous_yield(self):
        '''Aqueous product yield on the ash-free dry weight basis of the feedstock.'''
        return 0.2925

    @property
    def afdw_hydrochar_yield(self):
        '''Hydrochar product yield on the ash-free dry weight basis of the feedstock.'''
        return 0.01
    
    @property
    def afdw_gas_yield(self):
        '''Gas product yield on the ash-free dry weight basis of the feedstock.'''
        return 0.1756

    # @property
    # def biocrude_C_ratio(self):
    #     return (self.WWTP.AOSc*self.biocrude_C_slope + self.biocrude_C_intercept)/100 # [2]
    
    # @property
    # def biocrude_H_ratio(self):
    #     return (self.WWTP.AOSc*self.biocrude_H_slope + self.biocrude_H_intercept)/100 # [2]

    # @property
    # def biocrude_N_ratio(self):
    #     return self.biocrude_N_slope*self.WWTP.sludge_dw_protein # [2]
    
    # @property
    # def biocrude_C(self):
    #     return min(self.outs[2].F_mass*self.biocrude_C_ratio, self.WWTP.sludge_C)


    # @property
    # def HTLaqueous_C(self):
    #     return min(self.outs[1].F_vol*1000*self.HTLaqueous_C_slope*\
    #                self.WWTP.sludge_dw_protein*100/1000000/self.TOC_TC,
    #                self.WWTP.sludge_C - self.biocrude_C)

    # @property
    # def biocrude_H(self):
    #     return self.outs[2].F_mass*self.biocrude_H_ratio

    # @property
    # def biocrude_N(self):
    #     return min(self.outs[2].F_mass*self.biocrude_N_ratio, self.WWTP.sludge_N)
    
    # @property
    # def biocrude_HHV(self):
    #     return 30.74 - 8.52*self.WWTP.AOSc +\
    #            0.024*self.WWTP.sludge_dw_protein # [2]
               
    # @property
    # def energy_recovery(self):
    #     return self.biocrude_HHV*self.outs[2].imass['Biocrude']/\
    #            (self.WWTP.outs[0].F_mass -\
    #            self.WWTP.outs[0].imass['H2O'])/self.WWTP.sludge_HHV # [2]
        
    # @property
    # def offgas_C(self):
    #     carbon = sum(self.outs[3].imass[self.gas_composition]*
    #                  [cmp.i_C for cmp in self.components[self.gas_composition]])
    #     return min(carbon, self.WWTP.sludge_C - self.biocrude_C - self.HTLaqueous_C)
        
    # @property
    # def hydrochar_C_ratio(self):
    #     return min(self.hydrochar_C_slope*self.WWTP.sludge_dw_carbo, 0.65) # [2]

    # @property
    # def hydrochar_C(self):
    #     return min(self.outs[0].F_mass*self.hydrochar_C_ratio, self.WWTP.sludge_C -\
    #                self.biocrude_C - self.HTLaqueous_C - self.offgas_C)

    # @property
    # def hydrochar_P(self):
    #     return min(self.WWTP.sludge_P*self.hydrochar_P_recovery_ratio, self.outs[0].F_mass)

    # @property
    # def HTLaqueous_N(self):
    #     return self.WWTP.sludge_N - self.biocrude_N
        
    # @property
    # def HTLaqueous_P(self):
    #     return self.WWTP.sludge_P*(1 - self.hydrochar_P_recovery_ratio)

    def _design(self):
        Design = self.design_results
        Design['Feedstock dry flowrate'] = self.afdw_mass_in
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from((self.outs[1], self.outs[2], self.outs[3]))
        hx_outs0.copy_like(hx_ins0)
        hx_ins0.T = self.ins[0].T # temperature before/after HTL are similar
        hx_outs0.T = hx.T
        hx_ins0.P = hx_outs0.P = self.outs[0].P # cooling before depressurized, heating after pressurized
        # in other words, both heating and cooling are performed under relatively high pressure
        # hx_ins0.vle(T=hx_ins0.T, P=hx_ins0.P)
        # hx_outs0.vle(T=hx_outs0.T, P=hx_outs0.P)
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)

        self.P = self.ins[0].P
        # Reactor._design(self)
        # Design['Solid filter and separator weight'] = 0.2*Design['Weight']*Design['Number of reactors'] # assume stainless steel
        # # based on [6], case D design table, the purchase price of solid filter and separator to
        # # the purchase price of HTL reactor is around 0.2, therefore, assume the weight of solid filter
        # # and separator is 0.2*single HTL weight*number of HTL reactors
        # self.construction[0].quantity += Design['Solid filter and separator weight']*_lb_to_kg
        
        # self.kodrum.V = self.F_mass_out/_lb_to_kg/1225236*4230/_m3_to_gal
        # # in [6], when knockout drum influent is 1225236 lb/hr, single knockout
        # # drum volume is 4230 gal
        
        # self.kodrum.simulate()
        
    def _cost(self):
        # HydrothermalLiquefaction._cost(self)
        self.cost_items.clear() #!!! will not be needed if not inheriting from `HydrothermalLiquefaction`
        self._decorated_cost()
        
        purchase_costs = self.baseline_purchase_costs
        for item in purchase_costs.keys():
            purchase_costs[item] *= self.CAPEX_factor
        
        self.baseline_purchase_costs['Piping'] = self.baseline_purchase_cost*0.15
            
        # purchase_costs['Horizontal pressure vessel'] *= self.HTL_steel_cost_factor
        
        for aux_unit in self.auxiliary_units:
            purchase_costs = aux_unit.baseline_purchase_costs
            installed_costs = aux_unit.installed_costs
            for item in purchase_costs.keys():
                purchase_costs[item] *= self.CAPEX_factor
                installed_costs[item] *= self.CAPEX_factor

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

# @cost(basis='Feedstock dry flowrate', ID='Feedstock Tank', units='kg/h',
#       cost=4330, S=pilot_flowrate, CE=CEPCI_by_year[2011], n=0.77, BM=1.5)
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


ap_flowrate= 49.65 #kg/hr
@cost(basis='Aqueous flowrate', ID= 'Sand Filtration Unit', units='kg/h',
      cost=318, S=ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.7)
@cost(basis='Aqueous flowrate', ID= 'EC Oxidation Tank', units='kg/h',
      cost=1850, S=ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
@cost(basis='Aqueous flowrate', ID= 'Biological Treatment Tank', units='kg/h',
      cost=4330, S=ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
@cost(basis='Aqueous flowrate', ID= 'Liquid Fertilizer Storage', units='kg/h',
      cost=7549, S=ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
class AqueousFiltration(SanUnit):
    '''
    Placeholder for the aqueous filtration unit. All outs are copied from ins.
    '''
     
    _N_outs = 1
    _units= {
        'Aqueous flowrate': 'kg/h',
        }
    def _run(self):
        HTL_aqueous = self.ins[0]
        treated_aq = self.outs
        
        #treated_aq.copy_like(HTL_aqueous)
        
    def _design(self):
                aqueous = self.ins[0]
                self.design_results['Aqueous flowrate'] = aqueous.F_mass
                
biocrude_flowrate= 5.64 #kg/hr
@cost(basis='Biocrude flowrate', ID= 'Biocrude Storage Tank', units='kg/h',
      cost=7549, S=biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
@cost(basis='Biocrude flowrate', ID= 'Dewaering Tank', units='kg/h',
      cost=4330, S=biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
@cost(basis='Biocrude flowrate', ID= 'Deashing Tank', units='kg/h',
      cost=4330, S=biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
@cost(basis='Biocrude flowrate', ID= 'Fractional Distillation Column', units='kg/h',
      cost=63270, S=biocrude_flowrate, CE=CEPCI_by_year[2007],n=0.75, BM=2)
@cost(basis='Biocrude flowrate', ID= 'Heavy Fraction Tank', units='kg/h',
      cost=4330, S=biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
@cost(basis='Biocrude flowrate', ID= 'Medium Fraction Tank', units='kg/h',
      cost=4330, S=biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
@cost(basis='Biocrude flowrate', ID= 'Light Fraction Tank', units='kg/h',
      cost=4330, S=biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
class BiocrudeDeashing(SanUnit):
    '''
    Placeholder for the deashing unit.
    '''
    
    _N_outs = 2
    _units= {
        'Biocrude flowrate': 'kg/h',
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
            
    
    
# @cost(basis='Feedstock dry flowrate', ID='Feedstock Tank', units='kg/h',
#       cost=4330, S=pilot_flowrate, CE=CEPCI_by_year[2011], n=0.77, BM=1.5)
class BiocrudeDewatering(SanUnit):
    '''
    Placeholder for the dewatering unit.
    '''
    
    _N_outs = 2
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

            
class ShortcutColumn(bst.units.ShortcutColumn, qs.SanUnit):
    '''
    Similar to biosteam.units.ShortcutColumn.
    
    See Also
    --------
    `biosteam.units.ShortcutColumn <https://biosteam.readthedocs.io/en/latest/API/units/distillation.html>`_
    '''
            

# @cost(basis='Feedstock dry flowrate', ID='Feedstock Tank', units='kg/h',
#       cost=4330, S=pilot_flowrate, CE=CEPCI_by_year[2011], n=0.77, BM=1.5)
class Transportation(qsu.Copier):
    '''
    Placeholder for transportation. All outs are copied from ins.
    '''
    
    

