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

from math import ceil, log
from biosteam.units.decorators import cost
from qsdsan import SanUnit, Stream
from qsdsan.sanunits import HydrothermalLiquefaction, HXutility
# from qsdsan.sanunits._hydrothermal import KnockOutDrum
# from qsdsan.utils import auom
from exposan.biobinder import CEPCI_by_year

__all__ = (
    'PilotHTL',
    )

#!!! TO BE UPDATED THROUGHOUT
pilot_flowrate = 11.46 # kg/h
@cost(basis='Feedstock dry flowrate', ID='Feedstock Tank', units='kg/h',
      cost=4330, S=pilot_flowrate, CE=CEPCI_by_year[2023], n=0.77, BM=1.5)
class PilotHTL(HydrothermalLiquefaction):   
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
        **HydrothermalLiquefaction._F_BM_default,
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
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=self.eff_T, rigorous=True)
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
            
        # purchase_costs['Horizontal pressure vessel'] *= self.HTL_steel_cost_factor
        
        for aux_unit in self.auxiliary_units:
            purchase_costs = aux_unit.baseline_purchase_costs
            installed_costs = aux_unit.installed_costs
            for item in purchase_costs.keys():
                purchase_costs[item] *= self.CAPEX_factor
                installed_costs[item] *= self.CAPEX_factor
                
                
                