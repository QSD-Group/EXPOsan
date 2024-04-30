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
from qsdsan import SanUnit
from qsdsan.sanunits import HydrothermalLiquefaction
# from qsdsan.utils import auom
from exposan.biobinder import CEPCI_by_year

__all__ = (
    'PilotHTL',
    )

#!!! TO BE UPDATED THROUGHOUT
pilot_flowrate = 11.46 # kg/h
@cost(basis='Feedstock dry flowrate', ID='Feedstock Tank', units='kg/h',
      cost=4330, S=pilot_flowrate, CE=CEPCI_by_year[2011], n=0.77, BM=1.5)
class PilotHTL(HydrothermalLiquefaction):

    
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

    # def __init__(self, ID='', ins=None, outs=(), thermo=None,
    #              init_with='WasteStream',
    #              lipid_2_biocrude=0.846, # [1]
    #              protein_2_biocrude=0.445, # [1]
    #              carbo_2_biocrude=0.205, # [1]
    #              protein_2_gas=0.074, # [1]
    #              carbo_2_gas=0.418, # [1]
    #              biocrude_C_slope=-8.37, # [2]
    #              biocrude_C_intercept=68.55, # [2]
    #              biocrude_N_slope=0.133, # [2]
    #              biocrude_H_slope=-2.61, # [2]
    #              biocrude_H_intercept=8.20, # [2]
    #              HTLaqueous_C_slope=478, # [2]
    #              TOC_TC=0.764, # [3]
    #              hydrochar_C_slope=1.75, # [2]
    #              biocrude_moisture_content=0.063, # [4]
    #              hydrochar_P_recovery_ratio=0.86, # [5]
    #              gas_composition={'CH4':0.050, 'C2H6':0.032,
    #                               'CO2':0.918}, # [4]
    #              hydrochar_pre=3029.7*6894.76, # [4]
    #              HTLaqueous_pre=30*6894.76, # [4]
    #              biocrude_pre=30*6894.76, # [4]
    #              offgas_pre=30*6894.76, # [4]
    #              eff_T=60+273.15, # [4]
    #              P=None, tau=15/60, V_wf=0.45,
    #              length_to_diameter=None, diameter=6.875*_in_to_m,
    #              N=4, V=None, auxiliary=False,
    #              mixing_intensity=None, kW_per_m3=0,
    #              wall_thickness_factor=1,
    #              vessel_material='Stainless steel 316',
    #              vessel_type='Horizontal',
    #              CAPEX_factor=1,
    #              HTL_steel_cost_factor=2.7, # so the cost matches [6]
    #              mositure_adjustment_exist_in_the_system=False):
        
    #     SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
    #     self.lipid_2_biocrude = lipid_2_biocrude
    #     self.protein_2_biocrude = protein_2_biocrude
    #     self.carbo_2_biocrude = carbo_2_biocrude
    #     self.protein_2_gas = protein_2_gas
    #     self.carbo_2_gas = carbo_2_gas
    #     self.biocrude_C_slope = biocrude_C_slope
    #     self.biocrude_C_intercept = biocrude_C_intercept
    #     self.biocrude_N_slope = biocrude_N_slope
    #     self.biocrude_H_slope = biocrude_H_slope
    #     self.biocrude_H_intercept = biocrude_H_intercept
    #     self.HTLaqueous_C_slope = HTLaqueous_C_slope
    #     self.TOC_TC = TOC_TC
    #     self.hydrochar_C_slope = hydrochar_C_slope
    #     self.biocrude_moisture_content = biocrude_moisture_content
    #     self.hydrochar_P_recovery_ratio = hydrochar_P_recovery_ratio
    #     self.gas_composition = gas_composition
    #     self.hydrochar_pre = hydrochar_pre
    #     self.HTLaqueous_pre = HTLaqueous_pre
    #     self.biocrude_pre = biocrude_pre
    #     self.offgas_pre = offgas_pre
    #     hx_in = Stream(f'{ID}_hx_in')
    #     hx_out = Stream(f'{ID}_hx_out')
    #     self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=eff_T, rigorous=True)
    #     self.kodrum = KnockOutDrum(ID=f'.{ID}_KOdrum')
    #     self.P = P
    #     self.tau = tau
    #     self.V_wf = V_wf
    #     self.length_to_diameter = length_to_diameter
    #     self.diameter = diameter
    #     self.N = N
    #     self.V = V
    #     self.auxiliary = auxiliary
    #     self.mixing_intensity = mixing_intensity
    #     self.kW_per_m3 = kW_per_m3
    #     self.wall_thickness_factor = wall_thickness_factor
    #     self.vessel_material = vessel_material
    #     self.vessel_type = vessel_type
    #     self.CAPEX_factor = CAPEX_factor
    #     self.HTL_steel_cost_factor = HTL_steel_cost_factor
    #     self.mositure_adjustment_exist_in_the_system = mositure_adjustment_exist_in_the_system

    def _run(self):
        
        feedstock = self.ins[0]
        hydrochar, HTLaqueous, biocrude, offgas = self.outs
        
        
        #!!! Update so that it could be set by the users
        dewatered_sludge_afdw = feedstock.imass['Lipids'] +\
                                feedstock.imass['Proteins'] +\
                                feedstock.imass['Carbohydrates']
        # just use afdw in revised MCA model, other places use dw
        
        self.afdw_lipid_ratio = self.WWTP.sludge_afdw_lipid
        self.afdw_protein_ratio = self.WWTP.sludge_afdw_protein
        self.afdw_carbo_ratio = self.WWTP.sludge_afdw_carbo

        # the following calculations are based on revised MCA model
        hydrochar.imass['Hydrochar'] = 0.377*self.afdw_carbo_ratio*dewatered_sludge_afdw
        
        HTLaqueous.imass['HTLaqueous'] = (0.481*self.afdw_protein_ratio +\
                                          0.154*self.afdw_lipid_ratio)*\
                                          dewatered_sludge_afdw
        # HTLaqueous is TDS in aqueous phase
        # 0.377, 0.481, and 0.154 don't have uncertainties because they are calculated values
         
        gas_mass = (self.protein_2_gas*self.afdw_protein_ratio + self.carbo_2_gas*self.afdw_carbo_ratio)*\
                       dewatered_sludge_afdw
                       
        for name, ratio in self.gas_composition.items():
            offgas.imass[name] = gas_mass*ratio
            
        biocrude.imass['Biocrude'] = (self.protein_2_biocrude*self.afdw_protein_ratio +\
                                      self.lipid_2_biocrude*self.afdw_lipid_ratio +\
                                      self.carbo_2_biocrude*self.afdw_carbo_ratio)*\
                                      dewatered_sludge_afdw
        biocrude.imass['H2O'] = biocrude.imass['Biocrude']/(1 -\
                                self.biocrude_moisture_content) -\
                                biocrude.imass['Biocrude']
                                
        HTLaqueous.imass['H2O'] = feedstock.F_mass - hydrochar.F_mass -\
                                  biocrude.F_mass - gas_mass - HTLaqueous.imass['HTLaqueous']
        # assume ash (all soluble based on Jones) goes to water
        
        hydrochar.phase = 's'
        offgas.phase = 'g'
        HTLaqueous.phase = biocrude.phase = 'l'
        
        hydrochar.P = self.hydrochar_pre
        HTLaqueous.P = self.HTLaqueous_pre
        biocrude.P = self.biocrude_pre
        offgas.P = self.offgas_pre
        
        for stream in self.outs : stream.T = self.heat_exchanger.T
        
    # @property
    # def biocrude_yield(self):
    #     return self.protein_2_biocrude*self.afdw_protein_ratio +\
    #            self.lipid_2_biocrude*self.afdw_lipid_ratio +\
    #            self.carbo_2_biocrude*self.afdw_carbo_ratio

    # @property
    # def aqueous_yield(self):
    #     return 0.481*self.afdw_protein_ratio + 0.154*self.afdw_lipid_ratio
    
    # @property
    # def hydrochar_yield(self):
    #     return 0.377*self.afdw_carbo_ratio
    
    # @property
    # def gas_yield(self):
    #     return self.protein_2_gas*self.afdw_protein_ratio + self.carbo_2_gas*self.afdw_carbo_ratio

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

    # def _design(self):
        
    #     Design = self.design_results
    #     Design['Treatment capacity'] = self.ins[0].F_mass/_lb_to_kg
        
    #     hx = self.heat_exchanger
    #     hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
    #     hx_ins0.mix_from((self.outs[1], self.outs[2], self.outs[3]))
    #     hx_outs0.copy_like(hx_ins0)
    #     hx_ins0.T = self.ins[0].T # temperature before/after HTL are similar
    #     hx_outs0.T = hx.T
    #     hx_ins0.P = hx_outs0.P = self.outs[0].P # cooling before depressurized, heating after pressurized
    #     # in other words, both heating and cooling are performed under relatively high pressure
    #     # hx_ins0.vle(T=hx_ins0.T, P=hx_ins0.P)
    #     # hx_outs0.vle(T=hx_outs0.T, P=hx_outs0.P)
    #     hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)

    #     self.P = self.ins[0].P
    #     Reactor._design(self)
    #     Design['Solid filter and separator weight'] = 0.2*Design['Weight']*Design['Number of reactors'] # assume stainless steel
    #     # based on [6], case D design table, the purchase price of solid filter and separator to
    #     # the purchase price of HTL reactor is around 0.2, therefore, assume the weight of solid filter
    #     # and separator is 0.2*single HTL weight*number of HTL reactors
    #     self.construction[0].quantity += Design['Solid filter and separator weight']*_lb_to_kg
        
    #     self.kodrum.V = self.F_mass_out/_lb_to_kg/1225236*4230/_m3_to_gal
    #     # in [6], when knockout drum influent is 1225236 lb/hr, single knockout
    #     # drum volume is 4230 gal
        
    #     self.kodrum.simulate()
        
    def _cost(self):
        HydrothermalLiquefaction._cost(self)
        self._decorated_cost()