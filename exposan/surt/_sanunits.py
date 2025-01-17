#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import biosteam as bst
from warnings import warn
from math import ceil
import numpy as np
from qsdsan.sanunits._abstract import Mixer
from qsdsan.sanunits import IsothermalCompressor
from qsdsan.processes import Decay
from qsdsan import SanUnit,Construction, WasteStream
from qsdsan.sanunits import SludgeThickening, Copier
from biosteam.units.design_tools import flash_vessel_design
from qsdsan.utils import ospath, data_path, load_data, price_ratio
import CoolProp.CoolProp as CP
from exposan.g2rt._sanunits import VolumeReductionCombustor
surt_su_data_path = ospath.join(data_path, 'sanunit_data/surt')


__all__ = ('BiomassCombustion',
           )

#%%
biomass_combustion_path = ospath.join(surt_su_data_path, '_surt_biomass_combustion.csv')
@price_ratio()
class BiomassCombustion(VolumeReductionCombustor):
    '''
    Combust dry biomass (feces and wood pellets) to generate heat (in the form of power) 
    to offset the electrical heating demand for single unit reinvented toilet. This unit
    has automatic feeding system with pollution control.
    
    Parameters
    ----------
    ins : Iterable(WasteStream)
        Dewatered feces solid cakes, wood pellets, air
    outs : Iterable(WasteStream)
        Wood Ash, hot gas, fugitive N2O, fugitive CH4, fugitive NO, fugitive SO2
    energy_recovery_efficiency: float from 0 to 1
        The efficiency of the bioler to reuse the combustion heat to offset system
        heating demand
    biofuel_to_solids: float
        g-wood pellets/g-feces solids, the default value is 0
    if_sludge_service: bool
        if share combustion unit among 20 SURTs (i.e., 120 users), the default is False
    if_energy_recovery: bool
        if recovering combustion heat to offset system electricity, the default is True
    
    '''
    _N_ins = 3
    _N_outs = 7
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 energy_recovery_efficiency = 0.4,
                 biofuel_to_solids = 0,
                 if_sludge_service = False,
                 if_energy_recovery = True,
                 lifetime = 10,
                 **kwargs
                 ):
        VolumeReductionCombustor.__init__(self, ID=ID, ins=ins, outs=outs, 
                                          thermo=thermo, init_with=init_with, 
                                          if_sludge_service = if_sludge_service,
                                          lifetime = lifetime
                                          )
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.energy_recovery_efficiency = energy_recovery_efficiency
        self.biofuel_to_solids = biofuel_to_solids
        self.if_sludge_service = if_sludge_service
        self.if_energy_recovery = if_energy_recovery
        self.lifetime = lifetime
        
        data = load_data(path=biomass_combustion_path)
        for para in data.index:
            try:
                value = float(data.loc[para]['expected'])
            except:
                print(f"Error occurred at para={para}")
                breakpoint() 
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
    def _init_lca(self):
        self.construction = [
            Construction('Steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            Construction('Aluminum', linked_unit=self, item='Aluminum', quantity_unit='kg'),
            Construction('WoodPellet', linked_unit=self, item='WoodPellet', quantity_unit='kg'),
            Construction('StoneWool', linked_unit=self, item='StoneWool', quantity_unit='kg'),
            Construction('StainlessSteel', linked_unit=self, item='StainlessSteel', quantity_unit='kg'),
            ]
    
    def _design(self): #simulate energy in design to enture capturing system-wise utility
        solid_cakes, wood_pellets, air = self.ins
        wood_pellets.imass['WoodPellet'] = solid_cakes.F_mass * self.biofuel_to_solids
        ash, gas, CH4, N2O, NO, SO2, NH3 = self.outs
        # NO emissions
        NO.imass['NO'] = SO2.imass['SO2'] = NH3.imass['NH3'] = 0 #pollution control
        N = ceil(solid_cakes.F_mass/0.7) #assume 0.7 kg/hr combustion capacity
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.steel_weight * N
        design['Aluminum'] = constr[1].quantity = self.aluminum_weight* N
        design['WoodPellet'] = constr[2].quantity = self.lifetime*365*24* wood_pellets.imass['WoodPellet']* N
        design['StoneWool']= constr[3].quantity = self.stone_wool_weight * N
        design['StainlessSteel'] = constr[4].quantity = self.stainless_steel_weight * N
    
    
    def _cost(self):
        
        solid_cakes, wood_pellets, air = self.ins
        wood_pellets.imass['WoodPellet'] = solid_cakes.F_mass * self.biofuel_to_solids
        mixture = WasteStream()
        mixture.mix_from(self.ins)
        cmps = self.components
        thermal_energy = 0 #in kJ/hr
        for cmp in cmps:
            if not cmp.organic: #assume that only organic components release heat during combustion
                continue
            thermal_energy += mixture.imass[cmp.ID] * cmp.HHV #in kJ/hr
        heat_to_dry = mixture.imass['H2O'] * self.energy_required_to_evaporize_water #kWh/hr
        net_power_generated = (thermal_energy/3600-heat_to_dry) * self.energy_recovery_efficiency
        
        if self.if_energy_recovery:
            if net_power_generated > 0:
                self.power_utility.production = net_power_generated
                self.power_utility.consumption = 0
            else:
                self.power_utility.consumption = net_power_generated
                self.power_utility.production = 0
        N = ceil(solid_cakes.F_mass/0.7) #assume 0.7 kg/hr combustion capacity
        C = self.baseline_purchase_costs
        C["AirPreheating"] = self.air_preheating_unit_cost * self.price_ratio
        C["Combustor"] = self.material_combustion_unit_cost * self.price_ratio
        C["Meter&Boiler"] = self.metering_boiler_unit_cost * self.price_ratio
        service_factor = 0.05 if self.if_sludge_service else 1
        total_equipment = 0.
        for equipment, cost in C.items():
            C[equipment] = cost * service_factor * N
            total_equipment += cost * N
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_replacement_cost()+
                         self._calc_maintenance_labor_cost()*N) #USD/hr

            
    def _calc_replacement_cost(self):
        solid_cakes, wood_pellets, air = self.ins
        wood_pellet_cost = wood_pellets.imass['WoodPellet'] * self.wood_pellets_cost #USD/hr
        service_factor = 0.05 if self.if_sludge_service else 1
        return wood_pellet_cost* self.price_ratio * service_factor # USD/hr
            
    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.combustor_operation_maintenance * self.wages)
        service_factor = 0.05 if self.if_sludge_service else 1
        return maintenance_labor_cost / (365*24) * service_factor
        
        

        
        
        
    
    

