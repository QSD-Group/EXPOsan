#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>
    
    Joy Cheung

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import biosteam as bst
from warnings import warn
from math import ceil, exp, log
import numpy as np
from qsdsan.sanunits import (Mixer,
                             IsothermalCompressor,
                             SludgeThickening,
                             Copier,
                             HXutility
                             )
from qsdsan.processes import Decay
from qsdsan import SanUnit,Construction, WasteStream
from biosteam.units.design_tools import flash_vessel_design
from qsdsan.utils import ospath, data_path, load_data, price_ratio
import CoolProp.CoolProp as CP
from exposan.g2rt._sanunits import VolumeReductionCombustor
surt_su_data_path = ospath.join(data_path, 'sanunit_data/surt')

#demo for Siqi and Yuyao
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
    
    The following impact items should be pre-constructed for life cycle assessment:
    Steel, Aluminum, WoodPellet, StoneWool, StainlessSteel
    
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
        
        
#%%
htc_path = ospath.join(surt_su_data_path, '_surt_htc_reactor.csv')
@price_ratio()
class HydrothermalCarbonization(SanUnit):
    '''
    Hydrothermal carbonization unit that performs organic oxidation to CO2 by supercritical water oxidation.
    
    The following components should be included in system thermo object for simulation:
    H2O, OtherSS, N2O, NH3, CO2, O2, sCOD, xCOD, Tissue, N2, O2

    The following impact items should be pre-constructed for life cycle assessment:
    Aluminum, AluminumCasting, StainlessSteel, GlassFiber, Ceramic, Pump, Silicon, StainlessSteelMachining,
    Switch
    
    #TODO: add resource recovery into the unit and use self.hydrochar_selling_price as revenue
    If allowing resource recovery, this unit can produce hydrochar that can be 
    sold at a market price.

    Parameters
    ----------
    ins : Iterable(stream)
        homogenized human excreta, acetic acid catalysis
    outs : Iterable(stream)
        gas product mixture, process liquid, hydrochar, CH4

    References
    ----------
    The empirical model was developed using experimental data reported by:
    [1] Spitzer et al. Using hydrothermal carbonization for sustainable treatment 
    and reuse of human excreta. https://doi.org/10.1016/j.jclepro.2018.09.126
    
    [2] HTClean Technology Summary for BMGF in 2019

    '''
    _N_ins = 2
    _N_outs = 3
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
                
        data = load_data(path=htc_path)
        
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])

    def _run(self):
        mixed_waste, catalysis = self.ins
        gas, hydrochar_mixture, CH4 = self.outs
        solids = self.solids
        mixed_waste.phase = catalysis.phase = 'l'
        gas.phase = CH4.phase = 'g'
        mixture = WasteStream()
        mixture.mix_from(self.ins)
        hydrochar_mixture.mix_from(self.ins)
        catalysis.ivol['C2H4O2'] = mixed_waste.F_vol * self.catalyst_dosage /1000 #m3/hour
        severity_factor = 50*self.HTC_hold_time**0.2*exp(-3500/self.HTC_hold_temperature) #t, reaction time in [s], T, temperature [K]
        hydrochar_yield = 0.78952-0.74463*severity_factor #dry mass of hydrochar/dry mass of raw excreta based on severity factor.
        #C mass distribution in multiphase product
        H_hydrochar = 1.7474-.8712*severity_factor # molar ratio of N:H in hydrochar organic matters
        N_hydrochar = .0941-.1609*severity_factor # molar ratio of N:C in hydrochar organic matters
        COD_hydrochar = 1/((2+(H_hydrochar-N_hydrochar*3)*0.5)/2*32/12) # g carbon per g COD
        s_C = .7816-.3738*severity_factor
        l_C = -0.0112 + 2.907e-4*self.HTC_hold_temperature + 1.738e-5*self.HTC_hold_time
        g_C = 1- s_C - l_C - self.HTC_CH4_emission_factor
        hydrochar_mixture.imass['xCOD'] = (s_C * mixture.COD * mixture.F_vol/1000 * 
                                           self.carbon_COD_ratio/COD_hydrochar) #kg/hr
        hydrochar_mixture.imass['sCOD'] = (mixture.COD * mixture.F_vol/1000 * self.carbon_COD_ratio / COD_hydrochar * l_C) #kg/hr
        gas.imass['CO2'] = mixture.COD * mixture.F_vol/1000 * self.carbon_COD_ratio * g_C /12*44 #kg/hr
        #N mass distribution in multiphase product
        hydrochar_mixture.imass['NonNH3'] = N_hydrochar*hydrochar_mixture.imass['xCOD']*COD_hydrochar/12*14
        N_aq = -1.0992*severity_factor**2 + .5678*severity_factor + .1115 # molar ratio N:C in liquid phase
        hydrochar_mixture.imass['NH3'] = N_aq * hydrochar_mixture.imass['sCOD']*COD_hydrochar/12*14
        hydrochar_mixture.imass['Tissue'] = 0
        mc = mixture.imass['H2O']/mixture.F_mass
        if hydrochar_yield * mixed_waste.F_mass*(1-mc) > hydrochar_mixture.imass[solids]:
            hydrochar_mixture.imass['WoodAsh'] = (hydrochar_yield * mixed_waste.F_mass*(1-mc) - 
                                                  hydrochar_mixture.imass[solids])
        #CH4 emission factor
        CH4.imass['CH4']= (mixture.COD * mixture.F_vol/1000 * self.carbon_COD_ratio)/12*16 *self.HTC_CH4_emission_factor

    def _init_lca(self):
        self.construction = [
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            Construction('aluminum',linked_unit=self, 
                         item='Aluminum', 
                         quantity_unit='kg'),
            Construction('aluminum_casting',linked_unit=self, 
                         item='AluminumCasting', 
                         quantity_unit='kg'),
            Construction('glass_fiber',linked_unit=self, 
                         item='GlassFiber', 
                         quantity_unit='kg'),
            Construction('ceramic',linked_unit=self, 
                         item='Ceramic', 
                         quantity_unit='kg'),
            Construction('pump',linked_unit=self, 
                         item='Pump', 
                         quantity_unit='ea'),
            Construction('silicon',linked_unit=self, 
                         item='Silicon', 
                         quantity_unit='kg'),
            Construction('stainless_steel_machining',linked_unit=self, 
                         item='StainlessSteelMachining', 
                         quantity_unit='kg'),
            Construction('switch',linked_unit=self, 
                         item='Switch', 
                         quantity_unit='kg'),
            ]
        
    def _design(self):
        design=self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel
        design['Aluminum'] = constr[1].quantity = self.aluminum
        design['AluminumCasting'] = constr[2].quantity = self.aluminum_casting
        design['GlassFiber'] = constr[3].quantity = self.glass_fiber
        design['Ceramic'] = constr[4].quantity = self.ceramic
        design['Pump'] = constr[5].quantity = self.pump
        design['Silicon'] = constr[6].quantity = self.silicon
        design['StainlessSteelMachining'] = constr[7].quantity = self.stainless_steel_machining
        design['Switch'] = constr[8].quantity = self.switch
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Reactor'] = self.HTC_reactor_cost
        C['Preheater'] = self.HTC_preheater_cost 
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Reactor'] +
                               C['Preheater'])
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio
        mixed_waste, catalysis = self.ins
        gas, hydrochar_mixture, CH4 = self.outs
        mixed_waste.phase = catalysis.phase = 'l'
        gas.phase = CH4.phase = 'g'
        initial_enthalpy_flow = mixed_waste.H + catalysis.H #kJ/hour
        gas.T = hydrochar_mixture.T = CH4.T = self.HTC_hold_temperature
        final_enthalpy_flow = gas.H + hydrochar_mixture.H + CH4.H #kJ/hour
        self.power_input = ((final_enthalpy_flow-initial_enthalpy_flow)/3600*
                       (1+self.HTC_hold_phase_duty_cycle)/self.HTC_heating_efficiency) #kW
        self.power_recovered = self.power_input*self.HTC_heat_recovery #kW
        
        self.power_utility(self.power_input) if self.power_input>=0 else self.power_utility(0)  # kW
        self.power_utility.production = self.power_recovered #kW
        
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self.__calc_replacement_cost()+
                         self._calc_maintenance_labor_cost()) #USD/hr

    def _calc_replacement_cost(self):
        mixed_waste, catalysis = self.ins
        catalysis_cost = catalysis.ivol['C2H4O2'] * self.acetic_acid_cost *1000 #USD/hour
        return catalysis_cost* self.price_ratio # USD/hr
            
    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.HTC_operation_maintenance * self.wages)
        return maintenance_labor_cost / (365*24) #USD/hr 
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return self.power_input-self.power_recovered #kW
#%%
md_path = ospath.join(surt_su_data_path, '_surt_membrane_distillation.csv')
@price_ratio()
class MembraneDistillation(SanUnit):
    '''
    Membrane distillation unit that purifies water, reject organics, NH3, phosphate, 
    and other soluble salts. if energy recovery is True, heat recovered in the chiller will
    be reflected as electricity production.
    
    The following components should be included in system thermo object for simulation:
    H2O, OtherSS, N2O, NH3, CO2, O2, sCOD, xCOD, Tissue, N2, O2

    The following impact items should be pre-constructed for life cycle assessment:
    GFRPlastic, Steel, NylonGlassFilled, CastingBrass, Pump
    

    Parameters
    ----------
    ins : Iterable(stream)
        Mixed waste
    outs : Iterable(stream)
        Permeate, Brine
    if_energy_recovery: bool
        if recover the heat energy from the hot influent stream, if False, will not recover the heat

    References
    ----------
    [1] Pikaar, I.; Guest, J.; Ganigue, R.; Jensen, P.; Rabaey, K.; Seviour, T.;
        Trimmer, J.; van der Kolk, O.; Vaneeckhaute, C.; Verstraete, W.; Resource
        Recovery from Water: Principles and Applicaiton. IWA 2022.
    [2] C.J. Davey, P. Liu, F. Kamranvand, L. Williams,Y. Jiang, A. Parker, S. Tyrrel,
    E.J. McAdam. Membrane distillation for concentrated blackwater: Influence of configuration 
    (air gap, direct contact, vacuum) on selectivity and water productivity. 
    Separation and Purification Technology. 2021.
    [3] Khan, E.U.; Nordberg, Å.; Malmros, P. Waste Heat Driven Integrated 
    Membrane Distillation for Concentrating Nutrients and Process Water Recovery 
    at a Thermophilic Biogas Plant. Sustainability 2022

    '''
    _N_ins = 1
    _N_outs = 4
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 influent_temperature = 333.15, #K
                 if_energy_recovery = True,
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.influent_temperature = influent_temperature
        self.if_energy_recovery = if_energy_recovery
        data = load_data(path=md_path)
        
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])
        
    def _run(self):
        mixed_waste, = self.ins
        mixed_waste.T = self.influent_temperature
        md_permeate, md_retentate = self.outs
        md_retentate.copy_like(mixed_waste)
        solubles, solids = self.solubles, self.solids
        md_permeate.phase = md_retentate.phase = 'l'
        self.pKa = pKa = -log(exp(52.22/8.3145*1000*(1/298 - 1/self.ins[0].T))*10**(-9.252), 10)
        # calculation of pKa under different T follows van't Hoff relationship [1] page 292
        ammonia_to_ammonium = 10**(-pKa)/10**(-self.feed_pH)
        NH3_removal = 1/(ammonia_to_ammonium+1)
        water_flux = 0.2685*self.influent_temperature - 7.0297 #kg/m2/h
        # empirical relationship extracted from AGMD in [2] Fig.5
        self.water_recovery = (self.concentration_factor-1)/self.concentration_factor
        self.membrane_area = mixed_waste.F_mass/water_flux/self.effective_membrane_fraction #m2
        md_permeate.imass['H2O'] = mixed_waste.imass['H2O'] * self.water_recovery
        md_permeate.imass['NH3'] = mixed_waste.imass['NH3'] * (1-NH3_removal)
        md_permeate.imass['sCOD'] = mixed_waste.imass['sCOD'] * (1-self.COD_removal)
        md_permeate.imass[solubles] = mixed_waste.imass[solubles] * (1-self.TDS_removal)
        md_retentate.mass = mixed_waste.mass - md_permeate.mass

    def _init_lca(self):
        self.construction = [
            Construction(item='GFRPlastic', linked_unit=self, quantity_unit='kg'),
            Construction(item='Steel', linked_unit=self, quantity_unit='kg'),
            Construction(item='NylonGlassFilled', linked_unit=self, quantity_unit='kg'),
            Construction(item='CastingBrass', linked_unit=self, quantity_unit='kg'),
            Construction(item='Pump', linked_unit=self, quantity_unit='ea'),
            ]
        
    def _design(self):
        design=self.design_results
        design = self.design_results
        constr = self.construction
        design['GFRPlastic'] = constr[0].quantity = self.Plastic_weight
        design['Steel'] = constr[1].quantity = self.Steel_weight
        design['NylonGlassFilled'] = constr[2].quantity = self.Nylon_weight
        design['CastingBrass'] = constr[3].quantity = self.Brass_weight
        design['Pump'] = constr[4].quantity = 3 #inf, eff, flushing
        self.add_construction(add_cost=False)
        
    def _cost(self):#cost for MD is assumed based on [3]
        C = self.baseline_purchase_costs
        C['Membrane module'] = (self.membrane_area * self.membrane_cost + 
                                self.membrane_housing)
        C['Pumps'] = (self.flush_pump_cost + self.feed_pump_cost*2) #for both in and out
        C['Blower'] = self.air_scrubbing_blower_cost
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Membrane module'] +
                                                         C['Pumps']+
                                                         C['Blower'])
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio
        power_demand = (self.ins[0].F_mass+ self.outs[0].F_mass)* self.pump_energy/1000 #kW        
        self.power_utility(power_demand) # kW
        
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self.__calc_replacement_cost()+
                         self._calc_maintenance_labor_cost()) #USD/hr

    def _calc_replacement_cost(self):
        C = self.baseline_purchase_costs
        pump_replacement_cost = C['Pumps']/self.pump_lifetime/(365*24) #USD/hr
        membrane_replacement_cost = self.membrane_area * self.membrane_cost/self.MD_membrane_lifetime/(365*24) #USD/hr
        return (pump_replacement_cost+membrane_replacement_cost)* self.price_ratio # USD/hr
            
    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.membrane_distillation_maintenance * self.wages)
        return maintenance_labor_cost / (365*24) #USD/hr 
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.ins[0].F_mass+ self.outs[0].F_mass)* self.pump_energy/1000 #kW
#%%
hxu_path = ospath.join(surt_su_data_path, '_surt_heat_exchanger.csv')
@price_ratio()
class HeatExchanger(HXutility):
    '''
    Create a heat exchanger that changes temperature of the
    outlet stream using a heat utility. Cost are 
    specific for use in single unit toilet.
    
    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, HeatingUnit, Copper
    
    Parameters
    ----------
    ins : 
        Inlet.
    outs : 
        Outlet.
    T=None : float
        Temperature of outlet stream [K].
    V=None : float
        Vapor fraction of outlet stream.
    if_energy_recovery: bool
        If True, heat exchanger recovers a fraction of heat when cooling.
    heat_recovery_efficiency: float
        The fraction of heat recovered as electricity.
    rigorous=False : bool
        If true, calculate vapor liquid equilibrium
    U=None : float, optional
        Enforced overall heat transfer coefficent [kW/m^2/K].
    heat_exchanger_type : str, optional
        Heat exchanger type. Defaults to "Floating head".
    N_shells : int, optional
        Number of shells. Defaults to 2.
    ft : float, optional
        User imposed correction factor.
    heat_only : bool, optional
        If True, heat exchanger can only heat.
    cool_only : bool, optional
        If True, heat exchanger can only cool.
    heat_transfer_efficiency : bool, optional
        User enforced heat transfer efficiency. A value less than 1
        means that a fraction of heat transfered is lost to the environment.
        Defaults to the heat transfer efficiency of the utility agent.

    Notes
    -----
    Must specify either `T` or `V` when creating a HXutility object.

    '''
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None,
                 include_construction=True, if_energy_recovery = True, 
                 heat_recovery_efficiency = 0.4965,
                 heat_transfer_efficiency=0.7,
                 *, T=None, V=None, rigorous=False, U=None, H=None,
                 heat_exchanger_type="Floating head",
                 material="Carbon steel/carbon steel",
                 N_shells=2, ft=None, heat_only=None, cool_only=None,
                  **kwargs):
            HXutility.__init__(self,ID=ID, ins=ins, outs=outs, thermo=thermo,
                               init_with=init_with, F_BM_default=F_BM_default,
                               include_construction=include_construction,
                               T=T, V=V, rigorous=rigorous, U=U, H=H,
                               heat_exchanger_type=heat_exchanger_type,
                               material=material, N_shells=N_shells,
                               heat_only=cool_only,cool_only=cool_only,
                               heat_transfer_efficiency=heat_transfer_efficiency
                               )
            self.if_energy_recovery = if_energy_recovery
            self.heat_recovery_efficiency = heat_recovery_efficiency

            data = load_data(path=hxu_path)
            
            for para in data.index:
                value = float(data.loc[para]['expected'])
                setattr(self, para, value)
            del data
    
            for attr, value in kwargs.items():
                setattr(self, attr, value)
            cmps = self.components
            self.solids = tuple((cmp.ID for cmp in cmps.solids))
            self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])

    def _init_lca(self):
        self.construction = [
            Construction("stainless_steel", linked_unit=self,
                         item = "StainlessSteel", 
                         quantity_unit= "kg"),
            Construction("copper", linked_unit=self,
                         item = "Copper", 
                         quantity_unit= "kg"),
            Construction('heating_unit',linked_unit=self, 
                         item='HeatingUnit', 
                         quantity_unit='ea'),
            ]
        
    def _design(self):
        design=self.design_results
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.chiller_stainless_steel_weight
        design['Copper'] = constr[1].quantity = self.preheater_copper_weight
        design['HeatingUnit'] = constr[2].quantity = 1
        self.add_construction(add_cost=False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        if self.Hnet<0: #chiller
            C['Chiller'] = self.water_chiller_cost
            refer_cost = C['Chiller']
            power_demand = -self.Hnet/self.coefficient_of_performance/3600 #kW
            power_production = power_demand*self.heat_recovery_efficiency #kW
        else: #Heater
            C['Heaters'] = self.water_heater_cost + self.preheater_cost
            refer_cost = C['Heaters']
            power_demand = self.Hnet/self.heat_transfer_efficiency/3600 #kW
            power_production = 0
        C['Tank'] = self.liquids_tank_cost
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(refer_cost+C['Tank'])
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio     
        self.power_utility(power_demand) # kW
        if self.if_energy_recovery:
            self.power_utility.production = power_production
        
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
            
    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.Hnet<0:
            maintenance_labor_cost= (self.chiller_maintenance * self.wages)
        else:
            maintenance_labor_cost = (self.heaters_maintenance * self.wages)
        return maintenance_labor_cost / (365*24) #USD/hr 
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        if self.Hnet<0:
            return (-self.Hnet/self.coefficient_of_performance/3600) #kW
        else:
            return (self.Hnet/self.heat_transfer_efficiency/3600) #kW
