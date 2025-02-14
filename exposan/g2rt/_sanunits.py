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
g2rt_su_data_path = ospath.join(data_path, 'sanunit_data/g2rt')
surt_su_data_path = ospath.join(data_path, 'sanunit_data/surt')

__all__ = ('Excretion',
           'VolumeReductionFilterPress',
           'VRConcentrator',
           'VRdryingtunnel',
           'G2RThomogenizer',
           'VRpasteurization',
           'G2RTLiquidsTank',
           'G2RTSolidsTank',
           'G2RTControls',
           'G2RTHousing',
           'G2RTSolidsSeparation',
           'G2RTBeltSeparation',
           'G2RTUltrafiltration',
           'G2RTReverseOsmosis',
           'SURT',
           'FWMixer',
           'UFMixer',
           'VolumeReductionCombustor',
           'mSCWOGasModule',
           'mSCWOReactorModule',
           'mSCWOConcentratorModule',
           )

#%%
mscwo_gas_module_path = ospath.join(g2rt_su_data_path, '_mscwo_gas_handling_module.csv')
@price_ratio()
class mSCWOGasModule(IsothermalCompressor):
    '''
    Gas handling unit that consisted of compressors and an injection pressure vessel
    
    The following components should be included in system thermo object for simulation:
    O2, N2, CO2  #TODO: update

    The following impact items should be pre-constructed for life cycle assessment:
    Compressor_4kW, Compressor_300kW, StainlessSteel

    Parameters
    ----------
    ins : Iterable(stream)
        Air in
    outs : Iterable(stream)
        Compressed air out
    P : float
        Outlet pressure [Pa].
    eta : float
        Isothermal efficiency.
    vle : bool
        Whether to perform phase equilibrium calculations on
        the outflow. If False, the outlet will be assumed to be the same
        phase as the inlet.
    type: str
        Type of compressor : blower/centrifugal/reciprocating. If None, the type
        will be determined automatically.

    References
    ----------
    [1] YEE et al. Water oxidation non-sewered single unit toilet system. 
    https://patentimages.storage.googleapis.com/57/6a/81/72a168a92be44c/WO2023288331A1.pdf
    '''
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 P=1.55E7, eta =0.7,vle=False, compressor_type='Reciprocating', 
                 driver=None, material=None, driver_efficiency=None,include_construction= True,
                 **kwargs):
        IsothermalCompressor.__init__(self, ID=ID, ins=ins, outs=outs, 
                                      P = P, eta = eta, vle = vle, 
                                      compressor_type = compressor_type, driver= driver, 
                                      material = material, driver_efficiency = driver_efficiency 
                                      )
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.include_construction = include_construction
        data = load_data(path=mscwo_gas_module_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        
        self.ins.P = 101325
        self.outs.P = self.P
    
    def _init_lca(self):
        self.include_construction = True
        self.construction = [
            Construction('compressor_4kW', linked_unit=self, 
                         item='Compressor_4kW', 
                         quantity_unit='ea'),
            Construction('compressor_300kW', linked_unit=self, 
                         item='Compressor_300kW', 
                         quantity_unit='ea'),
            Construction("stainless_steel", linked_unit=self,
                         item = "StainlessSteel", 
                         quantity_unit= "kg"),
            ]
    
    def _run(self):
        air_in, = self.ins
        air_in.P = 101325
        air_out, = self.outs
        air_out.copy_like(air_in)
        air_out.P = self.P
        air_out.T = air_in.T
        if self.vle is True: air_out.vle(T=air_out.T, P=air_out.P)
        self.ideal_power, self.ideal_duty = self._calculate_ideal_power_and_duty()
    
    
    def _design(self):
        self._init_lca()
        air_in, = self.ins
        air_out, = self.outs
        air_in.imass['N2'] = air_out.imass['N2']
        air_in.imass['O2'] = air_out.imass['O2']
        air_out.P = self.P
        air_out.T = air_in.T
        super()._design()
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[2].quantity = self.stainless_steel_weight
        self.add_construction(add_cost = False)
        
    def _cost(self):
        super()._design()
        D = self.design_results
        C = self.baseline_purchase_costs
        C['Compressor'] = self.compressor_cost
        C['Valves'] = (self.injection_valve_cost +
                       self.dosing_valve_cost * self.dosing_valve_quantity
                       )
        C['Vessel_cost'] = self.vessel_cost
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Valves'] +
                               C['Vessel_cost'])
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio

        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
        self.power_utility(self.dosing_valve_power_demand * self.dosing_valve_daily_operation/24+
                           D['Ideal power']/self.eta/self.compressor_efficiency
                           ) # kWh/hr

    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.gas_module_maintenance * self.wages)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        super()._design()
        D = self.design_results
        return (self.dosing_valve_power_demand * self.dosing_valve_daily_operation/24+
                           D['Ideal power']/self.eta/self.compressor_efficiency) #kW

#%%
mscwo_reactor_module_path = ospath.join(g2rt_su_data_path, '_mscwo_reactor_module.csv')
@price_ratio()
class mSCWOReactorModule(SanUnit):
    '''
    Reactor unit that performs organic oxidation to CO2 by supercritical water oxidation.
    
    The following components should be included in system thermo object for simulation:
    H2O, OtherSS, N2O, NH3, CO2, O2, sCOD, xCOD, Tissue, N2, O2

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel

    Parameters
    ----------
    ins : Iterable(stream)
        compressed air (221 bar) , homogenized feces solids
    outs : Iterable(stream)
        gas product mixture, N2O, treated effluent, ash

    References
    ----------
    [1] YEE et al. Water oxidation non-sewered single unit toilet system.
    
    https://patentimages.storage.googleapis.com/57/6a/81/72a168a92be44c/WO2023288331A1.pdf
    '''
    _N_ins = 2
    _N_outs = 4
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', material_replacement_cost = None,
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        

        data = load_data(path=mscwo_reactor_module_path)
        for para in data.index:
            if material_replacement_cost is not None and para == "material_replacement_cost":
                self.material_replacement_cost = material_replacement_cost
                continue  # Skip the material_replacement_cost index if it's provided
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])

    def _run(self):
        compressed_air, feces = self.ins
        compressed_air.P = self.operating_pressure * 1e5
        compressed_air_vol = feces.F_vol * (self.reactor_volume-self.feces_batch_volume)/self.feces_batch_volume
        compressed_air.ivol['N2'] = compressed_air_vol * 0.79
        compressed_air.ivol['O2'] = compressed_air_vol* 0.21
        initial_enthalpy_flow = compressed_air.H + feces.H #kJ/hr
        
        gas_product, N2O, liquid_effluent, ash = self.outs
        gas_product.P = N2O.P = liquid_effluent.P = ash.P= 101325 #release into atmosphere
        N2O.phase = 'g'
        ash.phase = 's'
        N2O.imass['N2O'] = (feces.TN * feces.F_vol /1000/14*44
                            ) * self.N2O_nitrogen_fraction #kg/hr
        
        compressed_air.T = feces.T = self.operating_temperature
        final_enthalpy_flow = (compressed_air.H + 
                               CP.PropsSI('H', 'T', feces.T, 'P', self.operating_pressure * 1e5, 'Water')/1000*
                               feces.imass['H2O']) #kJ/hr
        #use CoolProp package for more accurate enthalpy for supercritical water
        mc = feces.imass['H2O']/ feces.F_mass
        chemical_energy = (feces.F_mass - feces.imass['H2O']) * self.HHV_feces_solids * 1000 #kJ/hr
        required_energy = final_enthalpy_flow - initial_enthalpy_flow #kJ/hr
        recovered_energy = chemical_energy * self.energy_recovery_efficiency * self.carbon_conversion_efficiency #kJ/hr
        self.power_input = (required_energy/self.heating_energy_efficiency - recovered_energy)/3600 #kW

        # self.power_input = required_energy/self.heating_energy_efficiency/3600 #kW
        # self.heat_output = recovered_energy #kJ/hr
        
        liquid_effluent.copy_like(feces)
        liquid_effluent.imass['NH3'] = (feces.TN * feces.F_vol /1000
                                        ) * self.ammonium_nitrogen_fraction #kg/hr
        liquid_effluent.imass['NonNH3'] = (feces.TN * feces.F_vol /1000
                                        ) * (1-self.ammonium_nitrogen_fraction-
                                             self.N2O_nitrogen_fraction) #kg/hr, 
        gas_product.copy_like(compressed_air)
        gas_product.imass['CO2'] = (feces.imass['sCOD'] + 
                                    feces.imass['xCOD'])* self.carbon_conversion_efficiency* self.carbon_COD_ratio #kg/hr
        gas_product.imol['O2'] -= gas_product.imol['CO2']
        ash.imass['xCOD'] = liquid_effluent.imass['xCOD'] * (1 - self.carbon_conversion_efficiency)
        ash.imass['WoodAsh'] = (feces.F_mass * (1-mc) * self.feces_ash_content) # kg ash /hr
        liquid_effluent.imass['sCOD'] *= (1-self.carbon_conversion_efficiency)
        liquid_effluent.imass['xCOD'] = 0
        liquid_effluent.imass['Tissue'] = 0
        gas_product.P = N2O.P = liquid_effluent.P = ash.P= 101325 #release into atmosphere
        gas_product.T = N2O.T = liquid_effluent.T = ash.T = self.operating_temperature

    def _init_lca(self):
        self.construction = [
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            ]
        
    def _design(self):
        design=self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Reactor'] = self.reactor_system_cost
        C['Injector'] = self.injector_cost 
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Reactor'] +
                               C['Injector'])
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio
        
        self.power_utility(self.power_input) if self.power_input>=0 else self.power_utility(0)  # kW
        # self.power_utility(self.power_input)    
        # self.heat_utilities(-self.heat_output)
        
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr

    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.reactor_maintenance * self.wages)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return self.power_input if self.power_input>=0 else 0  #kW

#%%
mscwo_concentrator_module_path = ospath.join(g2rt_su_data_path, '_mscwo_concentrator_module.csv')
@price_ratio()
class mSCWOConcentratorModule(SanUnit):
    '''
    Concentrator unit that evaporizes water by using heat from mSCWO effluent and external heat.
    
    The following components should be included in system thermo object for simulation:
    H2O, OtherSS, N2O, NH3, CO2, O2, sCOD, xCOD, Tissue

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, HeatingUnit, ElectricMotor, Pump, Fan, Polyethylene

    Parameters
    ----------
    ins : Iterable(stream)
        RO reject liquid, mscwo effluent liquid, mscwo effluent ash
    outs : Iterable(stream)
        Condensed effluent, recirculation water, fugitive N2O, fugitive CH4, fugitive NH3, water vapor

    References
    ----------
    [1] YEE et al. Water oxidation non-sewered single unit toilet system.
    
    https://patentimages.storage.googleapis.com/57/6a/81/72a168a92be44c/WO2023288331A1.pdf
    '''
    _N_ins = 3
    _N_outs = 6
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
                
        data = load_data(path=mscwo_concentrator_module_path)
        
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids)) + ('OtherSS',) #salts oversaturated and become solids
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])
    
    def _run(self):
        ro_waste_in, mscwo_liquid_in, mscwo_ash_in = self.ins
        # waste_out, N2O, CH4 = self.outs
        waste_out, excess_water, N2O, CH4, NH3_gas, water_vapor  = self.outs
        N2O.phase = CH4.phase = NH3_gas.phase = water_vapor.phase = 'g'
        recovered_energy_from_mscwo = (CP.PropsSI('H', 'T', mscwo_liquid_in.T, 'P', mscwo_liquid_in.P, 'Water')/1000*\
        mscwo_liquid_in.imass['H2O']) * self.heat_exchange_efficiency #kJ/hr
        waste_out.mix_from((self.ins[0],self.ins[1],self.ins[2]))
        
        solubles, solids = self.solubles, self.solids
        #Calculate water and solids in the condensed effluent
        mc_in = (ro_waste_in.imass['H2O'] + 
                 mscwo_liquid_in.imass['H2O'] + 
                 mscwo_ash_in.imass['H2O']) / (ro_waste_in.F_mass + 
                                                  mscwo_liquid_in. F_mass+
                                                  mscwo_ash_in.F_mass) # fraction
        mc_out = self.moisture_content_out/100 #convert to fraction
        if mc_in < mc_out*0.999:
            mc_out = mc_in
        TS_in = ro_waste_in.imass[solids].sum() + mscwo_liquid_in.imass[solids].sum() + mscwo_ash_in.imass[solids].sum() # kg TS dry/hr
        
        waste_out.imass['H2O'] = TS_in/(1-mc_out)*mc_out #kg water/hr
        
        #Calculate N2O and CH4 emissions
        CH4.imass['CH4'] = drying_CH4_to_air = \
            self.drying_CH4_emissions * self.carbon_COD_ratio * \
            (ro_waste_in.COD * ro_waste_in.F_vol + 
             mscwo_liquid_in.COD * mscwo_liquid_in.F_vol
             ) / 1000 # kg CH4 /hr
        drying_NH3_to_air = self.drying_NH3_emissions * (ro_waste_in.imass['NH3'] +
                                                         mscwo_liquid_in.imass['NH3']
                                                         ) # kg NH3 /hr
        
        N2O.imass['N2O'] = drying_NH3_to_air * self.NH3_to_N2O # kg N2O /hr
        waste_out.imass['NH3'] = ro_waste_in.imass['NH3'] + mscwo_liquid_in.imass['NH3'] - drying_NH3_to_air

        NH3_gas.imass['NH3'] = drying_NH3_to_air # kg NH3 /hr
        excess_water.imass['H2O'] = self.excess_water_recirculation * (ro_waste_in.imass['H2O'] + 
                                                                       mscwo_liquid_in.imass['H2O'] - 
                                                                       waste_out.imass['H2O'])
        excess_water.imass[solubles] = self.excess_water_recirculation * (ro_waste_in.imass[solubles]+
                                                                          mscwo_liquid_in.imass[solubles])
        water_vapor.imass['H2O'] = (ro_waste_in.imass['H2O'] + 
                                    mscwo_liquid_in.imass['H2O'] - 
                                    waste_out.imass['H2O'] - 
                                    excess_water.imass['H2O']) #kg H2O/hr
        # Store the calculated value of water_vapor.imass['H2O'] for use in the _cost function
        
        
        self.required_energy_input = water_vapor.imass['H2O'] * self.energy_required_to_evaporize_water - recovered_energy_from_mscwo/3600
        #Calculate COD
        drying_CO2_to_air = (self.drying_CO2_emissions * self.carbon_COD_ratio
                             * (ro_waste_in.COD * ro_waste_in.F_vol + 
                              mscwo_liquid_in.COD * mscwo_liquid_in.F_vol
                              ) / 1000) # kg CO2 /hr
        # 44/12/16 are the molecular weights of CO2, C, and CH4, respectively
        waste_out.imass['sCOD'] =  (waste_out.imass['sCOD'] + mscwo_liquid_in.imass['sCOD']) \
            -((drying_CO2_to_air/44*12+drying_CH4_to_air/16*12) / self.carbon_COD_ratio + 
              excess_water.imass['sCOD'])
    
    def _init_lca(self):
        self.construction = [
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            Construction('heating_unit',linked_unit=self, 
                         item='HeatingUnit', 
                         quantity_unit='ea'),
            Construction('electric_motor',linked_unit=self, 
                         item='ElectricMotor', 
                         quantity_unit='kg'),
            Construction('pump',linked_unit=self, 
                         item='Pump',
                         quantity_unit='ea'),
            Construction('fan',linked_unit=self, 
                         item='Fan',
                         quantity_unit='kg'),
            Construction('polyethylene',linked_unit=self,
                         item='Polyethylene', 
                         quantity_unit='kg'),
            ]
        
    def _design(self):
        design=self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        design['HeatingUnit'] = constr[1].quantity = 1
        design['ElectricMotor'] = constr[2].quantity = self.motor_weight
        design['Pump'] = constr[3].quantity = 1
        design['Fan'] = constr[4].quantity = self.fan_quantity*0.2 #one fan is 0.2 kg
        design['Polyethylene'] = constr[5].quantity = self.polyethylene_weight
        
    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        #C['Stainless steel'] = self.stainless_steel_cost * D['StainlessSteel']
        C['HeatingUnit'] = self.heating_unit_purchase_cost * D['HeatingUnit']
        C['ElectricMotor'] = self.motor_purchase_cost * D['ElectricMotor']
        C['Pump'] = self.pump_purchase_cost * D['Pump']
        #C['Polyethylene'] = self.polyethylene_cost * D['Polyethylene']
        C['Fan'] = self.fan_purchase_cost * self.fan_quantity
        C['Disc'] = self.disc_purchase_cost
        C['PhaseSeparator'] = self.phase_separator_cost
        C['Others'] = (self.enclosure_purchase_cost + 
                       self.open_concentrator_vessel_purchase_cost+
                       self.axle_purchase_cost+
                       self.thermistor_cost                       
                       )
        C["Misc.parts"] = self.miscellaneous_cost_ratio*(C['HeatingUnit'] +
                               C['ElectricMotor']+
                               C['Pump'] +
                               C['Fan'] +
                               C['Disc'] +
                               C['PhaseSeparator'] +
                               C['Others']
                               ) 
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio
        
        self.power_utility(self.pump_power_demand * self.pump_daily_operation/24+
                           self.required_energy_input) # kW
        
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
   
    def _calc_maintenance_labor_cost(self): #USD/hr
       maintenance_labor_cost= (self.concentrator_maintenance * self.wages) #USD/yr
       return maintenance_labor_cost / (365*24)
   
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.pump_power_demand * self.pump_daily_operation/24 + 
                self.required_energy_input) #kW

#%%
vr_combustor_path = ospath.join(g2rt_su_data_path, '_vr_combustor.csv')
@price_ratio()
class VolumeReductionCombustor(SanUnit):
    '''
    Combustor unit that uses wood pellets biofuels to burn dry feces solid cakes.
    There is no energy or heat recovery in this unit.
    
    The following components should be included in system thermo object for simulation:
    H2O, N, K, P OtherSS, N2O, CH4, NO, SO2, NH3

    The following impact items should be pre-constructed for life cycle assessment:
    Steel, Aluminum

    Parameters
    ----------
    ins : Iterable(WasteStream)
        Dewatered feces solid cakes, wood pellets, air
    outs : Iterable(WasteStream)
        Wood Ash, hot gas, fugitive N2O, fugitive CH4, fugitive NO, fugitive SO2
    if_sludge_service: bool
        If share combustor unit among multiple volume reduction toilets
        (assume 120 users per combustor unit,
        or 20 volume reduction toilets serving a population of 6 users per toilet).
    lifetime : 10 years
        This is used to estimate the amount of wood pellet used for TEA/LCA. Default to 10 years.

    References
    ----------
    [1] YEE et al. VOLUME REDUCTION SOLIDS TREATMENT SYSTEM. 
    https://patents.google.com/patent/WO2023288327A1/en?oq=WO2023288327A1
    '''
    _N_ins = 3
    _N_outs = 7
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_sludge_service = True, lifetime = 10,ppl =6, CH4_emission_factor=None,
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1, lifetime = lifetime)
        self.if_sludge_service = if_sludge_service
        self.lifetime = lifetime
        self.ppl = ppl
        data = load_data(path=vr_combustor_path)
        for para in data.index:
            if CH4_emission_factor is not None and para == "CH4_emission_factor":
                self.CH4_emission_factor = CH4_emission_factor
                continue  # Skip the material_replacement_cost index if it's provided
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _run(self):
        solid_cakes, wood_pellets, air = self.ins
        ash, gas, CH4, N2O, NO, SO2, NH3 = self.outs
        cmps = self.components
        ash.copy_like(self.ins[0])
        wood_pellets.imass['WoodPellet'] = solid_cakes.F_mass * self.biofuel_to_solids
        self.wood_pellets_kgphr = wood_pellets.imass['WoodPellet']
        wood_pellets.imass['H2O'] = wood_pellets.imass['WoodPellet'] * self.wood_moisture_content
        gas.phase = N2O.phase = CH4.phase = NO.phase = SO2.phase = NH3.phase = 'g'
        gas.T = self.combustion_temperature
        # Moisture content
        mc = solid_cakes.imass['H2O']/solid_cakes.F_mass if solid_cakes.F_mass!=0 else 0
        
        if mc > 0.351: # allow a small error
            warn(f'Moisture content of the influent is {mc:.1%}, '
                'larger than the maximum allowed level of 35%.')

        # ash.empty()
        
        ash_prcd = (solid_cakes.F_mass * (1-mc) * self.feces_ash_content + 
               wood_pellets.imass['WoodPellet']*(1-self.wood_moisture_content)*self.wood_ash_content) # kg ash /hr
        
        NPKCaMg = ('NH3','NonNH3','N', 'P', 'K','Ca','Mg')
        for cmp in cmps:
            ash.imass[cmp.ID] = 0 if cmp.ID not in NPKCaMg else ash.imass[cmp.ID]
        for element in NPKCaMg:
            if element == 'NH3':
                ash.imass[element] *= (1 - getattr(self, 'combustion_N_loss'))
            elif element in ('Ca','Mg','NonNH3'):
                continue
            else:
                ash.imass[element] *= (1 - getattr(self, f'combustion_{element}_loss'))
        ash.imass['WoodAsh'] = (ash_prcd - ash.imass[NPKCaMg].sum())
        ash.imass['H2O'] = 0.025 * ash.F_mass # kg H2O / hr with 2.5% moisture content
        #CH4 emissions
        CH4.imass['CH4'] = solid_cakes.COD * self.carbon_COD_ratio * solid_cakes.F_vol/1e3* self.fugitive_CH4_emission_factor #kg/hr
        # N2O emissions
        N2O.imass['N2O'] = solid_cakes.imass['N'] * self.N2O_emission_factor/28*44
        # NO emissions
        NO.imass['NO'] = solid_cakes.imass['N'] * self.NOx_emission_factor/14*30
        #SO2 emissions
        SO2.imass['SO2'] = solid_cakes.F_mass * self.SOx_emission_factor
        #NH3 emissions
        NH3.imass['NH3'] = solid_cakes.imass['N'] * self.NH3_emission_factor/14*17
        #gas
        gas.imass['CO2'] = (solid_cakes.COD * self.carbon_COD_ratio * solid_cakes.F_vol* self.combustion_C_loss/ 1e3 /12 *44 
                            + cmps['WoodPellet']._i_C * wood_pellets.imass['WoodPellet']
                            #assume all carbon in wood pellets convert to CO2
                            )
        air.imass['O2'] = gas.imass['CO2']/44*32
        air.imass['N2'] = air.imass['O2']/32/0.21*0.79*28 #assume air content 21% O2 and 79% N2 
        gas.imass['N2'] = (solid_cakes.imass['N']*self.combustion_N_loss - 
                           N2O.imass['N2O']/44*28 - NO.imass['NO']/30*14) + air.imass['N2']
        
        gas.imass['H2O'] = solid_cakes.imass['H2O'] + wood_pellets.imass['H2O'] - ash.imass['H2O']
        
        
    def _init_lca(self):
        self.construction = [
            Construction('Steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            Construction('Aluminum', linked_unit=self, item='Aluminum', quantity_unit='kg'),
            Construction('WoodPellet', linked_unit=self, item='WoodPellet', quantity_unit='kg'),
            ]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.steel_weight
        design['Aluminum'] = constr[1].quantity = self.aluminum_weight
        design['WoodPellet'] = constr[2].quantity = self.lifetime*365*24* self.wood_pellets_kgphr
        self.add_construction(add_cost = False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Combustor"] = self.combustor_cost * self.price_ratio #USD
        service_factor = self.ppl/120 if self.if_sludge_service else 1
        total_equipment = 0.
        for equipment, cost in C.items():
            C[equipment] = cost * service_factor
            total_equipment += cost
         
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_replacement_cost()+
                         self._calc_maintenance_labor_cost()) #USD/hr
        
    def _calc_replacement_cost(self):
        wood_pellet_cost = self.wood_pellets_kgphr * self.wood_pellets_cost #USD/hr
        service_factor = self.ppl/120 if self.if_sludge_service else 1
        return wood_pellet_cost* self.price_ratio * service_factor # USD/hr
            
    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.combustor_operation_maintenance * self.wages)
        service_factor = self.ppl/120 if self.if_sludge_service else 1
        return maintenance_labor_cost / (365*24) * service_factor
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return 0

#%%

class UFMixer(Mixer):
    '''
    Mixing ultrafiltration reject and liquid from solids separator before delivering to belt separator
    '''
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 rigorous=False, conserve_phases=False):
        Mixer.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        
        
    # def _run(self):
    #     separated_liquid, uf_retentate = self.ins
    #     waste_out, = self.outs
    #     waste_out.mix_from(self.ins)
    #     # print(f"The flushing water flow is {flushing.imass['H2O']} kg/h.")

# %%
toilet_path = ospath.join(g2rt_su_data_path, '_g2rt_toilet.csv')

class FWMixer(Mixer):
    '''
    Mixing tap water and recycled water from RO to meet the flushing water demand
    '''
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 rigorous=False, conserve_phases=False,N_user=None, N_toilet=None,
                 N_tot_user = None, if_flushing=True, flushing_water=None):
        Mixer.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        self.N_user = N_user
        self.N_toilet = N_toilet
        self.N_tot_user = N_tot_user
        self.if_flushing = if_flushing
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])
        
        data = load_data(path=toilet_path)    
        for para in data.index:
            if flushing_water is not None and para == "flushing_water":
                self.flushing_water = flushing_water
                continue  # Skip the material_replacement_cost index if it's provided
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
                               
    def _run(self):
        tap_water, ro_permeate = self.ins
        flushing, = self.outs
        solids = self.solids
        solubles = self.solubles
        N_tot_user = self.N_tot_user or self.N_toilet*self.N_user
        flushing.imass['H2O'] = self.flushing_water*N_tot_user
        tap_water.imass['H2O'] = flushing.imass['H2O'] - ro_permeate.imass['H2O']
        for i in solids+solubles:
            flushing.imass[i] = tap_water.imass[i] + ro_permeate.imass[i]

        # flushing.mix_from((tap_water,ro_permeate))
        # print(f"The flushing water flow is {flushing.imass['H2O']} kg/h.")

class Toilet(SanUnit, Decay, isabstract=True):
    '''
    Abstract class containing common parameters and design algorithms for toilets
    based on `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

    Parameters 
    ----------
    degraded_components : tuple
        IDs of components that will degrade (simulated by first-order decay).
    N_user : int, float
        Number of people per toilet.
        Note that this number can be a float when calculated from `N_tot_user` and `N_toilet`.
    N_toilet : int
        Number of parallel toilets.
        In calculation, `N_toilet` will be calculated as `ceil(N_tot_user/N_user)`.
    N_tot_user : int
        Total number of users.

        .. note::

            If `N_tot_user` is provided (i.e., not "None"),
            then updating `N_user` will recalculate `N_toilet`, and vice versa.

    if_toilet_paper : bool
        If toilet paper is used.
    if_flushing : bool
        If water is used for flushing.
    if_cleansing : bool
        If water is used for cleansing.
    if_desiccant : bool
        If desiccant is used for moisture and odor control.
    if_air_emission : bool
        If emission to air occurs
        (i.e., if the pit is completely sealed off from the atmosphere).
    if_ideal_emptying : bool
        If the toilet appropriately emptied to avoid contamination to the
        environmental.
    CAPEX : float
        Capital cost of a single toilet.
    OPEX_over_CAPEX : float
        Fraction of annual operating cost over total capital cost.
    price_ratio : float
        Calculated capital cost will be multiplied by this number
        to consider the effect in cost difference from different locations.

    References
    ----------
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
    Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
    Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
    https://doi.org/10.1021/acs.est.0c03296.

    See Also
    --------
    :ref:`qsdsan.processes.Decay <processes_Decay>`

    '''
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    density_dct = {
        'Sand': 1442,
        'Gravel': 1600,
        'Brick': 1750,
        'Plastic': 0.63,
        'Steel': 7900,
        'StainlessSteelSheet': 2.64
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), N_user=1, N_toilet=1, N_tot_user=None,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=None, OPEX_over_CAPEX=None, price_ratio=1.):

        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.degraded_components = tuple(degraded_components)
        self._N_user = self._N_toilet = self._N_tot_user = None
        self.N_user = N_user
        self.N_toilet = N_toilet
        self.N_tot_user = N_tot_user
        self.if_toilet_paper = if_toilet_paper
        self.if_flushing = if_flushing
        self.if_cleansing = if_cleansing
        self.if_desiccant = if_desiccant
        self.if_air_emission = if_air_emission
        self.if_ideal_emptying = if_ideal_emptying
        self.CAPEX = CAPEX
        self.OPEX_over_CAPEX = OPEX_over_CAPEX
        self.price_ratio = price_ratio

        data = load_data(path=toilet_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            if para in ('desiccant_V', 'desiccant_rho'):
                setattr(self, para, value)
            else:
                setattr(self, '_'+para, value)
        del data

        self._empty_ratio = 0.59


    def _run(self):
        ur, fec, tp = self.ins
        # if ro_permeate.imass['H2O'] <0.:
        #     print("RO permeate recycled water flow is less than zero!!")
        # elif ro_permeate.imass['H2O'] >1e10:
        #     print(f"RO permeate recycled water flow is too Large: {ro_permeate.imass['H2O']} kg/h!!")

        tp.imass['Tissue'] = int(self.if_toilet_paper)*self.toilet_paper
        # if fw.imass['H2O'] <0.:
        #     print(f"flush water flow is less than zero: {fw.imass['H2O']} kg/h")
        # print(ro_permeate.imass['H2O'])
        # ro_permeate.show()
        # print(f"treated RO influent to reuse as flushing water is {ro_permeate.imass['H2O']}") #TODO:debug
        # fw.show() #TODO:debug
        # - ro_permeate.imass['H2O']/N_tot_user
        # if fw.imass['H2O'] <0.:
        #     print("flushing water flow is less than zero!!")

    def _scale_up_outs(self):
        '''
        Scale up the effluent based on the number of user per toilet and
        toilet number.
        '''
        N_tot_user = self.N_tot_user or self.N_toilet*self.N_user
        for i in self.outs:
            if not i.F_mass == 0:
                i.F_mass *= N_tot_user
                
                # try:
                #     with np.errstate(over='raise'):
                #         i.F_mass *= N_tot_user  
                # except FloatingPointError:
                #     print("Overflow detected. Force RO permeate to be constant.")
                #     self.ins[5].imass['H2O'] = 0.302
                #     i.F_mass *= N_tot_user  
                #     breakpoint()

    def _cost(self):
        self.baseline_purchase_costs['Total toilets'] = self.CAPEX * self.N_toilet * self.price_ratio
        add_OPEX = self.baseline_purchase_costs['Total toilets']*self.OPEX_over_CAPEX/365/24
        self._add_OPEX = {'Additional OPEX': add_OPEX}


    @staticmethod
    def get_emptying_emission(waste, CH4, N2O, empty_ratio, CH4_factor, N2O_factor):
        '''
        Calculate emissions due to non-ideal emptying based on
        `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_,

        Parameters
        ----------
        stream : WasteStream
            Excreta stream that is not appropriately emptied (before emptying).
        CH4 : WasteStream
            Fugitive CH4 gas (before emptying).
        N2O : WasteStream
            Fugitive N2O gas (before emptying).
        empty_ratio : float
            Fraction of excreta that is appropriately emptied..
        CH4_factor : float
            Factor to convert COD removal to CH4 emission.
        N2O_factor : float
            Factor to convert COD removal to N2O emission.

        Returns
        -------
        stream : WasteStream
            Excreta stream that is not appropriately emptied (after emptying).
        CH4 : WasteStream
            Fugitive CH4 gas (after emptying).
        N2O : WasteStream
            Fugitive N2O gas (after emptying).
        '''
        COD_rmvd = waste.COD*(1-empty_ratio)/1e3*waste.F_vol
        CH4.imass['CH4'] += COD_rmvd * CH4_factor
        N2O.imass['N2O'] += COD_rmvd * N2O_factor
        waste.mass *= empty_ratio

    @property
    def N_user(self):
        '''[int, float] Number of people per toilet.'''
        return self._N_user or self.N_tot_user/self.N_toilet
    @N_user.setter
    def N_user(self, i):
        if i is not None:
            N_user = self._N_user = int(i)
            old_toilet = self._N_toilet
            if old_toilet and self.N_tot_user:
                new_toilet = ceil(self.N_tot_user/N_user)
                warn(f'With the provided `N_user`, the previous `N_toilet` of {old_toilet} '
                     f'is recalculated from `N_tot_user` and `N_user` as {new_toilet}.')
                self._N_toilet = None
        else:
            self._N_user = i

    @property
    def N_toilet(self):
        '''[int] Number of parallel toilets.'''
        return self._N_toilet or ceil(self.N_tot_user/self.N_user)
    @N_toilet.setter
    def N_toilet(self, i):
        if i is not None:
            N_toilet = self._N_toilet = ceil(i)
            old_user = self._N_user
            if old_user and self.N_tot_user:
                new_user = self.N_tot_user/N_toilet
                warn(f'With the provided `N_toilet`, the previous `N_user` of {old_user} '
                     f'is recalculated from `N_tot_user` and `N_toilet` as {new_user}.')
                self._N_user = None
        else:
            self._N_toilet = i

    @property
    def N_tot_user(self):
        '''[int] Number of total users.'''
        return self._N_tot_user
    @N_tot_user.setter
    def N_tot_user(self, i):
        if i is not None:
            self._N_tot_user = int(i)
        else:
            self._N_tot_user = None

    @property
    def toilet_paper(self):
        '''
        [float] Amount of toilet paper used
        (if ``if_toilet_paper`` is True), [kg/cap/hr].
        '''
        return self._toilet_paper
    @toilet_paper.setter
    def toilet_paper(self, i):
        self._toilet_paper = i

    @property
    def flushing_water(self):
        '''
        [float] Amount of water used for flushing
        (if ``if_flushing_water`` is True), [kg/cap/hr].
        '''
        return self._flushing_water
    @flushing_water.setter
    def flushing_water(self, i):
        self._flushing_water = i

    @property
    def cleansing_water(self):
        '''
        [float] Amount of water used for cleansing
        (if ``if_cleansing_water`` is True), [kg/cap/hr].
        '''
        return self._cleansing_water
    @cleansing_water.setter
    def cleansing_water(self, i):
        self._cleansing_water = i

    @property
    def desiccant(self):
        '''
        [float] Amount of desiccant used (if ``if_desiccant`` is True), [kg/cap/hr].

        .. note::

            Value set by ``desiccant_V`` and ``desiccant_rho``.

        '''
        return self.desiccant_V*self.desiccant_rho

    @property
    def N_volatilization(self):
        '''
        [float] Fraction of input N that volatilizes to the air
        (if ``if_air_emission`` is True).
        '''
        return self._N_volatilization
    @N_volatilization.setter
    def N_volatilization(self, i):
        self._N_volatilization = i

    @property
    def empty_ratio(self):
        '''
        [float] Fraction of excreta that is appropriately emptied.

        .. note::

            Will be 1 (i.e., 100%) if ``if_ideal_emptying`` is True.

        '''
        if self.if_ideal_emptying:
            return 1.
        return self._empty_ratio
    @empty_ratio.setter
    def empty_ratio(self, i):
        if self.if_ideal_emptying:
            warn(f'`if_ideal_emptying` is True, the set value {i} is ignored.')
        self._empty_ratio = i

    @property
    def MCF_aq(self):
        '''[float] Methane correction factor for COD lost due to inappropriate emptying.'''
        return self._MCF_aq
    @MCF_aq.setter
    def MCF_aq(self, i):
        self._MCF_aq = i

    @property
    def N2O_EF_aq(self):
        '''[float] Fraction of N emitted as N2O due to inappropriate emptying.'''
        return self._N2O_EF_aq
    @N2O_EF_aq.setter
    def N2O_EF_aq(self, i):
        self._N2O_EF_aq = i

    @property
    def if_N2O_emission(self):
        '''[bool] Whether to consider N degradation and fugitive N2O emission.'''
        return self.if_air_emission
    @if_N2O_emission.setter
    def if_N2O_emission(self, i):
        raise ValueError('Setting `if_N2O_emission` for `PitLatrine` is not supported, '
                         'please set `if_air_emission` instead.')


# %%
surt_path = ospath.join(surt_su_data_path, '_surt.csv')

@price_ratio()
class SURT(Toilet):
    '''
    Single-unit reinvented toilet.

    The following components should be included in system thermo object for simulation:
    Tissue, WoodAsh, H2O, NH3, NonNH3, P, K, Mg, CH4, N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    Ceramic, Fan.

    Parameters
    ----------
    ins : Iterable(stream)
        waste_in: mixed excreta.
    Outs : Iterable(stream)
        waste_out: degraded mixed excreta.
        CH4: fugitive CH4.
        N2O: fugitive N2O.
    N_squatting_pan_per_toilet : int
        The number of squatting pan per toilet.
    N_urinal_per_toilet : int
        The number of urinals per toilet.
    if_include_front_end : bool
        If False, will not consider the capital and operating costs of this unit.

    See Also
    --------
    :ref:`qsdsan.sanunits.Toilet <sanunits_toilet>`
    '''
    _N_outs = 3
    _units = {
        'Collection period': 'd',
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), N_user=1, N_tot_user=1,
                 N_toilet=None, if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=True, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=0, OPEX_over_CAPEX=0, lifetime=10,
                 N_squatting_pan_per_toilet=1, N_urinal_per_toilet=1,
                 if_include_front_end=True, **kwargs):

        Toilet.__init__(
            self, ID, ins, outs, thermo=thermo, init_with=init_with,
            degraded_components=degraded_components,
            N_user=N_user, N_tot_user=N_tot_user, N_toilet=N_toilet,
            if_toilet_paper=if_toilet_paper, if_flushing=if_flushing,
            if_cleansing=if_cleansing, if_desiccant=if_desiccant,
            if_air_emission=if_air_emission, if_ideal_emptying=if_ideal_emptying,
            CAPEX=CAPEX, OPEX_over_CAPEX=OPEX_over_CAPEX
            )
        self.N_squatting_pan_per_toilet = N_squatting_pan_per_toilet
        self.N_urinal_per_toilet = N_urinal_per_toilet
        self.if_include_front_end = if_include_front_end
        self._mixed_in = WasteStream(f'{self.ID}_mixed_in')

        data = load_data(path=surt_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [
            Construction(item='Ceramic', linked_unit=self, quantity_unit='kg'),
            Construction(item='Fan', linked_unit=self, quantity_unit='kg'),
        ]


    def _run(self):
        Toilet._run(self)
        mixed_out, CH4, N2O = self.outs
        CH4.phase = N2O.phase = 'g'

        mixed_in = self._mixed_in
        mixed_in.mix_from(self.ins)
        tot_COD_kg = sum(float(getattr(i, 'COD')) * i.F_vol for i in self.ins) / 1e3

        # Air emission
        if self.if_air_emission:
            # N loss due to ammonia volatilization
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(mixed_in.TN/1e3*mixed_in.F_vol*self.N_volatilization,
                                        mixed_in.imass['NH3'])
            mixed_in.imass ['NH3'] -= NH3_rmd
            mixed_in.imass['NonNH3'] -= NonNH3_rmd
            
            # Energy/N loss due to degradation
            mixed_in._COD = tot_COD_kg * 1e3 / mixed_in.F_vol # accounting for COD loss in leachate
            Decay._first_order_run(self, waste=mixed_in, treated=mixed_out, CH4=CH4, N2O=N2O)
        else:
            mixed_out.copy_like(mixed_in)
            CH4.empty()
            N2O.empty()
            
        # Aquatic emission when not ideally emptied
        if not self.if_ideal_emptying:
           self.get_emptying_emission(
                waste=mixed_out, CH4=CH4, N2O=N2O,
                empty_ratio=self.empty_ratio,
                CH4_factor=self.COD_max_decay*self.MCF_aq*self.max_CH4_emission,
                N2O_factor=self.N2O_EF_decay*44/28)
        self._scale_up_outs()


    def _design(self):
        design = self.design_results
        constr = self.construction
        if self.if_include_front_end:
            design['Number of users per toilet'] = self.N_user
            design['Parallel toilets'] = N = self.N_toilet
            design['Collection period'] = self.collection_period
            design['Ceramic'] = Ceramic_quant = (
                self.squatting_pan_weight * self.N_squatting_pan_per_toilet+
                self.urinal_weight * self.N_urinal_per_toilet
                )
            design['Fan'] = Fan_quant = 1  # assume fan quantity is 1
            constr[0].quantity = Ceramic_quant * N
            constr[1].quantity = Fan_quant * N
            self.add_construction(add_cost=False)
        else:
            design.clear()
            for i in constr: i.quantity = 0
            

    def _cost(self):
        C = self.baseline_purchase_costs
        if self.if_include_front_end:
            N_toilet = self.N_toilet
            C['Ceramic Toilets'] = (
                self.squatting_pan_cost * self.N_squatting_pan_per_toilet +
                self.urinal_cost * self.N_urinal_per_toilet
                ) * N_toilet
            C['Fan'] = self.fan_cost * N_toilet
            C['Misc. parts'] = (
                self.led_cost +
                self.anticor_floor_cost +
                self.circuit_change_cost +
                self.pipe_cost
                ) * N_toilet

            ratio = self.price_ratio
            for equipment, cost in C.items():
                C[equipment] = cost * ratio
        else:
            self.baseline_purchase_costs.clear()

        sum_purchase_costs = sum(v for v in C.values())
        self.add_OPEX = (
            self._calc_replacement_cost() +
            self._calc_maintenance_labor_cost() +
            sum_purchase_costs * self.OPEX_over_CAPEX / (365 * 24)
            )

    def _calc_replacement_cost(self):
        return 0

    def _calc_maintenance_labor_cost(self):
        return 0
    
    @property
    def power_kW(self):
        return 0.0 #kW
    
        
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day

    @property
    def collection_period(self):
        '''[float] Time interval between storage tank collection, [d].'''
        return self._collection_period
    @collection_period.setter
    def collection_period(self, i):
        self._collection_period = float(i)
        
    @property
    def tau(self):
        '''[float] Retention time of the unit, same as `collection_period`.'''
        return self.collection_period
    @tau.setter
    def tau(self, i):
        self.collection_period = i



#%%
excretion_path = ospath.join(data_path, 'sanunit_data/_excretion.tsv')
class Excretion(SanUnit):
    '''
    Estimation of N, P, K, and COD in urine and feces based on dietary intake
    for one person based on `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

    Parameters
    ----------
    waste_ratio : float
        A ratio in [0, 1] to indicate the amount of intake calories and nutrients
        (N, P, K) that is wasted.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_

    References
    ----------
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
    Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
    Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
    https://doi.org/10.1021/acs.est.0c03296.
    '''

    _N_ins = 0
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 waste_ratio=0, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.waste_ratio = waste_ratio

        data = load_data(path=excretion_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _run(self):
        ur, fec = self.outs
        ur.empty()
        fec.empty()

        not_wasted = 1 - self.waste_ratio
        factor = 24 * 1e3 # from g per person per day to kg per hour

        ur_N = (self.p_veg+self.p_anim)/factor*self.N_prot \
           * self.N_exc*self.N_ur*not_wasted
        ur.imass['NH3'] = ur_N * self.N_ur_NH3
        ur.imass['NonNH3'] = ur_N - ur.imass['NH3']

        ur.imass['P'] = (self.p_veg*self.P_prot_v+self.p_anim*self.P_prot_a)/factor \
            * self.P_exc*self.P_ur*not_wasted

        e_cal = self.e_cal / 24 * not_wasted
        ur.imass['K'] = e_cal/1e3 * self.K_cal/1e3 * self.K_exc*self.K_ur
        ur.imass['Mg'] = self.Mg_ur / factor
        ur.imass['Ca'] = self.Ca_ur / factor

        ur_exc = self.ur_exc / factor
        ur.imass['H2O'] = self.ur_moi * ur_exc
        ur.imass['OtherSS'] = ur_exc - ur.F_mass

        fec_exc = self.fec_exc / factor
        fec_N = (1-self.N_ur)/self.N_ur * ur_N
        fec.imass['NH3'] = fec_N * self.N_fec_NH3
        fec.imass['NonNH3'] = fec_N - fec.imass['NH3']
        fec.imass['P'] = (1-self.P_ur)/self.P_ur * ur.imass['P']
        fec.imass['K'] = (1-self.K_ur)/self.K_ur * ur.imass['K']
        fec.imass['Mg'] = self.Mg_fec / factor
        fec.imass['Ca'] = self.Ca_fec / factor
        fec.imass['H2O'] = self.fec_moi * fec_exc
        fec.imass['OtherSS'] = fec_exc - fec.F_mass

        # 14 kJ/g COD, the average lower heating value of excreta
        tot_COD = e_cal*self.e_exc*4.184/14/1e3 # in kg COD/hr
        ur.imass['sCOD'] = tot_COD*(1-self.e_fec) # in kg/hr
        fec.imass['xCOD'] = tot_COD*self.e_fec # in kg/hr
        
    @property
    def power_kW(self):
        return 0.0 #kW

    @property
    def e_cal(self):
        '''[float] Caloric intake, [kcal/cap/d].'''
        return self._e_cal
    @e_cal.setter
    def e_cal(self, i):
        self._e_cal = i

    @property
    def p_veg(self):
        '''[float] Vegetal protein intake, [g/cap/d].'''
        return self._p_veg
    @p_veg.setter
    def p_veg(self, i):
        self._p_veg = i

    @property
    def p_anim(self):
        '''[float] Animal protein intake, [g/cap/d].'''
        return self._p_anim
    @p_anim.setter
    def p_anim(self, i):
        self._p_anim = i

    @property
    def N_prot(self):
        '''[float] Nitrogen content in protein, [wt%].'''
        return self._N_prot
    @N_prot.setter
    def N_prot(self, i):
        self._N_prot = i

    @property
    def P_prot_v(self):
        '''[float] Phosphorus content in vegetal protein, [wt%].'''
        return self._P_prot_v
    @P_prot_v.setter
    def P_prot_v(self, i):
        self._P_prot_v = i

    @property
    def P_prot_a(self):
        '''[float] Phosphorus content in animal protein, [wt%].'''
        return self._P_prot_a
    @P_prot_a.setter
    def P_prot_a(self, i):
        self._P_prot_a = i

    @property
    def K_cal(self):
        '''[float] Potassium intake relative to caloric intake, [g K/1000 kcal].'''
        return self._K_cal
    @K_cal.setter
    def K_cal(self, i):
        self._K_cal = i

    @property
    def N_exc(self):
        '''[float] Nitrogen excretion factor, [% of intake].'''
        return self._N_exc
    @N_exc.setter
    def N_exc(self, i):
        self._N_exc = i

    @property
    def P_exc(self):
        '''[float] Phosphorus excretion factor, [% of intake].'''
        return self._P_exc
    @P_exc.setter
    def P_exc(self, i):
        self._P_exc = i

    @property
    def K_exc(self):
        '''[float] Potassium excretion factor, [% of intake].'''
        return self._K_exc
    @K_exc.setter
    def K_exc(self, i):
        self._K_exc = i

    @property
    def e_exc(self):
        '''[float] Energy excretion factor, [% of intake].'''
        return self._e_exc
    @e_exc.setter
    def e_exc(self, i):
        self._e_exc = i

    @property
    def N_ur(self):
        '''[float] Nitrogen recovered in urine, [wt%].'''
        return self._N_ur
    @N_ur.setter
    def N_ur(self, i):
        self._N_ur = i

    @property
    def P_ur(self):
        '''[float] Phosphorus recovered in urine, [wt%].'''
        return self._P_ur
    @P_ur.setter
    def P_ur(self, i):
        self._P_ur = i

    @property
    def K_ur(self):
        '''[float] Potassium recovered in urine, [wt%].'''
        return self._K_ur
    @K_ur.setter
    def K_ur(self, i):
        self._K_ur = i

    @property
    def e_fec(self):
        '''[float] Percent of excreted energy in feces, [%].'''
        return self._e_fec
    @e_fec.setter
    def e_fec(self, i):
        self._e_fec = i

    @property
    def N_ur_NH3(self):
        '''[float] Reduced inorganic nitrogen in urine, modeled as NH3, [% of total urine N].'''
        return self._N_ur_NH3
    @N_ur_NH3.setter
    def N_ur_NH3(self, i):
        self._N_ur_NH3 = i

    @property
    def N_fec_NH3(self):
        '''[float] Reduced inorganic nitrogen in feces, modeled as NH3, [% of total feces N].'''
        return self._N_fec_NH3
    @N_fec_NH3.setter
    def N_fec_NH3(self, i):
        self._N_fec_NH3 = i

    @property
    def ur_exc(self):
        '''[float] Urine generated per day, [g/cap/d].'''
        return self._ur_exc
    @ur_exc.setter
    def ur_exc(self, i):
        self._ur_exc = i

    @property
    def fec_exc(self):
        '''[float] Feces generated per day, [g/cap/d].'''
        return self._fec_exc
    @fec_exc.setter
    def fec_exc(self, i):
        self._fec_exc = i

    @property
    def ur_moi(self):
        '''[float] Moisture (water) content of urine, [wt%].'''
        return self._ur_moi
    @ur_moi.setter
    def ur_moi(self, i):
        self._ur_moi = i

    @property
    def fec_moi(self):
        '''[float] Moisture (water) content of feces, [wt%].'''
        return self._fec_moi
    @fec_moi.setter
    def fec_moi(self, i):
        self._fec_moi = i

    @property
    def Mg_ur(self):
        '''[float] Magnesium excreted in urine, [g Mg/cap/d].'''
        return self._Mg_ur
    @Mg_ur.setter
    def Mg_ur(self, i):
        self._Mg_ur = i

    @property
    def Mg_fec(self):
        '''[float] Magnesium excreted in feces, [g Mg/cap/d].'''
        return self._Mg_fec
    @Mg_fec.setter
    def Mg_fec(self, i):
        self._Mg_fec = i

    @property
    def Ca_ur(self):
        '''[float] Calcium excreted in urine, [g Ca/cap/d].'''
        return self._Ca_ur
    @Ca_ur.setter
    def Ca_ur(self, i):
        self._Ca_ur = i

    @property
    def Ca_fec(self):
        '''[float] Calcium excreted in feces, [g Ca/cap/d].'''
        return self._Ca_fec
    @Ca_fec.setter
    def Ca_fec(self, i):
        self._Ca_fec = i

    @property
    def waste_ratio(self):
        '''
        [float] The amount of intake calories and nutrients
        (N, P, K) that is wasted.

        .. note::
            Not considered for Mg and Ca.
        '''
        return self._waste_ratio
    @waste_ratio.setter
    def waste_ratio(self, i):
        self._waste_ratio = i

#%%
vr_pasteurization_path = ospath.join(g2rt_su_data_path, '_vr_pasteurization.csv')

@price_ratio()
class VRpasteurization(SanUnit):
    '''
    Pasteurizer in volume reduction toilet to remove pathogens in the homogenized 
    solids. Heating is from electricity (Joule heater).
    
    Parameters
    -----------
    ins: Iterable(stream)
        solids: solids produced from the homogenizer in the volume reduction toilet.

    outs: Iterable(stream)
        treated solids: solids treated from pasteurization.
    temp_pasteurization : float
        Pasteurization temperature is 70°C or 343.15 K.
    solids_inlet_temp : float
        Temperature of solids from the inlet.
    heat_loss : float
        Heat loss during pasteurization process is assumed to be 30%
    
    See Also
    --------
    :class:`~.sanunits.SludgePasteurization`
    '''
    
    _N_ins = 1
    _N_outs = 1
    
    # Specific Heat capacity of water
    Cp_w = 4.184 # kJ kg^-1 K^-1
    # Specific Heat capacity of dry matter (sludge)
    Cp_dm = 1.231 # kJ kg^-1 K^-1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 heat_loss=0.3, solids_inlet_temp=273.15+25.6,
                 temp_pasteurization= 273.15+90, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                          F_BM_default=1)
        self.heat_loss = heat_loss
        self.solids_inlet_temp = solids_inlet_temp
        self.temp_pasteurization = temp_pasteurization
        
        data = load_data(path=vr_pasteurization_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    def _run(self):
        solids_in, = self.ins
        solids_out, = self.outs
        solids_out.copy_like(solids_in)
        
    def _init_lca(self):
        self.construction = [Construction("stainless_steel", linked_unit=self,
                                          item = "StainlessSteel", 
                                          quantity_unit= "kg"),
                             Construction("polyethylene", linked_unit=self,
                                          item = "Polyethylene",
                                          quantity_unit= "kg")]
        
    def _cost(self):
        C=self.baseline_purchase_costs
        C["Heating"]= (self.joule_heater_cost+
                       self.thermocouple_port_cost* self.thermocouple_port_quantity+
                       self.temperature_control
                       ) * self.price_ratio
        C["Valves"] = self.valve_cost * self.valve_quantity * self.price_ratio
        C["Tubing"] = self.tubing_cost * self.price_ratio
        C["Misc.parts"] = self.miscellaneous_cost_ratio*(C["Heating"]+
                               C["Valves"]+
                               C["Tubing"])
        solids_in, = self.ins
        solids_in.T = self.solids_inlet_temp + 273.15
        solids_out, = self.outs
        solids_out.T = self.temp_pasteurization + 273.15
        # Overall heat required for pasteurization
        temp_diff = self.temp_pasteurization - self.solids_inlet_temp #K
        Q_d = (solids_in.imass['H2O']*self.Cp_w + 
               (solids_in.F_mass-solids_in.imass['H2O'])*self.Cp_dm)*temp_diff #kJ/hr
        Q_tot = Q_d/(1-self.heat_loss/100) # kJ/hr, 30% of the total generated is lost
        heating_electricity = Q_tot/3600 #kWh/hr
        # self.add_heat_utility(unit_duty=self.Hnet, T_in=solids_in.T, 
        #                       T_out = solids_out.T, heat_transfer_efficiency=(1-self.heat_loss/100),
        #                       hxn_ok= True)
        self.power_utility(heating_electricity) #kWh/hr
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + #USD/hr, assume 
                         #replacement cost a fraction of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
        
        # def _calc_replacement_cost(self): #USD/hr, assume 5% of CAPEX per year
        #     replacement_cost = 0.05*(
        #         C["Heating"]+ C["Valves"] + C["Tubing"] + C["Misc.parts"]
        #         ) #USD/yr
        #     return replacement_cost / (365*24)
    
    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.pasteurizer_maintenance * self.wages) #USD/yr
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return self.power_utility.rate #kW


#%%
G2RT_homogenizer_path = ospath.join(g2rt_su_data_path, '_g2rt_homogenizer.csv')

@price_ratio()
class G2RThomogenizer(Copier):
    '''
    Homogenizer and buffer tanks in generation II reinveted toilets to  break 
    up solids [1].
    
    .. note:

        This is a non-reactive unit (i.e., the effluent is copied from the influent)

    The following components should be included in system thermo object for simulation:
    H2O, OtherSS.

    The following impact items should be pre-constructed for life cycle assessment:
    Steel.

    Parameters
    ----------
    ins : Iterable(stream)
        Influent stream.
    outs : Iterable(stream)
        Effluent stream, is copied from the influent.
    moisture_content_out : float
        Moisture content of the effluent solids stream.

    References
    ----------
    [1] YEE et al. Buffer tank separation and homogenization system. 
    https://patents.google.com/patent/WO2023288114A1/en?oq=WO2023288114A1
    [2] YEE et al. Volume reduction non-sewered single unit toilet system.
    https://patents.google.com/patent/WO2023288326A1/en?oq=WO2023288326A1
    
    See Also
    ---------
    :class:`~.sanunits.BiogenicRefineryGrinder`
    '''
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',ppl=6,
                 **kwargs):
        Copier.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1) 
        self.ppl = ppl
        data = load_data(path=G2RT_homogenizer_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _init_lca(self):
        self.construction = [Construction("stainless_steel", linked_unit=self,
                                          item = "StainlessSteel", 
                                          quantity_unit= "kg"),
                             Construction("polyethylene", linked_unit=self,
                                          item = "Polyethylene",
                                          quantity_unit= "kg")]
    
    # def _run(self):
    #     waste_in = self.ins[0]
    #     waste_out = self.outs[0]
    #     waste_out.copy_like(waste_in)

    #     mc_in = waste_in.imass['H2O'] / waste_in.F_mass # fraction
    #     mc_out = self.moisture_content_out/100 #convert to fraction
    #     
    #     if mc_in < mc_out*0.999:
    #         raise RuntimeError(f'Moisture content of the influent stream ({mc_in:.2f}) '
    #                            f'is smaller than the desired moisture content ({mc_out:.2f}).')
    #     TS_in = waste_in.F_mass - waste_in.imass['H2O'] # kg TS dry/hr
    #     waste_out.imass['H2O'] = TS_in/(1-mc_out)*mc_out
    #     waste_out._COD = waste_in.COD * waste_in.F_vol / waste_out.F_vol
        
    def _design(self): #kg
        design = self.design_results
        constr = self.construction    
        self.design_results['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        self.design_results['Polyethylene'] = constr[1].quantity = self.polyethylene_weight
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Macerator"] = self.macerator_cost * self.price_ratio
        C["Misc.parts"] = self.miscellaneous_cost_ratio * C["Macerator"]
        
        self.power_utility(self.macerator_power_demand * self.macerator_daily_operation
                           *self.ppl/24) # kWh/hr
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost a fraction of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
    
    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.homogenizer_maintenance * self.wages) #USD/yr
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.macerator_power_demand * self.macerator_daily_operation
                           *self.ppl/24) #kW
    
#%%
vr_filter_press_path = ospath.join(g2rt_su_data_path, '_vr_filter_press.csv')

@price_ratio()
class VolumeReductionFilterPress(SanUnit): 
    
    '''
    A filter press unit for the dewatering of mixed excreta in volume reduction
    generation II reinveted toilet [1]
    
    The following componenet should be included in system thermo object for simulation:
    Water.
    
    The following impact items should be pre-constructed for life cycle assessment:
    Stainless steel, Pump.
    
    Parameters
    ----------
    ins: Iterable(stream)
      Pasteurized solids waste for dewatering treatment
    outs: Iterable(stream)
      Liquids and solids produced from filter press.
    sludge_moisture: float
      Moisture content of the solids cake after filter press [wt% water].
    
    References
    -----------
    [1] YEE et al. VOLUME REDUCTION SOLIDS TREATMENT SYSTEM. 
    https://patents.google.com/patent/WO2023288327A1/en?oq=WO2023288327A1
    
    [2] Bev Express multi-plate sheet filter (ME-10 model)
    https://417cb0.p3cdn1.secureserver.net/wp-content/Documents/Tech%20Sheets/Filtration%20Equipment/Plate%20&%20Frame/Beverage%20Express/ErtelAlsop_Bev_ExPRESS_Tech_Sheet.pdf?dl=1
    
    [3] https://www.earthshields.com/how-much-does-geotextile-fabric-cost/# accessed on yyyy-mm-dd
    
    [4] https://www.suezwaterhandbook.com/processes-and-technologies/liquid-sludge-treatment/filter-press/conventional-recessed-plate-filter-press#:~:text=The%20filter%20press%20consumes%20relatively,solids%20depending%20on%20sludge%20type.
    
    [5] https://multimedia.3m.com/mws/media/2113780O/3m-zeta-plus-vs-filter-press-for-beer-clarification-application-brief.pdf
        
    '''
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', 
                 ins=None, outs=(), thermo=None, 
                 init_with='WasteStream',
                 solids = (),
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                                solids=solids,
                                **kwargs)
                        
        data = load_data(path=vr_filter_press_path)
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
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            Construction('pump', linked_unit=self, 
                         item='Pump', 
                         quantity_unit='ea'),
            ]
    
    def _run(self):
        pasteurized_solids, = self.ins
        supernatant, solid_cakes = self.outs
        solid_cakes.phase = 's'
        solubles, solids = self.solubles, self.solids
        solid_cakes.copy_flow(pasteurized_solids,solids) #all solids go to sludge
        solid_cakes.imass[solids] = pasteurized_solids.imass[solids]*self.TSS_removal/100
        
        mc_in = pasteurized_solids.imass['H2O'] / pasteurized_solids.F_mass # fraction
        mc_out = self.moisture_content_out/100 #convert to fraction
        
        if mc_in < mc_out*0.999:
            mc_out = mc_in
        TS_in = pasteurized_solids.imass[solids].sum() # kg TS dry/hr
        TS_out = solid_cakes.imass[solids].sum()
        #calculate water and solid COD in the solid cakes
        solid_cakes.imass['H2O'] = TS_in/(1-mc_out)*mc_out
        solid_cakes.imass[solubles] = pasteurized_solids.imass[solubles]*\
            (TS_out/(1-mc_out)-TS_out)/(pasteurized_solids.F_mass-TS_in)
        supernatant.mass = pasteurized_solids.mass-solid_cakes.mass
        
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.filterpress_ss_weight
        design['Pump'] = constr[1].quantity = 1
        self.add_construction(add_cost = False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Filter Press"] = self.filterpress_purchase_cost * self.price_ratio #USD
        
        self.power_utility(self.filterpress_energy_persolids*self.outs[1].F_mass
                           /self.filterpress_energy_efficiency
                           ) # kW
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
            
    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.filterpress_maintenance * self.wages)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.filterpress_energy_persolids*self.outs[1].F_mass) #kW

#%%
vr_concentrator_path = ospath.join(g2rt_su_data_path, '_vr_concentrator.csv')

@price_ratio()
class VRConcentrator(SanUnit): 
    '''
    This concentrator unit is used in solids treatmnet in volume reduction 
    generation II reinveted toilet [1]. The heat source is from a heating coil.
#TODO: consider installing a heat exchange to receive heat from pasteurization.
    The following components should be included in system thermo object for simulation:
    H2O, N, CH4, N2O.
    
    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, HeatingUnit, ElectricMotor, Pump, Polyethylene.

    Parameters
    ----------
    ins : Iterable(stream)
        RO reject liquid.
    outs : Iterable(stream)
        Condensed effluent, recirculation water, fugitive N2O, fugitive CH4, fugitive NH3, water vapor
    
    Warnings
    --------
    Energy balance is not performed for this unit.

    References
    ----------
    [1] YEE et al. VOLUME REDUCTION SOLIDS TREATMENT SYSTEM. 
       https://patents.google.com/patent/WO2023288327A1/en?oq=WO2023288327A1

    See Also
    --------
    :class:`~.sanunits.BiogenicRefineryHHXdryer`
    :class:`~.sanunits.BiogenicRefineryHHX`
    '''
    _N_ins = 1
    _N_outs = 6
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
                
        data = load_data(path=vr_concentrator_path)
        
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids)) + ('OtherSS',) #dissolved salts dry out
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])
    
    def _run(self):
        waste_in, = self.ins
        waste_out, excess_water, N2O, CH4, NH3_gas, water_vapor  = self.outs
        N2O.phase = CH4.phase = NH3_gas.phase = water_vapor.phase = 'g'
        waste_out.copy_like(self.ins[0])
        solubles, solids = self.solubles, self.solids
        #Calculate water and solids in the condensed effluent
        mc_in = waste_in.imass['H2O'] / waste_in.F_mass # fraction
        mc_out = self.moisture_content_out/100 #convert to fraction
        if mc_in < mc_out*0.999:
            mc_out = mc_in
        TS_in = waste_in.imass[solids].sum() # kg TS dry/hr
        waste_out.imass['H2O'] = TS_in/(1-mc_out)*mc_out #kg water/hr
        #Calculate N2O and CH4 emissions
        CH4.imass['CH4'] = drying_CH4_to_air = \
            self.drying_CH4_emissions * self.carbon_COD_ratio * \
            waste_in.COD * waste_in.F_vol / 1000 # kg CH4 /hr
        drying_NH3_to_air = self.drying_NH3_emissions * waste_in.imass['NH3'] # kg NH3 /hr
        N2O.imass['N2O'] = drying_NH3_to_air * self.NH3_to_N2O # kg N2O /hr
        waste_out.imass['NH3'] = waste_in.imass['NH3'] - drying_NH3_to_air

        NH3_gas.imass['NH3'] = drying_NH3_to_air # kg NH3 /hr
        excess_water.imass['H2O'] = self.excess_water_recirculation * (waste_in.imass['H2O'] - waste_out.imass['H2O'])
        excess_water.imass[solubles] = self.excess_water_recirculation * waste_in.imass[solubles]
        water_vapor.imass['H2O'] = waste_in.imass['H2O'] - waste_out.imass['H2O'] -excess_water.imass['H2O']  #kg H2O/hr
        # Store the calculated value of water_vapor.imass['H2O'] for use in the _cost function
        self.water_vapor_H2O = water_vapor.imass['H2O'] 
        #Calculate COD
        drying_CO2_to_air = (self.drying_CO2_emissions * self.carbon_COD_ratio
                             * (waste_in.COD * waste_in.F_vol-
                                excess_water.COD * excess_water.F_vol
                                ) / 1000) # kg CO2 /hr
        # 44/12/16 are the molecular weights of CO2, C, and CH4, respectively
        waste_out.imass['sCOD'] -=  ((drying_CO2_to_air/44*12+drying_CH4_to_air/16*12)/ self.carbon_COD_ratio + 
                                     excess_water.imass['sCOD'])
    
    def _init_lca(self):
        self.construction = [
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            Construction('heating_unit',linked_unit=self, 
                         item='HeatingUnit', 
                         quantity_unit='ea'),
            Construction('electric_motor',linked_unit=self, 
                         item='ElectricMotor', 
                         quantity_unit='kg'),
            Construction('pump',linked_unit=self, 
                         item='Pump',
                         quantity_unit='ea'),
            Construction('fan',linked_unit=self, 
                         item='Fan',
                         quantity_unit='kg'),
            Construction('polyethylene',linked_unit=self,
                         item='Polyethylene', 
                         quantity_unit='kg'),
            ]
        
    def _design(self):
        design=self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        design['HeatingUnit'] = constr[1].quantity = 1
        design['ElectricMotor'] = constr[2].quantity = self.motor_weight
        design['Pump'] = constr[3].quantity = 1
        design['Fan'] = constr[4].quantity = self.fan_quantity*0.2 #one fan is 0.2 kg
        design['Polyethylene'] = constr[5].quantity = self.polyethylene_weight
        
    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        #C['Stainless steel'] = self.stainless_steel_cost * D['StainlessSteel']
        C['HeatingUnit'] = self.heating_unit_purchase_cost * D['HeatingUnit']
        C['ElectricMotor'] = self.motor_purchase_cost * D['ElectricMotor']
        C['Pump'] = self.pump_purchase_cost * D['Pump']
        #C['Polyethylene'] = self.polyethylene_cost * D['Polyethylene']
        C['Fan'] = self.fan_purchase_cost * self.fan_quantity
        C['Disc'] = self.disc_purchase_cost
        C['Others'] = (self.enclosure_purchase_cost + 
                       self.open_concentrator_vessel_purchase_cost+
                       self.axle_purchase_cost+
                       self.thermistor_cost                       
                       )
        C["Misc.parts"] = self.miscellaneous_cost_ratio*(C['HeatingUnit'] +
                               C['ElectricMotor']+
                               C['Pump'] +
                               C['Fan'] +
                               C['Disc'] + 
                               C['Others']
                               ) 
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio
        
        self.power_utility(self.pump_power_demand * self.pump_daily_operation/24+
                            self.water_vapor_H2O * self.energy_required_to_evaporize_water
                            ) # kW
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
   
    def _calc_maintenance_labor_cost(self): #USD/hr
       maintenance_labor_cost= (self.concentrator_maintenance * self.wages) #USD/yr
       return maintenance_labor_cost / (365*24)
   
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.pump_power_demand * self.pump_daily_operation/24+
                self.water_vapor_H2O * self.energy_required_to_evaporize_water) #kW
   
    
#%%
g2rt_liquids_tank_path = ospath.join(g2rt_su_data_path, '_g2rt_liquids_tank.csv')
@price_ratio()

class G2RTLiquidsTank(Mixer):
    '''
    Liquids storage unit for generation II reinveted toilets to accumulate enough
    liquid waste before ultrafiltration
    
    This is a non-reactive unit, (i.e., the effluent is copied from the mix of influent).
    
    The following impact items should be pre-constructed for life cycle assessment:
    Pump, Polyethylene, Polycarbonate
    
    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.

    See Also
    --------
    :class:`~.sanunits.Copier`
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 rigorous=False, conserve_phases=False):
        Mixer.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)

        data = load_data(path=g2rt_liquids_tank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

    def _init_lca(self):
        self.construction = [
            Construction(item='Polyethylene', linked_unit=self, quantity_unit='kg'),
            Construction(item='Polycarbonate', linked_unit=self, quantity_unit='kg'),
            Construction(item='Pump', linked_unit=self, quantity=1., quantity_unit='ea'),
            ]
    
    def _run(self):
        s_out, = self.outs
        s_out.mix_from(self.ins, vle=self.rigorous,
                       conserve_phases=getattr(self, 'conserve_phases', None))
        V = s_out.vapor_fraction
        if V == 0:
            self._B = 0
        elif V == 1:
            self._B = np.inf
        else:
            self._B = V / (1 - V)
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Polyethylene'] = constr[0].quantity = self.liquids_tank_polyethylene_weight
        design['Polycarbonate'] = constr[1].quantity = self.liquids_tank_polycarbonate_weight
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.liquids_tank_cost
        C['Tubing'] = self.liquids_tank_tubing
        C['Pump'] = self.liquids_tank_feed_pump_cost
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Tank']+
                                                         C['Tubing']+
                                                         C['Pump'])
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

        power_demand = self.pump_energy_per_cycle*self.pump_batch_cycle_per_day
        power_demand = power_demand / 24  # convert from kWh/d to kW
        self.power_utility(power_demand)  # kW
        
    def _calc_replacement_cost(self):
        pump_replacement_cost = self.liquids_tank_feed_pump_cost/ self.pump_lifetime
        return pump_replacement_cost/(365 * 24)* self.price_ratio # USD/hr

    def _calc_maintenance_labor_cost(self):
        liquids_tank_maintenance_labor_cost = (
            self.labor_pump_replacement +
            self.labor_tank_cleaning
            ) * self.wages
        return liquids_tank_maintenance_labor_cost/(365 * 24) # USD/hr
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return self.pump_energy_per_cycle*self.pump_batch_cycle_per_day/24 #kW

#%%
g2rt_solids_tank_path = ospath.join(g2rt_su_data_path, '_g2rt_solids_tank.csv')
@price_ratio()

class G2RTSolidsTank(Copier):
    '''
    Solids storage unit for generation II reinveted toilets to accumulate enough
    solid waste before homogenizer.
    
    This is a non-reactive unit (i.e., the effluent is copied from the influent).
    
    The following impact items should be pre-constructed for life cycle assessment:
    Polyethylene, Polycarbonate
    
    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.

    See Also
    --------
    :class:`~.sanunits.Copier`
    '''
    def __init__(self, ID='', ins=None, outs=(),  thermo=None, init_with='WasteStream',
                 **kwargs):
        Copier.__init__(self, ID, ins, outs, thermo, init_with)

        data = load_data(path=g2rt_solids_tank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [
            Construction(item='Polyethylene', linked_unit=self, quantity_unit='kg'),
            Construction(item='Polycarbonate', linked_unit=self, quantity_unit='kg'),
            ]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Polyethylene'] = constr[0].quantity = self.solids_tank_polyethylene_weight
        design['Polycarbonate'] = constr[1].quantity = self.solids_tank_polycarbonate_weight
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.solids_tank_cost
        C['Piping'] = self.solids_tank_piping
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Tank']+
                                                         C['Piping'])
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        self.add_OPEX =  self._calc_maintenance_labor_cost()

    def _calc_maintenance_labor_cost(self):
        solids_tank_maintenance_labor_cost = (
            self.labor_tank_cleaning
            ) * self.wages
        return solids_tank_maintenance_labor_cost/(365 * 24) # USD/hr
        #Macerator in the homogenizer functions as a pump. All the associated cost
        #should refere to sanunits.G2RThomogenizer
        
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return 0.0 #kW


#%%
vr_drying_tunnel_path = ospath.join(g2rt_su_data_path, '_vr_dryingtunnel.csv')
@price_ratio()
class VRdryingtunnel(SanUnit):
    '''
    This drying tunnel unit is used to produce solids cakes in volume reduction 
    generation II reinveted toilet [1]. The heat source is from a heating coil.
    #TODO: consider installing a heat exchange to receive heat from pasteurization.
    The following components should be included in system thermo object for simulation:
    H2O, N, CH4, N2O.
    
    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, ConveyorBelt, Fan, Polyethylene.

    Parameters
    ----------
    ins : Iterable(stream)
        condensed effluent from concentrator, filter press cakes.
    outs : Iterable(stream)
        solid cakes, fugitive N2O, fugitive CH4.
    moisture_content_out : float
        Desired moisture content of the Condensed effluent.

    Warnings
    --------
    Energy balance is not performed for this unit.

    References
    ----------
    [1] YEE et al. VOLUME REDUCTION SOLIDS TREATMENT SYSTEM. 
       https://patents.google.com/patent/WO2023288327A1/en?oq=WO2023288327A1

    See Also
    --------
    :class:`~.sanunits.BiogenicRefineryHHXdryer`
    :class:`~.sanunits.BiogenicRefineryHHX`
    '''
    
    _ins_size_is_fixed = False
    _N_outs = 5
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',ppl=6,
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
                
        data = load_data(path=vr_drying_tunnel_path)
        self.ppl = ppl
        
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids)) + ('OtherSS',)
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])

    def _run(self):
        mixture = WasteStream()
        mixture.mix_from(self.ins)
        # condensate, press_cakes = self.ins
        solid_cakes, N2O, CH4, NH3_gas, water_vapor  = self.outs
        solid_cakes.copy_like(mixture)
        solubles, solids = self.solubles, self.solids
        N2O.phase = CH4.phase = NH3_gas.phase = water_vapor.phase
        #Calculate water and solids in the solids cake
        mc_in = solid_cakes.imass['H2O'] / solid_cakes.F_mass # fraction
        mc_out = self.moisture_content_out/100
        if mc_in < mc_out*0.999:
            raise RuntimeError(f'Moisture content of the influent stream ({mc_in:.2f}) '
                               f'is smaller than the desired moisture content ({mc_out:.2f}).')
        TS_in = solid_cakes.imass[solids].sum()# kg TS dry/hr
        solid_cakes.imass['H2O'] = TS_in/(1-mc_out)*mc_out #kg water/hr
        
        #Calculate N2O and CH4 emissions
        CH4.imass['CH4'] = drying_CH4_to_air = \
            self.drying_CH4_emissions * self.carbon_COD_ratio * \
            solid_cakes.COD * solid_cakes.F_vol / 1000 # kg CH4 /hr
        drying_NH3_to_air = self.drying_NH3_emissions * solid_cakes.imass['NH3'] # kg NH3 /hr
        N2O.imass['N2O'] = drying_NH3_to_air * self.NH3_to_N2O # kg N2O /hr
        solid_cakes.imass['NH3'] -= drying_NH3_to_air
        NH3_gas.imass['NH3'] = drying_NH3_to_air # kg NH3 /hr
        water_vapor.imass['H2O'] = (mixture.imass['H2O'] - solid_cakes.imass['H2O']) #kg H2O/hr
        # Store the calculated value of water_vapor.imass['H2O'] for use in the _cost function
        self.water_vapor_H2O = water_vapor.imass['H2O']  #kg H2O/hr
        
        #Calculate COD
        drying_CO2_to_air = (self.drying_CO2_emissions * self.carbon_COD_ratio
                             * solid_cakes.COD * solid_cakes.F_vol / 1000) # kg CO2 /hr
        # 44/12/16 are the molecular weights of CO2, C, and CH4, respectively
        solid_cakes.imass['sCOD'] -=  (drying_CO2_to_air/44*12+drying_CH4_to_air/16*12) / self.carbon_COD_ratio
    
    def _init_lca(self):
        self.construction = [
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            Construction('conveyor_belt',linked_unit=self, 
                         item='ConveyorBelt', 
                         quantity_unit='m'),
            Construction('fan',linked_unit=self, 
                         item='Fan',
                         quantity_unit='kg'),
            Construction('polyethylene',linked_unit=self,
                         item='Polyethylene', 
                         quantity_unit='kg'),
            ]
        
    def _design(self):
        design=self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        design['ConveyorBelt'] = constr[1].quantity = 0.07 #0.07 m equivalent to 1m based on 0.2m equivalent to 3m width in ecoinvent
        design['Fan'] = constr[2].quantity = 3 # a fan roughly weigh 3 kg
        design['Polyethylene'] = constr[3].quantity = self.polyethylene_weight
        
    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        #C['Stainless steel'] = self.stainless_steel_cost * D['StainlessSteel']
        C['Housing'] = self.housing_cost
        C['AirDuct'] = self.duct_cost + self.ventilator_cost
        C['Conveyor'] = self.conveyor_cost
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Housing'] +
                               C['AirDuct']+
                               C['Conveyor'])
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio
        
        self.power_utility(self.water_vapor_H2O * self.energy_required_to_evaporize_water + 
                            self.conveyor_power_demand * self.conveyor_daily_operation*self.ppl/24
                            ) # kW
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr

    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.drying_tunnel_maintenance * self.wages)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.water_vapor_H2O * self.energy_required_to_evaporize_water + 
                self.conveyor_power_demand * self.conveyor_daily_operation*self.ppl/24) #kW
#%%
g2rt_housing_path = ospath.join(g2rt_su_data_path, '_g2rt_housing.csv')

@price_ratio()
class G2RTHousing(Copier):
    '''
    Housing of the generation II reinveted toilet.
    
    This is a non-reactive unit (i.e., the effluent is copied from the influent).
    
    The following impact items should be pre-constructed for life cycle assessment:
    Aluminum, Steel, ZincCoat, Polyethylene
    
    References
    ----------
    [1] Watabe et al. Advancing the Economic and Environmental Sustainability 
    of the NEWgenerator Nonsewered Sanitation System." ACS Environmental Au 
    3.4 (2023): 209-222.
    https://pubs.acs.org/doi/10.1021/acsenvironau.3c00001
    
    See Also
    --------
    :class:`~.sanunits.NEWgeneratorControls`
    '''
    def __init__(self, ID='', ins=None, outs=(),  thermo=None, init_with='WasteStream',
                 **kwargs):
        data = load_data(path=g2rt_housing_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        Copier.__init__(self, ID, ins, outs, thermo, init_with)
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    def _init_lca(self):
        self.construction = [
            Construction(item='Steel', linked_unit=self, quantity_unit='kg'),
            Construction(item='Aluminum', linked_unit=self, quantity_unit='kg'),
            Construction(item='ZincCoat', linked_unit=self, quantity_unit='m2'),
            Construction(item='Polyethylene', linked_unit=self, quantity_unit='kg'),
            ]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.housing_steel_weight
        design['Aluminum'] = constr[1].quantity = self.housing_zinc_coat
        design['ZincCoat'] = constr[2].quantity = self.housing_aluminum_weight
        design['Polyethylene'] = constr[3].quantity = self.housing_PE_weight
        self.add_construction(add_cost=False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Housing_frame"] = self.system_housing
        C["Misc.parts"] = self.miscellaneous_cost_ratio* C["Housing_frame"]
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
    
    @property
    def OPEX(self):
        return 0.0 #USD/day
    
    @property
    def labor_expense(self):
        return 0.0 #USD/day
    
    @property
    def power_kW(self):
        return 0.0 #kW

#%%
g2rt_controls_path = ospath.join(g2rt_su_data_path, '_g2rt_controls.csv')

@price_ratio()
class G2RTControls(Copier):
    '''
    Electronic control of the generation II reinveted toilet that interface with sensors,
    valves, pumps, and motors.
    
    This is a non-reactive unit (i.e., the effluent is copied from the influent).
    
    The following impact items should be pre-constructed for life cycle assessment:
    ControlUnits, Polycarbonate, Aluminum, ElectricCables, ElectronicsPassive, ElectronicsActive
    
    References
    ----------
    [1] Watabe et al. Advancing the Economic and Environmental Sustainability 
    of the NEWgenerator Nonsewered Sanitation System." ACS Environmental Au 
    3.4 (2023): 209-222.
    https://pubs.acs.org/doi/10.1021/acsenvironau.3c00001
    
    See Also
    --------
    :class:`~.sanunits.NEWgeneratorControls`
    '''
    def __init__(self, ID='', ins=None, outs=(),  thermo=None, init_with='WasteStream',
                 **kwargs):
        data = load_data(path=g2rt_controls_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        Copier.__init__(self, ID, ins, outs, thermo, init_with)

        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    def _init_lca(self):
        self.construction = [
            Construction(item='ControlUnits', linked_unit=self, quantity_unit='kg',
                         lifetime=self.control_system_PLC_lifetime),
            Construction(item='Polycarbonate', linked_unit=self, quantity_unit='kg'),
            Construction(item='Aluminum', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectricCables', linked_unit=self, quantity_unit='m'),
            Construction(item='ElectronicsPassive', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectronicsActive', linked_unit=self, quantity_unit='kg'),
            ]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['ControlUnits'] = constr[0].quantity = self.control_control_units_weight
        design['Polycarbonate'] = constr[1].quantity = self.control_polycarbonate_weight
        design['Aluminum'] = constr[2].quantity = self.control_aluminum_weight
        design['ElectricCables'] = constr[3].quantity = self.control_cable_length
        design['ElectronicsPassive'] = constr[4].quantity = self.control_electronics_passive_weight
        design['ElectronicsActive'] = constr[5].quantity = self.control_electronics_active_weight
        self.add_construction(add_cost=False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Controller"] = self.programmable_logic_controller_cost * self.PLC_quantities
        C["IO_relay_modules"] = self.IO_relay_module_cost * self.relay_module_quantity
        C["Sensors"] = self.sensors_cost
        C["Cables"] = self.cables
        C["Misc.parts"] = self.miscellaneous_cost_ratio*(C["Controller"]+
                                                         C["IO_relay_modules"]+
                                                         C["Sensors"]+
                                                         C["Cables"]
                                                         )
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        power_demand = (
            self.control_system_energy_percycle * self.control_batch_cycle_perday +
            self.control_background_runtime_energy_perday
            )
        power_demand = power_demand / 24  # convert from kWh/d to kW
        self.power_utility(power_demand) # kW
        
        self.add_OPEX =  self._calc_replacement_cost() + self._calc_maintenance_labor_cost()
        
    def _calc_replacement_cost(self):
        control_system_replacement_cost = (
            self.PLC_quantities + self.sensors_cost / self.control_sensor_lifetime
            )
        return control_system_replacement_cost / (365 * 24) * self.price_ratio # USD/hr
    
    def _calc_maintenance_labor_cost(self):
        control_system_maintenance_labor = self.control_labor_replacement_misc_repairs * self.wages
        return control_system_maintenance_labor / (365 * 24) # USD/hr
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.control_system_ORP_energy_percycle * self.control_batch_cycle_perday +
            self.control_background_runtime_energy_perday)/24 #kW

#%%
g2rt_solids_separation_path = ospath.join(g2rt_su_data_path, '_g2rt_solids_separation.csv')

@price_ratio()
class G2RTSolidsSeparation(SanUnit):
    '''
    Solids separation unit in generation II reinveted toilets as a frontend separator [1].
    
    .. note:

    Non-reactive. Moisture content of the effluent solid is adjusted to be 99% [2].

    The following components should be included in system thermo object for simulation:
    H2O

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, Polyethylene, ElectricMotor, Pump

    Parameters
    ----------
    ins : Iterable(stream)
        Influent stream.
    outs : Iterable(stream)
        liquid stream and solid stream after separation.
    moisture_content_out : float
        Moisture content of the effluent solids stream.

    References
    ----------
    [1] YEE et al. Buffer tank separation and homogenization system. 
    https://patents.google.com/patent/WO2023288114A1/en?oq=WO2023288114A1
    [2] YEE et al. Volume reduction non-sewered single unit toilet system.
    https://patents.google.com/patent/WO2023288326A1/en?oq=WO2023288326A1
    
    See Also
    ---------
    :class:`~.sanunits.BiogenicRefineryGrinder`
    '''
    _N_ins = 2
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', ppl=6,
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.ppl = ppl

        data = load_data(path=g2rt_solids_separation_path)
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
        self.construction = [Construction("stainless_steel", linked_unit=self,
                                          item = "StainlessSteel", 
                                          quantity_unit= "kg"),
                             Construction("polyethylene", linked_unit=self,
                                          item = "Polyethylene",
                                          quantity_unit= "kg"),
                             Construction("electric_motor", linked_unit=self,
                                          item = "ElectricMotor",
                                          quantity_unit= "kg"),
                             Construction("pump", linked_unit=self,
                                          item = "Pump",
                                          quantity_unit= "ea"),
                             ]
    
    def _run(self):
        waste_in,flushing_water = self.ins
        # waste_in.mix_from(self.ins)
        # print(waste_in.imass['H2O']) #TODO:debug
        liquid_stream, solid_stream = self.outs
        solubles, solids = self.solubles, self.solids
        solid_stream.mix_from(self.ins)
        solid_stream.imass[solids] = solid_stream.imass[solids] * self.solids_separator_TSS_removal/100
        mc_in = solid_stream.imass['H2O'] / (waste_in.F_mass + flushing_water.F_mass) # fraction
        mc_out = self.moisture_content_out/100 #convert to fraction
        TS_in = waste_in.imass[solids].sum()+ flushing_water.imass[solids].sum()# kg TS dry/hr
        TS_out = solid_stream.imass[solids].sum()
        if mc_in < mc_out*0.999:
            mc_out = mc_in
        # # solid_stream.copy_flow(waste_in,solids) #all solids go to sludge, remove from waste_in
        # solid_stream.imass[solids] = solid_stream.imass[solids] * self.solids_separator_TSS_removal/100 #add the removed solids back

        # if mc_in < mc_out*0.999:
        #     print(f"Moisture content of the influent stream ({mc_in:.4f})"
        #           f"is smaller than that of the desired effluent stream ({mc_out:.4f})."
        #           "High solids and low flushing event detected, adding more flushing water!")
        #     solid_stream.imass['H2O'] = (waste_in.F_mass + flushing_water.F_mass) * mc_out
        #     # raise RuntimeError(f'Moisture content of the influent stream ({mc_in:.4f}) '
        #     #                    f'is smaller than the desired moisture content ({mc_out:.4f}).')
        
        # TS_in = waste_in.imass[solids].sum()+ flushing_water.imass[solids].sum()# kg TS dry/hr
        # TS_out = solid_stream.imass[solids].sum()
        
        #calculate water in the solid cakes
        solid_stream.imass['H2O'] = TS_out/(1-mc_out)*mc_out

        solid_stream.imass[solubles] = (waste_in.imass[solubles]+flushing_water.imass[solubles])*\
            (TS_out/(1-mc_out)-TS_out)/(waste_in.F_mass+ flushing_water.F_mass-TS_in)
        liquid_stream.mass = flushing_water.mass + waste_in.mass-solid_stream.mass
        
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        design['Polyethylene'] = constr[1].quantity = self.polyethylene_weight
        design['ElectricMotor'] = constr[2].quantity = self.actuator_weight
        design['Pump'] = constr[1].quantity = 1.
        self.add_construction(add_cost = False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Vacuum tank"] = (self.vacuum_tank_cost + 
                            self.separating_filter_cost +
                            self.solid_outlet_cost + 
                            self.liquid_outlet_cost+
                            self.inlet_chamber_cost)
        C["Vacuum pump"] = self.vacuum_pump_cost
        C["Valves"] = self.valve_cost * self.valve_quantity
        C["Actuator"] = self.actuator_cost * self.actuator_quantity
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.power_utility(self.vacuum_pump_power_demand * 
                           self.vacuum_pump_operation*
                           self.ppl/24) # kW
        
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
    
    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.solids_separator_maintenance * self.wages)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.vacuum_pump_power_demand * 
                           self.vacuum_pump_operation*
                           self.ppl/24) #kW

#%%
g2rt_belt_separation_path = ospath.join(g2rt_su_data_path, '_g2rt_belt_separation.csv')

@price_ratio()
class G2RTBeltSeparation(SanUnit):
    '''
    Belt separation unit in generation II reinveted toilets before the buffer tank [1].
    
    .. note:

    Non-reactive. Moisture content of the effluent solid is adjusted to be 96-99% [2].
    The solids in the liquid stream influent is partially transferred to the solids effluent.

    The following components should be included in system thermo object for simulation:
    H2O, OtherSS.

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, ConveyorBelt, Polyethylene.

    Parameters
    ----------
    ins : Iterable(stream)
        liquid stream and solid stream from the solid separator, and retentate 
        from the ultrafiltration unit
    outs : Iterable(stream)
        liquid stream and solid stream are copied from the influent.
    moisture_content_out : float
        Moisture content of the effluent solids stream.

    References
    ----------
    [1] YEE et al. Buffer tank separation and homogenization system. 
    https://patents.google.com/patent/WO2023288114A1/en?oq=WO2023288114A1
    [2] YEE et al. Volume reduction non-sewered single unit toilet system.
    https://patents.google.com/patent/WO2023288326A1/en?oq=WO2023288326A1
    
    See Also
    ---------
    :class:`~.sanunits.BiogenicRefineryGrinder`
    '''
    _N_ins = 2
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', ppl=6,
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1) 
        self.ppl = ppl

        data = load_data(path=g2rt_belt_separation_path)
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
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            Construction('conveyor_belt',linked_unit=self, 
                         item='ConveyorBelt', 
                         quantity_unit='m'),
            Construction('fan',linked_unit=self, 
                         item='Fan',
                         quantity_unit='kg'),
            Construction('polyethylene',linked_unit=self,
                         item='Polyethylene', 
                         quantity_unit='kg'),
            ]
    
    def _run(self):
        liquid_in, solid_in = self.ins
        # uf_retentate.show() #TODO:debug
        liquid_out, solid_out = self.outs
        solid_out.phase = 's'
        solubles, solids = self.solubles, self.solids
        solid_out.mix_from(self.ins)
        TS_in = solid_out.imass[solids].sum()
        #solid_stream.copy_flow(solid_stream_in,solids) #all solids in go to solids out
        solid_out.imass[solids] = solid_in.imass[solids] + liquid_in.imass[solids] * self.belt_separator_TSS_removal/100 
        # the removed solids
        TS_out = solid_out.imass[solids].sum() # kg TS dry/hr
        
        # mc_in = solid_stream.imass['H2O'] / (solid_stream_in.F_mass + 
        #                                      liquid_stream_in.F_mass +
        #                                      uf_retentate.F_mass) # fraction
        mc_out = self.moisture_content_out/100 #convert to fraction; 
        mc_in = solid_out.imass['H2O'] / (liquid_in.F_mass + solid_in.F_mass)
        if mc_in < mc_out*0.999:
            mc_out = mc_in
        # if mc_in < mc_out*0.999:
        #     raise RuntimeError(f'Moisture content of the influent stream ({mc_in:.4f}) '
        #                        f'is smaller than the desired moisture content ({mc_out:.4f}).')
        # mix_in = WasteStream(ID='mix_in')
        # mix_in.mix_from((liquid_stream_in, uf_retentate,solid_stream_in))
        
        # TS_in = (solid_stream_in.imass[solids].sum()+ 
        #          liquid_stream_in.imass[solids].sum()+
        #          uf_retentate.imass[solids].sum()
        #          ) # kg TS dry/hr

        #calculate water and solid COD in the solid cakes
        solid_out.imass['H2O'] = TS_out/(1-mc_out)*mc_out

        # solid_stream._COD = (solid_stream_in.COD * solid_stream_in.F_vol +
        #                      (liquid_stream_in.COD * liquid_stream_in.F_vol 
        #                       + uf_retentate.COD * uf_retentate.F_vol) *
        #                      self.belt_separator_TSS_removal/100)/ solid_stream.F_vol

        solid_out.imass[solubles] = (liquid_in.imass[solubles] + 
                                        solid_in.imass[solubles]) *\
            (TS_out/(1-mc_out)-TS_out)/((solid_in.F_mass + liquid_in.F_mass)-TS_in)
        liquid_out.mass = solid_in.mass + liquid_in.mass -solid_out.mass

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        design['ConveyorBelt'] = constr[1].quantity = 0.07 #~0.07m equivalent to 0.7m based on 0.3m equivalent to 3m width in ecoinvent
        design['Polyethylene'] = constr[2].quantity = self.polyethylene_weight
        self.add_construction(add_cost=False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Belt_separator"] = (self.belt_conveyor_cost+
                               self.solid_inlet_chamber_cost+
                               self.liquid_inlet_chamber_cost+
                               self.housing_cost+
                               self.squeegee_cost+
                               self.splash_shield_cost
                               )
        C["Misc.parts"] = self.miscellaneous_cost_ratio * C["Belt_separator"]
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        power_demand = (
            self.belt_daily_operation * self.belt_power_demand * self.ppl/24
            ) # convert from kWh/d to kW
        self.power_utility(power_demand) # kW
                
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = (total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr

    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.belt_separator_maintenance * self.wages)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.belt_daily_operation * self.belt_power_demand * self.ppl/24) #kW

#%%

ultrafiltration_path = ospath.join(g2rt_su_data_path, '_g2rt_ultrafiltration.csv')

@price_ratio()
class G2RTUltrafiltration(SanUnit):
    '''
    Ultrafiltration in the generation II reinveted toilets is used for removing suspended solids
    with automated backwash.
    
    Modified from ultrafiltration unit in Duke Reclaimer system.

    The following impact items should be pre-constructed for life cycle assessment:
    GFRPlastic, Steel, NylonGlassFilled, CastingBrass

    Parameters
    ----------
    ins : Iterable(stream)
        waste: liquid waste stream to be treated by ultrafiltration unit.
    outs : Iterable(stream)
        treated: treated liquid leaving ultrafiltration unit.
        retentate: concentrated retentate leaving ultrafiltration unit.
    ppl: int
        Total number of users for scaling of costs.

    References
    ----------
    [1] Trotochaud et al., Laboratory Demonstration and Preliminary Techno-Economic Analysis of an Onsite
    Wastewater Treatment System Environ. Sci. Technol. 2020, 54, (24), 16147–16155.
    https://dx.doi.org/10.1021/acs.est.0c02755
    
    [2] Duke Center for WaSH-AID Reclaimer design team data and guidance
    https://washaid.pratt.duke.edu/work/water-sanitation/reinvent-toilet-challenge
    
    See Also
    ---------
    :class:`~.sanunits.ReclaimerUltrafiltration`
    
    '''
    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        data = load_data(path=ultrafiltration_path)
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
            Construction(item='GFRPlastic', linked_unit=self, quantity_unit='kg'),
            Construction(item='Steel', linked_unit=self, quantity_unit='kg'),
            Construction(item='NylonGlassFilled', linked_unit=self, quantity_unit='kg'),
            Construction(item='CastingBrass', linked_unit=self, quantity_unit='kg'),
            ]

    def _run(self):
        waste_in = self.ins[0]
        liquid_stream, solid_stream = self.outs
        solubles, solids = self.solubles, self.solids
        solid_stream.copy_flow(waste_in,solids) #all solids go to sludge, remove from waste_in
        liquid_stream.copy_flow(waste_in, solubles) #only copy solubles 
        # TS_in = waste_in.imass[solids].sum() # kg TS dry/hr
        # TS_out = solid_stream.imass[solids].sum()
        solid_stream.imass['H2O'] = waste_in.imass['H2O']*(1-self.water_recovery_rate/100)
        solid_stream.imass[solids] = waste_in.imass[solids] * self.TSS_removal/100 # the removed solids
        solid_stream.imass[solubles] = waste_in.imass[solubles]*\
            solid_stream.imass['H2O']/waste_in.imass['H2O']
        solid_stream.imass['sCOD'] += waste_in.imass['sCOD']*self.sCOD_removal/100
        liquid_stream.mass = waste_in.mass-solid_stream.mass
        # print(f"The recycled UF water flow is {solid_stream.imass['H2O']} kg/h.")

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['GFRPlastic'] = constr[0].quantity = self.Plastic_weight
        design['Steel'] = constr[1].quantity = self.Steel_weight
        design['NylonGlassFilled'] = constr[2].quantity = self.Nylon_weight
        design['CastingBrass'] = constr[3].quantity = self.Brass_weight
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs 
        C['Pipes'] = self.one_in_pipe_SCH40 + self.onehalf_in_pipe_SCH40 + self.three_in_pipe_SCH80
        C['Fittings'] = (
            self.one_in_elbow_SCH80 +
            self.one_in_tee_SCH80 +
            self.one_in_SCH80 +
            self.one_onehalf_in_SCH80 +
            self.onehalf_in_SCH80 +
            self.three_in_SCH80_endcap +
            self.one_one_NB_MTA +
            self.one_onehalf_NB_MTA +
            self.foot_valve +
            self.one_onehalf_in_SCH80_threadedtee +
            self.three_in_pipe_clamp +
            self.one_in_pipe_clamp +
            self.onehalf_in_pipe_clamp +
            self.two_way_valve +
            self.UF_brush
            )
        C['UF_unit'] = self.UF_unit
        C['Backwash_parts'] = self.air_scrubbing_blower_cost + self.back_flush_pump_cost
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio/3 #scaled down from 30 users to 6 users per day

        self.add_OPEX = self._calc_replacement_cost()
        # [W][1 kW/1000 W][hr/d][1 d/ 24 h] = [kW]
        self.back_wash_cycles = ceil(self.ins[0].get_TSS()*self.ins[0].F_vol/1000*24/self.accumulated_TSS_triggering_cleaning)
        power_demand = (self.ultrafiltration_energy_consumption * self.outs[0].F_mass+
                        (self.backwash_pump_energy_percycle+
                         self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24)  #kW
        self.power_utility(power_demand) #scaled down from 30 users to 6 users per day

    def _calc_replacement_cost(self):
        pipe_replacement_cost = (
            self.one_in_pipe_SCH40 / self.one_in_pipe_SCH40_lifetime +
            self.onehalf_in_pipe_SCH40 / self.onehalf_in_pipe_SCH40_lifetime +
            self.three_in_pipe_SCH80 / self.three_in_pipe_SCH80_lifetime
            )
        
        fittings_replacement_cost = (
            self.one_in_elbow_SCH80 / self.one_in_elbow_SCH80_lifetime +
            self.one_in_tee_SCH80 / self.one_in_tee_SCH80_lifetime +
            self.one_in_SCH80 / self.one_in_SCH80_lifetime +
            self.one_onehalf_in_SCH80 / self.one_onehalf_in_SCH80_lifetime +
            self.onehalf_in_SCH80 / self.onehalf_in_SCH80_lifetime +
            self.three_in_SCH80_endcap / self.three_in_SCH80_endcap_lifetime +
            self.one_one_NB_MTA / self.one_one_NB_MTA_lifetime +
            self.one_onehalf_NB_MTA / self.one_onehalf_NB_MTA_lifetime +
            self.foot_valve / self.foot_valve_lifetime +
            self.one_onehalf_in_SCH80_threadedtee / self.one_onehalf_in_SCH80_threadedtee_lifetime +
            self.three_in_pipe_clamp / self.three_in_pipe_clamp_lifetime +
            self.one_in_pipe_clamp / self.one_in_pipe_clamp_lifetime +
            self.onehalf_in_pipe_clamp / self.onehalf_in_pipe_clamp_lifetime +
            self.two_way_valve / self.two_way_valve_lifetime +
            self.UF_brush / self.UF_brush_lifetime
            )

        uf_replacement_cost = self.UF_unit/ self.UF_unit_lifetime 

        total_replacement_cost = self.price_ratio * (pipe_replacement_cost + fittings_replacement_cost + uf_replacement_cost)/3 # USD/year
        return total_replacement_cost / (365 * 24)  # USD/hr
    
    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.ultrafiltration_maintenance * self.wages)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.ultrafiltration_energy_consumption * self.outs[0].F_mass+
                        (self.backwash_pump_energy_percycle+
                         self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24) #kW

#%%
reverse_osmosis_path = ospath.join(g2rt_su_data_path, '_g2rt_reverse_osmosis.csv')

@price_ratio()
class G2RTReverseOsmosis(SanUnit):
    '''
    A reverse osmosis unit process to recover water and concentrate liquid stream.
    The model is based on a fraction of water recovered.
    
    Parameters
    ----------
    ins : 
        Inlet fluid to be split.
    outs : 
        * [0] Permeate
        * [1] Brine
    water_recovery : float, optional, 0 to 1
        Water recovered to 0th stream. Defaults to 0.6
    TDS_removal: float, optional, 0 to 1
        rejection rate to total dissolved salts. Defaults to 0.95
    permeate_recycle_ratio: float, optional, 0 to 1
        fraction of permeate that recirculates to front end as flushing water. Defaults to 1.
        
    The following impact items should be pre-constructed for life cycle assessment:
    GFRPlastic, Steel, ReverseOsmosisModule
    '''
    _N_ins = 1
    _N_outs = 3

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 water_recovery=0.6,TDS_removal = 0.95, permeate_recycle_ratio =1,**kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.water_recovery = water_recovery
        self.TDS_removal = TDS_removal
        self.permeate_recycle_ratio = permeate_recycle_ratio

        data = load_data(path=reverse_osmosis_path)
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
        waste_in = self.ins[0]
        permeate_recycle, brine, permeate_discharge = self.outs
        solubles, solids = self.solubles, self.solids
        permeate = WasteStream()
        for i in solubles:
            permeate.imass[i] = waste_in.imass[i]

        # permeate.copy_flow(waste_in, solubles) #copy_flow removes stream_impact_item 
        brine.copy_like(waste_in) #copy all the solids and solubles
        
        permeate.imass['H2O'] = waste_in.imass['H2O'] * self.water_recovery
        permeate.imass[solubles] = waste_in.imass[solubles]*(1-self.TDS_removal)
        permeate.imass['NH3'] = waste_in.imass['NH3']*(1-self.NH3_removal)
        permeate.split_to(permeate_recycle, permeate_discharge, self.permeate_recycle_ratio)
        # permeate.imass['P'] = waste_in.imass['P'] * 0.01 #higher rejection for P
        brine.mass = waste_in.mass - permeate.mass

    def _init_lca(self):
        self.construction = [
            Construction(item='GFRPlastic', linked_unit=self, quantity_unit='kg'),
            Construction(item='Steel', linked_unit=self, quantity_unit='kg'),
            Construction(item='ReverseOsmosisModule', linked_unit=self, quantity_unit='m2'),
            ]
        
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['GFRPlastic'] = constr[0].quantity = self.GFRPlastic_weight
        design['Steel'] = constr[1].quantity = self.Steel_weight
        design['ReverseOsmosisModule'] = constr[2].quantity = self.membrane_area
        self.add_construction(add_cost=False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Pipes'] = self.piping_cost
        C['RO_system'] = self.reverse_osmosis_system_cost
        C['Backwash_parts'] = self.air_scrubbing_blower_cost + self.back_flush_pump_cost
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        self.back_wash_cycles = ceil(self.ins[0].COD*self.ins[0].F_vol/1000*24/
                                     self.accumulated_COD_triggering_cleaning/self.membrane_area)
        
        power_demand = (self.RO_energy_consumption * self.outs[0].F_mass+
                        (self.backwash_pump_energy_percycle+
                         self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24)
        
        self.power_utility(power_demand) #kW
        
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
                  
        self.add_OPEX = ((total_equipment*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 4% of CAPEX per year
                         self._calc_membrane_replacement_cost() +
                         self._calc_maintenance_labor_cost())) #USD/hr
    def _calc_membrane_replacement_cost(self): #USD/hr
 
        membrane_replacement_cost = self.membrane_cost / self.membrane_life_time #USD/yr
        return membrane_replacement_cost/(365*24) #USD/hr
    
    def _calc_maintenance_labor_cost(self): #USD/hr
        maintenance_labor_cost= (self.reverse_osmosis_maintenance * self.wages)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.RO_energy_consumption * self.outs[0].F_mass+
                        (self.backwash_pump_energy_percycle+
                         self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24) #kW
        