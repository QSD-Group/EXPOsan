#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Yuyao Huang <yuyaoh20@gmail.com>
    Siqi Tang <siqit@outlook.com>

Part of this module is based on the qsdsan package:
https://github.com/QSD-Group/QSDsan/blob/main/qsdsan/sanunits/_pumping.py#L18

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.

Note: This module ignored the pump head, but there is an example in 'SelfPrimingPump'.
'''
# Only considering materials, cost, pump pressure, energy consumption, pump efficiency in this module because of a small system

import numpy as np
import os, pandas as pd
from qsdsan.sanunits import Pump
from qsdsan import Construction
from qsdsan.utils import price_ratio

__all__ = ('LiftPump', 'AgitationPump', 'DosingPump', 'ReturnPump', 'SelfPrimingPump', 'AirDissolvedPump', 'MicroBubblePump', 'ClearWaterPump')


# %%

@price_ratio()
class LiftPump(Pump):
    """
    LiftPump: Specialized pump for lifting waste water from collection tank to primary clarifiers
    """

    _N_ins = 1  # Number of input streams
    _N_outs = 1  # Number of output streams
    _ins_size_is_fixed = True  # Let the input interface be fixed
    _outs_size_is_fixed = True  # Let the output interface be fixed
    _pump_flow_rate = 250 * 60 / 1000  # The flow rate of each pump [m³/h]
    _CastIron_weight_per_pump = 12  # The weight of cast iron [kg]
    _pump_power = 0.45  # The power of each pump [kW]
    exponent_scale = 0.4

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', material='Cast iron',
                 F_BM_default=3.3, pump_type='Default', isdynamic=False, ignore_NPSH=True,
                 dP_design=None,  # The extra pressure of the pump [Pa]
                 P=101325,  # The pressure of the pump [Pa]
                 life_time=None, pump_cost=None, working_factor=None,
                 operation_time=None, ppl=100, baseline_ppl=100):
                 
        # Initialize the LiftPump with parent class logic
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with, F_BM_default=F_BM_default,
                        isdynamic=isdynamic, ignore_NPSH=ignore_NPSH, P=P, dP_design=dP_design, pump_type=pump_type,
                        material=material)

        self.life_time = life_time  # The life time of the pump [yrs]
        self.pump_cost = pump_cost  # The price of each pump [$/pump]
        self.working_factor = working_factor  # The working factor of the pump considering the head.
        self.operation_time = operation_time  # The operation time of the pump [h/d]
        self.ppl = ppl  # The number of all people used all toilets
        self.baseline_ppl = baseline_ppl  # The number of people per toilet
        
    
    @property
    def dP_factor(self):
        dP_factor = (self.P + self.dP_design) / self.P
        return dP_factor

    @property
    def actual_power(self):
        # Calculate the actual total power [kW]
        return self._pump_power / self.working_factor * self.dP_factor

    def _init_lca(self):
        self.construction = [Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Cast_iron'] = constr[0].quantity = self._CastIron_weight_per_pump
        self.add_construction(add_cost=False)
    
    def _run(self):
        
        # Input stream
        Water_in = self.ins[0]
        
        # Output stream
        Water_out = self.outs[0]
        
        # Inherite input stream
        Water_out.copy_like(Water_in)

    def _cost(self):

        C = self.baseline_purchase_costs
        C['Lift Pump'] = self.pump_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()

        power_demand = self.actual_power * self.operation_time
        self.power_utility(power_demand)
        
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        # Replacement cost: based on pump count
        replacement_cost = self.pump_cost / (self.life_time * 365 * 24) * scale
        return replacement_cost

# %%
 
@price_ratio()
class AgitationPump(Pump):
    """
    AgitationPump: Specialized pump for mixing chemicals with solutions or wastewater.
    """
    
    _N_ins = 1  # Number of input streams
    _N_outs = 1  # Number of output streams
    _ins_size_is_fixed = True  # Let the input interface be fixed
    _outs_size_is_fixed = True  # Let the output interface be fixed
    _CastIron_weight_per_pump = 9 # The weight of cast iron [kg]
    _pump_power = 0.18 # The power of each pump [kW]
    exponent_scale = 0.4
   
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream', material='Cast iron',
                 F_BM_default=3.3, pump_type='Default', isdynamic=False, ignore_NPSH=True,
                 dP_design=None,  # The extra pressure of the pump [Pa]
                 P=101325,  # The pressure of the pump [Pa]
                 life_time=None, pump_cost=None, working_factor=None,
                 operation_time=None, ppl=100, baseline_ppl=100):
                 
        # Initialize the LiftPump with parent class logic
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with, F_BM_default=F_BM_default,
                        isdynamic=isdynamic, ignore_NPSH=ignore_NPSH, P=P, dP_design=dP_design, pump_type=pump_type,
                        material=material)
        
        self.life_time = life_time  # The life time of the pump [yrs]
        self.pump_cost = pump_cost  # The price of each pump [$/pump]
        self.working_factor = working_factor  # The working factor of the pump considering the head.
        self.operation_time = operation_time  # The operation time of the pump [h/d] 
        self.ppl = ppl  # The number of all people used all toilets
        self.baseline_ppl = baseline_ppl  # The number of people per toilet

    @property
    def dP_factor(self):
        dP_factor = (self.P + self.dP_design) / self.P
        return dP_factor

    @property
    def actual_power(self):
        # Calculate the actual total power [kW]
        return self._pump_power / self.working_factor * self.dP_factor 
    
    def _init_lca(self):
        self.construction = [Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]
    
    def _run(self):
    
        # Input stream
        mass_in = self.ins[0]
        mass_in.phase = 'l'
        
        # Output stream
        mass_out = self.outs[0]
        
        # Inherite input stream
        mass_out.copy_like(mass_in)

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Cast_iron'] = constr[0].quantity = self._CastIron_weight_per_pump
        self.add_construction(add_cost=False)
    
    def _cost(self):

        C = self.baseline_purchase_costs
        C['Agitaion Pump'] = self.pump_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.actual_power * self.operation_time
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        # Replacement cost: based on pump count
        replacement_cost = self.pump_cost / (self.life_time * 365 * 24) *scale
        return replacement_cost

# %%

#@price_ratio
class DosingPump(Pump):
    """
    DosingPump: Specialized pump for input chemicals in water.
    """
    
    _N_ins = 1  # Number of input streams
    _N_outs = 1  # Number of output streams
    _ins_size_is_fixed = True  # Let the input interface be fixed
    _outs_size_is_fixed = True  # Let the output interface be fixed
    _CastIron_weight_per_pump = 2.5 # The weight of cast iron [kg]
    _pump_power = 0.5 # Assumption: The power of each pump [kW]
    exponent_scale = 0.4

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream', material='Cast iron',
                 F_BM_default=3.3, pump_type='Default', isdynamic=False, ignore_NPSH=True,
                 dP_design=None,  # The extra pressure of the pump [Pa]
                 P=101325,  # The pressure of the pump [Pa]
                 life_time=None, pump_cost=None, working_factor=None,
                 price_ratio=0.9,
                 operation_time=None, ppl=100, baseline_ppl=100):
                 
        # Initialize the LiftPump with parent class logic
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with, F_BM_default=F_BM_default,
                        isdynamic=isdynamic, ignore_NPSH=ignore_NPSH, P=P, dP_design=dP_design, pump_type=pump_type,
                        material=material)
        
        self.life_time = life_time  # The life time of the pump [yrs]
        self.pump_cost = pump_cost  # The price of each pump [$/pump]
        self.working_factor = working_factor  # The working factor of the pump considering the head.
        self.operation_time = operation_time  # The operation time of the pump [h/d]
        self.ppl = ppl  # The number of all people used all toilets
        self.baseline_ppl = baseline_ppl  # The number of people per toilet
        self.price_ratio = price_ratio

    @property
    def dP_factor(self):
        dP_factor = (self.P + self.dP_design) / self.P
        return dP_factor

    @property
    def actual_power(self):
        # Calculate the actual total power [kW]
        return self._pump_power / self.working_factor * self.dP_factor
    
    def _init_lca(self):
        self.construction = [Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]

    def _run(self):
        
        # Input stream
        mass_in = self.ins[0]
        #chemical_in = self.ins[1]
        
        # Output stream
        mass_out = self.outs[0]
        
        # Inherite input stream
        mass_out.copy_like(mass_in)

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Cast_iron'] = constr[0].quantity = self._CastIron_weight_per_pump
        self.add_construction(add_cost=False)
    
    def _cost(self):

        C = self.baseline_purchase_costs
        C['Dosing Pump'] = self.pump_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.actual_power * self.operation_time
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        # Replacement cost: based on pump count
        replacement_cost = self.pump_cost / (self.life_time * 365 * 24) *scale
        return replacement_cost

# %%

#@price_ratio
class ReturnPump(Pump):
    """
    ReturnPump: Specialized pump for returning sludge from membrane tank to primary clarifier 
    and anoxic tank or from primary clarifier to collection tank.
    """
    
    _N_ins = 1;  # Number of input streams
    _N_outs = 1;  # Number of output streams
    _ins_size_is_fixed = True;  # Let the input interface be fixed
    _outs_size_is_fixed = True;  # Let the output interface be fixed
    _CastIron_weight_per_pump = 18.2 # The weight of cast iron [kg]
    _pump_power = 0.45  # The power of each pump [kW]
    _pump_flow_rate = 200 * 60 / 1000  # The flow rate of each pump [m3/h]
    exponent_scale = 0.4
   
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', material='Cast iron',
                 F_BM_default=3.3, pump_type='Default', isdynamic=False, ignore_NPSH=True,
                 dP_design=None,  # The extra pressure of the pump [Pa]
                 P=101325,  # The pressure of the pump [Pa]
                 life_time=None, pump_cost=None, working_factor=None, 
                 price_ratio = 0.9, operation_time=None, ppl=100, baseline_ppl=100):
                 
        # Initialize the LiftPump with parent class logic
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with, F_BM_default=F_BM_default,
                        isdynamic=isdynamic, ignore_NPSH=ignore_NPSH, P=P, dP_design=dP_design, pump_type=pump_type,
                        material=material)
        
        self.life_time = life_time  # The life time of the pump [yrs]
        self.pump_cost = pump_cost  # The price of each pump [$/pump]
        self.working_factor = working_factor  # The working factor of the pump considering the head.
        self.operation_time = operation_time  # The operation time of the pump [h/d] 
        self.ppl = ppl  # The number of all people used all toilets
        self.baseline_ppl = baseline_ppl  # The number of people per toilet
        self.price_ratio = price_ratio

    @property
    def dP_factor(self):
        dP_factor = (self.P + self.dP_design) / self.P
        return dP_factor

    @property
    def actual_power(self):
        # Calculate the actual total power [kW]
        return self._pump_power / self.working_factor * self.dP_factor
    
    def _init_lca(self):
        self.construction = [Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Cast_iron'] = constr[0].quantity = self._CastIron_weight_per_pump
        self.add_construction(add_cost=False)
    
    def _run(self):
        
        # Input stream
        mass_in = self.ins[0]
        
        # Output stream
        mass_out = self.outs[0]
        
        # Inherite input stream
        mass_out.copy_like(mass_in)

    def _cost(self):

        C = self.baseline_purchase_costs
        C['Return Pump'] = self.pump_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.actual_power * self.operation_time
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        # Replacement cost: based on pump count
        replacement_cost = self.pump_cost / (self.life_time * 365 * 24) * scale
        return replacement_cost
            
# %%

#@price_ratio
class SelfPrimingPump(Pump):
    """
    SelfPrimingPump: Specialized pump for inputing water into the clear water tank
    """
    
    _N_ins = 1;  # Number of input streams
    _N_outs = 1;  # Number of output streams
    _ins_size_is_fixed = True;  # Let the input interface be fixed
    _outs_size_is_fixed = True;  # Let the output interface be fixed
    _CastIron_weight_per_pump = 9 # The weight of cast iron [kg]
    _pump_power = 0.6  # The power of each pump [kW]
    _pump_flow_rate = 30 * 60 / 1000  # The flow rate of each pump [m3/h]
    exponent_scale = 0.4

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', material='Cast iron',
                 F_BM_default=3.3, pump_type='Default', isdynamic=False, ignore_NPSH=True,
                 dP_design=None,  # The extra pressure of the pump [Pa]
                 P=101325,  # The pressure of the pump [Pa]
                 price_ratio=0.9,
                 life_time=None, pump_cost=None, working_factor=None,
                 operation_time=None, ppl=100, baseline_ppl=100):
                 
        # Initialize the LiftPump with parent class logic
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with, F_BM_default=F_BM_default,
                        isdynamic=isdynamic, ignore_NPSH=ignore_NPSH, P=P, dP_design=dP_design, pump_type=pump_type,
                        material=material)
        
        self.life_time = life_time  # The life time of the pump [yrs]
        self.pump_cost = pump_cost  # The price of each pump [$/pump]
        self.working_factor = working_factor  # The working factor of the pump considering the head.
        self.operation_time = operation_time  # The operation time of the pump [h/d] 
        self.ppl = ppl  # The number of all people used all toilets
        self.baseline_ppl = baseline_ppl  # The number of people per toilet
        self.price_ratio = price_ratio

    @property
    def dP_factor(self):
        dP_factor = (self.P + self.dP_design) / self.P
        return dP_factor

    @property
    def actual_power(self):
        # Calculate the actual total power [kW]
        return self._pump_power / self.working_factor * self.dP_factor 
    
    def _init_lca(self):
        self.construction = [Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Cast_iron'] = constr[0].quantity = self._CastIron_weight_per_pump
        self.add_construction(add_cost=False)

    def _run(self):
        
        # Input stream
        mass_in, = self.ins[0]
        
        # Output stream
        mass_out = self.outs[0]
        
        # Inherite input stream
        mass_out.copy_like(mass_in)
    
    def _cost(self):

        C = self.baseline_purchase_costs
        C['Self Priming Pump'] = self.pump_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.actual_power * self.operation_time
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        # Replacement cost: based on pump count
        replacement_cost = self.pump_cost / (self.life_time * 365 * 24) * scale
        return replacement_cost

# %%

#@price_ratio
class AirDissolvedPump(Pump):
    """
    In EnviroLoo system, Air-dissolved pump for extracting water from the clear water tank, injecting air into the water,
    and returning the water with dissolved air back to the clear water tank.
    """
    
    _N_ins = 2  # Number of input streams
    _N_outs = 1  # Number of output streams
    _ins_size_is_fixed = True  # Let the input interface be fixed
    _outs_size_is_fixed = True  # Let the output interface be fixed
    _CastIron_weight_per_pump = 6.8 # The weight of cast iron [kg]
    _pump_power = 0.55  # The power of each pump [kW]
    _pump_flow_rate = 30 * 60 / 1000  # The flow rate of each pump [m3/h]
    exponent_scale = 0.4
   
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', material='Cast iron',
                 F_BM_default=3.3, pump_type='Default', isdynamic=False, ignore_NPSH=True,
                 dP_design=50000,  # The extra pressure of the pump [Pa]
                 P=101325,  # The pressure of the pump [Pa]
                 price_ratio=0.9,
                 life_time=None, pump_cost=None, working_factor=None,
                 operation_time=None, ppl=100, baseline_ppl=100):
                 
        # Initialize the LiftPump with parent class logic
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with, F_BM_default=F_BM_default,
                        isdynamic=isdynamic, ignore_NPSH=ignore_NPSH, P=P, dP_design=dP_design, pump_type=pump_type,
                        material=material)
        
        self.life_time = life_time  # The life time of the pump [yrs]
        self.pump_cost = pump_cost  # The price of each pump [$/pump]
        self.working_factor = working_factor  # The working factor of the pump considering the head.
        self.operation_time = operation_time  # The operation time of the pump [h/d] 
        self.ppl = ppl  # The number of all people used all toilets
        self.baseline_ppl = baseline_ppl  # The number of people per toilet
        self.price_ratio = price_ratio

    @property
    def dP_factor(self):
        dP_factor = (self.P + self.dP_design) / self.P
        return dP_factor

    @property
    def actual_power(self):
        # Calculate the actual total power [kW]
        return self._pump_power / self.working_factor * self.dP_factor
    
    def _init_lca(self):
        self.construction = [Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Cast_iron'] = constr[0].quantity = self._CastIron_weight_per_pump
        self.add_construction(add_cost=False)

    def _run(self):
        
        # Input stream
        CWTwater = self.ins[0]
        air= self.ins[1]
        
        # Output stream
        Water_with_O2 = self.outs[0]
        
        # Inherite input stream
        Water_with_O2.mix_from([CWTwater, air])
    
    def _cost(self):

        C = self.baseline_purchase_costs
        C['Air Dissolved Pump'] = self.pump_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.actual_power * self.operation_time
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        # Replacement cost: based on pump count
        replacement_cost = self.pump_cost / (self.life_time * 365 * 24) * scale
        return replacement_cost

# %%

#@price_ratio
class MicroBubblePump(Pump):
    """
    Microbubble pump for injecting ozone (O₃) into clear water for disinfection.
    """

    _N_ins = 1;  # Number of input streams
    _N_outs = 1;  # Number of output streams
    _ins_size_is_fixed = True;  # Let the input interface be fixed
    _outs_size_is_fixed = True;  # Let the output interface be fixed
    _CastIron_weight_per_pump = 9 # The weight of cast iron [kg]
    _pump_power = 11.45  # The power of each pump [kW] #TODO: check pump power assumption
    _pump_flow_rate = 10  # Assumption: at 298K, 1.25atm=126656, 10g/L, rho=2.45kg/m3
    exponent_scale = 0.4
   
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream', material='Cast iron',
                 F_BM_default=3.3, pump_type='Default', isdynamic=False, ignore_NPSH=True,
                 dP_design=25331,  # The extra pressure of the pump [Pa]
                 P=101325,  # The pressure of the pump [Pa]
                 price_ratio=0.9,
                 life_time=None, pump_cost=None, working_factor=None,
                 operation_time=None, ppl=100, baseline_ppl=100):
                 
        # Initialize the LiftPump with parent class logic
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with, F_BM_default=F_BM_default,
                        isdynamic=isdynamic, ignore_NPSH=ignore_NPSH, P=P, dP_design=dP_design, pump_type=pump_type,
                        material=material)
        
        self.life_time = life_time  # The life time of the pump [yrs]
        self.pump_cost = pump_cost  # The price of each pump [$/pump]
        self.working_factor = working_factor  # The working factor of the pump considering the head.
        self.operation_time = operation_time  # The operation time of the pump [h/d] 
        self.ppl = ppl  # The number of all people used all toilets
        self.baseline_ppl = baseline_ppl  # The number of people per toilet
        self.price_ratio = price_ratio

    @property
    def dP_factor(self):
        dP_factor = (self.P + self.dP_design) / self.P
        return dP_factor

    @property
    def actual_power(self):
        # Calculate the actual total power [kW]
        return self._pump_power / self.working_factor * self.dP_factor
    
    def _init_lca(self):
        self.construction = [Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Cast_iron'] = constr[0].quantity = self._CastIron_weight_per_pump
        self.add_construction(add_cost=False)

    def _run(self):
        
        # Input stream
        gas_in = self.ins[0]
        
        # Output stream
        gas_out = self.outs[0]
        
        # Inherite input stream
        gas_out.copy_like(gas_in)
    
    def _cost(self):

        C = self.baseline_purchase_costs
        C['Micro Bubble Pump'] = self.pump_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.actual_power * self.operation_time
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        # Replacement cost: based on pump count
        replacement_cost = self.pump_cost / (self.life_time * 365 * 24) * scale
        return replacement_cost

# %%

#@price_ratio
class ClearWaterPump(Pump):
    """
    Clear Water Pump for reusing treated water in toilet flushing.
    """

    _N_ins = 1;  # Number of input streams
    _N_outs = 1;  # Number of output streams
    _ins_size_is_fixed = True;  # Let the input interface be fixed
    _outs_size_is_fixed = True;  # Let the output interface be fixed
    _CastIron_weight_per_pump = 9 # The weight of cast iron [kg]
    _pump_power = 0.75  # The power of each pump [kW]
    _pump_flow_rate = 30 * 60 / 1000  # The flow rate of each pump [m3/h]
    exponent_scale = 0.4
   
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', material='Cast iron',
                 F_BM_default=3.3, pump_type='Default', isdynamic=False, ignore_NPSH=True,
                 dP_design=202650,  # The extra pressure of the pump [Pa]
                 P=101325,  # The pressure of the pump [Pa]
                 price_ratio=0.9,
                 life_time=None, pump_cost=None, working_factor=None,
                 operation_time=None, ppl=100, baseline_ppl=100):
                 
        # Initialize the LiftPump with parent class logic
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with, F_BM_default=F_BM_default,
                        isdynamic=isdynamic, ignore_NPSH=ignore_NPSH, P=P, dP_design=dP_design, pump_type=pump_type,
                        material=material)
        
        self.life_time = life_time  # The life time of the pump [yrs]
        self.pump_cost = pump_cost  # The price of each pump [$/pump]
        self.working_factor = working_factor  # The working factor of the pump considering the head.
        self.operation_time = operation_time  # The operation time of the pump [h/d] 
        self.ppl = ppl  # The number of all people used all toilets
        self.baseline_ppl = baseline_ppl  # The number of people per toilet
        self.price_ratio = price_ratio

    @property
    def dP_factor(self):
        dP_factor = (self.P + self.dP_design) / self.P
        return dP_factor

    @property
    def actual_power(self):
        # Calculate the actual total power [kW]
        return self._pump_power / self.working_factor * self.dP_factor
    
    def _init_lca(self):
        self.construction = [Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Cast_iron'] = constr[0].quantity = self._CastIron_weight_per_pump
        self.add_construction(add_cost=False)

    def _run(self):
        
        # Input stream
        mass_in = self.ins[0]
        
        # Output stream
        mass_out = self.outs[0]
         
        # Inherite input stream
        mass_out.copy_like(mass_in)
    
    def _cost(self):

        C = self.baseline_purchase_costs
        C['Clear Water Pump'] = self.pump_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()

        power_demand = self.actual_power * self.operation_time
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        # Replacement cost: based on pump count
        replacement_cost = self.pump_cost / (self.life_time * 365 * 24) * scale
        return replacement_cost