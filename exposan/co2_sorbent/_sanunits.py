#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# import biosteam as bst
import biosteam as bst
from qsdsan import SanUnit
from math import sqrt, pi
from biosteam import Splitter
from biosteam.units.design_tools import PressureVessel

__all__ = (
    'ALFProduction',
    'ALFCrystallizer',
    'ALFTemperatureSwingAdsorption',
    'LowTemperatureElectrolysis'
    )

# =============================================================================
# ALFProduction
# =============================================================================

class ALFProduction(bst.CSTR):
    '''
    Reactor for ALF production.
    '''
    _N_ins = 1
    _N_outs = 1
    # TODO: confirm T and tau
    T_default = 60 + 273.15
    P_default = 101325
    tau_default = 8 # hr
    
    def _setup(self):
        super()._setup()
        self.ALF_production = bst.Reaction('AlH3O3,s + 3HCOOH,l -> C3H3AlO6,s + 3H2O,l', 'AlH3O3', 1)
    
    def _run(self):
        effluent = self.outs[0]
        effluent.copy_like(self.ins[0])
        self.ALF_production(effluent)
        effluent.T = self.T
        effluent.P = self.P

# =============================================================================
# ALFCrystallizer
# =============================================================================
class ALFCrystallizer(bst.BatchCrystallizer):
    '''
    Crystallier for ALF.
    
    Parameters
    ----------

    crystal_ALF_yield : float, optional
        ALF crystallization yield. Defacults to 1.
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 T, crystal_ALF_yield=1):
        bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                                       tau=5, V=1e6, T=T)
        self.crystal_ALF_yield = crystal_ALF_yield

    @property
    def Hnet(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        if 's' in feed.phases:
            H_in = - sum([i.Hfus * j for i,j in zip(self.chemicals, feed['s'].mol) if i.Hfus])
        else:
            H_in = 0.
        solids = effluent['s']
        H_out = - sum([i.Hfus * j for i,j in zip(self.chemicals, solids.mol) if i.Hfus])
        return H_out - H_in
        
    def _run(self):
        outlet = self.outs[0]
        outlet.phases = ('s', 'l')
        crystal_ALF_yield = self.crystal_ALF_yield
        feed = self.ins[0]
        ALF = feed.imass['C3H3AlO6']
        outlet.empty()
        
        outlet.imass['s', 'C3H3AlO6'] = ALF*crystal_ALF_yield
        outlet.imass['l', ('C3H3AlO6','HCOOH','H2O')] = [ALF*(1-crystal_ALF_yield), feed.imass['HCOOH'], feed.imass['H2O']]
        
        outlet.T = self.T

# =============================================================================
# ALFTemperatureSwingAdsorption
# =============================================================================
class ALFTemperatureSwingAdsorption(PressureVessel, Splitter):
    '''
    TSA using ALF as adsorbent for CO2 adsorption.
    
    Parameters
    ----------

    split : dict[str, float] or list[float], optional
        Component splits towards the effluent (0th outlet).
    superficial_velocity : float, optional
        Superficial velocity of the feed. The diameter of the receiving vessel adjusts
        accordingly. Defaults to 1080 m/h. Typical velocities are 540 to 2160 m/h for liquids [1]_.
    regeneration_velocity : float, optional
        Mean velocity of the fluid used to regenerate the bed. Defaults to 1080 m/h. 
        Common velocity range for gasses is 540 to 2160 m/h [1]_.
    cycle_time : float, optional
        Time at which the receiving vessel is switched. Defaults to 8 h [2]_.
    rho_adsorbent : float, optional
        The density of ALF. Defaults to 1441 [kg/m3] Table S1 [Y]_.
    adsorbent_capacity : float, optional
        Amount of CO2 that ALF can hold. Defaults to 2.7 mmol/g Table S7 [3]_.
    T_regeneration : float, optional
        Temperature during the regeneration phase. Defaults to 418 K Table S8 [3]_.
    vessel_material : float, optional
        Vessel material. Defaults to 'Stainless steel 316',
    vessel_type : float, optional
        Vessel type. Defaults to 'Vertical'.
    length_unused : float, optional
        Additional length of a column to account for mass transfer limitations (due to unused bed). Defaults to 2 ft per column.
    waste_ratio : float, optiona;
        Wasted ALF ratio per run. Defaults to 0.07.
    
    References
    ----------
    [1] Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
    [2] https://www.chemicalprocessing.com/processing-equipment/fluid-handling/article/11302111/select-the-right-valves-for-adsorption-processes-chemical-processing
    [3] Evans et al. Science Advances SI. # TODO: update the reference
    '''
    auxiliary_unit_names = ('heat_exchanger_regeneration','heat_exchanger_cooling')
    
    # in $/ft3 # TODO: this cost is just for the initial ALF in 3 columns, replace the value
    adsorbent_cost = 1
    
    # in year # TODO: confirm the value
    _default_equipment_lifetime = 10
    
    _N_ins = 3
    _N_outs = 4
    
    auxiliary_unit_names=('heat_exchanger_regeneration','heat_exchanger_cooling')
    
    def _init(self,
              split=dict(O2=0, N2=0, CO2=1),
              superficial_velocity=1080, # m/h
              regeneration_velocity=1080, # m/h # TODO: need to decide if this is needed since there is no carrier gas stream
              cycle_time=8, # h
              rho_adsorbent = 1441, # kg/m3
              adsorbent_capacity=2.7, # mmol/g # TODO: need to convert the unit
              T_regeneration=418, # K
              vessel_material='Stainless steel 316',
              vessel_type='Vertical',
              length_unused=1.219, # m
              waste_ratio=0.07, # ratio of ALF that is wasted (and is sent for CO2 storage) each cycle # TODO: confirm this is for each cycle
              ):
        bst.Splitter._init(self, split=split)
        self.superficial_velocity = superficial_velocity
        self.regeneration_velocity = regeneration_velocity
        self.cycle_time = cycle_time
        self.rho_adsorbent = rho_adsorbent
        self.adsorbent_capacity = adsorbent_capacity
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.T_regeneration = T_regeneration
        self.length_unused = length_unused
        self.waste_ratio = waste_ratio
        self.heat_exchanger_regeneration = bst.HXutility(None, None, None, thermo=self.thermo)
        self.heat_exchanger_cooling = bst.HXutility(None, None, None, thermo=self.thermo)
    
    def _run(self):
        feed, ALF, regen_in = self.ins
        offgas, CO2, used_ALF, regen_out = self.outs
        
        # TODO: flue gas temperature (may set high, and use a HX to cool first, then add HXN to offset)
        # TODO: ALF temperature, should be room temperature (same as the flue gas above)

        for i in self.outs: i.empty()
        
        # TODO: make sure the order of CO2 and offgas is correct
        feed.split_to(CO2, offgas, self.split)
        F_vol_feed = feed.F_vol
        superficial_velocity = self.superficial_velocity
        F_mass_adsorbate = CO2.imass['CO2']
        
        self.diameter = diameter = 2 * sqrt(F_vol_feed / (superficial_velocity * pi))
        self.area = area = pi * diameter * diameter / 4
        total_length = (
            self.cycle_time * F_mass_adsorbate / (self.adsorbent_capacity * self.rho_adsorbent * area)
        ) + self.length_unused # length of equilibrium section plus unused bed (LES + LUB)
        self.length = length = total_length / 2 # Size of each column
        self.vessel_volume = length * area
        
        ALF.phase = 'l'
        
        regen_in.copy_like(ALF)
        regen_in.imass['C3H3AlO6'] /= self.waste_ratio
        regen_in.T = ALF.T
        
        regen_out.copy_like(ALF)
        regen_out.imass['C3H3AlO6'] /= self.waste_ratio
        regen_out.T = self.T_regeneration
        
        used_ALF.copy_like(ALF)
    
    def _design(self):
        feed, ALF, regen_in = self.ins
        offgas, CO2, used_ALF, regen_out = self.outs 
        design_results = self.design_results
        diameter = self.diameter
        length = self.length
        design_results['Number of reactors'] = 3
        design_results.update(
            self._vessel_design(
                feed.P * 0.000145038, # Pa to psi
                diameter * 3.28084, # m to ft
                length * 3.28084, # m to ft
            )
        )
        hxr = self.heat_exchanger_regeneration
        hxr.ins.empty()
        hxr.outs.empty()
        hxr.ins[0] = regen_in.copy()
        hxr.outs[0] = regen_out.copy()
        hxr.T = regen_out.T
        hxr.simulate_as_auxiliary_exchanger(ins=hxr.ins, outs=hxr.outs)
        
        hxc = self.heat_exchanger_cooling
        hxc.ins.empty()
        hxc.outs.empty()
        hxc.ins[0] = regen_out.copy()
        hxc.outs[0] = regen_in.copy()
        hxc.T = regen_in.T
        hxc.simulate_as_auxiliary_exchanger(ins=hxc.ins, outs=hxc.outs)
    
    def _cost(self):
        design_results = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        baseline_purchase_costs.update(self._vessel_purchase_cost(
            design_results['Weight'], design_results['Diameter'], design_results['Length']))
        N_reactors = design_results['Number of reactors']
        for i, j in baseline_purchase_costs.items():
            baseline_purchase_costs[i] *= N_reactors
        baseline_purchase_costs['initial_ALF'] = N_reactors * 35.3147 * self.vessel_volume * self.adsorbent_cost





# =============================================================================
# LowTemperatureElectrolysis
# =============================================================================

class LowTemperatureElectrolysis(SanUnit):
    '''
    Electrochemical cell for nutrient recovery.

    This unit has the following equipment:
        - :class:`~.equipments.Column`
        - :class:`~.equipments.Machine`
        - :class:`~.equipments.Electrode`
        - :class:`~.equipments.Membrane`

    Parameters
    ----------
    recovery : dict
        Keys refer to chemical component IDs. Values refer to recovery fractions (with 1 being 100%) for the respective chemicals.
    removal : dict
        Keys refer to chemical component IDs. Values refer to removal fractions (with 1 being 100%) for the respective chemicals.
    equipment : list(obj)
        List of Equipment objects part of the Electrochemical Cell.
    OPEX_over_CAPEX : float
        Ratio with which operating costs are calculated as a fraction of capital costs

    Example
    -------
    >>> # Set components
    >>> import qsdsan as qs
    >>> kwargs = dict(particle_size='Soluble',
    ...               degradability='Undegradable',
    ...               organic=False)
    >>> H2O = qs.Component.from_chemical('H2O', phase='l', **kwargs)
    >>> NH3 = qs.Component.from_chemical('NH3', phase='g', **kwargs)
    >>> NH3.particle_size = 'Dissolved gas'
    >>> H2SO4 = qs.Component.from_chemical('H2SO4', phase='l', **kwargs)
    >>> AmmoniumSulfate = qs.Component.from_chemical('AmmoniumSulfate', phase='l',
    ...                                              **kwargs)
    >>> Na_ion = qs.Component.from_chemical('Na+', phase='s', **kwargs)
    >>> K_ion = qs.Component.from_chemical('K+', phase='s', **kwargs)
    >>> Cl_ion = qs.Component.from_chemical('Cl-', phase='s', **kwargs)
    >>> PO4_ion = qs.Component.from_chemical('Phosphate', phase='s', **kwargs)
    >>> SO4_ion = qs.Component.from_chemical('Sulfate', phase='s', **kwargs)
    >>> C = qs.Component.from_chemical('Carbon', phase='s', **kwargs)
    >>> COD = qs.Component.from_chemical('O2', phase='s', **kwargs)
    >>> cmps = qs.Components((H2O, NH3, H2SO4, AmmoniumSulfate, Na_ion, K_ion, Cl_ion, PO4_ion, SO4_ion, C, COD))
    >>> # Assuming all has the same molar volume as water for demonstration purpose
    >>> for cmp in cmps:
    ...     cmp.copy_models_from(H2O, names=['V'])
    ...     defaulted_properties = cmp.default()
    >>> qs.set_thermo(cmps)
    >>> # Set waste streams
    >>> influent = qs.WasteStream('influent')
    >>> influent.set_flow_by_concentration(flow_tot=0.5, concentrations={'NH3':3820,
    ...                                     'Na+':1620, 'K+':1470, 'Cl-':3060, 'Phosphate':169,
    ...                                     'Sulfate':1680, 'Carbon':1860, 'O2':3460}, units=('mL/min', 'mg/L'))
    >>> influent.show()
    WasteStream: influent
     phase: 'l', T: 298.15 K, P: 101325 Pa
     flow (g/hr): H2O        29.5
                  NH3        0.115
                  Na+        0.0486
                  K+         0.0441
                  Cl-        0.0918
                  Phosphate  0.00507
                  Sulfate    0.0504
                  Carbon     0.0558
                  O2         0.104
     WasteStream-specific properties:
      pH         : 7.0
      Alkalinity : 2.5 mg/L
      TC         : 1860.0 mg/L
      TP         : 31.9 mg/L
      TK         : 1470.0 mg/L
     Component concentrations (mg/L):
      H2O              984514.0
      NH3              3820.0
      Na+              1620.0
      K+               1470.0
      Cl-              3060.0
      Phosphate        169.0
      Sulfate          1680.0
      Carbon           1860.0
      O2               3460.0
    >>> catalysts = qs.WasteStream('catalysts', H2SO4=0.00054697, units=('kg/hr'))
    >>> # kg/hr imass value derived from 0.145 moles used on average for one, 26 hour experiment in the model cell
    >>> catalysts.price = 8.4 #in $/kg
    >>> # electricity price for ECS experiments in the default model tested by the Tarpeh Lab
    >>> # if want to change the price of electricity,
    >>> # use qs.PowerUtility, e.g., qs.PowerUtility.price = 0.1741
    
    >>> # Set the unit
    >>> U1 = qs.sanunits.ElectrochemicalCell('unit_1', ins=(influent, catalysts), outs=('recovered', 'removed', 'residual'))
    >>> # Simulate and look at the results
    >>> U1.simulate()
    >>> U1.results()
    Electrochemical cell                                      Units                              unit_1
    Electricity         Power                                    kW                            5.42e-05
                        Cost                                 USD/hr                            4.24e-06
    Design              Electrode unit_1_Main_Anode - Nu...    None                                   1
                        Electrode unit_1_Main_Anode - Ma...    None  titanium grid catalyst welded t...
                        Electrode unit_1_Main_Anode - Su...      m2                                   1
                        Electrode unit_1_Main_Cathode - ...    None                                   1
                        Electrode unit_1_Main_Cathode - ...    None  timesetl 3pcs stainless steel w...
                        Electrode unit_1_Main_Cathode - ...      m2                                30.2
                        Electrode unit_1_Current_Collect...    None                                   1
                        Electrode unit_1_Current_Collect...    None    stainless steel 26 gauge 5.5 x 7
                        Electrode unit_1_Current_Collect...      m2                                38.5
                        Electrode unit_1_Reference_Elect...    None                                   1
                        Electrode unit_1_Reference_Elect...    None  re-5b ag/agcl, 7.5 cm long, wit...
                        Electrode unit_1_Reference_Elect...      m2                                   1
                        Membrane unit_1_Cation_Exchange_...                                           1
                        Membrane unit_1_Cation_Exchange_...          CMI-7000S, polystyrene 0.45mm t...
                        Membrane unit_1_Cation_Exchange_...      m2                                30.2
                        Membrane unit_1_Gas_Permeable_Me...                                           1
                        Membrane unit_1_Gas_Permeable_Me...          Aquastill 0.3-micron polyethyle...
                        Membrane unit_1_Gas_Permeable_Me...      m2                                30.2
    Purchase cost       unit_1_Main_Anode                       USD                                 288
                        unit_1_Main_Cathode                     USD                               0.847
                        unit_1_Current_Collector_Cathode        USD                                39.9
                        unit_1_Reference_Electrode              USD                                  94
                        unit_1_Cation_Exchange_Membrane         USD                                 167
                        unit_1_Gas_Permeable_Membrane           USD                                29.3
                        Exterior                                USD                                60.5
    Total purchase cost                                         USD                                 679
    Utility cost                                             USD/hr                            4.24e-06
    Additional OPEX                                          USD/hr                                 136
    >>> U1.show() # doctest: +ELLIPSIS
    ElectrochemicalCell: unit_1
    ins...
    [0] influent
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): H2O        29.5
    ...
    '''

    _N_ins = 2
    _N_outs = 3


    def __init__(self, ID='', ins=(), outs=(),
                 recovery={'NH3':0.7}, removal={'NH3':0.83, 'K+':0.83, "Na+":0.8}, OPEX_over_CAPEX=0.2):
        SanUnit.__init__(self=self, ID=ID, ins=ins, outs=outs)
        self.recovery = recovery
        self.removal = removal
        self.OPEX_over_CAPEX = OPEX_over_CAPEX


        self.equipment = [
            Electrode('Main_Anode', linked_unit=self, N=1, electrode_type='anode',
                      material='Titanium grid catalyst welded to current collector tab both coated in iridium tantalum mixed metal oxide', surface_area=1, unit_cost=288), #288/unit, 1 unit
            Electrode('Main_Cathode', linked_unit=self, N=1, electrode_type='cathode',
                      material='TIMESETL 3pcs Stainless Steel Woven Wire 20 Mesh - 12"x8"(30x21cm) Metal Mesh Sheet 1mm Hole Great for Air Ventilation - A4', surface_area=30.25, unit_cost=0.847), #in in^2
            Electrode('Current_Collector_Cathode', linked_unit=self, N=1, electrode_type='cathode',
                      material='Stainless Steel 26 gauge 5.5'' x 7''', surface_area=38.5, unit_cost=39.9245), #in unknown units (94/unit, 1 unit)
            Electrode('Reference_Electrode', linked_unit=self, N=1, electrode_type='reference',
                      material='RE-5B Ag/AgCl, 7.5 cm long, with ceramic (MF-2056)', surface_area=1, unit_cost=94), #in unknown units (94/unit, 1 unit)
            Membrane('Cation_Exchange_Membrane', linked_unit=self, N=1,
                     material='CMI-7000S, polystyrene 0.45mm thick [48'' x 20'']',
                     unit_cost=5.5055, surface_area=30.25), # in in^2
            Membrane('Gas_Permeable_Membrane', linked_unit=self, N=1,
                     material='Aquastill 0.3-micron polyethylene membrane', unit_cost=0.968, surface_area=30.25), #in in^2
            ]


    def _run(self):
        influent, catalysts = self.ins
        recovered, removed, residual = self.outs[0], self.outs[1], self.outs[2]

        mixture = WasteStream()
        mixture.mix_from(self.ins)
        residual.copy_like(mixture)

        for chemical, recovery_ratio in self.recovery.items():
            recovered.imass[chemical] = mixture.imass[chemical]*recovery_ratio
            
        for chemical, removal_ratio in self.removal.items():
            recovery_ratio = 0 if recovered.imass[chemical] is None else recovered.imass[chemical]
            removed.imass[chemical] = mixture.imass[chemical]*removal_ratio - mixture.imass[chemical]*recovery_ratio
            residual.imass[chemical] = residual.imass[chemical]-residual.imass[chemical]*removal_ratio

    def _design(self):
        self.add_equipment_design()

    def _cost(self):
        self.add_equipment_cost()
        self.baseline_purchase_costs['Exterior'] = 60.52
        '''
        TOTAL CELL_EXTERIOR_COST = 60.52 USD
        Breakdown:
        Exterior frame (7'' x 7'' x 11/16'')	$21.46
        Interior half-cells (5.5'' x 5.5'' x 11/16'')	$19.87
        Rubber Sheets	$9.196
        Threaded Rods	$2.9696
        Wingnuts	$2.5504
        Flat Washers	$0.728
        Nylon Cable Glands	$3.744
        '''
        self.equip_costs = self.baseline_purchase_costs.values()
        add_OPEX = sum(self.equip_costs)*self.OPEX_over_CAPEX
        recovered, removed = self.outs[0], self.outs[1]

        self.power_utility.rate = recovered.imass['NH3']*0.67577
        # steady state value derived from 17.57 kWh used over 26 hrs
        self._add_OPEX = {'Additional OPEX': add_OPEX}