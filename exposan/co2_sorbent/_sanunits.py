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
import pandas as pd, biosteam as bst, qsdsan as qs
from qsdsan import SanUnit
from math import sqrt, pi
from biosteam import Splitter
from biosteam.units.design_tools import PressureVessel
from biosteam.units.decorators import cost

__all__ = (
    'ALFProduction',
    'ALFCrystallizer',
    'ALFTemperatureSwingAdsorption',
    'CO2ElectrolyzerSystem'
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
    [3] Evans, H. A.; Mullangi, D.; Deng, Z.; Wang, Y.; Peh, S. B.; Wei, F.;
        Wang, J.; Brown, C. M.; Zhao, D.; Canepa, P.; Cheetham, A. K.
        Aluminum Formate, Al(HCOO)3: An Earth-Abundant, Scalable, and Highly
        Selective Material for CO2 Capture. Science Advances 2022, 8 (44),
        eade1473. https://doi.org/10.1126/sciadv.ade1473.

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
        offgas, carbon_dioxide, used_ALF, regen_out = self.outs
        
        # TODO: flue gas temperature (may set high, and use a HX to cool first, then add HXN to offset)
        # TODO: ALF temperature, should be room temperature (same as the flue gas above)

        for i in self.outs: i.empty()
        
        # TODO: make sure the order of CO2 and offgas is correct
        feed.split_to(carbon_dioxide, offgas, self.split)
        F_vol_feed = feed.F_vol
        superficial_velocity = self.superficial_velocity
        F_mass_adsorbate = carbon_dioxide.imass['CO2']
        
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
        offgas, carbon_dioxide, used_ALF, regen_out = self.outs 
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
# CO2ElectrolyzerSystem
# =============================================================================

# TODO: do we need to update the cost data, e.g., H2A (see the referred paper) $/m2, e.g., 2018$ (the referred paper is published in 2018) to 2020$ (if correct and necessary)

# TODO: after determine CEPCI, update here, or find a way to match the CEPCI here with the CEPCI in the system (set CE here as bst.CE does not work)
@cost(basis='Electrolyzer area', ID='Electrolyzer', units='m^2',
      cost=250.25*1.75*175/1000*10*1.2*(1+1/0.65*0.35), S=1, CE=qs.CEPCI_by_year[2020], n=1)
@cost(basis='Electrolyte flow rate (ethanol)', ID='Distiller (ethanol)', units='L/min',
      cost=4162240, S=1000, CE=qs.CEPCI_by_year[2020], n=0.7)
@cost(basis='Electrolyte flow rate (formic acid)', ID='Distiller (formic acid)', units='L/min',
      cost=6896190, S=1000, CE=qs.CEPCI_by_year[2020], n=0.7)
@cost(basis='Electrolyte flow rate (methanol)', ID='Distiller (methanol)', units='L/min',
      cost=4514670, S=1000, CE=qs.CEPCI_by_year[2020], n=0.7)
@cost(basis='Electrolyte flow rate (propanol)', ID='Distiller (propanol)', units='L/min',
      cost=4687910, S=1000, CE=qs.CEPCI_by_year[2020], n=0.7)
@cost(basis='Total gas flow for PSA', ID='PSA', units='m^3/h',
      cost=1989043, S=1000, CE=qs.CEPCI_by_year[2020], n=0.7)
class CO2ElectrolyzerSystem(SanUnit):
    '''
    CO2 electrolyzer system that converts CO2 into reduced 1C, 2C, and nC products [1]_. 
    
    Parameters
    ----------
    target_product : str, optional
        target product of the CO2 reduction system, can only be 'carbon monoxide',
        'ethanol','ethylene','formic acid','methane','methanol', and 'propanol'.
        Defaults to 'formic acid'.
    current_density : float, optional
        Defaults to 0.2 A/cm2.
    cathodic_overpotential : float, optional
        Defaults to 0.454 V.
    cell_voltage : float, optional
        Defaults to 2.3 V.
    product_selectivity : float, optional
        Defaults to 0.9.
    converstion : float, optional
        Defaults to 0.5.
    PSA_operating_cost : float, optional
        Defaults to 0.25 kWh/m3.
    operating_days_per_year: float/int
        Same to the operating day per year of the system.
    
    References
    ----------
    [1] Jouny, M.; Luc, W.; Jiao, F. General Techno-Economic Analysis of CO2
        Electrolysis Systems. Ind. Eng. Chem. Res. 2018, 57 (6), 2165–2177.
        https://doi.org/10.1021/acs.iecr.7b03514.
    '''
    
    _N_ins = 2
    _N_outs = 2
    
    _units= {'Electrolyzer area': 'm^2',
             'Electrolyte flow rate (ethanol)': 'L/min',
             'Electrolyte flow rate (formic acid)': 'L/min',
             'Electrolyte flow rate (methanol)': 'L/min',
             'Electrolyte flow rate (propanol)': 'L/min',
             'Total gas flow for PSA': 'm^3/h'}
    
    def __init__(self, ID='', ins=(), outs=(), target_product='formic acid', current_density=0.2,
                 cell_voltage=2.3, cathodic_overpotential=0.454, product_selectivity=0.9,
                 converstion=0.5, PSA_operating_cost=0.25, operating_days_per_year=350):
        SanUnit.__init__(self=self, ID=ID, ins=ins, outs=outs)
        self.target_product = target_product
        self.current_density = current_density
        self.cathodic_overpotential = cathodic_overpotential
        self.cell_voltage = cell_voltage
        self.product_selectivity = product_selectivity
        self.converstion = converstion
        self.PSA_operating_cost = PSA_operating_cost
        self.operating_days_per_year = operating_days_per_year
    
    def _run(self):
        carbon_dioxide, water = self.ins
        product, mixed_offgas = self.outs
        
        if self.target_product not in ['carbon monoxide','ethanol','ethylene','formic acid','methane','methanol','propanol']:
            raise ValueError("target product must be in 'carbon monoxide', 'ethanol', 'ethylene', 'formic acid', 'methane', 'methanol', and 'propanol'")
        
        product_info = {'chemical': ['carbon monoxide','ethanol','ethylene','formic acid','methane','methanol','propanol'], # the propanol is n-propanol
                        'formula': ['CO','C2H6O','C2H4','HCOOH','CH4','CH4O','C3H8O'],
                        'market_price': [0.6, 1.003, 1.3, 0.735, 0.18, 0.577, 1.435], # $/kg
                        'elec_number': [2, 12, 12, 2, 8, 6, 18],
                        'elec_number_per_CO2': [2, 6, 6, 2, 8, 6, 6],
                        'mole_ratio': [1, 2, 2, 1, 1, 1, 3],
                        'MW': [28.01, 46.06, 28.05, 46.025, 16.04, 32.04, 60.06],
                        'density': [1.14, 789, 1.18, 1221, 0.656, 792, 803], # kg/m3
                        'state': ['gas','liq','gas','liq','gas','liq','liq'],
                        'potential': [-0.106, 0.084, 0.064, -0.25, 0.169, 0.016, 0.095],
                        'voltage_theory': [1.336, 1.146, 1.166, 1.48, 1.061, 1.214, 1.135],
                        'voltage_practical': [1.736, 1.546, 1.566, 1.88, 1.461, 1.614, 1.535],
                        'low_price': [i*0.85 for i in [0.6, 1.003, 1.3, 0.735, 0.18, 0.577, 1.435]],
                        'high_price': [i*1.15 for i in [0.6, 1.003, 1.3, 0.735, 0.18, 0.577, 1.435]]}
        
        product_info = pd.DataFrame.from_dict(product_info)
        
        chemical_info = self.chemical_info = product_info[product_info['chemical'] == self.target_product]
        
        # CO2 inlet flow rate [kg/h]
        CO2_inlet_flow_rate = carbon_dioxide.imass['CO2']
        
        # converted CO2 amount [kg/day]
        CO2_converted = CO2_inlet_flow_rate*24*self.converstion
        
        # CO2 outlet flow rate [kg/h]
        CO2_outlet_flow_rate_kg_per_h = CO2_inlet_flow_rate - CO2_converted/24
        
        # CO2 outlet flow rate [m3/h]
        self.CO2_outlet_flow_rate_m3_per_h = CO2_outlet_flow_rate_kg_per_h/1.98 # the density of CO2 is 1.98 kg/m3
        
        # production amount based on inlet CO2 [kg/day]
        product_production = CO2_converted/44/float(chemical_info['mole_ratio'])*float(chemical_info['MW'])
        
        # required current [A]
        current_needed = product_production/24/3600*1000/float(chemical_info['MW'])*float(chemical_info['elec_number'])*96485/self.product_selectivity
        
        # required electrolzer area [m2]
        self.electrolyzer_area = current_needed/self.current_density/10000
        
        # required power [MW]
        self.power_needed = current_needed*self.cell_voltage/1000000
        
        # gas product flow rate [m3/h]
        gas_product_flow_rate = 0 if chemical_info['state'].to_string(index=False) == 'liq' else product_production/float(chemical_info['density'])/24
        
        # liquid product flow rate [m3/h]
        liquid_product_flow_rate_m3_per_h = 0 if chemical_info['state'].to_string(index=False) == 'gas' else product_production/float(chemical_info['density'])/24
        
        # liquid product flow rate [l/min]
        liquid_product_flow_rate_l_per_min = liquid_product_flow_rate_m3_per_h*1000/60
        
        # electrolyte flow rate [l/min]
        # it is assumed that product-rich electrolyte is recycled until a steady-state volume concentration of 10% if reached
        self.electrolyte_flow_rate = liquid_product_flow_rate_l_per_min/0.1
        
        # hydrogen flow rate [mol/s]
        # 2 represents 2 e-
        hydrogen_flow_rate_mol_per_s = current_needed*(1-self.product_selectivity)/2/96485
        
        # hydrogen flow rate [m3/h]
        # hydrogen molar mass: 2.016, hydrogen density: 0.08375 kg/m3
        hydrogen_flow_rate_m3_per_h = hydrogen_flow_rate_mol_per_s*2.016/1000/0.08375*3600
        
        # total gas flow [m3/h]
        self.total_gas_flow = self.CO2_outlet_flow_rate_m3_per_h + gas_product_flow_rate + hydrogen_flow_rate_m3_per_h
        
        product.imass[chemical_info['formula'].to_string(index=False)] = product_production/24
        product.phase = 'g' if chemical_info['state'].to_string(index=False) == 'gas' else 'l'
        
        water.imass['H2O'] = current_needed/4/96485*18/1000*3600
        # TODO: check mass balance and energy balance
        mixed_offgas.imass['CO2'] = CO2_outlet_flow_rate_kg_per_h
        mixed_offgas.imass['H2'] = hydrogen_flow_rate_m3_per_h*0.08375 # hydrogen density: 0.08375 kg/m3
        mixed_offgas.imass['O2'] = carbon_dioxide.F_mass + water.F_mass - product.F_mass - mixed_offgas.imass['CO2'] - mixed_offgas.imass['H2']
        mixed_offgas.phase = 'g'
    
    def _design(self):
        D = self.design_results
        D['Electrolyzer area'] = self.electrolyzer_area
        
        for liquid_product in ['ethanol','formic acid','methanol','propanol']:
            D[f'Electrolyte flow rate ({liquid_product})'] = self.electrolyte_flow_rate if self.chemical_info['chemical'].to_string(index=False) == liquid_product else 0
        
        D['Total gas flow for PSA'] = 0 if self.CO2_outlet_flow_rate_m3_per_h == self.total_gas_flow else self.total_gas_flow
        
        self.add_power_utility(self.power_needed*1000 + D['Total gas flow for PSA']*self.PSA_operating_cost)
        
        # TODO: add construction for LCA
        # if self.include_construction:
        #     construction = getattr(self, 'construction', []) # would work for both biosteam/qsdsan units
        #     if construction: construction[0].quantity = pipe + pumps
        #     else:
        #         self.construction = [
        #             Construction('stainless_steel', linked_unit=self, item='Stainless_steel', 
        #                          quantity=pipe + pumps, quantity_unit='kg'),
        #             ]
    
    def _cost(self):
        self._decorated_cost()
        
        distillation_OPEX = {'ethanol': 3463310,
                             'formic acid': 11213200,
                             'methanol': 4027820,
                             'propanol': 5610420}
        if self.chemical_info['chemical'].to_string(index=False) in list(distillation_OPEX.keys()):
            self.add_OPEX = {f'{self.chemical_info["chemical"].to_string(index=False)}_distillation_OPEX':
                             self.electrolyte_flow_rate/1000*distillation_OPEX[self.chemical_info['chemical'].to_string(index=False)]/self.operating_days_per_year/24}