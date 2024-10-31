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
from qsdsan import (
    SanUnit,
    sanunits as qsu,
    Stream,
    )

__all__ = (
    'ElectrochemicalOxidation',
    'BiocrudeDeashing',
    'BiocrudeDewatering',
    'Disposal',
    'PilotHTL',
    'Scaler',
    )

CEPCI_by_year = qs.utils.tea_indices['CEPCI']

# %%

class Scaler(SanUnit):
    '''
    Scale up the influent or the effluent by a specified number.
    
    Parameters
    ----------
    ins : seq(obj)
        Stream before scaling.
    outs : seq(obj)
        Stream after scaling.
    scaling_factor : float
        Factor for which the effluent will be scaled.
    reverse : bool
        If True, will scale the influent based on the effluent.
        E.g., for a scaling factor of 2, when `reverse` is False, 
        all components in the effluent will have a mass flowrate that is 2X of the influent;
        when `reverse` is True,
        all components in the influent will have a mass flowrate that is 2X of the effluent.
    '''
    
    _N_ins = _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', F_BM_default=1,
                  scaling_factor=1, reverse=False, **kwargs,
                  ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
        self.scaling_factor = scaling_factor
        self.reverse = reverse
        for kw, arg in kwargs.items(): setattr(self, kw, arg)
        
    def _run(self):
        inf = self.ins[0]
        eff = self.outs[0]
        factor = self.scaling_factor
        if self.reverse is False:
            eff.copy_like(inf)
            eff.F_mass *= factor
        else:
            inf.copy_like(eff)
            inf.F_mass *= factor
        
        

# %%

base_feedstock_flowrate = 11.46 # kg/hr
salad_dressing_waste_yields = (0.5219, 0.2925, 0.1756)

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
    Pilot-scale reactor for hydrothermal liquefaction (HTL) of wet organics.
    Biocrude from mulitple pilot-scale reactors will be transported to a central plant
    for biocrude upgrading.
    
    Parameters
    ----------
    ins : obj
        Waste stream for HTL.
    outs : seq(obj)
        Hydrochar, aqueous, biocrude, offgas.
    tau : float
        Retention time, [hr].
    V_wf : float
        Reactor working volumne factor, volume of waste streams over total volume.
    N_unit : int
        Number of required HTL unit.
    afdw_yields : seq(float)
        Yields for biocrude, aqueous, and gas products on ash-free dry weight basis of the feedstock.
        Yield of the hydrochar product will be calculated by subtraction to close the mass balance.
        All ash assumed to go to the aqueous product.
    piping_cost_ratio : float
        Piping cost estimated as a ratio of the total reactor cost.
    accessory_cost_ratio : float
        Accessories (e.g., valves) cost estimated as a ratio of the total reactor cost.
    '''
    
    _N_ins = 1
    _N_outs = 4
    
    _units= {
        'Feedstock dry flowrate': 'kg/hr',
        'Non-scaling factor': 'ea',
        }
    
    
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
                  tau=15/60, V_wf=0.45,
                  N_unit=1,
                  afdw_yields=salad_dressing_waste_yields,
                  piping_cost_ratio=0.15,
                  accessory_cost_ratio=0.08,
                  **kwargs,
                  ):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        #!!! Need to compare the externally sourced HX cost and BioSTEAM default
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = qsu.HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=self.eff_T, rigorous=True)
        self.tau = tau
        self.V_wf = V_wf
        self.N_unit = N_unit
        self._afdw_yields = afdw_yields
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
        
        for stream in self.outs:
            stream.T = self.heat_exchanger.T
        

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

        
    def _cost(self):
        self.parallel['self'] = self.N_unit
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
        

    @property
    def dry_mass_in(self):
        '''[float] Total dry mass of the feedstock, kg/hr.'''
        feedstock = self.ins[0]
        return feedstock.F_mass-feedstock.imass[self.water_ID]

    @property
    def afdw_mass_in(self):
        '''[float] Total ash-free dry mass of the feedstock, kg/hr.'''
        feedstock = self.ins[0]
        return feedstock.F_mass-feedstock.imass[self.ash_ID]-feedstock.imass[self.water_ID]

    @property
    def afdw_biocrude_yield(self):
        '''[float] Biocrude product yield on the ash-free dry weight basis of the feedstock.'''
        return self._afdw_yields[0]

    @property
    def afdw_aqueous_yield(self):
        '''[float] Aquoues product yield on the ash-free dry weight basis of the feedstock.'''
        return self._afdw_yields[1]
    
    @property
    def afdw_gas_yield(self):
        '''[float] Gas product yield on the ash-free dry weight basis of the feedstock.'''
        return self._afdw_yields[2]

    @property
    def afdw_hydrochar_yield(self):
        '''[float] Hydrochar product yield on the ash-free dry weight basis of the feedstock.'''
        char_yield = 1-self.afdw_biocrude_yield-self.afdw_aqueous_yield-self.afdw_gas_yield
        if char_yield < 0:
            raise ValueError('Sum of biocrude, aqueous, and gas product exceeds 100%.')
        return char_yield
    
    @property
    def N_unit(self):
        '''
        [int] Number of HTL units.
        '''
        return self._N_unit
    @N_unit.setter
    def N_unit(self, i):
        self.parallel['self'] = self._N_unit = math.ceil(i)
    
        
# %%
        
base_biocrude_flowrate = 5.64 # kg/hr
@cost(basis='Biocrude flowrate', ID= 'Deashing Tank', units='kg/hr',
      cost=4330, S=base_biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
class BiocrudeDeashing(SanUnit):
    '''
    Biocrude deashing unit.
    
    Parameters
    ----------
    ins : obj
        HTL biocrude.
    outs : seq(obj)
        Deashed biocrude, ash for disposal.
    '''
    
    _N_outs = 2
    _units= {'Biocrude flowrate': 'kg/hr',}
    target_ash = 0.01 # dry weight basis
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', F_BM_default=1,
                  N_unit=1, **kwargs,
                  ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
        self.N_unit = N_unit
        for kw, arg in kwargs.items(): setattr(self, kw, arg)
    
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
        self.design_results['Biocrude flowrate'] = self.ins[0].F_mass
        self.parallel['self'] = self.N_unit

    @property
    def N_unit(self):
        '''
        [int] Number of deashing units.
        '''
        return self._N_unit
    @N_unit.setter
    def N_unit(self, i):
        self.parallel['self'] = self._N_unit = math.ceil(i)
            

@cost(basis='Biocrude flowrate', ID= 'Dewatering Tank', units='kg/hr',
      cost=4330, S=base_biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
class BiocrudeDewatering(SanUnit):
    '''
    Biocrude dewatering unit.
    
    Parameters
    ----------
    ins : obj
        HTL biocrude.
    outs : seq(obj)
        Dewatered biocrude, water for treatment.
    '''
    
    _N_outs = 2
    _units= {'Biocrude flowrate': 'kg/hr',}
    target_moisture = 0.01 # weight basis
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', F_BM_default=1,
                  N_unit=1, **kwargs,
                  ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
        self.N_unit = N_unit
        for kw, arg in kwargs.items(): setattr(self, kw, arg)
    
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
        self.design_results['Biocrude flowrate'] = self.ins[0].F_mass
        self.parallel['self'] = self.N_unit
        
    @property
    def N_unit(self):
        '''
        [int] Number of dewatering units.
        '''
        return self._N_unit
    @N_unit.setter
    def N_unit(self, i):
        self.parallel['self'] = self._N_unit = math.ceil(i)


# %%

base_ap_flowrate = 49.65 #kg/hr
# @cost(basis='Aqueous flowrate', ID= 'Sand Filtration Unit', units='kg/hr',
#       cost=318, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.7)
# # @cost(basis='Aqueous flowrate', ID= 'EC Oxidation Tank', units='kg/hr',
# #       cost=1850, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
# # @cost(basis='Aqueous flowrate', ID= 'Biological Treatment Tank', units='kg/hr',
# #       cost=4330, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
# @cost(basis='Aqueous flowrate', ID= 'Liquid Fertilizer Storage', units='kg/hr',
#       cost=7549, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
# class AqueousFiltration(SanUnit):
#     '''
#     HTL aqueous filtration unit.
    
#     Parameters
#     ----------
#     ins : seq(obj)
#         Any number of influent streams to be treated.
#     outs : seq(obj)
#         Fertilizer, recycled process water, waste.
#     N_unit : int
#         Number of required filtration unit.
#     '''
#     _ins_size_is_fixed = False
#     _N_outs = 3
 #    _units= {'Aqueous flowrate': 'kg/hr',}
    
#     def __init__(self, ID='', ins=None, outs=(), thermo=None,
#                   init_with='WasteStream', F_BM_default=1,
#                   N_unit=1, **kwargs,
#                   ):
#         SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
#         self._mixed = self.ins[0].copy(f'{self.ID}_mixed')
#         self.N_unit = N_unit
#         for kw, arg in kwargs.items(): setattr(self, kw, arg)
    
#     def _run(self):
#         mixed = self._mixed
#         mixed.mix_from(self.ins)
        
#         fertilizer, water, solids = self.outs
        
#         # Just to copy the conditions of the mixture
#         for i in self.outs:
#             i.copy_like(mixed)
#             i.empty()
        
#         water.imass['Water'] = mixed.imass['Water']
#         fertilizer.copy_flow(mixed, exclude=('Water', 'Ash'))
#         solids.copy_flow(mixed, IDs=('Ash',))

#     def _design(self):
#         self.design_results['Aqueous flowrate'] = self.F_mass_in
#         self.parallel['self'] = self.N_unit
        
#     @property
#     def N_unit(self):
#         '''
#         [int] Number of filtration units.
#         '''
#         return self._N_unit
#     @N_unit.setter
#     def N_unit(self, i):
#         self.parallel['self'] = self._N_unit = math.ceil(i)

from qsdsan.equipments import Electrode, Membrane
import thermosteam as tmo

@cost(basis='Aqueous flowrate', ID= 'Anode', units='kg/hr',
       cost=1649.95, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
@cost(basis='Aqueous flowrate', ID= 'Cathode', units='kg/hr',
       cost=18, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
@cost(basis='Aqueous flowrate', ID= 'Cell Exterior', units='kg/hr',
      cost=80, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
class ElectrochemicalOxidation(qs.SanUnit):
    _N_ins = 2
    _N_outs = 3
    _units = {'Aqueous flowrate': 'kg/hr',}

    def __init__(self, ID='', ins=(), outs=(),
                 recovery={'Carbon':0.7, 'Nitrogen':0.7, 'Phosphorus':0.7},  #consult Davidson group
                 removal={'Carbon':0.83, 'Nitrogen':0.83, 'Phosphorus':0.83},  #consult Davidson group
                 OPEX_over_CAPEX=0.2, N_unit=1, F_BM_default=1.0):
        super().__init__(ID, ins, outs)
        self.recovery = recovery
        self.removal = removal
        self.OPEX_over_CAPEX = OPEX_over_CAPEX
        self.N_unit = N_unit
        self.F_BM_default = F_BM_default

        # self.equipment = [
        #     Electrode('Anode', linked_unit=self, N=1, electrode_type='anode',
        #               material='Boron-doped diamond (BDD) on niobium', surface_area=1, unit_cost=1649.95),
        #     Electrode('Cathode', linked_unit=self, N=1, electrode_type='cathode',
        #               material='Stainless Steel 316', surface_area=1, unit_cost=18.00),
        #     Membrane('Proton_Exchange_Membrane', linked_unit=self, N=1,
        #              material='Nafion proton exchange membrane', unit_cost=50, surface_area=1)
        #]

    def _run(self): 
           HTL_aqueous, catalysts = self.ins
           #aqueous_flowrate = HTL_aqueous.imass['Aqueous flowrate']
           recovered, removed, residual = self.outs

           mixture = qs.WasteStream()
           mixture.mix_from(self.ins)
           residual.copy_like(mixture)
           #solids.copy_flow(mixture, IDs=('Ash',))

           # Check chemicals present in each stream
           #print("Available chemicals in mixture:", list(mixture.imass.chemicals))

           for chemical in set(self.recovery.keys()).union(set(self.removal.keys())):
             if chemical in mixture.imass.chemicals:
                recovery_amount = mixture.imass[chemical] * self.recovery.get(chemical, 0)
                recovered.imass[chemical] = recovery_amount
            
                removal_amount = mixture.imass[chemical] * self.removal.get(chemical, 0)
                removed.imass[chemical] = removal_amount - recovery_amount
                residual.imass[chemical] -= removal_amount
             else:
                print(f"Chemical '{chemical}' not found in mixture.imass")




    def _design(self):
        self.design_results['Aqueous flowrate'] = self.F_mass_in
        self.parallel['self'] = self.N_unit
        self.add_equipment_design()

  
    @property
    def N_unit(self):
        return self._N_unit
    
    @N_unit.setter
    def N_unit(self, i):
        self.parallel['self'] = self._N_unit = math.ceil(i)


# %%      

    
class Disposal(SanUnit):
    '''
    Mix any number of influents for waste disposal.
    Price for the disposal stream is given for dry weights.
    
    Parameters
    ----------
    ins : seq(obj)
        Any number of influent streams.
    outs : seq(obj)
        Waste, others. The "waste" stream is the disposal stream for price calculation,
        the "other" stream is a dummy stream for components excluded from disposal cost calculation
        (e.g., if the cost of a wastewater stream is given based on $/kg of organics,
         the "other" stream should contain the non-organics).
    disposal_price : float
        Price for the disposal stream.
    exclude_components : seq(str)
        IDs of the components to be excluded from disposal price calculation.
    '''
    
    _ins_size_is_fixed = False
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
        self._mixed = self.ins[0].copy(f'{self.ID}_mixed')
        for kw, arg in kwargs.items(): setattr(self, kw, arg)
    
    def _run(self):
        mixed = self._mixed
        mixed.mix_from(self.ins)
        waste, others = self.outs        
        
        waste.copy_like(mixed)
        waste.imass[self.exclude_components] = 0
        
        others.copy_like(mixed)
        others.imass[self.components.IDs] -= waste.imass[self.components.IDs]
        
    def _cost(self):
        self.outs[0].price = self.disposal_price
