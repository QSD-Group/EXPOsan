#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>
    
    Ali Ahmad <aa3056@scarletmail.rutgers.edu>
    
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
from exposan.saf import _units as safu
from exposan.biobinder._process_settings import dry_flowrate

__all__ = (
    'BiocrudeSplitter',
    'CentralizedHTL',
    'Conditioning',
    'Electrochemical',
    'Hydroprocessing',
    'PilotHTL',
    'ProcessWaterCenter',
    'Scaler',
    'Transportation',
    )

_psi_to_Pa = 6894.76
CEPCI_by_year = qs.utils.tea_indices['CEPCI']

BiocrudeSplitter = safu.BiocrudeSplitter
CentralizedHTL = safu.HydrothermalLiquefaction
Transportation = safu.Transportation
ProcessWaterCenter = safu.ProcessWaterCenter

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

# Modify needed units to allow scaling
class Conditioning(safu.Conditioning):

    _N_unit = 1
    
    def _cost(self):
        safu.Conditioning._cost(self)
        self.parallel['self'] = self._N_unit

    @property
    def N_unit(self):
        '''
        [int] Number of parallel units.
        '''
        return self._N_unit
    @N_unit.setter
    def N_unit(self, i):
        self.parallel['self'] = self._N_unit = int(i)
        
        
class Hydroprocessing(safu.Hydroprocessing):

    _N_unit = 1
    
    def _cost(self):
        safu.Hydroprocessing._cost(self)
        self.parallel['self'] = self._N_unit

    @property
    def N_unit(self):
        '''
        [int] Number of parallel units.
        '''
        return self._N_unit
    @N_unit.setter
    def N_unit(self, i):
        self.parallel['self'] = self._N_unit = int(i)


class Electrochemical(safu.Electrochemical):

    _N_unit = 1
    
    def _cost(self):
        safu.Electrochemical._cost(self)
        self.parallel['self'] = self._N_unit

    @property
    def N_unit(self):
        '''
        [int] Number of parallel units.
        '''
        return self._N_unit
    @N_unit.setter
    def N_unit(self, i):
        self.parallel['self'] = self._N_unit = int(i)


# %%

@cost(basis='Feedstock dry flowrate', ID='Feedstock Tank', units='kg/hr',
      cost=4330, S=dry_flowrate, CE=CEPCI_by_year[2023], n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Feedstock Pump', units='kg/hr',
      cost=6180, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2.3)
@cost(basis='Feedstock dry flowrate', ID= 'Inverter', units='kg/hr',
      cost=240, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'High Pressure Pump', units='kg/hr',
      cost=1634, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2.3)
@cost(basis='Feedstock dry flowrate', ID= 'Reactor Core', units='kg/hr',
      cost=30740, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2)
@cost(basis='Feedstock dry flowrate', ID= 'Reactor Vessel', units='kg/hr',
      cost=4330, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Heat Transfer Putty', units='kg/hr',
      cost=2723, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Electric Heaters', units='kg/hr',
      cost=8400, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'J Type Thermocouples', units='kg/hr',
      cost=497, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Ceramic Fiber', units='kg/hr',
      cost=5154, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Steel Jacket', units='kg/hr',
      cost=22515, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Feedstock dry flowrate', ID= 'Counterflow Heat Exchanger', units='kg/hr',
      cost=14355, S=dry_flowrate, CE=CEPCI_by_year[2013],n=0.77, BM=2.2)
@cost(basis='Feedstock dry flowrate', ID= 'Temperature Control and Data Logging Unit', units='kg/hr',
      cost=905, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'Pulsation Dampener', units='kg/hr',
      cost=3000, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'Fluid Accumulator', units='kg/hr',
      cost=995, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'Burst Rupture Discs', units='kg/hr',
      cost=1100, S=dry_flowrate, CE=CEPCI_by_year[2023], n=0.77, BM=1.6)
@cost(basis='Feedstock dry flowrate', ID= 'Pressure Relief Vessel', units='kg/hr',
      cost=4363, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=2)
@cost(basis='Feedstock dry flowrate', ID= 'Gas Scrubber', units='kg/hr',
      cost=1100, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.8)
@cost(basis='Feedstock dry flowrate', ID= 'BPR', units='kg/hr',
      cost=4900, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.6)
@cost(basis='Feedstock dry flowrate', ID= 'Primary Collection Vessel', units='kg/hr',
      cost=7549, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Belt Oil Skimmer', units='kg/hr',
      cost=2632, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Bag Filter', units='kg/hr',
      cost=8800, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.7)
@cost(basis='Feedstock dry flowrate', ID= 'Oil Vessel', units='kg/hr',
      cost=4330, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1.5)
@cost(basis='Feedstock dry flowrate', ID= 'Mobile HTL system', units='kg/hr',
      cost=23718, S=dry_flowrate, CE=CEPCI_by_year[2023],n=0.77, BM=1)
@cost(basis='Non-scaling factor', ID='Magnotrol Valves Set', units='ea',
      cost=343, S=1, CE=CEPCI_by_year[2023], n=1, BM=1)
class PilotHTL(safu.HydrothermalLiquefaction):
    '''
    Pilot-scale reactor for hydrothermal liquefaction (HTL) of wet organics.
    Biocrude from mulitple pilot-scale reactors will be transported to a central plant
    for biocrude upgrading.
    
    Parameters
    ----------
    ins : Iterable(stream)
        Feedstock into HTL.
    outs : Iterable(stream)
        Gas, aqueous, biocrude, char.
    N_unit : int
        Number of parallel units.
    piping_cost_ratio : float
        Piping cost estimated as a ratio of the total reactor cost.
    accessory_cost_ratio : float
        Accessories (e.g., valves) cost estimated as a ratio of the total reactor cost.
        
    See Also
    --------    
    :class:`qsdsan.sanunits.HydrothermalLiquefaction`
    '''

    _units= {
        'Feedstock dry flowrate': 'kg/hr',
        'Non-scaling factor': 'ea',
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', include_construction=False,
                 N_unit=1,
                 piping_cost_ratio=0.15,
                 accessory_cost_ratio=0.08,
                 **kwargs,
                 ):
        
        safu.HydrothermalLiquefaction.__init__(self, ID, ins, outs, thermo, init_with, include_construction, **kwargs)
        self.N_unit = N_unit
        self.piping_cost_ratio = piping_cost_ratio
        self.accessory_cost_ratio = accessory_cost_ratio
       

    def _design(self):
        safu.HydrothermalLiquefaction._design(self)
        Design = self.design_results
        feed = self.ins[0]
        Design.clear()
        Design['Feedstock dry flowrate'] = feed.F_mass-feed.imass['Water']
        Design['Non-scaling factor'] = 1

        
    def _cost(self):
        self.parallel['self'] = self.N_unit
        all_cost_items = self.cost_items.copy()
        HTL_cost_items = safu.HydrothermalLiquefaction.cost_items
        pilot_items = {k:v for k, v in all_cost_items.items() if k not in HTL_cost_items}
        self.cost_items = pilot_items
        self._decorated_cost()
        #!!! Need to compare the externally sourced HX cost and BioSTEAM default
        # also need to make sure the centralized HTL cost is not included
        baseline_purchase_cost = self.baseline_purchase_cost
        Cost = self.baseline_purchase_costs
        Cost['Piping'] = baseline_purchase_cost*self.piping_cost_ratio
        Cost['Accessories'] = baseline_purchase_cost*self.accessory_cost_ratio

    
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

# =============================================================================
# Legacy units for the record
# =============================================================================

# @cost(basis='Biocrude flowrate', ID= 'Deashing Tank', units='kg/hr',
#       cost=4330, S=base_biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
# class BiocrudeDeashing(SanUnit):
#     '''
#     Biocrude deashing unit.
    
#     Parameters
#     ----------
#     ins : obj
#         HTL biocrude.
#     outs : seq(obj)
#         Deashed biocrude, ash for disposal.
#     '''
    
#     _N_outs = 2
#     _units= {'Biocrude flowrate': 'kg/hr',}
#     target_ash = 0.01 # dry weight basis
    
#     def __init__(self, ID='', ins=None, outs=(), thermo=None,
#                   init_with='WasteStream', F_BM_default=1,
#                   N_unit=1, **kwargs,
#                   ):
#         SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
#         self.N_unit = N_unit
#         for kw, arg in kwargs.items(): setattr(self, kw, arg)
    
#     def _run(self):
#         biocrude = self.ins[0]
#         deashed, ash = self.outs
        
#         deashed.copy_like(biocrude)
#         ash.empty()
#         dw = deashed.F_mass - deashed.imass['Water']
#         excess_ash = deashed.imass['Ash'] - dw * self.target_ash
#         # Remove excess ash
#         if excess_ash >= 0:
#             deashed.imass['Ash'] -= excess_ash
#             ash.imass['Ash'] = excess_ash
            
#     def _design(self):
#         self.design_results['Biocrude flowrate'] = self.ins[0].F_mass
#         self.parallel['self'] = self.N_unit

#     @property
#     def N_unit(self):
#         '''
#         [int] Number of deashing units.
#         '''
#         return self._N_unit
#     @N_unit.setter
#     def N_unit(self, i):
#         self.parallel['self'] = self._N_unit = math.ceil(i)
            

# @cost(basis='Biocrude flowrate', ID= 'Dewatering Tank', units='kg/hr',
#       cost=4330, S=base_biocrude_flowrate, CE=CEPCI_by_year[2023],n=0.75, BM=1.5)
# class BiocrudeDewatering(SanUnit):
#     '''
#     Biocrude dewatering unit.
    
#     Parameters
#     ----------
#     ins : obj
#         HTL biocrude.
#     outs : seq(obj)
#         Dewatered biocrude, water for treatment.
#     '''
    
#     _N_outs = 2
#     _units= {'Biocrude flowrate': 'kg/hr',}
#     target_moisture = 0.01 # weight basis
    
#     def __init__(self, ID='', ins=None, outs=(), thermo=None,
#                   init_with='WasteStream', F_BM_default=1,
#                   N_unit=1, **kwargs,
#                   ):
#         SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
#         self.N_unit = N_unit
#         for kw, arg in kwargs.items(): setattr(self, kw, arg)
    
#     def _run(self):
#         biocrude = self.ins[0]
#         dewatered, water = self.outs
        
#         dewatered.copy_like(biocrude)
#         water.empty()
#         dw = dewatered.F_mass - dewatered.imass['Water']
#         excess_water = dw/(1-self.target_moisture) - dw
#         # Remove excess water
#         if excess_water >= 0:
#             dewatered.imass['Water'] -= excess_water
#             water.imass['Water'] = excess_water
            
#     def _design(self):
#         self.design_results['Biocrude flowrate'] = self.ins[0].F_mass
#         self.parallel['self'] = self.N_unit
        
#     @property
#     def N_unit(self):
#         '''
#         [int] Number of dewatering units.
#         '''
#         return self._N_unit
#     @N_unit.setter
#     def N_unit(self, i):
#         self.parallel['self'] = self._N_unit = math.ceil(i)


# @cost(basis='Aqueous flowrate', ID= 'Sand Filtration Unit', units='kg/hr',
#       cost=318, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.7)
# @cost(basis='Aqueous flowrate', ID= 'EC Oxidation Tank', units='kg/hr',
#       cost=1850, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
# @cost(basis='Aqueous flowrate', ID= 'Biological Treatment Tank', units='kg/hr',
#       cost=4330, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
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
#     _units= {'Aqueous flowrate': 'kg/hr',}
    
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


# #!!! need to add utility, etc.
# base_ap_flowrate = 49.65 #kg/hr

# @cost(basis='Aqueous flowrate', ID= 'Anode', units='kg/hr',
#        cost=1649.95, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
# @cost(basis='Aqueous flowrate', ID= 'Cathode', units='kg/hr',
#        cost=18, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
# @cost(basis='Aqueous flowrate', ID= 'Cell Exterior', units='kg/hr',
#       cost=80, S=base_ap_flowrate, CE=CEPCI_by_year[2023],n=0.65, BM=1.5)
# class ElectrochemicalOxidation(qs.SanUnit):
#     _N_ins = 2
#     _N_outs = 3
#     _units = {'Aqueous flowrate': 'kg/hr',}

#     def __init__(self, ID='', ins=(), outs=(),
#                  recovery={'Carbon':0.7, 'Nitrogen':0.7, 'Phosphorus':0.7},  #consult Davidson group
#                  removal={'Carbon':0.83, 'Nitrogen':0.83, 'Phosphorus':0.83},  #consult Davidson group
#                  OPEX_over_CAPEX=0.2, N_unit=1, F_BM_default=1.0):
#         super().__init__(ID, ins, outs)
#         self.recovery = recovery
#         self.removal = removal
#         self.OPEX_over_CAPEX = OPEX_over_CAPEX
#         self.N_unit = N_unit
#         self.F_BM_default = F_BM_default


#     def _run(self): 
#            HTL_aqueous, catalysts = self.ins
#            #aqueous_flowrate = HTL_aqueous.imass['Aqueous flowrate']
#            recovered, removed, residual = self.outs

#            mixture = qs.WasteStream()
#            mixture.mix_from(self.ins)
#            residual.copy_like(mixture)
#            #solids.copy_flow(mixture, IDs=('Ash',))

#            # Check chemicals present in each stream
#            #print("Available chemicals in mixture:", list(mixture.imass.chemicals))

#            for chemical in set(self.recovery.keys()).union(set(self.removal.keys())):
#              if chemical in mixture.imass.chemicals:
#                 recovery_amount = mixture.imass[chemical] * self.recovery.get(chemical, 0)
#                 recovered.imass[chemical] = recovery_amount
            
#                 removal_amount = mixture.imass[chemical] * self.removal.get(chemical, 0)
#                 removed.imass[chemical] = removal_amount - recovery_amount
#                 residual.imass[chemical] -= removal_amount
#              else:
#                 print(f"Chemical '{chemical}' not found in mixture.imass")


#     def _design(self):
#         self.design_results['Aqueous flowrate'] = self.F_mass_in
#         self.parallel['self'] = self.N_unit
#         self.add_equipment_design()

  
#     @property
#     def N_unit(self):
#         return self._N_unit
#     @N_unit.setter
#     def N_unit(self, i):
#         self.parallel['self'] = self._N_unit = math.ceil(i)