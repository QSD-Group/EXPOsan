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
import biosteam as bst, flexsolve as flx
from math import ceil, sqrt, pi
from .splitting import Splitter
from biosteam.units.design_tools import PressureVessel
import numpy as np
from numba import njit

__all__ = (
    'ALFProduction',
    'ALFCrystallizer',
    'ALFTSA'
    )

# =============================================================================
# CO2 absorber/stripper
# 
# electrochemical cells
# =============================================================================

# =============================================================================
# ALFProduction
# =============================================================================

class ALFProduction(bst.CSTR):
    
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
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 T, crystal_ALF_yield=1, order=None):
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
# ALFTSA
# =============================================================================
class ALFTSA(PressureVessel, Splitter):
    """
    
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
        Time at which the receiving vessel is switched. Defaults to 8 h [X]_.
    rho_adsorbent_solid : float, optional
        The density of ALF. Defaults to 1441 [kg/m3] Table S1 [Y]_.
    adsorbent_capacity : float, optional
        Amount of CO2 that ALF can hold. Defaults to 2.7 mmol/g Table S7 [Y]_.
    T_regeneration : float, optional
        Temperature during the regeneration phase. Defaults to 418 K Table S8 [Y]_.
    vessel_material : float, optional
        Vessel material. Defaults to 'Stainless steel 316',
    vessel_type : float, optional
        Vessel type. Defaults to 'Vertical'.
        
        
        
        
        
        
        
        
        
    regeneration_fluid : dict[str, float]
        Arguments to initialize fluid used to regenerate the bed.

    
    References
    ----------
    [1] Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
    [2] Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
    [3] 
    [4] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani,
        R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)
        
    [X] https://www.chemicalprocessing.com/processing-equipment/fluid-handling/article/11302111/select-the-right-valves-for-adsorption-processes-chemical-processing
    [Y] Evans et al. Science Advances SI.
    
    """
    auxiliary_unit_names = ('heat_exchanger_regeneration')
    
    # in $/ft3
    adsorbent_cost = XXX
    
    # in year
    _default_equipment_lifetime = XXX
    
    _N_ins = 2
    _N_outs = 3
    
    def _init(self,
              split=dict(O2=0, N2=0, CO2=1),
              superficial_velocity=1080, # m/h
              regeneration_velocity=1080, # m/h # TODO: need to decide if this is needed since there is no carrier gas stream
              cycle_time=8, # h
              rho_adsorbent_solid = 1441, # kg/m3
              adsorbent_capacity=2.7, # mmol/g
              T_regeneration=418, # K
              vessel_material='Stainless steel 316',
              vessel_type='Vertical',
              
              
              void_fraction=0.47, # Only matters when K given; 0.30 - 0.35 for activated carbon
              length_unused=1.219, # Additional length of a column to account for mass transfer limitations (due to unused bed). Defaults to +2 ft per column.
              
              
              
              waste_ratio=0.07, # ratio of ALF that is wasted (and is sent for CO2 storage) each cycle # TODO: confirm this is for each cycle
              

              ):
        bst.Splitter._init(self, split=split)
        self.superficial_velocity = superficial_velocity
        self.regeneration_velocity = regeneration_velocity
        self.cycle_time = cycle_time
        self.adsorbent_capacity = adsorbent_capacity
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.T_regeneration = T_regeneration

        self.void_fraction = void_fraction
        self.length_unused = length_unused


        self.rho_adsorbent_solid = rho_adsorbent_solid
        self.rho_adsorbent = rho_adsorbent_solid * (1-void_fraction)
        self.heat_exchanger_regeneration = bst.HXutility(None, None, None, thermo=self.thermo)
        
    @property
    def effluent(self):
        return self.outs[0]
        
    @property
    def regeneration_purge(self):
        return self.outs[1]
    
    def _run(self):
        feed, ALF = self.ins
        offgas, CO2, used_ALF = self.outs 
        
        
        
        regen.empty()
        dry_air.empty()
        for i in self.outs: i.empty()

        feed.split_to(effluent, purge, self.split)
        F_vol_feed = feed.F_vol
        superficial_velocity = self.superficial_velocity
        adsorbate_ID = self.adsorbate_ID
        F_mass_adsorbate = purge.imass[adsorbate_ID]

       
        
        self.diameter = diameter = 2 * sqrt(F_vol_feed / (superficial_velocity * pi))
        self.area = area = pi * diameter * diameter / 4
        total_length = (
            self.cycle_time * F_mass_adsorbate / (self.adsorbent_capacity * self.rho_adsorbent * area)
        ) + self.length_unused # length of equilibrium section plus unused bed (LES + LUB)
        self.length = length = total_length / 2 # Size of each column
        self.vessel_volume = length * area
        T_original = regen.T
        regen.reset_flow(**self.regeneration_fluid)
        purge.T = regen.T = self.T_regeneration
        regen.F_vol = area * self.regeneration_velocity
        regen.T = T_original
        
        
    
    def _design(self):
        feed, regen, dry_air = self.ins
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
        hxr.ins[0] = regen.copy()
        hxr.outs[0] = regen.copy()
        hxr.T = self.T_regeneration
        hxr.simulate()
    
    def _cost(self):
        design_results = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        baseline_purchase_costs.update(self._vessel_purchase_cost(
            design_results['Weight'], design_results['Diameter'], design_results['Length']))
        N_reactors = design_results['Number of reactors']
        for i, j in baseline_purchase_costs.items():
            baseline_purchase_costs[i] *= N_reactors
        baseline_purchase_costs[self.adsorbent] = N_reactors * 35.3147 * self.vessel_volume * self.adsorbent_cost[self.adsorbent]

# =============================================================================
# class XXX():
#     '''
#     XXX.
#     
#     Parameters
#     ----------
#     ins : Iterable(stream)
#         XXX, XXX.
#     outs : Iterable(stream)
#         XXX, XXX.
#         
#     References
#     ----------
# 
#     '''
#     _N_ins = 
#     _N_outs = 
#         
#     def __init__(self, ID='', ins=None, outs=(), thermo=None,
#                  init_with='WasteStream'):
#         
#         SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
#         self.XXX = XXX
# 
#     def _run(self):
#         
#         XXX, XXX = self.ins
#         XXX, XXX = self.outs
# 
#     def _design(self):
#         
#     def _cost(self):
#         
#     @property
#     def XXX(self):
#         return
# =============================================================================