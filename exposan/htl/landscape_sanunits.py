#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

The following units are fully or partially based on the BEAM*2024 model:
'Storage','Thickening','AerobicDigestion','AnaerobicDigestion','Dewatering',
'AlkalineStabilization','Composting','HeatDrying','Incineration','Pyrolysis',
'Landfilling','LandApplication'.

BEAM*2024 model:
North East Biosolids and Residuals Association (NEBRA), Northern Tilth LLC, and 
orthwest Biosolids. Estimating Greenhouse Gas Emissions from Biosolids Management.
BEAM*2024 Spreadsheet Model and Supporting Information, 2024.
https://www.BiosolidsGHGs.org (accessed 2025-06-16).

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, biosteam as bst
from qsdsan import SanUnit
from qsdsan.sanunits import SludgePump, WWTpump, HXutility
from qsdsan.utils import auom, calculate_excavation_volume
from qsdsan.equipments import Blower, GasPiping
from biosteam import Stream
from biosteam.units.decorators import cost
from biosteam.units.design_tools import CEPCI_by_year
from math import ceil, pi

# TODO: check natural gas CI, whether including upstream

# methane density, kg/m3
methane_density = 0.707

# methane heat content, BTU/m3
methane_heat = 35830

# natural gas heat content, BTU/m3
natural_gas_heat = 36263

_ton_to_tonne = auom('ton').conversion_factor('tonne')
_BTU_to_kWh = auom('BTU').conversion_factor('kWh')
_BTU_to_MJ = auom('BTU').conversion_factor('MJ')
_m_to_ft = auom('m').conversion_factor('ft')
_C_to_K = 273.15
_N_to_N2O = 44/28
_C_to_CH4 = 16/12
_C_to_CO2 = 44/12
_CaCO3_to_C = 12/100

folder = os.path.dirname(__file__)

# numbers without citations are likely from the BEAM*2024 model

# TODO: add pump for effluents of appropriate units

# TODO: add models
# TODO: HTL-based system
# TODO: HALT-based system
# TODO: SCWO
# TODO: gasification

# TODO: update
__all__ = (
    'Storage',
    'Thickening',
    'AerobicDigestion',
    'AnaerobicDigestion',
    'Dewatering',
    'AlkalineStabilization',
    'Composting',
    'HeatDrying',
    'Incineration',
    'Pyrolysis',
    'Landfilling',
    'LandApplication',
    )

# =============================================================================
# Storage
# =============================================================================

class Storage(SanUnit):
    '''
    Raw sludge long-term storage in a lagoon, with the assumption that the solids mass and
    moisture content does not change during this step.
    
    scope 1 emission: fugitive methane
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        raw_sludge.
    outs : iterable
        stored_sludge, fugitive_methane, vapor.
    BOD : float
        BOD of raw sludge, [mg/L].
    HRT : float
        hydraulic retention time, [day].
    aerated : bool
        aerated or anaerobic lagoon, [True, False].
    # TODO: confirm unit_electricity is based on the volume of the lagoon
    unit_electricity : float
        electricity for aeration (if aerated), [kW/m3].
    percentage_above_15_C : float
        percentage days per year with a temperature above 15 °C, [-].
    lagoon_depth : float
        the depth of the lagoon, [m].
    '''
    _N_ins = 1
    _N_outs = 3
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', BOD=15000, HRT=100,
                 aerated=False, unit_electricity=0.0056,
                 percentage_above_15_C=0.5, lagoon_depth=3.5):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.BOD = BOD
        self.HRT = HRT
        self.aerated = aerated
        self.unit_electricity = unit_electricity
        self.percentage_above_15_C = percentage_above_15_C
        self.lagoon_depth = lagoon_depth

    def _run(self):
        raw_sludge = self.ins[0]
        stored_sludge, fugitive_methane, vapor = self.outs
        
        stored_sludge.phase = 'l'
        fugitive_methane.phase = 'g'
        
        # kg BOD/day
        BOD_to_storage = raw_sludge.F_vol*1000*24*self.BOD/1000000
        
        if self.aerated:
            fugitive_methane.imass['CH4'] = 0
        else:
            if self.lagoon_depth <= 2:
                # kg CH4/kg BOD
                self.CH4_EF = 0.12
            else:
                self.CH4_EF = 0.48
            
            fugitive_methane.imass['CH4'] = BOD_to_storage*self.CH4_EF*self.percentage_above_15_C/24
        
        stored_sludge.copy_like(raw_sludge)
    
    def _design(self):
        if self.aerated:
            # kW
            self.add_power_utility(self.ins[0].F_vol*24*self.HRT*self.unit_electricity)
    
    @property
    def moisture(self):
        return self.outs[0].imass['H2O']/self.outs[0].F_mass
    
    @property
    def unit_process(self):
        return 'storage'

# =============================================================================
# Thickening
# =============================================================================

class Thickening(SanUnit):
    '''
    Thickening of sludge using either centrifugal or other types of thickener,
    with the assumption that polymer can be ignored in the mass balance.
    
    scope 1 emission: N/A
    scope 2 emission: electricity
    scope 3 emission: polymer
    
    Parameters
    ----------
    ins : iterable
        input_sludge, polymer.
    outs : iterable
        thickened_sludge, reject.
    target_moisture : float
        targeted moisture content, [-].
    polymer_addition : float
        polymer used for conditioning sludge, [kg polymer/dry tonne sludge].
    unit_electricity : float
        electricity for thickening (probably not including pumps), [kWh/dry tonne].
        typically, 33 for centrifugal thickening, 4.9 for others.
    max_capacity : float
        maximum hydraulic loading per belt thickener, [m3/h].
    '''
    _N_ins = 2
    _N_outs = 2
    auxiliary_unit_names = ('effluent_pump','sludge_pump')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', target_moisture=0.97,
                 polymer_addition=5, unit_electricity=4.9, max_capacity=100):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target_moisture = target_moisture
        self.polymer_addition = polymer_addition
        self.unit_electricity = unit_electricity
        self.max_capacity = max_capacity
        ID = self.ID
        sludge = self.outs[0].proxy(f'{ID}_sludge')
        eff = self.outs[1].proxy(f'{ID}_eff')
        self.sludge_pump = SludgePump(f'.{ID}_sludge_pump', ins=sludge, init_with=init_with)
        self.effluent_pump = SludgePump(f'.{ID}_eff_pump', ins=eff, init_with=init_with)
    
    def _run(self):
        input_sludge, polymer = self.ins
        thickened_sludge, reject = self.outs
        
        thickened_sludge.phase = 'l'
        reject.phase = 'l'
        
        # dry tonne/h
        self.dry_solids = (input_sludge.F_mass - input_sludge.imass['H2O'])/1000
        
        polymer.imass['PAM'] = self.dry_solids*self.polymer_addition
        
        # tonne/h
        thickened_sludge_mass_flow = self.dry_solids/(1 - self.target_moisture)
        
        thickened_sludge.copy_like(input_sludge)
       
        thickened_sludge.imass['H2O'] = thickened_sludge_mass_flow*1000*self.target_moisture
       
        reject.imass['H2O'] = input_sludge.F_mass - thickened_sludge.F_mass
    
    def _design(self):
        self._N_thickener = N = ceil(self.ins[0].F_vol/self.max_capacity)
        self.design_results['Number of thickeners'] = N
        self.F_BM['Thickeners'] = 1.7
        self.baseline_purchase_costs['Thickeners'] = 4000*N
        
        # kW
        self.add_power_utility(self.dry_solids*self.unit_electricity)
    
    def _cost(self):
        # add electricity consumption of pumps
        for p in (self.effluent_pump, self.sludge_pump): p.simulate()
    
    @property
    def moisture(self):
        return self.target_moisture
    
    @property
    def unit_process(self):
        return 'thickening'

# =============================================================================
# AerobicDigestion
# =============================================================================

class AerobicDigestion(SanUnit):
    '''
    Aerobic digestion of thickened (and conditioned) sludge, with the assumption
    that no water evaporizes and no natural gas is needed.
    
    scope 1 emission: N/A
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    Capital costs include tank, pump, and aeration system.
    
    Parameters
    ----------
    ins : iterable
        input_sludge.
    outs : iterable
        digested_sludge.
    O2_requirement : float
        O2 requirement for aerobic digestion, [kg O2/kg VS destroyed].
        2.3 kg O2/kg VS destroyed, [1].
    VS_reduction : float
        volatile solids reduction after aerobic digestion, [-].
    HRT : float
        hydraulic retention time of aerbic digestion, [day].
        14 days, [2].
    SRT : float
        solids retention time of aerobic digestion, [day].
    unit_electricity : float
        electricity for aerobic digestion, [kW/m3 reactor volume].
    depth : float
        side depth of the digester, [m].
    wall_concrete_unit_cost : float
        unit cost of the wall concrete, [$/ft3].
    slab_concrete_unit_cost : float
        unit cost of the slab concrete, [$/ft3].
    excavation_unit_cost : float
        unit cost of the excavation activity, [$/ft3].
    
    References
    ----------
    .. [1] Metcalf & Eddy, I. Wastewater Engineering : Treatment and Reuse;
           Fifth edition / revised by George Tchobanoglous, Franklin L.
           Burton, H. David Stensel. Boston : McGraw-Hill, 2014.
    .. [2] Wang, Q.; Yuan, Z. Enhancing Aerobic Digestion of Full-Scale Waste
           Activated Sludge Using Free Nitrous Acid Pre-Treatment. RSC Adv.
           2015, 5 (25), 19128–19134. https://doi.org/10.1039/C4RA17215A.
    '''
    _N_ins = 2
    _N_outs = 2
    
    # K
    _T_air = 17 + 273.15
    _T_earth = 10 + 273.15
    
    # ft
    _freeboard = 3
    _t_wall = 6/12
    _t_slab = 8/12
    
    # pump building, ft
    _L_PB = 50
    _W_PB = 30
    _D_PB = 10
    
    # excavation
    # horizontal/vertical
    _excav_slope = 1.5
    # ft
    _constr_access = 3
    
    auxiliary_unit_names = ('heat_exchanger','sludge_pump')
    
    F_BM_pump = 1.18*(1+0.007/100) # 0.007 is for miscellaneous costs
    # BioSTEAM: BM = 2.3 for all tanks
    default_F_BM = {'Pump': F_BM_pump,
                    'Pump building': F_BM_pump,
                    'Wall concrete': 2.3,
                    'Slab concrete': 2.3,
                    'Excavation': 2.3,
                    'Blower air piping': 1,
                    'Blowers': 2.22,
                    'Blower building': 1.11}
    # TODO: determine if this is needed
    default_equipment_lifetime = {'Pump': 15,
                                  'Pump pipe stainless steel': 15,
                                  'Pump stainless steel': 15}
    
    _units = {'HRT':'d',
              'SRT':'d',
              'Volume':'m3',
              'Surface area':'m2',
              'Depth':'m',
              'Diameter':'m',
              'Wall concrete':'ft3',
              'Slab concrete':'ft3',
              'Excavation':'ft3',
              'Pump pipe stainless steel':'kg',
              'Pump stainless steel':'kg'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', O2_requirement=2.3,
                 VS_reduction=0.475, HRT=14, SRT=40,
                 # TODO: this includes aeration, but there is a separate blower (check electricity)
                 unit_electricity=0.03, depth=10, wall_concrete_unit_cost=24,
                 slab_concrete_unit_cost=13, excavation_unit_cost=0.3,
                 F_BM=default_F_BM, lifetime=default_equipment_lifetime):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.O2_requirement = O2_requirement
        self.VS_reduction = VS_reduction
        self.HRT = HRT
        self.SRT = SRT
        self.unit_electricity = unit_electricity
        self.depth = depth
        self.wall_concrete_unit_cost = wall_concrete_unit_cost
        self.slab_concrete_unit_cost = slab_concrete_unit_cost
        self.excavation_unit_cost = excavation_unit_cost
        self.F_BM.update(F_BM)
        # TODO: determine if this is needed
        self._default_equipment_lifetime.update(lifetime)
        ID = self.ID
        self.sludge_pump = WWTpump(ID=f'{ID}_sludge', ins=self.ins[0].proxy(),
                                   pump_type='', Q_mgd=None, add_inputs=(1,),
                                   capacity_factor=1., include_pump_cost=True,
                                   include_building_cost=False,
                                   include_OM_cost=False)
        self.blower = blower = Blower('blower', linked_unit=self, N_reactor=1)
        self.air_piping = air_piping = GasPiping('air_piping', linked_unit=self, N_reactor=1)
        self.equipments = (blower, air_piping)
    
    def _run(self):
        input_sludge, air = self.ins
        digested_sludge, offgas = self.outs
        
        air.phase = 'g'
        digested_sludge.phase = 'l'
        offgas.phase = 'g'
        
        # kg/h
        VS_destroyed = (input_sludge.F_mass - input_sludge.imass['H2O'])*self.VS_reduction
        
        air.imass['O2'] = VS_destroyed*self.O2_requirement
        # 23.2% O2 in air by weight, [1]
        air.imass['N2'] = air.imass['O2']/0.232*(1 - 0.232)
        
        digested_sludge.imass['H2O'] = input_sludge.imass['H2O']
        digested_sludge.imass['Sludge_ash'] = input_sludge.imass['Sludge_ash']
        digested_sludge.imass['Sludge_lipid'] = input_sludge.imass['Sludge_lipid']*(1 - self.VS_reduction)
        digested_sludge.imass['Sludge_protein'] = input_sludge.imass['Sludge_protein']*(1 - self.VS_reduction)
        digested_sludge.imass['Sludge_carbo'] = input_sludge.imass['Sludge_carbo']*(1 - self.VS_reduction)
        
        # use CO2 to represent offgas
        offgas.imass['CO2'] = input_sludge.F_mass - digested_sludge.F_mass
    
    def _design(self):
        D = self.design_results
        Q = self.ins[0].F_vol*24
        
        # dimensions
        D['HRT'] = self.HRT
        D['SRT'] = self.SRT
        V = D['Volume'] = Q*self.HRT
        depth = D['Depth'] = self.depth
        A = D['Surface area'] = V/depth
        diameter = D['Diameter']= (A*4/pi)**0.5
        
        # concrete usage
        D['Wall concrete'] = self.t_wall*pi*(diameter*_m_to_ft)*(depth*_m_to_ft + self.freeboard)
        D['Slab concrete'] = 2*self.t_slab*A*(_m_to_ft**2)
        
        # excavation
        D['Excavation'] = calculate_excavation_volume(self.L_PB,
                                                           self.W_PB,
                                                           self.D_PB,
                                                           self.excav_slope,
                                                           self.constr_access)
        
        # TODO: is the sludge pump needed for other units?
        # pump
        # TODO: does this add the electricity from the pump (make sure it is far less than the unit_electricity, or make sure unit_electricity does not include pumps)
        sludge_pump = self.sludge_pump
        sludge_pump.simulate()
        D.update(sludge_pump.design_results)
        
        # TODO: just need their cost, confirm the parameter unit_electricty include electricity from aeration (AeD) and mixing (AD), but does not include pumping
        # blower and gas piping
        air_cfm = auom('m3/h').convert(self.ins[1].F_vol, 'cfm')
        blower, piping = self.equipments
        blower.N_reactor = piping.N_reactor = 1
        blower.gas_demand_per_reactor = piping.gas_demand_per_reactor = air_cfm
        self.add_equipment_design()
        
        D.update(blower.design_results)
        D.update(piping.design_results)
        
        # kW
        self.electricity_requirement = V*self.unit_electricity
        
        # kW
        self.add_power_utility(self.electricity_requirement)

    def _cost(self):
        D, C = self.design_results, self.baseline_purchase_costs
        blower, piping = self.equipments
        
        C['Wall concrete'] = D['Wall concrete']*self.wall_concrete_unit_cost
        C['Slab concrete'] = D['Slab concrete']*self.slab_concrete_unit_cost
        C['Excavation'] = D['Excavation']*self.excavation_unit_cost
        
        # note blower._cost() includes piping
        C.update(blower._cost())
        # C.update(piping._cost())
    
    @property
    def moisture(self):
        return self.outs[0].imass['H2O']/self.outs[0].F_mass
    
    @property
    def unit_process(self):
        return 'aerobic_digestion'
    
    @property
    def T_air(self):
        '''[float] temperature of the air, [K].'''
        return self._T_air
    @T_air.setter
    def T_air(self, i):
        self._T_air = i

    @property
    def T_earth(self):
        '''[float] temperature of the air, [K].'''
        return self._T_earth
    @T_earth.setter
    def T_earth(self, i):
        self._T_earth = i
    
    @property
    def t_wall(self):
        '''[float] concrete wall thickness, [ft].'''
        return self._t_wall
    @t_wall.setter
    def t_wall(self, i):
        self._t_wall = i
    
    @property
    def freeboard(self):
        '''[float] freeboard added to the depth of the reactor tank, [ft].'''
        return self._freeboard
    @freeboard.setter
    def freeboard(self, i):
        self._freeboard = i
    
    @property
    def t_slab(self):
        '''
        [float] concrete slab thickness, [ft],
        default to be 2 in thicker than the wall thickness.
        '''
        return self._t_slab or self.t_wall + 2/12
    @t_slab.setter
    def t_slab(self, i):
        self._t_slab = i
    
    @property
    def L_PB(self):
        '''[float] length of the pump building, [ft].'''
        return self._L_PB
    @L_PB.setter
    def L_PB(self, i):
        self._L_PB = i

    @property
    def W_PB(self):
        '''[float] width of the pump building, [ft].'''
        return self._W_PB
    @W_PB.setter
    def W_PB(self, i):
        self._W_PB = i

    @property
    def D_PB(self):
        '''[float] depth of the pump building, [ft].'''
        return self._D_PB
    @D_PB.setter
    def D_PB(self, i):
        self._D_PB = i

    @property
    def excav_slope(self):
        '''[float] slope for excavation (horizontal/vertical).'''
        return self._excav_slope
    @excav_slope.setter
    def excav_slope(self, i):
        self._excav_slope = i

    @property
    def constr_access(self):
        '''[float] extra room for construction access, [ft].'''
        return self._constr_access
    @constr_access.setter
    def constr_access(self, i):
        self._constr_access = i

# =============================================================================
# AnaerobicDigestion
# =============================================================================

class AnaerobicDigestion(SanUnit):
    '''
    Anaerobic digestion of thickened (and conditioned) sludge, with the assumption
    that no water evaporizes.
    
    scope 1 emission: fugitive methane, natural gas combustion
    scope 2 emission: electricity, natural gas upstream
    scope 3 emission: N/A
    
    Capital costs include tank, pump, mixer, heat exchanger, and gas collection system.
    
    Parameters
    ----------
    ins : iterable
        input_sludge, natural_gas.
    outs : iterable
        digested_sludge, fugitive_methane.
    VS_reduction : float
        volatile solids reduction after aerobic digestion, [-].
    biogas_yield : float
        biogas yield from destroyed VS, [m3 biogas/kg destroyed VS].
    methane_biogas : float
        methane content in biogas, [-].
    biogas_CHP_ratio : float
        ratio of biogas sent to CHP, [-].
    biogas_RNG_ratio : float
        ratio of biogas sent to RNG pipelines, [-].
    biogas_flare_ratio : float
        ratio of biogas flared, [-].
    flare_fugitive_ratio : float
        ratio of fugitive biogas during flaring, [-].
        typically, 0.01 for default, 0 for enclosed, 0.05 for candlestick.
    RNG_parasitic : float
        parasitic load from natural gas conditioning and injection, [-].
    # TODO: double check unit_heat is for m3 solids but not m3 reactor volume
    unit_heat : float
        heat for anaerobic digestion, [m3 natural gas/m3 solids].
    net_capacity_factor : float
        actual electricity generated to the theoretical maximum, [-].
    HRT : float
        hydraulic retention time of aerbic digestion, [day].
    SRT : float
        solids retention time of aerobic digestion, [day].
        typically, SRT = HRT in anaerobic digestion.
    unit_electricity : float
        electricity for anaerobic digestion, [kW/m3 reactor volume].
    depth : float
        side depth of the digester, [m].
    T : float
        temperature of anaerobic digestion, [°C].
    heat_transfer_coefficient : dict
        hat transfer coefficients for heat loss calculation, [W/m2/°C].
        keys should contain 'wall', 'floor', and 'ceiling'.
    wall_concrete_unit_cost : float
        unit cost of the wall concrete, [$/ft3].
    slab_concrete_unit_cost : float
        unit cost of the slab concrete, [$/ft3].
    excavation_unit_cost : float
        unit cost of the excavation activity, [$/ft3].
    gas_collection_cost_factor : float
        capital cost ratio of the gas collection system and the tank
        excluding execavation, [-].
    vertical_mixer_unit_power : float
        the power of one verticle mixer, [kW/unit].
    vertical_mixer_unit_price : float
        the price of one verticle mixer, [$/unit].
    '''
    _N_ins = 1
    _N_outs = 3
    
    # K
    _T_air = 17 + 273.15
    _T_earth = 10 + 273.15
    
    # ft
    _freeboard = 3
    _t_wall = 6/12
    _t_slab = 8/12
    
    # pump building, ft
    _L_PB = 50
    _W_PB = 30
    _D_PB = 10
    
    # excavation
    # horizontal/vertical
    _excav_slope = 1.5
    # ft
    _constr_access = 3
    
    auxiliary_unit_names = ('heat_exchanger','sludge_pump')
    
    F_BM_pump = 1.18*(1+0.007/100) # 0.007 is for miscellaneous costs
    # BioSTEAM: BM = 2.3 for all tanks
    default_F_BM = {'Pump': F_BM_pump,
                    'Pump building': F_BM_pump,
                    'Wall concrete': 2.3,
                    'Slab concrete': 2.3,
                    'Excavation': 2.3,
                    'Gas collection system': 2.3,
                    'Vertical mixer': 1}
    # TODO: determine if this is needed
    default_equipment_lifetime = {'Pump': 15,
                                  'Pump pipe stainless steel': 15,
                                  'Pump stainless steel': 15}
    
    _units = {'HRT':'d',
              'SRT':'d',
              'Volume':'m3',
              'Surface area':'m2',
              'Depth':'m',
              'Diameter':'m',
              'Wall concrete':'ft3',
              'Slab concrete':'ft3',
              'Excavation':'ft3',
              'Pump pipe stainless steel':'kg',
              'Pump stainless steel':'kg'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', VS_reduction=0.425, biogas_yield=0.9,
                 methane_biogas=0.65, biogas_CHP_ratio=0.55, biogas_RNG_ratio=0.35,
                 biogas_flare_ratio=0.09,
                 flare_fugitive_ratio=0.01, RNG_parasitic=0.169,
                 net_capacity_factor=0.85, HRT=22, SRT=22, unit_electricity=0.0065,
                 depth=10, T=35+273.15,
                 heat_transfer_coefficient=dict(wall=0.7, floor=1.7, ceiling=0.95),
                 wall_concrete_unit_cost=24, slab_concrete_unit_cost=13,
                 excavation_unit_cost=0.3,
                 # TODO: 0.2 is an arbitrary value, replace it if there is a number with citaions
                 gas_collection_cost_factor=0.2,
                 # TODO: check if vertical_mixer_unit_price include installation 
                 vertical_mixer_unit_power=3.7, vertical_mixer_unit_price=10200,
                 F_BM=default_F_BM, lifetime=default_equipment_lifetime):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.VS_reduction = VS_reduction
        self.biogas_yield = biogas_yield
        self.methane_biogas = methane_biogas
        self.biogas_CHP_ratio = biogas_CHP_ratio
        self.biogas_RNG_ratio = biogas_RNG_ratio
        self.biogas_flare_ratio = biogas_flare_ratio
        self.flare_fugitive_ratio = flare_fugitive_ratio
        self.RNG_parasitic = RNG_parasitic
        self.net_capacity_factor = net_capacity_factor
        self.HRT = HRT
        self.SRT = SRT
        self.unit_electricity = unit_electricity
        self.depth = depth
        self.T = T
        self.heat_transfer_coefficient = heat_transfer_coefficient
        self.wall_concrete_unit_cost = wall_concrete_unit_cost
        self.slab_concrete_unit_cost = slab_concrete_unit_cost
        self.excavation_unit_cost = excavation_unit_cost
        self.gas_collection_cost_factor = gas_collection_cost_factor
        self.vertical_mixer_unit_power = vertical_mixer_unit_power
        self.vertical_mixer_unit_price = vertical_mixer_unit_price
        self.F_BM.update(F_BM)
        # TODO: determine if this is needed
        self._default_equipment_lifetime.update(lifetime)
        ID = self.ID
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'{ID}_hx', ins=hx_in, outs=hx_out)
        self.sludge_pump = WWTpump(ID=f'{ID}_sludge', ins=self.ins[0].proxy(),
                                   pump_type='', Q_mgd=None, add_inputs=(1,),
                                   capacity_factor=1., include_pump_cost=True,
                                   include_building_cost=False,
                                   include_OM_cost=False)
    
    def _run(self):
        input_sludge = self.ins[0]
        digested_sludge, natural_gas, fugitive_methane = self.outs
        
        digested_sludge.phase = 'l'
        fugitive_methane.phase = 'g'
        
        digested_sludge.imass['H2O'] = input_sludge.imass['H2O']
        digested_sludge.imass['Sludge_ash'] = input_sludge.imass['Sludge_ash']
        digested_sludge.imass['Sludge_lipid'] = input_sludge.imass['Sludge_lipid']*(1 - self.VS_reduction)
        digested_sludge.imass['Sludge_protein'] = input_sludge.imass['Sludge_protein']*(1 - self.VS_reduction)
        digested_sludge.imass['Sludge_carbo'] = input_sludge.imass['Sludge_carbo']*(1 - self.VS_reduction)
        
        # kg/h
        VS_destroyed = (input_sludge.F_mass - input_sludge.imass['H2O'])*self.VS_reduction
        
        # m3/h
        biogas_volume = VS_destroyed*self.biogas_yield
        
        # m3/h
        methane_volume = biogas_volume*self.methane_biogas
        
        self.biogas_fugitive_ratio = 1 - self.biogas_CHP_ratio - self.biogas_RNG_ratio - self.biogas_flare_ratio
        
        assert 0 <= self.biogas_fugitive_ratio <= 1
        
        # from biogas
        fugitive_methane.imass['CH4'] = methane_volume*self.biogas_fugitive_ratio*methane_density
        
        # TODO: check += here
        # from flaring
        fugitive_methane.imass['CH4'] += methane_volume*self.biogas_flare_ratio*self.flare_fugitive_ratio*methane_density        
        
        # TODO: assume all generated natural gas is sent to CHP, use a heat exchanger to meet the heating need
        # m3/h
        natural_gas_generated = methane_volume*self.biogas_CHP_ratio +\
                                methane_volume*self.biogas_RNG_ratio*(1 - self.RNG_parasitic)
        
        # TODO: assume all generated natural gas is sent to CHP, use a heat exchanger to meet the heating need
        # TODO: need to add price and CI for natural_gas
        natural_gas.imass['CH4'] = natural_gas_generated*methane_density
    
    def _design(self):
        D = self.design_results
        Q = self.ins[0].F_vol*24
        
        # dimensions
        D['HRT'] = self.HRT
        D['SRT'] = self.SRT
        V = D['Volume'] = Q*self.HRT
        depth = D['Depth'] = self.depth
        A = D['Surface area'] = V/depth
        diameter = D['Diameter']= (A*4/pi)**0.5
        
        # calculate needed heating
        T = self.T
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.copy_flow(self.ins[0])
        hx_outs0.copy_flow(self.ins[0])
        hx_ins0.T = self.ins[0].T
        hx_outs0.T = T
        hx_ins0.P = hx_outs0.P = self.ins[0].P
        
        # heat loss
        coeff = self.heat_transfer_coefficient
        A_wall = pi*diameter*depth
        wall_loss = coeff['wall']*A_wall*(T - self.T_air)
        floor_loss = coeff['floor']*A*(T - self.T_earth)
        ceiling_loss = coeff['ceiling']*A*(T - self.T_air)
        duty = (wall_loss + floor_loss + ceiling_loss)*60*60/1e3
        hx.H = hx_ins0.H + duty
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        
        # concrete usage
        D['Wall concrete'] = self.t_wall*pi*(diameter*_m_to_ft)*(depth*_m_to_ft + self.freeboard)
        D['Slab concrete'] = 2*self.t_slab*A*(_m_to_ft**2)
        
        # excavation
        D['Excavation'] = calculate_excavation_volume(self.L_PB,
                                                           self.W_PB,
                                                           self.D_PB,
                                                           self.excav_slope,
                                                           self.constr_access)
        
        # TODO: is the sludge pump needed for other units?
        # pump
        # TODO: does this add the electricity from the pump (make sure it is far less than the unit_electricity, or make sure unit_electricity does not include pumps)
        sludge_pump = self.sludge_pump
        sludge_pump.simulate()
        D.update(sludge_pump.design_results)
        
        # kW
        self.electricity_requirement = V*self.unit_electricity
        
        # kW
        self.add_power_utility(self.electricity_requirement)

    def _cost(self):
        D, C = self.design_results, self.baseline_purchase_costs
        C['Wall concrete'] = D['Wall concrete']*self.wall_concrete_unit_cost
        C['Slab concrete'] = D['Slab concrete']*self.slab_concrete_unit_cost
        C['Excavation'] = D['Excavation']*self.excavation_unit_cost
        C['Gas collection system'] = (C['Wall concrete'] + C['Slab concrete'])*self.gas_collection_cost_factor
        C['Vertical mixer'] = self.electricity_requirement/self.vertical_mixer_unit_power*self.vertical_mixer_unit_price
    
    @property
    def moisture(self):
        return self.outs[0].imass['H2O']/self.outs[0].F_mass
    
    @property
    def unit_process(self):
        return 'anaerobic_digestion'
    
    @property
    def T_air(self):
        '''[float] temperature of the air, [K].'''
        return self._T_air
    @T_air.setter
    def T_air(self, i):
        self._T_air = i

    @property
    def T_earth(self):
        '''[float] temperature of the air, [K].'''
        return self._T_earth
    @T_earth.setter
    def T_earth(self, i):
        self._T_earth = i
    
    @property
    def t_wall(self):
        '''[float] concrete wall thickness, [ft].'''
        return self._t_wall
    @t_wall.setter
    def t_wall(self, i):
        self._t_wall = i
    
    @property
    def freeboard(self):
        '''[float] freeboard added to the depth of the reactor tank, [ft].'''
        return self._freeboard
    @freeboard.setter
    def freeboard(self, i):
        self._freeboard = i
    
    @property
    def t_slab(self):
        '''
        [float] concrete slab thickness, [ft],
        default to be 2 in thicker than the wall thickness.
        '''
        return self._t_slab or self.t_wall + 2/12
    @t_slab.setter
    def t_slab(self, i):
        self._t_slab = i
    
    @property
    def L_PB(self):
        '''[float] length of the pump building, [ft].'''
        return self._L_PB
    @L_PB.setter
    def L_PB(self, i):
        self._L_PB = i

    @property
    def W_PB(self):
        '''[float] width of the pump building, [ft].'''
        return self._W_PB
    @W_PB.setter
    def W_PB(self, i):
        self._W_PB = i

    @property
    def D_PB(self):
        '''[float] depth of the pump building, [ft].'''
        return self._D_PB
    @D_PB.setter
    def D_PB(self, i):
        self._D_PB = i

    @property
    def excav_slope(self):
        '''[float] slope for excavation (horizontal/vertical).'''
        return self._excav_slope
    @excav_slope.setter
    def excav_slope(self, i):
        self._excav_slope = i

    @property
    def constr_access(self):
        '''[float] extra room for construction access, [ft].'''
        return self._constr_access
    @constr_access.setter
    def constr_access(self, i):
        self._constr_access = i

# =============================================================================
# Dewatering
# =============================================================================

class Dewatering(SanUnit):
    '''
    Dewatering of sludge or biosolids, with the assumption that polymer can be
    ignored in the mass balance.
    
    scope 1 emission: fugitive methane (if from anaerobically digested biosolids)
    scope 2 emission: electricity
    scope 3 emission: polymer
    
    Parameters
    ----------
    ins : iterable
        input_sludge, polymer.
    outs : iterable
        dewatered_solids, reject, methane.
    target_moisture : float
        targeted moisture content, [-].
    polymer_addition : float
        polymer used for conditioning sludge, [kg polymer/dry tonne sludge].
    methane_concentration : float
        methane remains in anaerobically digested biosolids, [mg CH4/L].
    mathane_loss : float
        methane loss ratio during dewater, [-].
    unit_electricity : float
        electricity for dewatering, [kWh/m3].
        from BioSTEAM, 1.40.
    '''
    _N_ins = 2
    _N_outs = 3
    auxiliary_unit_names = ('effluent_pump','sludge_pump')
    solids_loading_range = (2, 40)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', target_moisture=0.8,
                 polymer_addition=5, methane_concentration=6,
                 mathane_loss=0.9, unit_electricity=1.40):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target_moisture = target_moisture
        self.polymer_addition = polymer_addition
        self.methane_concentration = methane_concentration
        self.mathane_loss = mathane_loss
        self.unit_electricity = unit_electricity
        ID = self.ID
        sludge = self.outs[0].proxy(f'{ID}_sludge')
        eff = self.outs[1].proxy(f'{ID}_eff')
        self.sludge_pump = SludgePump(f'.{ID}_sludge_pump', ins=sludge, init_with=init_with)
        self.effluent_pump = SludgePump(f'.{ID}_eff_pump', ins=eff, init_with=init_with)
    
    def _run(self):
        input_sludge, polymer = self.ins
        dewatered_solids, reject, methane = self.outs
        
        dewatered_solids.phase = 'l'
        reject.phase = 'l'
        methane.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (input_sludge.F_mass - input_sludge.imass['H2O'])/1000
        
        polymer.imass['PAM'] = self.dry_solids*self.polymer_addition
        
        # tonne/h
        dewatered_solids_mass_flow = self.dry_solids/(1 - self.target_moisture)
        
        dewatered_solids.copy_like(input_sludge)
       
        dewatered_solids.imass['H2O'] = dewatered_solids_mass_flow*1000*self.target_moisture
       
        reject.imass['H2O'] = input_sludge.F_mass - dewatered_solids.F_mass
        
        if self.ins[0]._source.unit_process == 'anaerobic_digestion':
            methane.imass['CH4'] = self.ins[0].F_vol*1000*self.methane_concentration*self.mathane_loss/1000000
        else:
            methane.imass['CH4'] = 0
    
    def _design(self):
        # totol solids, ton/h
        ts = self.dry_solids/_ton_to_tonne
        
        self.design_results['Solids loading'] = ts
        lb, ub = self.solids_loading_range
        self.design_results['Number of centrifuges'] = ceil(ts/ub)
        cost = 68040*(ts**0.5)
        cost *= bst.CE / 567
        self.baseline_purchase_costs['Centrifuges'] = cost
        self.F_BM['Centrifuges'] = 2.03
        self.design_results['Flow rate'] = self.ins[0].F_vol
        self.power_utility(self.ins[0].F_vol*self.unit_electricity)
    
    def _cost(self):
        # TODO: does this add the electricity from pumps
        for p in (self.effluent_pump, self.sludge_pump): p.simulate()
    
    @property
    def moisture(self):
        return self.target_moisture
    
    @property
    def unit_process(self):
        return 'dewatering'

# =============================================================================
# AlkalineStabilization
# =============================================================================

@cost(ID='Lime stabilization', basis='Dry solids flow',
      units='tonne/day', cost=311386, S=1, CE=CEPCI_by_year[2004], n=0.5623, BM=1)
class AlkalineStabilization(SanUnit):
    '''
    Alkaline stabilization of dewatered sludge, with the assumption
    that no water evaporizes.
    
    scope 1 emission: natural gas combustion
    scope 2 emission: electricity, natural gas upstream
    scope 3 emission: lime
    
    Annualized capital cost (which likely means installed cost) in 2004$
    (30 years, 7% discount rate [r]) from [1]:
        1  MGD (~1 tonne dry solids per day)  $  22,600
        4  MGD (~4 tonne dry solids per day)  $  64,700
        40 MGD (~40 tonne dry solids per day) $ 187,500
    
   The total installed cost can be calculated using the equation below:
        annualized installed cost = installed cost*r/(1 - (1 + r)^(-lifetime))
    
    The estimated total installed costs are:
        1  MGD (~1 tonne dry solids per day)  $   280,444
        4  MGD (~4 tonne dry solids per day)  $   802,865
        40 MGD (~40 tonne dry solids per day) $ 2,326,695
    
    The total installed cost follows a power function:
        total installed cost = 311386*X^0.5623
    
    Parameters
    ----------
    ins : iterable
        dewatered_sludge, lime.
    outs : iterable
        stabilized_solids, NG_as_CO2.
    lime_dose : float
        the dose of lime added to sludge, [tonne lime/dry tonne sludge].
        typically, 0.3 for class A biosolids, 0.2 for class B biosolids.
    # TODO: confirm this is just for combustion
    natural_gas_emission : float
        the CO2 eq emission from natural gas, [kg CO2 eq/dry tonne sludge].
        typically, 15.6 for class A biosolids, 0 for class B biosolids.
    # TODO: check, B>A?
    unit_electricity : float
        electricity for lime stabilization, [kWh/wet tonne sludge].
        typically, 3.7 for class A biosolids, 4.9 for class B biosolids.
    
    References
    ----------
    .. [1] Williford, C.; Chen, W.-Y.; Shamas, N. K.; Wang, L. K. Lime
        Stabilization. In Biosolids Treatment Processes; Wang, L. K.,
        Shammas, N. K., Hung, Y.-T., Eds.; Humana Press: Totowa, NJ,
        2007; pp 207–241. https://doi.org/10.1007/978-1-59259-996-7_7.
    '''
    _N_ins = 2
    _N_outs = 2
    _units = {'Dry solids flow':'tonne/day'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lime_dose=0.3,
                 natural_gas_emission=15.6, unit_electricity=3.7):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.lime_dose = lime_dose
        self.natural_gas_emission = natural_gas_emission
        self.unit_electricity = unit_electricity
    
    def _run(self):
        dewatered_sludge, lime = self.ins
        stabilized_solids, NG_as_CO2 = self.outs
        
        stabilized_solids.phase = 'l'
        
        # dry tonne/h
        self.dry_solids = (dewatered_sludge.F_mass - dewatered_sludge.imass['H2O'])/1000
        
        # TODO: add lime as a component or use CaO or CaCO3
        lime.imass['CaO'] = self.dry_solids*self.lime_dose*1000
        
        stabilized_solids.copy_like(dewatered_sludge)
        
        stabilized_solids.imass['CaO'] = lime.imass['lime']
        
        # TODO: need to add a price for NG_as_CO2
        # TODO: note this is not biogenic
        NG_as_CO2.imass['CO2'] = self.dry_solids*self.natural_gas_emission
        
    def _design(self):
        self.design_results['Dry solids flow'] = self.dry_solids*24
        
        # kW
        self.add_power_utility(self.ins[0].F_mass/1000*self.unit_electricity)
    
    @property
    def moisture(self):
        return self.outs[0].imass['H2O']/self.outs[0].F_mass
    
    @property
    def unit_process(self):
        return 'alkaline_stabilization'

# =============================================================================
# Composting
# =============================================================================

# TODO: BEAM*2024 model has min_solids_content_no_fugitive and C_N_cutoff as parameters for composting but does not use them
# instead, these two parameters are used in land application, considering add them and related algorithms here
# TODO: BEAM*2024 model has TVS as a parameter for sawdust, confirm this is not used
# TODO: the current code assume all bulking agents go to compost without reuse, which may not be true, can add a parameter to consider reuse
class Composting(SanUnit):
    '''
    # TODO: make sure this does not include the transportation of composted sludge from the composting facility (WRRF) to the land application site
    Composting of sludge, considering the land application of composted sludge
    but not considering the transportation of composted sludge from the composting
    facility (WRRF) to the land application site.
    
    scope 1 emission: diesel, fugitive methane, fugitive nitrous oxide, carbon sequestration
    scope 2 emission: electricity
    scope 3 emission: nitrogen, phosphorus
    
    Parameters
    ----------
    ins : iterable
        input_solids, bulking_agent, diesel.
    outs : iterable
        composted_solids, fugitive_methane, fugitive_nitrous_oxide.
    lipid_2_C : float
        lipid to C ratio, [-].
    protein_2_C : float
        protein to C ratio, [-].
    carbo_2_C : float
        carbohydrate to C ratio, [-].
    # TODO: what if there is N loss before composting, is protein_2_N here still accurate?
    protein_2_N, float
        N to protein ratio, [-].
    N_2_P, float
        P to N ratio, [-].
    bulking_agent_solids : float
        volumetric ratio of the bulking agent (sawdust) to solids, [-].
    bulking_agent_density : float
        density of the bulking agent (sawdust), [tonne/m3].
    bulking_agent_solids_ratio : float
        dry matter ratio of the bulking agent (sawdust), [-].
    bulking_agent_OC_ratio : float
        organic carbon ratio of the bulking agent (sawdust) based on dry weight, [-].
    bulking_agent_nitrogen_ratio : float
        nitrogen ratio of the bulking agent (sawdust) based on dry weight, [-].
    bulking_agent_grinding_onsite : bool
        whether the bulking agent (sawdust) grinded onsite, [True, False].
    unit_diesel_grinding : float
        diesel for grinding, [L diesel/wet tonne].
    unit_diesel_other_composting : float
        diesel for setting up and breaking down piles, [L diesel/wet tonne].
        typically, 5 for windrow, 2.5 for aerated static pile (ASP), 2.5 for in vessel.
    # TODO: this is not included in the BEAM*2024 model
    compost_solids_ratio : float
        compost volume to the dewatered solids volume, [-].
    load_size : float
        size of loads per truck, [m3/load].
    load_frequency : float
        load frequency, [load/hr].
    tractor_fuel : float
        tractor fuel, [L diesel/hr].
    # TODO: confirm this ratio is CH4:C, not CH4-C:C, and it just consider C in solids, not the bulking agent
    methane_fugitive_ratio : float
        CH4 to C ratio during composting, [-].
        typically, 0.0001 for well aerated conditions, 0.017 for inadequately aerated conditions.
    # TODO: confirm this ratio is N2O:N, not N2O-N:N, and it just consider N in solids, not the bulking agent
    nitrous_oxide_fugitive_ratio : flaot
        N2O to N ratio during composting, [-].
        typically, 0.00076 for digested solids, 0.018 for undigested solids.
    # TODO: confirm this ratio consider both solids and the bulking agent, not just solids
    unit_carbon_sequestration : ratio
        carbon sequestration potential, [tonne CO2 eq/dry tonne solids].
        typically, 0.15 for low-end, 0.4475 for mid-range, 0.745 for high-end.
    # TODO: double check whether the bulking agent should be included in the unit
    unit_electricity : float
        electricity for composting, [kWh/dry tonne solids (not including the bulking agent)].
        typically, 0 for windrow, 180 for ASP, 291 for in vessel.
    '''
    _N_ins = 3
    _N_outs = 3
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lipid_2_C=0.750,
                 protein_2_C=0.545, carbo_2_C=0.400, protein_2_N=0.159,
                 N_2_P=0.3927, bulking_agent_solids=3, bulking_agent_density=0.25,
                 bulking_agent_solids_ratio=0.61, bulking_agent_OC_ratio=0.518,
                 bulking_agent_nitrogen_ratio=0.0024, bulking_agent_grinding_onsite=True,
                 unit_diesel_grinding=3.3, unit_diesel_other_composting=5,
                 compost_solids_ratio=2, load_size=13, load_frequency=3,
                 tractor_fuel=25, methane_fugitive_ratio=0.0001,
                 nitrous_oxide_fugitive_ratio=0.00076, unit_carbon_sequestration=0.4475,
                 unit_electricity=0):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.lipid_2_C = lipid_2_C
        self.protein_2_C = protein_2_C
        self.carbo_2_C = carbo_2_C
        self.protein_2_N = protein_2_N
        self.N_2_P = N_2_P
        self.bulking_agent_solids = bulking_agent_solids
        self.bulking_agent_density = bulking_agent_density
        self.bulking_agent_solids_ratio = bulking_agent_solids_ratio
        self.bulking_agent_OC_ratio = bulking_agent_OC_ratio
        self.bulking_agent_nitrogen_ratio = bulking_agent_nitrogen_ratio
        self.bulking_agent_grinding_onsite = bulking_agent_grinding_onsite
        self.unit_diesel_grinding = unit_diesel_grinding
        self.unit_diesel_other_composting = unit_diesel_other_composting
        self.compost_solids_ratio = compost_solids_ratio
        self.load_size = load_size
        self.load_frequency = load_frequency
        self.tractor_fuel = tractor_fuel
        self.methane_fugitive_ratio = methane_fugitive_ratio
        self.nitrous_oxide_fugitive_ratio = nitrous_oxide_fugitive_ratio
        self.unit_carbon_sequestration = unit_carbon_sequestration
        self.unit_electricity = unit_electricity
    
    def _run(self):
        input_solids, bulking_agent, diesel = self.ins
        fugitive_methane, fugitive_nitrous_oxide, sequestered_carbon_dioxide = self.outs
        
        fugitive_methane.phase = 'g'
        fugitive_nitrous_oxide.phase = 'g'
        sequestered_carbon_dioxide.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (input_solids.F_mass - input_solids.imass['H2O'])/1000
        
        # TODO: add sawdust or something more general (e.g., bulking agent) or a representative chemical (if not aleady) as a component
        bulking_agent_mass_flow = self.ins[0].F_vol*self.bulking_agent_solids*self.bulking_agent_density*1000
        bulking_agent.imass['Sawdust'] = bulking_agent_mass_flow*self.bulking_agent_solids_ratio
        bulking_agent.imass['H2O'] = bulking_agent_mass_flow*(1 - self.bulking_agent_solids_ratio)
        
        # kg OC/h
        input_solids_OC_mass_flow = input_solids.imass['Sludge_lipid']*self.lipid_2_C +\
                                    input_solids.imass['Sludge_protein']*self.protein_2_C +\
                                    input_solids.imass['Sludge_carbo']*self.carbo_2_C
        
        # kg N/h
        input_solids_N_mass_flow = input_solids.imass['Sludge_protein']*self.protein_2_N
        
        # kg/h
        bulking_agent_OC_mass_flow = bulking_agent.imass['Sawdust']*self.bulking_agent_OC_ratio
        
        # kg/h
        bulking_agent_N_mass_flow = bulking_agent.imass['Sawdust']*self.bulking_agent_nitrogen_ratio
        
        # blended feedstock C:N ratio, -
        self.blended_C_N = (input_solids_OC_mass_flow + bulking_agent_OC_mass_flow)/\
                           (input_solids_N_mass_flow + bulking_agent_N_mass_flow)
        
        # blended feedstock solids content, -
        self.blended_solids_ratio = (self.dry_solids + bulking_agent.imass['Sawdust'])/\
                                    (input_solids.F_mass + bulking_agent.F_mass)    
        
        # L diesel/h
        diesel_grinding = bulking_agent_mass_flow/1000*self.unit_diesel_grinding
        
        # L diesel/h
        diesel_other_composting = (input_solids.F_mass + bulking_agent.F_mass)/1000*self.unit_diesel_other_composting
        
        # TODO: make sure this does not include the transporation from the composting facility (WRRF) to the land application site
        # L diesel/day
        diesel_land_application = (self.ins[0].F_vol*self.compost_solids_ratio/self.load_size/self.load_frequency*self.tractor_fuel)
        
        # TODO: add price and CI, note the unit here is actually L diesel/h, if data are for kg, convert the unit by considering the density of diesel
        diesel.imass['Diesel'] = diesel_grinding + diesel_other_composting + diesel_land_application
        
        fugitive_methane.imass['CH4'] = input_solids_OC_mass_flow*self.methane_fugitive_ratio
        
        fugitive_nitrous_oxide.imass['N2O'] = input_solids_N_mass_flow*self.nitrous_oxide_fugitive_ratio
        
        # TODO: add CI; any price?
        sequestered_carbon_dioxide.imass['CO2'] = -(self.dry_solids*1000 + bulking_agent.imass['Sawdust'])*self.unit_carbon_sequestration
        
        # kW
        self.add_power_utility(self.dry_solids*self.unit_electricity)
    
    # TODO: calculate the price and CI of ash based on the price and CI of nitrogen
    # TODO: may add a discount rate
    @property
    def nitrogen_mass_flow(self):
        # kg N/h
        return self.ins[0].imass['Sludge_protein']*self.protein_2_N
    
    # TODO: calculate the price and CI of ash based on the price and CI of phosphorus
    # TODO: may add a discount rate
    @property
    def phosphorus_mass_flow(self):
        # kg P/h
        return self.ins[0].imass['Sludge_protein']*self.protein_2_N*self.N_2_P
    
    @property
    def unit_process(self):
        return 'landfilling'

# =============================================================================
# HeatDrying
# =============================================================================

class HeatDrying(SanUnit):
    '''
    Heat drying of sludge or biosolids.
    
    scope 1 emission: natural gas combustion
    scope 2 emission: electricity, natural gas upstream
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        input_sludge.
    outs : iterable
        dried_solids, vapor.
    target_moisture : float
        targeted moisture content, [-].
    unit_heat : float
        energy for removing unit water from solids, [GJ/tonne water].
    T_in : float
        inlet process stream temperature, [K].
    unit_electricity : float
        electricity for heat drying, [kWh/dry tonne solids].
    '''
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', target_moisture=0.2,
                 unit_heat=4.5, T_in=650, unit_electricity=214):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target_moisture = target_moisture
        self.unit_heat = unit_heat
        self.T_in = T_in
        self.unit_electricity = unit_electricity
    
    def _run(self):
        input_sludge = self.ins[0]
        dried_solids, vapor = self.outs
        
        dried_solids.phase = 'l'
        vapor.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (input_sludge.F_mass - input_sludge.imass['H2O'])/1000
       
        dried_solids_mass_flow = self.dry_solids/(1 - self.target_moisture)
        
        dried_solids.copy_like(input_sludge)
       
        dried_solids.imass['H2O'] = dried_solids_mass_flow*1000*self.target_moisture
       
        vapor.imass['H2O'] = input_sludge.F_mass - dried_solids.F_mass
        
    def _design(self):
        # TODO: does this work, this might not be natural gas, but need to add GWP
        # TODO: if this does not work, can calculate NG amount manually based on the heating value of NG (and T_in will not be needed)
        self.add_heat_utility(unit_duty=self.outs[1].imass['H2O']/1000*self.unit_heat*1000000,
                              T_in=self.T_in, hxn_ok=True)
        
        # kW
        self.add_power_utility(self.dry_solids*self.unit_electricity)
        
    @property
    def moisture(self):
        return self.target_moisture
    
    @property
    def unit_process(self):
        return 'heat_drying'

# =============================================================================
# Incineration
# =============================================================================

class Incineration(SanUnit):
    '''
    Incineration of dried sludge.
    
    scope 1 emission: fugitive methane, fugitive nitrous oxide
    scope 2 emission: electricity
    scope 3 emission: phosphorus
    
    Parameters
    ----------
    ins : iterable
        dried_solids.
    outs : iterable
        ash, vapor, fugitive_methane, fugitive_nitrous_oxide.
    incineration_temperature : float
        the temperature in the incinerator, [°C].
    ash_moisture_content : float
        moisture content of ash, [-].
    heat_water_removal : float
        energy for removing unit water from solids, [GJ/tonne water].
    # TODO: note this value assumes 20% solids
    fugitive_methane_incineration : float
        fugitive methane during incineration, [tonne CH4/dry tonne solids].
    protein_2_N, float
        N to protein ratio, [-].
    N_2_P, float
        P to N ratio, [-].
    suzuki_constant_1 : float
        for nitrous oxide calculation, [-].
    suzuki_constant_2 : float
        for nitrous oxide calculation, [-].
    # TODO: check the equation name
    suzuki_lowest_temperature : float
        the lowest temperature to apply Suzuki equation, [°C].
    urea_catalyst : bool
        whether the urea catalyst is used during incineration, [True, False]
    N2O_adjustment_factor_urea : float
        N2O adjustment factor due to the use of the urea catalyst, [-].
    heat_solids_incineration : float
        energy from solids incineration, [MJ/dry tonne solids].
        typically, 12000 for digested soldis, 23000 for undigested solids.
    additional_fuel_ratio: float
        additional fuel ratio used in the incinerator, [-].
        typically, 0 for FBI, 0.2 for MHI.
    heat_recovery_ratio : float
        ratio of heat (potentially) recovered from solids incineration, [-].
    heat_recovery_efficiency : float
        efficiency of recovering heat from solids incineration, [-].
    electricity_recovery_ratio : float
        ratio of electricity (potentially) recovered from solids incineration, [-].
        typically, 0.
    BTU_to_kWh_with_efficiency : float
        convert BTU to kWh but considering loss, [kWh/BTU].
    net_capacity_factor : float
        actual electricity generated to the theoretical maximum, [-].
    unit_electricity : float
        electricity for incineration, [kWh/dry tonne solids].
        typically, 200 for FBI, 285 for MHI.
    '''
    _N_ins = 2
    _N_outs = 4
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', incineration_temperature=850,
                 ash_moisture_content=0.01, heat_water_removal=4.5,
                 fugitive_methane_incineration=0.0000097, protein_2_N=0.159,
                 N_2_P=0.3927, suzuki_constant_1=161.3, suzuki_constant_2=0.14,
                 suzuki_lowest_temperature=750, urea_catalyst=True,
                 N2O_adjustment_factor_urea=0.2, heat_solids_incineration=12000,
                 additional_fuel_ratio=0, heat_recovery_ratio=0.5,
                 heat_recovery_efficiency=0.8, electricity_recovery_ratio=0,
                 BTU_to_kWh_with_efficiency=0.0000854, net_capacity_factor=0.85,
                 unit_electricity=200):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.incineration_temperature = incineration_temperature
        self.ash_moisture_content = ash_moisture_content
        self.heat_water_removal = heat_water_removal
        self.fugitive_methane_incineration = fugitive_methane_incineration
        self.protein_2_N = protein_2_N
        self.N_2_P = N_2_P
        self.suzuki_constant_1 = suzuki_constant_1
        self.suzuki_constant_2 = suzuki_constant_2
        self.suzuki_lowest_temperature = suzuki_lowest_temperature
        self.urea_catalyst = urea_catalyst
        self.N2O_adjustment_factor_urea = N2O_adjustment_factor_urea
        self.heat_solids_incineration = heat_solids_incineration
        self.additional_fuel_ratio = additional_fuel_ratio
        self.heat_recovery_ratio = heat_recovery_ratio
        self.heat_recovery_efficiency = heat_recovery_efficiency
        self.electricity_recovery_ratio = electricity_recovery_ratio
        self.BTU_to_kWh_with_efficiency = BTU_to_kWh_with_efficiency
        self.net_capacity_factor = net_capacity_factor
        self.unit_electricity = unit_electricity
    
    def _run(self):
        dried_solids, natural_gas = self.ins
        ash, vapor, fugitive_methane, fugitive_nitrous_oxide = self.outs
        
        ash.phase = 's'
        vapor.phase = 'g'
        fugitive_methane.phase = 'g'
        fugitive_nitrous_oxide.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (dried_solids.F_mass - dried_solids.imass['H2O'])/1000
        
        # TODO: calculate the price and CI of ash based on the price and CI of phosphorus
        # TODO: may add a discount rate
        ash.imass['Sludge_ash'] = dried_solids.imass['Sludge_ash']
        ash.imass['H2O'] = ash.imass['Sludge_ash']/(1 - self.ash_moisture_content)*self.ash_moisture_content
        
        vapor.imass['H2O'] = dried_solids.imass['H2O'] - ash.imass['H2O']
        
        fugitive_methane.imass['CH4'] = self.dry_solids*self.fugitive_methane_incineration*1000
        
        # dry tonne N/day
        N_mass_flow = dried_solids.imass['Sludge_protein']*self.protein_2_N/1000*24
        
        if self.suzuki_constant_1 - (self.suzuki_constant_2*(self.incineration_temperature + _C_to_K))/100 < 0:
            # tonne N2O/day
            N2O_before_adjustment = 0
        else:
            # tonne N2O/day
            N2O_before_adjustment = (N_mass_flow*(self.suzuki_constant_1 - (self.suzuki_constant_2*(max(self.incineration_temperature, self.suzuki_lowest_temperature) + _C_to_K)))/100*_N_to_N2O)
        
        if self.urea_catalyst:
            # tonne N2O/day
            N2O_urea_catalyst = N2O_before_adjustment*self.N2O_adjustment_factor_urea
        else:
            # tonne N2O/day
            N2O_urea_catalyst = 0
        
        if 1 - self.ash_moisture_content < 0.24:
            # - 
            N2O_reduction_ratio = 0
        elif 1 - self.ash_moisture_content < 0.87:
            # - 
            N2O_reduction_ratio = 0.5
        else:
            # - 
            N2O_reduction_ratio = 0.6
        
        # tonne N2O/day
        N2O_reduction = -N2O_before_adjustment*N2O_reduction_ratio
        
        # tonne N2O/day
        N2O_fugitive = N2O_before_adjustment + N2O_urea_catalyst + N2O_reduction
        
        fugitive_nitrous_oxide.imass['N2O'] = N2O_fugitive*1000/24
        
        # MJ/h
        energy_solids_incineration = self.dry_solids*self.heat_solids_incineration
        
        # MJ/h
        energy_water_removal = vapor.imass['H2O']/1000*self.heat_water_removal*1000
        
        # m3 natural gas/h
        natural_gas_generated = energy_solids_incineration*self.heat_recovery_ratio*self.heat_recovery_efficiency/_BTU_to_MJ/natural_gas_heat
        
        # m3 natural gas/h
        natural_gas_requirement = (1 + self.additional_fuel_ratio)*energy_water_removal/_BTU_to_MJ/natural_gas_heat
        
        # m3 natural gas/h
        natural_gas_net = natural_gas_requirement - natural_gas_generated
        
        # TODO: assume excess heat, if any, is lost
        # TODO: need price and CI
        natural_gas.imass['CH4'] = max(0, natural_gas_net*methane_density)
        
        # kW
        electricity_generated = energy_solids_incineration*self.electricity_recovery_ratio/_BTU_to_MJ*self.BTU_to_kWh_with_efficiency*self.net_capacity_factor
        
        # kW
        electricity_requirement = self.dry_solids*self.unit_electricity
        
        # kW
        self.electricity_net = electricity_requirement - electricity_generated
        
    def _design(self):
        # kW
        self.add_power_utility(self.electricity_net)
        
    @property
    def moisture(self):
        return self.ash_moisture_content
    
    # TODO: calculate the price and CI of ash based on the price and CI of phosphorus
    # TODO: may add a discount rate
    @property
    def phosphorus_mass_flow(self):
        # kg P/h
        return self.ins[0].imass['Sludge_protein']*self.protein_2_N*self.N_2_P
    
    @property
    def unit_process(self):
        return 'incineration'

# =============================================================================
# Pyrolysis
# =============================================================================

class Pyrolysis(SanUnit):
    '''
    Pyrolysis of dried sludge, with the assumption that pyrolysis is autothermal.
    
    scope 1 emission: fugitive methane, fugitive nitrous oxide
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        dried_solids.
    outs : iterable
        stabilized_solids, fugitive_methane, fugitive_nitrous_oxide.
    mass_loss : float
        solids mass loss ratio during pyrolysis, [-].
    biochar_moisture_content : float
        moisture content of biochar, [-].
    unit_fugitive_methane : float
        fugitive methane during pyrolysis, [g CH4/dry tonne solids].
    unit_fugitive_nitrous_oxide : float
        fugitive nitrous oxide during pyrolysis, [g N2O/dry tonne solids].
    unit_electricity : float
        electricity for pyrolysis, [kWh/dry tonne solids].
    '''
    _N_ins = 1
    _N_outs = 3
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', mass_loss=0.5,
                 biochar_moisture_content=0.01, unit_fugitive_methane=2.65,
                 unit_fugitive_nitrous_oxide=5.23, unit_electricity=123.424):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.mass_loss = mass_loss
        self.biochar_moisture_content = biochar_moisture_content
    
    def _run(self):
        dried_solids = self.ins[0]
        biochar, fugitive_methane, fugitive_nitrous_oxide = self.outs
        
        biochar.phase = 's'
        fugitive_methane.phase = 'g'
        fugitive_nitrous_oxide.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (dried_solids.F_mass - dried_solids.imass['H2O'])/1000
        
        # TODO: add biochar or a representative chemical (if not already) as a component
        biochar.imass['biochar'] = self.dry_solids*(1 - self.mass_loss)
        biochar.imass['H2O'] = biochar.imass['biochar']/(1 - self.biochar_moisture_content)*self.biochar_moisture_content
        
        fugitive_methane.imass['CH4'] = self.dry_solids*self.unit_fugitive_methane/1000
        
        # TODO: check if N2O is already a component
        fugitive_nitrous_oxide.imass['N2O'] = self.dry_solids*self.unit_fugitive_nitrous_oxide/1000
        
        # kW
        self.add_power_utility(self.dry_solids*self.unit_electricity)
        
    @property
    def moisture(self):
        return self.biochar_moisture_content
    
    @property
    def unit_process(self):
        return 'pyrolysis'

# =============================================================================
# Landfilling
# =============================================================================

# TODO: in BEAM*2024, there is a parameter 'Landfill gas recovery, average U. S.' (value 0.75) not used in the code, double check this
class Landfilling(SanUnit):
    '''
    Landfilling of sludge or biosolids, without accounting for solids remaining
    the landfills in the mass balance.
    
    scope 1 emission: fugitive methane, fugitive nitrous oxide, carbon sequestration
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        input_solids.
    outs : iterable
        sequestered_carbon_dioxide, fugitive_methane, fugitive_nitrous_oxide.
    lipid_2_C : float
        lipid to C ratio, [-].
    protein_2_C : float
        protein to C ratio, [-].
    carbo_2_C : float
        carbohydrate to C ratio, [-].
    # TODO: what if there is N loss before landfilling, is protein_2_N here still accurate?
    protein_2_N, float
        N to protein ratio, [-].
    # TODO: this parameter is used for uncertainty, which can be done in Monte Carlo; remove this parameter
    landfilling_uncertainty : float
        landfilling uncertainty factor for methane emission, [-].
    methane_landfilling_gas : float
        methane content in landfilling gas, [-].
    MCF : float
        methane correction factor, [-].
    OC_decomposition_ratio : float
        organic carbon decomposition ratio, [-].
        typically, 0.5 for completely digested biosolids,
                   0.65 for partially digested biosolids,
                   0.8 for undigested sludge.
    # TODO: check the unit
    k_decay : float
        OC decay rate constant, [-].
        typically, 0.06 for cool & dry environment,
                   0.185 for cool & wet environment,
                   0.085 for warm & dry environment,
                   0.4 for warm & wet environment.
    flare_fugitive_ratio : float
        ratio of fugitive biogas during flaring, [-].
        typically, 0.01 for default, 0 for enclosed, 0.05 for candlestick.
    C_N_cutoff : float
        minimum C to N ratio to produce fugitive nitrous oxide, [-].
    N2O_N_landfilling : float
        N2O as N to the total N ratio in solids during landfilling, [-].
    methane_electricity : float
        percentage of captured methane used to generate electricity, [-].
    BTU_to_kWh_with_efficiency : float
        convert BTU to kWh but considering loss, [kWh/BTU].
    net_capacity_factor : float
        actual electricity generated to the theoretical maximum, [-].
    '''
    _N_ins = 1
    _N_outs = 3
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lipid_2_C=0.750,
                 protein_2_C=0.545, carbo_2_C=0.400, protein_2_N=0.159,
                 methane_landfilling_gas=0.5, MCF=1, OC_decomposition_ratio=0.5,
                 k_decay=0.06, flare_fugitive_ratio=0.01, C_N_cutoff=30,
                 N2O_N_landfilling=0.015, methane_electricity=0.5,
                 BTU_to_kWh_with_efficiency=0.0000854, net_capacity_factor=0.85):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.lipid_2_C = lipid_2_C
        self.protein_2_C = protein_2_C
        self.carbo_2_C = carbo_2_C
        self.protein_2_N = protein_2_N
        self.methane_landfilling_gas = methane_landfilling_gas
        self.MCF = MCF
        self.OC_decomposition_ratio = OC_decomposition_ratio
        self.k_decay = k_decay
        self.flare_fugitive_ratio = flare_fugitive_ratio
        self.C_N_cutoff = C_N_cutoff
        self.N2O_N_landfilling = N2O_N_landfilling
        self.methane_electricity = methane_electricity
        self.BTU_to_kWh_with_efficiency = BTU_to_kWh_with_efficiency
        self.net_capacity_factor = net_capacity_factor
    
    def _run(self):
        input_solids = self.ins[0]
        fugitive_methane, fugitive_nitrous_oxide, sequestered_carbon_dioxide = self.outs
        
        fugitive_methane.phase = 'g'
        fugitive_nitrous_oxide.phase = 'l'
        
        # dry tonne/h
        self.dry_solids = (input_solids.F_mass - input_solids.imass['H2O'])/1000
        
        # kg OC/h
        OC_mass_flow = input_solids.imass['Sludge_lipid']*self.lipid_2_C +\
                       input_solids.imass['Sludge_protein']*self.protein_2_C +\
                       input_solids.imass['Sludge_carbo']*self.carbo_2_C
        
        # decayed OC schedule
        decayed_OC_schedule = lambda x: self.k_decay*(1 - self.k_decay)**x
        
        # CH4 lost schedule
        CH4_lost_schedule = {'0_1': 1,
                             '2_4': 0.5,
                             '5_14': 0.25,
                             'after_capping': 0.1}
        
        # oxidized CH4 lost schedule
        oxidized_CH4_lost_schedule = {'0_1': 0.1,
                                      '2_4': 0.25,
                                      '5_14': 0.25,
                                      'after_capping': 0.35}
        
        # kg CH4/h
        methane_fugitive_years_0_1 = OC_mass_flow*self.OC_decomposition_ratio*\
                                     sum(decayed_OC_schedule(i) for i in range(0, 2))*_C_to_CH4*self.methane_landfilling_gas*\
                                     CH4_lost_schedule['0_1']*self.MCF*(1-oxidized_CH4_lost_schedule['0_1'])
        
        # kg CH4/h
        methane_fugitive_years_2_4 = OC_mass_flow*self.OC_decomposition_ratio*\
                                     sum(decayed_OC_schedule(i) for i in range(2, 5))*_C_to_CH4*self.methane_landfilling_gas*\
                                         CH4_lost_schedule['2_4']*self.MCF*(1-oxidized_CH4_lost_schedule['2_4'])
        
        # kg CH4/h
        methane_fugitive_years_5_14 = OC_mass_flow*self.OC_decomposition_ratio*\
                                      sum(decayed_OC_schedule(i) for i in range(5, 15))*_C_to_CH4*self.methane_landfilling_gas*\
                                          CH4_lost_schedule['5_14']*self.MCF*(1-oxidized_CH4_lost_schedule['5_14'])
        
        # kg CH4/h
        methane_fugitive_after_capping = OC_mass_flow*self.OC_decomposition_ratio*\
                                         sum(decayed_OC_schedule(i) for i in range(15, 30))*_C_to_CH4*self.methane_landfilling_gas*\
                                         CH4_lost_schedule['after_capping']*self.MCF*(1-oxidized_CH4_lost_schedule['after_capping'])
        
        # kg CH4/h
        methane_combustion_years_0_1 = OC_mass_flow*self.OC_decomposition_ratio*\
                                       sum(decayed_OC_schedule(i) for i in range(0, 2))*_C_to_CH4*self.methane_landfilling_gas*\
                                       (1 - CH4_lost_schedule['0_1'])*self.MCF
        
        # kg CH4/h
        methane_combustion_years_2_4 = OC_mass_flow*self.OC_decomposition_ratio*\
                                       sum(decayed_OC_schedule(i) for i in range(2, 5))*_C_to_CH4*self.methane_landfilling_gas*\
                                       (1 - CH4_lost_schedule['2_4'])*self.MCF
        
        # kg CH4/h
        methane_combustion_years_5_14 = OC_mass_flow*self.OC_decomposition_ratio*\
                                        sum(decayed_OC_schedule(i) for i in range(5, 15))*_C_to_CH4*self.methane_landfilling_gas*\
                                        (1 - CH4_lost_schedule['5_14'])*self.MCF
        
        # kg CH4/h
        methane_combustion_after_capping = OC_mass_flow*self.OC_decomposition_ratio*\
                                           sum(decayed_OC_schedule(i) for i in range(15, 30))*_C_to_CH4*self.methane_landfilling_gas*\
                                           (1 - CH4_lost_schedule['after_capping'])*self.MCF
        
        # kg CH4/h
        self.methane_combustion = methane_combustion_years_0_1 + methane_combustion_years_2_4 +\
                                  methane_combustion_years_5_14 + methane_combustion_after_capping
        
        fugitive_methane.imass['CH4'] = methane_fugitive_years_0_1 + methane_fugitive_years_2_4 + methane_fugitive_years_5_14 +\
                                        methane_fugitive_after_capping + self.methane_combustion*self.flare_fugitive_ratio
        
        # kg N/h
        N_mass_flow = input_solids.imass['Sludge_protein']*self.protein_2_N
        
        if OC_mass_flow/N_mass_flow < self.C_N_cutoff:
            fugitive_nitrous_oxide.imass['N2O'] = N_mass_flow*self.N2O_N_landfilling*_N_to_N2O
        else:
            fugitive_nitrous_oxide.imass['N2O'] = 0
        
        # TODO: add CI; any price?
        sequestered_carbon_dioxide.imass['CO2'] = -OC_mass_flow*(1 - self.OC_decomposition_ratio)*_C_to_CO2
        
        
    def _design(self):
        # TODO: confirm no methane is combusted to produce heat
        # kW
        self.add_power_utility(-self.methane_combustion*self.methane_electricity/methane_density*methane_heat*self.BTU_to_kWh_with_efficiency*self.net_capacity_factor)
        
    @property
    def unit_process(self):
        return 'landfilling'

# =============================================================================
# LandApplication
# =============================================================================

# TODO: shorten names of some parameters
class LandApplication(SanUnit):
    '''
    # TODO: make sure this does not include the transportation of biosolids from the WRRF to the land application site
    Land application of biosolids (not including composted sludge), not considering
    the transportation of biosolids from the WRRF to the land application site.
    
    scope 1 emission: diesel, fugitive methane, fugitive nitrous oxide, carbon sequestration, CaCO3
    scope 2 emission: N/A
    scope 3 emission: nitrogen, phosphorus
    
    Parameters
    ----------
    ins : iterable
        biosolids, diesel.
    outs : iterable
        fugitive_methane, fugitive_nitrous_oxide, carbon_dioxide.
    lipid_2_C : float
        lipid to C ratio, [-].
    protein_2_C : float
        protein to C ratio, [-].
    carbo_2_C : float
        carbohydrate to C ratio, [-].
    # TODO: what if there is N loss before landfilling, is protein_2_N here still accurate?
    protein_2_N, float
        N to protein ratio, [-].
    N_2_P, float
        P to N ratio, [-].
    load_size : float
        size of loads per truck, [m3/load].
    load_frequency : float
        load frequency, [load/hr].
    tractor_fuel : float
        tractor fuel, [L diesel/hr].
    min_solids_content_no_fugitive : float
        the solids ratio above which no fugitive emission druing storage, [-].
    storage_time : float
        biosolids storage time, [day].
    unit_fugitive_methane_strogae : float
        fugitive methane during storage before land application, [kg CH4/m3/day]
    C_N_cutoff : float
        minimum C to N ratio to produce fugitive nitrous oxide, [-].
    fine_textured_ratio : float
        ratio of soils with fine textures, [-].
    nitrous_oxide_N_TN_ratio : float
        ratio of N emitted as N2O, [-].
    min_solids_content_nitrous_oxide_reduction : float
        minimum solids content for N2O reduction, [-].
    nitrous_oxide_reduction_ratio : float
        fugitive nitrous oxide reduction ratio during land application when
        biosolids solids content higher than min_solids_content_nitrous_oxide_reduction.
    # TODO: check the unit
    nitrous_oxide_reduction_slope : float
        slope of the function calculating the nitrous oxide reduction when
        biosolids solids content is higher than min_solids_content_no_fugitive
        but less than min_solids_content_nitrous_oxide_reduction, [-].
    # TODO: check the unit
    nitrous_oxide_reduction_intercept : float
        intercept of the function calculating the nitrous oxide reduction when
        biosolids solids content is higher than min_solids_content_no_fugitive
        but less than min_solids_content_nitrous_oxide_reduction, [-].
    unit_fugitive_nitrous_oxide_strogae : float
        fugitive nitrous oxide during storage before land application, [kg N2O/m3/day]
    climate : str
        climate of the land application site, ['humid','arid']
    # TODO: add this here or in the Pyrolysis unit
    OC_sequestered_ratio : float
        organic carbon sequestered in biochar from pyrolysis, [-].
    unit_carbon_sequestration : float
        carbon sequestration potential, [tonne CO2 eq/dry tonne solids].
        typically, 0.15 for low-end, 0.4475 for mid-range, 0.745 for high-end.
    # TODO: update CaCO3 equivalence (%-dry weight) value (this might be related to lime stablization)
    calcium_carbonate_ratio : float
        ratio of calcium carbonate equivalent that can release CO2 during land application based on dry weight, [-].
    # calcium_carbonate_from_waste : bool
        whether calcium carbonate is from waste, [True, False].
    # calcium_carbonate_replacing_lime : bool
        whether calcium carbonate replaces purchased lime during land application, [True, False].
    '''
    _N_ins = 2
    _N_outs = 3
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lipid_2_C=0.750, protein_2_C=0.545,
                 carbo_2_C=0.400, protein_2_N=0.159, N_2_P=0.3927, load_size=13,
                 load_frequency=3, tractor_fuel=25,
                 min_solids_content_no_fugitive=0.55, storage_time=10,
                 unit_fugitive_methane_strogae=0.0091, C_N_cutoff=30,
                 fine_textured_ratio=0.5, nitrous_oxide_N_TN_ratio=0.0275,
                 min_solids_content_nitrous_oxide_reduction=0.8,
                 nitrous_oxide_reduction_ratio=1, 
                 nitrous_oxide_reduction_slope=0.276,
                 nitrous_oxide_reduction_intercept=0.1518,
                 unit_fugitive_nitrous_oxide_strogae=0.00043, climate='humid',
                 OC_sequestered_ratio=0.75, unit_carbon_sequestration=0.4475,
                 # TODO: update CaCO3 equivalence (%-dry weight) value (this might be related to lime stablization)
                 calcium_carbonate_ratio=0.1, calcium_carbonate_from_waste=False,
                 calcium_carbonate_replacing_lime=True):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.lipid_2_C = lipid_2_C
        self.protein_2_C = protein_2_C
        self.carbo_2_C = carbo_2_C
        self.protein_2_N = protein_2_N
        self.N_2_P = N_2_P
        self.load_size = load_size
        self.load_frequency = load_frequency
        self.tractor_fuel = tractor_fuel
        self.min_solids_content_no_fugitive = min_solids_content_no_fugitive
        self.storage_time = storage_time
        self.unit_fugitive_methane_strogae = unit_fugitive_methane_strogae
        self.C_N_cutoff = C_N_cutoff
        self.fine_textured_ratio = fine_textured_ratio
        self.nitrous_oxide_N_TN_ratio = nitrous_oxide_N_TN_ratio
        self.min_solids_content_nitrous_oxide_reduction = min_solids_content_nitrous_oxide_reduction
        self.nitrous_oxide_reduction_ratio = nitrous_oxide_reduction_ratio
        self.nitrous_oxide_reduction_slope = nitrous_oxide_reduction_slope
        self.nitrous_oxide_reduction_intercept = nitrous_oxide_reduction_intercept
        self.unit_fugitive_nitrous_oxide_strogae = unit_fugitive_nitrous_oxide_strogae
        self.climate = climate
        self.OC_sequestered_ratio = OC_sequestered_ratio
        self.unit_carbon_sequestration = unit_carbon_sequestration
        self.calcium_carbonate_ratio = calcium_carbonate_ratio
        self.calcium_carbonate_from_waste = calcium_carbonate_from_waste
        self.calcium_carbonate_replacing_lime = calcium_carbonate_replacing_lime
    
    def _run(self):
        biosolids, diesel = self.ins
        fugitive_methane, fugitive_nitrous_oxide, carbon_dioxide = self.outs
        
        fugitive_methane.phase = 'g'
        fugitive_nitrous_oxide.phase = 'g'
        carbon_dioxide.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (biosolids.F_mass - biosolids.imass['H2O'])/1000
        
        # kg OC/h
        input_solids_OC_mass_flow = biosolids.imass['Sludge_lipid']*self.lipid_2_C +\
                                    biosolids.imass['Sludge_protein']*self.protein_2_C +\
                                    biosolids.imass['Sludge_carbo']*self.carbo_2_C
        
        # kg N/h
        input_solids_N_mass_flow = biosolids.imass['Sludge_protein']*self.protein_2_N
        
        # TODO: make sure this does not include the transporation from the composting facility (WRRF) to the land application site
        # TODO: add price and CI, note the unit here is actually L diesel/h, if data are for kg, convert the unit by considering the density of diesel
        diesel.imass['Diesel'] = (self.ins[0].F_vol/self.load_size/self.load_frequency*self.tractor_fuel)
        
        if self.dry_solids/biosolids.F_mass < self.min_solids_content_no_fugitive:
            fugitive_methane.imass['CH4'] = self.ins[0].F_vol*self.storage_time*24*self.unit_fugitive_methane_strogae/24
            
            # kg N2O/h
            fugitive_nitrous_oxide_storage = self.ins[0].F_vol*self.storage_time*24*self.unit_fugitive_nitrous_oxide_strogae/24
        else:
            fugitive_methane.imass['CH4'] = 0
            
            # kg N2O/h
            fugitive_nitrous_oxide_storage = 0
        
        if input_solids_OC_mass_flow/input_solids_N_mass_flow < self.C_N_cutoff:
            # kg N2O/h
            fugitive_nitrous_oxide_land_application = input_solids_N_mass_flow*self.fine_textured_ratio*self.nitrous_oxide_N_TN_ratio*_N_to_N2O
        else:
            # kg N2O/h
            fugitive_nitrous_oxide_land_application = 0
        
        assert self.min_solids_content_nitrous_oxide_reduction > self.min_solids_content_no_fugitive
        
        if self.dry_solids/biosolids.F_mass >= self.min_solids_content_nitrous_oxide_reduction:
            fugitive_nitrous_oxide_land_application_reduction_ratio = self.nitrous_oxide_reduction_ratio
        elif self.dry_solids/biosolids.F_mass >= self.min_solids_content_no_fugitive:
            # TODO: double check the calculation here and in the BEAM*2024 model
            fugitive_nitrous_oxide_land_application_reduction_ratio = self.nitrous_oxide_reduction_slope*self.dry_solids/biosolids.F_mass - self.nitrous_oxide_reduction_intercept
        else:
            fugitive_nitrous_oxide_land_application_reduction_ratio = 0
        
        # kg N2O/h
        fugitive_nitrous_oxide_land_application_reduction = -fugitive_nitrous_oxide_land_application*fugitive_nitrous_oxide_land_application_reduction_ratio
        
        if self.climate == 'humid':
            fugitive_nitrous_oxide.imass['N2O'] = fugitive_nitrous_oxide_land_application +\
                                                  fugitive_nitrous_oxide_land_application_reduction +\
                                                  fugitive_nitrous_oxide_storage
        elif self.climate == 'arid':
            fugitive_nitrous_oxide.imass['N2O'] = fugitive_nitrous_oxide_storage
        else:
            raise ValueError('climate can only be `humid` or `arid`')
        
        # TODO: consider land application of biochar from pyrolysis here or in the Pyrolysis unit
        # if here, element calculation may not be accurate, and there might be some other problems
        # TODO: check this
        if self.ins[0]._source.unit_process == 'pyrolysis':
            # kg CO2 eq/h
            carbon_dioxide_credit = input_solids_OC_mass_flow*self.OC_sequestered_ratio*_C_to_CO2
        else:
            # kg CO2 eq/h
            carbon_dioxide_credit = self.dry_solids*self.unit_carbon_sequestration*1000
        
        # TODO: double check the relationship here
        if self.calcium_carbonate_from_waste and self.calcium_carbonate_replacing_lime:
            # kg CO2 eq/h
            carbon_dioxide_debit = 0
        else:
            # kg CO2 eq/h
            carbon_dioxide_debit = self.dry_solids*self.calcium_carbonate_ratio*_CaCO3_to_C*_C_to_CO2*1000
        
        carbon_dioxide.imass['CO2'] = carbon_dioxide_debit - carbon_dioxide_credit
        
    # TODO: calculate the price and CI of ash based on the price and CI of nitrogen
    # TODO: may add a discount rate
    @property
    def nitrogen_mass_flow(self):
        # kg N/h
        return self.ins[0].imass['Sludge_protein']*self.protein_2_N
    
    # TODO: calculate the price and CI of ash based on the price and CI of phosphorus
    # TODO: may add a discount rate
    @property
    def phosphorus_mass_flow(self):
        # kg P/h
        return self.ins[0].imass['Sludge_protein']*self.protein_2_N*self.N_2_P
    
    @property
    def unit_process(self):
        return 'land application'