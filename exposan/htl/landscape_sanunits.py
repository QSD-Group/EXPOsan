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

import os, qsdsan as qs, biosteam as bst
from qsdsan import SanUnit
from qsdsan.sanunits import SludgePump, HXutility, HXprocess, Pump
from qsdsan.utils import auom, calculate_excavation_volume
from qsdsan.equipments import Blower, GasPiping
from biosteam import Stream
from biosteam.units.decorators import cost
from math import ceil, pi

GDPCTPI = {2007: 86.352,
           2008: 87.977,
           2009: 88.557,
           2010: 89.619,
           2011: 91.466,
           2012: 93.176,
           2013: 94.786,
           2014: 96.436,
           2015: 97.277,
           2016: 98.208,
           2017: 100.000,
           2018: 102.290,
           2019: 104.008,
           2020: 105.407,
           2021: 110.220,
           2022: 117.995,
           2023: 122.284}

# use methane density for natural gas as well, probably consistent with BioSTEAM, kg/m3
# https://en.wikipedia.org/wiki/Methane (accessed 2025-02-10)
methane_density = 0.657

# diesel density, kg/m3
diesel_density = 850

# methane heat content, BTU/m3
methane_heat = 35830

# natural gas heat content, BTU/m3
natural_gas_heat = 36263

_m3perh_to_MGD = auom('m3/h').conversion_factor('MGD')
_ton_to_tonne = auom('ton').conversion_factor('tonne')
_BTU_to_kWh = auom('BTU').conversion_factor('kWh')
_BTU_to_MJ = auom('BTU').conversion_factor('MJ')
_m_to_ft = auom('m').conversion_factor('ft')
_m_to_in = auom('m').conversion_factor('in')
_psi_to_Pa = auom('psi').conversion_factor('Pa')
_kg_to_lb = auom('kg').conversion_factor('lb')
_m3_to_gal = auom('m3').conversion_factor('gal')
_C_to_K = 273.15
_N_to_N2O = 44/28
_C_to_CH4 = 16/12
_C_to_CO2 = 44/12

folder = os.path.dirname(__file__)

# numbers without citations are likely from the BEAM*2024 model

__all__ = (
    'WRRF',
    'Thickening',
    'AerobicDigestion',
    'AnaerobicDigestion',
    'Dewatering',
    'AlkalineStabilization',
    'Composting',
    'HeatDrying',
    'Landfilling',
    'LandApplication',
    'Incineration',
    'HydrothermalLiquefaction',
    'HydrothermalAlkalineTreatment',
    'Analyzer',
    'CatalyticHydrothermalGasification',
    'Pyrolysis',
    'Gasification',
    'SupercriticalWaterOxidation'
    )

# =============================================================================
# WRRF (water resource recovery facility)
# =============================================================================

class WRRF(SanUnit):
    '''
    WRRF is a fake unit that can set up sludge biochemical compositions
    and calculate sludge elemental compositions.
    
    Parameters
    ----------
    ins : iterable
        ww.
    outs : iterable
        sludge, treated.
    ww_2_dry_sludge : float
        Wastewater-to-dry-sludge conversion factor, [metric ton/day/MGD].
    sludge_moisture : float
        Sludge moisture content.
    sludge_dw_ash : float
        Sludge dry weight ash content.
    sludge_afdw_lipid : float
        Sludge ash free dry weight lipid content.
    sludge_afdw_protein : float
        Sludge ash free dry weight protein content.
    lipid_2_C : float
        Lipid to carbon factor.     
    protein_2_C : float
        Protein to carbon factor.
    carbo_2_C : float
        Carbohydrate to carbon factor.
    C_2_H : float
        Carbon to hydrogen factor.
    protein_2_N : float
        Protein to nitrogen factor.
    N_2_P : float
        Nitrogen to phosphorus factor.
    sludge_wet_density : float
        The density of sludge of 80% moisture content, [kg/m3].
    sludge_distance : float
        Normalized sludge transportation distance, [km].
    wage_adjustment : float
        A coefficient to adjust labor cost.
    PLI : float
        price level index as an approximation to adjust CAPEX, [-].
    
    References
    ----------
    [1] Metcalf and Eddy, Incorporated. 1991. Wastewater Engineering:
        Treatment Disposal and Reuse. New York: McGraw-Hill.
    [2] Cai, L.; Gao, D.; Chen, T.-B.; Liu, H.-T.; Zheng, G.-D.; Yang, Q.-W.
        Moisture Variation Associated with Water Input and Evaporation during
        Sewage Sludge Bio-Drying. Bioresource Technology 2012, 117, 13–19.
        https://doi.org/10.1016/j.biortech.2012.03.092.
    [3] Li, Y.; Leow, S.; Fedders, A. C.; Sharma, B. K.; Guest, J. S.;
        Strathmann, T. J. Quantitative Multiphase Model for Hydrothermal
        Liquefaction of Algal Biomass. Green Chem. 2017, 19 (4), 1163–1174.
        https://doi.org/10.1039/C6GC03294J.
    '''
    _N_ins = 1
    _N_outs = 2
    auxiliary_unit_names = ('sludge_pump',)

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', 
                 ww_2_dry_sludge=1, # [1]
                 sludge_moisture=0.99, sludge_dw_ash=0.231, 
                 sludge_afdw_lipid=0.206, sludge_afdw_protein=0.456,
                 sludge_wet_density=1040, # [2]
                 sludge_distance=100, wage_adjustment=1, PLI=1):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.ww_2_dry_sludge = ww_2_dry_sludge
        self.sludge_moisture = sludge_moisture
        self.sludge_dw_ash = sludge_dw_ash
        self.sludge_afdw_lipid = sludge_afdw_lipid
        self.sludge_afdw_protein = sludge_afdw_protein
        self.sludge_wet_density = sludge_wet_density
        self.sludge_distance = sludge_distance
        self.wage_adjustment = wage_adjustment
        self.PLI = PLI
        ID = self.ID
        sludge = self.outs[0].proxy(f'{ID}_sludge')
        self.sludge_pump = SludgePump(f'.{ID}_sludge_pump', ins=sludge, init_with=init_with)
    
    def _run(self):
        ww = self.ins[0]
        sludge, treated = self.outs

        self.sludge_afdw_carbo = round(1 - self.sludge_afdw_protein - self.sludge_afdw_lipid, 5)   
        
        if self.sludge_dw_ash >= 1:
            raise Exception ('ash can not be larger than or equal to 1')
        
        if self.sludge_afdw_protein + self.sludge_afdw_lipid > 1:
            raise Exception ('protein and lipid exceed 1')
            
        self.sludge_dw = ww.F_vol*_m3perh_to_MGD*self.ww_2_dry_sludge*1000/24

        sludge.imass['H2O'] = self.sludge_dw/(1-self.sludge_moisture)*self.sludge_moisture
        sludge.imass['Sludge_ash'] = self.sludge_dw*self.sludge_dw_ash

        sludge_afdw = self.sludge_dw*(1 - self.sludge_dw_ash)
        sludge.imass['Sludge_lipid'] = sludge_afdw*self.sludge_afdw_lipid
        sludge.imass['Sludge_protein'] = sludge_afdw*self.sludge_afdw_protein
        sludge.imass['Sludge_carbo'] = sludge_afdw*self.sludge_afdw_carbo

        treated.imass['H2O'] = ww.F_mass - sludge.F_mass
    
    def _design(self):
        self.sludge_pump.simulate()
    
    def _cost(self):
        for C in [self.baseline_purchase_costs] + [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]*self.PLI
    
    @property
    def unit_process(self):
        return 'WRRF'
    
    @property
    def sludge_dw_protein(self):
        return self.sludge_afdw_protein*(1-self.sludge_dw_ash)
    
    @property
    def sludge_dw_lipid(self):
        return self.sludge_afdw_lipid*(1-self.sludge_dw_ash)
    
    @property
    def sludge_dw_carbo(self):
        return self.sludge_afdw_carbo*(1-self.sludge_dw_ash)

# =============================================================================
# Thickening
# =============================================================================

class Thickening(SanUnit):
    '''
    Thickening of sludge using either centrifugal or other types of thickener,
    with the assumption that polymer can be ignored in the mass balance.
    
    scope 1 emission: N/A
    scope 2 emission: electricity
    scope 3 emission: polymer (polyacrylamide)
    
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
    PLI : float
        price level index as an approximation to adjust CAPEX, [-].
    '''
    _N_ins = 2
    _N_outs = 2
    auxiliary_unit_names = ('sludge_pump','effluent_pump')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=25, target_moisture=0.97,
                 polymer_addition=5, unit_electricity=4.9, max_capacity=100, PLI=1):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.target_moisture = target_moisture
        self.polymer_addition = polymer_addition
        self.unit_electricity = unit_electricity
        self.max_capacity = max_capacity
        self.PLI = PLI
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
        for p in (self.sludge_pump, self.effluent_pump): p.simulate()
        
        for C in [self.baseline_purchase_costs] + [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]*self.PLI
    
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
    PLI : float
        price level index as an approximation to adjust CAPEX, [-].
    
    References
    ----------
    [1] Metcalf & Eddy, I. Wastewater Engineering : Treatment and Reuse;
        Fifth edition / revised by George Tchobanoglous, Franklin L.
        Burton, H. David Stensel. Boston : McGraw-Hill, 2014.
    [2] Wang, Q.; Yuan, Z. Enhancing Aerobic Digestion of Full-Scale Waste
        Activated Sludge Using Free Nitrous Acid Pre-Treatment. RSC Adv.
        2015, 5 (25), 19128–19134. https://doi.org/10.1039/C4RA17215A.
    '''
    _N_ins = 2
    _N_outs = 2
    
    # K
    _T_air = 17 + _C_to_K
    _T_earth = 10 + _C_to_K
    
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
    
    auxiliary_unit_names = ('sludge_pump',)
    
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
                 init_with='WasteStream', lifetime=30, O2_requirement=2.3,
                 VS_reduction=0.475, HRT=14, SRT=40, unit_electricity=0.03,
                 depth=10, wall_concrete_unit_cost=24,
                 slab_concrete_unit_cost=13, excavation_unit_cost=0.3, PLI=1,
                 F_BM=default_F_BM):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.O2_requirement = O2_requirement
        self.VS_reduction = VS_reduction
        self.HRT = HRT
        self.SRT = SRT
        self.unit_electricity = unit_electricity
        self.depth = depth
        self.wall_concrete_unit_cost = wall_concrete_unit_cost
        self.slab_concrete_unit_cost = slab_concrete_unit_cost
        self.excavation_unit_cost = excavation_unit_cost
        self.PLI=PLI
        self.F_BM.update(F_BM)
        ID = self.ID
        eff = self.outs[0].proxy(f'{ID}_eff')
        self.sludge_pump = SludgePump(f'.{ID}_eff_pump', ins=eff, init_with=init_with)
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
        
        # the blower here does not include electricity
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
        
        self.sludge_pump.simulate()
        self.sludge_pump.power_utility.empty()
        
        for C in [self.baseline_purchase_costs] + [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]*self.PLI
    
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
    
    scope 1 emission: carbon dioxide (natural gas combustion), fugitive methane
    scope 2 emission: electricity, natural gas upstream
    scope 3 emission: N/A
    
    Capital costs include tank, pump, mixer, heat exchanger, and gas collection system.
    
    Parameters
    ----------
    ins : iterable
        input_sludge.
    outs : iterable
        digested_sludge, natural_gas, fugitive_methane.
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
    vertical_mixer_unit_power : float
        the power of one verticle mixer, [kW/unit].
    vertical_mixer_unit_price : float
        the price of one verticle mixer, [$/unit].
    gas_collection_cost_factor : float
        capital cost ratio of the gas collection system and the tank
        excluding execavation, [-].
    PLI : float
        price level index as an approximation to adjust CAPEX, [-].
    '''
    _N_ins = 1
    _N_outs = 4
    
    # K
    _T_air = 17 + _C_to_K
    _T_earth = 10 + _C_to_K
    
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
                 init_with='WasteStream', lifetime=30, VS_reduction=0.425,
                 biogas_yield=0.9, methane_biogas=0.65, biogas_CHP_ratio=0.55,
                 biogas_RNG_ratio=0.35, biogas_flare_ratio=0.09,
                 flare_fugitive_ratio=0.01, RNG_parasitic=0.169,
                 net_capacity_factor=0.85, HRT=22, SRT=22, unit_electricity=0.0065,
                 depth=10, T=35+_C_to_K,
                 heat_transfer_coefficient=dict(wall=0.7, floor=1.7, ceiling=0.95),
                 wall_concrete_unit_cost=24, slab_concrete_unit_cost=13,
                 excavation_unit_cost=0.3,
                 # refer to qsdsan/equipments/_vertical_mixer.py
                 vertical_mixer_unit_power=3.7, vertical_mixer_unit_price=10200,
                 # !!! 0.2 is an assumption, subject to wide uncertainties
                 gas_collection_cost_factor=0.2, PLI=1, F_BM=default_F_BM):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
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
        self.vertical_mixer_unit_power = vertical_mixer_unit_power
        self.vertical_mixer_unit_price = vertical_mixer_unit_price
        self.gas_collection_cost_factor = gas_collection_cost_factor
        self.PLI = PLI
        self.F_BM.update(F_BM)
        ID = self.ID
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'{ID}_hx', ins=hx_in, outs=hx_out)
        eff = self.outs[0].proxy(f'{ID}_eff')
        self.sludge_pump = SludgePump(f'.{ID}_eff_pump', ins=eff, init_with=init_with)
    
    def _run(self):
        input_sludge = self.ins[0]
        digested_sludge, natural_gas, fugitive_methane, carbon_dioxide = self.outs
        
        digested_sludge.phase = 'l'
        natural_gas.phase = 'g'
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
        fugitive_methane.ivol['CH4'] = methane_volume*self.biogas_fugitive_ratio
        
        # from flaring
        fugitive_methane.ivol['CH4'] += methane_volume*self.biogas_flare_ratio*self.flare_fugitive_ratio
        
        # m3/h
        natural_gas_generated = methane_volume*self.biogas_CHP_ratio +\
                                methane_volume*self.biogas_RNG_ratio*(1 - self.RNG_parasitic)
        
        natural_gas.ivol['CH4'] = natural_gas_generated
        
        carbon_dioxide.ivol['CO2'] = biogas_volume*(1 - self.methane_biogas)
    
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
        
        # kW
        self.electricity_requirement = V*self.unit_electricity
        
        # kW
        self.add_power_utility(self.electricity_requirement)

    def _cost(self):
        D, C = self.design_results, self.baseline_purchase_costs
        C['Wall concrete'] = D['Wall concrete']*self.wall_concrete_unit_cost
        C['Slab concrete'] = D['Slab concrete']*self.slab_concrete_unit_cost
        C['Excavation'] = D['Excavation']*self.excavation_unit_cost
        C['Vertical mixer'] = self.electricity_requirement/self.vertical_mixer_unit_power*self.vertical_mixer_unit_price
        # this just adds the cost, not indicating the gas collection system is made of concreate (do not affect LCA since construction is excluded)
        C['Gas collection system'] = (C['Wall concrete'] + C['Slab concrete'])*self.gas_collection_cost_factor
        
        self.sludge_pump.simulate()
        self.sludge_pump.power_utility.empty()
        
        for C in [self.baseline_purchase_costs] + [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]*self.PLI
    
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
    scope 3 emission: polymer (polyacrylamide)
    
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
    PLI : float
        price level index as an approximation to adjust CAPEX, [-].
    '''
    _N_ins = 2
    _N_outs = 3
    auxiliary_unit_names = ('effluent_pump',)
    solids_loading_range = (2, 40)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20, target_moisture=0.8,
                 polymer_addition=5, methane_concentration=6, mathane_loss=0.9,
                 unit_electricity=1.40, PLI=1):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.target_moisture = target_moisture
        self.polymer_addition = polymer_addition
        self.methane_concentration = methane_concentration
        self.mathane_loss = mathane_loss
        self.unit_electricity = unit_electricity
        self.PLI = PLI
        ID = self.ID
        eff = self.outs[1].proxy(f'{ID}_eff')
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
        self.effluent_pump.simulate()
        
        # self.baseline_purchase_costs has already been adjusted through bst.CE in the _design function
        for C in [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]*self.PLI
    
    @property
    def moisture(self):
        return self.target_moisture
    
    @property
    def unit_process(self):
        return 'dewatering'

# =============================================================================
# AlkalineStabilization
# =============================================================================

@cost(ID='Lime stabilization', basis='Dry solids flow', units='tonne/day',
      cost=311386, S=1, CE=qs.CEPCI_by_year[2004], n=0.5623, BM=1)
class AlkalineStabilization(SanUnit):
    '''
    Alkaline stabilization of dewatered sludge, with the assumption
    that no water evaporizes.
    The BEAM model assumes the use of fuel for class A bisolids, which is removed
    since it is atypical to use any fuel in lime stabilization.
    
    scope 1 emission: CO2 (natural gas combustion)
    scope 2 emission: electricity, natural gas upstream
    scope 3 emission: lime
    
    annualized capital cost (which likely means installed cost) in 2004$
    (30 years, 7% discount rate [r]) from [1]:
        1  MGD (~1 tonne dry solids per day)  $  22,600
        4  MGD (~4 tonne dry solids per day)  $  64,700
        40 MGD (~40 tonne dry solids per day) $ 187,500
    
   the total installed cost can be calculated using the equation below:
        annualized installed cost = installed cost*r/(1 - (1 + r)^(-lifetime))
    
    the estimated total installed costs are:
        1  MGD (~1 tonne dry solids per day)  $   280,444
        4  MGD (~4 tonne dry solids per day)  $   802,865
        40 MGD (~40 tonne dry solids per day) $ 2,326,695
    
    the total installed cost follows a power function:
        total installed cost = 311386*X^0.5623
    
    Parameters
    ----------
    ins : iterable
        dewatered_sludge, lime.
    outs : iterable
        stabilized_solids.
    lime_dose : float
        the dose of lime added to sludge, [tonne lime/dry tonne sludge].
        typically, 0.3 for class A biosolids, 0.2 for class B biosolids.
    unit_electricity : float
        electricity for lime stabilization, [kWh/wet tonne sludge].
        typically, 3.7 for class A biosolids, 4.9 for class B biosolids.
        B > A though, just use the larger value
    
    References
    ----------
    [1] Williford, C.; Chen, W.-Y.; Shamas, N. K.; Wang, L. K. Lime
        Stabilization. In Biosolids Treatment Processes; Wang, L. K.,
        Shammas, N. K., Hung, Y.-T., Eds.; Humana Press: Totowa, NJ,
        2007; pp 207–241. https://doi.org/10.1007/978-1-59259-996-7_7.
    '''
    _N_ins = 2
    _N_outs = 1
    
    _units = {'Dry solids flow':'tonne/day'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=25, lime_dose=0.3,
                 unit_electricity=4.9):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.lime_dose = lime_dose
        self.unit_electricity = unit_electricity
    
    def _run(self):
        dewatered_sludge, lime = self.ins
        stabilized_solids = self.outs[0]
        
        stabilized_solids.phase = 'l'
        
        # dry tonne/h
        self.dry_solids = (dewatered_sludge.F_mass - dewatered_sludge.imass['H2O'])/1000
        
        lime.imass['CaO'] = self.dry_solids*1000*self.lime_dose
        
        stabilized_solids.copy_like(dewatered_sludge)
        
        stabilized_solids.imass['CaO'] = lime.imass['lime']
        
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

class Composting(SanUnit):
    '''
    Composting of sludge. Diesel is for composting and land application of compost,
    not including the transportation of compost from the composting facility (WRRF)
    to the land application site. The equipment cost is based on [1].
    
    The current model assumes no reuse of compost (bulking agents) during the composting process.
    
    Fugitive emissions, carbon sequestration, electricity-related emissions,
    and fertilizers-related emissions are based on wastewater residual solids only (based on the BEAM model).
    
    scope 1 emission: CO2 (diesel combustion), fugitive methane, fugitive nitrous oxide
    scope 2 emission: electricity
    scope 3 emission: diesel, bulking agent (sawdust)
    
    Parameters
    ----------
    ins : iterable
        input_solids, bulking_agent, diesel.
    outs : iterable
        compost_cost, fugitive_methane, fugitive_nitrous_oxide, sequestered_carbon_dioxide.
    lipid_2_C : float
        lipid to C ratio, [-].
    protein_2_C : float
        protein to C ratio, [-].
    carbo_2_C : float
        carbohydrate to C ratio, [-].
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
    compost_moisture_content : float
        the moisture content of produced composts, [-].
    load_size : float
        size of loads per truck, [m3/load].
    load_frequency : float
        load frequency, [load/h].
    tractor_fuel : float
        tractor fuel, [L diesel/h].
    methane_fugitive_ratio : float
        CH4 to C ratio of wastewater residual solids, [-].
        typically, 0.0001 for well aerated conditions, 0.017 for inadequately aerated conditions.
    feedstock_digested : bool
        whether the feedstock is digested, [True, False].
    unit_carbon_sequestration : ratio
        carbon sequestration potential, [tonne CO2 eq/dry tonne solids].
        typically, 0.15 for low-end, 0.4475 for mid-range, 0.745 for high-end.
    unit_electricity : float
        electricity for composting, [kWh/dry tonne solids].
        typically, 0 for windrow, 180 for ASP, 291 for in vessel.
    solids_distance : float
        distance between WRRFs and land application sites, [km].
    PLI : float
        price level index as an approximation to adjust CAPEX, [-].
    
    References
    ----------
    [1] Stramer, Y.; Brenner, A.; Cohen, S. B.; Oron, G. Selection of a
        Multi-Stage System for Biosolids Management Applying Genetic Algorithm.
        Environ. Sci. Technol. 2010, 44 (14), 5503–5508.
        https://doi.org/10.1021/es902981t.
    '''
    _N_ins = 3
    _N_outs = 4
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=25, lipid_2_C=0.750,
                 protein_2_C=0.545, carbo_2_C=0.400, protein_2_N=0.159,
                 N_2_P=0.3927, bulking_agent_solids=3, bulking_agent_density=0.25,
                 bulking_agent_solids_ratio=0.61, bulking_agent_OC_ratio=0.518,
                 bulking_agent_nitrogen_ratio=0.0024, bulking_agent_grinding_onsite=True,
                 unit_diesel_grinding=3.3, unit_diesel_other_composting=5,
                 compost_moisture_content=0.5, load_size=13, load_frequency=3,
                 tractor_fuel=25, methane_fugitive_ratio=0.0001,
                 feedstock_digested=True, unit_carbon_sequestration=0.4475,
                 unit_electricity=0, solids_distance=100, PLI=1):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
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
        self.compost_moisture_content = compost_moisture_content
        self.load_size = load_size
        self.load_frequency = load_frequency
        self.tractor_fuel = tractor_fuel
        self.methane_fugitive_ratio = methane_fugitive_ratio
        self.feedstock_digested = feedstock_digested
        self.unit_carbon_sequestration = unit_carbon_sequestration
        self.unit_electricity = unit_electricity
        self.solids_distance = solids_distance
        self.PLI = PLI
    
    def _run(self):
        input_solids, bulking_agent, diesel = self.ins
        compost_cost, fugitive_methane, fugitive_nitrous_oxide, sequestered_carbon_dioxide = self.outs
        
        compost_cost.phase = input_solids.phase
        fugitive_methane.phase = 'g'
        fugitive_nitrous_oxide.phase = 'g'
        sequestered_carbon_dioxide.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (input_solids.F_mass - input_solids.imass['H2O'])/1000
        
        bulking_agent_mass_flow = self.ins[0].F_vol*self.bulking_agent_solids*self.bulking_agent_density*1000
        bulking_agent.imass['Sawdust'] = bulking_agent_mass_flow*self.bulking_agent_solids_ratio
        bulking_agent.imass['H2O'] = bulking_agent_mass_flow*(1 - self.bulking_agent_solids_ratio)
        
        # kg OC/h
        input_solids_OC_mass_flow = input_solids.imass['Sludge_lipid']*self.lipid_2_C +\
                                    input_solids.imass['Sludge_protein']*self.protein_2_C +\
                                    input_solids.imass['Sludge_carbo']*self.carbo_2_C
        
        compost_cost.imass['Compost'] = input_solids.F_mass - input_solids.imass['H2O'] + bulking_agent.F_mass - bulking_agent.imass['H2O']
        compost_cost.imass['H2O'] = compost_cost.imass['Compost']/(1 - self.compost_moisture_content)*self.compost_moisture_content
        
        # kg N/h
        input_solids_N_mass_flow = input_solids.imass['Sludge_protein']*self.protein_2_N
        
        # kg OC/h
        bulking_agent_OC_mass_flow = bulking_agent.imass['Sawdust']*self.bulking_agent_OC_ratio
        
        # kg N/h
        bulking_agent_N_mass_flow = bulking_agent.imass['Sawdust']*self.bulking_agent_nitrogen_ratio
        
        # blended feedstock C:N ratio, -
        self.blended_C_N = (input_solids_OC_mass_flow + bulking_agent_OC_mass_flow)/\
                           (input_solids_N_mass_flow + bulking_agent_N_mass_flow)
        
        # blended feedstock solids content, -
        self.blended_solids_ratio = (self.dry_solids*1000 + bulking_agent.imass['Sawdust'])/\
                                    (input_solids.F_mass + bulking_agent.F_mass)    
        
        # L diesel/h
        diesel_grinding = bulking_agent_mass_flow/1000*self.unit_diesel_grinding
        
        # L diesel/h
        diesel_other_composting = (input_solids.F_mass + bulking_agent.F_mass)/1000*self.unit_diesel_other_composting
        
        # this just includes the land application onsite spreading and does not include the transporation from composting facilities (WRRFs) to land application sites
        # L diesel/h
        diesel_land_application = self.ins[0].F_vol/self.load_size/self.load_frequency*self.tractor_fuel
        
        diesel.imass['Diesel'] = (diesel_grinding + diesel_other_composting + diesel_land_application)/1000*diesel_density
        
        fugitive_methane.imass['CH4'] = input_solids_OC_mass_flow*self.methane_fugitive_ratio
        
        if self.feedstock_digested:
            self.nitrous_oxide_fugitive_ratio = 0.00076
        else:
            self.nitrous_oxide_fugitive_ratio = 0.018
        
        fugitive_nitrous_oxide.imass['N2O'] = input_solids_N_mass_flow*self.nitrous_oxide_fugitive_ratio
        
        sequestered_carbon_dioxide.imass['CO2'] = -self.dry_solids*1000*self.unit_carbon_sequestration
        
        # kW
        self.add_power_utility(self.dry_solids*self.unit_electricity)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        # based on [1], this is likely to be the installed cost
        # based on [1], 1593.7 is marshall and swift equipment cost index 2017, 751 is the baseline marshall and swift equipment cost index, assuming this index has the same trend as CEPCI
        C['Composting equipment'] = (1560*self.ins[0].F_mass/1000*24 + 450000)*1593.7/qs.CEPCI_by_year[2017]*qs.CEPCI_by_year[2022]/751
        
        for C in [self.baseline_purchase_costs] + [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]*self.PLI
    
    @property
    def nitrogen_mass_flow(self):
        # kg N/h
        return self.ins[0].imass['Sludge_protein']*self.protein_2_N
    
    @property
    def phosphorus_mass_flow(self):
        # kg P/h
        return self.ins[0].imass['Sludge_protein']*self.protein_2_N*self.N_2_P
    
    @property
    def unit_process(self):
        return 'composting'

# =============================================================================
# HeatDrying
# =============================================================================

# n: 0.7, assumed
# BM: 3.17, from biosteam/units/heat_exchange.py
@cost(ID='Dryer 1', basis='Half dry solids flow', units='tonne/day',
      cost=220000*80/2/3.17, S=80, CE=qs.CEPCI_by_year[2018], n=0.7, BM=3.17)
@cost(ID='Dryer 2', basis='Half dry solids flow', units='tonne/day',
      cost=220000*80/2/3.17, S=80, CE=qs.CEPCI_by_year[2018], n=0.7, BM=3.17)
class HeatDrying(SanUnit):
    '''
    Heat drying of sludge or biosolids.
    
    scope 1 emission: CO2 (natural gas combustion)
    scope 2 emission: electricity, natural gas upstream
    scope 3 emission: N/A
    
    from [1]: 220000 2018$ CAPEX/dry tonne solids/day for a 80 dry tonne/day.
    # assume CAPEX / installed cost = 2 and installed cost / purchase cost (BM) = 3.17
    
    Parameters
    ----------
    ins : iterable
        input_sludge.
    outs : iterable
        dried_solids, vapor.
    target_moisture : float
        targeted moisture content, [-].
    T_out : float
        outlet solids temperature, [K].
    unit_heat : float
        energy for removing unit water from solids, [GJ/tonne water].
    unit_electricity : float
        electricity for heat drying, [kWh/dry tonne solids].
    natural_gas_HHV : float
        higher heating value of natural gas, [MJ/m3].
    
    References
    ----------
    [1] Hao, X.; Chen, Q.; van Loosdrecht, M. C. M.; Li, J.; Jiang, H.
        Sustainable Disposal of Excess Sludge: Incineration without
        Anaerobic Digestion. Water Research 2020, 170, 115298.
        https://doi.org/10.1016/j.watres.2019.115298.
    '''
    _N_ins = 2
    _N_outs = 2
    
    _units = {'Half dry solids flow':'tonne/day'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=25, target_moisture=0.2,
                 T_out=90 + _C_to_K, unit_heat=4.5, natural_gas_HHV=39,
                 unit_electricity=214):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.target_moisture = target_moisture
        self.T_out = T_out
        self.unit_heat = unit_heat
        self.unit_electricity = unit_electricity
        self.natural_gas_HHV = natural_gas_HHV
    
    def _run(self):
        input_sludge, natural_gas = self.ins
        dried_solids, vapor = self.outs
        
        natural_gas.phase = 'g'
        dried_solids.phase = 'l'
        vapor.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (input_sludge.F_mass - input_sludge.imass['H2O'])/1000
        
        # tonne/h
        dried_solids_mass_flow = self.dry_solids/(1 - self.target_moisture)
        
        dried_solids.copy_like(input_sludge)
        
        dried_solids.imass['H2O'] = dried_solids_mass_flow*1000*self.target_moisture
        
        vapor.imass['H2O'] = input_sludge.F_mass - dried_solids.F_mass
        
        # use natural gas for heat drying base on the BEAM model
        natural_gas.ivol['CH4'] = vapor.imass['H2O']/1000*self.unit_heat*1000/self.natural_gas_HHV
        
        dried_solids.T = vapor.T = self.T_out
    
    def _design(self):
        self.design_results['Half dry solids flow'] = self.dry_solids*24/2
        
        # kW
        self.add_power_utility(self.dry_solids*self.unit_electricity)
    
    @property
    def moisture(self):
        return self.target_moisture
    
    @property
    def unit_process(self):
        return 'heat_drying'

# =============================================================================
# Landfilling
# =============================================================================

class Landfilling(SanUnit):
    '''
    Landfilling of sludge or biosolids, without accounting for solids remaining
    in landfills in the mass balance.
    
    scope 1 emission: fugitive methane, fugitive nitrous oxide
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        input_solids.
    outs : iterable
        landfilled_solids, fugitive_methane, fugitive_nitrous_oxide, sequestered_carbon_dioxide.
    lipid_2_C : float
        lipid to C ratio, [-].
    protein_2_C : float
        protein to C ratio, [-].
    carbo_2_C : float
        carbohydrate to C ratio, [-].
    protein_2_N, float
        N to protein ratio, [-].
    methane_landfilling_gas : float
        methane content in landfilling gas, [-].
    MCF : float
        methane correction factor, [-].
    OC_decomposition_ratio : float
        organic carbon decomposition ratio, [-].
        typically, 0.5 for completely digested biosolids,
                   0.65 for partially digested biosolids,
                   0.8 for undigested sludge.
    k_decay : float
        OC decay rate constant, [year-1].
        typically, 0.06 for cool & dry environment,
                   0.185 for cool & wet environment,
                   0.085 for warm & dry environment,
                   0.4 for warm & wet environment.
        use the average (0.18) if no environment information is considered.
    flare_fugitive_ratio : float
        ratio of fugitive biogas during flaring, [-].
        typically, 0.01 for default, 0 for enclosed, 0.05 for candlestick.
    model_correction_factor : float
        model correction factor for methane emission, [-].
    C_N_cutoff : float
        minimum C to N ratio to produce fugitive nitrous oxide, [-].
    N2O_N_landfilling : float
        N2O as N to the total N ratio in solids during landfilling, [-].
    methane_electricity : float
        ratio of captured methane used to generate electricity, [-].
    BTU_to_kWh_with_efficiency : float
        convert BTU to kWh but considering loss, [kWh/BTU].
    net_capacity_factor : float
        actual electricity generated to the theoretical maximum, [-].
    solids_distance : float
        distance between WRRFs and landfills, [km].
    '''
    _N_ins = 1
    _N_outs = 4
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lipid_2_C=0.750,
                 protein_2_C=0.545, carbo_2_C=0.400, protein_2_N=0.159,
                 methane_landfilling_gas=0.5, MCF=1, OC_decomposition_ratio=0.5,
                 k_decay=0.06, flare_fugitive_ratio=0.01,
                 model_correction_factor=0.75, C_N_cutoff=30,
                 N2O_N_landfilling=0.015,
                 # assume all captured methane is used to generate electricity
                 methane_electricity=1,
                 BTU_to_kWh_with_efficiency=0.0000854, net_capacity_factor=0.85,
                 solids_distance=100):
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
        self.model_correction_factor = model_correction_factor
        self.C_N_cutoff = C_N_cutoff
        self.N2O_N_landfilling = N2O_N_landfilling
        self.methane_electricity = methane_electricity
        self.BTU_to_kWh_with_efficiency = BTU_to_kWh_with_efficiency
        self.net_capacity_factor = net_capacity_factor
        self.solids_distance = solids_distance
    
    def _run(self):
        input_solids = self.ins[0]
        landfilled_solids, fugitive_methane, fugitive_nitrous_oxide, sequestered_carbon_dioxide = self.outs
        
        fugitive_methane.phase = 'g'
        fugitive_nitrous_oxide.phase = 'g'
        sequestered_carbon_dioxide.phase = 'g'
        
        landfilled_solids.copy_like(input_solids)
        
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
        self.methane_fugitive = (methane_fugitive_years_0_1 + methane_fugitive_years_2_4 + methane_fugitive_years_5_14 +\
                                methane_fugitive_after_capping)*self.model_correction_factor
        
        # kg CH4/h
        self.methane_combustion = (methane_combustion_years_0_1 + methane_combustion_years_2_4 +\
                                  methane_combustion_years_5_14 + methane_combustion_after_capping)*self.model_correction_factor
        
        fugitive_methane.imass['CH4'] = self.methane_fugitive + self.methane_combustion*self.flare_fugitive_ratio
        
        # kg N/h
        N_mass_flow = input_solids.imass['Sludge_protein']*self.protein_2_N
        
        if OC_mass_flow/N_mass_flow < self.C_N_cutoff:
            fugitive_nitrous_oxide.imass['N2O'] = N_mass_flow*self.N2O_N_landfilling*_N_to_N2O
        else:
            fugitive_nitrous_oxide.imass['N2O'] = 0
        
        sequestered_carbon_dioxide.imass['CO2'] = -OC_mass_flow*(1 - self.OC_decomposition_ratio)*_C_to_CO2
    
    # add electricity production as a property since it is just included in the LCA (boundary: life cycle) but not TEA (boundary: WRRF)
    @property
    def electricity_kW(self):
        return self.methane_combustion*(1 - self.flare_fugitive_ratio)*self.methane_electricity/methane_density*methane_heat*self.BTU_to_kWh_with_efficiency*self.net_capacity_factor
    
    @property
    def unit_process(self):
        return 'landfilling'

# =============================================================================
# LandApplication
# =============================================================================

class LandApplication(SanUnit):
    '''
    Land application of biosolids (not including compost). Diesel is for land
    application of biosolids, not including the transportation of biosolids
    from the WRRF to the land application site.
    
    Carbon debit due to the lime in the wastewater residual solids is not
    included (CO2 emission during limestone decompoistion is considered in the
    chemical 'lime' in AlkalineStabilization and any CO2 combined with lime
    during other processes should be biogenic or directly from air), while the
    BEAM model does (though the contribution is limited).
    
    scope 1 emission: CO2 (diesel combustion), fugitive methane, fugitive nitrous oxide
    scope 2 emission: N/A
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        biosolids, diesel.
    outs : iterable
        biosolids_cost, fugitive_methane, fugitive_nitrous_oxide, carbon_dioxide.
    lipid_2_C : float
        lipid to C ratio, [-].
    protein_2_C : float
        protein to C ratio, [-].
    carbo_2_C : float
        carbohydrate to C ratio, [-].
    protein_2_N, float
        N to protein ratio, [-].
    N_2_P, float
        P to N ratio, [-].
    load_size : float
        size of loads per truck, [m3/load].
    load_frequency : float
        load frequency, [load/h].
    tractor_fuel : float
        tractor fuel, [L diesel/h].
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
        ratio of N emitted as N2O-N, [-].
    min_solids_content_nitrous_oxide_reduction : float
        minimum solids content for N2O reduction, [-].
    nitrous_oxide_reduction_ratio : float
        fugitive nitrous oxide reduction ratio during land application when
        biosolids solids content higher than min_solids_content_nitrous_oxide_reduction.
    nitrous_oxide_reduction_slope : float
        slope of the function calculating the nitrous oxide reduction when
        biosolids solids content is higher than min_solids_content_no_fugitive
        but less than min_solids_content_nitrous_oxide_reduction, [-].
    nitrous_oxide_reduction_intercept : float
        intercept of the function calculating the nitrous oxide reduction when
        biosolids solids content is higher than min_solids_content_no_fugitive
        but less than min_solids_content_nitrous_oxide_reduction, [-].
    unit_fugitive_nitrous_oxide_strogae : float
        fugitive nitrous oxide during storage before land application, [kg N2O/m3/day]
    climate : str
        climate of the land application site, ['humid','arid']
    unit_carbon_sequestration : float
        carbon sequestration potential, [tonne CO2 eq/dry tonne solids].
        typically, 0.15 for low-end, 0.4475 for mid-range, 0.745 for high-end.
    solids_distance : float
        distance between WRRFs and land application sites, [km].
    '''
    _N_ins = 2
    _N_outs = 4
    
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
                 nitrous_oxide_reduction_intercept=-0.1518,
                 unit_fugitive_nitrous_oxide_strogae=0.00043, climate='humid',
                 unit_carbon_sequestration=0.4475, solids_distance=100):
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
        self.unit_carbon_sequestration = unit_carbon_sequestration
        self.solids_distance = solids_distance
    
    def _run(self):
        biosolids, diesel = self.ins
        biosolids_cost, fugitive_methane, fugitive_nitrous_oxide, carbon_dioxide = self.outs
        
        biosolids_cost.phase = biosolids.phase
        fugitive_methane.phase = 'g'
        fugitive_nitrous_oxide.phase = 'g'
        carbon_dioxide.phase = 'g'
        
        biosolids_cost.copy_like(biosolids)
        
        # dry tonne/h
        self.dry_solids = (biosolids.F_mass - biosolids.imass['H2O'])/1000
        
        # kg OC/h
        input_solids_OC_mass_flow = biosolids.imass['Sludge_lipid']*self.lipid_2_C +\
                                    biosolids.imass['Sludge_protein']*self.protein_2_C +\
                                    biosolids.imass['Sludge_carbo']*self.carbo_2_C
        
        # kg N/h
        input_solids_N_mass_flow = biosolids.imass['Sludge_protein']*self.protein_2_N
        
        # this just includes the land application onsite spreading and does not include the transporation from WRRFs to land application sites
        diesel.imass['Diesel'] = (self.ins[0].F_vol/self.load_size/self.load_frequency*self.tractor_fuel)/1000*diesel_density
        
        if self.dry_solids*1000/biosolids.F_mass < self.min_solids_content_no_fugitive:
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
        
        if self.dry_solids*1000/biosolids.F_mass >= self.min_solids_content_nitrous_oxide_reduction:
            fugitive_nitrous_oxide_land_application_reduction_ratio = self.nitrous_oxide_reduction_ratio
        elif self.dry_solids*1000/biosolids.F_mass >= self.min_solids_content_no_fugitive:
            fugitive_nitrous_oxide_land_application_reduction_ratio = self.nitrous_oxide_reduction_slope*self.dry_solids*1000/biosolids.F_mass + self.nitrous_oxide_reduction_intercept
        else:
            fugitive_nitrous_oxide_land_application_reduction_ratio = 0
        
        # kg N2O/h
        fugitive_nitrous_oxide_land_application_reduction = -fugitive_nitrous_oxide_land_application*fugitive_nitrous_oxide_land_application_reduction_ratio
        
        if self.climate == 'humid':
            fugitive_nitrous_oxide.imass['N2O'] = fugitive_nitrous_oxide_storage +\
                                                  fugitive_nitrous_oxide_land_application +\
                                                  fugitive_nitrous_oxide_land_application_reduction
        elif self.climate == 'arid':
            fugitive_nitrous_oxide.imass['N2O'] = fugitive_nitrous_oxide_storage
        else:
            raise ValueError('climate can only be `humid` or `arid`')
        
        # kg CO2 eq/h
        carbon_dioxide.imass['CO2'] = -self.dry_solids*1000*self.unit_carbon_sequestration
    
    @property
    def nitrogen_mass_flow(self):
        # kg N/h
        return self.ins[0].imass['Sludge_protein']*self.protein_2_N
    
    @property
    def phosphorus_mass_flow(self):
        # kg P/h
        return self.ins[0].imass['Sludge_protein']*self.protein_2_N*self.N_2_P
    
    @property
    def unit_process(self):
        return 'land application'

# =============================================================================
# Incineration
# =============================================================================

# cost base year: 2014 (the original conference paper cited in [2] was published probably on Jan, 2015, so 1 year before that)
# n: 0.7753, [2]
# BM: 2, to be consistent with thermochemical units
@cost(ID='Incinerator 1', basis='Half wet mass flowrate', units='tonne/day',
      cost=2350700/2/2, S=1000/365, CE=qs.CEPCI_by_year[2014], n=0.7753, BM=2)
@cost(ID='Incinerator 2', basis='Half wet mass flowrate', units='tonne/day',
      cost=2350700/2/2, S=1000/365, CE=qs.CEPCI_by_year[2014], n=0.7753, BM=2)
class Incineration(SanUnit):
    '''
    Incineration of dried solids. N2O emission in the BEAM model is originally
    based on [1].
    
    scope 1 emission: fugitive methane, fugitive nitrous oxide
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    from [2]: I = 2.3507×C0.7753, where I is the investment cost in million
    dollars and C is the plant capacity (1000 metric tons of waste/year).
    # assume CAPEX / installed cost = 2 and installed cost / purchase cost (BM) = 2
    
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
    fugitive_methane_incineration : float
        fugitive methane during incineration, [tonne CH4/dry tonne solids].
        note this value assumes 20% solids (instead of 80% which is more common
        after heat drying), but this is based on the drt weight and therefore,
        can be used as an approximation.
    protein_2_N, float
        N to protein ratio, [-].
    N_2_P, float
        P to N ratio, [-].
    suzuki_constant_1 : float
        for nitrous oxide calculation, [-].
    suzuki_constant_2 : float
        for nitrous oxide calculation, [-].
    suzuki_lowest_temperature : float
        the lowest temperature to apply Suzuki equation, [°C].
    SNCR : bool
        whether the urea-based selective noncatalytic reduction emissions system is used, [True, False]
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
    
    References
    ----------
    [1] Suzuki, Y.; Ochi, S.-I.; Kawashima, Y.; Hiraide, R. Determination of
        Emission Factors of Nitrous Oxide from Fluidized Bed Sewage Sludge
        Incinerators by Long-Term Continuous Monitoring. Journal of Chemical
        Engineering of Japan 2003, 36 (4), 458–463.
        https://doi.org/10.1252/jcej.36.458.
    [2] Gergel, I. Cost of incineration plant. Waste To Energy International.
        https://wteinternational.com/news/cost-of-incineration-plant/
        (accessed 2025-08-06).
    '''
    _N_ins = 2
    _N_outs = 4
    
    _units= {'Half wet mass flowrate':'tonne/day'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20,
                 incineration_temperature=850, ash_moisture_content=0.01,
                 heat_water_removal=4.5, fugitive_methane_incineration=0.0000097,
                 protein_2_N=0.159, N_2_P=0.3927, suzuki_constant_1=161.3,
                 suzuki_constant_2=0.14, suzuki_lowest_temperature=750,
                 SNCR=False, N2O_adjustment_factor_urea=0.2,
                 heat_solids_incineration=12000, additional_fuel_ratio=0,
                 heat_recovery_ratio=0.5, heat_recovery_efficiency=0.8,
                 electricity_recovery_ratio=0,
                 BTU_to_kWh_with_efficiency=0.0000854, net_capacity_factor=0.85,
                 unit_electricity=200):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.incineration_temperature = incineration_temperature
        self.ash_moisture_content = ash_moisture_content
        self.heat_water_removal = heat_water_removal
        self.fugitive_methane_incineration = fugitive_methane_incineration
        self.protein_2_N = protein_2_N
        self.N_2_P = N_2_P
        self.suzuki_constant_1 = suzuki_constant_1
        self.suzuki_constant_2 = suzuki_constant_2
        self.suzuki_lowest_temperature = suzuki_lowest_temperature
        self.SNCR = SNCR
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
        
        natural_gas.phase = 'g'
        ash.phase = 's'
        vapor.phase = 'g'
        fugitive_methane.phase = 'g'
        fugitive_nitrous_oxide.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (dried_solids.F_mass - dried_solids.imass['H2O'])/1000
        
        ash.imass['Sludge_ash'] = dried_solids.imass['Sludge_ash']
        ash.imass['H2O'] = ash.imass['Sludge_ash']/(1 - self.ash_moisture_content)*self.ash_moisture_content
        
        vapor.imass['H2O'] = dried_solids.imass['H2O'] - ash.imass['H2O']
        
        fugitive_methane.imass['CH4'] = self.dry_solids*1000*self.fugitive_methane_incineration
        
        # dry tonne N/day
        N_mass_flow = dried_solids.imass['Sludge_protein']*self.protein_2_N/1000*24
        
        if self.suzuki_constant_1 - (self.suzuki_constant_2*(self.incineration_temperature + _C_to_K))/100 < 0:
            # tonne N2O/day
            N2O_before_adjustment = 0
        else:
            # tonne N2O/day
            N2O_before_adjustment = (N_mass_flow*(self.suzuki_constant_1 - (self.suzuki_constant_2*(max(self.incineration_temperature, self.suzuki_lowest_temperature) + _C_to_K)))/100*_N_to_N2O)
        
        if self.SNCR:
            # tonne N2O/day
            N2O_urea_catalyst = N2O_before_adjustment*self.N2O_adjustment_factor_urea
        else:
            # tonne N2O/day
            N2O_urea_catalyst = 0
        
        if self.dry_solids*1000/dried_solids.F_mass >= 0.87:
            # - 
            N2O_reduction_ratio = 0.6
        elif self.dry_solids*1000/dried_solids.F_mass >= 0.24:
            # - 
            N2O_reduction_ratio = 0.5
        else:
            # - 
            N2O_reduction_ratio = 0
        
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
        
        # assume excess heat, if any, is lost
        natural_gas.ivol['CH4'] = max(0, natural_gas_net)
        
        # kW
        electricity_generated = energy_solids_incineration*self.electricity_recovery_ratio/_BTU_to_MJ*self.BTU_to_kWh_with_efficiency*self.net_capacity_factor
        
        # kW
        electricity_requirement = self.dry_solids*self.unit_electricity
        
        # kW
        self.electricity_net = electricity_requirement - electricity_generated
        
    def _design(self):
        D = self.design_results
        D['Half wet mass flowrate'] = self.ins[0].F_mass/1000*24/2
        
        # kW
        self.add_power_utility(self.electricity_net)
        
    @property
    def moisture(self):
        return self.ash_moisture_content
    
    @property
    def phosphorus_mass_flow(self):
        # kg P/h
        return self.ins[0].imass['Sludge_protein']*self.protein_2_N*self.N_2_P
    
    @property
    def unit_process(self):
        return 'incineration'

# =============================================================================
# thermochemical units
# =============================================================================

# a brief comparison of thermochemical units: Table 2 in https://pubs.rsc.org/en/content/articlelanding/2024/ew/d4ew00278d/unauth

# =============================================================================
# HydrothermalLiquefaction
# =============================================================================

# based on [1]
@cost(ID='HTL system 1', basis='Half wet mass flowrate', units='lb/h',
      cost=18743378, S=306198, CE=qs.CEPCI_by_year[2011], n=0.77, BM=2.1)
@cost(ID='HTL system 2', basis='Half wet mass flowrate', units='lb/h',
      cost=18743378, S=306198, CE=qs.CEPCI_by_year[2011], n=0.77, BM=2.1)
@cost(ID='Solids filter oil/water separator', basis='Wet mass flowrate', units='lb/h',
      cost=3945523, S=1219765, CE=qs.CEPCI_by_year[2011], n=0.68, BM=1.9)
@cost(ID='Hot oil system', basis='Wet mass flowrate', units='lb/h',
      cost=4670532, S=306198, CE=qs.CEPCI_by_year[2011], n=0.6, BM=1.4)
class HydrothermalLiquefaction(SanUnit):
    '''
    HTL converts feedstock to gas, aqueous, biocrude, (hydro)char
    under elevated temperature and pressure.
    
    scope 1 emission: N/A
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        dewatered_solids.
    outs : iterable
        biocrude, HTLaqueous, hydrochar, offgas.
    lipid_2_C : float
        lipid to C ratio, [-].
    protein_2_C : float
        protein to C ratio, [-].
    carbo_2_C : float
        carbohydrate to C ratio, [-].
    lipid_2_H : float
        lipid to H ratio, [-].
    protein_2_H : float
        protein to H ratio, [-].
    carbo_2_H : float
        carbohydrate to H ratio, [-].
    protein_2_N : float
        N to protein ratio, [-].
    N_2_P : float
        P to N ratio, [-].
    eff_P : float
        HTL effluent pressure, [Pa].
    crude_oil_density : float
        density of crude oil, [kg/m3].
    crude_oil_HHV : float
        HHV of crude oil, [MJ/kg].
    biocrude_density : float
        density of biocrude, [kg/m3].
    biocrude_distance : float
        distance between WRRFs and oil refineries, [km].
    FOAK : bool
        whether the facility is first-of-a-king, [True, False].
    pctnew : float
        percentage of capital cost of commercially undemonstrated equipment, [%].
    impurities : int (float)
        impurities buildup and corrosion issues, ranging from 0 to 5, [-].
    complexity : int (float)
        number of continuously linked process steps, ranging from 1 to 11, [-].
    inclusiveness : float
        percentage of pre-startup personnel, inventory, and land purchase costs, [%].
    project_definition : int (float)
        levels of site-specific information and engineering included in estimate, ranging from 2 to 8, [-].
    newsteps : int (float)
        number of processes incorporating commercially unproven technologies, [-].
    baleqs : float
        percentage of heat and mass balance equations based on actual data from prior plants, [%].
    waste : int (float)
        difficulties with waste handling encountered during development, ranging from 0 to 5, [-].
    solids_handling : bool
        True for solids handling, False for no solids handling, [True, False].
    
    See Also
    --------
    :class:`qsdsan.sanunits.HydrothermalLiquefaction`
    
    :class:`exposan.saf._units.HydrothermalLiquefaction`
    
    References
    ----------
    [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    '''
    _N_ins = 1
    _N_outs = 4
    
    _units= {'Half wet mass flowrate':'lb/h',
             'Wet mass flowrate':'lb/h'}
    
    auxiliary_unit_names = ('pump','hx','inf_hx','eff_hx','kodrum')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20, lipid_2_C=0.750,
                 protein_2_C=0.545, carbo_2_C=0.400, lipid_2_H=0.125,
                 protein_2_H=0.068, carbo_2_H=0.067, protein_2_N=0.159,
                 N_2_P=0.3927, lipid_2_biocrude=0.846, protein_2_biocrude=0.445,
                 carbo_2_biocrude=0.205, protein_2_gas=0.074, carbo_2_gas=0.418,
                 biocrude_C_slope=-8.37, biocrude_C_intercept=68.55,
                 biocrude_N_slope=0.133, biocrude_H_slope=-2.61,
                 biocrude_H_intercept=8.20, HTLaqueous_C_slope=478,
                 TOC_TC=0.764, hydrochar_C_slope=1.75,
                 biocrude_moisture_content=0.063,
                 hydrochar_P_recovery_ratio=0.86,
                 gas_composition={'CH4': 0.050, 'C2H6': 0.032, 'CO2': 0.918},
                 T=350 + _C_to_K, P=3049.7*_psi_to_Pa, eff_T=60 + _C_to_K,
                 eff_P=30*_psi_to_Pa,
                 # https://www.transmountain.com/about-petroleum-liquids (accessed 2025-02-05)
                 crude_oil_density=850,
                 # crude oil HHV: https://world-nuclear.org/information-library/facts-and-figures/heat-values-of-various-fuels
                 crude_oil_HHV=44.5,
                 # https://doi.org/10.2172/1897670
                 biocrude_density=983, biocrude_distance=100, FOAK=True,
                 pctnew=41, impurities=4, complexity=6, inclusiveness=33,
                 project_definition=5, newsteps=3, baleqs=0, waste=2,
                 solids_handling=True):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.lipid_2_C = lipid_2_C
        self.protein_2_C = protein_2_C
        self.carbo_2_C = carbo_2_C
        self.lipid_2_H = lipid_2_H
        self.protein_2_H = protein_2_H
        self.carbo_2_H = carbo_2_H
        self.protein_2_N = protein_2_N
        self.N_2_P = N_2_P
        self.lipid_2_biocrude = lipid_2_biocrude
        self.protein_2_biocrude = protein_2_biocrude
        self.carbo_2_biocrude = carbo_2_biocrude
        self.protein_2_gas = protein_2_gas
        self.carbo_2_gas = carbo_2_gas
        self.biocrude_C_slope = biocrude_C_slope
        self.biocrude_C_intercept = biocrude_C_intercept
        self.biocrude_N_slope = biocrude_N_slope
        self.biocrude_H_slope = biocrude_H_slope
        self.biocrude_H_intercept = biocrude_H_intercept
        self.HTLaqueous_C_slope = HTLaqueous_C_slope
        self.TOC_TC = TOC_TC
        self.hydrochar_C_slope = hydrochar_C_slope
        self.biocrude_moisture_content = biocrude_moisture_content
        self.hydrochar_P_recovery_ratio = hydrochar_P_recovery_ratio
        self.gas_composition = gas_composition
        self.T = T
        self.P = P
        self.eff_T = eff_T
        self.eff_P = eff_P
        pump_in = Stream(f'{ID}_pump_in')
        pump_out = Stream(f'{ID}_pump_out')
        # purchase costs for Pump are related to bst.CE
        self.pump = Pump(ID=f'.{ID}_pump', ins=pump_in, outs=pump_out, P=P)
        inf_pre_hx = Stream(f'{ID}_inf_pre_hx')
        eff_pre_hx = Stream(f'{ID}_eff_pre_hx')
        inf_after_hx = Stream(f'{ID}_inf_after_hx')
        eff_after_hx = Stream(f'{ID}_eff_after_hx')
        self.hx = HXprocess(ID=f'.{ID}_hx', ins=(inf_pre_hx, eff_pre_hx), outs=(inf_after_hx, eff_after_hx))
        inf_hx_out = Stream(f'{ID}_inf_hx_out')
        self.inf_hx = HXutility(ID=f'.{ID}_inf_hx', ins=inf_after_hx, outs=inf_hx_out, T=T, rigorous=True)
        self._inf_at_temp = Stream(f'{ID}_inf_at_temp')
        self._eff_at_temp = Stream(f'{ID}_eff_at_temp')
        eff_hx_out = Stream(f'{ID}_eff_hx_out')
        self.eff_hx = HXutility(ID=f'.{ID}_eff_hx', ins=eff_after_hx, outs=eff_hx_out, T=eff_T, rigorous=True)
        self.crude_oil_density = crude_oil_density
        self.crude_oil_HHV = crude_oil_HHV
        self.biocrude_density = biocrude_density
        self.biocrude_distance = biocrude_distance
        self.FOAK = FOAK
        self.pctnew = pctnew
        self.impurities = impurities
        self.complexity = complexity
        self.inclusiveness = inclusiveness
        self.project_definition = project_definition
        self.newsteps = newsteps
        self.baleqs = baleqs
        self.waste = waste
        self.solids_handling = solids_handling
        if self.FOAK:
            self.plant_performance_factor = 85.77 - 9.69*self.newsteps + 0.33*self.baleqs -\
                                            4.12*self.waste - 17.91*self.solids_handling
        else:
            self.plant_performance_factor = 100
    
    def _run(self):
        dewatered_solids = self.ins[0]
        biocrude, HTLaqueous, hydrochar, offgas = self.outs
        
        biocrude.phase = 'l'
        HTLaqueous.phase = 'l'
        hydrochar.phase = 's'
        offgas.phase = 'g'
        
        dewatered_solids_afdw = dewatered_solids.imass['Sludge_lipid'] +\
                                dewatered_solids.imass['Sludge_protein'] +\
                                dewatered_solids.imass['Sludge_carbo']
        
        self.afdw_lipid_ratio = dewatered_solids.imass['Sludge_lipid']/dewatered_solids_afdw
        self.afdw_protein_ratio = dewatered_solids.imass['Sludge_lipid']/dewatered_solids_afdw
        self.afdw_carbo_ratio = 1 - self.afdw_lipid_ratio - self.afdw_protein_ratio
        
        # the following calculations are based on revised MCA model
        # HTLaqueous is TDS in aqueous phase
        # 0.377, 0.481, and 0.154 don't have uncertainties because they are calculated values
        
        biocrude.imass['Biocrude'] = (self.protein_2_biocrude*self.afdw_protein_ratio +\
                                      self.lipid_2_biocrude*self.afdw_lipid_ratio +\
                                      self.carbo_2_biocrude*self.afdw_carbo_ratio)*\
                                      dewatered_solids_afdw
        biocrude.imass['H2O'] = biocrude.imass['Biocrude']/(1 -\
                                self.biocrude_moisture_content) -\
                                biocrude.imass['Biocrude']
        
        HTLaqueous.imass['HTLaqueous'] = (0.481*self.afdw_protein_ratio +\
                                          0.154*self.afdw_lipid_ratio)*\
                                          dewatered_solids_afdw
        
        hydrochar.imass['Hydrochar'] = 0.377*self.afdw_carbo_ratio*dewatered_solids_afdw
        
        gas_mass_flow = (self.protein_2_gas*self.afdw_protein_ratio + self.carbo_2_gas*self.afdw_carbo_ratio)*\
                        dewatered_solids_afdw
        
        for name, ratio in self.gas_composition.items():
            offgas.imass[name] = gas_mass_flow*ratio
        
        # assume ash (all soluble based on [1]) goes to water
        HTLaqueous.imass['H2O'] = dewatered_solids.F_mass - hydrochar.F_mass -\
                                  biocrude.F_mass - gas_mass_flow - HTLaqueous.imass['HTLaqueous']
        
        for i in self.outs:
            i.T = self.T
            i.P = self.P
        
        self._eff_at_temp.mix_from([HTLaqueous, biocrude, offgas], vle=True)
        
        for attr, val in zip(('T','P'), (self.eff_T, self.eff_P)):
            if val:
                for i in self.outs: setattr(i, attr, val)
        
        # for properties
        self.solids_dw_protein = dewatered_solids.imass['Sludge_protein']/(dewatered_solids.F_mass -\
                                                                           dewatered_solids.imass['H2O'])
        self.solids_dw_carbo = dewatered_solids.imass['Sludge_carbo']/(dewatered_solids.F_mass -\
                                                                       dewatered_solids.imass['H2O'])
        self.solids_C = dewatered_solids.imass['Sludge_lipid']*self.lipid_2_C +\
                        dewatered_solids.imass['Sludge_protein']*self.protein_2_C +\
                        dewatered_solids.imass['Sludge_carbo']*self.carbo_2_C
        self.solids_H = dewatered_solids.imass['Sludge_lipid']*self.lipid_2_H +\
                        dewatered_solids.imass['Sludge_protein']*self.protein_2_H +\
                        dewatered_solids.imass['Sludge_carbo']*self.carbo_2_H
        self.solids_N = dewatered_solids.imass['Sludge_protein']*self.protein_2_N
        self.solids_P = self.solids_N*self.N_2_P
        self.solids_O = dewatered_solids.F_mass - dewatered_solids.imass['H2O'] -\
                        dewatered_solids.imass['Sludge_ash'] - self.solids_C -\
                        self.solids_H - self.solids_N
        self.AOSc = (3*self.solids_N/14.0067 + 2*self.solids_O/15.999 -\
                     self.solids_H/1.00784)/(self.solids_C/12.011)
    
    def _design(self):
        D = self.design_results
        D['Half wet mass flowrate'] = self.ins[0].F_mass*_kg_to_lb/2
        D['Wet mass flowrate'] = self.ins[0].F_mass*_kg_to_lb
        
        pump = self.pump
        pump.ins[0].copy_like(self.ins[0])
        pump.simulate()
        
        hx = self.hx
        inf_hx = self.inf_hx
        inf_hx_in, inf_hx_out = inf_hx.ins[0], inf_hx.outs[0]
        inf_pre_hx, eff_pre_hx = hx.ins
        inf_after_hx, eff_after_hx = hx.outs
        inf_pre_hx.copy_like(self.ins[0])
        eff_pre_hx.copy_like(self._eff_at_temp)
        
        # use product to heat up influent
        hx.T_lim1 = self.eff_T
        hx.simulate()
        for i in self.outs:
            i.T = eff_after_hx.T
        
        # additional inf HX
        inf_hx_in.copy_like(inf_after_hx)
        inf_hx_out.copy_flow(inf_hx_in)
        inf_hx_out.T = self.T
        inf_hx.simulate_as_auxiliary_exchanger(ins=inf_hx.ins, outs=inf_hx.outs)
        
        # additional eff HX
        eff_hx = self.eff_hx
        eff_hx_in, eff_hx_out = eff_hx.ins[0], eff_hx.outs[0]
        eff_hx_in.copy_like(eff_after_hx)
        eff_hx_out.mix_from(self.outs)
        eff_hx_out.T = self.eff_T
        eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs)
        
        for i in self.outs:
            i.T = self.eff_T
        
        # original code from exposan.saf._units, cooling duty very high
        # duty = self.Hnet + eff_hx.Hnet
        # eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs, duty=duty)
    
    def _cost(self):
        self.hx.baseline_purchase_costs.clear()
        self.inf_hx.baseline_purchase_costs.clear()
        self.eff_hx.baseline_purchase_costs.clear()
        self._decorated_cost()
        
        if self.FOAK:
            self.cost_growth_factor = 1.12196 - 0.00297*self.pctnew -\
                                      0.02125*self.impurities - 0.01137*self.complexity +\
                                      0.00111*self.inclusiveness - 0.06361*self.project_definition
        else:
            self.cost_growth_factor = 1
        
        for C in [self.baseline_purchase_costs] + [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]/self.cost_growth_factor
    
    @property
    def biocrude_yield(self):
        return self.protein_2_biocrude*self.afdw_protein_ratio +\
               self.lipid_2_biocrude*self.afdw_lipid_ratio +\
               self.carbo_2_biocrude*self.afdw_carbo_ratio
    
    @property
    def aqueous_yield(self):
        return 0.481*self.afdw_protein_ratio + 0.154*self.afdw_lipid_ratio
    
    @property
    def hydrochar_yield(self):
        return 0.377*self.afdw_carbo_ratio
    
    @property
    def gas_yield(self):
        return self.protein_2_gas*self.afdw_protein_ratio + self.carbo_2_gas*self.afdw_carbo_ratio
    
    @property
    def biocrude_C_ratio(self):
        return (self.AOSc*self.biocrude_C_slope + self.biocrude_C_intercept)/100
    
    @property
    def biocrude_H_ratio(self):
        return (self.AOSc*self.biocrude_H_slope + self.biocrude_H_intercept)/100
    
    @property
    def biocrude_N_ratio(self):
        return self.biocrude_N_slope*self.solids_dw_protein
    
    @property
    def biocrude_C(self):
        return min(self.outs[0].F_mass*self.biocrude_C_ratio, self.solids_C)
    
    @property
    def HTLaqueous_C(self):
        return min(self.outs[1].F_vol*1000*self.HTLaqueous_C_slope*\
                   self.solids_dw_protein*100/1000000/self.TOC_TC,
                   self.solids_C - self.biocrude_C)
    
    @property
    def biocrude_H(self):
        return self.outs[0].F_mass*self.biocrude_H_ratio
    
    @property
    def biocrude_N(self):
        return min(self.outs[0].F_mass*self.biocrude_N_ratio, self.solids_N)
    
    @property
    def biocrude_HHV(self):
        return 30.74 - 8.52*self.AOSc +\
               0.024*self.solids_dw_protein
    
    @property
    def offgas_C(self):
        carbon = sum(self.outs[3].imass[self.gas_composition]*
                     [cmp.i_C for cmp in self.components[self.gas_composition]])
        return min(carbon, self.solids_C - self.biocrude_C - self.HTLaqueous_C)
    
    @property
    def hydrochar_C_ratio(self):
        return min(self.hydrochar_C_slope*self.solids_dw_carbo, 0.65)
    
    @property
    def hydrochar_C(self):
        return min(self.outs[2].F_mass*self.hydrochar_C_ratio, self.solids_C -\
                   self.biocrude_C - self.HTLaqueous_C - self.offgas_C)
    
    @property
    def hydrochar_P(self):
        return min(self.solids_P*self.hydrochar_P_recovery_ratio, self.outs[2].F_mass)
    
    @property
    def HTLaqueous_N(self):
        return self.solids_N - self.biocrude_N
    
    @property
    def HTLaqueous_P(self):
        return self.solids_P*(1 - self.hydrochar_P_recovery_ratio)

# =============================================================================
# HydrothermalAlkalineTreatment
# =============================================================================

# assume to be the same as HTL
@cost(ID='HTL system 1', basis='Half wet mass flowrate', units='lb/h',
      cost=18743378, S=306198, CE=qs.CEPCI_by_year[2011], n=0.77, BM=2.1)
@cost(ID='HTL system 2', basis='Half wet mass flowrate', units='lb/h',
      cost=18743378, S=306198, CE=qs.CEPCI_by_year[2011], n=0.77, BM=2.1)
@cost(ID='Solids filter oil/water separator', basis='Wet mass flowrate', units='lb/h',
      cost=3945523, S=1219765, CE=qs.CEPCI_by_year[2011], n=0.68, BM=1.9)
@cost(ID='Hot oil system', basis='Wet mass flowrate', units='lb/h',
      cost=4670532, S=306198, CE=qs.CEPCI_by_year[2011], n=0.6, BM=1.4)
class HydrothermalAlkalineTreatment(SanUnit):
    '''
    HALT converts feedstock to gas, aqueous, biocrude, (hydro)char
    under elevated temperature and pressure. Land application of (hydro)char is
    included.
    
    Note diesel is used for land application and is therefore not listed here,
    though it is included in the model
    scope 1 emission: N/A
    scope 2 emission: electricity
    scope 3 emission: NaOH, HCl
    
    Parameters
    ----------
    ins : iterable
        dewatered_solids, sodium_hydroxide, hydrochloric_acid, diesel.
    outs : iterable
        biocrude, HALTaqueous, hydrochar, offgas.
    NaOH_conc : float
        the concentration of NaOH, [M].
    HALT_biocrude_adjustment_factor : float
        the biocrude yield in HALT compared to that in HTL, [-].
    HALT_hydrochar_adjustment_factor : float
        the hydrochar yield in HALT compared to that in HTL, [-].
    load_size : float
        size of loads per truck, [m3/load].
    load_frequency : float
        load frequency, [load/h].
    tractor_fuel : float
        tractor fuel, [L diesel/h].
    crude_oil_density : float
        density of crude oil, [kg/m3].
    crude_oil_HHV : float
        HHV of crude oil, [MJ/kg].
    biocrude_density : float
        density of biocrude, [kg/m3].
    biocrude_distance : float
        distance between WRRFs and oil refineries, [km].
    hydrochar_distance : float
        distance between WRRFs and land application sites, [km].
    FOAK : bool
        whether the facility is first-of-a-king, [True, False].
    pctnew : float
        percentage of capital cost of commercially undemonstrated equipment, [%].
    impurities : int (float)
        impurities buildup and corrosion issues, ranging from 0 to 5, [-].
    complexity : int (float)
        number of continuously linked process steps, ranging from 1 to 11, [-].
    inclusiveness : float
        percentage of pre-startup personnel, inventory, and land purchase costs, [%].
    project_definition : int (float)
        levels of site-specific information and engineering included in estimate, ranging from 2 to 8, [-].
    newsteps : int (float)
        number of processes incorporating commercially unproven technologies, [-].
    baleqs : float
        percentage of heat and mass balance equations based on actual data from prior plants, [%].
    waste : int (float)
        difficulties with waste handling encountered during development, ranging from 0 to 5, [-].
    solids_handling : bool
        True for solids handling, False for no solids handling, [True, False].
    
    See Also
    --------
    :class:`exposan.landscape_sanunits.HydrothermalLiquefaction`
    '''
    _N_ins = 4
    _N_outs = 4
    
    _units= {'Half wet mass flowrate':'lb/h',
             'Wet mass flowrate':'lb/h'}
    
    auxiliary_unit_names = ('pump','hx','inf_hx','eff_hx','kodrum')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=15, NaOH_conc=2,
                 lipid_2_C=0.750, protein_2_C=0.545, carbo_2_C=0.400,
                 lipid_2_H=0.125, protein_2_H=0.068, carbo_2_H=0.067,
                 protein_2_N=0.159, N_2_P=0.3927, lipid_2_biocrude=0.846,
                 protein_2_biocrude=0.445, carbo_2_biocrude=0.205,
                 protein_2_gas=0.074, carbo_2_gas=0.418, biocrude_C_slope=-8.37,
                 biocrude_C_intercept=68.55, biocrude_N_slope=0.133,
                 biocrude_H_slope=-2.61, biocrude_H_intercept=8.20,
                 HTLaqueous_C_slope=478, TOC_TC=0.764, hydrochar_C_slope=1.75,
                 biocrude_moisture_content=0.063,
                 hydrochar_P_recovery_ratio=0.86,
                 gas_composition={'CH4': 0.050, 'C2H6': 0.032, 'CO2': 0.918},
                 # Table 1, A.J.K. modeling paper
                 HALT_biocrude_adjustment_factor=0.85,
                 # Table 1, A.J.K. modeling paper
                 HALT_hydrochar_adjustment_factor=0.2,
                 T=350 + _C_to_K, P=3049.7*_psi_to_Pa, eff_T=60 + _C_to_K,
                 eff_P=30*_psi_to_Pa, load_size=13, load_frequency=3,
                 tractor_fuel=25, unit_carbon_sequestration=0.745,
                 # https://www.transmountain.com/about-petroleum-liquids (accessed 2025-02-05)
                 crude_oil_density=850,
                 # crude oil HHV: https://world-nuclear.org/information-library/facts-and-figures/heat-values-of-various-fuels
                 crude_oil_HHV=44.5,
                 # https://doi.org/10.2172/1897670
                 biocrude_density=983,  biocrude_distance=100,
                 hydrochar_distance=100, FOAK=True, pctnew=55, impurities=5,
                 complexity=7, inclusiveness=33, project_definition=6,
                 newsteps=5, baleqs=0, waste=2, solids_handling=True):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.NaOH_conc = NaOH_conc
        self.lipid_2_C = lipid_2_C
        self.protein_2_C = protein_2_C
        self.carbo_2_C = carbo_2_C
        self.lipid_2_H = lipid_2_H
        self.protein_2_H = protein_2_H
        self.carbo_2_H = carbo_2_H
        self.protein_2_N = protein_2_N
        self.N_2_P = N_2_P
        self.lipid_2_biocrude = lipid_2_biocrude
        self.protein_2_biocrude = protein_2_biocrude
        self.carbo_2_biocrude = carbo_2_biocrude
        self.protein_2_gas = protein_2_gas
        self.carbo_2_gas = carbo_2_gas
        self.biocrude_C_slope = biocrude_C_slope
        self.biocrude_C_intercept = biocrude_C_intercept
        self.biocrude_N_slope = biocrude_N_slope
        self.biocrude_H_slope = biocrude_H_slope
        self.biocrude_H_intercept = biocrude_H_intercept
        self.HTLaqueous_C_slope = HTLaqueous_C_slope
        self.TOC_TC = TOC_TC
        self.hydrochar_C_slope = hydrochar_C_slope
        self.biocrude_moisture_content = biocrude_moisture_content
        self.hydrochar_P_recovery_ratio = hydrochar_P_recovery_ratio
        self.HALT_biocrude_adjustment_factor = HALT_biocrude_adjustment_factor
        self.HALT_hydrochar_adjustment_factor = HALT_hydrochar_adjustment_factor
        self.gas_composition = gas_composition
        self.T = T
        self.P = P
        self.eff_T = eff_T
        self.eff_P = eff_P
        pump_in = Stream(f'{ID}_pump_in')
        pump_out = Stream(f'{ID}_pump_out')
        # purchase costs for Pump are related to bst.CE
        self.pump = Pump(ID=f'.{ID}_pump', ins=pump_in, outs=pump_out, P=P)
        inf_pre_hx = Stream(f'{ID}_inf_pre_hx')
        eff_pre_hx = Stream(f'{ID}_eff_pre_hx')
        inf_after_hx = Stream(f'{ID}_inf_after_hx')
        eff_after_hx = Stream(f'{ID}_eff_after_hx')
        self.hx = HXprocess(ID=f'.{ID}_hx', ins=(inf_pre_hx, eff_pre_hx), outs=(inf_after_hx, eff_after_hx))
        inf_hx_out = Stream(f'{ID}_inf_hx_out')
        self.inf_hx = HXutility(ID=f'.{ID}_inf_hx', ins=inf_after_hx, outs=inf_hx_out, T=T, rigorous=True)
        self._inf_at_temp = Stream(f'{ID}_inf_at_temp')
        self._eff_at_temp = Stream(f'{ID}_eff_at_temp')
        eff_hx_out = Stream(f'{ID}_eff_hx_out')
        self.eff_hx = HXutility(ID=f'.{ID}_eff_hx', ins=eff_after_hx, outs=eff_hx_out, T=eff_T, rigorous=True)
        self.load_size = load_size
        self.load_frequency = load_frequency
        self.tractor_fuel = tractor_fuel
        self.crude_oil_density = crude_oil_density
        self.crude_oil_HHV = crude_oil_HHV
        self.biocrude_density = biocrude_density
        self.biocrude_distance = biocrude_distance
        self.hydrochar_distance = hydrochar_distance
        self.FOAK = FOAK
        self.pctnew = pctnew
        self.impurities = impurities
        self.complexity = complexity
        self.inclusiveness = inclusiveness
        self.project_definition = project_definition
        self.newsteps = newsteps
        self.baleqs = baleqs
        self.waste = waste
        self.solids_handling = solids_handling
        if self.FOAK:
            self.plant_performance_factor = 85.77 - 9.69*self.newsteps + 0.33*self.baleqs -\
                                            4.12*self.waste - 17.91*self.solids_handling
        else:
            self.plant_performance_factor = 100
    
    def _run(self):
        dewatered_solids, sodium_hydroxide, hydrochloric_acid, diesel = self.ins
        biocrude, HALTaqueous, hydrochar, offgas = self.outs
        
        biocrude.phase = 'l'
        HALTaqueous.phase = 'l'
        hydrochar.phase = 's'
        offgas.phase = 'g'
        
        # this just includes the land application onsite spreading and does not include the transporation from WRRFs to land application sites
        diesel.imass['Diesel'] = (self.ins[0].F_vol/self.load_size/self.load_frequency*self.tractor_fuel)/1000*diesel_density
        
        dewatered_solids_afdw = dewatered_solids.imass['Sludge_lipid'] +\
                                dewatered_solids.imass['Sludge_protein'] +\
                                dewatered_solids.imass['Sludge_carbo']
        
        self.afdw_lipid_ratio = dewatered_solids.imass['Sludge_lipid']/dewatered_solids_afdw
        self.afdw_protein_ratio = dewatered_solids.imass['Sludge_lipid']/dewatered_solids_afdw
        self.afdw_carbo_ratio = 1 - self.afdw_lipid_ratio - self.afdw_protein_ratio
        
        sodium_hydroxide.imass['NaOH'] = dewatered_solids.ivol['H2O']*1000*self.NaOH_conc*40/1000
        hydrochloric_acid.imass['HCl'] = sodium_hydroxide.imass['NaOH']*36.5/40
        
        # change in gas composition with NaOH concentration
        # at NaOH 1.67 M and higher, all gases are H2; below, assumed proportional amounts of HTL mixture and H2
        if self.NaOH_conc <= 1.67:
            NaOH_factor = self.NaOH_conc/1.67
            self.gas_composition['CH4'] = self.gas_composition['CH4'] * (1 - NaOH_factor)
            self.gas_composition['C2H6'] = self.gas_composition['C2H6'] * (1 - NaOH_factor)
            self.gas_composition['CO2'] = self.gas_composition['CO2'] * (1 - NaOH_factor)
            self.gas_composition['H2'] = NaOH_factor
        else:
            # TODO: in the HTL PFAS code, update 'H2':100.00 to 'H2':1
            self.gas_composition={'CH4': 0, 'C2H6': 0,'CO2': 0, 'H2': 1}
        
        # the following calculations are based on revised MCA model
        # HALTaqueous is TDS in aqueous phase
        # 0.377, 0.481, and 0.154 don't have uncertainties because they are calculated values
        
        biocrude.imass['Biocrude'] = (self.protein_2_biocrude*self.afdw_protein_ratio +\
                                      self.lipid_2_biocrude*self.afdw_lipid_ratio +\
                                      self.carbo_2_biocrude*self.afdw_carbo_ratio)*\
                                      dewatered_solids_afdw*self.HALT_biocrude_adjustment_factor
        biocrude.imass['H2O'] = biocrude.imass['Biocrude']/(1 -\
                                self.biocrude_moisture_content) -\
                                biocrude.imass['Biocrude']
        
        hydrochar.imass['Hydrochar'] = 0.377*self.afdw_carbo_ratio*dewatered_solids_afdw*self.HALT_hydrochar_adjustment_factor
        
        gas_mass_flow = (self.protein_2_gas*self.afdw_protein_ratio + self.carbo_2_gas*self.afdw_carbo_ratio)*\
                        dewatered_solids_afdw
        
        for name, ratio in self.gas_composition.items():
            offgas.imass[name] = gas_mass_flow*ratio
        
        HALTaqueous.imass['HTLaqueous'] = dewatered_solids_afdw - biocrude.imass['Biocrude'] - hydrochar.imass['Hydrochar'] - gas_mass_flow
        # assume ash (all soluble based on [1]) goes to water
        HALTaqueous.imass['H2O'] = dewatered_solids.F_mass + sodium_hydroxide.F_mass + hydrochloric_acid.F_mass - hydrochar.F_mass -\
                                  biocrude.F_mass - gas_mass_flow - HALTaqueous.imass['HTLaqueous']
        
        for i in self.outs:
            i.T = self.T
            i.P = self.P
        
        self._eff_at_temp.mix_from([HALTaqueous, biocrude, offgas], vle=True)
        
        for attr, val in zip(('T','P'), (self.eff_T, self.eff_P)):
            if val:
                for i in self.outs: setattr(i, attr, val)
        
        # for properties
        self.solids_dw_protein = dewatered_solids.imass['Sludge_protein']/(dewatered_solids.F_mass -\
                                                                           dewatered_solids.imass['H2O'])
        self.solids_dw_carbo = dewatered_solids.imass['Sludge_carbo']/(dewatered_solids.F_mass -\
                                                                       dewatered_solids.imass['H2O'])
        self.solids_C = dewatered_solids.imass['Sludge_lipid']*self.lipid_2_C +\
                        dewatered_solids.imass['Sludge_protein']*self.protein_2_C +\
                        dewatered_solids.imass['Sludge_carbo']*self.carbo_2_C
        self.solids_H = dewatered_solids.imass['Sludge_lipid']*self.lipid_2_H +\
                        dewatered_solids.imass['Sludge_protein']*self.protein_2_H +\
                        dewatered_solids.imass['Sludge_carbo']*self.carbo_2_H
        self.solids_N = dewatered_solids.imass['Sludge_protein']*self.protein_2_N
        self.solids_P = self.solids_N*self.N_2_P
        self.solids_O = dewatered_solids.F_mass - dewatered_solids.imass['H2O'] -\
                        dewatered_solids.imass['Sludge_ash'] - self.solids_C -\
                        self.solids_H - self.solids_N
        self.AOSc = (3*self.solids_N/14.0067 + 2*self.solids_O/15.999 -\
                     self.solids_H/1.00784)/(self.solids_C/12.011)
    
    def _design(self):
        D = self.design_results
        D['Half wet mass flowrate'] = self.ins[0].F_mass*_kg_to_lb/2
        D['Wet mass flowrate'] = self.ins[0].F_mass*_kg_to_lb
        
        pump = self.pump
        pump.ins[0].copy_like(self.ins[0])
        pump.simulate()
        
        hx = self.hx
        inf_hx = self.inf_hx
        inf_hx_in, inf_hx_out = inf_hx.ins[0], inf_hx.outs[0]
        inf_pre_hx, eff_pre_hx = hx.ins
        inf_after_hx, eff_after_hx = hx.outs
        inf_pre_hx.copy_like(self.ins[0])
        eff_pre_hx.copy_like(self._eff_at_temp)
        
        # use product to heat up influent
        hx.T_lim1 = self.eff_T
        hx.simulate()
        for i in self.outs:
            i.T = eff_after_hx.T
        
        # additional inf HX
        inf_hx_in.copy_like(inf_after_hx)
        inf_hx_out.copy_flow(inf_hx_in)
        inf_hx_out.T = self.T
        inf_hx.simulate_as_auxiliary_exchanger(ins=inf_hx.ins, outs=inf_hx.outs)
        
        # additional eff HX
        eff_hx = self.eff_hx
        eff_hx_in, eff_hx_out = eff_hx.ins[0], eff_hx.outs[0]
        eff_hx_in.copy_like(eff_after_hx)
        eff_hx_out.mix_from(self.outs)
        eff_hx_out.T = self.eff_T
        eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs)
        
        for i in self.outs:
            i.T = self.eff_T
        
        # original code from exposan.saf._units, cooling duty very high
        # duty = self.Hnet + eff_hx.Hnet
        # eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs, duty=duty)
    
    def _cost(self):
        self.hx.baseline_purchase_costs.clear()
        self.inf_hx.baseline_purchase_costs.clear()
        self.eff_hx.baseline_purchase_costs.clear()
        self._decorated_cost()
        
        if self.FOAK:
            self.cost_growth_factor = 1.12196 - 0.00297*self.pctnew -\
                                      0.02125*self.impurities - 0.01137*self.complexity +\
                                      0.00111*self.inclusiveness - 0.06361*self.project_definition
        else:
            self.cost_growth_factor = 1
        
        for C in [self.baseline_purchase_costs] + [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]/self.cost_growth_factor
    
    @property
    def biocrude_yield(self):
        return self.protein_2_biocrude*self.afdw_protein_ratio +\
               self.lipid_2_biocrude*self.afdw_lipid_ratio +\
               self.carbo_2_biocrude*self.afdw_carbo_ratio
    
    @property
    def aqueous_yield(self):
        return 0.481*self.afdw_protein_ratio + 0.154*self.afdw_lipid_ratio
    
    @property
    def hydrochar_yield(self):
        return 0.377*self.afdw_carbo_ratio
    
    @property
    def gas_yield(self):
        return self.protein_2_gas*self.afdw_protein_ratio + self.carbo_2_gas*self.afdw_carbo_ratio
    
    @property
    def biocrude_C_ratio(self):
        return (self.AOSc*self.biocrude_C_slope + self.biocrude_C_intercept)/100
    
    @property
    def biocrude_H_ratio(self):
        return (self.AOSc*self.biocrude_H_slope + self.biocrude_H_intercept)/100
    
    @property
    def biocrude_N_ratio(self):
        return self.biocrude_N_slope*self.solids_dw_protein
    
    @property
    def biocrude_C(self):
        return min(self.outs[0].F_mass*self.biocrude_C_ratio, self.solids_C)
    
    @property
    def HTLaqueous_C(self):
        return min(self.outs[1].F_vol*1000*self.HTLaqueous_C_slope*\
                   self.solids_dw_protein*100/1000000/self.TOC_TC,
                   self.solids_C - self.biocrude_C)
    
    @property
    def biocrude_H(self):
        return self.outs[0].F_mass*self.biocrude_H_ratio
    
    @property
    def biocrude_N(self):
        return min(self.outs[0].F_mass*self.biocrude_N_ratio, self.solids_N)
    
    @property
    def biocrude_HHV(self):
        return 30.74 - 8.52*self.AOSc +\
               0.024*self.solids_dw_protein
    
    @property
    def offgas_C(self):
        carbon = sum(self.outs[3].imass[self.gas_composition]*
                     [cmp.i_C for cmp in self.components[self.gas_composition]])
        return min(carbon, self.solids_C - self.biocrude_C - self.HTLaqueous_C)
    
    @property
    def hydrochar_C_ratio(self):
        return min(self.hydrochar_C_slope*self.solids_dw_carbo, 0.65)
    
    @property
    def hydrochar_C(self):
        return min(self.outs[2].F_mass*self.hydrochar_C_ratio, self.solids_C -\
                   self.biocrude_C - self.HTLaqueous_C - self.offgas_C)
    
    @property
    def hydrochar_P(self):
        return min(self.solids_P*self.hydrochar_P_recovery_ratio, self.outs[2].F_mass)
    
    @property
    def HTLaqueous_N(self):
        return self.solids_N - self.biocrude_N
    
    @property
    def HTLaqueous_P(self):
        return self.solids_P*(1 - self.hydrochar_P_recovery_ratio)

# =============================================================================
# Analyzer
# =============================================================================

class Analyzer(SanUnit):
    '''
    A fake unit that calculates C, N, P, and H2O amount in the HTL/HALT aqueous effluent.
    
    Parameters
    ----------
    ins : iterable
        aqueous_undefined.
    outs : iterable
        aqueous_defined.
        
    References
    ----------
    [1] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
        https://doi.org/10.1021/acs.est.8b04035.
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
        
    def _run(self):
        aqueous_undefined = self.ins[0]
        aqueous_defined = self.outs[0]
        
        aqueous_defined.copy_like(aqueous_undefined)
        aqueous_defined.empty()
        
        aqueous_defined.imass['C'] = self.ins[0]._source.HTLaqueous_C
        aqueous_defined.imass['N'] = self.ins[0]._source.HTLaqueous_N
        aqueous_defined.imass['P'] = self.ins[0]._source.HTLaqueous_P
        # other compositions represented by H2O except C, N, P
        aqueous_defined.imass['H2O'] = aqueous_undefined.F_mass -\
                                            aqueous_defined.imass['C'] -\
                                            aqueous_defined.imass['N'] -\
                                            aqueous_defined.imass['P']
    
    @property
    def pH(self):
        return 7

# =============================================================================
# CatalyticHydrothermalGasification
# =============================================================================

# based on [1]
@cost(ID='Guard bed and CHG reactor 1', basis='Half wet mass flowrate', units='lb/h',
      cost=2041875*1.1, S=76235, CE=qs.CEPCI_by_year[2011], n=0.65, BM=2)
@cost(ID='Guard bed and CHG reactor 2', basis='Half wet mass flowrate', units='lb/h',
      cost=2041875*1.1, S=76235, CE=qs.CEPCI_by_year[2011], n=0.65, BM=2)
@cost(ID='Hydrocyclone', basis='Wet mass flowrate', units='lb/h',
      cost=5000000, S=968859, CE=qs.CEPCI_by_year[2009], n=0.65, BM=2.1)
class CatalyticHydrothermalGasification(SanUnit):
    '''
    CHG serves to reduce the COD content in the aqueous phase and produce fuel
    gas under elevated temperature (350°C) and pressure. The outlet will be
    cooled down and separated by a flash unit.
    
    scope 1 emission: N/A
    scope 2 emission: electricity
    scope 3 emission: 7.8% Ru/C
    
    Parameters
    ----------
    ins : iterable
        chg_in, catalyst_in.
    outs : iterable
        chg_out, catalyst_out.
    eff_T : float
        HTL effluent temperature, [K].
    eff_P : float
        HTL effluent pressure, [Pa].
    FOAK : bool
        whether the facility is first-of-a-king, [True, False].
    pctnew : float
        percentage of capital cost of commercially undemonstrated equipment, [%].
    impurities : int (float)
        impurities buildup and corrosion issues, ranging from 0 to 5, [-].
    complexity : int (float)
        number of continuously linked process steps, ranging from 1 to 11, [-].
    inclusiveness : float
        percentage of pre-startup personnel, inventory, and land purchase costs, [%].
    project_definition : int (float)
        levels of site-specific information and engineering included in estimate, ranging from 2 to 8, [-].
    newsteps : int (float)
        number of processes incorporating commercially unproven technologies, [-].
    baleqs : float
        percentage of heat and mass balance equations based on actual data from prior plants, [%].
    waste : int (float)
        difficulties with waste handling encountered during development, ranging from 0 to 5, [-].
    solids_handling : bool
        True for solids handling, False for no solids handling, [True, False].
    
    See Also
    --------
    :class:`qsdsan.sanunits.CatalyticHydrothermalGasification`
    
    References
    ----------
    [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    '''
    _N_ins = 2
    _N_outs = 2
    
    _units= {'Half wet mass flowrate':'lb/h',
             'Wet mass flowrate':'lb/h'}
    
    auxiliary_unit_names=('pump','hx','inf_hx','eff_hx')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20, gas_C_2_total_C=0.5981,
                 gas_composition={'CH4': 0.527, 'CO2': 0.432, 'C2H6': 0.011,
                                  'C3H8': 0.030, 'H2': 0.0001},
                 WHSV=3.562, catalyst_lifetime=7920, T=350 + _C_to_K,
                 P=3089.7*_psi_to_Pa, eff_T=60 + _C_to_K,
                 eff_P=3089.7*_psi_to_Pa, FOAK=True, pctnew=50, impurities=4,
                 complexity=6, inclusiveness=33, project_definition=5,
                 newsteps=4, baleqs=0, waste=2, solids_handling=False):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.gas_C_2_total_C = gas_C_2_total_C
        self.gas_composition = gas_composition
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.T = T
        self.P = P
        self.eff_T = eff_T
        self.eff_P = eff_P
        pump_in = Stream(f'{ID}_pump_in')
        pump_out = Stream(f'{ID}_pump_out')
        # purchase costs for Pump are related to bst.CE
        self.pump = Pump(ID=f'.{ID}_pump', ins=pump_in, outs=pump_out, P=P)
        inf_pre_hx = Stream(f'{ID}_inf_pre_hx')
        eff_pre_hx = Stream(f'{ID}_eff_pre_hx')
        inf_after_hx = Stream(f'{ID}_inf_after_hx')
        eff_after_hx = Stream(f'{ID}_eff_after_hx')
        self.hx = HXprocess(ID=f'.{ID}_hx', ins=(inf_pre_hx, eff_pre_hx), outs=(inf_after_hx, eff_after_hx))
        inf_hx_out = Stream(f'{ID}_inf_hx_out')
        self.inf_hx = HXutility(ID=f'.{ID}_inf_hx', ins=inf_after_hx, outs=inf_hx_out, T=T, rigorous=True)
        self._inf_at_temp = Stream(f'{ID}_inf_at_temp')
        self._eff_at_temp = Stream(f'{ID}_eff_at_temp')
        eff_hx_out = Stream(f'{ID}_eff_hx_out')
        self.eff_hx = HXutility(ID=f'.{ID}_eff_hx', ins=eff_after_hx, outs=eff_hx_out, T=eff_T, rigorous=True)
        self.FOAK = FOAK
        self.pctnew = pctnew
        self.impurities = impurities
        self.complexity = complexity
        self.inclusiveness = inclusiveness
        self.project_definition = project_definition
        self.newsteps = newsteps
        self.baleqs = baleqs
        self.waste = waste
        self.solids_handling = solids_handling
        if self.FOAK:
            self.plant_performance_factor = 85.77 - 9.69*self.newsteps + 0.33*self.baleqs -\
                                            4.12*self.waste - 17.91*self.solids_handling
        else:
            self.plant_performance_factor = 100
    
    def _run(self):
        chg_in, catalyst_in = self.ins
        chg_out, catalyst_out = self.outs
        
        chg_out.phase = 'l'
        
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts
        catalyst_in.imass['CHG_catalyst'] = chg_in.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        
        cmps = self.components
        gas_C_ratio = 0
        for name, ratio in self.gas_composition.items():
            gas_C_ratio += ratio*cmps[name].i_C
        
        gas_mass_flow = chg_in.imass['C']*self.gas_C_2_total_C/gas_C_ratio
        
        for name,ratio in self.gas_composition.items():
            chg_out.imass[name] = gas_mass_flow*ratio
        
        # all C, N, and P are accounted in H2O here, but will be calculated as properties
        chg_out.imass['H2O'] = chg_in.F_mass - gas_mass_flow
        
        for i in self.outs:
            i.T = self.T
            i.P = self.P
        
        self._eff_at_temp.mix_from([chg_out,], vle=True)
        
        for attr, val in zip(('T','P'), (self.eff_T, self.eff_P)):
            if val:
                for i in self.outs: setattr(i, attr, val)
    
    def _design(self):
        D = self.design_results
        D['Half wet mass flowrate'] = self.ins[0].F_mass*_kg_to_lb/2
        D['Wet mass flowrate'] = self.ins[0].F_mass*_kg_to_lb
        
        pump = self.pump
        pump.ins[0].copy_like(self.ins[0])
        pump.simulate()
        
        hx = self.hx
        inf_hx = self.inf_hx
        inf_hx_in, inf_hx_out = inf_hx.ins[0], inf_hx.outs[0]
        inf_pre_hx, eff_pre_hx = hx.ins
        inf_after_hx, eff_after_hx = hx.outs
        inf_pre_hx.copy_like(self.ins[0])
        eff_pre_hx.copy_like(self._eff_at_temp)
        
        # use product to heat up influent
        hx.T_lim1 = self.eff_T
        hx.simulate()
        for i in self.outs:
            i.T = eff_after_hx.T
        
        # additional inf HX
        inf_hx_in.copy_like(inf_after_hx)
        inf_hx_out.copy_flow(inf_hx_in)
        inf_hx_out.T = self.T
        inf_hx.simulate_as_auxiliary_exchanger(ins=inf_hx.ins, outs=inf_hx.outs)
        
        # additional eff HX
        eff_hx = self.eff_hx
        eff_hx_in, eff_hx_out = eff_hx.ins[0], eff_hx.outs[0]
        eff_hx_in.copy_like(eff_after_hx)
        eff_hx_out.mix_from(self.outs)
        eff_hx_out.T = self.eff_T
        eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs)
        
        # original code from exposan.saf._units, cooling duty very high
        # duty = self.Hnet + eff_hx.Hnet
        # eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs, duty=duty)
        
        for i in self.outs:
            i.T = self.eff_T
    
    def _cost(self):
        self._decorated_cost()
        
        if self.FOAK:
            self.cost_growth_factor = 1.12196 - 0.00297*self.pctnew -\
                                      0.02125*self.impurities - 0.01137*self.complexity +\
                                      0.00111*self.inclusiveness - 0.06361*self.project_definition
        else:
            self.cost_growth_factor = 1
        
        for C in [self.baseline_purchase_costs] + [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]/self.cost_growth_factor
    
    @property
    def CHGout_C(self):
        # not include carbon in the gas phase
        return self.ins[0].imass['C']*(1 - self.gas_C_2_total_C)
    
    @property
    def CHGout_N(self):
        return self.ins[0].imass['N']
    
    @property
    def CHGout_P(self):
        return self.ins[0].imass['P']

# =============================================================================
# Pyrolysis
# =============================================================================

# the cost here include heat drying (in addition to the heat drying unit in the system), HX, and other auxiliary units
# n=0.7 and BM=2 to be similar to HTL-based systems
@cost(ID='Pyrolysis system 1', basis='Half wet mass flowrate', units='tonne/day',
      cost=4011180, S=100, CE=qs.CEPCI_by_year[2021], n=0.7, BM=2)
@cost(ID='Pyrolysis system 2', basis='Half wet mass flowrate', units='tonne/day',
      cost=4011180, S=100, CE=qs.CEPCI_by_year[2021], n=0.7, BM=2)
class Pyrolysis(SanUnit):
    '''
    Pyrolysis of dried sludge or biosolids.
    
    Thermochemical units are sealed and fugitive emissions for pyrolysis are
    tested to contribute to <1% of stream life cycle impacts. Therefore, assume
    no fugitive CH4 and N2O emissions, which is also consistent with other
    thermochemical units (note the BEAM model includes fugitive emissions).
    
    Capital costs are based on [1].
    
    Biooil and pyrogas yields are based on [2].
    
    Note diesel is used for land application and is therefore not listed here,
    though it is included in the model
    scope 1 emission: N/A
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        dried_solids, diesel.
    outs : iterable
        biooil, biochar, pyrogas.
    biooil_yield : float
        biooil yield on a dry weight basis, [-].
    pyrogas_yield : float
        pyrogas yield on a dry weight basis, [-].
    pyrogas_composition : dict
        pyrogas composition, [-].
    unit_electricity : float
        electricity for pyrolysis, [kWh/dry tonne solids].
    load_size : float
        size of loads per truck, [m3/load].
    load_frequency : float
        load frequency, [load/h].
    tractor_fuel : float
        tractor fuel, [L diesel/h].
    crude_oil_density : float
        density of crude oil, [kg/m3].
    crude_oil_HHV : float
        HHV of crude oil, [MJ/kg].
    biooil_density : float
        density of biooil, [kg/m3].
    biooil_HHV : float
        HHV of biooil, [MJ/kg].
    biooil_distance : float
        distance between WRRFs and oil refineries, [km].
    biochar_distance : float
        distance between WRRFs and land application sites, [km].
    FOAK : bool
        whether the facility is first-of-a-king, [True, False].
    pctnew : float
        percentage of capital cost of commercially undemonstrated equipment, [%].
    impurities : int (float)
        impurities buildup and corrosion issues, ranging from 0 to 5, [-].
    complexity : int (float)
        number of continuously linked process steps, ranging from 1 to 11, [-].
    inclusiveness : float
        percentage of pre-startup personnel, inventory, and land purchase costs, [%].
    project_definition : int (float)
        levels of site-specific information and engineering included in estimate, ranging from 2 to 8, [-].
    newsteps : int (float)
        number of processes incorporating commercially unproven technologies, [-].
    baleqs : float
        percentage of heat and mass balance equations based on actual data from prior plants, [%].
    waste : int (float)
        difficulties with waste handling encountered during development, ranging from 0 to 5, [-].
    solids_handling : bool
        True for solids handling, False for no solids handling, [True, False].
    
    References
    ----------
    [1] Shoaib Ahmed Khan, M.; Grioui, N.; Halouani, K.; Benelmir, R. Techno-Economic
        Analysis of Production of Bio-Oil from Catalytic Pyrolysis of Olive Mill
        Wastewater Sludge with Two Different Cooling Mechanisms.
        Energy Conversion and Management: X 2022, 13, 100170.
        https://doi.org/10.1016/j.ecmx.2021.100170.
    [2] Pelagalli, V.; Langone, M.; Matassa, S.; Race, M.; Tuffi, R.; Papirio, S.;
        Lens, P. N. L.; Lazzazzara, M.; Frugis, A.; Petta, L.; Esposito, G. Pyrolysis
        of Municipal Sewage Sludge: Challenges, Opportunities and New Valorization
        Routes for Biochar, Bio-Oil, and Pyrolysis Gas. Environ. Sci.: Water Res. Technol.
        2024, 10 (10), 2282–2312. https://doi.org/10.1039/D4EW00278D.
    '''
    _N_ins = 2
    _N_outs = 3
    
    _units= {'Half wet mass flowrate':'tonne/day'}
    
    auxiliary_unit_names = ('hx','inf_hx','eff_hx')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20, biooil_yield=0.37,
                 pyrogas_yield=0.2,
                 pyrogas_composition={'H2': 0.019, 'CH4': 0.073, 'C2H6': 0.049,
                                      'C3H8': 0.072, 'CO': 0.238, 'CO2': 0.549},
                 unit_electricity=123.424, T=600 + _C_to_K, eff_T=60 + _C_to_K,
                 load_size=13, load_frequency=3, tractor_fuel=25,
                 # https://www.transmountain.com/about-petroleum-liquids (accessed 2025-02-05)
                 crude_oil_density=850,
                 # crude oil HHV: https://world-nuclear.org/information-library/facts-and-figures/heat-values-of-various-fuels
                 crude_oil_HHV=44.5,
                 # average, https://doi.org/10.1039/D4EW00278D
                 biooil_HHV=33.7,
                 # average, https://doi.org/10.1039/D4EW00278D
                 biooil_density=1072.5, biooil_distance=100,
                 biochar_distance=100, FOAK=True, pctnew=19, impurities=4,
                 complexity=6, inclusiveness=33, project_definition=3,
                 newsteps=2, baleqs=0, waste=2, solids_handling=True):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.biooil_yield = biooil_yield
        self.pyrogas_yield = pyrogas_yield
        self.pyrogas_composition = pyrogas_composition
        self.unit_electricity = unit_electricity
        self.T = T
        self.eff_T = eff_T
        inf_pre_hx = Stream(f'{ID}_inf_pre_hx')
        eff_pre_hx = Stream(f'{ID}_eff_pre_hx')
        inf_after_hx = Stream(f'{ID}_inf_after_hx')
        eff_after_hx = Stream(f'{ID}_eff_after_hx')
        self.hx = HXprocess(ID=f'.{ID}_hx', ins=(inf_pre_hx, eff_pre_hx), outs=(inf_after_hx, eff_after_hx))
        inf_hx_out = Stream(f'{ID}_inf_hx_out')
        self.inf_hx = HXutility(ID=f'.{ID}_inf_hx', ins=inf_after_hx, outs=inf_hx_out, T=T, rigorous=True)
        self._inf_at_temp = Stream(f'{ID}_inf_at_temp')
        self._eff_at_temp = Stream(f'{ID}_eff_at_temp')
        eff_hx_out = Stream(f'{ID}_eff_hx_out')
        self.eff_hx = HXutility(ID=f'.{ID}_eff_hx', ins=eff_after_hx, outs=eff_hx_out, T=eff_T, rigorous=True)
        self.load_size = load_size
        self.load_frequency = load_frequency
        self.tractor_fuel = tractor_fuel
        self.crude_oil_density = crude_oil_density
        self.crude_oil_HHV = crude_oil_HHV
        self.biooil_HHV = biooil_HHV
        self.biooil_density = biooil_density
        self.biooil_distance = biooil_distance
        self.biochar_distance = biochar_distance
        self.FOAK = FOAK
        self.pctnew = pctnew
        self.impurities = impurities
        self.complexity = complexity
        self.inclusiveness = inclusiveness
        self.project_definition = project_definition
        self.newsteps = newsteps
        self.baleqs = baleqs
        self.waste = waste
        self.solids_handling = solids_handling
        if self.FOAK:
            self.plant_performance_factor = 85.77 - 9.69*self.newsteps + 0.33*self.baleqs -\
                                            4.12*self.waste - 17.91*self.solids_handling
        else:
            self.plant_performance_factor = 100
    
    def _run(self):
        dried_solids, diesel = self.ins
        biooil, biochar, pyrogas = self.outs
        
        biooil.phase = 'l'
        biochar.phase = 's'
        pyrogas.phase = 'g'
        
        # this just includes the land application onsite spreading and does not include the transporation from WRRFs to land application sites
        diesel.imass['Diesel'] = (self.ins[0].F_vol/self.load_size/self.load_frequency*self.tractor_fuel)/1000*diesel_density
        
        # dry tonne/h
        self.dry_solids = (dried_solids.F_mass - dried_solids.imass['H2O'])/1000
        
        biooil.imass['Biooil'] = self.dry_solids*1000*self.biooil_yield
        
        gas_mass_flow = self.dry_solids*1000*self.pyrogas_yield
        
        pyrogas.imass['H2O'] = dried_solids.imass['H2O']
        
        for name, ratio in self.pyrogas_composition.items():
            pyrogas.imass[name] = gas_mass_flow*ratio
        
        biochar.imass['Biochar'] = self.dry_solids*1000*(1 - self.biooil_yield - self.pyrogas_yield)
        
        for i in self.outs:
            i.T = self.T
        
        self._eff_at_temp.mix_from([biooil, pyrogas], vle=True)
        
        for i in self.outs:
            i.T = self.eff_T
    
    def _design(self):
        D = self.design_results
        D['Half wet mass flowrate'] = self.ins[0].F_mass/1000*24/2
        
        hx = self.hx
        inf_hx = self.inf_hx
        inf_hx_in, inf_hx_out = inf_hx.ins[0], inf_hx.outs[0]
        inf_pre_hx, eff_pre_hx = hx.ins
        inf_after_hx, eff_after_hx = hx.outs
        inf_pre_hx.copy_like(self.ins[0])
        eff_pre_hx.copy_like(self._eff_at_temp)
        
        # use product to heat up influent
        hx.T_lim1 = self.eff_T
        hx.simulate()
        for i in self.outs:
            i.T = eff_after_hx.T
        
        # additional inf HX
        inf_hx_in.copy_like(inf_after_hx)
        inf_hx_out.copy_flow(inf_hx_in)
        inf_hx_out.T = self.T
        inf_hx.simulate_as_auxiliary_exchanger(ins=inf_hx.ins, outs=inf_hx.outs)
        
        # additional eff HX
        eff_hx = self.eff_hx
        eff_hx_in, eff_hx_out = eff_hx.ins[0], eff_hx.outs[0]
        eff_hx_in.copy_like(eff_after_hx)
        eff_hx_out.mix_from(self.outs)
        eff_hx_out.T = self.eff_T
        eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs)
        
        for i in self.outs:
            i.T = self.eff_T
        
        # kW
        self.add_power_utility(self.dry_solids*self.unit_electricity)
    
    def _cost(self):
        self.hx.baseline_purchase_costs.clear()
        self.inf_hx.baseline_purchase_costs.clear()
        self.eff_hx.baseline_purchase_costs.clear()
        self._decorated_cost()
        
        if self.FOAK:
            self.cost_growth_factor = 1.12196 - 0.00297*self.pctnew -\
                                      0.02125*self.impurities - 0.01137*self.complexity +\
                                      0.00111*self.inclusiveness - 0.04011*self.project_definition
        else:
            self.cost_growth_factor = 1
        
        # no auxiliary unit
        for C in [self.baseline_purchase_costs] + [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]/self.cost_growth_factor
    
    @property
    def moisture(self):
        return self.biochar_moisture_content
    
    @property
    def unit_process(self):
        return 'pyrolysis'

# =============================================================================
# Gasification
# =============================================================================

# the cost here include heat drying (in addition to the heat drying unit in the system), HX, and other auxiliary units
# n=0.7 and BM=2 to be similar to HTL-based systems
@cost(ID='Gasification system 1', basis='Half dry mass flowrate', units='tonne/day',
      cost=942500, S=5, CE=qs.CEPCI_by_year[2010], n=0.7, BM=2)
@cost(ID='Gasification system 2', basis='Half dry mass flowrate', units='tonne/day',
      cost=942500, S=5, CE=qs.CEPCI_by_year[2010], n=0.7, BM=2)
class Gasification(SanUnit):
    '''
    Gasification of dried sludge or biosolids with air as the process flow.
    
    Capital costs are based on [1].
    
    Note diesel is used for land application and is therefore not listed here,
    though it is included in the model
    scope 1 emission: N/A
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        dried_solids.
    outs : iterable
        tar, ash, syngas.
    VS_to_tar : float
        tar yield from volatile solids, [-].
        16 datapoints, min: 0.00396, average: 0.0190, max: 0.0319, [2]
    syngas_composition : dict
        syngas composition, [-].
    unit_electricity : float
        electricity for gasification, [kWh/dry tonne solids].
        based on [3], average of high-termperature (34.8 MW) and low-termperature
        (24.4 MW) gasification, 2000 dry tonne/day, 7446 hours/year.
    load_size : float
        size of loads per truck, [m3/load].
    load_frequency : float
        load frequency, [load/h].
    tractor_fuel : float
        tractor fuel, [L diesel/h].
    crude_oil_density : float
        density of crude oil, [kg/m3].
    crude_oil_HHV : float
        HHV of crude oil, [MJ/kg].
    tar_density : float
        density of tar, [kg/m3].
    tar_HHV : float
        HHV of tar, [MJ/kg].
    tar_distance : float
        distance between WRRFs and oil refineries, [km].
    biochar_distance : float
        distance between WRRFs and land application sites, [km].
    FOAK : bool
        whether the facility is first-of-a-king, [True, False].
    pctnew : float
        percentage of capital cost of commercially undemonstrated equipment, [%].
    impurities : int (float)
        impurities buildup and corrosion issues, ranging from 0 to 5, [-].
    complexity : int (float)
        number of continuously linked process steps, ranging from 1 to 11, [-].
    inclusiveness : float
        percentage of pre-startup personnel, inventory, and land purchase costs, [%].
    project_definition : int (float)
        levels of site-specific information and engineering included in estimate, ranging from 2 to 8, [-].
    newsteps : int (float)
        number of processes incorporating commercially unproven technologies, [-].
    baleqs : float
        percentage of heat and mass balance equations based on actual data from prior plants, [%].
    waste : int (float)
        difficulties with waste handling encountered during development, ranging from 0 to 5, [-].
    solids_handling : bool
        True for solids handling, False for no solids handling, [True, False].
    
    References
    ----------
    [1] Clack, K.; Rajagopal, D.; Hoek, E. M. V. Life Cycle and Techno-Economic
        Assessment of Bioresource Production from Wastewater. npj Clean Water
        2024, 7 (1), 1–17. https://doi.org/10.1038/s41545-024-00314-9.
    [2] Adegoroye, A.; Paterson, N.; Li, X.; Morgan, T.; Herod, A. A.; Dugwell, D. R.;
        Kandiyoti, R. The Characterisation of Tars Produced during the Gasification of
        Sewage Sludge in a Spouted Bed Reactor. Fuel 2004, 83 (14), 1949–1960.
        https://doi.org/10.1016/j.fuel.2004.04.006.
    [3] Swanson, R. M.; Platon, A.; Satrio, J. A.; Brown, R. C.; Hsu, D. D.
        Techno-Economic Analysis of Biofuels Production Based on Gasification; 2010.
    '''
    _N_ins = 1
    _N_outs = 3
    
    _units= {'Half dry mass flowrate':'tonne/day'}
    
    auxiliary_unit_names = ('hx','inf_hx','eff_hx')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20, VS_to_tar=0.0190,
                 syngas_composition={'H2': 0.04, 'CH4': 0.02, 'CO': 0.58, 'CO2': 0.36},
                 unit_electricity=301.92, T=900 + _C_to_K, eff_T=60 + _C_to_K,
                 FOAK=True, pctnew=19, impurities=4, complexity=6,
                 inclusiveness=33, project_definition=3, newsteps=2, baleqs=0,
                 waste=2, solids_handling=True):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.VS_to_tar = VS_to_tar
        self.syngas_composition = syngas_composition
        self.unit_electricity = unit_electricity
        self.T = T
        self.eff_T = eff_T
        inf_pre_hx = Stream(f'{ID}_inf_pre_hx')
        eff_pre_hx = Stream(f'{ID}_eff_pre_hx')
        inf_after_hx = Stream(f'{ID}_inf_after_hx')
        eff_after_hx = Stream(f'{ID}_eff_after_hx')
        self.hx = HXprocess(ID=f'.{ID}_hx', ins=(inf_pre_hx, eff_pre_hx), outs=(inf_after_hx, eff_after_hx))
        inf_hx_out = Stream(f'{ID}_inf_hx_out')
        self.inf_hx = HXutility(ID=f'.{ID}_inf_hx', ins=inf_after_hx, outs=inf_hx_out, T=T, rigorous=True)
        self._inf_at_temp = Stream(f'{ID}_inf_at_temp')
        self._eff_at_temp = Stream(f'{ID}_eff_at_temp')
        eff_hx_out = Stream(f'{ID}_eff_hx_out')
        self.eff_hx = HXutility(ID=f'.{ID}_eff_hx', ins=eff_after_hx, outs=eff_hx_out, T=eff_T, rigorous=True)
        self.FOAK = FOAK
        self.pctnew = pctnew
        self.impurities = impurities
        self.complexity = complexity
        self.inclusiveness = inclusiveness
        self.project_definition = project_definition
        self.newsteps = newsteps
        self.baleqs = baleqs
        self.waste = waste
        self.solids_handling = solids_handling
        if self.FOAK:
            self.plant_performance_factor = 85.77 - 9.69*self.newsteps + 0.33*self.baleqs -\
                                            4.12*self.waste - 17.91*self.solids_handling
        else:
            self.plant_performance_factor = 100
    
    def _run(self):
        dried_solids = self.ins[0]
        tar, ash, syngas = self.outs
        
        tar.phase = 'l'
        ash.phase = 's'
        syngas.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (dried_solids.F_mass - dried_solids.imass['H2O'])/1000
        
        ash.imass['Sludge_ash'] = dried_solids.imass['Sludge_ash']
        
        VS = dried_solids.F_mass - dried_solids.imass['H2O'] - dried_solids.imass['Sludge_ash']
        
        tar.imass['Tar'] = VS*self.VS_to_tar
        
        gas_mass_flow = VS*(1 - self.VS_to_tar)
        
        syngas.imass['H2O'] = dried_solids.imass['H2O']
        
        for name, ratio in self.syngas_composition.items():
            syngas.imass[name] = gas_mass_flow*ratio
        
        for i in self.outs:
            i.T = self.T
        
        self._eff_at_temp.mix_from([tar, syngas], vle=True)
        
        for i in self.outs:
            i.T = self.eff_T
    
    def _design(self):
        D = self.design_results
        D['Half dry mass flowrate'] = (self.ins[0].F_mass - self.ins[0].imass['H2O'])/1000*24/2
        
        hx = self.hx
        inf_hx = self.inf_hx
        inf_hx_in, inf_hx_out = inf_hx.ins[0], inf_hx.outs[0]
        inf_pre_hx, eff_pre_hx = hx.ins
        inf_after_hx, eff_after_hx = hx.outs
        inf_pre_hx.copy_like(self.ins[0])
        eff_pre_hx.copy_like(self._eff_at_temp)
        
        # use product to heat up influent
        hx.T_lim1 = self.eff_T
        hx.simulate()
        for i in self.outs:
            i.T = eff_after_hx.T
        
        # additional inf HX
        inf_hx_in.copy_like(inf_after_hx)
        inf_hx_out.copy_flow(inf_hx_in)
        inf_hx_out.T = self.T
        inf_hx.simulate_as_auxiliary_exchanger(ins=inf_hx.ins, outs=inf_hx.outs)
        
        # additional eff HX
        eff_hx = self.eff_hx
        eff_hx_in, eff_hx_out = eff_hx.ins[0], eff_hx.outs[0]
        eff_hx_in.copy_like(eff_after_hx)
        eff_hx_out.mix_from(self.outs)
        eff_hx_out.T = self.eff_T
        eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs)
        
        for i in self.outs:
            i.T = self.eff_T
        
        # kW
        self.add_power_utility(self.dry_solids*self.unit_electricity)
    
    def _cost(self):
        self.hx.baseline_purchase_costs.clear()
        self.inf_hx.baseline_purchase_costs.clear()
        self.eff_hx.baseline_purchase_costs.clear()
        self._decorated_cost()
        
        if self.FOAK:
            self.cost_growth_factor = 1.12196 - 0.00297*self.pctnew -\
                                      0.02125*self.impurities - 0.01137*self.complexity +\
                                      0.00111*self.inclusiveness - 0.04011*self.project_definition
        else:
            self.cost_growth_factor = 1
        
        # no auxiliary unit
        for C in [self.baseline_purchase_costs] + [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]/self.cost_growth_factor
    
    @property
    def moisture(self):
        return self.biochar_moisture_content
    
    @property
    def unit_process(self):
        return 'gasification'

# =============================================================================
# SupercriticalWaterOxidation
# =============================================================================

# the cost here include HX and other auxiliary units
# n=0.7 and BM=2 to be similar to HTL-based systems
@cost(ID='SCWO system 1', basis='Half wet mass flowrate', units='tonne/day',
      cost=1993701, S=24, CE=qs.CEPCI_by_year[2023], n=0.7, BM=2)
@cost(ID='SCWO system 2', basis='Half wet mass flowrate', units='tonne/day',
      cost=1993701, S=24, CE=qs.CEPCI_by_year[2023], n=0.7, BM=2)
class SupercriticalWaterOxidation(SanUnit):
    '''
    SCWO converts feedstock to ash miniral and gas (CO2 + H2O) under elevated
    temperature and pressure.
    
    Capital costs are based on [1].
    
    scope 1 emission: N/A
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        dewatered_solids.
    outs : iterable
        ash, offgas.
    lipid_2_C : float
        lipid to C ratio, [-].
    protein_2_C : float
        protein to C ratio, [-].
    carbo_2_C : float
        carbohydrate to C ratio, [-].
    T : float
        SCWO reactor termperature, [K].
    P : float
        SCWO reactor pressure, [Pa].
    eff_T : float
        SCWO effluent temperature, [K].
    eff_P : float
        SCWO effluent pressure, [Pa].
    FOAK : bool
        whether the facility is first-of-a-king, [True, False].
    pctnew : float
        percentage of capital cost of commercially undemonstrated equipment, [%].
    impurities : int (float)
        impurities buildup and corrosion issues, ranging from 0 to 5, [-].
    complexity : int (float)
        number of continuously linked process steps, ranging from 1 to 11, [-].
    inclusiveness : float
        percentage of pre-startup personnel, inventory, and land purchase costs, [%].
    project_definition : int (float)
        levels of site-specific information and engineering included in estimate, ranging from 2 to 8, [-].
    newsteps : int (float)
        number of processes incorporating commercially unproven technologies, [-].
    baleqs : float
        percentage of heat and mass balance equations based on actual data from prior plants, [%].
    waste : int (float)
        difficulties with waste handling encountered during development, ranging from 0 to 5, [-].
    solids_handling : bool
        True for solids handling, False for no solids handling, [True, False].
    
    References
    ----------
    [1] Qiu, Y.; Zhang, F.; Yuan, Y.; Zhao, Y.; Liu, Y.; Rong, W. Thermodynamic
        and Economic Comparisons of Supercritical Water Oxidation and Gasification
        of Oily Sludge under Hydrothermal Flames. International Journal of Hydrogen
        Energy 2024, 85, 571–585. https://doi.org/10.1016/j.ijhydene.2024.08.362.
    [2] Feng, S.; Guanghua, L. Chapter 4 - Hydrothermal and Solvothermal
        Syntheses. In Modern Inorganic Synthetic Chemistry; Xu, R., Pang, W.,
        Huo, Q., Eds.; Elsevier: Amsterdam, 2011; pp 63–95.
        https://doi.org/10.1016/B978-0-444-53599-3.10004-6.
    '''
    _N_ins = 1
    _N_outs = 2
    
    _units= {'Half wet mass flowrate':'tonne/day'}
    
    auxiliary_unit_names = ('pump','hx','inf_hx','eff_hx','kodrum')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=15, lipid_2_C=0.750,
                 protein_2_C=0.545, carbo_2_C=0.400,
                 # T and P, [2]
                 T=550 + _C_to_K, P=3e7,
                 # assume to be the same as HTL
                 eff_T=60 + _C_to_K, eff_P=30*_psi_to_Pa, FOAK=True, pctnew=50,
                 impurities=5, complexity=6, inclusiveness=33,
                 project_definition=5, newsteps=4, baleqs=0, waste=2,
                 solids_handling=True):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.lipid_2_C = lipid_2_C
        self.protein_2_C = protein_2_C
        self.carbo_2_C = carbo_2_C
        self.T = T
        self.P = P
        self.eff_T = eff_T
        self.eff_P = eff_P
        pump_in = Stream(f'{ID}_pump_in')
        pump_out = Stream(f'{ID}_pump_out')
        self.pump = Pump(ID=f'.{ID}_pump', ins=pump_in, outs=pump_out, P=P)
        inf_pre_hx = Stream(f'{ID}_inf_pre_hx')
        eff_pre_hx = Stream(f'{ID}_eff_pre_hx')
        inf_after_hx = Stream(f'{ID}_inf_after_hx')
        eff_after_hx = Stream(f'{ID}_eff_after_hx')
        self.hx = HXprocess(ID=f'.{ID}_hx', ins=(inf_pre_hx, eff_pre_hx), outs=(inf_after_hx, eff_after_hx))
        inf_hx_out = Stream(f'{ID}_inf_hx_out')
        self.inf_hx = HXutility(ID=f'.{ID}_inf_hx', ins=inf_after_hx, outs=inf_hx_out, T=T, rigorous=True)
        self._inf_at_temp = Stream(f'{ID}_inf_at_temp')
        self._eff_at_temp = Stream(f'{ID}_eff_at_temp')
        eff_hx_out = Stream(f'{ID}_eff_hx_out')
        self.eff_hx = HXutility(ID=f'.{ID}_eff_hx', ins=eff_after_hx, outs=eff_hx_out, T=eff_T, rigorous=True)
        self.FOAK = FOAK
        self.pctnew = pctnew
        self.impurities = impurities
        self.complexity = complexity
        self.inclusiveness = inclusiveness
        self.project_definition = project_definition
        self.newsteps = newsteps
        self.baleqs = baleqs
        self.waste = waste
        self.solids_handling = solids_handling
        if self.FOAK:
            self.plant_performance_factor = 85.77 - 9.69*self.newsteps + 0.33*self.baleqs -\
                                            4.12*self.waste - 17.91*self.solids_handling
        else:
            self.plant_performance_factor = 100
    
    def _run(self):
        dewatered_solids = self.ins[0]
        ash, offgas = self.outs
        
        ash.phase = 's'
        offgas.phase = 'g'
        
        ash.imass['Sludge_ash'] = dewatered_solids.imass['Sludge_ash']
        
        # kg OC/h
        input_solids_OC_mass_flow = dewatered_solids.imass['Sludge_lipid']*self.lipid_2_C +\
                                    dewatered_solids.imass['Sludge_protein']*self.protein_2_C +\
                                    dewatered_solids.imass['Sludge_carbo']*self.carbo_2_C
        
        offgas.imass['CO2'] = input_solids_OC_mass_flow*_C_to_CO2
        
        offgas.imass['H2O'] = dewatered_solids.F_mass - ash.imass['Sludge_ash'] - offgas.imass['CO2']
        
        for i in self.outs:
            i.T = self.T
            i.P = self.P
        
        self._eff_at_temp.mix_from([offgas,], vle=True)
        
        for attr, val in zip(('T','P'), (self.eff_T, self.eff_P)):
            if val:
                for i in self.outs: setattr(i, attr, val)
    
    def _design(self):
        D = self.design_results
        D['Half wet mass flowrate'] = self.ins[0].F_mass/1000*24/2
        
        pump = self.pump
        pump.ins[0].copy_like(self.ins[0])
        pump.simulate()
        
        hx = self.hx
        inf_hx = self.inf_hx
        inf_hx_in, inf_hx_out = inf_hx.ins[0], inf_hx.outs[0]
        inf_pre_hx, eff_pre_hx = hx.ins
        inf_after_hx, eff_after_hx = hx.outs
        inf_pre_hx.copy_like(self.ins[0])
        eff_pre_hx.copy_like(self._eff_at_temp)
        
        # use product to heat up influent
        hx.T_lim1 = self.eff_T
        hx.simulate()
        for i in self.outs:
            i.T = eff_after_hx.T
        
        # additional inf HX
        inf_hx_in.copy_like(inf_after_hx)
        inf_hx_out.copy_flow(inf_hx_in)
        inf_hx_out.T = self.T
        inf_hx.simulate_as_auxiliary_exchanger(ins=inf_hx.ins, outs=inf_hx.outs)
        
        # additional eff HX
        eff_hx = self.eff_hx
        eff_hx_in, eff_hx_out = eff_hx.ins[0], eff_hx.outs[0]
        eff_hx_in.copy_like(eff_after_hx)
        eff_hx_out.mix_from(self.outs)
        eff_hx_out.T = self.eff_T
        eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs)
        
        for i in self.outs:
            i.T = self.eff_T
        
        # original code from exposan.saf._units, cooling duty very high
        # duty = self.Hnet + eff_hx.Hnet
        # eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs, duty=duty)
    
    def _cost(self):
        self.pump.baseline_purchase_costs.clear()
        self.hx.baseline_purchase_costs.clear()
        self.inf_hx.baseline_purchase_costs.clear()
        self.eff_hx.baseline_purchase_costs.clear()
        self._decorated_cost()
        
        if self.FOAK:
            self.cost_growth_factor = 1.12196 - 0.00297*self.pctnew -\
                                      0.02125*self.impurities - 0.01137*self.complexity +\
                                      0.00111*self.inclusiveness - 0.06361*self.project_definition
        else:
            self.cost_growth_factor = 1
        
        # no auxiliary unit
        for C in [self.baseline_purchase_costs] + [unit.baseline_purchase_costs for unit in self.auxiliary_units]:
            for item in C.keys():
                C[item] = C[item]/self.cost_growth_factor