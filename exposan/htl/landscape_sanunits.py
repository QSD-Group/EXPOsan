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
from qsdsan.sanunits import SludgePump, WWTpump, HXutility, Reactor, HXprocess, Pump
from qsdsan.utils import auom, calculate_excavation_volume
from qsdsan.equipments import Blower, GasPiping
from thermosteam import indexer, equilibrium
from biosteam import Stream
from biosteam.units.decorators import cost
from math import ceil, pi, log, exp

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

# TODO: check the lifetime of each unit / equipment: see Table S1

# TODO: check natural gas CI, whether including upstream

# methane density, kg/m3
methane_density = 0.707

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
_CaCO3_to_C = 12/100

folder = os.path.dirname(__file__)

# numbers without citations are likely from the BEAM*2024 model

# TODO: add pump for effluents of appropriate units

__all__ = (
    'WRRF',
    'Thickening',
    'AerobicDigestion',
    'AnaerobicDigestion',
    'Dewatering',
    'AlkalineStabilization',
    'Composting',
    'HeatDrying',
    'Incineration',
    'Pyrolysis',
    'Gasification',
    'Landfilling',
    'LandApplication',
    'HydrothermalLiquefaction',
    'HydrothermalAlkalineTreatment',
    'HTLmixer',
    'CatalyticHydrothermalGasification',
    'AcidExtraction',
    'StruvitePrecipitation',
    'MembraneDistillation',
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
    operation_hours : float
        Plant yearly operation hour, [hr/yr].
    sludge_wet_density : float
        The density of sludge of 80% moisture content, [kg/m3].
    sludge_distance : float
        Normalized sludge transportation distance, [km].
    wage_adjustment : float
        A coefficient to adjust labor cost.
    
    References
    ----------
    .. [1] Metcalf and Eddy, Incorporated. 1991. Wastewater Engineering:
           Treatment Disposal and Reuse. New York: McGraw-Hill.
    .. [2] Cai, L.; Gao, D.; Chen, T.-B.; Liu, H.-T.; Zheng, G.-D.; Yang, Q.-W.
           Moisture Variation Associated with Water Input and Evaporation during
           Sewage Sludge Bio-Drying. Bioresource Technology 2012, 117, 13–19.
           https://doi.org/10.1016/j.biortech.2012.03.092.
    .. [3] Li, Y.; Leow, S.; Fedders, A. C.; Sharma, B. K.; Guest, J. S.;
           Strathmann, T. J. Quantitative Multiphase Model for Hydrothermal
           Liquefaction of Algal Biomass. Green Chem. 2017, 19 (4), 1163–1174.
           https://doi.org/10.1039/C6GC03294J.
    '''
    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', 
                 ww_2_dry_sludge=1, # [1]
                 sludge_moisture=0.99, sludge_dw_ash=0.257, 
                 sludge_afdw_lipid=0.204, sludge_afdw_protein=0.463, 
                 operation_hours=None,
                 sludge_wet_density=1040, # [2]
                 sludge_distance=100,
                 wage_adjustment=1):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.ww_2_dry_sludge = ww_2_dry_sludge
        self.sludge_moisture = sludge_moisture
        self.sludge_dw_ash = sludge_dw_ash
        self.sludge_afdw_lipid = sludge_afdw_lipid
        self.sludge_afdw_protein = sludge_afdw_protein
        self.operation_hours = operation_hours
        self.sludge_wet_density = sludge_wet_density
        self.sludge_distance = sludge_distance
        self.wage_adjustment = wage_adjustment
    
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
                 init_with='WasteStream', lifetime=25, target_moisture=0.97,
                 polymer_addition=5, unit_electricity=4.9, max_capacity=100):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
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
        # add pump power utility and costs
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
                 VS_reduction=0.475, HRT=14, SRT=40,
                 # TODO: this includes aeration, but there is a separate blower (check electricity)
                 unit_electricity=0.03, depth=10, wall_concrete_unit_cost=24,
                 slab_concrete_unit_cost=13, excavation_unit_cost=0.3,
                 F_BM=default_F_BM):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.O2_requirement = O2_requirement
        # TODO: use different reduction ratios for lipid, protein, and carbohydrate
        self.VS_reduction = VS_reduction
        self.HRT = HRT
        self.SRT = SRT
        self.unit_electricity = unit_electricity
        self.depth = depth
        self.wall_concrete_unit_cost = wall_concrete_unit_cost
        self.slab_concrete_unit_cost = slab_concrete_unit_cost
        self.excavation_unit_cost = excavation_unit_cost
        self.F_BM.update(F_BM)
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
    '''
    _N_ins = 1
    _N_outs = 3
    
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
                 # TODO: check if vertical_mixer_unit_price include installation 
                 vertical_mixer_unit_power=3.7, vertical_mixer_unit_price=10200,
                 # TODO: 0.2 is an arbitrary value, replace it if there is a number with citaions
                 gas_collection_cost_factor=0.2,
                 F_BM=default_F_BM):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        # TODO: use different reduction ratios for lipid, protein, and carbohydrate
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
        self.F_BM.update(F_BM)
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
        fugitive_methane.imass['CH4'] = methane_volume*self.biogas_fugitive_ratio*methane_density
        
        # from flaring
        fugitive_methane.imass['CH4'] += methane_volume*self.biogas_flare_ratio*self.flare_fugitive_ratio*methane_density        
        
        # TODO: consider adding an effluent RNG and add its own cost and CI (may not be necessary)
        # m3/h
        natural_gas_generated = methane_volume*self.biogas_CHP_ratio +\
                                methane_volume*self.biogas_RNG_ratio*(1 - self.RNG_parasitic)
        
        # TODO: should this be biogas (CH4 + CO2)?
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
        C['Vertical mixer'] = self.electricity_requirement/self.vertical_mixer_unit_power*self.vertical_mixer_unit_price
        C['Gas collection system'] = (C['Wall concrete'] + C['Slab concrete'])*self.gas_collection_cost_factor
    
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
                 init_with='WasteStream', lifetime=20, target_moisture=0.8,
                 polymer_addition=5, methane_concentration=6, mathane_loss=0.9,
                 unit_electricity=1.40):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
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
        # add pump power utility and costs
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

@cost(ID='Lime stabilization', basis='Dry solids flow', units='tonne/day',
      cost=311386, S=1, CE=qs.CEPCI_by_year[2004], n=0.5623, BM=1)
class AlkalineStabilization(SanUnit):
    '''
    Alkaline stabilization of dewatered sludge, with the assumption
    that no water evaporizes.
    
    scope 1 emission: natural gas combustion
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
                 init_with='WasteStream', lifetime=25, lime_dose=0.3,
                 natural_gas_emission=15.6, unit_electricity=3.7):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.lime_dose = lime_dose
        self.natural_gas_emission = natural_gas_emission
        self.unit_electricity = unit_electricity
    
    def _run(self):
        dewatered_sludge, lime = self.ins
        stabilized_solids, NG_as_CO2 = self.outs
        
        stabilized_solids.phase = 'l'
        NG_as_CO2.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (dewatered_sludge.F_mass - dewatered_sludge.imass['H2O'])/1000
        
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
    facility (WRRF) to the land application site. The equipment cost is based on [1].
    
    scope 1 emission: diesel, fugitive methane, fugitive nitrous oxide, carbon sequestration
    scope 2 emission: electricity
    scope 3 emission: nitrogen, phosphorus
    
    Parameters
    ----------
    ins : iterable
        input_solids, bulking_agent, diesel.
    outs : iterable
        compost_cost, composted_solids, fugitive_methane, fugitive_nitrous_oxide.
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
    compost_moisture_content : float
        the moisture content of produced composts, [-].
    load_size : float
        size of loads per truck, [m3/load].
    load_frequency : float
        load frequency, [load/h].
    tractor_fuel : float
        tractor fuel, [L diesel/h].
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
    solids_distance : float
        distance between WRRFs and land application sites, [km].
    
    References
    ----------
    .. [1] Stramer, Y.; Brenner, A.; Cohen, S. B.; Oron, G. Selection of a
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
                 compost_solids_ratio=2, compost_moisture_content=0.5, load_size=13,
                 load_frequency=3, tractor_fuel=25, methane_fugitive_ratio=0.0001,
                 nitrous_oxide_fugitive_ratio=0.00076, unit_carbon_sequestration=0.4475,
                 unit_electricity=0, solids_distance=100):
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
        self.compost_solids_ratio = compost_solids_ratio
        self.compost_moisture_content = compost_moisture_content
        self.load_size = load_size
        self.load_frequency = load_frequency
        self.tractor_fuel = tractor_fuel
        self.methane_fugitive_ratio = methane_fugitive_ratio
        self.nitrous_oxide_fugitive_ratio = nitrous_oxide_fugitive_ratio
        self.unit_carbon_sequestration = unit_carbon_sequestration
        self.unit_electricity = unit_electricity
        self.solids_distance = solids_distance
    
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
        
        # this just includes the land application onsite spreading and does not include the transporation from the composting facility (WRRF) to the land application site
        # L diesel/h
        diesel_land_application = (self.ins[0].F_vol*self.compost_solids_ratio/self.load_size/self.load_frequency*self.tractor_fuel)
        
        diesel.imass['Diesel'] = (diesel_grinding + diesel_other_composting + diesel_land_application)/1000*diesel_density
        
        fugitive_methane.imass['CH4'] = input_solids_OC_mass_flow*self.methane_fugitive_ratio
        
        fugitive_nitrous_oxide.imass['N2O'] = input_solids_N_mass_flow*self.nitrous_oxide_fugitive_ratio
        
        # TODO: based on wastewater solids mass flow only (not including sawdust)
        # TODO: any price?
        sequestered_carbon_dioxide.imass['CO2'] = -self.dry_solids*1000*self.unit_carbon_sequestration
        
        # kW
        self.add_power_utility(self.dry_solids*self.unit_electricity)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        # based on [1], 1593.7 is marshall and swift equipment cost index 2017, 751 is the baseline marshall and swift equipment cost index
        C['Composting equipment'] = (1560*self.ins[0].F_mass/1000*24 + 450000)*1593.7/qs.CEPCI_by_year[2017]*qs.CEPCI_by_year[2022]/751
    
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
        return 'composting'

# =============================================================================
# HeatDrying
# =============================================================================

# n: 0.7, assumed
# BM: 3.17, from biosteam/units/heat_exchange.py
@cost(ID='Dryer', basis='Dry solids flow', units='tonne/day',
      cost=220000*80/2/3.17, S=80, CE=qs.CEPCI_by_year[2018], n=0.7, BM=3.17)
class HeatDrying(SanUnit):
    '''
    Heat drying of sludge or biosolids.
    
    scope 1 emission: natural gas combustion
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
    .. [1] Hao, X.; Chen, Q.; van Loosdrecht, M. C. M.; Li, J.; Jiang, H.
           Sustainable Disposal of Excess Sludge: Incineration without
           Anaerobic Digestion. Water Research 2020, 170, 115298.
           https://doi.org/10.1016/j.watres.2019.115298.
    '''
    _N_ins = 2
    _N_outs = 2
    
    _units = {'Dry solids flow':'tonne/day'}
    
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
       
        dried_solids_mass_flow = self.dry_solids/(1 - self.target_moisture)
        
        dried_solids.copy_like(input_sludge)
       
        dried_solids.imass['H2O'] = dried_solids_mass_flow*1000*self.target_moisture
       
        vapor.imass['H2O'] = input_sludge.F_mass - dried_solids.F_mass
        
        # TODO: add cost and CI, CI can be based on MJ converted from ivol
        # use natural gas for heat drying base on the BEAM model
        natural_gas.ivol['CH4'] = vapor.imass['H2O']/1000*self.unit_heat*1000/self.natural_gas_HHV
        
        dried_solids.T = vapor.T = self.T_out
        
    def _design(self):
        self.design_results['Dry solids flow'] = self.dry_solids*24
        
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

# TODO: consider adding ash disposal

# cost base year: 2014 (the original conference paper cited in [1] was published probably on Jan, 2015, so 1 year before that)
# n: 0.7753, [1]
# BM: 2, to be consistent with thermochemical units
@cost(ID='Incinerator', basis='Wet mass flowrate', units='tonne/day',
      cost=2350700/2/2, S=1000/365, CE=qs.CEPCI_by_year[2014], n=0.7753, BM=2)
class Incineration(SanUnit):
    '''
    Incineration of dried sludge.
    
    scope 1 emission: fugitive methane, fugitive nitrous oxide
    scope 2 emission: electricity
    scope 3 emission: phosphorus
    
    from [1]: I = 2.3507×C0.7753, where I is the investment cost in million
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
    solids_distance : float
        distance between WRRFs and landfills, [km].
    
    References
    ----------
    .. [1] Gergel, I. Cost of incineration plant. Waste To Energy International.
           https://wteinternational.com/news/cost-of-incineration-plant/
           (accessed 2025-08-06).
    '''
    _N_ins = 2
    _N_outs = 4
    
    _units= {'Wet mass flowrate':'tonne/day'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20,
                 incineration_temperature=850, ash_moisture_content=0.01,
                 heat_water_removal=4.5, fugitive_methane_incineration=0.0000097,
                 protein_2_N=0.159, N_2_P=0.3927, suzuki_constant_1=161.3,
                 suzuki_constant_2=0.14, suzuki_lowest_temperature=750,
                 urea_catalyst=True, N2O_adjustment_factor_urea=0.2,
                 heat_solids_incineration=12000, additional_fuel_ratio=0,
                 heat_recovery_ratio=0.5, heat_recovery_efficiency=0.8,
                 electricity_recovery_ratio=0,
                 BTU_to_kWh_with_efficiency=0.0000854, net_capacity_factor=0.85,
                 unit_electricity=200, solids_distance=100):
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
        self.solids_distance = solids_distance
    
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
        
        # TODO: or just dispose of the ash, like what happens in the CHP
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
        
        # TODO: add natural gas in this unit and other units as a heating utility if it is used as a heating utility
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
        D = self.design_results
        D['Wet mass flowrate'] = self.ins[0].F_mass/1000*24
        
        # kW
        self.add_power_utility(self.electricity_net)
        
    @property
    def moisture(self):
        return self.ash_moisture_content
    
    # TODO: or just dispose of the ash, like what happens in the CHP
    # TODO: calculate the price and CI of ash based on the price and CI of phosphorus
    # TODO: may add a discount rate
    @property
    def phosphorus_mass_flow(self):
        # kg P/h
        return self.ins[0].imass['Sludge_protein']*self.protein_2_N*self.N_2_P
    
    @property
    def unit_process(self):
        return 'incineration'

# TODO: a brief comparison of thermochemical units: Table 2 in https://pubs.rsc.org/en/content/articlelanding/2024/ew/d4ew00278d/unauth

# =============================================================================
# Pyrolysis
# =============================================================================

# TODO: Aaron (Cusick Group) has an emperical pyrolysis model (check thesis if that is available)

# TODO: add citations
# TODO: yield: https://pubs.rsc.org/en/content/articlelanding/2024/ew/d4ew00278d/unauth

# TODO: carbon sequestration: see Table SI-12 in the page S17 of the SI of https://pubs.acs.org/doi/full/10.1021/acs.est.2c06083

# TODO: add citations from the excel (data_collection.xlsx)

# the cost here include heat drying (in addition to the heat drying unit in the system), HX, and other auxiliary units
# n=0.7 and BM=2 to be similar to HTL-based systems
@cost(ID='Pyrolysis system', basis='Wet mass flowrate', units='tonne/day',
      cost=4011180, S=100, CE=qs.CEPCI_by_year[2021], n=0.7, BM=2)
class Pyrolysis(SanUnit):
    '''
    Pyrolysis of dried sludge or biosolids.
    
    scope 1 emission: fugitive methane, fugitive nitrous oxide
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        dried_solids.
    outs : iterable
        biooil, biochar, pyrogas, fugitive_methane, fugitive_nitrous_oxide.
    biooil_yield : float
        biooil yield on a dry weight basis, [-].
    pyrogas_yield : float
        pyrogas yield on a dry weight basis, [-].
    pyrogas_composition : dict
        pyrogas composition, [-].
    unit_fugitive_methane : float
        fugitive methane during pyrolysis, [g CH4/dry tonne solids].
    unit_fugitive_nitrous_oxide : float
        fugitive nitrous oxide during pyrolysis, [g N2O/dry tonne solids].
    unit_electricity : float
        electricity for pyrolysis, [kWh/dry tonne solids].
    '''
    _N_ins = 1
    _N_outs = 5
    
    _units= {'Wet mass flowrate':'tonne/day'}
    
    auxiliary_unit_names = ('hx','inf_hx','eff_hx')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20, biooil_yield=0.37,
                 pyrogas_yield=0.2,
                 pyrogas_composition={'H2': 0.019, 'CH4': 0.073, 'C2H6': 0.049,
                                      'C3H8': 0.072, 'CO': 0.238, 'CO2': 0.549},
                 unit_fugitive_methane=2.65, unit_fugitive_nitrous_oxide=5.23,
                 unit_electricity=123.424,
                 T=600 + _C_to_K, eff_T=60 + _C_to_K):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.biooil_yield = biooil_yield
        self.pyrogas_yield = pyrogas_yield
        self.pyrogas_composition = pyrogas_composition
        self.unit_fugitive_methane = unit_fugitive_methane
        self.unit_fugitive_nitrous_oxide = unit_fugitive_nitrous_oxide
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
    
    def _run(self):
        dried_solids = self.ins[0]
        # TODO: for other units with fugitive emissions and streams to be combusted, considering writing fugitive emissions like this
        # TODO: or just assume 1% of produced CH4 is fugitive
        biooil, biochar, pyrogas, fugitive_methane, fugitive_nitrous_oxide = self.outs
        
        biooil.phase = 'l'
        biochar.phase = 's'
        pyrogas.phase = 'g'
        fugitive_methane.phase = 'g'
        fugitive_nitrous_oxide.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (dried_solids.F_mass - dried_solids.imass['H2O'])/1000
        
        biooil.imass['Biooil'] = self.dry_solids*1000*self.biooil_yield
        
        gas_mass_flow = self.dry_solids*1000*self.pyrogas_yield
        
        fugitive_methane.imass['CH4'] = self.dry_solids*self.unit_fugitive_methane/1000
        fugitive_nitrous_oxide.imass['N2O'] = self.dry_solids*self.unit_fugitive_nitrous_oxide/1000
        
        gas_mass_flow -= fugitive_methane.F_mass
        gas_mass_flow -= fugitive_nitrous_oxide.F_mass
        
        pyrogas.imass['H2O'] = dried_solids.imass['H2O']
        
        for name, ratio in self.pyrogas_composition.items():
            pyrogas.imass[name] = gas_mass_flow*ratio
        
        biochar.imass['Biochar'] = self.dry_solids*1000*(1 - self.biooil_yield - self.pyrogas_yield)
        
        for i in self.outs:
            i.T = self.T
        
        self._eff_at_temp.mix_from([biooil, pyrogas, fugitive_methane, fugitive_nitrous_oxide], vle=True)
        
        for i in self.outs:
            i.T = self.eff_T
    
    def _design(self):
        D = self.design_results
        D['Wet mass flowrate'] = self.ins[0].F_mass/1000*24
        
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
    
    @property
    def moisture(self):
        return self.biochar_moisture_content
    
    @property
    def unit_process(self):
        return 'pyrolysis'

# =============================================================================
# Gasification
# =============================================================================

# TODO: carbon sequestration: see Table SI-12 in the page S17 of the SI of https://pubs.acs.org/doi/full/10.1021/acs.est.2c06083

# TODO: add citations from the excel (data_collection.xlsx)

# the cost here include heat drying (in addition to the heat drying unit in the system), HX, and other auxiliary units
# n=0.7 and BM=2 to be similar to HTL-based systems
@cost(ID='Gasification system', basis='Dry mass flowrate', units='tonne/day',
      cost=942500, S=5, CE=qs.CEPCI_by_year[2010], n=0.7, BM=2)
class Gasification(SanUnit):
    '''
    Gasification of dried sludge or biosolids with air as the process flow.
    
    scope 1 emission: fugitive methane, fugitive nitrous oxide
    scope 2 emission: electricity
    scope 3 emission: N/A
    
    Parameters
    ----------
    ins : iterable
        dried_solids.
    outs : iterable
        biooil, biochar, pyrogas, fugitive_methane, fugitive_nitrous_oxide.
    VS_to_tar : float
        tar yield from volatile solids, [-].
    syngas_composition : dict
        syngas composition, [-].
    unit_fugitive_methane : float
        fugitive methane during pyrolysis, [g CH4/dry tonne solids].
    unit_fugitive_nitrous_oxide : float
        fugitive nitrous oxide during pyrolysis, [g N2O/dry tonne solids].
    unit_electricity : float
        electricity for pyrolysis, [kWh/dry tonne solids].
    '''
    _N_ins = 1
    _N_outs = 5
    
    _units= {'Dry mass flowrate':'tonne/day'}
    
    auxiliary_unit_names = ('hx','inf_hx','eff_hx')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20,
                 # TODO: if there is data, replace 0.1 with citations
                 VS_to_tar=0.1,
                 syngas_composition={'H2': 0.04, 'CH4': 0.02, 'CO': 0.58, 'CO2': 0.36},
                 # TODO: this is for pyrolysis, gasification may have a higher value, use this data only if there is no other data available
                 unit_fugitive_methane=2.65,
                 # TODO: this is for pyrolysis, gasification may have a higher value, use this data only if there is no other data available
                 unit_fugitive_nitrous_oxide=5.23,
                 # TODO: this is for pyrolysis, gasification may have a higher value
                 unit_electricity=123.424,
                 T=900 + _C_to_K, eff_T=60 + _C_to_K):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.VS_to_tar = VS_to_tar
        self.syngas_composition = syngas_composition
        self.unit_fugitive_methane = unit_fugitive_methane
        self.unit_fugitive_nitrous_oxide = unit_fugitive_nitrous_oxide
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
    
    def _run(self):
        # TODO: but may need a blower?
        # TODO: add an air stream (is using free air, may not need this stream)
        dried_solids = self.ins[0]
        # TODO: for other units with fugitive emissions and streams to be combusted, considering writing fugitive emissions like this
        tar, biochar, syngas, fugitive_methane, fugitive_nitrous_oxide = self.outs
        
        tar.phase = 'l'
        biochar.phase = 's'
        syngas.phase = 'g'
        fugitive_methane.phase = 'g'
        fugitive_nitrous_oxide.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (dried_solids.F_mass - dried_solids.imass['H2O'])/1000
        
        biochar.imass['Biochar'] = dried_solids.imass['Sludge_ash']
        
        VS = dried_solids.F_mass - dried_solids.imass['H2O'] - dried_solids.imass['Sludge_ash']
        tar.imass['Tar'] = VS*self.VS_to_tar
        
        gas_mass_flow = VS*(1 - self.VS_to_tar)
        
        fugitive_methane.imass['CH4'] = self.dry_solids*self.unit_fugitive_methane/1000
        fugitive_nitrous_oxide.imass['N2O'] = self.dry_solids*self.unit_fugitive_nitrous_oxide/1000
        
        gas_mass_flow -= fugitive_methane.F_mass
        gas_mass_flow -= fugitive_nitrous_oxide.F_mass
        
        syngas.imass['H2O'] = dried_solids.imass['H2O']
        
        for name, ratio in self.syngas_composition.items():
            syngas.imass[name] = gas_mass_flow*ratio
        
        for i in self.outs:
            i.T = self.T
        
        self._eff_at_temp.mix_from([tar, syngas, fugitive_methane, fugitive_nitrous_oxide], vle=True)
        
        for i in self.outs:
            i.T = self.eff_T
    
    def _design(self):
        D = self.design_results
        D['Dry mass flowrate'] = (self.ins[0].F_mass - self.ins[0].imass['H2O'])/1000*24
        
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
    
    @property
    def moisture(self):
        return self.biochar_moisture_content
    
    @property
    def unit_process(self):
        return 'gasification'

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
        landfilled_solids, fugitive_methane, fugitive_nitrous_oxide, sequestered_carbon_dioxide.
    lipid_2_C : float
        lipid to C ratio, [-].
    protein_2_C : float
        protein to C ratio, [-].
    carbo_2_C : float
        carbohydrate to C ratio, [-].
    # TODO: what if there is N loss before landfilling, is protein_2_N here still accurate?
    protein_2_N, float
        N to protein ratio, [-].
    # TODO: or this is actually a model correction factor?
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
    tipping_fee : float
        tipping fee, [$·wet tonne-1].
        see exposan/htl/data/landfilling_tipping_fee.xlsx.
        U.S. average: 56.8 $·wet ton-1 ~ 62.6 $·wet tonne-1.
    solids_distance : float
        distance between WRRFs and landfills, [km].
    '''
    _N_ins = 1
    _N_outs = 4
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lipid_2_C=0.750,
                 protein_2_C=0.545, carbo_2_C=0.400, protein_2_N=0.159,
                 methane_landfilling_gas=0.5, MCF=1, OC_decomposition_ratio=0.5,
                 k_decay=0.06, flare_fugitive_ratio=0.01, C_N_cutoff=30,
                 N2O_N_landfilling=0.015, methane_electricity=0.5,
                 BTU_to_kWh_with_efficiency=0.0000854, net_capacity_factor=0.85,
                 tipping_fee=62.6, solids_distance=100):
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
        self.tipping_fee = tipping_fee
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
    
    def _cost(self):
        landfilled_solids = self.outs[0]
        landfilled_solids.price = -self.tipping_fee/1000
    
    @property
    def electricity_kW(self):
        # TODO: confirm no methane is combusted to produce heat
        return self.methane_combustion*self.methane_electricity/methane_density*methane_heat*self.BTU_to_kWh_with_efficiency*self.net_capacity_factor
    
    @property
    def unit_process(self):
        return 'landfilling'

# =============================================================================
# LandApplication
# =============================================================================

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
        biosolids_cost, fugitive_methane, fugitive_nitrous_oxide, carbon_dioxide.
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
                 nitrous_oxide_reduction_intercept=0.1518,
                 unit_fugitive_nitrous_oxide_strogae=0.00043, climate='humid',
                 OC_sequestered_ratio=0.75, unit_carbon_sequestration=0.4475,
                 # TODO: update CaCO3 equivalence (%-dry weight) value (this might be related to lime stablization)
                 calcium_carbonate_ratio=0.1, calcium_carbonate_from_waste=False,
                 calcium_carbonate_replacing_lime=True, solids_distance=100):
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
        
        # this just includes the land application onsite spreading and does not include the transporation from the composting facility (WRRF) to the land application site
        # TODO: add price and CI, note the unit here is actually L diesel/h, if data are for kg, convert the unit by considering the density of diesel
        diesel.imass['Diesel'] = (self.ins[0].F_vol/self.load_size/self.load_frequency*self.tractor_fuel)/1000*diesel_density
        
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

# =============================================================================
# thermochemical units
# =============================================================================

# TODO: add notes - KnockOutDrum and HydrothermalLiquefaction are based on exposan.saf._units and the original HTL model
# TODO: see if the Transportation unit can be used
# TODO: to be removed
# __all__ = (
#     'HydrothermalLiquefaction',
#     'Hydroprocessing',
#     'PressureSwingAdsorption',
#     'Electrochemical',
#     'HydrogenCenter',
#     'ProcessWaterCenter',
#     'BiocrudeSplitter',
#     'Conditioning',
#     'Transportation',
#     )

# =============================================================================
# HydrothermalLiquefaction
# =============================================================================

# TODO: need to decide whether to recover N and P or not (to the same complexity level as other thermochemical units)

# TODO: the cost here is based on the nth plant assumption, change to the 1st plant and use the learning curve, same for other (especially thermochemical) units / systems

# TODO: these comments apply to HALT as well
# TODO: check this note
# original [1] is 1339 dry-ash free ton per day (tpd), ash content is 13%,
# which is 58176 wet lb/h, scaling basis in the equipment table in [1]
# (flow for S100) does not seem right
# TODO: check data with [1]
# TODO: try to add N=2 for the HTL system, based on [1]
# based on [1]
@cost(ID='HTL system', basis='Wet mass flowrate', units='lb/h',
      cost=18743378, S=306198, CE=qs.CEPCI_by_year[2011], n=0.77, BM=2.1)
@cost(ID='Solids filter oil/water separator', basis='Wet mass flowrate', units='lb/h',
      cost=3945523, S=1219765, CE=qs.CEPCI_by_year[2011], n=0.68, BM=1.9)
@cost(ID='Hot oil system', basis='Wet mass flowrate', units='lb/h',
      cost=4670532, S=306198, CE=qs.CEPCI_by_year[2011], n=0.6, BM=1.4)
class HydrothermalLiquefaction(SanUnit):
    '''
    HTL converts feedstock to gas, aqueous, biocrude, (hydro)char
    under elevated temperature and pressure.
    
    Parameters
    ----------
    ins : iterable
        dewatered_solids.
    outs : iterable
        hydrochar, HTLaqueous, biocrude, offgas.
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
    # TODO: what if there is N loss before composting, is protein_2_N here still accurate?
    protein_2_N : float
        N to protein ratio, [-].
    N_2_P, float
        P to N ratio, [-].
    eff_P : float
        HTL effluent pressure, [Pa].
    
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
    
    _units= {'Wet mass flowrate':'lb/h'}
    
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
                 # TODO: count CH4 as a fugitive GHG emission, assume the offgas is flared (send or not send it to the CHP unit) and use the fugitive ratio the same as other units
                 gas_composition={'CH4': 0.050, 'C2H6': 0.032, 'CO2': 0.918},
                 T=350 + _C_to_K, P=3049.7*_psi_to_Pa, eff_T=60 + _C_to_K,
                 eff_P=30*_psi_to_Pa):
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
    
    def _run(self):
        dewatered_solids = self.ins[0]
        hydrochar, HTLaqueous, biocrude, offgas = self.outs
        
        hydrochar.phase = 's'
        offgas.phase = 'g'
        HTLaqueous.phase = biocrude.phase = 'l'
        
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
        self._decorated_cost()
    
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
        return min(self.outs[2].F_mass*self.biocrude_C_ratio, self.solids_C)
    
    @property
    def HTLaqueous_C(self):
        return min(self.outs[1].F_vol*1000*self.HTLaqueous_C_slope*\
                   self.solids_dw_protein*100/1000000/self.TOC_TC,
                   self.solids_C - self.biocrude_C)
    
    @property
    def biocrude_H(self):
        return self.outs[2].F_mass*self.biocrude_H_ratio
    
    @property
    def biocrude_N(self):
        return min(self.outs[2].F_mass*self.biocrude_N_ratio, self.solids_N)
    
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
        return min(self.outs[0].F_mass*self.hydrochar_C_ratio, self.solids_C -\
                   self.biocrude_C - self.HTLaqueous_C - self.offgas_C)
    
    @property
    def hydrochar_P(self):
        return min(self.solids_P*self.hydrochar_P_recovery_ratio, self.outs[0].F_mass)
    
    @property
    def HTLaqueous_N(self):
        return self.solids_N - self.biocrude_N
    
    @property
    def HTLaqueous_P(self):
        return self.solids_P*(1 - self.hydrochar_P_recovery_ratio)

# =============================================================================
# HydrothermalAlkalineTreatment
# =============================================================================

# TODO: need to decide whether to recover N and P or not (to the same complexity level as other thermochemical units)

# TODO: same as HTL for now (HTL itself may get updated)
@cost(ID='HTL system', basis='Wet mass flowrate', units='lb/h',
      cost=18743378, S=306198, CE=qs.CEPCI_by_year[2011], n=0.77, BM=2.1)
@cost(ID='Solids filter oil/water separator', basis='Wet mass flowrate', units='lb/h',
      cost=3945523, S=1219765, CE=qs.CEPCI_by_year[2011], n=0.68, BM=1.9)
@cost(ID='Hot oil system', basis='Wet mass flowrate', units='lb/h',
      cost=4670532, S=306198, CE=qs.CEPCI_by_year[2011], n=0.6, BM=1.4)
class HydrothermalAlkalineTreatment(SanUnit):
    '''
    HTL converts feedstock to gas, aqueous, biocrude, (hydro)char
    under elevated temperature and pressure.
    
    Parameters
    ----------
    ins : iterable
        dewatered_solids, sodium_hydroxide, hydrochloric_acid.
    outs : iterable
        hydrochar, HTLaqueous, biocrude, offgas.
    NaOH_conc : float
        the concentration of NaOH, [M].
    HALT_biocrude_adjustment_factor : float
        the biocrude yield in HALT compared to that in HTL, [-].
    HALT_hydrochar_adjustment_factor : float
        the hydrochar yield in HALT compared to that in HTL, [-].
    
    See Also
    --------
    :class:`exposan.landscape_sanunits.HydrothermalLiquefaction`
    '''
    _N_ins = 3
    _N_outs = 4
    
    _units= {'Wet mass flowrate':'lb/h'}
    
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
                 # TODO: count CH4 as a fugitive GHG emission, assume the offgas is flared (send or not send it to the CHP unit) and use the fugitive ratio the same as other units
                 gas_composition={'CH4': 0.050, 'C2H6': 0.032, 'CO2': 0.918},
                 # Table 1, A.J.K. modeling paper
                 HALT_biocrude_adjustment_factor=0.85,
                 # Table 1, A.J.K. modeling paper
                 HALT_hydrochar_adjustment_factor=0.2,
                 T=350 + _C_to_K, P=3049.7*_psi_to_Pa, eff_T=60 + _C_to_K,
                 eff_P=30*_psi_to_Pa):
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
    
    def _run(self):
        dewatered_solids, sodium_hydroxide, hydrochloric_acid = self.ins
        hydrochar, HTLaqueous, biocrude, offgas = self.outs
        
        hydrochar.phase = 's'
        offgas.phase = 'g'
        HTLaqueous.phase = biocrude.phase = 'l'
        
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
        # HTLaqueous is TDS in aqueous phase
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
        
        HTLaqueous.imass['HTLaqueous'] = dewatered_solids_afdw - biocrude.imass['Biocrude'] - hydrochar.imass['Hydrochar'] - gas_mass_flow
        # assume ash (all soluble based on [1]) goes to water
        HTLaqueous.imass['H2O'] = dewatered_solids.F_mass + sodium_hydroxide.F_mass + hydrochloric_acid.F_mass - hydrochar.F_mass -\
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
        self._decorated_cost()
    
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
        return min(self.outs[2].F_mass*self.biocrude_C_ratio, self.solids_C)
    
    @property
    def HTLaqueous_C(self):
        return min(self.outs[1].F_vol*1000*self.HTLaqueous_C_slope*\
                   self.solids_dw_protein*100/1000000/self.TOC_TC,
                   self.solids_C - self.biocrude_C)
    
    @property
    def biocrude_H(self):
        return self.outs[2].F_mass*self.biocrude_H_ratio
    
    @property
    def biocrude_N(self):
        return min(self.outs[2].F_mass*self.biocrude_N_ratio, self.solids_N)
    
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
        return min(self.outs[0].F_mass*self.hydrochar_C_ratio, self.solids_C -\
                   self.biocrude_C - self.HTLaqueous_C - self.offgas_C)
    
    @property
    def hydrochar_P(self):
        return min(self.solids_P*self.hydrochar_P_recovery_ratio, self.outs[0].F_mass)
    
    @property
    def HTLaqueous_N(self):
        return self.solids_N - self.biocrude_N
    
    @property
    def HTLaqueous_P(self):
        return self.solids_P*(1 - self.hydrochar_P_recovery_ratio)

# =============================================================================
# HTLmixer
# =============================================================================

class HTLmixer(SanUnit):
    '''
    A fake unit that calculates C, N, P, and H2O amount in the mixture of HTL
    aqueous and AcidEx effluent.
    
    Parameters
    ----------
    ins : iterable
        HTLaqueous, extracted.
    outs : iterable
        mixture.
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)

    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
        
    def _run(self):
        HTLaqueous, extracted = self.ins
        mixture = self.outs[0]
        
        mixture.mix_from(self.ins)
        mixture.empty()
        
        mixture.imass['C'] = self.ins[0]._source.HTLaqueous_C
        mixture.imass['N'] = self.ins[0]._source.HTLaqueous_N
        mixture.imass['P'] = self.ins[0]._source.HTLaqueous_P +\
                             extracted.imass['P']
        # other compositions represented by H2O except C, N, P
        mixture.imass['H2O'] = HTLaqueous.F_mass + extracted.F_mass -\
                               mixture.imass['C'] - mixture.imass['N'] -\
                               mixture.imass['P']
    
    @property
    def pH(self):
        HTLaqueous, extracted = self.ins
        mixture = self.outs[0]

        base_mol_per_h = HTLaqueous.imass['NaOH']*1000/40
        acid_mol_per_h = extracted.imass['H2SO4']*1000/98*2
        
        volume_L_per_h = mixture.F_vol*1000
        
        base_M = base_mol_per_h/volume_L_per_h
        acid_M = acid_mol_per_h/volume_L_per_h
        
        if base_M > acid_M:
            hydrogen_ion_M = 10**-14/(base_M-acid_M)
        elif base_M == acid_M:
            hydrogen_ion_M = 10**(-7)
        else:
            hydrogen_ion_M = acid_M - base_M

        return -log(hydrogen_ion_M, 10)

# =============================================================================
# CatalyticHydrothermalGasification
# =============================================================================

# TODO: any fugitive GHG emissions?

# TODO: check data with [1]
# TODO: try to add N=6 for the CHG reactor, based on [1]
# TODO: check Guard bed
# based on [1]
@cost(ID='Hydrocyclone', basis='Wet mass flowrate', units='lb/h',
      cost=5000000, S=968859, CE=qs.CEPCI_by_year[2009], n=0.65, BM=2.1)
@cost(ID='Guard bed and CHG reactor', basis='Wet mass flowrate', units='lb/h',
      cost=2041875*1.1, S=76235, CE=qs.CEPCI_by_year[2011], n=0.65, BM=2)
class CatalyticHydrothermalGasification(SanUnit):
    '''
    CHG serves to reduce the COD content in the aqueous phase and produce fuel
    gas under elevated temperature (350°C) and pressure. The outlet will be
    cooled down and separated by a flash unit.
    
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
    
    _units= {'Wet mass flowrate':'lb/h',
              'Hydrocyclone weight':'lb'}
    
    auxiliary_unit_names=('pump','hx','inf_hx','eff_hx')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20, gas_C_2_total_C=0.5981,
                 # TODO: count CH4 as a fugitive GHG emission, assume the offgas is flared (send or not send it to the CHP unit) and use the fugitive ratio the same as other units
                 gas_composition={'CH4': 0.527, 'CO2': 0.432, 'C2H6': 0.011,
                                  'C3H8': 0.030, 'H2': 0.0001},
                 WHSV=3.562, catalyst_lifetime=7920, T=350 + _C_to_K,
                 P=3089.7*_psi_to_Pa, eff_T=60 + _C_to_K,
                 eff_P=3089.7*_psi_to_Pa):
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
# Acid Extraction
# =============================================================================

class AcidExtraction(Reactor):
    '''
    H2SO4 is added to hydrochar from HTL to extract phosphorus.
    
    Parameters
    ----------
    ins : iterable
        hydrochar, acid.
    outs : iterable
        residual, extracted.
    
    See Also
    --------
    :class:`exposan.htl._sanunits.AcidExtraction`
    
    References
    ----------
    .. [1] Zhu, Y.; Schmidt, A.; Valdez, P.; Snowden-Swan, L.; Edmundson, S.
           Hydrothermal Liquefaction and Upgrading of Wastewater-Grown Microalgae:
           2021 State of Technology; PNNL-32695, 1855835; 2022; p PNNL-32695, 1855835.
           https://doi.org/10.2172/1855835.
    '''
    _N_ins = 2
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=20, acid_vol=7,
                 P_acid_recovery_ratio=0.8, P=None, tau=2, V_wf=0.8, # tau: [1]
                 length_to_diameter=2, N=1, V=10, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0, # use MixTank default value
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 304', # acid condition
                 vessel_type='Vertical'):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.acid_vol = acid_vol
        self.P_acid_recovery_ratio = P_acid_recovery_ratio
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
    
    def _run(self):
        hydrochar, acid = self.ins
        residual, extracted = self.outs
        
        self.HTL = self.ins[0]._source
        
        if hydrochar.F_mass <= 0:
            pass
        else:
            if self.HTL.hydrochar_P <= 0:
                residual.copy_like(hydrochar)
            else: 
                acid.imass['H2SO4'] = hydrochar.F_mass*self.acid_vol*0.5*98.079/1000
                # 0.5 M H2SO4 acid_vol (10 mL/1 g) hydrochar
                # 0.5 M H2SO4 density: 1.03 kg/L 
                # https://www.fishersci.se/shop/products/sulfuric-acid-0-5m-4/11943303 (accessed 2024-08-03)
                acid.imass['H2O'] = hydrochar.F_mass*self.acid_vol*1.03 -\
                                    acid.imass['H2SO4']
                
                residual.imass['Residual'] = hydrochar.F_mass - self.ins[0]._source.\
                                             hydrochar_P*self.P_acid_recovery_ratio
                
                extracted.copy_like(acid)
                extracted.imass['P'] = hydrochar.F_mass - residual.F_mass
                # assume just P can be extracted
                
                residual.phase = 's'
                
                residual.T = extracted.T = hydrochar.T
                residual.P = hydrochar.P
                # H2SO4 reacts with hydrochar to release heat and temperature will increase
    
    @property
    def residual_C(self):
        return self.ins[0]._source.hydrochar_C
    
    @property
    def residual_P(self):
        return self.ins[0]._source.hydrochar_P - self.outs[1].imass['P']
    
    def _design(self):
        # conservatively based on the orginal HTL model
        # assuming 1 unit for 250 kg/h of hydrochar
        self.N = ceil(self.ins[0].F_mass/250)
        self.P = self.ins[1].P
        Reactor._design(self)

# =============================================================================
# StruvitePrecipitation
# =============================================================================

class StruvitePrecipitation(Reactor):
    '''
    Extracted and HTL aqueous are mixed together before adding MgCl2 for struvite precipitation.
    If mol(N)<mol(P), add NH4Cl to mol(N):mol(P)=1:1.
    
    Parameters
    ----------
    ins : iterable
        mixture, supply_MgCl2, supply_NH4Cl, base.
    outs : iterable
        struvite, effluent.
    
    See Also
    --------
    :class:`exposan.htl._sanunits.StruvitePrecipitation`
    
    References
    ----------
    .. [1] Zhu, Y.; Schmidt, A.; Valdez, P.; Snowden-Swan, L.; Edmundson, S.
           Hydrothermal Liquefaction and Upgrading of Wastewater-Grown Microalgae:
           2021 State of Technology; PNNL-32695, 1855835; 2022; p PNNL-32695, 1855835.
           https://doi.org/10.2172/1855835.
    .. [2] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
           Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
           Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
           Process Design and Economics for the Conversion of Algal Biomass to
           Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
           PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    '''
    _N_ins = 4
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', lifetime=20, target_pH = 9,
                  Mg_P_ratio=1, P_pre_recovery_ratio=0.828, # [1]
                  HTLaqueous_NH3_N_2_total_N = 0.853, # [2]
                  P=None, tau=1, # [1]
                  V_wf=0.8, length_to_diameter=2, N=1, V=20, auxiliary=False,
                  mixing_intensity=None, kW_per_m3=0, # use MixTank default value
                  wall_thickness_factor=1, vessel_material='Carbon steel', # basic condition
                  vessel_type='Vertical'):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.target_pH = target_pH
        self.Mg_P_ratio = Mg_P_ratio
        self.P_pre_recovery_ratio = P_pre_recovery_ratio
        self.HTLaqueous_NH3_N_2_total_N = HTLaqueous_NH3_N_2_total_N
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        
    def _run(self):
        mixture, supply_MgCl2, supply_NH4Cl, base = self.ins
        struvite, effluent = self.outs
        
        self.HTLmixer = self.ins[0]._source
        
        if self.HTLmixer.outs[0].imass['P'] == 0:
            effluent.copy_like(mixture)
        else:
            old_pH = self.HTLmixer.pH
            if old_pH > 9:
                base.imass['MgO'] = 0
            elif old_pH >= 7:
                OH_M = 10**(old_pH-14)
                OH_M_needed = 10**(self.target_pH-14) - OH_M
                base.imass['MgO'] = OH_M_needed/2 * 40.3044/1000*self.ins[0].F_vol*1000
            else:
                neutral_OH_M = 10**(-old_pH)
                to_target_OH_M = 10**(self.target_pH - 14)
                OH_M_needed = neutral_OH_M + to_target_OH_M
                base.imass['MgO'] = OH_M_needed/2 * 40.3044/1000*self.ins[0].F_vol*1000
            
            supply_MgCl2.imass['MgCl2'] = max((mixture.imass['P']/30.973762*self.Mg_P_ratio -\
                                            base.imass['MgO']/40.3044)*95.211, 0)
            
            if mixture.imass['P']/30.973762 > mixture.imass['N']*self.HTLaqueous_NH3_N_2_total_N/14.0067:
            # if P > N, add NH4Cl to make sure N ≥ P
                supply_NH4Cl.imass['NH4Cl'] = (mixture.imass['P']/30.973762 - mixture.imass['N']*\
                                                self.HTLaqueous_NH3_N_2_total_N/14.0067)*53.491
            
            struvite.imass['Struvite'] = mixture.imass['P']*\
                                          self.P_pre_recovery_ratio/\
                                          30.973762*245.41
            supply_MgCl2.phase = supply_NH4Cl.phase = base.phase = 's'
            
            effluent.copy_like(mixture)
            effluent.imass['P'] -= struvite.imass['Struvite']*30.973762/245.41
            effluent.imass['N'] += supply_NH4Cl.imass['NH4Cl']*14.0067/53.491 -\
                                    struvite.imass['Struvite']*14.0067/245.41
            effluent.imass['H2O'] = self.F_mass_in - struvite.F_mass -\
                                    effluent.imass['C'] - effluent.imass['N'] -\
                                    effluent.imass['P']
            struvite.phase = 's'    
            
            struvite.T = mixture.T
            effluent.T = mixture.T
    
    @property
    def struvite_P(self):
        return self.outs[0].imass['Struvite']*30.973762/245.41
    
    @property
    def struvite_N(self):
        return self.struvite_P*14.0067/30.973762
    
    def _design(self):
        # conservatively based on the orginal HTL model
        # assuming 1 unit for 10 m3/h of hydrochar
        self.N = ceil(self.ins[0].F_vol/10)
        self.P = self.ins[0].P
        Reactor._design(self)

# =============================================================================
# MembraneDistillation
# =============================================================================

class MembraneDistillation(SanUnit):
    '''
    Membrane distillation recovers nitrogen as ammonia sulfate based on vapor
    pressure difference across the hydrophobic membrane. Water flux across
    membrane is ignored.
    
    Parameters
    ----------
    ins : Iterable(stream)
        influent, acid, base, mem_in.
    outs : Iterable(stream)
        ammonium sulfate, ww, mem_out.
    
    See Also
    --------
    :class:`qsdsan.sanunits.MembraneDistillation`
        
    References
    ----------
    [1] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
        https://doi.org/10.1021/acs.est.8b04035.
    [2] Doran, P. M. Chapter 11 - Unit Operations. In Bioprocess Engineering
        Principles (Second Edition); Doran, P. M., Ed.; Academic Press: London,
        2013; pp 445–595. https://doi.org/10.1016/B978-0-12-220851-5.00011-3.
    [3] Spiller, L. L. Determination of Ammonia/Air Diffusion Coefficient Using
        Nafion Lined Tube. Analytical Letters 1989, 22 (11–12), 2561–2573.
        https://doi.org/10.1080/00032718908052375.
    [4] Scheepers, D. M.; Tahir, A. J.; Brunner, C.; Guillen-Burrieza, E.
        Vacuum Membrane Distillation Multi-Component Numerical Model for Ammonia
        Recovery from Liquid Streams. Journal of Membrane Science
        2020, 614, 118399. https://doi.org/10.1016/j.memsci.2020.118399.
    [5] Ding, Z.; Liu, L.; Li, Z.; Ma, R.; Yang, Z. Experimental Study of Ammonia
        Removal from Water by Membrane Distillation (MD): The Comparison of Three
        Configurations. Journal of Membrane Science 2006, 286 (1), 93–103.
        https://doi.org/10.1016/j.memsci.2006.09.015.
    [6] Al-Obaidani, S.; Curcio, E.; Macedonio, F.; Di Profio, G.; Al-Hinai, H.;
        Drioli, E. Potential of Membrane Distillation in Seawater Desalination:
        Thermal Efficiency, Sensitivity Study and Cost Estimation.
        Journal of Membrane Science 2008, 323 (1), 85–98.
        https://doi.org/10.1016/j.memsci.2008.06.006.
    [7] Kogler, A.; Farmer, M.; Simon, J. A.; Tilmans, S.; Wells, G. F.;
        Tarpeh, W. A. Systematic Evaluation of Emerging Wastewater Nutrient Removal
        and Recovery Technologies to Inform Practice and Advance Resource
        Efficiency. ACS EST Eng. 2021, 1 (4), 662–684.
        https://doi.org/10.1021/acsestengg.0c00253.
    [8] Pikaar, I.; Guest, J.; Ganigue, R.; Jensen, P.; Rabaey, K.; Seviour, T.;
        Trimmer, J.; van der Kolk, O.; Vaneeckhaute, C.; Verstraete, W.; Resource
        Recovery from Water: Principles and Applicaiton. IWA 2022.
    '''
    _N_ins = 3
    _N_outs = 3
    _F_BM_default = {'Membrane': 1}
    
    _units = {'Area':'m2',
              'Total volume':'m3'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=5,
                 influent_pH=8.16, # CHG effluent pH: 8.16 ± 0.25 [1]
                 target_pH=10,
                 N_S_ratio=2,
                 # S is excess since not all N can be transferred to form ammonia sulfate
                 # for now, assume N_S_ratio = 2 is ok
                 m2_2_m3=1/1200, # specific surface area, for hollow fiber membrane [2]
                 Dm=2.28*10**(-5), # (2.28 ± 0.12)*10^-5 m^2/s NH3 molecular diffusivity in air [3]
                 # (underestimate, this value may be at 15 or 25 C, our feed is 60 C, should be higher)
                 porosity=0.9, # [4]
                 thickness=7*10**(-5), # m [4]
                 tortuosity=1.2, # [4]
                 Henry=1.61*10**(-5), # atm*m3/mol
                 # https://webwiser.nlm.nih.gov/substance?substanceId=315&identifier=\
                 # Ammonia&identifierType=name&menuItemId=81&catId=120#:~:text=The%20\
                 # Henry's%20Law%20constant%20for,m%2Fmole(2). (accessed 2022-11-11)
                 Ka=1.75*10**(-5), # overall mass transfer coefficient 1.2~2.3*10^-5 m/s [5]
                 capacity=6.01, # kg/m2/h [6]
                 # permeate flux, similar values can be found in many other papers ([176], [222] ,[223] in [7])
                 # [177], [223] in [7] show high nitrogen recovery ratio (>85% under optimal conditions)
                 membrane_price=90/GDPCTPI[2008]*GDPCTPI[2022]): #[6]
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.influent_pH = influent_pH
        self.target_pH = target_pH
        self.N_S_ratio = N_S_ratio
        self.m2_2_m3 = m2_2_m3
        self.Dm = Dm
        self.porosity = porosity
        self.thickness = thickness
        self.tortuosity = tortuosity
        self.Henry = Henry
        self.Ka = Ka
        self.capacity = capacity
        self.membrane_price = membrane_price
    
    def _run(self):
        influent, acid, base = self.ins
        ammoniumsulfate, ww, ammoniumsulfatesolution = self.outs
        
        self.CHG = self.ins[0]._source.ins[0]._source.ins[0]._source
        
        if self.CHG.CHGout_N == 0:
            ww.copy_like(influent)
            self.membrane_area = self.F_vol_in*1000/self.capacity
        else:
            NaOH_conc = 10**(self.target_pH - 14) - 10**(self.influent_pH - 14)
            NaOH_mol = NaOH_conc*self.ins[0].F_mass
            base.imass['NaOH'] = NaOH_mol*39.997/1000
            
            # 0.5 M H2SO4 solution, density: 1.03 g/mL (https://www.fishersci.se/shop/products/sulfuric-acid-0-5m-4/11943303, accessed 2024-08-03)
            acid.imass['H2SO4'] = self.CHG.CHGout_N/14.0067/self.N_S_ratio*98.079
            acid.imass['H2O'] = acid.imass['H2SO4']*1000/98.079/0.5*1.03 -\
                                acid.imass['H2SO4']
            
            # calculation of pKa under different T follows van't Hoff relationship [8] page 292
            self.pKa = pKa = -log(exp(52.22/8.3145*1000*(1/298 - 1/self.ins[0].T))*10**(-9.252), 10)

            ammonia_to_ammonium = 10**(-pKa)/10**(-self.target_pH)
            ammonia = self.CHG.CHGout_N*ammonia_to_ammonium/(1 +\
                       ammonia_to_ammonium)*17.031/14.0067
            others = influent.F_mass - ammonia
            
            self.membrane_area = self.F_vol_in*1000/self.capacity
            # N2:O2 = 0.79:0.21 in the air, air density is 1.204 kg/m3
            # https://en.wikipedia.org/wiki/Density_of_air#:~:text=It%20also%20\
            # changes%20with%20variation,International%20Standard%20Atmosphere%2\
            # 0(ISA). (accessed 2022-11-14)
            N2_in_air = self.membrane_area*self.m2_2_m3*self.porosity*0.79*1.204
            O2_in_air = self.membrane_area*self.m2_2_m3*self.porosity*0.21*1.204
            
            # N2 amount will be changed based on design, maybe also add O2 (N2:O2 = 4:1)
            imass = indexer.MassFlowIndexer(l=[('H2O', others),
                                               ('NH3', ammonia),
                                               ('N2', 0),
                                               ('O2', 0)],
                                            g=[('H2O', 0),
                                               ('NH3', 0), 
                                               ('N2', N2_in_air),
                                               ('O2', O2_in_air)])
            
            vle = equilibrium.VLE(imass)
            vle(T=influent.T, P=influent.P)
            X_NH3_f_m = vle.imol['g','NH3']/(vle.imol['g','H2O'] + vle.imol['g','NH3'])
            X_NH3_f = vle.imol['l','NH3']/(vle.imol['l','H2O'] + vle.imol['l','NH3'])
    
            km = self.Dm*self.porosity/self.tortuosity/self.thickness
            
            # https://www.sciencedirect.com/topics/chemistry/henrys-law#:~:text=\
            # Values%20for%20Henry's%20law%20constants,gas%20constant%20(8.20575%\
            # 20%C3%97%2010 (accessed 2022-11-11)
            dimensionless_Henry = self.Henry/8.20575/(10**(-5))/influent.T # H' = H/RT
            
            kf = 1/(1/self.Ka - 1/dimensionless_Henry/km*(1 + (10**-14/10**(-pKa))*10**\
                  (-self.target_pH)/(10**(-14))))
            
            J = kf*ammonia/influent.F_mass*1000*log(X_NH3_f_m/X_NH3_f)*3600 # in kg/m2/h
            
            NH3_mass_flow = J*self.membrane_area
            
            ammonia_transfer_ratio = min(1, NH3_mass_flow/ammonia)
            
            ammoniumsulfate.imass['NH42SO4'] = ammonia*ammonia_transfer_ratio/34.062*132.14
            ammoniumsulfatesolution.imass['H2O'] = acid.imass['H2O']
            ammoniumsulfatesolution.imass['H2SO4'] = acid.imass['H2SO4'] -\
                                             ammoniumsulfate.imass['NH42SO4']/\
                                             132.14*98.079
            
            # ww has the same T and P as influent          
            ww.copy_like(influent)
            
            ww.imass['N'] = self.CHG.CHGout_N*(1 - ammonia_to_ammonium/(1 +\
                             ammonia_to_ammonium)*ammonia_transfer_ratio)
                                                              
            ww.imass['C'] = self.CHG.CHGout_C
                             
            ww.imass['P'] = self.CHG.CHGout_P
                       
            ww.imass['H2O'] -= (ww.imass['C'] + ww.imass['N'] + ww.imass['P'])
            
            ww.imass['H2O'] += self.ins[2].F_mass
            
            # ammonium sulfate has the same T and P as acid
            ammoniumsulfate.T = ammoniumsulfatesolution.T = acid.T
            ammoniumsulfate.P = ammoniumsulfatesolution.P = acid.P
    
    @property
    def N_recovery_ratio(self):
        return 1 - self.outs[1].imass['N']/self.CHG.CHGout_N
    
    def _design(self):
        Design = self.design_results
        Design['Area'] = self.membrane_area
        Design['Total volume'] = Design['Area']*self.m2_2_m3
    
    def _cost(self):
        Design = self.design_results
        purchase_costs = self.baseline_purchase_costs
        purchase_costs['Membrane'] = Design['Area']*self.membrane_price

# =============================================================================
# SupercriticalWaterOxidation
# =============================================================================

# TODO: add citations from the excel (data_collection.xlsx)
# the cost here include HX and other auxiliary units
# n=0.7 and BM=2 to be similar to HTL-based systems
@cost(ID='SCWO system', basis='Wet mass flowrate', units='tonne/day',
      cost=1993701, S=24, CE=qs.CEPCI_by_year[2023], n=0.7, BM=2)
class SupercriticalWaterOxidation(SanUnit):
    '''
    SCWO converts feedstock to ash miniral and gas (CO2 + H2O) under elevated
    temperature and pressure.
    
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
        
    References
    ----------
    .. [1] Feng, S.; Guanghua, L. Chapter 4 - Hydrothermal and Solvothermal
           Syntheses. In Modern Inorganic Synthetic Chemistry; Xu, R., Pang, W.,
           Huo, Q., Eds.; Elsevier: Amsterdam, 2011; pp 63–95.
           https://doi.org/10.1016/B978-0-444-53599-3.10004-6.
    '''
    _N_ins = 1
    _N_outs = 2
    
    _units= {'Wet mass flowrate':'tonne/day'}
    
    auxiliary_unit_names = ('pump','hx','inf_hx','eff_hx','kodrum')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=15, lipid_2_C=0.750,
                 protein_2_C=0.545, carbo_2_C=0.400,
                 # T and P, [1]
                 T=550 + _C_to_K, P=3e7,
                 # assume to be the same as HTL
                 eff_T=60 + _C_to_K, eff_P=30*_psi_to_Pa):
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
        D['Wet mass flowrate'] = self.ins[0].F_mass/1000*24
        
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