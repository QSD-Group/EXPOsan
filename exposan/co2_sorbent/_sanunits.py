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
from biosteam import Unit
from qsdsan.utils import auom
from biosteam.units.decorators import cost

__all__ = (
    'BauxiteHammerMill',
    'ALFProduction',
    'PhaseChanger',
    'SolidPressureFilter',
    'ALFCrystallizer',
    'ALFPressureFilter',
    'ReverseOsmosis',
    'S2WS',
    'ALFTSA',
    'CO2ElectrolyzerSystem'
    )

_mgd_to_cmh = auom('gallon').conversion_factor('m3')*1e6/24
CEPCI = bst.units.design_tools.CEPCI_by_year

# =============================================================================
# BauxiteHammerMill
# =============================================================================

@cost(basis='Flow rate', ID='Hammer mill', units='kg/hr',
           kW=13.23*4.54, cost=105225, S=4540, CE=CEPCI[2011], n=1, BM=1,
           lifetime=40000/(365*24))
class BauxiteHammerMill(SanUnit):
    '''
    See biorefineries/ethanol_adipic/_preprocessing.py
    Note there are two types of hammer mill: CPP (conventional pelleting process)
    hammer mill and HMPP (high-moisture pelleting process) hammer mill.
    CPP hammer mill is used for this system.
    '''

# =============================================================================
# ALFProduction
# =============================================================================
class ALFProduction(bst.CSTR):
    '''
    Reactor for ALF production. See biosteam/units/stirred_tank_reactor.py.
    '''
    _N_ins = 1
    _N_outs = 1
    
    # TODO: add a parameter for X (conversion rate), but it seems that X=1 makes sense since Al(OH)3 is a solid and HCOOH is over amount
    # TODO: check the reaction between Fe2O3 and HCOOH
    def _setup(self):
            self.ALF_production_AlH3O3 = bst.Reaction('AlH3O3,s + 3HCOOH,l -> C3H3AlO6,s + 3H2O,l', 'AlH3O3', 1)
            self.ALF_production_bauxite = bst.Reaction('Al2O3,s + 6HCOOH,l -> 2C3H3AlO6,s + 3H2O,l', 'Al2O3', 1)
            self.Fe_side_reaciton = bst.Reaction('Fe2O3,s + 6HCOOH,l -> 2Fe,s + 3H2O,l + 3CO2,l', 'Fe2O3', 1)
    
    def _run(self):
        effluent = self.outs[0]
        effluent.copy_like(self.ins[0])
        self.ALF_production_AlH3O3(effluent)
        self.ALF_production_bauxite(effluent)
        self.Fe_side_reaciton(effluent)
        effluent.T = self.T
        effluent.P = self.P

# =============================================================================
# PhaseChanger        
# =============================================================================
class PhaseChanger(Unit):
    '''
    A fake unit that change the phase of streams (specifically developed for the system B).
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    ins : iterable
        Inlet streams.
    outs : iterable
        Outlet streams.
    '''
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=(), outs=(), thermo=None):
        super().__init__(ID, ins, outs, thermo)
        
    def _run(self):
        inlet = self.ins[0]
        outlet, carbon_dioxide = self.outs
        
        carbon_dioxide.phase='g'
        carbon_dioxide.imass['CO2'] = inlet.imass['CO2']
        outlet.copy_like(inlet)
        outlet.imass['l','CO2'] = 0
        outlet.imass['l','C3H3AlO6'] = inlet.imass['C3H3AlO6']
        outlet.imass['s','C3H3AlO6'] = 0

# =============================================================================
# SolidPressureFilter
# =============================================================================
_hp2kW = 0.7457
@cost(basis='Retentate flow rate', ID='Flitrate tank agitator',
      cost=26e3, CE=551, kW=7.5*_hp2kW, S=31815, n=0.5, BM=1.5)
@cost(basis='Retentate flow rate', ID='Discharge pump',
      cost=13040, CE=551, S=31815, n=0.8, BM=2.3)
@cost(basis='Retentate flow rate', ID='Filtrate tank',
      cost=103e3, S=31815, CE=551, BM=2.0, n=0.7)
@cost(basis='Retentate flow rate', ID='Feed pump', kW=74.57,
      cost= 18173, S=31815, CE=551, n=0.8, BM=2.3)
@cost(basis='Retentate flow rate', ID='Stillage tank 531',
      cost=174800, CE=551, S=31815, n=0.7, BM=2.0)
@cost(basis='Retentate flow rate', ID='Mafifold flush pump', kW=74.57,
      cost=17057, CE=551, S=31815, n=0.8, BM=2.3)
@cost(basis='Retentate flow rate', ID='Recycled water tank',
      cost=1520, CE=551, S=31815, n=0.7, BM=3.0)
@cost(basis='Retentate flow rate', ID='Wet cake screw',  kW=15*_hp2kW,
      cost=2e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost(basis='Retentate flow rate', ID='Wet cake conveyor', kW=10*_hp2kW,
      cost=7e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost(basis='Retentate flow rate', ID='Pressure filter',
      cost=3294700, CE=551, S=31815, n=0.8, BM=1.7)
@cost(basis='Retentate flow rate', ID='Pressing air compressor receiver tank',
      cost=8e3, CE=551, S=31815, n=0.7, BM=3.1)
@cost(basis='Retentate flow rate', ID='Cloth wash pump', kW=150*_hp2kW,
      cost=29154, CE=551, S=31815, n=0.8, BM=2.3)
@cost(basis='Retentate flow rate', ID='Dry air compressor receiver tank',
      cost=17e3, CE=551, S=31815, n=0.7, BM=3.1)
@cost(basis='Retentate flow rate', ID='Pressing air pressure filter',
      cost=75200, CE=521.9, S=31815, n=0.6, kW=112, BM=1.6)
@cost(basis='Retentate flow rate', ID='Dry air pressure filter (2)',
      cost=405000, CE=521.9, S=31815, n=0.6, kW=1044, BM=1.6)
class SolidPressureFilter(Unit):
    """
    Create a pressure filter for the separation of ALF. See biosteam/units/solids_separation.py.
    Capital costs are based on [1]_.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    ins : iterable
        Inlet streams.
    outs : iterable
        Outlet streams.
    moisture_content : float, optional
        Moisture content of retentate. Defaults to 0.35.
    split : array_like or dict[str, float]
        Splits of chemicals to the retentate.
        Assume completely separation for ALF.
        From biosteam/units/solids_separation.py: soluble chemicals~0.036.
    
    References
    ----------
    [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    """
    _units = {'Retentate flow rate': 'kg/hr'}
    
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=(), outs=(), thermo=None,
                 moisture_content=0.35, split={'SiO2':1,
                                               'Fe':1,
                                               'C3H3AlO6':0.036,
                                               'HCOOH':0.036}):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo)
        self.moisture_content = moisture_content
        self.split = split
    
    def _run(self):
        inlet = self.ins[0]
        retentate, permeate = self.outs
        
        retentate.empty()
        retentate.phases = ('s','l')
        
        retentate.imass['s','SiO2'] = inlet.imass['s','SiO2']*self.split['SiO2']
        retentate.imass['s','Fe'] = inlet.imass['s','Fe']*self.split['Fe']
        retentate.imass['l','C3H3AlO6'] = inlet.imass['C3H3AlO6']*self.split['C3H3AlO6']
        retentate.imass['l','HCOOH'] = inlet.imass['l','HCOOH']*self.split['HCOOH']
        retentate.imass['l','H2O'] = retentate.F_mass/(1-self.moisture_content)*self.moisture_content
        
        permeate.phases = ('s','l')
        
        permeate.imass['s','SiO2'] = inlet.imass['s','SiO2']*(1-self.split['SiO2'])
        permeate.imass['s','Fe'] = inlet.imass['s','Fe']*(1-self.split['Fe'])
        permeate.imass['l','C3H3AlO6'] = inlet.imass['C3H3AlO6']*(1-self.split['C3H3AlO6'])
        permeate.imass['l','HCOOH'] = inlet.imass['l','HCOOH']*(1-self.split['HCOOH'])
        permeate.imass['l','H2O'] = inlet.imass['l','H2O'] - retentate.imass['l','H2O']
        
        retentate.T = inlet.T
        permeate.T = inlet.T
        
    def _design(self):
        self.design_results['Retentate flow rate'] = self.outs[0].F_mass

# =============================================================================
# ALFCrystallizer
# =============================================================================
class ALFCrystallizer(bst.BatchCrystallizer):
    '''
    Crystallier for ALF. See biosteam/units/_batch_crystallizer.py.
    
    Parameters
    ----------
    crystal_ALF_yield : float, optional
        ALF crystallization yield. Defacults to 1.
    
    References
    ----------
    [1] Evans, H. A.; Mullangi, D.; Deng, Z.; Wang, Y.; Peh, S. B.; Wei, F.;
        Wang, J.; Brown, C. M.; Zhao, D.; Canepa, P.; Cheetham, A. K.
        Aluminum Formate, Al(HCOO)3: An Earth-Abundant, Scalable, and Highly
        Selective Material for CO2 Capture. Science Advances 2022, 8 (44),
        eade1473. https://doi.org/10.1126/sciadv.ade1473.
    '''
    def __init__(self, ID='', ins=(), outs=(), thermo=None,
                 tau=5, N=3, T=298.15, crystal_ALF_yield=1):
        bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                                       tau=tau, N=N, T=T)
        self.crystal_ALF_yield = crystal_ALF_yield
        
    def _run(self):
        inlet = self.ins[0]
        outlet = self.outs[0]
        outlet.phases = ('s','l')
        crystal_ALF_yield = self.crystal_ALF_yield
        ALF = inlet.imass['C3H3AlO6']
        outlet.empty()
        
        outlet.copy_like(inlet)
        outlet.imass['s','C3H3AlO6'] = ALF*crystal_ALF_yield
        outlet.imass['l','C3H3AlO6'] = ALF*(1-crystal_ALF_yield)
        
        outlet.T = self.T

# =============================================================================
# ALFPressureFilter
# =============================================================================
_hp2kW = 0.7457
@cost(basis='Retentate flow rate', ID='Flitrate tank agitator',
      cost=26e3, CE=551, kW=7.5*_hp2kW, S=31815, n=0.5, BM=1.5)
@cost(basis='Retentate flow rate', ID='Discharge pump',
      cost=13040, CE=551, S=31815, n=0.8, BM=2.3)
@cost(basis='Retentate flow rate', ID='Filtrate tank',
      cost=103e3, S=31815, CE=551, BM=2.0, n=0.7)
@cost(basis='Retentate flow rate', ID='Feed pump', kW=74.57,
      cost= 18173, S=31815, CE=551, n=0.8, BM=2.3)
@cost(basis='Retentate flow rate', ID='Stillage tank 531',
      cost=174800, CE=551, S=31815, n=0.7, BM=2.0)
@cost(basis='Retentate flow rate', ID='Mafifold flush pump', kW=74.57,
      cost=17057, CE=551, S=31815, n=0.8, BM=2.3)
@cost(basis='Retentate flow rate', ID='Recycled water tank',
      cost=1520, CE=551, S=31815, n=0.7, BM=3.0)
@cost(basis='Retentate flow rate', ID='Wet cake screw',  kW=15*_hp2kW,
      cost=2e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost(basis='Retentate flow rate', ID='Wet cake conveyor', kW=10*_hp2kW,
      cost=7e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost(basis='Retentate flow rate', ID='Pressure filter',
      cost=3294700, CE=551, S=31815, n=0.8, BM=1.7)
@cost(basis='Retentate flow rate', ID='Pressing air compressor receiver tank',
      cost=8e3, CE=551, S=31815, n=0.7, BM=3.1)
@cost(basis='Retentate flow rate', ID='Cloth wash pump', kW=150*_hp2kW,
      cost=29154, CE=551, S=31815, n=0.8, BM=2.3)
@cost(basis='Retentate flow rate', ID='Dry air compressor receiver tank',
      cost=17e3, CE=551, S=31815, n=0.7, BM=3.1)
@cost(basis='Retentate flow rate', ID='Pressing air pressure filter',
      cost=75200, CE=521.9, S=31815, n=0.6, kW=112, BM=1.6)
@cost(basis='Retentate flow rate', ID='Dry air pressure filter (2)',
      cost=405000, CE=521.9, S=31815, n=0.6, kW=1044, BM=1.6)
class ALFPressureFilter(Unit):
    """
    Create a pressure filter for the separation of ALF. See biosteam/units/solids_separation.py.
    Capital costs are based on [1]_.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    ins : iterable
        Inlet streams.
    outs : iterable
        Outlet streams.
    moisture_content : float, optional
        Moisture content of retentate. Defaults to 0.35.
    split : array_like or dict[str, float]
        Splits of chemicals to the retentate.
        Assume completely separation for ALF.
        From biosteam/units/solids_separation.py: soluble chemicals~0.036.
    
    References
    ----------
    [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    """
    _units = {'Retentate flow rate': 'kg/hr'}
    
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=(), outs=(), thermo=None,
                 moisture_content=0.35, split={'C3H3AlO6_s':1,
                                               'C3H3AlO6_l':0,
                                               'HCOOH':0.036}):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo)
        self.moisture_content = moisture_content
        self.split = split
    
    def _run(self):
        inlet = self.ins[0]
        retentate, permeate = self.outs
        
        retentate.empty()
        retentate.phases = ('s','l')
        
        retentate.imass['s','C3H3AlO6'] = inlet.imass['s','C3H3AlO6']*self.split['C3H3AlO6_s']
        retentate.imass['l','C3H3AlO6'] = inlet.imass['l','C3H3AlO6']*self.split['C3H3AlO6_l']
        retentate.imass['l','HCOOH'] = inlet.imass['l','HCOOH']*self.split['HCOOH']
        retentate.imass['l','H2O'] = retentate.F_mass/(1-self.moisture_content)*self.moisture_content
        
        permeate.phases = ('s','l')
        
        permeate.imass['s','C3H3AlO6'] = inlet.imass['s','C3H3AlO6']*(1-self.split['C3H3AlO6_s'])
        permeate.imass['l','C3H3AlO6'] = inlet.imass['l','C3H3AlO6']*(1-self.split['C3H3AlO6_l'])
        permeate.imass['l','HCOOH'] = inlet.imass['l','HCOOH']*(1-self.split['HCOOH'])
        permeate.imass['l','H2O'] = inlet.imass['l','H2O'] - retentate.imass['l','H2O']
    
    def _design(self):
        self.design_results['Retentate flow rate'] = self.outs[0].F_mass

# =============================================================================
# ReverseOsmosis
# =============================================================================
@cost(basis='Volumetric flow', ID='Reactor', units='m3/hr',
      cost=2450000, S=2.7*_mgd_to_cmh, CE=CEPCI[2012], n=1, BM=1.8)
@cost(basis='Volumetric flow', ID='Evaporator', units='m3/hr',
      kW=1103.636, cost=5000000, S=2.7*_mgd_to_cmh, CE=CEPCI[2012], n=0.6, BM=1.6)
class ReverseOsmosis(SanUnit):
    """
    See biorefineries/wwt/_wwt_process.py for cost.
    See biosteam/wastewater/conventional.py for the default water recovery rate.
    """
    _N_ins = 1
    _N_outs = 2
    _units = {'Volumetric flow': 'm3/hr'}
    
    def __init__(self, ID='', ins=(), outs=(), thermo=None,
                 init_with='WasteStream', water_recovery=0.987):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with)
        self.water_recovery = water_recovery

    def _run(self):
        influent = self.ins[0]
        water, brine = self.outs

        self.design_results['Volumetric flow'] = self.F_vol_in
        
        water.imass['Water'] = influent.imass['Water']*self.water_recovery
        brine.mol = influent.mol - water.mol
        water.T = brine.T = influent.T

# =============================================================================
# S2WS
# =============================================================================
class S2WS(SanUnit):
    """
    S2WS: Stream to WasteStream.
    A fake unit that enables the calculation of LCA for 'Stream'
    by converting 'Stream' to 'WasteStream'.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    ins : iterable
        Inlet streams.
    outs : iterable
        Outlet streams.
    """
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream'):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with)

    def _run(self):
        inlet = self.ins[0]
        outlet = self.outs[0]
        # if inlet is a 'MultiStream', for now, just assume the outlet
        # contains water that has the same weight as the inlet
        try:
            outlet.copy_like(inlet)
        except TypeError:
            outlet.imass['H2O'] = inlet.F_mass

# =============================================================================
# ALFTSA
# =============================================================================
class ALFTSA(bst.AdsorptionColumnTSA):
    '''
    TSA using ALF as adsorbent for CO2 adsorption.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    ins : iterable
        Inlet streams.
    outs : iterable
        Outlet streams.
    adsorbent : str
        Defaults to 'ALF'.
    adsorbent_cost : dict
        Defaults to {'ALF': 60} [$/ft3].
    # TODO: what does the lifetime mean here, for ALF or for the entire unit?
    equipment_lifetime : dict
        The lifetime of adsorbents. Defaults to {'ALF': 10} [year].
    adsorbate_ID : str
        Defaults to 'CO2'.
    split : dict[str, float] or list[float], optional
        Component splits towards the effluent (0th outlet).
    superficial_velocity : float, optional
        Superficial velocity of the feed. The diameter of the receiving vessel adjusts
        accordingly. Defaults to 1080 m/h. Typical velocities are 540 to 2160 m/h for liquids [1]_.
    regeneration_velocity : float, optional
        Mean velocity of the fluid used to regenerate the bed. Defaults to 540 m/h. 
        Common velocity range for gasses is 504-2160 m / hr [1]_.
    cycle_time : float, optional
        Time at which the receiving vessel is switched. Defaults to 8 h [2]_.
    rho_adsorbent : float, optional
        The density of ALF. Defaults to 1441 kg/m3, Table S1 [3]_.
    adsorbent_capacity : float, optional
        Amount of CO2 that ALF can hold. Defaults to 0.1188 g/g (2.7 mmol/g), Table S7 [3]_.
    T_regeneration : float, optional
        Temperature during the regeneration phase. Defaults to 418 K, Table S8 [3]_.
    vessel_material : float, optional
        Vessel material. Defaults to 'Stainless steel 316',
    vessel_type : float, optional
        Vessel type. Defaults to 'Vertical'.
    length_unused : float, optional
        Additional length of a column to account for mass transfer limitations (due to unused bed). Defaults to 1.219 m.
    
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
    def _init(self, ID='ALF_TSA', ins=(), outs=(),
              adsorbent='ALF',
              adsorbent_cost={'ALF': 1},
              equipment_lifetime={'ALF': 10},
              adsorbate_ID='CO2',
              split=dict(O2=0.0036, N2=0.0036, CO2=0.945),
              superficial_velocity=1080,
              regeneration_velocity=540,
              cycle_time=8,
              rho_adsorbent = 1441,
              adsorbent_capacity=0.1188,
              T_regeneration=418,
              vessel_material='Stainless steel 316',
              vessel_type='Vertical',
              length_unused=1.219):
        super()._init(split=split,
                      adsorbent=adsorbent,
                      adsorbate_ID=adsorbate_ID,
                      superficial_velocity=superficial_velocity,
                      regeneration_velocity=regeneration_velocity,
                      cycle_time=cycle_time,
                      rho_adsorbent = rho_adsorbent,
                      adsorbent_capacity=adsorbent_capacity,
                      T_regeneration=T_regeneration,
                      vessel_material=vessel_material,
                      vessel_type=vessel_type,
                      length_unused=length_unused)
        self.adsorbent_cost = adsorbent_cost
        self._default_equipment_lifetime = self.equipment_lifetime = equipment_lifetime
        
# =============================================================================
# CO2ElectrolyzerSystem
# =============================================================================
# cost basis year: 2010, see Jouny et al. 2018 SI Excel - tab Data - cell C4
@cost(basis='Electrolyzer area', ID='Electrolyzer', units='m^2',
      cost=250.25*1.75*175/1000*10*(1+1/0.65*0.35), S=1, CE=qs.CEPCI_by_year[2010], n=1, BM=1.2)
@cost(basis='Electrolyte flow rate (ethanol)', ID='Distiller (ethanol)', units='L/min',
      cost=4162240, S=1000, CE=qs.CEPCI_by_year[2010], n=0.7)
@cost(basis='Electrolyte flow rate (formic acid)', ID='Distiller (formic acid)', units='L/min',
      cost=6896190, S=1000, CE=qs.CEPCI_by_year[2010], n=0.7)
@cost(basis='Electrolyte flow rate (methanol)', ID='Distiller (methanol)', units='L/min',
      cost=4514670, S=1000, CE=qs.CEPCI_by_year[2010], n=0.7)
@cost(basis='Electrolyte flow rate (propanol)', ID='Distiller (propanol)', units='L/min',
      cost=4687910, S=1000, CE=qs.CEPCI_by_year[2010], n=0.7)
@cost(basis='Total gas flow for PSA', ID='PSA', units='m^3/h',
      cost=1989043, S=1000, CE=qs.CEPCI_by_year[2010], n=0.7)
class CO2ElectrolyzerSystem(SanUnit):
    '''
    CO2 electrolyzer system that converts CO2 into reduced 1C, 2C, and nC products [1]_. 
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    ins : iterable
        Inlet streams.
    outs : iterable
        Outlet streams.
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
    
    def __init__(self, ID='', ins=(), outs=(), thermo=None, target_product='formic acid',
                 current_density=0.2, cathodic_overpotential=0.454, cell_voltage=2.3,
                 product_selectivity=0.9, converstion=0.5, operating_days_per_year=350):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo)
        self.target_product = target_product
        self.current_density = current_density
        self.cathodic_overpotential = cathodic_overpotential
        self.cell_voltage = cell_voltage
        self.product_selectivity = product_selectivity
        self.converstion = converstion
        self.operating_days_per_year = operating_days_per_year
    
    def _run(self):
        carbon_dioxide, water = self.ins
        product, mixed_offgas = self.outs
        
        if self.target_product not in ['carbon monoxide','ethanol','ethylene',
                                       'formic acid','methane','methanol','propanol']:
            raise ValueError("target product must be in 'carbon monoxide', 'ethanol', " + 
                             "'ethylene', 'formic acid', 'methane', 'methanol', and 'propanol'")
        
        product_info = {# the propanol is n-propanol
                        'chemical': ['carbon monoxide','ethanol','ethylene',
                                     'formic acid','methane','methanol','propanol'],
                        'formula': ['CO','C2H6O','C2H4','HCOOH','CH4','CH4O','C3H8O'],
                        'elec_number': [2, 12, 12, 2, 8, 6, 18],
                        'mole_ratio': [1, 2, 2, 1, 1, 1, 3],
                        'MW': [28, 46, 28, 46, 16, 32, 60],
                        # kg/m3
                        'density': [1.14, 789, 1.18, 1221, 0.656, 792, 803],
                        'state': ['gas','liq','gas','liq','gas','liq','liq'],
                        'potential': [-0.106, 0.084, 0.064, -0.25, 0.169, 0.016, 0.095],
                        'cathode_water_production': [1, 3, 4, 0, 2, 1, 5]}
        
        product_info = pd.DataFrame.from_dict(product_info)
        
        chemical_info = self.chemical_info =\
            product_info[product_info['chemical'] == self.target_product]
        
        # CO2 inlet flow rate [kg/h]
        CO2_inlet_flow_rate = carbon_dioxide.imass['CO2']
        
        # converted CO2 amount [kg/day]
        CO2_converted = CO2_inlet_flow_rate*24*self.converstion
        
        # CO2 outlet flow rate [kg/h]
        CO2_outlet_flow_rate_kg_per_h = CO2_inlet_flow_rate - CO2_converted/24
        
        # CO2 outlet flow rate [m3/h]
        # the density of CO2 is 1.98 kg/m3, Jouny et al. 2018
        self.CO2_outlet_flow_rate_m3_per_h = CO2_outlet_flow_rate_kg_per_h/1.98
        
        # production amount based on inlet CO2 [kg/day]
        product_production = CO2_converted/44/float(chemical_info['mole_ratio'])*\
            float(chemical_info['MW'])
        
        # required current [A]
        current_needed = product_production/24/3600*1000/float(chemical_info['MW'])*\
            float(chemical_info['elec_number'])*96485/self.product_selectivity
        
        # required electrolzer area [m2]
        self.electrolyzer_area = current_needed/self.current_density/10000
        
        # required power [MW]
        self.power_needed = current_needed*self.cell_voltage/1000000
        
        # gas product flow rate [m3/h]
        gas_product_flow_rate = 0\
            if chemical_info['state'].to_string(index=False) == 'liq'\
            else product_production/float(chemical_info['density'])/24
        
        # liquid product flow rate [m3/h]
        liquid_product_flow_rate_m3_per_h = 0\
            if chemical_info['state'].to_string(index=False) == 'gas'\
            else product_production/float(chemical_info['density'])/24
        
        # liquid product flow rate [l/min]
        liquid_product_flow_rate_l_per_min = liquid_product_flow_rate_m3_per_h*1000/60
        
        # electrolyte flow rate [l/min]
        # product-rich electrolyte is recycled until a steady-state volume concentration
        # of 10% if reached, Jouny et al. 2018
        self.electrolyte_flow_rate = liquid_product_flow_rate_l_per_min/0.1
        
        # hydrogen flow rate [mol/s]
        # 2 represents 2 e-
        hydrogen_flow_rate_mol_per_s = current_needed*(1-self.product_selectivity)/2/96485
        
        # hydrogen flow rate [m3/h]
        # hydrogen molar mass: 2, hydrogen density: 0.08375 kg/m3, Jouny et al. 2018
        hydrogen_flow_rate_m3_per_h = hydrogen_flow_rate_mol_per_s*2/1000/0.08375*3600
        
        # total gas flow [m3/h]
        self.total_gas_flow = self.CO2_outlet_flow_rate_m3_per_h +\
            gas_product_flow_rate + hydrogen_flow_rate_m3_per_h
        
        product.imass[chemical_info['formula'].to_string(index=False)] = product_production/24
        product.phase = 'g' if chemical_info['state'].to_string(index=False) == 'gas' else 'l'
        
        # TODO: explain why for formic acid, the cathode does not produce water but it still needs distillation
        # note Jouny et al. 2018 calculated process water amount as 'current_needed/4/96485*18/1000*3600'
        # TODO: confirm it is correct to change it to 'current_needed/2/96485*18/1000*3600'
        # as every 18 g H2O can transfer 2 e-
        water.imass['H2O'] = current_needed/2/96485*18/1000*3600
        # TODO: 
        mixed_offgas.empty()
        mixed_offgas.imass['CO2'] = CO2_outlet_flow_rate_kg_per_h
        # hydrogen density: 0.08375 kg/m3, Jouny et al. 2018
        # TODO: hydrogen is produced in the cathode. Should it be in the stream 'product'?
        # Jouny et al. 2018 kind of ingored this
        # TODO: decide how to deal with it. Two cases: gas product (with mixed hydrogen) and liquid product (selling hydrogen?)
        mixed_offgas.imass['H2'] = hydrogen_flow_rate_m3_per_h*0.08375
        mixed_offgas.imass['N2'] = carbon_dioxide.imass['N2']
        
        # H2O after distillation or PSA
        # TODO: clarify why formic acid is produced in the cathode with no water produced but still needs distillation?
        mixed_offgas.imass['H2O'] = product.F_mass/float(chemical_info['MW'])*float(chemical_info['cathode_water_production'])*18
        
        mixed_offgas.imass['O2'] = carbon_dioxide.F_mass + water.F_mass -\
            product.F_mass - mixed_offgas.F_mass
        mixed_offgas.phase = 'g'
    
    def _design(self):
        D = self.design_results
        D['Electrolyzer area'] = self.electrolyzer_area
        
        for liquid_product in ['ethanol','formic acid','methanol','propanol']:
            D[f'Electrolyte flow rate ({liquid_product})'] = self.electrolyte_flow_rate\
                if self.chemical_info['chemical'].to_string(index=False) == liquid_product\
                    else 0
        
        D['Total gas flow for PSA'] = 0\
            if self.CO2_outlet_flow_rate_m3_per_h == self.total_gas_flow\
                else self.total_gas_flow
        
        # PSA operating power: 0.25 kWh/m3, Jouny et al. 2018
        self.add_power_utility(self.power_needed*1000 + D['Total gas flow for PSA']*0.25)
    
    def _cost(self):
        self._decorated_cost()
        
        distillation_OPEX = {'ethanol': 3463310,
                             'formic acid': 11213200,
                             'methanol': 4027820,
                             'propanol': 5610420}
        if self.chemical_info['chemical'].to_string(index=False) in list(distillation_OPEX.keys()):
            self.add_OPEX =\
                {f'{self.chemical_info["chemical"].to_string(index=False)}_distillation_OPEX':
                 self.electrolyte_flow_rate/1000*
                 distillation_OPEX[self.chemical_info['chemical'].to_string(index=False)]/
                 self.operating_days_per_year/24}