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
from biosteam.units.splitting import Splitter
from biosteam.units.design_tools import PressureVessel
from math import sqrt, pi, ceil

__all__ = (
    'BauxiteHammerMill',
    'ALFProduction',
    'PhaseChanger',
    'SolidPressureFilter',
    'ALFCrystallizer',
    'ReverseOsmosis',
    'S2WS',
    'ALFTSA',
    'CO2ElectrolyzerSystem'
    )

_mgd_to_cmh = auom('gallon').conversion_factor('m3')*1e6/24
_Pa_to_psi = auom('Pa').conversion_factor('psi')
_m_to_ft = auom('m').conversion_factor('ft')

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
# TODO: from Ben: the reactants are slurry, should use a batch reactor and do not need a crystallizer
class ALFProduction(bst.CSTR):
    '''
    Reactor for ALF production. See biosteam/units/stirred_tank_reactor.py.
    '''
    _N_ins = 1
    _N_outs = 2
    T_default = 60 + 273.15
    P_default = 101325
    tau_default = 2
    
    # TODO: add a parameter for X (conversion rate), but it seems that X=1 makes sense since Al(OH)3 is a solid and HCOOH is over amount
    # TODO: check the reaction between Fe2O3 and HCOOH
    def _setup(self):
        super()._setup()
        chemicals = self.chemicals
        self.ALF_production_AlH3O3 = bst.Reaction('AlH3O3 + 3HCOOH -> C3H3AlO6 + 3H2O', 'AlH3O3', 1, chemicals)
        self.ALF_production_bauxite = bst.Reaction('Al2O3 + 6HCOOH -> 2C3H3AlO6 + 3H2O', 'Al2O3', 1, chemicals)
        self.Fe_side_reaciton = bst.Reaction('Fe2O3 + 6HCOOH -> 2Fe + 3H2O + 3CO2', 'Fe2O3', 1, chemicals)
    
    def _run(self):
        vent, effluent = self.outs
        effluent.mix_from(self.ins, energy_balance=False)
        self.ALF_production_AlH3O3(effluent)
        self.ALF_production_bauxite(effluent)
        self.Fe_side_reaciton(effluent)
        effluent.T = vent.T = self.T
        effluent.P = vent.P = self.P
        vent.phase = 'g'
        vent.empty()
        vent.receive_vent(effluent, energy_balance=False)

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
    '''
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
    '''
    _units = {'Retentate flow rate':'kg/hr'}
    
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
        retentate.phases = permeate.phases = inlet.phase
        
        if len(inlet.phase) == 1:
            for chemical in self.split.keys():
                retentate.imass[chemical] = inlet.imass[chemical]*self.split[chemical]
                permeate.imass[chemical] = inlet.imass[chemical]*(1-self.split[chemical])
            retentate.imass['H2O'] = retentate.F_mass/(1-self.moisture_content)*self.moisture_content
            permeate.imass['H2O'] = inlet.imass['H2O'] - retentate.imass['H2O']
        else:
            for chemical in self.split.keys():
                chem = chemical if chemical.rfind('_') == -1 else chemical[:chemical.rfind('_')]
                if chemical[-2:] == '_s':
                    retentate.imass['s', chem] = inlet.imass['s', chem]*self.split[chemical]
                    permeate.imass['s', chem] = inlet.imass['s', chem]*(1-self.split[chemical])
                else:
                    retentate.imass['l', chem] = inlet.imass['l', chem]*self.split[chemical]
                    permeate.imass['l', chem] = inlet.imass['l', chem]*(1-self.split[chemical])
            retentate.imass['l','H2O'] = retentate.F_mass/(1-self.moisture_content)*self.moisture_content
            permeate.imass['l','H2O'] = inlet.imass['l','H2O'] - retentate.imass['l','H2O']
        
        retentate.T = permeate.T = inlet.T
        
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
# ReverseOsmosis
# =============================================================================
@cost(basis='Volumetric flow', ID='Reactor', units='m3/hr',
      cost=2450000, S=2.7*_mgd_to_cmh, CE=CEPCI[2012], n=1, BM=1.8)
@cost(basis='Volumetric flow', ID='Evaporator', units='m3/hr',
      kW=1103.636, cost=5000000, S=2.7*_mgd_to_cmh, CE=CEPCI[2012], n=0.6, BM=1.6)
class ReverseOsmosis(SanUnit):
    '''
    See biorefineries/wwt/_wwt_process.py for cost.
    See biosteam/wastewater/conventional.py for the default water recovery rate.
    '''
    _N_ins = 1
    _N_outs = 2
    _units = {'Volumetric flow':'m3/hr'}
    
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
    '''
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
    '''
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
class ALFTSA(PressureVessel, Splitter):
    '''
    TSA using ALF as adsorbent for CO2 adsorption.
    See biosteam/units/adsorption.py.
    
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
        Defaults to 60 $/ft3.
    adsorbent_lifetime : float
        The lifetime of adsorbents. Defaults to 20 years.
    adsorbate_ID : str
        Defaults to 'CO2'.
    split : dict[str, float] or list[float], optional
        Component splits towards the effluent (0th outlet).
        Defaults to dict(O2=0.9964, N2=0.9964, CO2=0.055).
    superficial_velocity : float, optional
        Superficial velocity of the feed. The diameter of the receiving vessel adjusts
        accordingly. Defaults to 2160 m/h. Typical velocities are 540 to 2160 m/h for liquids [1]_.
    regeneration_velocity : float, optional
        Mean velocity of the fluid used to regenerate the bed. Defaults to 1332 m/h.
    cycle_time : float, optional
        Time at which the receiving vessel is switched. Defaults to 5 h.
        Typical 2 h [biosteam/units/adsorption.py] - 8 h [2]_.
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
    treatment_capacity : float
        Flow rate of flue gas to a single TSA unit. Defaults to 10000 m3/h.
        10000 m3/h is an assumed value, at which the diameter of TSA columns is about the half of its length.
        1.74 m3/s for a temperature swing adsorption unit [4]_.
        100000 m3/h at 50 kPa vacuum for a vacuum swing adsorption unit [5]_.
    regen_CO2 : float
        CO2 ratio in the regeneration gas. Defaults to 0.975 (CO2 purity in carbon_dioxide).
    regen_N2 : float
        N2 ratio in the regeneration gas. Defaults to 0.025*0.82/0.87 (flue gas N2:O2 = 82:5).
    regen_O2 : float
        O2 ratio in the regeneration gas. Defaults to 0.025*0.05/0.87 (flue gas N2:O2 = 82:5).
    
    References
    ----------
    [1] Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
    [2] https://www.chemicalprocessing.com/processing-equipment/fluid-handling/\
        article/11302111/select-the-right-valves-for-adsorption-processes-chemical-processing
    [3] Evans, H. A.; Mullangi, D.; Deng, Z.; Wang, Y.; Peh, S. B.; Wei, F.;
        Wang, J.; Brown, C. M.; Zhao, D.; Canepa, P.; Cheetham, A. K.
        Aluminum Formate, Al(HCOO)3: An Earth-Abundant, Scalable, and Highly
        Selective Material for CO2 Capture. Science Advances 2022, 8 (44),
        eade1473. https://doi.org/10.1126/sciadv.ade1473.
    [4] Jareteg, A.; Maggiolo, D.; Larsson, A.; Thunman, H.; Sasic, S.; Ström,
        H. Industrial-Scale Benzene Adsorption: Assessment of a Baseline
        One-Dimensional Temperature Swing Model against Online Industrial Data.
        Ind. Eng. Chem. Res. 2020, 59 (26), 12239–12249.
        https://doi.org/10.1021/acs.iecr.0c01590.
    [5] Webley, P. A. Adsorption Technology for CO2 Separation and Capture:
        A Perspective. Adsorption 2014, 20 (2), 225–231.
        https://doi.org/10.1007/s10450-014-9603-2.
    [6] Ntiamoah, A.; Ling, J.; Xiao, P.; Webley, P. A.; Zhai, Y. CO2 Capture by
        Temperature Swing Adsorption: Use of Hot CO2-Rich Gas for Regeneration.
        Ind. Eng. Chem. Res. 2016, 55 (3), 703–713.
        https://doi.org/10.1021/acs.iecr.5b01384.
    '''
    auxiliary_unit_names = ('heat_exchanger_regeneration',)
    
    _N_ins = 2
    _N_outs = 2
    
    def _init(self, ID='ALF_TSA', ins=(), outs=(),
              adsorbent='ALF',
              adsorbent_cost=60,
              adsorbent_lifetime=20,
              adsorbate_ID='CO2',
              split=dict(O2=0.0036, N2=0.0036, CO2=0.945),
              superficial_velocity=2160,
              regeneration_velocity=1332,
              cycle_time=5,
              rho_adsorbent=1441,
              adsorbent_capacity=0.1188,
              T_regeneration=418,
              vessel_material='Stainless steel 316',
              vessel_type='Vertical',
              length_unused=1.219,
              treatment_capacity=10000,
              regen_CO2=0.975,
              regen_N2=0.025*0.82/0.87,
              regen_O2=0.025*0.05/0.87):
        bst.Splitter._init(self, split=split)
        self.adsorbent = adsorbent
        self.adsorbent_cost = adsorbent_cost
        self.adsorbent_lifetime = adsorbent_lifetime
        self._default_equipment_lifetime = self.equipment_lifetime = {'ALF': adsorbent_lifetime}
        self.adsorbate_ID = adsorbate_ID
        self.superficial_velocity = superficial_velocity
        self.regeneration_velocity = regeneration_velocity
        self.cycle_time = cycle_time
        self.rho_adsorbent = rho_adsorbent
        self.adsorbent_capacity = adsorbent_capacity
        self.T_regeneration = T_regeneration
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.length_unused = length_unused
        self.treatment_capacity=treatment_capacity
        self.regen_CO2 = regen_CO2
        self.regen_N2 = regen_N2
        self.regen_O2 = regen_O2
        self.heat_exchanger_regeneration = bst.HXutility(None, None, None, thermo=self.thermo)
    
    def _run(self):
        feed, regen = self.ins
        carbon_dioxide, effluent = self.outs 
        regen.empty()
        for i in self.outs: i.empty()
        
        self.N = ceil(feed.F_vol/self.treatment_capacity)

        feed.split_to(carbon_dioxide, effluent, self.split)
        F_vol_feed = feed.F_vol/self.N
        
        superficial_velocity = self.superficial_velocity
        adsorbate_ID = self.adsorbate_ID
        
        F_mass_adsorbate = carbon_dioxide.imass[adsorbate_ID]/self.N
        
        self.diameter = diameter = 2*sqrt(F_vol_feed/(superficial_velocity*pi))
        self.area = area = pi*diameter*diameter/4
        # length of equilibrium section plus unused bed
        total_length = self.cycle_time*F_mass_adsorbate/\
            (self.adsorbent_capacity*self.rho_adsorbent*area) + self.length_unused
        # size of each column
        self.length = length = total_length/2
        self.vessel_volume = length*area
        T_original = regen.T
        # the regen should have the same gas composition as carbon_dioxide, see [6]_
        # carbon dioxide: 97.5% CO2, the rest should be 0.82/0.87 N2 and 0.05/0.87 O2 (see flue gas composition)
        regen.reset_flow(CO2=self.regen_CO2, N2=self.regen_N2, O2=self.regen_O2, phase='g', units='kg/hr')
        carbon_dioxide.T = regen.T = self.T_regeneration
        regen.F_vol = area*self.regeneration_velocity*self.N
        regen.T = T_original
    
    def _design(self):
        feed, regen = self.ins
        design_results = self.design_results
        diameter = self.diameter
        length = self.length
        design_results['Number of sets'] = self.N
        design_results['Number of reactors'] = 3
        design_results.update(self._vessel_design(feed.P*_Pa_to_psi,
                                                  diameter*_m_to_ft,
                                                  length*_m_to_ft))
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
        baseline_purchase_costs.update(self._vessel_purchase_cost(design_results['Weight'],
                                                                  design_results['Diameter'],
                                                                  design_results['Length']))
        N_sets = design_results['Number of sets']
        N_reactors = design_results['Number of reactors']
        for i, j in baseline_purchase_costs.items():
            baseline_purchase_costs[i] *= N_sets*N_reactors
        # 35.3147 ft3/m3
        baseline_purchase_costs[self.adsorbent] = N_sets*N_reactors*35.3147*\
            self.vessel_volume*self.adsorbent_cost

# =============================================================================
# CO2ElectrolyzerSystem
# =============================================================================
# cost basis year: 2010 (may be later years, but to be conservative)
# see Jouny et al. 2018 SI Excel - tab Data - cell C4
@cost(basis='Electrolyzer area', ID='Electrolyzer', units='m2',
      cost=250.25*1.75*175/1000*10*(1+1/0.65*0.35), S=1, CE=qs.CEPCI_by_year[2010], n=1, BM=1.2)
@cost(basis='Electrolyte flow rate (ethanol)', ID='Distiller (ethanol)', units='L/min',
      cost=4162240, S=1000, CE=qs.CEPCI_by_year[2010], n=0.7)
@cost(basis='Electrolyte flow rate (formic acid)', ID='Distiller (formic acid)', units='L/min',
      cost=6896190, S=1000, CE=qs.CEPCI_by_year[2010], n=0.7)
@cost(basis='Electrolyte flow rate (methanol)', ID='Distiller (methanol)', units='L/min',
      cost=4514670, S=1000, CE=qs.CEPCI_by_year[2010], n=0.7)
@cost(basis='Electrolyte flow rate (propanol)', ID='Distiller (propanol)', units='L/min',
      cost=4687910, S=1000, CE=qs.CEPCI_by_year[2010], n=0.7)
@cost(basis='Total gas flow for PSA', ID='PSA', units='m3/h',
      cost=1989043, S=1000, CE=qs.CEPCI_by_year[2010], n=0.7)
@cost(basis='Gas product and hydrogen flow for PSA', ID='H2_PSA', units='m3/h',
      cost=1989043, S=1000, CE=qs.CEPCI_by_year[2010], n=0.7)
class CO2ElectrolyzerSystem(SanUnit):
    '''
    CO2 electrolyzer system that converts CO2 into reduced 1C, 2C, and 3C products [1]_.
    Compared to [1]_, this system adds a PSA (if gas products) to separate H2.
    Assume the CAPEX and OPEX of the PSA for H2 separation are the same as the PSA for
    CO2 separation.
    
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
    conversion : float, optional
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
    _N_outs = 3
    
    _units= {'Electrolyzer area':'m2',
             'Electrolyte flow rate (ethanol)':'L/min',
             'Electrolyte flow rate (formic acid)':'L/min',
             'Electrolyte flow rate (methanol)':'L/min',
             'Electrolyte flow rate (propanol)':'L/min',
             'Total gas flow for PSA':'m3/h',
             'Gas product and hydrogen flow for PSA':'m3/h'}
    
    def __init__(self, ID='', ins=(), outs=(), thermo=None, target_product='formic acid',
                 current_density=0.2, cathodic_overpotential=0.454, cell_voltage=2.3,
                 product_selectivity=0.9, conversion=0.5, operating_days_per_year=350):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo)
        self.target_product = target_product
        self.current_density = current_density
        self.cathodic_overpotential = cathodic_overpotential
        self.cell_voltage = cell_voltage
        self.product_selectivity = product_selectivity
        self.conversion = conversion
        self.operating_days_per_year = operating_days_per_year
    
    def _run(self):
        carbon_dioxide, water = self.ins
        product, hydrogen, offgas = self.outs
        
        # propanol is n-propanol
        if self.target_product not in ['carbon monoxide','ethanol','ethylene',
                                       'formic acid','methane','methanol','propanol']:
            raise ValueError("target product must be in 'carbon monoxide', 'ethanol', " + 
                             "'ethylene', 'formic acid', 'methane', 'methanol', and 'propanol'")
        
        # density in kg/m3
        product_info = {'chemical': ['carbon monoxide','ethanol','ethylene',
                                     'formic acid','methane','methanol','propanol'],
                        'formula': ['CO','C2H6O','C2H4','HCOOH','CH4','CH4O','C3H8O'],
                        'elec_number': [2, 12, 12, 2, 8, 6, 18],
                        'mole_ratio': [1, 2, 2, 1, 1, 1, 3],
                        'MW': [28, 46, 28, 46, 16, 32, 60],
                        'density': [1.14, 789, 1.18, 1221, 0.656, 792, 803],
                        'state': ['gas','liq','gas','liq','gas','liq','liq'],
                        'potential': [-0.106, 0.084, 0.064, -0.25, 0.169, 0.016, 0.095],
                        'cathode_water_production': [1, 3, 4, 0, 2, 1, 5]}
        
        product_info = pd.DataFrame.from_dict(product_info)
        
        chemical_info = self.chemical_info =\
            product_info[product_info['chemical'] == self.target_product]
        
        # CO2 inlet flow rate [kg/h]
        CO2_inlet_flow_rate = carbon_dioxide.imass['CO2']
        
        # unconverted CO2 will be recycled
        CO2_converted = CO2_inlet_flow_rate*24
        
        # CO2 recycle flow rate [kg/h]
        CO2_recycle_flow_rate_kg_per_h = CO2_inlet_flow_rate*(1-self.conversion)/self.conversion
        
        # CO2 recycle flow rate [m3/h]
        # the density of CO2 is 1.98 kg/m3, Jouny et al. 2018
        self.CO2_recycle_flow_rate_m3_per_h = CO2_recycle_flow_rate_kg_per_h/1.98
        
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
        self.gas_product_flow_rate = 0\
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
        self.total_gas_flow = self.CO2_recycle_flow_rate_m3_per_h +\
            self.gas_product_flow_rate + hydrogen_flow_rate_m3_per_h
        
        product.empty()
        product.imass[chemical_info['formula'].to_string(index=False)] = product_production/24
        product.phase = 'g' if chemical_info['state'].to_string(index=False) == 'gas' else 'l'
        
        # note Jouny et al. 2018 calculated process water amount as 'current_needed/4/96485*18/1000*3600'
        # change it to 'current_needed/2/96485*18/1000*3600' as every 18 g H2O can transfer 2 e-
        water.imass['H2O'] = current_needed/2/96485*18/1000*3600
        
        # we added a PSA to separate H2
        # hydrogen density: 0.08375 kg/m3, Jouny et al. 2018
        hydrogen.empty()
        hydrogen.imass['H2'] = hydrogen_flow_rate_m3_per_h*0.08375
        
        offgas.empty()
        offgas.imass['N2'] = carbon_dioxide.imass['N2']
        # H2O after distillation or PSA
        # note formic acid is produced in the cathode with no water produced
        # but still needs distillation since there must be carried-out water
        offgas.imass['H2O'] = product.F_mass/float(chemical_info['MW'])*\
            float(chemical_info['cathode_water_production'])*18
        
        offgas.imass['O2'] = carbon_dioxide.F_mass + water.F_mass -\
            product.F_mass - hydrogen.F_mass - offgas.F_mass
        offgas.phase = 'g'
    
    def _design(self):
        D = self.design_results
        D['Electrolyzer area'] = self.electrolyzer_area
        
        for liquid_product in ['ethanol','formic acid','methanol','propanol']:
            D[f'Electrolyte flow rate ({liquid_product})'] = self.electrolyte_flow_rate\
                if self.chemical_info['chemical'].to_string(index=False) == liquid_product\
                    else 0
        
        D['Total gas flow for PSA'] = 0\
            if self.CO2_recycle_flow_rate_m3_per_h == self.total_gas_flow\
                else self.total_gas_flow
        
        D['Gas product and hydrogen flow for PSA'] = 0\
            if self.gas_product_flow_rate == 0\
            else self.total_gas_flow - self.CO2_recycle_flow_rate_m3_per_h
        
        # PSA operating power: 0.25 kWh/m3, Jouny et al. 2018
        self.add_power_utility(self.power_needed*1000 +\
                               D['Total gas flow for PSA']*0.25 +\
                               D['Gas product and hydrogen flow for PSA']*0.25)
    
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