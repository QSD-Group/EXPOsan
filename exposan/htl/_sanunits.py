#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs
from math import ceil, pi, log
from warnings import warn
from qsdsan import SanUnit, Stream, Construction
from qsdsan.sanunits import Reactor, HXutility, Tank
from qsdsan.utils import auom
from biosteam.units import StorageTank as BSTStorageTank
from biosteam.units.decorators import cost
from biosteam.units.design_tools import CEPCI_by_year, size_batch, flash_vessel_design
from biosteam.units.design_tools.specification_factors import material_densities_lb_per_ft3
from biosteam.exceptions import bounds_warning, DesignWarning

__all__ = (
    'AmineAbsorption',
    'AcidExtraction',
    'BiocrudeTank',
    'DAPSynthesis',
    'FuelMixer',
    'HTLaqueous',
    'HTLmixer',
    'Humidifier',
    'PreStripper',
    'StreamTypeConverter',
    'StruvitePrecipitation',
    'UANSynthesis',
    'UreaSynthesis',
    'WWmixer',
    'WWTP',
    )

_hp2kW = 0.7457
_lb_to_kg = auom('lb').conversion_factor('kg')
_m3perh_to_MGD = auom('m3/h').conversion_factor('MGD')
_Pa_to_psi = auom('Pa').conversion_factor('psi')
_m_to_ft = auom('m').conversion_factor('ft')

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
    acid_vol : float
        0.5 M H2SO4 to hydrochar ratio: mL/g.
    P_acid_recovery_ratio : float
        The ratio of phosphorus that can be extracted.
        
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
                 init_with='WasteStream', acid_vol=7, P_acid_recovery_ratio=0.8,
                 P=None, tau=2, V_wf=0.8, # tau: [1]
                 length_to_diameter=2, N=1, V=10, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0, # use MixTank default value
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 304', # acid condition
                 vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
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
        self.N = ceil(self.HTL.WWTP.ins[0].F_vol/788.627455/self.V)
        # 1/788.627455 m3 reactor/m3 wastewater/h (50 MGD ~ 10 m3)
        self.P = self.ins[1].P
        Reactor._design(self)

# =============================================================================
# AmineAbsorption 
# =============================================================================

@cost(basis='CO2 flow', ID='Reactor', units='kmol/hr',
      cost=3.063e6, S=760*1000000/24/44/1000, CE=CEPCI_by_year[2019], n=0.67, BM=1.47)
@cost(basis='CO2 flow', ID='Pumps electricity', units='kmol/hr',
      kW=55595.96/(613*(1000/44))*(24123*0.1186), S=24123, CE=CEPCI_by_year[2009])
class AmineAbsorption(SanUnit):
    '''
    Similar to biosteam.units.AmineAbsorption, but can init with 'WasteStream'.
    Electricity use is from biosteam.units.AmineAbsorption.
    Cost function from [1], only keep pumps electricity from biosteam.units.AmineAbsorption.
    
    Parameters
    ----------
    ins : iterable
        flue_gas, MEA, water.
    outs : iterable
        vent, CO2.
    CO2_recovery : float
        Percentage of CO2 that can be captured.
    MEA_to_CO2 : float
        Net usage of MEA (kg pure MEA/metric tonne CO2 captured).
        The default is 1.5 based on [1]_ and [3]_.
    
    References
    ----------
    .. [1] Devkota, S.; Karmacharya, P.; Maharjan, S.; Khatiwada, D.;
        Uprety, B. Decarbonizing Urea: Techno-Economic and Environmental
        Analysis of a Model Hydroelectricity and Carbon Capture Based
        Green Urea Production. Applied Energy 2024, 372, 123789.
        https://doi.org/10.1016/j.apenergy.2024.123789.
    '''
    
    auxiliary_unit_names=('heat_ex_heating_1','heat_ex_heating_2')
    
    _N_ins = 3
    _N_outs = 2
    _units = {'Total flow':'kmol/hr',
              'CO2 flow':'kmol/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', CO2_recovery=0.9, MEA_to_CO2=1.5):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.CO2_recovery = CO2_recovery
        self.MEA_to_CO2 = MEA_to_CO2
        hx_ht_in_1 = Stream(f'{ID}_hx_ht_in_1')
        hx_ht_out_1 = Stream(f'{ID}_hx_ht_out_1')
        self.heat_ex_heating_1 = HXutility(ID=f'.{ID}_hx_ht_1', ins=hx_ht_in_1, outs=hx_ht_out_1, T=273.15+40, rigorous=True)
        hx_ht_in_2 = Stream(f'{ID}_hx_ht_in_2')
        hx_ht_out_2 = Stream(f'{ID}_hx_ht_out_2')
        self.heat_ex_heating_2 = HXutility(ID=f'.{ID}_hx_ht_2', ins=hx_ht_in_2, outs=hx_ht_out_2, T=273.15+40, rigorous=True)
    
    def _run(self):
        flue_gas, MEA, water = self.ins
        vent, CO2 = self.outs
        vent.copy_like(flue_gas)
        CO2.imol['CO2'] = self.CO2_recovery * flue_gas.imol['CO2']
        vent.imol['CO2'] = flue_gas.imol['CO2'] - CO2.imol['CO2']
        MEA.imass['MEA'] = self.MEA_to_CO2 * CO2.imass['CO2']/1000
        water.imass['Water'] = MEA.imass['MEA']/0.3*(1-0.3)
        vent.T = CO2.T = 273.15 + 40
        CO2.phase = 'g'
        
    def _design(self):
        self.design_results['Total flow'] = self.ins[0].F_mol
        self.design_results['CO2 flow'] = self.outs[1].F_mol
        
        hx_ht_1 = self.heat_ex_heating_1
        hx_ht_1_ins0, hx_ht_1_outs0 = hx_ht_1.ins[0], hx_ht_1.outs[0]
        hx_ht_1_ins0.copy_like(self.outs[0])
        hx_ht_1_outs0.copy_like(hx_ht_1_ins0)
        # the tempeture of CHP exhaust is probably to be higher, assume the temperature to be 298.15 K to be conservative
        # set hxn_ok=False to be conservative
        hx_ht_1_ins0.T = 298.15
        hx_ht_1_outs0.T = hx_ht_1.T
        
        hx_ht_1_ins0.vle(T=hx_ht_1_ins0.T, P=hx_ht_1_ins0.P)
        hx_ht_1_outs0.vle(T=hx_ht_1_outs0.T, P=hx_ht_1_outs0.P)
        # also set hxn_ok=False to be conservative
        hx_ht_1.simulate_as_auxiliary_exchanger(ins=hx_ht_1.ins, outs=hx_ht_1.outs, vle=True, hxn_ok=False)
            
        hx_ht_2 = self.heat_ex_heating_2
        hx_ht_2_ins0, hx_ht_2_outs0 = hx_ht_2.ins[0], hx_ht_2.outs[0]
        hx_ht_2_ins0.copy_like(self.outs[1])
        hx_ht_2_outs0.copy_like(hx_ht_2_ins0)
        # the tempeture of CHP exhaust is probably to be higher, assume the temperature to be 298.15 K to be conservative
        hx_ht_2_ins0.T = 298.15
        hx_ht_2_outs0.T = hx_ht_2.T
        
        hx_ht_2_ins0.vle(T=hx_ht_2_ins0.T, P=hx_ht_2_ins0.P)
        hx_ht_2_outs0.vle(T=hx_ht_2_outs0.T, P=hx_ht_2_outs0.P)
        # set hxn_ok=False to be conservative
        hx_ht_2.simulate_as_auxiliary_exchanger(ins=hx_ht_2.ins, outs=hx_ht_2.outs, vle=True, hxn_ok=False)

# =============================================================================
# BiocrudeTank 
# =============================================================================

class BiocrudeTank(Tank, BSTStorageTank):
    '''
    Similar to the :class:`biosteam.units.MixTank`, but can calculate material usage.
    Also add three parameters (crude_oil_density, crude_oil_HHV, and biocrude_wet_density)
    to correct biocrude price and CI.
    
    See Also
    --------
    :class:`biosteam.units.StorageTank`
    
    Parameters
    ----------
    crude_oil_density : float
        Density of crude oil, kg/m3.
    crude_oil_HHV : float
        HHV of crude oil, MJ/kg.
    biocrude_wet_density : float
        Density of biocrude, kg/m3.
    biocrude_distance : float
        Distance between WRRFs and oil refineries, [km].
    
    References
    ----------
    .. [1] https://www.transmountain.com/about-petroleum-liquids (accessed 2025-02-05).
    .. [2] https://world-nuclear.org/information-library/facts-and-figures/heat-values-of-various-fuels
    .. [3] Snowden-Swan, L. J.; Li, S.; Thorson, M. R.; Schmidt, A. J.; Cronin, D. J.;
        Zhu, Y.; Hart, T. R.; Santosa, D. M.; Fox, S. P.; Lemmon, T. L.; Swita, M. S.
        Wet Waste Hydrothermal Liquefaction and Biocrude Upgrading to Hydrocarbon Fuels:
        2022 State of Technology; PNNL-33622; Pacific Northwest National Lab. (PNNL),
        Richland, WA (United States), 2022. https://doi.org/10.2172/1897670.        
    '''
    
    _units = {'Diameter': 'ft',
              'Length': 'ft',
              'Wall thickness': 'in',
              'Weight': 'lb'}
    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Horizontal vessel diameter': (3, 21),
               'Vertical vessel length': (12, 40)}
    _vessel_material = 'Stainless steel'
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 # [1]
                 crude_oil_density=850,
                 # [2]
                 crude_oil_HHV=44.5,
                 # [3]
                 biocrude_wet_density=983,
                 biocrude_distance=100,
                 vessel_type=None, tau=None, V_wf=None,
                 vessel_material=None, kW_per_m3=0.,
                 init_with='WasteStream', F_BM_default=None,
                 include_construction=True, length_to_diameter=2):
        Tank.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo,
                      init_with=init_with, F_BM_default=F_BM_default,
                      include_construction=include_construction,
                      vessel_type=vessel_type, tau=tau, V_wf=V_wf,
                      vessel_material=vessel_material, kW_per_m3=kW_per_m3)
        self.crude_oil_density = crude_oil_density
        self.crude_oil_HHV = crude_oil_HHV
        self.biocrude_wet_density = biocrude_wet_density
        self.biocrude_distance = biocrude_distance
        self.length_to_diameter = length_to_diameter
    
    def _init_lca(self):
        item_name = self.vessel_material.replace(' ', '_')
        self.construction = [
            Construction(item_name.lower(), linked_unit=self, item=item_name, quantity_unit='kg'),
            ]
    
    def _design(self):
        BSTStorageTank._design(self)
        D = self.design_results
                
        Diameter = (4*D['Total volume']/pi/self.length_to_diameter)**(1/3)
        Diameter *= _m_to_ft # convert from m to ft
        L = Diameter * self.length_to_diameter # ft
        D.update(self._horizontal_vessel_design(self.ins[0].P*_Pa_to_psi, Diameter, L))
        D['Material'] = self.vessel_material
        if self.include_construction: self.construction[0].quantity = D['Weight']*_lb_to_kg
    
    def _horizontal_vessel_design(self, pressure, diameter, length) -> dict:
        pressure = pressure
        diameter = diameter
        length = length
        # Calculate vessel weight and wall thickness
        if self.vessel_material == 'Carbon steel':
            rho_M = material_densities_lb_per_ft3[self.vessel_material]
        else:
            rho_M = material_densities_lb_per_ft3['Stainless steel 304']
        if pressure < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=DesignWarning)
        VW, VWT = flash_vessel_design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M)
        bounds_warning(self, 'Horizontal vessel weight', VW, 'lb',
                       self._bounds['Horizontal vessel weight'], 'cost')
        bounds_warning(self, 'Horizontal vessel diameter', diameter, 'ft',
                       self._bounds['Horizontal vessel diameter'], 'cost')
        Design = {}
        Design['Vessel type'] = 'Horizontal'
        Design['Length'] = length # ft
        Design['Diameter'] = diameter # ft
        Design['Weight'] = VW # lb
        Design['Wall thickness'] = VWT # in
        return Design
    
    @property
    def vessel_material(self):
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, i):
        exist_material = getattr(self, '_vessel_material', None)
        BSTStorageTank.vessel_material.fset(self, i)
        if i and exist_material == i: return # type doesn't change, no need to reload construction items
        self._init_lca()

# =============================================================================
# DAPSynthesis
# =============================================================================

@cost('Crystallizer volume', 'Crystallizer',
      CE=444., S=0.003785411784, # originally 1 gal
      BM=2.0, N='Number of crystallizers',
      f=lambda S: 222.4 * S**0.71 + 35150)
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
class DAPSynthesis(Reactor):
    '''
    Synthesize DAP followed by crystallization and filtration.
    If ammonia is excess, additional ammonia is calculated as excess_ammonia
    (as if it never enters the reactor).
    
    Parameters
    ----------
    ins : iterable
        P_solution, ammonia.
    outs : iterable
        DAP, excess_ammonia, effluent.
    P_pre_recovery_ratio : float
        Ratio of phosphorus that can be precipitated out.
    crystallizer_electricity : float
        Electricity usage per volume in kW/gal. Defaults to 0.00746, a 
        heuristic value for suspension of solids. Note the unit should be
        kW/m3, not a driver of the cost or CI, keep consistent with
        BioSTEAM for now
    '''
    _N_ins = 3
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default}
    
    _units= {'Crystallizer volume':'m3',
             'Batch time':'hr',
             'Loading time':'hr',
             'Retentate flow rate':'kg/h'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', 
                  # due to lack of data, assume P_syn_recovery_ratio to be the same as
                  # P_pre_recovery_ratio in StruvitePrecipitation
                  P_syn_recovery_ratio=0.828,
                  crystallizer_electricity=0.00746,
                  # due to lack of data, assume tau to be the same as
                  # tau in StruvitePrecipitation
                  P=None, tau=1,
                  V_wf=0.8,length_to_diameter=2, N=1, V=None,
                  auxiliary=False, mixing_intensity=None,
                  kW_per_m3=0, wall_thickness_factor=1,
                  vessel_material='Stainless steel 316',
                  vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.P_syn_recovery_ratio = P_syn_recovery_ratio
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
        self.crystallizer_electricity = crystallizer_electricity
    
    def _run(self):
        
        P_solution, ammonia = self.ins
        DAP, excess_ammonia, effluent = self.outs
        
        DAP.imass['DAP'] = P_solution.imass['P']*self.P_syn_recovery_ratio/30.97*132.06
        DAP.phase = 's'
        
        excess_ammonia.imass['NH3'] = ammonia.imass['NH3'] - DAP.imass['DAP']/132.06*2*17.031
        excess_ammonia.imass['H2O'] = ammonia.imass['H2O']
        excess_ammonia.T = ammonia.T
        excess_ammonia.phase = 'g'
        
        effluent.imass['H2O'] = P_solution.F_mass + ammonia.F_mass - DAP.F_mass - excess_ammonia.F_mass
    
    def _design(self):
        # reactor
        self.N = 1
        self.V = (self.outs[0].F_vol + self.outs[2].F_vol)/self.V_wf
        self.P = self.ins[0].P
        Reactor._design(self)
        
        Design = self.design_results
        
        # crystallizer
        crystallizer_total_volume = self.outs[0].F_vol + self.outs[2].F_vol
        # from TAL.units in BioSTEAM: assumed 8 h; uncertainty range is 2-14 h
        crystallizer_tau = 8
        # cleaning and unloading time
        crystallizer_tau_0 = 1
        # assume the same number of DAP synthesizer and crystallizer
        crystallizer_N = max(2, self.N)
        # fraction of filled tank to total tank volume
        crystallizer_V_wf = 0.9
        
        dct = size_batch(crystallizer_total_volume,
                         crystallizer_tau,
                         crystallizer_tau_0,
                         crystallizer_N,
                         crystallizer_V_wf)
        
        Design['Crystallizer volume'] = volume = dct.pop('Reactor volume')
        Design.update(dct)
        Design['Number of crystallizers'] = crystallizer_N
        
        self.add_heat_utility(self.Hnet, self.outs[2].T)
        self.add_power_utility(self.crystallizer_electricity * crystallizer_V_wf * volume * crystallizer_N)
        
        # pressure filter
        Design['Retentate flow rate'] = self.outs[0].F_mass

# =============================================================================
# FuelMixer
# =============================================================================

class FuelMixer(SanUnit):
    '''
    Convert gasoline to diesel or diesel to gasoline based on LHV.
    
    Parameters
    ----------
    ins: iterable
        gasoline, diesel.
    outs: iterable
        fuel.
    target : str
        The target can only be 'gasoline' or 'diesel'.
    gasoline_price : float
        Gasoline price, [$/kg].
    diesel_price : float
        Diesel price, [$/kg].
    '''
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', target='diesel',
                 gasoline_gal_2_kg=2.834894885,
                 diesel_gal_2_kg=3.220628346,
                 gasoline_price=0.9388,
                 diesel_price=0.9722):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target = target
        self.gasoline_gal_2_kg = gasoline_gal_2_kg
        self.diesel_gal_2_kg = diesel_gal_2_kg
        self.gasoline_price = gasoline_price
        self.diesel_price = diesel_price

    def _run(self):
        
        gasoline, diesel = self.ins
        fuel = self.outs[0]
        target = self.target
        
        gasoline_LHV_2_diesel_LHV = (gasoline.LHV/gasoline.F_mass)/(diesel.LHV/diesel.F_mass)
        # KJ/kg gasoline:KJ/kg diesel
        
        if target == 'gasoline':
            fuel.imass['Gasoline'] = gasoline.F_mass + diesel.F_mass/gasoline_LHV_2_diesel_LHV
            fuel.T = gasoline.T
            fuel.P = gasoline.P
        elif target == 'diesel':
            fuel.imass['Diesel'] = diesel.F_mass + gasoline.F_mass*gasoline_LHV_2_diesel_LHV
            fuel.T = diesel.T
            fuel.P = diesel.P
    
    def _cost(self):
        if self.target == 'gasoline':
            self.outs[0].price = self.gasoline_price
        elif self.target == 'diesel':
            self.outs[0].price = self.diesel_price
            
    @property
    def target(self):
        return self._target
    @target.setter
    def target(self, i):
        if i not in ('gasoline', 'diesel'):
            raise ValueError('`target` must be either "gasoline" or "diesel" ',
                             f'the provided "{i}" is not valid.')
        self._target = i

# =============================================================================
# HTLaqueous
# =============================================================================

class HTLaqueous(SanUnit):
    '''
    A fake unit that calculates C, N, P, and H2O amount in the HTL aqueous effluent.
    
    Parameters
    ----------
    ins : iterable
        HTL_effluent_undefined.
    outs : iterable
        HTL_effluent_defined.
        
    References
    ----------
    .. [1] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
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
        
        HTL_effluent_undefined = self.ins[0]
        HTL_effluent_defined = self.outs[0]
        
        HTL_effluent_defined.copy_like(HTL_effluent_undefined)
        HTL_effluent_defined.empty()
        
        HTL_effluent_defined.imass['C'] = self.ins[0]._source.HTLaqueous_C
        HTL_effluent_defined.imass['N'] = self.ins[0]._source.HTLaqueous_N
        HTL_effluent_defined.imass['P'] = self.ins[0]._source.HTLaqueous_P
        # other compositions represented by H2O except C, N, P
        HTL_effluent_defined.imass['H2O'] = HTL_effluent_undefined.F_mass -\
                                            HTL_effluent_defined.imass['C'] -\
                                            HTL_effluent_defined.imass['N'] -\
                                            HTL_effluent_defined.imass['P']
    
    @property
    def pH(self):
        return 7

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
        
    References
    ----------
    .. [1] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
        https://doi.org/10.1021/acs.est.8b04035.
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

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
# Humidifier
# =============================================================================

class Humidifier(SanUnit):
    '''
    A fake unit increases the moisture content of HTL feedstocks to 80%.
    Assume 80% of the water can be recovered from the end. Calcalute makeup water.
        feedstock: X H2O + Y dry matter
        target:    4Y H2O + Y dry matter
        water recovery = 0.8 = (4Y - X - makeup water)/4Y
        makeup water = 0.8Y - X
    
    Parameters
    ----------
    ins : iterable
        feedstock, makeup, recycle.
    outs : iterable
        mixture.
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

    _N_ins = 3
    _N_outs = 1
        
    def _run(self):
        
        feedstock, makeup, recycle = self.ins
        mixture = self.outs[0]
        
        makeup.imass['H2O'] = max(0, 0.8*(feedstock.F_mass - feedstock.imass['H2O']) - feedstock.imass['H2O'])
        
        recycle.imass['H2O'] = (feedstock.F_mass - feedstock.imass['H2O'])/0.2 - feedstock.F_mass - makeup.imass['H2O']

        mixture.mix_from(self.ins)

# =============================================================================
# PreStripper
# =============================================================================

class PreStripper(SanUnit):
    '''
    Calculate the NH3 concentration in the influent to the stripper.
    
    Parameters
    ----------
    ins : iterable
        influent.
    outs : iterable
        effluent.
    influent_pH : float
        CHG effluent pH: 8.16 ± 0.25 [1]
    target_pH : float
        2 unit higher than pKa (9.25 + 2 = 11.25)
    
    References
    ----------
    [1] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
        https://doi.org/10.1021/acs.est.8b04035.
    '''
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 influent_pH=8.16, # CHG effluent pH: 8.16 ± 0.25 [1]
                 target_pH=11.25): # 2 unit higher than pKa (9.25)
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.influent_pH = influent_pH
        self.target_pH = target_pH
    
    def _run(self):
        influent, base = self.ins
        effluent = self.outs[0]
        
        NaOH_conc = 10**(self.target_pH - 14) - 10**(self.influent_pH - 14)
        NaOH_mol = NaOH_conc*self.ins[0].F_mass
        base.imass['NaOH'] = NaOH_mol*39.997/1000
        
        self.CHG = self.ins[0]._source.ins[0]._source.ins[0]._source
        
        effluent.imass['NH3'] = self.CHG.CHGout_N/14.0067*17.031
        effluent.imass['H2O'] = influent.F_mass + base.F_mass - effluent.imass['NH3']

# =============================================================================
# StreamTypeConverter
# =============================================================================
class StreamTypeConverter(SanUnit):
    '''
    A fake unit that converts Stream or MultiStream to WasteStream to enable LCA.
    
    Parameters
    ----------
    ins : iterable
        inlet.
    outs : iterable
        outlet.
    '''
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream'):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with)

    def _run(self):
        inlet = self.ins[0]
        outlet = self.outs[0]
        
        try:
            outlet.copy_like(inlet)
        except TypeError:     
            for item in qs.get_components():
                if inlet.imass[item.ID] != 0:
                    outlet.imass[item.ID] = inlet.imass[item.ID]
                    outlet.T = inlet.T
                    outlet.P = inlet.P
                    outlet.phase = inlet.phase

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
    target_pH : float
        Target pH for struvite precipitation.
    Mg_P_ratio : float
        mol(Mg) to mol(P) ratio.   
    P_pre_recovery_ratio : float
        Ratio of phosphorus that can be precipitated out.
    HTLaqueous_NH3_N_2_total_N : float
        Ratio of NH3-N to TN in HTL aqueous phase.
        
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
                  init_with='WasteStream', 
                  target_pH = 9,
                  Mg_P_ratio=1,
                  P_pre_recovery_ratio=0.828, # [1]
                  HTLaqueous_NH3_N_2_total_N = 0.853, # [2]
                  P=None, tau=1, # [1]
                  V_wf=0.8, length_to_diameter=2, N=1, V=20,
                  auxiliary=False, mixing_intensity=None,
                  kW_per_m3=0, # use MixTank default value
                  wall_thickness_factor=1,
                  vessel_material='Carbon steel', # basic condition
                  vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
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
        # 2/788.627455 m3 reactor/m3 wastewater/h (50 MGD ~ 20 m3)
        self.N = ceil(self.HTLmixer.ins[0]._source.WWTP.ins[0].F_vol*2/788.627455/self.V)
        self.P = self.ins[0].P
        Reactor._design(self)

# =============================================================================
# UANSynthesis
# =============================================================================

# EURO to USD: [1]
@cost(basis='Urea production capacity', ID='Urea synthesizer', units='kg/h',
      cost=28/0.951*1000000, S=500*1000/24,
      CE=CEPCI_by_year[2016], n=0.67, BM=1.47)
class UANSynthesis(Reactor):
    '''
    Urea ammonium nitrate 30 (UAN30) synthesis.
    Cost of urea systhesis from [2].
    
    Parameters
    ----------
    ins : iterable
        ammonia, carbon_dioxide, additional_carbon_dioxide, nitric_acid, water.
    outs : iterable
        UAN_solution, vapor, waste.
    ratio : float
        The overall ratio between NH3 and CO2 as reactants (after considering
        recycling of unconverted reactants).
    efficiency : float
        The overall conversion efficiency (after considering recycling of
        unconverted reactants) of CO2 to urea.
    loss : float
        The loss ratio of unconverted reactants before recycling.
    
    References
    ----------
    .. [1] https://www.irs.gov/individuals/international-taxpayers/yearly-\
        average-currency-exchange-rates (accessed 2025-02-09)
    .. [2] Devkota, S.; Karmacharya, P.; Maharjan, S.; Khatiwada, D.;
        Uprety, B. Decarbonizing Urea: Techno-Economic and Environmental
        Analysis of a Model Hydroelectricity and Carbon Capture Based
        Green Urea Production. Applied Energy 2024, 372, 123789.
        https://doi.org/10.1016/j.apenergy.2024.123789.
    .. [3] https://www.cropnutrition.com/resource-library/urea-ammonium-nitrate/
        (accessed 2025-02-09)
    .. [4] Palys, M. J.; Daoutidis, P. Techno-Economic Optimization of
        Renewable Urea Production for Sustainable Agriculture and CO2 
        tilization. J. Phys. Energy 2023, 6 (1), 015013.
        https://doi.org/10.1088/2515-7655/ad0ee6.
    '''
    _N_ins = 5
    _N_outs = 3
    _F_BM_default = {**Reactor._F_BM_default}
    _units= {'Urea production capacity':'kg/h'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 ratio=3.5, efficiency=0.8, loss=0.02,
                 P=None, tau=1, V_wf=0.8, # assume to be the same as StruvitePrecipitation
                 length_to_diameter=2, N=1, V=20, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316', # basic condition
                 vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.ratio = ratio
        self.efficiency = efficiency
        self.loss = loss
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
        
        ammonia, carbon_dioxide, additional_carbon_dioxide, nitric_acid, water = self.ins
        UAN_solution, vapor, waste = self.outs
        
        NH3_CO2_molar_ratio = (self.ratio - (self.ratio - 2*self.efficiency)*(1-self.loss))/\
                              (1 - (1 - self.efficiency)*(1-self.loss))
        
        # UAN30: 33 wt% urea, 42 wt% ammonium nitrate, [2]
        urea_ratio = 33/60.06*2/(33/60.06*2 + 42/80.043)
        ammonium_nitrate_ratio = 1 - urea_ratio
        
        ammonia_to_urea = ammonia.imass['NH3']*urea_ratio
        ammonia_to_ammonium_nitrate = ammonia.imass['NH3']*ammonium_nitrate_ratio
        
        nitric_acid.imass['HNO3'] = ammonia_to_ammonium_nitrate/17.031*63.01
        
        required_CO2 = ammonia_to_urea/17.031/NH3_CO2_molar_ratio*44.009
        
        if required_CO2 > carbon_dioxide.imass['CO2']:
            additional_carbon_dioxide.imass['CO2'] = required_CO2 - carbon_dioxide.imass['CO2']
        else:
            additional_carbon_dioxide.imass['CO2'] = 0
        
        urea_amount = self.urea_amount = required_CO2/44.009*self.efficiency/\
                                         (1 - (1 - self.efficiency)*(1-self.loss))*60.06
        vapor.imass['H2O'] = required_CO2/44.009*self.efficiency/\
                             (1 - (1 - self.efficiency)*(1-self.loss))*18.01528
        waste.imass['NH3'] = required_CO2/44.009*(self.ratio - 2*self.efficiency)*\
                             self.loss/(1 - (1 - self.efficiency)*(1-self.loss))*17.031
        waste.imass['CO2'] = required_CO2/44.009*(1 - self.efficiency)*\
                             self.loss/(1 - (1 - self.efficiency)*(1-self.loss))*44.009
        waste.imass['CO2'] += max(0, carbon_dioxide.imass['CO2'] - required_CO2)
        
        ammonium_nitrate_amount = ammonia_to_ammonium_nitrate/17.031*80.043
        
        UAN_solution.imass['UAN'] = urea_amount + ammonium_nitrate_amount
        
        N_amount = urea_amount/60.06*2*14.0067 + ammonium_nitrate_amount/80.043*2*14.0067
        
        UAN_total_amount = N_amount/0.3
        
        UAN_solution.imass['H2O'] = UAN_total_amount - UAN_solution.imass['UAN']
        
        water.imass['H2O'] = UAN_solution.imass['H2O']
        
        vapor.phase = 'g'
        waste.phase = 'g'
        
        # convert 0.18 MWh/tonne-urea and 0.95 MWh/tonne-urea to kW, [4]
        self.power_utility.consumption = (0.18 + 0.95)*1000/1000*urea_amount
    
    def _design(self):
        
        # reactor
        self.N = 1
        self.V = self.outs[0].F_vol
        self.P = self.ins[0].P
        Reactor._design(self)
        
        Design = self.design_results
        Design['Urea production capacity'] = self.urea_amount

# =============================================================================
# UreaSynthesis
# =============================================================================

# EURO to USD: [1]
@cost(basis='Production capacity', ID='Urea synthesizer', units='kg/h',
      cost=28/0.951*1000000, S=500*1000/24,
      CE=CEPCI_by_year[2016], n=0.67, BM=1.47)
class UreaSynthesis(SanUnit):
    '''
    Urea Synthesis.
    Cost of urea systhesis from [2].
    
    Parameters
    ----------
    ins : iterable
        ammonia, carbon_dioxide, additional_carbon_dioxide.
    outs : iterable
        urea, vapor, waste.
    ratio : float
        The overall ratio between NH3 and CO2 as reactants (after considering
        recycling of unconverted reactants).
    efficiency : float
        The overall conversion efficiency (after considering recycling of
        unconverted reactants) of CO2 to urea.
    loss : float
        The loss ratio of unconverted reactants before recycling.
    
    References
    ----------
    .. [1] https://www.irs.gov/individuals/international-taxpayers/yearly-\
        average-currency-exchange-rates (accessed 2025-02-09)
    .. [2] Devkota, S.; Karmacharya, P.; Maharjan, S.; Khatiwada, D.;
        Uprety, B. Decarbonizing Urea: Techno-Economic and Environmental
        Analysis of a Model Hydroelectricity and Carbon Capture Based
        Green Urea Production. Applied Energy 2024, 372, 123789.
        https://doi.org/10.1016/j.apenergy.2024.123789.
    .. [3] Palys, M. J.; Daoutidis, P. Techno-Economic Optimization of
        Renewable Urea Production for Sustainable Agriculture and CO2 
        tilization. J. Phys. Energy 2023, 6 (1), 015013.
        https://doi.org/10.1088/2515-7655/ad0ee6.
    '''
    _N_ins = 3
    _N_outs = 3
    _units= {'Production capacity':'kg/h'}

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 ratio=3.5, efficiency=0.8, loss=0.02):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.ratio = ratio
        self.efficiency = efficiency
        self.loss = loss

    def _run(self):
        
        ammonia, carbon_dioxide, additional_carbon_dioxide = self.ins
        urea, vapor, waste = self.outs
        
        NH3_CO2_molar_ratio = (self.ratio - (self.ratio - 2*self.efficiency)*(1-self.loss))/\
                              (1 - (1 - self.efficiency)*(1-self.loss))
        
        required_CO2 = ammonia.imass['NH3']/17.031/NH3_CO2_molar_ratio*44.009
        
        if required_CO2 > carbon_dioxide.imass['CO2']:
            additional_carbon_dioxide.imass['CO2'] = required_CO2 - carbon_dioxide.imass['CO2']
        else:
            additional_carbon_dioxide.imass['CO2'] = 0
        
        urea.imass['Urea'] = required_CO2/44.009*self.efficiency/\
                             (1 - (1 - self.efficiency)*(1-self.loss))*60.06
        vapor.imass['H2O'] = required_CO2/44.009*self.efficiency/\
                             (1 - (1 - self.efficiency)*(1-self.loss))*18.01528
        waste.imass['NH3'] = required_CO2/44.009*(self.ratio - 2*self.efficiency)*\
                             self.loss/(1 - (1 - self.efficiency)*(1-self.loss))*17.031
        waste.imass['CO2'] = required_CO2/44.009*(1 - self.efficiency)*\
                             self.loss/(1 - (1 - self.efficiency)*(1-self.loss))*44.009
        waste.imass['CO2'] += max(0, carbon_dioxide.imass['CO2'] - required_CO2)
        
        vapor.phase = 'g'
        waste.phase = 'g'
        
        # convert 0.18 MWh/tonne-urea and 0.95 MWh/tonne-urea to kW, [3]
        self.power_utility.consumption = (0.18 + 0.95)*1000/1000*urea.imass['urea']
    
    def _design(self):
        
        Design = self.design_results
        Design['Production capacity'] = self.outs[0].F_mass

# =============================================================================
# WWmixer
# =============================================================================

class WWmixer(SanUnit):
    '''
    A fake unit that mixes all wastewater streams and calculates C, N, P, and H2O
    amount.
    Parameters
    ----------
    ins : iterable
        supernatant, memdis_ww, ht_ww.
    outs : iterable(stream)
        mixture.
    '''
    _N_ins = 3
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        
    def _run(self):
        
        supernatant, memdis_ww, ht_ww = self.ins
        mixture = self.outs[0]
        
        mixture.mix_from(self.ins)
        
        HT = self.ins[2]._source.ins[0]._source.ins[0]._source.ins[0]._source.\
             ins[0]._source.ins[0]._source
        
        # only account for C and N from HT if they are not less than 0
        if HT.HTaqueous_C >= 0:
            mixture.imass['C'] += HT.HTaqueous_C
            mixture.imass['H2O'] -= HT.HTaqueous_C
        if HT.HTaqueous_N >=0:
            mixture.imass['N'] += HT.HTaqueous_N
            mixture.imass['H2O'] -= HT.HTaqueous_N

# =============================================================================
# WWTP
# =============================================================================

class WWTP(SanUnit):
    '''
    WWTP is a fake unit that can set up sludge biochemical compositions
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
    operation_hour : float
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
                 ww_2_dry_sludge=0.94, # [1]
                 sludge_moisture=0.99, sludge_dw_ash=0.257, 
                 sludge_afdw_lipid=0.204, sludge_afdw_protein=0.463, 
                 lipid_2_C=0.750,
                 protein_2_C=0.545,
                 carbo_2_C=0.400, 
                 lipid_2_H=0.125,
                 protein_2_H=0.068,
                 carbo_2_H=0.067, 
                 protein_2_N=0.159,
                 N_2_P=0.3927,
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
        self.lipid_2_C = lipid_2_C
        self.protein_2_C = protein_2_C
        self.carbo_2_C = carbo_2_C
        self.lipid_2_H = lipid_2_H
        self.protein_2_H = protein_2_H
        self.carbo_2_H = carbo_2_H
        self.protein_2_N = protein_2_N
        self.N_2_P = N_2_P
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
    
    @property
    def sludge_C_ratio(self):
       return self.sludge_dw_protein*self.protein_2_C + self.sludge_dw_lipid*self.lipid_2_C +\
           self.sludge_dw_carbo*self.carbo_2_C
    
    @property
    def sludge_H_ratio(self):
       return self.sludge_dw_protein*self.protein_2_H + self.sludge_dw_lipid*self.lipid_2_H +\
           self.sludge_dw_carbo*self.carbo_2_H
    
    @property
    def sludge_N_ratio(self):
       return self.sludge_dw_protein*self.protein_2_N
    
    @property
    def sludge_P_ratio(self):
       return self.sludge_N_ratio*self.N_2_P
    
    @property
    def sludge_O_ratio(self):
       return 1 - self.sludge_C_ratio - self.sludge_H_ratio -\
           self.sludge_N_ratio - self.sludge_dw_ash
    
    @property
    def sludge_C(self):
       return self.sludge_C_ratio*self.sludge_dw
    
    @property
    def sludge_H(self):
       return self.sludge_H_ratio*self.sludge_dw
    
    @property
    def sludge_N(self):
       return self.sludge_N_ratio*self.sludge_dw
    
    @property
    def sludge_P(self):
       return self.sludge_P_ratio*self.sludge_dw
    
    @property
    def sludge_O(self):
       return self.sludge_O_ratio*self.sludge_dw
    
    @property
    def AOSc(self):
       return (3*self.sludge_N_ratio/14.0067 + 2*self.sludge_O_ratio/15.999 -\
               self.sludge_H_ratio/1.00784)/(self.sludge_C_ratio/12.011)
    
    @property
    def sludge_HHV(self):
       return 100*(0.338*self.sludge_C_ratio + 1.428*(self.sludge_H_ratio -\
              self.sludge_O_ratio/8)) # [3]
    
    @property
    def H_C_eff(self):
        return (self.sludge_H/1.00784-2*self.sludge_O/15.999)/self.sludge_C*12.011