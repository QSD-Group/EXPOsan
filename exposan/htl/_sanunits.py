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
from qsdsan.unit_operations import Reactor, HXutility, Tank, IsothermalCompressor
from qsdsan import unit_operations as qsu
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
    'HydrothermalLiquefaction',
    'Hydrocracking',
    'Hydrotreating',
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
_m3_to_gal = auom('m3').conversion_factor('gallon')
_in_to_m = auom('inch').conversion_factor('m')
_m3perh_to_mmscfd = 1/1177.17 # H2

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
        residue, extracted.
    acid_vol : float
        0.5 M H2SO4 to hydrochar ratio, [mL/g].
    P_acid_recovery_ratio : float
        The ratio of phosphorus that can be extracted.
    hard_coal_HHV : float
        HHV of hard coal, MJ/kg.
    hydrochar_distance : float
        Distance between WRRFs and coal-based power plants, [km].
        
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
                 hard_coal_HHV=27.91, # ecoinvent
                 hydrochar_distance=100, P=None, tau=2, V_wf=0.8, # tau: [1]
                 length_to_diameter=2, N=1, V=10, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0, # use MixTank default value
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 304', # acid condition
                 vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.acid_vol = acid_vol
        self.P_acid_recovery_ratio = P_acid_recovery_ratio
        self.hard_coal_HHV = hard_coal_HHV
        self.hydrochar_distance = hydrochar_distance
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
        residue, extracted = self.outs
        
        self.HTL = self.ins[0]._source
        
        if hydrochar.F_mass <= 0:
            pass
        else:
            if self.HTL.hydrochar_P <= 0:
                residue.copy_like(hydrochar)
            else: 
                acid.imass['H2SO4'] = hydrochar.F_mass*self.acid_vol*0.5*98.079/1000
                # 0.5 M H2SO4 acid_vol (10 mL/1 g) hydrochar
                # 0.5 M H2SO4 density: 1.03 kg/L 
                # https://www.fishersci.se/shop/products/sulfuric-acid-0-5m-4/11943303 (accessed 2024-08-03)
                acid.imass['H2O'] = hydrochar.F_mass*self.acid_vol*1.03 -\
                                    acid.imass['H2SO4']
                
                residue.imass['Residue'] = hydrochar.F_mass - self.ins[0]._source.\
                                           hydrochar_P*self.P_acid_recovery_ratio
                
                extracted.copy_like(acid)
                extracted.imass['P'] = hydrochar.F_mass - residue.F_mass
                # assume just P can be extracted
                
                residue.phase = 's'
                
                residue.T = extracted.T = hydrochar.T
                residue.P = hydrochar.P
                # H2SO4 reacts with hydrochar to release heat and temperature will increase
            
    @property
    def residue_C(self):
        return self.ins[0]._source.hydrochar_C
    
    @property
    def residue_P(self):
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
    # TODO: add access date
    .. [2] https://world-nuclear.org/information-library/facts-and-figures/heat-values-of-various-fuels.
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

@cost(basis='Crystallizer volume', ID='Crystallizer',
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
        BioSTEAM for now.
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
# HydrothermalLiquefaction (sludge-specific)
# =============================================================================

class KnockOutDrum(Reactor):
    '''
    Knockout drum is an auxiliary unit for :class:`HydrothermalLiquefaction`.

    References
    ----------
    [1] Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels via
        Liquefaction - Hydrothermal Liquefaction Reactor Design: April 5, 2013;
        NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462, 1111191.
        https://doi.org/10.2172/1111191.

    See Also
    --------
    :class:`HydrothermalLiquefaction`
    '''
    _N_ins = 3
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 P=3049.7*6894.76, tau=0, V_wf=0,
                 length_to_diameter=2, diameter=None,
                 N=4, V=None,
                 auxiliary=True,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 drum_steel_cost_factor=1.5):
        # drum_steel_cost_factor: so the cost matches [1]
        # when do comparison, if fully consider scaling factor (2000 tons/day to 100 tons/day),
        # drum_steel_cost_factor should be around 3
        # but that is too high, we use 1.5 instead.

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.diameter = diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.drum_steel_cost_factor = drum_steel_cost_factor

    def _run(self):
        pass

    def _cost(self):
        Reactor._cost(self)

        purchase_costs = self.baseline_purchase_costs
        purchase_costs['Vertical pressure vessel'] *= self.drum_steel_cost_factor


@cost(basis='Treatment capacity', ID='Solids filter oil/water separator', units='lb/h',
      cost=3945523, S=1219765,
      CE=CEPCI_by_year[2011], n=0.68, BM=1.9)
class HydrothermalLiquefaction(Reactor):
    '''
    HTL converts dewatered sludge to biocrude, aqueous, off-gas, and hydrochar
    under elevated temperature (350°C) and pressure. The products
    percentage (wt%) is evaluated using the revised MCA model (Li et al.,
    2017; Leow et al., 2018) from known sludge composition (protein%, lipid%,
    and carbohydrate%, all afdw%).

    Notice that for HTL we just calculate each phases' total mass (except gas)
    and calculate C, N, and P amount in each phase as properties. We don't
    specify components for oil/char since we want to use MCA model to
    calculate C and N amount and it is not necessary to calculate every
    possible components since they will be treated in HT/AcidEx anyway. We
    also don't specify components for aqueous since we want to calculate
    aqueous C, N, and P based on mass balance closure. But later for CHG, HT, and HC, we specify
    each components (except aqueous phase) for the application of flash,
    distillation column, and CHP units.

    Parameters
    ----------
    ins : Iterable(stream)
        dewatered_sludge.
    outs : Iterable(stream)
        hydrochar, HTLaqueous, biocrude, offgas.
    lipid_2_biocrude : float
        Lipid to biocrude factor.
    protein_2_biocrude : float
        Protein to biocrude factor.
    carbo_2_biocrude : float
        Carbohydrate to biocrude factor.
    protein_2_gas : float
        Protein to gas factor.
    carbo_2_gas : float
        Carbohydrate to gas factor.
    biocrude_C_slope : float
        Biocrude carbon content slope.
    biocrude_C_intercept : float
        Biocrude carbon content intercept.
    biocrude_N_slope : float
        Biocrude nitrogen content slope.
    biocrude_H_slope : float
        Biocrude hydrogen content slope.
    biocrude_H_intercept : float
        Biocrude hydrogen content intercept.
    HTLaqueous_C_slope : float
        HTLaqueous carbon content slope.
    TOC_TC : float
        HTL TOC/TC.
    hydrochar_C_slope : float
        Hydrochar carbon content slope.
    hydrochar_H_slope : float
        Hydrochar hydrogen content slope.
    biocrude_moisture_content : float
        Biocrude moisture content.
    hydrochar_P_recovery_ratio : float
        Hydrochar phosphorus to total phosphorus ratio.
    gas_composition : dict
        HTL offgas compositions.
    hydrochar_pre : float
        Hydrochar pressure, [Pa].
    HTLaqueous_pre : float
        HTL aqueous phase pressure, [Pa].
    biocrude_pre : float
        Biocrude pressure, [Pa].
    offgas_pre : float
        Offgas pressure, [Pa].
    eff_T : float
        HTL effluent temperature, [K].
    CAPEX_factor : float
        Factor used to adjust CAPEX.
    HTL_steel_cost_factor : float
        Factor used to adjust the cost of stainless steel.
    mositure_adjustment_exist_in_the_system : bool
        If a moisture adjustment unit exists, set to true.

    References
    ----------
    [1] Leow, S.; Witter, J. R.; Vardon, D. R.; Sharma, B. K.;
        Guest, J. S.; Strathmann, T. J. Prediction of Microalgae Hydrothermal
        Liquefaction Products from Feedstock Biochemical Composition.
        Green Chem. 2015, 17 (6), 3584-3599. https://doi.org/10.1039/C5GC00574D.
    [2] Li, Y.; Leow, S.; Fedders, A. C.; Sharma, B. K.; Guest, J. S.;
        Strathmann, T. J. Quantitative Multiphase Model for Hydrothermal
        Liquefaction of Algal Biomass. Green Chem. 2017, 19 (4), 1163-1174.
        https://doi.org/10.1039/C6GC03294J.
    [3] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J.
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products.
        Environ. Sci. Technol. 2018, 52 (21), 12717-12727.
        https://doi.org/10.1021/acs.est.8b04035.
    [4] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.;
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.;
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C.
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    [5] Matayeva, A.; Rasmussen, S. R.; Biller, P. Distribution of Nutrients and
        Phosphorus Recovery in Hydrothermal Liquefaction of Waste Streams.
        BiomassBioenergy 2022, 156, 106323.
        https://doi.org/10.1016/j.biombioe.2021.106323.
    [6] Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels
        via Liquefaction - Hydrothermal Liquefaction Reactor Design:
        April 5, 2013; NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462,
        1111191. https://doi.org/10.2172/1111191.
    '''
    _N_ins = 1
    _N_outs = 4
    _units = {'Treatment capacity': 'lb/h',
              'Solid filter and separator weight': 'lb'}

    auxiliary_unit_names = ('heat_exchanger', 'kodrum')

    _F_BM_default = {**Reactor._F_BM_default,
                      'Heat exchanger': 3.17}

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 lipid_2_biocrude=0.846, # [1]
                 protein_2_biocrude=0.445, # [1]
                 carbo_2_biocrude=0.205, # [1]
                 protein_2_gas=0.074, # [1]
                 carbo_2_gas=0.418, # [1]
                 biocrude_C_slope=-8.37, # [2]
                 biocrude_C_intercept=68.55, # [2]
                 biocrude_N_slope=0.133, # [2]
                 biocrude_H_slope=-2.61, # [2]
                 biocrude_H_intercept=8.20, # [2]
                 HTLaqueous_C_slope=478, # [2]
                 TOC_TC=0.764, # [3]
                 hydrochar_C_slope=1.75, # [2]
                 hydrochar_H_slope=0.141, # [2]
                 biocrude_moisture_content=0.063, # [4]
                 hydrochar_P_recovery_ratio=0.86, # [5]
                 gas_composition={'CH4': 0.050, 'C2H6': 0.032,
                                  'CO2': 0.918}, # [4]
                 hydrochar_pre=3029.7*6894.76, # [4]
                 HTLaqueous_pre=30*6894.76, # [4]
                 biocrude_pre=30*6894.76, # [4]
                 offgas_pre=30*6894.76, # [4]
                 eff_T=60+273.15, # [4]
                 P=None, tau=15/60, V_wf=0.45,
                 length_to_diameter=None, diameter=6.875*_in_to_m,
                 N=4, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Horizontal',
                 CAPEX_factor=1,
                 # this is equivalent to F_M = 2.1*2.7, which is a reasonable value for high-temperature, high-pressure, sludge-fed HTL reactors
                 HTL_steel_cost_factor=2.7, # so the cost matches [6]
                 mositure_adjustment_exist_in_the_system=False):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
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
        self.hydrochar_H_slope = hydrochar_H_slope
        self.biocrude_moisture_content = biocrude_moisture_content
        self.hydrochar_P_recovery_ratio = hydrochar_P_recovery_ratio
        self.gas_composition = gas_composition
        self.hydrochar_pre = hydrochar_pre
        self.HTLaqueous_pre = HTLaqueous_pre
        self.biocrude_pre = biocrude_pre
        self.offgas_pre = offgas_pre
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=eff_T, rigorous=True)
        self.kodrum = KnockOutDrum(ID=f'.{ID}_KOdrum')
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.diameter = diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.CAPEX_factor = CAPEX_factor
        self.HTL_steel_cost_factor = HTL_steel_cost_factor
        self.mositure_adjustment_exist_in_the_system = mositure_adjustment_exist_in_the_system

    def _run(self):

        dewatered_sludge = self.ins[0]
        hydrochar, HTLaqueous, biocrude, offgas = self.outs

        if self.mositure_adjustment_exist_in_the_system == True:
            self.WWTP = self.ins[0]._source.ins[0]._source.ins[0].\
                             _source.ins[0]._source
        else:
            self.WWTP = self.ins[0]._source.ins[0]._source.ins[0]._source

        dewatered_sludge_afdw = dewatered_sludge.imass['Sludge_lipid'] +\
                                dewatered_sludge.imass['Sludge_protein'] +\
                                dewatered_sludge.imass['Sludge_carbo']
        # just use afdw in revised MCA model, other places use dw

        self.afdw_lipid_ratio = self.WWTP.sludge_afdw_lipid
        self.afdw_protein_ratio = self.WWTP.sludge_afdw_protein
        self.afdw_carbo_ratio = self.WWTP.sludge_afdw_carbo

        # the following calculations are based on revised MCA model
        hydrochar.imass['Hydrochar'] = 0.377*self.afdw_carbo_ratio*dewatered_sludge_afdw

        HTLaqueous.imass['HTLaqueous'] = (0.481*self.afdw_protein_ratio +\
                                          0.154*self.afdw_lipid_ratio)*\
                                          dewatered_sludge_afdw
        # HTLaqueous is TDS in aqueous phase
        # 0.377, 0.481, and 0.154 don't have uncertainties because they are calculated values

        gas_mass = (self.protein_2_gas*self.afdw_protein_ratio + self.carbo_2_gas*self.afdw_carbo_ratio)*\
                       dewatered_sludge_afdw

        for name, ratio in self.gas_composition.items():
            offgas.imass[name] = gas_mass*ratio

        biocrude.imass['Biocrude'] = (self.protein_2_biocrude*self.afdw_protein_ratio +\
                                      self.lipid_2_biocrude*self.afdw_lipid_ratio +\
                                      self.carbo_2_biocrude*self.afdw_carbo_ratio)*\
                                      dewatered_sludge_afdw
        biocrude.imass['H2O'] = biocrude.imass['Biocrude']/(1 -\
                                self.biocrude_moisture_content) -\
                                biocrude.imass['Biocrude']

        HTLaqueous.imass['H2O'] = dewatered_sludge.F_mass - hydrochar.F_mass -\
                                  biocrude.F_mass - gas_mass - HTLaqueous.imass['HTLaqueous']
        # assume ash (all soluble based on Jones) goes to water

        hydrochar.phase = 's'
        offgas.phase = 'g'
        HTLaqueous.phase = biocrude.phase = 'l'

        hydrochar.P = self.hydrochar_pre
        HTLaqueous.P = self.HTLaqueous_pre
        biocrude.P = self.biocrude_pre
        offgas.P = self.offgas_pre

        for stream in self.outs : stream.T = self.heat_exchanger.T

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
        return (self.WWTP.AOSc*self.biocrude_C_slope + self.biocrude_C_intercept)/100 # [2]

    @property
    def biocrude_H_ratio(self):
        return (self.WWTP.AOSc*self.biocrude_H_slope + self.biocrude_H_intercept)/100 # [2]

    @property
    def biocrude_N_ratio(self):
        return self.biocrude_N_slope*self.WWTP.sludge_dw_protein # [2]

    @property
    def biocrude_C(self):
        return min(self.outs[2].F_mass*self.biocrude_C_ratio, self.WWTP.sludge_C)

    @property
    def HTLaqueous_C(self):
        return min(self.outs[1].F_vol*1000*self.HTLaqueous_C_slope*\
                   self.WWTP.sludge_dw_protein*100/1000000/self.TOC_TC,
                   self.WWTP.sludge_C - self.biocrude_C)

    @property
    def biocrude_H(self):
        return self.outs[2].F_mass*self.biocrude_H_ratio

    @property
    def biocrude_N(self):
        return min(self.outs[2].F_mass*self.biocrude_N_ratio, self.WWTP.sludge_N)

    # MJ/kg
    @property
    def biocrude_HHV(self):
        return 30.74 - 8.52*self.WWTP.AOSc +\
               0.024*self.WWTP.sludge_dw_protein # [2]

    @property
    def energy_recovery(self):
        return self.biocrude_HHV*self.outs[2].imass['Biocrude']/\
               (self.WWTP.outs[0].F_mass -\
               self.WWTP.outs[0].imass['H2O'])/self.WWTP.sludge_HHV # [2]

    @property
    def offgas_C(self):
        carbon = sum(self.outs[3].imass[self.gas_composition]*
                     [cmp.i_C for cmp in self.components[self.gas_composition]])
        return min(carbon, self.WWTP.sludge_C - self.biocrude_C - self.HTLaqueous_C)

    @property
    def hydrochar_C_ratio(self):
        return min(self.hydrochar_C_slope*self.WWTP.sludge_dw_carbo, 0.65) # [2]

    @property
    def hydrochar_C(self):
        return min(self.outs[0].F_mass*self.hydrochar_C_ratio, self.WWTP.sludge_C -\
                   self.biocrude_C - self.HTLaqueous_C - self.offgas_C)

    @property
    def hydrochar_P(self):
        return min(self.WWTP.sludge_P*self.hydrochar_P_recovery_ratio, self.outs[0].F_mass)

    @property
    def hydrochar_P_ratio(self):
        return min(self.WWTP.sludge_P*self.hydrochar_P_recovery_ratio, self.outs[0].F_mass)/self.outs[0].F_mass

    @property
    def hydrochar_H_ratio(self):
        return self.hydrochar_H_slope*self.WWTP.sludge_dw_carbo

    @property
    def hydrochar_O_ratio(self):
        return 1 - self.hydrochar_C_ratio - self.hydrochar_P_ratio - self.hydrochar_H_ratio

    # assume no N and no ash in hydrochar
    # MJ/kg
    @property
    def hydrochar_HHV(self):
        return (0.338*self.hydrochar_C_ratio + 1.428*(self.hydrochar_H_ratio - self.hydrochar_O_ratio/8))*100

    @property
    def HTLaqueous_N(self):
        return self.WWTP.sludge_N - self.biocrude_N

    @property
    def HTLaqueous_P(self):
        return self.WWTP.sludge_P*(1 - self.hydrochar_P_recovery_ratio)

    def _design(self):

        Design = self.design_results
        Design['Treatment capacity'] = self.ins[0].F_mass/_lb_to_kg

        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from((self.outs[1], self.outs[2], self.outs[3]))
        hx_outs0.copy_like(hx_ins0)
        hx_ins0.T = self.ins[0].T # temperature before/after HTL are similar
        hx_outs0.T = hx.T
        hx_ins0.P = hx_outs0.P = self.outs[0].P # cooling before depressurized, heating after pressurized
        # in other words, both heating and cooling are performed under relatively high pressure
        # hx_ins0.vle(T=hx_ins0.T, P=hx_ins0.P)
        # hx_outs0.vle(T=hx_outs0.T, P=hx_outs0.P)
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)

        self.P = self.ins[0].P
        Reactor._design(self)
        Design['Solid filter and separator weight'] = 0.2*Design['Weight']*Design['Number of reactors'] # assume stainless steel
        # based on [6], case D design table, the purchase price of solid filter and separator to
        # the purchase price of HTL reactor is around 0.2, therefore, assume the weight of solid filter
        # and separator is 0.2*single HTL weight*number of HTL reactors
        if self.include_construction:
            self.construction[0].quantity += Design['Solid filter and separator weight']*_lb_to_kg

        self.kodrum.V = self.F_mass_out/_lb_to_kg/1225236*4230/_m3_to_gal
        # in [6], when knockout drum influent is 1225236 lb/hr, single knockout
        # drum volume is 4230 gal

        self.kodrum.simulate()

    def _cost(self):
        Reactor._cost(self)
        self._decorated_cost()

        purchase_costs = self.baseline_purchase_costs
        for item in purchase_costs.keys():
            purchase_costs[item] *= self.CAPEX_factor

        purchase_costs['Horizontal pressure vessel'] *= self.HTL_steel_cost_factor

        for aux_unit in self.auxiliary_units:
            purchase_costs = aux_unit.baseline_purchase_costs
            installed_costs = aux_unit.installed_costs
            for item in purchase_costs.keys():
                purchase_costs[item] *= self.CAPEX_factor
                installed_costs[item] *= self.CAPEX_factor

# =============================================================================
# Hydrocracking (sludge-specific)
# =============================================================================

class Hydrocracking(Reactor):
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature (405 degC)
    and pressure to produce upgraded biooil. Co-product includes fuel gas.

    Parameters
    ----------
    ins : Iterable(stream)
        heavy_oil, hydrogen, catalyst_in.
    outs : Iterable(stream)
        hc_out, catalyst_out.
    WHSV : float
        Weight Hourly Space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime : float
        HC catalyst lifetime, [hr].
    catalyst_ID : str
        ID of the catalyst.
    hydrogen_P : float
        Hydrogen pressure, [Pa].
    hydrogen_rxned_to_inf_oil : float
        Reacted H2 to influent oil mass ratio.
    hydrogen_excess : float
        Actual hydrogen amount = hydrogen_rxned_to_biocrude*hydrogen_excess
    oil_yield : float
        Mass ratio of cracked oil to the sum of heavy oil and reacted H2,
        gas yield is calculated as 1-oil_yield (about 100% conversion as in [1]).
    HCin_T : float
        HC influent temperature, [K].
    HCrxn_T : float
        HC effluent (after reaction) temperature, [K].
    gas_composition : dict
        Composition of the gas products, will be normalized to 100% sum.
    oil_composition : dict
        Composition of the cracked oil, will be normalized to 100% sum.

    References
    ----------
    [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.;
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.;
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C.
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    '''
    _N_ins = 3
    _N_outs = 2

    auxiliary_unit_names = ('compressor', 'heat_exchanger')

    _F_BM_default = {**Reactor._F_BM_default,
                      'Heat exchanger': 3.17,
                      'Compressor': 1.1}

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 include_construction=False,
                 WHSV=0.625, # wt./hr per wt. catalyst [1]
                 catalyst_lifetime=5*7920, # 5 years [1]
                 catalyst_ID='HC_catalyst',
                 hydrogen_P=1039.7*6894.76,
                 hydrogen_rxned_to_inf_oil=0.01125,
                 hydrogen_excess=5.556,
                 oil_yield=1-0.03880-0.00630,
                 HCin_T=394+273.15,
                 HCrxn_T=451+273.15,
                 gas_composition={'CO2': 0.03880, 'CH4': 0.00630,},
                 oil_composition={
                    'CYCHEX': 0.03714, 'HEXANE': 0.01111,
                    'HEPTANE': 0.11474, 'OCTANE': 0.08125,
                    'C9H20': 0.09086, 'C10H22': 0.11756,
                    'C11H24': 0.16846, 'C12H26': 0.13198,
                    'C13H28': 0.09302, 'C14H30': 0.04643,
                    'C15H32': 0.03250, 'C16H34': 0.01923,
                    'C17H36': 0.00431, 'C18H38': 0.00099,
                    'C19H40': 0.00497, 'C20H42': 0.00033,
                    },
                 #combine C20H42 and PHYTANE as C20H42
                 # will not be a variable in uncertainty/sensitivity analysis
                 P=None, tau=5, void_fraciton=0.4, # Towler
                 length_to_diameter=2, diameter=None,
                 N=None, V=None, auxiliary=False, mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1.5,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, include_construction=include_construction)
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.catalyst_ID = catalyst_ID
        self.hydrogen_P = hydrogen_P
        self.hydrogen_rxned_to_inf_oil = hydrogen_rxned_to_inf_oil
        self.hydrogen_excess = hydrogen_excess
        self.oil_yield = oil_yield
        self.HCin_T = HCin_T
        self._mixed_in = Stream(f'{ID}_mixed_in')
        self.HCrxn_T = HCrxn_T
        self.gas_composition = gas_composition
        self.oil_composition = oil_composition
        IC_in = Stream(f'{ID}_IC_in')
        IC_out = Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in,
                                               outs=IC_out, P=None)
        hx_H2_in = Stream(f'{ID}_hx_H2_in')
        hx_H2_out = Stream(f'{ID}_hx_H2_out')
        self.heat_exchanger_H2 = HXutility(ID=f'.{ID}_hx_H2', ins=hx_H2_in, outs=hx_H2_out)
        hx_oil_in = Stream(f'{ID}_hx_oil_in')
        hx_oil_out = Stream(f'{ID}_hx_oil_out')
        self.heat_exchanger_oil = HXutility(ID=f'.{ID}_hx_oil', ins=hx_oil_in, outs=hx_oil_out)
        self.P = P
        self.tau = tau
        self.void_fraciton = void_fraciton
        self.length_to_diameter = length_to_diameter
        self.diameter = diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

    def _run(self):

        heavy_oil, hydrogen, catalyst_in = self.ins
        hc_out, catalyst_out = self.outs

        catalyst_in.imass[self.catalyst_ID] = heavy_oil.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts

        hydrogen_rxned_to_inf_oil = self.hydrogen_rxned_to_inf_oil
        hydrogen.imass['H2'] = heavy_oil.F_mass*hydrogen_rxned_to_inf_oil*self.hydrogen_excess
        hydrogen.phase = 'g'

        hydrocarbon_mass = heavy_oil.F_mass*(1 + hydrogen_rxned_to_inf_oil)
        # 100 wt% of heavy oil and reacted H2
        # nearly all input heavy oils and H2 will be converted to products [1]
        # spreadsheet HC calculation
        hc_out.phase = 'g'

        for name, ratio in self.HC_composition.items():
            hc_out.imass[name] = hydrocarbon_mass*ratio

        hc_out.imass['H2'] = heavy_oil.F_mass*hydrogen_rxned_to_inf_oil*(self.hydrogen_excess - 1)

        hc_out.P = heavy_oil.P
        hc_out.T = self.HCrxn_T

        hc_out.vle(T=hc_out.T, P=hc_out.P)

        cmps = self.components
        C_in = 0
        total_num = len(list(cmps))
        for num in range(total_num):
            C_in += heavy_oil.imass[str(list(cmps)[num])]*list(cmps)[num].i_C

        C_out = self.hydrocarbon_C

        if C_out < 0.95*C_in or C_out > 1.05*C_out :
            raise Exception('carbon mass balance is out of +/- 5% for HC')
        # make sure that carbon mass balance is within +/- 5%. Otherwise, an
        # exception will be raised.

    def _normalize_composition(self, dct):
        total = sum(dct.values())
        if total <=0: raise ValueError(f'Sum of total yields/composition should be positive, not {total}.')
        return {k:v/total for k, v in dct.items()}

    @property
    def gas_composition(self):
        return self._gas_composition
    @gas_composition.setter
    def gas_composition(self, comp_dct):
        self._gas_composition = self._normalize_composition(comp_dct)

    @property
    def oil_composition(self):
        return self._oil_composition
    @oil_composition.setter
    def oil_composition(self, comp_dct):
        self._oil_composition = self._normalize_composition(comp_dct)

    @property
    def HC_composition(self):
        '''Composition of gas and oil products, normalized to 100%.'''
        gas_composition = self.gas_composition
        oil_composition = self.oil_composition
        oil_yield = self.oil_yield
        gas_yield = 1 - oil_yield
        HC_composition = {k:v*gas_yield for k, v in gas_composition.items()}
        HC_composition.update({k:v*oil_yield for k, v in oil_composition.items()})
        return self._normalize_composition(HC_composition)

    @property
    def hydrocarbon_C(self):
        return sum(self.outs[0].imass[self.HC_composition]*
                   [cmp.i_C for cmp in self.components[self.HC_composition]])

    def _design(self):
        IC = self.compressor
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]
        IC_ins0.copy_like(self.ins[1])
        IC_outs0.copy_like(self.ins[1])
        IC_outs0.P = IC.P = self.hydrogen_P
        IC_ins0.phase = IC_outs0.phase = 'g'
        IC.simulate()

        hx_H2 = self.heat_exchanger_H2
        hx_H2_ins0, hx_H2_outs0 = hx_H2.ins[0], hx_H2.outs[0]
        hx_H2_ins0.copy_like(self.ins[1])
        hx_H2_outs0.copy_like(hx_H2_ins0)
        hx_H2_ins0.phase = hx_H2_outs0.phase = 'g'
        self._mixed_in.mix_from(self.ins)
        if not self.HCin_T: self.HCin_T = self._mixed_in.T
        hx_H2_outs0.T = self.HCin_T
        hx_H2_ins0.P = hx_H2_outs0.P = IC_outs0.P
        hx_H2.simulate_as_auxiliary_exchanger(ins=hx_H2.ins, outs=hx_H2.outs)

        hx_oil = self.heat_exchanger_oil
        hx_oil_ins0, hx_oil_outs0 = hx_oil.ins[0], hx_oil.outs[0]
        hx_oil_ins0.copy_like(self.ins[0])
        hx_oil_outs0.copy_like(hx_oil_ins0)
        hx_oil_outs0.T = self.HCin_T
        hx_oil_ins0.P = hx_oil_outs0.P = self.ins[0].P
        hx_oil.simulate_as_auxiliary_exchanger(ins=hx_oil.ins, outs=hx_oil.outs)

        self.P = min(IC_outs0.P, self.ins[0].P)

        V_H2 = self.ins[1].F_vol/self.hydrogen_excess*101325/self.hydrogen_P
        # just account for reacted H2
        V_biocrude = self.ins[0].F_vol
        self.V_wf = self.void_fraciton*V_biocrude/(V_biocrude + V_H2)
        Reactor._design(self)


# =============================================================================
# Hydrotreating (sludge-specific)
# =============================================================================

@cost(basis='Hydrogen_PSA', ID='PSA', units='mmscfd',
      cost=1750000, S=10,
      CE=CEPCI_by_year[2004], n=0.8, BM=2.47)
class Hydrotreating(Reactor):
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature (405 degC)
    and pressure to produce upgraded biooil. Co-product includes fuel gas.
    A pressure swing adsorption (PSA) process can be optionally included
    for H2 recovery.

    Parameters
    ----------
    ins : Iterable(stream)
        biocrude, hydrogen, catalyst_in.
    outs : Iterable(stream)
        ht_out, catalyst_out.
    WHSV : float
        Weight Hourly Space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime : float
        HT catalyst lifetime, [hr].
    catalyst_ID : str
        ID of the catalyst.
    hydrogen_P : float
        Hydrogen pressure, [Pa].
    hydrogen_rxned_to_inf_oil : float
        Reacted H2 to influent oil mass ratio.
    hydrogen_excess : float
        Actual hydrogen amount = hydrogen_rxned_to_biocrude*hydrogen_excess
    hydrocarbon_ratio : float
        Mass ratio of produced hydrocarbon to the sum of biocrude and reacted H2.
    HTin_T : float
        HT influent temperature, [K].
    HTrxn_T : float
        HT effluent (after reaction) temperature, [K].
    HT_composition : dict
        HT effluent composition.
    CAPEX_factor : float
        Factor used to adjust CAPEX.
    include_PSA : bool
        Whether to include pressure swing adsorption for H2 recovery.

    References
    ----------
    [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.;
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.;
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C.
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    [2] Towler, G.; Sinnott, R. Chapter 14 - Design of Pressure Vessels.
        In Chemical Engineering Design (Second Edition); Towler, G., Sinnott, R.,
        Eds.; Butterworth-Heinemann: Boston, 2013; pp 563-629.
        https://doi.org/10.1016/B978-0-08-096659-5.00014-6.
    '''
    _N_ins = 3
    _N_outs = 2
    auxiliary_unit_names = ('compressor', 'heat_exchanger')

    _F_BM_default = {**Reactor._F_BM_default,
                      'Heat exchanger': 3.17,
                      'Compressor': 1.1}

    _units = {'Hydrogen': 'mmscfd',
              'Hydrogen_PSA': 'mmscfd'}

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 WHSV=0.625, # wt./hr per wt. catalyst [1]
                 catalyst_lifetime=2*7920, # 2 years [1]
                 catalyst_ID='HT_catalyst',
                 hydrogen_P=1530*6894.76,
                 hydrogen_rxned_to_inf_oil=0.046,
                 hydrogen_excess=3,
                 # spreadsheet HT calculation
                 hydrocarbon_ratio=0.875, # 87.5 wt% of biocrude and reacted H2 [1]
                 # Tin = 174 C (345 F) based on Jones PNNL report
                 # however, the reaction releases heat and increase the temperature of effluent to 402 C (755.5 F)
                 HTin_T=174+273.15,
                 HTrxn_T=402+273.15, # [1]
                 HT_composition={'CH4': 0.02280, 'C2H6': 0.02923,
                                 'C3H8': 0.01650, 'C4H10': 0.00870,
                                 'TWOMBUTAN': 0.00408, 'NPENTAN': 0.00678,
                                 'TWOMPENTA': 0.00408, 'HEXANE': 0.00401,
                                 'TWOMHEXAN': 0.00408, 'HEPTANE': 0.00401,
                                 'CC6METH': 0.01020, 'PIPERDIN': 0.00408,
                                 'TOLUENE': 0.01013, 'THREEMHEPTA': 0.01020,
                                 'OCTANE': 0.01013, 'ETHCYC6': 0.00408,
                                 'ETHYLBEN': 0.02040, 'OXYLENE': 0.01020,
                                 'C9H20': 0.00408, 'PROCYC6': 0.00408,
                                 'C3BENZ': 0.01020, 'FOURMONAN': 0,
                                 'C10H22': 0.00240, 'C4BENZ': 0.01223,
                                 # C10H22 was originally 0.00203, but it is not
                                 # good for distillation column, the excess amount
                                 # is substracted from HEXANE, HEPTANE, TOLUENE,
                                 # OCTANE, and C9H20, which were originally 0.00408,
                                 # 0.00408, 0.01020, 0.01020, and 0.00408
                                 'C11H24': 0.02040, 'C10H12': 0.02040,
                                 'C12H26': 0.02040, 'OTTFNA': 0.01020,
                                 'C6BENZ': 0.02040, 'OTTFSN': 0.02040,
                                 'C7BENZ': 0.02040, 'C8BENZ': 0.02040,
                                 'C10H16O4': 0.01837, 'C15H32': 0.06120,
                                 'C16H34': 0.18360, 'C17H36': 0.08160,
                                 'C18H38': 0.04080, 'C19H40': 0.04080,
                                 'C20H42': 0.10200, 'C21H44': 0.04080,
                                 'TRICOSANE': 0.04080, 'C24H38O4': 0.00817,
                                 'C26H42O4': 0.01020, 'C30H62': 0.00203}, # [1]
                 # spreadsheet HT calculation
                 # will not be a variable in uncertainty/sensitivity analysis
                 P=None, tau=0.5, void_fraciton=0.4, # [2]
                 length_to_diameter=2, diameter=None,
                 N=None, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 CAPEX_factor=1,
                 include_PSA=False,
                 PSA_pre=717.4*6894.76,
                 PSA_efficiency=0.9,):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.catalyst_ID = catalyst_ID
        self.hydrogen_P = hydrogen_P
        self.hydrogen_rxned_to_inf_oil = hydrogen_rxned_to_inf_oil
        self.hydrogen_excess = hydrogen_excess
        self.hydrocarbon_ratio = hydrocarbon_ratio
        self.HTin_T = HTin_T
        self.HTrxn_T = HTrxn_T
        self.HT_composition = HT_composition
        IC_in = Stream(f'{ID}_IC_in')
        IC_out = Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in,
                                               outs=IC_out, P=None)
        hx_H2_in = Stream(f'{ID}_hx_H2_in')
        hx_H2_out = Stream(f'{ID}_hx_H2_out')
        self.heat_exchanger_H2 = HXutility(ID=f'.{ID}_hx_H2', ins=hx_H2_in, outs=hx_H2_out)
        hx_oil_in = Stream(f'{ID}_hx_oil_in')
        hx_oil_out = Stream(f'{ID}_hx_oil_out')
        self.heat_exchanger_oil = HXutility(ID=f'.{ID}_hx_oil', ins=hx_oil_in, outs=hx_oil_out)
        self.P = P
        self.tau = tau
        self.void_fraciton = void_fraciton
        self.length_to_diameter = length_to_diameter
        self.diameter = diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.CAPEX_factor = CAPEX_factor
        self.include_PSA = include_PSA
        self.PSA_pre = PSA_pre
        self.PSA_efficiency = PSA_efficiency

    def _run(self):

        biocrude, hydrogen, catalyst_in = self.ins
        ht_out, catalyst_out = self.outs

        self.HTL = self.ins[0]._source.ins[0]._source

        HT_composition = self.HT_composition
        if self.HTL.biocrude_N == 0:
            remove = HT_composition['PIPERDIN']
            for chemical in self.HT_composition.keys():
                HT_composition[chemical] /= (1-remove)
            HT_composition['PIPERDIN'] = 0

        catalyst_in.imass[self.catalyst_ID] = biocrude.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts

        hydrogen_excess = self.hydrogen_excess
        H2_rxned =  biocrude.imass['Biocrude']*self.hydrogen_rxned_to_inf_oil
        recovered_frac =  (hydrogen_excess - 1)*self.PSA_efficiency*float(self.include_PSA)
        hydrogen.imass['H2'] = H2_rxned*(hydrogen_excess - recovered_frac)
        hydrogen.phase = 'g'

        hydrocarbon_mass = biocrude.imass['Biocrude']*\
                           (1 + self.hydrogen_rxned_to_inf_oil)*\
                           self.hydrocarbon_ratio

        ht_out.phase = 'g'

        for name, ratio in self.HT_composition.items():
            ht_out.imass[name] = hydrocarbon_mass*ratio

        ht_out.imass['H2'] = H2_rxned*(self.hydrogen_excess - recovered_frac - 1)

        ht_out.imass['H2O'] = biocrude.F_mass + hydrogen.F_mass -\
                              hydrocarbon_mass - ht_out.imass['H2']
        # use water to represent HT aqueous phase,
        # C and N can be calculated base on MB closure.

        ht_out.P = biocrude.P

        ht_out.T = self.HTrxn_T

        ht_out.vle(T=ht_out.T, P=ht_out.P)

        if self.HTaqueous_C < -0.1*self.HTL.WWTP.sludge_C:
            raise Exception('carbon mass balance is out of +/- 10% for the whole system')
        # allow +/- 10% out of mass balance
        # should be no C in the aqueous phase, the calculation here is just for MB

        if self.HTaqueous_N < -0.1*self.HTL.WWTP.sludge_N:
            warn('nitrogen mass balance is out of +/- 10% for the whole system')
        # allow +/- 10% out of mass balance

        # possibility exist that more carbon is in biooil and gas than in
        # biocrude because we use the biooil/gas compositions to calculate
        # carbon. In this case, the C in HT aqueous phase will be negative.
        # It's OK if the mass balance is within +/- 10% of total carbon in
        # sludge. Otherwise, an exception will be raised.

    @property
    def hydrocarbon_C(self):
        return sum(self.outs[0].imass[self.HT_composition]*
                   [cmp.i_C for cmp in self.components[self.HT_composition]])

    @property
    def hydrocarbon_N(self):
        return sum(self.outs[0].imass[self.HT_composition]*
                   [cmp.i_N for cmp in self.components[self.HT_composition]])

    @property
    def HTaqueous_C(self):
        return self.HTL.biocrude_C - self.hydrocarbon_C
    # should be no C in the aqueous phase, the calculation here is just for MB

    @property
    def HTaqueous_N(self):
        return self.HTL.biocrude_N - self.hydrocarbon_N

    @property
    def PSA_efficiency(self):
        return self._PSA_efficiency
    @PSA_efficiency.setter
    def PSA_efficiency(self, i):
        if i > 1: raise Exception('PSA efficiency cannot be larger than 1.')
        self._PSA_efficiency = i

    def _design(self):

        IC = self.compressor
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]
        IC_ins0.copy_like(self.ins[1])
        IC_outs0.copy_like(self.ins[1])
        IC_outs0.P = IC.P = self.hydrogen_P
        IC_ins0.phase = IC_outs0.phase = 'g'
        IC.simulate()

        hx_H2 = self.heat_exchanger_H2
        hx_H2_ins0, hx_H2_outs0 = hx_H2.ins[0], hx_H2.outs[0]
        hx_H2_ins0.copy_like(self.ins[1])
        hx_H2_outs0.copy_like(hx_H2_ins0)
        hx_H2_ins0.phase = hx_H2_outs0.phase = 'g'
        hx_H2_outs0.T = self.HTin_T
        hx_H2_ins0.P = hx_H2_outs0.P = IC_outs0.P
        hx_H2.simulate_as_auxiliary_exchanger(ins=hx_H2.ins, outs=hx_H2.outs)

        hx_oil = self.heat_exchanger_oil
        hx_oil_ins0, hx_oil_outs0 = hx_oil.ins[0], hx_oil.outs[0]
        hx_oil_ins0.copy_like(self.ins[0])
        hx_oil_outs0.copy_like(hx_oil_ins0)
        hx_oil_outs0.T = self.HTin_T
        hx_oil_ins0.P = hx_oil_outs0.P = self.ins[0].P
        hx_oil.simulate_as_auxiliary_exchanger(ins=hx_oil.ins, outs=hx_oil.outs)

        self.P = min(IC_outs0.P, self.ins[0].P)

        V_H2 = self.ins[1].F_vol/self.hydrogen_excess*101325/self.hydrogen_P
        # just account for reacted H2
        V_biocrude = self.ins[0].F_vol
        self.V_wf = self.void_fraciton*V_biocrude/(V_biocrude + V_H2)
        Reactor._design(self)

        Design = self.design_results
        factor = float(self.include_PSA)
        Design['Hydrogen_PSA'] = self.ins[1].F_vol*_m3perh_to_mmscfd*101325/self.PSA_pre*factor
        Design['PSA'] = 0.5*Design['Weight']*Design['Number of reactors']*factor # assume stainless steel
        # based on [1], page 54, the purchase price of PSA to the purchase price of
        # HT reactor (excluding vessels and columns) is around 0.5,
        # therefore, assume the weight of PSA is 0.5*single HT weight*number of HT reactors
        if self.include_construction:
            self.construction[0].quantity += Design['PSA']*_lb_to_kg*factor

    def _cost(self):
        Reactor._cost(self)

        purchase_costs = self.baseline_purchase_costs
        CAPEX_factor = self.CAPEX_factor
        if self.include_PSA: self._decorated_cost()
        else: purchase_costs['Hydrogen_PSA'] = 0

        for item in purchase_costs.keys():
            purchase_costs[item] *= CAPEX_factor

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
       return (0.338*self.sludge_C_ratio + 1.428*(self.sludge_H_ratio -\
              self.sludge_O_ratio/8))*100 # [3]
    
    @property
    def H_C_eff(self):
        return (self.sludge_H/1.00784-2*self.sludge_O/15.999)/self.sludge_C*12.011