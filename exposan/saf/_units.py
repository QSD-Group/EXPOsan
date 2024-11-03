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

from biosteam import Facility, ProcessWaterCenter as PWC
from biosteam.units.decorators import cost
from biosteam.units.design_tools import CEPCI_by_year
from qsdsan import SanUnit, Stream, WasteStream
from qsdsan.sanunits import Reactor, IsothermalCompressor, HXutility, HXprocess, MixTank

__all__ = (
    # To be moved to QSDsan
    'HydrothermalLiquefaction',
    'Hydroprocessing',
    'PressureSwingAdsorption',
    'Electrochemical',
    'HydrogenCenter',
    'ProcessWaterCenter',
    'BiocrudeSplitter',
    'Conditioning',
    'Transportation',
    )

_lb_to_kg = 0.453592
_m3_to_gal = 264.172
_barrel_to_m3 = 42/_m3_to_gal # 1 barrel is 42 gallon
_in_to_m = 0.0254
_psi_to_Pa = 6894.76
_m3perh_to_mmscfd = 1/1177.17 # H2


# %%

# =============================================================================
# Knock-out Drum
# =============================================================================

class KnockOutDrum(Reactor):
    '''
    Knockout drum is an auxiliary unit for :class:`HydrothermalLiquefaction`,
    when the cost is calculated using generic pressure vessel algorithms.
    
    Parameters
    ----------
    F_M : dict
        Material factors used to adjust cost (only used `use_decorated_cost` is False).
    
    See Also
    --------
    :class:`qsdsan.sanunits.HydrothermalLiquefaction`
    
    :class:`qsdsan.sanunits.Reactor`
    
    :class:`biosteam.units.design_tools.PressureVessel`
    
    References
    ----------
    [1] Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels via
        Liquefaction - Hydrothermal Liquefaction Reactor Design: April 5, 2013;
        NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462, 1111191.
        https://doi.org/10.2172/1111191.
    '''
    _N_ins = 3
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', include_construction=False,
                 P=3049.7*_psi_to_Pa, tau=0, V_wf=0,
                 length_to_diameter=2, diameter=None,
                 N=4, V=None,
                 auxiliary=True,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 F_M={
                     'Horizontal pressure vessel': 1.5,
                     'Vertical pressure vessel': 1.5,
                     },
                 ):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                         include_construction=include_construction)
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
        self.F_M = F_M
    
    def _run(self):
        pass
    
    def _cost(self):
        Reactor._cost(self)


# =============================================================================
# Hydrothermal Liquefaction
# =============================================================================

# Original [1] is 1339 dry-ash free ton per day (tpd), ash content is 13%,
# which is 58176 wet lb/hr, scaling basis in the equipment table in [1]
# (flow for S100) does not seem right.
@cost(basis='Wet mass flowrate', ID='HTL system', units='lb/hr',
      cost=37486757, S=574476,
      CE=CEPCI_by_year[2011], n=0.77, BM=2.1)
@cost(basis='Wet mass flowrate', ID='Solids filter oil/water separator', units='lb/hr',
      cost=3945523, S=574476,
      CE=CEPCI_by_year[2011], n=0.68, BM=1.9)
@cost(basis='Wet mass flowrate', ID='Hot oil system', units='lb/hr',
      cost=4670532, S=574476,
      CE=CEPCI_by_year[2011], n=0.6, BM=1.4)
class HydrothermalLiquefaction(Reactor):
    '''
    HTL converts feedstock to gas, aqueous, biocrude, (hydro)char
    under elevated temperature and pressure.
    
    Parameters
    ----------
    ins : Iterable(stream)
        Feedstock into HTL.
    outs : Iterable(stream)
        Gas, aqueous, biocrude, char.
    T : float
        Temperature of the HTL reaction, [K].
    P : float
        Pressure when the reaction is at temperature, [Pa].
    dw_yields : dict
        Dry weight percentage yields of the four products (gas, aqueous, biocrude, char),
        will be normalized to 100% sum.
        Keys must be 'gas', 'aqueous', 'biocrude', and 'char'.
    gas_composition : dict
        Composition of the gaseous products INCLUDING water, will be normalized to 100% sum.
    aqueous_composition : dict
        Composition of the aqueous products EXCLUDING water, will be normalized to 100% sum.
        Water not allocated to other products will all go to aqueous.
    biocrude_composition : dict
        Composition of the biocrude products INCLUDING water, will be normalized to 100% sum.
    char_composition : dict
        Composition of the char products INCLUDING water, will be normalized to 100% sum.
        "Ash" in the feedstock will be simulated based on the setting of `adjust_char_by_ash`.
    adjust_char_by_ash : bool
        If True, all ash in 
        but EXCLUDING ash (all ash will remain as ash),
    internal_heat_exchanging : bool
        If to use product to preheat feedstock.
    eff_T: float
        HTL effluent temperature [K],
        if provided, will use an additional HX to control effluent temperature.
    eff_P: float
        HTL effluent pressure [Pa].
    use_decorated_cost : bool
        If True, will use cost scaled based on [1], otherwise will use generic
        algorithms for ``Reactor`` (``PressureVessel``).
    F_M : dict
        Material factors used to adjust cost (only used `use_decorated_cost` is False).
        

    See Also
    --------
    :class:`qsdsan.sanunits.KnockOutDrum`
    
    :class:`qsdsan.sanunits.Reactor`
    
    :class:`biosteam.units.design_tools.PressureVessel`

        
    References
    ----------
    [1] Leow, S.; Witter, J. R.; Vardon, D. R.; Sharma, B. K.;
        Guest, J. S.; Strathmann, T. J. Prediction of Microalgae Hydrothermal
        Liquefaction Products from Feedstock Biochemical Composition.
        Green Chem. 2015, 17 (6), 3584–3599. https://doi.org/10.1039/C5GC00574D.
    
    [2] Li, Y.; Leow, S.; Fedders, A. C.; Sharma, B. K.; Guest, J. S.;
        Strathmann, T. J. Quantitative Multiphase Model for Hydrothermal
        Liquefaction of Algal Biomass. Green Chem. 2017, 19 (4), 1163–1174.
        https://doi.org/10.1039/C6GC03294J.
    
    [3] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J.
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products.
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727.
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
    _units= {
        'Wet mass flowrate': 'lb/hr',
        'Solid filter and separator weight': 'lb',
        }
    
    auxiliary_unit_names=('hx', 'inf_heating_hx', 'eff_cooling_hx','kodrum')

    _F_BM_default = {
        **Reactor._F_BM_default,
        'Heat exchanger': 3.17,
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', include_construction=False,
                 T=280+273.15,
                 P=101325,
                 dw_yields={
                     'gas': 0,
                     'aqueous': 0,
                     'biocrude': 0,
                     'char': 1,
                     },
                 gas_composition={'HTLgas': 1},
                 aqueous_composition={'HTLaqueous': 1},
                 biocrude_composition={'HTLbiocrude': 1},
                 char_composition={'HTLchar': 1},
                 internal_heat_exchanging=True,
                 eff_T=60+273.15, # 140.7°F
                 eff_P=30*_psi_to_Pa,
                 use_decorated_cost=True,
                 tau=15/60, V_wf=0.45,
                 length_to_diameter=None,
                 diameter=6.875*_in_to_m,
                 N=4, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Horizontal',
                 # Use material factors so that the calculated reactor cost matches [6]
                 F_M={
                    'Horizontal pressure vessel': 2.7,
                    'Vertical pressure vessel': 2.7,
                    }
                 ):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                         include_construction=include_construction)
        self.T = T
        self.P = P
        self.dw_yields = dw_yields
        self.gas_composition = gas_composition
        self.aqueous_composition = aqueous_composition
        self.biocrude_composition = biocrude_composition
        self.char_composition = char_composition
        self.internal_heat_exchanging = internal_heat_exchanging
        inf_pre_hx = Stream(f'{ID}_inf_pre_hx')
        eff_pre_hx = Stream(f'{ID}_eff_pre_hx')
        inf_after_hx = Stream(f'{ID}_inf_after_hx')
        eff_after_hx = Stream(f'{ID}_eff_after_hx')
        self.hx = HXprocess(ID=f'.{ID}_hx',
                            ins=(inf_pre_hx, eff_pre_hx),
                            outs=(inf_after_hx, eff_after_hx))
        inf_heating_hx_out = Stream(f'{ID}_inf_heating_hx_out')
        self.inf_heating_hx = HXutility(ID=f'.{ID}_inf_heating_hx', ins=inf_after_hx, outs=inf_heating_hx_out, T=T, rigorous=True)
        self._inf_at_temp = Stream(f'{ID}_inf_at_temp')
        self._eff_at_temp = Stream(f'{ID}_eff_at_temp')
        eff_cooling_hx_out = Stream(f'{ID}eff_cooling_hx_out')
        self.eff_T = eff_T
        self.eff_P = eff_P
        self.eff_cooling_hx = HXutility(ID=f'.{ID}_eff_cooling_hx', ins=eff_after_hx, outs=eff_cooling_hx_out, T=eff_T, rigorous=True)
        self.use_decorated_cost = use_decorated_cost
        self.kodrum = KnockOutDrum(ID=f'.{ID}_KOdrum', include_construction=include_construction)
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
        self.F_M = F_M
        

    def _run(self):
        feed = self.ins[0]
        gas, aq, crude, char = outs = self.outs        
        tot_dw = feed.F_mass - feed.imass['Water']
        comps = (
            self.gas_composition,
            self.aqueous_composition,
            self.biocrude_composition,
            self.char_composition,
            )
        for out, comp in zip(outs, comps):
            out.empty()
            for k, v in comp.items():
                out.imass[k] = v
        
        dw_yields = self.dw_yields
        gas.F_mass = tot_dw * dw_yields['gas']
        aq.F_mass = tot_dw * dw_yields['aqueous']
        crude.F_mass = tot_dw * dw_yields['biocrude']
        char.F_mass = tot_dw * dw_yields['char']

        aq.imass['Water'] = feed.imass['Water'] - sum(i.imass['Water'] for i in (gas, crude, char))
        
        for i in outs:
            i.T = self.T
            i.P = self.P
        
        self._eff_at_temp.mix_from(outs)
        
        gas.phase = 'g'
        char.phase = 's'
        aq.phase = crude.phase = 'l'
        
        for attr, val in zip(('T', 'P'), (self.eff_T, self.eff_P)):
            if val:
                for i in self.outs: setattr(i, attr, val)

    
    def _design(self):
        hx = self.hx
        inf_heating_hx = self.inf_heating_hx
        inf_hx_in, inf_hx_out = inf_heating_hx.ins[0], inf_heating_hx.outs[0]
        
        if self.internal_heat_exchanging:
            # Use HTL product to heat up influent
            inf_pre_hx, eff_pre_hx = hx.ins
            inf_pre_hx.copy_like(self.ins[0])
            eff_pre_hx.copy_like(self._eff_at_temp)
            hx.simulate()
    
            # Additional heating, if needed
            inf_hx_in.copy_like(hx.outs[0])
            inf_hx_out.copy_flow(inf_hx_in)
        else:
            hx.empty()
            # Influent heating to HTL conditions
            inf_hx_in.copy_like(self.ins[0])

        inf_hx_out.copy_flow(inf_hx_in)
        inf_hx_out.T = self.T
        inf_hx_out.P = self.P # this may lead to HXN error, when at pressure
        inf_heating_hx.simulate_as_auxiliary_exchanger(ins=inf_heating_hx.ins, outs=inf_heating_hx.outs)
            
        # Additional cooling, if needed
        eff_cooling_hx = self.eff_cooling_hx        
        if self.eff_T:
            eff_hx_in, eff_hx_out = eff_cooling_hx.ins[0], eff_cooling_hx.outs[0]
            eff_hx_in.copy_like(self._eff_at_temp)
            eff_hx_out.mix_from(self.outs)
            eff_cooling_hx.simulate_as_auxiliary_exchanger(ins=eff_cooling_hx.ins, outs=eff_cooling_hx.outs)
        else:
            eff_cooling_hx.empty()

        Reactor._design(self)
        
        Design = self.design_results
        Design['Solid filter and separator weight'] = 0.2*Design['Weight']*Design['Number of reactors'] # assume stainless steel
        # based on [6], case D design table, the purchase price of solid filter and separator to
        # the purchase price of HTL reactor is around 0.2, therefore, assume the weight of solid filter
        # and separator is 0.2*single HTL weight*number of HTL reactors
        if self.include_construction:
            self.construction[0].quantity += Design['Solid filter and separator weight']*_lb_to_kg
        
        kodrum = self.kodrum
        if self.use_decorated_cost is True:
            kodrum.empty()
        else:
            kodrum.V = self.F_mass_out/_lb_to_kg/1225236*4230/_m3_to_gal
            # in [6], when knockout drum influent is 1225236 lb/hr, single knockout
            # drum volume is 4230 gal
            
            kodrum.simulate()
        
    def _cost(self):
        Design = self.design_results
        Design.clear()
        self.baseline_purchase_costs.clear()
        if self.use_decorated_cost:
            ins0 = self.ins[0]
            Design['Wet mass flowrate'] = ins0.F_mass/_lb_to_kg         
            self._decorated_cost()
        else: Reactor._cost(self)
        
                
    def _normalize_composition(self, dct):
        total = sum(dct.values())
        if total <=0: raise ValueError(f'Sum of total yields/compositions should be positive, not {total}.')
        return {k:v/total for k, v in dct.items()}
    
    @property
    def dw_yields(self):
        return self._dw_yields
    @dw_yields.setter
    def dw_yields(self, comp_dct):
        self._dw_yields = self._normalize_composition(comp_dct)

    @property
    def gas_composition(self):
        return self._gas_composition
    @gas_composition.setter
    def gas_composition(self, comp_dct):
        self._gas_composition = self._normalize_composition(comp_dct)
        
    @property
    def aqueous_composition(self):
        return self._aqueous_composition
    @aqueous_composition.setter
    def aqueous_composition(self, comp_dct):
        self._aqueous_composition = self._normalize_composition(comp_dct)
        
    @property
    def biocrude_composition(self):
        return self._biocrude_composition
    @biocrude_composition.setter
    def biocrude_composition(self, comp_dct):
        self._biocrude_composition = self._normalize_composition(comp_dct)
        
    @property
    def char_composition(self):
        return self._char_composition
    @char_composition.setter
    def char_composition(self, comp_dct):
        self._char_composition = self._normalize_composition(comp_dct)
    
    @property
    def biocrude_HHV(self):
        '''Higher heating value of the biocrude, MJ/kg.'''
        crude = self.outs[2]
        return crude.HHV/crude.F_mass/1e3
    
    @property
    def energy_recovery(self):
        '''Energy recovery calculated as the HHV of the biocrude over the HHV of the feedstock.'''
        feed = self.ins[0]
        return self.biocrude_HHV/(feed.HHV/feed.F_mass/1e3)
    


# %%

# =============================================================================
# Hydroprocessing
# =============================================================================

@cost(basis='Oil lb flowrate', ID='Hydrocracker', units='lb/hr',
      cost=25e6, # installed cost
      S=5963, # S338 in [1]
      CE=CEPCI_by_year[2007], n=0.75, BM=1)
@cost(basis='Oil lb flowrate', ID='Hydrotreater', units='lb/hr',
      cost=27e6, # installed cost
      S=69637, # S135 in [1]
      CE=CEPCI_by_year[2007], n=0.68, BM=1)
@cost(basis='PSA H2 lb flowrate', ID='PSA', units='lb/hr', # changed scaling basis
      cost=1750000, S=5402, # S135 in [1]
      CE=CEPCI_by_year[2004], n=0.8, BM=2.47)
class Hydroprocessing(Reactor):
    '''
    For fuel upgrading processes such as hydrocracking and hydrotreating.
    Co-product includes fuel gas and aqueous stream.
    
    Note that addition units are needed to fractionate the product
    into the gas, aqueous, and oil streams.
    
    Parameters
    ----------
    ins : Iterable(stream)
        Influent crude oil, makeup H2, catalyst_in.
        Note that the amount of makeup H2 will be back-calculated upon simulation.
    outs : Iterable(stream)
        Mixed products (oil and excess hydrogen, fuel gas, as well as the aqueous stream), catalyst_out.
    T: float
        Operating temperature, [K].
    P : float
        Operating pressure, [Pa].
    WHSV: float
        Weight hourly space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime: float
        Catalyst lifetime, [hr].
    catalyst_ID : str
        ID of the catalyst.
    hydrogen_rxned_to_inf_oil: float
        Reacted H2 to influent oil mass ratio.
    hydrogen_ratio : float
        Total hydrogen amount = hydrogen_rxned * hydrogen_ratio,
        excess hydrogen will be included in the fuel gas.
    include_PSA : bool
        Whether to include a pressure swing adsorption unit to recover H2.
    PSA_efficiency : float
        H2 recovery efficiency of the PSA unit,
        will be set to 0 if `include_PSA` is False.
    gas_yield : float
        Mass ratio of fuel gas to the sum of influent oil and reacted H2.
    oil_yield : float
        Mass ratio of treated oil to the sum of influent oil and reacted H2.
    gas_composition: dict
        Composition of the gas products (excluding excess H2), will be normalized to 100% sum.
    oil_composition: dict
        Composition of the treated oil, will be normalized to 100% sum.
    aqueous_composition: dict
        Composition of the aqueous product, yield will be calculated as 1-gas-oil.
    internal_heat_exchanging : bool
        If to use effluent to preheat influent.
    use_decorated_cost : str
        Either 'Hydrotreater' or 'Hydrotreater' to use the corresponding
        decorated cost, otherwise, will use generic
        algorithms for ``Reactor`` (``PressureVessel``).

        
    See Also
    --------    
    :class:`qsdsan.sanunits.Reactor`
    
    :class:`biosteam.units.design_tools.PressureVessel`
        
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
    _units= {
        'Oil lb flowrate': 'lb/hr',
        'PSA H2 lb flowrate': 'lb/hr',
        }  
    _F_BM_default = {
        **Reactor._F_BM_default,
        'Heat exchanger': 3.17,
        'Compressor': 1.1,
        }
    auxiliary_unit_names=('compressor','hx', 'hx_inf_heating',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 include_construction=False,
                 T=451+273.15,
                 P=1039.7*_psi_to_Pa,
                 WHSV=0.625, # wt./hr per wt. catalyst [1]
                 catalyst_lifetime=5*7920, # 5 years [1]
                 catalyst_ID='HC_catalyst',
                 hydrogen_rxned_to_inf_oil=0.01125,
                 hydrogen_ratio=5.556,
                 include_PSA=False,
                 PSA_efficiency=0.9,
                 gas_yield=0.03880-0.00630,
                 oil_yield=1-0.03880-0.00630,
                 gas_composition={'CO2':0.03880, 'CH4':0.00630,},
                 oil_composition={
                    'CYCHEX':0.03714, 'HEXANE':0.01111,
                    'HEPTANE':0.11474, 'OCTANE':0.08125,
                    'C9H20':0.09086, 'C10H22':0.11756,
                    'C11H24':0.16846, 'C12H26':0.13198,
                    'C13H28':0.09302, 'C14H30':0.04643,
                    'C15H32':0.03250, 'C16H34':0.01923,
                    'C17H36':0.00431, 'C18H38':0.00099,
                    'C19H40':0.00497, 'C20H42':0.00033, # combine C20H42 and PHYTANE as C20H42
                    },
                 aqueous_composition={'Water':1},
                 internal_heat_exchanging=True,
                 use_decorated_cost=True,
                 tau=15/60, # set to the same as HTL as in [1]
                 V_wf=0.4, # void_fraciton=0.4, # Towler
                 length_to_diameter=2, diameter=None,
                 N=None, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1.5,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 ):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, include_construction=include_construction)
        self.T = T
        self.P = P
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.catalyst_ID = catalyst_ID
        self.hydrogen_rxned_to_inf_oil = hydrogen_rxned_to_inf_oil
        self.hydrogen_ratio = hydrogen_ratio
        self.include_PSA = include_PSA
        self.PSA_efficiency = PSA_efficiency
        self.gas_yield = gas_yield
        self.oil_yield = oil_yield
        self.gas_composition = gas_composition
        self.oil_composition = oil_composition
        self.aqueous_composition = aqueous_composition
        self.internal_heat_exchanging = internal_heat_exchanging
        # For H2 compressing
        IC_in = Stream(f'{ID}_IC_in')
        IC_out = Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in,
                                               outs=IC_out, P=P)
        # For influent heating      
        inf_pre_hx = Stream(f'{ID}_inf_pre_hx')
        eff_pre_hx = Stream(f'{ID}_eff_pre_hx')
        inf_after_hx = Stream(f'{ID}_inf_after_hx')
        eff_after_hx = Stream(f'{ID}_eff_after_hx')
        self.hx = HXprocess(ID=f'.{ID}_hx',
                            ins=(inf_pre_hx, eff_pre_hx),
                            outs=(inf_after_hx, eff_after_hx))
        inf_heating_hx_out = Stream(f'{ID}_inf_heating_hx_out')
        self.inf_heating_hx = HXutility(ID=f'.{ID}_inf_heating_hx', ins=inf_after_hx, outs=inf_heating_hx_out, T=T, rigorous=True)
        self.use_decorated_cost = use_decorated_cost
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
        
    def _run(self):
        inf_oil, makeup_hydrogen, catalyst_in = self.ins
        eff_oil, catalyst_out = self.outs
        
        catalyst_in.imass[self.catalyst_ID] = inf_oil.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        hydrogen_rxned_to_inf_oil = self.hydrogen_rxned_to_inf_oil
        hydrogen_ratio = self.hydrogen_ratio
        H2_rxned =  inf_oil.F_mass * hydrogen_rxned_to_inf_oil
        H2_tot = H2_rxned * hydrogen_ratio
        H2_residual = self._H2_residual = H2_tot - H2_rxned
        H2_recycled = self._H2_recycled = H2_residual * self.PSA_efficiency
        H2_wasted = H2_residual - H2_recycled

        eff_oil.copy_like(inf_oil)
        eff_oil.phase = inf_oil.phase
        eff_oil.empty()
        eff_oil.imass[self.eff_composition.keys()] = self.eff_composition.values()
        eff_oil.F_mass = inf_oil.F_mass*(1 + hydrogen_rxned_to_inf_oil)
        eff_oil.imass['H2'] = H2_wasted
        eff_oil.P = self.P
        eff_oil.T = self.T
        eff_oil.vle(T=eff_oil.T, P=eff_oil.P)
        
        makeup_hydrogen.imass['H2'] = H2_rxned + H2_wasted
        makeup_hydrogen.phase = 'g'

    def _design(self):
        Design = self.design_results
        Design.clear()
        Design['PSA H2 lb flowrate'] = self._H2_residual / _lb_to_kg
        Design['H2 recycled'] = recycled = self._H2_recycled / _lb_to_kg
        Design['Oil lb flowrate'] = self.ins[0].F_mass/_lb_to_kg
        
        IC = self.compressor # for H2 compressing
        H2 = self.ins[1]
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]
        IC_ins0.copy_like(H2)
        IC_ins0.F_mass += recycled # including the compressing needs for the recycled H2
        IC_outs0.copy_like(IC_ins0)
        IC_outs0.P = IC.P = self.P
        IC_ins0.phase = IC_outs0.phase = 'g'
        IC.simulate()
        
        hx = self.hx
        inf_heating_hx = self.inf_heating_hx
        inf_hx_in, inf_hx_out = inf_heating_hx.ins[0], inf_heating_hx.outs[0]
        
        if self.internal_heat_exchanging:
            # Use HTL product to heat up influent
            inf_pre_hx, eff_pre_hx = hx.ins
            inf_pre_hx.copy_like(self.ins[0])
            eff_pre_hx.copy_like(self.outs[0])
            hx.simulate()
    
            # Additional heating, if needed
            inf_hx_in.copy_like(hx.outs[0])
            inf_hx_out.copy_flow(inf_hx_in)
        else:
            hx.empty()
            # Influent heating to HTL conditions
            inf_hx_in.copy_like(self.ins[0])

        inf_hx_out.copy_flow(inf_hx_in)
        inf_hx_out.T = self.T
        inf_hx_out.P = self.P
        inf_heating_hx.simulate_as_auxiliary_exchanger(ins=inf_heating_hx.ins, outs=inf_heating_hx.outs)

        Reactor._design(self)

    def _cost(self):
        Cost = self.baseline_purchase_costs
        Cost.clear()
        
        use_decorated_cost = self.use_decorated_cost
        include_PSA = self.include_PSA
        self._decorated_cost()
        
        if use_decorated_cost == 'Hydrocracker':
            Cost.pop('Hydrotreater')
        elif use_decorated_cost == 'Hydrotreater':
            Cost.pop('Hydrocracker')
        else:
            Cost.pop('Hydrocracker')
            Cost.pop('Hydrotreater')
            Reactor._cost(self)
        
        if not include_PSA: Cost.pop('PSA')

        
    def _normalize_yields(self):
        gas = self._gas_yield
        oil = self._oil_yield
        gas_oil = gas + oil
        aq = 0
        if gas_oil > 1:
            gas /=  gas_oil
            oil /= gas_oil
        else:
            aq = 1 - gas_oil
        self._gas_yield = gas
        self._oil_yield = oil
        self._aq_yield = aq
        
    _normalize_composition = HydrothermalLiquefaction._normalize_composition

    
    @property
    def gas_yield(self):
        return self._gas_yield
    @gas_yield.setter
    def gas_yield(self, gas):
        self._gas_yield = gas
        if hasattr(self, '_oil_yield'):
            self._normalize_yields()

    @property
    def oil_yield(self):
        return self._oil_yield
    @oil_yield.setter
    def oil_yield(self, oil):
        self._oil_yield = oil
        if hasattr(self, '_gas_yield'):
            self._normalize_yields()
            
    @property
    def aq_yield(self):
        return self._aq_yield

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
    def aqueous_composition(self):
        return self._aqueous_composition
    @aqueous_composition.setter
    def aqueous_composition(self, comp_dct):
        self._aqueous_composition = self._normalize_composition(comp_dct)
        
    @property
    def eff_composition(self):
        '''Composition of products, normalized to 100% sum.'''
        gas_composition = self.gas_composition
        oil_composition = self.oil_composition
        aqueous_composition = self.aqueous_composition
        oil_yield = self.oil_yield
        gas_yield = self.gas_yield
        aq_yield = self.aq_yield
        eff_composition = {k:v*gas_yield for k, v in gas_composition.items()}
        eff_composition.update({k:v*oil_yield for k, v in oil_composition.items()})
        eff_composition.update({k:v*aq_yield for k, v in aqueous_composition.items()})
        return self._normalize_composition(eff_composition)

    @property
    def PSA_efficiency(self):
        '''
        [float] H2 recovery efficiency of the PSA unit,
        will be set to 0 if `include_PSA` is False.
        '''
        if self.include_PSA: return self._PSA_efficiency
        return 0
    @PSA_efficiency.setter
    def PSA_efficiency(self, i):
        if i > 1: raise ValueError('PSA_efficiency cannot be larger than 1.')
        self._PSA_efficiency  = i


# %%

# =============================================================================
# Pressure Swing Adsorption
# =============================================================================

@cost(basis='PSA H2 lb flowrate', ID='PSA', units='lb/hr', # changed scaling basis
      cost=1750000, S=5402, # S135 in [1]
      CE=CEPCI_by_year[2004], n=0.8, BM=2.47)
class PressureSwingAdsorption(SanUnit):
    '''
    A pressure swing adsorption (PSA) process can be optionally included
    for H2 recovery.
    
    Parameters
    ----------
    ins : Iterable(stream)
        Mixed gas streams for H2 recovery.
    outs : Iterable(stream)
        Hydrogen, other gases.
    efficiency : float
        H2 recovery efficiency.
    PSA_compressor_P : float
        Pressure to compressed the generated H2 to, if desired, [Pa].
        
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
    _N_outs = 2  
    _units= {'PSA H2 lb flowrate': 'lb/hr',}
    _F_BM_default = {'Compressor': 1.1,}
    auxiliary_unit_names=('compressor',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 efficiency=0.9,
                 PSA_compressor_P=101325,
                 ):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.efficiency = efficiency
        # For H2 compressing
        P = self.PSA_compressor_P = PSA_compressor_P
        IC_in = Stream(f'{ID}_IC_in')
        IC_out = Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in, outs=IC_out, P=P)
        
    def _run(self):
        H2, others = self.outs       
        others.mix_from(self.ins)
        H2.imass['H2'] = recovered = others.imass['H2'] * self.efficiency
        others.imass['H2'] -= recovered
        
    def _design(self):
        self.design_results['PSA H2 lb flowrate'] = self.F_mass_in/_lb_to_kg
        IC = self.compressor # for H2 compressing
        H2 = self.ins[0]
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]
        IC_ins0.copy_like(H2)
        IC_outs0.copy_like(IC_ins0)
        IC_outs0.P = IC.P = self.PSA_compressor_P
        IC_ins0.phase = IC_outs0.phase = 'g'
        IC.simulate()
        
    @property
    def efficiency (self):
        return self._efficiency 
    @efficiency.setter
    def efficiency(self, i):
        if i > 1: raise Exception('Efficiency cannot be larger than 1.')
        self._efficiency  = i


# %%

@cost(basis='PSA H2 lb flowrate', ID='PSA', units='lb/hr', # changed scaling basis
      cost=1750000, S=5402, # S135 in [1]
      CE=CEPCI_by_year[2004], n=0.8, BM=2.47)
class Electrochemical(SanUnit):
    '''
    An electrochemical unit alternatively operated in 
    electrochemical oxidation (EO) and electrodialysis (ED) modes.
    
    The original design was costed for a 50,000 kg H2/d system.
    
    The `replacement_surrogate` stream is used to represent the annual replacement cost,
    its price is set based on `annual_replacement_ratio`.
    
    A pressure swing adsorption unit is included to clean up the recycled H2, if desired.
    
    Parameters
    ----------
    ins : Iterable(stream)
        Influent water, replacement_surrogate.
    outs : Iterable(stream)
        Mixed gas, recycled H2, recovered N, recovered P, treated water.
    removal : float or dict
        Removal of non-water components either by a universal factor when given as a float,
        or as indicated by the dict.
    gas_yield : float
        Dry mass yield of the gas products.
    N_IDs : Iterable(str)
        IDs of the components for nitrogen recovery.
    P_IDs : Iterable(str)
        IDs of the components for phosphorus recovery.
    K_IDs : Iterable(str)
        IDs of the components for potassium recovery.
    N_recovery : float
        Recovery efficiency for nitrogen components (set by `N_IDs`).
    P_recovery : float
        Recovery efficiency for phosphorus components (set by `P_IDs`).
    K_recovery : float
        Recovery efficiency for potassium components (set by `K_IDs`).
    EO_current_density : float
        Currenty density when operating in the electrochemical oxidation, [A/m2].
    ED_current_density : float
        Currenty density when operating in the electrodialysis mode, [A/m2].
    EO_voltage : float
        Voltage when operating in the electrochemical oxidation mode, [V].
    ED_voltage : float
        Voltage when operating in the electrodialysis mode, [V].
    EO_online_time_ratio : float
        Ratio of time operated in the electrochemical oxidation model,
        ED_online_time_ratio is calculated as 1 - EO_online_time_ratio.
    chamber_thickness : float
        Thickness of the unit chamber, [m].
    electrode_cost : float
        Unit cost of the electrodes, [$/m2].
    anion_exchange_membrane_cost : float
        Unit cost of the anion exchange membrane, [$/m2].
    cation_exchange_membrane_cost : float
        Unit cost of the cation exchange membrane, [$/m2].
    electrolyte_load : float
        Load of the electrolyte per unit volume of the unit, [kg/m3].
    electrolyte_price : float
        Unit price of the electrolyte, [$/kg].
        Note that the electrolyte is calculated as a capital cost because
        theoretically it is not consumed during operation
        (replacement cost calculated through `annual_replacement_ratio`).
    annual_replacement_ratio : float
        Annual replacement cost as a ratio of the total purchase cost.
    include_PSA : bool
        Whether to include a pressure swing adsorption unit to recover H2.
    PSA_efficiency : float
        H2 recovery efficiency of the PSA unit,
        will be set to 0 if `include_PSA` is False.
    PSA_compressor_P : float
        Pressure to compressed the generated H2 to, if desired, [Pa].


    References
    ----------
    [1] Jiang et al., 2024.
    '''
    
    _N_ins = 2
    _N_outs = 6
    _units= {'PSA H2 lb flowrate': 'lb/hr',}
    _F_BM_default = {'Compressor': 1.1,}
    auxiliary_unit_names=('compressor',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=1,
                 removal=0.75,
                 gas_yield=0.056546425,
                 gas_composition={
                     'N2': 0.000795785,
                     'H2': 0.116180614,
                     'O2': 0.48430472,
                     'CO2': 0.343804756,
                     'CO': 0.054914124,
                     },
                 N_IDs=('N',),
                 P_IDs=('P',),
                 K_IDs=('K',),
                 N_recovery=0.8,
                 P_recovery=0.99,
                 K_recovery=0.8,
                 EO_current_density=1500, # A/m2
                 ED_current_density=100, # A/m2
                 EO_voltage=5, # V
                 ED_voltage=30, # V
                 EO_online_time_ratio=8/(8+1.5),
                 chamber_thickness=0.02, # m
                 electrode_cost=40000, # $/m2
                 anion_exchange_membrane_cost=170, # $/m2
                 cation_exchange_membrane_cost=190, # $/m2
                 electrolyte_load=3*136, # kg/m3, 3 M of KH2PO4 (MW=136 k/mole)
                 electrolyte_price=30, # $/kg
                 annual_replacement_ratio=0.02,
                 include_PSA=False,
                 PSA_efficiency=0.9,
                 PSA_compressor_P=101325,
                 ):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
        self.removal = removal
        self.gas_yield = gas_yield
        self.gas_composition = gas_composition
        self.N_recovery = N_recovery
        self.P_recovery = P_recovery
        self.K_recovery = K_recovery
        self.N_IDs = N_IDs
        self.P_IDs = P_IDs
        self.K_IDs = K_IDs
        self.EO_current_density = EO_current_density
        self.ED_current_density = ED_current_density
        self.EO_voltage = EO_voltage
        self.ED_voltage = ED_voltage
        self.EO_online_time_ratio = EO_online_time_ratio
        self.chamber_thickness = chamber_thickness
        self.electrode_cost = electrode_cost
        self.anion_exchange_membrane_cost = anion_exchange_membrane_cost
        self.cation_exchange_membrane_cost = cation_exchange_membrane_cost
        self.electrolyte_price = electrolyte_price # costing like a CAPEX due to minimal replacement requirements
        self.electrolyte_load = electrolyte_load
        self.annual_replacement_ratio = annual_replacement_ratio
        self.include_PSA = include_PSA
        self.PSA_efficiency = PSA_efficiency
        # For H2 compressing
        P = self.PSA_compressor_P = PSA_compressor_P
        IC_in = Stream(f'{ID}_IC_in')
        IC_out = Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in, outs=IC_out, P=P)
        
    
    def _run(self):
        gas, H2, N, P, K, eff = self.outs
        eff.copy_like(self.ins[0])
        water_in = eff.imass['Water']
        eff.imass['Water'] = 0
        dry_in = eff.F_mass
               
        fert_IDs = self.N_IDs, self.P_IDs, self.K_IDs
        recoveries = self.N_recovery, self.P_recovery, self.K_recovery
        for IDs, out, recovery in zip(fert_IDs, (N, P, K), recoveries):
            out.imass[IDs] = eff.imass[IDs] * recovery
            eff.imass[IDs] -= out.imass[IDs]
        
        gas.empty()
        comp = self.gas_composition
        gas.imass[list(comp.keys())] = list(comp.values())
        gas.F_mass = dry_in * self.gas_yield
        self._PSA_H2_lb_flowrate = gas.F_mass / _lb_to_kg
        gas.phase = 'g'
        
        removal = self.removal
        if type(removal) is float: eff.F_mass *= (1-removal)
        else:
            for k, v in removal.items():
                eff.imass[k] *= (1-v)

        eff.imass['Water'] = water_in
        
        H2_tot = gas.imass['H2']
        H2.imass['H2'] = H2_recycled = H2_tot * self.PSA_efficiency
        gas.imass['H2'] = H2_tot - H2_recycled
        
        
    def _design(self):
        Design = self.design_results
        
        # 96485 is the Faraday constant C/mol (A·s/mol e)
        # MW of H2 is 2 g/mol, 2 electrons per mole of H2
        factor = 2/(2/1e3) * 96485 # (A·s/kg H2)
        H2_production = self.outs[1].imass['H2'] / 3600 # kg/s
        current_eq = factor * H2_production # A
        area = current_eq / self.average_current_density
        Design['Area'] = area
        Design['Volume'] = area * self.chamber_thickness
        
        EO_power = self.EO_current_density * self.EO_voltage # W/m2, when online
        EO_electricity_per_area = EO_power/1e3 * self.EO_online_time_ratio # kWh/h/m2
        Design['EO electricity'] = EO_electricity = area * EO_electricity_per_area # kWh/h
        
        ED_power = self.ED_current_density * self.ED_voltage # W/m2, when online
        ED_electricity_per_area = ED_power/1e3 * self.ED_online_time_ratio # kWh/h/m2
        Design['ED electricity'] = ED_electricity = area * ED_electricity_per_area # kWh/h
        total_power = EO_electricity + ED_electricity
        self.power_utility.consumption = total_power
        #!!! unsure of this calculation
        self._FE = current_eq/(total_power*1e3) if total_power else 0
        
        Design['PSA H2 lb flowrate'] = self._PSA_H2_lb_flowrate
        IC = self.compressor # for H2 compressing
        H2 = self.outs[1]
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]
        IC_ins0.copy_like(H2)
        IC_outs0.copy_like(IC_ins0)
        IC_outs0.P = IC.P = self.PSA_compressor_P
        IC_ins0.phase = IC_outs0.phase = 'g'
        IC.simulate()
        

    def _cost(self):
        Design = self.design_results
        Cost = self.baseline_purchase_costs
        Cost.clear()
        
        stack_cost = self.electrode_cost+self.anion_exchange_membrane_cost+self.cation_exchange_membrane_cost
        Cost['Stack'] = stack_cost * Design['Area']
        Cost['Electrolyte'] = self.electrolyte_load*Design['Volume']*self.electrolyte_price
        cell_cost = sum(Cost.values())
        
        self._decorated_cost()
        
        replacement = self.ins[1]
        replacement.imass['Water'] = 1
        try: hours = self.system.operating_hours
        except: hours = 365*24
        replacement.price = cell_cost*self.annual_replacement_ratio/replacement.F_mass/hours
        
        
    def _normalize_composition(self, dct):
        total = sum(dct.values())
        if total <=0: raise ValueError(f'Sum of total compositions should be positive, not {total}.')
        return {k:v/total for k, v in dct.items()}
    
    @property
    def gas_composition(self):
        return self._gas_composition
    @gas_composition.setter
    def gas_composition(self, comp_dct):
        self._gas_composition = self._normalize_composition(comp_dct)
        
    @property
    def ED_online_time_ratio(self):
        '''Ratio of electrodialysis in operation.'''
        return 1 - self.EO_online_time_ratio
        
    @property
    def average_current_density(self):
        '''Currenty density of EO/ED averaged by online hours, [A/m2].'''
        return (self.EO_current_density*self.EO_online_time_ratio +
                self.ED_current_density*self.ED_online_time_ratio)
    
    @property
    def EO_electricity_ratio(self):
        '''Ratio of electricity used by electrochemical oxidation.'''
        EO = self.EO_current_density * self.EO_online_time_ratio
        ED = self.ED_current_density * self.ED_online_time_ratio
        return EO/(EO+ED)
    
    @property
    def ED_electricity_ratio(self):
        '''Ratio of electricity used by electrodialysis.'''
        return 1 - self.EO_electricity_ratio
    
    @property
    def FE(self):
        '''Faradaic efficiency of the combined EO and ED unit.'''
        return self._FE

    @property
    def PSA_efficiency(self):
        '''
        [float] H2 recovery efficiency of the PSA unit,
        will be set to 0 if `include_PSA` is False.
        '''
        if self.include_PSA: return self._PSA_efficiency
        return 0
    @PSA_efficiency.setter
    def PSA_efficiency(self, i):
        if i > 1: raise ValueError('PSA_efficiency cannot be larger than 1.')
        self._PSA_efficiency  = i


# %%

class HydrogenCenter(Facility):
    '''
    Calculate the amount of needed makeup hydrogen based on recycles and demands.
    
    This unit is for mass balance purpose only, not design/capital cost is included.
    
    ins and outs will be automatically created based on provided
    process and recycled H2 streams.
    
    ins: makeup H2, recycled H2.
    
    outs: process H2, excess H2.
    
    Notes
    -----
    When creating the unit, no ins and outs should be give (will be automatically created),
    rather, recycled and process H2 streams should be provided.
    
    
    Parameters
    ----------
    process_H2_streams : Iterable(stream)
        Process H2 streams (i.e., H2 demand) across the system.
    recycled_H2_streams : Iterable(stream)
        Recycled H2 streams across the system.
    makeup_H2_price : float
        Price of the makeup H2 (cost).
    excess_H2_price : float
        Price of the excess H2 (revenue).
        
    See Also
    --------
    :class:`biosteam.facilities.ProcessWaterCenter`
    '''
    
    ticket_name = 'H2C'
    network_priority = 2
    _N_ins = 2
    _N_outs = 2

    def __init__(self, ID='',
                 process_H2_streams=(), recycled_H2_streams=(),
                 makeup_H2_price=0, excess_H2_price=0):
        ins = (WasteStream('makeup_H2'), WasteStream('recycled_H2'))
        outs = (WasteStream('process_H2'), WasteStream('excess_H2'))
        Facility.__init__(self, ID, ins=ins, outs=outs)
        self.process_H2_streams = process_H2_streams
        self.recycled_H2_streams = recycled_H2_streams
        self.makeup_H2_price = makeup_H2_price
        self.excess_H2_price = excess_H2_price

    def _run(self):
        makeup, recycled = self.ins
        process, excess = self.outs
        
        for i in self.ins+self.outs:
            if i.F_mass != i.imass['H2']:
                raise RuntimeError(f'Streams in `{self.ID}` should only include H2, '
                                   f'the stream {i.ID} contains other components.')
        
        process_streams = self.process_H2_streams
        if process_streams:
            process.mix_from(process_streams)
        else:
            process.empty()
        
        recycled_streams = self.recycled_H2_streams
        if recycled_streams:
            recycled.mix_from(process_streams)
        else:
            recycled.empty()
        
        demand = process.F_mass - recycled.F_mass
        if demand >= 0:
            excess.empty()
            makeup.imass['H2'] = demand
        else:
            makeup.empty()
            excess.imass['H2'] = -demand
        
    @property
    def makeup_H2_price(self):
        '''[float] Price of the makeup H2, will be used to set the price of ins[0].'''
        return self.ins[0].price
    @makeup_H2_price.setter
    def makeup_H2_price(self, i):
        self.ins[0].price = i

    @property
    def excess_H2_price(self):
        '''[float] Price of the excess H2, will be used to set the price of outs[1].'''
        return self.outs[1].price
    @excess_H2_price.setter
    def excess_H2_price(self, i):
        self.outs[1].price = i


class ProcessWaterCenter(PWC, SanUnit):
    '''
    biosteam.facilities.ProcessWaterCenter with QSDsan properties.
    
    See Also
    --------
    `biosteam.facilities.ProcessWaterCenter <https://biosteam.readthedocs.io/en/latest/API/facilities/ProcessWaterCenter.html>`_
    '''


class Conditioning(MixTank):
    '''
    Adjust the composition and moisture content of the feedstock.
    
    Parameters
    ----------
    ins : seq(obj)
        Raw feedstock, process water for moisture adjustment.
    outs : obj
        Conditioned feedstock with appropriate composition and moisture for conversion.
    feedstock_composition : dict
        Composition of the influent feedstock,
        note that water in the feedstock will be adjusted using `target_HTL_solid_loading`.
    feedstock_dry_flowrate : float
        Feedstock dry mass flowrate for 1 reactor.
    target_HTL_solid_loading : float
        Target solid loading.
    tau : float
        Retention time for the mix tank.
    add_mixtank_kwargs : dict
        Additional keyword arguments for MixTank unit.
    '''
    _N_ins = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', F_BM_default=1,
                  feedstock_composition={ # salad dressing waste
                      'Water': 0.7566,
                      'Lipids': 0.2434*0.6245,
                      'Proteins': 0.2434*0.0238,
                      'Carbohydrates': 0.2434*0.2946,
                      'Ash': 0.2434*0.0571,
                      },
                  feedstock_dry_flowrate=1,
                  target_HTL_solid_loading=0.2,
                  tau=1, **add_mixtank_kwargs,
                  ):
        mixtank_kwargs = add_mixtank_kwargs.copy()
        mixtank_kwargs['tau'] = tau
        MixTank.__init__(self, ID, ins, outs, thermo, 
                         init_with=init_with, F_BM_default=F_BM_default, **mixtank_kwargs)
        self.feedstock_composition = feedstock_composition
        self.feedstock_dry_flowrate = feedstock_dry_flowrate
        self.target_HTL_solid_loading = target_HTL_solid_loading
        
    
    def _run(self):
        feedstock_in, htl_process_water = self.ins
        feedstock_out = self.outs[0]
        
        feedstock_composition = self.feedstock_composition
        if feedstock_composition is not None:
            for i, j in feedstock_composition.items():
                feedstock_in.imass[i] = j
        
        feedstock_dry_flowrate = self.feedstock_dry_flowrate
        feedstock_dw = 1 - feedstock_in.imass['Water']/feedstock_in.F_mass
        feedstock_in.imass['Water'] = 0
        feedstock_in.F_mass = feedstock_dry_flowrate # scale flowrate
        feedstock_in.imass['Water'] = feedstock_dry_flowrate/feedstock_dw - feedstock_dry_flowrate
              
        feedstock_out.copy_like(feedstock_in)
        total_wet = feedstock_dry_flowrate/self.target_HTL_solid_loading
        required_water = total_wet - feedstock_dry_flowrate - feedstock_in.imass['Water']
        htl_process_water.imass['Water'] = max(0, required_water)
        
        MixTank._run(self)


class Transportation(SanUnit):    
    '''
    To account for transportation cost using the price of the surrogate stream.
    The surrogate stream total mass is set to the total feedstock mass (accounting for `N_unit`),
    the price is set to `transportation_distance*transportation_distance`.
    
    Parameters
    ----------
    ins : seq(obj)
        Influent streams to be transported,
        with a surrogate flow to account for the transportation cost.
    outs : obj
        Mixture of the influent streams to be transported.        
    transportation_unit_cost : float
        Transportation cost in $/kg/km.
    transportation_distance : float
        Transportation distance in km.
    N_unit : int
        Number of parallel units.
    copy_ins_from_outs : bool
        If True, will copy influent from effluent, otherwise,
        effluent will be copied from influent.
    '''
    
    _N_ins = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', F_BM_default=1,
                  transportation_unit_cost=0,
                  transportation_distance=0,
                  N_unit=1,
                  copy_ins_from_outs=False,
                  **kwargs,
                  ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
        self.transportation_distance = transportation_distance
        self.transportation_unit_cost = transportation_unit_cost
        self.N_unit = N_unit
        self.copy_ins_from_outs = copy_ins_from_outs
        for kw, arg in kwargs.items(): setattr(self, kw, arg)
    
    def _run(self):
        inf, surrogate = self.ins
        eff = self.outs[0]
        
        if self.copy_ins_from_outs is False:
            eff.copy_like(inf)
        else:
            inf.copy_like(eff)
        
        surrogate.copy_like(inf)
        surrogate.F_mass *= self.N_unit

    def _cost(self):
        # Use the surrogate price to account for transportation cost
        self.ins[1].price = self.transportation_unit_cost * self.transportation_distance


# %%

# Jone et al., Table C-1
default_biocrude_ratios = {
    '1E2PYDIN':     0.067912,
    # 'C5H9NS':       0.010257,
    'ETHYLBEN':     0.025467,
    '4M-PHYNO':     0.050934,
    '4EPHYNOL':     0.050934,
    'INDOLE':       0.050934,
    '7MINDOLE':     0.033956,
    'C14AMIDE':     0.033956,
    'C16AMIDE':     0.152801,
    'C18AMIDE':     0.067912,
    'C16:1FA':      0.135823,
    'C16:0FA':      0.101868,
    'C18FACID':     0.016978,
    'NAPHATH':      0.050934,
    'CHOLESOL':     0.016978,
    'AROAMINE':     0.081424,
    'C30DICAD':     0.050934,
    }

class BiocrudeSplitter(SanUnit):
    '''
    Split biocrude into the respective components that meet specific boiling point
    and faction specifics.
    
    Parameters
    ----------
    ins : obj
        HTL biocrude containing the gross components.
    outs : obj
        HTL biocrude split into specific components.
    biocrude_IDs : seq(str)
        IDs of the gross components used to represent biocrude in the influent,
        will be normalized to 100% sum.
    cutoff_Tbs : Iterable(float)
        Cutoff boiling points of different fractions.
    cutoff_fracs : Iterable(float)
        Mass fractions of the different cuts, will be normalized to 100% sum.
        If there is N cutoff_Tbs, then there should be N+1 fractions.
    biocrude_ratios : dict(str, float)
        Ratios of all the components in the biocrude.
    '''
    _N_ins = _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', F_BM_default=1,
                  biocrude_IDs=('Biocrude',),
                  cutoff_Tbs=(273.15+343,), cutoff_fracs=(0.5316, 0.4684),
                  biocrude_ratios=default_biocrude_ratios,
                   **kwargs,
                  ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
        self.cutoff_Tbs = cutoff_Tbs
        self.cutoff_fracs = cutoff_fracs
        self._update_component_ratios()
        self.biocrude_IDs = biocrude_IDs
        self.biocrude_ratios = biocrude_ratios
        for kw, arg in kwargs.items(): setattr(self, kw, arg)
        
    def _update_component_ratios(self):
        '''Update the light and heavy ratios of the biocrude components.'''
        if not hasattr(self, 'cutoff_Tbs'): return
        if not hasattr(self, 'biocrude_ratios'): return

        cmps = self.components
        Tbs = self.cutoff_Tbs
        fracs = self.cutoff_fracs
        if not len(fracs)-len(Tbs) == 1:
            raise ValueError(f'Based on the number of `cutoff_Tbs` ({len(Tbs)})), '
                             f'there should be {len(Tbs)+1} `cutoff_fracs`,' 
                             f'currently there is {len(fracs)}.')
        ratios = self.biocrude_ratios.copy()

        keys = []
        frac_dcts = dict.fromkeys(fracs)
        lighter_IDs = []
        for n, Tb in enumerate(Tbs):
            frac_dct = {}
            for ID, ratio in ratios.items():
                if ID in lighter_IDs: continue
                if cmps[ID].Tb <= Tb:
                    frac_dct[ID] = ratio
                    light_key = ID
                else: 
                    keys.append((light_key, ID))
                    lighter_IDs.extend(list(frac_dct.keys()))
                    break
                    
            frac_tot = sum(frac_dct.values())
            frac_dcts[fracs[n]] = {k: v/frac_tot for k, v in frac_dct.items()}

        frac_dct_last = {k:v for k,v in ratios.items() if k not in lighter_IDs}
        frac_last_tot = sum(frac_dct_last.values())
        frac_dcts[fracs[n+1]] = {k: v/frac_last_tot for k, v in frac_dct_last.items()}
        
        self._keys = keys # light and heavy key pairs
        self._frac_dcts = frac_dcts # fractions for each cut
        
        
    def _run(self):
        biocrude_in = self.ins[0]
        biocrude_out = self.outs[0]
        
        biocrude_IDs = self.biocrude_IDs
        biocrude_out.copy_like(biocrude_in) # for the non-biocrude part, biocrude will be updated later
        
        total_crude = biocrude_in.imass[self.biocrude_IDs].sum()
        frac_dcts = self.frac_dcts
        
        for frac, dct in frac_dcts.items():
            frac_mass = frac * total_crude
            for ID, ratio in dct.items():
                biocrude_out.imass[ID] = frac_mass * ratio
        
        biocrude_out.imass[biocrude_IDs] = 0 # clear out biocrude


    @property
    def cutoff_Tbs(self):
        '''[Iterable] Boiling point cutoffs for different fractions.'''
        return self._cutoff_Tbs
    @cutoff_Tbs.setter
    def cutoff_Tbs(self, Tbs):
        try: iter(Tbs)
        except: Tbs = [Tbs]
        self._cutoff_Tbs = Tbs
        if hasattr(self, '_cutoff_fracs'):
            self._update_component_ratios()
        
    @property
    def cutoff_fracs(self):
        '''
        [Iterable] Mass fractions of the different cuts, will be normalized to 100% sum.
        If there is N cutoff_Tbs, then there should be N+1 fractions.
        '''
        return self._cutoff_fracs
    @cutoff_fracs.setter
    def cutoff_fracs(self, fracs):
        try: iter(fracs)
        except: fracs = [fracs]
        tot = sum(fracs)
        self._cutoff_fracs = [i/tot for i in fracs]
        if hasattr(self, '_cutoff_Tbs'):
            self._update_component_ratios()

    @property
    def frac_dcts(self):
        '''Fractions of the different cuts.'''
        return self._frac_dcts

    @property
    def keys(self):
        '''Light and heavy key pairs.'''
        return self._keys

    @property
    def light_component_ratios(self):
        '''Mass ratios of the components in the light fraction of the biocrude.'''
        return self._light_component_ratios

    @property
    def heavy_component_ratios(self):
        '''Mass ratios of the components in the heavy fraction of the biocrude.'''
        return self._heavy_component_ratios
    
    @property
    def light_key(self):
        '''ID of the component that has the highest boiling point in the light fraction of the biocrude.'''
        return self._light_key
    
    @property
    def heavy_key(self):
        '''ID of the component that has the lowest boiling point in the heavy fraction of the biocrude.'''
        return self._heavy_key
    
    @property
    def biocrude_ratios(self):
        '''[dict] Mass ratios of the components used to model the biocrude.'''
        return self._biocrude_ratios
    @biocrude_ratios.setter
    def biocrude_ratios(self, ratios):
        cmps = self.components
        # Sort the biocrude ratios by the boiling point
        tot = sum(ratios.values())
        ratios = {ID: ratio/tot for ID, ratio in 
                  sorted(ratios.items(), key=lambda item: cmps[item[0]].Tb)}
        self._biocrude_ratios = ratios
        self._update_component_ratios()