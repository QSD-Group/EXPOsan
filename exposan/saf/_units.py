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

from biosteam.units.decorators import cost
from biosteam.units.design_tools import CEPCI_by_year
from qsdsan import SanUnit, Stream
from qsdsan.sanunits import Reactor, IsothermalCompressor, HXutility

__all__ = (
    'HydrothermalLiquefaction',
    'Hydroprocessing',
    )

_lb_to_kg = 0.453592
_m3_to_gal = 264.172
_in_to_m = 0.0254
_psi_to_Pa = 6894.76
_m3perh_to_mmscfd = 1/1177.17 # H2


# %%

# =============================================================================
# KOdrum
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
    :class:`qsdsan.sanunits.HydrothermalLiquefaction`
    '''
    _N_ins = 3
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    _F_BM_default = {
        **Reactor._F_BM_default,
        'Horizontal pressure vessel': 1.5,
        'Vertical pressure vessel': 1.5,
        }
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', include_construction=False,
                 P=3049.7*_psi_to_Pa, tau=0, V_wf=0,
                 length_to_diameter=2, diameter=None,
                 N=4, V=None,
                 auxiliary=True,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',):
        # drum_steel_cost_factor: so the cost matches [1]
        # when do comparison, if fully consider scaling factor (2000 tons/day to 100 tons/day),
        # drum_steel_cost_factor should be around 3
        # but that is too high, we use 1.5 instead.
        
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
    
    def _run(self):
        pass
    
    def _cost(self):
        Reactor._cost(self)



# =============================================================================
# HTL
# =============================================================================

# separator
@cost(basis='Mass flow', ID='Solids filter oil/water separator', units='lb/h',
      cost=3945523, S=1219765,
      CE=CEPCI_by_year[2011], n=0.68, BM=1.9)
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
        Temperature of the HTL reaction, K.
    P : float
        Pressure of the HTL reaction, Pa.
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
    eff_T: float
        HTL effluent temperature, K.
    eff_P: float
        HTL effluent pressure, Pa.
    CAPEX_factor: float
        Factor used to adjust the total installed cost,
        this is on top of all other factors to individual equipment of this unit
        (e.g., bare module, material factors).
        
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
    _units= {'Mass flow': 'lb/h',
             'Solid filter and separator weight': 'lb'}
    
    auxiliary_unit_names=('inf_hx', 'eff_hx','kodrum')

    _F_BM_default = {
        **Reactor._F_BM_default,
        'Heat exchanger': 3.17,
        'Horizontal pressure vessel': 2.7, # so the cost matches [6]
        'Vertical pressure vessel': 2.7, # so the cost matches [6]
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
                 eff_T=60+273.15, # [4]
                 eff_P=30*6894.76, # [4], all set to 30 psi
                 tau=15/60, V_wf=0.45,
                 length_to_diameter=None,
                 diameter=6.875*_in_to_m,
                 N=4, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Horizontal',
                 CAPEX_factor=1,
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
        inf_hx_in = Stream(f'{ID}_inf_hx_in')
        inf_hx_out = Stream(f'{ID}_inf_hx_out')
        self.inf_hx = HXutility(ID=f'.{ID}_inf_hx', ins=inf_hx_in, outs=inf_hx_out, T=T, rigorous=True)
        eff_hx_in = Stream(f'{ID}_eff_hx_in')
        eff_hx_out = Stream(f'{ID}_eff_hx_out')
        self._eff_at_temp = Stream(f'{ID}_eff_at_temp')
        self.eff_hx = HXutility(ID=f'.{ID}_eff_hx', ins=eff_hx_in, outs=eff_hx_out, T=eff_T, rigorous=True)
        self.kodrum = KnockOutDrum(ID=f'.{ID}_KOdrum')
        self.eff_T = eff_T
        self.eff_P = eff_P
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
            for k, v in comp.items():
                out.imass[k] = v
                
        dw_yields = self.dw_yields
        gas.F_mass = tot_dw * dw_yields['gas']
        aq.F_mass = tot_dw * dw_yields['aqueous']
        crude.F_mass = tot_dw * dw_yields['biocrude']
        char.F_mass = tot_dw * dw_yields['char']
        aq.imass['Water'] = feed.imass['Water'] - sum(i.imass['Water'] for i in (gas, crude, char))
        
        self._eff_at_temp.mix_from(outs)
        
        gas.phase = 'g'
        char.phase = 's'
        aq.phase = crude.phase = 'l'
        
        for i in outs:
            i.T = self.eff_T
            i.P = self.eff_P
    
    def _design(self):
        Design = self.design_results
        Design['Mass flow'] = self.ins[0].F_mass/_lb_to_kg
        
        feed = self.ins[0]
        
        # Influent heating to HTL conditions
        inf_hx = self.inf_hx
        inf_hx_in, inf_hx_out = inf_hx.ins[0], inf_hx.outs[0]
        inf_hx_in.copy_like(feed)
        inf_hx_out.copy_flow(inf_hx_in)
        inf_hx_out.T = self.T
        inf_hx_out.P = self.P # this may lead to HXN error, when at pressure
        inf_hx.simulate_as_auxiliary_exchanger(ins=inf_hx.ins, outs=inf_hx.outs)
        
        # Effluent cooling to near ambient conditions
        eff_hx = self.eff_hx
        eff_hx_in, eff_hx_out = eff_hx.ins[0], eff_hx.outs[0]
        eff_hx_in.copy_like(self._eff_at_temp)
        eff_hx_out.mix_from(self.outs)
        eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs)
        
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
        if hasattr(self, '_decorated_cost'): self._decorated_cost()
        
        purchase_costs = self.baseline_purchase_costs
        for item in purchase_costs.keys():
            purchase_costs[item] *= self.CAPEX_factor
        
        for aux_unit in self.auxiliary_units:
            purchase_costs = aux_unit.baseline_purchase_costs
            installed_costs = aux_unit.installed_costs
            for item in purchase_costs.keys():
                purchase_costs[item] *= self.CAPEX_factor
                installed_costs[item] *= self.CAPEX_factor
                
    def _normalize_composition(self, dct):
        total = sum(dct.values())
        if total <=0: raise ValueError(f'Sum of total yields/compositions should be positive, not {total}.')
        return {k:v/total for k, v in dct.items()}
    
    @property
    def yields(self):
        return self._yields
    @yields.setter
    def yields(self, comp_dct):
        self._yields = self._normalize_composition(comp_dct)

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
        """Higher heating value of the biocrude, MJ/kg."""
        crude = self.outs[2]
        return crude.HHV/crude.F_mass/1e3
    
    @property
    def energy_recovery(self):
        """Energy recovery calculated as the HHV of the biocrude over the HHV of the feedstock."""
        feed = self.ins[0]
        return self.biocrude_HHV/(feed.HHV/feed.F_mass/1e3)



# =============================================================================
# Hydroprocessing
# =============================================================================

class Hydroprocessing(Reactor):
    '''
    For fuel upgrading processes such as hydrocracking and hydrotreating.
    Co-product includes fuel gas and aqueous stream.
    
    Note that addition units are needed to fractionate the product
    into the gas, aqueous, and oil streams.
    
    Parameters
    ----------
    ins : Iterable(stream)
        Influent crude oil, hydrogen, catalyst_in.
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
    CAPEX_factor: float
        Factor used to adjust the total installed cost,
        this is on top of all other factors to individual equipment of this unit
        (e.g., bare module, material factors).
        
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
    
    auxiliary_unit_names=('compressor','heat_exchanger',)
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Compressor': 1.1}
    
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
                 tau=15/60, # set to the same as HTL as in [1]
                 V_wf=0.4, # void_fraciton=0.4, # Towler
                 length_to_diameter=2, diameter=None,
                 N=None, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1.5,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 CAPEX_factor=1.,):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, include_construction=include_construction)
        self.T = T
        self.P = P
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.catalyst_ID = catalyst_ID
        self.hydrogen_rxned_to_inf_oil = hydrogen_rxned_to_inf_oil
        self.hydrogen_ratio = hydrogen_ratio
        self.gas_yield = gas_yield
        self.oil_yield = oil_yield
        self.gas_composition = gas_composition
        self.oil_composition = oil_composition
        self.aqueous_composition = aqueous_composition
        self.CAPEX_factor = CAPEX_factor
        # For H2 compressing
        IC_in = Stream(f'{ID}_IC_in')
        IC_out = Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in,
                                               outs=IC_out, P=P)
        # For influent heating
        self._mixed_in = Stream(f'{ID}_mixed_in')
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out)
        self.tau = tau
        self._V_wf = V_wf # will be adjusted later
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
        inf_oil, hydrogen, catalyst_in = self.ins
        eff_oil, catalyst_out = self.outs
        
        catalyst_in.imass[self.catalyst_ID] = inf_oil.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        
        hydrogen_rxned_to_inf_oil = self.hydrogen_rxned_to_inf_oil
        hydrogen_ratio = self.hydrogen_ratio
        H2_rxned =  inf_oil.F_mass*hydrogen_rxned_to_inf_oil
        hydrogen.imass['H2'] = H2_rxned*hydrogen_ratio
        hydrogen.phase = 'g'

        eff_oil.copy_like(inf_oil)
        eff_oil.empty()
        eff_oil.imass[self.eff_composition.keys()] = self.eff_composition.values()
        eff_oil.F_mass = inf_oil.F_mass*(1 + hydrogen_rxned_to_inf_oil)
        eff_oil.imass['H2'] = H2_rxned*(hydrogen_ratio - 1)
        
        eff_oil.P = self.P
        eff_oil.T = self.T
        eff_oil.vle(T=eff_oil.T, P=eff_oil.P)


    def _design(self):
        IC = self.compressor # for H2 compressing
        H2 = self.ins[1]
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]
        IC_ins0.copy_like(H2)
        IC_outs0.copy_like(H2)
        IC_outs0.P = IC.P = self.P
        IC_ins0.phase = IC_outs0.phase = 'g'
        IC.simulate()
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from(self.ins)
        hx_outs0.copy_like(hx_ins0)
        hx_outs0.T = self.T
        hx_ins0.P = hx_outs0.P = IC_outs0.P
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        
        V_oil = self.ins[0].F_vol
        V_H2 = H2.F_vol*(H2.P/self.P)
        self.V_wf = self._V_wf*V_oil/(V_oil + V_H2) # account for the added H2
        Reactor._design(self)

    _cost = HydrothermalLiquefaction._cost

        
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

    # @property
    # def C_balance(self):
    #     '''Total carbon in the outs over total in the ins.'''
    #     cmps = self.components
    #     C_in = sum(self.ins[0].imass[cmp.ID]*cmp.i_C for cmp in cmps)
    #     C_out = sum(self.outs[0].imass[cmp.ID]*cmp.i_C for cmp in cmps)
    #     return C_out/C_in


# =============================================================================
# Pressure Swing Adsorption
# =============================================================================

@cost(basis='H2 flowrate', ID='PSA', units='mmscfd',
      cost=1750000, S=10,
      CE=CEPCI_by_year[2004], n=0.8, BM=2.47)
class PressureSwingAdsorption:
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
    _ins_size_is_fixed = False
    _units = {'H2 flowrate': 'mmscfd',}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='WasteStream', efficiency=0.9,):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.efficiency = efficiency
    
    @property
    def efficiency (self):
        return self._efficiency 
    @efficiency.setter
    def efficiency(self, i):
        if i > 1: raise Exception('Efficiency cannot be larger than 1.')
        self._efficiency  = i
        
    def _run(self):
        H2, others = self.outs       
        others.mix_from(self.ins)
        
        H2.imass['H2'] = recovered = others.imass['H2'] * self.efficiency
        others.imass['H2'] -= recovered
        
    def _design(self):
        self.design_results['Hydrogen_PSA'] = self.F_vol_in*_m3perh_to_mmscfd*(101325/self.outs[0].P)