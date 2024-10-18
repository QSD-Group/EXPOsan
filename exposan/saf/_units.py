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
    'Hydrotreating',
    )

_lb_to_kg = 0.453592
_m3_to_gal = 264.172
_in_to_m = 0.0254
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
                 P=3049.7*6894.76, tau=0, V_wf=0,
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

    def _normalize_composition(self, dct):
        total = sum(dct.values())
        if total <=0: raise ValueError(f'Sum of total yields/composition should be positive, not {total}.')
        return {k:v/total for k, v in dct.items()}

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', include_construction=False,
                 T=350+273.15,
                 P=None,
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
        self.eff_at_temp = Stream(f'{ID}_eff_at_temp')
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
        
        gas.phase = 'g'
        char.phase = 's'
        aq.phase = crude.phase = 'l'
        
        self.eff_at_temp.mix_from(outs)
        
        for i in outs:
            i.T = self.eff_T
            i.P = self.eff_P
    
    def _design(self):
        Design = self.design_results
        Design['Mass flow'] = self.ins[0].F_mass/_lb_to_kg
        
        feed = self.ins[0]
        self.P = self.P or feed.P
        
        # Influent heating to HTL conditions
        inf_hx = self.inf_hx
        inf_hx_in, inf_hx_out = inf_hx.ins[0], inf_hx.outs[0]
        inf_hx_in.copy_like(feed)
        inf_hx_out.copy_flow(inf_hx_in)
        inf_hx_out.T = self.T
        inf_hx_out.P = self.P
        inf_hx.simulate_as_auxiliary_exchanger(ins=inf_hx.ins, outs=inf_hx.outs)
        
        # Effluent cooling to near ambient conditions
        eff_hx = self.eff_hx
        eff_hx_in, eff_hx_out = eff_hx.ins[0], eff_hx.outs[0]
        eff_hx_in.copy_like(self.eff_at_temp)
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
        self._decorated_cost()
        
        purchase_costs = self.baseline_purchase_costs
        for item in purchase_costs.keys():
            purchase_costs[item] *= self.CAPEX_factor
        
        for aux_unit in self.auxiliary_units:
            purchase_costs = aux_unit.baseline_purchase_costs
            installed_costs = aux_unit.installed_costs
            for item in purchase_costs.keys():
                purchase_costs[item] *= self.CAPEX_factor
                installed_costs[item] *= self.CAPEX_factor
                

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
# Hydrocracking
# =============================================================================

#!!! Hydrocracking and hydrotreating can be potentially combined
class Hydrocracking(Reactor):
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature (405°C)
    and pressure to produce upgraded biooil. Co-product includes fuel gas.
    
    Parameters
    ----------
    ins : Iterable(stream)
        Influent heavy oil, hydrogen, catalyst_in.
    outs : Iterable(stream)
        Crakced oil, catalyst_out.
    WHSV: float
        Weight hourly space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime: float
        HC catalyst lifetime, [hr].
    catalyst_ID : str
        ID of the catalyst.
    hydrogen_P: float
        Hydrogen pressure, [Pa].
    hydrogen_rxned_to_inf_oil: float
        Reacted H2 to influent oil mass ratio.
    hydrogen_ratio : float
        Actual hydrogen amount = hydrogen_rxned_to_biocrude*hydrogen_ratio
    gas : float
        Mass ratio of fuel gas to the sum of influent oil and reacted H2.
    oil_yield : float
        Mass ratio of cracked oil to the sum of influent oil and reacted H2.
    HCin_T: float
        HC influent temperature, [K].
    HCrxn_T: float
        HC effluent (after reaction) temperature, [K].
    gas_composition: dict
        Composition of the gas products (excluding excess H2), will be normalized to 100% sum.
    oil_composition: dict
        Composition of the cracked oil, will be normalized to 100% sum.
    aq_composition: dict
        Composition of the aqueous product, yield will be calculated as 1-gas-oil.
        
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
                 WHSV=0.625, # wt./hr per wt. catalyst [1]
                 catalyst_lifetime=5*7920, # 5 years [1]
                 catalyst_ID='HC_catalyst',
                 hydrogen_P=1039.7*6894.76,
                 hydrogen_rxned_to_inf_oil=0.01125,
                 hydrogen_ratio=5.556,
                # 100 wt% of heavy oil and reacted H2
                # nearly all input heavy oils and H2 will be converted to products [1]
                # spreadsheet HC calculation
                 gas_yield=0.03880-0.00630,
                 oil_yield=1-0.03880-0.00630,
                 HCin_T=394+273.15,
                 HCrxn_T=451+273.15,
                 gas_composition={'CO2':0.03880, 'CH4':0.00630,},
                 oil_composition={
                    'CYCHEX':0.03714, 'HEXANE':0.01111,
                    'HEPTANE':0.11474, 'OCTANE':0.08125,
                    'C9H20':0.09086, 'C10H22':0.11756,
                    'C11H24':0.16846, 'C12H26':0.13198,
                    'C13H28':0.09302, 'C14H30':0.04643,
                    'C15H32':0.03250, 'C16H34':0.01923,
                    'C17H36':0.00431, 'C18H38':0.00099,
                    'C19H40':0.00497, 'C20H42':0.00033,
                    },
                 aq_composition={'Water':1},
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
        self.hydrogen_ratio = hydrogen_ratio
        self.gas_yield = gas_yield
        self.oil_yield = oil_yield
        self.HCin_T = HCin_T
        self._mixed_in = Stream(f'{ID}_mixed_in')
        self.HCrxn_T = HCrxn_T
        self.gas_composition = gas_composition
        self.oil_composition = oil_composition
        self.aq_composition = aq_composition
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
        inf_oil, hydrogen, catalyst_in = self.ins
        cracked_oil, catalyst_out = self.outs
        
        catalyst_in.imass[self.catalyst_ID] = inf_oil.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts
        
        hydrogen_rxned_to_inf_oil = self.hydrogen_rxned_to_inf_oil
        hydrogen_ratio = self.hydrogen_ratio
        H2_rxned =  inf_oil.F_mass*hydrogen_rxned_to_inf_oil
        hydrogen.imass['H2'] = inf_oil.F_mass*hydrogen_ratio
        hydrogen.phase = 'g'

        cracked_oil.empty()
        cracked_oil.imass[self.eff_composition.keys()] = self.eff_composition.values()
        cracked_oil.F_mass = inf_oil.F_mass*(1 + hydrogen_rxned_to_inf_oil)
        
        cracked_oil.imass['H2'] = H2_rxned*(hydrogen_ratio - 1)
        
        cracked_oil.P = inf_oil.P
        cracked_oil.T = self.HCrxn_T
        
        cracked_oil.vle(T=cracked_oil.T, P=cracked_oil.P)
        
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
        
    def _normalize_composition(self, dct):
        total = sum(dct.values())
        if total <=0: raise ValueError(f'Sum of total yields/composition should be positive, not {total}.')
        return {k:v/total for k, v in dct.items()}
    
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
    def aq_composition(self):
        return self._aq_composition
    @aq_composition.setter
    def aq_composition(self, comp_dct):
        self._aq_composition = self._normalize_composition(comp_dct)
        
        
    @property
    def eff_composition(self):
        '''Composition of products, normalized to 100% sum.'''
        gas_composition = self.gas_composition
        oil_composition = self.oil_composition
        aq_composition = self.aq_composition
        oil_yield = self.oil_yield
        gas_yield = self.gas_yield
        aq_yield = self.aq_yield
        eff_composition = {k:v*gas_yield for k, v in gas_composition.items()}
        eff_composition.update({k:v*oil_yield for k, v in oil_composition.items()})
        eff_composition.update({k:v*aq_yield for k, v in aq_composition.items()})
        return self._normalize_composition(eff_composition)

    # @property
    # def C_balance(self):
    #     '''Total carbon in the outs over total in the ins.'''
    #     cmps = self.components
    #     C_in = sum(self.ins[0].imass[cmp.ID]*cmp.i_C for cmp in cmps)
    #     C_out = sum(self.outs[0].imass[cmp.ID]*cmp.i_C for cmp in cmps)
    #     return C_out/C_in

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
        
        V_H2 = self.ins[1].F_vol/self.hydrogen_ratio*101325/self.hydrogen_P
        # just account for reacted H2
        V_biocrude = self.ins[0].F_vol
        self.V_wf = self.void_fraciton*V_biocrude/(V_biocrude + V_H2)
        Reactor._design(self)
    
# =============================================================================
# Hydrotreating
# =============================================================================

class Hydrotreating(Reactor):
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature
    and pressure to produce upgraded biooil. Co-product includes fuel gas.
    
    Parameters
    ----------
    ins : Iterable(stream)
        Influent oil, hydrogen, catalyst_in.
    outs : Iterable(stream)
        Treated oil, catalyst_out.
    WHSV: float
        Weight hourly space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime: float
        HT catalyst lifetime, [hr].
    hydrogen_P: float
        Hydrogen pressure, [Pa].
    hydrogen_rxned_to_inf_oil: float
        Reacted H2 to influent oil mass ratio.
    hydrogen_ratio: float
        Actual hydrogen amount = hydrogen_rxned_to_biocrude*hydrogen_ratio
    gas_yield: float
        Mass ratio of gas to the sum of influent oil and reacted H2.
    oil_yield: float
        Mass ratio of treated oil to the sum of influent oil and reacted H2.
    HTin_T: float
        HT influent temperature, [K].
    HTrxn_T: float
        HT effluent (after reaction) temperature, [K].
    gas_composition: dict
        Composition of the gas products (excluding excess H2), will be normalized to 100% sum.
    oil_composition: dict
        Composition of the cracked oil, will be normalized to 100% sum.
    CAPEX_factor: float
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
        Eds.; Butterworth-Heinemann: Boston, 2013; pp 563–629.
        https://doi.org/10.1016/B978-0-08-096659-5.00014-6.
    '''
    _N_ins = 3
    _N_outs = 2
    auxiliary_unit_names=('compressor','heat_exchanger',)
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Compressor': 1.1}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 WHSV=0.625, # wt./hr per wt. catalyst [1]
                 catalyst_lifetime=2*7920, # 2 years [1]
                 catalyst_ID='HT_catalyst',
                 hydrogen_P=1530*6894.76,
                 hydrogen_rxned_to_inf_oil=0.046,
                 hydrogen_ratio=3,
                 # gas and oil combined is 87.5 wt% of biocrude and reacted H2 [1]
                 # spreadsheet HT calculation
                 gas_yield=0.07707875,
                 oil_yield=0.875-0.07707875,
                 HTin_T=174+273.15,
                 HTrxn_T=402+273.15, # [1]
                 gas_composition={
                     'CH4':0.02280, 'C2H6':0.02923,
                     'C3H8':0.01650, 'C4H10':0.00870,
                     'TWOMBUTAN':0.00408, 'NPENTAN':0.00678,
                     },
                 oil_composition={
                    'TWOMPENTA':0.00408, 'HEXANE':0.00401,
                    'TWOMHEXAN':0.00408, 'HEPTANE':0.00401,
                    'CC6METH':0.01020, 'PIPERDIN':0.00408,
                    'TOLUENE':0.01013, 'THREEMHEPTA':0.01020,
                    'OCTANE':0.01013, 'ETHCYC6':0.00408,
                    'ETHYLBEN':0.02040, 'OXYLENE':0.01020,
                    'C9H20':0.00408, 'PROCYC6':0.00408,
                    'C3BENZ':0.01020, 'FOURMONAN':0,
                    'C10H22':0.00240, 'C4BENZ':0.01223,
                    # C10H22 was originally 0.00203, but it is not
                    # good for distillation column, the excess amount
                    # is substracted from HEXANE, HEPTANE, TOLUENE,
                    # OCTANE, and C9H20, which were originally 0.00408,
                    # 0.00408, 0.01020, 0.01020, and 0.00408
                    'C11H24':0.02040, 'C10H12':0.02040,
                    'C12H26':0.02040, 'OTTFNA':0.01020,
                    'C6BENZ':0.02040, 'OTTFSN':0.02040,
                    'C7BENZ':0.02040, 'C8BENZ':0.02040,
                    'C10H16O4':0.01837, 'C15H32':0.06120,
                    'C16H34':0.18360, 'C17H36':0.08160, 
                    'C18H38':0.04080, 'C19H40':0.04080,
                    'C20H42':0.10200, 'C21H44':0.04080,
                    'TRICOSANE':0.04080, 'C24H38O4':0.00817,
                    'C26H42O4':0.01020, 'C30H62':0.00203, # [1]              
                    },
                 aq_composition={'Water':1},
                 # spreadsheet HT calculation
                 # will not be a variable in uncertainty/sensitivity analysis
                 P=None, tau=0.5, void_fraciton=0.4, # [2]
                 length_to_diameter=2, diameter=None,
                 N=None, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 CAPEX_factor=1,):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.catalyst_ID = catalyst_ID
        self.hydrogen_P = hydrogen_P
        self.hydrogen_rxned_to_inf_oil = hydrogen_rxned_to_inf_oil
        self.hydrogen_ratio = hydrogen_ratio
        self.gas_yield = gas_yield
        self.oil_yield = oil_yield
        self.HTin_T = HTin_T
        self._mixed_in = Stream(f'{ID}_mixed_in')
        self.HTrxn_T = HTrxn_T
        self.gas_composition = gas_composition
        self.oil_composition = oil_composition
        self.aq_composition = aq_composition
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
        
    def _run(self):
        inf_oil, hydrogen, catalyst_in = self.ins
        treated_oil, catalyst_out = self.outs
        
        catalyst_in.imass[self.catalyst_ID] = inf_oil.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts
        
        hydrogen_rxned_to_inf_oil = self.hydrogen_rxned_to_inf_oil
        hydrogen_ratio = self.hydrogen_ratio
        H2_rxned =  inf_oil.F_mass*hydrogen_rxned_to_inf_oil
        hydrogen.imass['H2'] = H2_rxned*hydrogen_ratio
        hydrogen.phase = 'g'
        
        treated_oil.empty()
        treated_oil.imass[self.eff_composition.keys()] = self.eff_composition.values()
        treated_oil.F_mass = inf_oil.F_mass*(1 + hydrogen_rxned_to_inf_oil)
            
        treated_oil.imass['H2'] = H2_rxned*(hydrogen_ratio - 1)
        
        treated_oil.P = inf_oil.P        
        treated_oil.T = self.HTrxn_T
        
        treated_oil.vle(T=treated_oil.T, P=treated_oil.P)
        

    _normalize_yields = Hydrocracking._normalize_yields
    _normalize_composition = Hydrocracking._normalize_composition
    gas_yield = Hydrocracking.gas_yield
    oil_yield = Hydrocracking.oil_yield
    aq_yield = Hydrocracking.aq_yield
    gas_composition = Hydrocracking.gas_composition
    oil_composition = Hydrocracking.oil_composition
    aq_composition = Hydrocracking.aq_composition
    eff_composition = Hydrocracking.eff_composition


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
        if not self.HTin_T: self.HTin_T = self._mixed_in.T
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
        
        V_H2 = self.ins[1].F_vol/self.hydrogen_ratio*101325/self.hydrogen_P
        # just account for reacted H2
        V_biocrude = self.ins[0].F_vol
        self.V_wf = self.void_fraciton*V_biocrude/(V_biocrude + V_H2)
        Reactor._design(self)
        
    
    def _cost(self):
        Reactor._cost(self)
        purchase_costs = self.baseline_purchase_costs
        CAPEX_factor = self.CAPEX_factor        
        for item in purchase_costs.keys():
            purchase_costs[item] *= CAPEX_factor


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