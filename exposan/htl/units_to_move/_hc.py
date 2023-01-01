#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from qsdsan import SanUnit, Stream
from . import Reactor, IsothermalCompressor, HXutility

__all__ = ('HC',)


# =============================================================================
# HC
# =============================================================================

class HC(Reactor):
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature (405°C)
    and pressure to produce upgraded biooil. Co-product includes fuel gas.
    
    Parameters
    ----------
    ins : Iterable(stream)
        heavy_oil, hydrogen, catalyst_in.
    outs : Iterable(stream)
        hc_out, catalyst_out.
    WHSV: float
        Weight Hourly Space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime: float
        CHG catalyst lifetime, [hr].
    hydrogen_P: float
        Hydrogen pressure, [Pa].
    hydrogen_rxned_to_heavy_oil: float
        Reacted H2 to heavy oil mass ratio.
    hydrogen_excess: float
        Actual hydrogen amount = hydrogen_rxned_to_biocrude*hydrogen_excess
    hydrocarbon_ratio: float
        Mass ratio of produced hydrocarbon to the sum of heavy oil and reacted H2.
    HCin_T: float
        HC influent temperature, [K].
    HCrxn_T: float
        HC effluent (after reaction) temperature, [K].
    HC_composition: dict
        HC effluent composition.
        
    References
    ----------
    .. [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
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
                 WHSV=0.625, # wt./hr per wt. catalyst [1]
                 catalyst_lifetime=5*7920, # 5 years [1]
                 hydrogen_P=1039.7*6894.76,
                 hydrogen_rxned_to_heavy_oil=0.01125,
                 hydrogen_excess=5.556,
                 hydrocarbon_ratio=1, # 100 wt% of heavy oil and reacted H2
                 # nearly all input heavy oils and H2 will be converted to
                 # products [1]
                 # spreadsheet HC calculation
                 HCin_T=394+273.15,
                 HCrxn_T=451+273.15,
                 HC_composition={'CO2':0.03880, 'CH4':0.00630,
                                 'CYCHEX':0.03714, 'HEXANE':0.01111,
                                 'HEPTANE':0.11474, 'OCTANE':0.08125,
                                 'C9H20':0.09086, 'C10H22':0.11756,
                                 'C11H24':0.16846, 'C12H26':0.13198,
                                 'C13H28':0.09302, 'C14H30':0.04643,
                                 'C15H32':0.03250, 'C16H34':0.01923,
                                 'C17H36':0.00431, 'C18H38':0.00099,
                                 'C19H40':0.00497, 'C20H42':0.00033},
                 #combine C20H42 and PHYTANE as C20H42
                 # will not be a variable in uncertainty/sensitivity analysis
                 P=None, tau=5, void_fraciton=0.4, # Towler
                 length_to_diameter=2, N=None, V=None, auxiliary=False, mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1.5,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.hydrogen_P = hydrogen_P
        self.hydrogen_rxned_to_heavy_oil = hydrogen_rxned_to_heavy_oil
        self.hydrogen_excess = hydrogen_excess
        self.hydrocarbon_ratio = hydrocarbon_ratio
        self.HCin_T = HCin_T
        self.HCrxn_T = HCrxn_T
        self.HC_composition = HC_composition
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
        
        catalyst_in.imass['HC_catalyst'] = heavy_oil.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts
        
        hydrogen.imass['H2'] = heavy_oil.F_mass*self.hydrogen_rxned_to_heavy_oil*self.hydrogen_excess
        hydrogen.phase = 'g'

        hydrocarbon_mass = heavy_oil.F_mass*(1 +\
                           self.hydrogen_rxned_to_heavy_oil)*\
                           self.hydrocarbon_ratio

        hc_out.phase = 'g'

        for name, ratio in self.HC_composition.items():
            hc_out.imass[name] = hydrocarbon_mass*ratio
        
        hc_out.imass['H2'] = heavy_oil.F_mass*self.hydrogen_rxned_to_heavy_oil*(self.hydrogen_excess - 1)
        
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