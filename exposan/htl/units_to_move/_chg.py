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

from biosteam.units.decorators import cost
from qsdsan import SanUnit, Stream, CEPCI_by_year
from qsdsan.utils import auom
from . import Reactor, HXutility, Pump

__all__ = ('CHG',)

_lb_to_kg = auom('lb').conversion_factor('kg')

# =============================================================================
# CHG
# =============================================================================

# hydrocyclone
@cost(basis='Treatment capacity', ID='Hydrocyclone', units='lb/h',
      cost=5000000, S=968859,
      CE=CEPCI_by_year[2009], n=0.65, BM=2.1)
class CHG(Reactor):
    '''
    CHG serves to reduce the COD content in the aqueous phase and produce fuel
    gas under elevated temperature (350°C) and pressure. The outlet will be
    cooled down and separated by a flash unit.
    
    Parameters
    ----------
    ins : Iterable(stream)
        chg_in, catalyst_in.
    outs : Iterable(stream)
        chg_out, catalyst_out.
    pump_pressure: float
        CHG influent pressure, [Pa].
    heat_temp: float
        CHG influent temperature, [K].
    cool_temp: float
        CHG effluent temperature, [K].
    WHSV: float
        Weight Hourly Space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime: float
        CHG catalyst lifetime, [hr].
    gas_composition: dict
        CHG gas composition.
    gas_C_2_total_C: dict
        CHG gas carbon content to feed carbon content.
    CAPEX_factor: float
        Factor used to adjust CAPEX.
        
    References
    ----------
    .. [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    .. [2] Davis, R. E.; Grundl, N. J.; Tao, L.; Biddy, M. J.; Tan, E. C.;
        Beckham, G. T.; Humbird, D.; Thompson, D. N.; Roni, M. S. Process
        Design and Economics for the Conversion of Lignocellulosic Biomass
        to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case
        Update; Biochemical Deconstruction and Conversion of Biomass to Fuels
        and Products via Integrated Biorefinery Pathways; NREL/TP--5100-71949,
        1483234; 2018; p NREL/TP--5100-71949, 1483234.
        https://doi.org/10.2172/1483234.
    .. [3] Elliott, D. C.; Neuenschwander, G. G.; Hart, T. R.; Rotness, L. J.;
        Zacher, A. H.; Santosa, D. M.; Valkenburg, C.; Jones, S. B.;
        Rahardjo, S. A. T. Catalytic Hydrothermal Gasification of Lignin-Rich
        Biorefinery Residues and Algae Final Report. 87.
    '''
    _N_ins = 2
    _N_outs = 2
    
    _F_BM_default = {**Reactor._F_BM_default,
                      'Heat exchanger': 3.17,
                      'Sulfur guard': 2.0}
    _units= {'Treatment capacity': 'lb/h', # hydrocyclone
              'Hydrocyclone weight': 'lb'}
    
    auxiliary_unit_names=('pump','heat_ex_heating','heat_ex_cooling')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='Stream',
                  pump_pressure=3089.7*6894.76,
                  heat_temp=350+273.15,
                  cool_temp=60+273.15,
                  WHSV=3.562,
                  catalyst_lifetime=7920, # 1 year [1]
                  gas_composition={'CH4':0.527,
                                  'CO2':0.432,
                                  'C2H6':0.011,
                                  'C3H8':0.030,
                                  'H2':0.0001}, # [1]
                  gas_C_2_total_C=0.5981, # [1]
                  P=None, tau=20/60, void_fraction=0.5, # [2, 3]
                  length_to_diameter=2, N=6, V=None, auxiliary=False,
                  mixing_intensity=None, kW_per_m3=0,
                  wall_thickness_factor=1,
                  vessel_material='Stainless steel 316',
                  vessel_type='Vertical',
                  CAPEX_factor=1):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        
        self.pump_pressure = pump_pressure
        self.heat_temp = heat_temp
        self.cool_temp = cool_temp
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.gas_composition = gas_composition
        self.gas_C_2_total_C = gas_C_2_total_C
        pump_in = Stream(f'{ID}_pump_in')
        pump_out = Stream(f'{ID}_pump_out')
        self.pump = Pump(ID=f'.{ID}_pump', ins=pump_in, outs=pump_out, P=pump_pressure)
        hx_ht_in = Stream(f'{ID}_hx_ht_in')
        hx_ht_out = Stream(f'{ID}_hx_ht_out')
        self.heat_ex_heating = HXutility(ID=f'.{ID}_hx_ht', ins=hx_ht_in, outs=hx_ht_out, T=heat_temp, rigorous=True)
        hx_cl_in = Stream(f'{ID}_hx_cl_in')
        hx_cl_out = Stream(f'{ID}_hx_cl_out')
        self.heat_ex_cooling = HXutility(ID=f'.{ID}_hx_cl', ins=hx_cl_in, outs=hx_cl_out, T=cool_temp, rigorous=True)
        self.P = P
        self.tau = tau
        self.V_wf = void_fraction
        # no headspace, gases produced will be vented, so V_wf = void fraction [2, 3]
        self.length_to_diameter = length_to_diameter
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
        
        chg_in, catalyst_in = self.ins
        chg_out, catalyst_out = self.outs
        
        catalyst_in.imass['CHG_catalyst'] = chg_in.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts
            
        cmps = self.components
        gas_C_ratio = 0
        for name, ratio in self.gas_composition.items():
            gas_C_ratio += ratio*cmps[name].i_C
            
        gas_mass = chg_in.imass['C']*self.gas_C_2_total_C/gas_C_ratio
        
        for name,ratio in self.gas_composition.items():
            chg_out.imass[name] = gas_mass*ratio
                
        chg_out.imass['H2O'] = chg_in.F_mass - gas_mass
        # all C, N, and P are accounted in H2O here, but will be calculated as properties.
        
        chg_out.T = self.cool_temp
        chg_out.P = self.pump_pressure
        
    @property
    def CHGout_C(self):
        # not include carbon in gas phase
        return self.ins[0].imass['C']*(1 - self.gas_C_2_total_C)
    
    @property
    def CHGout_N(self):
        return self.ins[0].imass['N']
    
    @property
    def CHGout_P(self):
        return self.ins[0].imass['P']
        
    def _design(self):
        Design = self.design_results
        Design['Treatment capacity'] = self.ins[0].F_mass/_lb_to_kg
        
        pump = self.pump
        pump.ins[0].copy_like(self.ins[0])
        pump.simulate()
        
        hx_ht = self.heat_ex_heating
        hx_ht_ins0, hx_ht_outs0 = hx_ht.ins[0], hx_ht.outs[0]
        hx_ht_ins0.copy_like(self.ins[0])
        hx_ht_outs0.copy_like(hx_ht_ins0)
        hx_ht_ins0.T = self.ins[0].T
        hx_ht_outs0.T = hx_ht.T
        hx_ht_ins0.P = hx_ht_outs0.P = pump.P
        
        hx_ht_ins0.vle(T=hx_ht_ins0.T, P=hx_ht_ins0.P)
        hx_ht_outs0.vle(T=hx_ht_outs0.T, P=hx_ht_outs0.P)
        
        hx_ht.simulate_as_auxiliary_exchanger(ins=hx_ht.ins, outs=hx_ht.outs)
            
        hx_cl = self.heat_ex_cooling
        hx_cl_ins0, hx_cl_outs0 = hx_cl.ins[0], hx_cl.outs[0]
        hx_cl_ins0.copy_like(self.outs[0])
        hx_cl_outs0.copy_like(hx_cl_ins0)
        hx_cl_ins0.T = hx_ht.T
        hx_cl_outs0.T = hx_cl.T
        hx_cl_ins0.P = hx_cl_outs0.P = self.outs[0].P

        hx_cl_ins0.vle(T=hx_cl_ins0.T, P=hx_cl_ins0.P)
        hx_cl_outs0.vle(T=hx_cl_outs0.T, P=hx_cl_outs0.P)        
        
        hx_cl.simulate_as_auxiliary_exchanger(ins=hx_cl.ins, outs=hx_cl.outs)

        self.P = self.pump_pressure
        Reactor._design(self)
        Design['Hydrocyclone weight'] = 0.3*Design['Weight']*Design['Number of reactors'] # assume stainless steel
        # based on [1], page 54, the purchase price of hydrocyclone to the purchase price of CHG
        # reactor is around 0.3, therefore, assume the weight of hydrocyclone is 0.3*single CHG weight*number of CHG reactors
        self.construction[0].quantity += Design['Hydrocyclone weight']*_lb_to_kg
    
    def _cost(self):
        Reactor._cost(self)
        purchase_costs = self.baseline_purchase_costs
        current_cost = 0 # cost w/o sulfur guard
        for item in purchase_costs.keys():
            current_cost += purchase_costs[item]
        purchase_costs['Sulfur guard'] = current_cost*0.05
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