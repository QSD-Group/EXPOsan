# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 08:54:00 2024

@author: joy_c
"""
import qsdsan as qs
from qsdsan import SanUnit, Component, Components
from biosteam.units import BatchBioreactor
from biorefineries.cane import create_sugarcane_chemicals

__all__ = ('HApFermenter',)

def create_hap_cmps(set_thermo=True):
    cmps_df = Components.load_default()
    chems = create_sugarcane_chemicals()
    Yeast = Component.from_chemical('Yeast', chemical=chems.Yeast,
                                    particle_size='Particulate', organic=True,
                                    degradability='Slowly')
    NH3 = Component('NH3', particle_size='Dissolved_gas', 
                    degradability='Undegradable', organic=False)
    CO2 = Component('CO2', particle_size='Dissolved_gas', 
                    degradability='Undegradable', organic=False)
    HAP = Component('HAP', search_ID='hydroxyapatite', 
                    particle_size='particulate', 
                    degradability='Undegradable', organic=False)
    CaCl2 = Component('CaCl2', particle_size='Soluble', 
                    degradability='Undegradable', organic=False)
    
    # common urine consitituents
    org_kwargs = dict(particle_size='Soluble', degradability='Readily', organic=True)
    urea = Component('Urea', **org_kwargs)
    creatinine = Component('Creatinine', **org_kwargs)
    Hhip = Component('Hippuric acid', **org_kwargs)
    Hcit = Component('Citric acid', **org_kwargs)
    Hglu = Component.from_chemical('Glucuronic_acid', chemical='C6H10O7', **org_kwargs)
    Huric = Component('Uric acid', **org_kwargs)
    other_COD = cmps_df.S_F.copy('Other_COD')
    
    ig_kwargs = dict(particle_size='Soluble', degradability='Undegradable', organic=False)
    chloride = Component('Cl-', **ig_kwargs)
    sodium = Component('Na+', **ig_kwargs)
    potassium = Component('K+', **ig_kwargs)
    IS = Component('SO4-2', measured_as='S', **ig_kwargs)
    IP = Component('PO4-3', measured_as='P', **ig_kwargs)
    
    ash = Component.from_chemical('Ash', chemical=chems.Ash, **ig_kwargs)
    other_SS = ash.copy('other_SS')
    
    cmps = Components([cmps_df.H2O, Yeast, NH3, CO2, HAP, CaCl2, 
                       urea, creatinine, Hhip, Hcit, Hglu, Huric, other_COD,
                       chloride, sodium, potassium, IS, IP, ash, other_SS])
    cmps.default_compile()
    if set_thermo: qs.set_thermo(cmps)
    return cmps

class HApFermenter(qs.SanUnit, BatchBioreactor):
    
    _units = BatchBioreactor._units
    _N_ins = 3      # [0] feed fresh urine, [1] osteoyeast, [2] CaCl2 
    _N_outs = 3     # [0] vent, [1] liquid effluent, [2] precipitates (yeast cells + HAP)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 tau=5, N=2, T=273.15+37, P=101325,
                 f_maximum_hap_yield=0.66, precipitate_moisture=90,
                 biomass_yield=0.1, inoculum_concentration=0.5):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self._init(tau=tau*24, N=N, T=T, P=P)
        self.f_maximum_hap_yield=f_maximum_hap_yield
        self.precipitate_moisture=precipitate_moisture
        self.biomass_yield=biomass_yield
        self.inoculum_concentration=inoculum_concentration
    
    @property
    def precipitate(self):
        return self.outs[2]
    
    @property
    def precipitate_moisture(self):
        '''[float] Moisture content [w%] of the precipitated product.'''
        return self._p_h2o
    @precipitate_moisture.setter
    def precipitate_moisture(self, p):
        if p >= 100 or p <= 0:
            raise ValueError(f'Moisture content must be between 0-100%, not {p}%')
        self._p_h2o = p
        
    @property
    def f_maximum_hap_yield(self):
        '''[float] Achieved fraction of maximum theoretical yield of hydroxyapatite, unitless.'''
        return self._f_max_y
    @f_maximum_hap_yield.setter
    def f_maximum_hap_yield(self, f):
        if f <= 0 or f > 1:
            raise ValueError(f'fraction of maximum yield must be in (0,1], not {f}')
        self._f_max_y = f
    
    @property
    def biomass_yield(self):
        '''[float] Overall yeast biomass yield in g yeast COD/g urine COD.'''
        return self._y_bio
    @biomass_yield.setter
    def biomass_yield(self, y):
        if y > 1:
            raise ValueError(f'overall yeast biomass yield must be less than 1, not {y}')
        self._y_bio = y
    
    @property
    def inoculum_concentration(self):
        '''[float] Overall concentration of osteoyeast inoculation, in kg/m3'''
        return self._cinocu
    @inoculum_concentration.setter
    def inoculum_concentration(self, c):
        self._cinocu = c
    
    def _setup(self):
        SanUnit._setup()
        vent, effluent, precip = self.outs
        vent.phase = 'g'
        vent.T = effluent.T = precip.T = self.T
        vent.P = effluent.P = precip.P = self.P
    
    def _run(self):
        urine, inocu, Ca = self.ins
        vent, eff, precip = self.outs
        f_hap = self.f_maximum_hap_yield
        y_bio = self.biomass_yield
        c_inocu = self.inoculum_concentration
        p_h2o = self.precipitate_moisture / 100
        Q = urine.F_vol   # m3/hr
        yeast_i_COD = inocu.components.Yeast.i_COD
        inocu.imass['Yeast'] = Qin_yeast = Q * c_inocu # kg/hr
        sol_IP = urine.composite('P', flow=True, particle_size='s', organic=False) # kg/hr
        sol_Ca = urine.composite('Ca', flow=True, particle_size='s', organic=False) # kg/hr
        COD = urine.composite('COD', flow=True)
        need_Ca = max(0, sol_IP * 2.156557 - sol_Ca)
        Ca.imass['CaCl2'] = need_Ca * 2.7693
        vent.empty()
        precip.imass['HAP'] = m_hap = sol_IP * f_hap / 0.184987
        precip.imass['Yeast'] = m_bio = max(Qin_yeast + COD * y_bio / yeast_i_COD, 0)
        precip.imass['H2O'] = m_h2o = (m_hap + m_bio) * p_h2o / (1-p_h2o)
        cmps = eff.components
        eff.copy_flow(urine)
        eff.mass[(cmps.i_P * cmps.s) > 0] *= 1-f_hap
        eff.imass['Yeast'] = eff.imass['Urea'] = 0
        eff.mass[cmps.i_COD > 0] *= 1-y_bio # conservative estimation of COD degradation
        eff.imol['NH3'] += urine.imol['Urea'] * 2  # ignore gaseous NH3 in vent
        eff.imass['H2O'] -= m_h2o
        
        
    