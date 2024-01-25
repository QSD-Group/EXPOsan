# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References
----------
.. [1] Zheng, G., Yu, X., Li, Z., Yu, M., Yao, J., Chi, B., & Wang, J. (2013). 
    Boulardii active dry yeasts and production method thereof 
    (Patent No. CN103374531A).

'''
import qsdsan as qs, numpy as np
from qsdsan import SanUnit
from biosteam import Stream
from biosteam.units import BatchBioreactor
from exposan.hap import SimpleBlower

__all__ = ('HApFermenter', 'SBoulardiiFermenter', )

#%%
class HApFermenter(qs.SanUnit, BatchBioreactor):
    
    _units = BatchBioreactor._units
    _N_ins = 3      # [0] feed fresh urine, [1] osteoyeast, [2] CaCl2 
    _N_outs = 3     # [0] vent, [1] liquid effluent, [2] precipitates (yeast cells + HAP)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 tau=120, N=2, T=273.15+37, P=101325,
                 f_maximum_hap_yield=0.66, precipitate_moisture=90,
                 biomass_yield=0.1, inoculum_concentration=0.5,
                 CaCl2_price=0.3):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self._init(tau=tau, N=N, T=T, P=P)
        self.f_maximum_hap_yield=f_maximum_hap_yield
        self.precipitate_moisture=precipitate_moisture
        self.biomass_yield=biomass_yield
        self.inoculum_concentration=inoculum_concentration
        self.ins[2].price=CaCl2_price
    
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
        SanUnit._setup(self)
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
        eff.mass[cmps.i_COD > 0] *= 1-y_bio        # conservative estimation of COD degradation
        eff.imol['NH3'] += urine.imol['Urea'] * 2  # ignore gaseous NH3 in vent
        eff.imass['H2O'] -= m_h2o

#%%

class SBoulardiiFermenter(qs.SanUnit, BatchBioreactor):
    
    '''
    Fermenter for the production of S. Boulardii active yeast, process design
    follows [1]_.
    
    '''
    
    _units = {
        **BatchBioreactor._units,
        'Aeration duty': 'm3/hr',
              }
    _N_ins = 3      # [0] carbon source (e.g., molasses + water), [1] nutrients (N & P), [2] minerals, vitamins etc. + seed
    _N_outs = 2     # [0] vent, [1] fermentation broth
    V_wf = 0.7
    auxiliary_unit_names = ('blower',)
    
    # USD/kg
    _prices = {
        'molasses': 2.5,
        'glucose': 3.5,
        'ammonium sulfate': 0.18,
        'H3PO4': 1.7,
        'MgSO4': 0.35,
        'VB1': 52.0,
        'VB2': 35.0,
        'VB5': 35.0,
        'VB6': 30.0,
        'VB7': 300,
        'NaCl': 0.1,
        }
    
    _micronutrient_mixture = (
        150,  # MgSO4, mass-based
        370,  # VB1
        800,  # VB2
        2300, # VB5
        360,  # VB6
        1200  # VB7
        )
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 tau=12, N=2, T=273.15+35, P=101325,
                 yield_on_sugar=1.45, design_production_rate=80, 
                 sugar_concentration=275.95,
                 N_to_sugar_ratio=0.01245,
                 P_to_sugar_ratio=0.00522,
                 minerals_vitamins_to_sugar_ratio=0.005349,
                 aeration_duty=60, sugar_price=None,
                 nutrient_price=None, mineral_vitamin_price=None,
                 reactor_height_to_diameter_ratio=2,
                 ):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self._init(tau=tau, N=N, T=T, P=P)
        self.blower = SimpleBlower(self.ID+'_blower', 
                                   ins=Stream(self.ID+'blower_in', phase='g'),
                                   outs=Stream(self.ID+'blower_out', phase='g'))
        self.yield_on_sugar = yield_on_sugar                    # kg yeast / kg sugar
        self.design_production_rate = design_production_rate    # kg yeast / hr
        self.sugar_concentration = sugar_concentration          # kg / m3 fermentation broth
        self.N_to_sugar_ratio = N_to_sugar_ratio
        self.P_to_sugar_ratio = P_to_sugar_ratio
        self.minerals_vitamins_to_sugar_ratio = minerals_vitamins_to_sugar_ratio
        self.aeration_duty = aeration_duty                      # m3 sanitized air / m3 broth / hr
        self.sugar_price = sugar_price
        self.nutrient_price = nutrient_price
        self.mineral_vitamin_price = mineral_vitamin_price
        self.reactor_height_to_diameter_ratio = reactor_height_to_diameter_ratio
        
    @property
    def design_production_rate(self):
        '''[float] Design rate of production of fresh yeast biomass, in kg/hr.'''
        return self._r_yeast
    @design_production_rate.setter
    def design_production_rate(self, r_yeast):
        self._r_yeast = r_yeast
        
    @property
    def yield_on_sugar(self):
        '''[float] Yeast yield in kg fresh yeast biomass / kg sugar fed.'''
        return self._y_yeast
    @yield_on_sugar.setter
    def yield_on_sugar(self, y_yeast):
        self._y_yeast = y_yeast
    
    @property
    def sugar_concentration(self):
        '''[float] Overall concentration of sugar, in kg / m3 fermentation broth.'''
        return self._c_sugar
    @sugar_concentration.setter
    def sugar_concentration(self, c_sugar):
        self._c_sugar = c_sugar    
    
    @property
    def N_to_sugar_ratio(self):
        '''[float] Nitrogen-to-sugar ratio, in kg-N/kg sugar fed.'''
        return self._n2s
    @N_to_sugar_ratio.setter
    def N_to_sugar_ratio(self, n2s):
        self._n2s = n2s

    @property
    def P_to_sugar_ratio(self):
        '''[float] Phosphorus-to-sugarratio, in kg-P/kg sugar fed.'''
        return self._p2s
    @P_to_sugar_ratio.setter
    def P_to_sugar_ratio(self, p2s):
        self._p2s = p2s
        
    @property
    def minerals_vitamins_to_sugar_ratio(self):
        '''[float] Weight ratio of added minerals and vitamins relative to sugar,
        in kg mixture / kg sugar fed.'''
        return self._mv2s
    @minerals_vitamins_to_sugar_ratio.setter
    def minerals_vitamins_to_sugar_ratio(self, mv2s):
        self._mv2s = mv2s
        
    @property
    def aeration_duty(self):
        '''[float] In m3 sanitized air / m3 fermentation broth / hr.'''
        return self._aeration_duty
    @aeration_duty.setter
    def aeration_duty(self, duty):
        self._aeration_duty = duty
    
    @property
    def sugar_price(self):
        '''[float] in USD / kg sugar fed'''
        if self._sugar_price is None:
            pm = self._prices['molasses']
            pg = self._prices['glucose']
            sugars = self.ins[0]
            mm, mg = sugars.imass['Molasses', 'Glucose']
            if mm + mg == 0: 
                return pm*2/3 + pg*1/3
            else:
                return (pm*mm + pg*mg) / (mm + mg)
        return self._sugar_price
    @sugar_price.setter
    def sugar_price(self, p):
        self._sugar_price = p
    
    @property
    def nutrient_price(self):
        '''[float] in USD / kg nutrient fed'''
        if self._nutrient_price is None:
            pn = self._prices['ammonium sulfate']
            pp = self._prices['H3PO4']
            nutrients = self.ins[1]
            mn, mp = nutrients.imass['Ammonium_sulfate', 'H3PO4']
            if mn + mp == 0: 
                return pn*17/21 + pp*4/21
            else:
                return (pn*mn + pp*mp) / (mn + mp)
        return self._nutrient_price
    @nutrient_price.setter
    def nutrient_price(self, p):
        self._nutrient_price = p
    
    @property
    def mineral_vitamin_price(self):
        '''[float] in USD / kg mineral-vitamin mixture fed'''
        if self._nutrient_price is None:
            cmps = ['MgSO4', 'VB1', 'VB2', 'VB5', 'VB6', 'VB7']
            prices = np.array([self._prices[k] for k in cmps])
            seed = self.ins[2]
            ms = seed.imass[cmps]
            if sum(ms) == 0:
                composition = np.array(self._micronutrient_mixture)
                return sum(prices * composition)/sum(composition)
            else:
                return sum(prices * ms)/sum(ms)
        return self._mineral_vitamin_price
    @mineral_vitamin_price.setter
    def mineral_vitamin_price(self, p):
        self._mineral_vitamin_price = p

    @property
    def reactor_height_to_diameter_ratio(self):
        '''[float] Height-to-diameter ratio of an individual cylindrical reactor.'''
        return self._reactor_height_to_diameter_ratio
    @reactor_height_to_diameter_ratio.setter
    def reactor_height_to_diameter_ratio(self, ratio):
        self._reactor_height_to_diameter_ratio = ratio
        
    def _setup(self):
        super()._setup()

    
    def _run(self):
        sugars, nutrients, seed = self.ins
        y_yeast = self.yield_on_sugar
        if sugars.isempty():
            r_yeast = self.design_production_rate
            c_sugar = self.sugar_concentration # kg/m3
            n2s = self.N_to_sugar_ratio
            p2s = self.P_to_sugar_ratio
            mv2s = self.minerals_vitamins_to_sugar_ratio
            _micro_nutrients = np.array(self._micronutrient_mixture)
            Fm_molasses = r_yeast / y_yeast  # kg/hr
            Q = Fm_molasses / c_sugar  # m3/hr
            sugars.set_flow_by_concentration(Q, {'Molasses': c_sugar}, 
                                             units=('m3/hr', 'kg/m3'))
            nutrients.imass['Ammonium_sulfate'] = Fm_molasses * n2s / 0.212
            nutrients.imass['H3PO4'] = Fm_molasses * p2s / 0.316
            seed.imass['Yeast'] = r_yeast * 0.04
            seed.imass['MgSO4', 'VB1', 'VB2', 'VB5', 'VB6', 'VB7'] = \
                _micro_nutrients / sum(_micro_nutrients) * Fm_molasses * mv2s
        vent, effluent = self.outs
        effluent.imass['Yeast'] = yeast = (sugars.F_mass - sugars.imass['H2O']) * y_yeast
        effluent.imass['Ethanol'] = ethanol = sugars.F_mass * 3.5e-4
        effluent.imass['H2O'] = (nutrients.F_mass + seed.F_mass + sugars.F_mass) - ethanol - yeast

    def _design(self):
        super()._design()
        D = self.design_results
        V_i =  D['Reactor volume']
        N = D['Number of reactors']
        D['Aeration duty'] = Q_air = V_i * N * self.V_wf * self.aeration_duty
        r_dim = self.reactor_height_to_diameter_ratio
        d_diffuser = r_dim * self.V_wf * (4*V_i/np.pi/r_dim) ** (1/3)
        blower = self.blower
        blower.P = d_diffuser * 9804.14 * 1.05 + 101325
        air, = blower.ins
        air.imol['O2'] = 21
        air.imol['N2'] = 78
        air.set_total_flow(Q_air, 'm3/hr')
        blower.simulate()
    
    def _cost(self):
        super()._cost()
        sugars, nutrients, seed = self.ins
        sugars.price = self.sugar_price * (1 - sugars.imass['H2O']/sugars.F_mass)
        nutrients.price = self.nutrient_price
        seed.price = self.mineral_vitamin_price * (1 - seed.imass['Yeast']/seed.F_mass)
        self.add_OPEX['NaCl'] = self.effluent.F_vol * 0.762 * self._prices['NaCl']  # final concentration of NaCl in kg/m3 fermentation broth
    