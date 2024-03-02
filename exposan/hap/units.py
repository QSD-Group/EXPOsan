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
from biosteam import Stream, Facility
from biosteam.units import BatchBioreactor
from biosteam.units.decorators import cost
from exposan.hap import SimpleBlower, Locations, SimpleCVRP
# from math import ceil

__all__ = ('HApFermenter', 
           'YeastProduction', 
           'CollectionDistribution',
           'PrecipitateProcessing',)

#%%
@cost('Recirculation flow rate', 'Recirculation pumps', kW=30, S=77.22216,
      cost=47200, n=0.8, BM=2.3, CE=522, N='Total number of reactors')
@cost('Reactor volume', 'Cleaning in place', CE=801, cost=2000, S=1, n=0.6, 
      BM=1.8, N='Sites in parallel')  # assume ~1/2 reactor cost
@cost('Reactor volume', 'Agitators', kW=2.2, CE=801, cost=245, S=1, n=0.54, 
      lb=0.1, BM=1.5, N='Total number of reactors')
@cost('Reactor volume', 'Reactors', CE=801, cost=3936, S=1, n=0.54, 
      lb=0.1, BM=1.5, N='Total number of reactors')
@cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900,
      S=20920000.0, n=0.7, BM=2.2, N='Total number of reactors',
      magnitude=True) # Based on a similar heat exchanger
class HApFermenter(qs.SanUnit, BatchBioreactor):
    
    cost_items = {}
    _units = BatchBioreactor._units
    _N_ins = 3      # [0] feed fresh urine, [1] osteoyeast, [2] CaCl2 
    _N_outs = 3     # [0] vent, [1] liquid effluent, [2] precipitates (yeast cells + HAP)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 N_parallel_HApFermenter=10, tau=60, T=273.15+37, P=101325,
                 f_maximum_hap_yield=0.66, precipitate_moisture=90,
                 biomass_yield=0.1, inoculum_concentration=0.5,
                 CaCl2_price=0.3,
                 # labor_wage=21.68,  # assume 20% over San Francisco minimum wage by default
                 ):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self._init(tau=tau, N=2, T=T, P=P)
        self.N_parallel_HApFermenter = N_parallel_HApFermenter
        self.f_maximum_hap_yield=f_maximum_hap_yield
        self.precipitate_moisture=precipitate_moisture
        self.biomass_yield=biomass_yield
        self.inoculum_concentration=inoculum_concentration
        self.ins[2].price=CaCl2_price
        # self.labor_wage = labor_wage
    
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
        yeast_i_N = inocu.components.Yeast.i_N
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
        eff.imass['NH3'] -= (m_bio - Qin_yeast) * yeast_i_N
        eff.imass['H2O'] -= m_h2o
    
    def _design(self):
        super()._design()
        D = self.design_results
        D['Sites in parallel'] = self.N_parallel_HApFermenter
        D['Total number of reactors'] = D['Number of reactors'] * self.N_parallel_HApFermenter
        

#%%
@cost('Recirculation flow rate', 'Recirculation pumps', kW=30, S=77.22216,
      cost=47200, n=0.8, BM=2.3, CE=522, N='Number of reactors')
@cost('Reactor volume', 'Cleaning in place', CE=801, cost=2000, S=1, n=0.6, BM=1.8)  # assume ~1/2 reactor cost
@cost('Reactor volume', 'Agitators', kW=2.2, CE=801, cost=245, S=1, n=0.54, 
      lb=0.1, BM=1.5, N='Number of reactors')
@cost('Reactor volume', 'Reactors', CE=801, cost=3936, S=1, n=0.54, 
      lb=0.1, BM=1.5, N='Number of reactors')
@cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900,
      S=20920000.0, n=0.7, BM=2.2, N='Number of reactors',
      magnitude=True) # Based on a similar heat exchanger
class YeastProduction(Facility, qs.SanUnit, BatchBioreactor):
    
    '''
    Fermenter for the production of S. Boulardii active yeast, process design
    follows [1]_.
    
    '''
    network_priority = 0
    cost_items = {}
    _units = {
        **BatchBioreactor._units,
        'Effluent flowrate': 'gal/hr',
        'Aeration duty': 'm3/hr',
              }
    _N_ins = 3      # [0] carbon source (e.g., molasses + water), [1] nutrients (N & P), [2] minerals, vitamins etc. + seed
    _N_outs = 2     # [0] vent, [1] fermentation broth
    V_wf = 0.7
    auxiliary_unit_names = ('blower', )
    
    _centrifuge_cost_factor = 1.0
    
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
    
    def __init__(self, ID='', N_parallel_HApFermenter=10,
                 tau=12, T=273.15+35, P=101325,
                 yield_on_sugar=1.45, 
                 sugar_concentration=275.95,
                 N_to_sugar_ratio=0.01245,
                 P_to_sugar_ratio=0.00522,
                 minerals_vitamins_to_sugar_ratio=0.005349,
                 aeration_duty=60, sugar_price=None,
                 nutrient_price=None, mineral_vitamin_price=None,
                 reactor_height_to_diameter_ratio=2,
                 labor_wage=21.68,  # assume 20% over San Francisco minimum wage by default
                 ):

        SanUnit.__init__(self, ID, ins=None, outs=['vent', 'yeast'], thermo=None, init_with='WasteStream')
        self._system = None
        self._other_units = None
        self._init(tau=tau, N=2, T=T, P=P)
        self.N_parallel_HApFermenter = N_parallel_HApFermenter
        self.blower = SimpleBlower(self.ID+'_blower', 
                                   ins=Stream(self.ID+'blower_in', phase='g'),
                                   outs=Stream(self.ID+'blower_out', phase='g'))
        self.yield_on_sugar = yield_on_sugar                    # kg yeast / kg sugar
        self.sugar_concentration = sugar_concentration          # kg / m3 fermentation broth
        self.N_to_sugar_ratio = N_to_sugar_ratio
        self.P_to_sugar_ratio = P_to_sugar_ratio
        self.minerals_vitamins_to_sugar_ratio = minerals_vitamins_to_sugar_ratio
        self.aeration_duty = aeration_duty                      # m3 sanitized air / m3 broth / hr
        self.sugar_price = sugar_price
        self.nutrient_price = nutrient_price
        self.mineral_vitamin_price = mineral_vitamin_price
        self.reactor_height_to_diameter_ratio = reactor_height_to_diameter_ratio
        self.labor_wage = labor_wage

    @property
    def N_parallel_HApFermenter(self):
        return self._np
    @N_parallel_HApFermenter.setter
    def N_parallel_HApFermenter(self, np):
        self._np = int(np)

    @property
    def design_production_rate(self):
        '''[float] Design rate of production of fresh yeast biomass, in kg/hr.'''
        isa = isinstance
        units = [u for u in self._system.path if isa(u, HApFermenter)]
        r_yeast = sum(u.ins[1].imass['Yeast']/0.8 for u in units) * self.N_parallel_HApFermenter
        return r_yeast

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
    
    @property
    def labor_wage(self):
        '''[float] Hourly labor wage for fermenter operation and maintenance, in USD/hr.'''
        return self._wage
    @labor_wage.setter
    def labor_wage(self, wage):
        self._wage = wage
    
    def _run(self):
        sugars, nutrients, seed = self.ins
        vent, effluent = self.outs
        y_yeast = self.yield_on_sugar
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
        effluent.imass['Yeast'] = yeast = (sugars.F_mass - sugars.imass['H2O']) * y_yeast
        effluent.imass['Ethanol'] = ethanol = sugars.F_mass * 3.5e-4
        effluent.imass['H2O'] = (nutrients.F_mass + seed.F_mass + sugars.F_mass) - ethanol - yeast

    def _design(self):
        super()._design()
        D = self.design_results
        V_i =  D['Reactor volume']
        D['Effluent flowrate'] = self.effluent.get_total_flow('gal/hr')
        N = D['Number of reactors']
        D['Aeration duty'] = Q_air = V_i * N * self.V_wf * self.aeration_duty
        r_dim = self.reactor_height_to_diameter_ratio
        d_diffuser = r_dim * self.V_wf * (4*V_i/np.pi/r_dim) ** (1/3)
        blower = self.blower
        blower.P = max(d_diffuser * 9804.14 * 1.05, 10132.5) + 101325
        air, = blower.ins
        air.imol['O2'] = 21
        air.imol['N2'] = 78
        air.set_total_flow(Q_air, 'm3/hr')
        blower.simulate()
    
    def _cost(self):
        super()._cost()
        D = self.design_results
        C = self.baseline_purchase_costs
        q = D['Effluent flowrate']
        C['Centrifuge'] = (894 + 211*q) * self._centrifuge_cost_factor
        self.F_BM['Centrifuge'] = 1.5
        self.power_utility.consumption += q*0.035
        sugars, nutrients, seed = self.ins
        sugars.price = self.sugar_price * (1 - sugars.imass['H2O']/sugars.F_mass)
        nutrients.price = self.nutrient_price
        seed.price = self.mineral_vitamin_price * (1 - seed.imass['Yeast']/seed.F_mass)
        opex = self.add_OPEX = {}
        opex['NaCl'] = self.effluent.F_vol * 0.762 * self._prices['NaCl']  # final concentration of NaCl in kg/m3 fermentation broth
        opex['Operator'] = self.labor_wage * 8/24

#%%
class CollectionDistribution(Facility, SanUnit):
    
    network_priority = -1
    _N_ins = 0     
    _N_outs = 0
    
    _units = {
        'Collection interval': 'hr',
        'Total travel distance': 'mile/trip'
        }
    
    _default_area = {
        'x': (0, 13.36e3),  # in meter
        'y': (0, 11.31e3)
        }
    
    def __init__(self, ID='', N_parallel_HApFermenter=10, 
                 locations=None, vehicle_vol_capacity=20, 
                 route_max_duration=8*60, duration_per_location=20,
                 vehicle_fuel_cost=0.75, # USD/mile
                 vehicle_rental_price=200, # USD/day
                 labor_wage=21.68,  # assume 20% over San Francisco minimum wage by default
                 ):
        super().__init__(ID, ins=None, outs=(), thermo=None)
        self.cvr = SimpleCVRP()
        self.N_parallel_HApFermenter = N_parallel_HApFermenter
        if locations is not None: self.locations = locations
        self.vehicle_vol_capacity = vehicle_vol_capacity
        self.route_max_duration = route_max_duration
        self.duration_per_location = duration_per_location
        self.vehicle_fuel_cost = vehicle_fuel_cost
        self.vehicle_rental_price = vehicle_rental_price
        self.labor_wage=labor_wage

    @property
    def N_parallel_HApFermenter(self):
        return self._np
    @N_parallel_HApFermenter.setter
    def N_parallel_HApFermenter(self, np):
        self._np = int(np)
        self.locations = None

    @property
    def locations(self):
        return self._locs
    @locations.setter
    def locations(self, locs):
        if locs is None:
            area = self._default_area
            locs = Locations.random_within_area(
                self._np+1, distance_metric='cityblock',
                x_range=area['x'], # in meter
                y_range=area['y'],
                demands=1,
                )
        else:
            assert isinstance(locs, Locations)
            self._np = locs.size - 1
        self._locs = self.cvr.locations = locs

    @property
    def vehicle_vol_capacity(self) -> int:
        '''[int] Vehicle volumetric capacity, in m3.'''
        return self._vol_max
    @vehicle_vol_capacity.setter
    def vehicle_vol_capacity(self, vol):
        self._capacity = None
        self._vol_max = int(vol)
    
    @property
    def route_max_duration(self) -> int:
        return self._t_max
    @route_max_duration.setter
    def route_max_duration(self, t):
        '''[int] Maximum duration allocated for each route, in min.'''
        self._capacity = None
        self._t_max = int(t)
        
    @property
    def duration_per_location(self):
        return self._ti
    @duration_per_location.setter
    def duration_per_location(self, ti):
        '''[int] Amount of time needed to collect from and distribute to each 
        location on average, in min.'''
        self._capacity = None
        self._ti = int(ti)
    
    @property
    def vol_per_location(self):
        if self._system is None: return
        else:
            isa = isinstance
            for u in self._system.path:
                if isa(u, HApFermenter):
                    vol = u.design_results['Batch time'] * u.outs[2].F_vol # m3
                    return vol
            raise RuntimeError('no HApFermenter in system path')                    
    
    @property
    def capacity(self):
        if self._capacity is None:
            cap_t = int(self.route_max_duration / self.duration_per_location)
            cap_vol = int(self.vehicle_vol_capacity / self.vol_per_location)
            self._capacity = min(cap_t, cap_vol)
        return self._capacity
    
    def _run(self):
        if self.cvr.vehicle_capacity != self.capacity:
            self.cvr.vehicle_capacity = self.capacity
        self.cvr.register()
        self.cvr.solve(100, True)

    def _design(self):
        D = self.design_results
        locs = self.locations
        isa = isinstance
        for u in self._system.path:
            if isa(u, HApFermenter):
                tau = u.design_results['Batch time'] # hr
                break
        D['Collection interval'] = tau
        D['Total travel distance'] = dis = sum(locs.path_cost(r) for r in self.cvr.routes) * 6.2137e-4 # mile
        D['Number of routes'] = n_routes = len(self.cvr.routes)
        opex = self.add_OPEX = {}
        opex['Fuel'] = self.vehicle_fuel_cost * dis / tau
        opex['Vehicle rental'] = self.vehicle_rental_price * n_routes / tau
        opex['Labor'] = self.labor_wage * n_routes * self.route_max_duration/60 / tau
        

#%%
class PrecipitateProcessing(Facility, SanUnit):
    
    '''Simple precipitate post processing facility involving sizing and costing 
    of storage tank, recessed plate filter press dryer, incinerator, and relevant equipment.'''
    
    network_priority = 1
    _units = {
        'Dry solid loading rate': 'lb/hr',
        'Dryer cake volume': 'ft3/hr',
        'Total dryer chamber volume': 'ft^3',
        'Dryer power': 'kW',
        'Incinerator burn rate': 'lb/hr',
        'Incinerator power': 'kW',
        'Auxiliary fuel demand': 'gal/hr',
              }
    _N_ins = 1      
    _N_outs = 2     # [0] product, [1] wastewater

    _dryer_cost_factor = 1.0
    _incinerator_cost_factor = 1.0
    _F_BM_default = {
        'Recessed plate filter press': 3,
        'Incinerator': 1.7,
        }

    def __init__(self, ID='', N_parallel_HApFermenter=10, dryer_operating_hours=8,
                 dryer_cycle_time=None, dryer_cake_moisture=60, dryer_cake_density=71,
                 auxiliary_fuel_HV=137381, auxiliary_fuel_price=5.23,
                 labor_wage=21.68,  # assume 20% over San Francisco minimum wage by default
                 ):
        
        super().__init__(ID, ins=None, outs=['product', 'ww'], thermo=None)
        self.N_parallel_HApFermenter = N_parallel_HApFermenter
        self.dryer_operating_hours = dryer_operating_hours
        self.dryer_cycle_time = dryer_cycle_time
        self.dryer_cake_moisture = dryer_cake_moisture
        self.dryer_cake_density = dryer_cake_density
        self.auxiliary_fuel_HV = auxiliary_fuel_HV
        self.auxiliary_fuel_price = auxiliary_fuel_price
        self.labor_wage = labor_wage
        
    @property
    def N_parallel_HApFermenter(self):
        return self._np
    @N_parallel_HApFermenter.setter
    def N_parallel_HApFermenter(self, np):
        self._np = int(np)
    
    @property
    def dryer_operating_hours(self):
        '''[float] hours per day recessed plate filter press is operated, in hr.'''
        return self._dryer_hpd
    @dryer_operating_hours.setter
    def dryer_operating_hours(self, hpd):
        if hpd <= 0 or hpd > 24:
            raise ValueError('invalid value of operating hours per day')
        self._dryer_hpd = hpd
    
    @property
    def dryer_cycle_time(self):
        '''[float] recessed plate filter cycle time, in hr.'''
        if self._fct is None:
            feed, = self.ins
            ss = sum(feed.mass * feed.chemicals.x) / feed.F_mass * 100
            return 2.9 - 0.225 * ss + 0.0125 * ss**2
        else:
            return self._fct
    @dryer_cycle_time.setter
    def dryer_cycle_time(self, fct):
        self._fct = fct
    
    @property
    def dryer_cake_moisture(self):
        '''[float] recessed plate filter press cake moisture content, in %.'''
        return self._cm
    @dryer_cake_moisture.setter
    def dryer_cake_moisture(self, cm):
        self._cm = cm
        
    @property
    def dryer_cake_density(self):
        '''[float] filter cake density, in lb/ft^3.'''
        return self._crho
    @dryer_cake_density.setter
    def dryer_cake_density(self, rho):
        self._crho = rho
        
    @property
    def auxiliary_fuel_HV(self):
        '''[float] heating value of auxiliary fuel, in Btu/gal.'''
        return self._fuel_HV
    @auxiliary_fuel_HV.setter
    def auxiliary_fuel_HV(self, hv):
        self._fuel_HV = hv
    
    @property
    def auxiliary_fuel_price(self):
        '''[float] price of auxiliary fuel, in USD/gal.'''
        return self._fp
    @auxiliary_fuel_price.setter
    def auxiliary_fuel_price(self, p):
        self._fp = p
        
    @property
    def labor_wage(self):
        '''[float] Hourly labor wage for facility operation and maintenance, in USD/hr.'''
        return self._wage
    @labor_wage.setter
    def labor_wage(self, wage):
        self._wage = wage
    
    def _run(self):
        cake, ww = self.outs
        feed, = self.ins
        n = self.N_parallel_HApFermenter
        isa = isinstance
        precipitates = [u.outs[2] for u in self._system.path if isa(u, HApFermenter)]
        feed.mix_from(precipitates)
        feed.scale(n)
        ww.copy_like(feed)
        cake.copy_like(feed)
        cm = self.dryer_cake_moisture
        cake.imass['H2O'] = (cake.imass['Yeast'] + cake.imass['HAP'])/(100-cm)*cm
        ww.separate_out(cake)
        product = cake
        yeast = cake.chemicals.Yeast
        product.imass['Ash'] = cake.imass['Yeast'] * (1-yeast.f_Vmass_Totmass)
        product.imass['Yeast'] = product.imass['H2O'] = 0.
        
    def _design(self):
        feed, = self.ins
        product, ww = self.outs
        cmps = feed.chemicals
        D = self.design_results
        D['Dry solid loading rate'] = LR = sum(feed.mass * cmps.x) * 2.204623
        fct = self.dryer_cycle_time
        hpd = self.dryer_operating_hours
        cm = self.dryer_cake_moisture
        csg = self.dryer_cake_density
        D['Dryer cake volume'] = CV = LR * 100 / (100-cm) / csg
        D['Total dryer chamber volume'] = TCV = CV * 24 / fct / hpd
        D['Dryer power'] = dkW = TCV / 10 * 19.863
        D['Incinerator burn rate'] = BR = LR / (1-cm/100)
        D['Incinerator power'] = ikW = BR / 87.5 * 3.5
        HIR = 10**(3.247 + 0.0126*cm) * LR
        D['Auxiliary fuel demand'] = (HIR - feed.HHV*0.94782)/self.auxiliary_fuel_HV
        self.power_utility.consumption = dkW * hpd/24 + ikW
    
    def _cost(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        fbm_dry, fbm_inc = self.F_BM.values()
        TCV = D['Total dryer chamber volume']
        BR = D['Incinerator burn rate']
        C['Recessed plate filter press'] = capex_dry = \
            (TCV * 829.7 + 4227.7)*self._dryer_cost_factor
        C['Incinerator'] = capex_inc = 894.4*BR**0.8 *self._incinerator_cost_factor
        C['Storage tank'] = (capex_dry + capex_inc) * 0.1
        opex = self.add_OPEX = {}
        opex['Dryer parts & maintenance'] = 0.0013 * fbm_dry * capex_dry /365/24
        opex['Incinerater parts & maintenance'] = 0.0045 * fbm_inc * capex_inc /365/24
        opex['Labor'] = 3 * self.labor_wage * self.dryer_cycle_time / 24
        opex['Auxiliary fuel'] = D['Auxiliary fuel demand'] * self.auxiliary_fuel_price
        