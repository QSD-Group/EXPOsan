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
from qsdsan.unit_operations import (
    Reactor, IsothermalCompressor, HXutility, HXprocess, MixTank,
    HydrothermalLiquefaction, Hydroprocessing, KnockOutDrum,
    )

__all__ = (
    'HydrothermalLiquefaction',
    'Hydroprocessing',
    'KnockOutDrum',
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
    COD_removal : float or dict
        Removal of influent COD.
    H2_yield : float
        H2 yield as in g H2/g COD removed.
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
    N_chamber : int
        Number of cell chambers.
    chamber_thickness : float
        Thickness of a single chamber, [m].
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
                 COD_removal=0.95, # assumed
                 H2_yield=0.157888654,
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
                 N_chamber=3,
                 chamber_thickness=0.02, # m
                 electrode_cost=40000, # $/m2
                 anion_exchange_membrane_cost=170, # $/m2
                 cation_exchange_membrane_cost=190, # $/m2
                 electrolyte_load=13.6, # kg/m3, 0.1 M of KH2PO4 (MW=136 k/mole)
                 electrolyte_price=30, # $/kg
                 annual_replacement_ratio=0, # Jiang assumed 2%, but 3% of maintenance already considered in TEA
                 include_PSA=True,
                 PSA_efficiency=0.95,
                 PSA_compressor_P=101325,
                 ):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=F_BM_default)
        self.COD_removal = COD_removal
        self.H2_yield = H2_yield
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
        self.N_chamber = N_chamber
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
        inf = self.ins[0]
        gas, H2, N, P, K, eff = self.outs
        eff.copy_like(inf)
        water_in = eff.imass['Water']
        eff.imass['Water'] = 0
               
        fert_IDs = self.N_IDs, self.P_IDs, self.K_IDs
        recoveries = self.N_recovery, self.P_recovery, self.K_recovery
        for IDs, out, recovery in zip(fert_IDs, (N, P, K), recoveries):
            out.imass[IDs] = eff.imass[IDs] * recovery
            eff.imass[IDs] -= out.imass[IDs]
        
        gas.empty()
        comp = self.gas_composition
        gas.imass[list(comp.keys())] = list(comp.values())
        cmps = self.components
        COD_removal = self.COD_removal
        COD_in = sum(inf.imass[i.ID]*i.i_COD for i in cmps)       
        H2_mass = COD_in * COD_removal * self.H2_yield
        scale_factor = H2_mass/gas.imass['H2']
        gas.F_mass *= scale_factor
        self._PSA_H2_lb_flowrate = gas.F_mass / _lb_to_kg
        gas.phase = 'g'
        
        for i in cmps:
            if i.ID not in ('Water', *fert_IDs):
                eff.imass[i.ID] *= (1-COD_removal)
        eff.imass['Water'] = water_in
        
        H2_tot = gas.imass['H2']
        H2.imass['H2'] = H2_recycled = H2_tot * self.PSA_efficiency
        gas.imass['H2'] = H2_tot - H2_recycled

        
    def _design(self):
        Design = self.design_results
        
        # 96485 is the Faraday constant C/mol (A·s/mol e)
        # MW of H2 is 2 g/mol, 2 electrons per mole of H2
        factor = 2/(2/1e3) * 96485 # (A·s/kg H2)
        H2_production = (self.outs[0].imass['H2']+self.outs[1].imass['H2']) / 3600 # kg/s
        current_eq = factor * H2_production # A
        area = current_eq / self.average_current_density
        Design['Area'] = area
        Design['Volume'] = volume = area * self.N_chamber * self.chamber_thickness
        try: hours = self.system.operating_hours
        except: hours = 365*24
        Design['Total Electrolyte'] = tot_ec = self.electrolyte_load * volume
        Design['Annual Electrolyte'] = annual_ec = tot_ec * self.annual_replacement_ratio
        self.ins[1].imass['Electrolyte'] = annual_ec / hours
        
        EO_power = self.EO_current_density * self.EO_voltage # W/m2, when online
        EO_electricity_per_area = EO_power/1e3 * self.EO_online_time_ratio # kWh/h/m2
        Design['EO electricity'] = EO_electricity = area * EO_electricity_per_area # kWh/h
        
        ED_power = self.ED_current_density * self.ED_voltage # W/m2, when online
        ED_electricity_per_area = ED_power/1e3 * self.ED_online_time_ratio # kWh/h/m2
        Design['ED electricity'] = ED_electricity = area * ED_electricity_per_area # kWh/h
        total_power = EO_electricity + ED_electricity
        self.power_utility.consumption = total_power
        
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
        Cost['Electrolyte'] = Design['Total Electrolyte']*self.electrolyte_price # initial capital cost
        cell_cost = Cost['Stack'] + Cost['Electrolyte']
        
        self._decorated_cost()
        
        # Cost is based on all replacement costs, mass just considering the electrolyte
        replacement = self.ins[1]
        annual_replacement_ratio = self.annual_replacement_ratio
        if annual_replacement_ratio:
            replacement.price = cell_cost*annual_replacement_ratio/Design['Annual Electrolyte']
        else:
            replacement.price = 0
        
        
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
    def normalized_CAPEX(self):
        '''Installed equipment cost per kg/hr of H2, [$/(kg H2/hr)].'''
        return self.installed_cost/(self.outs[0].imass['H2']+self.outs[1].imass['H2'])

    @property
    def normalized_energy_consumption(self):
        '''Electricity consumption per kg of H2, [kWh/kg H2].'''
        return self.power_utility.rate/(self.outs[0].imass['H2']+self.outs[1].imass['H2'])


class SAFElectrochemical(Electrochemical):
    '''To allow skipping unit simulation for different configurations.'''
    
    skip = False
    include_PSA_cost = True
    
    def _run(self):
        if self.skip:
            self.ins[1].empty() # replacement_surrogate
            for i in self.outs: i.empty()
            self.outs[-1].copy_like(self.ins[0])
        else: Electrochemical._run(self)

    def _design(self):
        if self.skip: self.design_results.clear()
        else: Electrochemical._design(self)

    def _cost(self):
        if self.skip: self.baseline_purchase_costs.clear()
        else: 
            Electrochemical._cost(self)
            if not self.include_PSA_cost:
                self.baseline_purchase_costs.pop('PSA')


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
            recycled.mix_from(recycled_streams)
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