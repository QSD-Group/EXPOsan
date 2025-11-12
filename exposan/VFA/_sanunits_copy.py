#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:27:03 2025

@author: blues
"""

from qsdsan import SanUnit, Construction, WasteStream, Stream, Chemical, Component, \
Components, set_thermo as qs_set_thermo
from qsdsan.utils import ospath, data_path, load_data, price_ratio
# from exposan.g2rt._sanunits import G2RTSolidsSeparation
import biosteam as bst
import qsdsan as qs
from biosteam import Splitter
from biosteam.units.decorators import cost
from qsdsan import CEPCI_by_year
from qsdsan.sanunits import SludgePump
import numpy as np, flexsolve as flx
from warnings import warn
from biosteam.exceptions import lb_warning
from math import ceil, pi, sqrt
from qsdsan.sanunits import Pump

__all__ = ('SolidsSeparation', 'SludgeThickening', 'BaseDosing', 'RedoxED',)


#%% Pretreatment: solids separation using centrifugation
centrifuge_path = ospath.join(data_path, 'sanunit_data/VFA/_centrifuge.csv')

hazen_william_C = 110 # dimensionless
pipe_length = 150 # ft
min_velocity = 2.5 # ft, based on minimum flowrate suggestion from manual of practice No. 8 6th ed page 173
gas_pressure = 10/12 # ft, positive gas pressure in anaerobic digestor

@cost(basis='Solids loading', ID='SolidsSeparation', units='kg/d',
      cost=215000, S=750*0.453592, N = 'Number of centrifuges',
      CE=CEPCI_by_year[2000], n=0.6, BM=650000/215000)
    # the cost is for one centrifuge unit, baseline cost and solids loading is from EPA fact sheet 2000
    #!!! check what exponential scaling factor to use? Normal values are 2.1 or 2.3, but the EPA design manual provided 3.023
    #!!! exponent n: 0.6 is considered default in chemical engineering and 0.3 is considered stronge economy of scale
class SolidsSeparation(SanUnit):
    '''
    Function of unit:
    A solids separation class for simulation of a solid bowl centrifuge 
    to separate solids from mixed stream.
    Calculation:
    Separation split is calculated based on moisture content 
    in the sludge and solids separation rate, both defined by user.
    Assume solubles and water share the same split.
    If the moiture content in the feed is smaller than the targeted moiture
    moisture content in sludge, the target moiture content will be ignored.
    Material in cost and design:
    References:
    Reference units are G2RTSolidsSeparation developed by Zixuan Wang, 
    the SolidCetrifuge developed by Yoel, and Thickener by Yalin.
    
    Parameters
    ----------
    Ins: Iterable (stream)
        Incoming waste stream.
    Outs: Iterable (stream)
        outs[0] is supernatent effluent
        outs[1] is sludge
    sludge_moisture : float
        Moisture content of the sludge, [wt% water].
    removal_rate: float
        Percent removal of solids from the mixed influent.
    solids : Iterable(stream)
        IDs of the solid components.
    solubles: Iterable(stream)
        IDs of the soluble components, NOT including water.
    
    '''
# !!! This is the first unit of the system, the in flow would be AD effluent   
    _N_ins = 1
    _N_outs = 2
    _units = {'Flow rate': 'm3/hr'}
    flow_rate_range = (0.95, 3.42) 
    # solids loading rate for solid boal centrifuge determined from reported cases in the EPA design manual
    m3_per_gal = 0.00378541
    min_per_hr = 60
    auxiliary_unit_names = ('feed_pump', 'dosing_pump', 'centrate_pump')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 removal_rate = 0.82, sludge_moisture = 0.926, polymer_dose = 0.0054, 
                 disposal_cost = 125/907.18474, kWhr_per_m3_per_s = 85.74, 
                 operating_hour = 5127.86, labor_hour = 1.088, auger_cost = 500,
                 conveyor_belt_cost = 229.67, sludge_density = 1.03*998.2, 
                 feed_flowrate = 1.50, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.tss_removal = removal_rate
        self.sludge_moisture = sludge_moisture
        self.polymer_dose = polymer_dose
        self.disposal_cost = disposal_cost
        self.kWhr_per_m3_per_s = kWhr_per_m3_per_s # power consumption per gpm of sludge into the centrifuge
        self.operating_hour = operating_hour # this is number of operations per year
        self.labor_hour = labor_hour
        self.auger_cost = auger_cost
        self.conveyor_belt_cost = conveyor_belt_cost
        self.sludge_density = sludge_density
        self.feed_flowrate = feed_flowrate
        # sludge density is estimated by specific gravity of primary slidge + WAS and water density at 20 C from MtCalf&Eddy; density of AD effluent was not foudn in the book
# !!! Disposal cost is a place holder value taken from Yalin's sludge thickening unit, check recent value
# !!! Find recent polymer costs or adjust inflation for EPA design manual numbers
# !!! Check the default value of pump_pressure
# !!! Check what F_BM_default = 1 means, deleted scale up, ppl, and estreme arguments from parent class

# !!! in the line below, what components is it referring to? Should import from _components.py?        
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids])

    def _init_lca(self):
        self.construction = [Construction("stainless_steel", linked_unit=self,
                                          item = "StainlessSteel", 
                                          quantity_unit= "kg"),
                             Construction("polyethylene", linked_unit=self,
                                          item = "Polyethylene",
                                          quantity_unit= "kg"),
                             Construction("electric_motor", linked_unit=self,
                                          item = "ElectricMotor",
                                          quantity_unit= "kg"),
                             Construction("pump", linked_unit=self,
                                          item = "Pump",
                                          quantity_unit= "ea"),
                             ]        
        
    def _run(self):
        inf, polymer = self.ins # index or comma needed when there is only one stream in inlets
        liquid_stream, solid_stream = self.outs
# This following is defining the subset of components that later could be accessed through []
        solubles, solids = self.solubles, self.solids
        
        TL_in = inf.F_mass - inf.imass[solids].sum()
        mc_in = TL_in / inf.F_mass # including both water and solubles in the moisture content
        mc_out = self.sludge_moisture # this moisture data should assume the total mass of water and solubles
# !!! Below the logic check is for moisture content in the solids but the focus here is the liquid stream
        # if mc_in < mc_out*0.999:
        if mc_in < mc_out:
            mc_out = mc_in
        
        solid_stream.imass[solids] = inf.imass[solids] * self.tss_removal
        TS_out = solid_stream.imass[solids].sum() # total soilds in the solids stream
        TL_out = TS_out / (1 - mc_out) * mc_out # total solubles and water mass flowrate in the solids stream
        solid_stream.imass[solubles] = inf.imass[solubles] * TL_out / TL_in 
        solid_stream.imass['H2O'] = TL_out - solid_stream.imass[solubles].sum()
        # the above line assumes that the solubles (including water and soluble chemicals) partition together and by the same ratio
        #liquid_stream.imass['H2O'] = TS_out*(1-mc_out) * mc_out - liquid_stream.imass[solubles].sum()
        
        liquid_stream.imass[solids] = inf.imass[solids] * (1 - self.tss_removal)
        liquid_stream.imass[solubles] = inf.imass[solubles] - solid_stream.imass[solubles]
        liquid_stream.imass['H2O'] = inf.imass['H2O'] - solid_stream.imass['H2O']
        
        polymer.imass['Polymer'] = inf.imass[solids] * self.polymer_dose
        
        pipe_diameter = sqrt(4 * self.feed_flowrate * 35.3147/3600 / pi / min_velocity)
        friction_head = 3.02 * pipe_length * (min_velocity ** 1.85) * \
        (hazen_william_C ** (-1.85)) * (pipe_diameter ** (-1.17)) # Using equation ESI-7 in the SI of Brian Shoener's 2016 paper
        feed_pump_p = friction_head + gas_pressure + 10 # include 10ft water column for pressure head
        centrate_pump_p = friction_head
        feed_pump = self.auxiliary('feed_pump', cls = SludgePump, ins = self.ins[0].copy('feed_pump_in'), P = feed_pump_p)
        dosing_pump = self.auxiliary('dosing_pump', cls = Pump, ins = self.ins[1].copy('dose_pump_in'))
        centrate_pump = self.auxiliary('centrate_pump', cls = SludgePump, outs = self.outs[0].copy('centrate_pump_out'), P = centrate_pump_p)

        # !!! the pump_pressure argument is the pressure of the output stream, assume atmospheric pressure for the feeding
        # and the dosing pumps.
        #!!! estimate the pressure of the centrate pu,p considering 3m (10ft) of elevation lift and the frictional loss in the pipes
        # calculate using the energy balance equation in fluid mechanics

    def _design(self):
        self.design_results['Feeding rate'] = feeding_rate = self.F_vol_in
        lower_bound, upper_bound = self.flow_rate_range
        if feeding_rate < lower_bound:
            lb_warning(self, 'Feeding rate', feeding_rate, 'm3/hr', lower_bound)
        self.design_results['Number of centrifuges'] = ceil(feeding_rate/upper_bound) 
        self.power_utility(self.F_vol_in * self.kWhr_per_m3_hr + self.conveyor_power * 
                           self.operating_hour) # this is electricity consumption not including pumping
    
        
    def _cost(self):
        # capital cost for auger and conveyor belt for sludge removal
        ts_out = self.outs[-1].F_mass
        self.baseline_purchase_costs['Auger'] = self.auger_cost * (ts_out ** 0.6) * current_CE / CE #!!! scale with total solids and multiply with the number of auger?
        self.F_BM['Auger'] = 2.03
        self.baseline_purchase_costs['Conveyor belt'] = self.conveyor_belt_cost * current_CE / CE #!!! scale with total solids?
        self.F_BM['Conveyor belt'] = 2.03
        # !!! find the auger and conveyor belt costs, also find the corresponding CE valeus and bare module factors
        #OPEX = sludge disposal, material replacement, labor, electricity, and polymer dosage
        self.add_OPEX = {'Sludge disposal': ts_out * self.disposal_cost, 
                         'Labor': self.wages * self.operating_hours * self.labor_hour} # two operators per shift
        # !!! find costs in the add_OPEX dictionary
        for p in (self.feed_pump, self.dose_pump, self.centrate_pump): p.simulate()

#%% Pretreatment: solids separation using microfiltration       
class SludgeThickening(SanUnit, Splitter):
    '''
    A generic class for concentrating (i.e., thickening) of sludge
    from wastewater treatment processes based on
    `Shoener et al. <https://doi.org/10.1039/C5EE03715H>`_

    The 0th outs is the water-rich supernatant (effluent) and
    the 1st outs is the solid-rich sludge.

    Two pumps (one for the supernatant and one for sludge) are included.

    Separation split is determined by the moisture (i.e., water)
    content of the sludge, soluble components will have the same split as water,
    insolubles components will all go to the retentate.

    Note that if the moisture content of the incoming feeds are smaller than
    the target moisture content, the target moisture content will be ignored.

    The following components should be included in system thermo object for simulation:
    Water.

    Parameters
    ----------
    ins : Iterable(stream)
        Dilute sludge stream.
    outs : Iterable(stream)
        Water/bulk-liquid-rich stream, concentrated sludge.
    sludge_moisture : float
        Moisture content of the sludge, [wt% water].
    solids : Iterable(stream)
        IDs of the solid components.
        If not provided, will be set to the default `solids` attribute of the components.
    disposal_cost : float
        Disposal cost of the dewatered solids. [$/kg].

    References
    ----------
    [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.
    '''

    SKIPPED = False
    _graphics = Splitter._graphics
    _ins_size_is_fixed = False
    _N_outs = 2
    auxiliary_unit_names = ('effluent_pump', 'sludge_pump')

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 sludge_moisture=0.96, solids=(),
                 disposal_cost=125/907.18474): # converting from $/U.S. ton
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with=init_with)
        self.sludge_moisture = sludge_moisture
        cmps = self.components
        self.solids = solids or tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids])
        self.disposal_cost = disposal_cost
        ID = self.ID
        eff = self.outs[0].proxy(f'{ID}_eff')
        sludge = self.outs[1].proxy(f'{ID}_sludge')
        # Add '.' in ID for auxiliary units
        self.effluent_pump = SludgePump(f'.{ID}_eff_pump', ins=eff, init_with=init_with)
        self.sludge_pump = SludgePump(f'.{ID}_sludge_pump', ins=sludge, init_with=init_with)
        self._mixed = self.ins[0].copy(f'{ID}_mixed')
        self._set_split_at_mc()

# _set_split_at_mc calculates the solubles and water content in the sludge using
# the moisture content in the sludege and assuming all solids go into the sludge
    def _set_split_at_mc(self):
        mixed = self._mixed
        eff, sludge = self.outs
        mixed.mix_from(self.ins)
        mc = self.sludge_moisture
        mixed_F_mass = mixed.F_mass
        if mixed_F_mass == 0: # empty streams
            split = 0
        else:
            mixed_mc = mixed.imass['Water']/mixed_F_mass
            if mixed_mc < mc: # not enough water in the feeds
                warn(f'The set `sludge_moisture` {mc} is smaller than the '
                     f'moisture content in the influent ({mixed_mc:.3f}) and is ignored.')
                sludge.copy_like(mixed)
                eff.empty()
                split = 0
                self.SKIPPED = True
            else:
                solubles, solids = self.solubles, self.solids
                sludge.copy_flow(mixed)
                solids_mass = sludge.imass[solids].sum()
                sludge.imass[('Water', *solubles)] *= \
                    (solids_mass/(1-mc)-solids_mass)/(mixed_F_mass-solids_mass)
                eff.mass = mixed.mass - sludge.mass
                split = mixed.mass.value.copy()
                idx = np.where(mixed.mass!=0)
                split[idx] = eff.mass[idx]/mixed.mass[idx]
                self.SKIPPED = False
        self._isplit = self.thermo.chemicals.isplit(split)

# _mc_at_split calculates moisture content from the solubles separation, and
# then compare the calculated mc with the target mc
    @staticmethod
    def _mc_at_split(split, solubles, mixed, eff, sludge, target_mc):
        eff.imass[solubles] = mixed.imass[solubles] * split
        sludge.imass[solubles] = mixed.imass[solubles] - eff.imass[solubles]
        mc = sludge.imass['Water'] / sludge.F_mass
        return mc-target_mc


    def _run(self):
        eff, sludge = self.outs
        solubles, solids = self.solubles, self.solids

        mixed = self._mixed
        mixed.mix_from(self.ins)
        eff.T = sludge.T = mixed.T
        sludge.copy_flow(mixed, solids, remove=True) # all solids go to sludge
        eff.copy_flow(mixed, solubles)

        self._set_split_at_mc()
        flx.IQ_interpolation(
            f=self._mc_at_split, x0=1e-3, x1=1.-1e-3,
            args=(solubles, mixed, eff, sludge, self.sludge_moisture),
            checkbounds=False)
        self._set_split_at_mc() #!!! not sure if still needs this

# the IQ_interpolation finds split in a specified range, and then calls _mc_at_split
# to calculate the mc based on solubles split. It then compare the calculated mc with
# the target mc. Through many iterations, it finds the split that yields the targeted
# mc. 
# !!! however, the IQ_interpolation does not talk to _set_split_at_mc to set the same
# split

# =============================================================================
#     def _cost(self):
#         if self.SKIPPED == False:
#             m_solids = self.outs[-1].F_mass
#             self.add_OPEX = {'Sludge disposal': m_solids*self.disposal_cost}
#             for p in (self.effluent_pump, self.sludge_pump): p.simulate()
#         else:
#             self.add_OPEX = {}
#             self.baseline_purchase_costs.clear()
#             self.power_utility.rate = 0
#         
# =============================================================================
 
#%% Pretreatment: pH adjustment by adding 1 M NaOH
base_dosing_path = ospath.join(data_path, 'sanunit_data/VFA/base_dosing_path.csv')
class BaseDosing(SanUnit):
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 base_conc = 1, titr_factor = 0.01225):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        
        self.base_conc = base_conc
        self.titr_factor = titr_factor
        # titr_factor calculated based on the experimental data from Wangsuk: 1.98 ml of 1 M NaOH into 160 ml of AD effluent.
        
    def _run(self):
        inf, base = self.ins
        eff, = self.outs
        
        base.ivol['H2O'] = inf.F_vol * self.titr_factor * 1000
        base.imass['Na'] = base.ivol['H2O'] * self.base_conc * 23 # MW of NaOH = 40 g/mol
        # by defining the base stream here, no need to initiate the stream in system.py
        eff = eff.mix_from(self.ins)
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
   
#%% Separation: redox-ED
redox_ED_path = ospath.join(data_path, 'sanunit_data/VFA/_redox_ED.csv')
class RedoxED(SanUnit):
    '''
    Redox-electrodialysis unit to upconcentrate volatile fatty acids from waste.
    
    Parameters
    ----------
    ins: Iterable(stream)
        fc_in is feeding channel inflow
        ac_in is accumulating channel inflow
    outs: Iterable(stream)
        fc_out is feeding channel effluent
        ac_out is accumulating channel effluent
    voltage: float
        cell voltage (V)
    op_time: float
        operational time (hr)
    fc_c3: float
        feeding channel propionate inflow concentration (M)
    fc_c4: float
        feeding channel butyrate inflow concentration (M)
    fc_c6: float
        feeding channel hexanoate inflow concentration (M)
    
    Unit conventions:
    -----------------
    concentration: mol/L or M
    volume: L
    voltage: V
    current: A
    time: hr
    
    Concentration calculations:
    ---------------------------
    Assume a balck box model for the ionic flux, continuously-stirred tank reactors
    for the feeding and accumulating channels. Model membrane transfer as 1st order
    reaction and solve a steady state mass balance for concentration.
    
    '''
    _N_ins = 2
    _N_outs = 2
    
    # !!!
    # Linear regression coefficients to extrapolate flux, calculated in excel
    # and imported as instance attributes. Linear regression coefficients
    # should be updated when new experimental data is used for future instances
    # !!! Check units for the mass balance and concentrations!
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 voltage=None, op_time=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        
        # Defining instance attributes, operation related for redox-ED
        # !!! When do the users input values? because this script is only run 
        # implicitly when a system is created
        
        data = load_data(path=redox_ED_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
        self.voltage = voltage
        self.area = self.compartment_length * self.compartment_width
        
    def _run(self):
        fc_in, ac_in = self.ins
        # !!! need to edit the inputs, fc_in is from the wastestream, ac_in is
        # a supporting electrolyte solution that only inlcudes sodium, potassium,
        # and chloride
        fc_out, ac_out = self.outs
        # Calculate the average flux from cell voltage
        # !!! how to reduce repetitiveness? when flux extrapolation constants
        # are imported in other ways instead of as class attributes in future,
        # need to change ways to retrieve this data as well
        c3_flux = self.c3_slope*self.voltage + self.c3_const
        c4_flux = self.c4_slope*self.voltage + self.c4_const
        c6_flux = self.c6_slope*self.voltage + self.c6_const
        
        # Calculate accumulating channel effluent concentration
        ac_out.imass['Propionate'] = ac_in.imass['Propionate'] + \
        c3_flux * self.area / self.flowrate 
        ac_out.imass['Butyrate'] = ac_in.imass['Butyrate'] + \
        c4_flux * self.area / self.flowrate
        ac_out.imass['Hexanoate'] = ac_in.imass['Hexanoate'] + \
        c6_flux * self.area / self.flowrate
        
        # Calcualte feeding channel effluent concentration        
        fc_out.imass['Propionate'] = fc_in.imass['Propionate'] - \
        c3_flux * self.area / self.flowrate
        fc_out.imass['Butyrate'] = fc_in.imass['Butyrate'] - \
        c4_flux * self.area / self.flowrate
        fc_out.imass['Hexanoate'] = fc_in.imass['Hexanoate'] - \
        c6_flux * self.area / self.flowrate
        
# !!! need to consider the mass transport of water?
        
# =============================================================================
#     def _init_lca(self):
#         self.construction = [Construction('carbon_cloth', linked_unit=self, 
#                                            item='CarbonCloth', 
#                                            quantity_unit='kg'),
#                              Construction('current_collector', linked_unit=self,
#                                            item='CurrentCollector', 
#                                            quantity_unit='kg'),
#                              Construction('silicone_spacer', linked_unit=self,
#                                            item='SiliconeSpacer',
#                                            quantity_unit='kg'),
#                              Construction('IEM', linked_unit=self,
#                                            item='IEM',
#                                            quantity_unit='m2'),    #is there a common format for units?
#                              Construction('housing', linked_unit=self,
#                                            item='Housing',
#                                            quantity_unit='kg'),
#                              Construction('piping', linked_unit=self,
#                                            item='Piping',
#                                            quantity_unit='m')]
# 
# # !!! change construction items when scaling up system. Now the materials are
# # based on bench top experiment.
#     
#     def _design(self):
#         design = self.design_results
#         constr = self.construction
#         design['CarbonCloth'] = constr[0].quantity = self.area
#         design['CurrentCollector'] = constr[1].quantity = self.area*self.num_current_collector
#         design['SiliconSpacer'] = constr[2].quantity = self.area*self.num_spacer
#         design['IEM'] = constr[3].quantity = self.area*self.num_IEM
#         design['Housing'] = constr[4].quantity # where would I input quantity value?
#         design['Piping'] = constr[5].quantity
#     
#     def _cost(self):
#         pass
# 
# =============================================================================
# =============================================================================
#     def _init_lca(self): 
#         self.construction = [Construction('carbon_cloth', linked_unit=self, 
#                                           item='CarbonCloth', 
#                                           quantity_unit='kg'),
#                              Construction('current_collector', linked_unit=self,
#                                           item='CurrentCollector', 
#                                           quantity_unit='kg'),
#                              Construction('silicon_spacer', linked_unit=self,
#                                           item='SiliconSpacer',
#                                           quantity_unit='kg'),
#                              Construction('electrolyte_membrane', linked_unit=self,
#                                           item='ElectrolyteMembrane',
#                                           quantity_unit='m2'),    #is there a common format for units?
#                              Construction('housing', linked_unit=self,
#                                           item='Housing',
#                                           quantity_unit='kg'),
#                              Construction('piping', linked_unit=self,
#                                           item='Piping',
#                                           quantity_unit='m')] 
#             
#     def _design(self):
#         design = self.design_results
#         constr = self.construction
# =============================================================================
        
        
        
        
        