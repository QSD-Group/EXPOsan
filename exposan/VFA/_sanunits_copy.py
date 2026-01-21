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

#%%
## some global constants, might include uncertainty range in the future

# conversion factors
m3_per_gal = 0.00378541
min_per_hr = 60

# assumption for the vonnecting pipes at WWTP
hazen_william_C = 110 # dimensionless
pipe_length = 150 # ft
min_velocity = 2.5 # ft, based on minimum flowrate suggestion from manual of practice No.8, 6th ed page 173

# assume positive gas pressure in anaerobic digestor
gas_pressure = 10/12 # ft

# densities
ss_weight = 8000 # density of stainless steel 316 in kg/m3


#%% Pretreatment: solids separation using centrifugation
centrifuge_path = ospath.join(data_path, 'sanunit_data/VFA/_centrifuge.csv')


@cost(basis='Power required', ID='SolidsSeparation', units='hp',
      cost=1730000, S=150, N = 'Number of centrifuges',
      CE=CEPCI_by_year[2014], n=0.8023)
# parameters in @cost are obtained from CapdetWorks, exponential is calculated by fitting unit costs of different power requirement, 
# the default BM of 1 is used because CapdetWorks lumps equipment and construction costs.
# the construction and equipment costs of a single centrifuge is dependent on the total required power in an exponential relationship
# the total power required can be calculated by the feeding rate multiplied by the power of centrifuge in hp/gpm
class SolidsSeparation(SanUnit):
    '''
    -Function of unit:
        A solids separation class for simulation of a solid bowl centrifuge 
        to separate solids from mixed stream.
    -Mass balance calculation:
        Separation split is calculated based on moisture content 
        in the sludge and solids separation rate, both defined by user.
        Assume solubles and water share the same split.
        If the moiture content in the feed is smaller than the targeted moiture 
        content in sludge, the target moiture content will be ignored.
    -Material in cost and design:
        
    -References:
        Reference units are G2RTSolidsSeparation developed by Zixuan Wang, 
        the SolidCetrifuge developed by Yoel, and Thickener by Yalin.
    
    Parameters
    ----------
    Ins: Iterable (stream)
        ins[0] is incoming waste stream.
        ins[1] is polymer addition.
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
    _N_ins = 2
    _N_outs = 2
    _units = {'Power required': 'hp'}
    power_range = (0, 200) 
    # solids loading rate for solid boal centrifuge determined from reported cases in the EPA design manual

    auxiliary_unit_names = ('feed_pump', 'dosing_pump', 'centrate_pump')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 removal_rate = 0.82, sludge_moisture = 0.926, polymer_dose = 0.0054, 
                 disposal_cost = 60, kW_per_m3_per_hr = 85.74, 
                 operating_hours = 24 * 365, labor_hour = 1.088, auger_cost = 500,
                 conveyor_belt_cost = 229.67, sludge_density = 1.03*998.2, 
                 feed_flowrate = 182 , material_repair_and_replacement_cost = 0.035, 
                 conveyer_power = 30.08, centrifuge_power = 180, auger_conveyer_cost = 0.0004, 
                 wages = 29.32, wall_thickness = 0.015, bowl_diameter = 0.88, bowl_length = 2.914, 
                 motor_power = 156.25, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.tss_removal = removal_rate
        self.sludge_moisture = sludge_moisture
        self.polymer_dose = polymer_dose
        self.disposal_cost = disposal_cost # 450 $/dry ton
        self.kW_per_m3_per_hr = kW_per_m3_per_hr # power consumption per gpm of sludge into the centrifuge
        self.operating_hours = operating_hours # this is number of operations per year
        self.labor_hour = labor_hour # num of labor hour required per operating hour
        self.auger_cost = auger_cost
        self.conveyor_belt_cost = conveyor_belt_cost
        self.sludge_density = sludge_density
        self.feed_flowrate = feed_flowrate # in gpm
        self.material_repair_and_replacement_cost = material_repair_and_replacement_cost
        self.conveyer_power = conveyer_power # HP
        # !!! replace place holder value for conveyer_power
        self.centrifuge_power = centrifuge_power # hp, typical values can be 150 - 250 hp based on case studies in EPA design manual
        self.auger_conveyer_cost = auger_conveyer_cost
        self.wages = wages
        self.wall_thickness = wall_thickness
        self.bowl_diameter = bowl_diameter
        self.bowl_length = bowl_length
        # sludge density is estimated by specific gravity of primary slidge + WAS and water density at 20 C from MtCalf&Eddy; density of AD effluent was not foudn in the book


# !!! Check what F_BM_default = 1 means, deleted scale up, ppl, and estreme arguments from parent class

# !!! in the line below, what components is it referring to? Should import from _components.py?        
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids])

    def _init_lca(self):
        self.construction = [Construction("stainless_steel", linked_unit=self,
                                          item = "StainlessSteel", 
                                          quantity_unit= "kg"),
                             Construction("electric_motor", linked_unit=self,
                                          item = "ElectricMotor",
                                          quantity_unit= "ea"),
                             Construction("auguer_conveyor", linked_unit=self,
                                          item = "AugerConveyor", 
                                          quantity_unit= "m"),
                             ]
    #!!! unit of conveyor belt should be m but is not calculated in LCA       
        
    def _run(self):
        AD_effluent, polymer = self.ins # index or comma needed when there is only one stream in inlets
        liquid_stream, solid_stream = self.outs
# This following is defining the subset of components that later could be accessed through []
        solubles, solids = self.solubles, self.solids
        TL_in = AD_effluent.F_mass - AD_effluent.imass[solids].sum()
        mc_in = TL_in / AD_effluent.F_mass # including both water and solubles in the moisture content
        mc_out = self.sludge_moisture # this moisture data should assume the total mass of water and solubles
# !!! Below the logic check is for moisture content in the solids but the focus here is the liquid stream
        # if mc_in < mc_out*0.999:
        if mc_in < mc_out:
            mc_out = mc_in
        
        solid_stream.imass[solids] = AD_effluent.imass[solids] * self.tss_removal
        TS_out = solid_stream.imass[solids].sum() # total soilds in the solids stream
        TL_out = TS_out / (1 - mc_out) * mc_out # total solubles and water mass flowrate in the solids stream
        solid_stream.imass[solubles] = AD_effluent.imass[solubles] * TL_out / TL_in
        # solid_stream.imass['H2O'] = TL_out - solid_stream.imass[solubles].sum()
        # the above line assumes that the solubles (including water and soluble chemicals) partition together and by the same ratio
        #liquid_stream.imass['H2O'] = TS_out*(1-mc_out) * mc_out - liquid_stream.imass[solubles].sum()
        
        liquid_stream.imass[solids] = AD_effluent.imass[solids] * (1 - self.tss_removal)
        liquid_stream.imass[solubles] = AD_effluent.imass[solubles] - solid_stream.imass[solubles]
        #liquid_stream.imass['H2O'] = AD_effluent.imass['H2O'] - solid_stream.imass['H2O']
        
        polymer.imass['Polymer'] = AD_effluent.imass[solids].sum() * self.polymer_dose
        
        pipe_diameter = sqrt(4 * self.feed_flowrate * 35.3147/3600 / pi / min_velocity)
        friction_head = 3.02 * pipe_length * (min_velocity ** 1.85) * \
        (hazen_william_C ** (-1.85)) * (pipe_diameter ** (-1.17)) # Using equation ESI-7 in the SI of Brian Shoener's 2016 paper
        feed_pump_p = friction_head + gas_pressure + 10 # include 10ft water column for pressure head
        centrate_pump_p = friction_head + 101325 # Pa
        feed_pump = self.auxiliary('feed_pump', cls = SludgePump, ins = self.ins[0].copy('feed_pump_in'), P = feed_pump_p)
        dosing_pump = self.auxiliary('dosing_pump', cls = Pump, ins = self.ins[1].copy('dose_pump_in'))
        centrate_pump = self.auxiliary('centrate_pump', cls = SludgePump, ins = self.outs[0].copy('centrate_pump_out'), P = centrate_pump_p)

        # !!! the pump_pressure argument is the pressure of the output stream, assume atmospheric pressure for the feeding
        # and the dosing pumps.
        #!!! estimate the pressure of the centrate pu,p considering 3m (10ft) of elevation lift and the frictional loss in the pipes
        # calculate using the energy balance equation in fluid mechanics
        self.ins[0].P = centrate_pump_p
        
    def _design(self):
        design = self.design_results
        design['Feeding_rate'] = feeding_rate = self.F_vol_in
        design['Power_required'] = power_required = feeding_rate * self.centrifuge_power
        lower_bound, upper_bound = self.power_range
        if power_required < lower_bound:
            lb_warning(self, 'Power_required', power_required, 'hp', lower_bound)
        design['Number_of_centrifuges'] = ceil(feeding_rate/upper_bound)
        # assuming the power requirement per gpm is 1 hp for each centrifuge, and using the total feeding rate and the typical total power that one centrifuge can provide from CapdetWorks,
        # we calculate the number of centrifuges needed.
        design['Number_of_motor'] = self.construction[1].quantity \
            = design['Number_of_centrifuges']
        design['Auger_conveyor'] = self.construction[2].quantity = 10 #!!! assume 10-m auger conveyor
        
        out_d = self.bowl_diameter
        in_d = out_d - self.wall_thickness
        design['Stainless_Steel'] = self.construction[0].quantity = (out_d**2 - in_d**2) / 4 * pi * self.bowl_length * ss_weight \
            + out_d * out_d / 2 * self.bowl_length * self.wall_thickness
        # centrifuge bowl and hosuing stianless steel requirement
        
        
        self.power_utility(self.F_vol_in * self.kW_per_m3_per_hr + self.conveyer_power * 0.7457 *
                           self.operating_hours)
        #!!! check calculation
        # this is electricity consumption not including pumping, 1 hp = 0.7457 kW
        # this uses the energy consumption per m3/hr from the EPA design manual, which might include more that the energy requirement of the centrifuge pump
        # on top of that add the energy consumption of the conveyor
        
        
    def _cost(self):
        # capital cost for auger and conveyor belt for sludge removal
        ts_out = self.outs[-1].F_mass
        self.baseline_purchase_costs['Centrifuge'] = 31188 * self.centrifuge_power ** 0.8023 * CEPCI_by_year[2023] / CEPCI_by_year[2014]
        self.F_BM['Centrifuge']= 1
        self.baseline_purchase_costs['Auger_conveyer'] = self.baseline_purchase_costs['Centrifuge'] * self.auger_conveyer_cost 
        # using a ratio here becasue it's time consuming to find accurate price that statisfied the solids loading rate as well as the scaling factors
        self.F_BM['Auger_conveyer'] = 1
        # using 1 because we icluded the other construction costs in the purchase cost
        tot_equip_and_constr = sum(self.baseline_purchase_costs.values())

        # OPEX = sludge disposal, material replacement (include as a ratio), labor, 
        # electricity (scale with flow in the design fxn), and polymer dosage (scale with flow int he design fxn)
        self.add_OPEX = {'Sludge_disposal': ts_out * 24 * self.disposal_cost * 0.00110231, # $/day
                         'Labor': self.wages * self.operating_hours / 365 * self.labor_hour, # $/day
                         'Maintenance': tot_equip_and_constr * self.material_repair_and_replacement_cost / 365} # two operators per shift, units in $/day
        # 0.00110231 short ton per kg
        # units all converted to $ per day: sludge disposal - kg/hr * hr/day * $/kg; labor - $/hr * hr/year / (day/year) * hr/hr
        # !!! find costs in the add_OPEX dictionary
        for p in (self.feed_pump, self.dosing_pump, self.centrate_pump): p.simulate()
        # !!! should include multiple feed, centrate, and dosing pumps matching with the number of centrifuge units?
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
    all insolubles components will all go to the retentate.

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
                 disposal_cost=450): # converting from $/U.S. ton
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
# the moisture content in the sludge and assuming all solids go into the sludge
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
        # might be an inssue with line 346 calculating self._isplit   

# _mc_at_split calculates moisture content from the solubles separation, and
# then compare the calculated mc with the target mc
    @staticmethod
    def _mc_at_split(split, solubles, mixed, eff, sludge, target_mc):
        eff.imass[solubles] = mixed.imass[solubles] * split
        sludge.imass[solubles] = mixed.imass[solubles] - eff.imass[solubles]
        mc = sludge.imass['Water'] / sludge.F_mass
        return mc-target_mc

    def _run(self):
        effluent, sludge = self.outs
        solubles, solids = self.solubles, self.solids
        mixed = self._mixed
        mixed.mix_from(self.ins)
        effluent.T = sludge.T = mixed.T
        sludge.copy_flow(mixed, solids, remove=True) # all solids go to sludge
        effluent.copy_flow(mixed, solubles)
        self._set_split_at_mc()
        flx.IQ_interpolation(
            f=self._mc_at_split, x0=1e-3, x1=1.-1e-3,
            args=(solubles, mixed, effluent, sludge, self.sludge_moisture),
            checkbounds=False)
        self._set_split_at_mc() #!!! not sure if still needs this
        
    def _design(self):
        design = self.design_results
        pass

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
    _N_ins = 3
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 base_conc = 1, titr_factor = 0.01225):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        
        self.base_conc = base_conc
        self.titr_factor = titr_factor
        # titr_factor calculated based on the experimental data from Wangsuk: 1.98 ml of 1 M NaOH into 160 ml of AD effluent.
        # AD_red only used to calculate base addition and not participate in mass balance
    def _run(self):
        inf, ad_ref, base = self.ins
        effluent, = self.outs
        
        base_vol = ad_ref.F_vol * self.titr_factor
        base.ivol['H2O'] = base_vol
        base.imol['Na'] = base_vol * self.base_conc * 1000 # kmol/hr
        # by defining the base stream here, no need to initiate the stream in system.py
        effluent = effluent.mix_from((inf, base))
        
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
        
        fc_out.copy_like(fc_in)
        ac_out.copy_like(ac_in)

        fc_out.imass['H2O'] = fc_in.imass['H2O']
        ac_out.imass['H2O'] = ac_in.imass['H2O']
        
        # Calculate the average flux from cell voltage
        # !!! how to reduce repetitiveness? when flux extrapolation constants
        # are imported in other ways instead of as class attributes in future,
        # need to change ways to retrieve this data as well
        c3_flux = self.c3_slope*self.voltage + self.c3_const
        c4_flux = self.c4_slope*self.voltage + self.c4_const
        c6_flux = self.c6_slope*self.voltage + self.c6_const
        
        # Calculate accumulating channel effluent concentration
        ac_out.imol['Propionate'] = ac_in.imol['Propionate'] + \
        c3_flux * self.area / self.flowrate 
        ac_out.imol['Butyrate'] = ac_in.imol['Butyrate'] + \
        c4_flux * self.area / self.flowrate
        ac_out.imol['Hexanoate'] = ac_in.imol['Hexanoate'] + \
        c6_flux * self.area / self.flowrate
        
        # Calcualte feeding channel effluent concentration        
        fc_out.imol['Propionate'] = fc_in.imol['Propionate'] - \
        c3_flux * self.area / self.flowrate
        fc_out.imol['Butyrate'] = fc_in.imol['Butyrate'] - \
        c4_flux * self.area / self.flowrate
        fc_out.imol['Hexanoate'] = fc_in.imol['Hexanoate'] - \
        c6_flux * self.area / self.flowrate
        
# !!! need to consider the mass transport of water?
# =============================================================================
#         J3 = self.c3_slope * self.voltage + self.c3_const  # mol/m2/s
#         J4 = self.c4_slope * self.voltage + self.c4_const
#         J6 = self.c6_slope * self.voltage + self.c6_const
# 
#         A = self.area  # m2
#         n3_tr = J3 * A * 3600  # mol/hr
#         n4_tr = J4 * A * 3600
#         n6_tr = J6 * A * 3600
# 
#         def transfer(solute_id, n_tr):
#             n_avail = fc_in.imol[solute_id]
#             n_move = min(max(n_tr, 0.0), n_avail)  # no negative transfer, cannot exceed available
#             breakpoint()
#             fc_out.imol[solute_id] = n_avail - n_move
#             ac_out.imol[solute_id] = ac_in.imol[solute_id] + n_move
# 
#         transfer('Propionate', n3_tr)
#         transfer('Butyrate', n4_tr)
#         transfer('Hexanoate', n6_tr)
# =============================================================================
        
        # Set Na from electroneutrality (Cl conservative)
        # Na+ = Cl- + sum(VFA-)
        fc_vfa = sum(fc_out.imol[i] for i in ('Propionate', 'Butyrate', 'Hexanoate'))
        ac_vfa = sum(ac_out.imol[i] for i in ('Propionate', 'Butyrate', 'Hexanoate'))

        fc_out.imol['Na'] = max(fc_out.imol['Cl'] + fc_vfa, 0.0)
        ac_out.imol['Na'] = max(ac_out.imol['Cl'] + ac_vfa, 0.0)


        
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
        
        
        
        
        