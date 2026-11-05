#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:27:03 2025

@author: blues
"""

from qsdsan import SanUnit, Construction, WasteStream, Stream, Chemical, Component, \
Components, set_thermo as qs_set_thermo, CEPCI_by_year
#from qsdsan.sanunits import pumping
from qsdsan.utils import ospath, data_path, load_data, price_ratio, components, auom
# from exposan.g2rt._sanunits import G2RTSolidsSeparation
from biosteam import Splitter
from biosteam.units.decorators import cost
from qsdsan.sanunits import SludgePump, WWTpump, Pump, AnMBR
import numpy as np, flexsolve as flx
from warnings import warn
from biosteam.exceptions import lb_warning
from math import ceil, pi, sqrt

__all__ = ('SolidsSeparation', 'Microfiltration', 'BaseDosing', 'RedoxED',)

#%% Pretreatment: solids separation using centrifugation
## some global constants, might include uncertainty range in the future

# conversion factors
m3_per_gal = 0.00378541
min_per_hr = 60
_ft_to_m = 0.3048
_kg_to_short_ton = 0.00110231

# assumption for the vonnecting pipes at WWTP
hazen_william_C = 110 # dimensionless for preassure head in ft4
pipe_length = 150 # ft
min_velocity = 2.5 # ft, based on minimum flowrate suggestion from manual of practice No.8, 6th ed page 173
_g = 9.81 # m/s2

# assume positive gas pressure in anaerobic digestor
gas_pressure = 10/12 # ft
water_column = 10 # ft

# densities
ss_density = 8000 # density of stainless steel 316 in kg/m3
water_density = 998.2 # kg/m3 at 20 degrees C from MtCalf & Eddy


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
    power_range = (0, 200) # hp
    # solids loading rate for solid boal centrifuge determined from reported cases in the EPA design manual

    auxiliary_unit_names = ('feed_pump', 'dosing_pump', 'centrate_pump')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 tss_removal = 0.926, sludge_moiture = 0.822, 
                 polymer_dose = 0.0054, kW_per_m3_per_hr = 85.743, 
                 operating_hours = 5127.855263, 
                 centrifuge_power = 1, 
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
# =============================================================================
#         self.tss_removal = removal_rate
#         self.sludge_moisture = sludge_moisture
#         self.polymer_dose = polymer_dose
#         self.disposal_cost = disposal_cost # 450 $/dry ton
#         self.kW_per_m3_per_hr = kW_per_m3_per_hr # power consumption per gpm of sludge into the centrifuge
#         self.operating_hours = operating_hours # this is number of operations per year
#         self.labor_hour = labor_hour # num of labor hour required per operating hour
#         self.auger_conveyer_cost = auger_conveyer_cost
#         self.material_repair_and_replacement_cost = material_repair_and_replacement_cost
#         self.conveyer_power = conveyer_power # HP
#         # !!! replace place holder value for conveyer_power
#         self.centrifuge_power = centrifuge_power # hp/gpm, typical values can be 0.5 - 2 hp/gpm based on case studies in EPA design manual
#         self.auger_conveyer_cost = auger_conveyer_cost
#         self.wages = wages
#         self.wall_thickness = wall_thickness # m
#         self.bowl_diameter = bowl_diameter # m 
#         self.bowl_length = bowl_length # m
#         # sludge density is estimated by specific gravity of primary slidge + WAS and water density at 20 C from MtCalf&Eddy; density of AD effluent was not foudn in the book
# 
# =============================================================================
        data = load_data(path=centrifuge_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

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
        
        # pipe_diameter = sqrt(4 * self.feed_flowrate * 35.3147/3600 / pi / min_velocity)
        # friction_head = 3.02 * pipe_length * (min_velocity ** 1.85) * \
        # (hazen_william_C ** (-1.85)) * (pipe_diameter ** (-1.17)) # Using equation ESI-7 in the SI of Brian Shoener's 2016 paper
        # feed_pump_p = (friction_head + gas_pressure + 10) * _ft_to_m * water_density * _g # include 10ft water column for pressure head
        # centrate_pump_p = friction_head # Pa
        
        feed_pump = self.auxiliary('feed_pump', cls = SludgePump, ins = self.ins[0].copy('feed_pump_in'))  # pipe length is not modifiable when using the self.
        # feeding from the anaerobic digestor at 10 ft deep of the water column, explains the minus sign
        dosing_pump = self.auxiliary('dosing_pump', cls = SludgePump, ins = self.ins[1].copy('dose_pump_in')) # tried ot use WWTPump but got a string index out of range error
        centrate_pump = self.auxiliary('centrate_pump', cls = SludgePump, ins = self.outs[0].copy('centrate_pump_out')) # L_s = 0, L_d = pipe_length)
        #!!! might consider adding pumps as individual units because modifications with auxiliary units are limited
        
        # !!! the pump_pressure argument is the pressure of the output stream, assume atmospheric pressure for the feeding
        # and the dosing pumps.
        #!!! estimate the pressure of the centrate pump considering 3m (10ft) of elevation lift and the frictional loss in the pipes
        # calculate using the energy balance equation in fluid mechanics
        
    def _design(self):
        design = self.design_results
        design['Feeding_rate'] = feeding_rate = self.F_vol_in / (24*60) / m3_per_gal
        design['Power_required'] = power_required = feeding_rate * self.centrifuge_power
        lower_bound, upper_bound = self.power_range
        if power_required < lower_bound:
            lb_warning(self, 'Power_required', power_required, 'hp', lower_bound)
        design['Number_of_centrifuges'] = ceil(power_required/upper_bound)
        # assuming the power requirement per gpm is 1 hp for each centrifuge, and using the total feeding rate and the typical total power that one centrifuge can provide from CapdetWorks,
        # we calculate the number of centrifuges needed.
        design['Number_of_motor'] = self.construction[1].quantity \
            = design['Number_of_centrifuges']
        design['Auger_conveyor'] = self.construction[2].quantity = 10 #!!! assume 10-m auger conveyor
        
        out_d = self.bowl_diameter
        in_d = out_d - self.wall_thickness
        design['Stainless_Steel'] = self.construction[0].quantity = (out_d**2 - in_d**2) / 4 * pi * self.bowl_length * ss_density \
            + out_d * out_d / 2 * self.bowl_length * self.wall_thickness
        # centrifuge bowl and hosuing stianless steel requirement
        
        
        self.power_utility((self.F_vol_in / 24 * self.kW_per_m3_per_hr + self.conveyer_power * 0.7457))
        # here only power in kW is needed. total energy consumption will be calculated using operating hours in system
        #!!! how to use different operating hours for different units?
        #!!! check calculation
        # this is electricity consumption not including pumping, 1 hp = 0.7457 kW
        # this uses the energy consumption per m3/hr from the EPA design manual, which might include more that the energy requirement of the centrifuge pump
        # on top of that add the energy consumption of the conveyor
        
        
    def _cost(self):
        # capital cost for auger and conveyor belt for sludge removal
        ts_out = self.outs[-1].F_mass - self.outs[-1].imass['Water']
        design = self.design_results
        self.baseline_purchase_costs['Centrifuge'] = 31188 * design['Power_required'] ** 0.8023 * CEPCI_by_year[2023] / CEPCI_by_year[2014] # the cost curve was obtained based on 2014 data in Capdet Works
        self.F_BM['Centrifuge']= 1
        self.baseline_purchase_costs['Auger_conveyer'] = self.baseline_purchase_costs['Centrifuge'] * self.auger_conveyer_cost 
        # using a ratio here becasue it's time consuming to find accurate price that statisfied the solids loading rate as well as the scaling factors
        self.F_BM['Auger_conveyer'] = 1
        # using 1 because we icluded the other construction costs in the purchase cost
        tot_equip_and_constr = sum(self.baseline_purchase_costs.values())

        # OPEX = sludge disposal, material replacement (include as a ratio), labor, 
        # electricity (scale with flow in the design fxn), and polymer dosage (scale with flow int he design fxn)
        self.add_OPEX = {'Sludge_disposal': ts_out * 24 * self.disposal_cost * _kg_to_short_ton, # $/day 
                         'Labor': self.wages * 8 * self.labor_hour, # assume operaters work for 8 hours per day, [$/day]
                         'Maintenance': tot_equip_and_constr * self.material_repair_and_replacement_cost / 365} # two operators per shift, units in $/day
        # 0.00110231 short ton per kg
        # units all converted to $ per day: sludge disposal - kg/hr * hr/day * $/kg; labor - $/hr * hr/year / (day/year) * hr/hr
        # !!! find costs in the add_OPEX dictionary
        for p in (self.feed_pump, self.dosing_pump, self.centrate_pump): p.simulate()
        # !!! should include multiple feed, centrate, and dosing pumps matching with the number of centrifuge units?
#%% Pretreatment: solids separation using microfiltration       
_ft2_to_m2 = auom('ft2').conversion_factor('m2')
_m3_to_gal = auom('m3').conversion_factor('gal')
_m_to_ft = 3.28084
_H2O_spec_grav = 9.81 #[kN/m3]
_kPa_to_psi = 0.145038
_sludge_spec_grav = 1.005

microfiltration_path = ospath.join(data_path, 'sanunit_data/VFA/_microfiltration.csv')

class Microfiltration(AnMBR):
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
    
    
    A cross-flow microfiltration unit with multi-tube membranes.
    '''

    SKIPPED = False
    _graphics = Splitter._graphics
    _ins_size_is_fixed = False
    _N_outs = 2
    auxiliary_unit_names = ('effluent_pump', 'sludge_pump')
    
    # Equipment-related parameters, values all for cross-flow, multitube configuration
    _cas_per_tank_spare = 2

    _mod_surface_area = {
        'hollow fiber': 370,
        'flat sheet': 1.45/_ft2_to_m2,
        'multi-tube': 32/_ft2_to_m2
        }

    _mod_per_cas = None
    _mod_per_cas_range = {
        'hollow fiber': (30, 48), # min, max
        'flat sheet': (150, 200),
        'multi-tube': (44, 48)
        }


    _cas_per_tank = None
    _cas_per_tank_range = (16, 22)
    
    # Operation related parameters for cross-flow, multitube configuration
    _recir_ratio = 2.25 # from the 0.5-4 uniform range in ref [1]
    _SGD = 0.625 # from the 0.05-1.2 uniform range in ref [1]
    
    # !!! update parameter values combining other sources
    
    _membrane_configuration = 'cross-flow'
    _membrane_type = 'multi-tube'
    _membrane_material = 'PET' # membrane material could be plastics or ceramic, according to Choener et al. for multi-tube membranes
    # allowed plastic material include: PES, PVDF, PET, PTFE. Using plastic membrane material here to be consistant with product information
    # of the life cycle inventory item.
    
    pumps = ('naocl', 'perm', 'backwash', 'retent', 'sludge', 'citric', 'bisulfite' )

    _discharge_pipe_length = 10 #ft #!!! check pipe length

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',solids=(), solubles=(), organics=(), N_train=2, **kwargs): # converting from $/U.S. ton
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with=init_with)
        
# =============================================================================
#         self.sludge_solids_conc = sludge_solids_conc # kg/m3 from metcalf&eddie
#         #!!! update sludge moisture data, which data is more easily available: sludge moisture or solids concentration?
#         self.disposal_cost = disposal_cost
#         self.organics_split = organics_split
#         # approximations based on rejection data from MtCalf&Eddy Table 11-31,
#         # it means the portion of compounds retained in the sludge
#         # !!! check if the splitting here causes a big portion of organic acid wasted
#         self.bw_flux_ratio = bw_flux_ratio
#         self.bw_time_ratio = bw_time_ratio
#         self.m_surface_area = m_surface_area # from Pentair manufacture website in units of m2
#         self.TMP = TMP
#         self.TMP_back = TMP_back
#         
# =============================================================================
        
        data = load_data(path=microfiltration_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)        
        
        cmps = self.components
        self.solids = solids or tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids])
        self.organics = tuple(('Propionate', 'Butyrate', 'Hexanoate'))
        self.salts = tuple(('Na', 'Cl', 'K'))
        ID = self.ID
         
        self._inf = inf = self.ins[0].copy() # this stream will be preserved (i.e., no reaction)
        self._retent = inf.copy(f'{ID}_retent') # for pump design, this will be the retentate which get recirculated
        #self.Q_mgd = self.ins[0].F_vol*_m3_to_gal*24/1e6
        # Add '.' in ID for auxiliary units
        # self.effluent_pump = SludgePump(f'.{ID}_eff_pump', ins=eff, init_with=init_with)
        # self.sludge_pump = SludgePump(f'.{ID}_sludge_pump', ins=sludge, init_with=init_with)
    
        self.N_train = N_train    
    
    def _init_lca(self):
        self.construction = [Construction("stainless_steel", linked_unit=self,
                                          item = "StainlessSteel",
                                          quantity_unit= "kg"),
                             Construction("electric_motor", linked_unit=self,
                                          item = "ElectricMotor",
                                          quantity_unit= "ea"),
                             Construction("HDPE", linked_unit=self,
                                          item = "Hdpe",
                                          quantity_unit= "ea"),
                             Construction("membrane_module", linked_unit = self,
                                          item = "MembraneModule",
                                          quantity_unit = "ea"),
                             ]
    
    def _run(self):
        inf, naocl, citric, bisulfite = self.ins    
        # inf is centrifuge effluent, not including recycled in the mass balance, treating as a black box model. The recycled will be internally calculated for pumps.
        perm, sludge = self.outs
        
        # TODO: JF added this
        self._inf = inf
        
        # Chemicals for cleaning, assume all chemicals will be used up, cost of chemicals are calculated in systems.py
        # 2.2 L/yr/cmd of 12.5 wt% solution (15% vol)
        naocl.empty()
        naocl.imass['NaOCl', 'Water'] = [0.125, 1-0.125]
        naocl.F_vol = (2.2/1e3/365/24) * (inf.F_vol*24) # m3/hr solution

        # 0.6 L/yr/cmd of 100 wt% solution, 13.8 lb/kg
        citric.empty()
        citric.ivol['CitricAcid'] = (0.6/1e3/365/24) * (inf.F_vol*24) # m3/hr pure

        # 0.35 L/yr/cmd of 38% solution, 3.5 lb/gal
        bisulfite.empty()
        bisulfite.imass['Bisulfite', 'Water'] = [0.38, 1-0.38]
        bisulfite.F_vol = (0.35/1e3/365/24) * (inf.F_vol*24) # m3/hr solution
        
        ## mass balance for separation into the permeate and sludge streams
        # inf.split_to(perm, sludge, self._isplit.data)
        solids, solubles, salts, organics = self.solids, self.solubles, self.salts, self.organics
        org_split = self.organics_split
        sludge.imass[solubles] = inf.imass[solubles] * 0
        sludge.imass[organics] = inf.imass[organics] * org_split
        sludge.imass[solids] = inf.imass[solids] * 1
        # assume all solids rejected and all salts pass through
        
        perm.imass['Water'] = inf.imass['Water'] * self.water_recovery
        perm.imass[solubles] = inf.imass[solubles] - sludge.imass[solubles]
        perm.imass[organics] = inf.imass[organics] * (1-org_split)
        perm.imass[solids] = inf.imass[solids] - sludge.imass[solids]
        
        sludge.imass['Water'] = inf.imass['Water'] - perm.imass['Water']

        # flowrate calculation for retentate
        self._compute_mod_case_tank_N()
        Q_retent = self._compute_liq_flows()
        retent = self._retent
        retent.F_mass *= Q_retent / inf.F_vol
        self._retent = retent
        #!!! is it okay to assume the retentate has similar composition to raw influent?
        #This is at steady state so we assume there has been retentate circulating in the system from the start?
        
    def _compute_liq_flows(self):
        inf = self.ins[0]
        self.cross_section_area = (self.m_diameter / 2) ** 2 * pi # membrane diameter in m
        Q_cross_flow = self.cross_section_area * self.cross_section_velocity
        Q_retent = self.N_mod_tot * Q_cross_flow - inf.F_vol
        # check if the recirculation ratio matches with numbers in Shoener 2016
        return Q_retent
    
    def _compute_mod_case_tank_N(self):
        N_mod_min, N_mod_max = self.mod_per_cas_range[self.membrane_type]
        N_cas_min, N_cas_max = self.cas_per_tank_range

        mod_per_cas, cas_per_tank = N_mod_min, N_cas_min
        
        # flux adjustment
        SA = (self.N_train-1) * self.cas_per_tank * self.mod_per_cas * self.m_surface_area
        J = self._inf.F_vol*1e3/SA
        J_max = self.forward_flux
        J_adj = J_max * (1-self.bw_time_ratio) + self.bw_flux_ratio * J_max * self.bw_time_ratio
        # Total flux of the microfiltration unit is adjusted based on the backwash duration and flux
        while J > J_adj:
            mod_per_cas += 1
            if mod_per_cas == N_mod_max + 1:
                if cas_per_tank == N_cas_max + 1:
                    self.N_train += 1
                    mod_per_cas, cas_per_tank = N_mod_min, N_cas_min
                else:
                    cas_per_tank += 1
                    mod_per_cas = N_mod_min

        self._mod_per_cas, self._cas_per_tank = mod_per_cas, cas_per_tank
        
    def _design(self):
        constr = self.construction
        # design membranes
        D = self.design_results
        D['Treatment train'] = self.N_train
        D['Cassette per train'] = self.cas_per_tank
        D['Module per cassette'] = self.mod_per_cas
        D['Total membrane modules'] = constr[3].quantity = self.N_mod_tot
        D['Total membrane area'] = self._design_multi_tube()
        
        # design pumps and pipes
        pipe, pumps, hdpe = self._design_pump()
        D['Pump_pipe_stainless_steel'] = pipe
        D['Pump stainless steel'] = pumps
        D['Pump chemical storage HDPE'] =constr[2].quantity = hdpe * 950 # volume * density; HDPE density average of 930 -970 kg/m3
        D['Pump motor'] = constr[1].quantity = 2*self.cas_per_tank + 5 # permeate, backwash, retentate, sludge, naocl, citric, and bisulfate pumps, refer to _design_pump, N_pump
        
        constr[0].quantity = pipe + pumps
        
    def _design_multi_tube(self):
        return self.N_mod_tot * self.m_surface_area
        
        # !!! Add a feeding tank or a reception tank for the permeate?
        
    def _design_pump(self):
        ID, ins, outs = self.ID, self.ins, self.outs
        m_config, pumps = \
            self._membrane_configuration, self.pumps
        
        ins_dct = {
            'perm': outs[0].proxy(f'{ID}_perm'),
            'backwash': self.bw_flux_ratio * ins[0].F_vol/outs[0].F_vol * \
                outs[0].proxy(f'{ID}_perm'),
            'retent': self._retent,
            'sludge': outs[1].proxy(f'{ID}_sludge'),
            'naocl': ins[1].proxy(f'{ID}_NaOCl'),
            'citric': ins[2].proxy(f'{ID}_citric'),
            'bisulfite': ins[3].proxy(f'{ID}_bisulfite'),
            }
        # influent for the backwash pump is scaled to higher than the influent 
        # flowrate but using the permeate
        type_dct = {
            'perm': f'permeate_{m_config}',
            'backwash':f'permeate_{m_config}',
            'retent': 'retentate_CSTR',
            'sludge': 'sludge',
            'naocl': 'chemical',
            'citric': 'chemical',
            'bisulfite': 'chemical',
            }
        inputs_dct = {
            'perm': (self.cas_per_tank, self._discharge_pipe_length, self.TMP * _kPa_to_psi, False),
            'backwash': (self.cas_per_tank, self._discharge_pipe_length, self.TMP_back * _kPa_to_psi, False),
            'retent': (1,),
            'sludge': (1,),
            'naocl': (1,),
            'citric': (1,),
            'bisulfite': (1,),
            }
        # D (or depth of tank in WWTPump will not be used when aerobic_filter is False)
        # TMP is expressed in ft, backwash TMP is higher than that of forward TMP
        # !!! There is no depth of the reactor needed in microfiltration,
        # however the WWTPump class was written for AnMBR which include 
        # bioreactor tank. For this unit, can we just not pass this argument? 
        # Same thing applies to self.include_aerobic_filter
        
        for i in pumps:
            if hasattr(self, f'{i}_pump'):
                p = getattr(self, f'{i}_pump')
                setattr(p, 'add_inputs', inputs_dct[i])
            else:
                # Add '.' in ID for auxiliary units
                ID = f'.{ID}_{i}'
                capacity_factor=2. if i=='perm' else self.recir_ratio if i=='retent' else 1.
                pump = WWTpump(
                    ID=ID, ins=ins_dct[i], pump_type=type_dct[i],
                    Q_mgd=None, add_inputs=inputs_dct[i],
                    capacity_factor=capacity_factor,
                    include_pump_cost=True,
                    include_building_cost=False,
                    include_OM_cost=True,
                    )
                setattr(self, f'{i}_pump', pump)
# !!! should pass self.Q_mgd to argument in WWTPump?        
        pipe_ss, pump_ss, chemical_hdpe = 0., 0., 0.
        for i in (*pumps,):
            p = getattr(self, f'{i}_pump')
            if p == None: continue
            p.simulate()
            p_design = p.design_results               
            pipe_ss += p_design['Pump pipe stainless steel']
            pump_ss += p_design['Pump stainless steel']
            if 'Pump chemical storage HDPE' in p_design.keys():
                chemical_hdpe += p_design['Pump chemical storage HDPE']
        return pipe_ss, pump_ss, chemical_hdpe
    
    def _cost(self):
        D, C = self.design_results, self.baseline_purchase_costs
        ### Capital ###
        # Membrane
        C['Membrane'] = self.membrane_unit_cost * D['Total membrane area'] / _ft2_to_m2 * CEPCI_by_year[2023] / CEPCI_by_year[2014]
        # !!! update membrane unit cost, might be in per module, per area or per volume
        # Pump
        pumps, add_OPEX = self.pumps, self.add_OPEX
        pump_cost, building_cost, opex_o, opex_m = 0., 0., 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            if p == None:
                continue
            p_cost, p_add_opex = p.baseline_purchase_costs, p.add_OPEX
            pump_cost += p_cost['Pump']
            building_cost += p_cost['Pump building']
            opex_o += p_add_opex['Pump operating']
            opex_m += p_add_opex['Pump maintenance']

        C['Pumps'] = pump_cost
        C['Pump building'] = building_cost
        C['Storage HDPE'] = D['Pump chemical storage HDPE'] * self.hdpe_cost
        add_OPEX['Pump operating'] = opex_o
        add_OPEX['Pump maintenance'] = opex_m
        
        # Power for pumping
        pumping_power = 0.
        for ID in self.pumps: #!!! check if cost/power of AF_pump/AeF_pump included in AF/AeF
            p = getattr(self, f'{ID}_pump')
            if p is None: continue
            pumping_power += p.power_utility.rate * self.bw_time_ratio if ID == 'backwash' else p.power_utility.rate
        self.power_utility.rate = pumping_power
        # note: power_utility.rate is the aggregate of power_utility.consumption and 
        # power_utility.production, it has the same unit as power consumption and production
        
        self.add_OPEX = {'Sludge_disposal': (self.outs[1].F_mass - self.outs[1].imass['Water']) * 24 * self.disposal_cost * _kg_to_short_ton, # $/day
                         'Labor': self.wages * 8, # assume 1 operater work for 8 hours per day, [$/day]
                         'Maintenance': C['Membrane'] * self.material_repair_and_replacement_cost / 365} # two operators per shift, units in $/day
#!!! MtCalf&Eddy Table 11-30 has typical energy consumption for microfiltration, might use that for sanity check
#!!! compare WWTPump with biosteam pump, this class does not include the pump motor?
 
    @property
    def Q_mgd(self):
        '''
        [float] Influent volumetric flow rate in million gallon per day, [mgd].
        '''
        return self.ins[0].F_vol*_m3_to_gal*24/1e6
    
#%% Pretreatment: pH adjustment by adding 1 M NaOH
_ft_to_m = auom('ft').conversion_factor('m')
_ft3_to_gal = auom('ft3').conversion_factor('gallon')
_m3_to_gal = auom('m3').conversion_factor('gallon')

base_dosing_path = ospath.join(data_path, 'sanunit_data/VFA/_basedosing.xlsx')

class BaseDosing(WWTpump):
    '''
    Perform the duty of a pump and a mixer. Base addition amount is based on
    the titration ratio provided by Su Group which is the ratio of base to
    AD sludge. Design and cost function calculate for the base dosing pump
    and is based on biosteam _pump.py. Material requirement for the base 
    storage tank was included referencing the chemical pump within the qsdsan
    WWTPump class. Tank material was assumed to be HDPE.
    '''
    
    _N_ins = 3
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 base_conc = 1, titr_factor = 0.01225, **kwargs):
        WWTpump.__init__(self, ins=ins, outs=outs, pump_type='chemical', N_pump = 1, include_OM_cost=True)
        
        self.base_conc = base_conc
        self.titr_factor = titr_factor
        # titr_factor calculated based on the experimental data from Wangsuk: 1.98 ml of 1 M NaOH into 160 ml of AD effluent.
        # AD_red only used to calculate base addition and not participate in mass balance
        
        data = load_data(path=base_dosing_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _run(self):
        inf, base, ad_ref = self.ins
        effluent, = self.outs
        
        base_vol = ad_ref.F_vol * self.titr_factor
        base.ivol['H2O'] = base_vol
        base.imol['Na'] = base_vol * self.base_conc * 1000 # kmol/hr
        # by defining the base stream here, no need to initiate the stream in system.py
        effluent = effluent.mix_from((inf, base))
    
    def _init_lca(self):
        D = self.design_results
        self.construction = [Construction("stainless_steel", linked_unit=self,
                                          item = "StainlessSteel",
                                          quantity_unit= "kg"),
                             Construction("electric_motor", linked_unit=self,
                                          item = "ElectricMotor",
                                          quantity_unit= "ea"),
                             Construction("HDPE", linked_unit=self,
                                          item = "Hdpe",
                                          quantity_unit= "ea"),
                             ]
        
    def _design(self):
        constr = self.construction
        # design_func = getattr(WWTPump, f'design_{pump_type}')
        # pipe, pump_ss, hdpe = design_func() # since pump_type is set to 'chemical' in _init_, design_chemical will be called within _design of WWTpump
        pipe, pump_ss, hdpe = super().design_chemical()
        D = self.design_results
        D['Pump pipe stainless steel'] = pipe
        D['Pump stainless steel'] = pump_ss
        D['Storage HDPE'] = constr[2].quantity = hdpe
        D['Pump motor'] = constr[1].quantity = self.N_pump
       
        constr[0].quantity = pipe + pump_ss
            
    def _cost(self):
        D = self.design_results
        super()._cost() #!!! do we assume WWTpump _cost includes the price of entire pump including the motor?
        self.baseline_purchase_costs['Storage HDPE'] = D['Storage HDPE'] * self.hdpe_cost
        #!!! add HDPE cost

    
    @property
    def Q_mgd(self):
        '''
        [float] Volumetric flow rate in million gallon per day, [mgd].
        Will use total volumetric flow through the unit if not provided.
        '''
        return self.F_vol_in*_m3_to_gal*24/1e6
    @Q_mgd.setter
    def Q_mgd(self, i):
        self._Q_mgd = i
        
    @property
    def Q_cfs(self):
        '''[float] Volumetric flow rate in cubic feet per second, [cfs].'''
        return self.Q_mgd*1e6/24/60/60/_ft3_to_gal
        
    @property
    def H_ts(self):
        '''[float] Total static head, [ft].'''
        return self._H_ts
    
    @property
    def H_p(self):
        '''[float] Pressure head, [ft].'''
        return self._H_p
    
    @property
    def v(self):
        '''[float] Fluid velocity, [ft/s].'''
        return self._v
    @v.setter
    def v(self, i):
        self._v = i
        
    @property
    def C(self):
        '''[float] Hazen-Williams coefficient to calculate fluid friction.'''
        return self._C
    @C.setter
    def C(self, i):
        self._C = i
        
    @property
    def SS_per_pump(self):
        '''[float] Quantity of stainless steel per pump, [kg/ea].'''
        return self._SS_per_pump
    @SS_per_pump.setter
    def SS_per_pump(self, i):
        self._SS_per_pump = i
    
        
#%% Separation: redox-ED
redox_ED_path = ospath.join(data_path, 'sanunit_data/VFA/_redox_ED.csv')
_MW_RC = 484.06 + 329.24 # molecular weight of sodium ferrocyanide decahydrate and potassium ferricyanide

class RedoxED(SanUnit):
    '''
    Redox-electrodialysis unit to separate volatile fatty acids from waste.
    
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
    Assume a black box model for the ionic flux, continuously-stirred tank reactors
    for the feeding and accumulating channels. Model membrane transfer as 1st order
    reaction and solve a steady state mass balance for concentration.
    
    Reference:
    Nutrient recovery from wastewater through pilot scale electrodialysis. Ward et al. 
    2018 Water Research. 
    https://www.sciencedirect.com/science/article/pii/S0043135418301209?via%3Dihub
    
    '''
    _N_ins = 2
    _N_outs = 2
    
    # !!!
    # Linear regression coefficients to extrapolate flux, calculated in excel
    # and imported as instance attributes. Linear regression coefficients
    # should be updated when new experimental data is used for future instances
    # !!! Check units for the mass balance and concentrations!
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 voltage=None, op_time=None, m_area = 0.12, cell_pair_number = 30, 
                 **kwargs):
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
        self.m_area = m_area #m2 from Ward et al.
        self.cell_pair_number = cell_pair_number
        #!!! rc_voume is a place holder value now, adjust to overall size and flow rate
        #self.area = self.compartment_length * self.compartment_width
    
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
        #!!! are these two lines needed after copy_like? or should they be moved?
        # the purpose was to keep the volume constant in feeding and accumulating
        # channels.
        
        # Calculate the average flux from cell voltage
        # !!! should calculate the distribution of these linear regression parameters for their uncertainty range?
        # !!! how to reduce repetitiveness? when flux extrapolation constants
        # are imported in other ways instead of as class attributes in future,
        # need to change ways to retrieve this data as well
        c3_flux = (self.c3_slope*self.voltage + self.c3_const) * 10**(-7) # [mol m-2 s-1]
        c4_flux = (self.c4_slope*self.voltage + self.c4_const) * 10**(-7)
        c6_flux = (self.c6_slope*self.voltage + self.c6_const) * 10**(-7)
        self.c3_flux = c3_flux
        self.c4_flux = c4_flux
        self.c6_flux = c6_flux
        
# =============================================================================
#         # Determination of membrane area based on typical recovery of propionate
#         m_area = fc_in.imol['Propionate'] * (1000/3600) / c3_flux * self.c3_recovery #wastestream.imol is in kmol/hr, converted to mol/s
#         
# =============================================================================
        active_area = self.m_area * self.cell_pair_number
        # Calculate accumulating channel effluent concentration
        ac_out.imol['Propionate'] = ac_in.imol['Propionate'] + \
        c3_flux * active_area * 3600 * 10**(-3)
        ac_out.imol['Butyrate'] = ac_in.imol['Butyrate'] + \
        c4_flux * active_area * 3600 * 10**(-3)
        ac_out.imol['Hexanoate'] = ac_in.imol['Hexanoate'] + \
        c6_flux * active_area * 3600 * 10**(-3)
        
        # Calcualte feeding channel effluent concentration        
        fc_out.imol['Propionate'] = fc_in.imol['Propionate'] - \
        c3_flux * active_area * 3600 * 10**(-3)
        fc_out.imol['Butyrate'] = fc_in.imol['Butyrate'] - \
        c4_flux * active_area * 3600 * 10**(-3)
        fc_out.imol['Hexanoate'] = fc_in.imol['Hexanoate'] - \
        c6_flux * active_area * 3600 * 10**(-3)
        
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

        recirc_pump = self.auxiliary('recirc_pump', cls = SludgePump, ins = self.outs[1].copy('accum_channel_out')) # L_s = 1, L_d = 1; in the self.auxiliary method, cannot modify the length of pipes anciliary to the pump
        # assume 1 ft for the suction and discharge pipe
        
    def _init_lca(self):
        self.construction = [Construction('carbon_cloth', linked_unit=self, 
                                           item='CarbonCloth', 
                                           quantity_unit='kg'),
                             Construction('electrode', linked_unit=self,
                                           item='TiElectrode', 
                                           quantity_unit='kg'),
                             Construction('silicone_spacer', linked_unit=self,
                                           item='SiliconeSpacer',
                                           quantity_unit='kg'),
                             Construction('IEM', linked_unit=self,
                                           item='IEM',
                                           quantity_unit='kg'),    #is there a common format for units?
                             Construction('redox_couple', linked_unit=self,
                                          item='RedoxCouple',
                                          quantity_unit='kg'),
                             Construction('piping', linked_unit=self,
                                           item='PvcPiping',
                                           quantity_unit='m'),
                             Construction('housing', linked_unit=self,
                                           item='Housing',
                                           quantity_unit='kg'),]

# !!! change construction items when scaling up system. Now the materials are
# based on bench top experiment.
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        area = self.m_area / self.active_area_ratio
        #dimension = area ** 0.5
        cp_number = self.cell_pair_number
        tot_area = area * cp_number * 2
        design['CarbonCloth'] = constr[0].quantity = tot_area * self.carbon_cloth_thickness/1000  # carbon cloth area, one carbon cloth at each electrode
        design['TiElectrode'] = constr[1].quantity = tot_area * self.ti_thickness/1000 # electrode volume, two electrodes
        design['IEM'] = constr[3].quantity = area * (cp_number * 2 + 1) * 0.1/1000 # IEM area multiplied by thickness (asume thickness to be 0.1 mm), based on Nayeong's paper https://pubs.acs.org/doi/suppl/10.1021/acsenergylett.3c00482/suppl_file/nz3c00482_si_001.pdf
        design['SiliconeSpacer'] = constr[2].quantity = (area - self.m_area) * (cp_number * 6) * self.spacer_thickness/1000 # silicone spacer volume, based on stacked cells configuration, including for rc channel, 6 spacers per cell
        # design['HDPE'] = constr[4].quantity = cp_number * 2 * (self.rc_thickness / area * dimension * 4 + area * 2) * 0.003 # approximate for HDPE volume; assume square shape of membrane and 0.003 thickness of HDPE
        design['RedoxCouple'] = constr[4].quantity = self.rc_thickness /100 * area * self.rc_conc / 1000 * _MW_RC / 1000 # kg, required mass of sodium ferrocyanide decahydrate, based on concentration used by Oh et al. converted
        design['Piping'] = constr[5].quantity = 10 * cp_number # ft, assume a total of 10 ft piping of the same material for the redox-ED unit
        design['AcrylicHousing'] = constr[6].quantity = tot_area * self.acrylic_thickness/100 # in m3, assume 1 cm thickness for acrylic board
            # one cell pair is estimated to be 25 mm based on cell pair configuration and material thickness
            # one electrode is estimated to be 2 cm based on cell pair configuration and material thickness
            # assume square shape of membranes
            # assume stacking ccell pairs which include electrodes instead of just repeating feeding and acucmulating channels with single electrodes at both end of the stacks. This assumption seems inconsistent in Nayeong's papers
        # !!! redox channel volume to be calculated and scaled based on membrane area in the future, rc_chickness stay constant?
        # !!! what is not included is the amount of cleaning chemicals, however data is available for the cost of cleaning chemicals per water treated, so cost of cleaning chemicals accounted in system.py
        #electricity requirement
        specific_energy = self.e_slope * self.voltage # linear regression of spec energy wrp to voltage, forcing intercept to zero
        mass_transferred = sum((self.c3_flux, self.c4_flux, self.c6_flux)) * self.m_area * self.cell_pair_number * 8 * 3600 # mass transferred in a day for 8 hrs of operation per day
        self.power_utility((specific_energy * mass_transferred))


        
        #!!! any storage tank needed?
        #!!! calculate using decorators instead of typing equations here?
    def _cost(self):
        area = self.m_area / self.active_area_ratio
        cp_number = self.cell_pair_number
        tot_area = area * cp_number * 2        
        Cost = self.baseline_purchase_costs
        Design = self.design_results
        Cost['CarbonCloth'] = tot_area * self.carbon_cloth_cost 
        Cost['TiElectrode'] = tot_area * self.ti_electrode_cost 
        Cost['SiliconeSpacer'] = (area - self.m_area) * (cp_number * 6) * self.si_spacer_cost 
        Cost['IEM'] = area * (cp_number * 2 + 1) * self.iem_cost 
        Cost['Piping'] = Design['Piping'] * self.pipe_cost
        Cost['AcrylicHousing'] = tot_area * self.acrylic_cost
        
        tot_equip_and_constr = self.baseline_purchase_cost
        self.add_OPEX = {'RedoxCouple': Design['RedoxCouple'] * self.rc_replacement * self.rc_cost /365, # 2 channels. assume redox couple is replaced every two years
                         'Labor': self.wages * 8, # assume 1 operator needed who work for 8 hours per day, [$/day]
                         'Maintenance': tot_equip_and_constr * self.material_repair_and_replacement_cost / 365}
        #!!! check the O&M assumption from CapdetWorks
        # maintenance OpEx does not include chemical costs Generous et al. 2021 https://www.sciencedirect.com/science/article/pii/S1383586621005864?via%3Dihub
        # by dividing by 365, we assume continuous operation everyday for the entire year
        self.recirc_pump.simulate()        
        
        
        