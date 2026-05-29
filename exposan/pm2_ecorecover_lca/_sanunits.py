#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""
import os, numpy as np
from qsdsan import SanUnit, Construction, WasteStream, System, unit_operations as su
from qsdsan import processes as pc
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
import math

__all__ = ('Photobioreactor'
           )
#%%
CSTR = su.CSTR
lb_to_kg = 0.453592
acre_to_sq_m = 4046.86
sq_feet_to_sq_m = 10.7639
CEPCI_by_year.update({
    2022: 816.0,
    2023: 797.9,
    2024: 798.0,   # https://reg.lub.lu.se/luur/download?func=downloadFile&recordOId=9209157&fileOId=9209158
    2025: 811.0,   # estimated as 2024 * 1.016
})
@cost(basis = 'Aerial footage-volume-to-area ratio', ID='Glass tube & fittings', units='m3/m2',
      cost = 233240/acre_to_sq_m, S=0.029499829299047615, CE=CEPCI_by_year[2014], n=1, BM=1.1) #ref:https://docs.nrel.gov/docs/fy19osti/72716.pdf
@cost(basis='Aerial footage', ID='Greenhouse',units='m2',cost= 13*sq_feet_to_sq_m, S=1, CE=CEPCI_by_year[1994], n=1, BM=1) #conventional greenhouse in ref page 38: https://blog.uvm.edu/cwcallah/files/2021/03/NRAES-33_Web.pdf
@cost(basis='Aerial footage-tube area per support area', ID='Support structure', units='m2',
      cost=118937/acre_to_sq_m, S=0.00983, CE=CEPCI_by_year[2014], n=1, BM=1.1) #ref:https://docs.nrel.gov/docs/fy19osti/72716.pdf
@cost(basis='Aerial footage', ID= 'LED system', units = 'm2', cost= 11.27*sq_feet_to_sq_m, S=1, CE=CEPCI_by_year[2025],n=1, BM=1.1)

class Photobioreactor(CSTR):
    '''
    Helical fence-type tubular photobioreactor.

    Parameters
    ----------
    ID : str
        ID for the reactor.
    ins : :class:`WasteStream`
        Influents to the reactor. Can be an array of up to 3 WasteStream objects by
        default, typically wastewater to be treated, recycled effluent, recycled
        activated sludge.
    outs : :class:`WasteStream`
        Treated effluent.
    split : iterable of float
        Volumetric splits of effluent flows if there are more than one effluent.
        The default is None.
    InD: float
        The designed inner diameter [m] of the tube. The default is 0.1 m.
    OD: float
        The designed outer diameter [m] of the tube. The default is 0.105 m.
    spt: float
        The designed single pass time [hr] for an micralgae cell to travel through one set of photobioreactor. The default is 1.67 hrs.
    hrt: float
        The designed hydraulic retention time [hr] for influent. The default is 4 hrs.
    spacing: float
        The designed average space [m] between each set of photobioreactors. The default is 1.05 m.
    support_spacing: float
        The designed average space [m] between each support structure for glass tubes. The default is 2.0 m.
    length: float
        The designed maximum length [m] for each row of photobioreactor. The default is 80 m.
    v: float
        The designed maximum flow velocity [m/s] in the tube. The default is 0.4 m/s. Note this flow is the combined flow from influent and recycled flow.
    rr: float
        The designed recirculation rate proportional to the influent flow rate. The default is 1.
    light_intensity: int
        The designed light intensity [umol/(m2·s)] over the illuminated area. The default is 100.
    annual_heating_days: int
        The average winter days per year that requires heaters, default to 150 days.
    heater_up_time_ratio: float
        The average daily up time ratio [0-1] of heaters, default to 0.15.
    ===============================
    V_max : float
        Designed volume, in [m^3]. The default is 1000.
        
    # W_tank : float
    #     The design width of the tank, in [m]. The default is 6.4 m (21 ft). [1, Yalin's adaptation of code]
    # D_tank : float
    #     The design depth of the tank in [m]. The default is 3.65 m (12 ft). [1, Yalin's adaptation of code]
    # freeboard : float
    #     Freeboard added to the depth of the reactor tank, [m]. The default is 0.61 m (2 ft). [1, Yalin's adaptation of code]
        
    aeration : float or :class:`Process`, optional
        Aeration setting. Either specify a targeted dissolved oxygen concentration
        in [mg O2/L] or provide a :class:`Process` object to represent aeration,
        or None for no aeration. The default is 2.0.
    DO_ID : str, optional
        The :class:`Component` ID for dissolved oxygen, only relevant when the
        reactor is aerated. The default is 'S_O2'.
    suspended_growth_model : :class:`Processes`, optional
        The suspended growth biokinetic model. The default is None.
    exogenous_var : iterable[:class:`ExogenousDynamicVariable`], optional
        Any exogenous dynamic variables that affect the process mass balance,
        e.g., temperature, sunlight irradiance. Must be independent of state 
        variables of the suspended_growth_model (if has one).
    
    References
    ----------
     [1] #TODO
    
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                split=None,V_max=1000, W_tank = 6.4, D_tank = 3.65,
                freeboard = 0.61, t_wall = None, t_slab = None, aeration=2.0, 
                DO_ID='S_O2', suspended_growth_model=None, 
                gas_stripping=False, gas_IDs=None, stripping_kLa_min=None, 
                K_Henry=None, D_gas=None, p_gas_atm=None,
                isdynamic=True, exogenous_vars=(), 
                InD=0.1, 
                OD =0.105,
                spt=1.67,
                hrt=4,
                spacing=1.05,
                support_spacing=2,
                length = 80,
                v= 0.4,
                rr=1,
                light_intensity = 100,
                annual_heating_days = 150,
                heater_up_time_ratio = 0.15,
                include_construction= True,
                **kwargs):
        CSTR.__init__(self,ID=ID,ins=ins,outs=outs,split=None,V_max=V_max, W_tank = W_tank, D_tank = D_tank,
                freeboard = freeboard, t_wall = t_wall, t_slab = t_slab, aeration=aeration, 
                DO_ID=DO_ID, suspended_growth_model=suspended_growth_model, 
                gas_stripping=gas_stripping, gas_IDs=gas_IDs, stripping_kLa_min=stripping_kLa_min, 
                K_Henry=K_Henry, D_gas=D_gas, p_gas_atm=p_gas_atm,
                isdynamic=isdynamic, exogenous_vars=exogenous_vars, )
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,F_BM_default=None)
        self.InD = InD
        self.OD = OD
        self.spt = spt
        self.hrt = hrt
        self.spacing = spacing
        self.support_spacing = support_spacing
        self.length = length
        self.v = v
        self.rr = rr
        self.light_intensity = light_intensity
        self.annual_heating_days = annual_heating_days
        self.heater_up_time_ratio = heater_up_time_ratio
    
    def _init_lca(self):#TODO
        self.include_construction = True
        self.construction = [
            Construction('glasses', linked_unit=self, 
                         item='XXXX', 
                         quantity_unit='kg'),
            Construction('PVC', linked_unit=self, 
                         item='XXXX', 
                         quantity_unit='kg'),
            Construction('HDPE', linked_unit=self, 
                         item='XXXX', 
                         quantity_unit='kg'),
            Construction("stainless_steel", linked_unit=self,
                         item = "StainlessSteel", 
                         quantity_unit= "kg"),
            Construction("aluminum", linked_unit=self,
                         item = "XXXX", 
                         quantity_unit= "kg"),
            Construction("LED", linked_unit=self,
                         item = "XXXX", 
                         quantity_unit= "kg"),
            Construction("FRP", linked_unit=self,
                         item = "XXXX", 
                         quantity_unit= "kg"),
            Construction("polycarbonate_glass", linked_unit=self,
                         item = "XXXX", 
                         quantity_unit= "kg"),
            Construction('fan',linked_unit=self, 
                         item='Fan',
                         quantity_unit='kg'),
            ]
    def _design(self):
        self.design_results['Total flow'] = self.ins[0].F_vol #m3/hr
        ###Glass tube
        #vworking olume of glass tube
        self.design_results['Working volume'] = self.design_results['Total flow']*self.hrt #m3
        #total length
        self.design_results['Total length'] = self.design_results['Working volume']/math.pi/(self.InD/2)**2 #m
        #the number of sets required
        self.design_results['sets of PBR blocks'] = math.ceil(self.design_results['Total flow']/3600*
                                                              (1+self.rr)/math.pi/(self.InD/2)**2/self.v)
        self.design_results['number of rows per set'] = math.ceil(self.design_results['Total length']/self.design_results['sets of PBR blocks']/self.length)
        self.design_results['Aerial footage'] = self.length*(self.design_results['sets of PBR blocks']+1)*self.spacing #m2
        self.design_results['volume-to-area ratio'] = self.design_results['Working volume'] / self.design_results['Aerial footage'] #m3/m2
        self.design_results['Aerial footage-volume-to-area ratio'] = self.design_results['volume-to-area ratio']*self.design_results['Aerial footage']
        #total glass material weight
        rho_glass = 2230 #kg/m3, reference: https://www.imetra.com/borosilicate-glass-material-properties/
        self.design_results['Total glass weight'] = (
        math.pi * ((self.OD / 2) ** 2 - (self.InD / 2) ** 2)
        * self.design_results['Total length'])*rho_glass #kg
        self.construction[0].quantity = self.design_results['Total glass weight']
        self.design_results['U-bends'] = self.design_results['Total length']/self.length*2 #ea
        self.design_results['Plates'] = self.design_results['anchor bolts'] = self.design_results['U-bends']/2 #ea,based on Clearas estimation
        self.design_results['Plate mounting bolts'] =  self.design_results['nuts'] = self.design_results['washers']=self.design_results['plates']*30 #ea,based on Clearas estimation
        self.design_results['Coupling'] = self.design_results['Total length']/3.1 #ea,based on Clearas estimation
        self.design_results['Hose clamps'] = self.design_results['Total length']/0.76 #ea,based on Clearas estimation
        
        self.construction[1].quantity = (self.design_results['U-bends']*1.4+ #U-bend, PVC
                                         self.design_results['coupling']*0.91 #coupling, PVC
                                         )
        self.construction[4].quantity = self.design_results['Plates'] * 0.2 # aluminum plates
        self.construction[3].quantity = (self.design_results['Plate mounting bolts'] * 0.03 + 
                                         self.design_results['Hose clamps'] * 0.06 + #based on Clearas estimation
                                         self.design_results['anchor bolts'] * 0.7 +  #based on Clearas estimation
                                         self.design_results['nuts'] * 0.02 +
                                         self.design_results['washers'] * 0.005
                                         )
        self.construction[2].quantity = 424*lb_to_kg/79*self.design_results['sets of PBR blocks'] #kg HDPE pipe
        
        ###Support structure
        self.design_results['Support structure'] = math.ceil((self.length/self.support_spacing+1))*self.design_results['sets of PBR blocks'] #ea
        self.construction[3].quantity += self.design_results['Support structure']*50 #assume 50 kg stainless steel per support structure
        self.design_results['Aerial footage-tube area per support area'] = (self.design_results['Aerial footage']*
                                                                        (math.pi*(self.InD/2)**2*
                                                                     self.design_results['number of rows per set'] 
                                                                     /self.support_spacing/self.spacing)) #m2
        
        ###Greenhouse
        self.design_results['Greenhouse panel area'] = 1.02*self.design_results['Aerial footage']+15.5*self.design_results['Aerial footage']**0.5 #m2, gable roof with 3m sidewall and 4m total height, aerial footage L:W= 4 
        panel_thickness = 0.008 #m,  http://www.unitedgreenhouse.com/accessories/multilayered-poly-coverings.php
        polycarbonate_density = 1200 #kg/m3
        self.construction[7].quantity = self.design_results['Greenhouse panel area']*panel_thickness*polycarbonate_density #kg
        self.construction[8].quantity = self.design_results['Fans weight']= self.design_results['Aerial footage']/504*722.98*lb_to_kg #fan in kg based on clearas estimation, linear scale up
        self.design_results['250K BTU heater'] = math.ceil(self.design_results['Aerial footage']/6720*20) #ea, Modine HD125 model, number based on Clearas estimation
        self.construction[3].quantity += self.design_results['250K BTU heater']*143* lb_to_kg #kg stainless steel
        self.construction[6].quantity = self.design_results['Aerial footage']/6720*294*lb_to_kg #kg FRP
        
        #TODO
        #need to convert natural gas into a stream
        self.design_results['natural gas'] = (self.design_results['250K BTU heater']*
                                              self.heater_up_time_ratio*
                                              self.annual_heating_days*2.5*1.9) #1.9 kg natural gas per therm, 2.5 therm per hr fuel input
        
        ###LED lighting
        self.design_results['LED light'] = self.light_intensity*self.design_results['Aerial footage']/1100 #ea, 1100 PPF per light based on Clearas estimation
        self.self.construction[4].quantity += self.design_results['LED light'] * 10 * lb_to_kg #kg, aluminum as the main component, https://californialightworks.com/megadrive/megadrive-linear-400/
        self.design_results['LED mounting structure'] = 32.3*self.design_results['LED light'] #kg, 32.3 kg each based on Clearas estimation
        self.construction[3].quantity += self.design_results['LED mounting structure'] #kg, stainless steel
        self.construction[5].quantity = self.design_results['LED light'] * 0.012 #kg LED, 12 g LED per fixture
        
        ###Pigging & interconnects
        self.design_results['PIG assemblies PVC'] = 19907/6720*self.design_results['Aerial footage']*lb_to_kg #kg of PVC
        self.construction[1].quantity += self.design_results['PIG assemblies PVC']
        self.design_results['PIG assemblies stainless steel'] =1934/6720*self.design_results['Aerial footage']*lb_to_kg #kg of stainless steel
        self.construction[3].quantity += self.design_results['PIG assemblies stainless steel']
        
        ###Pumps
        #TODO
        
        
        self._init_lca()
        inf, RAA=self.ins
        eff, = self.outs
        