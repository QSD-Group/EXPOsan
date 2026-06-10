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
from qsdsan.unit_operations import WWTpump
from qsdsan.unit_operations.bst._pumping import Pump
from qsdsan import processes as pc
from qsdsan.utils import auom, select_pipe, format_str
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from warnings import warn

import math
from math import pi, ceil
__all__ = ('Photobioreactor'
           )
#%%
CSTR = su.CSTR
lb_to_kg = 0.453592
acre_to_sq_m = 4046.86
sq_feet_to_sq_m = 10.7639
m_to_feet = 3.28084
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
    rho: float
        The density of the mixed liquor [kg/m3]. The default is 1000.
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
                rho = 1000,
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
        self.rho = rho
    
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

        self._init_lca()
        inf, RAA=self.ins
        eff, = self.outs

#%%
class Ecorecoverypump(WWTpump):
    '''
    Pump subclass for Ecorecovery systems.

    Adds a feed_PBR pump type for the hydraulic design of fence-type tubular photobioreactors.

    Additional inputs for pump_type='feed_PBR'
    ------------------------------------------------
    InD : float
        Inner diameter of the PBR tube, [ft].
        Default is 0.1 m converted to ft.

    OD : float
        Outer diameter, or vertical spacing basis, of the PBR tube, [ft].

    hrt : float
        Hydraulic retention time based on influent flow, [hr].
        Default is 4 hr.

    length : float
        Maximum length of each PBR row, [ft].
        Default is 80 m converted to ft.

    v : float
        Designed flow velocity inside the PBR tube, [ft/s].
        This should represent the combined influent plus recirculation flow.
        Default is 0.4 m/s converted to ft/s.

    rr : float
        Recirculation ratio relative to influent flow, dimensionless.
        Default is 1.

    H_p : float
        Pressure head, [ft].
        Default is 0.

    L_s : float
        Suction pipe length, [ft].
        Default is 10 m converted to ft.

    N_ubends : int or None
        Number of U-bends per PBR set. If None, it is estimated as
        number of rows per set minus 1.
    '''
    _valid_pump_types = WWTpump._valid_pump_types + ('feed_PBR',)
    _ft_to_m = auom('ft').conversion_factor('m')
    _m_to_ft = 1 / _ft_to_m
    _g = 32.174  # gravitational acceleration, [ft/s2]
    F_BM_pump = 1.18*(1+0.007/100)
    _lb_to_kg = auom('lb').conversion_factor('kg')
    default_F_BM = {
            'Pump': F_BM_pump,
            'Pump building': F_BM_pump,
            }
    default_equipment_lifetime = {
        'Pump': 15,
        'Pump pipe stainless steel': 15,
        'Pump stainless steel': 15,
        'Pump chemical storage HDPE': 30,
        }
    _feed_PBR_input_order = (
        'InD',      # [ft]
        'OD',       # [ft]
        'hrt',      # [hr]
        'length',   # [ft]
        'v',        # [ft/s]
        'rr',       # dimensionless
        'H_p',      # [ft]
        'L_s',      # [ft]
        'N_ubends', # dimensionless, optional
    )

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 prefix='', pump_type='', Q_mgd=None, add_inputs=(),
                 InD=None, OD=None, hrt=4.,
                 length=None, v=None, rr=1.,
                 H_p=0., L_s=None, N_ubends=None,
                 capacity_factor=1.,
                 include_pump_cost=True, include_building_cost=False,
                 include_OM_cost=False,
                 F_BM=default_F_BM,
                 lifetime=default_equipment_lifetime,
                 **kwargs):

        super().__init__(
            ID=ID, ins=ins, outs=outs, thermo=thermo,
            init_with=init_with,
            prefix=prefix,
            pump_type=pump_type,
            Q_mgd=Q_mgd,
            add_inputs=add_inputs,
            capacity_factor=capacity_factor,
            include_pump_cost=include_pump_cost,
            include_building_cost=include_building_cost,
            include_OM_cost=include_OM_cost,
            F_BM=F_BM,
            lifetime=lifetime,
            **kwargs
        )

        # Store feed_PBR specific values separately so other WWTpump types
        # can still use WWTpump's default v, H_p, etc.
        self.feed_PBR_InD = 0.1 * self._m_to_ft if InD is None else InD
        self.feed_PBR_OD = 0.105 * self._m_to_ft if OD is None else OD
        self.feed_PBR_hrt = hrt
        self.feed_PBR_length = 80 * self._m_to_ft if length is None else length
        self.feed_PBR_v = 0.4 * self._m_to_ft if v is None else v
        self.feed_PBR_rr = rr
        self.feed_PBR_H_p = H_p
        self.feed_PBR_L_s = 10 * self._m_to_ft if L_s is None else L_s
        self.feed_PBR_N_ubends = N_ubends

    @property
    def pump_type(self):
        return self._pump_type

    @pump_type.setter
    def pump_type(self, i):
        i = i or ''
        i_lower = i.lower()
        i_lower = i_lower.replace('cstr', 'CSTR')
        i_lower = i_lower.replace('af', 'AF')
        i_lower = i_lower.replace('pbr', 'PBR')

        if i_lower not in self.valid_pump_types:
            raise ValueError(
                f'The given `pump_type` "{i}" is not valid, '
                'check `valid_pump_types` for acceptable pump types.'
            )

        if i_lower == 'feed_PBR':
            warn(
                "`pump_type='feed_PBR'` requires PBR-specific inputs. "
                "Pass them as keyword arguments or through `add_inputs` in this order: "
                "InD [ft], OD [ft], hrt [hr], length [ft], v [ft/s], "
                "rr [-], H_p [ft], L_s [ft], N_ubends [-]. ",
                RuntimeWarning,
                stacklevel=3,
            )

        self._pump_type = i_lower

    def _get_feed_PBR_inputs(self, InD=None, OD=None, hrt=None,
                             length=None, v=None, rr=None,
                             H_p=None, L_s=None, N_ubends=None):
        values = {
            'InD': self.feed_PBR_InD,
            'OD': self.feed_PBR_OD,
            'hrt': self.feed_PBR_hrt,
            'length': self.feed_PBR_length,
            'v': self.feed_PBR_v,
            'rr': self.feed_PBR_rr,
            'H_p': self.feed_PBR_H_p,
            'L_s': self.feed_PBR_L_s,
            'N_ubends': self.feed_PBR_N_ubends,
        }

        # Allow add_inputs to follow the WWTpump style.
        for name, value in zip(self._feed_PBR_input_order, self.add_inputs):
            if value is not None:
                values[name] = value

        explicit_values = {
            'InD': InD,
            'OD': OD,
            'hrt': hrt,
            'length': length,
            'v': v,
            'rr': rr,
            'H_p': H_p,
            'L_s': L_s,
            'N_ubends': N_ubends,
        }

        for name, value in explicit_values.items():
            if value is not None:
                values[name] = value

        missing = [
            name for name in ('InD', 'OD', 'hrt', 'length', 'v', 'rr', 'H_p', 'L_s')
            if values[name] is None
        ]
        if missing:
            raise ValueError(
                "`pump_type='feed_PBR'` is missing required input(s): "
                f"{', '.join(missing)}."
            )

        if values['InD'] <= 0:
            raise ValueError('`InD` must be positive, [ft].')
        if values['OD'] <= 0:
            raise ValueError('`OD` must be positive, [ft].')
        if values['hrt'] <= 0:
            raise ValueError('`hrt` must be positive, [hr].')
        if values['length'] <= 0:
            raise ValueError('`length` must be positive, [ft].')
        if values['v'] <= 0:
            raise ValueError('`v` must be positive, [ft/s].')
        if values['rr'] < 0:
            raise ValueError('`rr` must be non-negative.')
        if values['L_s'] < 0:
            raise ValueError('`L_s` must be non-negative, [ft].')

        return values

    def design_feed_PBR(self, Q_mgd=None, InD=None, OD=None, hrt=None,
                        length=None, v=None, rr=None,
                        H_p=None, L_s=None, N_ubends=None,
                        **kwargs):
        '''
        Design the feed pump for a fence-type tubular PBR.

        This method intentionally does not call WWTpump._design_generic().
        It follows the PBRpump hydraulic design structure, but all length,
        volume, and velocity calculations are in US units.
        '''

        # Allow additional overrides through kwargs.
        if kwargs:
            for k, val in kwargs.items():
                if k in self._feed_PBR_input_order:
                    locals()[k] = val
                else:
                    setattr(self, k, val)

        if Q_mgd is None:
            Q_mgd = self.Q_mgd
        self.Q_mgd = Q_mgd

        vals = self._get_feed_PBR_inputs(
            InD=InD, OD=OD, hrt=hrt, length=length,
            v=v, rr=rr, H_p=H_p, L_s=L_s, N_ubends=N_ubends,
        )

        InD = vals['InD']          # [ft]
        OD = vals['OD']            # [ft]
        hrt = vals['hrt']          # [hr]
        length = vals['length']    # [ft]
        v = vals['v']              # [ft/s]
        rr = vals['rr']            # [-]
        H_p = vals['H_p']          # [ft]
        L_s = vals['L_s']          # [ft]
        N_ubends = vals['N_ubends']

        D = self.design_results

        # Flow and PBR sizing, all in US units.
        Q_cfs = self.Q_cfs                         # [ft3/s], influent flow
        Q_cfh = Q_cfs * 3600                       # [ft3/hr]
        tube_area = pi / 4 * InD**2                # [ft2]

        if Q_cfs <= 0:
            self.N_pump = 0
            self._H_ts = self._H_sf = self._H_df = self._H_p = 0.
            D['Total flow'] = 0.
            D['Working volume'] = 0.
            D['Total length'] = 0.
            D['sets of PBR blocks'] = 0
            D['number of rows per set'] = 0
            D['U-bends'] = 0
            return 0., 0., 0.


        # Working volume = influent flow * HRT
        # Total length = working volume / tube cross sectional area
        working_volume = Q_cfh * hrt               # [ft3]
        total_length = working_volume / tube_area  # [ft]

        # Number of PBR sets, equivalent to N_pump in WWTpump logic.
        # Uses combined influent plus recirculation flow.
        N_pump = ceil(Q_cfs * (1 + rr) / tube_area / v)
        N_pump = max(N_pump, 1)
        self.N_pump = N_pump

        rows_per_set = ceil(total_length / N_pump / length)
        rows_per_set = max(rows_per_set, 1)

        if N_ubends is None:
            N_ubends = max(rows_per_set, 0) * 2
        else:
            N_ubends = ceil(N_ubends)

        D['Total flow'] = Q_cfh                         # [ft3/hr]
        D['Working volume'] = working_volume            # [ft3]
        D['Total length'] = total_length                # [ft]
        D['sets of PBR blocks'] = N_pump                # [-]
        D['number of rows per set'] = rows_per_set      # [-]
        D['U-bends'] = N_ubends                         # [-]

        # Store units for the additional feed_PBR design results.
        self._units['Total flow'] = 'ft3/hr'
        self._units['Working volume'] = 'ft3'
        self._units['Total length'] = 'ft'
        self._units['sets of PBR blocks'] = ''
        self._units['number of rows per set'] = ''
        self._units['U-bends'] = ''

        self._v = v
        C = self.C

        # Suction side pipe sizing.
        OD_s, t_s, ID_s = select_pipe(Q_cfs / N_pump, v)  # [in]

        # Static head.
        # Height simplifies to PBR vertical height

        self._H_ts = rows_per_set * OD  # [ft]

        # Suction friction head, Hazen-Williams form.
        self._H_sf = (
            3.02 * L_s * v**1.85 * C**(-1.85) * InD**(-1.17)
        )
        self._H_sf *= self.headloss_multiplication_factor

        # Discharge friction head along one PBR set.
        length_per_set = total_length / N_pump  # [ft]

        # Minor loss from U-bends.
        #https://www.engineeringtoolbox.com/minor-loss-coefficients-pipes-d_626.html
        H_bends = N_ubends * 0.2 * v**2 / (2 * self._g)  # [ft]

        self._H_df = (
            3.02 * length_per_set * v**1.85 * C**(-1.85) * InD**(-1.17)
            + H_bends
        )
        self._H_df *= self.headloss_multiplication_factor

        self._H_p = H_p

        # Stainless steel pipe mass.
        # only counts suction-side pipe SS.
        V_s = N_pump * pi / 4 * (OD_s**2 - ID_s**2) * (L_s * 12)  # [in3]
        M_SS_pipe = 0.29 * V_s * _lb_to_kg          # [kg]

        # Stainless steel pump mass follows WWTpump's N_pump handling.
        M_SS_pump = N_pump * self.SS_per_pump                       # [kg]

        return M_SS_pipe, M_SS_pump, 0.
        