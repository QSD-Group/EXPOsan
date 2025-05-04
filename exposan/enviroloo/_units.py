#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''

This module is developed by:
    Yuyao Huang <yuyaoh20@gmail.com>
    Siqi Tang <siqit@outlook.com>
for Enviroloo (EL) Clear Toilet system

'''
# %%
import numpy as np
from math import ceil
from warnings import warn
from qsdsan import WasteStream
from qsdsan import SanUnit, Construction, Process
from qsdsan.sanunits import IdealClarifier, Mixer, dydt_cstr, _non_reactive, Copier
from biosteam.units import Mixer as BSTMixer
from qsdsan.sanunits._tank import StorageTank
from qsdsan.processes._decay import Decay
#from qsdsan.utils import ospath, load_data, data_path, price_ratio
from qsdsan.utils import load_data, price_ratio
import os
from exposan.utils import _init_modules
from qsdsan.sanunits._toilet import Toilet, MURT
from qsdsan.sanunits._excretion import ExcretionmASM2d
from qsdsan.sanunits._suspended_growth_bioreactor import CSTR
from qsdsan.sanunits._membrane_bioreactor import CompletelyMixedMBR, dydt_mbr
from qsdsan import Processes, CompiledProcesses
from qsdsan.utils import auom

from warnings import warn
from math import pi
from biosteam.units import StorageTank as BSTStorageTank
from biosteam.exceptions import bounds_warning, DesignWarning
from biosteam.units.design_tools import flash_vessel_design
from biosteam.units.design_tools.specification_factors import material_densities_lb_per_ft3
from collections.abc import Iterable
from biosteam._graphics import UnitGraphics
from exposan import enviroloo as el

# %% This callable file will be reposited to qsdsan.SanUnit subbranch with the name of _enviroloo
__all__ = (
    'EL_CT', # Collection tank
    'EL_PC', # Primary clarifier
    'EL_Anoxic', # Anoxic tank
    'EL_Aerobic', # Aerobic tank
    'EL_CMMBR', # Membrane filter
    'EL_CWT', # Clear water tank
    'EL_System', # System-level summary
    )
el_path = os.path.dirname(__file__)
module = os.path.split(el_path)[-1]
data_path, results_path = _init_modules(module, include_data_path = True)
EL_su_data_path = os.path.join(data_path, 'units_data')


# %%

CollectionTank_path = os.path.join(EL_su_data_path, '_EL_CT.tsv')

@price_ratio()
class EL_CT(CSTR):
    
    '''
    Name
    ----
    Collection tank in the Enviroloo (EL) Clear Toilet system.
    
    Parameters
    ----------
    Ins: 
    (1) Mixed wastewater
    (2) Primary clarifier effluent spill
    (3) Clear water tank effluent spill 
    (4) Primary clarifier sludge return

    Outs:
    (1) Treated water
    (2) Methane (CH4)
    (3) Nitrous oxide (N2O)

    Attributes
    ----------
    length_to_diameter : float
        Ratio of the tank length to diameter.
    vessel_material : str
        Material used for constructing the vessel.
    sludge_moisture_content : float
        Moisture content of the sludge (mass of water/total mass).
    COD_removal : float
        Fraction of COD removed in the collection tank.

    References
    ----------
    Similar to the :class:`biosteam.units.MixTank`, but can calculate material usage.

    See Also
    ----------
    class:`biosteam.units.StorageTank`
    '''
    
    '''
    Similar to :class:`biosteam.units.Mixer`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`,
    and allows dynamic simulation.

    See Also
    --------
    `biosteam.units.Mixer <https://biosteam.readthedocs.io/en/latest/units/mixing.html>`_
    '''
    _N_ins = 1 # treated water, PAC, blower
    _N_outs = 1  # treated water, CH4, N2O
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    exponent_scale = 0.1
    ppl = el.ppl
    baseline_ppl = el.baseline_ppl
    
    
    _D_O2 = 2.10e-9   # m2/s

    def __init__(self, 
                 ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=10, W_tank = None, 
                 # D_tank = 3.65,
                 # freeboard = 0.61, 
                 t_wall = None, t_slab = None, aeration=None, 
                 DO_ID='S_O2', suspended_growth_model=None, 
                 gas_stripping=False, gas_IDs=None, stripping_kLa_min=None, 
                 K_Henry=None, D_gas=None, p_gas_atm=None, 
                 tank_steel_volume=19.5, steel_density= 7850, ppl=ppl, baseline_ppl = baseline_ppl,
                 isdynamic=True, exogenous_vars=(), **kwargs):
        super().__init__(
            ID=ID, ins=ins, outs=outs, split=split, thermo=thermo,
            init_with=init_with, V_max=V_max, W_tank = W_tank, 
            # D_tank = D_tank,
            # freeboard = freeboard, 
            t_wall = t_wall, t_slab = t_slab, aeration=aeration, 
            DO_ID=DO_ID, suspended_growth_model=suspended_growth_model, 
            gas_stripping=gas_stripping, gas_IDs=gas_IDs, 
            stripping_kLa_min=stripping_kLa_min, 
            K_Henry=K_Henry, D_gas=D_gas, p_gas_atm=p_gas_atm, tank_steel_volume= tank_steel_volume, steel_density= steel_density,
            ppl=ppl, baseline_ppl = baseline_ppl,
            isdynamic=isdynamic, exogenous_vars=exogenous_vars, **kwargs
            )
        
        self.tank_steel_volume=tank_steel_volume
        self.steel_density=steel_density
        self.ppl=ppl
        self.baseline_ppl=baseline_ppl
        
        # # Design parameters 
        # self._W_tank = W_tank
        # self._D_tank = D_tank
        # self._freeboard = freeboard
        # self._t_wall = t_wall
        # self._t_slab = t_slab

 
    
        data = load_data(path = CollectionTank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)    
        ###############################################

    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),
                             Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_weight # * (self.ppl / self.baseline_ppl)
        design['Cast_iron'] = constr[1].quantity = self.lift_pump_cast_iron
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.collection_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
        C['Pumps'] = self.lift_pump_cost
    
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_CT
        self.power_utility(power_demand)
    
    def _calc_replacement_cost(self):
        scale  = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        CT_replacement_cost = (
            self.collection_tank_cost / self.collection_tank_lifetime +               
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
            self.lift_pump_cost / self.pump_lifetime) * scale
        CT_replacement_cost = CT_replacement_cost / (365 * 24)  # convert to USD/hr
        return CT_replacement_cost

# %%


PrimaryClarifier_path = os.path.join(EL_su_data_path, '_EL_PC.tsv')

@price_ratio()
class EL_PC(IdealClarifier):
    """
    Primary clarifier in the Enviroloo (EL) Clear Toilet system for COD and suspended solids removal.

    Parameters
    ----------
    ID : str
        Unique identifier for the unit.
    ins : tuple
        Input streams: (0) Wastewater from lift pump, (1) Nitrate return from membrane tank.
    outs : tuple
        Output streams: (0) Treated water, (1) Spill return, (2) Sludge return.
    sludge_flow_rate : float
        Sludge flow rate (m³/d). Required for sludge return calculation.
    solids_removal_efficiency : float
        Fraction of suspended solids removed (0 to 1). Default is 0.5.
    max_overflow : float
        Maximum allowable overflow rate (m³/h). Default is 15.
    ppl : float
        Current population served.
    baseline_ppl : float
        Baseline population for scaling design and cost.

    Notes
    -----
    - COD and suspended solids are tracked via the `components` object.
    - Spill return occurs if treated water flow exceeds `max_overflow`.
    - Inherits from `SanUnit`, not `IdealClarifier`, for flexibility.
    """

    _N_ins = 2
    _N_outs = 2  # [0] effluent overflow, [1] sludge underflow
    _outs_size_is_fixed = True
    exponent_scale = 0.1
    ppl=1000,
    baseline_ppl= 1000 

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 sludge_flow_rate=10, 
                 solids_removal_efficiency=.85,
                 sludge_MLSS=None, isdynamic=False, init_with='WasteStream', tank_steel_volume=3.34, steel_density=7850,
                 ppl=ppl, baseline_ppl=baseline_ppl,
                 F_BM_default=None, **kwargs):

        super().__init__(
            ID=ID, ins=ins, outs=outs, thermo=thermo,
            sludge_flow_rate=sludge_flow_rate, 
            solids_removal_efficiency=solids_removal_efficiency,
            sludge_MLSS=sludge_MLSS, isdynamic=isdynamic, init_with=init_with, tank_steel_volume=tank_steel_volume, steel_density=steel_density,
            ppl=ppl, baseline_ppl=baseline_ppl,
            F_BM_default=F_BM_default, **kwargs
            )
        
        self.tank_steel_volume=tank_steel_volume
        self.steel_density=steel_density
        self.ppl=ppl
        self.baseline_ppl=baseline_ppl
        
        data = load_data(path = PrimaryClarifier_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)    
            
    def _init_state(self):
        inf = self._mixed
        C_in = inf.conc
        Q_in = inf.F_vol * 24
        self._state = np.append(C_in, Q_in)
        self._dstate = self._state * 0.
        
    def _update_state(self):
        arr = self._state
        Cs = arr[:-1]
        Qi = arr[-1]
        Qs, e_rmv, mlss = self._Qs, self._e_rmv, self._MLSS
        x = self.components.x
        i_tss = x * self.components.i_mass

        of, uf = self.outs
        if uf.state is None: uf.state = np.zeros(len(x)+1)
        if of.state is None: of.state = np.zeros(len(x)+1)

        if Qs:
            Qe = Qi - Qs
            if e_rmv:
                fuf = e_rmv * Qi/Qs + (1-e_rmv)
                fof = 1-e_rmv
            elif mlss:
                tss_in = sum(Cs * i_tss)
                tss_e = (Qi * tss_in - Qs * mlss)/Qe
                fuf = mlss/tss_in
                fof = tss_e/tss_in
        elif e_rmv and mlss:
            tss_in = sum(Cs * i_tss)
            Qs = Qi * e_rmv / (mlss/tss_in - (1-e_rmv))
            Qe = Qi - Qs
            fuf = mlss/tss_in
            fof = 1-e_rmv
        else:
            raise RuntimeError('missing parameter')
            
        if Qs >= Qi: 
            uf.state[:] = arr
            of.state[:] = 0.
        else:
            self._f_uf = fuf
            self._f_of = fof
            uf.state[:-1] = Cs * ((1-x) + x*fuf)
            uf.state[-1] = Qs
            of.state[:-1] = Cs * ((1-x) + x*fof)
            of.state[-1] = Qe

    def _update_dstate(self):
        of, uf = self.outs
        x = self.components.x
        if uf.dstate is None: uf.dstate = np.zeros(len(x)+1)
        if of.dstate is None: of.dstate = np.zeros(len(x)+1)
    
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):        
        _state = self._state
        # _dstate = self._dstate
        _update_state = self._update_state
        # _update_dstate = self._update_dstate
        def yt(t, QC_ins, dQC_ins):
            Q_ins = QC_ins[:, -1]
            C_ins = QC_ins[:, :-1]
            # dQ_ins = dQC_ins[:, -1]
            # dC_ins = dQC_ins[:, :-1]
            Q = Q_ins.sum()
            C = Q_ins @ C_ins / Q
            _state[-1] = Q
            _state[:-1] = C
            # Q_dot = dQ_ins.sum()
            # C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
            # _dstate[-1] = Q_dot
            # _dstate[:-1] = C_dot
            _update_state()
            # _update_dstate()
        self._AE = yt
        ###############################################
        
        # self.sludge_flow_rate = sludge_flow_rate
        # self.solids_removal_efficiency = solids_removal_efficiency
        # self.sludge_MLSS = sludge_MLSS
        # self._mixed = WasteStream()
        # self._f_uf = None
        # self._f_of = None
    
    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),
                             Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]


    def _design(self):
        """Calculate design parameters."""
        
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_weight  + self.primary_return_pump_stainless_steel # * (self.ppl / self.baseline_ppl)  # assume linear scale
        design['Cast_iron'] = constr[1].quantity = self.primary_return_pump_cast_iron
        self.add_construction(add_cost=False)

    def _cost(self):
        """Calculate capital and operating costs."""
        C = self.baseline_purchase_costs
        C['Tank'] = self.PC_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
        C['Pumps'] = self.primary_return_pump_cost

        ratio = self.price_ratio  # Assume price_ratio decorator sets this
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()
        power_demand = self.power_demand_PC  # Default to 0 if not set
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        """Calculate replacement cost in USD/hr."""
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        replacement_cost = (
            self.PC_tank_cost / self.PC_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
            self.primary_return_pump_cost / self.pump_lifetime) * scale
        return replacement_cost / (365 * 24)  # Convert to USD/hr



# %%


Anoxic_path = os.path.join(EL_su_data_path, '_EL_Anoxic.tsv')
@price_ratio()

class EL_Anoxic(CSTR):
    '''
    A single continuous stirred tank reactor.

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
     [1] Shoener, B. D.; Zhong, C.; Greiner, A. D.; Khunjar, W. O.; Hong, P.-Y.; Guest, J. S.
         Design of Anaerobic Membrane Bioreactors for the Valorization
         of Dilute Organic Carbon Waste Streams.
         Energy Environ. Sci. 2016, 9 (3), 1102–1112.
         https://doi.org/10.1039/C5EE03715H.
    
    '''
    _N_ins = 3
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    exponent_scale=0.1
    ppl = el.ppl
    baseline_ppl = el.baseline_ppl
    dosing_flow = el.dosing_flow
    scale_factor = el.scale_factor  
    
    # _D_O2 = 2.10e-9   # m2/s
    #TODO: check mixing value 
    
    def __init__(self, ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=7.3, W_tank = 2.09, mixing_ratio = .035,#scale_factor = 1, dosing_flow = 1,
                 # D_tank = 3.65,
                 # freeboard = 0.61, 
                 t_wall = None, t_slab = None, aeration=None, 
                 DO_ID=None, suspended_growth_model=None, 
                 gas_stripping=False, gas_IDs=None, stripping_kLa_min=None, 
                 K_Henry=None, D_gas=None, p_gas_atm=None,
                 isdynamic=True, exogenous_vars=(), **kwargs):
        super().__init__(
            ID=ID, ins=ins, outs=outs, split=split, thermo=thermo,
            init_with=init_with, V_max=V_max, W_tank = W_tank, 
            # D_tank = D_tank,
            # freeboard = freeboard, 
            t_wall = t_wall, t_slab = t_slab, aeration=aeration, 
            DO_ID=DO_ID, suspended_growth_model=suspended_growth_model, 
            gas_stripping=gas_stripping, gas_IDs=gas_IDs, 
            stripping_kLa_min=stripping_kLa_min, 
            K_Henry=K_Henry, D_gas=D_gas, p_gas_atm=p_gas_atm,
            isdynamic=isdynamic, exogenous_vars=exogenous_vars, **kwargs
            )
        
    
        # # Design parameters 
        self._W_tank = W_tank
        self.mixing_ratio = mixing_ratio
        # self._D_tank = D_tank
        # self._freeboard = freeboard
        # self._t_wall = t_wall
        # self._t_slab = t_slab
    
        data = load_data(path = Anoxic_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)    
        ###############################################
    
    def _run(self):
        glucose = self.ins[2]
        glucose.imass['S_F'] = self.chemical_glucose_mixing_ratio * self.dosing_flow * self.scale_factor * 1.067 # kg/L and L/h    to kg/h TODO: add source for 1.067 kg COD/kg glucose
        
        super()._run()
            
    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),
                             Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),
                             Construction(item='PVC_generic', linked_unit=self, quantity_unit='kg'),]
        

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_weight # * (self.ppl / self.baseline_ppl)  # assume linear scale
        design['Cast_iron'] = constr[1].quantity = self.anoxic_mixer_cast_iron + self.glucose_mixer_cast_iron
        design['PVC_generic'] = constr[2].quantity = self.dosing_pump_PVC
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        # massflow_anoxic = self.ins[0].mass
        C['Tank'] = self.anoxic_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
        # C['Chemcial_glucose'] = self.chemical_glucose_dosage * massflow_anoxic * self.chemical_glucose_price  # make sense the unit of treated water flow
        # Glucose cost is already accounted for in the WasteStream
        C['Mixers'] = self.anoxic_mixer_cost + self.glucose_mixer_cost
        C['Pumps'] = self.dosing_pump_cost
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = (self._calc_replacement_cost() + self.chemical_glucose_mixing_ratio * self.dosing_flow * self.scale_factor *
            self.chemical_glucose_price)
        
        power_demand = self.power_demand_AnoxicTank
        self.power_utility(power_demand)
    
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        Anoxic_tank_replacement_cost = (self.anoxic_tank_cost /self.anoxic_tank_lifetime +
                                        self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
                                        self.pipeline_connectors / self.pipeline_connectors_lifetime +
                                        self.anoxic_mixer_cost / self.anoxic_mixer_lifetime + 
                                        self.glucose_mixer_cost / self.glucose_mixer_lifetime + 
                                        self.dosing_pump_cost / self.pump_lifetime) * scale
        Anoxic_tank_replacement_cost = Anoxic_tank_replacement_cost / (365 * 24)  # convert to USD/hr
        return Anoxic_tank_replacement_cost

# %%

Aerobic_path = os.path.join(EL_su_data_path, '_EL_Aerobic.tsv')
@price_ratio()

class EL_Aerobic(CSTR):
    '''
    A single continuous stirred tank reactor.

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
     [1] Shoener, B. D.; Zhong, C.; Greiner, A. D.; Khunjar, W. O.; Hong, P.-Y.; Guest, J. S.
         Design of Anaerobic Membrane Bioreactors for the Valorization
         of Dilute Organic Carbon Waste Streams.
         Energy Environ. Sci. 2016, 9 (3), 1102–1112.
         https://doi.org/10.1039/C5EE03715H.
    
    '''
    _N_ins = 2 # treated water, PAC, blower
    _N_outs = 3  # treated water, CH4, N2O
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    exponent_scale = 0.1
    ppl = el.ppl
    baseline_ppl = el.baseline_ppl
    dosing_flow = el.dosing_flow
    scale_factor = el.scale_factor  
    
    _D_O2 = 2.10e-9   # m2/s

    def __init__(self, 
                 ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=7.3, W_tank = 2.09, mixing_ratio = .035, #scale_factor = 1, dosing_flow = 1,
                 # D_tank = 3.65,
                 # freeboard = 0.61, 
                 t_wall = None, t_slab = None, aeration=2.0, 
                 DO_ID='S_O2', suspended_growth_model=None, 
                 gas_stripping=False, gas_IDs=None, stripping_kLa_min=None, 
                 K_Henry=None, D_gas=None, p_gas_atm=None,
                 isdynamic=True, exogenous_vars=(), **kwargs):
        super().__init__(
            ID=ID, ins=ins, outs=outs, split=split, thermo=thermo,
            init_with=init_with, V_max=V_max, W_tank = W_tank, 
            # D_tank = D_tank,
            # freeboard = freeboard, 
            t_wall = t_wall, t_slab = t_slab, aeration=aeration, 
            DO_ID=DO_ID, suspended_growth_model=suspended_growth_model, 
            gas_stripping=gas_stripping, gas_IDs=gas_IDs, 
            stripping_kLa_min=stripping_kLa_min, 
            K_Henry=K_Henry, D_gas=D_gas, p_gas_atm=p_gas_atm,
            isdynamic=isdynamic, exogenous_vars=exogenous_vars, **kwargs
            )
       

        # # Design parameters 
        self._W_tank = W_tank
        # self._D_tank = D_tank
        # self._freeboard = freeboard
        # self._t_wall = t_wall
        # self._t_slab = t_slab

 
    
        data = load_data(path = Aerobic_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)    
        ############################################### 
        
    def _run(self):
     
        PAC = self.ins[1]
        PAC.imass['X_AlOH'] = self.chemical_PAC_mixing_ratio * self.dosing_flow * self.scale_factor * 0.2886 # kg/L and L/h    to kg/h TODO: add source for 0.2886 kg AlOH/kg PAC
        super()._run()
        
    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),
                             Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),
                             Construction(item='PVC_generic', linked_unit=self, quantity_unit='kg'),]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_weight + self.blower_steel_weight # (self.ppl / self.baseline_ppl)
        design['Cast_iron'] = constr[1].quantity = self.PAC_mixer_cast_iron 
        design['PVC_generic'] = constr[2].quantity = self.dosing_pump_PVC
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        massflow_aerobic = self.ins[0].mass #kg/h
        C['Tank'] = self.aerobic_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
        C['Mixers'] = self.PAC_mixer_cost
        C['Pumps'] = self.dosing_pump_cost
        C['Blower'] = self.blower_cost
        
       # C['Chemical_PAC'] = 
        #PAC cost is already accounted in WasteStream

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
            
        self.add_OPEX = (self._calc_replacement_cost() + 
                         self.chemical_PAC_mixing_ratio * self.dosing_flow * self.scale_factor *
                         self.chemical_PAC_price / 24)
        
        # self.add_OPEX = (self._calc_replacement_cost() + 
        #                  self.chemical_PAC_mixing_ratio * self.PAC_pump_flow * 
        #                  self.chemical_PAC_price) #USD/hr
        
        power_demand = self.power_demand_AerobicTank
        self.power_utility(power_demand) # kWh
    
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) * self.exponent_scale
        Aerobic_tank_replacement_cost = (self.aerobic_tank_cost / self.aerobic_tank_lifetime +
                                        self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
                                        self.pipeline_connectors / self.pipeline_connectors_lifetime + 
                                        self.PAC_mixer_cost / self.PAC_mixer_lifetime + 
                                        self.blower_cost / self.blower_lifetime + 
                                        self.dosing_pump_cost / self.pump_lifetime) * scale
                                        
        Aerobic_tank_replacement_cost = Aerobic_tank_replacement_cost / (365 * 24)  # convert to USD/hr
        return Aerobic_tank_replacement_cost

# %%

MBR_path = os.path.join(EL_su_data_path, '_EL_MBR.tsv')
@price_ratio()

class EL_CMMBR(CompletelyMixedMBR):
    '''
    Completely mixed membrane bioreactor, equivalent to a CSTR with ideal 
    membrane filtration at the outlet.

    See Also
    --------
    :class:`qsdsan.processes.DiffusedAeration`
    :class:`qsdsan.sanunits.CSTR`

    
    '''
    _N_ins = 1
    _N_outs = 2  # [0] filtrate, [1] pumped flow
    _outs_size_is_fixed = True
    ppl = el.ppl
    baseline_ppl = el.baseline_ppl
    exponent_scale=0.1
    
    
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', isdynamic=True, 
                 aeration=2, DO_ID='S_O2', suspended_growth_model=None,
                 pumped_flow=0.103, 
                 solids_capture_rate=0.999, 
                 V_max=3.3, 
                 crossflow_air=None,
                 tank_steel_volume=3.34,
                 steel_density=7850,                
                 **kwargs):
        super().__init__(
            ID=ID, ins=ins, outs=outs, thermo=thermo,
            init_with=init_with, V_max=V_max, aeration=aeration, DO_ID=DO_ID,
            suspended_growth_model=suspended_growth_model, isdynamic=isdynamic, 
            pumped_flow=pumped_flow, 
            solids_capture_rate=solids_capture_rate,  
            crossflow_air=crossflow_air, tank_steel_volume=tank_steel_volume, steel_density=steel_density,
            **kwargs
            )
        
        self.tank_steel_volume=tank_steel_volume
        self.steel_density=steel_density
        
        data = load_data(path = MBR_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)    
        ############################################### 
        
    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),
                             Construction(item ='PVDF_membrane', linked_unit=self, quantity_unit='kg'),
                             Construction(item ='Cast_iron', linked_unit=self, quantity_unit='kg'),]
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_weight + self.nitrate_return_pump_stainless_steel # * (self.ppl / self.baseline_ppl)  # assume linear scaling
        design['PVDF_membrane'] = constr[1].quantity = self.membrane_material_weight
        design['Cast_iron'] = constr[2].quantity = self.nitrate_return_pump_cast_iron + self.self_priming_pump_cast_iron
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['MBR_tank'] = self.MBR_tank_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings
        C['Membrane_material'] = self.membrane_material_price * self.membrane_material_weight
        C['Membrane_cleaning'] = self.membrane_cleaning_fee
        C['Pumps'] = self.nitrate_return_pump_cost + self.self_priming_pump_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio 
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_MBR #TODO: power_demand unit should be kW, but the tsv file is in kWh/day
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        MBR_replacement_cost = (
            self.MBR_tank_cost / self.MBR_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
            self.membrane_material_price * self.membrane_material_weight / self.membrane_material_lifetime +
            self.nitrate_return_pump_cost / self.pump_lifetime + 
            self.self_priming_pump_cost / self.pump_lifetime
            ) * scale
        MBR_replacement_cost = MBR_replacement_cost / (365 * 24) * self.price_ratio # USD/hr
        return MBR_replacement_cost

# %%
ClearWaterTank_path = os.path.join(EL_su_data_path, '_EL_CWT.tsv')

@price_ratio()
class EL_CWT(CSTR):

    '''
    Introduction
    ------------
    To only collect the treated water

    Parameters
    ----------
    Ins:
    (1) influent of treated wastetwater from self-priming pump
    (2) O3 dosing
    (3) air-dissolve influent

    Outs:
    (1) effluent of treated wastewater (clear water)
    (2) spill flow to collection tank


    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.sanunits.storagetank module

    '''
    _N_ins = 1 # treated water, PAC, blower
    _N_outs = 1  # treated water, CH4, N2O
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    exponent_scale = 0.1
    ppl = el.ppl
    baseline_ppl = el.baseline_ppl

    
    _D_O2 = 2.10e-9   # m2/s

    def __init__(self, 
                 ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=12, W_tank = None, 
                 # D_tank = 3.65,
                 # freeboard = 0.61, 
                 t_wall = None, t_slab = None, aeration=None, 
                 DO_ID='S_O2', suspended_growth_model=None, 
                 gas_stripping=False, gas_IDs=None, stripping_kLa_min=None, 
                 K_Henry=None, D_gas=None, p_gas_atm=None,
                 isdynamic=True, exogenous_vars=(), **kwargs):
        super().__init__(
            ID=ID, ins=ins, outs=outs, split=split, thermo=thermo,
            init_with=init_with, V_max=V_max, W_tank = W_tank, 
            # D_tank = D_tank,
            # freeboard = freeboard, 
            t_wall = t_wall, t_slab = t_slab, aeration=aeration, 
            DO_ID=DO_ID, suspended_growth_model=suspended_growth_model, 
            gas_stripping=gas_stripping, gas_IDs=gas_IDs, 
            stripping_kLa_min=stripping_kLa_min, 
            K_Henry=K_Henry, D_gas=D_gas, p_gas_atm=p_gas_atm,
            isdynamic=isdynamic, exogenous_vars=exogenous_vars, **kwargs
            )
       

        # # Design parameters 
        # self._W_tank = W_tank
        # self._D_tank = D_tank
        # self._freeboard = freeboard
        # self._t_wall = t_wall
        # self._t_slab = t_slab

 
    
        data = load_data(path = ClearWaterTank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)    
        ###############################################
        
    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),
                             Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]
       
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_weight # * (self.ppl / self.baseline_ppl)  # to be defined in .tsv file
        design['Cast_iron'] = constr[1].quantity = self.clear_water_air_dissolved_pump_cast_iron + self.clear_water_pump_cast_iron
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['Clear water tank'] = self.clear_water_tank_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings
        C['O3 generator'] = self.O3_generation_machine_cost
        C['Pumps'] = self.clear_water_air_dissolved_pump_cost + self.clear_water_pump_cost
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_CWT
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        CWT_replacement_cost = (
            self.clear_water_tank_cost / self.clear_water_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
            self.O3_generation_machine_cost / self.O3_generation_machine_lifetime + 
            self.clear_water_air_dissolved_pump_cast_iron / self.pump_lifetime + 
            self.clear_water_pump_cast_iron / self.pump_lifetime
            ) * scale
        CWT_replacement_cost = CWT_replacement_cost / (365 * 24) # convert to USD/hr
        return CWT_replacement_cost



# %%
housing_path = os.path.join(EL_su_data_path, '_EL_housing.tsv')

@price_ratio()
class EL_Housing(CSTR):
    '''
     non_reactive unit for the Enviroloo Clear system
    '''
    _N_ins = 1 # treated water, PAC, blower
    _N_outs = 1  # treated water, CH4, N2O
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    exponent_scale = 0.1
    ppl = el.ppl
    baseline_ppl = el.baseline_ppl
 
    
    _D_O2 = 2.10e-9   # m2/s

    def __init__(self, 
                 ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=12, W_tank = None, 
                 # D_tank = 3.65,
                 # freeboard = 0.61, 
                 t_wall = None, t_slab = None, aeration=None, 
                 DO_ID='S_O2', suspended_growth_model=None, 
                 gas_stripping=False, gas_IDs=None, stripping_kLa_min=None, 
                 K_Henry=None, D_gas=None, p_gas_atm=None,
                 isdynamic=True, exogenous_vars=(), **kwargs):
        super().__init__(
            ID=ID, ins=ins, outs=outs, split=split, thermo=thermo,
            init_with=init_with, V_max=V_max, W_tank = W_tank, 
            # D_tank = D_tank,
            # freeboard = freeboard, 
            t_wall = t_wall, t_slab = t_slab, aeration=aeration, 
            DO_ID=DO_ID, suspended_growth_model=suspended_growth_model, 
            gas_stripping=gas_stripping, gas_IDs=gas_IDs, 
            stripping_kLa_min=stripping_kLa_min, 
            K_Henry=K_Henry, D_gas=D_gas, p_gas_atm=p_gas_atm,
            isdynamic=isdynamic, exogenous_vars=exogenous_vars, **kwargs
            )
       
        data = load_data(path = housing_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
    
        for attr, value in kwargs.items():
            setattr(self, attr, value)    
        ############################################### 

    def _init_lca(self): # replace the actual materials used in the EL
        self.construction = [
            Construction(item = 'StainlessSteel', linked_unit= self, quantity_unit= 'kg'),
            Construction(item = 'Plastic', linked_unit= self, quantity_unit= 'kg'),]

    def _design(self): # replace the actual materials used in the EL
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = (self.steel_weight + self.steel_framework_weight + self.steel_fittings_weight) #* (self.ppl / self.baseline_ppl)  # assume linear scaling
        design['Plastic'] = constr[1].quantity = (self.LLDPE_weight) #* (self.ppl / self.baseline_ppl)   # assume linear scaling
        self.add_construction(add_cost= False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Housing'] = (self.frame + self.extrusion + 
                        self.angle_frame + self.angle +
                        self.door_sheet + self.plate +
                        self.powder_coating) * (1 + 0.1 * (self.N_EL -1))
        C['Frontend'] = self.frontend_1000_ppl_cost
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        self.add_OPEX = self._calc_replacement_cost()

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        housing_replacement_cost = (
            self.frame / self.frame_lifetime +
            self.angle_frame / self.angle_frame_lifetime +
            self.door_sheet / self.door_sheet_lifetime +
            self.powder_coating / self.powder_coating_lifetime + 
            self.frontend_1000_ppl_cost / self.frontend_lifetime
            ) * scale
        housing_replacement_cost = housing_replacement_cost / (365 * 24) # convert to USD/hr
        return housing_replacement_cost
    
    @property
    def N_EL(self): # determine the number of EL system needed
        return ceil(self.ppl / self.baseline_ppl)
    
    @property
    def N_toilets(self): # determine the number of toilets needed
        return ceil(self.ppl / self.ppl_per_MURT)

# %%
system_path = os.path.join(EL_su_data_path, '_EL_system.tsv')

@price_ratio()
class EL_System(CSTR):
    '''
    Relate to connection components in the EL system
    '''
    _N_ins = 1 # treated water, PAC, blower
    _N_outs = 1  # treated water, CH4, N2O
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    exponent_scale = 0.1
    ppl = el.ppl
    baseline_ppl = el.baseline_ppl
 
    _D_O2 = 2.10e-9   # m2/s

    def __init__(self, 
                     ID='', ins=None, outs=(), split=None, thermo=None,
                     init_with='WasteStream', V_max=12, W_tank = None, 
                     # D_tank = 3.65,
                     # freeboard = 0.61, 
                     t_wall = None, t_slab = None, aeration=None, 
                     DO_ID='S_O2', suspended_growth_model=None, 
                     gas_stripping=False, gas_IDs=None, stripping_kLa_min=None, 
                     K_Henry=None, D_gas=None, p_gas_atm=None,
                     isdynamic=True, exogenous_vars=(), **kwargs):
        super().__init__(
                ID=ID, ins=ins, outs=outs, split=split, thermo=thermo,
                init_with=init_with, V_max=V_max, W_tank = W_tank, 
                # D_tank = D_tank,
                # freeboard = freeboard, 
                t_wall = t_wall, t_slab = t_slab, aeration=aeration, 
                DO_ID=DO_ID, suspended_growth_model=suspended_growth_model, 
                gas_stripping=gas_stripping, gas_IDs=gas_IDs, 
                stripping_kLa_min=stripping_kLa_min, 
                K_Henry=K_Henry, D_gas=D_gas, p_gas_atm=p_gas_atm,
                isdynamic=isdynamic, exogenous_vars=exogenous_vars, **kwargs
                )
           
        data = load_data(path = windSolar_path)
        for para in data.index:
                value = float(data.loc[para]['expected'])
                setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
                setattr(self, attr, value)    
            ###############################################

    def _init_lca(self):
        self.construction = [
            Construction(item = 'PVC_generic', linked_unit= self, quantity_unit= 'kg'),
            Construction(item = 'HDPE', linked_unit= self, quantity_unit= 'kg'),
            Construction(item='Cast_iron', linked_unit=self, quantity_unit='kg'),]

    def _design(self):
        design = self.design_results
        design['PVC_generic'] = self.construction[0].quantity = self.PVC_weight * (self.ppl / self.baseline_ppl)
        design['HDPE'] = self.construction[1].quantity = self.HDPE_weight * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost= False)
    
    def _cost(self): # replace these items below by that listed in the _EL_system.tsv
        C = self.baseline_purchase_costs
        C['System'] = (
            self.membrane_filters_M +
            self.membrane_filters_size +
            self.membrane_filters_pause_size +
            self.membrane_filters_chassis_M +
            self.membrane_filters_air_diffuser +
            self.membrane_filters_air_diffuser_chassis +
            self.overflow_membrane2collection +
            self.overflow_clear_water2collection +
            self.overflow_primary_clarifier2anoixc +
            self.overflow_anoxic2aerobic +
            self.overflow_aerobic2membrane +
            self.ozone_pipeline_200m +
            self.pipeline_fittings_32mm +
            self.pipeline_fittings_40mm +
            self.pipeline_fittings_60mm +
            self.ball_valves_50mm +
            self.aerobic_air_diffuser)

        ratio = self.price_ratio # ratio of the price of the new system to the baseline system
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()  # add the cost of replacement

        if self.if_gridtied:
            power_demand = (self.power_demand_system / 1000) * self.N_EL  # in W/d
        else:
            power_demand = 0
        
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        system_replacement_cost = (
            self.membrane_filters_M / self.lifetime_membrane_filters_M +
            self.membrane_filters_size / self.lifetime_membrane_filters_size +
            # self.aerobic_basin / self.lifetime_aerobic_basin +
            self.membrane_filters_air_diffuser / self.lifetime_membrane_filters_air_diffuser +
            self.membrane_filters_air_diffuser_chassis / self.lifetime_membrane_filters_air_diffuser_chassis +
            self.overflow_membrane2collection / self.lifetime_overflow_membrane2collection +
            self.overflow_clear_water2collection / self.lifetime_overflow_clear_water2collection +
            self.overflow_primary_clarifier2anoixc / self.lifetime_overflow_primary_clarifier2anoixc +
            self.overflow_anoxic2aerobic / self.lifetime_overflow_anoxic2aerobic +
            self.overflow_aerobic2membrane / self.lifetime_overflow_aerobic2membrane +
            self.ozone_pipeline_200m / self.lifetime_200m_ozone_pipeline +
            self.pipeline_fittings_32mm / self.lifetime_32mm_pipeline_fittings +
            self.pipeline_fittings_40mm / self.lifetime_40mm_pipeline_fittings +
            self.pipeline_fittings_60mm / self.lifetime_60mm_pipeline_fittings +
            self.ball_valves_50mm / self.lifetime_50mm_ball_valves) * scale
        system_replacement_cost = system_replacement_cost / (365 * 24)   # convert from USD/year to USD/hour
        return system_replacement_cost
    @property
    def N_EL(self): # determine the number of EL system needed
       return ceil(self.ppl / self.baseline_ppl)
# %%

windSolar_path = os.path.join(EL_su_data_path, '_EL_photovoltaic_wind.tsv')

@price_ratio()
class EL_WindSolar(CSTR):
    
    '''
    
    The photovoltaic and wind power configuration of the EnviroLoo System,
    which is composed of solar panels, batteries, inverters, cables, carport,
    charger, turbine, backplates, and other parts. Costs include sales tax, licensing, 
    consulting, and installation fees.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    Parameters
    ----------
    user_scale_up : float
        Scaling up of consumables, electricity demand, capital parts,
        and replacement parts based on number of input users compared to the
        baseline design number of users

    References
    ----------
   '''

    _N_ins = 1 # treated water, PAC, blower
    _N_outs = 1  # treated water, CH4, N2O
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    exponent_scale = 0.1
    pl = el.ppl
    baseline_ppl = el.baseline_ppl
    
    _D_O2 = 2.10e-9   # m2/s

    def __init__(self, 
                 ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=12, W_tank = None, 
                 # D_tank = 3.65,
                 # freeboard = 0.61, 
                 t_wall = None, t_slab = None, aeration=None, 
                 DO_ID='S_O2', suspended_growth_model=None, 
                 gas_stripping=False, gas_IDs=None, stripping_kLa_min=None, 
                 K_Henry=None, D_gas=None, p_gas_atm=None,
                 isdynamic=True, exogenous_vars=(), **kwargs):
        super().__init__(
            ID=ID, ins=ins, outs=outs, split=split, thermo=thermo,
            init_with=init_with, V_max=V_max, W_tank = W_tank, 
            # D_tank = D_tank,
            # freeboard = freeboard, 
            t_wall = t_wall, t_slab = t_slab, aeration=aeration, 
            DO_ID=DO_ID, suspended_growth_model=suspended_growth_model, 
            gas_stripping=gas_stripping, gas_IDs=gas_IDs, 
            stripping_kLa_min=stripping_kLa_min, 
            K_Henry=K_Henry, D_gas=D_gas, p_gas_atm=p_gas_atm,
            isdynamic=isdynamic, exogenous_vars=exogenous_vars, **kwargs
            )
       
        data = load_data(path = windSolar_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
    
        for attr, value in kwargs.items():
            setattr(self, attr, value)    
        ############################################### 
        
    def _init_lca(self):
        self.construction = [
            Construction(item='PhotovoltaicPanel', linked_unit=self, quantity_unit='m2'),
            Construction(item='Battery', linked_unit=self, quantity_unit='kg',),
            Construction(item='ElectricCables', linked_unit=self, quantity_unit='m'),
            Construction(item='Aluminum', linked_unit=self, quantity_unit='kg'),
            Construction(item='Steel', linked_unit=self, quantity_unit='kg')            
            ]

    def _design(self):
        design = self.design_results
        constr = self.construction
        # user_scale_up = self.user_scale_up
        design['PhotovoltaicPanel'] = self.construction[0].quantity = self.pv_photovoltaic_panel_area # * (user_scale_up**0.6)
        design['Battery'] = constr[1].quantity = self.el_battery_weight        #* (user_scale_up**0.6)
        design['ElectricCables'] = constr[2].quantity = self.pv_cable_length
        # design['Aluminum'] = constr[3].quantity = self.pv_aluminum_weight
        # design['Aluminum'] = constr[4].quantity = self.pv_aluminum_weight
        design['Steel'] = constr[4].quantity = self.el_wind_galvanized_steel_weight 
        + self.pv_carport_weight_galvanized_metal + self.pv_inverter_weight_steel 
        + self.pv_charger_weight_steel

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Battery'] = (
            self.el_battery
            )
        #C['Battery'] += self.pv_battery * (self.user_scale_up**0.6) optional scale-up ratio

        C['PhotovoltaicPanel'] = (
            #self.el_panels * (self.user_scale_up**0.6) +
            self.el_panels + 
            self.el_carport +
            self.el_backplate
            )
        C['Misc'] = (
            self.el_template +
            self.el_controller +
            self.el_certificate +
            self.el_consulting +
            self.el_installation +
            self.el_inverter + 
            self.el_charger +
            self.el_meter
            )
        C['Wind'] = (
            self.el_wind_turbine
            )
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

    def _calc_replacement_cost(self):
        battery_replacememt_parts_annual_cost = self.el_battery / self.el_battery_lifetime
        wind_replacement_cost = self.el_wind_turbine / self.wind_turbine_lifetime
        photovoltaic_replacement_cost = self.el_panels / self.el_panel_lifetime
        misc_replacement_cost = self.el_inverter / self.el_inverter_lifetime
        
        total_replacement_cost = (
            battery_replacememt_parts_annual_cost +
            wind_replacement_cost +
            photovoltaic_replacement_cost +
            misc_replacement_cost
            )
        return total_replacement_cost/ (365 * 24) * self.price_ratio # USD/hr

    def _calc_maintenance_labor_cost(self):
        photovoltaic_maintenance_labor = (
            self.pv_labor_replacement_misc_repairs
            ) * self.wages
        return photovoltaic_maintenance_labor/ (365 * 24) # USD/hr
        