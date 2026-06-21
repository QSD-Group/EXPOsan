#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''

This module is developed by:
    Aaron Marszewski <aaronpm3@illinois.edu>
    Rishabh Puri <rp34@illinois.edu>
    Yuyao Huang <yuyaoh20@gmail.com>
    Siqi Tang <siqit@outlook.com>
for Enviroloo (EL) Clear Toilet system

'''
# %%
import numpy as np
import qsdsan
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
    'StruviteReactor', # Struvite precipitation reactor (upstream P removal)
    'StruviteRedissolution', # Struvite redissolution (unsettled fines)
    )
el_path = os.path.dirname(__file__)
module = os.path.split(el_path)[-1]
data_path, results_path = _init_modules(module, include_data_path = True)
EL_su_data_path = os.path.join(data_path, 'units_data')
SR_data_path = os.path.join(EL_su_data_path, '_EL_SR.tsv')


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
        C['Membrane_material'] = self.membrane_material_price
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
            self.membrane_material_price / self.membrane_material_lifetime +
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
# =============================================================================
# StruviteReactor
#
# Continuous flow reactor for upstream phosphorus removal from
# source-separated urine in the EnviroLoo Clear MURT system.
#
# Design basis:
#     - Continuous flow reactor (not batch); HRT = 1 hr (Etter et al. 2011)
#     - precip_yield = 0.70 (harvest efficiency, Triangular(0.55, 0.70, 0.85)).
#       Expected value confirmed by EnviroLoo field reactor trial harvest
#       efficiency of 68.7% (IC_P_Recovery_Overview_Illinois_University_
#       Updated.pdf). Low bound (0.55) consistent with pilot-scale struvite
#       harvest losses from crystal washout (Antonini et al. 2011 pilot:
#       94% precipitation but only 55% recovery). High bound (0.85)
#       reflects potential for improved harvest under controlled conditions.
#     - eff_PO4_mgL = 6.55 mg/L (fixed-target effluent PO4, distinct from
#       precip_yield). Triangular(1.31, 6.55, 27.50). Calculated from the
#       corrected influent PO4 to SR (~65.5 mg/L, see system file for the
#       urine diversion mass-balance correction) at 90% P removal (Etter
#       et al. 2011; Maurer et al. 2006). Low bound (1.31 mg/L) uses 98%
#       removal (Antonini et al. 2011, validated on diluted urine). High
#       bound (27.50 mg/L) uses 58% removal (EnviroLoo field reactor data).
#     - Mg:P molar ratio = 1.1:1 (Etter et al. 2011); Mg source: MgCl2
#     - Mixing achieved via air blower (aeration), consistent with field
#       reactor design (IC_P_Recovery_Overview: air blower component 3) —
#       no mechanical stirrer required
#
# CAPEX items (from _EL_SR.tsv):
#     - Conical plastic reactor tank ($2,600/m3, Lohman et al. 2020 SI
#       Table S17)
#     - MgCl2 dosing pump ($150, AliExpress BOM)
#     - Air blower ($40, BOM — same source/value as _EL_blower.tsv)
#     - Pipeline connectors ($13.16, eco-filtration.co.uk BOM)
#     - Weld female adapter fittings ($6.87, thedrainagesource.com BOM)
#
# OPEX items (from _EL_SR.tsv):
#     - MgCl2 chemical cost ($0.50/kg expected, South Africa-specific
#       range $0.19-$1.50/kg) — dominant cost driver
#     - Air blower electricity (2 kWh/day — consistent with _EL_blower.tsv
#       and _EL_Aerobic.tsv power_demand_blower)
#     - Component replacement costs (tank, pump, blower, piping)
#
# PARAMETER PRECEDENCE:
#     All cost/design parameters are loaded from _EL_SR.tsv FIRST as the
#     baseline/default. Any value explicitly passed as a constructor
#     keyword argument (e.g., eff_PO4_mgL=6.55, precip_yield=0.70) is
#     applied AFTER the TSV load and therefore correctly takes precedence.
#     Monte Carlo uncertainty sampling (via batch_setting_unit_params in
#     Enviroloo_model.py) sets attributes directly on the already-
#     constructed instance per trial, unaffected by this ordering.
#
# References:
#     Etter et al. 2011. Water Res. 45, 852-862.
#     Maurer, M., Pronk, W., Larsen, T.A. 2006. Treatment processes for
#         source-separated urine. Water Res. 40, 3151-3166.
#     Antonini, S., et al. 2011. Nitrogen and Phosphorus Recovery from
#         Human Urine by Struvite Precipitation and Air Stripping in
#         Vietnam. CLEAN-Soil Air Water 39(12), 1099-1104.
#     Lohman et al. 2020. Environ. Sci. Technol. 54, 9217-9227. SI Table S17.
#     IC_P_Recovery_Overview_Illinois_University_Updated.pdf (EnviroLoo
#         field reactor report — confirms air blower for mixing, harvest
#         efficiency 68.7%, P removal from 32.8 to 13.8 mg/L TP = 58%).
# =============================================================================

@price_ratio()
class StruviteReactor(SanUnit):
    """
    Struvite precipitation reactor — dynamic-capable (AE unit).

    Continuous flow reactor for upstream phosphorus removal from
    source-separated urine in the EnviroLoo Clear MURT system.

    Mixing is achieved via an air blower (aeration), consistent with the
    EnviroLoo field reactor design (IC_P_Recovery_Overview). No mechanical
    stirrer is used.

    Two distinct, independently uncertain performance parameters:
        eff_PO4_mgL  : fixed-target effluent PO4 concentration (mg/L),
                       representing P removal from solution.
        precip_yield : fraction of precipitated struvite captured as solid
                       (harvest efficiency), representing how much of the
                       precipitated mass is actually recovered vs. carried
                       forward as unsettled fines (handled downstream by
                       StruviteRedissolution).

    Parameters loaded from _EL_SR.tsv (same approach as other EnviroLoo units):
        tank_cost_per_m3, tank_lifetime
        dosing_pump_cost, pump_lifetime
        blower_cost, blower_lifetime, power_demand_blower
        pipeline_connectors, pipeline_connectors_lifetime
        weld_female_adapter_fittings, weld_female_adapter_fittings_lifetime
        price_MgCl2_per_kg
        eff_PO4_mgL, precip_yield

    Inputs
    ------
    ins[0] : liquid influent (urine + MgCl2 dose, pre-mixed via Mixer)

    Outputs
    -------
    outs[0] : struvite_recovered  — solid struvite product  (phase 's')
    outs[1] : loss                — reserved, always zero   (phase 'l')
    outs[2] : effluent            — liquid effluent          (phase 'l')
    """

    _N_ins  = 1
    _N_outs = 3
    _ins_size_is_fixed  = True
    _outs_size_is_fixed = True
    _neg_tol = -1e-9

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 component_ID_NH3='S_NH4',
                 component_ID_P='S_PO4',
                 component_ID_Mg='S_Mg',
                 component_ID_struvite='X_struv',
                 precip_yield=None,
                 eff_PO4_mgL=None,
                 HRT_hr=1.0,
                 init_with='WasteStream',
                 F_BM_default=None,
                 isdynamic=True,
                 dose_MgCl2_kg_d=None,
                 pH_ctrl=None,
                 **kwargs):

        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo,
                         init_with=init_with, isdynamic=isdynamic,
                         F_BM_default=F_BM_default, **kwargs)

        self._component_ID_NH3      = component_ID_NH3
        self._component_ID_P        = component_ID_P
        self._component_ID_Mg       = component_ID_Mg
        self._component_ID_struvite = component_ID_struvite

        self.HRT_hr          = float(HRT_hr)
        self.pH_ctrl         = pH_ctrl
        self.dose_MgCl2_kg_d = dose_MgCl2_kg_d

        # -------------------------------------------------------------------
        # STEP 1: Load ALL parameters from _EL_SR.tsv FIRST as defaults.
        # This includes eff_PO4_mgL and precip_yield among all other rows.
        # Same pattern as EL_CT, EL_PC, EL_Anoxic etc.
        # -------------------------------------------------------------------
        data = load_data(path=SR_data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        # -------------------------------------------------------------------
        # STEP 2: Apply explicit constructor arguments AFTER the TSV load,
        # so they correctly take precedence over the TSV default whenever
        # they are explicitly provided (not None).
        # -------------------------------------------------------------------
        if precip_yield is not None:
            self.precip_yield = float(precip_yield)
        if eff_PO4_mgL is not None:
            self.eff_PO4_mgL = float(eff_PO4_mgL)

        # -------------------------------------------------------------------
        # STEP 3: Apply any remaining **kwargs LAST, so they take final
        # precedence over both the TSV and the explicit named arguments
        # above (matches the override behavior of every other EL unit).
        # -------------------------------------------------------------------
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    # --- Component ID properties ---
    @property
    def component_ID_NH3(self):      return self._component_ID_NH3
    @property
    def component_ID_P(self):        return self._component_ID_P
    @property
    def component_ID_Mg(self):       return self._component_ID_Mg
    @property
    def component_ID_struvite(self): return self._component_ID_struvite

    # -------------------------------------------------------------------------
    # Static mass balance
    # -------------------------------------------------------------------------
    def _run(self):
        influent, = self.ins
        recovered, loss, effluent = self.outs

        recovered.phase = 's'
        loss.phase      = 'l'
        effluent.phase  = 'l'

        cmps = influent.components
        IDs  = cmps.IDs

        P_id     = self.component_ID_P
        NH4_id   = self.component_ID_NH3
        Mg_id    = self.component_ID_Mg
        struv_id = self.component_ID_struvite

        MW_P   = cmps[P_id].MW
        MW_NH4 = cmps[NH4_id].MW
        MW_Mg  = cmps[Mg_id].MW
        MW_STR = cmps[struv_id].MW

        effluent.mass[:]  = influent.mass
        recovered.mass[:] = 0.0
        loss.mass[:]      = 0.0

        if influent.state is not None:
            Q = max(float(influent.state[-1]), 1e-12)
        else:
            Q = max(float(influent.F_vol) * 24.0, 1e-12)

        C_PO4_in   = float(influent.iconc[P_id])   if P_id   in IDs else 0.0
        C_NH4_in   = float(influent.iconc[NH4_id]) if NH4_id in IDs else 0.0
        C_Mg_in    = float(influent.iconc[Mg_id])  if Mg_id  in IDs else 0.0
        C_struv_in = float(influent.iconc[struv_id]) if struv_id in IDs else 0.0

        if (C_PO4_in < self._neg_tol or
            C_NH4_in < self._neg_tol or
            C_Mg_in  < self._neg_tol):
            raise ValueError(
                f'{self.ID}: negative inlet concentration '
                f'(PO4={C_PO4_in:g}, NH4={C_NH4_in:g}, Mg={C_Mg_in:g} mg/L).')

        M_PO4_in   = C_PO4_in   * Q
        M_NH4_in   = C_NH4_in   * Q
        M_Mg_in    = C_Mg_in    * Q
        M_struv_in = C_struv_in * Q

        M_PO4_target_eff = self.eff_PO4_mgL * Q
        M_PO4_req        = max(M_PO4_in - M_PO4_target_eff, 0.0)

        mol_form = max(min(
            M_PO4_req / MW_P,
            M_NH4_in  / MW_NH4,
            M_Mg_in   / MW_Mg,
        ), 0.0)

        M_PO4_eff = max(M_PO4_in - mol_form * MW_P,   0.0)
        M_NH4_eff = max(M_NH4_in - mol_form * MW_NH4, 0.0)
        M_Mg_eff  = max(M_Mg_in  - mol_form * MW_Mg,  0.0)

        M_struv_total = M_struv_in + mol_form * MW_STR
        frac_rec      = float(np.clip(self.precip_yield, 0.0, 1.0))
        M_struv_rec   = M_struv_total * frac_rec
        M_struv_unrec = M_struv_total * (1.0 - frac_rec)

        effluent.copy_like(influent)
        recovered.mass[:] = 0.0
        loss.mass[:]      = 0.0

        effluent.imass[P_id]     = M_PO4_eff    / 24.0 / 1000.0
        effluent.imass[NH4_id]   = M_NH4_eff    / 24.0 / 1000.0
        effluent.imass[Mg_id]    = M_Mg_eff     / 24.0 / 1000.0
        effluent.imass[struv_id] = M_struv_unrec / 24.0 / 1000.0

        recovered.imass[struv_id] = M_struv_rec / 24.0 / 1000.0
        for cmp in IDs:
            if cmp != struv_id:
                recovered.imass[cmp] = 0.0

        if 'H2O' in IDs:
            w = recovered.imass['H2O']
            if w > 0:
                recovered.imass['H2O'] = 0.0
                effluent.imass['H2O'] += w

        if self.pH_ctrl is not None:
            for ws in (recovered, loss, effluent):
                try:
                    ws.pH  = self.pH_ctrl
                    ws._pH = self.pH_ctrl
                except AttributeError:
                    pass

    # -------------------------------------------------------------------------
    # Dynamic methods — AE
    # -------------------------------------------------------------------------
    def _init_state(self):
        eff = self.outs[2]
        n   = len(self.components)
        self._state      = np.empty(n + 1)
        self._state[:-1] = eff.conc
        self._state[-1]  = max(eff.F_vol * 24, 1e-12)
        self._dstate     = np.zeros(n + 1)

        dummy = np.zeros(n + 1)
        dummy[-1] = 1e-12
        for ws in (self.outs[0], self.outs[1]):
            if ws.state  is None: ws.state  = dummy.copy()
            else:                 ws.state[:]  = dummy
            if ws.dstate is None: ws.dstate = np.zeros(n + 1)
            else:                 ws.dstate[:] = 0.0

        if eff.state  is None: eff.state  = self._state.copy()
        else:                  eff.state[:]  = self._state
        if eff.dstate is None: eff.dstate = self._dstate.copy()
        else:                  eff.dstate[:] = self._dstate

    def _update_state(self):
        dummy = np.zeros_like(self._state)
        dummy[-1] = 1e-12
        for ws in (self.outs[0], self.outs[1]):
            if ws.state is None: ws.state = dummy.copy()
            else:                ws.state[:] = dummy
        eff = self.outs[2]
        if eff.state is None: eff.state = self._state.copy()
        else:                 eff.state[:] = self._state

    def _update_dstate(self):
        for ws in (self.outs[0], self.outs[1]):
            if ws.dstate is None: ws.dstate = np.zeros_like(self._dstate)
            else:                 ws.dstate[:] = 0.0
        eff = self.outs[2]
        if eff.dstate is None: eff.dstate = self._dstate.copy()
        else:                  eff.dstate[:] = self._dstate

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _state         = self._state
        _dstate        = self._dstate
        _update_state  = self._update_state
        _update_dstate = self._update_dstate

        cmps = self.components
        IDs  = list(cmps.IDs)

        idx_P     = IDs.index(self.component_ID_P)
        idx_NH4   = IDs.index(self.component_ID_NH3)
        idx_Mg    = IDs.index(self.component_ID_Mg)
        idx_struv = IDs.index(self.component_ID_struvite)

        MW_P    = cmps[self.component_ID_P].MW
        MW_NH4  = cmps[self.component_ID_NH3].MW
        MW_Mg   = cmps[self.component_ID_Mg].MW
        MW_STR  = cmps[self.component_ID_struvite].MW

        eff_PO4_gm3 = self.eff_PO4_mgL
        frac_rec    = float(np.clip(self.precip_yield, 0.0, 1.0))
        frac_unrec  = 1.0 - frac_rec

        def y_t(t, y_ins, dy_ins):
            Q_in  = max(y_ins[0, -1], 1e-12)
            C_in  = y_ins[0, :-1]
            dQ_in = dy_ins[0, -1]

            M_in = C_in * Q_in

            m_PO4_target_eff = eff_PO4_gm3 * Q_in
            m_PO4_req = max(M_in[idx_P] - m_PO4_target_eff, 0.0)

            mol_form = max(min(
                m_PO4_req        / MW_P,
                M_in[idx_NH4]    / MW_NH4,
                M_in[idx_Mg]     / MW_Mg,
            ), 0.0)

            m_PO4_eff  = max(M_in[idx_P]   - mol_form * MW_P,   0.0)
            m_NH4_eff  = max(M_in[idx_NH4] - mol_form * MW_NH4, 0.0)
            m_Mg_eff   = max(M_in[idx_Mg]  - mol_form * MW_Mg,  0.0)

            m_struv_total = M_in[idx_struv] + mol_form * MW_STR
            m_struv_un    = m_struv_total * frac_unrec

            C_eff            = C_in.copy()
            C_eff[C_eff < 1e-12] = 1e-12
            C_eff[idx_P]     = m_PO4_eff  / Q_in
            C_eff[idx_NH4]   = m_NH4_eff  / Q_in
            C_eff[idx_Mg]    = m_Mg_eff   / Q_in
            C_eff[idx_struv] = m_struv_un / Q_in

            _state[:-1] = C_eff
            _state[-1]  = max(Q_in, 1e-9)
            _dstate[:-1] = 0.0
            _dstate[-1]  = dQ_in

            _update_state()
            _update_dstate()

        self._AE = y_t

    # -------------------------------------------------------------------------
    # Design — loads values from _EL_SR.tsv
    # -------------------------------------------------------------------------
    def _design(self):
        d        = self.design_results
        influent = self.ins[0]

        # Tank volume based on HRT = 1 hr (continuous flow)
        # F_vol in m3/hr × HRT_hr = tank volume in m3
        d['Tank_volume_m3'] = float(influent.F_vol) * self.HRT_hr
        d['MgCl2_kg_d']     = float(self.dose_MgCl2_kg_d or 0.0)

    # -------------------------------------------------------------------------
    # Cost — all values loaded from _EL_SR.tsv via __init__
    # Blower included for aeration/mixing, consistent with field reactor design
    # and _EL_blower.tsv / _EL_Aerobic.tsv precedent in other EL units
    # -------------------------------------------------------------------------
    def _cost(self):
        C = self.baseline_purchase_costs
        d = self.design_results

        # --- CAPEX ---
        # Tank cost: tank_cost_per_m3 × volume (Lohman 2020 SI Table S17)
        C['Tank']     = self.tank_cost_per_m3 * d['Tank_volume_m3']

        # Dosing pump (AliExpress BOM, same as EL_Anoxic dosing_pump_cost)
        C['Dosing_pump'] = self.dosing_pump_cost

        # Air blower for mixing/CO2 stripping (BOM, same source/value as
        # _EL_blower.tsv blower_cost; confirmed by field reactor report)
        C['Blower'] = self.blower_cost

        # Piping and fittings (same references and values as all other EL units)
        C['Pipeline_connectors']          = self.pipeline_connectors
        C['Weld_female_adapter_fittings'] = self.weld_female_adapter_fittings

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        # --- OPEX ---
        # MgCl2 chemical cost — dominant cost driver
        mgcl2_kg_d = float(d.get('MgCl2_kg_d', 0) or 0)
        mgcl2_cost_per_hr = (mgcl2_kg_d * self.price_MgCl2_per_kg) / 24.0

        # Blower electricity cost
        # power_demand_blower in kWh/day → converted to kW for power_utility
        blower_power_kW = self.power_demand_blower / 24.0   # kWh/day ÷ 24 = kW
        self.power_utility(blower_power_kW)

        # Total OPEX = MgCl2 + component replacement
        self.add_OPEX = mgcl2_cost_per_hr + self._calc_replacement_cost()

    # -------------------------------------------------------------------------
    # Replacement cost — includes blower consistent with _EL_Aerobic pattern
    # -------------------------------------------------------------------------
    def _calc_replacement_cost(self):
        d = self.design_results
        replacement = (
            (self.tank_cost_per_m3 * d['Tank_volume_m3']) / (self.tank_lifetime     * 365 * 24) +
            self.dosing_pump_cost                          / (self.pump_lifetime      * 365 * 24) +
            self.blower_cost                               / (self.blower_lifetime    * 365 * 24) +
            self.pipeline_connectors                       / (self.pipeline_connectors_lifetime * 365 * 24) +
            self.weld_female_adapter_fittings              / (self.weld_female_adapter_fittings_lifetime * 365 * 24)
        )
        return replacement

    # -------------------------------------------------------------------------
    # LCA — plastic tank + blower steel material
    # Blower steel weight from _EL_blower.tsv: 24 kg (same blower)
    # -------------------------------------------------------------------------
    def _init_lca(self):
        vol = self.design_results.get('Tank_volume_m3', 0.175)
        self.construction = [
            Construction(item='Plastic',      linked_unit=self,
                         quantity_unit='kg',
                         quantity=vol * 950.0 * 0.01),
            Construction(item='Steel',        linked_unit=self,
                         quantity_unit='kg',
                         quantity=24.0),   # blower steel weight from _EL_blower.tsv
        ]

    # -------------------------------------------------------------------------
    # Utility methods
    # -------------------------------------------------------------------------
    def report_recovery(unit):
        ws_in  = unit.ins[0]
        ws_eff = unit.outs[2]

        cmps     = unit.components
        P_id     = unit.component_ID_P
        NH4_id   = unit.component_ID_NH3
        Mg_id    = unit.component_ID_Mg
        struv_id = unit.component_ID_struvite

        MW_P   = cmps[P_id].MW
        MW_NH4 = cmps[NH4_id].MW
        MW_Mg  = cmps[Mg_id].MW
        MW_STR = cmps[struv_id].MW

        Q_m3_d = (max(float(ws_in.state[-1]), 1e-12)
                  if ws_in.state is not None
                  else max(float(ws_in.F_vol) * 24.0, 1e-12))

        C_PO4_in  = float(ws_in.iconc[P_id])
        C_NH4_in  = float(ws_in.iconc[NH4_id])
        C_Mg_in   = float(ws_in.iconc[Mg_id])
        C_PO4_eff = float(ws_eff.iconc[P_id])
        C_NH4_eff = float(ws_eff.iconc[NH4_id])
        C_Mg_eff  = float(ws_eff.iconc[Mg_id])

        removal_P = (max((C_PO4_in - C_PO4_eff) / C_PO4_in, 0.0)
                     if C_PO4_in > 1e-12 else 0.0)

        M_PO4_req   = max(C_PO4_in*Q_m3_d - unit.eff_PO4_mgL*Q_m3_d, 0.0)
        mol_form    = max(min(M_PO4_req/MW_P,
                              C_NH4_in*Q_m3_d/MW_NH4,
                              C_Mg_in*Q_m3_d/MW_Mg), 0.0)
        C_STR_in    = float(ws_in.iconc[struv_id])
        M_struv_rec = (C_STR_in*Q_m3_d + mol_form*MW_STR) * float(
            np.clip(unit.precip_yield, 0, 1))

        tank_vol = unit.design_results.get(
            'Tank_volume_m3', unit.HRT_hr * float(unit.ins[0].F_vol))

        print("\n===== STRUVITE REACTOR SUMMARY =====")
        print(f"Q urine (m3/d)          : {Q_m3_d:.3f}")
        print(f"Tank volume (m3)        : {tank_vol:.4f}  [HRT={unit.HRT_hr} hr]")
        print(f"Inlet  PO4  (mg/L)      : {C_PO4_in:.3f}")
        print(f"Inlet  NH4  (mg/L)      : {C_NH4_in:.3f}")
        print(f"Inlet  Mg   (mg/L)      : {C_Mg_in:.3f}")
        print(f"Effluent PO4 (mg/L)     : {C_PO4_eff:.3f}")
        print(f"Effluent NH4 (mg/L)     : {C_NH4_eff:.3f}")
        print(f"Effluent Mg  (mg/L)     : {C_Mg_eff:.3f}")
        print(f"P removal               : {removal_P*100:.1f}%")
        print(f"Struvite recovered(g/hr): {M_struv_rec/24:.3f}")
        print(f"MgCl2 dose (kg/d)       : {unit.dose_MgCl2_kg_d:.4f}")
        print(f"eff_PO4_mgL (target)    : {unit.eff_PO4_mgL}")
        print(f"precip_yield            : {unit.precip_yield}")
        print(f"tank_cost_per_m3 (USD)  : {unit.tank_cost_per_m3}")
        print(f"dosing_pump_cost (USD)  : {unit.dosing_pump_cost}")
        print(f"blower_cost (USD)       : {unit.blower_cost}")
        print(f"power_demand_blower     : {unit.power_demand_blower} kWh/day")
        print(f"price_MgCl2 (USD/kg)    : {unit.price_MgCl2_per_kg}")
        print("=====================================\n")

    def update_product_stream(unit):
        ws_in = unit.ins[0]
        recovered, loss, effluent = unit.outs

        cmps     = unit.components
        P_id     = unit.component_ID_P
        NH4_id   = unit.component_ID_NH3
        Mg_id    = unit.component_ID_Mg
        struv_id = unit.component_ID_struvite

        MW_P   = cmps[P_id].MW
        MW_NH4 = cmps[NH4_id].MW
        MW_Mg  = cmps[Mg_id].MW
        MW_STR = cmps[struv_id].MW

        Q_m3_d = (max(float(ws_in.state[-1]), 1e-12)
                  if ws_in.state is not None
                  else max(float(ws_in.F_vol) * 24.0, 1e-12))

        C_PO4_in  = float(ws_in.iconc[P_id])
        C_NH4_in  = float(ws_in.iconc[NH4_id])
        C_Mg_in   = float(ws_in.iconc[Mg_id])
        C_STR_in  = float(ws_in.iconc[struv_id])

        M_PO4_req   = max(C_PO4_in*Q_m3_d - unit.eff_PO4_mgL*Q_m3_d, 0.0)
        mol_form    = max(min(M_PO4_req/MW_P,
                              C_NH4_in*Q_m3_d/MW_NH4,
                              C_Mg_in*Q_m3_d/MW_Mg), 0.0)
        M_struv_rec = (C_STR_in*Q_m3_d + mol_form*MW_STR) * float(
            np.clip(unit.precip_yield, 0, 1))

        recovered.mass[:] = 0.0
        if struv_id in recovered.components.IDs:
            recovered.imass[struv_id] = M_struv_rec / 24.0 / 1000.0


# %%
# =============================================================================
# StruviteRedissolution
#
# Dissolves residual unsettled X_struv particles (from SR effluent) back
# into S_PO4, S_NH4, S_Mg using Shrinking Object model kinetics.
# No capital or operating cost — minor unit, no separate TSV file needed.
#
# References:
#     Aguiar 2019 (Shrinking Object model kinetics)
# =============================================================================

class StruviteRedissolution(SanUnit):
    """
    Struvite redissolution unit — dynamic-capable (AE unit).

    Dissolves residual unsettled X_struv particles (from SR effluent)
    back into S_PO4, S_NH4, S_Mg using Shrinking Object model kinetics.

    No capital or operating cost — minor unit, no separate TSV file needed.

    References:
        Aguiar 2019 (Shrinking Object model kinetics)
    """

    _N_ins  = 1
    _N_outs = 1
    _ins_size_is_fixed  = True
    _outs_size_is_fixed = True
    _neg_tol = -1e-9

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 component_ID_P='S_PO4',
                 component_ID_NH3='S_NH4',
                 component_ID_Mg='S_Mg',
                 component_ID_struvite='X_struv',
                 k_max=2.61,
                 d_p=0.3,
                 HRT_min=5.0,
                 isdynamic=True,
                 pH_ctrl=None,
                 **kwargs):

        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo,
                         init_with='WasteStream',
                         isdynamic=isdynamic,
                         **kwargs)

        self._component_ID_P        = component_ID_P
        self._component_ID_NH3      = component_ID_NH3
        self._component_ID_Mg       = component_ID_Mg
        self._component_ID_struvite = component_ID_struvite

        self.k_max   = float(k_max)
        self.d_p     = float(d_p)
        self.HRT_min = float(HRT_min)
        self.pH_ctrl = pH_ctrl

        R0           = self.d_p / 2.0
        t_diss       = R0 / self.k_max
        self._f_diss = float(np.clip(1.0 - np.exp(-self.HRT_min / t_diss), 0.0, 1.0))

    @property
    def component_ID_P(self):        return self._component_ID_P
    @property
    def component_ID_NH3(self):      return self._component_ID_NH3
    @property
    def component_ID_Mg(self):       return self._component_ID_Mg
    @property
    def component_ID_struvite(self): return self._component_ID_struvite

    def _run(self):
        influent, = self.ins
        effluent, = self.outs
        effluent.copy_like(influent)
        if 'H2O' in influent.components.IDs:
            effluent.imass['H2O'] = max(float(effluent.imass['H2O']), 1e-9)

        cmps   = influent.components
        MW_P   = cmps[self.component_ID_P].MW
        MW_NH4 = cmps[self.component_ID_NH3].MW
        MW_Mg  = cmps[self.component_ID_Mg].MW
        MW_STR = cmps[self.component_ID_struvite].MW

        m_struv_in = float(influent.imass[self.component_ID_struvite])
        if m_struv_in < self._neg_tol:
            raise ValueError(
                f'{self.ID}: negative influent struvite ({m_struv_in:g} kg/hr).')

        m_struv_diss = m_struv_in * self._f_diss
        mol_diss     = m_struv_diss / MW_STR

        effluent.imass[self.component_ID_P]        += mol_diss * MW_P
        effluent.imass[self.component_ID_NH3]      += mol_diss * MW_NH4
        effluent.imass[self.component_ID_Mg]       += mol_diss * MW_Mg
        effluent.imass[self.component_ID_struvite]  = max(m_struv_in - m_struv_diss, 0.0)

        if self.pH_ctrl is not None:
            try:
                effluent.pH  = self.pH_ctrl
                effluent._pH = self.pH_ctrl
            except AttributeError:
                pass

    def _init_state(self):
        eff = self.outs[0]
        n   = len(self.components)
        self._state      = np.empty(n + 1)
        self._state[:-1] = eff.conc
        self._state[-1]  = max(eff.F_vol * 24, 1e-12)
        self._dstate     = np.zeros(n + 1)
        if eff.state  is None: eff.state  = self._state.copy()
        if eff.dstate is None: eff.dstate = self._dstate.copy()
        eff.state[:]  = self._state
        eff.dstate[:] = self._dstate

    def _update_state(self):
        eff = self.outs[0]
        if eff.state is None: eff.state = self._state.copy()
        else:                 eff.state[:] = self._state

    def _update_dstate(self):
        eff = self.outs[0]
        if eff.dstate is None: eff.dstate = self._dstate.copy()
        else:                  eff.dstate[:] = self._dstate

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _state         = self._state
        _dstate        = self._dstate
        _update_state  = self._update_state
        _update_dstate = self._update_dstate

        cmps = self.components
        IDs  = list(cmps.IDs)

        idx_P     = IDs.index(self.component_ID_P)
        idx_NH4   = IDs.index(self.component_ID_NH3)
        idx_Mg    = IDs.index(self.component_ID_Mg)
        idx_struv = IDs.index(self.component_ID_struvite)

        MW_P   = cmps[self.component_ID_P].MW
        MW_NH4 = cmps[self.component_ID_NH3].MW
        MW_Mg  = cmps[self.component_ID_Mg].MW
        MW_STR = cmps[self.component_ID_struvite].MW

        f_diss    = self._f_diss
        f_rem     = 1.0 - f_diss
        ratio_P   = f_diss * MW_P   / MW_STR
        ratio_NH4 = f_diss * MW_NH4 / MW_STR
        ratio_Mg  = f_diss * MW_Mg  / MW_STR

        def y_t(t, y_ins, dy_ins):
            Q_in  = max(y_ins[0, -1],  1e-12)
            C_in  = y_ins[0, :-1]
            dQ_in = dy_ins[0, -1]
            dC_in = dy_ins[0, :-1]

            c_struv  = C_in[idx_struv]
            dc_struv = dC_in[idx_struv]

            C_eff             = C_in.copy()
            C_eff[idx_P]     += ratio_P   * c_struv
            C_eff[idx_NH4]   += ratio_NH4 * c_struv
            C_eff[idx_Mg]    += ratio_Mg  * c_struv
            C_eff[idx_struv]  = f_rem * c_struv

            _state[:-1] = C_eff
            _state[-1]  = Q_in

            dC_eff             = dC_in.copy()
            dC_eff[idx_P]     += ratio_P   * dc_struv
            dC_eff[idx_NH4]   += ratio_NH4 * dc_struv
            dC_eff[idx_Mg]    += ratio_Mg  * dc_struv
            dC_eff[idx_struv]  = f_rem * dc_struv

            _dstate[:-1] = dC_eff
            _dstate[-1]  = dQ_in

            _update_state()
            _update_dstate()

        self._AE = y_t

    def _design(self):
        pass  # no capital cost


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
        C['Housing'] = (self.catwalk + self.control_room_frame + 
                        self.control_room_board + self.misc_electronics)
        C['Frontend'] = self.frontend_1000_ppl_cost
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        self.add_OPEX = self._calc_replacement_cost()

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        housing_replacement_cost = (
            self.catwalk / self.catwalk_lifetime +
            self.control_room_board / self.control_room_board_lifetime +
            self.control_room_frame / self.control_room_frame_lifetime +
            self.misc_electronics / self.misc_electronics_lifetime + 
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
    """
    Hybrid photovoltaic, wind, and grid power configuration of the EnviroLoo System.

    Toggled by el.PV_ONLY (see exposan/enviroloo/__init__.py):
        False (default) -- PV+grid hybrid: 80% wind/solar, 20% grid.
            Grid cost (energy + fixed monthly charge) added to OPEX.
        True -- PV-only: 100% wind/solar, 0% grid. No grid cost in OPEX.

    - Grid cost and contribution are added to OPEX (no LCA emissions),
      only when PV_ONLY = False.
    - Electricity cost is calculated using Overstrand Municipality tariffs
      (Prepaid 60 A Single Phase, 2024–2025).

    Assumptions (PV+grid hybrid mode)
    -----------
    • Annual electricity demand = 911 kWh × 12 months = 10 932 kWh / year  
    • Grid share = 20 % → 2 186 kWh / year ≈ 182 kWh / month  
    • Fixed monthly charge = 632.15 ZAR (≈ \$37 USD @ 17.1 ZAR/USD)  
    • Variable tariff = 2.0033 ZAR / kWh (≈ \$0.117 USD / kWh)  
    • Effective blended electricity price ≈ \$0.32 USD / kWh
    """

    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    exponent_scale = 0.1
    pl = el.ppl
    baseline_ppl = el.baseline_ppl
    _D_O2 = 2.10e-9   # m2/s

    def __init__(self,
                 ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=12, W_tank=None,
                 t_wall=None, t_slab=None, aeration=None,
                 DO_ID='S_O2', suspended_growth_model=None,
                 gas_stripping=False, gas_IDs=None, stripping_kLa_min=None,
                 K_Henry=None, D_gas=None, p_gas_atm=None,
                 isdynamic=True, exogenous_vars=(), **kwargs):
        super().__init__(
            ID=ID, ins=ins, outs=outs, split=split, thermo=thermo,
            init_with=init_with, V_max=V_max, W_tank=W_tank,
            t_wall=t_wall, t_slab=t_slab, aeration=aeration,
            DO_ID=DO_ID, suspended_growth_model=suspended_growth_model,
            gas_stripping=gas_stripping, gas_IDs=gas_IDs,
            stripping_kLa_min=stripping_kLa_min,
            K_Henry=K_Henry, D_gas=D_gas, p_gas_atm=p_gas_atm,
            isdynamic=isdynamic, exogenous_vars=exogenous_vars, **kwargs
        )

        data = load_data(path=windSolar_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

        # ---------------- HYBRID / PV-ONLY TOGGLE ----------------
        # Controlled by el.PV_ONLY (see exposan/enviroloo/__init__.py).
        # True  -> 100% wind/solar, 0% grid (no grid cost added to OPEX)
        # False -> PV+grid hybrid, 80% wind/solar / 20% grid (default)
        self.total_energy_demand_kWh = 911 * 12   # Annual total demand (kWh)
        if el.PV_ONLY:
            self.grid_share       = 0.0
            self.wind_solar_share = 1.0
        else:
            self.grid_share       = 0.2
            self.wind_solar_share = 0.8

        # ---- Electricity cost assumptions (updated from Overstrand tariffs) ----
        self.grid_cost_per_kWh = 0.32              # USD per kWh (effective cost incl. fixed)
        self.grid_fixed_monthly_cost = 37          # USD/month (≈ 632 ZAR/month)
        # ------------------------------------------------------------------------

    def _init_lca(self):
        self.construction = [
            Construction(item='PhotovoltaicPanel', linked_unit=self, quantity_unit='m2'),
            Construction(item='Battery', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectricCables', linked_unit=self, quantity_unit='m'),
            Construction(item='Aluminum', linked_unit=self, quantity_unit='kg'),
            Construction(item='Steel', linked_unit=self, quantity_unit='kg')
        ]

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['PhotovoltaicPanel'] = constr[0].quantity = self.pv_photovoltaic_panel_area
        design['Battery'] = constr[1].quantity = self.el_battery_weight
        design['ElectricCables'] = constr[2].quantity = self.pv_cable_length
        design['Steel'] = constr[4].quantity = (
            self.el_wind_galvanized_steel_weight +
            self.pv_carport_weight_galvanized_metal +
            self.pv_inverter_weight_steel +
            self.pv_charger_weight_steel
        )

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Battery'] = self.el_battery
        C['PhotovoltaicPanel'] = (
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
        C['Wind'] = self.el_wind_turbine

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        # ---------------- GRID COST ADDITION ----------------
        # Skipped entirely when el.PV_ONLY = True (no grid connection at
        # all -- self.grid_share = 0.0 already, but the fixed monthly
        # charge must also be excluded since a PV-only system has no
        # grid connection to be charged for).
        if el.PV_ONLY:
            grid_cost_per_hour = 0.0
        else:
            annual_grid_energy = self.grid_share * self.total_energy_demand_kWh
            annual_grid_cost = (annual_grid_energy * self.grid_cost_per_kWh) + (self.grid_fixed_monthly_cost * 12)
            grid_cost_per_hour = annual_grid_cost / (365 * 24)  # USD/hr
        # ----------------------------------------------------

        self.add_OPEX = (
            self._calc_replacement_cost() +
            self._calc_maintenance_labor_cost() +
            grid_cost_per_hour
        )

    def _calc_replacement_cost(self):
        battery_replacement_parts_annual_cost = self.el_battery / self.el_battery_lifetime
        wind_replacement_cost = self.el_wind_turbine / self.wind_turbine_lifetime
        photovoltaic_replacement_cost = self.el_panels / self.el_panel_lifetime
        misc_replacement_cost = self.el_inverter / self.el_inverter_lifetime

        total_replacement_cost = (
            battery_replacement_parts_annual_cost +
            wind_replacement_cost +
            photovoltaic_replacement_cost +
            misc_replacement_cost
        )
        return total_replacement_cost / (365 * 24) * self.price_ratio  # USD/hr

    def _calc_maintenance_labor_cost(self):
        photovoltaic_maintenance_labor = self.pv_labor_replacement_misc_repairs * self.wages
        return photovoltaic_maintenance_labor / (365 * 24)  # USD/hr