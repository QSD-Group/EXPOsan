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
from qsdsan.sanunits import IdealClarifier, Mixer, dydt_cstr
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
# %% This callable file will be reposited to qsdsan.SanUnit subbranch with the name of _enviroloo
__all__ = (
    'EL_Excretion', # excretion
    'EL_Toilet', # toilet
    'EL_CT', # Collection tank
    'EL_PC', # Primary clarifier
    'EL_Anoxic', # Anoxic tank
    'EL_Aerobic', # Aerobic tank
    'EL_CMMBR', # Membrane filter
    'EL_CWT', # Clear water tank
    'EL_PT', # Pressure tank
    'EL_System', # System-level summary
    # 'EL_Housing', # Housing of EL_System, such as equipment's armor
    )
el_path = os.path.dirname(__file__)
module = os.path.split(el_path)[-1]
data_path, results_path = _init_modules(module, include_data_path = True)
EL_su_data_path = os.path.join(data_path, 'units_data')

# %%

excretion_path = os.path.join(EL_su_data_path, '_EL_excretion.tsv')

class EL_Excretion(ExcretionmASM2d):
    '''
    Estimation of N, P, K, and COD in urine and feces based on dietary intake
    for one person based on `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

    Parameters
    ----------
    waste_ratio : float
        A ratio in [0, 1] to indicate the amount of intake calories and nutrients
        (N, P, K) that is wasted.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_

    References
    ----------
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
    Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
    Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
    https://doi.org/10.1021/acs.est.0c03296
    '''
    
    _N_ins = 0
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 waste_ratio=0, **kwargs):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with,
                     waste_ratio=waste_ratio, **kwargs)
         # if isdyn: self._init_dynamic()
    
    # def _run(self):
    #     ur, fec = self.outs
    #     ur.empty()
    #     fec.empty()
    #     cmps = ur.components
    #     sf_iN = cmps.S_F.i_N
    #     xs_iN = cmps.X_S.i_N
    #     xb_iN = cmps.X_H.i_N
    #     sxi_iN = cmps.S_I.i_N
    #     i_mass = cmps.i_mass
    #     i_P = cmps.i_P
    #     hco3_imass = cmps.S_IC.i_mass
        
    #     # breakpoint()

    #     not_wasted = 1 - self.waste_ratio
    #     factor = 24 * 1e3 # from g/cap/d to kg/hr(/cap)
    #     e_cal = self.e_cal / 24 * not_wasted # kcal/cap/d --> kcal/cap/hr
    #     ur_exc = self.ur_exc / factor
    #     fec_exc = self.fec_exc / factor

    #     # 14 kJ/g COD, the average lower heating value of excreta
    #     tot_COD = e_cal*self.e_exc*4.184/14/1e3 # in kg COD/hr
    #     fec_COD = tot_COD*self.e_fec
    #     ur_COD = tot_COD - fec_COD
        
    #     tot_N = (self.p_veg+self.p_anim)*self.N_prot/factor \
    #         * self.N_exc*not_wasted
    #     ur_N = tot_N*self.N_ur
    #     fec_N = tot_N - ur_N
        
    #     tot_P = (self.p_veg*self.P_prot_v+self.p_anim*self.P_prot_a)/factor \
    #         * self.P_exc*not_wasted
    #     ur_P = tot_P*self.P_ur
    #     fec_P = tot_P - ur_P
        
    #     # breakpoint()
    #     ur.imass['S_NH4'] = ur_nh4 = ur_N * self.N_ur_NH3
    #     req_sf_cod = (ur_N - ur_nh4) / sf_iN
    #     if req_sf_cod <= ur_COD:
    #         ur.imass['S_F'] = sf = req_sf_cod
    #         ur.imass['S_A'] = ur_COD - sf  # contains no N or P
    #     else:
    #         req_si_cod = (ur_N - ur_nh4) / sxi_iN
    #         if req_si_cod <= ur_COD:
    #             ur.imass['S_F'] = sf = (sxi_iN * ur_COD - (ur_N - ur_nh4))/(sxi_iN - sf_iN)
    #             ur.imass['S_I'] = ur_COD - sf
    #         else:
    #             ur.imass['S_F'] = sf = ur_COD
    #             ur_other_n = ur_N - ur_nh4 - sf * sf_iN
    #             warn(f"Excess non-NH3 nitrogen cannot be accounted for by organics "
    #                  f"in urine: {ur_other_n} kg/hr. Added to NH3-N.")
    #             ur.imass['S_NH4'] += ur_other_n # debatable, has negative COD # raise warning/error
        
    #     ur.imass['S_PO4'] = ur_P - sum(ur.mass * i_P)
    #     ur.imass['S_K'] = e_cal/1e3 * self.K_cal/1e3 * self.K_exc*self.K_ur
    #     ur.imass['S_Mg'] = self.Mg_ur / factor
    #     ur.imass['S_Ca'] = self.Ca_ur / factor

    #     ur.imass['H2O'] = self.ur_moi * ur_exc
    #     ur_others = ur_exc - sum(ur.mass * i_mass)
    #     ur.imass['S_IC'] = ur_others * 0.34 / hco3_imass
    #     ur.imass['S_Na'] = ur_others * 0.35
    #     ur.imass['S_Cl'] = ur_others * 0.31

    #     fec.imass['S_NH4'] = fec_nh4 = fec_N * self.N_fec_NH3
    #     req_xs_cod = (fec_N - fec_nh4) / xs_iN
    #     if req_xs_cod <= fec_COD:
    #         fec.imass['X_S'] = xs = req_xs_cod
    #         fec.imass['S_A'] = fec_COD - xs
    #     else:
    #         req_xi_cod = (fec_N - fec_nh4) / sxi_iN
    #         if req_xi_cod <= fec_COD:
    #             fec.imass['X_S'] = xs = (sxi_iN * fec_COD - (fec_N - fec_nh4))/(sxi_iN - xs_iN)
    #             fec.imass['X_I'] = fec_COD - xs
    #         else:
    #             req_xb_cod = (fec_N - fec_nh4) / xb_iN
    #             if req_xb_cod <= fec_COD:
    #                 fec.imass['X_S'] = xs = (xb_iN * fec_COD - (fec_N - fec_nh4))/(xb_iN - xs_iN)
    #                 fec.imass['X_H'] = fec_COD - xs
    #             else:
    #                 fec.imass['X_S'] = xs = fec_COD
    #                 fec_other_n = fec_N - fec_nh4 - xs * xs_iN
    #                 warn(f"Excess non-NH3 nitrogen cannot be accounted for by organics "
    #                      f"in feces: {fec_other_n} kg/hr. Added to NH3-N.")
    #                 fec.imass['S_NH4'] += fec_other_n # debatable, has negative COD
        
    #     fec.imass['S_PO4'] = fec_P - sum(fec.mass * i_P)
    #     fec.imass['S_K'] = (1-self.K_ur)/self.K_ur * ur.imass['S_K']
    #     fec.imass['S_Mg'] = self.Mg_fec / factor
    #     fec.imass['S_Ca'] = self.Ca_fec / factor
    #     fec.imass['H2O'] = self.fec_moi * fec_exc
        
    #     fec_others = fec_exc - sum(fec.mass * i_mass)
    #     fec.imass['S_IC'] = fec_others * 0.34 / hco3_imass
    #     fec.imass['S_Na'] = fec_others * 0.35
    #     fec.imass['S_Cl'] = fec_others * 0.31
        
    
    # @property
    # def AE(self):
    #     if self._AE is None:
    #         self._compile_AE()
    #     return self._AE
    
    # def _compile_AE(self):
    #     def yt(t, QC_ins, dQC_ins):
    #         pass
    #     self._AE = yt
        
    # def _init_state(self):
    #     ur, fec = self.outs
    #     self._state = np.append(ur.mass, fec.mass)
    #     for ws in self.outs:
    #         ws.state = np.append(ws.conc, ws.F_vol * 24)
    #         ws.dstate = np.zeros_like(ws.state)

    # def _update_state(self):
    #     pass

    # def _update_dstate(self):
    #     pass


# %%

toilet_path = os.path.join(EL_su_data_path, '_EL_Toilet.tsv')

@price_ratio()
class EL_Toilet(Toilet):
    _N_ins = 6
    _outs_size_is_fixed = False
    density_dct = {
        'Sand': 1442,
        'Gravel': 1600,
        'Brick': 1750,
        'Plastic': 0.63,
        'Steel': 7900,
        'StainlessSteelSheet': 2.64
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), N_user=100, N_toilet=1, N_tot_user=None,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=True,
                 if_desiccant=True, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=None, OPEX_over_CAPEX=None, price_ratio=1.):
                 # F_BM_default=1):
        super().__init__(
            ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with,
            degraded_components=degraded_components, N_user=N_user, N_toilet=N_toilet, N_tot_user=N_tot_user,
            if_toilet_paper=if_toilet_paper, if_flushing=if_flushing, if_cleansing=if_cleansing,
            if_desiccant=if_desiccant, if_air_emission=if_air_emission, if_ideal_emptying=if_ideal_emptying,
            CAPEX=CAPEX, OPEX_over_CAPEX=OPEX_over_CAPEX, price_ratio=price_ratio, 
            )
                     # F_BM_default=F_BM_default,)

    #     Toilet.__init__(self, ID, ins, outs, thermo, init_with)
    #     self.degraded_components = tuple(degraded_components)
    #     self._N_user = self._N_toilet = self._N_tot_user = None
    #     self.N_user = N_user
    #     self.N_toilet = N_toilet
    #     self.N_tot_user = N_tot_user
    #     self.if_toilet_paper = if_toilet_paper
    #     self.if_flushing = if_flushing
    #     self.if_cleansing = if_cleansing
    #     self.if_desiccant = if_desiccant
    #     self.if_air_emission = if_air_emission
    #     self.if_ideal_emptying = if_ideal_emptying
    #     self.CAPEX = CAPEX
    #     self.OPEX_over_CAPEX = OPEX_over_CAPEX
    #     self.price_ratio = price_ratio

    #     data = load_data(path=toilet_path)
    #     for para in data.index:
    #         value = float(data.loc[para]['expected'])
    #         if para in ('desiccant_V', 'desiccant_rho'):
    #             setattr(self, para, value)
    #         else:
    #             setattr(self, '_'+para, value)
    #     del data

    #     self._empty_ratio = 0.59


    # # def _run(self):
    # #     ur, fec, tp, fw, cw, des = self.ins
    # #     tp.imass['Tissue'] = int(self.if_toilet_paper)*self.toilet_paper
    # #     fw.imass['H2O'] = int(self.if_flushing)*self.flushing_water
    # #     cw.imass['H2O'] = int(self.if_cleansing)*self.cleansing_water
    # #     des.imass['WoodAsh'] = int(self.if_desiccant)*self.desiccant
        
        
    # def _run(self):
    #     ur, fec, tp, fw, cw, des = self.ins
    #     tp.imass['Tissue'] = int(self.if_toilet_paper)*self.toilet_paper
    #     fw.imass['H2O'] = int(self.if_flushing)*self.flushing_water
    #     cw.imass['H2O'] = int(self.if_cleansing)*self.cleansing_water
    #     des.imass['WoodAsh'] = int(self.if_desiccant)*self.desiccant

    #     mw = self.outs[0]

    #     # Mix all relevant inputs into the mixed_waste output
    #     mw.mix_from([ur, fec, tp, fw, cw, des])

    #     # if self.if_air_emission:
    #     #   self.get_emptying_emission(
    #     #       mw,
    #     #       self.empty_ratio,
    #     #       self.MCF_aq,
    #     #       self.N2O_EF_aq
    #     # )
        
    # def _init_state(self):
    #     pass

    # def _update_state(self):
    #     pass

    # def _update_dstate(self):
    #     pass

    # def _scale_up_outs(self):
    #     '''
    #     Scale up the effluent based on the number of user per toilet and
    #     toilet number.
    #     '''
    #     N_tot_user = self.N_tot_user or self.N_toilet*self.N_user
    #     for i in self.outs:
    #         if not i.F_mass == 0:
    #             i.F_mass *= N_tot_user


    # def _cost(self):
    #     self.baseline_purchase_costs['Total toilets'] = self.CAPEX * self.N_toilet * self.price_ratio
    #     add_OPEX = self.baseline_purchase_costs['Total toilets']*self.OPEX_over_CAPEX/365/24
    #     self._add_OPEX = {'Additional OPEX': add_OPEX}


    # @staticmethod
    # def get_emptying_emission(waste, CH4, N2O, empty_ratio, CH4_factor, N2O_factor):
    #     '''
    #     Calculate emissions due to non-ideal emptying based on
    #     `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_,

    #     Parameters
    #     ----------
    #     stream : WasteStream
    #         Excreta stream that is not appropriately emptied (before emptying).
    #     CH4 : WasteStream
    #         Fugitive CH4 gas (before emptying).
    #     N2O : WasteStream
    #         Fugitive N2O gas (before emptying).
    #     empty_ratio : float
    #         Fraction of excreta that is appropriately emptied..
    #     CH4_factor : float
    #         Factor to convert COD removal to CH4 emission.
    #     N2O_factor : float
    #         Factor to convert COD removal to N2O emission.

    #     Returns
    #     -------
    #     stream : WasteStream
    #         Excreta stream that is not appropriately emptied (after emptying).
    #     CH4 : WasteStream
    #         Fugitive CH4 gas (after emptying).
    #     N2O : WasteStream
    #         Fugitive N2O gas (after emptying).
    #     '''
    #     COD_rmvd = waste.COD*(1-empty_ratio)/1e3*waste.F_vol
    #     CH4.imass['CH4'] += COD_rmvd * CH4_factor
    #     N2O.imass['N2O'] += COD_rmvd * N2O_factor
    #     waste.mass *= empty_ratio

    # @property
    # def N_user(self):
    #     '''[int, float] Number of people per toilet.'''
    #     return self._N_user or self.N_tot_user/self.N_toilet
    # @N_user.setter
    # def N_user(self, i):
    #     if i is not None:
    #         N_user = self._N_user = int(i)
    #         old_toilet = self._N_toilet
    #         if old_toilet and self.N_tot_user:
    #             new_toilet = ceil(self.N_tot_user/N_user)
    #             warn(f'With the provided `N_user`, the previous `N_toilet` of {old_toilet} '
    #                  f'is recalculated from `N_tot_user` and `N_user` as {new_toilet}.')
    #             self._N_toilet = None
    #     else:
    #         self._N_user = i

    # @property
    # def N_toilet(self):
    #     '''[int] Number of parallel toilets.'''
    #     return self._N_toilet or ceil(self.N_tot_user/self.N_user)
    # @N_toilet.setter
    # def N_toilet(self, i):
    #     if i is not None:
    #         N_toilet = self._N_toilet = ceil(i)
    #         old_user = self._N_user
    #         if old_user and self.N_tot_user:
    #             new_user = self.N_tot_user/N_toilet
    #             warn(f'With the provided `N_toilet`, the previous `N_user` of {old_user} '
    #                  f'is recalculated from `N_tot_user` and `N_toilet` as {new_user}.')
    #             self._N_user = None
    #     else:
    #         self._N_toilet = i

    # @property
    # def N_tot_user(self):
    #     '''[int] Number of total users.'''
    #     return self._N_tot_user
    # @N_tot_user.setter
    # def N_tot_user(self, i):
    #     if i is not None:
    #         self._N_tot_user = int(i)
    #     else:
    #         self._N_tot_user = None

    # @property
    # def toilet_paper(self):
    #     '''
    #     [float] Amount of toilet paper used
    #     (if ``if_toilet_paper`` is True), [kg/cap/hr].
    #     '''
    #     return self._toilet_paper
    # @toilet_paper.setter
    # def toilet_paper(self, i):
    #     self._toilet_paper = i

    # @property
    # def flushing_water(self):
    #     '''
    #     [float] Amount of water used for flushing
    #     (if ``if_flushing_water`` is True), [kg/cap/hr].
    #     '''
    #     return self._flushing_water
    # @flushing_water.setter
    # def flushing_water(self, i):
    #     self._flushing_water = i

    # @property
    # def cleansing_water(self):
    #     '''
    #     [float] Amount of water used for cleansing
    #     (if ``if_cleansing_water`` is True), [kg/cap/hr].
    #     '''
    #     return self._cleansing_water
    # @cleansing_water.setter
    # def cleansing_water(self, i):
    #     self._cleansing_water = i

    # @property
    # def desiccant(self):
    #     '''
    #     [float] Amount of desiccant used (if ``if_desiccant`` is True), [kg/cap/hr].

    #     .. note::

    #         Value set by ``desiccant_V`` and ``desiccant_rho``.

    #     '''
    #     return self.desiccant_V*self.desiccant_rho

    # @property
    # def N_volatilization(self):
    #     '''
    #     [float] Fraction of input N that volatilizes to the air
    #     (if ``if_air_emission`` is True).
    #     '''
    #     return self._N_volatilization
    # @N_volatilization.setter
    # def N_volatilization(self, i):
    #     self._N_volatilization = i

    # @property
    # def empty_ratio(self):
    #     '''
    #     [float] Fraction of excreta that is appropriately emptied.

    #     .. note::

    #         Will be 1 (i.e., 100%) if ``if_ideal_emptying`` is True.

    #     '''
    #     if self.if_ideal_emptying:
    #         return 1.
    #     return self._empty_ratio
    # @empty_ratio.setter
    # def empty_ratio(self, i):
    #     if self.if_ideal_emptying:
    #         warn(f'`if_ideal_emptying` is True, the set value {i} is ignored.')
    #     self._empty_ratio = i

    # @property
    # def MCF_aq(self):
    #     '''[float] Methane correction factor for COD lost due to inappropriate emptying.'''
    #     return self._MCF_aq
    # @MCF_aq.setter
    # def MCF_aq(self, i):
    #     self._MCF_aq = i

    # @property
    # def N2O_EF_aq(self):
    #     '''[float] Fraction of N emitted as N2O due to inappropriate emptying.'''
    #     return self._N2O_EF_aq
    # @N2O_EF_aq.setter
    # def N2O_EF_aq(self, i):
    #     self._N2O_EF_aq = i

    # @property
    # def if_N2O_emission(self):
    #     '''[bool] Whether to consider N degradation and fugitive N2O emission.'''
    #     return self.if_air_emission
    # @if_N2O_emission.setter
    # def if_N2O_emission(self, i):
    #     raise ValueError('Setting `if_N2O_emission` for `PitLatrine` is not supported, '
    #                      'please set `if_air_emission` instead.')





# %%

CollectionTank_path = os.path.join(EL_su_data_path, '_EL_CT.tsv')

@price_ratio()
class EL_CT(Mixer):
    
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
    # _graphics = BSTMixer._graphics
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 rigorous=False, conserve_phases=False):
        super().__init__(
            ID=ID, ins=ins, outs=outs, thermo=thermo,
            init_with=init_with, F_BM_default=F_BM_default, isdynamic=isdynamic,
            rigorous=rigorous, conserve_phases=conserve_phases
            )
         # self.rigorous = rigorous
         # self.conserve_phases = conserve_phases


    # @property
    # def state(self):
    #     '''The state of the Mixer, including component concentrations [mg/L] and flow rate [m^3/d].'''
    #     if self._state is None: return None
    #     else:
    #         return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    # def _init_state(self):
    #     '''initialize state by specifying or calculating component concentrations
    #     based on influents. Total flow rate is always initialized as the sum of
    #     influent wastestream flows.'''
    #     QCs = self._ins_QC
        
    #     if QCs.shape[0] <= 1: self._state = QCs[0]
    #     else:
    #         Qs = QCs[:,-1]
    #         # breakpoint()
    #         Cs = QCs[:,:-1]
            
    #         self._state = np.append(Qs @ Cs / Qs.sum(), Qs.sum())
    #     self._dstate = self._state * 0.

    # def _update_state(self):
    #     '''updates conditions of output stream based on conditions of the Mixer'''
    #     self._outs[0].state = self._state

    # def _update_dstate(self):
    #     '''updates rates of change of output stream from rates of change of the Mixer'''
    #     self._outs[0].dstate = self._dstate

    # @property
    # def AE(self):
    #     if self._AE is None:
    #         self._compile_AE()
    #     return self._AE

    # def _compile_AE(self):
    #     _n_ins = len(self.ins)
    #     _state = self._state
    #     _dstate = self._dstate
    #     _update_state = self._update_state
    #     _update_dstate = self._update_dstate
    #     def yt(t, QC_ins, dQC_ins):
    #         if _n_ins > 1:
    #             Q_ins = QC_ins[:, -1]
    #             C_ins = QC_ins[:, :-1]
    #             dQ_ins = dQC_ins[:, -1]
    #             dC_ins = dQC_ins[:, :-1]
    #             Q = Q_ins.sum()
    #             C = Q_ins @ C_ins / Q
    #             _state[-1] = Q
    #             _state[:-1] = C
    #             Q_dot = dQ_ins.sum()
    #             C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
    #             _dstate[-1] = Q_dot
    #             _dstate[:-1] = C_dot
    #         else:
    #             _state[:] = QC_ins[0]
    #             _dstate[:] = dQC_ins[0]
    #         _update_state()
    #         _update_dstate()
    #     self._AE = yt

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.collection_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
    
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
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime) * scale
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

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 sludge_flow_rate=10, 
                 solids_removal_efficiency=.85,
                 sludge_MLSS=None, isdynamic=False, init_with='WasteStream',
                 F_BM_default=None, **kwargs):

        super().__init__(
            ID=ID, ins=ins, outs=outs, thermo=thermo,
            sludge_flow_rate=sludge_flow_rate, 
            solids_removal_efficiency=solids_removal_efficiency,
            sludge_MLSS=sludge_MLSS, isdynamic=isdynamic, init_with=init_with,
            F_BM_default=F_BM_default, **kwargs
            )
        
        # self.sludge_flow_rate = sludge_flow_rate
        # self.solids_removal_efficiency = solids_removal_efficiency
        # self.sludge_MLSS = sludge_MLSS
        # self._mixed = WasteStream()
        # self._f_uf = None
        # self._f_of = None

    # @property
    # def sludge_flow_rate(self):
    #     '''[float] The designed sludge flow rate (wasted + recycled) in m3/d.'''
    #     return self._Qs

    # @sludge_flow_rate.setter
    # def sludge_flow_rate(self, Qs):
    #     self._Qs = Qs

    # @property
    # def solids_removal_efficiency(self):
    #     return self._e_rmv

    # @solids_removal_efficiency.setter
    # def solids_removal_efficiency(self, f):
    #     if f is not None and (f > 1 or f < 0):
    #         raise ValueError(f'solids removal efficiency must be within [0, 1], not {f}')
    #     self._e_rmv = f

    # @property
    # def sludge_MLSS(self):
    #     return self._MLSS

    # @sludge_MLSS.setter
    # def sludge_MLSS(self, MLSS):
    #     if MLSS is not None:
    #         warn(f'sludge MLSS {MLSS} mg/L is only used to estimate '
    #              f'sludge flowrate or solids removal efficiency, when either '
    #              f'one of them is unspecified.')
    #     self._MLSS = MLSS

    # def _run(self):
    #     inf = self._mixed
    #     inf.mix_from(self.ins)
    #     of, uf = self.outs
    #     TSS_in = inf.get_TSS()
    #     if TSS_in <= 0:
    #         uf.empty()
    #         of.copy_like(inf)
    #     else:
    #         Q_in = inf.F_vol * 24 # m3/d
    #         x = inf.components.x
    #         Qs, e_rmv, mlss = self._Qs, self._e_rmv, self._MLSS
    #         if Qs and e_rmv:
    #             f_Qu = Qs/Q_in
    #             f_Xu = e_rmv + (1-e_rmv) * f_Qu
    #         elif Qs and mlss:
    #             f_Qu = Qs/Q_in
    #             f_Xu = f_Qu*mlss/TSS_in
    #         elif e_rmv and mlss:
    #             f_Qu = e_rmv / (mlss/TSS_in - (1-e_rmv))
    #             f_Xu = e_rmv + (1-e_rmv) * f_Qu
    #         split_to_uf = (1-x)*f_Qu + x*f_Xu
    #         if any(split_to_uf > 1): split_to_uf = 1
    #         inf.split_to(uf, of, split_to_uf)

    # def _init_state(self):
    #     inf = self._mixed
    #     C_in = inf.conc
    #     Q_in = inf.F_vol * 24
    #     self._state = np.append(C_in, Q_in)
    #     self._dstate = self._state * 0.
        
    # def _update_state(self):
    #     arr = self._state
    #     Cs = arr[:-1]
    #     Qi = arr[-1]
    #     Qs, e_rmv, mlss = self._Qs, self._e_rmv, self._MLSS
    #     x = self.components.x
    #     i_tss = x * self.components.i_mass

    #     of, uf = self.outs
    #     if uf.state is None: uf.state = np.zeros(len(x)+1)
    #     if of.state is None: of.state = np.zeros(len(x)+1)

    #     if Qs:
    #         Qe = Qi - Qs
    #         if e_rmv:
    #             fuf = e_rmv * Qi/Qs + (1-e_rmv)
    #             fof = 1-e_rmv
    #         elif mlss:
    #             tss_in = sum(Cs * i_tss)
    #             tss_e = (Qi * tss_in - Qs * mlss)/Qe
    #             fuf = mlss/tss_in
    #             fof = tss_e/tss_in
    #     elif e_rmv and mlss:
    #         tss_in = sum(Cs * i_tss)
    #         Qs = Qi * e_rmv / (mlss/tss_in - (1-e_rmv))
    #         Qe = Qi - Qs
    #         fuf = mlss/tss_in
    #         fof = 1-e_rmv
    #     else:
    #         raise RuntimeError('missing parameter')
            
    #     if Qs >= Qi: 
    #         uf.state[:] = arr
    #         of.state[:] = 0.
    #     else:
    #         self._f_uf = fuf
    #         self._f_of = fof
    #         uf.state[:-1] = Cs * ((1-x) + x*fuf)
    #         uf.state[-1] = Qs
    #         of.state[:-1] = Cs * ((1-x) + x*fof)
    #         of.state[-1] = Qe

    # def _update_dstate(self):
    #     of, uf = self.outs
    #     x = self.components.x
    #     if uf.dstate is None: uf.dstate = np.zeros(len(x)+1)
    #     if of.dstate is None: of.dstate = np.zeros(len(x)+1)
    
    # @property
    # def AE(self):
    #     if self._AE is None:
    #         self._compile_AE()
    #     return self._AE

    # def _compile_AE(self):        
    #     _state = self._state
    #     # _dstate = self._dstate
    #     _update_state = self._update_state
    #     # _update_dstate = self._update_dstate
    #     def yt(t, QC_ins, dQC_ins):
    #         Q_ins = QC_ins[:, -1]
    #         C_ins = QC_ins[:, :-1]
    #         # dQ_ins = dQC_ins[:, -1]
    #         # dC_ins = dQC_ins[:, :-1]
    #         Q = Q_ins.sum()
    #         C = Q_ins @ C_ins / Q
    #         _state[-1] = Q
    #         _state[:-1] = C
    #         # Q_dot = dQ_ins.sum()
    #         # C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
    #         # _dstate[-1] = Q_dot
    #         # _dstate[:-1] = C_dot
    #         _update_state()
    #         # _update_dstate()
    #     self._AE = yt
 

    # def _design(self):
    #     """Calculate design parameters."""
    #     self.design_results['StainlessSteel'] = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
    #     self.construction = [
    #         Construction(item='StainlessSteel', quantity=self.design_results['StainlessSteel'], quantity_unit='kg'),
    #     ]
    #     self.add_construction(add_cost=False)
        ####### These below codes are essential to perform LCA and TEA #########
        data = load_data(path = PrimaryClarifier_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)    
        ###############################################
    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]  
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost=False)

    def _cost(self):
        """Calculate capital and operating costs."""
        C = self.baseline_purchase_costs
        C['Tank'] = self.PC_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings

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
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime
        ) * scale
        return replacement_cost / (365 * 24)  # Convert to USD/hr

# %%


Anoxic_path = os.path.join(EL_su_data_path, '_EL_Anoxic.tsv')

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
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    ppl = 100
    baseline_ppl =30
    
    # _D_O2 = 2.10e-9   # m2/s

    def __init__(self, ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=7.3, W_tank = 2.09, 
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
        # self._D_tank = D_tank
        # self._freeboard = freeboard
        # self._t_wall = t_wall
        # self._t_slab = t_slab
    
        # for attr, value in kwargs.items():
        #     setattr(self, attr, value)
        data = load_data(path = Anoxic_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)    
        ###############################################
    
    
         
    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]      
        

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # assume linear scale
        self.add_construction(add_cost=False)
    
    # def _cost(self):
    #     C = self.baseline_purchase_costs
    #     massflow_anoxic = self.ins[0].mass
    #     C['Tank'] = self.anoxic_tank_cost
    #     C['Pipes'] = self.pipeline_connectors
    #     C['Fittings'] = self.weld_female_adapter_fittings
    #     C['Chemcial_glucose'] = self.chemical_glucose_dosage * massflow_anoxic * self.chemical_glucose_price  # make sense the unit of treated water flow

    #     ratio = self.price_ratio
    #     for equipment, cost in C.items():
    #         C[equipment] = cost * ratio
        
    #     self.add_OPEX = self._calc_replacement_cost()
        
    #     power_demand = self.power_demand_AnoxicTank
    #     self.power_utility(power_demand)
    
    # def _calc_replacement_cost(self):
    #     scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
    #     Anoxic_tank_replacement_cost = (self.anoxic_tank_cost /self.anoxic_tank_lifetime +
    #                                     self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
    #                                     self.pipeline_connectors / self.pipeline_connectors_lifetime) * scale
    #     Anoxic_tank_replacement_cost = Anoxic_tank_replacement_cost / (365 * 24)  # convert to USD/hr
    #     return Anoxic_tank_replacement_cost

# %%


Aerobic_path = os.path.join(EL_su_data_path, '_EL_Aerobic.tsv')

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
    _N_ins = 1 # treated water, PAC, blower
    _N_outs = 3  # treated water, CH4, N2O
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    # exponent_scale = 0.1
    ppl = 100
    baseline_ppl =30
    
    _D_O2 = 2.10e-9   # m2/s

    def __init__(self, 
                 ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=7.3, W_tank = 2.09, 
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
    
    
        

    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]  


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost=False)
    
    # def _cost(self):
    #     C = self.baseline_purchase_costs
    #     massflow_aerobic = self.ins[0].mass
    #     C['Tank'] = self.aerobic_tank_cost
    #     C['Pipes'] = self.pipeline_connectors
    #     C['Fittings'] = self.weld_female_adapter_fittings
    #     C['Chemical_PAC'] = self.chemical_PAC_dosage * massflow_aerobic * self.chemical_PAC_price

    #     ratio = self.price_ratio
    #     for equipment, cost in C.items():
    #         C[equipment] = cost * ratio
        
    #     self.add_OPEX = self._calc_replacement_cost()
        
    #     power_demand = self.power_demand_AerobicTank
    #     self.power_utility(power_demand) # kWh
    
    # def _calc_replacement_cost(self):
    #     scale = (self.ppl / self.baseline_ppl) * self.exponent_scale
    #     Aerobic_tank_replacement_cost = (self.aerobic_tank_cost / self.aerobic_tank_lifetime +
    #                                     self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
    #                                     self.pipeline_connectors / self.pipeline_connectors_lifetime) * scale
    #     Aerobic_tank_replacement_cost = Aerobic_tank_replacement_cost / (365 * 24)  # convert to USD/hr
    #     return Aerobic_tank_replacement_cost

# %%

MBR_path = os.path.join(EL_su_data_path, '_EL_MBR.tsv')

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
    ppl = 100
    baseline_ppl =30
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', isdynamic=True, 
                 aeration=2, DO_ID='S_O2', suspended_growth_model=None,
                 pumped_flow=0.103, 
                 solids_capture_rate=0.999, 
                 V_max=3.3, 
                 crossflow_air=None,
                 **kwargs):
        super().__init__(
            ID=ID, ins=ins, outs=outs, thermo=thermo,
            init_with=init_with, V_max=V_max, aeration=aeration, DO_ID=DO_ID,
            suspended_growth_model=suspended_growth_model, isdynamic=isdynamic, 
            pumped_flow=pumped_flow, 
            solids_capture_rate=solids_capture_rate,  
            crossflow_air=crossflow_air,
            **kwargs
            )
    
    # @property
    # def pumped_flow(self):
    #     '''[float] Pumped flow rate, in m3/d'''
    #     return self._Q_pump
    # @pumped_flow.setter
    # def pumped_flow(self, Q):
    #     self._Q_pump = Q
    
    # @property
    # def solids_capture_rate(self):
    #     '''[float] Membrane solid capture rate, i.e., 
    #     filtrate-to-internal solids concentration ratio, unitless.'''
    #     return self._f_rtn
    # @solids_capture_rate.setter
    # def solids_capture_rate(self, f):
    #     if f < 0 or f > 1:
    #         raise ValueError(f'membrane solids capture rate must be within [0,1], not {f}')
    #     self._f_rtn = f
    #     cmps = self._mixed.components
    #     self._flt2in_conc_ratio = (1-cmps.x) + cmps.x * (1-f)
    
    # @property
    # def crossflow_air(self):
    #     '''[:class:`qsdsan.Process` or NoneType]
    #     Membrane cross flow air specified for process modeling, such as `qsdsan.processes.DiffusedAeration`. 
    #     Ignored if DO setpoint is specified by the `aeration` attribute.
    #     '''
    #     return self._cfa
    # @crossflow_air.setter
    # def crossflow_air(self, cfa):
    #     if cfa is None or isinstance(cfa, Process):
    #         self._cfa = cfa
    #     else:
    #         raise TypeError('crossflow_air must be a `Process` object or None, '
    #                         f'not {type(cfa)}')

    # split = None
        
    # def _run(self):
    #     '''Only to converge volumetric flows.'''
    #     mixed = self._mixed
    #     mixed.mix_from(self.ins)
    #     cmps = mixed.components
    #     Q = mixed.F_vol*24 # m3/d
    #     Qp = self._Q_pump
    #     f_rtn = self._f_rtn
    #     xsplit = Qp / ((1-f_rtn)*(Q-Qp) + Qp) # mass split of solids to pumped flow
    #     qsplit = Qp / Q
    #     flt, rtn = self.outs
    #     mixed.split_to(rtn, flt, xsplit*cmps.x + qsplit*(1-cmps.x))
    
    # def _compile_ODE(self):
    #     aer = self._aeration
    #     cfa = self._cfa
    #     isa = isinstance
    #     cmps = self.components
    #     if self._model is None:
    #         warn(f'{self.ID} was initialized without a suspended growth model, '
    #              f'and thus run as a non-reactive unit')
    #         r = lambda state_arr: np.zeros(cmps.size)
    #     else:
    #         r = self._model.production_rates_eval

    #     _dstate = self._dstate
    #     _update_dstate = self._update_dstate
    #     V = self._V_max
    #     Qp = self.pumped_flow
    #     f_rtn = self.solids_capture_rate
    #     xarr = cmps.x
    #     gstrip = self.gas_stripping
    #     if gstrip:
    #         gas_idx = self.components.indices(self.gas_IDs)
    #         if isa(aer, Process): kLa = aer.kLa
    #         else: kLa = 0.
    #         if cfa: kLa += cfa.kLa
    #         S_gas_air = np.asarray(self.K_Henry)*np.asarray(self.p_gas_atm)
    #         kLa_stripping = np.maximum(kLa*self.D_gas/self._D_O2, self.stripping_kLa_min)
    #     hasexo = bool(len(self._exovars))
    #     f_exovars = self.eval_exo_dynamic_vars
         
    #     if isa(aer, (float, int)):
    #         i = cmps.index(self._DO_ID)
    #         def dy_dt(t, QC_ins, QC, dQC_ins):
    #             QC[i] = aer
    #             dydt_mbr(QC_ins, QC, V, Qp, f_rtn, xarr, _dstate)
    #             if hasexo: QC = np.append(QC, f_exovars(t))
    #             _dstate[:-1] += r(QC)
    #             if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
    #             _dstate[i] = 0
    #             _update_dstate()
    #     else:        
    #         if cfa:
    #             cfa_stoi = cfa._stoichiometry
    #             cfa_frho = cfa.rate_function
    #             dy_cfa = lambda QC: cfa_stoi * cfa_frho(QC)
    #         else:
    #             dy_cfa = lambda QC: 0.
            
    #         if isa(aer, Process):
    #             aer_stoi = aer._stoichiometry
    #             aer_frho = aer.rate_function
    #             dy_aer = lambda QC: aer_stoi * aer_frho(QC)
    #         else:
    #             dy_aer = lambda QC: 0.
                
    #         def dy_dt(t, QC_ins, QC, dQC_ins):
    #             dydt_mbr(QC_ins, QC, V, Qp, f_rtn, xarr, _dstate)
    #             if hasexo: QC = np.append(QC, f_exovars(t))
    #             _dstate[:-1] += r(QC) + dy_aer(QC) + dy_cfa(QC)
    #             if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
    #             _update_dstate()
    #     self._ODE = dy_dt
    
    # def _update_state(self):
    #     arr = self._state
    #     arr[arr < 1e-16] = 0.
    #     arr[-1] = sum(ws.state[-1] for ws in self.ins)
    #     for ws in self.outs:
    #         if ws.state is None: 
    #             ws.state = np.zeros_like(arr)
    #             ws.dstate = np.zeros_like(arr)
    #     flt, rtn = self.outs
    #     Qp = self.pumped_flow
    #     flt.state[:-1] = arr[:-1] * self._flt2in_conc_ratio
    #     flt.state[-1] = arr[-1] - Qp
    #     rtn.state[:-1] = arr[:-1]
    #     rtn.state[-1] = Qp
        
    # def _update_dstate(self):
    #     arr = self._dstate
    #     arr[-1] = sum(ws.dstate[-1] for ws in self.ins)
    #     flt, rtn = self.outs
    #     flt.dstate[:-1] = arr[:-1] * self._flt2in_conc_ratio
    #     flt.dstate[-1] = arr[-1]
    #     rtn.dstate[:-1] = arr[:-1]
    #     rtn.dstate[-1] = 0
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
                             Construction(item ='PVDF_membrane', linked_unit=self, quantity_unit='kg'),]
    

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # assume linear scaling TODO: scale to mass not volume
        design['PVDF_membrane'] = constr[1].quantity = self.membrane_material_weight
        self.add_construction(add_cost=False)

    # def _cost(self):
    #     C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
    #     C['MBR_tank'] = self.MBR_tank_cost
    #     C['pipeline'] = self.pipeline_connectors
    #     C['fittings'] = self.weld_female_adapter_fittings
    #     C['Membrane_material'] = self.membrane_material_price * self.membrane_material_weight
    #     C['Membrane_cleaning'] = self.membrane_cleaning_fee

    #     ratio = self.price_ratio
    #     for equipment, cost in C.items():
    #         C[equipment] = cost * ratio 
        
    #     self.add_OPEX = self._calc_replacement_cost()
        
    #     power_demand = self.power_demand_MBR
    #     self.power_utility(power_demand)

    # def _calc_replacement_cost(self):
    #     scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
    #     MBR_replacement_cost = (
    #         self.MBR_tank_cost / self.MBR_tank_lifetime +
    #         self.pipeline_connectors / self.pipeline_connectors_lifetime +
    #         self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
    #         self.membrane_material_price * self.membrane_material_weight / self.membrane_material_lifetime
    #         ) * scale
    #     MBR_replacement_cost = MBR_replacement_cost / (365 * 24) * self.price_ratio # USD/hr
    #     return MBR_replacement_cost

# %%
ClearWaterTank_path = os.path.join(EL_su_data_path, '_EL_CWT.tsv')

@price_ratio()
class EL_CWT(Mixer):

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
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 rigorous=False, conserve_phases=False):
        super().__init__(
            ID=ID, ins=ins, outs=outs, thermo=thermo,
            init_with=init_with, F_BM_default=F_BM_default, isdynamic=isdynamic,
            rigorous=rigorous, conserve_phases=conserve_phases)
         # self.rigorous = rigorous
         # self.conserve_phases = conserve_phases


    # @property
    # def state(self):
    #     '''The state of the Mixer, including component concentrations [mg/L] and flow rate [m^3/d].'''
    #     if self._state is None: return None
    #     else:
    #         return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    # def _init_state(self):
    #     '''initialize state by specifying or calculating component concentrations
    #     based on influents. Total flow rate is always initialized as the sum of
    #     influent wastestream flows.'''
    #     QCs = self._ins_QC
        
    #     if QCs.shape[0] <= 1: self._state = QCs[0]
    #     else:
    #         Qs = QCs[:,-1]
    #         # breakpoint()
    #         Cs = QCs[:,:-1]
            
    #         self._state = np.append(Qs @ Cs / Qs.sum(), Qs.sum())
    #     self._dstate = self._state * 0.

    # def _update_state(self):
    #     '''updates conditions of output stream based on conditions of the Mixer'''
    #     self._outs[0].state = self._state

    # def _update_dstate(self):
    #     '''updates rates of change of output stream from rates of change of the Mixer'''
    #     self._outs[0].dstate = self._dstate

    # @property
    # def AE(self):
    #     if self._AE is None:
    #         self._compile_AE()
    #     return self._AE

    # def _compile_AE(self):
    #     _n_ins = len(self.ins)
    #     _state = self._state
    #     _dstate = self._dstate
    #     _update_state = self._update_state
    #     _update_dstate = self._update_dstate
    #     def yt(t, QC_ins, dQC_ins):
    #         if _n_ins > 1:
    #             Q_ins = QC_ins[:, -1]
    #             C_ins = QC_ins[:, :-1]
    #             dQ_ins = dQC_ins[:, -1]
    #             dC_ins = dQC_ins[:, :-1]
    #             Q = Q_ins.sum()
    #             C = Q_ins @ C_ins / Q
    #             _state[-1] = Q
    #             _state[:-1] = C
    #             Q_dot = dQ_ins.sum()
    #             C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
    #             _dstate[-1] = Q_dot
    #             _dstate[:-1] = C_dot
    #         else:
    #             _state[:] = QC_ins[0]
    #             _dstate[:] = dQC_ins[0]
    #         _update_state()
    #         _update_dstate()
    #     self._AE = yt

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.collection_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
    
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
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime) * scale
        CT_replacement_cost = CT_replacement_cost / (365 * 24)  # convert to USD/hr
        return CT_replacement_cost

# %%
PressureTank_path = os.path.join(EL_su_data_path, '_EL_PT.tsv')

@price_ratio()
class EL_PT(StorageTank):
    
    '''
    Introduction
    ------------
    To only collect the clear water within pressure tank

    Parameters
    ----------
    Ins:
    (1) influent of pressurized clear water from clear water tank

    Outs:
    (1) recycle to flush toilet


    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.sanunits.storagetank module

    '''
    '''
    Similar to the :class:`biosteam.units.MixTank`, but can calculate material usage.

    See Also
    --------
    :class:`biosteam.units.StorageTank`
    '''
    _N_ins = 1
    _N_outs = 2
    _units = {'Diameter': 'ft',
              'Length': 'ft',
              'Wall thickness': 'in',
              'Weight': 'lb'}
    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Horizontal vessel diameter': (3, 21),
               'Vertical vessel length': (12, 40)}
    _vessel_material = 'Stainless steel'
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 vessel_type=None, tau=None, V_wf=None,
                 vessel_material=None, kW_per_m3=0.,ppl=100, baseline_ppl = 100,
                 init_with='WasteStream', F_BM_default=None,
                 include_construction=True, length_to_diameter=2):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo,
                      init_with=init_with, F_BM_default=F_BM_default,
                      include_construction=include_construction,
                      vessel_type=vessel_type, tau=tau, V_wf=V_wf,
                      vessel_material=vessel_material, kW_per_m3=kW_per_m3,)
        self.length_to_diameter = length_to_diameter
        
        
    # def _init_state(self):
    #     pass

    # def _update_state(self):
    #     pass

    # def _update_dstate(self):
    #     pass
        
    def _init_lca(self):
        item_name = self.vessel_material.replace(' ', '_')
        self.construction = [
            Construction(item_name.lower(), linked_unit=self, item=item_name, quantity_unit='kg'),
            ]
        
    
    def _design(self):
        BSTStorageTank._design(self)
        D = self.design_results
        
        _lb_to_kg = auom('lb').conversion_factor('kg')
        _m_to_ft = auom('m').conversion_factor('ft')
        _Pa_to_psi = auom('Pa').conversion_factor('psi')
        
        Diameter = (4*D['Total volume']/pi/self.length_to_diameter)**(1/3)
        Diameter *= _m_to_ft # convert from m to ft
        L = Diameter * self.length_to_diameter # ft
        D.update(self._horizontal_vessel_design(self.ins[0].P*_Pa_to_psi, Diameter, L))
        D['Material'] = self.vessel_material
        if self.include_construction: self.construction[0].quantity = D['Weight']*_lb_to_kg


    def _horizontal_vessel_design(self, pressure, diameter, length) -> dict:
        pressure = pressure
        diameter = diameter
        length = length
        # Calculate vessel weight and wall thickness
        if self.vessel_material == 'Carbon steel':
            rho_M = material_densities_lb_per_ft3[self.vessel_material]
        else:
            rho_M = material_densities_lb_per_ft3['Stainless steel 304']
        if pressure < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=DesignWarning)
        VW, VWT = flash_vessel_design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M)
        bounds_warning(self, 'Horizontal vessel weight', VW, 'lb',
                       self._bounds['Horizontal vessel weight'], 'cost')
        bounds_warning(self, 'Horizontal vessel diameter', diameter, 'ft',
                       self._bounds['Horizontal vessel diameter'], 'cost')
        Design = {}
        Design['Vessel type'] = 'Horizontal'
        Design['Length'] = length # ft
        Design['Diameter'] = diameter # ft
        Design['Weight'] = VW # lb
        Design['Wall thickness'] = VWT # in
        return Design

    @property
    def vessel_material(self):
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, i):
        exist_material = getattr(self, '_vessel_material', None)
        StorageTank.vessel_material.fset(self, i)
        if i and exist_material == i: return # type doesn't change, no need to reload construction items
        self._init_lca()
        
        
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # to be defined in .tsv file TODO: scale to mass not volume
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['Pressure water tank'] = self.pressure_water_tank_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()

        power_demand = self.power_demand_PT
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        PW_replacement_cost = (
            self.pressure_water_tank_cost / self.pressure_water_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime) * scale
        PW_replacement_cost = PW_replacement_cost / (365 * 24) # convert to USD/hr
        return PW_replacement_cost

# %%
# blower_path = ospath.join(EL_su_data_path, '_EL_blower.tsv')

# @price_ratio()
# class EL_blower(SanUnit):
#     '''
#     Introduction
#     ------------
#     To areate air for aerobic tank and membrane tank

#     Parameters
#     ----------
#     Ins:
#     (1) air

#     Outs:
#     (1) air


#     Attributes
#     ----------

    
#     References
#     ----------
#      refer to the qsdsan.equipments.Blower module

#     '''
#     _N_ins = 1  # number of ins
#     _N_outs = 1  # number of outs
#     _ins_size_is_fixed = True
#     _outs_size_is_fixed = True
#     exponent_scale = 0.1

#     def __init__(self, ID = '', ins = None, outs = (), init_with = 'WasteStream',
#                  # F_BM={
#                  #     'Blowers': 2.22,
#                  #     'Blower piping': 1,
#                  #     'Blower building': 1.11,
#                  #     },
#                  F_BM = 2.22,
#                  lifetime=15, lifetime_unit='yr',
#                  # units={
#                  #     'Total gas flow': 'CFM',
#                  #     'Blower capacity': 'CFM',
#                  #     'Number of blowers': '',
#                  #     'Total blower power': 'kW',
#                  #     },
#                  # N_reactor=2, # the number of the reactors where the gas sparging modules will be installed
#                  # gas_demand_per_reactor=1, # gas demand per reactor
#                  # TDH=6, # total dynamic head for rhe blower, in psi
#                  # eff_blower=0.7, # efficiency of the blower in fraction
#                  # eff_motor=0.7, # efficiency of the motor in fraction
#                  # AFF=3.33, # air flow fraction
#                  # building_unit_cost=9, # unit cost of the building, in USD/ft2
#                  thermo = None, ppl = None, baseline_ppl = None, **kwargs):
#         # super().__init__(ID=ID, lifetime = lifetime, lifetime_unit = lifetime_unit, F_BM=F_BM,
#         #                 units=units, N_reactor=N_reactor, gas_demand_per_reactor=gas_demand_per_reactor,
#         #                 TDH=TDH, eff_blower=eff_blower, eff_motor=eff_motor, AFF=AFF, building_unit_cost=building_unit_cost,)
#         SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=F_BM)

#         self.ppl = ppl
#         self.baseline_ppl = baseline_ppl
#         self.lifetime = lifetime
#         self.lifetime_unit = lifetime_unit

#         data = load_data(path = blower_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             setattr(self, para, value)
#         del data

#         for attr, value in kwargs.items():
#             setattr(self, attr, value)

#     def _init_lca(self):
#         self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]
    
#     def _design(self):
#         design = self.design_results
#         constr = self.construction
#         design['StainlessSteel'] = constr[0].quantity = self.blower_steel_weight  # to be defined in .tsv file
#         self.add_construction(add_cost=False)
    

#     def _cost(self):
#         C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
#         C['Blower'] = self.blower_cost
#         C['pipeline'] = self.pipeline_connectors
#         C['fittings'] = self.weld_female_adapter_fittings

#         ratio = self.price_ratio
#         for equipment, cost in C.items():
#             C[equipment] = cost * ratio

#         self.add_OPEX = self._calc_replacement_cost()
        
#         power_demand = self.power_demand_blower
#         self.power_utility(power_demand)
    
#     def _calc_replacement_cost(self):
#         scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
#         Blower_replacement_cost = (
#             self.blower_cost * self.blower_lifetime +
#             self.pipeline_connectors / self.pipeline_connectors_lifetime +
#             self.weld_female_adapter_fittings / self.fittings_lifetime) * scale
#         Blower_replacement_cost = Blower_replacement_cost / (365 * 24) # convert to USD/hr
#         return Blower_replacement_cost

# %%
# housing_path = ospath.join(EL_su_data_path, '_EL_housing.tsv')

# #@price_ratio()
# class EL_Housing(SanUnit):
#     '''
#      non_reactive unit for the Enviroloo Clear system
#     '''
#     _N_ins = 1  # number of ins
#     _N_outs = 1  # number of outs
#     _ins_size_is_fixed = True
#     _outs_size_is_fixed = True
#     ppl_per_MURT = 30  # number of people per MURT

#     def __init__(self, ID = '', ins = None, outs = (), thermo = None, init_with = None,
#                  price_ratio=0.9,
#                  ppl = 1000, baseline_ppl = 30, F_BM_default= 1, **kwargs):
#         init_with = init_with or {}
#         super().__init__(ID=ID, ins=ins, outs=outs, thermo = thermo, 
#                          init_with = init_with, F_BM_default=F_BM_default)
        
#         self.ppl = ppl
#         self.baseline_ppl = baseline_ppl
#         self.price_ratio = price_ratio

#         data = load_data(path = housing_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             setattr(self, para, value)
#         del data

#         for attr, value in kwargs.items():
#             setattr(self, attr, value)

#     def _init_lca(self): # replace the actual materials used in the EL
#         self.construction = [
#             Construction(item = 'StainlessSteel', linked_unit= self, quantity_unit= 'kg'),
#             Construction(item = 'Plastic', linked_unit= self, quantity_unit= 'kg'),]

#     def _design(self): # replace the actual materials used in the EL
#         design = self.design_results
#         constr = self.construction
#         design['StainlessSteel'] = constr[0].quantity = (self.steel_weight + self.steel_framework_weight + self.steel_fittings_weight) * (self.ppl / self.baseline_ppl)  # assume linear scaling
#         design['Plastic'] = constr[1].quantity = (self.LLDPE_weight) * (self.ppl / self.baseline_ppl)   # assume linear scaling
#         self.add_construction(add_cost= False)
    
#     def _cost(self):
#         C = self.baseline_purchase_costs
#         C['Housing'] = (self.frame + self.extrusion + 
#                         self.angle_frame + self.angle +
#                         self.door_sheet + self.plate +
#                         self.powder_coating) * (1 + 0.1 * (self.N_EL -1))
        
#         ratio = self.price_ratio
#         for equipment, cost in C.items():
#             C[equipment] = cost * ratio
    
#     @property
#     def N_EL(self): # determine the number of EL system needed
#         return ceil(self.ppl / self.baseline_ppl)
    
#     @property
#     def N_toilets(self): # determine the number of toilets needed
#         return ceil(self.ppl / self.ppl_per_MURT)

# %%
system_path = os.path.join(EL_su_data_path, '_EL_system.tsv')

@price_ratio()
class EL_System(SanUnit, isabstract=True):
    '''
    Relate to connection components in the EL system
    '''
    _N_ins = 1
    _N_outs = 0
    exponent_scale = 0.1

    def __init__(self, ID='', ins=(), outs=None, thermo = None, 
                 init_with = 'WasteStream',
                 # init_with=None,
                 if_gridtied = True, ppl = None, baseline_ppl = None, F_BM_default = 1, **kwargs):
        SanUnit.__init__(self, ID=ID, ins=ins, outs=outs, thermo = thermo, init_with = init_with,  F_BM_default = F_BM_default)
        
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self.if_gridtied = if_gridtied

        data = load_data(path = system_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [
            Construction(item = 'PVC_generic', linked_unit= self, quantity_unit= 'kg'),
            Construction(item = 'HDPE', linked_unit= self, quantity_unit= 'kg'),
            ]

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
