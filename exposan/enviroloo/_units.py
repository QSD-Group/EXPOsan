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
from qsdsan import SanUnit, Construction
from qsdsan.sanunits import IdealClarifier, Mixer
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
from qsdsan.sanunits._membrane_bioreactor import CompletelyMixedMBR
from qsdsan import Processes, CompiledProcesses
# %% This callable file will be reposited to qsdsan.SanUnit subbranch with the name of _enviroloo
__all__ = (
    'EL_Excretion', # excretion
    # 'EL_Toilet', # toilet
    'EL_MURT', # toilet
    'EL_CT', # Collection tank
    'EL_PC', # Primary clarifier
    'EL_Anoxic', # Anoxic tank
    'EL_Aerobic', # Aerobic tank
    'EL_CMMBR', # Membrane filter
    'EL_CWT', # Clear water tank
    'EL_PT', # Pressure tank
    # 'EL_blower', # blower
    'EL_System', # System-level summary
    # 'EL_Housing', # Housing of EL_System, such as equipment's armor
    )
el_path = os.path.dirname(__file__)
module = os.path.split(el_path)[-1]
data_path, results_path = _init_modules(module, include_data_path = True)
EL_su_data_path = os.path.join(data_path, 'units_data')

# %%
# excretion_path = ospath.join(EL_su_data_path, '_EL_excretion.tsv')

# class EL_Excretion(SanUnit):
#     '''
#     Estimation of N, P, K, and COD in urine and feces based on dietary intake
#     for one person based on `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

#     Parameters
#     ----------
#     waste_ratio : float
#         A ratio in [0, 1] to indicate the amount of intake calories and nutrients
#         (N, P, K) that is wasted.

#     Examples
#     --------
#     `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_

#     References
#     ----------
#     [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
#     Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
#     Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
#     https://doi.org/10.1021/acs.est.0c03296
#     '''

#     _N_ins = 0
#     _N_outs = 2

#     def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
#                  waste_ratio=0, **kwargs):
#         SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
#         self.waste_ratio = waste_ratio

#         data = load_data(path=excretion_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             # value = float(data.loc[para]['low'])
#             # value = float(data.loc[para]['high'])
#             setattr(self, '_'+para, value)
#         del data

#         for attr, value in kwargs.items():
#             setattr(self, attr, value)

#     def _run(self):
#         ur, fec = self.outs
#         ur.empty()
#         fec.empty()

#         not_wasted = 1 - self.waste_ratio
#         factor = 24 * 1e3 # from g per person per day to kg per hour

#         ur_N = (self.p_veg+self.p_anim)/factor*self.N_prot \
#            * self.N_exc*self.N_ur*not_wasted
#         ur.imass['NH3'] = ur_N * self.N_ur_NH3
#         ur.imass['NonNH3'] = ur_N - ur.imass['NH3']

#         ur.imass['P'] = (self.p_veg*self.P_prot_v+self.p_anim*self.P_prot_a)/factor \
#             * self.P_exc*self.P_ur*not_wasted

#         e_cal = self.e_cal / 24 * not_wasted
#         ur.imass['K'] = e_cal/1e3 * self.K_cal/1e3 * self.K_exc*self.K_ur
#         ur.imass['Mg'] = self.Mg_ur / factor
#         ur.imass['Ca'] = self.Ca_ur / factor

#         ur_exc = self.ur_exc / factor
#         ur.imass['H2O'] = self.ur_moi * ur_exc
#         ur.imass['OtherSS'] = ur_exc - ur.F_mass

#         fec_exc = self.fec_exc / factor
#         fec_N = (1-self.N_ur)/self.N_ur * ur_N
#         fec.imass['NH3'] = fec_N * self.N_fec_NH3
#         fec.imass['NonNH3'] = fec_N - fec.imass['NH3']
#         fec.imass['P'] = (1-self.P_ur)/self.P_ur * ur.imass['P']
#         fec.imass['K'] = (1-self.K_ur)/self.K_ur * ur.imass['K']
#         fec.imass['Mg'] = self.Mg_fec / factor
#         fec.imass['Ca'] = self.Ca_fec / factor
#         fec.imass['H2O'] = self.fec_moi * fec_exc
#         fec.imass['OtherSS'] = fec_exc - fec.F_mass

#         # 14 kJ/g COD, the average lower heating value of excreta
#         tot_COD = e_cal*self.e_exc*4.184/14/1e3 # in kg COD/hr
#         ur._COD = tot_COD*(1-self.e_fec) / (ur.F_vol/1e3) # in mg/L
#         fec._COD = tot_COD*self.e_fec / (fec.F_vol/1e3) # in mg/L

#     @property
#     def e_cal(self):
#         '''[float] Caloric intake, [kcal/cap/d].'''
#         return self._e_cal
#     @e_cal.setter
#     def e_cal(self, i):
#         self._e_cal = i

#     @property
#     def p_veg(self):
#         '''[float] Vegetal protein intake, [g/cap/d].'''
#         return self._p_veg
#     @p_veg.setter
#     def p_veg(self, i):
#         self._p_veg = i

#     @property
#     def p_anim(self):
#         '''[float] Animal protein intake, [g/cap/d].'''
#         return self._p_anim
#     @p_anim.setter
#     def p_anim(self, i):
#         self._p_anim = i

#     @property
#     def N_prot(self):
#         '''[float] Nitrogen content in protein, [wt%].'''
#         return self._N_prot
#     @N_prot.setter
#     def N_prot(self, i):
#         self._N_prot = i

#     @property
#     def P_prot_v(self):
#         '''[float] Phosphorus content in vegetal protein, [wt%].'''
#         return self._P_prot_v
#     @P_prot_v.setter
#     def P_prot_v(self, i):
#         self._P_prot_v = i

#     @property
#     def P_prot_a(self):
#         '''[float] Phosphorus content in animal protein, [wt%].'''
#         return self._P_prot_a
#     @P_prot_a.setter
#     def P_prot_a(self, i):
#         self._P_prot_a = i

#     @property
#     def K_cal(self):
#         '''[float] Potassium intake relative to caloric intake, [g K/1000 kcal].'''
#         return self._K_cal
#     @K_cal.setter
#     def K_cal(self, i):
#         self._K_cal = i

#     @property
#     def N_exc(self):
#         '''[float] Nitrogen excretion factor, [% of intake].'''
#         return self._N_exc
#     @N_exc.setter
#     def N_exc(self, i):
#         self._N_exc = i

#     @property
#     def P_exc(self):
#         '''[float] Phosphorus excretion factor, [% of intake].'''
#         return self._P_exc
#     @P_exc.setter
#     def P_exc(self, i):
#         self._P_exc = i

#     @property
#     def K_exc(self):
#         '''[float] Potassium excretion factor, [% of intake].'''
#         return self._K_exc
#     @K_exc.setter
#     def K_exc(self, i):
#         self._K_exc = i

#     @property
#     def e_exc(self):
#         '''[float] Energy excretion factor, [% of intake].'''
#         return self._e_exc
#     @e_exc.setter
#     def e_exc(self, i):
#         self._e_exc = i

#     @property
#     def N_ur(self):
#         '''[float] Nitrogen recovered in urine, [wt%].'''
#         return self._N_ur
#     @N_ur.setter
#     def N_ur(self, i):
#         self._N_ur = i

#     @property
#     def P_ur(self):
#         '''[float] Phosphorus recovered in urine, [wt%].'''
#         return self._P_ur
#     @P_ur.setter
#     def P_ur(self, i):
#         self._P_ur = i

#     @property
#     def K_ur(self):
#         '''[float] Potassium recovered in urine, [wt%].'''
#         return self._K_ur
#     @K_ur.setter
#     def K_ur(self, i):
#         self._K_ur = i

#     @property
#     def e_fec(self):
#         '''[float] Percent of excreted energy in feces, [%].'''
#         return self._e_fec
#     @e_fec.setter
#     def e_fec(self, i):
#         self._e_fec = i

#     @property
#     def N_ur_NH3(self):
#         '''[float] Reduced inorganic nitrogen in urine, modeled as NH3, [% of total urine N].'''
#         return self._N_ur_NH3
#     @N_ur_NH3.setter
#     def N_ur_NH3(self, i):
#         self._N_ur_NH3 = i

#     @property
#     def N_fec_NH3(self):
#         '''[float] Reduced inorganic nitrogen in feces, modeled as NH3, [% of total feces N].'''
#         return self._N_fec_NH3
#     @N_fec_NH3.setter
#     def N_fec_NH3(self, i):
#         self._N_fec_NH3 = i

#     @property
#     def ur_exc(self):
#         '''[float] Urine generated per day, [g/cap/d].'''
#         return self._ur_exc
#     @ur_exc.setter
#     def ur_exc(self, i):
#         self._ur_exc = i

#     @property
#     def fec_exc(self):
#         '''[float] Feces generated per day, [g/cap/d].'''
#         return self._fec_exc
#     @fec_exc.setter
#     def fec_exc(self, i):
#         self._fec_exc = i

#     @property
#     def ur_moi(self):
#         '''[float] Moisture (water) content of urine, [wt%].'''
#         return self._ur_moi
#     @ur_moi.setter
#     def ur_moi(self, i):
#         self._ur_moi = i

#     @property
#     def fec_moi(self):
#         '''[float] Moisture (water) content of feces, [wt%].'''
#         return self._fec_moi
#     @fec_moi.setter
#     def fec_moi(self, i):
#         self._fec_moi = i

#     @property
#     def Mg_ur(self):
#         '''[float] Magnesium excreted in urine, [g Mg/cap/d].'''
#         return self._Mg_ur
#     @Mg_ur.setter
#     def Mg_ur(self, i):
#         self._Mg_ur = i

#     @property
#     def Mg_fec(self):
#         '''[float] Magnesium excreted in feces, [g Mg/cap/d].'''
#         return self._Mg_fec
#     @Mg_fec.setter
#     def Mg_fec(self, i):
#         self._Mg_fec = i

#     @property
#     def Ca_ur(self):
#         '''[float] Calcium excreted in urine, [g Ca/cap/d].'''
#         return self._Ca_ur
#     @Ca_ur.setter
#     def Ca_ur(self, i):
#         self._Ca_ur = i

#     @property
#     def Ca_fec(self):
#         '''[float] Calcium excreted in feces, [g Ca/cap/d].'''
#         return self._Ca_fec
#     @Ca_fec.setter
#     def Ca_fec(self, i):
#         self._Ca_fec = i

#     @property
#     def waste_ratio(self):
#         '''
#         [float] The amount of intake calories and nutrients
#         (N, P, K) that is wasted.

#         .. note::
#             Not considered for Mg and Ca.
#         '''
#         return self._waste_ratio
#     @waste_ratio.setter
#     def waste_ratio(self, i):
#         self._waste_ratio = i

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

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 waste_ratio=0, **kwargs):
        super().__init__(ID, ins, outs, thermo, init_with, waste_ratio, **kwargs)
        isdyn = kwargs.pop('isdynamic', False)
        if isdyn: self._init_dynamic()
    
    def _run(self):
        ur, fec = self.outs
        ur.empty()
        fec.empty()
        cmps = ur.components
        sf_iN = cmps.S_F.i_N
        xs_iN = cmps.X_S.i_N
        xb_iN = cmps.X_H.i_N
        sxi_iN = cmps.S_I.i_N
        i_mass = cmps.i_mass
        i_P = cmps.i_P
        hco3_imass = cmps.S_IC.i_mass

        not_wasted = 1 - self.waste_ratio
        factor = 24 * 1e3 # from g/cap/d to kg/hr(/cap)
        e_cal = self.e_cal / 24 * not_wasted # kcal/cap/d --> kcal/cap/hr
        ur_exc = self.ur_exc / factor
        fec_exc = self.fec_exc / factor

        # 14 kJ/g COD, the average lower heating value of excreta
        tot_COD = e_cal*self.e_exc*4.184/14/1e3 # in kg COD/hr
        fec_COD = tot_COD*self.e_fec
        ur_COD = tot_COD - fec_COD
        
        tot_N = (self.p_veg+self.p_anim)*self.N_prot/factor \
            * self.N_exc*not_wasted
        ur_N = tot_N*self.N_ur
        fec_N = tot_N - ur_N
        
        tot_P = (self.p_veg*self.P_prot_v+self.p_anim*self.P_prot_a)/factor \
            * self.P_exc*not_wasted
        ur_P = tot_P*self.P_ur
        fec_P = tot_P - ur_P
        
        # breakpoint()
        ur.imass['S_NH4'] = ur_nh4 = ur_N * self.N_ur_NH3
        req_sf_cod = (ur_N - ur_nh4) / sf_iN
        if req_sf_cod <= ur_COD:
            ur.imass['S_F'] = sf = req_sf_cod
            ur.imass['S_A'] = ur_COD - sf  # contains no N or P
        else:
            req_si_cod = (ur_N - ur_nh4) / sxi_iN
            if req_si_cod <= ur_COD:
                ur.imass['S_F'] = sf = (sxi_iN * ur_COD - (ur_N - ur_nh4))/(sxi_iN - sf_iN)
                ur.imass['S_I'] = ur_COD - sf
            else:
                ur.imass['S_F'] = sf = ur_COD
                ur_other_n = ur_N - ur_nh4 - sf * sf_iN
                warn(f"Excess non-NH3 nitrogen cannot be accounted for by organics "
                     f"in urine: {ur_other_n} kg/hr. Added to NH3-N.")
                ur.imass['S_NH4'] += ur_other_n # debatable, has negative COD # raise warning/error
        
        ur.imass['S_PO4'] = ur_P - sum(ur.mass * i_P)
        ur.imass['S_K'] = e_cal/1e3 * self.K_cal/1e3 * self.K_exc*self.K_ur
        ur.imass['S_Mg'] = self.Mg_ur / factor
        ur.imass['S_Ca'] = self.Ca_ur / factor

        ur.imass['H2O'] = self.ur_moi * ur_exc
        ur_others = ur_exc - sum(ur.mass * i_mass)
        ur.imass['S_IC'] = ur_others * 0.34 / hco3_imass
        ur.imass['S_Na'] = ur_others * 0.35
        ur.imass['S_Cl'] = ur_others * 0.31

        fec.imass['S_NH4'] = fec_nh4 = fec_N * self.N_fec_NH3
        req_xs_cod = (fec_N - fec_nh4) / xs_iN
        if req_xs_cod <= fec_COD:
            fec.imass['X_S'] = xs = req_xs_cod
            fec.imass['S_A'] = fec_COD - xs
        else:
            req_xi_cod = (fec_N - fec_nh4) / sxi_iN
            if req_xi_cod <= fec_COD:
                fec.imass['X_S'] = xs = (sxi_iN * fec_COD - (fec_N - fec_nh4))/(sxi_iN - xs_iN)
                fec.imass['X_I'] = fec_COD - xs
            else:
                req_xb_cod = (fec_N - fec_nh4) / xb_iN
                if req_xb_cod <= fec_COD:
                    fec.imass['X_S'] = xs = (xb_iN * fec_COD - (fec_N - fec_nh4))/(xb_iN - xs_iN)
                    fec.imass['X_H'] = fec_COD - xs
                else:
                    fec.imass['X_S'] = xs = fec_COD
                    fec_other_n = fec_N - fec_nh4 - xs * xs_iN
                    warn(f"Excess non-NH3 nitrogen cannot be accounted for by organics "
                         f"in feces: {fec_other_n} kg/hr. Added to NH3-N.")
                    fec.imass['S_NH4'] += fec_other_n # debatable, has negative COD
        
        fec.imass['S_PO4'] = fec_P - sum(fec.mass * i_P)
        fec.imass['S_K'] = (1-self.K_ur)/self.K_ur * ur.imass['S_K']
        fec.imass['S_Mg'] = self.Mg_fec / factor
        fec.imass['S_Ca'] = self.Ca_fec / factor
        fec.imass['H2O'] = self.fec_moi * fec_exc
        
        fec_others = fec_exc - sum(fec.mass * i_mass)
        fec.imass['S_IC'] = fec_others * 0.34 / hco3_imass
        fec.imass['S_Na'] = fec_others * 0.35
        fec.imass['S_Cl'] = fec_others * 0.31
        
    
    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE
    
    def _compile_AE(self):
        def yt(t, QC_ins, dQC_ins):
            pass
        self._AE = yt
        
    def _init_state(self):
        ur, fec = self.outs
        self._state = np.append(ur.mass, fec.mass)
        for ws in self.outs:
            ws.state = np.append(ws.conc, ws.F_vol * 24)
            ws.dstate = np.zeros_like(ws.state)

    def _update_state(self):
        pass

    def _update_dstate(self):
        pass

# %%
# toilet_path = ospath.join(EL_su_data_path, '_EL_toilet.tsv')

# class EL_Toilet(SanUnit, Decay, isabstract=True):
#     '''
#     Abstract class containing common parameters and design algorithms for toilets
#     based on `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

#     Parameters
#     ----------
#     degraded_components : tuple
#         IDs of components that will degrade (simulated by first-order decay).
#     N_user : int, float
#         Number of people per toilet.
#         Note that this number can be a float when calculated from `N_tot_user` and `N_toilet`.
#     N_toilet : int
#         Number of parallel toilets.
#         In calculation, `N_toilet` will be calculated as `ceil(N_tot_user/N_user)`.
#     N_tot_user : int
#         Total number of users.

#         .. note::

#             If `N_tot_user` is provided (i.e., not "None"),
#             then updating `N_user` will recalculate `N_toilet`, and vice versa.

#     if_toilet_paper : bool
#         If toilet paper is used.
#     if_flushing : bool
#         If water is used for flushing.
#     if_cleansing : bool
#         If water is used for cleansing.
#     if_desiccant : bool
#         If desiccant is used for moisture and odor control.
#     if_air_emission : bool
#         If emission to air occurs
#         (i.e., if the pit is completely sealed off from the atmosphere).
#     if_ideal_emptying : bool
#         If the toilet appropriately emptied to avoid contamination to the
#         environmental.
#     CAPEX : float
#         Capital cost of a single toilet.
#     OPEX_over_CAPEX : float
#         Fraction of annual operating cost over total capital cost.
#     price_ratio : float
#         Calculated capital cost will be multiplied by this number
#         to consider the effect in cost difference from different locations.

#     References
#     ----------
#     [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
#     Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
#     Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
#     https://doi.org/10.1021/acs.est.0c03296.

#     See Also
#     --------
#     :ref:`qsdsan.processes.Decay <processes_Decay>`

#     '''
#     _N_ins = 6
#     _outs_size_is_fixed = False
#     density_dct = {
#         'Sand': 1442,
#         'Gravel': 1600,
#         'Brick': 1750,
#         'Plastic': 0.63,
#         'Steel': 7900,
#         'StainlessSteelSheet': 2.64
#         }

#     def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
#                  degraded_components=('OtherSS',), N_user=1, N_toilet=1, N_tot_user=None,
#                  if_toilet_paper=True, if_flushing=True, if_cleansing=False,
#                  if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
#                  CAPEX=None, OPEX_over_CAPEX=None, price_ratio=1.):

#         SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=1)
#         self.degraded_components = tuple(degraded_components)
#         self._N_user = self._N_toilet = self._N_tot_user = None
#         self.N_user = N_user
#         self.N_toilet = N_toilet
#         self.N_tot_user = N_tot_user
#         self.if_toilet_paper = if_toilet_paper
#         self.if_flushing = if_flushing
#         self.if_cleansing = if_cleansing
#         self.if_desiccant = if_desiccant
#         self.if_air_emission = if_air_emission
#         self.if_ideal_emptying = if_ideal_emptying
#         self.CAPEX = CAPEX
#         self.OPEX_over_CAPEX = OPEX_over_CAPEX
#         self.price_ratio = price_ratio

#         data = load_data(path=toilet_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             if para in ('desiccant_V', 'desiccant_rho'):
#                 setattr(self, para, value)
#             else:
#                 setattr(self, '_'+para, value)
#         del data

#         # self._empty_ratio = 0.59 default
#         self._empty_ratio = 0


#     def _run(self):
#         ur, fec, tp, fw, cw, des = self.ins
#         tp.imass['Tissue'] = int(self.if_toilet_paper)*self.toilet_paper
#         fw.imass['H2O'] = int(self.if_flushing)*self.flushing_water
#         cw.imass['H2O'] = int(self.if_cleansing)*self.cleansing_water
#         des.imass['WoodAsh'] = int(self.if_desiccant)*self.desiccant

#     def _scale_up_outs(self):
#         '''
#         Scale up the effluent based on the number of user per toilet and
#         toilet number.
#         '''
#         N_tot_user = self.N_tot_user or self.N_toilet*self.N_user
#         for i in self.outs:
#             if not i.F_mass == 0:
#                 i.F_mass *= N_tot_user


#     def _cost(self):
#         self.baseline_purchase_costs['Total toilets'] = self.CAPEX * self.N_toilet * self.price_ratio
#         add_OPEX = self.baseline_purchase_costs['Total toilets']*self.OPEX_over_CAPEX/365/24
#         self._add_OPEX = {'Additional OPEX': add_OPEX}


#     @staticmethod
#     def get_emptying_emission(waste, CH4, N2O, empty_ratio, CH4_factor, N2O_factor):
#         '''
#         Calculate emissions due to non-ideal emptying based on
#         `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_,

#         Parameters
#         ----------
#         stream : WasteStream
#             Excreta stream that is not appropriately emptied (before emptying).
#         CH4 : WasteStream
#             Fugitive CH4 gas (before emptying).
#         N2O : WasteStream
#             Fugitive N2O gas (before emptying).
#         empty_ratio : float
#             Fraction of excreta that is appropriately emptied..
#         CH4_factor : float
#             Factor to convert COD removal to CH4 emission.
#         N2O_factor : float
#             Factor to convert COD removal to N2O emission.

#         Returns
#         -------
#         stream : WasteStream
#             Excreta stream that is not appropriately emptied (after emptying).
#         CH4 : WasteStream
#             Fugitive CH4 gas (after emptying).
#         N2O : WasteStream
#             Fugitive N2O gas (after emptying).
#         '''
#         COD_rmvd = waste.COD*(1-empty_ratio)/1e3*waste.F_vol
#         CH4.imass['CH4'] += COD_rmvd * CH4_factor
#         N2O.imass['N2O'] += COD_rmvd * N2O_factor
#         waste.mass *= empty_ratio

#     @property
#     def N_user(self):
#         '''[int, float] Number of people per toilet.'''
#         return self._N_user or self.N_tot_user/self.N_toilet
#     @N_user.setter
#     def N_user(self, i):
#         if i is not None:
#             N_user = self._N_user = int(i)
#             old_toilet = self._N_toilet
#             if old_toilet and self.N_tot_user:
#                 new_toilet = ceil(self.N_tot_user/N_user)
#                 warn(f'With the provided `N_user`, the previous `N_toilet` of {old_toilet} '
#                      f'is recalculated from `N_tot_user` and `N_user` as {new_toilet}.')
#                 self._N_toilet = None
#         else:
#             self._N_user = i

#     @property
#     def N_toilet(self):
#         '''[int] Number of parallel toilets.'''
#         return self._N_toilet or ceil(self.N_tot_user/self.N_user)
#     @N_toilet.setter
#     def N_toilet(self, i):
#         if i is not None:
#             N_toilet = self._N_toilet = ceil(i)
#             old_user = self._N_user
#             if old_user and self.N_tot_user:
#                 new_user = self.N_tot_user/N_toilet
#                 warn(f'With the provided `N_toilet`, the previous `N_user` of {old_user} '
#                      f'is recalculated from `N_tot_user` and `N_toilet` as {new_user}.')
#                 self._N_user = None
#         else:
#             self._N_toilet = i

#     @property
#     def N_tot_user(self):
#         '''[int] Number of total users.'''
#         return self._N_tot_user
#     @N_tot_user.setter
#     def N_tot_user(self, i):
#         if i is not None:
#             self._N_tot_user = int(i)
#         else:
#             self._N_tot_user = None

#     @property
#     def toilet_paper(self):
#         '''
#         [float] Amount of toilet paper used
#         (if ``if_toilet_paper`` is True), [kg/cap/hr].
#         '''
#         return self._toilet_paper
#     @toilet_paper.setter
#     def toilet_paper(self, i):
#         self._toilet_paper = i

#     @property
#     def flushing_water(self):
#         '''
#         [float] Amount of water used for flushing
#         (if ``if_flushing_water`` is True), [kg/cap/hr].
#         '''
#         return self._flushing_water
#     @flushing_water.setter
#     def flushing_water(self, i):
#         self._flushing_water = i

#     @property
#     def cleansing_water(self):
#         '''
#         [float] Amount of water used for cleansing
#         (if ``if_cleansing_water`` is True), [kg/cap/hr].
#         '''
#         return self._cleansing_water
#     @cleansing_water.setter
#     def cleansing_water(self, i):
#         self._cleansing_water = i

#     @property
#     def desiccant(self):
#         '''
#         [float] Amount of desiccant used (if ``if_desiccant`` is True), [kg/cap/hr].

#         .. note::

#             Value set by ``desiccant_V`` and ``desiccant_rho``.

#         '''
#         return self.desiccant_V*self.desiccant_rho

#     @property
#     def N_volatilization(self):
#         '''
#         [float] Fraction of input N that volatilizes to the air
#         (if ``if_air_emission`` is True).
#         '''
#         return self._N_volatilization
#     @N_volatilization.setter
#     def N_volatilization(self, i):
#         self._N_volatilization = i

#     @property
#     def empty_ratio(self):
#         '''
#         [float] Fraction of excreta that is appropriately emptied.

#         .. note::

#             Will be 1 (i.e., 100%) if ``if_ideal_emptying`` is True.

#         '''
#         if self.if_ideal_emptying:
#             return 1.
#         return self._empty_ratio
#     @empty_ratio.setter
#     def empty_ratio(self, i):
#         if self.if_ideal_emptying:
#             warn(f'`if_ideal_emptying` is True, the set value {i} is ignored.')
#         self._empty_ratio = i

#     @property
#     def MCF_aq(self):
#         '''[float] Methane correction factor for COD lost due to inappropriate emptying.'''
#         return self._MCF_aq
#     @MCF_aq.setter
#     def MCF_aq(self, i):
#         self._MCF_aq = i

#     @property
#     def N2O_EF_aq(self):
#         '''[float] Fraction of N emitted as N2O due to inappropriate emptying.'''
#         return self._N2O_EF_aq
#     @N2O_EF_aq.setter
#     def N2O_EF_aq(self, i):
#         self._N2O_EF_aq = i

#     @property
#     def if_N2O_emission(self):
#         '''[bool] Whether to consider N degradation and fugitive N2O emission.'''
#         return self.if_air_emission
#     @if_N2O_emission.setter
#     def if_N2O_emission(self, i):
#         raise ValueError('Setting `if_N2O_emission` for `PitLatrine` is not supported, '
#                          'please set `if_air_emission` instead.')

# %%
# murt_path = ospath.join(EL_su_data_path, '_EL_murt.tsv')

# @price_ratio()
# class EL_MURT(EL_Toilet):
#     '''
#     Multi-unit reinvented toilet.

#     The following components should be included in system thermo object for simulation:
#     Tissue, WoodAsh, H2O, NH3, NonNH3, P, K, Mg, CH4, N2O.

#     The following impact items should be pre-constructed for life cycle assessment:
#     Ceramic, Fan.

#     Parameters
#     ----------
#     ins : Iterable(stream)
#         waste_in: mixed excreta.
#     Outs : Iterable(stream)
#         waste_out: degraded mixed excreta.
#         CH4: fugitive CH4.
#         N2O: fugitive N2O.
#     N_squatting_pan_per_toilet : int
#         The number of squatting pan per toilet.
#     N_urinal_per_toilet : int
#         The number of urinals per toilet.
#     if_include_front_end : bool
#         If False, will not consider the capital and operating costs of this unit.

#     See Also
#     --------
#     :ref:`qsdsan.sanunits.Toilet <sanunits_toilet>`
#     '''
#     _N_outs = 3
#     _units = {
#         'Collection period': 'd',
#         }

#     def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
#                  degraded_components=('OtherSS',), N_user=1, N_tot_user=1,
#                  N_toilet=None, if_toilet_paper=True, if_flushing=True, if_cleansing=False,
#                  if_desiccant=True, if_air_emission=True, if_ideal_emptying=True,
#                  CAPEX=0, OPEX_over_CAPEX=0, lifetime=10,
#                  N_squatting_pan_per_toilet=1, N_urinal_per_toilet=1,
#                  if_include_front_end=True, **kwargs):

#         EL_Toilet.__init__(
#             self, ID, ins, outs, thermo=thermo, init_with=init_with,
#             degraded_components=degraded_components,
#             N_user=N_user, N_tot_user=N_tot_user, N_toilet=N_toilet,
#             if_toilet_paper=if_toilet_paper, if_flushing=if_flushing,
#             if_cleansing=if_cleansing, if_desiccant=if_desiccant,
#             if_air_emission=if_air_emission, if_ideal_emptying=if_ideal_emptying,
#             CAPEX=CAPEX, OPEX_over_CAPEX=OPEX_over_CAPEX
#             )
#         self.N_squatting_pan_per_toilet = N_squatting_pan_per_toilet
#         self.N_urinal_per_toilet = N_urinal_per_toilet
#         self.if_include_front_end = if_include_front_end
#         self._mixed_in = WasteStream(f'{self.ID}_mixed_in')

#         data = load_data(path=murt_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             setattr(self, para, value)
#         del data

#         for attr, value in kwargs.items():
#             setattr(self, attr, value)

#     def _init_lca(self):
#         self.construction = [
#             Construction(item='Ceramic', linked_unit=self, quantity_unit='kg'),
#             Construction(item='Fan', linked_unit=self, quantity_unit='ea'),
#         ]

#     def _run(self):
#         EL_Toilet._run(self)
#         mixed_out, CH4, N2O = self.outs
#         CH4.phase = N2O.phase = 'g'

#         mixed_in = self._mixed_in
#         mixed_in.mix_from(self.ins)
#         #tot_COD_kg = sum(float(getattr(i, 'COD', 0)) * i.F_vol for i in self.ins) / 1e3
#         tot_COD_kg = sum(float(getattr(i, 'COD')) * i.F_vol for i in self.ins) / 1e3
        
#         # Air emission
#         if self.if_air_emission:
#             # N loss due to ammonia volatilization
#             NH3_rmd, NonNH3_rmd = \
#                 self.allocate_N_removal(mixed_in.TN/1e3*mixed_in.F_vol*self.N_volatilization,
#                                         mixed_in.imass['NH3'])
#             mixed_in.imass ['NH3'] -= NH3_rmd
#             mixed_in.imass['NonNH3'] -= NonNH3_rmd
            
#             # Energy/N loss due to degradation
#             mixed_in._COD = tot_COD_kg * 1e3 / mixed_in.F_vol # accounting for COD loss in leachate
#             Decay._first_order_run(self, waste=mixed_in, treated=mixed_out, CH4=CH4, N2O=N2O)
#         else:
#             mixed_out.copy_like(mixed_in)
#             CH4.empty()
#             N2O.empty()
            
#         # Aquatic emission when not ideally emptied
#         if not self.if_ideal_emptying:
#            self.get_emptying_emission(
#                 waste=mixed_out, CH4=CH4, N2O=N2O,
#                 empty_ratio=self.empty_ratio,
#                 CH4_factor=self.COD_max_decay*self.MCF_aq*self.max_CH4_emission,
#                 N2O_factor=self.N2O_EF_decay*44/28)
        
#         self._scale_up_outs()

#     def _design(self):
#         design = self.design_results
#         constr = self.construction
#         if self.if_include_front_end:
#             design['Number of users per toilet'] = self.N_user
#             design['Parallel toilets'] = N = self.N_toilet
#             design['Collection period'] = self.collection_period
#             design['Ceramic'] = Ceramic_quant = (
#                 self.squatting_pan_weight * self.N_squatting_pan_per_toilet+
#                 self.urinal_weight * self.N_urinal_per_toilet
#                 )
#             design['Fan'] = Fan_quant = 1  # assume fan quantity is 1
#             constr[0].quantity = Ceramic_quant * N
#             constr[1].quantity = Fan_quant * N
#             self.add_construction(add_cost=False)
#         else:
#             design.clear()
#             for i in constr: i.quantity = 0

#     def _cost(self):
#         C = self.baseline_purchase_costs
#         if self.if_include_front_end:
#             N_toilet = self.N_toilet
#             C['Ceramic Toilets'] = (
#                 self.squatting_pan_cost * self.N_squatting_pan_per_toilet +
#                 self.urinal_cost * self.N_urinal_per_toilet
#                 ) * N_toilet
#             C['Fan'] = self.fan_cost * N_toilet
#             C['Misc. parts'] = (
#                 self.led_cost +
#                 self.anticor_floor_cost +
#                 self.circuit_change_cost +
#                 self.pipe_cost
#                 ) * N_toilet

#             ratio = self.price_ratio
#             for equipment, cost in C.items():
#                 C[equipment] = cost * ratio
#         else:
#             self.baseline_purchase_costs.clear()

#         sum_purchase_costs = sum(v for v in C.values())
#         self.add_OPEX = (
#             self._calc_replacement_cost() +
#             self._calc_maintenance_labor_cost() +
#             sum_purchase_costs * self.OPEX_over_CAPEX / (365 * 24)
#             )

#     def _calc_replacement_cost(self):
#         return 0

#     def _calc_maintenance_labor_cost(self):
#         return 0

#     @property
#     def collection_period(self):
#         '''[float] Time interval between storage tank collection, [d].'''
#         return self._collection_period
#     @collection_period.setter
#     def collection_period(self, i):
#         self._collection_period = float(i)
        
#     @property
#     def tau(self):
#         '''[float] Retention time of the unit, same as `collection_period`.'''
#         return self.collection_period
#     @tau.setter
#     def tau(self, i):
#         self.collection_period = i

murt_path = os.path.join(EL_su_data_path, '_EL_murt.tsv')

@price_ratio()
class EL_MURT(SanUnit):
    '''
    Multi-unit reinvented toilet.

    The following components should be included in system thermo object for simulation:
    Tissue, WoodAsh, H2O, NH3, NonNH3, P, K, Mg, CH4, N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    Ceramic, Fan.

    Parameters
    ----------
    ins : Iterable(stream)
        waste_in: mixed excreta.
    Outs : Iterable(stream)
        waste_out: degraded mixed excreta.
        CH4: fugitive CH4.
        N2O: fugitive N2O.
    N_squatting_pan_per_toilet : int
        The number of squatting pan per toilet.
    N_urinal_per_toilet : int
        The number of urinals per toilet.
    if_include_front_end : bool
        If False, will not consider the capital and operating costs of this unit.

    See Also
    --------
    :ref:`qsdsan.sanunits.Toilet <sanunits_toilet>`
    '''
    _N_ins = 6
    _N_outs = 3
    _units = {
        'Collection period': 'd',
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), N_user=1, N_tot_user=1,
                 N_toilet=None, 
                 # if_toilet_paper=False,
                 if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=0, OPEX_over_CAPEX=0, lifetime=10,
                 N_squatting_pan_per_toilet=1, N_urinal_per_toilet=1,
                 if_include_front_end=True, **kwargs):

        SanUnit.__init__(
            self, ID, ins, outs, thermo=thermo, init_with=init_with,
            degraded_components=degraded_components,
            N_user=N_user, N_tot_user=N_tot_user, N_toilet=N_toilet,
            # if_toilet_paper=if_toilet_paper, 
            if_flushing=if_flushing,
            if_cleansing=if_cleansing, if_desiccant=if_desiccant,
            if_air_emission=if_air_emission, if_ideal_emptying=if_ideal_emptying,
            CAPEX=CAPEX, OPEX_over_CAPEX=OPEX_over_CAPEX
            )
        self.N_squatting_pan_per_toilet = N_squatting_pan_per_toilet
        self.N_urinal_per_toilet = N_urinal_per_toilet
        self.if_include_front_end = if_include_front_end
        self._mixed_in = WasteStream(f'{self.ID}_mixed_in')

        data = load_data(path=murt_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [
            Construction(item='Ceramic', linked_unit=self, quantity_unit='kg'),
            Construction(item='Fan', linked_unit=self, quantity_unit='ea'),
        ]

    def _run(self):
        # Toilet._run(self)
        mixed_out, CH4, N2O = self.outs
        CH4.phase = N2O.phase = 'g'

        mixed_in = self._mixed_in
        mixed_in.mix_from(self.ins)
        #tot_COD_kg = sum(float(getattr(i, 'COD', 0)) * i.F_vol for i in self.ins) / 1e3
        tot_COD_kg = sum(float(getattr(i, 'COD')) * i.F_vol for i in self.ins) / 1e3
        
        # # Air emission
        # if self.if_air_emission:
        #     # N loss due to ammonia volatilization
        #     NH3_rmd, NonNH3_rmd = \
        #         self.allocate_N_removal(mixed_in.TN/1e3*mixed_in.F_vol*self.N_volatilization,
        #                                 mixed_in.imass['NH3'])
        #     mixed_in.imass ['NH3'] -= NH3_rmd
        #     mixed_in.imass['NonNH3'] -= NonNH3_rmd
            
        #     # Energy/N loss due to degradation
        #     mixed_in._COD = tot_COD_kg * 1e3 / mixed_in.F_vol # accounting for COD loss in leachate
        #     Decay._first_order_run(self, waste=mixed_in, treated=mixed_out, CH4=CH4, N2O=N2O)
        # else:
        #     mixed_out.copy_like(mixed_in)
        #     CH4.empty()
        #     N2O.empty()
            
        # Aquatic emission when not ideally emptied
        if not self.if_ideal_emptying:
           self.get_emptying_emission(
                waste=mixed_out, CH4=CH4, N2O=N2O,
                empty_ratio=self.empty_ratio,
                CH4_factor=self.COD_max_decay*self.MCF_aq*self.max_CH4_emission,
                N2O_factor=self.N2O_EF_decay*44/28)
        
        # self._scale_up_outs()

    def _design(self):
        design = self.design_results
        constr = self.construction
        if self.if_include_front_end:
            design['Number of users per toilet'] = self.N_user
            design['Parallel toilets'] = N = self.N_toilet
            design['Collection period'] = self.collection_period
            design['Ceramic'] = Ceramic_quant = (
                self.squatting_pan_weight * self.N_squatting_pan_per_toilet+
                self.urinal_weight * self.N_urinal_per_toilet
                )
            design['Fan'] = Fan_quant = 1  # assume fan quantity is 1
            constr[0].quantity = Ceramic_quant * N
            constr[1].quantity = Fan_quant * N
            self.add_construction(add_cost=False)
        else:
            design.clear()
            for i in constr: i.quantity = 0

    def _cost(self):
        C = self.baseline_purchase_costs
        if self.if_include_front_end:
            N_toilet = self.N_toilet
            C['Ceramic Toilets'] = (
                self.squatting_pan_cost * self.N_squatting_pan_per_toilet +
                self.urinal_cost * self.N_urinal_per_toilet
                ) * N_toilet
            C['Fan'] = self.fan_cost * N_toilet
            C['Misc. parts'] = (
                self.led_cost +
                self.anticor_floor_cost +
                self.circuit_change_cost +
                self.pipe_cost
                ) * N_toilet

            ratio = self.price_ratio
            for equipment, cost in C.items():
                C[equipment] = cost * ratio
        else:
            self.baseline_purchase_costs.clear()

        sum_purchase_costs = sum(v for v in C.values())
        self.add_OPEX = (
            self._calc_replacement_cost() +
            self._calc_maintenance_labor_cost() +
            sum_purchase_costs * self.OPEX_over_CAPEX / (365 * 24)
            )

    def _calc_replacement_cost(self):
        return 0

    def _calc_maintenance_labor_cost(self):
        return 0

    @property
    def collection_period(self):
        '''[float] Time interval between storage tank collection, [d].'''
        return self._collection_period
    @collection_period.setter
    def collection_period(self, i):
        self._collection_period = float(i)
        
    @property
    def tau(self):
        '''[float] Retention time of the unit, same as `collection_period`.'''
        return self.collection_period
    @tau.setter
    def tau(self, i):
        self.collection_period = i

# %%
# CollectionTank_path = ospath.join(EL_su_data_path, '_EL_CT.tsv')

# @price_ratio()
# class EL_CT(StorageTank):
    
#     '''
#     Name
#     ----
#     Collection tank in the Enviroloo (EL) Clear Toilet system.
    
#     Parameters
#     ----------
#     Ins: 
#     (1) Mixed wastewater
#     (2) Primary clarifier effluent spill
#     (3) Clear water tank effluent spill 
#     (4) Primary clarifier sludge return

#     Outs:
#     (1) Treated water
#     (2) Methane (CH4)
#     (3) Nitrous oxide (N2O)

#     Attributes
#     ----------
#     length_to_diameter : float
#         Ratio of the tank length to diameter.
#     vessel_material : str
#         Material used for constructing the vessel.
#     sludge_moisture_content : float
#         Moisture content of the sludge (mass of water/total mass).
#     COD_removal : float
#         Fraction of COD removed in the collection tank.

#     References
#     ----------
#     Similar to the :class:`biosteam.units.MixTank`, but can calculate material usage.

#     See Also
#     ----------
#     class:`biosteam.units.StorageTank`
#     '''
    
#     _N_ins = 4 
#     _N_outs = 1
#     _ins_size_is_fixed = True
#     _outs_size_is_fixed = True
#     exponent_scale = 0.1
    
#     def __init__(self, ID='', ins=None, outs=(), thermo=None, ppl=None, baseline_ppl=None,
#                  vessel_type= 'Field erected', tau=24, V_wf=None, vessel_material='Stainless steel', kW_per_m3=0.1,
#                  init_with='WasteStream', F_BM_default=1,
#                  include_construction=True, **kwargs):
#         StorageTank.__init__(self, 
#                       # Basic parameters
#                       ID=ID, # The unique identifier of the tank
#                       ins=ins, # The input stream to the tank
#                       outs=outs, # The output streams from the tank
#                       thermo=thermo, # The thermodynamic property package for simulating physical or chemical processes
#                       # Other control parameters
#                       init_with=init_with, # The method to initialize the tank contents.
#                       F_BM_default=F_BM_default, # The default bare module factor for the tank cost estimation
#                       include_construction=include_construction, # A Boolean value indicating whether the tank's construction material is considered in the cost analysis or life cycle assessment.
#                       # Design parameters
#                       vessel_type=vessel_type, # The type of the tank
#                       tau=tau, # The retention time of the tank contents, the important parameters involved in mixing, reaction or separation processes.
#                       V_wf=V_wf, # The volume working fraction of the tank, the ratio of the volume of the tank contents to the total volume of the tank.
#                       vessel_material=vessel_material, # The material of the tank
#                       kW_per_m3=kW_per_m3, # The power consumption per unit volume of the tank
#                       )
#         self.ppl = ppl
#         self.baseline_ppl = baseline_ppl
#         self._mixed_in = WasteStream(f'{self.ID}_mixed_in')

#         data = load_data(path=CollectionTank_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             setattr(self, para, value)
#         del data

#         for attr, value in kwargs.items():
#             setattr(self, attr, value)

#     def _init_lca(self):
#         self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]

#     def _run(self):
#         # Input stream
#         WasteWater = self.ins[0]
#         sludge_return = self.ins[1]
#         PC_spill_return = self.ins[2]
#         CWT_spill_return = self.ins[3]

#         mixed_in = self._mixed_in
#         # Define input streams
#         input_streams = [WasteWater, sludge_return, PC_spill_return, CWT_spill_return]
#         mixed_in.mix_from(input_streams)

#         # Output stream
#         TreatedWater = self.outs[0]
        
#         # Copy the mixed result to the outflow
#         TreatedWater.copy_like(mixed_in)

#     def _design(self):
#         design = self.design_results
#         constr = self.construction
#         design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
#         self.add_construction(add_cost=False)
    
#     def _cost(self):
#         C = self.baseline_purchase_costs
#         C['Tank'] = self.collection_tank_cost
#         C['Pipes'] = self.pipeline_connectors
#         C['Fittings'] = self.weld_female_adapter_fittings
    
#         ratio = self.price_ratio
#         for equipment, cost in C.items():
#             C[equipment] = cost * ratio
        
#         self.add_OPEX = self._calc_replacement_cost()
        
#         power_demand = self.power_demand_CT
#         self.power_utility(power_demand)
    
#     def _calc_replacement_cost(self):
#         scale  = (self.ppl / self.baseline_ppl) ** self.exponent_scale
#         CT_replacement_cost = (
#             self.collection_tank_cost / self.collection_tank_lifetime +               
#             self.pipeline_connectors / self.pipeline_connectors_lifetime +
#             self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime) * scale
#         CT_replacement_cost = CT_replacement_cost / (365 * 24)  # convert to USD/hr
#         return CT_replacement_cost

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
        Mixer.__init__(self, ID, ins, outs, thermo, init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        self.rigorous = rigorous
        self.conserve_phases = conserve_phases


    @property
    def state(self):
        '''The state of the Mixer, including component concentrations [mg/L] and flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    def _init_state(self):
        '''initialize state by specifying or calculating component concentrations
        based on influents. Total flow rate is always initialized as the sum of
        influent wastestream flows.'''
        QCs = self._ins_QC
        if QCs.shape[0] <= 1: self._state = QCs[0]
        else:
            Qs = QCs[:,-1]
            Cs = QCs[:,:-1]
            self._state = np.append(Qs @ Cs / Qs.sum(), Qs.sum())
        self._dstate = self._state * 0.

    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Mixer'''
        self._outs[0].state = self._state

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Mixer'''
        self._outs[0].dstate = self._dstate

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _n_ins = len(self.ins)
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        def yt(t, QC_ins, dQC_ins):
            if _n_ins > 1:
                Q_ins = QC_ins[:, -1]
                C_ins = QC_ins[:, :-1]
                dQ_ins = dQC_ins[:, -1]
                dC_ins = dQC_ins[:, :-1]
                Q = Q_ins.sum()
                C = Q_ins @ C_ins / Q
                _state[-1] = Q
                _state[:-1] = C
                Q_dot = dQ_ins.sum()
                C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
                _dstate[-1] = Q_dot
                _dstate[:-1] = C_dot
            else:
                _state[:] = QC_ins[0]
                _dstate[:] = dQC_ins[0]
            _update_state()
            _update_dstate()
        self._AE = yt

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
# PrimaryClarifier_path = ospath.join(EL_su_data_path, '_EL_PC.tsv')

# @price_ratio()
# class EL_PC(SanUnit):
#     """
#     Primary clarifier in the Enviroloo (EL) Clear Toilet system for COD and suspended solids removal.

#     Parameters
#     ----------
#     ID : str
#         Unique identifier for the unit.
#     ins : tuple
#         Input streams: (0) Wastewater from lift pump, (1) Nitrate return from membrane tank.
#     outs : tuple
#         Output streams: (0) Treated water, (1) Spill return, (2) Sludge return.
#     sludge_flow_rate : float
#         Sludge flow rate (m³/d). Required for sludge return calculation.
#     solids_removal_efficiency : float
#         Fraction of suspended solids removed (0 to 1). Default is 0.5.
#     max_overflow : float
#         Maximum allowable overflow rate (m³/h). Default is 15.
#     ppl : float
#         Current population served.
#     baseline_ppl : float
#         Baseline population for scaling design and cost.

#     Notes
#     -----
#     - COD and suspended solids are tracked via the `components` object.
#     - Spill return occurs if treated water flow exceeds `max_overflow`.
#     - Inherits from `SanUnit`, not `IdealClarifier`, for flexibility.
#     """

#     _N_ins = 2
#     _N_outs = 3
#     _ins_size_is_fixed = True
#     _outs_size_is_fixed = True
#     exponent_scale = 0.1

#     def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
#                  sludge_flow_rate=None, solids_removal_efficiency=0.5, max_overflow=15,
#                  F_BM=1, transportation=[],
#                  ppl=None, baseline_ppl=None, **kwargs):
#         SanUnit.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with, 
#                          F_BM_default=F_BM, transportation=transportation)
#         self.sludge_flow_rate = sludge_flow_rate  # m³/d
#         self.solids_removal_efficiency = solids_removal_efficiency  # 0 to 1
#         self.max_overflow = max_overflow  # m³/h
#         self.ppl = ppl
#         self.baseline_ppl = baseline_ppl

#         # Load default parameters from data file
#         data = load_data(path=PrimaryClarifier_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             setattr(self, para, value)
#         del data

#         for attr, value in kwargs.items(): 
#             setattr(self, attr, value)

#     def _init_lca(self):
#         self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]

#     def _run(self):
#         """Simulate the primary clarifier operation."""
#         # Inputs
#         wastewater, nitrate_return = self.ins
#         # Outputs
#         treated_water, sludge_return, spill_return = self.outs

#         # Mix input streams
#         mixed = wastewater.copy()  # Temporary stream for mixing
#         mixed.mix_from([wastewater, nitrate_return])

#         # Get total inflow (m³/d) and TSS (mg/L)
#         Q_in = mixed.F_vol * 24  # Convert m³/h to m³/d
#         TSS_in = sludge_return.F_mass # assume all return flow only including 
#                                       # sewage sludge that was treated as TSS
#         TSS_in = mixed.get_TSS()  # Total suspended solids in mg/L

#         # Handle case with no solids
#         if TSS_in <= 0:
#             treated_water.copy_like(mixed)
#             sludge_return.empty()
#             spill_return.empty()
#             return

#         # Validate parameters
#         if not self.sludge_flow_rate or not self.solids_removal_efficiency:
#             raise ValueError("Must specify 'sludge_flow_rate' and 'solids_removal_efficiency'.")

#         # Calculate sludge split
#         Qs = self.sludge_flow_rate  # m³/d
#         e_rmv = self.solids_removal_efficiency
#         f_Qu = Qs / Q_in  # Fraction of flow to sludge
#         f_Xu = e_rmv + (1 - e_rmv) * f_Qu  # Fraction of solids to sludge

#         # Component-wise split (simplified for COD and SS)
#         cmps = self.components
#         split_to_sludge = np.zeros(len(cmps))
#         for i, cmp in enumerate(cmps):
#             if cmp.ID == 'OtherSS':  # Suspended solids
#                 split_to_sludge[i] = f_Xu
#             elif cmp.ID == 'S_COD':  # Soluble COD, assume some removal with solids
#                 split_to_sludge[i] = f_Qu * 0.1  # Arbitrary small fraction
#             else:
#                 split_to_sludge[i] = f_Qu  # Other components follow flow split
#         split_to_sludge = np.clip(split_to_sludge, 0, 1)

#         # Split mixed stream
#         mixed.split_to(sludge_return, treated_water, split_to_sludge)

#         # Handle spill return
#         # if self.max_overflow is not None:
#         #     if treated_water.F_vol > self.max_overflow:
#         #         # Spill return exists
#         #         spill_vol = treated_water.F_vol - self.max_overflow  
#         #         if not hasattr(self, '_f_spill'):
#         #             self._f_spill = None
#         #         if not hasattr(self, '_f_overflow'):
#         #             self._f_overflow = None
#         #             self._f_spill = spill_vol / treated_water.F_vol  
#         #             self._f_overflow = 1 - self._f_spill  

#         #             spill_return.copy_like(treated_water)  
#         #             spill_return.F_mass *= self._f_spill  

#         #             TreatedWater.F_mass *= self._f_overflow  
#         #     else:
#         #         # max_overflow is not none, but TreatedWater < max_overflow
#         #         spill_return.empty()
#         #         if hasattr(self, '_f_spill'):
#         #             del self._f_spill  
#         #         if hasattr(self, '_f_overflow'):
#         #             del self._f_overflow  
#         # else:
#         #     # max_overflow is none, no spill return
#         #     spill_return.empty()
        
#         max_overflow_m3d = self.max_overflow * 24  # Convert m³/h to m³/d
#         Q_treated = treated_water.F_vol * 24  # m³/d
#         if Q_treated > max_overflow_m3d:
#             spill_vol = Q_treated - max_overflow_m3d
#             f_spill = spill_vol / Q_treated
#             spill_return.copy_like(treated_water)
#             spill_return.F_mass *= f_spill
#             treated_water.F_mass *= (1 - f_spill)
#         else:
#             spill_return.empty()

#     def _design(self):
#         """Calculate design parameters."""
#         self.design_results['StainlessSteel'] = (
#             self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
#         )
#         self.construction = [
#             Construction(item='StainlessSteel', quantity=self.design_results['StainlessSteel'], quantity_unit='kg')
#         ]
#         self.add_construction(add_cost=False)

#     def _cost(self):
#         """Calculate capital and operating costs."""
#         C = self.baseline_purchase_costs
#         C['Tank'] = self.PC_tank_cost
#         C['Pipes'] = self.pipeline_connectors
#         C['Fittings'] = self.weld_female_adapter_fittings

#         ratio = self.price_ratio  # Assume price_ratio decorator sets this
#         for equipment, cost in C.items():
#             C[equipment] = cost * ratio

#         self.add_OPEX = self._calc_replacement_cost()
#         power_demand = self.power_demand_PC  # Default to 0 if not set
#         self.power_utility(power_demand)

#     def _calc_replacement_cost(self):
#         """Calculate replacement cost in USD/hr."""
#         scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
#         replacement_cost = (
#             self.PC_tank_cost / self.PC_tank_lifetime +
#             self.pipeline_connectors / self.pipeline_connectors_lifetime +
#             self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime
#         ) * scale
#         return replacement_cost / (365 * 24)  # Convert to USD/hr

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
    _N_outs = 3  # [0] effluent overflow, [1] sludge underflow
    _outs_size_is_fixed = True

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 sludge_flow_rate=2000, solids_removal_efficiency=0.995,
                 sludge_MLSS=None, isdynamic=False, init_with='WasteStream',
                 F_BM_default=None, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic,
                         init_with=init_with, F_BM_default=F_BM_default)
        self.sludge_flow_rate = sludge_flow_rate
        self.solids_removal_efficiency = solids_removal_efficiency
        self.sludge_MLSS = sludge_MLSS
        self._mixed = WasteStream()
        self._f_uf = None
        self._f_of = None

    @property
    def sludge_flow_rate(self):
        '''[float] The designed sludge flow rate (wasted + recycled) in m3/d.'''
        return self._Qs

    @sludge_flow_rate.setter
    def sludge_flow_rate(self, Qs):
        self._Qs = Qs

    @property
    def solids_removal_efficiency(self):
        return self._e_rmv

    @solids_removal_efficiency.setter
    def solids_removal_efficiency(self, f):
        if f is not None and (f > 1 or f < 0):
            raise ValueError(f'solids removal efficiency must be within [0, 1], not {f}')
        self._e_rmv = f

    @property
    def sludge_MLSS(self):
        return self._MLSS

    @sludge_MLSS.setter
    def sludge_MLSS(self, MLSS):
        if MLSS is not None:
            warn(f'sludge MLSS {MLSS} mg/L is only used to estimate '
                 f'sludge flowrate or solids removal efficiency, when either '
                 f'one of them is unspecified.')
        self._MLSS = MLSS

    def _run(self):
        inf = self._mixed
        inf.mix_from(self.ins)
        of, uf, spill_PC = self.outs
        TSS_in = inf.get_TSS()
        if TSS_in <= 0:
            uf.empty()
            of.copy_like(inf)
        else:
            Q_in = inf.F_vol * 24 # m3/d
            x = inf.components.x
            Qs, e_rmv, mlss = self._Qs, self._e_rmv, self._MLSS
            if Qs and e_rmv:
                f_Qu = Qs/Q_in
                f_Xu = e_rmv + (1-e_rmv) * f_Qu
            elif Qs and mlss:
                f_Qu = Qs/Q_in
                f_Xu = f_Qu*mlss/TSS_in
            elif e_rmv and mlss:
                f_Qu = e_rmv / (mlss/TSS_in - (1-e_rmv))
                f_Xu = e_rmv + (1-e_rmv) * f_Qu
            split_to_uf = (1-x)*f_Qu + x*f_Xu
            if any(split_to_uf > 1): split_to_uf = 1
            inf.split_to(uf, of, split_to_uf)

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
   
# Total COD removal efficiency
nCOD = lambda f_corr, fx, HRT: f_corr*(2.88*fx - 0.118)*(1.45 + 6.15*np.log(HRT*24*60))

def calc_f_i(fx, f_corr, HRT):
    '''calculates the effluent-to-influent ratio of solid concentrations'''
    nX = nCOD(f_corr, fx, HRT)/fx
    if nX > 100: nX = 100
    if nX < 0: nX = 0
    return 1-(nX/100)

    def _design(self):
        """Calculate design parameters."""
        self.design_results['StainlessSteel'] = (
            self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        )
        self.construction = [
            Construction(item='StainlessSteel', quantity=self.design_results['StainlessSteel'], quantity_unit='kg')
        ]
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
# Anoxic_path = ospath.join(EL_su_data_path, '_EL_Anoxic.tsv')

# @price_ratio()
# class EL_Anoxic(SanUnit, Decay):
#     '''
#     Anoxic treatment unit in the EL system.

#     Parameters
#     ----------
#     Ins:
#     (1) influent of treated wastetwater from primary clarifier
#     (2) nitrate return flow from membrane tank
#     (3) glucose addition
#     (4) agitation pump

#     Outs:
#     (1) effluent of treated wastewater
#     (2) fugitive CH4 emission
#     (3) fugtivie N2O emission

#     Attributes
#     ----------

    
#     References
#     ----------
#      refer to the qsdsan.sanunits.PrimaryClarifier module
    
    
#     '''
#     _N_ins = 4
#     _N_outs = 3
#     _ins_size_is_fixed = True
#     _outs_size_is_fixed = True
#     baseline_ppl = 30
#     exponent_scale = 0.1

#     def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', methane_density=0.657,
#                  degraded_components=('OtherSS',),  F_BM_default=1, ppl = None, baseline_ppl = None,
#                  glucose_density=1560, if_capture_biogas=False, if_N2O_emission=True, **kwargs):
#         Decay.__init__(self, ID, ins, outs, thermo=thermo,
#                        init_with=init_with, F_BM_default=F_BM_default,
#                        degraded_components=degraded_components,
#                        if_capture_biogas=if_capture_biogas,
#                        if_N2O_emission=if_N2O_emission,)
        
#         self.ppl = ppl
#         self.baseline_ppl = baseline_ppl
#         self.methane_density = methane_density  # kg/m^3
#         self.glucose_density = glucose_density  # kg/m^3
#         self._mixed_in = WasteStream(f'{self.ID}_mixed_in')

#         data = load_data(path=Anoxic_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             setattr(self, para, value)
#         del data

#         for attr, value in kwargs.items():
#             setattr(self, attr, value)
#     def _init_lca(self):
#         self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]      
        
#     def _run(self):

#           # Input stream
#           WasteWater = self.ins[0]
#           nitrate_return = self.ins[1]  # Nitrate from membrane tank over return pump
          
#           glucose = self.ins[2]  # Extra carbon source
#           glucose_consumed = WasteWater.F_mass * self.chemical_glucose_dosage 
#           #TODO: this dosage is the same as PAC. Check data source
#           glucose.imass['Glucose'] = glucose_consumed #kg/h
#           # glucose.mass = glucose_consumed
          
#           agitation = self.ins[3]  # Agitation pump works here, but it does not attend mass flow balance calculation
#           agitation.empty()
          
#           mixed_in = self._mixed_in
#           # Define input streams
#           input_streams = [WasteWater, nitrate_return, glucose, agitation]
#           mixed_in.mix_from(input_streams)
          
#           # Output stream
#           TreatedWater = self.outs[0]
#           CH4_emission = self.outs[1]
#           N2O_emission = self.outs[2]

#           # Output stream
#           # TreatedWater = self.outs[0]
#           # Copy the mixed result to the outflow
#           TreatedWater.copy_like(mixed_in)
#           # Inherited input stream
#           # TreatedWater.copy_like(WasteWater)
          
#           # Manually add return nitrate (NO3)
#           TreatedWater.imass['NO3'] += nitrate_return.imass['NO3']
          
#           # Glucose consumption
#           # glucose_consumed = WasteWater.F_vol * self.chemical_glucose_dosage * self.glucose_density
#           # glucose.imass['Glucose'] = glucose_consumed
#           TreatedWater.imass['Glucose'] += glucose.imass['Glucose']
          
#           # glucose.mass = glucose_consumed

#           # COD removal
#           COD_removal = self.EL_anoT_COD_removal
#           removed_COD = (WasteWater.COD + glucose.COD) / 1e3 * WasteWater.F_vol * COD_removal  # kg/hr
          
#           # Now we explicitly remove Glucose from TreatedWater
#           TreatedWater.imass['Glucose'] = 0  # All glucose is consumed in reaction
#           glucose.empty()
          
#           # Sludge produced
#           sludge_prcd = self.EL_anoT_sludge_yield * removed_COD  # produced biomass
          
#           for component in ('Mg', 'Ca', 'OtherSS', 'Tissue', 'WoodAsh'):
#               TreatedWater.imass[component] += sludge_prcd  # New sludge produced
          
#           # CH4 produced
#           CH4_prcd = removed_COD * self.EL_anoT_methane_yield * self.methane_density  # kg CH4 produced/hr
#           CH4_soluble = self.EL_anoT_soluble_methane_fraction * CH4_prcd
#           CH4_emitted = CH4_prcd - CH4_soluble
#           CH4_emission.imass['SolubleCH4'] = CH4_emitted
#           TreatedWater.imass['SolubleCH4'] = CH4_soluble
          
#           # N2O produced
#           N_removal = self.EL_anoT_TN_removal
#           if self.if_N2O_emission:
#           # Assume the removal part covers both degradation loss
#           # and other unspecified removal as well
#               N_loss = self.first_order_decay(k = self.decay_k_N,
#                                     t = self.EL_anoT_tau / 365,
#                                     max_decay = self.N_max_decay)
#               if N_loss > N_removal:
#                   warn(f'Nitrogen degradation loss ({N_loss:.2f}) larger than the given removal ({N_removal:.2f})), '
#                         'the degradation loss is assumed to be the same as the removal.')
#                   N_loss = N_removal
            
#               # N2O only from the degraded part
#               N_loss_tot = N_loss * WasteWater.TN / 1e3 * WasteWater.F_vol
#               N2O_emission.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44 / 28
#           else:
#               N2O_emission.empty()
#               N2O_emission.imass['N2O'] = 0  # All N2O is emitted
              
#           # Assume all NO3 is consumed and does not appear in TreatedWater
#           TreatedWater.imass['NO3']  = 0  
          
#           # NH3 & NonNH3, P, K calculation
#           total_solubles = np.array([
#               CH4_soluble,
#               WasteWater.imass['NH3'] * (1 - N_loss),  # In water or sludge
#               WasteWater.imass['NonNH3'] * (1 - N_loss),  # In water or sludge
#               0,  # NO3 tolly consumed
#               WasteWater.imass['P'],  # In water or sludge
#               WasteWater.imass['K'],  # In water or sludge
#               ])
          
#           # Removed solubles in the sludge, assume minimal used for growth
#           sludge_solubles = total_solubles * np.array([
#                             0,  # CH₄ will not go into sludge
#                             N_removal,  # Part of it in water and in sludge
#                             N_removal,  # Part of it in water and in sludge
#                             0,  # NO3 tolly consumed
#                             self.EL_anoT_TP_removal,  # Part of it in water and in sludge
#                             0,  # K will not go into sludge
#                             ])
          
#           # Calculate solutes entering TreatedWater
#           liquid_solubles = total_solubles - sludge_solubles
          
#           # Assign to outputs (TreatedWater now includes sludge)
#           TreatedWater.imass['SolubleCH4', 'NH3', 'NonNH3', 'NO3', 'P', 'K'] = liquid_solubles
#           TreatedWater.imass['SolubleCH4', 'NH3', 'NonNH3', 'NO3', 'P', 'K'] += sludge_solubles  # Here is only one output stream (TreatedWater + sludge)
          
#           # Final COD
#           TreatedWater._COD = WasteWater.COD * (1 - self.EL_anoT_COD_removal)
          
    # def _run(self):
        
    #     # Input stream
    #     WasteWater = self.ins[0]
    #     sludge_return = self.ins[1]
    #     glucose = self.ins[2]
    #     agiation = self.ins[3]
        
    #     # Output stream
    #     TreatedWater = self.outs[0]
    #     CH4_emission = self.outs[1]
    #     N2O_emission = self.outs[2]
        
    #     input_streams = [WasteWater, sludge_return, glucose]
            
    #     # Mix all inputs into a single stream
    #     self._mixed.empty()
    #     self._mixed.mix_from(input_streams)

    #     # Copy the mixed result to the outflow
    #     TreatedWater.copy_like(self._mixed)
    
   
    # def _run(self):
        
    #     # Input stream
    #     WasteWater = self.ins[0]
    #     sludge_return = self.ins[1]
    #     glucose = self.ins[2]

    #     # Output stream
    #     TreatedWater = self.outs[0]
    #     CH4_emission = self.outs[1]
    #     N2O_emission = self.outs[2]

    #     # Combined sludge return and influent
    #     WasteWater.F_mass += sludge_return.F_mass
    #     WasteWater.imass += sludge_return.imass

    #     # Glucose in
    #     WasteWater.imass['Glucose'] += glucose.imass['Glucose']
    #     WasteWater.F_mass += glucose.F_mass

    #     # Simulate complete glucose consumption
    #     consumed_glucose = WasteWater.imass['Glucose']
    #     WasteWater.imass['Glucose'] -= consumed_glucose
    #     WasteWater.F_mass -= consumed_glucose

    #     # Adjust N2O emission factor based on reduction ratio
    #     original_N2O_EF = self.N2O_EF_decay
    #     self.N2O_EF_decay *= 0.25

    #     # Call Decay._first_order_run with the updated N2O_EF_decay
    #     super()._first_order_run(waste=WasteWater,
    #                              treated=TreatedWater,
    #                              CH4=CH4_emission,
    #                              N2O=N2O_emission
    #                              )

    #     # Restore the original N2O emission factor
    #     self.N2O_EF_decay = original_N2O_EF

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
    _N_ins = 3
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    _D_O2 = 2.10e-9   # m2/s

    def __init__(self, ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=1000, W_tank = 6.4, D_tank = 3.65,
                 freeboard = 0.61, t_wall = None, t_slab = None, aeration=None, 
                 DO_ID='S_O2', suspended_growth_model=None, 
                 gas_stripping=False, gas_IDs=None, stripping_kLa_min=None, 
                 K_Henry=None, D_gas=None, p_gas_atm=None,
                 isdynamic=True, exogenous_vars=(), **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, isdynamic=isdynamic,
                         exogenous_vars=exogenous_vars, **kwargs)
        self._V_max = V_max
        self._aeration = aeration
        self._DO_ID = DO_ID
        self._model = suspended_growth_model
        self.gas_IDs = gas_IDs
        self.stripping_kLa_min = stripping_kLa_min
        self.K_Henry = K_Henry
        self.D_gas = D_gas
        self.p_gas_atm = p_gas_atm
        self.gas_stripping = gas_stripping
        self._concs = None
        self._mixed = WasteStream()
        self.split = split

        # # Design parameters 
        # self._W_tank = W_tank
        # self._D_tank = D_tank
        # self._freeboard = freeboard
        # self._t_wall = t_wall
        # self._t_slab = t_slab
    
        # for attr, value in kwargs.items():
        #     setattr(self, attr, value)
    
    @property
    def V_max(self):
        '''[float] The designed maximum liquid volume, not accounting for increased volume due to aeration, in m^3.'''
        return self._V_max

    @V_max.setter
    def V_max(self, Vm):
        self._V_max = Vm
        
    # @property
    # def W_tank(self):
    #     '''[float] The design width of the tank, in m.'''
    #     return self._W_tank

    # @W_tank.setter
    # def W_tank(self, W_tank):
    #     self._W_tank = W_tank
        
    # @property
    # def D_tank(self):
    #     '''[float] The design depth of the tank, in m.'''
    #     return self._D_tank

    # @D_tank.setter
    # def D_tank(self, D_tank):
    #     self._D_tank = D_tank
        
    # @property
    # def freeboard(self):
    #     '''[float] Freeboard added to the depth of the reactor tank, [m].'''
    #     return self._freeboard
    
    # @freeboard.setter
    # def freeboard(self, i):
    #     self._freeboard = i
        
    # @property
    # def t_wall(self):
    #     '''
    #     [float] Thickness of the wall concrete, [m].
    #     default to be minimum of 1 ft with 1 in added for every ft of depth over 12 ft.
    #     '''
    #     D_tank = self.D_tank*39.37 # m to inches 
    #     return self._t_wall or (1 + max(D_tank - 12, 0)/12)*0.3048 # from feet to m
    
    # @t_wall.setter
    # def t_wall(self, i):
    #     self._t_wall = i

    # @property
    # def t_slab(self):
    #     '''
    #     [float] Concrete slab thickness, [m],
    #     default to be 2 in thicker than the wall thickness.
    #     '''
    #     return self._t_slab or (self.t_wall + 2/12)*0.3048 # from feet to m
    
    # @t_slab.setter
    # def t_slab(self, i):
    #     self._t_slab = i
     
    @property
    def aeration(self):
        '''[:class:`Process` or float or NoneType] Aeration model.'''
        return self._aeration

    @aeration.setter
    def aeration(self, ae):
        if ae is None or isinstance(ae, Process): self._aeration = ae
        elif isinstance(ae, (float, int)):
            if ae < 0:
                raise ValueError('targeted dissolved oxygen concentration for aeration must be non-negative.')
            else:
                if ae > 14:
                    warn(f'targeted dissolved oxygen concentration for {self.ID} might exceed the saturated level.')
                self._aeration = ae
        else:
            raise TypeError(f'aeration must be one of the following types: float, '
                            f'int, Process, NoneType. Not {type(ae)}')

    @property
    def suspended_growth_model(self):
        '''[:class:`CompiledProcesses` or NoneType] Suspended growth model.'''
        return self._model

    @suspended_growth_model.setter
    def suspended_growth_model(self, model):
        if isinstance(model, CompiledProcesses) or model is None: self._model = model
        else: raise TypeError(f'suspended_growth_model must be one of the following '
                              f'types: CompiledProesses, NoneType. Not {type(model)}')

    @property
    def DO_ID(self):
        '''[str] The `Component` ID for dissolved oxygen used in the suspended growth model and the aeration model.'''
        return self._DO_ID

    @DO_ID.setter
    def DO_ID(self, doid):
        if doid not in self.components.IDs:
            raise ValueError(f'DO_ID must be in the set of `CompiledComponents` used to set thermo, '
                             f'i.e., one of {self.components.IDs}.')
        self._DO_ID = doid

    @property
    def gas_stripping(self):
        return self._gstrip
    @gas_stripping.setter
    def gas_stripping(self, strip):
        self._gstrip = strip = bool(strip)
        if strip:
            if self.gas_IDs:
                cmps = self.components.IDs
                for i in self.gas_IDs:
                    if i not in cmps:
                        raise RuntimeError((f'gas ID {i} not in component set: {cmps}.'))
            else:
                mdl = self._model
                self.gas_IDs = mdl.gas_IDs
                self.stripping_kLa_min = np.array(mdl.kLa_min)
                self.D_gas = np.array(mdl.D_gas)
                self.K_Henry = np.array(mdl.K_Henry)
                self.p_gas_atm = np.array(mdl.p_gas_atm)
                
    @property
    def split(self):
        '''[numpy.1darray or NoneType] The volumetric split of outs.'''
        return self._split

    @split.setter
    def split(self, split):
        if split is None: self._split = split
        else:
            if len(split) != len(self._outs):
                raise ValueError('split and outs must have the same size')
            self._split = np.array(split)/sum(split)

    @property
    def state(self):
        '''The state of the CSTR, including component concentrations [mg/L] and flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    @state.setter
    def state(self, QCs):
        QCs = np.asarray(QCs)
        if QCs.shape != (len(self.components)+1, ):
            raise ValueError(f'state must be a 1D array of length {len(self.components) + 1},'
                              'indicating component concentrations [mg/L] and total flow rate [m^3/d]')
        self._state = QCs

    def set_init_conc(self, **kwargs):
        '''set the initial concentrations [mg/L] of the CSTR.'''
        self._concs = self.components.kwarray(kwargs)

    def _init_state(self):
        mixed = self._mixed
        Q = mixed.get_total_flow('m3/d')
        if self._concs is not None: Cs = self._concs
        else: Cs = mixed.conc
        self._state = np.append(Cs, Q).astype('float64')
        self._dstate = self._state * 0.

    def _update_state(self):
        arr = self._state
        arr[arr < 1e-16] = 0.
        arr[-1] = sum(ws.state[-1] for ws in self.ins)
        if self.split is None: self._outs[0].state = arr
        else:
            for ws, spl in zip(self._outs, self.split):
                y = arr.copy()
                y[-1] *= spl
                ws.state = y

    def _update_dstate(self):
        arr = self._dstate
        if self.split is None: self._outs[0].dstate = arr
        else:
            for ws, spl in zip(self._outs, self.split):
                y = arr.copy()
                y[-1] *= spl
                ws.dstate = y

    def _run(self):
        '''Only to converge volumetric flows.'''
        mixed = self._mixed # avoid creating multiple new streams
        mixed.mix_from(self.ins)
        Q = mixed.F_vol # m3/hr
        if self.split is None: self.outs[0].copy_like(mixed)
        else:
            for ws, spl in zip(self._outs, self.split):
                ws.copy_like(mixed)
                ws.set_total_flow(Q*spl, 'm3/hr')

    def get_retained_mass(self, biomass_IDs):
        cmps = self.components
        mass = cmps.i_mass * self._state[:len(cmps)]
        return self._V_max * mass[cmps.indices(biomass_IDs)].sum()

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE

    def _init_model(self):
        if self._model is None:
            warn(f'{self.ID} was initialized without a suspended growth model, '
                 f'and thus run as a non-reactive unit')
            r = lambda state_arr: np.zeros(self.components.size)
        else:
            # processes = _add_aeration_to_growth_model(aer, self._model)
            r = self._model.production_rates_eval
        return r

    def _compile_ODE(self):
        isa = isinstance
        aer = self._aeration
        r = self._init_model()

        _dstate = self._dstate
        _update_dstate = self._update_dstate
        V = self._V_max
        gstrip = self.gas_stripping
        if gstrip:
            gas_idx = self.components.indices(self.gas_IDs)
            if isa(aer, Process): kLa = aer.kLa
            else: kLa = 0.
            S_gas_air = np.asarray(self.K_Henry)*np.asarray(self.p_gas_atm)
            kLa_stripping = np.maximum(kLa*self.D_gas/self._D_O2, self.stripping_kLa_min)
        hasexo = bool(len(self._exovars))
        f_exovars = self.eval_exo_dynamic_vars
        
        if isa(aer, (float, int)):
            i = self.components.index(self._DO_ID)
            fixed_DO = self._aeration
            def dy_dt(t, QC_ins, QC, dQC_ins):
                QC[i] = fixed_DO
                dydt_cstr(QC_ins, QC, V, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC)
                if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
                _dstate[i] = 0
                _update_dstate()
        elif isa(aer, Process):
            aer_stoi = aer._stoichiometry
            aer_frho = aer.rate_function
            def dy_dt(t, QC_ins, QC, dQC_ins):
                dydt_cstr(QC_ins, QC, V, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC) + aer_stoi * aer_frho(QC)
                if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
                _update_dstate()
        else:
            def dy_dt(t, QC_ins, QC, dQC_ins):
                dydt_cstr(QC_ins, QC, V, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC)
                if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
                _update_dstate()

        self._ODE = dy_dt
         
    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]      
        

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # assume linear scale
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        massflow_anoxic = self.ins[0].mass
        C['Tank'] = self.anoxic_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
        C['Chemcial_glucose'] = self.chemical_glucose_dosage * massflow_anoxic * self.chemical_glucose_price  # make sense the unit of treated water flow

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_AnoxicTank
        self.power_utility(power_demand)
    
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        Anoxic_tank_replacement_cost = (self.anoxic_tank_cost /self.anoxic_tank_lifetime +
                                        self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
                                        self.pipeline_connectors / self.pipeline_connectors_lifetime) * scale
        Anoxic_tank_replacement_cost = Anoxic_tank_replacement_cost / (365 * 24)  # convert to USD/hr
        return Anoxic_tank_replacement_cost

# %%
# Aerobic_path = ospath.join(EL_su_data_path, '_EL_Aerobic.tsv')

# @price_ratio()
# class EL_Aerobic(SanUnit, Decay):

#     '''
#     Aerobic treatment unit in the EL system.

#     Parameters
#     ----------
#     Ins:
#     (1) influent of treated wastetwater from anoxic tank
#     (2) PAC addition
#     (3) blower

#     Outs:
#     (1) effluent of treated wastewater
#     (2) fugitive CH4 emission
#     (3) fugtivie N2O emission

#     Attributes
#     ----------

    
#     References
#     ----------
#      refer to the qsdsan.sanunits.PrimaryClarifier module
#     '''

#     _N_ins = 3  # treated water, PAC, blower
#     _N_outs = 3  # treated water, CH4, N2O
#     _ins_size_is_fixed = True
#     _outs_size_is_fixed = True
#     exponent_scale = 0.1
    
#     def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
#                  degraded_components=('OtherSS',), ppl = None, baseline_ppl = None, F_BM_default=1,
#                  if_capture_biogas=False, if_N2O_emission=True,**kwargs):
#         Decay.__init__(self, ID, ins=ins, outs=outs, thermo=thermo,
#                        init_with=init_with, F_BM_default=F_BM_default,
#                        degraded_components=degraded_components,
#                        if_capture_biogas=if_capture_biogas,
#                        if_N2O_emission=if_N2O_emission,)

#         self.ppl = ppl
#         self.baseline_ppl = baseline_ppl

#         data = load_data(path=Aerobic_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             setattr(self, para, value)
#         del data

#         for attr, value in kwargs.items():
#             setattr(self, attr, value)

#     def _init_lca(self):
#         self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]  
    
#     def _run(self):

#         # Input streams
#         WasteWater = self.ins[0]
#         PAC = self.ins[1]
#         air = self.ins[2]
        
#         # Output streams
#         TreatedWater = self.outs[0]
#         CH4_emission = self.outs[1] 
#         N2O_emission = self.outs[2]
        
#         # Inherited input stream
#         TreatedWater.copy_like(WasteWater)
        
#         # Manually add return nitrate (NO3)
#         # TreatedWater.imass['PAC'] += PAC.imass['PAC']
#         PAC.imass['PAC'] = WasteWater.imass['H2O'] * self.chemical_PAC_dosage/100
        
#         # COD removal
#         COD_removal = self.EL_aeroT_COD_removal
#         removed_COD = WasteWater.COD / 1e3 * WasteWater.F_vol * COD_removal  # kg/hr
        
#         # Sludge produced
#         sludge_prcd = self.EL_aeroT_sludge_yield * removed_COD  # produced biomass
          
#         for component in ('Mg', 'Ca', 'OtherSS', 'Tissue', 'WoodAsh'):
#             TreatedWater.imass[component] += sludge_prcd  # New sludge produced
            
#         # CH4 emission
#         CH4_emission.imass['CH4'] += TreatedWater.imass['SolubleCH4']  # Let all soluble CH4 transfer from solution phase to gas phase
#         TreatedWater.imass['SolubleCH4'] = 0  # Ensure that treated water will not include soluble CH4
        
#         # N2O produced
#         N_removal = self.EL_aeroT_TN_removal
#         if self.if_N2O_emission:
#           # Assume the removal part covers both degradation loss
#           # and other unspecified removal as well
#               N_loss = self.first_order_decay(k = self.decay_k_N,
#                                     t = self.EL_aeroT_tau / 365,
#                                     max_decay = self.N_max_decay)
#               if N_loss > N_removal:
#                   warn(f'Nitrogen degradation loss ({N_loss:.2f}) larger than the given removal ({N_removal:.2f})), '
#                         'the degradation loss is assumed to be the same as the removal.')
#                   N_loss = N_removal
            
#               # N2O only from the degraded part
#               N_loss_tot = N_loss * WasteWater.TN / 1e3 * WasteWater.F_vol
#               N2O_emission.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44 / 28
#         else:
#               N2O_emission.empty()
#               N2O_emission.imass['N2O'] = 0  # All N2O is emitted
        
#         # NO3 conversion
#         NH3_mass = TreatedWater.imass['NH3']  # Inherite NH4 property from anoxic tank
#         NO3_mass_generated = NH3_mass * self.NO3_produced_ratio
#         TreatedWater.imass['NO3'] += NO3_mass_generated
#         TreatedWater.imass['NH3'] = 0
        
#         # N2O emission
#         N2O_mass_generated = NH3_mass * (1 - self.NO3_produced_ratio)
#         N2O_emission.imass['N2O'] += N2O_mass_generated
    
#     # def _run(self):

#     #     # Input streams
#     #     WasteWater = self.ins[0]
#     #     PAC = self.ins[1]
#     #     air = self.ins[2]

#     #     # Output stream
#     #     TreatedWater = self.outs[0] 
#     #     CH4_emission = self.outs[1] 
#     #     N2O_emission = self.outs[2]
        
#     #     input_streams = [WasteWater, PAC, air]
            
#     #     # Mix all inputs into a single stream
#     #     self._mixed.empty()
#     #     self._mixed.mix_from(input_streams)

#     #     # Copy the mixed result to the outflow
#     #     TreatedWater.copy_like(self._mixed)


#     def _design(self):
#         design = self.design_results
#         constr = self.construction
#         design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
#         self.add_construction(add_cost=False)
    
#     def _cost(self):
#         C = self.baseline_purchase_costs
#         # massflow_aerobic = self.ins[0].mass #kg/h
#         C['Tank'] = self.aerobic_tank_cost
#         C['Pipes'] = self.pipeline_connectors
#         C['Fittings'] = self.weld_female_adapter_fittings
#         # C['Chemical_PAC'] = self.chemical_PAC_dosage * massflow_aerobic * self.chemical_PAC_price
#         #PAC cost is already accounted in WasteStream

#         ratio = self.price_ratio
#         for equipment, cost in C.items():
#             C[equipment] = cost * ratio
        
#         self.add_OPEX = self._calc_replacement_cost()
        
#         power_demand = self.power_demand_AerobicTank
#         self.power_utility(power_demand) # kWh
    
#     def _calc_replacement_cost(self):
#         scale = (self.ppl / self.baseline_ppl) * self.exponent_scale
#         Aerobic_tank_replacement_cost = (self.aerobic_tank_cost / self.aerobic_tank_lifetime +
#                                         self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
#                                         self.pipeline_connectors / self.pipeline_connectors_lifetime) * scale
#         Aerobic_tank_replacement_cost = Aerobic_tank_replacement_cost / (365 * 24)  # convert to USD/hr
#         return Aerobic_tank_replacement_cost

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
    _N_ins = 2 # treated water, PAC, blower
    _N_outs = 3  # treated water, CH4, N2O
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    # exponent_scale = 0.1
    
    _D_O2 = 2.10e-9   # m2/s

    def __init__(self, ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=1000, W_tank = 6.4, D_tank = 3.65,
                 freeboard = 0.61, t_wall = None, t_slab = None, aeration=2.0, 
                 DO_ID='S_O2', suspended_growth_model=None, 
                 gas_stripping=False, gas_IDs=None, stripping_kLa_min=None, 
                 K_Henry=None, D_gas=None, p_gas_atm=None,
                 isdynamic=True, exogenous_vars=(), **kwargs):
        CSTR.__init__(self, ID, ins, outs, thermo, init_with, isdynamic=isdynamic,
                         exogenous_vars=exogenous_vars, **kwargs)
        self._V_max = V_max
        self._aeration = aeration
        self._DO_ID = DO_ID
        self._model = suspended_growth_model
        self.gas_IDs = gas_IDs
        self.stripping_kLa_min = stripping_kLa_min
        self.K_Henry = K_Henry
        self.D_gas = D_gas
        self.p_gas_atm = p_gas_atm
        self.gas_stripping = gas_stripping
        self._concs = None
        self._mixed = WasteStream()
        self.split = split

        # # Design parameters 
        # self._W_tank = W_tank
        # self._D_tank = D_tank
        # self._freeboard = freeboard
        # self._t_wall = t_wall
        # self._t_slab = t_slab
    
        # for attr, value in kwargs.items():
        #     setattr(self, attr, value)
    
    @property
    def V_max(self):
        '''[float] The designed maximum liquid volume, not accounting for increased volume due to aeration, in m^3.'''
        return self._V_max

    @V_max.setter
    def V_max(self, Vm):
        self._V_max = Vm
        
    # @property
    # def W_tank(self):
    #     '''[float] The design width of the tank, in m.'''
    #     return self._W_tank

    # @W_tank.setter
    # def W_tank(self, W_tank):
    #     self._W_tank = W_tank
        
    # @property
    # def D_tank(self):
    #     '''[float] The design depth of the tank, in m.'''
    #     return self._D_tank

    # @D_tank.setter
    # def D_tank(self, D_tank):
    #     self._D_tank = D_tank
        
    # @property
    # def freeboard(self):
    #     '''[float] Freeboard added to the depth of the reactor tank, [m].'''
    #     return self._freeboard
    
    # @freeboard.setter
    # def freeboard(self, i):
    #     self._freeboard = i
        
    # @property
    # def t_wall(self):
    #     '''
    #     [float] Thickness of the wall concrete, [m].
    #     default to be minimum of 1 ft with 1 in added for every ft of depth over 12 ft.
    #     '''
    #     D_tank = self.D_tank*39.37 # m to inches 
    #     return self._t_wall or (1 + max(D_tank - 12, 0)/12)*0.3048 # from feet to m
    
    # @t_wall.setter
    # def t_wall(self, i):
    #     self._t_wall = i

    # @property
    # def t_slab(self):
    #     '''
    #     [float] Concrete slab thickness, [m],
    #     default to be 2 in thicker than the wall thickness.
    #     '''
    #     return self._t_slab or (self.t_wall + 2/12)*0.3048 # from feet to m
    
    # @t_slab.setter
    # def t_slab(self, i):
    #     self._t_slab = i
     
    @property
    def aeration(self):
        '''[:class:`Process` or float or NoneType] Aeration model.'''
        return self._aeration

    @aeration.setter
    def aeration(self, ae):
        if ae is None or isinstance(ae, Process): self._aeration = ae
        elif isinstance(ae, (float, int)):
            if ae < 0:
                raise ValueError('targeted dissolved oxygen concentration for aeration must be non-negative.')
            else:
                if ae > 14:
                    warn(f'targeted dissolved oxygen concentration for {self.ID} might exceed the saturated level.')
                self._aeration = ae
        else:
            raise TypeError(f'aeration must be one of the following types: float, '
                            f'int, Process, NoneType. Not {type(ae)}')

    @property
    def suspended_growth_model(self):
        '''[:class:`CompiledProcesses` or NoneType] Suspended growth model.'''
        return self._model

    @suspended_growth_model.setter
    def suspended_growth_model(self, model):
        if isinstance(model, CompiledProcesses) or model is None: self._model = model
        else: raise TypeError(f'suspended_growth_model must be one of the following '
                              f'types: CompiledProesses, NoneType. Not {type(model)}')

    @property
    def DO_ID(self):
        '''[str] The `Component` ID for dissolved oxygen used in the suspended growth model and the aeration model.'''
        return self._DO_ID

    @DO_ID.setter
    def DO_ID(self, doid):
        if doid not in self.components.IDs:
            raise ValueError(f'DO_ID must be in the set of `CompiledComponents` used to set thermo, '
                             f'i.e., one of {self.components.IDs}.')
        self._DO_ID = doid

    @property
    def gas_stripping(self):
        return self._gstrip
    @gas_stripping.setter
    def gas_stripping(self, strip):
        self._gstrip = strip = bool(strip)
        if strip:
            if self.gas_IDs:
                cmps = self.components.IDs
                for i in self.gas_IDs:
                    if i not in cmps:
                        raise RuntimeError((f'gas ID {i} not in component set: {cmps}.'))
            else:
                mdl = self._model
                self.gas_IDs = mdl.gas_IDs
                self.stripping_kLa_min = np.array(mdl.kLa_min)
                self.D_gas = np.array(mdl.D_gas)
                self.K_Henry = np.array(mdl.K_Henry)
                self.p_gas_atm = np.array(mdl.p_gas_atm)
                
    @property
    def split(self):
        '''[numpy.1darray or NoneType] The volumetric split of outs.'''
        return self._split

    @split.setter
    def split(self, split):
        if split is None: self._split = split
        else:
            if len(split) != len(self._outs):
                raise ValueError('split and outs must have the same size')
            self._split = np.array(split)/sum(split)

    @property
    def state(self):
        '''The state of the CSTR, including component concentrations [mg/L] and flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    @state.setter
    def state(self, QCs):
        QCs = np.asarray(QCs)
        if QCs.shape != (len(self.components)+1, ):
            raise ValueError(f'state must be a 1D array of length {len(self.components) + 1},'
                              'indicating component concentrations [mg/L] and total flow rate [m^3/d]')
        self._state = QCs

    def set_init_conc(self, **kwargs):
        '''set the initial concentrations [mg/L] of the CSTR.'''
        self._concs = self.components.kwarray(kwargs)

    def _init_state(self):
        mixed = self._mixed
        Q = mixed.get_total_flow('m3/d')
        if self._concs is not None: Cs = self._concs
        else: Cs = mixed.conc
        self._state = np.append(Cs, Q).astype('float64')
        self._dstate = self._state * 0.

    def _update_state(self):
        arr = self._state
        arr[arr < 1e-16] = 0.
        arr[-1] = sum(ws.state[-1] for ws in self.ins)
        if self.split is None: self._outs[0].state = arr
        else:
            for ws, spl in zip(self._outs, self.split):
                y = arr.copy()
                y[-1] *= spl
                ws.state = y

    def _update_dstate(self):
        arr = self._dstate
        if self.split is None: self._outs[0].dstate = arr
        else:
            for ws, spl in zip(self._outs, self.split):
                y = arr.copy()
                y[-1] *= spl
                ws.dstate = y

    def _run(self):
        '''Only to converge volumetric flows.'''
        mixed = self._mixed # avoid creating multiple new streams
        mixed.mix_from(self.ins)
        Q = mixed.F_vol # m3/hr
        if self.split is None: self.outs[0].copy_like(mixed)
        else:
            for ws, spl in zip(self._outs, self.split):
                ws.copy_like(mixed)
                ws.set_total_flow(Q*spl, 'm3/hr')

    def get_retained_mass(self, biomass_IDs):
        cmps = self.components
        mass = cmps.i_mass * self._state[:len(cmps)]
        return self._V_max * mass[cmps.indices(biomass_IDs)].sum()

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE

    def _init_model(self):
        if self._model is None:
            warn(f'{self.ID} was initialized without a suspended growth model, '
                 f'and thus run as a non-reactive unit')
            r = lambda state_arr: np.zeros(self.components.size)
        else:
            # processes = _add_aeration_to_growth_model(aer, self._model)
            r = self._model.production_rates_eval
        return r

    def _compile_ODE(self):
        isa = isinstance
        aer = self._aeration
        r = self._init_model()

        _dstate = self._dstate
        _update_dstate = self._update_dstate
        V = self._V_max
        gstrip = self.gas_stripping
        if gstrip:
            gas_idx = self.components.indices(self.gas_IDs)
            if isa(aer, Process): kLa = aer.kLa
            else: kLa = 0.
            S_gas_air = np.asarray(self.K_Henry)*np.asarray(self.p_gas_atm)
            kLa_stripping = np.maximum(kLa*self.D_gas/self._D_O2, self.stripping_kLa_min)
        hasexo = bool(len(self._exovars))
        f_exovars = self.eval_exo_dynamic_vars
        
        if isa(aer, (float, int)):
            i = self.components.index(self._DO_ID)
            fixed_DO = self._aeration
            def dy_dt(t, QC_ins, QC, dQC_ins):
                QC[i] = fixed_DO
                dydt_cstr(QC_ins, QC, V, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC)
                if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
                _dstate[i] = 0
                _update_dstate()
        elif isa(aer, Process):
            aer_stoi = aer._stoichiometry
            aer_frho = aer.rate_function
            def dy_dt(t, QC_ins, QC, dQC_ins):
                dydt_cstr(QC_ins, QC, V, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC) + aer_stoi * aer_frho(QC)
                if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
                _update_dstate()
        else:
            def dy_dt(t, QC_ins, QC, dQC_ins):
                dydt_cstr(QC_ins, QC, V, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC)
                if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
                _update_dstate()

        self._ODE = dy_dt
        
    # _units = {
    #     'Tank volume': 'm3',
    #     'Tank width': 'm',
    #     'Tank depth': 'm',
    #     'Tank length': 'm',
    #     'Volume of concrete wall': 'm3',
    #     'Volume of concrete slab': 'm3' 
    #     }
        

    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]  


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        massflow_aerobic = self.ins[0].mass
        C['Tank'] = self.aerobic_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
        C['Chemical_PAC'] = self.chemical_PAC_dosage * massflow_aerobic * self.chemical_PAC_price

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_AerobicTank
        self.power_utility(power_demand) # kWh
    
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) * self.exponent_scale
        Aerobic_tank_replacement_cost = (self.aerobic_tank_cost / self.aerobic_tank_lifetime +
                                        self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
                                        self.pipeline_connectors / self.pipeline_connectors_lifetime) * scale
        Aerobic_tank_replacement_cost = Aerobic_tank_replacement_cost / (365 * 24)  # convert to USD/hr
        return Aerobic_tank_replacement_cost

# %%
# MBR_path = ospath.join(EL_su_data_path, '_EL_MBR.tsv')

# @price_ratio()
# class EL_MBR(SanUnit, Decay):

#     '''
#     Aerobic treatment unit in the EL system.

#     Parameters
#     ----------
#     Ins:
#     (1) influent of treated wastetwater from aerobic tank
#     (2) blower

#     Outs:
#     (1) effluent of treated wastewater
#     (2) nitrate return flow to anoxic tank
#     (3) nitrate return flow to primary clarifier
#     (4) fugitive CH4 emission
#     (5) fugtivie N2O emission

#     Attributes
#     ----------

    
#     References
#     ----------
#      refer to the exposan.eco-san.MBR module
#     '''
#     _N_ins = 2  # treated water, Blower
#     _N_outs = 6  # treated water, CH4, N2O, Nitrate return to Primary Clarifier, Nitrate return to Anoxic Tank
#     _ins_size_is_fixed = True
#     _outs_size_is_fixed = True
#     exponent_scale = 0.1

#     def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
#                  degraded_components=('OtherSS',), ppl = None, baseline_ppl = None, F_BM_default=1,
#                  if_capture_biogas=False, if_N2O_emission=True, **kwargs):
#         Decay.__init__(self, ID, ins=ins, outs=outs, thermo=thermo,
#                        init_with=init_with, F_BM_default=F_BM_default,
#                        degraded_components=degraded_components,
#                        if_capture_biogas=if_capture_biogas,
#                        if_N2O_emission=if_N2O_emission,)

#         self.ppl = ppl
#         self.baseline_ppl = baseline_ppl

#         data = load_data(path=MBR_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             setattr(self, para, value)
#         del data

#         for attr, value in kwargs.items():
#             setattr(self, attr, value)

#     def _init_lca(self):
#         self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),
#                              Construction(item ='PVDF_membrane', linked_unit=self, quantity_unit='kg'),]
    
#     def _run(self):
        
#         # Input stream
#         WasteWater = self.ins[0]
#         air = self.ins[1]
        
#         # Output stream
#         TreatedWater = self.outs[0]
#         PC_nitrate_return = self.outs[1]
#         AnoT_nitrate_return = self.outs[2]
#         CH4_emission = self.outs[3]
#         N2O_emission = self.outs[4]
#         sludge = self.outs[5]
        
#         # Inherited input stream
#         TreatedWater.copy_like(WasteWater)
#         # TreatedWater.empty()
    
#         # # Transfer 99% of components to sludge
#         sludge.empty()
#         for component in ('P','K','NH3','NonNH3','Mg', 'Ca', 'OtherSS', 'Tissue', 'WoodAsh'):
#             mass_in_treated = TreatedWater.imass[component]  # Obtain every components' property
#             mass_to_sludge = 0.99 * mass_in_treated          # Transfer 99% components content to sludge
#             sludge.imass[component] = mass_to_sludge
#             TreatedWater.imass[component] -= mass_to_sludge  # The last components content to treated water
        
#         # COD removal
#         COD_removal = self.EL_mbrT_COD_removal
#         removed_COD = WasteWater.COD / 1e3 * WasteWater.F_vol * COD_removal  # kg/hr
        
#         # Sludge produced
#         sludge_prcd = removed_COD * self.EL_mbrT_sludge_yield  # Typical yield: 0.5 kg VSS/kg COD
            
#         # Sludge handling
#         sludge.copy_flow(TreatedWater, ('P', 'K', 'NH3','NonNH3', 'Mg', 'Ca', 'OtherSS', 'Tissue', 'WoodAsh'), remove=True)
#         sludge.imass['OtherSS'] += sludge_prcd
#         sludge.imass['H2O'] = sludge.imass['OtherSS']/(1-self.sludge_moisture_content) - sludge.imass['OtherSS']
        
#         # Water balance
#         TreatedWater.imass['H2O'] -= sludge.imass['H2O']
#         if TreatedWater.imass['H2O'] <= 0:
#             sludge.imass['H2O'] = WasteWater.imass['H2O']
#             TreatedWater.imass['H2O'] = 0
        
#         # CH4 emission
#         CH4_emission.imass['CH4'] = WasteWater.imass['SolubleCH4']  # Let all soluble CH4 transfer from solution phase to gas phase
#         TreatedWater.imass['SolubleCH4'] = 0  # Ensure that treated water will not include soluble CH4
        
#         # N2O produced
#         N_removal = self.EL_mbrT_TN_removal
#         if self.if_N2O_emission:
#           # Assume the removal part covers both degradation loss
#           # and other unspecified removal as well
#               N_loss = self.first_order_decay(k = self.decay_k_N,
#                                     t = self.EL_mbrT_tau / 365,
#                                     max_decay = self.N_max_decay)
#               if N_loss > N_removal:
#                   warn(f'Nitrogen degradation loss ({N_loss:.2f}) larger than the given removal ({N_removal:.2f})), '
#                         'the degradation loss is assumed to be the same as the removal.')
#                   N_loss = N_removal
            
#               # N2O only from the degraded part
#               N_loss_tot = N_loss * WasteWater.TN / 1e3 * WasteWater.F_vol
#               N2O_emission.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44 / 28
#         else:
#               N2O_emission.empty()
        
#         # NO3 conversion
#         NH3_mass = TreatedWater.imass['NH3']  # Inherite NH4 property from anoxic tank
#         NO3_mass_generated = NH3_mass * self.NO3_produced_ratio
#         TreatedWater.imass['NO3'] += NO3_mass_generated
#         TreatedWater.imass['NH3'] = 0
        
#         # N2O emission
#         N2O_mass_generated = NH3_mass * (1 - self.NO3_produced_ratio)
#         N2O_emission.imass['N2O'] += N2O_mass_generated
        
#         # Split NO3 into PC_nitrate_return and AnoT_nitrate_return
#         NO3_return_ratio = self.NO3_split_ratio  # Ratio of NO3 going to PC vs. AnoT
#         PC_nitrate_return.imass['NO3'] = TreatedWater.imass['NO3'] * NO3_return_ratio
#         AnoT_nitrate_return.imass['NO3'] = TreatedWater.imass['NO3'] * (1 - NO3_return_ratio)

#         # Remove NO3 from TreatedWater since it is returned to the previous tanks
#         TreatedWater.imass['NO3'] -= (PC_nitrate_return.imass['NO3'] + AnoT_nitrate_return.imass['NO3'])
        
#         # Update COD
#         TreatedWater._COD = WasteWater.COD * (1 - self.EL_mbrT_COD_removal)
#         sludge._COD = WasteWater.COD * self.EL_mbrT_COD_removal
    
#     # def _run(self):

#     #     # Input streams
#     #     WasteWater = self.ins[0]
#     #     air = self.ins[1]

#     #     # Output stream
#     #     TreatedWater = self.outs[0] 
#     #     PC_sludge_return = self.outs[1]
#     #     AnoT_sludge_return = self.outs[2]
#     #     CH4_emission = self.outs[3] 
#     #     N2O_emission = self.outs[4]
        
#     #     input_streams = [WasteWater, air]
            
#     #     # Mix all inputs into a single stream
#     #     self._mixed.empty()
#     #     self._mixed.mix_from(input_streams)

#     #     # Copy the mixed result to the outflow
#     #     TreatedWater.copy_like(self._mixed)
    
#     # def _run(self):
        
#     #     # Input streams
#     #     WasteWater = self.ins[0]
#     #     air = self.ins[1]

#     #     # Output streams
#     #     TreatedWater = self.outs[0]
#     #     PC_sludge_return = self.outs[1]
#     #     AnoT_sludge_return = self.outs[2]
#     #     CH4_emission = self.outs[3]
#     #     N2O_emission = self.outs[4]

#     #     # Air does not affect mass balance
#     #     # Air is used for aeration, so we do not modify WasteWater mass

#     #     # Step 1: Copy WasteWater to TreatedWater
#     #     TreatedWater.copy_like(WasteWater)

#     #     # Step 2: Split sludge for returns (1% to PC, 99% to AnoT)
#     #     total_f_sludge = 0.05  # Total sludge accounts for 5% of the wastewater mass
#     #     f_PC_sludge = 0.01     # Primary clarifier accounts for 1% of the total sludge
#     #     f_AnoT_sludge = 0.99   # Anoxic tank accounts for 99% of the total sludge

#     #     # Sludge mass
#     #     total_sludge_mass = WasteWater.F_mass * total_f_sludge
#     #     PC_sludge_return.mass = total_sludge_mass * f_PC_sludge
#     #     AnoT_sludge_return.mass = total_sludge_mass * f_AnoT_sludge

#     #     # Update TreatedWater after sludge removal
#     #     TreatedWater.F_mass -= total_sludge_mass
#     #     for component in TreatedWater.components:
#     #         TreatedWater.imass[component] -= (
#     #             PC_sludge_return.imass[component] + AnoT_sludge_return.imass[component]
#     #         )

#     #     # Step 3: Calculate emissions using Decay._first_order_run
#     #     super()._first_order_run(waste=WasteWater,
#     #                              treated=TreatedWater,
#     #                              biogas=None,  # Assume no biogas is captured
#     #                              CH4=CH4_emission,
#     #                              N2O=N2O_emission
#     #                              )

#     def _design(self):
#         design = self.design_results
#         constr = self.construction
#         design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # assume linear scaling
#         design['PVDF_membrane'] = constr[1].quantity = self.membrane_material_weight
#         self.add_construction(add_cost=False)

#     def _cost(self):
#         C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
#         C['MBR_tank'] = self.MBR_tank_cost
#         C['pipeline'] = self.pipeline_connectors
#         C['fittings'] = self.weld_female_adapter_fittings
#         C['Membrane_material'] = self.membrane_material_price * self.membrane_material_weight
#         C['Membrane_cleaning'] = self.membrane_cleaning_fee

#         ratio = self.price_ratio
#         for equipment, cost in C.items():
#             C[equipment] = cost * ratio 
        
#         self.add_OPEX = self._calc_replacement_cost()
        
#         power_demand = self.power_demand_MBR #TODO: power_demand unit should be kW, but the tsv file is in kWh/day
#         self.power_utility(power_demand)

#     def _calc_replacement_cost(self):
#         scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
#         MBR_replacement_cost = (
#             self.MBR_tank_cost / self.MBR_tank_lifetime +
#             self.pipeline_connectors / self.pipeline_connectors_lifetime +
#             self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
#             self.membrane_material_price * self.membrane_material_weight / self.membrane_material_lifetime
#             ) * scale
#         MBR_replacement_cost = MBR_replacement_cost / (365 * 24) * self.price_ratio # USD/hr
#         return MBR_replacement_cost

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
    _N_outs = 6 # [0] filtrate, [1] pumped flow
    _outs_size_is_fixed = True
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', isdynamic=True, 
                 pumped_flow=50, solids_capture_rate=0.999, 
                 V_max=1000, crossflow_air=None,
                 **kwargs):
        super().__init__(ID, ins, outs, 
                         # split=None, 
                         thermo=thermo,
                         init_with=init_with, V_max=V_max, isdynamic=isdynamic, 
                         **kwargs)
        self.pumped_flow = pumped_flow
        self.solids_capture_rate = solids_capture_rate
        self.crossflow_air = crossflow_air
    
    @property
    def pumped_flow(self):
        '''[float] Pumped flow rate, in m3/d'''
        return self._Q_pump
    @pumped_flow.setter
    def pumped_flow(self, Q):
        self._Q_pump = Q
    
    @property
    def solids_capture_rate(self):
        '''[float] Membrane solid capture rate, i.e., 
        filtrate-to-internal solids concentration ratio, unitless.'''
        return self._f_rtn
    @solids_capture_rate.setter
    def solids_capture_rate(self, f):
        if f < 0 or f > 1:
            raise ValueError(f'membrane solids capture rate must be within [0,1], not {f}')
        self._f_rtn = f
        cmps = self._mixed.components
        self._flt2in_conc_ratio = (1-cmps.x) + cmps.x * (1-f)
    
    @property
    def crossflow_air(self):
        '''[:class:`qsdsan.Process` or NoneType]
        Membrane cross flow air specified for process modeling, such as `qsdsan.processes.DiffusedAeration`. 
        Ignored if DO setpoint is specified by the `aeration` attribute.
        '''
        return self._cfa
    @crossflow_air.setter
    def crossflow_air(self, cfa):
        if cfa is None or isinstance(cfa, Process):
            self._cfa = cfa
        else:
            raise TypeError('crossflow_air must be a `Process` object or None, '
                            f'not {type(cfa)}')

    split = None
        
    def _run(self):
        '''Only to converge volumetric flows.'''
        mixed = self._mixed
        mixed.mix_from(self.ins)
        cmps = mixed.components
        Q = mixed.F_vol*24 # m3/d
        Qp = self._Q_pump
        f_rtn = self._f_rtn
        xsplit = Qp / ((1-f_rtn)*(Q-Qp) + Qp) # mass split of solids to pumped flow
        qsplit = Qp / Q
        flt, rtn = self.outs
        mixed.split_to(rtn, flt, xsplit*cmps.x + qsplit*(1-cmps.x))
    
    def _compile_ODE(self):
        aer = self._aeration
        cfa = self._cfa
        isa = isinstance
        cmps = self.components
        if self._model is None:
            warn(f'{self.ID} was initialized without a suspended growth model, '
                 f'and thus run as a non-reactive unit')
            r = lambda state_arr: np.zeros(cmps.size)
        else:
            r = self._model.production_rates_eval

        _dstate = self._dstate
        _update_dstate = self._update_dstate
        V = self._V_max
        Qp = self.pumped_flow
        f_rtn = self.solids_capture_rate
        xarr = cmps.x
        gstrip = self.gas_stripping
        if gstrip:
            gas_idx = self.components.indices(self.gas_IDs)
            if isa(aer, Process): kLa = aer.kLa
            else: kLa = 0.
            if cfa: kLa += cfa.kLa
            S_gas_air = np.asarray(self.K_Henry)*np.asarray(self.p_gas_atm)
            kLa_stripping = np.maximum(kLa*self.D_gas/self._D_O2, self.stripping_kLa_min)
        hasexo = bool(len(self._exovars))
        f_exovars = self.eval_exo_dynamic_vars
         
        if isa(aer, (float, int)):
            i = cmps.index(self._DO_ID)
            def dy_dt(t, QC_ins, QC, dQC_ins):
                QC[i] = aer
                dydt_mbr(QC_ins, QC, V, Qp, f_rtn, xarr, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC)
                if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
                _dstate[i] = 0
                _update_dstate()
        else:        
            if cfa:
                cfa_stoi = cfa._stoichiometry
                cfa_frho = cfa.rate_function
                dy_cfa = lambda QC: cfa_stoi * cfa_frho(QC)
            else:
                dy_cfa = lambda QC: 0.
            
            if isa(aer, Process):
                aer_stoi = aer._stoichiometry
                aer_frho = aer.rate_function
                dy_aer = lambda QC: aer_stoi * aer_frho(QC)
            else:
                dy_aer = lambda QC: 0.
                
            def dy_dt(t, QC_ins, QC, dQC_ins):
                dydt_mbr(QC_ins, QC, V, Qp, f_rtn, xarr, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC) + dy_aer(QC) + dy_cfa(QC)
                if gstrip: _dstate[gas_idx] -= kLa_stripping * (QC[gas_idx] - S_gas_air)
                _update_dstate()
        self._ODE = dy_dt
    
    def _update_state(self):
        arr = self._state
        arr[arr < 1e-16] = 0.
        arr[-1] = sum(ws.state[-1] for ws in self.ins)
        for ws in self.outs:
            if ws.state is None: 
                ws.state = np.zeros_like(arr)
                ws.dstate = np.zeros_like(arr)
        flt, rtn = self.outs
        Qp = self.pumped_flow
        flt.state[:-1] = arr[:-1] * self._flt2in_conc_ratio
        flt.state[-1] = arr[-1] - Qp
        rtn.state[:-1] = arr[:-1]
        rtn.state[-1] = Qp
        
    def _update_dstate(self):
        arr = self._dstate
        arr[-1] = sum(ws.dstate[-1] for ws in self.ins)
        flt, rtn = self.outs
        flt.dstate[:-1] = arr[:-1] * self._flt2in_conc_ratio
        flt.dstate[-1] = arr[-1]
        rtn.dstate[:-1] = arr[:-1]
        rtn.dstate[-1] = 0

    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),
                             Construction(item ='PVDF_membrane', linked_unit=self, quantity_unit='kg'),]
    

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # assume linear scaling
        design['PVDF_membrane'] = constr[1].quantity = self.membrane_material_weight
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['MBR_tank'] = self.MBR_tank_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings
        C['Membrane_material'] = self.membrane_material_price * self.membrane_material_weight
        C['Membrane_cleaning'] = self.membrane_cleaning_fee

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio 
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_MBR
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        MBR_replacement_cost = (
            self.MBR_tank_cost / self.MBR_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
            self.membrane_material_price * self.membrane_material_weight / self.membrane_material_lifetime
            ) * scale
        MBR_replacement_cost = MBR_replacement_cost / (365 * 24) * self.price_ratio # USD/hr
        return MBR_replacement_cost

# %%
# ClearWaterTank_path = ospath.join(EL_su_data_path, '_EL_CWT.tsv')

# @price_ratio()
# class EL_CWT(StorageTank):

#     '''
#     Introduction
#     ------------
#     To only collect the treated water

#     Parameters
#     ----------
#     Ins:
#     (1) influent of treated wastetwater from self-priming pump
#     (2) O3 dosing
#     (3) air-dissolve influent

#     Outs:
#     (1) effluent of treated wastewater (clear water)
#     (2) spill flow to collection tank


#     Attributes
#     ----------

    
#     References
#     ----------
#      refer to the qsdsan.sanunits.storagetank module

#     '''
#     _N_ins = 3  # number of ins
#     _N_outs = 2  # number of outs
#     _ins_size_is_fixed = True
#     _outs_size_is_fixed = True
#     exponent_scale = 0.1

#     def __init__(self, ID = '', ins = None, outs = (), thermo = None, ppl = None, baseline_ppl = None,
#                  vessel_type='Field erected',tau = 24, V_wf = None, vessel_material='Stainless steel', kW_per_m3 = 0.1,
#                  init_with = 'WasteStream', F_BM_default = 1, max_overflow=15, include_construction = True, **kwargs):
#         StorageTank.__init__(self, ID=ID, ins=ins, outs=outs, thermo = thermo, init_with = init_with, 
#                              F_BM_default = F_BM_default, include_construction = include_construction,
#                              kW_per_m3= kW_per_m3, vessel_type= vessel_type, tau= tau,
#                              vessel_material= vessel_material, V_wf= V_wf)

#         self.ppl = ppl
#         self.baseline_ppl = baseline_ppl
#         self.max_overflow = max_overflow

#         data = load_data(path = ClearWaterTank_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             setattr(self, para, value)
#         del data

#         for attr, value in kwargs.items():
#             setattr(self, attr, value)

#     def _init_lca(self):
#         self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]
    
#     def _run(self):
        
#         # Input streams
#         WasteWater = self.ins[0]
#         ozone = self.ins[1]
#         air = self.ins[2]
        
#         # Output streams
#         TreatedWater = self.outs[0]
#         spill_return = self.outs[1]
#         # CWT_cycle = self.outs[2]
        
#         # Inherited input stream
#         TreatedWater.copy_like(WasteWater)
        
#         # Spill water return
#         max_overflow_m3d = self.max_overflow * 24  # Convert m³/h to m³/d
#         Q_treated = TreatedWater.F_vol * 24  # m³/d
#         if Q_treated > max_overflow_m3d:
#             spill_vol = Q_treated - max_overflow_m3d
#             f_spill = spill_vol / Q_treated
#             spill_return.copy_like(TreatedWater)
#             spill_return.F_mass *= f_spill
#             TreatedWater.F_mass *= (1 - f_spill)
#         else:
#             # max_overflow is none, no spill return
#             spill_return.empty()
        

#         # # Input streams
#         # WasteWater = self.ins[0]
#         # ozone = self.ins[1]
#         # air = self.ins[2]

#         # # Output streams
#         # TreatedWater = self.outs[0]
#         # spill_return = self.outs[1]
#         # CWT_cycle = self.outs[2]
        
#         # input_streams = [WasteWater, ozone, air]
            
#         # # Mix all inputs into a single stream
#         # self._mixed.empty()
#         # self._mixed.mix_from(input_streams)
        

#         # # Copy the mixed result to the outflow
#         # TreatedWater.copy_like(self._mixed)
        
    
    
#     # def _run(self):
        
#     #     # Input streams
#     #     WasteWater = self.ins[0]
#     #     ozone = self.ins[1]
#     #     air = self.ins[2]

#     #     # Output streams
#     #     TreatedWater = self.outs[0]
#     #     CT_spill_return = self.outs[1]
#     #     PC_spill_return = self.outs[2]
 
#     #     # Constants for overflow conditions
#     #     CT_spill_fraction = 0.6
#     #     PC_spill_fraction = 0.4

#     #     # Step 1: Process inputs
#     #     # Air and ozone are balanced, not affecting mass balance
#     #     TreatedWater.copy_like(WasteWater)

#     #     # Step 2: Handle overflow if necessary
#     #     if TreatedWater.F_vol > self.max_overflow:
#     #         # Calculate excess volume
#     #         spill_vol = TreatedWater.F_vol - self.max_overflow

#     #         # Calculate spill fractions
#     #         self._f_spill = spill_vol / TreatedWater.F_vol
#     #         self._f_treated = 1 - self._f_spill

#     #         # Assign spill fractions to CT_spill_return and PC_spill_return
#     #         CT_spill_return.F_mass[:] = TreatedWater.F_mass[:] * self._f_spill * CT_spill_fraction
#     #         PC_spill_return.F_mass[:] = TreatedWater.F_mass[:] * self._f_spill * PC_spill_fraction

#     #         # Adjust the normal treated water to within max capacity
#     #         TreatedWater.F_mass[:] *= self._f_treated
#     #     else:
#     #         # No spill return if overflow is within capacity
#     #         CT_spill_return.empty()
#     #         PC_spill_return.empty()
#     #         self._f_spill = 0.0
#     #         self._f_treated = 1.0

#     def _design(self):
#         design = self.design_results
#         constr = self.construction
#         design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # to be defined in .tsv file
#         self.add_construction(add_cost=False)

#     def _cost(self):
#         C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
#         C['Clear water tank'] = self.clear_water_tank_cost
#         C['pipeline'] = self.pipeline_connectors
#         C['fittings'] = self.weld_female_adapter_fittings
#         C['O3 generator'] = self.O3_generation_machine_cost

#         ratio = self.price_ratio
#         for equipment, cost in C.items():
#             C[equipment] = cost * ratio

#         self.add_OPEX = self._calc_replacement_cost()
        
#         power_demand = self.power_demand_CWT
#         self.power_utility(power_demand)

#     def _calc_replacement_cost(self):
#         scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
#         CWR_replacement_cost = (
#             self.clear_water_tank_cost / self.clear_water_tank_lifetime +
#             self.pipeline_connectors / self.pipeline_connectors_lifetime +
#             self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
#             self.O3_generation_machine_cost / self.O3_generation_machine_lifetime
#             ) * scale
#         CWR_replacement_cost = CWR_replacement_cost / (365 * 24) # convert to USD/hr
#         return CWR_replacement_cost

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
    '''
    Similar to :class:`biosteam.units.Mixer`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`,
    and allows dynamic simulation.

    See Also
    --------
    `biosteam.units.Mixer <https://biosteam.readthedocs.io/en/latest/units/mixing.html>`_
    '''
    _N_ins = 3
    _N_outs = 2
    _graphics = BSTMixer._graphics
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 rigorous=False, conserve_phases=False):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        self.rigorous = rigorous
        self.conserve_phases = conserve_phases


    @property
    def state(self):
        '''The state of the Mixer, including component concentrations [mg/L] and flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    def _init_state(self):
        '''initialize state by specifying or calculating component concentrations
        based on influents. Total flow rate is always initialized as the sum of
        influent wastestream flows.'''
        QCs = self._ins_QC
        if QCs.shape[0] <= 1: self._state = QCs[0]
        else:
            Qs = QCs[:,-1]
            Cs = QCs[:,:-1]
            self._state = np.append(Qs @ Cs / Qs.sum(), Qs.sum())
        self._dstate = self._state * 0.

    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Mixer'''
        self._outs[0].state = self._state

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Mixer'''
        self._outs[0].dstate = self._dstate

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _n_ins = len(self.ins)
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        def yt(t, QC_ins, dQC_ins):
            if _n_ins > 1:
                Q_ins = QC_ins[:, -1]
                C_ins = QC_ins[:, :-1]
                dQ_ins = dQC_ins[:, -1]
                dC_ins = dQC_ins[:, :-1]
                Q = Q_ins.sum()
                C = Q_ins @ C_ins / Q
                _state[-1] = Q
                _state[:-1] = C
                Q_dot = dQ_ins.sum()
                C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
                _dstate[-1] = Q_dot
                _dstate[:-1] = C_dot
            else:
                _state[:] = QC_ins[0]
                _dstate[:] = dQC_ins[0]
            _update_state()
            _update_dstate()
        self._AE = yt


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # to be defined in .tsv file
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['Clear water tank'] = self.clear_water_tank_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings
        C['O3 generator'] = self.O3_generation_machine_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_CWT
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        CWR_replacement_cost = (
            self.clear_water_tank_cost / self.clear_water_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
            self.O3_generation_machine_cost / self.O3_generation_machine_lifetime
            ) * scale
        CWR_replacement_cost = CWR_replacement_cost / (365 * 24) # convert to USD/hr
        return CWR_replacement_cost

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
        StorageTank.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo,
                      init_with=init_with, F_BM_default=F_BM_default,
                      include_construction=include_construction,
                      vessel_type=vessel_type, tau=tau, V_wf=V_wf,
                      vessel_material=vessel_material, kW_per_m3=kW_per_m3,)
        self.length_to_diameter = length_to_diameter
        
    def _init_lca(self):
        item_name = self.vessel_material.replace(' ', '_')
        self.construction = [
            Construction(item_name.lower(), linked_unit=self, item=item_name, quantity_unit='kg'),
            ]
        
    
    def _design(self):
        BSTStorageTank._design(self)
        D = self.design_results
        
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
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # to be defined in .tsv file
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
