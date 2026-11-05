# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

Silicon Valley Clean Water (SVCW) treatment train based on WERF F1, with a
trickling filter before activated sludge and a final sand filter.

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from math import isfinite

import numpy as np

import qsdsan as qs
from qsdsan import process_models as pc, unit_operations as su
from qsdsan.utils import load_data, ospath

from exposan.werf import data_path, default_fctss_init
from exposan.werf._units import TricklingFilter
from exposan.werf.systems.D1 import _tf_tolerance_solver

from .TricklingFilter_validation import (
    _create_validation_components,
    _create_validation_model,
    _set_tf_init_conc,
    )

__all__ = ('create_svcw_system',)

ID = 'SVCW'


# %% Source data and base constants
_dfs = load_data(
    ospath.join(data_path, 'initial_conditions.xlsx'),
    sheet=None,
    )
_asinit = _dfs['E2P']
_fcinit = _asinit.iloc[-1].to_dict()
_adinit = _dfs['adm'].loc['F1'].to_dict()

MGD2cmd = 3785.412
T_ad = 273.15 + 35


# %% SVCW calibration
# Flow, COD, and NH4-N are June 2026 SVCW means. PO4-P is a calibration
# assumption because the plant export has no contemporaneous orthophosphate
# measurement. Influent ions are inferred from measured digestate samples.
_SVCW_SUMMER_2026_FLOW_MGD = 13.685
_SVCW_SUMMER_2026_WASTEWATER = dict(
    COD=600, NH4_N=44, PO4_P=2.0,
    S_Na=278.5, S_K=66.5, S_Mg=70, S_Ca=650, S_Cl=643.9,
    )

# Moving COD from X_I to S_I matches June TSS without changing BOD or total N.
_SVCW_SUMMER_2026_ORGANIC_FRACTIONS = dict(
    fr_SI=0.14, fr_SF=0.24, fr_SA=0.05, fr_XI=0.235,
    )

_SVCW_PARTICULATE_N_FACTOR = 1

# The TF uses base mASM2d decay/kinetics with the calibrated rates below. These
# settings are intentionally independent from the activated-sludge model.
_SVCW_TF_MODEL_KWARGS = dict(
    gpsx_decay=False,
    gpsx_kinetics=False,
    gpsx_aob_proxy=False,
    anaerobic_decay_factor=0.6,
    )
# _SVCW_TF_ASM_KWARGS = dict(q_PHA=6.0, q_PP=2.0, b_PP=0.05, mu_H=3.2) #GPS-X values
_SVCW_TF_ASM_KWARGS = dict(q_PHA=3.0,q_PP=1.5, b_PP=0.2,mu_H=6.0) #QSDsan default values
_SVCW_AD_LIQUID_VOLUME_MGD = 0.85*1.3685
_SVCW_AD_INITIAL_S_IP = 24.5
_SVCW_AD_INITIAL_S_MG = 10.0
_SVCW_AD_INITIAL_S_CA = 80.0
_SVCW_AD_PH = 7.81

# Apparent equilibrium values are calibrated against the measured AD sample at
# pH 7.81 and 35 deg C; CaCO3 and MgCO3 also assume 270 mM inorganic carbon.
_SVCW_AD_CACO3_PKSP = 5.5 #default value 8.5
_SVCW_AD_STRUVITE_PKSP = 11.6 #default value 13.7, apparent value 11.6
_SVCW_AD_ACP_PKSP = 23.3 #default value 28.6, apparent value 23.3
_SVCW_AD_MGCO3_PKSP = 5.6 #default value 7.6

# ADM1p order: CaCO3, struvite, newberyite, ACP, MgCO3, AlPO4, FePO4.
_SVCW_AD_K_MMP = (0.24, 240.0, 0.24, 72.0, 0.24, 0.24, 0.24)


# %% Initialization safeguard: mineral saturation projection
def _equilibrate_ad_mineral_init(
        model, concentrations, pH, T, rtol=1e-8, maxiter=25,
        saturation_target=0.99):
    '''Return precipitation-equilibrated ADM1p mineral initial conditions.

    The projection follows the KM mineral model used during simulation. It
    sequentially advances the currently fastest seeded mineral to saturation,
    conserves every component through ADM1p's own stoichiometric matrix, and
    leaves undersaturated existing minerals intact because the KM rate function
    does not model their dissolution.

    ``saturation_target`` defaults to 0.99 so the initial state remains just
    below ADM1p's piecewise precipitation threshold at SI = 1. This small
    numerical buffer avoids Jacobian perturbations repeatedly toggling mineral
    rates without materially changing the projected mineral inventory.
    A target slightly above one may instead initialize the active precipitation
    branch of a continuously fed reactor whose influent supplies mineral ions.
    '''
    pH = float(pH)
    T = float(T)
    if not isfinite(pH) or not 0. < pH < 14.:
        raise ValueError('ad_pH must be finite and between 0 and 14')
    if not isfinite(T) or T <= 0.:
        raise ValueError('AD temperature must be finite and positive')
    saturation_target = float(saturation_target)
    if (not isfinite(saturation_target)
            or not 0. < saturation_target <= 1.1):
        raise ValueError('saturation_target must be in the interval (0, 1.1]')

    adjusted = dict(concentrations)
    cmps = model._components
    values = np.array(
        [adjusted.get(ID, 0.) for ID in cmps.IDs], dtype=float,
        ) * 1e-3  # ADM initial concentrations are mg/L; rates use kg/m3.
    params = model.rate_function.params
    Ka = params['Ka_base'] * pc.T_correction_factor(
        params['T_base'], T, params['Ka_dH'],
        )
    Knh, Kc1, Kc2, Kp1, Kp2, Kp3 = Ka[1:7]
    h = 10.**(-pH)
    nh4_fraction = h / (Knh + h)
    co3_fraction = pc.ion_speciation(h, Kc1, Kc2)[2]
    phosphate_fractions = pc.ion_speciation(h, Kp1, Kp2, Kp3)
    hpo4_fraction, po4_fraction = phosphate_fractions[2:]
    Ksp = np.asarray(params['Ksp'], dtype=float)

    component_indices = {
        ID: cmps.index(ID) for ID in (
            'S_IC', 'S_IN', 'S_IP', 'S_Mg', 'S_Ca',
            'X_CaCO3', 'X_struv', 'X_newb', 'X_ACP', 'X_MgCO3',
            )
        }
    process_IDs = (
        'CaCO3_precipitation_dissolution',
        'struvite_precipitation_dissolution',
        'newberyite_precipitation_dissolution',
        'ACP_precipitation_dissolution',
        'MgCO3_precipitation_dissolution',
        )
    mineral_IDs = ('X_CaCO3', 'X_struv', 'X_newb', 'X_ACP', 'X_MgCO3')
    mineral_indices = np.array([component_indices[ID] for ID in mineral_IDs])
    stoichiometry = model.stoichio_eval()
    mineral_stoichiometry = np.array(
        [stoichiometry[model.index(ID)] for ID in process_IDs], dtype=float,
        )
    rate_constants = np.asarray(params['k_mmp'][:5], dtype=float)
    rate_powers = np.array((2., 3., 2., 2., 2.))

    def saturation_indices(state):
        S_IC = max(0., state[component_indices['S_IC']])
        S_IN = max(0., state[component_indices['S_IN']])
        S_IP = max(0., state[component_indices['S_IP']])
        S_Mg = max(0., state[component_indices['S_Mg']])
        S_Ca = max(0., state[component_indices['S_Ca']])
        nh4 = nh4_fraction * S_IN
        co3 = co3_fraction * S_IC
        hpo4 = hpo4_fraction * S_IP
        po4 = po4_fraction * S_IP
        return np.array((
            np.sqrt(S_Ca * co3 / Ksp[0]) if S_Ca and co3 else 0.,
            np.cbrt(S_Mg * nh4 * po4 / Ksp[1])
                if S_Mg and nh4 and po4 else 0.,
            np.sqrt(S_Mg * hpo4 / Ksp[2]) if S_Mg and hpo4 else 0.,
            (S_Ca**3 * po4**2 / Ksp[3])**0.2
                if S_Ca and po4 else 0.,
            np.sqrt(S_Mg * co3 / Ksp[4]) if S_Mg and co3 else 0.,
            ))

    for _ in range(maxiter):
        saturation = saturation_indices(values)
        active = ((saturation > saturation_target + rtol)
                  & (values[mineral_indices] > 0.))
        if not active.any():
            break
        rates = (rate_constants * values[mineral_indices]
                 * np.maximum(saturation - 1., 0.)**rate_powers)
        mineral = int(np.argmax(np.where(active, rates, -np.inf)))
        nu = mineral_stoichiometry[mineral]
        reactants = np.flatnonzero(nu < 0.)
        upper = min(values[i] / -nu[i] for i in reactants)
        lower = 0.
        for _ in range(80):
            midpoint = 0.5 * (lower + upper)
            trial = values + nu * midpoint
            if saturation_indices(trial)[mineral] > saturation_target:
                lower = midpoint
            else:
                upper = midpoint
            if upper - lower <= max(1e-15, rtol * upper):
                break
        values += nu * upper
        values[abs(values) < 1e-15] = 0.

    for ID in component_indices:
        adjusted[ID] = max(0., values[component_indices[ID]] * 1e3)
    return adjusted


# %% Runtime safeguards: kinetic trial states and solver bookkeeping
class _SVCWAnaerobicCSTR(su.AnaerobicCSTR):
    '''SVCW AD with a nonnegative projection for kinetic evaluation.

    QSDsan's base implementation clips the shared system state in
    ``_update_state``. This subclass retains the clipped values only for the
    physical outlet streams, restores the raw solver state before evaluating
    the ODE, and projects a copy solely for ADM1p kinetics. The hydraulic term
    can therefore drive a negative solver trial back toward zero without
    evaluating concentration-dependent rates at an impossible state.
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._raw_flex_rhos = self.model.flex_rhos
        self._raw_rate_function = self.model.rate_function.function
        self.model.__dict__['flex_rhos'] = self._projected_flex_rhos
        self.model.rate_function._function = self._projected_rate_function

    def _project_kinetic_state(self, state):
        projected = np.array(state, dtype=float, copy=True)
        stop = len(self.components) + self._n_gas
        projected[:stop] = np.maximum(projected[:stop], 0.)
        return projected

    def _projected_flex_rhos(self, state, params, h=None):
        return self._raw_flex_rhos(
            self._project_kinetic_state(state), params, h=h)

    def _projected_rate_function(self, state, params):
        return self._raw_rate_function(
            self._project_kinetic_state(state), params)

    def _update_state(self):
        raw = self._state.copy()
        super()._update_state()
        # Retain the updated algebraic flow, but do not let outlet formatting
        # silently overwrite the concentrations seen by the solver and ODE.
        self._state[:-1] = raw[:-1]


class _SVCWSystem(qs.System):
    '''SVCW system with TF-specific solver tolerances.'''

    def dynamic_run(self, **kwargs):
        '''Run the requested time span with targeted trickling-filter tolerances.'''
        method = kwargs.get('method', self.dynsim_kwargs.get('method', 'BDF'))
        kwargs['method'] = _tf_tolerance_solver(method, self)
        try:
            result = super().dynamic_run(**kwargs)
            self._restore_initial_solution()
            return result
        finally:
            self.dynsim_kwargs['method'] = method

    def _restore_initial_solution(self):
        '''Restore the true initial sample after BioSTEAM's aliased t0 state.'''
        solution = getattr(self.scope, 'sol', None)
        if solution is None or not getattr(solution, 'y', np.empty(0)).size:
            return
        initial = self.dynsim_kwargs.get('y0')
        if initial is not None and solution.y.shape[0] == len(initial):
            # BioSTEAM's result list aliases System._state at t0; later RHS
            # evaluations otherwise replace the first saved state in place.
            solution.y[:, 0] = initial


# %% Component and interface setup
def _create_svcw_components(particulate_n_factor):
    '''Create SVCW components with adjustable particulate N-to-COD ratios.'''
    factor = float(particulate_n_factor)
    if not isfinite(factor) or factor <= 0.:
        raise ValueError('particulate_n_factor must be finite and positive')
    cmps = _create_validation_components('mid')
    for component_ID in ('X_S', 'X_I'):
        component = getattr(cmps, component_ID)
        component.i_N *= factor
    cmps.refresh_constants()
    return cmps


def _synchronize_adm1p_components(cmps_asm, cmps_adm):
    '''Copy shared ASM component indices before ADM1p is compiled.'''
    specifications = {
        'X_PHA': ('measured_as', 'i_COD', 'i_C', 'i_N', 'i_P'),
        'X_PP': ('measured_as', 'i_COD', 'i_C', 'i_N', 'i_P'),
        'X_PAO': ('measured_as', 'i_COD', 'i_C', 'i_N', 'i_P'),
        'S_I': ('measured_as', 'i_C', 'i_N', 'i_P'),
        'X_I': ('measured_as', 'i_C', 'i_N', 'i_P'),
        }
    changes = {}
    for component_id, attributes in specifications.items():
        source = getattr(cmps_asm, component_id)
        target = getattr(cmps_adm, component_id)
        component_changes = {}
        for attribute in attributes:
            old = getattr(target, attribute)
            new = getattr(source, attribute)
            if old != new:
                setattr(target, attribute, new)
                component_changes[attribute] = (old, new)
        if component_changes:
            changes[component_id] = component_changes
    cmps_adm.refresh_constants()
    return changes


def _create_svcw_wastewater(
        flow_mgd=_SVCW_SUMMER_2026_FLOW_MGD,
        composition=None, organic_fractions=None):
    '''Return the adjustable summer-2026 SVCW wastewater calibration profile.'''
    flow_mgd = float(flow_mgd)
    if not isfinite(flow_mgd) or flow_mgd <= 0.:
        raise ValueError('wastewater_flow_mgd must be finite and positive')
    calibrated_composition = dict(_SVCW_SUMMER_2026_WASTEWATER)
    calibrated_composition.update(composition or ())
    calibrated_fractions = dict(_SVCW_SUMMER_2026_ORGANIC_FRACTIONS)
    calibrated_fractions.update(organic_fractions or ())
    return pc.create_masm2d_inf(
        'RWW', flow_mgd, 'MGD', T=293.15,
        **calibrated_composition, **calibrated_fractions)


# %% SVCW construction
def create_svcw_system(
        flowsheet=None, default_init_conds=True,
        tf_kwargs=None, sf_kwargs=None,
        wastewater_flow_mgd=_SVCW_SUMMER_2026_FLOW_MGD,
        wastewater_kwargs=None, organic_fractions=None,
        particulate_n_factor=_SVCW_PARTICULATE_N_FACTOR,
        ad_liquid_volume_mgd=_SVCW_AD_LIQUID_VOLUME_MGD,
        ad_initial_s_ip=_SVCW_AD_INITIAL_S_IP,
        ad_pH=_SVCW_AD_PH,
        use_adm1p_default_pKsp=False,
        ad_caco3_pKsp=_SVCW_AD_CACO3_PKSP,
        ad_struvite_pKsp=_SVCW_AD_STRUVITE_PKSP,
        ad_acp_pKsp=_SVCW_AD_ACP_PKSP,
        ad_mgco3_pKsp=_SVCW_AD_MGCO3_PKSP,
        equilibrate_ad_minerals=True,
        ):
    '''Create the Silicon Valley Clean Water treatment system.

    The process is WERF F1 with ``PC -> TF -> ASR -> FC -> SF`` in the liquid
    train. Primary sludge and WAS retain F1's separate thickening, anaerobic
    digestion, and dewatering train. DW separates ADM1p digestate directly;
    only its liquid outlet passes through J2 to mASM2d before joining M2.
    Cake retains ADM1p biomass and mineral components without J2's
    instantaneous biomass conversion. Sand-filter backwash joins gravity-
    thickener overflow, mechanical-thickener overflow, and dewatering centrate
    in M2 before returning to PC through HD.

    Parameters
    ----------
    flowsheet : qsdsan.Flowsheet, optional
        Flowsheet to populate. A new SVCW flowsheet is created by default.
    default_init_conds : bool, optional
        Apply F1 initial conditions to ASR, FC, and AD and the validated
        attached-growth initialization to TF. AD dissolved Mg and Ca are
        fixed at 10 and 80 mg/L after optional mineral projection, independent
        of the selected pKsp profile.
    tf_kwargs : dict, optional
        Overrides for the trickling filter. Defaults are six horizontal
        sections, five biofilm layers, 9464 m2 bed surface, recirculation ratio
        one, dynamic biofilm thickness, and 4e-4 m initial total thickness.
    sf_kwargs : dict, optional
        Overrides for the constant-splitting sand filter. Defaults are 0.60
        particulate mass removal to backwash and 0.02 liquid backwash fraction.
        The model cannot be changed because empirical fouling has a different
        stream topology.
    wastewater_flow_mgd : float, optional
        Raw wastewater flow in MGD. The default is the June 2026 mean of
        13.685 MGD from the SVCW plant export.
    wastewater_kwargs : dict, optional
        Overrides for ``create_masm2d_inf`` composition arguments. Default COD
        and NH4-N are June 2026 means. Default Na, K, Mg, Ca, and Cl are
        digestate-targeted calibration estimates, not measured influent data.
    organic_fractions : dict, optional
        Overrides for influent COD fractions. The calibrated default transfers
        0.065 of total COD from particulate X_I to soluble S_I so predicted
        influent TSS approaches the June 2026 mean while BOD is unchanged.
    particulate_n_factor : float, optional
        Multiplier applied to the Table-S1 ``i_N`` values of X_S and X_I. The
        calibrated default is 0.9.
    ad_liquid_volume_mgd : float, optional
        Anaerobic-digester liquid volume in million gallons. The F1 value of
        0.85 is retained by default; gas volume is 10 percent of liquid volume.
    ad_initial_s_ip : float, optional
        Initial ADM1p dissolved inorganic phosphorus concentration [mg P/L]
        before optional mineral-equilibrium projection. The SVCW default is
        24.5 mg P/L.
    ad_pH : float or None, optional
        Fixed anaerobic-digester pH. The default is 7.81. Use ``None`` to let
        ADM1p calculate pH from its charge balance.
    use_adm1p_default_pKsp : bool, optional
        If true, use the ADM1p KM defaults for CaCO3, struvite, ACP, and MgCO3
        (8.5, 13.7, 28.6, and 7.6) instead of the SVCW values. This switch takes
        precedence over the four individual ``ad_*_pKsp`` arguments.
    ad_caco3_pKsp : float, optional
        Apparent molar calcium-carbonate pKsp used only by this SVCW system.
        The calibrated default is 5.8.
    ad_struvite_pKsp : float, optional
        Apparent molar struvite pKsp used only by this SVCW system. The
        calibrated default is 13.0.
    ad_acp_pKsp : float, optional
        Apparent molar amorphous-calcium-phosphate pKsp used only by this SVCW
        system. The calibrated default is 25.8.
    ad_mgco3_pKsp : float, optional
        Apparent molar MgCO3 pKsp. The calibrated default is 5.9.
    equilibrate_ad_minerals : bool, optional
        If true, project the F1 AD inorganic C/N/P and mineral pools along
        precipitation reaction rows at ``ad_pH`` before imposing the fixed
        initial dissolved Mg and Ca concentrations. The final startup state
        can therefore be supersaturated. Existing undersaturated minerals are
        not dissolved; this option is ignored when ``ad_pH`` is ``None``.

    Returns
    -------
    qsdsan.System
        Dynamic SVCW system. ``FC_SE`` is the unfiltered secondary effluent;
        ``SE`` is the final sand-filter effluent.
    '''
    flowsheet = flowsheet or qs.Flowsheet(ID)
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and process models
    cmps_asm = _create_svcw_components(particulate_n_factor)
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)
    tf_asm = _create_validation_model(
        **_SVCW_TF_MODEL_KWARGS,
        asm_kwargs=dict(_SVCW_TF_ASM_KWARGS),
        )
    rww = _create_svcw_wastewater(
        wastewater_flow_mgd, wastewater_kwargs, organic_fractions,
        )
    thermo_asm = qs.get_thermo()

    # Liquid treatment train
    PC = su.PrimaryClarifier(
        'PC', ins=[rww, 'reject'], outs=('PE', 'PS'),
        isdynamic=True,
        sludge_flow_rate=0.074*MGD2cmd*1.3685*0.8, #scale based on 13.685/10, scale down to match S_IN in AD with actual sample
        solids_removal_efficiency=0.6,
        )

    GT = su.IdealClarifier(
        'GT', PC-1, outs=('', 'thickened_PS'),
        sludge_flow_rate=0.026*MGD2cmd*1.3685*0.8, #scale based on 13.685/10, scale down to match S_IN in AD with actual sample
        solids_removal_efficiency=0.9,
        )

    tf_kwargs = dict(tf_kwargs or ())
    n_sections = int(tf_kwargs.get('N_section', 6))
    tf_defaults = dict(
        N_section=n_sections,
        N_biofilm_layer=4,
        influent_fractions=[[1] + [0]*(n_sections-1)],
        recirculation_ratio=1.,
        suspended_growth_model=tf_asm,
        DO_ID='S_O2',
        max_liquid_film_thickness=5e-5,
        max_biofilm_thickness=1e-3,
        filter_bed_depth=2.,
        filter_bed_surface=9464.,
        media_specific_surface=100.,
        dynamic_biofilm_thickness=True,
        initial_biofilm_thickness=4e-4,
        detachment_model='capacity',
        detachment_exponent=1.,
        detachment_rate=0.065,
        )
    tf_defaults.update(tf_kwargs)
    TF = TricklingFilter('TF', ins=PC-0, outs='TF_eff', **tf_defaults)

    n_zones = 6
    V_tot = 2.4*MGD2cmd
    ASR = su.PFR(
        'ASR', ins=[TF-0, 'RAS'], outs='treated',
        N_tanks_in_series=n_zones,
        V_tanks=[V_tot/n_zones]*n_zones,
        influent_fractions=[[1] + [0]*5]*2,
        internal_recycles=[], DO_ID='S_O2',
        kLa=[50]*n_zones,
        DO_setpoints=[2.5]*n_zones,
        suspended_growth_model=asm,
        gas_stripping=True,
        )

    FC = su.FlatBottomCircularClarifier(
        'FC', ins=ASR-0, outs=['FC_SE', 1-ASR, 'WAS'],
        underflow=0.67*10*MGD2cmd, wastage=0.163*MGD2cmd*1.3685,
        surface_area=1579.352, height=3.6576,
        N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=410, v_max_practical=274,
        rh=4e-4, rp=0.1, fns=0.01,
        maximum_nonsettleable_solids=8.0,
        )

    sf_kwargs = dict(sf_kwargs or ())
    sf_defaults = dict(
        model='constant_splitting',
        backwash_solids_mass_fraction=0.6,
        backwash_flow_fraction=0.02,
        )
    sf_defaults.update(sf_kwargs)
    sf_model = su.SandFiltration._normalize_model(sf_defaults['model'])
    if sf_model != 'constant_splitting':
        raise ValueError(
            "SVCW requires model='constant_splitting' for its two-outlet "
            'sand-filter topology')
    sf_defaults['model'] = sf_model
    SF = su.SandFiltration(
        'SF', ins=FC-0, outs=('SE', 'SF_backwash'), **sf_defaults)

    # Sludge treatment and anaerobic digestion
    MT = su.IdealClarifier(
        'MT', FC-2, outs=('', 'thickened_WAS'),
        sludge_flow_rate=0.0196*MGD2cmd*1.3685*0.8, #scale based on 13.685/10, scale down to match S_IN in AD with actual sample
        solids_removal_efficiency=0.95,
        )
    M1 = su.Mixer('M1', ins=[GT-1, MT-1])

    cmps_adm = pc.create_adm1p_cmps()
    thermo_adm = qs.get_thermo()
    adm1p_component_sync = _synchronize_adm1p_components(
        cmps_asm, cmps_adm)
    if use_adm1p_default_pKsp:
        adm1p_defaults = pc.ADM1p._pKsp['KM']
        ad_caco3_pKsp = adm1p_defaults[0]
        ad_struvite_pKsp = adm1p_defaults[1]
        ad_acp_pKsp = adm1p_defaults[3]
        ad_mgco3_pKsp = adm1p_defaults[4]

    pKsp_names = (
        'ad_caco3_pKsp', 'ad_struvite_pKsp',
        'ad_acp_pKsp', 'ad_mgco3_pKsp',
        )
    selected_mineral_pKsps = np.array(
        (ad_caco3_pKsp, ad_struvite_pKsp, ad_acp_pKsp, ad_mgco3_pKsp),
        dtype=float,
        )
    for name, value in zip(pKsp_names, selected_mineral_pKsps):
        if not isfinite(value) or value <= 0.:
            raise ValueError(f'{name} must be finite and positive')
    ad_pKsps = list(pc.ADM1p._pKsp['KM'])
    # ADM1p order: CaCO3, struvite, newberyite, ACP, MgCO3, AlPO4, FePO4.
    ad_pKsps[0], ad_pKsps[1], ad_pKsps[3], ad_pKsps[4] = (
        selected_mineral_pKsps)
    adm = pc.ADM1p(
        kLa=10.0,
        pKsp=ad_pKsps,
        k_mmp=_SVCW_AD_K_MMP,
        )
    adm1p_p_residual = (
        np.asarray(adm.stoichio_eval(), dtype=float)
        @ np.asarray(cmps_adm.i_P, dtype=float)
        )
    maximum_p_residual = float(np.max(np.abs(adm1p_p_residual)))
    if maximum_p_residual > 1e-10:
        index = int(np.argmax(np.abs(adm1p_p_residual)))
        raise RuntimeError(
            'ADM1p phosphorus indices are inconsistent with its compiled '
            f'stoichiometry; largest residual is in {adm.IDs[index]!r}: '
            f'{adm1p_p_residual[index]:.6g}.')

    J1 = su.mASM2dtoADM1p(
        'J1', upstream=M1-0, thermo=thermo_adm, isdynamic=True,
        adm1_model=adm, asm2d_model=asm)
    ad_liquid_volume_mgd = float(ad_liquid_volume_mgd)
    if not isfinite(ad_liquid_volume_mgd) or ad_liquid_volume_mgd <= 0.:
        raise ValueError('ad_liquid_volume_mgd must be finite and positive')
    ad_initial_s_ip = float(ad_initial_s_ip)
    if not isfinite(ad_initial_s_ip) or ad_initial_s_ip < 0.:
        raise ValueError('ad_initial_s_ip must be finite and nonnegative')
    if ad_pH is not None:
        ad_pH = float(ad_pH)
        if not isfinite(ad_pH) or not 0. < ad_pH < 14.:
            raise ValueError('ad_pH must be finite and between 0 and 14')
    AD = _SVCWAnaerobicCSTR(
        'AD', ins=J1-0, outs=('biogas', 'digestate'),
        V_liq=ad_liquid_volume_mgd*MGD2cmd,
        V_gas=0.1*ad_liquid_volume_mgd*MGD2cmd,
        fixed_headspace_P=False, fraction_retain=0,
        T=T_ad, model=adm, pH_ctrl=ad_pH,
        )
    AD.algebraic_h2 = False
    # Remove solids in the ADM1p basis before J2 converts residual biomass.
    DW = su.IdealClarifier(
        'DW', AD-1, outs=('', 'cake'), thermo=thermo_adm,
        sludge_flow_rate=5.77e-3*MGD2cmd,
        solids_removal_efficiency=0.9,
        )
    J2 = su.ADM1ptomASM2d(
        'J2', upstream=DW-0, thermo=thermo_asm, isdynamic=True,
        adm1_model=adm, asm2d_model=asm)
    qs.set_thermo(thermo_asm)
    M2 = su.Mixer('M2', ins=[GT-0, MT-0, J2-0, SF-1])
    HD = su.HydraulicDelay('HD', ins=M2-0, outs=1-PC)

    # Initial conditions
    if default_init_conds:
        _set_tf_init_conc(TF, _asinit)
        ASR.set_init_conc(concentrations=_asinit)
        FC.set_init_solubles(**_fcinit)
        FC.set_init_sludge_solids(**_fcinit)
        FC.set_init_TSS(default_fctss_init)
        ad_init = dict(_adinit)
        ad_init['S_IP'] = ad_initial_s_ip
        if equilibrate_ad_minerals and ad_pH is not None:
            ad_init = _equilibrate_ad_mineral_init(
                adm, ad_init, ad_pH, T_ad,
                saturation_target=0.99,
                )
        # These are imposed t=0 concentrations, not pre-equilibration guesses.
        ad_init['S_Mg'] = _SVCW_AD_INITIAL_S_MG
        ad_init['S_Ca'] = _SVCW_AD_INITIAL_S_CA
        AD.set_init_conc(**ad_init)

    # System assembly and default tracking
    system = _SVCWSystem(
        ID,
        path=(PC, GT, TF, ASR, FC, SF, MT, M1, J1, AD, DW, J2, M2, HD),
        recycle=(FC-1, HD-0),
        )
    system.adm1p_component_sync = adm1p_component_sync
    system.adm1p_p_balance_residual = adm1p_p_residual.copy()
    system.set_dynamic_tracker(
        PC-0, TF, TF-0, ASR-0, SF, SF-0, M1-0, AD-1, DW-0, J2-0,
        )
    return system


# %% Interactive entry point
if __name__ == '__main__':
    sys = create_svcw_system()
    globals().update(sys.flowsheet.to_dict())
