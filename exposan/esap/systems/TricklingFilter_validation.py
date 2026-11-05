# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

Standalone RWW -> TricklingFilter -> effluent system for ESAP validation.

Usage example
-------------

The following example creates and simulates the validation system, then reports
the influent and effluent organic concentrations, BOD removal, and biofilm
thickness::

    from exposan.esap.systems.TricklingFilter_validation import (
        create_trickling_filter_validation_system,
        create_trickling_filter_validation_model,
    )

    sys = create_trickling_filter_validation_system(
        wastewater_strength='mid',
        gpsx_aob_proxy=True,
        tf_kwargs={
            'N_section': 6,
            'N_biofilm_layer': 4,
            'filter_bed_surface': 9464,
            'recirculation_ratio': 1,
            'dynamic_biofilm_thickness': True,
            'initial_biofilm_thickness': 4e-4,
        },
    )

    sys.simulate(
        method='BDF', t_span=(0, 300),
        rtol=1e-3, atol=1e-2,
        stop_at_steady_state=True,
        steady_state_basis='inventory',
        steady_state_window=5,
        steady_state_rtol=0.01,
    )

    TF = sys.flowsheet.unit.TF
    for label, stream in [('Influent', TF.ins[0]), ('Effluent', TF.outs[0])]:
        print(f'\nTF {label}')
        print(f"Soluble COD:     {stream.composite('COD', particle_size='s'):.3f} mg/L")
        print(f"Particulate COD: {stream.composite('COD', particle_size='x'):.3f} mg/L")
        print(f"Total COD:       {stream.COD:.3f} mg/L")
        print(f"Total BOD:       {stream.BOD:.3f} mg/L")

    BOD_in, BOD_out = TF.ins[0].BOD, TF.outs[0].BOD
    if BOD_in > 0:
        print(f'\nBOD removal efficiency: {100 * (BOD_in - BOD_out) / BOD_in:.2f}%')
    else:
        print('\nBOD removal efficiency: undefined (influent BOD is zero)')

    import numpy as np
    np.set_printoptions(precision=7, suppress=True)
    TF.biofilm_thickness
    TF.outs[0]

'''

from collections import deque
from math import isfinite

import numpy as np
from scipy.integrate import LSODA
import qsdsan as qs
from qsdsan import process_models as pc
from qsdsan.utils import load_data, ospath
from exposan.werf import data_path, default_rww
from exposan.werf._units import TricklingFilter
from exposan.werf.systems.D1 import (
    _ValidationBDF, _tf_tolerance_solver,
    )

__all__ = (
    'create_trickling_filter_validation_system',
    'create_trickling_filter_validation_model',
    )

ID = 'TricklingFilter_validation'


# %% Programmable startup for validation and site-specific component sets
asinit = load_data(ospath.join(data_path, 'initial_conditions.xlsx'), sheet='E2P')


def _set_tf_init_conc(TF, asinit, *, pao_fraction=0.01,
                      pp_saturation=0.5, pha_ratio=0.05):
    # CODEX removes the following 4 lines
    # '''
    # Initialize D1 trickling filter with a reproducible, physically motivated
    # profile instead of copying one suspended-growth state everywhere.
    # '''
    # CODEX adds the followings 7 lines
    '''
    Initialize SVCW and validation liquid/biofilm profiles from E2P composition.

    Biofilm particulate states are normalized in each layer to
    biofilm_density * biofilm_dry_content [mg dry solids/L biofilm], retaining
    species ratios after seeding PAOs. Liquid, soluble profiles and prescribed
    thickness are unchanged. The startup assumptions (not measured composition)
    allocate pao_fraction of X_H + X_PAO COD to PAOs, set PP/PAO to
    pp_saturation*K_MAX, and set PHA/PAO to pha_ratio. Initial elemental
    inventories may change; no material is added during integration.
    Set pao_fraction=None to retain the legacy E2P storage composition.
    Non-reactive filters also retain that composition.
    '''
    cmps = TF.components
    IDs = cmps.IDs
    N = TF.N_section
    L = TF.N_biofilm_layer
    n_cmps = len(cmps)
    init = np.zeros((N, L+1, n_cmps))
    n_rows = len(asinit)
    # CODEX adds the followings 17 lines
    seed_pao = L and TF._model is not None and pao_fraction is not None
    if seed_pao:
        for name, value, upper in (('pao_fraction', pao_fraction, 1.),
                                    ('pp_saturation', pp_saturation, 1.),
                                    ('pha_ratio', pha_ratio, np.inf)):
            if not np.isfinite(value) or not 0. < value < upper:
                raise ValueError(f'{name} must be finite and between 0 and {upper}.')
        required = ('X_H', 'X_PAO', 'X_PP', 'X_PHA')
        if not all(ID in IDs for ID in required):
            raise ValueError('PAO initialization requires X_H, X_PAO, X_PP and X_PHA.')
        params = getattr(getattr(TF._model, 'rate_function', None), 'params', {})
        kmax = params.get('K_MAX')
        if kmax is None or not np.isfinite(kmax) or kmax <= 0.:
            raise ValueError('PAO initialization requires a positive model K_MAX.')
        xh, xpao, xpp, xpha = cmps.indices(required)
        pp_ratio = pp_saturation * kmax

    def get_idx(ID):
        return cmps.index(ID) if ID in IDs else None

    DO = get_idx('S_O2')
    soluble_organics = [get_idx(i) for i in ('S_F', 'S_A')]
    soluble_organics = [i for i in soluble_organics if i is not None]
    NH4 = get_idx('S_NH4')
    XH = get_idx('X_H')
    XAUT = get_idx('X_AUT')

    if N > 1:
        section_factors = np.linspace(1.0, 0.45, N)
    else:
        section_factors = np.ones(N)

    for i in range(N):
        # Borrow F1/E2P section states as the baseline so the initial
        # component magnitudes remain traceable to existing WERF runs.
        row = asinit.iloc[min(i, n_rows-1)].to_dict()
        base = cmps.kwarray(row)
        section_factor = section_factors[i]

        # The upper trickling-filter sections receive fresher influent, so
        # soluble organics and attached heterotrophs start higher near the top.
        liquid = base.copy()
        if DO is not None:
            liquid[DO] = 2.0
        liquid[soluble_organics] *= section_factor
        if XH is not None:
            liquid[XH] *= section_factor
        init[i, 0, :] = liquid

        for k in range(L):
            # Oxygen enters from the liquid film, so inner biofilm layers start
            # oxygen-limited while outer layers remain mildly aerobic.
            if L > 1:
                depth_fraction = k/(L-1)
            else:
                depth_fraction = 0.
            depth_factor = max(0.25, 1.0 - 0.55*depth_fraction)
            bio = base.copy()
            if DO is not None:
                bio[DO] = max(0.05, 0.35 - 0.30*depth_fraction)
            bio[soluble_organics] *= section_factor * depth_factor
            if NH4 is not None:
                bio[NH4] *= max(section_factor, 0.4)
            if XH is not None:
                bio[XH] *= 1.5 * section_factor * depth_factor
            if XAUT is not None:
                bio[XAUT] *= 0.8 * depth_factor
            # CODEX adds the followings 8 lines
            if seed_pao:
                # Reallocate heterotrophic COD instead of adding active biomass.
                heterotrophs = bio[xh] + bio[xpao]
                bio[xpao] = pao_fraction * heterotrophs
                bio[xh] = (1. - pao_fraction) * heterotrophs
                bio[xpp] = pp_ratio * bio[xpao]
                bio[xpha] = pha_ratio * bio[xpao]
                # The common dry-mass normalization below preserves these ratios.
            init[i, k+1, :] = bio

    # CODEX adds the followings 19 lines
    if L:
        # Convert COD-based particulate states to dry mass before scaling.
        # Each layer represents packed biofilm at wet density times dry fraction;
        # E2P supplies composition, not the attached-solids inventory.
        particulate = np.array([
            cmp.particle_size == 'Particulate' for cmp in cmps
            ])
        biofilm = init[:, 1:, :]
        mass_factors = np.asarray(cmps.i_mass)[particulate]
        target = TF.biofilm_density * TF.biofilm_dry_content  # mg dry solids/L
        dry_solids = biofilm[:, :, particulate] @ mass_factors
        if not np.isfinite(target) or target <= 0.:
            raise ValueError('Initial biofilm dry density must be finite and positive.')
        if not np.isfinite(dry_solids).all() or (dry_solids <= 0.).any():
            raise ValueError('Initial biofilm profiles require positive finite dry solids.')
        # Scale all particulate species together, preserving composition ratios.
        # Liquid and biofilm solubles and the user-specified thickness are unchanged.
        biofilm[:, :, particulate] = biofilm[:, :, particulate] \
            * (target / dry_solids)[:, :, None]

    # CODEX removes the following 1 line
    # TF.set_init_conc(init)
    # CODEX adds the followings 4 lines
    if L:
        TF.set_init_conc(init)
    else:
        TF.set_init_conc(init[:, 0, :])

    # Start with thin, equal biofilm layers rather than full layers. This keeps
    # the model away from the GPS-X layer-filling switch at t=0 and makes later
    # convergence tests easier to reproduce.
    if L and TF._initial_biofilm_thickness is None:
        layer_H = min(1e-4, 0.5*TF.max_biofilm_thickness/L)
        TF.initial_biofilm_thickness = np.full((N, L), layer_H)


# %% Optional benchmark termination and reporting
class _SteadyStateWindow:
    """Check accepted states over a complete rolling window, including extrema."""

    def __init__(self, window, rtol, floor):
        if not np.isfinite(window) or window <= 0:
            raise ValueError('steady_state_window must be finite and positive')
        if not np.isfinite(rtol) or not 0 < rtol < 1:
            raise ValueError('steady_state_rtol must be between zero and one')
        self.window, self.rtol = window, rtol
        self.floor = np.asarray(floor)
        self.history = deque()
        self.variation = np.inf

    def update(self, t, y):
        if self.history and t <= self.history[-1][0]:
            raise ValueError('Steady-state checks require increasing accepted times')
        self.history.append((float(t), np.array(y, copy=True)))
        cutoff = t - self.window
        while len(self.history) > 2 and self.history[1][0] <= cutoff:
            self.history.popleft()
        if self.history[0][0] > cutoff or len(self.history) < 2:
            return False
        t0, y0 = self.history[0]
        t1, y1 = self.history[1]
        boundary = y0 + (y1-y0)*((cutoff-t0)/(t1-t0))
        values = np.vstack([boundary] + [v for _, v in list(self.history)[1:]])
        scale = np.maximum(np.max(np.abs(values), axis=0), self.floor)
        self.variation = float(np.max(np.ptp(values, axis=0)/scale))
        return np.isfinite(self.variation) and self.variation < self.rtol


class TricklingFilterValidationSystem(qs.System):
    """Single-filter benchmark with BDF/LSODA steady-state termination.

    simulate(..., stop_at_steady_state=True, steady_state_window=5.,
    steady_state_rtol=0.01) stops after all liquid/biofilm component states and
    section thicknesses vary by less than 1% over five days. The range of
    sampled accepted states is checked, not just endpoints. Near-zero concentrations
    use a 1 mg/L relative scale floor; thickness uses 1e-6 m. These are stopping
    criteria, not integration tolerances. Use stop_at_steady_state=False for
    time-varying forcing or a prescribed full-duration trajectory.

    BDF and LSODA accept ordinary scalar rtol/atol. Other integration methods
    run for the requested interval. steady_state_report records whether
    convergence was reached; reaching t_span[1] alone does not imply convergence.

    steady_state_basis='inventory' checks all section liquid concentrations,
    component inventories per media area, and thicknesses. Internal layer
    redistribution may continue at termination; use 'states' (default) to
    require each layer to settle. Checks occur at most every min(0.25, window/20)
    days, on accepted steps. A five-day window is a local stationarity test,
    not proof that slow changes cannot occur later.
    """

    def nitrogen_summary(self):
        """Return influent/effluent N species in mg N/L after simulation.

        NOx is the lumped S_NO3 state; no separate NO2 prediction is available.
        Organic/other N is TN minus NH4-N and NOx-N, excluding dissolved N2.
        """
        tf = self.flowsheet.unit.TF
        result = {}
        for label, stream in (('influent', tf.ins[0]), ('effluent', tf.outs[0])):
            nh4 = float(stream.iconc['S_NH4'])
            nox = float(stream.iconc['S_NO3'])
            total = float(stream.TN)
            result[label] = {
                'NH4-N': nh4, 'NOx-N': nox, 'TN': total,
                'organic/other N': total-nh4-nox,
                'particulate N': float(stream.composite('N', particle_size='x')),
            }
        return result

    def get_influent_pumping_energy(self):
        """Return raw-influent pumping power in kW."""
        from exposan.esap.models import get_influent_pumping_power
        return get_influent_pumping_power(self.flowsheet.stream.RWW)

    def get_recirculation_pumping_energy(self):
        """Return trickling-filter recirculation pumping power in kW."""
        from exposan.esap.models import (
            get_trickling_filter_recirculation_pumping_power,
            )
        return get_trickling_filter_recirculation_pumping_power(
            self.flowsheet.unit.TF)

    def get_rotating_arm_energy(self):
        """Return rotary-distributor power in kW."""
        from exposan.esap.models import get_trickling_filter_rotating_arm_power
        return get_trickling_filter_rotating_arm_power(
            self.flowsheet.unit.TF)

    def get_pumping_energy(self):
        """Return influent plus recirculation pumping power in kW."""
        return (self.get_influent_pumping_energy()
                + self.get_recirculation_pumping_energy())

    def get_energy_summary(self):
        """Return separated validation electricity demands in kW."""
        return {
            'influent pumping energy [kW]': self.get_influent_pumping_energy(),
            'recirculation pumping energy [kW]':
                self.get_recirculation_pumping_energy(),
            'rotating arm energy consumption [kW]':
                self.get_rotating_arm_energy(),
            }

    def dynamic_run(self, **kwargs):
        """Integrate with targeted tolerances and optional accepted-step monitoring."""
        # Inventory mode permits redistribution within a section at steady output.
        # It is an engineering stopping criterion, not full layer equilibrium.
        basis = kwargs.pop('steady_state_basis', 'states')
        if basis not in ('states', 'inventory'):
            raise ValueError("steady_state_basis must be 'states' or 'inventory'")
        self.steady_state_basis = basis
        stop = kwargs.pop('stop_at_steady_state', True)
        window = float(kwargs.pop('steady_state_window', 5.))
        tolerance = float(kwargs.pop('steady_state_rtol', 0.01))
        method = kwargs.get('method', self.dynsim_kwargs.get('method', 'BDF'))
        self.steady_state_report = dict(converged=False, window_days=window,
                                       relative_tolerance=tolerance)
        if not stop or method not in ('BDF', 'LSODA'):
            kwargs['method'] = _tf_tolerance_solver(method, self)
            try:
                return super().dynamic_run(**kwargs)
            finally:
                self.dynsim_kwargs['method'] = method
        tf = self.flowsheet.unit.TF
        N, L, n = tf.N_section, tf.N_biofilm_layer, len(tf.components)
        floor = np.ones(N*tf._ncol + N*L*n + (N if L else 0))
        if L:
            floor[-N:] = 1e-6
        def monitored_state(y):
            if basis == 'states' or not L:
                return y
            liquid = y[:N*tf._ncol].reshape(N, tf._ncol)[:, :-1]
            bio_stop = N*tf._ncol + N*L*n
            biofilm = y[N*tf._ncol:bio_stop].reshape(N, L, n)
            thickness = y[bio_stop:bio_stop+N]
            # Equal layer widths: mean(C)*H is inventory per media area.
            # Divide by Hmax to retain mg/L-scale floors for the monitor.
            inventory = biofilm.mean(axis=1)*(
                thickness[:, None]/tf.max_biofilm_thickness)
            return np.concatenate((liquid.ravel(), inventory.ravel(), thickness))

        if basis == 'inventory' and L:
            floor = np.ones(2*N*n + N)
            floor[-N:] = 1e-6
        # Check at most four times per day; avoid scanning a growing history
        # at every tiny startup solver step. Solver max_step still bounds gaps.
        check_interval = min(window/20., 0.25)
        last_check = [None]
        self.steady_state_report['basis'] = basis
        self.steady_state_report['criterion'] = (
            'all states' if basis == 'states' else
            'liquid concentrations, section component inventories, thickness')
        monitor = _SteadyStateWindow(window, tolerance, floor)
        report = self.steady_state_report

        base_solver = _ValidationBDF if method == 'BDF' else LSODA
        class SteadySolver(base_solver):
            def __init__(solver, fun, t0, y0, t_bound, **options):
                if t_bound <= t0:
                    raise ValueError('Steady-state stopping requires forward integration')
                # Several accepted samples per window prevent an endpoint-only check.
                options['max_step'] = min(options.get('max_step', np.inf), window/5)
                # System._state may alias y0; preserve the solver's initial record.
                super().__init__(fun, t0, y0.copy(), t_bound, **options)
                monitor.update(t0, monitored_state(solver.y))
                last_check[0] = t0

            def step(solver):
                message = super().step()
                if solver.status != 'failed':
                    if solver.t-last_check[0] < check_interval and solver.status != 'finished':
                        return message
                    last_check[0] = solver.t
                    converged = monitor.update(solver.t, monitored_state(solver.y))
                    report['monitor_checks'] = report.get('monitor_checks', 0) + 1
                    report.update(time=solver.t, max_variation=monitor.variation)
                    if converged:
                        report.update(converged=True, state=solver.y.copy())
                        solver.status = 'finished'
                return message

        kwargs['method'] = _tf_tolerance_solver(SteadySolver, self)
        try:
            result = super().dynamic_run(**kwargs)
            sol = self.scope.sol
            if report['converged'] and sol.status == 0:
                final = report.pop('state')
                # Preserve the actual stopping point even with user-provided t_eval.
                if not len(sol.t) or sol.t[-1] < report['time']:
                    sol.t = np.append(sol.t, report['time'])
                    sol.y = np.column_stack((sol.y, final))
                sol.message = (f"Steady state reached: <{100*tolerance:g}% variation "
                               f"over {window:g} days at t={report['time']:.6g} d.")
                # Publish the accepted endpoint, not the last numerical-Jacobian probe.
                self.DAE(report['time'], final.copy())
                self._write_state()
            return result
        finally:
            self.dynsim_kwargs['method'] = method
            report.pop('state', None)


# %% Calibration profiles: kinetic matches, wastewater, and component indices
# Direct matches transcribed from the user's kinetic_screenshot.jpg.
# Only mismatches are listed; units and previous mASM2d defaults are recorded.
_GPSX_KINETIC_OVERRIDES = dict(
    mu_H=3.2,             # d^-1; heterotrophic growth (was 6.0)
    K_F=5.0,              # mgCOD/L; Ss / S_F saturation (was 4.0)
    K_A_H=5.0,            # mgCOD/L; acetate saturation (was 4.0)
    eta_NO3=0.28,         # -; anoxic hydrolysis reduction (was 0.6)
    q_PHA=6.0,            # gCOD/(gCOD-PAO*d); PHA storage (was 3.0)
    q_PP=2.0,             # gP/(gCOD-PAO*d); poly-P storage (was 1.5)
    b_PP=0.05,            # d^-1; poly-P lysis (was 0.2)
    eta_NO3_PAOl=0.9,     # -; anoxic PAO biomass decay reduction (was 0.33)
    b_AUT=0.17,           # d^-1; same for BOTH GPS-X AOB and NOB (was 0.15)
    eta_NO3_AUTl=0.5,     # -; same for BOTH GPS-X AOB and NOB (was 0.33)
    )
# Influent profiles from analyses/F1_HA_uncertainty.py (Metcalf & Eddy Table 3-18).
_WASTEWATER_ORGANIC_FRACTIONS = dict(
    fr_SI=0.075, fr_SF=0.24, fr_SA=0.05, fr_XI=0.3,
    )
_WASTEWATER_STRENGTHS = {
    'low': dict(COD=339., NH4_N=14., PO4_P=1.6, S_K=11., S_Cl=39.),
    'mid': dict(COD=508., NH4_N=20., PO4_P=2.4, S_K=16., S_Cl=59.),
    'high': dict(COD=1016., NH4_N=41., PO4_P=4.7, S_K=32., S_Cl=118.),
    }
_NONDEFAULT_COMPONENT_PROPERTIES = {
    'S_F': dict(i_N=0.023464, i_P=0.006149),
    'S_I': dict(i_N=0.042021, i_P=0.007139),
    'X_S': dict(i_N=0.023464, i_P=0.006149, i_mass=0.63),
    'X_I': dict(i_N=0.042021, i_P=0.007139, i_mass=0.63),
    'X_H': dict(i_N=0.060305, i_P=0.02369),
    'X_PAO': dict(i_N=0.060305, i_P=0.02369),
    'X_AUT': dict(i_N=0.060305, i_P=0.02369),
    }


# %% Shared component and influent construction (also used by SF and SVCW)
def _normalize_wastewater_strength(strength):
    """Return a validated, lower-case validation wastewater profile name."""
    strength = str(strength).lower()
    if strength != 'default' and strength not in _WASTEWATER_STRENGTHS:
        raise ValueError(
            "wastewater_strength must be 'default', 'low', 'mid', or 'high', "
            f"not {strength!r}")
    return strength


def _create_validation_components(strength):
    """Create local mASM2d components with Table S1 non-default indices."""
    cmps = pc.create_masm2d_cmps()
    if strength == 'default':
        return cmps
    for component_ID, properties in _NONDEFAULT_COMPONENT_PROPERTIES.items():
        component = getattr(cmps, component_ID)
        for name, value in properties.items():
            setattr(component, name, value)
    cmps.refresh_constants()
    return cmps


def _create_validation_wastewater(strength):
    """Return a 10-MGD, 20-degC influent for the selected strength profile."""
    strength = _normalize_wastewater_strength(strength)
    if strength == 'default':
        return default_rww()
    composition = _WASTEWATER_STRENGTHS[strength]
    return pc.create_masm2d_inf(
        'RWW', 10., 'MGD', T=293.15,
        **composition, **_WASTEWATER_ORGANIC_FRACTIONS)


# %% Process-model construction and explicit GPS-X approximations
# Do not map separate NOx reduction steps to the lumped denitrification factor,
# or select one guild's growth/O2 constants for mASM2d's lumped nitrifiers.
# GPS-X inert residue (0.001/d) and inert organics (0.03/d) hydrolysis have no
# direct mASM2d rate: X_I pools both sources and is not hydrolyzed.
# The visible mineral constants match the existing Musvoto preset; leave
# mmp_kinetics=None (precipitation disabled) unless the caller enables it.


def _create_validation_model(gpsx_decay, asm_kwargs, anaerobic_decay_factor,
                             gpsx_kinetics=False, gpsx_aob_proxy=False):
    """Apply independently switchable GPS-X decay and direct kinetic matches."""
    parameters = dict(electron_acceptor_dependent_decay=True)
    if gpsx_decay:
        # Values supplied from the user's GPS-X benchmark, not mASM2d defaults.
        # k_h/K_X match; the separate kinetic profile adjusts anoxic hydrolysis.
        parameters.update(b_H=0.62, eta_NO3_Hl=0.9, f_XI_H=0.08, f_SI=0.)
    if gpsx_kinetics:
        parameters.update(_GPSX_KINETIC_OVERRIDES)
    if gpsx_aob_proxy:
        # kinetic_screenshot.jpg: AOB mu [1/d], NH4 [mgN/L], O2 [mgO2/L].
        # Approximation for lumped X_AUT, not a separate nitrite pathway.
        # Keep user overrides last; do not map multistep denitrification factors.
        parameters.update(mu_AUT=0.9, K_NH4_AUT=0.7, K_O2_AUT=0.25, b_AUT=0.17)
    parameters.update(asm_kwargs or {})
    if gpsx_decay and not parameters['electron_acceptor_dependent_decay']:
        raise ValueError('gpsx_decay=True requires electron_acceptor_dependent_decay=True')
    asm = pc.mASM2d(**parameters)
    if not gpsx_decay:
        return asm
    factor = float(anaerobic_decay_factor)
    if not isfinite(factor) or not 0. <= factor <= 1.:
        raise ValueError('anaerobic_decay_factor must be finite and between 0 and 1')

    kinetics = asm.rate_function
    base_rate_function = kinetics.function
    decay_idx = asm.index('hetero_lysis')
    oxygen_idx, nitrate_idx, biomass_idx = kinetics.params['cmps'].indices(
        ('S_O2', 'S_NO3', 'X_H'))

    def rates_with_anaerobic_decay(state, params):
        # Preserve mASM2d hydrolysis, growth, mineral reactions, and other decay.
        # Copy because the base rate function reuses its output buffer.
        rates = base_rate_function(state, params).copy()
        oxygen_inhibition = params['K_O2_H'] / (
            params['K_O2_H'] + state[oxygen_idx])
        nitrate_inhibition = params['K_NO3_H'] / (
            params['K_NO3_H'] + state[nitrate_idx])
        # Reuse mASM2d's smooth redox gates as an explicit approximation:
        # b_H*X_H*[aerobic + eta_anoxic*anoxic + eta_anaerobic*anaerobic].
        # Only the last contribution is missing from the base rate function.
        rates[decay_idx] += (
            params['b_H'] * state[biomass_idx] * factor
            * oxygen_inhibition * nitrate_inhibition)
        return rates

    kinetics.function = rates_with_anaerobic_decay
    return asm


# %% Standalone system and metric-model factories
def create_trickling_filter_validation_system(
        flowsheet=None, default_init_conds=True, tf_kwargs=None,
        gpsx_decay=True, asm_kwargs=None, anaerobic_decay_factor=0.6,
        gpsx_kinetics=True, wastewater_strength='default',
        gpsx_aob_proxy=False, tf_init_kwargs=None):
    '''
    Create a single-unit trickling-filter system fed directly by D1's raw
    wastewater (10 MGD, 20 degC, from default_rww).

    Parameters
    ----------
    flowsheet : qsdsan.Flowsheet, optional
        Flowsheet to populate; a new validation flowsheet is used by default.
    default_init_conds : bool, optional
        Initialize liquid and biofilm concentrations with D1's E2P-based
        profiles. Disabling this permits custom concentrations before simulation;
        the specified initial biofilm thickness still applies.
    tf_kwargs : dict, optional
        TricklingFilter keyword overrides, as in create_d1_system. Defaults
        are six sections, five biofilm layers, 200 m2 plan area,
        2 m bed depth, 100 m2/m3 media surface, and no recirculation.
        D1 instead defaults to 9464 m2 area and recirculation ratio one.
        Initial total biofilm thickness defaults to 3e-4 m (0.3 mm) in every
        section. Thickness dynamics and passive oxygen transfer follow the
        TricklingFilter defaults; both remain adjustable through this dict.

        The validation profile uses detachment_model='capacity' with exponent
        one: J_det = k_det * (H/(Hmax-H)), scaled by surface dry density.
        This is an approximation to GPS-X's increasing sloughing near Hmax,
        not its unpublished exact equation. Set detachment_model='linear'
        to select the previous capped linear law. The validation detachment
        scale is calibrated to 0.065 kg dry solids/(m2*d), versus the GPS-X
        input of 0.070. This modest adjustment improves thickness and BOD in
        paired 10-day tests; it is not a transcribed GPS-X parameter. Override
        tf_kwargs['detachment_rate']=0.07 to restore the previous scale.
        No section-specific thickness targets are imposed.

    wastewater_strength : {'default', 'low', 'mid', 'high'}, optional
        Influent profile. ``'default'`` preserves D1's existing raw wastewater.
        ``'low'``, ``'mid'``, and ``'high'`` use the corresponding 10-MGD
        COD, NH4-N, PO4-P, K, Cl, and
        organic-COD fractions from analyses/F1_HA_uncertainty.py while keeping
        the validation temperature at 20 degC. Their COD values are 339, 508,
        and 1016 mg/L, respectively. Matching is case-insensitive.
        For any non-default profile, component N and P contents and the
        X_S/X_I COD-to-mass indices use Table S1 values from the validation
        source. These changes are made on a fresh component set owned by this
        system; they do not modify default-strength or previously created
        systems. The corrected project biomass i_P is 0.02369 gP/gCOD;
        D1 and F1 retain the stock mASM2d component definitions.

    gpsx_decay : bool, optional
        Default True applies the user's GPS-X heterotrophic decay settings:
        b_H=0.62 d^-1, eta_NO3_Hl=0.9, f_XI_H=0.08 gCOD/gCOD, and f_SI=0.
        Adds anaerobic heterotrophic decay using mASM2d's O2/NO3 inhibition
        gates; this is not a full MANTIS3 implementation. False restores the
        original mASM2d decay settings, unless overridden by asm_kwargs.
        D1 and other WERF configurations are not changed.
    tf_init_kwargs : dict, optional
        Startup biofilm settings: pao_fraction=0.01 of X_H+X_PAO COD,
        pp_saturation=0.5 of K_MAX, and pha_ratio=0.05 gCOD/gCOD-PAO.
        These assumptions preserve layer dry density, not the legacy elemental
        inventory. Liquid profiles and thickness are unchanged. Use
        {'pao_fraction': None} for legacy initialization. Only used when
        default_init_conds=True. The same default seed is shared by D1/SVCW.
    asm_kwargs : dict, optional
        Overrides passed to mASM2d, including k_h [d^-1], K_X [gCOD/gCOD],
        b_H [d^-1], eta_NO3_Hl [-], and f_XI_H [gCOD/gCOD]. Overrides take
        precedence over both GPS-X parameter options. Growth yields, influent
        composition, and transport are not changed by those options. gpsx_decay=True
        requires electron_acceptor_dependent_decay=True. A model supplied
        through tf_kwargs['suspended_growth_model'] takes precedence.
    anaerobic_decay_factor : float, optional
        Heterotrophic anaerobic decay reduction factor [-], between 0 and 1.
        Default 0.6 gives b_H*0.6=0.372 d^-1 when O2 and NO3 are both zero.
        Ignored when gpsx_decay=False. Recreate the model to change this factor.

    gpsx_kinetics : bool, optional
        Default True applies direct kinetic matches from kinetic_screenshot.jpg:
        mu_H=3.2 d^-1, K_F=K_A_H=5 mgCOD/L, eta_NO3=0.28,
        q_PHA=6 gCOD/(gCOD-PAO*d), q_PP=2 gP/(gCOD-PAO*d),
        b_PP=0.05 d^-1, eta_NO3_PAOl=0.9, b_AUT=0.17 d^-1,
        eta_NO3_AUTl=0.5. Independent of gpsx_decay; set both flags False
        to reproduce the original validation model before GPS-X adjustments.
        AOB/NOB growth and O2 constants, multistep denitrification, colloidal
        adsorption, and inert-material hydrolysis are not direct matches.
        Precipitation remains disabled unless enabled through asm_kwargs.

    gpsx_aob_proxy : bool, optional
        Default False. True approximates lumped autotrophic kinetics using
        GPS-X AOB values from kinetic_screenshot.jpg: mu_AUT=0.9 d^-1,
        K_NH4_AUT=0.7 mg N/L, K_O2_AUT=0.25 mg O2/L. This is an AOB proxy,
        not separate AOB/NOB kinetics or a fitted high-strength calibration.
        asm_kwargs takes precedence. S_NO3 must be compared with GPS-X's
        nitrate plus nitrite, never nitrite alone. Heterotrophic denitrification
        can be explored with asm_kwargs['eta_NO3_H']; multistep GPS-X factors
        and methylotrophic kinetics have no direct single-parameter mapping.

    Returns
    -------
    qsdsan.System
        One TF unit connecting RWW to effluent, with all three tracked.
        Raw wastewater enters directly, without D1's primary clarification
        or sludge-treatment return flows.
        With method='BDF' or 'LSODA', simulation stops once all state variations
        stay below 1% over five days. See TricklingFilterValidationSystem for
        stopping controls and near-zero scale floors. For changing influent,
        pass stop_at_steady_state=False. Scalar integration tolerances suffice.
    '''
    flowsheet = flowsheet or qs.Flowsheet(ID)
    qs.main_flowsheet.set_flowsheet(flowsheet)

    wastewater_strength = _normalize_wastewater_strength(wastewater_strength)
    _create_validation_components(wastewater_strength)
    asm = _create_validation_model(
        gpsx_decay, asm_kwargs, anaerobic_decay_factor, gpsx_kinetics,
        gpsx_aob_proxy)
    rww = _create_validation_wastewater(wastewater_strength)

    tf_kwargs = dict(tf_kwargs or ())
    n_zones = int(tf_kwargs.get('N_section', 6))
    tf_defaults = dict(
        N_section=n_zones,
        N_biofilm_layer=5,
        influent_fractions=[[1] + [0]*(n_zones-1)],
        recirculation_ratio=0.0,
        suspended_growth_model=asm,
        DO_ID='S_O2',
        max_liquid_film_thickness=5e-5,
        max_biofilm_thickness=1e-3,
        filter_bed_depth=2.0,
        filter_bed_surface=200.0,
        media_specific_surface=100.0,
        initial_biofilm_thickness=3e-4,
        detachment_model='capacity',
        detachment_exponent=1.,
        # Validation calibration: 7.1% below the GPS-X input of 0.070.
        # Paired 10-day runs improved BOD and thickness with little nitrate change.
        # The general TricklingFilter default remains at 0.070; D1 uses 0.065.
        detachment_rate=0.065,
        )
    tf_defaults.update(tf_kwargs)

    TF = TricklingFilter('TF', ins=rww, outs='effluent', **tf_defaults)
    if default_init_conds:
        _set_tf_init_conc(TF, asinit, **(tf_init_kwargs or {}))

    sys = TricklingFilterValidationSystem(ID, path=(TF,))
    sys.set_dynamic_tracker(rww, TF, TF-0)
    return sys


def create_trickling_filter_validation_model(system=None, **system_kwargs):
    """Return an energy model for a standalone validation system.

    The model includes TF pumping power, rotary-distributor power, and their
    electricity cost. Metric evaluation reads the current system state and does
    not rerun the dynamic simulation.
    Influent and recirculation pumping are reported separately; the former uses
    50 ft lift plus 3.28 ft loss at 60% efficiency.
    Recirculation also uses 60% efficiency, with bed depth plus 5.5 ft and
    an additional 3.28 ft pump loss.

    Parameters
    ----------
    system : TricklingFilterValidationSystem, optional
        Existing validation system. When omitted, ``system_kwargs`` are passed
        to :func:`create_trickling_filter_validation_system`.
    **system_kwargs
        Validation-system creation arguments; invalid when ``system`` is given.
    """
    if system is None:
        system = create_trickling_filter_validation_system(**system_kwargs)
    elif system_kwargs:
        raise TypeError('system_kwargs cannot be used with an existing system')
    from exposan.esap.models import add_trickling_filter_energy_metrics
    model = qs.Model(system, specification=lambda: None)
    add_trickling_filter_energy_metrics(model)
    return model


# %% Interactive benchmark entry point
if __name__ == '__main__':
    sys = create_trickling_filter_validation_system()
    RWW = sys.flowsheet.stream.RWW
    TF = sys.flowsheet.unit.TF
    effluent = sys.flowsheet.stream.effluent
    sys.simulate(
        method='BDF', t_span=(0, 300),
        rtol=1e-3, atol=1e-2, stop_at_steady_state=True,
        )
    print(sys.steady_state_report)
