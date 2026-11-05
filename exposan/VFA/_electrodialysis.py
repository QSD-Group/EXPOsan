"""Steady-state, well-mixed, single-cell-pair VFA electrodialysis.

Sources: WaterTAP electrodialysis_1D.py, eq_mass_transfer_term_diluate,
eq_conc_polarization_ratio, eq_current_dens_lim_x and resistance equations:
https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/electrodialysis_1D.py
https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/electrodialysis_1D.html
This is a well-mixed adaptation, not WaterTAP's axial discretization.
"""
import numpy as np
from scipy.optimize import minimize, NonlinearConstraint, Bounds
from qsdsan import SanUnit

__all__ = ('Electrodialysis',)
F = 96485.33212  # C/mol, Faraday constant


class Electrodialysis(SanUnit):
    """Size length and ion transfers at fixed current density for one cell pair.

    Parameters
    ----------
    ins, outs : iterable of streams
        Diluate and concentrate, in that order; liquid, charge-balanced feeds.
    charges : dict[str, int]
        All ionic component IDs and signed charges; water must be excluded.
    target : str
        Ion whose net transfer / diluate inlet molar flow defines recovery.
    recovery : float
        Target in (0, 1). Concentrate inlet inventory is excluded.
    current_density : float
        Positive current density [A/m2].
    width, depth : float
        Positive channel width and depth [m], identical for both channels.
    transport_numbers : array_like
        Shape (2, number of ions), ordered CEM, AEM and as charges keys.
        Nonnegative current fractions, each membrane row summing to one.
    permeabilities : array_like
        Same shape; membrane diffusivity / thickness [m/s].
        Back diffusion must jointly conserve charge; incompatible independent
        permeabilities are rejected, rather than corrected by creating ions.
    ionic_conductivities : array_like
        One positive ionic molar conductivity [S m2/mol] per ion.
    limiting_coefficient : float
        i_lim / diluate equivalent concentration [A m/mol]. The equivalent
        concentration is sum(z*c for cations), in mol charge/m3.
        This empirical calibration must match the specified channel hydraulics.
    membrane_resistance : array_like
        Two nonnegative base areal resistances [ohm m2], CEM then AEM.
    membrane_resistance_coefficient : array_like
        Two nonnegative inverse-concentration coefficients [ohm mol/m].
        r_m = r_base + coefficient / mean_bulk_equivalent_concentration.
        This is an explicit provisional closure, not a fitted WaterTAP law.
    polarization : bool, optional
        Default True. Surface concentrations use cD*(1-i/i_lim) and
        cC*(1+i/i_lim). This symmetric-electrolyte approximation for mixed
        VFAs needs validation. Operation at/above limiting current is rejected; a hard nonlinear
        constraint maintains a 1e-6 relative margin, even if polarization is off.
    max_length : float, optional
        Positive upper sizing bound [m], default 1000.

    Notes
    -----
    All transport calculations use mol/s, mol/m3, m, A and V. Stream imol is
    kmol/hr (multiply by 1000/3600). Water and nonionic solutes stay in their
    original channels; outlet volumes are recalculated from stream properties.
    Isothermal, no pressure drop, water splitting, pH adjustment, electrode
    chemistry, nonohmic potentials, pump power, or capital cost is modeled.
    Voltage/power cover the cell pair's ohmic contribution only, not the
    complete redox/electrode apparatus. Liquid resistance uses full depth;
    no additional diffusion-layer ohmic resistance is added (no double count).
    TODO: obtain mixed-VFA transport, limiting-current and resistance data,
    electrolyte activity/property models and configuration-specific electrodes.
    """
    _N_ins = _N_outs = 2
    _units = {'Length': 'm', 'Projected area': 'm2', 'Installed membrane area': 'm2',
              'Current': 'A', 'Voltage': 'V', 'Areal resistance': 'ohm*m2',
              'Limiting current density': 'A/m2', 'Cell-pair power': 'kW'}

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, charges,
                 target, recovery, current_density, width, depth,
                 transport_numbers, permeabilities, ionic_conductivities,
                 limiting_coefficient, membrane_resistance,
                 membrane_resistance_coefficient, polarization=True,
                 max_length=1000., init_with='WasteStream'):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with)
        self.charges = dict(charges)
        self.ions = tuple(charges)
        self.z = np.asarray(tuple(charges.values()), dtype=float)
        self.target = target
        self.recovery = recovery
        self.current_density = current_density
        self.width, self.depth = width, depth
        self.t = np.array(transport_numbers, dtype=float, copy=True)
        self.permeabilities = np.array(permeabilities, dtype=float, copy=True)
        self.conductivities = np.array(ionic_conductivities, dtype=float, copy=True)
        self.limiting_coefficient = limiting_coefficient
        self.membrane_resistance = np.array(membrane_resistance, dtype=float, copy=True)
        self.membrane_resistance_coefficient = np.array(membrane_resistance_coefficient, dtype=float, copy=True)
        self.polarization, self.max_length = polarization, max_length
        self._validate()

    def _validate(self):
        n = len(self.ions)
        if self.target not in self.ions or 'H2O' in self.ions:
            raise ValueError('Target must be an ion, and water cannot be an ion.')
        if not np.all(np.isfinite(self.z)) or np.any(self.z == 0):
            raise ValueError('Charges must be finite and nonzero.')
        if not (np.any(self.z > 0) and np.any(self.z < 0)):
            raise ValueError('Both cations and anions must be tracked.')
        if not 0 < self.recovery < 1:
            raise ValueError('Recovery must be strictly between zero and one.')
        for name in ('current_density', 'width', 'depth', 'limiting_coefficient', 'max_length'):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f'{name} must be finite and positive.')
        for name, shape in (('t', (2, n)), ('permeabilities', (2, n)),
                            ('conductivities', (n,)), ('membrane_resistance', (2,)),
                            ('membrane_resistance_coefficient', (2,))):
            a = getattr(self, name)
            if a.shape != shape or not np.all(np.isfinite(a)) or np.any(a < 0):
                raise ValueError(f'{name} must be nonnegative, finite, shape {shape}.')
        if np.any(self.conductivities <= 0) or not np.allclose(self.t.sum(axis=1), 1., atol=1e-10, rtol=0):
            raise ValueError('Conductivities must be positive; each transport-number row must sum to one.')
        for ID in self.ions:
            self.components.index(ID)

    def _bulk_properties(self, transfers):
        # Trial streams are detached copies; failed solves do not alter outlets.
        d, c = (s.copy(None) for s in self.ins)
        # Inlet ion flows and signed D -> C transfers are all in mol/s.
        # mol/s -> kmol/hr: 3600 s/hr / 1000 mol/kmol = 3.6.
        d.imol[self.ions] = self._diluate_inlet_mol_s * 3.6 - transfers * 3.6
        c.imol[self.ions] = self._concentrate_inlet_mol_s * 3.6 + transfers * 3.6
        # Stream F_vol is m3/hr; divide by 3600 s/hr to obtain m3/s.
        q = np.array([d.F_vol, c.F_vol]) / 3600.
        if np.any(q <= 0) or not np.all(np.isfinite(q)):
            raise ValueError('Both trial channels require positive liquid volume.')
        bulk = np.array([self._diluate_inlet_mol_s - transfers, self._concentrate_inlet_mol_s + transfers]) / q[:, None]
        return (d, c), bulk

    def _polarization(self, bulk):
        # TODO: validate the limiting-current law, concentration basis and
        # coefficient against mixed-VFA data for the chosen membranes/hydraulics;
        # this charge-equivalent closure is not WaterTAP's exact InitialValue law.
        equivalent = bulk[:, self.z > 0] @ self.z[self.z > 0]
        i_lim = self.limiting_coefficient * equivalent[0]
        # ASSUMPTION: all ions share one mixture-level i_lim and therefore the
        # same fractional polarization. Extending the symmetric 1:1-electrolyte
        # relation to mixed VFAs requires validation; individual ion currents
        # still differ through their membrane transport numbers.
        ratio = self.current_density / max(i_lim, 1e-30) if self.polarization else 0.
        # Clipping only extends residuals into infeasible trial regions; a
        # nonlinear constraint and final strict check reject the limiting regime.
        r = min(ratio, 1. - 1e-12)
        surface = bulk * np.array([1. - r, 1. + r])[:, None]
        return surface, i_lim, ratio

    def _flux(self, surface):
        migration = (self.t[0] - self.t[1]) * self.current_density / (self.z * F)
        diffusion = self.permeabilities.sum(axis=0) * (surface[1] - surface[0])
        return migration - diffusion

    def _resistance(self, bulk):
        conductivity = bulk @ self.conductivities
        if np.any(conductivity <= 0):
            raise ValueError('Both outlet conductivities must be positive.')
        eq = bulk[:, self.z > 0] @ self.z[self.z > 0]
        mean_eq = eq.mean()
        if mean_eq <= 0:
            raise ValueError('Mean equivalent concentration must be positive.')
        return (self.depth / conductivity).sum() + np.sum(
            self.membrane_resistance + self.membrane_resistance_coefficient / mean_eq)

    def _run(self):
        self._validate()
        self.design_results.clear()
        self.power_utility.empty()
        for stream in self.ins:
            if stream.phase != 'l' or np.any(np.asarray(stream.mol) < 0) or not np.all(np.isfinite(np.asarray(stream.mol))):
                raise ValueError('ED requires nonnegative liquid feeds.')
        # imol [kmol/hr] -> mol/s: 1000 mol/kmol / 3600 s/hr = 1/3.6.
        self._diluate_inlet_mol_s, self._concentrate_inlet_mol_s = (np.asarray(s.imol[self.ions], dtype=float) / 3.6 for s in self.ins)
        for inlet in (self._diluate_inlet_mol_s, self._concentrate_inlet_mol_s):
            if abs(inlet @ self.z) > 1e-9 * max(np.abs(inlet * self.z).sum(), 1e-12):
                raise ValueError('Inlet ions must be charge balanced; no counterions are silently added.')
        k = self.ions.index(self.target)
        target = self.recovery * self._diluate_inlet_mol_s[k]
        if target <= 0:
            raise ValueError('Target species must be present in diluate feed.')
        active = self._diluate_inlet_mol_s + self._concentrate_inlet_mol_s > 0
        scale = np.maximum(self._diluate_inlet_mol_s + self._concentrate_inlet_mol_s, target * 1e-4)
        lower, upper = -self._concentrate_inlet_mol_s[active] / scale[active], self._diluate_inlet_mol_s[active] / scale[active]
        def unpack(x):
            # Recover trial ion transfers [mol/s] from scaled solver unknowns.
            transfer = np.zeros(len(self.ions))
            transfer[active] = x[:-1] * scale[active]
            return transfer, x[-1] * self.max_length
        def residual(x):
            transfer, length = unpack(x)
            _, bulk = self._bulk_properties(transfer)
            surface, _, ratio = self._polarization(bulk)
            flux = self._flux(surface)
            # Include absent ions too: any nonzero predicted flux is infeasible.
            return np.r_[(transfer - self.width * length * flux) / scale,
                         (transfer[k] - target) / target]
        x0 = np.r_[0.1 * upper + 0.9 * lower, min(1., 1. / self.max_length)]
        # A strict inequality needs a numerical margin: i <= (1-1e-6)*i_lim.
        # Recompute i_lim from each trial outlet, not the inlet concentration.
        current_margin = 1e-6
        def current_slack(x):
            transfer, _ = unpack(x)
            _, bulk = self._bulk_properties(transfer)
            _, i_lim, _ = self._polarization(bulk)
            return (1. - current_margin) * i_lim / self.current_density - 1.
        def objective(x):
            errors = residual(x)
            return 0.5 * np.dot(errors, errors)
        fit = minimize(objective, x0, method='SLSQP', jac='3-point',
                       bounds=Bounds(np.r_[lower, 0.], np.r_[upper, 1.]),
                       constraints=[NonlinearConstraint(current_slack, 0., np.inf)],
                       options={'ftol': 1e-22, 'maxiter': 3000})
        transfer, length = unpack(fit.x)
        streams, bulk = self._bulk_properties(transfer)
        surface, i_lim, ratio = self._polarization(bulk)
        # Independently verify both the balances and the hard current constraint.
        if not fit.success or np.max(np.abs(residual(fit.x))) > 1e-7 or current_slack(fit.x) < -1e-9 or self.current_density >= i_lim:
            raise RuntimeError('ED sizing infeasible or unconverged: check recovery, limiting current, back diffusion and max_length; '
                               f'max scaled residual={np.max(np.abs(residual(fit.x))):.3g}, i/i_lim={ratio:.4g}.')
        if abs(transfer @ self.z) > 1e-8 * max(np.abs(transfer * self.z).sum(), 1e-12):
            raise ValueError('Transport parameters cause net charge transfer. Supply charge-consistent coupled permeabilities/transport numbers.')
        resistance = self._resistance(bulk)
        area = self.width * length
        voltage = self.current_density * resistance
        current = self.current_density * area
        for out, trial in zip(self.outs, streams):
            out.copy_like(trial)
        self.transfers = dict(zip(self.ions, transfer))  # mol/s, signed D -> C
        self.solver_residual = float(np.max(np.abs(residual(fit.x))))
        # Voltage [V] * current [A] gives W; divide by 1000 W/kW for utility power.
        self.design_results.update({'Length': length, 'Projected area': area,
            'Installed membrane area': 2 * area, 'Current': current,
            'Voltage': voltage, 'Areal resistance': resistance,
            'Limiting current density': i_lim, 'Cell-pair power': voltage * current / 1000.})

    def _design(self):
        self.add_power_utility(self.design_results['Cell-pair power'])
