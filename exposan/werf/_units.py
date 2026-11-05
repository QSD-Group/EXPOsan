# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>
    Zixuan Wang <wyatt4428@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, pandas as pd
from warnings import warn
from qsdsan import CompiledProcesses, SanUnit
from qsdsan.unit_operations import Splitter

__all__ = ('SelectiveRecovery',
           'TricklingFilter')

class SelectiveRecovery(Splitter):
    '''
    Selective removal and recovery of certain wastewater constituents (e.g., Ammonia). 
    Specify recovery efficiency with `split` between two outlet streams: 
    [0] recovered constituents and [1] wastewater.
    
    '''
    def _init_state(self):
        self._state = self._ins_QC[0]
        self._dstate = self._state * 0.

    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Splitter'''
        arr = self._state
        Cs = arr[:-1]
        Q = arr[-1]
        sp = self.split
        recv, eff = self.outs
        if recv.state is None: recv.state = 0 * arr
        if eff.state is None: eff.state = 1 * arr
        recv.state[-1] = 1 # assume recovered in a liquid solution of 1 m3/d
        recv.state[:-1] = sp * Cs * Q  # mass flows in g/d, equivalent to g/m3 concentrations in 1m3/d solution
        eff.state[:-1] = (1-sp) * Cs
        eff.state[-1] = Q

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Splitter'''
        pass

class TricklingFilter(SanUnit):
    # '''
    # Trickling filter scaffold discretized as horizontal sections in series.
    #
    # This first implementation is a bulk-liquid surrogate for the GPS-X-style
    # trickling filter: each horizontal section behaves as a CSTR running the
    # supplied suspended-growth model. Biofilm geometry and transport parameters
    # are stored as user-adjustable attributes for the next implementation stage.
    # '''

    '''
    Trickling filter discretized as horizontal sections in series.

    Each horizontal section includes one liquid-film CSTR and a customizable
    number of biofilm layers. Biological reactions are evaluated only in the
    biofilm layers through the supplied mASM2d-compatible process model, while
    soluble components move by liquid-film and biofilm diffusion and particulate
    components move by attachment, detachment, and internal solids exchange.
    Set ``N_biofilm_layer=0`` to recover the liquid-only CSTR-chain behavior.
    '''
    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = True

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', N_section=6, N_biofilm_layer=5,
                 V_sections=None, influent_fractions=None,
                 recirculation_ratio=0.0,
                 suspended_growth_model=None, DO_ID='S_O2',
                 DO_setpoints=(), kLa=(), DO_sat=8.0,
                 air_O2_transfer_coefficient=None,
                 filter_bed_depth=2.0, filter_bed_surface=200.0,
                 media_specific_surface=100.0,
                 liquid_retention_time=10/1440,
                 max_liquid_film_thickness=5e-5,
                 max_biofilm_thickness=1e-3,
                 initial_biofilm_thickness=None,
                 biofilm_density=1020000.0,
                 biofilm_dry_content=0.1,
                 f_diff=0.5, diffusivity=None,
                 attachment_rate=0.5,
                 detachment_rate=0.07,
                 detachment_exponent=1.0,
                 detachment_model='linear',
                 internal_solids_exchange_rate=2e-5,
                 dynamic_biofilm_thickness=True,
                 isdynamic=True, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                         isdynamic=isdynamic, **kwargs)
        self.N_section = N_section
        self.N_biofilm_layer = N_biofilm_layer
        self.V_sections = V_sections
        self.influent_fractions = influent_fractions
        self.recirculation_ratio = recirculation_ratio
        self.DO_ID = DO_ID
        self.DO_setpoints = DO_setpoints
        self.kLa = kLa
        self.DO_sat = DO_sat
        self.air_O2_transfer_coefficient = air_O2_transfer_coefficient
        self.filter_bed_depth = filter_bed_depth
        self.filter_bed_surface = filter_bed_surface
        self.media_specific_surface = media_specific_surface
        self.liquid_retention_time = liquid_retention_time
        self.max_liquid_film_thickness = max_liquid_film_thickness
        self.max_biofilm_thickness = max_biofilm_thickness
        self.initial_biofilm_thickness = initial_biofilm_thickness
        self.biofilm_density = biofilm_density
        self.biofilm_dry_content = biofilm_dry_content
        self.f_diff = f_diff
        self.diffusivity = diffusivity
        self.attachment_rate = attachment_rate
        self.detachment_rate = detachment_rate
        self.detachment_exponent = detachment_exponent
        self.detachment_model = detachment_model
        self.internal_solids_exchange_rate = internal_solids_exchange_rate
        self.dynamic_biofilm_thickness = dynamic_biofilm_thickness
        self.suspended_growth_model = suspended_growth_model
        self._concs = None
        self._Qs = np.zeros(self.N_section)

    @property
    def N_section(self):
        '''[int] Number of horizontal sections modeled as CSTRs in series.'''
        return self._N_section
    @N_section.setter
    def N_section(self, n):
        n = int(n)
        if n <= 0:
            raise ValueError(f'N_section must be positive, not {n}')
        self._N_section = n
        if hasattr(self, '_Qs') and self._Qs.size != n:
            self._Qs = np.zeros(n)
        if hasattr(self, '_ncol'):
            self._Qs_idx = list(self._ncol * np.arange(1, 1+n) - 1)
            self._refresh_state_header()
        self._ODE = None

    @property
    def N_biofilm_layer(self):
        '''[int] Number of biofilm layers; zero keeps the liquid-only CSTR chain.'''
        return self._N_biofilm_layer
    @N_biofilm_layer.setter
    def N_biofilm_layer(self, n):
        n = int(n)
        if n < 0:
            raise ValueError(f'N_biofilm_layer cannot be negative, not {n}')
        self._N_biofilm_layer = n
        if hasattr(self, '_ncol'):
            self._refresh_state_header()
        self._ODE = None

    @property
    def V_sections(self):
        '''[numpy.ndarray] Liquid volume of each horizontal section, in m3.'''
        return self._V_sections
    @V_sections.setter
    def V_sections(self, Vs):
        if Vs is None:
            self._V_sections = None
            return
        if np.isscalar(Vs):
            Vs = [float(Vs)] * self.N_section
        elif len(Vs) != self.N_section:
            raise ValueError(f'V_sections must have {self.N_section} values, '
                             f'not {len(Vs)}')
        self._V_sections = np.asarray(Vs, dtype=float)
        if (self._V_sections <= 0).any():
            raise ValueError('all V_sections values must be positive')
        self._ODE = None

    @property
    def influent_fractions(self):
        '''[numpy.ndarray] Fraction of each influent entering each section.'''
        return self._f_ins
    @influent_fractions.setter
    def influent_fractions(self, fs):
        if fs is None:
            fs = np.zeros((len(self.ins), self.N_section))
            fs[:, 0] = 1.0
        fs = np.asarray(fs, dtype=float)
        if fs.ndim == 1:
            fs = fs.reshape(1, -1)
        if fs.shape != (len(self.ins), self.N_section):
            raise ValueError('influent_fractions must have shape '
                             f'({len(self.ins)}, {self.N_section}), not {fs.shape}')
        if (fs < 0).any():
            raise ValueError('influent_fractions cannot contain negative values')
        row_sums = fs.sum(axis=1)
        if (row_sums == 0).any():
            raise ValueError('each influent must enter at least one section')
        self._f_ins = fs / row_sums[:, None]
        self._ODE = None

    @property
    def recirculation_ratio(self):
        '''[float] Recirculation ratio, defined as recycle flow divided by influent flow.'''
        return self._recirculation_ratio
    @recirculation_ratio.setter
    def recirculation_ratio(self, value):
        value = float(value)
        if value < 0:
            raise ValueError('recirculation_ratio cannot be negative')
        self._recirculation_ratio = value
        self._reset_auto_filter_bed_depth()
        self._ODE = None

    @property
    def DO_setpoints(self):
        '''[numpy.ndarray] Dissolved oxygen setpoint in each section, in mg/L.'''
        return self._DOs
    @DO_setpoints.setter
    def DO_setpoints(self, DOs):
        if DOs is None or len(DOs) == 0:
            self._DOs = np.zeros(self.N_section)
        elif len(DOs) == self.N_section:
            self._DOs = np.asarray(DOs, dtype=float)
        else:
            raise ValueError(f'DO_setpoints must be empty or have '
                             f'{self.N_section} values')
        self._ODE = None

    @property
    def kLa(self):
        '''[numpy.ndarray] Oxygen transfer coefficient in each section, in 1/d.'''
        return self._kLa
    @kLa.setter
    def kLa(self, ks):
        if ks is None or len(ks) == 0:
            self._kLa = np.zeros(self.N_section)
        elif len(ks) == self.N_section:
            self._kLa = np.asarray(ks, dtype=float)
        else:
            raise ValueError(f'kLa must be empty or have {self.N_section} values')
        self._ODE = None

    @property
    def DO_sat(self):
        '''[float] Saturation dissolved oxygen concentration, in mg/L.'''
        return self._DO_sat
    @DO_sat.setter
    def DO_sat(self, value):
        self._DO_sat = float(value)
        self._ODE = None

    @property
    def air_O2_transfer_coefficient(self):
        '''[numpy.ndarray] Passive air-to-liquid-film O2 transfer coefficient in each section, in m/d.'''
        if self._k_O2_air is None:
            D = self.diffusivity[self.components.index(self.DO_ID)]
            l_liq = max(self.max_liquid_film_thickness, 1e-12)
            k = D / l_liq
            return np.full(self.N_section, k)
        return self._k_O2_air
    @air_O2_transfer_coefficient.setter
    def air_O2_transfer_coefficient(self, ks):
        if ks is None:
            self._k_O2_air = None
        elif np.isscalar(ks):
            k = float(ks)
            if k < 0:
                raise ValueError('air_O2_transfer_coefficient cannot be negative')
            self._k_O2_air = np.full(self.N_section, k)
        elif len(ks) == self.N_section:
            self._k_O2_air = np.asarray(ks, dtype=float)
            if (self._k_O2_air < 0).any():
                raise ValueError('air_O2_transfer_coefficient cannot contain negative values')
        else:
            raise ValueError(f'air_O2_transfer_coefficient must be a scalar, None, '
                             f'or have {self.N_section} values')
        self._ODE = None

    @property
    def filter_bed_depth(self):
        '''[float] Filter bed depth, in m.'''
        return self._filter_bed_depth
    @filter_bed_depth.setter
    def filter_bed_depth(self, value):
        if value is None:
            raise ValueError('filter_bed_depth must be positive')
        value = float(value)
        if value <= 0:
            raise ValueError('filter_bed_depth must be positive')
        self._filter_bed_depth = value
        self._ODE = None

    @property
    def filter_bed_surface(self):
        '''[float] Filter bed surface area, in m2.'''
        return self._filter_bed_surface
    @filter_bed_surface.setter
    def filter_bed_surface(self, value):
        value = float(value)
        if value <= 0:
            raise ValueError('filter_bed_surface must be positive')
        self._filter_bed_surface = value
        self._ODE = None

    @property
    def hydraulic_loading_rate(self):
        '''[float] Applied hydraulic loading rate on plan filter-bed area, in m3/(m2*d).'''
        Q_m3_d = self._get_applied_flow_rate()
        if Q_m3_d > 0:
            A_plan_m2 = self.filter_bed_surface
            if A_plan_m2 > 0:
                return Q_m3_d / A_plan_m2
        return 0.

    def _get_applied_flow_rate(self):
        Q_external = 0.
        for ws in self.ins:
            state = getattr(ws, 'state', None)
            if state is not None and state[-1] > 0:
                Q = state[-1]
            else:
                Q = ws.F_vol * 24
            Q_external += Q
        if Q_external <= 0:
            return 0.
        R = getattr(self, '_recirculation_ratio', 0.)
        return Q_external * (1. + R)

    @property
    def media_specific_surface(self):
        '''[float] Specific surface area of filter media, in 1/m.'''
        return self._media_specific_surface
    @media_specific_surface.setter
    def media_specific_surface(self, value):
        value = float(value)
        if value <= 0:
            raise ValueError('media_specific_surface must be positive')
        self._media_specific_surface = value
        self._ODE = None

    def _reset_auto_filter_bed_depth(self):
        return

    def _calculate_filter_bed_depth(self, Q_m3_d):
        warn('filter_bed_depth is not calculated from hydraulic_loading_rate; '
             'hydraulic_loading_rate is a derived value based on filter_bed_surface plan area.')
        return self.filter_bed_depth

    def _update_filter_bed_depth_from_flow(self, Q_m3_d):
        return

    @property
    def liquid_retention_time(self):
        '''[float] Liquid retention time in the filter, in d.'''
        return self._liquid_retention_time
    @liquid_retention_time.setter
    def liquid_retention_time(self, value):
        self._liquid_retention_time = float(value)
        self._ODE = None

    @property
    def max_liquid_film_thickness(self):
        '''[float] Maximum attached liquid film thickness, in m.'''
        return self._max_liquid_film_thickness
    @max_liquid_film_thickness.setter
    def max_liquid_film_thickness(self, value):
        self._max_liquid_film_thickness = float(value)
        self._ODE = None

    @property
    def max_biofilm_thickness(self):
        '''[float] Maximum biofilm thickness, in m.'''
        return self._max_biofilm_thickness
    @max_biofilm_thickness.setter
    def max_biofilm_thickness(self, value):
        self._max_biofilm_thickness = float(value)
        self._ODE = None

    @property
    def initial_biofilm_thickness(self):
        '''[numpy.ndarray] Initial biofilm thickness by section and layer, in m.'''
        if self._initial_biofilm_thickness is None:
            if not self.N_biofilm_layer:
                return np.zeros((self.N_section, 0))
            return np.full((self.N_section, self.N_biofilm_layer),
                           self.max_biofilm_thickness/self.N_biofilm_layer)
        return self._initial_biofilm_thickness
    @initial_biofilm_thickness.setter
    def initial_biofilm_thickness(self, value):
        if value is None:
            self._initial_biofilm_thickness = None
            self._ODE = None
            return
        N = self.N_section
        L = self.N_biofilm_layer
        if not L:
            raise ValueError('initial_biofilm_thickness can only be set when N_biofilm_layer > 0')
        layer_max = self.max_biofilm_thickness / L
        if np.isscalar(value):
            value = np.full(N, float(value))
        value = np.asarray(value, dtype=float)
        if value.shape == (N,):
            totals = value
            value = np.zeros((N, L))
            for i, total in enumerate(totals):
                remaining = total
                for j in range(L):
                    filled = min(max(remaining, 0.), layer_max)
                    value[i, j] = filled
                    remaining -= filled
        elif value.shape != (N, L):
            raise ValueError('initial_biofilm_thickness must be a scalar, '
                             f'have shape ({N},), or have shape ({N}, {L}), '
                             f'not {value.shape}')
        if (value < 0).any():
            raise ValueError('initial_biofilm_thickness cannot contain negative values')
        if (value > layer_max).any():
            raise ValueError('layer initial_biofilm_thickness cannot exceed layer capacity')
        self._initial_biofilm_thickness = value
        self._ODE = None

    @property
    def biofilm_density(self):
        '''[float] Density of biofilm, in mg/L.'''
        return self._biofilm_density
    @biofilm_density.setter
    def biofilm_density(self, value):
        self._biofilm_density = float(value)
        self._ODE = None

    @property
    def biofilm_dry_content(self):
        '''[float] Dry material content of biofilm, dimensionless.'''
        return self._biofilm_dry_content
    @biofilm_dry_content.setter
    def biofilm_dry_content(self, value):
        self._biofilm_dry_content = float(value)
        self._ODE = None

    @property
    def f_diff(self):
        '''[float] Diffusion reduction factor in biofilm, dimensionless.'''
        return self._f_diff
    @f_diff.setter
    def f_diff(self, value):
        self._f_diff = float(value)
        self._ODE = None

    @property
    def diffusivity(self):
        '''[numpy.ndarray] Diffusivities of components in water, in m2/d.'''
        return self._diffusivity
    @diffusivity.setter
    def diffusivity(self, values):
        cmps = self.components
        n_cmps = len(cmps)
        D = np.zeros(n_cmps)
        if values is None:
            for ID, value in self._default_diffusivities.items():
                if ID in cmps.IDs:
                    D[cmps.index(ID)] = value
        elif np.isscalar(values):
            for i, cmp in enumerate(cmps):
                if getattr(cmp, 'particle_size', '') == 'Soluble':
                    D[i] = float(values)
        elif isinstance(values, dict):
            self.diffusivity = None
            D = self._diffusivity.copy()
            idx = cmps.indices(values.keys())
            D[idx] = list(values.values())
        else:
            D = np.asarray(values, dtype=float)
            if D.shape != (n_cmps,):
                raise ValueError(f'diffusivity must have {n_cmps} values, '
                                 f'not shape {D.shape}')
        self._diffusivity = D
        self._ODE = None

    _cm2_s_to_m2_d = 8.64
    _default_diffusivities = {
        'S_O2': 2.5e-5 * _cm2_s_to_m2_d,
        'S_h2': 5.84e-5 * _cm2_s_to_m2_d,
        'S_N2': 1.9e-5 * _cm2_s_to_m2_d,
        'S_ch4': 1.0e-5 * _cm2_s_to_m2_d,
        'S_I': 1.0e-5 * _cm2_s_to_m2_d,
        'S_F': 6.9e-6 * _cm2_s_to_m2_d,
        'S_A': 1.24e-5 * _cm2_s_to_m2_d,
        'S_NH4': 2.0e-5 * _cm2_s_to_m2_d,
        'S_NO2': 1.23e-5 * _cm2_s_to_m2_d,
        'S_NO3': 1.23e-5 * _cm2_s_to_m2_d,
        'S_PO4': 1.0e-5 * _cm2_s_to_m2_d,
        'S_IC': 1.96e-5 * _cm2_s_to_m2_d,
        'S_ALK': 1.96e-5 * _cm2_s_to_m2_d,
        'S_Ca': 1.0e-5 * _cm2_s_to_m2_d,
        'S_K': 1.0e-5 * _cm2_s_to_m2_d,
        'S_Mg': 1.0e-5 * _cm2_s_to_m2_d,
        'S_Na': 1.0e-5 * _cm2_s_to_m2_d,
        'S_Cl': 3.1e-5 * _cm2_s_to_m2_d,
        }

    @property
    def attachment_rate(self):
        '''[float] Solids attachment rate to biofilm, in m/d.'''
        return self._attachment_rate
    @attachment_rate.setter
    def attachment_rate(self, value):
        self._attachment_rate = float(value)
        self._ODE = None

    @property
    def detachment_rate(self):
        '''[float] Detachment scale in kg dry solids/(m2*d), at the configured
        dry-solids density. 
        The capacity law reaches it at Hmax/2 and increases toward Hmax.
        Both laws scale by actual outer-layer dry density / configured density.
        Net excess solids at the maximum-thickness constraint also slough.
        '''
        return self._detachment_rate
    @detachment_rate.setter
    def detachment_rate(self, value):
        self._detachment_rate = float(value)
        self._ODE = None

    @property
    def detachment_exponent(self):
        '''[float] Dimensionless exponent applied to relative thickness in
        detachment; default 1.0, finite and at least one.'''
        return self._detachment_exponent
    @detachment_exponent.setter
    def detachment_exponent(self, value):
        value = float(value)
        if not np.isfinite(value) or value < 1.:
            raise ValueError('detachment_exponent must be finite and at least one')
        self._detachment_exponent = value
        self._ODE = None

    @property
    def detachment_model(self):
        '''[str] Thickness dependence of detachment: ``linear`` (legacy) or
        ``capacity``. With r=H/Hmax and n=detachment_exponent, use r**n or
        (r/(1-r))**n, respectively. The capacity law approximates the GPS-X
        description of increasing sloughing approaching maximum thickness;
        it is not a published GPS-X equation. A 0.001 minimum free fraction
        regularizes solver trial states at the cap. At n=1 and half-full
        thickness, capacity-law flux equals detachment_rate at design density.
        '''
        return self._detachment_model

    @detachment_model.setter
    def detachment_model(self, value):
        if value not in ('linear', 'capacity'):
            raise ValueError("detachment_model must be 'linear' or 'capacity'")
        self._detachment_model = value
        self._ODE = None

    @property
    def internal_solids_exchange_rate(self):
        '''[float] Internal solids exchange rate between biofilm layers, in m/d.'''
        return self._internal_solids_exchange_rate
    @internal_solids_exchange_rate.setter
    def internal_solids_exchange_rate(self, value):
        self._internal_solids_exchange_rate = float(value)
        self._ODE = None

    @property
    def dynamic_biofilm_thickness(self):

        '''[bool] Evolve one total biofilm thickness per section (default True).

        Equal-width layers use -C*(dV/dt)/V to conserve component inventory.
        Conservative upwind particulate transfers redistribute local solids
        production among layers, preserving initially uniform dry-solids density
        while allowing component fractions to differ. Solubles retain diffusion.
        At maximum thickness, net excess dry solids are sloughed to the liquid
        using the outer-layer composition. Numerical volume floors can change
        mean density. Custom unequal initial densities are not renormalized.
        Fixed-thickness mode omits this moving-layer transport and overflow.
        '''
        return self._dynamic_biofilm_thickness
    @dynamic_biofilm_thickness.setter
    def dynamic_biofilm_thickness(self, value):
        self._dynamic_biofilm_thickness = bool(value)
        self._ODE = None

    @property
    def suspended_growth_model(self):
        '''[:class:`CompiledProcesses`] Suspended-growth process model.'''
        return self._model
    @suspended_growth_model.setter
    def suspended_growth_model(self, model):
        if not isinstance(model, CompiledProcesses) and model is not None:
            raise TypeError('suspended_growth_model must be a CompiledProcesses '
                            f'object or None, not {type(model)}')
        self._model = model
        self._state_keys = list(self.components.IDs) + ['Q']
        self._ncol = len(self._state_keys)
        self._Qs_idx = list(self._ncol * np.arange(1, 1+self.N_section) - 1)
        self._refresh_state_header()
        self._ODE = None

    def _refresh_state_header(self):
        cmp_names = [f'{i} [mg/L]' for i in self.components.IDs]
        liq_names = cmp_names + ['Q [m3/d]']
        header = []
        for i in range(self.N_section):
            header.extend(f'section {i} liquid {name}' for name in liq_names)
        if self.N_biofilm_layer:
            for i in range(self.N_section):
                for j in range(self.N_biofilm_layer):
                    header.extend(f'section {i} biofilm {j} {name}'
                                  for name in cmp_names)
            for i in range(self.N_section):
                header.append(f'section {i} biofilm thickness [m]')
        self._state_header = header

    @property
    def biofilm_thickness(self):
        '''[numpy.ndarray] Current biofilm thickness by section and layer, in m.'''
        if self._state is None:
            return self.initial_biofilm_thickness.copy()
        if not self.N_biofilm_layer:
            return np.zeros(self.N_section)
        N = self.N_section
        L = self.N_biofilm_layer
        n_cmps = len(self.components)
        start = N*self._ncol + N*L*n_cmps
        H_total = np.clip(self._state[start:start+N], 0.,
                          self.max_biofilm_thickness)
        return np.repeat((H_total/L)[:, None], L, axis=1)

    @property
    def total_biofilm_thickness(self):
        '''[numpy.ndarray] Total biofilm thickness in each section, in m.'''
        if not self.N_biofilm_layer:
            return np.zeros(self.N_section)
        if self._state is None:
            return self.initial_biofilm_thickness.sum(axis=1)
        N = self.N_section
        L = self.N_biofilm_layer
        n_cmps = len(self.components)
        start = N*self._ncol + N*L*n_cmps
        return np.clip(self._state[start:start+N], 0., self.max_biofilm_thickness)

    @property
    def DO_ID(self):
        '''[str] Component ID used as dissolved oxygen in the process model.'''
        return self._DO_ID
    @DO_ID.setter
    def DO_ID(self, ID):
        if ID not in self.components.IDs:
            raise ValueError(f'DO_ID must be one of {self.components.IDs}, not {ID}')
        self._DO_ID = ID
        self._ODE = None

    @property
    def state(self):
        '''[pandas.DataFrame] Liquid-film and biofilm states for each section.'''
        if self._state is None:
            return None
        n_cmps = len(self.components)
        N = self.N_section
        L = self.N_biofilm_layer
        y_liq = self._state[:N*self._ncol].reshape((N, self._ncol))
        records = [*y_liq]
        index = [f'section {i} liquid' for i in range(N)]
        if L:
            start = N*self._ncol
            stop = start + N*L*n_cmps
            y_bio = self._state[start:stop].reshape((N, L, n_cmps))
            records.extend(np.insert(y_bio.reshape((N*L, n_cmps)), n_cmps, np.nan, axis=1))
            index.extend(f'section {i} biofilm {j}' for i in range(N) for j in range(L))
        return pd.DataFrame(records, index=index, columns=self._state_keys)

    def get_retained_mass(self, biomass_IDs):
        """Return selected biomass inventory in g dry solids for SRT reporting.

        Includes liquid and attached biomass using component i_mass factors.
        Biofilm volumes follow the same clipped thickness and numerical floors
        as the concentration balances; this method does not modify any state.
        The dynamic state and liquid section volumes must be initialized.
        """
        if self._state is None or self.V_sections is None:
            raise RuntimeError('Initialize the TF dynamic state before reporting biomass.')
        cmps = self.components
        indices = cmps.indices(biomass_IDs)
        factors = cmps.i_mass[indices]
        N, L, n = self.N_section, self.N_biofilm_layer, len(cmps)
        liquid_stop = N*self._ncol
        liquid = self._state[:liquid_stop].reshape(N, self._ncol)
        retained = (liquid[:, indices] @ factors) @ self.V_sections
        if L:
            biofilm = self._state[liquid_stop:liquid_stop+N*L*n].reshape(N, L, n)
            area = max(self.filter_bed_surface*self.filter_bed_depth
                       * self.media_specific_surface/N, 1e-12)
            minimum = max(self.max_biofilm_thickness/L*1e-6, 1e-12)
            volumes = np.maximum(area*np.maximum(self.biofilm_thickness, minimum),
                                 1e-12)
            retained += np.sum((biofilm[:, :, indices] @ factors)*volumes)
        return float(retained)

    def set_init_conc(self, concentrations=None, i_section=None, **kwargs):
        """Set liquid and biofilm concentrations without changing their values.

        A two-level DataFrame index describes zero-based (section, layer),
        with layer -1 reserved for liquid and layers 0, 1, ... for biofilm.
        All components must be provided. Profiles on a different uniform grid
        are conservatively averaged over overlapping section/depth cells;
        identical grids load directly. Concentrations are not renormalized to
        the configured dry density. Other DataFrames retain liquid-profile
        broadcasting, and arrays retain the existing shape conventions.
        """
        cmps = self.components
        n_cmps = len(cmps)
        N = self.N_section
        L = self.N_biofilm_layer
        shape = (N, n_cmps) if not L else (N, L+1, n_cmps)
        if self._concs is None:
            self._concs = np.zeros(shape)
        if concentrations is None:
            concs = cmps.kwarray(kwargs)
            if i_section is None:
                self._concs[...] = concs
            else:
                self._concs[i_section, ...] = concs
        elif isinstance(concentrations, pd.DataFrame) and isinstance(concentrations.index, pd.MultiIndex):
            if concentrations.index.nlevels != 2 or not concentrations.index.is_unique:
                raise ValueError('TF profiles require unique (section, layer) rows.')
            sections = concentrations.index.get_level_values(0).unique().sort_values()
            layers = concentrations.index.get_level_values(1).unique().sort_values()
            source_N, source_L = len(sections), len(layers)-1
            if not source_N or source_L < 0:
                raise ValueError('TF initial profiles cannot be empty.')
            if list(sections) != list(range(source_N)) or list(layers) != list(range(-1, source_L)):
                raise ValueError('TF section/layer indices must be contiguous from 0/-1.')
            index = pd.MultiIndex.from_product((sections, layers))
            if len(index) != len(concentrations):
                raise ValueError('Every TF section must include the same layer rows.')
            values = concentrations.loc[index, list(cmps.IDs)].to_numpy(dtype=float)
            if not np.isfinite(values).all() or (values < 0.).any():
                raise ValueError('TF initial concentrations must be finite and nonnegative.')
            source = values.reshape(source_N, source_L+1, n_cmps)

            def overlap_weights(old, new):
                if old == new:
                    return np.eye(old)
                source_edges = np.linspace(0., 1., old+1)
                target_edges = np.linspace(0., 1., new+1)
                return np.maximum(0., np.minimum(target_edges[1:, None], source_edges[None, 1:])
                                  - np.maximum(target_edges[:-1, None], source_edges[None, :-1]))*new

            section_profile = np.einsum('ij,jkc->ikc', overlap_weights(source_N, N), source)
            if L:
                if not source_L:
                    raise ValueError('The workbook profile contains no biofilm layers.')
                self._concs = np.empty(shape)
                self._concs[:, 0, :] = section_profile[:, 0, :]
                self._concs[:, 1:, :] = np.einsum(
                    'ij,njc->nic', overlap_weights(source_L, L), section_profile[:, 1:, :])
            else:
                self._concs = section_profile[:, 0, :].copy()
        elif isinstance(concentrations, pd.DataFrame):
            concentrations.index = range(self.N_section)
            for i, concs in concentrations.to_dict('index').items():
                self._concs[i, ...] = cmps.kwarray(concs)
        else:
            concentrations = np.asarray(concentrations, dtype=float)
            if L and concentrations.shape == (N, n_cmps):
                self._concs = np.zeros(shape)
                self._concs[...] = concentrations[:, None, :]
            elif L and concentrations.shape == (N*(L+1), n_cmps):
                self._concs = concentrations.reshape(shape)
            elif concentrations.shape == shape:
                self._concs = concentrations
            else:
                raise ValueError(f'concentrations must have shape {shape}, '
                                 f'not {concentrations.shape}')

    def _run(self):
        out, = self.outs
        out.mix_from(self.ins)

    def _init_state(self):
        out, = self.outs
        y = np.zeros((self.N_section, self._ncol))
        if self._concs is not None:
            if self.N_biofilm_layer:
                y[:, :-1] = self._concs[:, 0, :]
            else:
                y[:, :-1] = self._concs
        else:
            y[:, :-1] = out.conc
        Q_ins = np.asarray([ws.F_vol*24 for ws in self.ins])
        Q_external = sum(Q_ins)
        Q_recycle = self.recirculation_ratio * Q_external
        self._Qs[:] = np.dot(Q_ins, self.influent_fractions).cumsum() + Q_recycle
        if self.V_sections is None:
            V_total = max(Q_external + Q_recycle, 1e-12) * self.liquid_retention_time
            self._V_sections = np.full(self.N_section, V_total/self.N_section)
        y[:, -1] = self._Qs

        if self.N_biofilm_layer:
            if self._concs is not None:
                y_bio = self._concs[:, 1:, :]
            else:
                y_bio = np.tile(out.conc, (self.N_section, self.N_biofilm_layer, 1))
            H_total = self.initial_biofilm_thickness.sum(axis=1)
            self._state = np.append(np.append(y.flatten(), y_bio.flatten()),
                                    H_total)
        else:
            self._state = y.flatten()
        self._dstate = self._state * 0.

    def _update_state(self):
        out, = self.outs
        if out.state is None:
            out.state = np.zeros(self._ncol)
        concentrations = self._state[stop-self._ncol:stop-1]
        out.state[:-1] = np.where(concentrations < 1e-16, 0., concentrations)
        DO = self.DO_setpoints
        if len(DO) and DO[-1] > 0.:
            out.state[self.components.index(self.DO_ID)] = DO[-1]
        out.state[-1] = sum(ws.state[-1] for ws in self.ins)

    def _update_dstate(self):
        out, = self.outs
        if out.dstate is None:
            out.dstate = np.zeros(self._ncol)
        stop = self.N_section * self._ncol
        out.dstate[:-1] = self._dstate[stop-self._ncol:stop-1]
        out.dstate[-1] = sum(ws.dstate[-1] for ws in self.ins)

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE

    def _compile_ODE(self):
        _Qs = self._Qs
        _dstate = self._dstate
        _update_dstate = self._update_dstate
        cmps = self.components
        N = self.N_section
        L = self.N_biofilm_layer
        ncol = self._ncol
        n_cmps = len(cmps)
        _1_ov_V = np.diag(1/self.V_sections)
        f_in = self.influent_fractions
        DO = self.DO_setpoints
        kLa = self.kLa
        DO_idx = cmps.index(self.DO_ID)
        DOsat = self.DO_sat
        k_O2_air = self.air_O2_transfer_coefficient
        R = self.recirculation_ratio
        dynamic_H = self.dynamic_biofilm_thickness

        if self._model is None:
            warn(f'{self.ID} was initialized without a suspended growth model, '
                 'and thus run as a non-reactive unit')
            def Rs(Cs):
                return np.zeros_like(Cs)
        else:
            f_rho = self._model.rate_function
            M_stoi = self._model.stoichio_eval()
            def Rs(Cs):
                rhos = np.apply_along_axis(f_rho, 1, np.maximum(Cs, 0.))
                return rhos @ M_stoi

        if not L:
            if any(DO):
                aerated = DO > 0
                aerated_DO = DO[aerated]
                def dy_dt(t, QC_ins, QC, dQC_ins):
                    y = QC.reshape((N, ncol))
                    Cs = y[:, :-1].copy()
                    Cs[aerated, DO_idx] = aerated_DO
                    Q_recycle = R * sum(QC_ins[:, -1])
                    _Qs[:] = Qs = np.dot(QC_ins[:, -1], f_in).cumsum() + Q_recycle
                    M_ins = f_in.T @ np.diag(QC_ins[:, -1]) @ QC_ins[:, :-1]
                    M_ins[0] += Q_recycle * Cs[-1]
                    M_outs = np.diag(Qs) @ Cs
                    M_ins[1:] += M_outs[:-1]
                    dy = np.zeros_like(y)
                    dy[:, :-1] = _1_ov_V @ (M_ins - M_outs) + Rs(Cs)
                    dy[aerated, DO_idx] = 0.
                    _dstate[:] = dy.flatten()
                    _update_dstate()
            else:
                def dy_dt(t, QC_ins, QC, dQC_ins):
                    y = QC.reshape((N, ncol))
                    Cs = y[:, :-1].copy()
                    Q_recycle = R * sum(QC_ins[:, -1])
                    _Qs[:] = Qs = np.dot(QC_ins[:, -1], f_in).cumsum() + Q_recycle
                    M_ins = f_in.T @ np.diag(QC_ins[:, -1]) @ QC_ins[:, :-1]
                    M_ins[0] += Q_recycle * Cs[-1]
                    M_outs = np.diag(Qs) @ Cs
                    M_ins[1:] += M_outs[:-1]
                    dy = np.zeros_like(y)
                    dy[:, :-1] = _1_ov_V @ (M_ins - M_outs) + Rs(Cs)
                    dy[:, DO_idx] += kLa * (DOsat - Cs[:, DO_idx])
                    _dstate[:] = dy.flatten()
                    _update_dstate()
            self._ODE = dy_dt
            return

        V_liq = self.V_sections
        A_total = self.filter_bed_surface * self.filter_bed_depth \
            * self.media_specific_surface
        A = max(A_total / N, 1e-12)
        h_layer_max = self.max_biofilm_thickness / L
        h_min = max(h_layer_max * 1e-6, 1e-12)
        h_atol = max(h_layer_max * 1e-9, 1e-15)
        l_liq = max(self.max_liquid_film_thickness, 1e-12)
        D_water = self.diffusivity
        D_bio = D_water * self.f_diff
        soluble_idx = np.flatnonzero(D_water > 0)
        particulate_idx = np.asarray([
            i for i, cmp in enumerate(cmps)
            if getattr(cmp, 'particle_size', '') == 'Particulate'
            ], dtype=int)
        k_liq = D_water[soluble_idx] / l_liq
        v_att = self.attachment_rate
        dry_density = max(self.biofilm_density * self.biofilm_dry_content * 1e-3,
                          1e-12)
        v_det = self.detachment_rate / dry_density
        v_int = self.internal_solids_exchange_rate
        liq_stop = N * ncol
        bio_stop = liq_stop + N * L * n_cmps
        dry_density_g_m3 = dry_density * 1e3
        solid_mass_factors = np.asarray(cmps.i_mass)[particulate_idx]
        solid_mass_factors = np.nan_to_num(solid_mass_factors, nan=0.)
        solid_mass_factors[solid_mass_factors < 0.] = 0.
        has_solid_mass = particulate_idx.size and solid_mass_factors.any()
        if any(DO):
            aerated = DO > 0
            aerated_DO = DO[aerated]
        else:
            aerated = None

        def dy_dt(t, QC_ins, QC, dQC_ins):
            y_liq = QC[:liq_stop].reshape((N, ncol))
            C_liq = y_liq[:, :-1].copy()
            C_bio = QC[liq_stop:bio_stop].reshape((N, L, n_cmps)).copy()
            H_total_raw = QC[bio_stop:bio_stop+N]
            H_total = np.clip(H_total_raw, 0., self.max_biofilm_thickness)
            H_below = H_total_raw < 0.
            H_above = H_total_raw > self.max_biofilm_thickness
            bound_recovery_rate = 1.0
            H = np.repeat((H_total/L)[:, None], L, axis=1)
            H_eff = np.maximum(H, h_min)
            V_bio = np.maximum(A * H_eff, 1e-12)
            dz_interface = np.maximum((H_eff[:, :-1] + H_eff[:, 1:])/2, h_min)
            k_bio_ov_dz = D_bio[soluble_idx] / dz_interface[:, :, None]
            A_ov_V_bio = A / V_bio
            if aerated is not None:
                C_liq[aerated, DO_idx] = aerated_DO
            Q_recycle = R * sum(QC_ins[:, -1])
            _Qs[:] = Qs = np.dot(QC_ins[:, -1], f_in).cumsum() + Q_recycle
            M_ins = f_in.T @ np.diag(QC_ins[:, -1]) @ QC_ins[:, :-1]
            M_ins[0] += Q_recycle * C_liq[-1]
            M_outs = np.diag(Qs) @ C_liq
            M_ins[1:] += M_outs[:-1]
            dy_liq = np.zeros_like(y_liq)
            dy_bio = np.zeros_like(C_bio)
            dy_liq[:, :-1] = _1_ov_V @ (M_ins - M_outs)

            if soluble_idx.size:
                J_liq_bio = k_liq * (C_liq[:, soluble_idx]
                                      - C_bio[:, 0, soluble_idx])
                dy_liq[:, soluble_idx] -= (A / V_liq[:, None]) * J_liq_bio
                dy_bio[:, 0, soluble_idx] += A_ov_V_bio[:, 0, None] * J_liq_bio
                if L > 1:
                    S_bio = C_bio[:, :, soluble_idx]
                    J_bio = k_bio_ov_dz * (S_bio[:, :-1, :] - S_bio[:, 1:, :])
                    dy_bio[:, :-1, soluble_idx] -= A_ov_V_bio[:, :-1, None] * J_bio
                    dy_bio[:, 1:, soluble_idx] += A_ov_V_bio[:, 1:, None] * J_bio

            if particulate_idx.size:
                attachment = v_att * C_liq[:, particulate_idx]
                attachment_layer = np.zeros((N, L, len(particulate_idx)))
                attachment_layer[:, 0, :] = attachment
                detachment_layer = np.zeros((N, L, len(particulate_idx)))
                relative_thickness = H_total / self.max_biofilm_thickness
                if self.detachment_model == 'capacity':
                    fullness = relative_thickness / np.maximum(1. - relative_thickness, 1e-3)
                else:
                    fullness = relative_thickness
                detachment_fraction = fullness ** self.detachment_exponent
                # Both liquid and biofilm balances use exactly the same flux.
                detachment_layer[:, 0, :] = v_det \
                    * C_bio[:, 0, particulate_idx] * detachment_fraction[:, None]
                dy_liq[:, particulate_idx] += A / V_liq[:, None] \
                    * np.sum(detachment_layer - attachment_layer, axis=1)
                dy_bio[:, :, particulate_idx] += A_ov_V_bio[:, :, None] \
                    * (attachment_layer - detachment_layer)
                if L > 1 and v_int:
                    exchange = v_int * C_bio[:, 1:, particulate_idx]
                    dy_bio[:, :-1, particulate_idx] += A_ov_V_bio[:, :-1, None] * exchange
                    dy_bio[:, 1:, particulate_idx] -= A_ov_V_bio[:, 1:, None] * exchange
            Rs_bio = Rs(C_bio.reshape((N*L, n_cmps))).reshape((N, L, n_cmps))
            dy_bio += Rs_bio
            dh_total = np.zeros(N)
            if has_solid_mass and dynamic_H:
                net_surface_solid = (attachment_layer - detachment_layer) \
                    @ solid_mass_factors
                net_rxn_solid = np.sum(Rs_bio[:, :, particulate_idx]
                                       * solid_mass_factors * H[:, :, None],
                                       axis=2)
                dh_total = np.sum(net_surface_solid + net_rxn_solid,
                                  axis=1) / dry_density_g_m3
                dh_total[(H_total >= self.max_biofilm_thickness - h_atol)
                         & (dh_total > 0.)] = 0.
                dh_total[(H_total <= h_min) & (dh_total < 0.)] = 0.
            if H_above.any():
                dh_total[H_above] = np.minimum(dh_total[H_above],
                    -bound_recovery_rate*(H_total_raw[H_above] - self.max_biofilm_thickness))
            if H_below.any():
                dh_total[H_below] = np.maximum(dh_total[H_below],
                    -bound_recovery_rate*H_total_raw[H_below])
            if dynamic_H:
                dh_layer = np.broadcast_to((dh_total/L)[:, None], H.shape).copy()
                dh_layer[H_below | H_above] = 0.
                dh_layer[(H < h_min)
                         | ((H == h_min) & (dh_layer < 0.))] = 0.
                dV_bio = A * dh_layer
                V_raw = A * H_eff
                dV_bio[(V_raw < 1e-12)
                       | ((V_raw == 1e-12) & (dV_bio < 0.))] = 0.
                dy_bio -= C_bio * (dV_bio / V_bio)[:, :, None]
                at_capacity = H_total >= self.max_biofilm_thickness - h_atol
                if has_solid_mass and at_capacity.any():
                    # A capped volume cannot store further dry solids at fixed
                    # density. Discharge its net excess through the outer surface.
                    particles = np.maximum(C_bio[:, :, particulate_idx], 0.)
                    density_rate = dy_bio[:, :, particulate_idx] @ solid_mass_factors
                    excess = np.maximum(np.sum(V_bio * density_rate, axis=1), 0.)
                    excess *= at_capacity
                    surface_density = particles[:, 0, :] @ solid_mass_factors
                    slough = np.divide(excess[:, None] * particles[:, 0, :],
                                       surface_density[:, None],
                                       out=np.zeros_like(particles[:, 0, :]),
                                       where=surface_density[:, None] > 0.)
                    dy_bio[:, 0, particulate_idx] -= slough / V_bio[:, 0, None]
                    dy_liq[:, particulate_idx] += slough / V_liq[:, None]
                if L > 1 and has_solid_mass:
                    # Moving equal-width cells need particulate advection as
                    # local growth, attachment and erosion change their inventory.
                    # Use the density derivative after the moving-volume term.
                    # Its volume-weighted mean is the closed-section constraint;
                    # interior face fluxes redistribute only the deviations.
                    particles = np.maximum(C_bio[:, :, particulate_idx], 0.)
                    layer_density = particles @ solid_mass_factors
                    density_rate = dy_bio[:, :, particulate_idx] @ solid_mass_factors
                    mean_rate = np.sum(V_bio * density_rate, axis=1) \
                        / np.sum(V_bio, axis=1)
                    imbalance = V_bio * (density_rate - mean_rate[:, None])
                    solid_flux = np.cumsum(imbalance, axis=1)[:, :-1]  # g dry/d
                    # Positive flux is outer-to-inner, relative to the moving grid.
                    # No added flux crosses the liquid surface or media boundary.
                    donors = np.where(solid_flux[:, :, None] >= 0.,
                                      particles[:, :-1, :], particles[:, 1:, :])
                    donor_density = np.where(solid_flux >= 0.,
                                             layer_density[:, :-1], layer_density[:, 1:])
                    fractions = np.divide(donors, donor_density[:, :, None],
                                          out=np.zeros_like(donors),
                                          where=donor_density[:, :, None] > 0.)
                    transfer = solid_flux[:, :, None] * fractions
                    dy_bio[:, :-1, particulate_idx] -= transfer / V_bio[:, :-1, None]
                    dy_bio[:, 1:, particulate_idx] += transfer / V_bio[:, 1:, None]
            if aerated is not None:
                dy_liq[aerated, DO_idx] = 0.
            else:
                K_O2_air = k_O2_air * A / V_liq
                dy_liq[:, DO_idx] += K_O2_air \
                    * (DOsat - C_liq[:, DO_idx])
            _dstate[:liq_stop] = dy_liq.flatten()
            _dstate[liq_stop:bio_stop] = dy_bio.flatten()
            _dstate[bio_stop:bio_stop+N] = dh_total
            _update_dstate()

        self._ODE = dy_dt

    def _design(self):
        pass
