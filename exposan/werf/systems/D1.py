# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''
# CODEX adds the followings 250 lines
# CODEX adds the followings 1 line
import numpy as np
from scipy.integrate import BDF, LSODA
import qsdsan as qs
from qsdsan import (
    process_models as pc,
    unit_operations as su,
    )
from qsdsan.utils import ospath, time_printer, load_data
from exposan.werf import data_path, default_fctss_init, default_rww
from exposan.werf._units import TricklingFilter

__all__ = ('create_d1_system',)

ID = 'D1'
#%%
dfs = load_data(
    ospath.join(data_path, 'initial_conditions.xlsx'),
    sheet=None,
    )
asinit = dfs[ID].set_index('Layer (-1=liquid)', append=True)
fcinit = dfs['E2P'].iloc[-1].to_dict()
adinit = dfs['adm'].loc[ID].to_dict()

MGD2cmd = 3785.412
Temp = 273.15+20 # temperature [K]
T_ad = 273.15+35

# %% TF biology: base mASM2d rates, matching the current SVCW TF profile
_D1_TF_ASM_PARAMETERS = dict(
    electron_acceptor_dependent_decay=True,
    mu_H=6., q_PHA=3., q_PP=1.5, b_PP=0.2,
    )


def _create_d1_tf_model():
    """Match SVCW's base TF rates without changing WERF component indices."""
    return pc.mASM2d(**_D1_TF_ASM_PARAMETERS)


# %% Shared TF numerical safeguards (also used by ESAP validation and SVCW)
class _ValidationBDF(BDF):
    """Bound numerical-Jacobian perturbations for inactive TF states."""

    def _validate_jac(self, jac, sparsity):
        evaluate, initial = super()._validate_jac(jac, sparsity)
        if jac is not None:
            return evaluate, initial

        def bounded_jacobian(t, y):
            # Zero-sensitivity columns can grow SciPy's persistent factors
            # without bound. Limit their relative perturbation to 0.1%;
            # this changes neither integration tolerances nor model rates.
            factors = self.jac_factor
            if factors is not None:
                invalid = ~np.isfinite(factors) | (factors <= 0.)
                factors[invalid] = np.sqrt(np.finfo(float).eps)
                np.minimum(factors, 1e-3, out=factors)
            return evaluate(t, y)

        return bounded_jacobian, initial


def _tf_tolerance_solver(base, system):
    """Apply TF-specific absolute tolerances after system state assembly."""
    if isinstance(base, str):
        base = {'BDF': _ValidationBDF, 'LSODA': LSODA}.get(base, base)
    if isinstance(base, str):
        return base

    class TFToleranceSolver(base):
        def __init__(solver, fun, t0, y0, t_bound, **options):
            tf = system.flowsheet.unit.TF
            start, stop = system._state_idx[tf.ID]
            atol = np.broadcast_to(np.asarray(options.get('atol', 1e-6), dtype=float), y0.shape).copy()
            local = atol[start:stop]
            N, L, n = tf.N_section, tf.N_biofilm_layer, len(tf.components)
            ids = [ID for ID in ('S_O2', 'S_NH4', 'S_NO3', 'S_PO4',
                                'X_PHA', 'X_PP', 'X_PAO')
                   if ID in tf.components.IDs]
            indices = tf.components.indices(ids)
            liquid_stop = N * tf._ncol
            liquid = local[:liquid_stop].reshape(N, tf._ncol)
            liquid[:, indices] = np.minimum(liquid[:, indices], 1e-6)
            if L:
                bio_stop = liquid_stop + N * L * n
                bio = local[liquid_stop:bio_stop].reshape(N, L, n)
                bio[:, :, indices] = np.minimum(bio[:, :, indices], 1e-6)
                local[bio_stop:bio_stop+N] = np.minimum(local[bio_stop:bio_stop+N], 1e-9)
            options['atol'] = atol
            system.tf_effective_atol = local.copy()
            super().__init__(fun, t0, y0.copy(), t_bound, **options)

    return TFToleranceSolver


class _D1System(qs.System):
    """WERF D1 with the TF-specific integration safeguards used by SVCW."""

    def dynamic_run(self, **kwargs):
        """Tighten TF state tolerances without changing other unit tolerances."""
        method = kwargs.get('method', self.dynsim_kwargs.get('method', 'BDF'))
        kwargs['method'] = _tf_tolerance_solver(method, self)
        try:
            return super().dynamic_run(**kwargs)
        finally:
            self.dynsim_kwargs['method'] = method


def create_d1_system(flowsheet=None, default_init_conds=True, tf_kwargs=None):
    """Create WERF D1 with SVCW-matched TF kinetics and targeted accuracy.

    Parameters
    ----------
    flowsheet : qsdsan.Flowsheet, optional
        Flowsheet to populate; a new D1 flowsheet is used by default.
    default_init_conds : bool, optional
        Read TF liquid/biofilm profiles from sheet D1 and AD from row D1 of
        sheet adm in initial_conditions.xlsx. FC uses F1's E2P-derived
        soluble/solids composition and default ten-layer TSS profile.
        Defaults to True.
    tf_kwargs : dict, optional
        TricklingFilter overrides. Physical defaults remain six sections,
        five biofilm layers, 9464 m2 bed area, recirculation ratio one,
        initial thickness 4e-4 m, and capacity detachment at 0.065 kg/(m2*d).
        A supplied suspended_growth_model overrides the default TF biology.

    Notes
    -----
    Retains F1's stock mASM2d components and default_rww feed: 10 MGD,
    358 mg COD/L, 25.91 mg NH4-N/L and 5 mg PO4-P/L at 20 degC.
    Only the biological rate profile matches SVCW; no site-specific component
    indices or influent composition are copied. The model uses base
    electron-acceptor-dependent decay, mu_H=6, q_PHA=3, q_PP=1.5 and b_PP=0.2,
    without the former additional anaerobic heterotrophic-decay term.

    The workbook stores the previous six-section/five-layer TF startup values
    explicitly, including the meaningful PAO/PHA/PP seed and dry-solids density.
    Different uniform grids use conservative profile averaging in set_init_conc.
    No seeding or density renormalization is performed while constructing D1;
    edit the workbook if changing its initial concentrations or dry density.

    FC uses F1's ten-layer FlatBottomCircularClarifier, geometry, settling
    parameters and 0.163 MGD wastage. Its RAS underflow is zero: attached-growth
    D1 does not return sludge to the TF. All settled sludge goes through WAS
    (FC outlet 2) to the joint primary/secondary sludge thickener GT.

    With BDF or LSODA, selected TF nutrient/storage absolute tolerances are
    at most 1e-6 mg/L and thickness tolerances at most 1e-9 m. Tighter caller
    values are preserved. Other unit states retain the caller's tolerances;
    rtol and max_step are not changed. Effective TF tolerances are available
    as system.tf_effective_atol after state assembly.
    """
    flowsheet = flowsheet or qs.Flowsheet(ID)
    qs.main_flowsheet.set_flowsheet(flowsheet)

    pc.create_masm2d_cmps()
    asm = _create_d1_tf_model()
    rww = default_rww()
    thermo_asm = qs.get_thermo()

    PC = su.PrimaryClarifier(
        'PC', ins=[rww, 'reject'],
        outs=('PE', 'PS'),
        isdynamic=True,
        sludge_flow_rate=0.074*MGD2cmd,
        solids_removal_efficiency=0.6
        )

    tf_kwargs = dict(tf_kwargs or ())
    n_zones = int(tf_kwargs.get('N_section', 6))
    tf_defaults = dict(
        N_section=n_zones,
        N_biofilm_layer=5,
        influent_fractions=[[1] + [0]*(n_zones-1)],
        recirculation_ratio=1.0,
        suspended_growth_model=asm,
        DO_ID='S_O2',
        max_liquid_film_thickness=5e-5,
        max_biofilm_thickness=1e-3,
        filter_bed_depth=2.0,
        filter_bed_surface=9464.0,
        media_specific_surface=100.0,
        initial_biofilm_thickness=4e-4,
        dynamic_biofilm_thickness=True,
        detachment_model='capacity',
        detachment_exponent=1.,
        detachment_rate=0.065,
        )
    if tf_kwargs:
        tf_defaults.update(tf_kwargs)

    TF = TricklingFilter(
        'TF', ins=PC-0, outs='TF_eff',
        **tf_defaults
        )

    FC = su.FlatBottomCircularClarifier(
        'FC', ins=TF-0, outs=['SE', '', 'WAS'],
        underflow=0., wastage=0.33*MGD2cmd,
        surface_area=1579.352, height=3.6576, N_layer=10, feed_layer=6,
        X_threshold=3000, v_max=410, v_max_practical=274,
        rh=3e-4, rp=5.2e-3, fns=0.001, 
        maximum_nonsettleable_solids=8.0,
        )

    GT = su.IdealClarifier(
        'GT', ins=[PC-1, FC-2], outs=['', 'thickened_sludge'],
        sludge_flow_rate=0.0475*MGD2cmd,     # aim for 5% TS
        solids_removal_efficiency=0.85
        )

    pc.create_adm1p_cmps()
    thermo_adm = qs.get_thermo()
    adm = pc.ADM1p(kLa=10.0)

    J1 = su.mASM2dtoADM1p('J1', upstream=GT-1, thermo=thermo_adm, isdynamic=True,
                          adm1_model=adm, asm2d_model=asm)
    AD = su.AnaerobicCSTR(
        'AD', ins=J1-0, outs=('biogas', 'digestate'),
        V_liq=1.2*MGD2cmd, V_gas=0.12*MGD2cmd,
        fixed_headspace_P=False, fraction_retain=0,
        T=T_ad, model=adm,
        pH_ctrl=7.0,
        )
    AD.algebraic_h2 = False
    J2 = su.ADM1ptomASM2d('J2', upstream=AD-1, thermo=thermo_asm, isdynamic=True,
                          adm1_model=adm, asm2d_model=asm)
    qs.set_thermo(thermo_asm)

    DW = su.IdealClarifier(
        'DW', J2-0, outs=('', 'cake'),
        sludge_flow_rate=0.0069*MGD2cmd,   # aim for 18% TS
        solids_removal_efficiency=0.9
        )
    MX = su.Mixer('MX', ins=[GT-0, DW-0])

    HD = su.HydraulicDelay('HD', ins=MX-0, outs=1-PC)

    if default_init_conds:
        TF.set_init_conc(concentrations=asinit)
        FC.set_init_solubles(**fcinit)
        FC.set_init_sludge_solids(**fcinit)
        FC.set_init_TSS(default_fctss_init)
        AD.set_init_conc(**adinit)

    sys = _D1System(
        ID,
        path=(PC, TF, FC, GT, J1, AD, J2, DW, MX, HD),
        recycle=(HD-0,)
        )

    sys.set_dynamic_tracker(TF, TF-0, FC-0, AD)

    return sys

#%%
@time_printer
def run(sys, t, t_step, method=None, **kwargs):
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}')
    print(f'Time span 0-{t}d \n')

    sys.simulate(
        # state_reset_hook='reset_cache',
        t_span=(0,t),
        # t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        print_t=True,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)


#%%
if __name__ == '__main__':
    sys = create_d1_system()
    dct = globals()
    dct.update(sys.flowsheet.to_dict())

    t = 300
    t_step = 1
    # method = 'RK45'
    # method = 'RK23'
    # method = 'DOP853'
    # method = 'Radau'
    method = 'BDF'
    # method = 'LSODA'

    run(sys, t, t_step, method=method)
    # biomass_IDs = ('X_H', 'X_PAO', 'X_AUT')
    # srt = get_SRT(sys, biomass_IDs,
    #               wastage=[WAS],
    #               active_unit_IDs=('TF',))
    # if srt: print(f'Estimated SRT assuming at steady state is {round(srt, 2)} days')

    # from exposan.werf import figures_path
    # sys.diagram(format='png', file=ospath.join(figures_path, f'{ID}'))
