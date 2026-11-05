# -*- coding: utf-8 -*-
"""Run the WERF baseline-performance workflow for original and new configurations.

Edit BASELINE_CONFIGURATIONS to include future registered systems. This script
uses the same TSS targets, underflow adjustment and reporting as
baseline_performance.py, with the original configurations plus D1 and D1E.

This module is under the University of Illinois/NCSA Open Source License.
See https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt.
"""

from pathlib import Path
import time as tm

import numpy as np
import pandas as pd
from biosteam.evaluation._utils import var_columns
from qsdsan import Model

from exposan.werf import (
    SYSTEM_CREATORS, add_OPEX_metrics, add_performance_metrics, create_system,
    results_path,
    )
from exposan.werf.utils import (
    cache_state, plantwide_N_mass_flows, plantwide_P_mass_flows,
    )

__all__ = ('run_baseline_performance_extension',)


# %% Editable configuration and simulation settings
BASELINE_CONFIGURATIONS = (
    'B1', 'B1E', 'B2', 'B3',
    'C1', 'C1E', 'C2', 'C3',
    'D1', 'D1E',
    'E2', 'E2P',
    'F1', 'F1E',
    'G1', 'G1E', 'G2', 'G3',
    'H1', 'H1E',
    'I1', 'I1E', 'I2', 'I3',
    'N1', 'N1E', 'N2',
    )
SIMULATE_KWARGS = dict(t_span=(0, 300), method='BDF')
THICKENED_SLUDGE_TSS = 5e4  # mg/L; same 5% TSS proxy as the original analysis.
TSS_RTOL = 0.01
MAX_TSS_ADJUSTMENTS = 30


# %% Simulation and baseline evaluation
def _simulate(system, options):
    """Retry once after a cache reset, preserving both errors if that fails."""
    try:
        system.simulate(**options)
    except Exception as initial_error:
        try:
            system.simulate(**{**options, 'state_reset_hook': 'reset_cache'})
        except Exception as retry_error:
            raise RuntimeError(
                f'{system.ID}: initial simulation failed '
                f'({type(initial_error).__name__}: {initial_error}); '
                f'cache-reset retry failed '
                f'({type(retry_error).__name__}: {retry_error}).'
                ) from retry_error
    solution = system.scope.sol
    if not solution.success or solution.status != 0:
        raise RuntimeError(f'{system.ID}: {solution.message}')
    if not np.isfinite(solution.y).all():
        raise RuntimeError(f'{system.ID}: simulation returned nonfinite states.')


def _simulate_to_tss(system, options, thickened, thickener, cake_tss,
                     max_tss_adjustments=MAX_TSS_ADJUSTMENTS):
    """Simulate and adjust thickener/DW underflows to the WERF TSS targets.

    Shared by baseline generation and manual operating-condition comparisons,
    including configurations with no published baseline_underflows entry.
    Return the number of underflow updates. Reject failed/nonfinite runs and
    invalid solids concentrations; stop after the configured update limit.
    """
    if not isinstance(max_tss_adjustments, int) or max_tss_adjustments < 0:
        raise ValueError('max_tss_adjustments must be a nonnegative integer.')
    streams, units = system.flowsheet.stream, system.flowsheet.unit
    for adjustment in range(max_tss_adjustments + 1):
        _simulate(system, options)
        ratios = np.array((thickened.get_TSS()/THICKENED_SLUDGE_TSS,
                           streams.cake.get_TSS()/cake_tss))
        if not np.isfinite(ratios).all() or (ratios <= 0.).any():
            raise RuntimeError('Sludge TSS must be finite and positive.')
        if (np.abs(ratios - 1.) <= TSS_RTOL).all():
            return adjustment
        if adjustment == max_tss_adjustments:
            raise RuntimeError(
                f'TSS targets not reached after {adjustment} adjustments; '
                f'target ratios: {ratios}.')
        print(f'Adjusting TS% ... {ratios[0]:.3f}  {ratios[1]:.3f}')
        thickener.sludge_flow_rate *= ratios[0]
        units.DW.sludge_flow_rate *= ratios[1]


def run_baseline_performance_extension(
        configurations=None, *, simulate_kwargs=None, output_dir=None,
        max_tss_adjustments=MAX_TSS_ADJUSTMENTS):
    """Simulate, adjust sludge underflows and report WERF baselines.

    Parameters
    ----------
    configurations : iterable of str, optional
        Registered configuration IDs. Defaults to BASELINE_CONFIGURATIONS,
        containing the original baseline list plus D1 and D1E.
        Future additions need only be registered and added to that tuple.
    simulate_kwargs : dict, optional
        Overrides for SIMULATE_KWARGS, which defaults to 300 days with BDF.
        D1 and D1E automatically retain their TF-specific tolerance guards.
    output_dir : path-like, optional
        Defaults to the WERF results directory. The three workbooks have an
        extension suffix to avoid overwriting the original baseline workbooks.
        States use steady_states/baseline_unopt/<ID>.npy beneath this directory.
    max_tss_adjustments : int, optional
        Maximum underflow updates after the initial simulation; defaults to 30.
        Unconverged cases are reported as failures rather than cached.

    Returns
    -------
    pandas.DataFrame
        Successful configuration rows with the original metric-column format.
        attrs['failures'] contains any omitted configuration/error pairs.

    Notes
    -----
    Uses factory underflows, not another configuration's calibrated values.
    As in the original script, targets thickened sludge at 5e4 mg TSS/L and
    cake at 18e4 (primary sludge plus AD), 17e4 (AED/no primary sludge), or
    20e4 mg TSS/L (primary sludge without digestion), all within 1%.
    These are empirical TSS proxies, not a separate total-solids model.
    Each AD case is simulated once per adjustment and its feed is reported
    in the same pass, rather than repeated in a separate AD-only loop.
    Importing this module does not run simulations or write outputs.
    """
    selected = BASELINE_CONFIGURATIONS if configurations is None else configurations
    if isinstance(selected, str):
        selected = (selected,)
    selected = tuple(ID.upper() for ID in selected)
    if not selected or len(set(selected)) != len(selected):
        raise ValueError('Provide at least one unique configuration ID.')
    unknown = set(selected).difference(SYSTEM_CREATORS)
    if unknown:
        raise ValueError(f'Unregistered configurations: {sorted(unknown)}')
    if not isinstance(max_tss_adjustments, int) or max_tss_adjustments < 0:
        raise ValueError('max_tss_adjustments must be a nonnegative integer.')

    options = {**SIMULATE_KWARGS, **(simulate_kwargs or {})}
    output_dir = Path(results_path if output_dir is None else output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics, nitrogen, phosphorus, failures = {}, {}, {}, {}

    for ID in selected:
        print(f'System {ID}\n{"="*30}')
        system = None
        try:
            system = create_system(ID)
            streams, units = system.flowsheet.stream, system.flowsheet.unit
            model = Model(system)
            add_performance_metrics(model)
            add_OPEX_metrics(model)

            if 'PS' not in streams:
                cake_tss = 17e4
            elif 'AD' in units:
                cake_tss = 18e4
            elif 'AED' in units:
                cake_tss = 17e4
            else:
                cake_tss = 20e4
            if 'thickened_WAS' in streams:
                thickened, thickener = streams.thickened_WAS, units.MT
            else:
                thickened, thickener = streams.thickened_sludge, units.GT

            start = tm.perf_counter()
            print('Start time:', tm.strftime('%H:%M:%S', tm.localtime()))
            _simulate_to_tss(system, options, thickened, thickener, cake_tss,
                             max_tss_adjustments)

            ndf = plantwide_N_mass_flows(system)
            pdf = plantwide_P_mass_flows(system)
            values = pd.Series([metric() for metric in model.metrics],
                               index=var_columns(model.metrics))
            cache_state(system, folder=str(output_dir/'steady_states'/'baseline_unopt'))
            nitrogen[ID], phosphorus[ID], metrics[ID] = ndf, pdf, values
            print('Final underflows:',
                  f'{thickener.sludge_flow_rate:.2f}, {units.DW.sludge_flow_rate:.2f}')
            print(f'Duration: {tm.perf_counter()-start:.1f} s\n')
            if 'AD' in units:
                feed = units.AD.ins[0]
                print(f'System {ID}-baseline AD influent flow ~ {feed.F_mass} '
                      f'kg/hour, AD influent temperature ~ {feed.T} K')
        except Exception as error:
            failures[ID] = str(error)
            print(f'{ID} failed: {error}\n')
        finally:
            if system is not None:
                system.flowsheet.clear()

    if not metrics:
        raise RuntimeError(f'No extension configurations completed: {failures}')

    #  Export the same table structures under extension-specific file names.
    for filename, tables in (
            ('N_mass_unopt_extension.xlsx', nitrogen),
            ('P_mass_unopt_extension.xlsx', phosphorus)):
        with pd.ExcelWriter(output_dir/filename) as writer:
            for ID, table in tables.items():
                table.to_excel(writer, sheet_name=ID)
    performance = pd.DataFrame(metrics).T
    performance.attrs['failures'] = failures
    performance.to_excel(output_dir/'_baseline_unopt_performance_extension.xlsx')
    return performance


# %% Run the selected baseline configurations when launched as a script.
if __name__ == '__main__':
    metrics = run_baseline_performance_extension()
