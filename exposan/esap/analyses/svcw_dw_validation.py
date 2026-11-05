#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot the two DW liquid-outlet panels against the SVCW-AD grab sample.

Run this file to create/simulate SVCW and save a 300-d comparison figure, or
plot an existing simulation without rerunning it:

    from exposan.esap.analyses.svcw_dw_validation import plot_svcw_dw_validation
    fig, axes = plot_svcw_dw_validation(sys)

The style follows batch_time_resolved.py: Arial, Guest group colors, 3-pt
lines/spines, inward/outward primary ticks and inward mirrored ticks, bold
25-pt labels, and mg per liter written as mg.L^-1 with a centered dot.
"""

from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.mathtext import _mathtext as mathtext
from matplotlib.patches import Patch
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, MultipleLocator
import numpy as np
from qsdsan.utils import palettes

__all__ = ('plot_svcw_dw_validation', 'measured_concentration_ranges')

# %% Editable simulation and figure settings
SYSTEM_KWARGS = {}
SIMULATE_KWARGS = {
    'method': 'BDF',
    't_span': (0, 300),
    'rtol': 1e-3,
    'atol': 1e-2,
}
# SVCW integrates the full interval; it has no ESAP early-stop monitor.
PLOT_WINDOW = (0, 300)
CONCENTRATION_BASIS = 'element'  # 'element': NH4-N/PO4-P; 'ion': NH4+/PO4(3-).
OUTPUT_DIR = Path(
    '/Users/zixuanwang/Library/CloudStorage/GoogleDrive-'
    'wyatt4428@gmail.com/My Drive/4- Research/430- UIUC/434. Projects/'
    '43412. ARPA_E_N_Mg_recovery/process_models/SVCW_configuration/'
)
FIGURE_PATH = OUTPUT_DIR / 'SVCW_DW_liquid_outlet.jpg'

# %% Measured mean and supplied +/- range, in mg of the stated ion per liter
MEASURED_ION_MG_L = {
    'S_NH4': (1380.06, 137.9),
    'S_PO4': (131.1, 1.3),
    'S_Mg': (19.3, 1.6),
    'S_Ca': (64.9, 2.1),
}
# Match the reference's PO4 -> P conversion (31/95); use 14/18 for NH4 -> N.
ION_TO_ELEMENT = {'S_NH4': 14/18, 'S_PO4': 31/95, 'S_Mg': 1., 'S_Ca': 1.}

# %% Reference figure style and Guest group colors
GUEST = palettes['Guest']
COLORS = {
    'S_NH4': GUEST.blue.HEX,
    'S_Mg': GUEST.green.HEX,
    'S_Ca': GUEST.purple.HEX,
    'S_PO4': GUEST.red.HEX,
}
STYLE = {
    'axes.linewidth': 3,
    'hatch.linewidth': 3,
    'xtick.labelsize': 25,
    'ytick.labelsize': 21,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'mathtext.fontset': 'custom',
    'mathtext.default': 'regular',
    'mathtext.rm': 'Arial',
    'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
    'figure.max_open_warning': 100,
}


# %% Measurement basis conversion and recorded stream extraction
def measured_concentration_ranges(concentration_basis='element', measured=None):
    """Return measured mean/lower/upper bounds in the selected plotting basis.

    Input values always use ion mass, as in the provided SVCW-AD table. The
    supplied +/- quantities are preserved without assuming SD or confidence
    intervals. In element mode NH4 and PO4 are converted to N and P; Mg and
    Ca retain their elemental metal concentrations in either mode.
    """
    if concentration_basis not in ('element', 'ion'):
        raise ValueError("concentration_basis must be 'element' or 'ion'.")
    measurements = MEASURED_ION_MG_L if measured is None else measured
    ranges = {}
    for component in ION_TO_ELEMENT:
        mean, uncertainty = np.asarray(measurements[component], dtype=float)
        if not np.isfinite([mean, uncertainty]).all() or not 0 <= uncertainty <= mean:
            raise ValueError(f'Invalid measured concentration/range for {component}.')
        factor = ION_TO_ELEMENT[component] if concentration_basis == 'element' else 1.
        ranges[component] = {
            'mean': mean*factor,
            'lower': (mean-uncertainty)*factor,
            'upper': (mean+uncertainty)*factor,
        }
    return ranges


def _tracked_dw_series(system, start, stop, concentration_basis):
    """Read ADM1p DW records, also accepting older mASM2d DW histories."""
    frame = system.scope.export()
    time = frame[('-', 't [d]')].to_numpy(dtype=float)
    mask = np.isfinite(time) & (time >= start) & (time <= stop)
    if not mask.any():
        raise ValueError(f'No recorded simulation times within {start:g}-{stop:g} d.')
    time = time[mask]
    if np.any(np.diff(time) < 0):
        raise ValueError('The recorded simulation times are not ordered.')
    stream = system.flowsheet.unit.DW.outs[0]
    series = {}
    source_ids = {
        'S_NH4': ('S_NH4', 'S_IN'),
        'S_PO4': ('S_PO4', 'S_IP'),
    }
    for component, factor in ION_TO_ELEMENT.items():
        aliases = source_ids.get(component, (component,))
        matches = [column for column in frame.columns
                   if column[0] == stream.ID
                   and column[1].split(' [', 1)[0] in aliases]
        if len(matches) != 1:
            raise KeyError(f'Expected one tracked {stream.ID}.{component} column; '
                           f'accepted IDs: {aliases}; found {matches}. '
                           'Track DW.outs[0] before simulation.')
        values = frame[matches[0]].to_numpy(dtype=float)[mask]
        if not np.isfinite(values).all():
            raise ValueError(f'Non-finite DW outlet concentrations for {component}.')
        series[component] = values if concentration_basis == 'element' else values/factor
    return time, series


# %% Axis formatting shared by the two DW panels
def _style_axis(ax, start, stop, limits):
    ax.set_xlim(start, stop)
    locator = MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10])
    ticks = locator.tick_values(*limits)
    ax.set_ylim(ticks[0], ticks[-1])
    ax.yaxis.set_major_locator(MultipleLocator(ticks[1]-ticks[0]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator((stop-start)/6))
    ax.xaxis.set_minor_locator(MultipleLocator((stop-start)/30))
    for axis in ('x', 'y'):
        ax.tick_params(axis=axis, which='major', direction='inout',
                       length=16, width=3, color='black', pad=5)
        ax.tick_params(axis=axis, which='minor', direction='inout',
                       length=8, width=3, color='black')
    top = ax.secondary_xaxis('top')
    right = ax.secondary_yaxis('right')
    for mirror, axis in ((top, 'x'), (right, 'y')):
        mirror.tick_params(axis=axis, which='major', direction='in',
                           length=8, width=3, color='black',
                           labeltop=False, labelright=False)
        mirror.tick_params(axis=axis, which='minor', direction='in',
                           length=4, width=3, color='black')
    top.xaxis.set_major_locator(MultipleLocator((stop-start)/6))
    top.xaxis.set_minor_locator(MultipleLocator((stop-start)/30))
    right.yaxis.set_major_locator(MultipleLocator(ticks[1]-ticks[0]))
    right.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.set_xlabel(r'$\mathbf{Duration}$'+' [d]', fontname='Arial', fontsize=25)
    ax.grid(False)


# %% Two-panel figure using only the DW liquid outlet
def plot_svcw_dw_validation(system, start=0, stop=300,
                            concentration_basis='element', measured=None,
                            save_as=None, show=True):
    """Plot a recorded SVCW run and the measured SVCW-AD concentration ranges.

    Parameters
    ----------
    system : qsdsan.System
        Already-simulated SVCW system with DW.outs[0] tracked. This function
        does not run or change the system. DW now uses ADM1p S_IN/S_IP;
        older mASM2d S_NH4/S_PO4 histories are also accepted on the same
        elemental N/P concentration basis.
    start, stop : float, optional
        Requested x-axis limits in days; only recorded data are plotted.
    concentration_basis : {'element', 'ion'}, optional
        Default 'element' uses model mg N/L and mg P/L, converting measured
        NH4+/PO4(3-) masses. 'ion' instead converts the model to mg NH4+/L
        and mg PO4/L and uses the original measured 1242.16-1517.96 mg NH4/L
        range. Metal concentrations remain mg Mg/L and mg Ca/L.
    measured : dict, optional
        Component ID -> (mean, supplied +/- value), always in ion mg/L.
        Defaults to the four values in the user's SVCW-AD table.
    save_as : path-like, optional
        JPG output path; saved at 300 dpi like batch_time_resolved.py.
        Other filename extensions supported by Matplotlib are also accepted.
    show : bool, optional
        Display the figure. Set False for an unattended run.

    Returns
    -------
    fig, axes
        Matplotlib figure and its two primary axes. Measurement bands span
        the requested x range at alpha=0.5, beneath the simulation curves.
    """
    start, stop = float(start), float(stop)
    if not np.isfinite([start, stop]).all() or stop <= start:
        raise ValueError('Plot limits must be finite with stop > start.')
    ranges = measured_concentration_ranges(concentration_basis, measured)
    time, series = _tracked_dw_series(system, start, stop, concentration_basis)
    labels = {
        'S_NH4': r'NH$_4^+$-N' if concentration_basis == 'element' else r'NH$_4^+$',
        'S_Mg': r'Mg$^{2+}$',
        'S_Ca': r'Ca$^{2+}$',
        'S_PO4': r'PO$_4^{3-}$-P' if concentration_basis == 'element' else r'PO$_4^{3-}$',
    }
    original_sup = mathtext.FontConstantsBase.sup1
    try:
        mathtext.FontConstantsBase.sup1 = 0.35
        with plt.rc_context(STYLE):
            fig, axes = plt.subplots(1, 2, figsize=(16, 6.4))
            for ax, components in zip(axes, (('S_NH4',), ('S_Mg', 'S_Ca', 'S_PO4'))):
                minimum, maximum = 0., 0.
                for component in components:
                    sample = ranges[component]
                    ax.axhspan(sample['lower'], sample['upper'], color=COLORS[component],
                               alpha=0.5, linewidth=0, zorder=1)
                    ax.plot(time, series[component], color=COLORS[component],
                            linewidth=3, label=labels[component], zorder=3)
                    minimum = min(minimum, float(np.min(series[component])), sample['lower'])
                    maximum = max(maximum, float(np.max(series[component])), sample['upper'])
                padding = max((maximum-minimum)*0.08, 1.)
                _style_axis(ax, start, stop, (minimum, maximum+padding))
                ax.legend(frameon=False, fontsize=17, loc='upper right')

            nutrient = r'$\mathbf{NH_4^+-N}$' if concentration_basis == 'element' else r'$\mathbf{NH_4^+}$'
            units = '\n[mg\u00b7L$^{-1}$]'
            axes[0].set_ylabel(nutrient+units, fontname='Arial', fontsize=25)
            axes[1].set_ylabel(r'$\mathbf{Concentration}$'+units, fontname='Arial', fontsize=25)
            fig.legend(handles=[
                Line2D([], [], color='black', linewidth=3, label='Simulation'),
                Patch(facecolor=GUEST.gray.HEX, alpha=0.5, label='SVCW-AD measured range'),
            ], loc='lower center', ncol=2, frameon=False, fontsize=16,
                bbox_to_anchor=(0.5, 0.005))
            fig.tight_layout(rect=(0, 0.09, 1, 1), w_pad=3)
            if save_as is not None:
                output = Path(save_as).expanduser()
                output.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(output, dpi=300, bbox_inches='tight', facecolor='white')
                print(f'Saved: {output.resolve()}')
            if show:
                plt.show()
    finally:
        mathtext.FontConstantsBase.sup1 = original_sup
    return fig, axes


# %% Standalone simulation and export (also runnable with %run in a notebook)
if __name__ == '__main__':
    from exposan import esap
    from exposan.esap.provenance import snapshot_sources

    stamp = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%fZ')
    figure_path = FIGURE_PATH.with_stem(f'{FIGURE_PATH.stem}_{stamp}')
    source_path = snapshot_sources(
        figure_path.with_suffix('.sources.zip'),
        settings={'system_kwargs': SYSTEM_KWARGS, 'simulate_kwargs': SIMULATE_KWARGS,
                  'concentration_basis': CONCENTRATION_BASIS,
                  'measured_ion_mg_L': MEASURED_ION_MG_L,
                  'plot_window': PLOT_WINDOW, 'figure': str(figure_path)},
    )
    print(f'Source archive: {source_path}')
    sys = esap.create_svcw_system(**deepcopy(SYSTEM_KWARGS))
    sys.simulate(**deepcopy(SIMULATE_KWARGS))
    solution = sys.scope.sol
    if not solution.success:
        raise RuntimeError(f'SVCW simulation failed: {solution.message}')
    if not np.isclose(solution.t[-1], SIMULATE_KWARGS['t_span'][1]):
        raise RuntimeError(f'SVCW simulation ended early at {solution.t[-1]:g} d.')
    print(f'Final recorded time: {solution.t[-1]:g} d')
    print(f'Measured bands ({CONCENTRATION_BASIS} basis), mg/L:')
    for component, values in measured_concentration_ranges(CONCENTRATION_BASIS).items():
        print(f'  {component}: {values["mean"]:.4f} '
              f'[{values["lower"]:.4f}, {values["upper"]:.4f}]')
    fig, axes = plot_svcw_dw_validation(
        sys, start=PLOT_WINDOW[0], stop=PLOT_WINDOW[1],
        concentration_basis=CONCENTRATION_BASIS,
        save_as=figure_path,
    )
