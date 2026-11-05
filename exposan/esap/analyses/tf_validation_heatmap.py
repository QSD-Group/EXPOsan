# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed and tested by:

    Zixuan Wang <wyatt4428@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from pathlib import Path
import re

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.mathtext import _mathtext as mathtext
import numpy as np
import pandas as pd


# %% Shared figure style
plt.rcParams.update({
    'font.sans-serif': 'Arial',
    'mathtext.fontset': 'custom',
    'mathtext.default': 'regular',
    'mathtext.rm': 'Arial',
    'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
    'figure.max_open_warning': 100,
    })
mathtext.FontConstantsBase.sup1 = 0.35
mathtext.SHRINK_FACTOR = 0.60

# %% Workbook inputs and metric definitions
SOURCE_PATH = Path(
    '/Users/zixuanwang/Library/CloudStorage/GoogleDrive-'
    'wyatt4428@gmail.com/My Drive/4- Research/430- UIUC/434. Projects/'
    '43412. ARPA_E_N_Mg_recovery/process_models/TF_validation/'
    'TF_specific_GPS-X_validation.xlsx'
    )
FIGURES_PATH = Path(__file__).resolve().parents[1] / 'figures'

STRENGTHS = ('low', 'mid', 'high')
LAYERS = (3, 4, 5)
ROW_PREFIX = {'low': 'L', 'mid': 'M', 'high': 'H'}
QSD_SAN_ROW = re.compile(r'QSDsan\s*-\s*(\d+)\s+biofilm\s+layers', re.I)

# Each tuple contains the heatmap abbreviation and source columns to sum before
# comparison. Single-column tuples are retained for one consistent calculation.
METRICS = (
    ('DO', ('out_DO',)),
    ('BOD', ('out_cBOD',)),
    ('COD', ('out_COD',)),
    ('NH4+-N', ('out_Ammonia N',)),
    ('NOx--N', ('out_Nitrate N', 'out_Nitrite N')),
    ('TN', ('out_TN',)),
    ('PO43--P', ('out_PO4-P',)),
    ('TP', ('out_TP',)),
    ('BT1', ('biofilm_1_thickness',)),
    ('BT2', ('biofilm_2_thickness',)),
    ('BT3', ('biofilm_3_thickness',)),
    ('BT4', ('biofilm_4_thickness',)),
    ('BT5', ('biofilm_5_thickness',)),
    ('BT6', ('biofilm_6_thickness',)),
    ('E', ('influent_pumping', 'recirculation_pumping', 'rotating_arm')),
    )

COLUMN_LABELS = (
    'DO', 'BOD', 'COD', 'NH$_4^+$-N', 'NO$_X^-$-N', 'TN', 'H$_x$PO$_4^{3-x}$-P','TP',
    'BT1', 'BT2', 'BT3', 'BT4', 'BT5', 'BT6', 'E',
    )


# %% Workbook parsing and max-normalized comparison
def _metric_values(row):
    """Return the selected and pre-summed validation metrics for one result."""
    values = {}
    for abbreviation, columns in METRICS:
        values[abbreviation] = row.loc[list(columns)].sum(
            min_count=len(columns))
    values = pd.Series(values, dtype=float)
    if not np.isfinite(values).all():
        invalid = values.index[~np.isfinite(values)].tolist()
        raise ValueError(f'Non-finite validation metrics: {invalid}')
    return values


def calculate_percent_variation(reference, evaluated):
    """Return the signed max-normalized difference in the interval [-1, 1]
    for nonnegative validation metrics without mutating inputs."""
    reference = pd.Series(reference, dtype=float)
    evaluated = pd.Series(evaluated, dtype=float)
    if (reference < 0.).any() or (evaluated < 0.).any():
        raise ValueError(
            'Max-normalized validation metrics must be nonnegative')
    scale = pd.concat((reference.abs(), evaluated.abs()), axis=1).max(axis=1)
    variation = evaluated.sub(reference).div(scale)
    variation.loc[scale == 0.] = 0.
    if (variation.abs() > 1. + 1e-12).any():
        raise ValueError('Max-normalized variation is outside [-1, 1]')
    return variation


def compile_percent_variation(source_path=SOURCE_PATH):
    """Compile the 9-by-15 QSDsan variation matrix from the validation workbook.

    Rows combine wastewater strength and biofilm layers (L3-L5, M3-M5,
    H3-H5). Columns follow ``METRICS``. NOx and energy are summed before the
    percentage comparison is made.
    """
    source_path = Path(source_path).expanduser()
    if not source_path.is_file():
        raise FileNotFoundError(source_path)

    sheets = pd.read_excel(
        source_path, sheet_name=list(STRENGTHS), header=0,
        skiprows=[1], index_col=0,
        )
    required_columns = {
        column for _, columns in METRICS for column in columns
        }
    records = []
    labels = []

    for strength in STRENGTHS:
        frame = sheets[strength].copy()
        frame.index = frame.index.map(lambda value: str(value).strip())
        frame.columns = frame.columns.map(lambda value: str(value).strip())
        missing = sorted(required_columns.difference(frame.columns))
        if missing:
            raise KeyError(f'{strength!r} sheet is missing columns: {missing}')
        frame.loc[:, list(required_columns)] = frame.loc[
            :, list(required_columns)].apply(pd.to_numeric, errors='raise')
        if 'GPS-X' not in frame.index:
            raise KeyError(f'{strength!r} sheet is missing the GPS-X row')

        qsd_rows = {}
        for row_name in frame.index:
            match = QSD_SAN_ROW.fullmatch(row_name)
            if match:
                qsd_rows[int(match.group(1))] = row_name

        reference = _metric_values(frame.loc['GPS-X'])
        for layers in LAYERS:
            try:
                row_name = qsd_rows[layers]
            except KeyError:
                raise KeyError(
                    f'{strength!r} sheet is missing the QSDsan {layers}-layer row'
                    ) from None
            evaluated = _metric_values(frame.loc[row_name])
            records.append(calculate_percent_variation(reference, evaluated))
            labels.append(f'{ROW_PREFIX[strength]}{layers}')

    return pd.DataFrame(records, index=labels,
                        columns=[name for name, _ in METRICS])


# %% Figure rendering and export
def _format_percent(value):
    if not np.isfinite(value):
        return ''
    if abs(value) >= 0.095:
        return f'{value:.0%}'
    return f'{value:.1%}'


def heatmap(data, annotate=None, save_as=None):
    """Draw a WERF-style diverging percentage heatmap."""
    fig, ax = plt.subplots(figsize=(14, data.shape[0]*0.8))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'validation_difference', ('#60c1cf', '#ffffff', '#ED586F'), N=257)
    im = ax.imshow(data, cmap=cmap, vmin=-1., vmax=1.)
    cbar = ax.figure.colorbar(
        im, ax=ax, orientation='horizontal', pad=0.05, fraction=0.05,
        )
    cbar.ax.set_xticks(
        ticks=[-1, 0, 1], labels=['-100%', '0%', '+100%'], fontsize=12,
        )

    ax.set_xticks(range(data.shape[1]), labels=COLUMN_LABELS, fontsize=12)
    ax.set_yticks(range(data.shape[0]), labels=data.index, fontsize=12)
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    ax.spines[:].set_visible(False)
    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which='minor', color='w', linestyle='-', linewidth=3)
    ax.tick_params(which='minor', bottom=False, left=False)

    if annotate is not None:
        for i in range(annotate.shape[0]):
            for j in range(annotate.shape[1]):
                value = annotate.iat[i, j]
                fill = data.iat[i, j]
                color = 'white' if abs(fill) > 0.6 else 'black'
                ax.text(
                    j, i, _format_percent(value), color=color,
                    horizontalalignment='center', verticalalignment='center',
                    size=9,
                    )

    fig.subplots_adjust(bottom=0.12, top=0.98, left=0.06, right=0.99)
    if save_as:
        output = Path(save_as)
        if not output.is_absolute():
            output = FIGURES_PATH / output
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=300, transparent=True, bbox_inches='tight')
    return fig, ax, im


def plot_validation_heatmap(source_path=SOURCE_PATH,
                            save_as='validation_heatmap.tif'):
    """Compile percentage differences and draw the validation heatmap."""
    variation = compile_percent_variation(source_path)
    fill = variation.copy()
    fig, ax, im = heatmap(fill, annotate=variation, save_as=save_as)
    return variation, fig, ax, im


# %% Plot an existing workbook without simulation
if __name__ == '__main__':
    percentage_variation, *_ = plot_validation_heatmap()
    print(percentage_variation.to_string(
        float_format=lambda value: f'{value:.2%}'))
