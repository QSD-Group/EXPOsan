# -*- coding: utf-8 -*-
"""Run the nine TF validation cases, preserve the GPS-X workbook, and plot.

Example
-------
    from exposan.esap.analyses.tf_validation_sweep import run_validation_sweep
    run = run_validation_sweep(
        system_kwargs={'asm_kwargs': {'mu_H': 3.3, 'q_PHA': 6.0}},
        show_plot=True,
    )

Nested tf_kwargs/asm_kwargs/tf_init_kwargs overrides merge with the settings
from DEFAULT_SYSTEM_KWARGS below. Pass None to clear a nested override dict.
Concentrations use mg/L, section-total thicknesses use m, and electricity
uses kWh/day. mASM2d's S_NO3 is lumped NOx-N; it fills out_Nitrate N while
out_Nitrite N is zero, so the existing heatmap compares total NOx correctly.
"""

from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime
from hashlib import sha256
from inspect import signature
from io import BytesIO
import json
from math import ceil
from pathlib import Path
import posixpath
from xml.etree import ElementTree as ET
from zipfile import ZIP_DEFLATED, ZipFile
from zoneinfo import ZoneInfo

import matplotlib.pyplot as plt
import numpy as np
import qsdsan as qs

from exposan.esap.systems.TricklingFilter_validation import (
    create_trickling_filter_validation_system,
)
from exposan.esap.provenance import snapshot_sources
from .tf_validation_heatmap import (
    LAYERS, METRICS, QSD_SAN_ROW, SOURCE_PATH, STRENGTHS,
    plot_validation_heatmap,
)

__all__ = ('run_validation_sweep', 'ValidationRun')

# %% Calibrated library defaults (the launcher can select a different profile)
DEFAULT_SYSTEM_KWARGS = {
    'gpsx_kinetics': False,
    'gpsx_decay': False,
    'gpsx_aob_proxy': False,
    'tf_kwargs': {
        'N_section': 6,
        'filter_bed_surface': 9464,
        'recirculation_ratio': 1,
        'dynamic_biofilm_thickness': True,
        'initial_biofilm_thickness': 4e-4,
    },
    'asm_kwargs': {'mu_H': 3.3, 'q_PHA': 6.0, 'q_PP': 2.0, 'b_PP': 0.05},
}
DEFAULT_SIMULATE_KWARGS = {
    'method': 'BDF', 't_span': (0, 300), 'rtol': 1e-3, 'atol': 1e-2,
    'stop_at_steady_state': False,
}

_NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
_REL = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
_PACKAGE_REL = 'http://schemas.openxmlformats.org/package/2006/relationships'
_CONTENT = 'http://schemas.openxmlformats.org/package/2006/content-types'
def _TAG(name):
    return f'{{{_NS}}}{name}'


@dataclass
class ValidationRun:
    """Saved outputs, numeric case results, and the existing heatmap figure."""

    workbook_path: Path
    heatmap_path: Path
    results: dict
    variation: object
    figure: object
    source_path: Path | None = None


# %% Configuration and workbook layout
def _json_value(value):
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f'Cannot record a reproducible setting of type {type(value).__name__}')


def _prepare_kwargs(overrides):
    kwargs = deepcopy(DEFAULT_SYSTEM_KWARGS)
    for name, value in (overrides or {}).items():
        if name in ('tf_kwargs', 'asm_kwargs', 'tf_init_kwargs') and value is not None:
            kwargs[name] = {**(kwargs.get(name) or {}), **deepcopy(value)}
        else:
            kwargs[name] = deepcopy(value)
    if 'wastewater_strength' in kwargs or 'N_biofilm_layer' in (kwargs.get('tf_kwargs') or {}):
        raise ValueError('The sweep controls wastewater_strength and N_biofilm_layer.')
    if kwargs.get('flowsheet') is not None:
        raise ValueError('The sweep creates a fresh flowsheet for each case.')
    if (kwargs.get('tf_kwargs') or {}).get('N_section', 6) != 6:
        raise ValueError('The six-section workbook requires N_section=6.')
    signature(create_trickling_filter_validation_system).bind_partial(**kwargs)
    json.dumps(_json_value(kwargs), allow_nan=False)
    return kwargs


def _worksheet_parts(parts):
    workbook = ET.fromstring(parts['xl/workbook.xml'])
    relations = ET.fromstring(parts['xl/_rels/workbook.xml.rels'])
    targets = {r.get('Id'): r.get('Target') for r in relations}
    sheets = {}
    for sheet in workbook.find(_TAG('sheets')):
        target = targets[sheet.get(f'{{{_REL}}}id')]
        sheets[sheet.get('name')] = posixpath.normpath(
            target.lstrip('/') if target.startswith('/') else f'xl/{target}')
    return workbook, relations, sheets


def _cell_value(cell, strings):
    kind = cell.get('t')
    if kind == 'inlineStr':
        return ''.join(t.text or '' for t in cell.iter(_TAG('t')))
    value = cell.find(_TAG('v'))
    if value is None or value.text is None:
        return None
    if kind == 's':
        return strings[int(value.text)]
    if kind == 'b':
        return value.text == '1'
    if kind in ('str', 'e'):
        return value.text
    return float(value.text)


def _load_template(path):
    with ZipFile(path) as archive:
        parts = {name: archive.read(name) for name in archive.namelist()}
    _, _, sheet_parts = _worksheet_parts(parts)
    strings = []
    if 'xl/sharedStrings.xml' in parts:
        root = ET.fromstring(parts['xl/sharedStrings.xml'])
        strings = [''.join(t.text or '' for t in item.iter(_TAG('t'))) for item in root]
    required = {c for _, columns in METRICS for c in columns}
    required.update(f'{prefix}_{name}' for prefix in ('in', 'out')
                    for name in ('DO', 'cBOD', 'COD', 'Ammonia N',
                                 'Nitrate N', 'TN', 'PO4-P', 'TP'))
    layout = {}
    for strength in STRENGTHS:
        if strength not in sheet_parts:
            raise ValueError(f'Template is missing the {strength!r} sheet.')
        root = ET.fromstring(parts[sheet_parts[strength]])
        rows = list(root.find(_TAG('sheetData')))
        headers = {str(_cell_value(c, strings)).strip(): c.get('r').rstrip('0123456789')
                   for c in rows[0] if _cell_value(c, strings) is not None}
        missing = sorted(required.difference(headers))
        if missing:
            raise ValueError(f'{strength} sheet is missing columns: {missing}')
        cases, reference = {}, None
        for row in rows[1:]:
            cell = row.find(f"{_TAG('c')}[@r='A{row.get('r')}']")
            label = str(_cell_value(cell, strings)).strip() if cell is not None else ''
            if label == 'GPS-X':
                if reference is not None:
                    raise ValueError(f'Duplicate GPS-X row in {strength}.')
                reference = row
            match = QSD_SAN_ROW.fullmatch(label)
            if match:
                layers = int(match.group(1))
                if layers in cases:
                    raise ValueError(f'Duplicate {layers}-layer row in {strength}.')
                cases[layers] = row
        if reference is None or not set(LAYERS).issubset(cases):
            raise ValueError(f'{strength} requires GPS-X and all three QSDsan rows.')
        cells = {c.get('r'): c for c in reference}
        for _, columns in METRICS:
            for name in columns:
                value = _cell_value(cells[f'{headers[name]}{reference.get("r")}'], strings)
                if value is None or not np.isfinite(float(value)) or float(value) < 0:
                    raise ValueError(f'GPS-X {strength}.{name} has no finite, nonnegative value. '
                                     'Recalculate and save the template in Excel first.')
        layout[strength] = (root, headers, cases)
    return parts, layout


# %% Numeric results and reproducibility settings
def _case_results(system):
    tf = system.flowsheet.unit.TF
    values = {}
    for prefix, stream in (('in', tf.ins[0]), ('out', tf.outs[0])):
        for name, value in {
            'DO': stream.iconc['S_O2'], 'cBOD': stream.BOD, 'COD': stream.COD,
            'Ammonia N': stream.iconc['S_NH4'], 'Nitrate N': stream.iconc['S_NO3'],
            'TN': stream.TN, 'PO4-P': stream.iconc['S_PO4'], 'TP': stream.TP,
        }.items():
            values[f'{prefix}_{name}'] = float(value)
    values['out_Nitrite N'] = 0.
    for section, thickness in enumerate(tf.total_biofilm_thickness, 1):
        values[f'biofilm_{section}_thickness'] = float(thickness)
    values.update(
        influent_pumping=24.*system.get_influent_pumping_energy(),
        recirculation_pumping=24.*system.get_recirculation_pumping_energy(),
        rotating_arm=24.*system.get_rotating_arm_energy(),
    )
    invalid = [name for name, value in values.items() if not np.isfinite(value) or value < 0]
    if invalid:
        raise ValueError(f'Cannot export non-finite or negative final metrics: {invalid}')
    return values


def _setting_rows(kwargs, simulate_kwargs, generated, template, source_digest, results,
                  provenance=None):
    bound = signature(create_trickling_filter_validation_system).bind_partial(**kwargs)
    bound.apply_defaults()
    rows = [('Category', 'Argument / item', 'Value')]

    def add(category, mapping, prefix=''):
        for name, value in mapping.items():
            if name in ('wastewater_strength', 'N_biofilm_layer'):
                continue
            key = f'{prefix}.{name}' if prefix else name
            if isinstance(value, dict) and value:
                add(category, value, key)
            else:
                value = _json_value(value)
                if value is None or isinstance(value, (dict, list)):
                    value = json.dumps(value, allow_nan=False)
                rows.append((category, key, value))

    add('System creation', bound.arguments)
    add('Simulation', simulate_kwargs)
    add('Run metadata', {
        'Generated (America/Chicago)': generated.isoformat(),
        'Template': str(template), 'Template SHA256': source_digest,
        'QSDsan version': qs.__version__,
        'NOx convention': 'S_NO3 is total NOx-N; out_Nitrite N=0 (no separate NO2 state).',
        'BOD convention': 'in_cBOD/out_cBOD use WasteStream.BOD (carbonaceous BOD5).',
        'Energy convention': 'Power [kW] multiplied by 24 to obtain kWh/day.',
        'Thickness convention': 'Total thickness per section, not per computational layer.',
    })
    add('Source provenance', provenance or {})
    for (strength, layers), case in results.items():
        rows.append(('Completed cases', f'{strength} / {layers} layers: final time [d]',
                     case['final_time_d']))
    return rows


# %% Preserve the template's cells, formula caches, styles, and other ZIP parts
def _xml_bytes(root, source=None):
    """Retain prefixes named by Excel's markup-compatibility attributes."""
    root = deepcopy(root)
    namespaces = dict(item for _, item in ET.iterparse(
        BytesIO(source), events=('start-ns',))) if source else {'': _NS}
    for prefix, uri in namespaces.items():
        ET.register_namespace(prefix, uri)
    payload = ET.tostring(root, encoding='utf-8', xml_declaration=True)
    emitted = dict(item for _, item in ET.iterparse(BytesIO(payload), events=('start-ns',)))
    # ElementTree otherwise discards namespaces used only in mc:Ignorable/QNames.
    for prefix, uri in namespaces.items():
        if emitted.get(prefix) != uri:
            root.set(f'xmlns:{prefix}' if prefix else 'xmlns', uri)
    return ET.tostring(root, encoding='utf-8', xml_declaration=True)


def _set_cell(cell, value):
    for child in list(cell):
        cell.remove(child)
    if isinstance(value, bool):
        cell.set('t', 'b')
        ET.SubElement(cell, _TAG('v')).text = str(int(value))
    elif isinstance(value, (float, int)):
        cell.set('t', 'n')
        ET.SubElement(cell, _TAG('v')).text = str(value)
    else:
        cell.set('t', 'inlineStr')
        ET.SubElement(ET.SubElement(cell, _TAG('is')), _TAG('t')).text = str(value)


def _settings_xml(rows, body_style, header_style):
    root = ET.Element(_TAG('worksheet'))
    ET.SubElement(root, _TAG('dimension'), ref=f'A1:C{len(rows)}')
    views = ET.SubElement(root, _TAG('sheetViews'))
    view = ET.SubElement(views, _TAG('sheetView'), workbookViewId='0', showGridLines='0')
    ET.SubElement(view, _TAG('pane'), ySplit='1', topLeftCell='A2',
                  activePane='bottomLeft', state='frozen')
    ET.SubElement(root, _TAG('sheetFormatPr'), defaultRowHeight='18')
    columns = ET.SubElement(root, _TAG('cols'))
    for index, width in enumerate((24, 52, 90), 1):
        ET.SubElement(columns, _TAG('col'), min=str(index), max=str(index),
                      width=str(width), customWidth='1')
    data = ET.SubElement(root, _TAG('sheetData'))
    for index, values in enumerate(rows, 1):
        lines = max(ceil(len(str(value))/width) for value, width in zip(values, (24, 52, 90)))
        row = ET.SubElement(data, _TAG('row'), r=str(index),
                            ht=str(18*lines+4), customHeight='1')
        for letter, value in zip('ABC', values):
            _set_cell(ET.SubElement(row, _TAG('c'), r=f'{letter}{index}',
                                   s=str(header_style if index == 1 else body_style)), value)
    ET.SubElement(root, _TAG('autoFilter'), ref=f'A1:C{len(rows)}')
    return _xml_bytes(root)


def _save_workbook(parts, layout, results, settings, path):
    workbook, relations, sheet_parts = _worksheet_parts(parts)
    output = dict(parts)
    styles = ET.fromstring(parts['xl/styles.xml'])
    formats = styles.find(_TAG('numFmts'))
    if formats is None:
        formats = ET.Element(_TAG('numFmts'), count='0')
        styles.insert(0, formats)
    xfs = styles.find(_TAG('cellXfs'))
    format_ids = {f.get('formatCode'): int(f.get('numFmtId')) for f in formats}
    next_format = max([163, *format_ids.values()])+1
    style_ids = {}

    def result_style(base, code, wrap=False, bold=False):
        nonlocal next_format
        key = (base, code, wrap, bold)
        if key not in style_ids:
            if code not in format_ids:
                format_ids[code] = next_format
                ET.SubElement(formats, _TAG('numFmt'), numFmtId=str(next_format), formatCode=code)
                next_format += 1
                formats.set('count', str(len(formats)))
            xf = deepcopy(xfs[base])
            xf.set('numFmtId', str(format_ids[code]))
            xf.set('applyNumberFormat', '1')
            if wrap:
                alignment = xf.find(_TAG('alignment'))
                if alignment is None:
                    alignment = ET.SubElement(xf, _TAG('alignment'))
                alignment.set('wrapText', '1')
                alignment.set('vertical', 'center')
                xf.set('applyAlignment', '1')
            if bold:
                fonts = styles.find(_TAG('fonts'))
                font = deepcopy(fonts[int(xf.get('fontId', '0'))])
                if font.find(_TAG('b')) is None:
                    ET.SubElement(font, _TAG('b'))
                xf.set('fontId', str(len(fonts)))
                fonts.append(font)
                fonts.set('count', str(len(fonts)))
                xf.set('applyFont', '1')
            style_ids[key] = len(xfs)
            xfs.append(xf)
            xfs.set('count', str(len(xfs)))
        return style_ids[key]

    for strength, (root, headers, cases) in layout.items():
        for layers in LAYERS:
            row = cases[layers]
            cells = {c.get('r'): c for c in row}
            for name, value in results[strength, layers]['values'].items():
                coordinate = f'{headers[name]}{row.get("r")}'
                cell = cells.get(coordinate)
                if cell is None:
                    cell = ET.SubElement(row, _TAG('c'), r=coordinate)
                _set_cell(cell, value)
                code = ('0.0000000' if name.startswith('biofilm_') else
                        '0.00' if name in ('influent_pumping', 'recirculation_pumping', 'rotating_arm')
                        else '0.0000')
                cell.set('s', str(result_style(int(cell.get('s', '0')), code)))
        part = sheet_parts[strength]
        output[part] = _xml_bytes(root, parts[part])

    if 'run_settings' in sheet_parts:
        settings_part = sheet_parts['run_settings']
    else:
        settings_part = 'xl/worksheets/run_settings.xml'
        if settings_part in parts:
            raise ValueError('Template already uses the run_settings.xml package path.')
        ids = {relation.get('Id') for relation in relations}
        index = 1
        while f'rId{index}' in ids:
            index += 1
        relation_id = f'rId{index}'
        ET.SubElement(relations, f'{{{_PACKAGE_REL}}}Relationship',
                      Id=relation_id, Type=f'{_REL}/worksheet', Target='worksheets/run_settings.xml')
        sheets = workbook.find(_TAG('sheets'))
        sheet_id = 1 + max(int(sheet.get('sheetId')) for sheet in sheets)
        ET.SubElement(sheets, _TAG('sheet'), name='run_settings', sheetId=str(sheet_id),
                      attrib={f'{{{_REL}}}id': relation_id})
        content = ET.fromstring(parts['[Content_Types].xml'])
        ET.SubElement(content, f'{{{_CONTENT}}}Override', PartName=f'/{settings_part}',
                      ContentType='application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml')
        output['xl/workbook.xml'] = _xml_bytes(workbook, parts['xl/workbook.xml'])
        output['xl/_rels/workbook.xml.rels'] = _xml_bytes(relations, parts['xl/_rels/workbook.xml.rels'])
        output['[Content_Types].xml'] = _xml_bytes(content, parts['[Content_Types].xml'])
    output[settings_part] = _settings_xml(
        settings, result_style(0, 'General', wrap=True),
        result_style(0, 'General', wrap=True, bold=True))
    output['xl/styles.xml'] = _xml_bytes(styles, parts['xl/styles.xml'])
    with ZipFile(path, 'x', compression=ZIP_DEFLATED) as archive:
        for name, payload in output.items():
            archive.writestr(name, payload)


# %% Nine simulations followed by workbook-driven heatmap generation
def run_validation_sweep(system_kwargs=None, simulate_kwargs=None,
                         template_path=SOURCE_PATH, output_dir=None, show_plot=True):
    """Simulate low/mid/high wastewater at 3/4/5 layers and save one report.

    Parameters
    ----------
    system_kwargs : dict, optional
        Overrides for create_trickling_filter_validation_system. Nested
        parameter dictionaries merge with DEFAULT_SYSTEM_KWARGS. The sweep
        controls wastewater_strength, N_biofilm_layer, and fresh flowsheets;
        N_section must remain 6 to match the workbook. Use serializable values.
    simulate_kwargs : dict, optional
        Overrides for System.simulate; defaults are BDF, 0-300 d, rtol=1e-3,
        atol=1e-2, stop_at_steady_state=False. Each run starts from fresh initial
        conditions. Explicit early steady-state termination is permitted and
        its actual endpoint is recorded, never labelled as a 300-day result.
    template_path : path-like, optional
        Workbook containing the immutable GPS-X reference rows and case layout.
    output_dir : path-like, optional
        Destination for the timestamped .xlsx, .tif, and source .zip; defaults to the
        template folder. Existing reports and the original template are retained.
    show_plot : bool, optional
        Display the heatmap after saving it. Set False for unattended runs.

    Returns
    -------
    ValidationRun
        Saved paths (including source_path), case metrics, variation, and figure.
        Source/input capture occurs before simulation. A failed/incomplete
        simulation may leave that archive, but does not write a workbook.
    """
    kwargs = _prepare_kwargs(system_kwargs)
    simulation = {**deepcopy(DEFAULT_SIMULATE_KWARGS), **deepcopy(simulate_kwargs or {})}
    json.dumps(_json_value(simulation), allow_nan=False)
    template = Path(template_path).expanduser().resolve()
    parts, layout = _load_template(template)
    output_dir = Path(output_dir).expanduser().resolve() if output_dir else template.parent
    generated = datetime.now(ZoneInfo('America/Chicago'))
    stem = f'TF_specific_GPS-X_validation_{generated:%Y%m%d_%H%M%S_%f}'
    output_dir.mkdir(parents=True, exist_ok=True)
    source_digest = sha256(template.read_bytes()).hexdigest()
    source_path = snapshot_sources(
        output_dir / f'{stem}_sources.zip',
        settings={'system_kwargs': _json_value(kwargs),
                  'simulate_kwargs': _json_value(simulation),
                  'strengths': STRENGTHS, 'layers': LAYERS},
        extra_files=(template,),
    )
    provenance = {'Archive': str(source_path),
                  'Archive SHA256': sha256(source_path.read_bytes()).hexdigest()}
    results = {}
    for strength in STRENGTHS:
        for layers in LAYERS:
            print(f'Running {strength} wastewater / {layers} layers...', flush=True)
            case_kwargs = deepcopy(kwargs)
            case_kwargs['tf_kwargs'] = {**(case_kwargs.get('tf_kwargs') or {}),
                                        'N_biofilm_layer': layers}
            with qs.Flowsheet(f'TF_validation_{strength}_{layers}') as flowsheet:
                case_kwargs.update(flowsheet=flowsheet, wastewater_strength=strength)
                system = create_trickling_filter_validation_system(**case_kwargs)
                system.simulate(**deepcopy(simulation))
                sol = system.scope.sol
                if not sol.success or not len(sol.t) or not np.isfinite(sol.y).all():
                    raise RuntimeError(f'{strength}/{layers} simulation failed: {sol.message}')
                endpoint = float(sol.t[-1])
                if not simulation['stop_at_steady_state'] and not np.isclose(
                        endpoint, simulation['t_span'][1], rtol=0., atol=1e-8):
                    raise RuntimeError(f'{strength}/{layers} stopped at {endpoint:g} d, '
                                       f'before {simulation["t_span"][1]:g} d.')
                results[strength, layers] = {'values': _case_results(system),
                                            'final_time_d': endpoint}
            values = results[strength, layers]['values']
            print(f'  t={endpoint:g} d; COD={values["out_COD"]:.4f}; '
                  f'NH4-N={values["out_Ammonia N"]:.4f}; '
                  f'PO4-P={values["out_PO4-P"]:.4f} mg/L', flush=True)

    workbook_path = output_dir / f'{stem}.xlsx'
    heatmap_path = output_dir / f'{stem}_heatmap.tif'
    settings = _setting_rows(kwargs, simulation, generated, template,
                             source_digest, results, provenance)
    _save_workbook(parts, layout, results, settings, workbook_path)
    variation, figure, _, _ = plot_validation_heatmap(workbook_path, save_as=heatmap_path)
    print(f'Workbook: {workbook_path}\nHeatmap:  {heatmap_path}\nSources:  {source_path}', flush=True)
    if show_plot:
        plt.show()
    return ValidationRun(workbook_path, heatmap_path, results, variation, figure,
                         source_path=source_path)


# %% Command-line entry point with the calibrated library profile
if __name__ == '__main__':
    run_validation_sweep()
