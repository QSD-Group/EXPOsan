# -*- coding: utf-8 -*-
"""Trickling-filter calibration, startup, solver, reporting and unit-physics tests.

This module is under the University of Illinois/NCSA Open Source License.
See https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt.
"""

from copy import copy
from hashlib import sha256
from inspect import signature
from io import BytesIO
import json
from types import SimpleNamespace
from xml.etree import ElementTree as ET
from zipfile import ZIP_DEFLATED, ZipFile

import matplotlib.pyplot as plt
import numpy as np
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font
import pytest
from qsdsan import process_models as pc

import exposan.esap.models as esap_models
from exposan.esap.analyses import tf_validation_sweep as sweep
from exposan.esap.models import get_trickling_filter_rotating_arm_power
from exposan.esap.systems.TricklingFilter_validation import (
    _create_validation_model,
    create_trickling_filter_validation_model,
    create_trickling_filter_validation_system,
)
from exposan.werf import default_as_init, default_rww


# %% Sweep fixtures, workbook/heatmap export and numerical configuration
@pytest.fixture
def template(tmp_path):
    headers = [None, 'in_DO', 'out_DO', 'in_cBOD', 'out_cBOD',
               'in_COD', 'out_COD', 'in_Ammonia N', 'out_Ammonia N',
               'in_Nitrate N', 'out_Nitrate N', 'out_Nitrite N',
               'in_TN', 'out_TN', 'in_PO4-P', 'out_PO4-P', 'in_TP', 'out_TP']
    headers += [f'biofilm_{i}_thickness' for i in range(1, 7)]
    headers += ['influent_pumping', 'recirculation_pumping', 'rotating_arm']
    workbook = Workbook()
    workbook.remove(workbook.active)
    for strength in sweep.STRENGTHS:
        sheet = workbook.create_sheet(strength)
        sheet.append(headers)
        sheet.append([None] + ['mg/L']*17 + ['m']*6 + ['kwh/day']*3)
        sheet.append(['GPS-X'] + [10.]*17 + [4e-4]*6 + ['=2*24', '=3*24', 96.])
        for cell in sheet[3]:
            cell.font = Font(bold=True, color='225577')
        for layers in sweep.LAYERS:
            sheet.append([f'QSDsan - {layers} biofilm layers'] + [999.]*26)
        sheet.column_dimensions['A'].width = 30
    path = tmp_path / 'template.xlsx'
    workbook.save(path)
    buffer = BytesIO()
    with ZipFile(path) as source, ZipFile(buffer, 'w', ZIP_DEFLATED) as output:
        for item in source.infolist():
            payload = source.read(item.filename)
            if item.filename.startswith('xl/worksheets/sheet'):
                root = ET.fromstring(payload)
                for coordinate, value in [('Y3', '48'), ('Z3', '72')]:
                    cell = root.find(f'.//{sweep._TAG("c")}[@r="{coordinate}"]')
                    cell.find(sweep._TAG('v')).text = value
                payload = ET.tostring(root, encoding='utf-8', xml_declaration=True)
            output.writestr(item, payload)
    path.write_bytes(buffer.getvalue())
    return path


@pytest.fixture
def fake_cases(monkeypatch):
    calls = []

    def stream(cod):
        return SimpleNamespace(
            iconc={'S_O2': 2., 'S_NH4': 0.5, 'S_NO3': 6., 'S_PO4': 1.5},
            BOD=20., COD=cod, TN=8., TP=2.,
        )

    def create(**kwargs):
        calls.append(kwargs)
        layers = kwargs['tf_kwargs']['N_biofilm_layer']
        tf = SimpleNamespace(ins=[stream(339.)], outs=[stream(100.+layers)],
                             total_biofilm_thickness=np.full(6, 4e-4))
        system = SimpleNamespace(
            flowsheet=SimpleNamespace(unit=SimpleNamespace(TF=tf)),
            scope=SimpleNamespace(sol=SimpleNamespace(
                success=True, t=np.array([0., 300.]), y=np.ones((2, 2)), message='ok')),
            get_influent_pumping_energy=lambda: 1.,
            get_recirculation_pumping_energy=lambda: 2.,
            get_rotating_arm_energy=lambda: 3.,
        )
        system.simulate = lambda **options: None
        return system

    create.__signature__ = signature(sweep.create_trickling_filter_validation_system)
    monkeypatch.setattr(sweep, 'create_trickling_filter_validation_system', create)
    return calls


def test_sweep_preserves_reference_and_plots_saved_workbook(template, fake_cases, tmp_path):
    original_bytes = template.read_bytes()
    source = load_workbook(template, data_only=False)
    run = sweep.run_validation_sweep(
        system_kwargs={'asm_kwargs': {'mu_H': 4.6}},
        template_path=template, output_dir=tmp_path, show_plot=False,
    )
    try:
        assert template.read_bytes() == original_bytes
        assert len(fake_cases) == 9
        assert {(c['wastewater_strength'], c['tf_kwargs']['N_biofilm_layer'])
                for c in fake_cases} == {(s, n) for s in sweep.STRENGTHS for n in sweep.LAYERS}
        assert len({id(c['flowsheet']) for c in fake_cases}) == 9
        assert all(c['asm_kwargs']['mu_H'] == 4.6 for c in fake_cases)
        assert all(c['asm_kwargs']['q_PHA'] == 6. for c in fake_cases)
        report = load_workbook(run.workbook_path, data_only=False)
        numeric = load_workbook(run.workbook_path, data_only=True)
        assert report.sheetnames == [*source.sheetnames, 'run_settings']
        for strength in sweep.STRENGTHS:
            for row in (1, 2, 3):
                assert [c.value for c in source[strength][row]] == [c.value for c in report[strength][row]]
                assert [copy(c._style) for c in source[strength][row]] == [copy(c._style) for c in report[strength][row]]
            assert numeric[strength]['Y3'].value == 48.
            assert numeric[strength]['Z3'].value == 72.
            assert numeric[strength]['G4'].value == 103.
            assert report[strength]['G4'].number_format == '0.0000'
            assert report[strength]['S4'].number_format == '0.0000000'
            assert report[strength]['Y4'].number_format == '0.00'
            assert numeric[strength]['L4'].value == 0.
            assert [numeric[strength].cell(4, c).value for c in range(25, 28)] == [24., 48., 72.]
            assert all(numeric[strength].cell(4, c).value == 4e-4 for c in range(19, 25))
        settings = {row[1]: row[2] for row in report['run_settings'].iter_rows(min_row=2, values_only=True)}
        assert settings['asm_kwargs.mu_H'] == 4.6
        assert settings['stop_at_steady_state'] is False
        assert settings['Archive'] == str(run.source_path)
        assert settings['Archive SHA256'] == sha256(run.source_path.read_bytes()).hexdigest()
        with ZipFile(run.source_path) as archive:
            manifest = json.loads(archive.read('manifest.json'))
            assert manifest['settings']['system_kwargs']['asm_kwargs']['mu_H'] == 4.6
            assert archive.read('inputs/01_template.xlsx') == original_bytes
        assert report['run_settings']['A1'].font.bold
        assert report['run_settings']['C2'].alignment.wrap_text
        assert not any('wastewater_strength' in key or 'N_biofilm_layer' in key for key in settings)
        assert run.variation.shape == (9, 15)
        assert run.variation.loc['L3', 'COD'] == pytest.approx((103.-10.)/103.)
        assert run.variation.loc['L3', 'E'] == pytest.approx((144.-216.)/216.)
        assert run.heatmap_path.is_file()
    finally:
        plt.close(run.figure)


@pytest.mark.parametrize('kwargs', [
    {'wastewater_strength': 'low'},
    {'tf_kwargs': {'N_biofilm_layer': 4}},
    {'tf_kwargs': {'N_section': 5}},
])
def test_sweep_rejects_settings_incompatible_with_case_layout(kwargs, template, fake_cases, tmp_path):
    with pytest.raises(ValueError):
        sweep.run_validation_sweep(kwargs, template_path=template, output_dir=tmp_path, show_plot=False)
    assert not fake_cases
    assert list(tmp_path.glob('TF_specific*.xlsx')) == []


def test_sweep_does_not_publish_a_failed_simulation(template, monkeypatch, tmp_path):
    def create(**kwargs):
        return SimpleNamespace(
            simulate=lambda **options: None,
            scope=SimpleNamespace(sol=SimpleNamespace(
                success=False, t=np.array([0., 10.]), y=np.ones((2, 2)), message='failed solve')),
        )
    create.__signature__ = signature(sweep.create_trickling_filter_validation_system)
    monkeypatch.setattr(sweep, 'create_trickling_filter_validation_system', create)
    with pytest.raises(RuntimeError, match='failed solve'):
        sweep.run_validation_sweep(template_path=template, output_dir=tmp_path, show_plot=False)
    assert list(tmp_path.glob('TF_specific*.xlsx')) == []


def test_excel_compatibility_prefixes_survive_serialization():
    source = b'''<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
        xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006"
        xmlns:x14ac="http://schemas.microsoft.com/office/spreadsheetml/2009/9/ac"
        mc:Ignorable="x14ac"><sheetData/></worksheet>'''
    payload = sweep._xml_bytes(ET.fromstring(source), source)
    namespaces = dict(item for _, item in ET.iterparse(BytesIO(payload), events=('start-ns',)))
    assert namespaces['x14ac'] == 'http://schemas.microsoft.com/office/spreadsheetml/2009/9/ac'
    assert ET.fromstring(payload).get(
        '{http://schemas.openxmlformats.org/markup-compatibility/2006}Ignorable') == 'x14ac'


def test_launcher_uses_current_calibration_profile():
    from exposan.esap.analyses.tf_validation_script import SYSTEM_KWARGS
    effective = sweep._prepare_kwargs(SYSTEM_KWARGS)
    assert effective['asm_kwargs'] == sweep.DEFAULT_SYSTEM_KWARGS['asm_kwargs']
    assert effective['asm_kwargs'] == {
        'mu_H': 3.3, 'q_PHA': 6.0, 'q_PP': 2.0, 'b_PP': 0.05,
    }
    assert effective['gpsx_kinetics'] is False
    assert effective['gpsx_decay'] is False
    assert effective['gpsx_aob_proxy'] is False
    assert effective['tf_kwargs']['detachment_rate'] == pytest.approx(0.07)
    cleared = sweep._prepare_kwargs({**SYSTEM_KWARGS, 'asm_kwargs': None})
    assert cleared['asm_kwargs'] is None


def test_launcher_selects_tighter_simulation_profile():
    from exposan.esap.analyses.tf_validation_script import SIMULATE_KWARGS

    assert SIMULATE_KWARGS == {
        'method': 'BDF', 't_span': (0, 300),
        'rtol': 1e-5, 'atol': 1e-6, 'max_step': 0.25,
        'stop_at_steady_state': False,
    }
    assert sweep.DEFAULT_SIMULATE_KWARGS['rtol'] == 1e-3
    assert sweep.DEFAULT_SIMULATE_KWARGS['atol'] == 1e-2


@pytest.mark.parametrize('layers', [3, 4, 5])
@pytest.mark.parametrize('requested_atol', [1e-2, 1e-6, 1e-10])
@pytest.mark.parametrize('vector_atol', [False, True])
def test_tf_tolerances_scale_after_state_assembly(
        layers, requested_atol, vector_atol):
    from exposan.esap.systems.TricklingFilter_validation import (
        _tf_tolerance_solver,
    )

    ids = ('S_F', 'S_O2', 'S_NH4', 'S_NO3', 'S_PO4',
           'X_PHA', 'X_PP', 'X_PAO', 'X_I')

    class Components(tuple):
        IDs = ids

        def indices(self, names):
            return np.array([self.index(ID) for ID in names])

    components = Components(ids)
    n, sections = len(ids), 6
    columns = n + 1
    liquid_size = sections * columns
    biofilm_size = sections * layers * n
    start = 2
    stop = start + liquid_size + biofilm_size + sections
    tf = SimpleNamespace(
        ID='TF', N_section=sections, N_biofilm_layer=layers,
        components=components, _ncol=columns,
    )
    system = SimpleNamespace(
        flowsheet=SimpleNamespace(unit=SimpleNamespace(TF=tf)),
        _state_idx={'TF': (start, stop)},
    )
    y0 = np.ones(stop + 3)
    original_atol = (np.full(y0.size, requested_atol)
                     if vector_atol else requested_atol)

    class CaptureSolver:
        def __init__(self, fun, t0, y0, t_bound, **options):
            self.y0 = y0
            self.atol = options['atol']

    solver_class = _tf_tolerance_solver(CaptureSolver, system)
    solver = solver_class(None, 0., y0, 1., atol=original_atol)

    expected = np.full(y0.size, requested_atol)
    selected = components.indices(ids[1:-1])
    for section in range(sections):
        expected[start + section * columns + selected] = min(requested_atol, 1e-6)
        for layer in range(layers):
            offset = start + liquid_size + (section * layers + layer) * n
            expected[offset + selected] = min(requested_atol, 1e-6)
    expected[stop-sections:stop] = min(requested_atol, 1e-9)

    np.testing.assert_array_equal(solver.atol, expected)
    np.testing.assert_array_equal(system.tf_effective_atol, expected[start:stop])
    np.testing.assert_array_equal(original_atol, np.full(y0.size, requested_atol)
                                  if vector_atol else requested_atol)
    np.testing.assert_array_equal(y0, np.ones(y0.size))
    assert not np.shares_memory(solver.y0, y0)
    if vector_atol:
        assert not np.shares_memory(solver.atol, original_atol)


# %% Standalone validation profiles, kinetics, energy and solver safeguards
def test_validation_bdf_inactive_column_perturbations():
    # from exposan.werf.systems.TricklingFilter_validation import _ValidationBDF
    from exposan.esap.systems.TricklingFilter_validation import _ValidationBDF
    solver = _ValidationBDF(
        lambda t, y: np.array([-y[0], 1.]),
        0., np.array([1., 0.]), 1., rtol=1e-3, atol=1e-2)
    # Simulate accumulated factor growth on a derivative-independent column.
    for factor in (1e300, np.inf, np.nan):
        solver.jac_factor[:] = factor
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            jac = solver.jac(0., solver.y)
        np.testing.assert_allclose(jac, [[-1., 0.], [0., 0.]], atol=1e-7)
        assert np.isfinite(solver.jac_factor).all()
    while solver.status == 'running':
        solver.step()
    assert solver.status == 'finished'
    np.testing.assert_allclose(solver.y, [np.exp(-1.), 1.], atol=0.02)


@pytest.mark.parametrize('oxygen,nitrate', [
    (0., 0.), (1e10, 0.), (0., 1e10), (0.2, 0.5), (2., 1.),
    ])
def test_validation_heterotrophic_decay_redox_limits(oxygen, nitrate):
    cmps = pc.create_masm2d_cmps()
    model = _create_validation_model(True, None, 0.6)
    reference = pc.mASM2d(
        electron_acceptor_dependent_decay=True,
        b_H=0.62, eta_NO3_Hl=0.9, f_XI_H=0.08, f_SI=0.)
    state = cmps.kwarray(default_as_init)
    state[cmps.index('S_O2')] = oxygen
    state[cmps.index('S_NO3')] = nitrate
    params = model.rate_function.params
    aerobic = oxygen / (params['K_O2_H'] + oxygen)
    nitrate_gate = nitrate / (params['K_NO3_H'] + nitrate)
    expected = 0.62 * state[cmps.index('X_H')] * (
        aerobic + 0.9*(1-aerobic)*nitrate_gate
        + 0.6*(1-aerobic)*(1-nitrate_gate))
    rates = model.rate_function(state).copy()
    reference_rates = reference.rate_function(state).copy()
    idx = model.index('hetero_lysis')
    np.testing.assert_allclose(rates[idx], expected, rtol=1e-12)
    other = np.arange(len(rates)) != idx
    # Hydrolysis, nitrification, PAO reactions and all other rates are unchanged.
    np.testing.assert_array_equal(rates[other], reference_rates[other])
    np.testing.assert_allclose(model.stoichio_eval(), reference.stoichio_eval())
    assert model.rate_function.params['k_h'] == reference.rate_function.params['k_h']
    assert model.rate_function.params['K_X'] == reference.rate_function.params['K_X']


def test_validation_decay_stoichiometry_and_legacy_mode():
    cmps = pc.create_masm2d_cmps()
    calibrated = _create_validation_model(True, None, 0.6)
    legacy = _create_validation_model(False, None, 0.6)
    stock = pc.mASM2d(electron_acceptor_dependent_decay=True)
    state = cmps.kwarray(default_as_init)
    for oxygen, nitrate in ((2., 1.), (0., 0.)):
        state[cmps.index('S_O2')] = oxygen
        state[cmps.index('S_NO3')] = nitrate
        expected = stock.rate_function(state).copy()
        calibrated.rate_function(state)
        np.testing.assert_array_equal(legacy.rate_function(state), expected)
    matrix = calibrated.stoichio_eval()
    idx = calibrated.index('hetero_lysis')
    row = matrix[idx]
    np.testing.assert_allclose(row[cmps.indices(('X_H', 'X_I', 'X_S'))],
                               [-1., 0.08, 0.92])
    for conversion in (cmps.i_COD, cmps.i_C, cmps.i_N, cmps.i_P):
        np.testing.assert_allclose(row @ conversion, 0., atol=1e-12)
    other = np.arange(len(matrix)) != idx
    # In particular, do not change the growth yields during this decay comparison.
    np.testing.assert_allclose(matrix[other], stock.stoichio_eval()[other])
    np.testing.assert_allclose(legacy.stoichio_eval(), stock.stoichio_eval())


def test_validation_decay_parameter_overrides():
    cmps = pc.create_masm2d_cmps()
    overrides = dict(b_H=0.7, eta_NO3_Hl=0.8, f_XI_H=0.12, k_h=4., K_X=0.2)
    model = _create_validation_model(True, overrides, 0.4)
    assert overrides == dict(b_H=0.7, eta_NO3_Hl=0.8, f_XI_H=0.12, k_h=4., K_X=0.2)
    state = cmps.kwarray(default_as_init)
    state[cmps.indices(('S_O2', 'S_NO3'))] = 0.
    idx = model.index('hetero_lysis')
    biomass = state[cmps.index('X_H')]
    np.testing.assert_allclose(model.rate_function(state)[idx], 0.7*0.4*biomass)
    # Parameter updates must affect the anaerobic addition as well as base kinetics.
    model.set_parameters(b_H=0.5)
    np.testing.assert_allclose(model.rate_function(state)[idx], 0.5*0.4*biomass)
    assert model.rate_function.params['k_h'] == 4.
    assert model.rate_function.params['K_X'] == 0.2


@pytest.mark.parametrize('factor', [-0.1, 1.1, float('nan'), float('inf')])
def test_validation_decay_rejects_invalid_factors(factor):
    pc.create_masm2d_cmps()
    with pytest.raises(ValueError, match='anaerobic_decay_factor'):
        _create_validation_model(True, None, factor)


def test_validation_decay_requires_redox_dependence():
    pc.create_masm2d_cmps()
    with pytest.raises(ValueError, match='electron_acceptor_dependent_decay'):
        _create_validation_model(
            True, dict(electron_acceptor_dependent_decay=False), 0.6)


@pytest.mark.parametrize('gpsx_kinetics', [False, True])
@pytest.mark.parametrize('gpsx_decay', [False, True])
# def test_validation_decay_system_simulation(gpsx_decay):
def test_validation_decay_system_simulation(gpsx_decay, gpsx_kinetics):
    system = create_trickling_filter_validation_system(
        gpsx_decay=gpsx_decay,
        gpsx_kinetics=gpsx_kinetics,
        tf_kwargs=dict(N_section=1, N_biofilm_layer=2, filter_bed_surface=9464.))
    tf = system.flowsheet.unit.TF
    # assert tf.detachment_exponent == 2.
    assert tf.detachment_exponent == 1.
    reference_feed = default_rww()
    np.testing.assert_array_equal(tf.ins[0].conc, reference_feed.conc)
    np.testing.assert_allclose(tf.ins[0].F_vol, reference_feed.F_vol)
    assert len(system.units) == 1
    np.testing.assert_allclose(tf.total_biofilm_thickness, 3e-4)
    assert tf.suspended_growth_model.rate_function.params['b_H'] == (
        0.62 if gpsx_decay else 0.4)
    system.simulate(method='BDF', t_span=(0., 0.02), rtol=1e-5, atol=1e-8)
    assert system.scope.sol.success, system.scope.sol.message
    assert np.isfinite(system.scope.sol.y).all()
    assert (tf.total_biofilm_thickness >= 0.).all()


@pytest.mark.parametrize(('strength', 'composition'), [
    ('low', dict(COD=339., NH4_N=14., PO4_P=1.6, S_K=11., S_Cl=39.)),
    ('mid', dict(COD=508., NH4_N=20., PO4_P=2.4, S_K=16., S_Cl=59.)),
    ('HIGH', dict(COD=1016., NH4_N=41., PO4_P=4.7, S_K=32., S_Cl=118.)),
    ])
def test_validation_wastewater_strength_profiles(strength, composition):
    system = create_trickling_filter_validation_system(
        wastewater_strength=strength, default_init_conds=False,
        tf_kwargs=dict(N_section=1, N_biofilm_layer=1))
    rww = system.flowsheet.stream.RWW
    reference = pc.create_masm2d_inf(
        'strength_reference', 10., 'MGD', T=293.15,
        **composition,
        fr_SI=0.075, fr_SF=0.24, fr_SA=0.05, fr_XI=0.3)

    assert system.flowsheet.unit.TF.ins[0] is rww
    np.testing.assert_allclose(rww.F_vol, reference.F_vol)
    np.testing.assert_allclose(rww.conc, reference.conc)
    assert rww.T == pytest.approx(293.15)


_TABLE_S1_NONDEFAULT_PROPERTIES = {
    'S_F': dict(i_N=0.023464, i_P=0.006149),
    'S_I': dict(i_N=0.042021, i_P=0.007139),
    'X_S': dict(i_N=0.023464, i_P=0.006149, i_mass=0.63),
    'X_I': dict(i_N=0.042021, i_P=0.007139, i_mass=0.63),
    'X_H': dict(i_N=0.060305, i_P=0.02369),
    'X_PAO': dict(i_N=0.060305, i_P=0.02369),
    'X_AUT': dict(i_N=0.060305, i_P=0.02369),
    }


@pytest.mark.parametrize(('strength', 'expected_bulk'), [
    ('low', dict(BOD=136.4, TN=23.9, TP=3.7, TSS=135.6)),
    # ('mid', dict(BOD=204.3, TN=35.6, TP=5.5, TSS=203.2)),
    ('mid', dict(BOD=204.3, TN=34.9, TP=5.6, TSS=203.2)),
    ('high', dict(BOD=408.7, TN=70.7, TP=11.0, TSS=406.5)),
    ])
def test_validation_nondefault_table_s1_properties(strength, expected_bulk):
    system = create_trickling_filter_validation_system(
        wastewater_strength=strength, default_init_conds=False,
        tf_kwargs=dict(N_section=1, N_biofilm_layer=1))
    tf = system.flowsheet.unit.TF
    cmps, rww = tf.components, tf.ins[0]

    for component_ID, expected in _TABLE_S1_NONDEFAULT_PROPERTIES.items():
        component = getattr(cmps, component_ID)
        for name, value in expected.items():
            assert getattr(component, name) == pytest.approx(value)

    observed_bulk = dict(
        BOD=rww.BOD, TN=rww.TN, TP=rww.TP, TSS=rww.get_TSS())
    for name, value in expected_bulk.items():
        assert observed_bulk[name] == pytest.approx(value, abs=0.05)


def test_validation_nondefault_component_properties_are_local():
    high = create_trickling_filter_validation_system(
        wastewater_strength='high', default_init_conds=False,
        tf_kwargs=dict(N_section=1, N_biofilm_layer=1))
    high_cmps = high.flowsheet.unit.TF.components
    default = create_trickling_filter_validation_system(
        wastewater_strength='default', default_init_conds=False,
        tf_kwargs=dict(N_section=1, N_biofilm_layer=1))
    default_cmps = default.flowsheet.unit.TF.components

    assert high_cmps is not default_cmps
    assert high_cmps.S_F.i_N == pytest.approx(0.023464)
    assert high_cmps.X_H.i_P == pytest.approx(0.02369)
    assert high_cmps.X_S.i_mass == pytest.approx(0.63)
    assert default_cmps.S_F.i_N == pytest.approx(0.03352)
    assert default_cmps.X_H.i_P == pytest.approx(0.02154)
    assert default_cmps.X_S.i_mass == pytest.approx(0.75)


def test_validation_trickling_filter_energy_metrics(monkeypatch):
    base = create_trickling_filter_validation_system(
        default_init_conds=False,
        tf_kwargs=dict(
            N_section=1, N_biofilm_layer=1,
            filter_bed_surface=9464., recirculation_ratio=0.))
    recycled = create_trickling_filter_validation_system(
        default_init_conds=False,
        tf_kwargs=dict(
            N_section=1, N_biofilm_layer=1,
            filter_bed_surface=9464., recirculation_ratio=1.))
    tf = recycled.flowsheet.unit.TF
    # base_power = base.get_pumping_energy()
    # recycled_power = recycled.get_pumping_energy()
    # assert base_power > 0.
    # assert recycled_power > base_power
    base_influent = base.get_influent_pumping_energy()
    recycled_influent = recycled.get_influent_pumping_energy()
    base_recirculation = base.get_recirculation_pumping_energy()
    recycled_recirculation = recycled.get_recirculation_pumping_energy()
    expected_influent = (
        1000. * 9.81 * (base.flowsheet.stream.RWW.F_vol/3600.)
        * ((50. + 3.28)*0.3048) / 0.60 / 1000.)
    assert base_influent == pytest.approx(expected_influent)
    assert recycled_influent == pytest.approx(base_influent)
    assert base_recirculation == 0.
    assert recycled_recirculation > 0.
    expected_recirculation = (
        1000. * 9.81 * (tf.ins[0].F_vol*tf.recirculation_ratio/3600.)
        * (tf.filter_bed_depth + (5.5 + 3.28)*0.3048)
        / 0.60 / 1000.)
    assert recycled_recirculation == pytest.approx(
        expected_recirculation)
    with monkeypatch.context() as patch:
        # patch.setattr(werf_models, 'TF_recirculation_pump_headloss', 0.)
        patch.setattr(esap_models, 'TF_recirculation_pump_headloss', 0.)
        no_pump_loss = recycled.get_recirculation_pumping_energy()
    assert recycled_recirculation > no_pump_loss
    assert base.get_pumping_energy() == pytest.approx(base_influent)
    assert recycled.get_pumping_energy() == pytest.approx(
        recycled_influent + recycled_recirculation)

    rotating_power = get_trickling_filter_rotating_arm_power(tf)
    assert rotating_power == pytest.approx(
        tf.filter_bed_surface * tf.filter_bed_depth * 3e-3)

    model = create_trickling_filter_validation_model(system=recycled)
    # values = {metric.name: metric() for metric in model.metrics}
    values = {metric.name.lower(): metric() for metric in model.metrics}
    # assert values['trickling filter pumping energy'] == pytest.approx(
    #     recycled_power)
    assert values['influent pumping energy'] == pytest.approx(
        recycled_influent)
    assert values['recirculation pumping energy'] == pytest.approx(
        recycled_recirculation)
    assert values['rotating arm energy consumption'] == pytest.approx(
        rotating_power)
    # assert values['trickling filter electricity cost'] == pytest.approx(
    #     (recycled_power + rotating_power) * 24.
    #     * recycled.power_utility.price)
    assert values['trickling filter electricity cost'] == pytest.approx(
        (recycled_influent + recycled_recirculation + rotating_power)
        * 24. * recycled.power_utility.price)
    assert recycled.get_energy_summary()['influent pumping energy [kW]'] \
        == pytest.approx(recycled_influent)

    with pytest.raises(TypeError, match='system_kwargs'):
        create_trickling_filter_validation_model(
            recycled, wastewater_strength='high')


@pytest.mark.parametrize('strength', ['medium', '', None])
def test_validation_rejects_unknown_wastewater_strength(strength):
    with pytest.raises(ValueError, match='wastewater_strength'):
        create_trickling_filter_validation_system(
            wastewater_strength=strength, default_init_conds=False,
            tf_kwargs=dict(N_section=1, N_biofilm_layer=1))


@pytest.mark.parametrize('gpsx_decay', [False, True])
def test_validation_direct_kinetic_matches(gpsx_decay):
    pc.create_masm2d_cmps()
    model = _create_validation_model(gpsx_decay, None, 0.6, gpsx_kinetics=True)
    params = model.rate_function.params
    expected = dict(
        mu_H=3.2, K_F=5., K_A_H=5., eta_NO3=0.28,
        q_PHA=6., q_PP=2., b_PP=0.05, eta_NO3_PAOl=0.9,
        b_AUT=0.17, eta_NO3_AUTl=0.5)
    for key, value in expected.items():
        assert params[key] == value, key

    # These directly comparable screenshot values already match mASM2d.
    matching = dict(
        k_h=3., K_X=0.1, eta_fe=0.4, K_NO3=0.5,
        K_O2_H=0.2, K_NH4_H=0.05, K_P_H=0.01, K_NO3_H=0.5,
        mu_PAO=1., K_A_PAO=4., K_PP=0.01, K_PHA=0.01, K_O2_PAO=0.2,
        K_MAX=0.34, K_IPP=0.02, K_P_S=0.2, b_PAO=0.2, b_PHA=0.2,
        K_NO3_PAO=0.5)
    for key, value in matching.items():
        assert params[key] == value, key

    # Do not impose separate-guild or individual NOx-step parameters on
    # lumped mASM2d processes; the PAO biomass factor is not a PP/PHA factor.
    untouched = dict(
        mu_AUT=1., K_O2_AUT=0.5, K_NH4_AUT=1.,
        eta_NO3_H=0.8, eta_NO3_PAO=0.6,
        eta_NO3_PPl=0.33, eta_NO3_PHAl=0.33)
    for key, value in untouched.items():
        assert params[key] == value, key
    assert params['mmp_kinetics'] is None
    assert params['b_H'] == (0.62 if gpsx_decay else 0.4)
    assert params['eta_NO3_Hl'] == (0.9 if gpsx_decay else 0.5)


def test_validation_kinetic_overrides_and_switches():
    cmps = pc.create_masm2d_cmps()
    overrides = dict(mu_H=4.2, eta_NO3=0.45, q_PHA=4.5, b_AUT=0.19)
    model = _create_validation_model(True, overrides, 0.6, gpsx_kinetics=True)
    for key, value in overrides.items():
        assert model.rate_function.params[key] == value
    assert overrides == dict(mu_H=4.2, eta_NO3=0.45, q_PHA=4.5, b_AUT=0.19)
    legacy = _create_validation_model(False, None, 0.6, gpsx_kinetics=False)
    stock = pc.mASM2d(electron_acceptor_dependent_decay=True)
    state = cmps.kwarray(default_as_init)
    expected = stock.rate_function(state).copy()
    np.testing.assert_array_equal(legacy.rate_function(state), expected)


def test_validation_kinetic_rate_effects():
    cmps = pc.create_masm2d_cmps()
    old = _create_validation_model(True, None, 0.6, gpsx_kinetics=False)
    new = _create_validation_model(True, None, 0.6, gpsx_kinetics=True)
    state = cmps.kwarray(default_as_init)
    state[cmps.index('S_O2')] = 0.
    state[cmps.index('S_NO3')] = 1.
    old_rates = old.rate_function(state).copy()
    new_rates = new.rate_function(state).copy()
    for process, ratio in (
            ('anox_hydrolysis', 0.28/0.6),
            ('storage_PHA', 6./3.),
            ('PAO_lysis', 0.9/0.33),
            ('PP_lysis', 0.05/0.2),
            ('auto_lysis', (0.17/0.15)*(0.5/0.33))):
        i = new.index(process)
        assert old_rates[i] > 0., process
        np.testing.assert_allclose(new_rates[i], old_rates[i]*ratio)
    # Kinetic changes must not alter stoichiometry or the existing decay wrapper.
    np.testing.assert_array_equal(new.stoichio_eval(), old.stoichio_eval())
    i = new.index('hetero_lysis')
    np.testing.assert_allclose(new_rates[i], old_rates[i])


@pytest.mark.parametrize('N,L,density,dry_fraction', [
    (6, 5, 1020000., 0.1), (2, 2, 900000., 0.08), (1, 1, 1100000., 0.12),
    ])
def test_validation_initial_biofilm_inventory(N, L, density, dry_fraction):
    from exposan.esap.systems.TricklingFilter_validation import asinit
    system = create_trickling_filter_validation_system(tf_kwargs=dict(
        N_section=N, N_biofilm_layer=L, biofilm_density=density,
        biofilm_dry_content=dry_fraction, initial_biofilm_thickness=3e-4))
    tf = system.flowsheet.unit.TF
    cmps = tf.components
    particulate = np.array([c.particle_size == 'Particulate' for c in cmps])
    target = density * dry_fraction
    init = tf._concs
    dry = init[:, 1:, particulate] @ cmps.i_mass[particulate]
    np.testing.assert_allclose(dry, target, rtol=1e-12)
    np.testing.assert_allclose(tf.total_biofilm_thickness, 3e-4)
    # Independent reconstruction of the borrowed profiles before normalization.
    section_factors = np.linspace(1., 0.45, N) if N > 1 else np.ones(N)
    for i, sf in enumerate(section_factors):
        base = cmps.kwarray(asinit.iloc[min(i, len(asinit)-1)].to_dict())
        liquid = base.copy()
        liquid[cmps.index('S_O2')] = 2.
        liquid[cmps.indices(('S_F', 'S_A', 'X_H'))] *= sf
        np.testing.assert_array_equal(init[i, 0], liquid)
        for k in range(L):
            depth_fraction = k/(L-1) if L > 1 else 0.
            depth_factor = max(0.25, 1.-0.55*depth_fraction)
            original = base.copy()
            original[cmps.index('S_O2')] = max(0.05, 0.35-0.30*depth_fraction)
            original[cmps.indices(('S_F', 'S_A'))] *= sf*depth_factor
            original[cmps.index('S_NH4')] *= max(sf, 0.4)
            original[cmps.index('X_H')] *= 1.5*sf*depth_factor
            original[cmps.index('X_AUT')] *= 0.8*depth_factor
            xh, pao, pp, pha = cmps.indices(('X_H', 'X_PAO', 'X_PP', 'X_PHA'))
            heterotrophs = original[xh] + original[pao]
            original[xh] = 0.99 * heterotrophs
            original[pao] = 0.01 * heterotrophs
            original[pp] = 0.5 * tf._model.rate_function.params['K_MAX'] * original[pao]
            original[pha] = 0.05 * original[pao]
            scale = target / (original[particulate] @ cmps.i_mass[particulate])
            np.testing.assert_allclose(init[i, k+1, particulate],
                                       original[particulate]*scale)
            np.testing.assert_array_equal(init[i, k+1, ~particulate],
                                          original[~particulate])
    # Initialize the actual ODE state too: inventory must equal density times volume.
    tf.ins[0]._init_state()
    tf._run()
    tf._init_state()
    stop = N*tf._ncol
    bio = tf._state[stop:stop+N*L*len(cmps)].reshape(N, L, len(cmps))
    H = tf._state[-N:]
    A = tf.filter_bed_surface*tf.filter_bed_depth*tf.media_specific_surface/N
    inventory = (bio[:, :, particulate] @ cmps.i_mass[particulate]).sum(axis=1) \
        * A*H/L
    np.testing.assert_allclose(inventory, target*A*H, rtol=1e-12)


def test_validation_initialization_without_biofilm():
    system = create_trickling_filter_validation_system(tf_kwargs=dict(
        # N_section=2, N_biofilm_layer=0))
        N_section=2, N_biofilm_layer=0, initial_biofilm_thickness=None))
    tf = system.flowsheet.unit.TF
    assert np.isfinite(tf._concs).all()


def test_validation_initialization_rejects_empty_solids():
    from exposan.esap.systems.TricklingFilter_validation import _set_tf_init_conc, asinit
    system = create_trickling_filter_validation_system(default_init_conds=False)
    tf = system.flowsheet.unit.TF
    with pytest.raises(ValueError, match='positive finite dry solids'):
        _set_tf_init_conc(tf, asinit*0.)


@pytest.mark.parametrize('thickness', [0., 1e-4, 3e-4, 9e-4, 9.99e-4])
def test_validation_detachment_flux_and_inventory(thickness):
    # # Isolate sloughing: no reactions, attachment, internal exchange or diffusion.
    # # This audits the current constant-velocity law, not GPS-X equivalence.
    # Isolate sloughing to verify that flux increases continuously with available
    # biofilm inventory and reaches the specified rate only at maximum thickness.
    # This avoids the old discontinuous near-zero detachment switch.
    system = create_trickling_filter_validation_system(tf_kwargs=dict(
        N_section=1, N_biofilm_layer=5, suspended_growth_model=None,
        initial_biofilm_thickness=thickness, attachment_rate=0.,
        internal_solids_exchange_rate=0., air_O2_transfer_coefficient=0.,
        dynamic_biofilm_thickness=True))
    tf = system.flowsheet.unit.TF
    cmps = tf.components
    particulate = np.array([c.particle_size == 'Particulate' for c in cmps])
    tf.diffusivity = np.zeros(len(cmps))
    tf.ins[0]._init_state()
    tf._run()
    tf._init_state()
    state = tf._state
    # Equal inlet/liquid concentrations remove the advective derivative.
    state[:len(cmps)] = tf.ins[0].state[:-1]
    tf.ODE(0., np.array([tf.ins[0].state]), state,
           np.zeros((1, len(cmps)+1)))
    deriv = tf._dstate
    A = tf.filter_bed_surface*tf.filter_bed_depth*tf.media_specific_surface
    liquid_gain = deriv[:len(cmps)][particulate] @ cmps.i_mass[particulate] \
        * tf.V_sections[0]
    # np.testing.assert_allclose(liquid_gain/A/1000., tf.detachment_rate, rtol=1e-12)
    # expected_flux = tf.detachment_rate \
    #     * thickness/tf.max_biofilm_thickness
    # expected_flux = tf.detachment_rate \
    #     * (thickness/tf.max_biofilm_thickness)**tf.detachment_exponent
    relative = thickness/tf.max_biofilm_thickness
    expected_flux = tf.detachment_rate \
        * (relative/max(1-relative, 1e-3))**tf.detachment_exponent
    np.testing.assert_allclose(liquid_gain/A/1000., expected_flux, rtol=1e-12, atol=1e-15)
    start = tf._ncol
    bio = state[start:-1].reshape(5, len(cmps))
    dbio = deriv[start:-1].reshape(5, len(cmps))
    # Moving-volume term is essential when checking mass, not just concentrations.
    bio_mass_change = ((dbio[:, particulate] @ cmps.i_mass[particulate])*thickness
                       + (bio[:, particulate] @ cmps.i_mass[particulate])*deriv[-1]) \
        .sum()*A/5
    np.testing.assert_allclose(bio_mass_change, -liquid_gain, rtol=1e-12)
    # np.testing.assert_allclose(deriv[-1],
    #     -tf.detachment_rate/(tf.biofilm_density*tf.biofilm_dry_content*1e-3))
    np.testing.assert_allclose(deriv[-1],
        -expected_flux/(tf.biofilm_density*tf.biofilm_dry_content*1e-3),
        atol=1e-15)


@pytest.mark.parametrize('exponent', [0., 0.5, np.nan, np.inf])
def test_validation_rejects_invalid_detachment_exponent(exponent):
    with pytest.raises(ValueError, match='detachment_exponent'):
        create_trickling_filter_validation_system(tf_kwargs=dict(
            N_section=1,
            N_biofilm_layer=1,
            detachment_exponent=exponent))


def test_validation_layer_density_during_simulation():
    system = create_trickling_filter_validation_system(tf_kwargs=dict(
        N_section=2, N_biofilm_layer=3, filter_bed_surface=9464.,
        recirculation_ratio=1., detachment_exponent=1.))
    tf = system.flowsheet.unit.TF
    cmps = tf.components
    p = np.array([c.particle_size == 'Particulate' for c in cmps])
    system.simulate(method='BDF', t_span=(0., 0.5), rtol=1e-6, atol=1e-9)
    sol = system.scope.sol
    assert sol.success, sol.message
    start = 2*tf._ncol
    bio = sol.y[start:start+2*3*len(cmps)].T.reshape(-1, 2, 3, len(cmps))
    dry_density = bio[:, :, :, p] @ cmps.i_mass[p]
    target = tf.biofilm_density*tf.biofilm_dry_content
    np.testing.assert_allclose(dry_density/target, 1., rtol=1e-5, atol=1e-7)
    assert (sol.y[-2:] > 0.).all()
    assert (sol.y[-2:] < tf.max_biofilm_thickness).all()


@pytest.mark.parametrize('law', ['linear', 'capacity'])
def test_detachment_capacity_flux_conserves_mass(law):
    system = create_trickling_filter_validation_system(tf_kwargs=dict(
        N_section=1, N_biofilm_layer=2, suspended_growth_model=None,
        detachment_model=law, initial_biofilm_thickness=5e-4,
        attachment_rate=0., internal_solids_exchange_rate=0.,
        air_O2_transfer_coefficient=0.))
    tf = system.flowsheet.unit.TF
    tf.diffusivity = np.zeros(len(tf.components))
    tf.ins[0]._init_state()
    tf._run()
    tf._init_state()
    cmps = tf.components
    p = np.array([c.particle_size == 'Particulate' for c in cmps])
    state = tf._state.copy()
    state[:len(cmps)] = tf.ins[0].state[:-1]
    tf.ODE(0., np.array([tf.ins[0].state]), state,
           np.zeros((1, len(cmps)+1)))
    A = tf.filter_bed_surface*tf.filter_bed_depth*tf.media_specific_surface
    gain = tf._dstate[:len(cmps)][p] @ cmps.i_mass[p]*tf.V_sections[0]
    expected = tf.detachment_rate*(1. if law == 'capacity' else 0.5)
    np.testing.assert_allclose(gain/A/1000., expected, rtol=1e-12)
    np.testing.assert_allclose(tf._dstate[-1],
                              -expected/(tf.biofilm_density*tf.biofilm_dry_content/1000.))


def test_steady_window_rejects_drift_and_matching_oscillatory_endpoints():
    # from exposan.werf.systems.TricklingFilter_validation import _SteadyStateWindow
    from exposan.esap.systems.TricklingFilter_validation import _SteadyStateWindow
    flat = _SteadyStateWindow(5., 0.01, np.array([1., 1e-6]))
    for t in range(5):
        assert not flat.update(t, [10., 5e-4])
    assert flat.update(5., [10., 5e-4])
    oscillating = _SteadyStateWindow(5., 0.01, np.array([1., 1e-6]))
    for t in range(6):
        # Same endpoints must not hide an intermediate thickness excursion.
        assert not oscillating.update(t, [10., 6e-4 if t == 3 else 5e-4])
    drifting = _SteadyStateWindow(5., 0.01, np.array([1., 1e-6]))
    for t in range(11):
        assert not drifting.update(t, [10., 5e-4*(1.+0.01*t)])


@pytest.mark.parametrize('stop', [True, False])
@pytest.mark.parametrize('method', ['BDF', 'LSODA'])
# def test_validation_steady_stop_preserves_endpoint_and_opt_out(stop, method):
@pytest.mark.parametrize('basis', ['states', 'inventory'])
def test_validation_steady_stop_preserves_endpoint_and_opt_out(stop, method, basis):
    system = create_trickling_filter_validation_system(
        default_init_conds=False, tf_kwargs=dict(
            N_section=1, N_biofilm_layer=0, initial_biofilm_thickness=None,
            suspended_growth_model=None, kLa=[0.]))
    tf = system.flowsheet.unit.TF
    # tf.set_init_conc(tf.ins[0].conc[None, :].copy())
    # system.simulate(method='BDF', t_span=(0., 10.), t_eval=[0., 10.],
    tf.set_init_conc(np.asarray(tf.ins[0].conc)[None, :].copy())
    system.simulate(method=method, t_span=(0., 10.), t_eval=[0., 10.],
                    steady_state_basis=basis,
                    stop_at_steady_state=stop, rtol=1e-3, atol=1e-2)
    sol = system.scope.sol
    assert sol.success
    # assert system.dynsim_kwargs['method'] == 'BDF'
    assert system.dynsim_kwargs['method'] == method
    assert system.steady_state_report['converged'] == stop
    if stop:
        assert 5. <= sol.t[-1] < 7.
        assert 'Steady state reached' in sol.message
    else:
        assert sol.t[-1] == 10.
    np.testing.assert_allclose(sol.y[:, -1], tf._state, atol=1e-8)


def test_aob_proxy_preserves_nitrogen_stoichiometry_and_overrides():
    pc.create_masm2d_cmps()
    base = _create_validation_model(True, None, 0.6, True)
    proxy = _create_validation_model(True, None, 0.6, True, True)
    params = proxy.rate_function.params
    assert (params['mu_AUT'], params['K_NH4_AUT'], params['K_O2_AUT']) == (0.9, 0.7, 0.25)
    np.testing.assert_array_equal(proxy.stoichio_eval(), base.stoichio_eval())
    overridden = _create_validation_model(True, {'K_O2_AUT': 0.4}, 0.6, True, True)
    assert overridden.rate_function.params['K_O2_AUT'] == 0.4
    # At low oxygen and the reported NH4, the proxy increases the growth gate.
    old_gate = 1.0*0.1/(0.5+0.1)*8.3/(1.0+8.3)
    new_gate = 0.9*0.1/(0.25+0.1)*8.3/(0.7+8.3)
    assert new_gate > old_gate


def test_steady_state_rejects_unknown_basis():
    with pytest.raises(ValueError, match='steady_state_basis'):
        create_trickling_filter_validation_system().dynamic_run(steady_state_basis='invalid')


# %% WERF D1 biology, feed compatibility and TF-specific integration
def test_d1_workbook_startup_preserves_previous_default():
    import qsdsan as qs
    from exposan.werf.systems import D1
    from exposan.werf._units import TricklingFilter
    from exposan.esap.systems.TricklingFilter_validation import _set_tf_init_conc, asinit

    system = D1.create_d1_system()
    tf = system.flowsheet.unit.TF
    reference = TricklingFilter(
        'reference', ins=qs.WasteStream('reference_feed'), outs='reference_effluent',
        N_section=6, N_biofilm_layer=5, initial_biofilm_thickness=4e-4,
        suspended_growth_model=tf.suspended_growth_model,
        biofilm_density=tf.biofilm_density, biofilm_dry_content=tf.biofilm_dry_content,
        )
    _set_tf_init_conc(reference, asinit)
    np.testing.assert_array_equal(tf._concs, reference._concs)
    np.testing.assert_array_equal(
        system.flowsheet.unit.AD._concs,
        system.flowsheet.unit.AD.components.kwarray(D1.dfs['adm'].loc['F1'].to_dict()))
    assert '_set_tf_init_conc' not in D1.__dict__
    assert D1.asinit.index[0] == (0, -1)
    assert D1.asinit.index[-1] == (5, 4)


@pytest.mark.parametrize('sections,layers', [(6, 5), (1, 1), (2, 3), (3, 4), (6, 0), (8, 7)])
def test_tf_workbook_profile_grid_averaging(sections, layers):
    import pandas as pd
    import qsdsan as qs
    from exposan.werf.systems.D1 import asinit
    from exposan.werf._units import TricklingFilter

    cmps = pc.create_masm2d_cmps()
    original = asinit.copy(deep=True)
    source = asinit.loc[:, list(cmps.IDs)].to_numpy().reshape(6, 6, len(cmps))
    tf = TricklingFilter(
        'TF', ins=qs.WasteStream('feed'), outs='effluent', N_section=sections,
        N_biofilm_layer=layers, initial_biofilm_thickness=4e-4 if layers else None,
        suspended_growth_model=None,
        )
    tf.set_init_conc(asinit)
    if layers:
        assert tf._concs.shape == (sections, layers+1, len(cmps))
        np.testing.assert_allclose(tf._concs[:, 0, :].mean(axis=0), source[:, 0, :].mean(axis=0),
                                   rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(tf._concs[:, 1:, :].mean(axis=(0, 1)), source[:, 1:, :].mean(axis=(0, 1)),
                                   rtol=1e-12, atol=1e-12)
        particulate = np.asarray([c.particle_size=='Particulate' for c in cmps])
        np.testing.assert_allclose(tf._concs[:, 1:, particulate] @ cmps.i_mass[particulate], 102000.)
        pao, pp, pha = cmps.indices(('X_PAO', 'X_PP', 'X_PHA'))
        np.testing.assert_allclose(tf._concs[:, 1:, pp]/tf._concs[:, 1:, pao], .17)
        np.testing.assert_allclose(tf._concs[:, 1:, pha]/tf._concs[:, 1:, pao], .05)
        if (sections, layers) == (6, 5):
            np.testing.assert_array_equal(tf._concs, source)
    else:
        assert tf._concs.shape == (sections, len(cmps))
        np.testing.assert_allclose(tf._concs.mean(axis=0), source[:, 0, :].mean(axis=0))
    pd.testing.assert_frame_equal(asinit, original)


@pytest.mark.parametrize('invalid', ['duplicate', 'missing_row', 'missing_liquid',
                                   'negative', 'nonfinite', 'empty', 'liquid_only', 'missing_component'])
def test_tf_workbook_profile_rejects_invalid_inputs(invalid):
    import pandas as pd
    import qsdsan as qs
    from exposan.werf.systems.D1 import asinit
    from exposan.werf._units import TricklingFilter

    pc.create_masm2d_cmps()
    frame = asinit.copy(deep=True)
    if invalid == 'duplicate':
        frame = pd.concat((frame, frame.iloc[[0]]))
    elif invalid == 'missing_row':
        frame = frame.drop(index=(1, 0))
    elif invalid == 'missing_liquid':
        frame = frame[frame.index.get_level_values(1) != -1]
    elif invalid == 'negative':
        frame.iloc[0, 0] = -1.
    elif invalid == 'nonfinite':
        frame.iloc[0, 0] = np.nan
    elif invalid == 'empty':
        frame = frame.iloc[:0]
    elif invalid == 'liquid_only':
        frame = frame[frame.index.get_level_values(1) == -1]
    elif invalid == 'missing_component':
        frame = frame.drop(columns='X_PAO')
    tf = TricklingFilter(
        'TF', ins=qs.WasteStream('feed'), outs='effluent', N_section=2,
        N_biofilm_layer=3, initial_biofilm_thickness=4e-4, suspended_growth_model=None,
        )
    with pytest.raises(KeyError if invalid == 'missing_component' else ValueError):
        tf.set_init_conc(frame)


def test_d1_tf_rates_match_svcw_and_base_decay():
    from exposan.esap.systems import create_svcw_system
    from exposan.werf.systems.D1 import create_d1_system

    d1 = create_d1_system(default_init_conds=False)
    tf = d1.flowsheet.unit.TF
    model = tf.suspended_growth_model
    svcw = create_svcw_system(default_init_conds=False)
    reference = svcw.flowsheet.unit.TF.suspended_growth_model
    parameters = model.rate_function.params
    reference_parameters = reference.rate_function.params
    assert parameters.keys() == reference_parameters.keys()
    for name, value in parameters.items():
        if isinstance(value, (bool, int, float, np.number)):
            assert value == pytest.approx(reference_parameters[name]), name

    # Component indices remain WERF-specific, so compare stoichiometry and
    # rates against base mASM2d compiled with D1's own components.
    stock = pc.mASM2d(components=tf.components,
                      electron_acceptor_dependent_decay=True)
    np.testing.assert_allclose(model.stoichio_eval(), stock.stoichio_eval())
    state = tf.components.kwarray(default_as_init)
    for oxygen, nitrate in ((2., 1.), (0., 1.), (0., 0.)):
        state[tf.components.indices(('S_O2', 'S_NO3'))] = oxygen, nitrate
        expected = stock.rate_function(state).copy()
        actual = model.rate_function(state).copy()
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    assert actual[model.index('hetero_lysis')] == pytest.approx(0.)


def test_d1_preserves_f1_influent_and_components():
    from exposan.werf.systems.D1 import create_d1_system
    from exposan.werf.systems.F1 import create_f1_system

    d1 = create_d1_system(default_init_conds=False)
    f1 = create_f1_system(default_init_conds=False)
    d1_components = d1.flowsheet.unit.TF.components
    f1_components = f1.flowsheet.unit.ASR.components
    assert d1_components.IDs == f1_components.IDs
    for attribute in ('i_N', 'i_P', 'i_C', 'i_COD', 'i_mass'):
        np.testing.assert_array_equal(
            getattr(d1_components, attribute), getattr(f1_components, attribute))
    for component_id in ('X_H', 'X_PAO', 'X_AUT'):
        assert getattr(d1_components, component_id).i_P == pytest.approx(0.02154)
        assert getattr(d1_components, component_id).i_N == pytest.approx(0.08615)

    rww = d1.flowsheet.stream.RWW
    reference = f1.flowsheet.stream.RWW
    np.testing.assert_array_equal(rww.conc, reference.conc)
    assert rww.F_vol == pytest.approx(reference.F_vol)
    assert rww.T == pytest.approx(reference.T)
    assert rww.get_total_flow('MGD') == pytest.approx(10.)
    assert rww.COD == pytest.approx(358.)
    assert rww.iconc['S_NH4'] == pytest.approx(25.91)
    assert rww.iconc['S_PO4'] == pytest.approx(5.)
    assert d1.flowsheet.unit.TF.N_biofilm_layer == 5


@pytest.mark.parametrize('initialize', [True, False])
def test_d1_final_clarifier_matches_f1_without_ras(initialize):
    from qsdsan import unit_operations as su
    from exposan.werf.systems.D1 import create_d1_system, MGD2cmd
    from exposan.werf.systems.F1 import create_f1_system

    d1 = create_d1_system(default_init_conds=initialize)
    f1 = create_f1_system(default_init_conds=initialize)
    units = d1.flowsheet.unit
    actual, reference = units.FC, f1.flowsheet.unit.FC
    assert type(actual) is type(reference) is su.FlatBottomCircularClarifier
    for name in ('A_settle', 'height', 'V_settle', 'h_layer', 'N_layer',
                 'feed_layer', 'X_t', 'v_max', 'v_max_p', 'rh', 'rp', 'fns',
                 'maximum_nonsettleable_solids', 'wastage'):
        assert getattr(actual, name) == pytest.approx(getattr(reference, name)), name
    assert actual.underflow == 0.
    assert reference.underflow == pytest.approx(6.7*MGD2cmd)
    assert actual.wastage == pytest.approx(0.163*MGD2cmd)
    assert actual.ins[0] is units.TF.outs[0]
    assert actual.outs[0] is d1.flowsheet.stream.SE
    assert actual.outs[1].sink is None
    assert actual.outs[2] is d1.flowsheet.stream.WAS is units.GT.ins[1]
    assert 'RAS' not in d1.flowsheet.stream
    for name in ('_solubles', '_solids', '_X_comp'):
        np.testing.assert_array_equal(getattr(actual, name), getattr(reference, name))


@pytest.mark.parametrize('method', ['BDF', 'LSODA'])
def test_d1_applies_tf_specific_tolerances(method):
    from exposan.werf.systems.D1 import create_d1_system

    system = create_d1_system(tf_kwargs={'N_section': 1, 'N_biofilm_layer': 1})
    system.simulate(method=method, t_span=(0., 1e-5), max_step=2e-6,
                    rtol=1e-3, atol=1e-2)
    solution = system.scope.sol
    assert solution.success, solution.message
    assert solution.t[-1] == pytest.approx(1e-5)
    assert np.isfinite(solution.y).all()
    assert system.dynsim_kwargs['method'] == method

    tf = system.flowsheet.unit.TF
    start, stop = system._state_idx[tf.ID]
    expected = np.full(stop-start, 1e-2)
    indices = np.asarray(tf.components.indices((
        'S_O2', 'S_NH4', 'S_NO3', 'S_PO4', 'X_PHA', 'X_PP', 'X_PAO')))
    expected[indices] = 1e-6
    expected[tf._ncol + indices] = 1e-6
    expected[-1] = 1e-9
    np.testing.assert_array_equal(system.tf_effective_atol, expected)

    # The dynamic settler stores solids; flux imbalance equals inventory change.
    fc = system.flowsheet.unit.FC
    feed = fc.ins[0].state.copy()
    fc.ODE(solution.t[-1], feed[None, :], fc._state.copy(), np.zeros((1, feed.size)))
    tss_in = feed[:-1] @ (fc.components.x*fc.components.i_mass)
    net_in = (feed[-1]*tss_in - (feed[-1]-fc.wastage)*fc._state[-fc.N_layer]
              - fc.wastage*fc._state[-1])
    accumulation = fc.V_settle/fc.N_layer*fc._dstate[-fc.N_layer:].sum()
    assert accumulation == pytest.approx(net_in, rel=1e-10, abs=1e-6)
    assert fc.outs[1].state[-1] == 0.


def test_d1_restores_solver_method_after_failure(monkeypatch):
    import qsdsan as qs
    from exposan.werf.systems.D1 import create_d1_system, _ValidationBDF

    system = create_d1_system(default_init_conds=False)
    captured = []

    def fail(self, **options):
        captured.append(options['method'])
        raise RuntimeError('test solver failure')

    monkeypatch.setattr(qs.System, 'dynamic_run', fail)
    with pytest.raises(RuntimeError, match='test solver failure'):
        system.dynamic_run(method='BDF', t_span=(0., 1.))
    assert issubclass(captured[0], _ValidationBDF)
    assert system.dynsim_kwargs['method'] == 'BDF'


def test_d1_preserves_tf_overrides_and_shared_solver_imports():
    from exposan.werf.systems import D1
    from exposan.esap.systems import TricklingFilter_validation as validation

    assert validation._ValidationBDF is D1._ValidationBDF
    assert validation._tf_tolerance_solver is D1._tf_tolerance_solver
    options = {'N_section': 2, 'N_biofilm_layer': 3, 'detachment_rate': 0.08}
    system = D1.create_d1_system(default_init_conds=False, tf_kwargs=options)
    tf = system.flowsheet.unit.TF
    assert tf.N_section == 2
    assert tf.N_biofilm_layer == 3
    assert tf.detachment_rate == pytest.approx(0.08)
    assert options == {'N_section': 2, 'N_biofilm_layer': 3, 'detachment_rate': 0.08}


# %% WERF D1E CHP wrapper and accounting
@pytest.mark.parametrize('initialize', [True, False])
def test_d1e_uses_shared_chp_wrapper(monkeypatch, initialize):
    from exposan.werf.systems import D1, D1E

    flowsheet, expected = object(), object()
    calls = []

    def create(*args):
        calls.append(args)
        return expected

    monkeypatch.setattr(D1E, 'create_chp_system', create)
    actual = D1E.create_d1e_system(flowsheet, default_init_conds=initialize)
    assert actual is expected
    assert calls == [('D1E', D1.create_d1_system, flowsheet, initialize)]


def test_d1e_registry_preserves_d1_configuration_and_initialization():
    from exposan.werf import systems

    assert systems.SYSTEM_CREATORS['D1E'] is systems.D1E.create_d1e_system
    assert systems.create_d1e_system is systems.D1E.create_d1e_system
    reference = systems.create_d1_system()
    system = systems.create_system('d1e')
    assert system.ID == system.flowsheet.ID == 'D1E'
    assert type(system) is type(reference)
    assert tuple(unit.ID for unit in system.path) == tuple(
        unit.ID for unit in reference.path)
    np.testing.assert_array_equal(
        system.flowsheet.stream.RWW.conc, reference.flowsheet.stream.RWW.conc)
    for unit_id in ('TF', 'AD'):
        actual = getattr(system.flowsheet.unit, unit_id)
        expected = getattr(reference.flowsheet.unit, unit_id)
        np.testing.assert_array_equal(actual._concs, expected._concs)
    for name in ('_solubles', '_solids', '_X_comp'):
        np.testing.assert_array_equal(
            getattr(system.flowsheet.unit.FC, name),
            getattr(reference.flowsheet.unit.FC, name))
    assert system.flowsheet.unit.FC.underflow == 0.
    assert system.flowsheet.unit.FC.wastage == reference.flowsheet.unit.FC.wastage
    tf = system.flowsheet.unit.TF
    expected_tf = reference.flowsheet.unit.TF
    assert tf.N_biofilm_layer == expected_tf.N_biofilm_layer == 5
    np.testing.assert_array_equal(
        tf.initial_biofilm_thickness, expected_tf.initial_biofilm_thickness)
    for name, value in tf.suspended_growth_model.rate_function.params.items():
        if isinstance(value, (bool, int, float, np.number)):
            assert value == pytest.approx(
                expected_tf.suspended_growth_model.rate_function.params[name])


def test_d1e_enables_chp_recovery_and_opex_credits():
    import qsdsan as qs
    from exposan.werf import models, systems

    for configuration in ('D1', 'D1E'):
        system = systems.create_system(configuration)
        model = qs.Model(system)
        models.add_OPEX_metrics(model)
        metrics = {metric.name.lower(): metric for metric in model.metrics}
        assert 'ad heating demand' in metrics
        assert 'ad heating cost' in metrics
        chp_names = {name for name in metrics if name.startswith('chp ')}
        if configuration == 'D1':
            assert not models._has_CHP(system)
            assert not chp_names
            continue

        assert models._has_CHP(system)
        assert chp_names == {
            'chp electricity recovery', 'chp heat recovery', 'chp total recovery',
            'chp electricity credit', 'chp heat credit', 'chp total credit',
            }
        biogas = system.flowsheet.stream.biogas
        biogas.imass['S_ch4'] = 100. / biogas.components.S_ch4.i_mass
        electricity, heat = models.get_CHP_outputs(
            100., system.flowsheet.stream.RWW.F_vol*24.)
        flow = system.flowsheet.stream.RWW.F_vol
        expected = {
            'chp electricity recovery': -electricity*flow,
            'chp heat recovery': -heat*flow,
            'chp total recovery': -(electricity+heat)*flow,
            'chp electricity credit': -electricity*flow*24.*system.power_utility.price,
            'chp heat credit': -heat*flow*24.*models.heat_price,
            }
        expected['chp total credit'] = (
            expected['chp electricity credit'] + expected['chp heat credit'])
        for name, value in expected.items():
            assert metrics[name]() == pytest.approx(value), name
        assert metrics['total opex']() == pytest.approx(expected['chp total credit'])


# %% Baseline-performance extension workflow and exports
@pytest.fixture
def cached_state_system():
    class Unit:
        hasode = True

        def __init__(self, ID, size):
            self.ID = ID
            self.size = size
            self.updates = 0

        def _init_dynamic(self):
            self._state = self._dstate = None

        def _update_state(self):
            assert self._state.shape == (self.size,)
            self.updates += 1

        def _update_dstate(self):
            self._dstate = np.zeros_like(self._state)

    units = (Unit('first', 2), Unit('second', 3))
    system = SimpleNamespace(ID='D1', units=units, convergences=0,
                             _state=np.full(5, -1.), _state_idx=None, _DAE=object())

    def converge():
        system.convergences += 1

    def assemble():
        system._state = np.arange(5, dtype=float)
        system._state_idx = {'first': (0, 2), 'second': (2, 5)}
        for unit in units:
            start, stop = system._state_idx[unit.ID]
            unit._state = system._state[start:stop]
        return system._state, system._state_idx, 0

    system.converge = converge
    system._load_state = assemble
    return system


def test_cached_state_load_preserves_views_and_input(cached_state_system):
    from exposan.werf.utils import load_state

    system = cached_state_system
    saved = np.arange(5, dtype=float) + 10.
    original = saved.copy()
    assert load_state(system, state_arr=saved) is True
    np.testing.assert_array_equal(system._state, original)
    assert system._DAE is None
    assert not np.shares_memory(system._state, saved)
    for unit in system.units:
        assert np.shares_memory(unit._state, system._state)
        assert unit.updates == 1
    saved[:] = -100.
    np.testing.assert_array_equal(system._state, original)


@pytest.mark.parametrize('strict', [True, False])
@pytest.mark.parametrize('saved', [np.ones(4), np.ones((1, 5)),
                                  np.array([1., 2., np.nan, 4., 5.]),
                                  np.array([1., 2., np.inf, 4., 5.])])
def test_cached_state_rejects_incompatible_before_copy(cached_state_system, strict, saved):
    from exposan.werf.utils import load_state

    system = cached_state_system
    if strict:
        with pytest.raises(ValueError, match='cached state'):
            load_state(system, state_arr=saved, strict=True)
    else:
        with pytest.warns(RuntimeWarning, match='Using fresh initialization'):
            assert load_state(system, state_arr=saved, strict=False) is False
    if system.convergences:
        np.testing.assert_array_equal(system._state, np.arange(5))
        for unit in system.units:
            assert unit._state.shape == (unit.size,)
            assert np.shares_memory(unit._state, system._state)


@pytest.mark.parametrize('strict', [True, False])
def test_cached_state_missing_file(cached_state_system, strict, tmp_path):
    from exposan.werf.utils import load_state

    system = cached_state_system
    if strict:
        with pytest.raises(FileNotFoundError):
            load_state(system, folder=str(tmp_path))
    else:
        with pytest.warns(RuntimeWarning, match='state cache not found'):
            assert load_state(system, folder=str(tmp_path), strict=False) is False
        np.testing.assert_array_equal(system._state, np.arange(5))


@pytest.mark.parametrize('exact', [True, False])
def test_cached_state_chp_filename_precedence(cached_state_system, tmp_path, exact):
    from exposan.werf.utils import load_state

    np.save(tmp_path/'D1.npy', np.ones(5))
    if exact:
        np.save(tmp_path/'D1E.npy', np.full(5, 2.))
    system = cached_state_system
    system.ID = 'D1E'
    assert load_state(system, folder=str(tmp_path))
    np.testing.assert_array_equal(system._state, np.full(5, 2. if exact else 1.))


def test_cached_state_does_not_hide_initialization_errors(cached_state_system):
    from exposan.werf.utils import load_state

    def fail():
        raise RuntimeError('model initialization failed')

    cached_state_system.converge = fail
    with pytest.raises(RuntimeError, match='model initialization failed'):
        load_state(cached_state_system, state_arr=np.ones(4), strict=False)


def test_d1_cached_state_stale_layout_leaves_valid_initialization():
    from exposan.werf import create_system
    from exposan.werf.utils import load_state

    system = create_system('D1')
    with pytest.warns(RuntimeWarning, match='1206 entries.*1248'):
        assert load_state(system, state_arr=np.zeros(1206), strict=False) is False
    assert system._state.shape == (1248,)
    assert system._state_idx == {
        'TF': (0, 1128), 'FC': (1128, 1170),
        'AD': (1170, 1216), 'HD': (1216, 1248),
        }
    for unit in system.units:
        if unit.hasode:
            assert np.shares_memory(unit._state, system._state)
    assert system.flowsheet.unit.AD._state.shape == (46,)
    assert np.isfinite(system._state).all()


def test_d1_cached_state_roundtrip_and_restart(tmp_path):
    from exposan.werf.systems.D1 import create_d1_system
    from exposan.werf.utils import cache_state, load_state

    system = create_d1_system(tf_kwargs={'N_section': 1, 'N_biofilm_layer': 1})
    options = dict(method='BDF', t_span=(0., 1e-5), max_step=2e-6)
    system.simulate(**options)
    assert system.scope.sol.success
    saved = system._state.copy()
    cache_state(system, folder=str(tmp_path))
    assert load_state(system, folder=str(tmp_path))
    material = np.ones(saved.size, dtype=bool)
    for unit in system.units:
        if not unit.hasode:
            continue
        start, stop = system._state_idx[unit.ID]
        if unit.ID == 'TF':
            material[start+unit._ncol-1:start+unit.N_section*unit._ncol:unit._ncol] = False
        elif unit.ID == 'FC':
            material[start+len(unit.components)] = False
        else:
            material[stop-1] = False
        assert np.shares_memory(unit._state, system._state)
    np.testing.assert_array_equal(system._state[material], saved[material])
    np.testing.assert_array_equal(np.load(tmp_path/'D1.npy'), saved)
    system.simulate(**options)
    assert system.scope.sol.success
    assert system.scope.sol.t[-1] == pytest.approx(1e-5)


@pytest.mark.parametrize('registered', [True, False])
def test_baseline_operating_cell_handles_unregistered_underflows(registered):
    import ast
    from pathlib import Path

    source = Path(__file__).resolve().parents[1]/'exposan/werf/analyses/baseline_optimize_extension.py'
    tree = ast.parse(source.read_text())
    cell = next(n for n in tree.body if isinstance(n, ast.Try)
                and any(isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                        and node.func.id == 'load_state' for node in ast.walk(n)))
    ID = 'F1' if registered else 'D1'
    thickener, dw = SimpleNamespace(sludge_flow_rate=10.), SimpleNamespace(sludge_flow_rate=2.)
    system = SimpleNamespace(ID=ID)
    calls = []

    def load(actual, **options):
        assert actual is system
        assert options == {'folder': 'steady_states/baseline_unopt', 'strict': False}
        return False

    def simulate(actual, options, thickened, actual_thickener, target):
        assert actual is system and actual_thickener is thickener
        assert options == {'t_span': (0, 300), 'method': 'BDF'}
        assert target == 18e4
        calls.append((thickener.sludge_flow_rate, dw.sludge_flow_rate))

    clock = SimpleNamespace(time=lambda: 0., strftime=lambda *args: '',
                             localtime=lambda: None, gmtime=lambda *args: None)
    env = dict(sys=system, ID=ID, tm=clock, load_state=load,
               baseline_underflows={'F1': (62.5, 19.35)},
               thickener=thickener, u=SimpleNamespace(DW=dw),
               thickened=object(), cake_tss=18e4, _simulate_to_tss=simulate,
               mdl=object(), display_metrics=lambda model: None)
    exec(compile(ast.Module(body=[cell], type_ignores=[]), str(source), 'exec'), env)
    assert calls == [(62.5, 19.35) if registered else (10., 2.)]


@pytest.mark.parametrize('configuration', ['D1', 'D1E'])
def test_d1_restarts_after_underflow_adjustment(configuration):
    from exposan.werf import create_system

    system = create_system(configuration)
    options = dict(method='BDF', t_span=(0., 1e-5), max_step=2e-6)
    system.simulate(**options)
    assert system.scope.sol.success
    scopes = [(subject, subject.scope) for subject in system.scope.subjects]
    units = system.flowsheet.unit
    units.GT.sludge_flow_rate *= 0.699
    units.DW.sludge_flow_rate *= 0.706
    system.simulate(**options)
    assert system.scope.sol.success
    assert system.scope.sol.t[-1] == pytest.approx(options['t_span'][1])
    for subject, tracker in scopes:
        assert subject.scope is tracker
        assert len(tracker.time_series) == len(system.scope.time_series)
    system.simulate(**options, state_reset_hook='reset_cache')
    assert system.scope.sol.success
    assert system.scope.sol.t[-1] == pytest.approx(options['t_span'][1])


def test_baseline_extension_reports_initial_and_retry_errors():
    from exposan.werf.analyses import baseline_performance_extension as extension

    initial = ValueError('original failure')
    retry = AttributeError('missing tracker')
    errors = iter((initial, retry))
    calls = []

    def simulate(**options):
        calls.append(options)
        raise next(errors)

    system = SimpleNamespace(ID='D1', simulate=simulate)
    with pytest.raises(RuntimeError) as caught:
        extension._simulate(system, extension.SIMULATE_KWARGS)
    assert 'initial simulation failed (ValueError: original failure)' in str(caught.value)
    assert 'cache-reset retry failed (AttributeError: missing tracker)' in str(caught.value)
    assert caught.value.__cause__ is retry
    assert 'state_reset_hook' not in calls[0]
    assert calls[1]['state_reset_hook'] == 'reset_cache'


@pytest.fixture
def extension_cases(monkeypatch):
    import pandas as pd
    from exposan.werf.analyses import baseline_performance_extension as extension

    class Registry(dict):
        __getattr__ = dict.__getitem__

    context = SimpleNamespace(systems=[], calls=[], cleared=[], failed=set())

    def create(ID):
        thickener = SimpleNamespace(sludge_flow_rate=10.)
        dw = SimpleNamespace(sludge_flow_rate=2.)
        streams = Registry(
            PS=object(),
            thickened_sludge=SimpleNamespace(
                get_TSS=lambda: 1e5*10./thickener.sludge_flow_rate),
            cake=SimpleNamespace(get_TSS=lambda: 54e4*2./dw.sludge_flow_rate),
            )
        units = Registry(GT=thickener, DW=dw, AD=SimpleNamespace(
            ins=[SimpleNamespace(F_mass=1000., T=293.15)]))
        system = SimpleNamespace(
            ID=ID, _state=np.array([1., 2.]),
            flowsheet=SimpleNamespace(stream=streams, unit=units,
                                      clear=lambda: context.cleared.append(ID)),
            scope=SimpleNamespace(sol=SimpleNamespace(
                success=True, status=0, message='ok', y=np.ones((2, 2)))),
            )

        def simulate(**options):
            context.calls.append((ID, options.copy()))
            if ID in context.failed:
                raise RuntimeError('test simulation failure')

        system.simulate = simulate
        context.systems.append(system)
        return system

    def add_metric(model, element, name, value):
        def metric():
            return value
        metric.column = element, name
        model.metrics.append(metric)

    def opex(model):
        add_metric(model, 'OPEX', 'total OPEX [USD/d]', 20.)
        if model.system.ID.endswith('E'):
            add_metric(model, 'CHP', 'electricity [kW]', -5.)

    monkeypatch.setattr(extension, 'create_system', create)
    monkeypatch.setattr(extension, 'Model',
                        lambda system: SimpleNamespace(system=system, metrics=[]))
    monkeypatch.setattr(extension, 'add_performance_metrics',
                        lambda model: add_metric(model, 'Effluent', 'COD [mg/L]', 10.))
    monkeypatch.setattr(extension, 'add_OPEX_metrics', opex)
    monkeypatch.setattr(extension, 'var_columns',
                        lambda metrics: pd.MultiIndex.from_tuples(
                            [metric.column for metric in metrics]))
    monkeypatch.setattr(extension, 'plantwide_N_mass_flows',
                        lambda system: pd.DataFrame({'TN': [1.]}))
    monkeypatch.setattr(extension, 'plantwide_P_mass_flows',
                        lambda system: pd.DataFrame({'TP': [2.]}))
    return extension, context


def test_baseline_extension_adjusts_underflows_and_preserves_original_outputs(
        extension_cases, tmp_path):
    import pandas as pd

    extension, context = extension_cases
    expected = (
        'B1', 'B1E', 'B2', 'B3', 'C1', 'C1E', 'C2', 'C3', 'D1', 'D1E',
        'E2', 'E2P', 'F1', 'F1E', 'G1', 'G1E', 'G2', 'G3', 'H1', 'H1E',
        'I1', 'I1E', 'I2', 'I3', 'N1', 'N1E', 'N2',
        )
    assert extension.BASELINE_CONFIGURATIONS == expected
    original = tmp_path/'_baseline_unopt_performance.xlsx'
    original.write_bytes(b'original baseline results')
    results = extension.run_baseline_performance_extension(output_dir=tmp_path)
    assert tuple(results.index) == expected
    assert results.attrs['failures'] == {}
    assert np.isnan(results.loc['D1', ('CHP', 'electricity [kW]')])
    assert results.loc['D1E', ('CHP', 'electricity [kW]')] == -5.
    assert original.read_bytes() == b'original baseline results'
    assert [ID for ID, options in context.calls] == [ID for ID in expected for _ in range(2)]
    assert tuple(context.cleared) == expected
    for system in context.systems:
        assert system.flowsheet.unit.GT.sludge_flow_rate == pytest.approx(20.)
        assert system.flowsheet.unit.DW.sludge_flow_rate == pytest.approx(6.)
        np.testing.assert_array_equal(
            np.load(tmp_path/'steady_states'/'baseline_unopt'/f'{system.ID}.npy'),
            system._state)
    for filename in ('N_mass_unopt_extension.xlsx', 'P_mass_unopt_extension.xlsx'):
        with pd.ExcelFile(tmp_path/filename) as workbook:
            assert tuple(workbook.sheet_names) == expected
    saved = pd.read_excel(tmp_path/'_baseline_unopt_performance_extension.xlsx',
                          header=[0, 1], index_col=0)
    pd.testing.assert_frame_equal(saved, results, check_dtype=False)


def test_baseline_extension_supports_future_cases_and_solver_overrides(
        extension_cases, monkeypatch, tmp_path):
    extension, context = extension_cases
    monkeypatch.setattr(extension, 'SYSTEM_CREATORS', {'FUTURE': object()})
    options = dict(t_span=(0, 400), rtol=1e-5, atol=1e-6)
    results = extension.run_baseline_performance_extension(
        'future', simulate_kwargs=options, output_dir=tmp_path)
    assert list(results.index) == ['FUTURE']
    assert context.calls[0][1] == {'method': 'BDF', **options}
    assert options == dict(t_span=(0, 400), rtol=1e-5, atol=1e-6)
    assert extension.SIMULATE_KWARGS == dict(t_span=(0, 300), method='BDF')


def test_baseline_extension_omits_failed_cases(extension_cases, tmp_path):
    import pandas as pd

    extension, context = extension_cases
    context.failed.add('D1')
    results = extension.run_baseline_performance_extension(
        ('D1', 'D1E'), output_dir=tmp_path)
    assert list(results.index) == ['D1E']
    failure = results.attrs['failures']['D1']
    assert 'initial simulation failed (RuntimeError: test simulation failure)' in failure
    assert 'cache-reset retry failed (RuntimeError: test simulation failure)' in failure
    assert context.calls[1][1]['state_reset_hook'] == 'reset_cache'
    assert context.cleared == ['D1', 'D1E']
    assert not (tmp_path/'steady_states'/'baseline_unopt'/'D1.npy').exists()
    with pd.ExcelFile(tmp_path/'N_mass_unopt_extension.xlsx') as workbook:
        assert workbook.sheet_names == ['D1E']


def test_baseline_extension_bounds_tss_adjustment(extension_cases, tmp_path):
    extension, context = extension_cases
    with pytest.raises(RuntimeError, match='TSS targets not reached'):
        extension.run_baseline_performance_extension(
            ('D1', 'D1E'), output_dir=tmp_path, max_tss_adjustments=0)
    assert context.cleared == ['D1', 'D1E']
    assert len(context.calls) == 2
    assert not list(tmp_path.rglob('*.xlsx'))
    assert not list(tmp_path.rglob('*.npy'))


@pytest.mark.parametrize('configurations', [(), ('D1', 'd1'), ('unknown',)])
def test_baseline_extension_rejects_invalid_selection(
        extension_cases, tmp_path, configurations):
    extension, context = extension_cases
    with pytest.raises(ValueError):
        extension.run_baseline_performance_extension(
            configurations, output_dir=tmp_path)
    assert not context.systems


@pytest.mark.parametrize('success,status,states', [
    (False, -1, [[1.]]), (True, 1, [[1.]]), (True, 0, [[np.nan]]),
    ])
def test_baseline_extension_rejects_incomplete_or_nonfinite_simulation(
        success, status, states):
    from exposan.werf.analyses import baseline_performance_extension as extension

    system = SimpleNamespace(
        ID='D1', simulate=lambda **options: None,
        scope=SimpleNamespace(sol=SimpleNamespace(
            success=success, status=status, message='incomplete', y=np.array(states))),
        )
    with pytest.raises(RuntimeError):
        extension._simulate(system, extension.SIMULATE_KWARGS)


@pytest.mark.parametrize('layers,total_height', [
    (0, 0.), (1, 4e-4), (3, 4e-4), (3, 0.), (3, 2e-3),
    ])
def test_tf_retained_mass_counts_liquid_and_biofilm(layers, total_height):
    import qsdsan as qs
    from qsdsan.utils import get_SRT
    from exposan.werf._units import TricklingFilter

    cmps = pc.create_masm2d_cmps()
    tf = TricklingFilter(
        'TF', ins=qs.WasteStream('feed'), outs='effluent', N_section=2,
        N_biofilm_layer=layers, V_sections=[2., 3.], filter_bed_surface=10.,
        filter_bed_depth=2., media_specific_surface=4.,
        initial_biofilm_thickness=0. if layers else None, suspended_growth_model=None,
        )
    indices = cmps.indices(('X_H', 'X_AUT'))
    factors = cmps.i_mass[indices]
    liquid = np.zeros((2, len(cmps)+1))
    liquid[:, indices] = [[10., 20.], [30., 40.]]
    liquid[:, -1] = 1e6  # Flow entries are not biomass concentrations.
    expected = 2.*(10.*factors[0]+20.*factors[1])
    expected += 3.*(30.*factors[0]+40.*factors[1])
    tf._state = liquid.ravel()
    if layers:
        biofilm = np.zeros((2, layers, len(cmps)))
        biofilm[:, :, indices] = [100., 200.]
        biofilm[:, :, cmps.index('X_S')] = 1e6
        heights = np.array([total_height, total_height/2.])
        tf._state = np.concatenate((liquid.ravel(), biofilm.ravel(), heights))
        for height in heights:
            effective = max(min(height, 1e-3)/layers, 1e-3/layers*1e-6)
            expected += layers*40.*effective*(100.*factors[0]+200.*factors[1])
    saved_state = tf._state.copy()
    assert tf.get_retained_mass(('X_H', 'X_AUT')) == pytest.approx(expected)
    np.testing.assert_array_equal(tf._state, saved_state)

    waste = qs.WasteStream('waste')
    waste.set_flow_by_concentration(
        100., {'X_H': 10., 'X_AUT': 20.}, units=('m3/d', 'mg/L'))
    system = qs.System('inventory_test', path=(tf,))
    assert get_SRT(system, ('X_H', 'X_AUT'), wastage=[waste]) == pytest.approx(
        expected/(100.*(10.*factors[0]+20.*factors[1])))


# %% PAO, PP and PHA initialization
@pytest.mark.parametrize('strength', ['low', 'mid', 'high'])
def test_seed_ratios_density_and_unchanged_liquid(strength):
    options = dict(wastewater_strength=strength, gpsx_aob_proxy=True,
                   tf_kwargs=dict(N_section=6, N_biofilm_layer=4,
                                  initial_biofilm_thickness=4e-4))
    legacy = create_trickling_filter_validation_system(**options, tf_init_kwargs={'pao_fraction': None}).flowsheet.unit.TF
    tf = create_trickling_filter_validation_system(**options).flowsheet.unit.TF
    cmps = tf.components
    h, pao, pp, pha = cmps.indices(('X_H', 'X_PAO', 'X_PP', 'X_PHA'))
    bio = tf._concs[:, 1:]
    np.testing.assert_allclose(bio[..., pao]/(bio[..., pao]+bio[..., h]), .01)
    np.testing.assert_allclose(bio[..., pp]/bio[..., pao], .17)
    np.testing.assert_allclose(bio[..., pha]/bio[..., pao], .05)
    assert bio[..., pao].min() > 1.
    particulate = np.array([c.particle_size == 'Particulate' for c in cmps])
    np.testing.assert_allclose(bio[..., particulate] @ cmps.i_mass[particulate],
                               tf.biofilm_density * tf.biofilm_dry_content)
    np.testing.assert_array_equal(tf._concs[:, 0], legacy._concs[:, 0])
    np.testing.assert_array_equal(bio[..., ~particulate], legacy._concs[:, 1:, ~particulate])
    np.testing.assert_allclose(tf.total_biofilm_thickness, 4e-4)


def test_custom_seed_and_capacity():
    tf = create_trickling_filter_validation_system(asm_kwargs={'K_MAX': .2}, tf_init_kwargs=dict(
        pao_fraction=.02, pp_saturation=.4, pha_ratio=.1)).flowsheet.unit.TF
    bio = tf._concs[:, 1:]
    h, pao, pp, pha = tf.components.indices(('X_H', 'X_PAO', 'X_PP', 'X_PHA'))
    np.testing.assert_allclose(bio[..., pao]/(bio[..., pao]+bio[..., h]), .02)
    np.testing.assert_allclose(bio[..., pp]/bio[..., pao], .08)
    np.testing.assert_allclose(bio[..., pha]/bio[..., pao], .1)


@pytest.mark.parametrize('options', [dict(pao_fraction=-.1), dict(pao_fraction=1),
                                    dict(pp_saturation=1.1), dict(pha_ratio=np.nan)])
def test_invalid_seed(options):
    with pytest.raises(ValueError):
        create_trickling_filter_validation_system(tf_init_kwargs=options)


def test_short_seeded_run():
    system = create_trickling_filter_validation_system(wastewater_strength='high', gpsx_aob_proxy=True,
                    tf_kwargs={'N_section': 6, 'N_biofilm_layer': 4})
    tf = system.flowsheet.unit.TF
    pao = tf._concs[:, 1:, tf.components.index('X_PAO')]
    print('Initial biofilm PAO range [mg COD/L]:', pao.min(), pao.max())
    system.simulate(method='BDF', t_span=(0, .05), rtol=1e-3, atol=1e-2,
                    stop_at_steady_state=False)
    assert system.scope.sol.success
    assert np.isfinite(system.scope.sol.y).all()


# %% Core TF transport, reactions and dynamic-thickness conservation
def test_trickling_filter_matches_pfr_ode():
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc, unit_operations as su
    from exposan.werf import default_as_init
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    cmps = qs.get_components()
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)
    N = 4
    Q = 120.0
    Vs = np.array([18., 19., 20., 21.])
    f_ins = [[1., 0., 0., 0.]]

    inf_tf = WasteStream('inf_tf', T=293.15)
    inf_pfr = WasteStream('inf_pfr', T=293.15)
    for inf in (inf_tf, inf_pfr):
        inf.set_flow_by_concentration(
            Q,
            concentrations=default_as_init,
            units=('m3/d', 'mg/L'),
            )

    init = np.tile(cmps.kwarray(default_as_init), (N, 1))
    tf = TricklingFilter(
        'TF', ins=inf_tf, outs='tf_eff',
        N_section=N,
        N_biofilm_layer=0,
        V_sections=Vs,
        influent_fractions=f_ins,
        suspended_growth_model=asm,
        DO_ID='S_O2',
        )
    pfr = su.PFR(
        'PFR', ins=inf_pfr, outs='pfr_eff',
        N_tanks_in_series=N,
        V_tanks=Vs,
        influent_fractions=f_ins,
        internal_recycles=[],
        DO_ID='S_O2',
        kLa=[0.]*N,
        suspended_growth_model=asm,
        gas_stripping=False,
        )
    tf.set_init_conc(concentrations=init)
    pfr.set_init_conc(concentrations=init)

    for unit in (tf, pfr):
        unit._run()
        unit._init_state()
        QC_ins = np.vstack([stream.state for stream in unit.ins])
        dQC_ins = np.zeros_like(QC_ins)
        unit.ODE(0., QC_ins, unit._state.copy(), dQC_ins)

    np.testing.assert_allclose(tf._dstate, pfr._dstate, rtol=1e-12, atol=1e-12)


def test_trickling_filter_biofilm_ode_smoke():
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf import default_as_init
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_biofilm_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    cmps = qs.get_components()
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)

    inf = WasteStream('inf_tf_biofilm', T=293.15)
    inf.set_flow_by_concentration(
        120.0,
        concentrations=default_as_init,
        units=('m3/d', 'mg/L'),
        )
    tf = TricklingFilter(
        'TF_biofilm', ins=inf, outs='tf_biofilm_eff',
        N_section=2,
        N_biofilm_layer=3,
        initial_biofilm_thickness=0.,
        # diffusivity=1e-4,
        suspended_growth_model=asm,
        DO_ID='S_O2',
        )
    # tf.set_init_conc(np.tile(cmps.kwarray(default_as_init), (2, 1)))
    init = np.zeros((2, 4, len(cmps)))
    init[:, 0, :] = cmps.kwarray(default_as_init)
    tf.set_init_conc(init)
    tf._run()
    tf._init_state()
    QC_ins = np.vstack([stream.state for stream in tf.ins])
    dQC_ins = np.zeros_like(QC_ins)
    tf.ODE(0., QC_ins, tf._state.copy(), dQC_ins)
    assert tf._dstate.shape == tf._state.shape
    assert np.isfinite(tf._dstate).all()
    liq_stop = tf.N_section * tf._ncol
    bio_stop = liq_stop + tf.N_section * tf.N_biofilm_layer * len(cmps)
    assert np.allclose(tf.biofilm_thickness, 0.)
    # assert (tf._dstate[bio_stop:] > 0.).all()
    # dh = tf._dstate[bio_stop:].reshape((tf.N_section, tf.N_biofilm_layer))
    # assert (dh[:, 0] > 0.).all()
    # assert np.allclose(dh[:, 1:], 0.)
    dh = tf._dstate[bio_stop:]
    assert (dh > 0.).all()
    # H = tf._state[bio_stop:].reshape((tf.N_section, tf.N_biofilm_layer))
    # H[:, 0] = tf.max_biofilm_thickness / tf.N_biofilm_layer
    tf._state[bio_stop:] = tf.max_biofilm_thickness
    tf.ODE(0., QC_ins, tf._state.copy(), dQC_ins)
    # dh = tf._dstate[bio_stop:].reshape((tf.N_section, tf.N_biofilm_layer))
    # assert np.allclose(dh[:, 0], 0.)
    # assert (dh[:, 1] > 0.).all()
    # assert np.allclose(dh[:, 2:], 0.)
    dh = tf._dstate[bio_stop:]
    assert (dh <= 0.).all()


def test_trickling_filter_transport_only_mass_conservation():
    import warnings
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_transport_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    cmps = qs.get_components()
    idx = cmps.index('S_A')
    N = 1
    L = 3
    V_liq = 500.

    inf = WasteStream('inf_tf_transport', T=293.15)
    inf.set_flow_by_concentration(
        # 0.,
        # concentrations={'S_A': 0.},
        1.,
        concentrations={'S_A': 100.},
        units=('m3/d', 'mg/L'),
        )
    tf = TricklingFilter(
        'TF_transport', ins=inf, outs='tf_transport_eff',
        N_section=N,
        N_biofilm_layer=L,
        V_sections=[V_liq],
        suspended_growth_model=None,
        DO_ID='S_O2',
        kLa=[0.],
        attachment_rate=0.,
        detachment_rate=0.,
        internal_solids_exchange_rate=0.,
        )
    init = np.zeros((N, L+1, len(cmps)))
    init[0, 0, idx] = 100.
    init[0, 1, idx] = 10.
    init[0, 2, idx] = 30.
    init[0, 3, idx] = 60.
    tf.set_init_conc(init)
    tf._run()
    tf._init_state()
    QC_ins = np.vstack([stream.state for stream in tf.ins])
    dQC_ins = np.zeros_like(QC_ins)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tf.ODE(0., QC_ins, tf._state.copy(), dQC_ins)

    liq_stop = tf.N_section * tf._ncol
    bio_stop = liq_stop + tf.N_section * tf.N_biofilm_layer * len(cmps)
    dC_liq = tf._dstate[:liq_stop].reshape((N, tf._ncol))[:, :-1]
    dC_bio = tf._dstate[liq_stop:bio_stop].reshape((N, L, len(cmps)))
    # dH = tf._dstate[bio_stop:].reshape((N, L))
    dH = tf._dstate[bio_stop:]
    C_bio = tf._state[liq_stop:bio_stop].reshape((N, L, len(cmps)))
    H = tf.biofilm_thickness
    A = tf.filter_bed_surface * tf.filter_bed_depth \
        * tf.media_specific_surface / N
    # dM = V_liq * dC_liq[0, idx] \
    #     + np.sum(A * (H[0] * dC_bio[0, :, idx] + C_bio[0, :, idx] * dH[0]))
    dM = (
        V_liq * dC_liq[0, idx]
        + np.sum(A * (H[0] * dC_bio[0, :, idx] + C_bio[0, :, idx] * dH[0]/L)))
    assert np.isclose(dM, 0., atol=1e-8)
    assert np.allclose(dH, 0.)


def test_trickling_filter_initial_biofilm_thickness_fills_layers():
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_thickness_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    inf = WasteStream('inf_tf_thickness', T=293.15)
    inf.set_flow_by_concentration(
        1.,
        concentrations={'X_S': 10.},
        units=('m3/d', 'mg/L'),
        )
    tf = TricklingFilter(
        'TF_thickness', ins=inf, outs='tf_thickness_eff',
        N_section=2,
        N_biofilm_layer=3,
        initial_biofilm_thickness=[4e-4, 5e-4],
        suspended_growth_model=None,
        kLa=[0., 0.],
        )
    layer_max = tf.max_biofilm_thickness / tf.N_biofilm_layer
    expected = np.array([
        [layer_max, 4e-4-layer_max, 0.],
        [layer_max, 5e-4-layer_max, 0.],
        ])
    np.testing.assert_allclose(tf.initial_biofilm_thickness, expected)
    explicit = np.full((2, 3), layer_max/2)
    tf.initial_biofilm_thickness = explicit
    np.testing.assert_allclose(tf.initial_biofilm_thickness, explicit)
    with np.testing.assert_raises(ValueError):
        tf.initial_biofilm_thickness = np.full((2, 3), layer_max*2)


def test_trickling_filter_HLR_from_filter_bed_surface():
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_depth_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    Q = 1000.
    surface = 200.
    media_surface = 100.
    depth = 2.
    # HLR = 13.
    # ft_to_m = 0.3048
    # HLR = 70.

    inf = WasteStream('inf_tf_depth', T=293.15)
    inf.set_flow_by_concentration(
        Q,
        concentrations={'S_A': 10.},
        units=('m3/d', 'mg/L'),
        )
    tf = TricklingFilter(
        'TF_depth', ins=inf, outs='tf_depth_eff',
        N_section=1,
        N_biofilm_layer=0,
        filter_bed_depth=depth,
        filter_bed_surface=surface,
        media_specific_surface=media_surface,
        # hydraulic_loading_rate=HLR,
        suspended_growth_model=None,
        )
    tf._run()
    tf._init_state()
    # expected_depth = Q / (HLR*ft_to_m*surface*media_surface)
    # expected_contact_area = Q / (HLR*ft_to_m)
    # expected_depth = Q / (HLR*surface*media_surface)
    # expected_contact_area = Q / HLR
    expected_HLR = Q / surface
    expected_contact_area = surface * depth * media_surface
    assert np.isclose(tf.filter_bed_depth, depth)
    assert np.isclose(tf.hydraulic_loading_rate, expected_HLR)
    assert np.isclose(
        tf.filter_bed_surface * tf.filter_bed_depth
        * tf.media_specific_surface,
        expected_contact_area,
        )

    tf_manual = TricklingFilter(
        'TF_depth_manual', ins=inf, outs='tf_depth_manual_eff',
        N_section=1,
        N_biofilm_layer=0,
        filter_bed_depth=2.0,
        filter_bed_surface=surface,
        media_specific_surface=media_surface,
        # hydraulic_loading_rate=HLR,
        suspended_growth_model=None,
        )
    tf_manual._run()
    tf_manual._init_state()
    assert np.isclose(tf_manual.filter_bed_depth, 2.0)
    # expected_manual_HLR = Q / (surface * tf_manual.filter_bed_depth * media_surface)
    expected_manual_HLR = Q / surface
    assert np.isclose(tf_manual.hydraulic_loading_rate, expected_manual_HLR)
    tf_manual_recycle = TricklingFilter(
        'TF_depth_manual_recycle', ins=inf, outs='tf_depth_manual_recycle_eff',
        N_section=1,
        N_biofilm_layer=0,
        filter_bed_depth=2.0,
        filter_bed_surface=surface,
        media_specific_surface=media_surface,
        # hydraulic_loading_rate=HLR,
        recirculation_ratio=1.0,
        suspended_growth_model=None,
        )
    tf_manual_recycle._run()
    tf_manual_recycle._init_state()
    # expected_recycle_HLR = 2 * Q / (
    #     surface * tf_manual_recycle.filter_bed_depth * media_surface)
    expected_recycle_HLR = 2 * Q / surface
    assert np.isclose(tf_manual_recycle.hydraulic_loading_rate, expected_recycle_HLR)


def test_trickling_filter_recirculation_mass_balance():
    import warnings
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_recycle_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    cmps = qs.get_components()
    idx = cmps.index('S_A')
    Q = 100.
    R = 0.5
    V = np.array([10., 20.])
    C_in = 10.
    C0 = 0.
    C1 = 100.

    inf = WasteStream('inf_tf_recycle', T=293.15)
    inf.set_flow_by_concentration(
        Q,
        concentrations={'S_A': C_in},
        units=('m3/d', 'mg/L'),
        )
    # inf.state = np.append(cmps.kwarray({'S_A': C_in}), Q)
    tf = TricklingFilter(
        'TF_recycle', ins=inf, outs='tf_recycle_eff',
        N_section=2,
        N_biofilm_layer=0,
        V_sections=V,
        recirculation_ratio=R,
        suspended_growth_model=None,
        kLa=[0., 0.],
        )
    inf.state = np.append(cmps.kwarray({'S_A': C_in}), Q)
    init = np.zeros((2, len(cmps)))
    init[0, idx] = C0
    init[1, idx] = C1
    tf.set_init_conc(init)
    tf._run()
    tf._init_state()
    QC_ins = np.vstack([stream.state for stream in tf.ins])
    dQC_ins = np.zeros_like(QC_ins)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tf.ODE(0., QC_ins, tf._state.copy(), dQC_ins)

    Q_recycle = R * Q
    Q_applied = Q + Q_recycle
    dC = tf._dstate.reshape((2, tf._ncol))[:, :-1]
    expected = np.array([
        (Q*C_in + Q_recycle*C1 - Q_applied*C0) / V[0],
        (Q_applied*C0 - Q_applied*C1) / V[1],
        ])
    np.testing.assert_allclose(dC[:, idx], expected)
    np.testing.assert_allclose(tf._Qs, np.full(2, Q_applied))
    dM = np.sum(V * dC[:, idx])
    assert np.isclose(dM, Q*C_in - Q*C1)


def test_trickling_filter_particulate_mass_conservation():
    import warnings
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_particulate_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    cmps = qs.get_components()
    idx = cmps.index('X_S')
    N = 1
    L = 3
    V_liq = 500.

    inf = WasteStream('inf_tf_particulate', T=293.15)
    inf.set_flow_by_concentration(
        1.,
        concentrations={'X_S': 100.},
        units=('m3/d', 'mg/L'),
        )
    tf = TricklingFilter(
        'TF_particulate', ins=inf, outs='tf_particulate_eff',
        N_section=N,
        N_biofilm_layer=L,
        V_sections=[V_liq],
        # initial_biofilm_thickness=1e-4,
        initial_biofilm_thickness=np.full((N, L), 1e-4),
        suspended_growth_model=None,
        DO_ID='S_O2',
        kLa=[0.],
        attachment_rate=0.5,
        detachment_rate=0.04,
        internal_solids_exchange_rate=0.,
        dynamic_biofilm_thickness=False,
        )
    init = np.zeros((N, L+1, len(cmps)))
    init[0, 0, idx] = 100.
    init[0, 1, idx] = 20.
    tf.set_init_conc(init)
    tf._run()
    tf._init_state()
    QC_ins = np.vstack([stream.state for stream in tf.ins])
    dQC_ins = np.zeros_like(QC_ins)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tf.ODE(0., QC_ins, tf._state.copy(), dQC_ins)

    liq_stop = tf.N_section * tf._ncol
    bio_stop = liq_stop + tf.N_section * tf.N_biofilm_layer * len(cmps)
    dC_liq = tf._dstate[:liq_stop].reshape((N, tf._ncol))[:, :-1]
    dC_bio = tf._dstate[liq_stop:bio_stop].reshape((N, L, len(cmps)))
    # dH = tf._dstate[bio_stop:].reshape((N, L))
    dH = tf._dstate[bio_stop:]
    C_bio = tf._state[liq_stop:bio_stop].reshape((N, L, len(cmps)))
    H = tf.biofilm_thickness
    A = tf.filter_bed_surface * tf.filter_bed_depth \
        * tf.media_specific_surface / N
    # dM = V_liq * dC_liq[0, idx] \
    #     + np.sum(A * (H[0] * dC_bio[0, :, idx] + C_bio[0, :, idx] * dH[0]))
    dM = (
        V_liq * dC_liq[0, idx]
        + np.sum(A * (H[0] * dC_bio[0, :, idx] + C_bio[0, :, idx] * dH[0]/L)))
    assert np.isclose(dM, 0., atol=1e-8)
    # assert dH[0, 0] > 0.
    # # assert np.allclose(dH[0, 1:], 0.)
    # assert dH[0, 0] > dH[0, 1] > dH[0, 2] >= 0.
    assert np.allclose(dH, 0.)


def test_trickling_filter_masm2d_reaction_in_biofilm_only():
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf import default_as_init
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_reaction_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    cmps = qs.get_components()
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)
    N = 1
    L = 2
    Q = 1.
    init_conc = cmps.kwarray(default_as_init)

    inf = WasteStream('inf_tf_reaction', T=293.15)
    inf.set_flow_by_concentration(
        Q,
        concentrations=default_as_init,
        units=('m3/d', 'mg/L'),
        )
    tf = TricklingFilter(
        'TF_reaction', ins=inf, outs='tf_reaction_eff',
        N_section=N,
        N_biofilm_layer=L,
        V_sections=[500.],
        initial_biofilm_thickness=1e-4,
        suspended_growth_model=asm,
        DO_ID='S_O2',
        kLa=[0.],
        air_O2_transfer_coefficient=0.,
        diffusivity=np.zeros(len(cmps)),
        attachment_rate=0.,
        detachment_rate=0.,
        internal_solids_exchange_rate=0.,
        )
    init = np.zeros((N, L+1, len(cmps)))
    init[0, :, :] = init_conc
    tf.set_init_conc(init)
    tf._run()
    tf._init_state()
    QC_ins = np.vstack([stream.state for stream in tf.ins])
    dQC_ins = np.zeros_like(QC_ins)
    tf.ODE(0., QC_ins, tf._state.copy(), dQC_ins)

    liq_stop = tf.N_section * tf._ncol
    bio_stop = liq_stop + tf.N_section * tf.N_biofilm_layer * len(cmps)
    dC_liq = tf._dstate[:liq_stop].reshape((N, tf._ncol))[:, :-1]
    dC_bio = tf._dstate[liq_stop:bio_stop].reshape((N, L, len(cmps)))
    assert np.allclose(dC_liq, 0., atol=1e-12)
    assert np.any(np.abs(dC_bio) > 1e-12)
    assert np.isfinite(tf._dstate).all()


def test_trickling_filter_DO_aeration_and_biofilm_diffusion():
    import warnings
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_do_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    cmps = qs.get_components()
    DO_idx = cmps.index('S_O2')
    N = 1
    L = 2
    V_liq = 500.
    Q = 1.
    DO_liq = 1.
    DO_sat = 8.
    # kLa = 2.
    k_air = 2.

    inf = WasteStream('inf_tf_do', T=293.15)
    inf.set_flow_by_concentration(
        Q,
        concentrations={'S_O2': DO_liq},
        units=('m3/d', 'mg/L'),
        )
    D = np.zeros(len(cmps))
    D[DO_idx] = 2.5e-5 * 8.64
    tf = TricklingFilter(
        'TF_do', ins=inf, outs='tf_do_eff',
        N_section=N,
        N_biofilm_layer=L,
        V_sections=[V_liq],
        # initial_biofilm_thickness=1e-4,
        initial_biofilm_thickness=np.full((N, L), 1e-4),
        suspended_growth_model=None,
        DO_ID='S_O2',
        # kLa=[kLa],
        air_O2_transfer_coefficient=k_air,
        DO_sat=DO_sat,
        diffusivity=D,
        attachment_rate=0.,
        detachment_rate=0.,
        internal_solids_exchange_rate=0.,
        )
    init = np.zeros((N, L+1, len(cmps)))
    init[0, 0, DO_idx] = DO_liq
    tf.set_init_conc(init)
    tf._run()
    tf._init_state()
    QC_ins = np.vstack([stream.state for stream in tf.ins])
    dQC_ins = np.zeros_like(QC_ins)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tf.ODE(0., QC_ins, tf._state.copy(), dQC_ins)

    liq_stop = tf.N_section * tf._ncol
    bio_stop = liq_stop + tf.N_section * tf.N_biofilm_layer * len(cmps)
    dC_liq = tf._dstate[:liq_stop].reshape((N, tf._ncol))[:, :-1]
    dC_bio = tf._dstate[liq_stop:bio_stop].reshape((N, L, len(cmps)))
    A = tf.filter_bed_surface * tf.filter_bed_depth \
        * tf.media_specific_surface / N
    k_liq = D[DO_idx] / tf.max_liquid_film_thickness
    # expected_liq = kLa * (DO_sat - DO_liq) - A / V_liq * k_liq * DO_liq
    expected_liq = k_air * A / V_liq * (DO_sat - DO_liq) - A / V_liq * k_liq * DO_liq
    assert np.isclose(dC_liq[0, DO_idx], expected_liq)
    assert dC_bio[0, 0, DO_idx] > 0.
    assert np.isclose(dC_bio[0, 1, DO_idx], 0.)


def test_trickling_filter_dynamic_simulation_smoke():
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf import default_as_init
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_dynamic_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    cmps = qs.get_components()
    asm = pc.mASM2d(electron_acceptor_dependent_decay=True)

    inf = WasteStream('inf_tf_dynamic', T=293.15)
    inf.set_flow_by_concentration(
        1.,
        concentrations=default_as_init,
        units=('m3/d', 'mg/L'),
        )
    tf = TricklingFilter(
        'TF_dynamic', ins=inf, outs='tf_dynamic_eff',
        N_section=1,
        N_biofilm_layer=1,
        V_sections=[10.],
        initial_biofilm_thickness=1e-4,
        suspended_growth_model=asm,
        DO_ID='S_O2',
        kLa=[0.],
        diffusivity=np.zeros(len(cmps)),
        attachment_rate=0.,
        detachment_rate=0.,
        internal_solids_exchange_rate=0.,
        )
    init = np.zeros((1, 2, len(cmps)))
    init[0, :, :] = cmps.kwarray(default_as_init)
    tf.set_init_conc(init)
    tf.simulate(t_span=(0., 0.01), method='BDF')

    eff, = tf.outs
    assert eff.state.shape == (len(cmps)+1,)
    assert np.isfinite(eff.state).all()
    assert tf._state.shape == tf._dstate.shape
    assert tf.biofilm_thickness.shape == (1, 1)
    assert np.all(tf.biofilm_thickness >= 0.)


def test_trickling_filter_biofilm_thickness_bounds():
    import warnings
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_thickness_bounds_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    cmps = qs.get_components()
    N = 1
    L = 2

    inf = WasteStream('inf_tf_thickness_bounds', T=293.15)
    inf.set_flow_by_concentration(
        1.,
        concentrations={'S_A': 0.},
        units=('m3/d', 'mg/L'),
        )
    tf = TricklingFilter(
        'TF_thickness_bounds', ins=inf, outs='tf_thickness_bounds_eff',
        N_section=N,
        N_biofilm_layer=L,
        V_sections=[10.],
        initial_biofilm_thickness=np.full((N, L), 1e-4),
        suspended_growth_model=None,
        DO_ID='S_O2',
        kLa=[0.],
        diffusivity=np.zeros(len(cmps)),
        attachment_rate=0.,
        detachment_rate=0.,
        internal_solids_exchange_rate=0.,
        )
    init = np.zeros((N, L+1, len(cmps)))
    tf.set_init_conc(init)
    tf._run()
    tf._init_state()
    QC_ins = np.vstack([stream.state for stream in tf.ins])
    dQC_ins = np.zeros_like(QC_ins)
    QC = tf._state.copy()
    liq_stop = tf.N_section * tf._ncol
    bio_stop = liq_stop + tf.N_section * tf.N_biofilm_layer * len(cmps)
    # layer_max = tf.max_biofilm_thickness / L
    # QC[bio_stop] = 2 * layer_max
    # QC[bio_stop+1] = -layer_max
    QC[bio_stop] = 2 * tf.max_biofilm_thickness
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # tf.ODE(0., QC, dQC_ins, dQC_ins)
        tf.ODE(0., QC_ins, QC, dQC_ins)
    # dH = tf._dstate[bio_stop:].reshape((N, L))
    # assert dH[0, 0] < 0.
    # assert dH[0, 1] > 0.
    dH = tf._dstate[bio_stop:]
    assert dH[0] < 0.
    QC[bio_stop] = -tf.max_biofilm_thickness
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tf.ODE(0., QC_ins, QC, dQC_ins)
    assert tf._dstate[bio_stop] > 0.


def test_trickling_filter_frozen_biofilm_thickness():
    import warnings
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_frozen_thickness_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    cmps = qs.get_components()
    idx = cmps.index('X_S')
    N = 1
    L = 2

    inf = WasteStream('inf_tf_frozen_thickness', T=293.15)
    inf.set_flow_by_concentration(
        1.,
        concentrations={'X_S': 100.},
        units=('m3/d', 'mg/L'),
        )
    tf = TricklingFilter(
        'TF_frozen_thickness', ins=inf, outs='tf_frozen_thickness_eff',
        N_section=N,
        N_biofilm_layer=L,
        V_sections=[10.],
        initial_biofilm_thickness=np.full((N, L), 1e-4),
        suspended_growth_model=None,
        DO_ID='S_O2',
        kLa=[0.],
        diffusivity=np.zeros(len(cmps)),
        attachment_rate=0.5,
        detachment_rate=0.04,
        internal_solids_exchange_rate=0.,
        dynamic_biofilm_thickness=False,
        )
    init = np.zeros((N, L+1, len(cmps)))
    init[0, 0, idx] = 100.
    init[0, 1, idx] = 20.
    tf.set_init_conc(init)
    tf._run()
    tf._init_state()
    QC_ins = np.vstack([stream.state for stream in tf.ins])
    dQC_ins = np.zeros_like(QC_ins)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tf.ODE(0., QC_ins, tf._state.copy(), dQC_ins)
    liq_stop = tf.N_section * tf._ncol
    bio_stop = liq_stop + tf.N_section * tf.N_biofilm_layer * len(cmps)
    # dH = tf._dstate[bio_stop:].reshape((N, L))
    dH = tf._dstate[bio_stop:]
    assert np.allclose(dH, 0.)


def test_trickling_filter_internal_solids_exchange_conservation():
    import warnings
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf._units import TricklingFilter

    flowsheet = qs.Flowsheet('tf_internal_exchange_test')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    pc.create_masm2d_cmps()
    cmps = qs.get_components()
    idx = cmps.index('X_S')
    N = 1
    L = 3
    V_liq = 500.

    inf = WasteStream('inf_tf_internal_exchange', T=293.15)
    inf.set_flow_by_concentration(
        1.,
        concentrations={'X_S': 0.},
        units=('m3/d', 'mg/L'),
        )
    tf = TricklingFilter(
        'TF_internal_exchange', ins=inf, outs='tf_internal_exchange_eff',
        N_section=N,
        N_biofilm_layer=L,
        # Isolate the specified one-way exchange from moving-layer packing
        # transport, which is tested separately with dynamic thickness.
        dynamic_biofilm_thickness=False,
        V_sections=[V_liq],
        # initial_biofilm_thickness=1e-4,
        initial_biofilm_thickness=np.full((N, L), 1e-4),
        suspended_growth_model=None,
        DO_ID='S_O2',
        kLa=[0.],
        attachment_rate=0.,
        detachment_rate=0.,
        internal_solids_exchange_rate=2e-5,
        )
    init = np.zeros((N, L+1, len(cmps)))
    init[0, 1, idx] = 0.
    init[0, 2, idx] = 100.
    init[0, 3, idx] = 50.
    tf.set_init_conc(init)
    tf._run()
    tf._init_state()
    QC_ins = np.vstack([stream.state for stream in tf.ins])
    dQC_ins = np.zeros_like(QC_ins)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tf.ODE(0., QC_ins, tf._state.copy(), dQC_ins)

    liq_stop = tf.N_section * tf._ncol
    bio_stop = liq_stop + tf.N_section * tf.N_biofilm_layer * len(cmps)
    dC_bio = tf._dstate[liq_stop:bio_stop].reshape((N, L, len(cmps)))
    # dH = tf._dstate[bio_stop:].reshape((N, L))
    dH = tf._dstate[bio_stop:]
    H = tf.biofilm_thickness
    A = tf.filter_bed_surface * tf.filter_bed_depth \
        * tf.media_specific_surface / N
    dM = np.sum(A * H[0] * dC_bio[0, :, idx])
    assert np.isclose(dM, 0., atol=1e-8)
    assert np.allclose(dH, 0.)
    assert dC_bio[0, 0, idx] > 0.
    assert dC_bio[0, 2, idx] < 0.


def _make_thickness_balance_filter(N, L, case):
    import numpy as np
    import qsdsan as qs
    from qsdsan import WasteStream, process_models as pc
    from exposan.werf import default_as_init
    from exposan.werf._units import TricklingFilter

    qs.main_flowsheet.set_flowsheet(qs.Flowsheet('tf_inventory_test'))
    cmps = pc.create_masm2d_cmps()
    model = pc.mASM2d(electron_acceptor_dependent_decay=True) \
        if case == 'reaction' else None
    inf = WasteStream('inf_tf_inventory', T=293.15)
    inf.set_flow_by_concentration(
        100., concentrations={'S_A': 50., 'S_I': 30., 'X_I': 100.},
        units=('m3/d', 'mg/L'),
        )
    shrinking = case in ('shrink', 'max_shrink')
    # Direct ODE tests must initialize stream states as System.simulate does.
    # inf._init_state()
    tf = TricklingFilter(
        'TF_inventory', ins=inf, outs='tf_inventory_eff',
        N_section=N, N_biofilm_layer=L,
        V_sections=np.linspace(100., 200., N),
        initial_biofilm_thickness=3e-4,
        suspended_growth_model=model,
        recirculation_ratio=0.7,
        air_O2_transfer_coefficient=0.,
        attachment_rate=0. if shrinking or model is not None else 0.5,
        detachment_rate=0.04 if shrinking else 0.,
        internal_solids_exchange_rate=2e-5,
        dynamic_biofilm_thickness=case != 'frozen',
        )
    # Connect the inlet first: its state is stored in the receiving unit.
    inf._init_state()
    if case == 'volume_floor':
        tf.filter_bed_surface = 1e-15
    if case in ('floor', 'floor_growth', 'below_zero'):
        # Isolate the volume floor from very large diffusion gradients across
        # sub-nanometre layers; other cases exercise the full diffusion balance.
        # Soluble tracers still test dilution/concentration as volume changes.
        tf.diffusivity = np.zeros(len(cmps))
    init = np.zeros((N, L+1, len(cmps)))
    init[:, 0, :] = cmps.kwarray({'S_A': 50., 'S_I': 30., 'X_I': 100.})
    init[:, 0, :] *= np.linspace(1., 0.5, N)[:, None]
    init[:, 1:, :] = cmps.kwarray(default_as_init)
    init[:, 1:, :] *= np.linspace(1., 0.5, L)[None, :, None]
    tf.set_init_conc(init)
    tf._run()
    tf._init_state()
    h_min = max(tf.max_biofilm_thickness/L * 1e-6, 1e-12)
    if case in ('max_growth', 'max_shrink'):
        tf._state[-N:] = tf.max_biofilm_thickness
    elif case == 'floor':
        tf._state[-N:] = 0.5 * L * h_min
    elif case == 'floor_growth':
        # tf._state[-N:] = L * h_min
        # Multiplication/division by L can round below the floor; test the
        # growing side explicitly without changing the production threshold.
        tf._state[-N:] = np.nextafter(L * h_min, np.inf)
    elif case == 'above_max':
        tf._state[-N:] = 1.1 * tf.max_biofilm_thickness
    elif case == 'below_zero':
        tf._state[-N:] = -tf.max_biofilm_thickness
    return tf


@pytest.mark.parametrize('N,L', [(1, 1), (2, 3), (6, 5)])
@pytest.mark.parametrize('case', [
    'growth', 'shrink', 'reaction', 'frozen', 'max_growth', 'max_shrink',
    'floor', 'floor_growth', 'above_max', 'below_zero', 'volume_floor',
    ])
def test_trickling_filter_dynamic_thickness_component_balance(N, L, case):
    """Changing biofilm volume must not introduce a component source or sink."""
    import warnings
    import numpy as np

    tf = _make_thickness_balance_filter(N, L, case)
    cmps = tf.components
    n = len(cmps)
    liq_stop = N * tf._ncol
    bio_stop = liq_stop + N * L * n
    QC = tf._state.copy()
    QC_ins = np.vstack([stream.state for stream in tf.ins])
    assert (QC_ins[:, -1] > 0.).all()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='.*non-reactive unit')
        tf.ODE(0., QC_ins, QC, np.zeros_like(QC_ins))

    liquid = QC[:liq_stop].reshape(N, tf._ncol)[:, :-1]
    biofilm = QC[liq_stop:bio_stop].reshape(N, L, n)
    d_liquid = tf._dstate[:liq_stop].reshape(N, tf._ncol)[:, :-1]
    d_biofilm = tf._dstate[liq_stop:bio_stop].reshape(N, L, n)
    dH = tf._dstate[bio_stop:]
    A = max(tf.filter_bed_surface * tf.filter_bed_depth
            * tf.media_specific_surface / N, 1e-12)
    h_min = max(tf.max_biofilm_thickness/L * 1e-6, 1e-12)
    H = np.clip(QC[bio_stop:], 0., tf.max_biofilm_thickness)
    V = np.repeat(np.maximum(A * np.maximum(H/L, h_min), 1e-12)[:, None],
                  L, axis=1)
    if case in ('floor', 'above_max', 'below_zero', 'volume_floor'):
        dV = np.zeros_like(V)
    else:
        dV = np.repeat((A * dH/L)[:, None], L, axis=1)
    if case in ('growth', 'floor', 'floor_growth', 'volume_floor'):
        assert (dH > 0.).all()
    elif case in ('shrink', 'max_shrink'):
        assert (dH < 0.).all()
    elif case in ('frozen', 'max_growth'):
        np.testing.assert_array_equal(dH, 0.)

    # Concentrations in mg/L are numerically g/m3, so inventories are in g.
    accumulation = (tf.V_sections[:, None] * d_liquid).sum(axis=0)
    accumulation += (V[:, :, None] * d_biofilm
                     + dV[:, :, None] * biofilm).sum(axis=(0, 1))
    external = (QC_ins[:, -1, None] * QC_ins[:, :-1]).sum(axis=0)
    external -= QC_ins[:, -1].sum() * liquid[-1]
    reaction = np.zeros(n)
    if case == 'reaction':
        model = tf.suspended_growth_model
        rates = np.apply_along_axis(model.rate_function, 2, biofilm)
        reaction = (V[:, :, None] * (rates @ model.stoichio_eval())).sum(axis=(0, 1))
        assert np.any(np.abs(dH) > 1e-12)
    np.testing.assert_allclose(accumulation, external + reaction,
                               rtol=1e-10, atol=1e-6)
    assert np.isfinite(tf._dstate).all()


@pytest.mark.parametrize('case', ['growth', 'shrink'])
def test_trickling_filter_dynamic_thickness_integrated_inventory(case):
    """BDF inventory change equals cumulative net input during growth/shrinkage."""
    import warnings
    import numpy as np
    from scipy.integrate import solve_ivp

    tf = _make_thickness_balance_filter(2, 3, case)
    N, L, n = tf.N_section, tf.N_biofilm_layer, len(tf.components)
    liq_stop = N * tf._ncol
    bio_stop = liq_stop + N * L * n
    QC_ins = np.vstack([stream.state for stream in tf.ins])
    dQC_ins = np.zeros_like(QC_ins)
    A = tf.filter_bed_surface * tf.filter_bed_depth * tf.media_specific_surface / N
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='.*non-reactive unit')
        ode = tf.ODE

    def rhs(t, y):
        state = y[:bio_stop+N].copy()
        ode(t, QC_ins, state, dQC_ins)
        liquid = state[:liq_stop].reshape(N, tf._ncol)[:, :-1]
        net_input = (QC_ins[:, -1, None] * QC_ins[:, :-1]).sum(axis=0)
        net_input -= QC_ins[:, -1].sum() * liquid[-1]
        return np.concatenate((tf._dstate, net_input))

    y0 = np.concatenate((tf._state, np.zeros(n)))
    sol = solve_ivp(rhs, (0., 0.05), y0, method='BDF', rtol=1e-8, atol=1e-11)
    assert sol.success, sol.message
    liquid = sol.y[:liq_stop].T.reshape(-1, N, tf._ncol)[:, :, :-1]
    biofilm = sol.y[liq_stop:bio_stop].T.reshape(-1, N, L, n)
    H = sol.y[bio_stop:bio_stop+N].T
    assert (H > 0.).all() and (H < tf.max_biofilm_thickness).all()
    inventory = (tf.V_sections[None, :, None] * liquid).sum(axis=1)
    inventory += (A * H[:, :, None, None]/L * biofilm).sum(axis=(1, 2))
    cumulative_input = sol.y[bio_stop+N:].T
    scale = np.maximum(np.abs(inventory[0]), 1.)
    np.testing.assert_allclose(
        (inventory - inventory[0] - cumulative_input)/scale, 0., atol=1e-7)
    assert np.all((H[-1] - H[0]) * (1 if case == 'growth' else -1) > 0.)


@pytest.mark.parametrize('N,L,case', [
    (1, 1, 'growth'), (2, 3, 'growth'), (2, 3, 'shrink'),
    (2, 3, 'reaction'), (6, 5, 'growth'),
    (1, 1, 'max_growth'), (2, 3, 'max_growth'),
    ])
def test_trickling_filter_preserves_layer_dry_density(N, L, case):
    """Normalized layers retain their density under unequal local sources."""
    import numpy as np
    tf = _make_thickness_balance_filter(N, L, case)
    cmps = tf.components
    p = np.array([c.particle_size == 'Particulate' for c in cmps])
    weights = cmps.i_mass[p]
    start = N*tf._ncol
    bio = tf._state[start:-N].reshape(N, L, len(cmps))
    rho = tf.biofilm_density*tf.biofilm_dry_content
    bio[:, :, p] = bio[:, :, p] * (rho/(bio[:, :, p] @ weights))[:, :, None]
    inlet = np.vstack([s.state for s in tf.ins])
    tf.ODE(0., inlet, tf._state.copy(), np.zeros_like(inlet))
    d_bio = tf._dstate[start:-N].reshape(N, L, len(cmps))
    np.testing.assert_allclose(d_bio[:, :, p] @ weights, 0., atol=1e-7)
    # Redistribution must not erase the component composition gradients.
    if L > 1:
        assert np.any(np.abs(d_bio[:, :, p]) > 0.)


def test_trickling_filter_redistribution_uses_donor_composition():
    """Attachment displaces outer-layer solids inward without mixing species."""
    import numpy as np
    tf = _make_thickness_balance_filter(1, 2, 'growth')
    tf.internal_solids_exchange_rate = 0.
    tf.diffusivity = np.zeros(len(tf.components))
    cmps = tf.components
    xi, xs = cmps.indices(('X_I', 'X_S'))
    n = len(cmps)
    rho = tf.biofilm_density*tf.biofilm_dry_content
    liquid = tf._state[:tf._ncol]
    liquid[:-1] = 0.
    liquid[xs] = 100.
    bio = tf._state[tf._ncol:-1].reshape(2, n)
    bio[:] = 0.
    bio[0, xi] = rho/cmps.i_mass[xi]
    bio[1, xs] = rho/cmps.i_mass[xs]
    inlet = np.vstack([s.state for s in tf.ins])
    tf.ODE(0., inlet, tf._state.copy(), np.zeros_like(inlet))
    d_bio = tf._dstate[tf._ncol:-1].reshape(2, n)
    # Half the attached dry mass expands each equal-width layer. The inner
    # layer receives the displaced X_I of the outer layer, not attached X_S.
    dry_attachment = tf.attachment_rate*100.*cmps.i_mass[xs]
    layer_depth = tf._state[-1]/2.
    np.testing.assert_allclose(d_bio[1, xi],
                               dry_attachment/2./cmps.i_mass[xi]/layer_depth)
    np.testing.assert_allclose(d_bio[0, xs], tf.attachment_rate*100./layer_depth)


def test_trickling_filter_density_at_capacity_during_integration():
    """Excess solids leave a full film without increasing its layer density."""
    import numpy as np
    from scipy.integrate import solve_ivp
    tf = _make_thickness_balance_filter(1, 2, 'growth')
    cmps = tf.components
    p = np.array([c.particle_size == 'Particulate' for c in cmps])
    start = tf._ncol
    n = len(cmps)
    rho = tf.biofilm_density*tf.biofilm_dry_content
    bio = tf._state[start:-1].reshape(2, n)
    bio[:, p] = bio[:, p] * (rho/(bio[:, p] @ cmps.i_mass[p]))[:, None]
    tf._state[-1] = 0.9999*tf.max_biofilm_thickness
    inlet = np.vstack([s.state for s in tf.ins])
    ode = tf.ODE
    def rhs(t, y):
        ode(t, inlet, y.copy(), np.zeros_like(inlet))
        return tf._dstate.copy()
    atol = np.full(tf._state.size, 1e-7)
    atol[-1] = 1e-11
    sol = solve_ivp(rhs, (0., 0.02), tf._state.copy(),
                    method='BDF', rtol=1e-7, atol=atol)
    assert sol.success, sol.message
    assert sol.y[-1, -1] >= tf.max_biofilm_thickness*(1.-1e-6)
    bio = sol.y[start:-1].T.reshape(-1, 2, n)
    np.testing.assert_allclose((bio[:, :, p] @ cmps.i_mass[p])/rho,
                               1., rtol=1e-5, atol=1e-7)
