# -*- coding: utf-8 -*-
"""SVCW construction, calibration, initialization and DW reporting tests.

This module is under the University of Illinois/NCSA Open Source License.
See https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt.
"""

from types import SimpleNamespace

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import numpy as np
import pandas as pd
import pytest
import qsdsan as qs

from exposan.esap.analyses.svcw_dw_validation import (
    COLORS, measured_concentration_ranges, plot_svcw_dw_validation,
)
from exposan.esap.systems import create_svcw_system


# %% Plant construction, calibration and AD/DW consistency
def test_svcw_tf_kinetics_are_independent_of_asr():
    system = create_svcw_system()
    units = system.flowsheet.unit
    tf_model = units.TF.suspended_growth_model
    asr_model = units.ASR.suspended_growth_model
    assert tf_model is not asr_model
    tf = tf_model.rate_function.params
    asr = asr_model.rate_function.params
    for name, value in dict(q_PHA=3.0, q_PP=1.5, b_PP=0.2, mu_H=6.0,
                            b_H=0.4, eta_NO3_Hl=0.5,
                            eta_NO3_PAOl=0.33, mu_AUT=1.0,
                            K_NH4_AUT=1.0, K_O2_AUT=0.5).items():
        assert tf[name] == pytest.approx(value)
    for name, value in dict(q_PHA=3.0, q_PP=1.5, b_PP=0.2, mu_H=6.0).items():
        assert asr[name] == pytest.approx(value)

    # Equal default rates must not conceal a shared mutable process model.
    tf_model.set_parameters(q_PHA=4.5)
    assert tf_model.rate_function.params['q_PHA'] == pytest.approx(4.5)
    assert asr_model.rate_function.params['q_PHA'] == pytest.approx(3.0)


def test_svcw_topology_and_defaults():
    system = create_svcw_system()
    unit = system.flowsheet.unit
    stream = system.flowsheet.stream

    assert system.ID == 'SVCW'
    assert unit.TF.ins[0] is unit.PC.outs[0]
    assert unit.TF.ins[0].ID == 'PE'
    assert unit.ASR.ins[0] is unit.TF.outs[0]
    assert unit.FC.ins[0] is unit.ASR.outs[0]
    assert unit.SF.ins[0] is unit.FC.outs[0]
    assert unit.SF.ins[0].ID == 'FC_SE'
    assert unit.SF.outs[0] is stream.SE
    assert unit.M2.ins[-1] is unit.SF.outs[1]
    assert unit.HD.outs[0] is unit.PC.ins[1]
    assert unit.FC.outs[1] is unit.ASR.ins[1]

    TF = unit.TF
    assert TF.N_section == 6
    assert TF.N_biofilm_layer == 4
    assert TF.filter_bed_surface == pytest.approx(9464.)
    assert TF.recirculation_ratio == pytest.approx(1.)
    assert TF.dynamic_biofilm_thickness is True
    assert TF.detachment_model == 'capacity'
    assert TF.detachment_exponent == pytest.approx(1.)
    assert TF.detachment_rate == pytest.approx(0.065)
    np.testing.assert_allclose(
        TF.initial_biofilm_thickness.sum(axis=1), 4e-4)

    SF = unit.SF
    assert SF.model == 'constant_splitting'
    assert SF.backwash_solids_mass_fraction == pytest.approx(0.6)
    assert SF.backwash_flow_fraction == pytest.approx(0.02)
    assert [i.ID for i in SF.outs] == ['SE', 'SF_backwash']

    rww = stream.RWW
    assert rww.get_total_flow('MGD') == pytest.approx(13.685)
    assert rww.COD == pytest.approx(600.)
    assert rww.iconc['S_NH4'] == pytest.approx(44.)
    assert rww.iconc['S_PO4'] == pytest.approx(2.0)
    assert rww.iconc['S_Na'] == pytest.approx(278.5)
    assert rww.iconc['S_K'] == pytest.approx(66.5)
    assert rww.iconc['S_Mg'] == pytest.approx(70.)
    assert rww.iconc['S_Ca'] == pytest.approx(650.)
    assert rww.iconc['S_Cl'] == pytest.approx(643.9)
    assert rww.iconc['S_I']/rww.COD == pytest.approx(0.14)
    assert rww.iconc['S_F']/rww.COD == pytest.approx(0.24)
    assert rww.iconc['S_A']/rww.COD == pytest.approx(0.05)
    assert rww.iconc['X_I']/rww.COD == pytest.approx(0.235)
    assert rww.BOD == pytest.approx(241.338, rel=1e-6)
    assert rww.get_TSS() == pytest.approx(215.46, rel=1e-6)
    assert rww.components.X_S.i_N == pytest.approx(0.023464)
    assert rww.components.X_I.i_N == pytest.approx(0.042021)
    assert rww.components.S_F.i_P == pytest.approx(0.006149)
    assert rww.components.X_S.i_P == pytest.approx(0.006149)
    assert rww.components.S_I.i_P == pytest.approx(0.007139)
    assert rww.components.X_I.i_P == pytest.approx(0.007139)
    for component_id in ('X_H', 'X_PAO', 'X_AUT'):
        assert getattr(rww.components, component_id).i_P == pytest.approx(0.02369)
    assert rww.TP == pytest.approx(5.72768)
    assert unit.AD.V_liq == pytest.approx(0.85*1.3685*3785.412)
    assert unit.AD.V_gas == pytest.approx(0.085*1.3685*3785.412)


def test_svcw_calibration_inputs_are_overridable():
    system = create_svcw_system(
        default_init_conds=False,
        wastewater_flow_mgd=12.,
        wastewater_kwargs=dict(COD=600., S_Na=200.),
        organic_fractions=dict(fr_SI=0.10, fr_XI=0.20),
        particulate_n_factor=1.,
        ad_liquid_volume_mgd=0.75,
        )
    unit = system.flowsheet.unit
    rww = system.flowsheet.stream.RWW

    assert rww.get_total_flow('MGD') == pytest.approx(12.)
    assert rww.COD == pytest.approx(600.)
    assert rww.iconc['S_Na'] == pytest.approx(200.)
    assert rww.iconc['S_I']/rww.COD == pytest.approx(0.10)
    assert rww.iconc['X_I']/rww.COD == pytest.approx(0.20)
    assert rww.components.X_S.i_N == pytest.approx(0.023464)
    assert rww.components.X_I.i_N == pytest.approx(0.042021)
    assert unit.AD.V_liq == pytest.approx(0.75*3785.412)
    assert unit.AD.V_gas == pytest.approx(0.075*3785.412)


def test_svcw_dewaters_adm1p_before_the_return_junction():
    system = create_svcw_system()
    unit = system.flowsheet.unit
    adm_ids = unit.AD.components.IDs
    asm_ids = unit.M1.components.IDs

    assert unit.DW.ins[0] is unit.AD.outs[1]
    assert unit.J2.ins[0] is unit.DW.outs[0]
    assert unit.M2.ins[2] is unit.J2.outs[0]
    assert unit.DW.components.IDs == adm_ids
    assert unit.DW._mixed.components.IDs == adm_ids
    assert all(stream.components.IDs == adm_ids for stream in unit.DW.outs)
    assert unit.J2.outs[0].components.IDs == asm_ids
    assert unit.M2.components.IDs == asm_ids


def test_svcw_dw_preserves_dissolved_nutrients_and_component_balances():
    system = create_svcw_system()
    unit = system.flowsheet.unit
    dw = unit.DW
    feed = dw.ins[0]
    feed.set_flow_by_concentration(
        flow_tot=200., units=('m3/d', 'mg/L'),
        concentrations=dict(S_IN=1000., S_IP=40., S_IC=3000.,
                            S_Mg=20., S_Ca=65., X_su=1000., X_PP=100.),
        )
    dw._run()
    np.testing.assert_allclose(
        dw.outs[0].mass + dw.outs[1].mass, feed.mass, atol=1e-10,
        )

    # Exercise the algebraic split used inside the integrated dynamic run.
    unit.J2._init_dynamic()
    dw._init_dynamic()
    dw._init_state()
    initial = dw._state.copy()
    dw.AE(0., initial[None, :], np.zeros((1, len(initial))))
    liquid, cake = (stream.state.copy() for stream in dw.outs)
    cmps = dw.components
    dissolved = cmps.indices(('S_IN', 'S_IP', 'S_Mg', 'S_Ca'))
    np.testing.assert_array_equal(liquid[dissolved], initial[dissolved])
    np.testing.assert_array_equal(cake[dissolved], initial[dissolved])
    np.testing.assert_allclose(
        liquid[:-1]*liquid[-1] + cake[:-1]*cake[-1],
        initial[:-1]*initial[-1], atol=1e-9,
        )
    particulate = cmps.indices(('X_su', 'X_PP'))
    np.testing.assert_allclose(
        liquid[particulate], initial[particulate]*(1-dw.solids_removal_efficiency),
        )
    assert cake[-1] == pytest.approx(dw.sludge_flow_rate)

    mapped = unit.J2.reactions(liquid[:-1].copy())
    np.testing.assert_allclose(
        np.dot(mapped, unit.J2.outs[0].components.i_P),
        np.dot(liquid[:-1], cmps.i_P), atol=1e-10,
        )
    assert mapped[unit.J2.outs[0].components.index('S_PO4')] > liquid[cmps.index('S_IP')]


def test_svcw_rejects_empirical_sand_filter_topology():
    with pytest.raises(ValueError, match='constant_splitting'):
        create_svcw_system(
            default_init_conds=False,
            sf_kwargs=dict(model='empirical_fouling'))


def test_svcw_ad_initial_phosphate_and_component_indices():
    system = create_svcw_system()
    units = system.flowsheet.unit
    ad = units.AD
    adm_components = ad.components
    asm_components = units.J1.ins[0].components

    assert ad._concs[adm_components.index('S_IP')] == pytest.approx(24.5)
    assert np.max(np.abs(system.adm1p_p_balance_residual)) < 1e-10
    np.testing.assert_allclose(
        ad.model.rate_function.params['k_mmp'],
        (0.24, 240.0, 0.24, 72.0, 0.24, 0.24, 0.24),
        )

    for component_id in ('S_I', 'X_I', 'X_PHA', 'X_PP', 'X_PAO'):
        assert getattr(adm_components, component_id).i_P == pytest.approx(
            getattr(asm_components, component_id).i_P)

    assert system.adm1p_component_sync['X_PAO']['i_P'][1] == pytest.approx(
        asm_components.X_PAO.i_P)
    qs.main_flowsheet.clear()


def test_svcw_ad_initial_phosphate_is_overridable():
    system = create_svcw_system(
        ad_initial_s_ip=18.75,
        equilibrate_ad_minerals=False,
        )
    ad = system.flowsheet.unit.AD
    assert ad._concs[ad.components.index('S_IP')] == pytest.approx(18.75)
    qs.main_flowsheet.clear()


@pytest.mark.parametrize('mineral_settings', [
    {},
    dict(use_adm1p_default_pKsp=True),
    dict(ad_caco3_pKsp=8.5, ad_struvite_pKsp=13.7,
         ad_acp_pKsp=28.6, ad_mgco3_pKsp=7.6),
    dict(ad_caco3_pKsp=8.500001, ad_struvite_pKsp=13.7,
         ad_acp_pKsp=28.6, ad_mgco3_pKsp=7.6),
    dict(ad_caco3_pKsp=5.5, ad_struvite_pKsp=11.6,
         ad_acp_pKsp=23.3, ad_mgco3_pKsp=5.6),
])
@pytest.mark.parametrize('initialization_settings', [
    {}, dict(equilibrate_ad_minerals=False), dict(ad_pH=None),
])
def test_svcw_ad_initial_metals_are_independent_of_pKsp(
        mineral_settings, initialization_settings):
    with qs.Flowsheet('svcw_initial_metals'):
        system = create_svcw_system(
            **mineral_settings, **initialization_settings,
            wastewater_kwargs=dict(S_Mg=23., S_Ca=45.),
            )
        ad = system.flowsheet.unit.AD
        indices = ad.components.indices(('S_Mg', 'S_Ca'))
        np.testing.assert_array_equal(ad._concs[indices], (10., 80.))

        ad._init_state()
        np.testing.assert_allclose(
            ad._state[indices] * 1e3, (10., 80.), rtol=0., atol=1e-12,
            )
        rww = system.flowsheet.stream.RWW
        assert rww.iconc['S_Mg'] == pytest.approx(23.)
        assert rww.iconc['S_Ca'] == pytest.approx(45.)


@pytest.mark.parametrize('value', [-1., np.nan, np.inf])
def test_svcw_rejects_invalid_ad_initial_phosphate(value):
    with pytest.raises(ValueError, match='ad_initial_s_ip'):
        create_svcw_system(
            default_init_conds=False,
            ad_initial_s_ip=value,
            )
    qs.main_flowsheet.clear()


# %% DW plotting fixtures, nutrient mass bases and measured ranges
@pytest.fixture(params=['asm', 'adm'])
def dw_plot_system(request):
    nitrogen, phosphorus = (
        ('S_IN', 'S_IP') if request.param == 'adm' else ('S_NH4', 'S_PO4'))
    columns = pd.MultiIndex.from_tuples([
        ('-', 't [d]'),
        ('DW_liquid', f'{nitrogen} [mg/L]'),
        ('DW_liquid', 'S_Mg [mg/L]'),
        ('DW_liquid', 'S_Ca [mg/L]'),
        ('DW_liquid', f'{phosphorus} [mg/L]'),
    ])
    frame = pd.DataFrame([
        [-1, 9999., 999., 999., 999.],
        [0, 1000., 10., 60., 30.],
        [100, 1100., 15., 63., 40.],
        [300, 1200., 20., 65., 43.],
        [301, 9999., 999., 999., 999.],
    ], columns=columns)
    return SimpleNamespace(
        scope=SimpleNamespace(export=lambda: frame),
        flowsheet=SimpleNamespace(unit=SimpleNamespace(
            DW=SimpleNamespace(outs=[SimpleNamespace(ID='DW_liquid')]))),
    )


def test_measured_ranges_use_element_mass():
    ranges = measured_concentration_ranges()
    assert ranges['S_NH4']['mean'] == pytest.approx(1380.06*14/18)
    assert ranges['S_NH4']['lower'] == pytest.approx(1242.16*14/18)
    assert ranges['S_NH4']['upper'] == pytest.approx(1517.96*14/18)
    assert ranges['S_PO4']['mean'] == pytest.approx(42.78)
    assert ranges['S_PO4']['lower'] == pytest.approx(129.8*31/95)
    assert ranges['S_PO4']['upper'] == pytest.approx(132.4*31/95)
    assert ranges['S_Mg']['lower'] == pytest.approx(17.7)
    assert ranges['S_Mg']['upper'] == pytest.approx(20.9)
    assert ranges['S_Ca']['lower'] == pytest.approx(62.8)
    assert ranges['S_Ca']['upper'] == pytest.approx(67.)


@pytest.mark.parametrize('basis', ['element', 'ion'])
def test_dw_panels_and_bands_share_one_mass_basis(dw_plot_system, basis):
    original_font = plt.rcParams['font.sans-serif'].copy()
    fig, axes = plot_svcw_dw_validation(dw_plot_system, concentration_basis=basis, show=False)
    try:
        assert len(fig.axes) == 2
        assert [len(ax.lines) for ax in axes] == [1, 3]
        assert [len(ax.patches) for ax in axes] == [1, 3]
        np.testing.assert_array_equal(axes[0].lines[0].get_xdata(), [0, 100, 300])
        expected_n = np.array([1000., 1100., 1200.])
        expected_p = np.array([30., 40., 43.])
        if basis == 'ion':
            expected_n *= 18/14
            expected_p *= 95/31
        np.testing.assert_allclose(axes[0].lines[0].get_ydata(), expected_n)
        np.testing.assert_allclose(axes[1].lines[2].get_ydata(), expected_p)
        ranges = measured_concentration_ranges(basis)
        for ax, components in zip(axes, [('S_NH4',), ('S_Mg', 'S_Ca', 'S_PO4')]):
            for patch, line, component in zip(ax.patches, ax.lines, components):
                bbox = patch.get_path().transformed(patch.get_patch_transform()).get_extents()
                assert bbox.ymin == pytest.approx(ranges[component]['lower'])
                assert bbox.ymax == pytest.approx(ranges[component]['upper'])
                assert patch.get_alpha() == 0.5
                assert patch.get_facecolor() == pytest.approx(to_rgba(COLORS[component], 0.5))
                assert line.get_color() == COLORS[component]
                assert line.get_linewidth() == 3
            assert ax.get_xlim() == pytest.approx((0, 300))
            assert ax.get_ylim()[1] >= max(ranges[c]['upper'] for c in components)
        assert plt.rcParams['font.sans-serif'] == original_font
    finally:
        plt.close(fig)


def test_dw_plot_rejects_missing_records(dw_plot_system):
    with pytest.raises(ValueError, match='No recorded simulation times'):
        plot_svcw_dw_validation(dw_plot_system, start=400, stop=500, show=False)


def test_dw_plot_rejects_a_missing_component(dw_plot_system):
    frame = dw_plot_system.scope.export()
    phosphorus = [column for column in frame.columns
                  if column[1].split(' [', 1)[0] in ('S_IP', 'S_PO4')]
    frame = frame.drop(columns=phosphorus)
    dw_plot_system.scope.export = lambda: frame
    with pytest.raises(KeyError, match='S_PO4'):
        plot_svcw_dw_validation(dw_plot_system, show=False)


def test_dw_plot_rejects_ambiguous_nutrient_aliases(dw_plot_system):
    frame = dw_plot_system.scope.export().copy()
    nitrogen = next(column for column in frame.columns
                    if column[1].split(' [', 1)[0] in ('S_IN', 'S_NH4'))
    other = 'S_NH4' if nitrogen[1].startswith('S_IN ') else 'S_IN'
    frame[('DW_liquid', f'{other} [mg/L]')] = frame[nitrogen]
    dw_plot_system.scope.export = lambda: frame
    with pytest.raises(KeyError, match='S_NH4'):
        plot_svcw_dw_validation(dw_plot_system, show=False)
