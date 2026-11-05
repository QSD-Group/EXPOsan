# -*- coding: utf-8 -*-
"""Sand-filter feed definitions, separation and dynamic recycle tests.

This module is under the University of Illinois/NCSA Open Source License.
See https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt.
"""

import numpy as np
import pytest

import exposan.esap.systems.SandFilter_validation as sfv
import exposan.esap.systems.TricklingFilter_validation as tfv


# %% Shared wastewater definitions and filtration behavior
def test_sand_filter_validation_reuses_wastewater_definitions():
    assert sfv._WASTEWATER_ORGANIC_FRACTIONS \
        is tfv._WASTEWATER_ORGANIC_FRACTIONS
    assert sfv._WASTEWATER_STRENGTHS is tfv._WASTEWATER_STRENGTHS
    assert sfv._NONDEFAULT_COMPONENT_PROPERTIES \
        is tfv._NONDEFAULT_COMPONENT_PROPERTIES

    system = sfv.create_sand_filter_validation_system(
        wastewater_strength='high', sf_kwargs=dict(isdynamic=False))
    SF = system.flowsheet.unit.SF
    rww = system.flowsheet.stream.RWW

    assert rww.COD == pytest.approx(
        sfv._WASTEWATER_STRENGTHS['high']['COD'])
    for component_ID, expected in sfv._NONDEFAULT_COMPONENT_PROPERTIES.items():
        component = getattr(SF.components, component_ID)
        for property_name, value in expected.items():
            assert getattr(component, property_name) == pytest.approx(value)


def test_sand_filter_validation_constant_splitting():
    system = sfv.create_sand_filter_validation_system(
        wastewater_strength='mid', sf_kwargs=dict(isdynamic=False))
    SF = system.flowsheet.unit.SF
    SF.simulate()

    assert SF.model == 'constant_splitting'
    assert len(SF.ins) == 1
    assert [stream.ID for stream in SF.outs] == ['effluent', 'backwash']
    np.testing.assert_allclose(
        sum(stream.mass for stream in SF.outs), SF.ins[0].mass,
        rtol=1e-12, atol=1e-9)
    # assert SF.outs[1].F_vol/SF.ins[0].F_vol == pytest.approx(0.02)
    assert SF.outs[1].imass['H2O']/SF.ins[0].imass['H2O'] \
        == pytest.approx(0.02)


def test_sand_filter_validation_empirical_recycle():
    system = sfv.create_sand_filter_validation_system(
        wastewater_strength='low', model='empirical_fouling')
    SF = system.flowsheet.unit.SF
    recycle = system.flowsheet.stream.filtered_water_recycle

    assert SF.ins[1] is recycle
    assert SF.outs[2] is recycle
    assert [stream.ID for stream in SF.outs] == [
        'effluent', 'backwash', 'filtered_water_recycle']

    system.simulate(
        method='BDF', t_span=(0., 0.1), t_eval=(0., 0.1),
        rtol=1e-5, atol=1e-8)
    assert system.scope.sol.success, system.scope.sol.message
    assert recycle.state[-1] == pytest.approx(
        SF.ins[0].state[-1]*SF.average_backwash_fraction)
    assert np.isfinite(SF.captured_solids_inventory).all()
    assert (SF.captured_solids_inventory >= -1e-8).all()
