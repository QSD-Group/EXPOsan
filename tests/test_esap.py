#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = (
    'test_esap_acid_consumption',
    'test_esap_crystallization_and_separation',
    'test_esap_service_prices',
    'test_esap_equipment_cost_index',
    'test_esap_evaporation_fraction',
    )


def test_esap_evaporation_fraction():
    """The overall evaporation fraction uses feed water, not total feed moles."""
    import math
    from types import SimpleNamespace
    from exposan.esap._systems import _get_evaporation_fraction

    components = SimpleNamespace(Salt=SimpleNamespace(MW=100.))
    # 1 M at 1000 L/hr and 1 kg/L: 100 kg salt + 900 kg water.
    # A 50 wt% target retains 100 kg water, removing 800/900 of feed water.
    def fraction(target, flow=1000., density=1.):
        return _get_evaporation_fraction('Salt', 1., flow, target, density, components)

    assert math.isclose(fraction(0.5), 8 / 9)
    assert fraction(0.1) == 0.  # Already at target.
    assert fraction(0.05) == 0.  # Feed exceeds target concentration.
    assert fraction(1.) == 0.95  # Retain the existing evaporation cap.
    assert fraction(0.5, flow=0.) == 0.
    assert fraction(0.5, density=0.1) == 0.  # No feed water.


def test_esap_equipment_cost_index():
    """All nonzero equipment purchase items follow the shared cost index."""
    import math
    import qsdsan as qs
    import biosteam as bst
    from exposan.esap._process_settings import equipment_cost_index
    from exposan.esap._systems import create_system

    for product in ('ammonium_sulfate', 'ammonium_citrate',
                    'ammonium_maleate', 'ammonium_lactate'):
        system = create_system(product, simulate=True)
        assert qs.CEPCI == bst.CE == equipment_cost_index
        costs = {(u.ID, item): cost for u in system.units
                 for item, cost in u.purchase_costs.items() if cost > 0}
        assert costs
        try:
            qs.CEPCI = equipment_cost_index * 1.1
            system.simulate()
            for unit in system.units:
                for item, cost in unit.purchase_costs.items():
                    if (unit.ID, item) in costs:
                        assert math.isclose(cost, costs[unit.ID, item] * 1.1,
                                            rel_tol=1e-7), (product, unit.ID, item)
        finally:
            qs.CEPCI = equipment_cost_index


def test_esap_service_prices():
    """Service rates are assigned before simulation without repeated escalation."""
    import math
    import qsdsan as qs
    from exposan.esap._systems import create_system

    expected = {
        'low_pressure_steam': (0., 0.2378 * 4.58 / 3.51),
        'cooling_water': (0., 4.8785e-4 * 0.1429 / 0.0676),
        'chilled_water': (5e-6 * 0.1429 / 0.0676, 0.),
    }
    # A different product on the second build must not compound shared rates.
    for product in ('ammonium_sulfate', 'ammonium_lactate'):
        create_system(product, simulate=False)
        for ID, (heat_price, regeneration_price) in expected.items():
            agent = qs.HeatUtility.get_agent(ID)
            assert math.isclose(agent.heat_transfer_price, heat_price)
            assert math.isclose(agent.regeneration_price, regeneration_price)


def test_esap_acid_consumption():
    import numpy as np
    from exposan.esap._systems import (
        _acid_components,
        _ammonium_per_product,
        _acid_per_ammonium,
        _products,
        create_system,
        )

    purchased_product_price = 2.5
    purchased_purity = 0.8
    products = (
        'ammonium_sulfate',
        'ammonium_citrate',
        'ammonium_lactate',
        'ammonium_maleate',
        )

    for product in products:
        sys = create_system(
            product,
            acid_prices={product: purchased_product_price},
            acid_purities={product: purchased_purity},
            simulate=False,
            )
        solute_ID = _products[product]['solute']
        acid_ID = _acid_components[product]
        feed = getattr(sys.flowsheet.stream, _products[product]['feed'])
        acid = getattr(sys.flowsheet.stream, f'{acid_ID}_solution')

        expected_acid_flow = (
            feed.imol[solute_ID]
            * _ammonium_per_product[product]
            * _acid_per_ammonium[product]
            )
        assert np.isclose(acid.imol[acid_ID], expected_acid_flow)
        nominal_solution_volume = acid.F_mass  # L/hr at the default 1 kg/L
        assert np.isclose(
            acid.imol[acid_ID] * 1000 / nominal_solution_volume,
            1,
            )
        assert np.isclose(
            acid.cost,
            acid.imass[acid_ID] / purchased_purity * purchased_product_price,
            )
        assert acid in sys.feeds


def test_esap_crystallization_and_separation():
    import numpy as np
    from biosteam.units.design_tools import size_batch
    from exposan.esap._sanunits import ESAPBatchCrystallizer
    from exposan.esap._systems import (
        _default_crystallizer_parameters,
        _products,
        create_system,
        )

    m3_to_gal = 1 / 0.003785411784
    for product in (
            'ammonium_sulfate',
            'ammonium_citrate',
            'ammonium_maleate',
            ):
        sys = create_system(product, simulate=True)
        crystallizer = sys.flowsheet.unit.Crystallizer
        assert isinstance(crystallizer, ESAPBatchCrystallizer)

        feed = crystallizer.ins[0]
        cake, mother_liquor = crystallizer.outs
        final_product = sys.flowsheet.unit.Dryer.outs[0]
        solute_ID = _products[product]['solute']
        crystal_ID = _products[product]['crystal']
        parameters = _default_crystallizer_parameters[product]
        formed = (
            feed.imol[solute_ID]
            * parameters['crystallization_recovery']
            )
        captured = formed * parameters['crystal_capture_efficiency']

        assert np.isclose(cake.imol[crystal_ID], captured)
        assert np.isclose(
            mother_liquor.imol[solute_ID],
            feed.imol[solute_ID] - formed,
            )
        assert np.isclose(
            mother_liquor.imol[crystal_ID], formed - captured,
            )
        assert np.isclose(
            cake.imass['H2O'] / cake.F_mass,
            parameters['wet_cake_moisture'],
            )
        assert np.isclose(final_product.imass[crystal_ID], cake.imass[crystal_ID])
        assert np.isclose(final_product.imass['H2O'] / final_product.F_mass, 0.01)
        assert np.isclose(crystallizer.mass_balance_error(), 0, atol=1e-9)
        assert all(
            np.isclose(i, 0, atol=1e-9)
            for i in crystallizer.atomic_balance_error().values()
            )

        flow = sum(i.F_vol for i in crystallizer.outs)
        expected = size_batch(
            flow,
            crystallizer.tau,
            crystallizer.tau_0,
            crystallizer.V_wf,
            N_reactors=crystallizer.N,
            )
        design = crystallizer.design_results
        assert design['Crystallizer volume'] > 0
        assert np.isclose(
            design['Crystallizer volume'], expected['Reactor volume'],
            )
        expected_agitation = (
            crystallizer.kW * crystallizer.V_wf
            * design['Crystallizer volume']
            * design['Number of crystallizers'] * m3_to_gal
            )
        expected_separation = parameters['separator_kWhr_per_m3'] * flow
        assert np.isclose(design['Agitator power'], expected_agitation)
        assert np.isclose(design['Separator power'], expected_separation)
        assert np.isclose(
            crystallizer.power_utility.rate,
            expected_agitation + expected_separation,
            )
        assert crystallizer.baseline_purchase_costs['Crystallizer'] > 0
        assert crystallizer.baseline_purchase_costs['Crystal separator'] > 0
        assert mother_liquor in sys.products


if __name__ == '__main__':
    test_esap_acid_consumption()
    test_esap_crystallization_and_separation()
