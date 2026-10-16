"""Check ESAP cash-flow closure and product/N reporting bases."""


def test_esap_msp():
    from math import isclose
    from exposan.esap._systems import create_system
    from exposan.esap.model import calculate_msp

    for name in ('ammonium_sulfate', 'ammonium_citrate',
                 'ammonium_lactate', 'ammonium_maleate'):
        system = create_system(name, simulate=False)
        result = calculate_msp(system, discount_rate=0.1, lifetime=20,
                               operating_hours=8000)
        assert abs(result['NPV_USD']) < 1
        assert isclose(result['tea'].operating_hours, 8000)
        assert result['annual_product_kg'] > result['annual_N_kg'] > 0
        assert isclose(
            result['MSP_USD_per_kg_product'] * result['annual_product_kg'],
            result['MSP_USD_per_kg_N'] * result['annual_N_kg'],
        )
        if name == 'ammonium_lactate':
            # Use actual C3H9NO3 concentration, not the nominal 60 wt% target.
            product = system.flowsheet.unit.ProductConverter.outs[0]
            salt_fraction = product.imass['NH4Lactate'] / product.F_mass
            assert 0 < salt_fraction < 1
            expected_N_fraction = salt_fraction * 14.0067 / 107.10946
            assert isclose(result['annual_N_kg'] / result['annual_product_kg'],
                           expected_N_fraction, rel_tol=1e-4)
