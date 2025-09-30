# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Ali Ahmad <aa3056@scarletmail.rutgers.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''

__all__ = ('test_biobinder',)

from numpy.testing import assert_allclose
from qsdsan.utils import clear_lca_registries
from exposan.biobinder import create_system

EXPECTED = {
    'c-HTL_No_EC': {'MSP': -0.42, 'GWP': -6.2911},
    'c-HTL_EC':    {'MSP': 5.78, 'GWP': 2.3823},
    'd-HTL_No_EC': {'MSP': -0.68, 'GWP': -6.3264},
    'd-HTL_EC':    {'MSP': 5.80, 'GWP': 2.347},
}

def run_test(config_name, config_kwargs, rtol=0.15):
    clear_lca_registries()
    sys = create_system(**config_kwargs)
    sys.simulate()

    biobinder = sys.flowsheet.stream.biobinder
    tea = sys.TEA
    lca = sys.LCA

    MSP = tea.solve_price(biobinder)
    GWP = lca.get_allocated_impacts((biobinder,), operation_only=True, annual=True)['GWP']
    GWP /= biobinder.F_mass * lca.system.operating_hours

    print(f"{config_name} → MSP: {MSP:.2f}, GWP: {GWP:.4f}")

    ref = EXPECTED[config_name]
    assert_allclose(MSP, ref['MSP'], rtol=rtol)
    assert_allclose(GWP, ref['GWP'], rtol=rtol)

def test_biobinder():
    # CHCU no EC
    run_test('c-HTL_No_EC', dict(
        decentralized_HTL=False, decentralized_upgrading=False,
        skip_EC=True, generate_H2=False, EC_config=None,
    ))

    # CHCU with EC
    run_test('c-HTL_EC', dict(
        decentralized_HTL=False, decentralized_upgrading=False,
        skip_EC=False, generate_H2=False, EC_config=None,
    ))

    # DHCU no EC
    run_test('d-HTL_No_EC', dict(
        decentralized_HTL=True, decentralized_upgrading=False,
        skip_EC=True, generate_H2=False, EC_config=None,
    ))

    # DHCU with EC
    run_test('d-HTL_EC', dict(
        decentralized_HTL=True, decentralized_upgrading=False,
        skip_EC=False, generate_H2=False, EC_config=None,
    ))

if __name__ == '__main__':
    # test_biobinder() # temporarily remove test for distillation to be fixed
    from exposan import biobinder
