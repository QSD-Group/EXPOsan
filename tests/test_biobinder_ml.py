# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

"""
Regression test for representative feedstocks in the biobinder_ml module.
"""

__all__ = ('test_biobinder_ml',)

import numpy as np
from numpy.testing import assert_allclose
from qsdsan.utils import clear_lca_registries
from exposan.biobinder_ml.Dist_flex import create_system, get_EBITDA


# IRR [-], NPV [USD], Revenue [USD/yr], EBITDA [USD/yr], GWP [kg CO2e/kg biobinder]
EXPECTED = {
    'food':   {'IRR': 0.2282213933, 'NPV': 42348257.21, 'Revenue': 11811350.36, 'EBITDA': 11769308.31, 'GWP': -6.404441273},
    'green':  {'IRR': 0.1423643610, 'NPV': 14207664.50, 'Revenue': 11811350.36, 'EBITDA':  8112550.925, 'GWP': -0.1568124352},
    'manure': {'IRR': 0.1014802320, 'NPV':   513459.6925, 'Revenue': 11811350.36, 'EBITDA': 6515878.854, 'GWP': -1.088927266},
    'sludge': {'IRR': 0.2936024966, 'NPV': 47415221.38, 'Revenue': 10991927.73, 'EBITDA': 10813439.99, 'GWP': -4.100078838},
}


def run_test(feedstock_id, rtol=0.01):
    clear_lca_registries()

    sys = create_system(
        feedstock_id=feedstock_id,
        decentralized_HTL=False,
        decentralized_upgrading=False,
        skip_EC=True,
        generate_H2=False,
        EC_config=None,
    )
    sys.simulate()

    tea = sys.TEA
    lca = sys.LCA
    biobinder = sys.flowsheet.stream.biobinder

    biobinder.price = 0.10

    GWP = lca.get_allocated_impacts(
        streams=(biobinder,),
        operation_only=True,
        annual=True,
    )['GWP']
    GWP /= biobinder.F_mass * sys.operating_hours

    results = {
        'IRR': tea.solve_IRR(),
        'NPV': tea.NPV,
        'Revenue': tea.sales,
        'EBITDA': get_EBITDA(tea),
        'GWP': GWP,
    }

    print(
    f"{feedstock_id}: "
    f"IRR={results['IRR'] * 100:.2f}%, "
    f"NPV=${results['NPV'] / 1e6:.2f} MM, "
    f"Revenue=${results['Revenue'] / 1e6:.2f} MM/yr, "
    f"EBITDA=${results['EBITDA'] / 1e6:.2f} MM/yr, "
    f"GWP={results['GWP']:.4f} kg CO2e/kg"
)

    for name, value in results.items():
        assert np.isfinite(value), f"{feedstock_id}: non-finite {name} ({value})"
        assert_allclose(
            value,
            EXPECTED[feedstock_id][name],
            rtol=rtol,
            err_msg=f"{feedstock_id}: {name}",
        )


def test_biobinder_ml():
    for feedstock_id in EXPECTED:
        run_test(feedstock_id)


if __name__ == '__main__':
    pass
    # test_biobinder_ml()