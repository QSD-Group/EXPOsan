# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_hap',)

def test_hap():
    from numpy.testing import assert_allclose as ac
    from exposan.hap import create_model

    rtol = 0.01

    mdl = create_model()
    for p in mdl.parameters: 
        p.setter(p.baseline)
    mdl.system.simulate()
    out = [m() for m in mdl.metrics]
    transport_duty = out.pop(1)
    
    #!!! disabled test w.r.t. CVRP solution until figure out random seed setting
    # import numpy as np
    # assert np.isclose(transport_duty, 0.966498304330359, rtol=0.1)
    
    ac(out, 
        [63172.32562288159,
        # 0.966498304330359,
        14.103592623653629,
        -866121.2322857676,
        32.364184632972425,
        67.63581536702758,
        97.27280534057341,
        0.3445350538637175,
        2.382659605562873,
        31.36148888930313,
        20.189553959652866,
        22.19030301875137,
        20.284001927385358,
        5.974652204907268], 
        rtol=rtol)

if __name__ == '__main__':
    test_hap()