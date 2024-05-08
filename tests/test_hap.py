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
    import numpy as np
    from numpy.testing import assert_allclose as ac
    from exposan.hap import create_model

    rtol = 0.01

    mdl = create_model()
    sys = mdl.system
    for p in mdl.parameters: 
        p.setter(p.baseline)
    sys.simulate()
    out = np.array([m() for m in mdl.metrics])
    ac(out, 
       [63172.32562288159,
        1.0195815674295157,
        14.107632178545428,
        -866376.4203627636,
        32.3546518781029,
        67.6453481218971,
        97.27280534057341,
        0.3445350538637175,
        2.382659605562873,
        31.347833235692512,
        20.180762873349565,
        22.224183497872545,
        20.27516971584288,
        5.972050677242491], 
       rtol=rtol)

if __name__ == '__main__':
    test_hap()