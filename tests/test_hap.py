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
         # 1.048288575203588,
         17.597048423328086,
         -1086811.1628681163,
         25.792040993957844,
         74.20795900604216,
         97.2731791145954,
         0.3441417797725412,
         2.38267910563208,
         49.52777806264908,
         15.221503548445778,
         16.166861157662943,
         14.744119672295634,
         4.339737558946575], 
        rtol=rtol)

if __name__ == '__main__':
    test_hap()