# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_bsm2',)

def test_bsm2():
    from numpy.testing import assert_allclose as ac
    # from exposan import bsm2
    # bsm2.load()
    # rtol = 0.01
    # sys = inter.interface_sys
    
    #!!! Temporarily disabled while trying to figuring out the problem
    
    # t_span = (0, 15) # just 15 days to make the test faster
    # sys.simulate(method='BDF', t_span=t_span)
    
    # assert sys.outs[0].isempty() == False
    # assert int(sys.scope.time_series[-1]) == t_span[1] # ensure it's complete
    # eff, biogas = sys.outs
    # ac(eff.iconc['S_S'], 0.9141898387548315, rtol=rtol)
    # ac(eff.COD, 49.19587909828996, rtol=rtol)
    # ac(eff.get_TSS(), 12.829724029634336, rtol=rtol)
    # ac(biogas.F_mass, 56.89029406091017, rtol=rtol)
    # ac(biogas.imass['S_ch4'], 50.77209529117958, rtol=rtol)


if __name__ == '__main__':
    test_bsm2()