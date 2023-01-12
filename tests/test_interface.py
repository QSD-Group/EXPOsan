# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_interface',)

def test_interface():
    from numpy.testing import assert_allclose as ac
    from exposan import interface as inter
    inter.load()
    sys = inter.interface_sys
    t_span = (0, 15) # just 15 days to make the test faster
    sys.simulate(method='BDF', t_span=t_span)
    
    assert sys.outs[0].isempty() == False
    assert int(sys.scope.time_series[-1]) == t_span[1] # ensure it's complete
    eff, biogas = sys.outs
    ac(eff.iconc['S_S'], 0.9141898387548315, rtol=1e-2)
    ac(eff.COD, 49.19587909828996, rtol=1e-2)
    ac(eff.get_TSS(), 12.829724029634336, rtol=1e-2)
    ac(biogas.F_mass, 56.89029406091017, rtol=1e-2)
    ac(biogas.imass['S_ch4'], 50.77209529117958, rtol=1e-2)


if __name__ == '__main__':
    test_interface()