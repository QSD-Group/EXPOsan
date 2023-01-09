# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_bsm1',)

def test_bsm1():
    from numpy.testing import assert_allclose as ac
    from numpy import arange, array, vstack
    from exposan import bsm1
    bsm1_ss_concs = array([
        [30, 2.81, 1149., 82.1, 2552., 148., 449., 0.0043, 5.37, 7.92, 1.22, 5.28, 4.93*12],
        [30, 1.46, 1149., 76.4, 2553., 148., 450., 6.31e-5, 3.66, 8.34, 0.882, 5.03, 5.08*12],
        [30, 1.15, 1149., 64.9, 2557., 149., 450., 1.72, 6.54, 5.55, 0.829, 4.39, 4.67*12],
        [30, 0.995, 1149., 55.7, 2559., 150., 451., 2.43, 9.3, 2.97, 0.767, 3.88, 4.29*12],
        [30, 0.889, 1149., 49.3, 2559., 150., 452., 0.491, 10.4, 1.73, 0.688, 3.53, 4.13*12]
        ])
    bsm1_ss_c1tss = array([12.5, 18.1, 29.5, 69.0, 356., 356., 356., 356., 356., 6394])
    bsm1.load()
    sys = bsm1.sys
    sys.simulate(t_span=(0,100), method='BDF', state_reset_hook='reset_cache')
    assert sys.outs[0].isempty() == False
    # ac(float(sys.outs[0].iconc['S_S']), 0.895, rtol=1e-2)
    # ac(float(sys.outs[1].iconc['X_BH']), 4994.3, rtol=1e-2)
    # ac(sys.outs[0].COD, 47.5, rtol=1e-2)
    # ac(sys.outs[1].get_TSS(), 6377.9, rtol=1e-2)
    u = sys.flowsheet.unit
    concs = vstack((u.A1._state[:13], u.A2._state[:13], u.O1._state[:13], u.O2._state[:13], u.O3._state[:13]))
    ac(concs, bsm1_ss_concs, rtol=2e-2)
    ac(u.C1._state[-10:], bsm1_ss_c1tss, rtol=2e-2)

if __name__ == '__main__':
    test_bsm1()