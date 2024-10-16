#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Ga-Yeong Kim <gayeong1225@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_pm2',)

def test_pm2():
    # Batch mode
    from exposan import pm2_batch
    pm2_batch.load()
    sys = pm2_batch.sys
    sys.simulate(t_span=(0,1), method='RK23')
    PBR = pm2_batch.PBR
    fig, ax = PBR.scope.plot_time_series(('S_P'))
    
    # EcoRecover
    from exposan import pm2_ecorecover
    pm2_ecorecover.load()
    sys = pm2_ecorecover.sys
    sys.simulate(t_span=(0,1), method='RK23')
    PBR = pm2_ecorecover.PBR20
    fig, ax = PBR.scope.plot_time_series(('S_P'))


if __name__ == '__main__':
    test_pm2()