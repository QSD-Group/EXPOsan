#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os
from exposan.pou_disinfection import (
    create_model,
    results_path,
    run_uncertainty,
    )

sys_ID = 'A'
water_source = 'GW'
ppl = 1000
N = 100

modelA = create_model(sys_ID, water_source=water_source, ppl=ppl)

if __name__ == '__main__':
    path = os.path.join(results_path, f'{sys_ID}_{water_source}_ppl{ppl}_{N}.xlsx')
    run_uncertainty(modelA, path=path, N=N)
