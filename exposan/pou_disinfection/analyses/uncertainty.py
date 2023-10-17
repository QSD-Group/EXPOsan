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

import os, numpy as np
from exposan.pou_disinfection import (
    create_model,
    results_path,
    run_uncertainty,
    )

water_sources = ['GW', 'SW']
# water_sources = ['GW',]

sys_IDs = ['A', 'B', 'C', 'D']
# sys_IDs = ['A',]

ppls = np.arange(500, 1100, 100)
# ppls = [500, 1000]

N = 1000

if __name__ == '__main__':
    for sys_ID in sys_IDs:
        for water_source in water_sources:
            for ppl in ppls:
                model = create_model(
                    model_ID=sys_ID,
                    water_source=water_source,
                    ppl=ppl,
                    )
                path = os.path.join(results_path, f'{sys_ID}_{water_source}_ppl{ppl}_{N}.xlsx')
                run_uncertainty(model, path=path, N=N)
