#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>
    Shion Watabe <shionwatabe@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np
from exposan import new_generator as ng
from exposan.new_generator import results_path
from exposan.new_generator.models import create_model, run_uncertainty

# ppls = np.arange(100, 600, 100) # 100 to 600 with a step size of 100
ppls = (200,) # if you just want to run one number
N_uncertainty = 100

if __name__ == '__main__':
    ng.INCLUDE_RESOURCE_RECOVERY = True
    for ppl in ppls:
        modelA = create_model('A', ppl=ppl)
        modelB = create_model('B', ppl=ppl)
        run_uncertainty(modelA, N=N_uncertainty, path=os.path.join(results_path, f'uncertaintyA_{ppl}users.xlsx'))
        run_uncertainty(modelB, N=N_uncertainty, path=os.path.join(results_path, f'uncertaintyB_{ppl}users.xlsx'))