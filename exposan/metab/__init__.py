# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os

folder = os.path.dirname(__file__)
data_path = os.path.join(folder, 'data')
results_path = os.path.join(folder, 'results')
figures_path = os.path.join(folder, 'figures')
# To save simulation results and generated figures
if not os.path.isdir(results_path): os.mkdir(results_path)
if not os.path.isdir(figures_path): os.mkdir(figures_path)

from qsdsan import (
    ImpactIndicator as IInd, 
    ImpactItem as IItm
    )

_impact_item_loaded = False
def load_lca_data():
    global _impact_item_loaded
    if _impact_item_loaded:
        IInd.clear_registry()
        IItm.clear_registry()
    ind_path = os.path.join(data_path, 'TRACI_indicators.xlsx')
    itm_path = os.path.join(data_path, '_impact_items.xlsx')
    IInd.load_from_file(ind_path, sheet=0)
    IItm.load_from_file(itm_path)
    IItm('Stainless_steel', source='stainless_steel')
    _impact_item_loaded = True

from . import process
from .process import *

from . import equipment
from .equipment import *

from . import units
from .units import *

from . import systems
from .systems import *

# from . import models
# from .models import *

__all__ = (
    'folder',
    'data_path',
    'results_path',
    'figures_path',
    *process.__all__,
    *equipment.__all__,
    *units.__all__,
    *systems.__all__,
    # *models.__all__,
	)


