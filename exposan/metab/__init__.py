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

from qsdsan import (
    ImpactIndicator as IInd, 
    ImpactItem as IItm
    )

_impact_item_loaded = False
def _load_lca_data():
    global _impact_item_loaded
    if _impact_item_loaded:
        IItm.clear_registry()  # items are per-flowsheet; indicators are global
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

from . import models
from .models import *


_loaded_config = None
def load(reload=False, n_stages=1, reactor_type='UASB', gas_extraction='P',
         tot_HRT=12, simulate=False, t_span=(0, 400), method='BDF', **kwargs):
    '''
    Construct a METAB configuration and expose its units and streams at the
    module level (e.g. ``metab.sys``).

    Parameters
    ----------
    reload : bool
        Rebuild even if an identical configuration is already loaded.
    n_stages : int
        Number of reactor stages.
    reactor_type : str
        Reactor type, e.g. 'UASB', 'FB', or 'PB'.
    gas_extraction : str
        Biogas extraction mode, e.g. 'P' (passive), 'M' (membrane), or 'H'
        (headspace).
    tot_HRT : float
        Total hydraulic retention time [h].
    simulate : bool
        If True, also run the dynamic simulation. Defaults to False (construct
        only) because the run is long (t_span=(0, 400), BDF); pass simulation
        kwargs via ``t_span``/``method``.
    **kwargs
        Forwarded to ``create_system`` (e.g. lifetime, T, Q, inf_concs).
    '''
    global sys, _loaded_config
    config = (n_stages, reactor_type, gas_extraction, tot_HRT)
    if reload or _loaded_config != config:
        sys = create_system(n_stages=n_stages, reactor_type=reactor_type,
                            gas_extraction=gas_extraction, tot_HRT=tot_HRT, **kwargs)
        _loaded_config = config
    if simulate:
        sys.simulate(t_span=t_span, method=method, state_reset_hook='reset_cache')
    globals().update(sys.flowsheet.to_dict())
    return sys


__all__ = (
    'folder',
    'data_path',
    'results_path',
    'figures_path',
    *process.__all__,
    *equipment.__all__,
    *units.__all__,
    *systems.__all__,
    *models.__all__,
	)


