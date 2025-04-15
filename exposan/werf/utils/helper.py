# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, os
from .. import results_path

__all__ = ('cache_state', 'load_state')


def cache_state(sys, folder=''):
    path = os.path.join(results_path, f'{folder}/{sys.ID}.npy')
    np.save(path, sys._state) 
    
def load_state(sys, state_arr=None, folder='steady_states'):
    """Load a previously saved state to a system object for faster model evaluation."""
    if state_arr is None:
        path = os.path.join(results_path, f'{folder}/{sys.ID}.npy')
        state_arr = np.load(path)
    nr = sys._n_rotate
    units = sys.units[nr:] + sys.units[:nr]
    sys.converge()
    for unit in sys.units: 
        unit._init_dynamic()
    for ws in sys.feeds:
        if not ws.state.all(): ws._init_state()
    if sys.recycle:
        for ws in sys.recycle:
            if not ws.state.all(): ws._init_state()
    for inf in units[0].ins:
        if not inf.state.all(): inf._init_state()
        y = np.array([])
        idx = {}
        for unit in units: 
            unit._init_state()
            if unit.hasode:
                start = len(y)
                y = np.append(y, unit._state)
                stop = len(y)
                idx[unit._ID] = (start, stop)
                unit._state = state_arr[start: stop]
            unit._update_state()
            unit._update_dstate()
        assert len(y) == len(state_arr)
        sys._state = state_arr
        sys._state_idx = idx

