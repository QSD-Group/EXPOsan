# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>
    Zixuan Wang <wyatt4428@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, os
from warnings import warn
from .. import results_path

__all__ = ('cache_state', 'load_state')


def cache_state(sys, folder='', filename=''):
    filename = filename or sys.ID
    path = os.path.join(results_path, f'{folder}/{filename}.npy')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.save(path, sys._state)

def load_state(sys, state_arr=None, folder='steady_states', *, strict=True):
    """Load a dimension-compatible dynamic-state snapshot; return whether loaded.

    Check a finite, one-dimensional vector against the current assembled system
    before applying any cached values. With strict=True (default), missing or
    incompatible caches raise an informative error. With strict=False, warn and
    retain fresh initialization instead. Initialization/model errors still raise.

    Prefer the exact system ID's file; E-suffix configurations may fall back to
    their base ID when an exact file is absent. A matching shape is not proof of
    matching component order or calibration: regenerate caches after structural
    changes. Sludge-flow settings are not stored in the state vector.
    Unit state views are preserved; instrumental flow states are refreshed
    from the current connections by each unit's update routine.
    """
    problem = None
    if state_arr is None:
        filename = sys.ID
        path = os.path.join(results_path, f'{folder}/{filename}.npy')
        if not os.path.isfile(path) and filename.upper().endswith('E'):
            path = os.path.join(results_path, f'{folder}/{filename[:-1]}.npy')
        if not os.path.isfile(path) and not strict:
            problem = f'{sys.ID}: state cache not found: {path}.'
        else:
            state_arr = np.load(path)
    if problem is None:
        state_arr = np.asarray(state_arr, dtype=float)
        if state_arr.ndim != 1 or not np.isfinite(state_arr).all():
            problem = (f'{sys.ID}: cached state must be a finite one-dimensional vector; '
                       f'got shape {state_arr.shape}.')
            if strict:
                raise ValueError(problem)

    sys.converge()
    for unit in sys.units:
        unit._init_dynamic()
    sys._state = sys._state_idx = sys._DAE = None
    current, idx, nr = sys._load_state()
    if problem is None and state_arr.shape != current.shape:
        sizes = ', '.join(f'{ID}={stop-start}' for ID, (start, stop) in idx.items())
        problem = (
            f'{sys.ID}: cached state has {state_arr.size} entries, but the current '
            f'system requires {current.size} ({sizes}). Regenerate the cache with '
            'the current unit types, component sets and discretization.')
    if problem is not None:
        if strict:
            raise ValueError(problem)
        warn(problem + ' Using fresh initialization.', RuntimeWarning, stacklevel=2)
        return False

    # Preserve the state views established by the engine's own assembly routine.
    current[:] = state_arr
    units = sys.units[nr:] + sys.units[:nr]
    for unit in units:
        unit._update_state()
        unit._update_dstate()
    return True
