#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Jianan Feng <jiananf2@illinois.edu>
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, pandas as pd, qsdsan as qs
from qsdsan.utils import auom
from exposan.utils import _init_modules

htl_path = os.path.dirname(__file__)
module = os.path.split(htl_path)[-1]
data_path, results_path = _init_modules(module, include_data_path=True)

_kg_to_g = auom('kg').conversion_factor('g')
_m3perh_to_MGD = auom('m3/h').conversion_factor('MGD')
_MJ_to_MMBTU = auom('MJ').conversion_factor('MMBTU')
_MMgal_to_L = auom('gal').conversion_factor('L')*1000000


# %%

# =============================================================================
# Load components and systems
# =============================================================================

from . import _components
from ._components import *
_components_loaded = False
def _load_components(reload=False):
    global components, _components_loaded
    if not _components_loaded or reload:
        components = create_components()
        qs.set_thermo(components)
        _components_loaded = True

from . import _process_settings
from ._process_settings import *

from . import _tea
from ._tea import *

from . import systems
from .systems import *

_system_loaded = False
def load(configuration='baseline'):
    global sys, tea, lca, flowsheet, _system_loaded
    sys = create_system(configuration)
    tea = sys.TEA
    lca = sys.LCA
    flowsheet = sys.flowsheet
    _system_loaded = True
    dct = globals()
    dct.update(sys.flowsheet.to_dict())


def __getattr__(name):
    if not _components_loaded or not _system_loaded:
        raise AttributeError(
            f'Module {__name__} does not have the attribute "{name}" '
            'and the module has not been loaded, '
            f'loading the module with `{__name__}.load()` may solve the issue.')
        
from . import models
from .models import *

def simulate_and_save(model, N, rule='L', seed=None, path='',
                      include_spearman=True, **spearman_kwargs):
    samples = model.sample(N=N, rule=rule, seed=seed)
    model.load_samples(samples)
    model.evaluate()
    idx = len(model.parameters)
    parameters = model.table.iloc[:, :idx]
    results = model.table.iloc[:, idx:]
    percentiles = results.quantile([0, 0.05, 0.25, 0.5, 0.75, 0.95, 1])
    if include_spearman:
        kwargs = {'nan_policy': 'omit'}
        kwargs.update(spearman_kwargs)
        r_df, p_df = qs.stats.get_correlations(model, kind='Spearman', **kwargs)
    
    ID = model.system.flowsheet.ID
    N = model.table.shape[0]
    path = path or os.path.join(results_path, f'_{ID}_{N}.xlsx')
    with pd.ExcelWriter(path) as writer:
        parameters.to_excel(writer, sheet_name='Parameters')
        results.to_excel(writer, sheet_name='Results')
        percentiles.to_excel(writer, sheet_name='Percentiles')
        if include_spearman:
            r_df.to_excel(writer, sheet_name='Spearman_r')
            p_df.to_excel(writer, sheet_name='Spearman_p')
    writer.save()


__all__ = (
    'htl_path',
    'data_path',
    'results_path',
    'simulate_and_save',
    *_components.__all__,
    *_process_settings.__all__,
    *_tea.__all__,
    *systems.__all__,
    *models.__all__,
)



#%%
fig, ax = qs.stats.plot_uncertainties(model)
fig

#%%
fig, ax = qs.stats.plot_uncertainties(model, x_axis=model.metrics[0], y_axis=model.metrics[1],
                                      kind='kde-kde', center_kws={'fill': True})
fig