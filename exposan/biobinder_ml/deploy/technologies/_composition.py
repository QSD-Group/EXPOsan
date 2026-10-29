# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Access deployment-local feedstocks without initializing the HTL/ML stack."""
import importlib.util
from pathlib import Path
import sys


def feedstock_module():
    name = 'exposan.biobinder_ml.deploy.feedstocks'
    if name not in sys.modules:
        path = Path(__file__).resolve().parents[1]/'feedstocks.py'
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


def characterize(name, overrides=None):
    module = feedstock_module()
    record = dict(module.BIOBINDER_FEEDSTOCKS[name])
    record.update(overrides or {})
    if overrides and 'ad_biodegradability' in overrides:
        record['ad_biodegradability_status'] = 'explicit characterization override'
        record['ad_biodegradability_source'] = overrides.get('ad_biodegradability_source','')
        record['ad_biodegradability_substrate'] = overrides.get('ad_biodegradability_substrate','')
        record['ad_biodegradability_basis'] = overrides.get('ad_biodegradability_basis','explicit fraction')
    return record, module.get_recovery_composition(record)
