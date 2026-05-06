#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com> with Codex assistance
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from pathlib import Path
from urllib.request import urlretrieve

from . import data_path

__all__ = ('get_htl_data_path', 'download_htl_data')


_MISSING_OPTIONAL_DATA = (
    'Large HTL optional data files are not included in the PyPI package. '
    'Download the HTL data bundle, or use exposan.htl._data.download_htl_data '
    'if your workflow provides a data-bundle URL, and place the files in: '
    '{data_path}\n'
    'Missing file: {filename}'
)


def get_htl_data_path(filename):
    path = Path(data_path) / filename
    if not path.exists():
        raise FileNotFoundError(
            _MISSING_OPTIONAL_DATA.format(
                data_path=data_path,
                filename=filename,
            )
        )
    return path


def download_htl_data(url, filename=None):
    if not url:
        raise ValueError('A URL is required to download the HTL data bundle.')
    path = Path(data_path)
    path.mkdir(parents=True, exist_ok=True)
    target = path / (filename or Path(url).name)
    urlretrieve(url, target)
    return target
