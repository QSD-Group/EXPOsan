#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest


def test_missing_optional_htl_data_has_actionable_message():
    from exposan.htl._data import get_htl_data_path

    missing_file = 'missing_optional_dataset.xlsx'
    with pytest.raises(FileNotFoundError) as error:
        get_htl_data_path(missing_file)

    message = str(error.value)
    assert missing_file in message
    assert 'Large HTL optional data files are not included in the PyPI package' in message
    assert 'download_htl_data' in message


def test_download_htl_data_requires_url():
    from exposan.htl._data import download_htl_data

    with pytest.raises(ValueError) as error:
        download_htl_data('')

    assert 'URL' in str(error.value)
