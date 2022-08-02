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

# =============================================================================
# Note: keep getting an error about:
#     E       fixture 'filename' not found
# thus avoid using "test" in naming to avoid being auto-discovered by `pytest`
#
# Run this file directly for testing, no need for frequent test
# =============================================================================

import os
from doctest import testfile


def run_bsm1():
    from exposan import bsm1 as b1
    b1.bsm1.reset_cache()
    readme = os.path.join(b1.bsm1_path, 'README.rst')
    testfile(readme, module_relative=False)


def run_biogenic_refinery():
    from exposan.biogenic_refinery import br_path
    readme = os.path.join(br_path, 'README.rst')
    testfile(readme, module_relative=False)


def run_bwaise():
    from exposan.bwaise import bw_path
    readme = os.path.join(bw_path, 'README.rst')
    testfile(readme, module_relative=False)


def run_eco_san():
    from exposan.eco_san import es_path
    readme = os.path.join(es_path, 'README.rst')
    testfile(readme, module_relative=False)


def run_reclaimer():
    from exposan.reclaimer import re_path
    readme = os.path.join(re_path, 'README.rst')
    testfile(readme, module_relative=False)


if __name__ == '__main__':
    run_bsm1()
    run_biogenic_refinery()
    run_bwaise()
    run_eco_san()
    run_reclaimer()