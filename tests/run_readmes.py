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
__all__ = (
    'run_bsm1_readme',
    'run_bwaise_readme',
    )

def run_bsm1_readme():
    from exposan import bsm1 as b1
    b1.bsm1.reset_cache()
    readme = os.path.join(b1.bsm1_path, 'README.rst')
    testfile(readme, module_relative=False)
    del b1

def run_bwaise_readme():
    from exposan import bwaise as bw
    readme = os.path.join(bw.bwaise_path, 'README.rst')
    testfile(readme, module_relative=False)
    del bw


if __name__ == '__main__':
    run_bsm1_readme()
    run_bwaise_readme()