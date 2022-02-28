#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_cas',)

#TODO: for now just make sure it can run, add more costing, etc
def test_cas():
    from qsdsan import set_thermo
    from exposan import activated_sludge as cas

    set_thermo(cas.cmps)
    cas.cas.simulate()


if __name__ == '__main__':
    test_cas()