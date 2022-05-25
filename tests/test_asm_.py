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

__all__ = ('test_asm',)

def test_asm():
    import numpy as np
    from exposan import asm
    asm1_an = asm.create_asm_validation_system(process_model='ASM1', aerated=False)
    asm1_aer = asm.create_asm_validation_system(process_model='ASM1', aerated=True)
    asm2d_an = asm.create_asm_validation_system(process_model='ASM2d', aerated=False)
    asm2d_aer = asm.create_asm_validation_system(process_model='ASM2d', aerated=True)

    t = 10
    t_step = 0.5
    syses = (asm1_an, asm1_aer, asm2d_an, asm2d_aer)
    states = (
        dict(S_S=69.50000000001846, X_S=202.7347508802057, X_BH=27.719183825843597, S_NH=32.04570170429094),
        dict(S_S=64.32684709384519, X_S=196.16142052970068, X_BH=35.46730553472383, S_NH=31.745090457273598),
        dict(S_A=4.114072165874767, X_S=123.55649395125677, X_H=29.519897578137343, S_NH4=16.157541424757596),
        dict(S_A=0.578711484217068, X_S=118.32521252898434, X_H=37.24707568599702, S_NH4=15.914049030156331),
        )
    for sys, state in zip(syses, states):
        sys.simulate(state_reset_hook='reset_cache',
                     t_span=(0, t),
                     t_eval=np.arange(0, t+t_step, t_step),
                     method='BDF')
        CSTR = sys.flowsheet.unit.CSTR
        for k, v in state.items():
            assert np.isclose(v, CSTR.state[k], rtol=0.01)

if __name__ == '__main__':
    test_asm()