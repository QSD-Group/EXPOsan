# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np, qsdsan as qs
from qsdsan import System
from qsdsan.sanunits import ADMtoASM, ASMtoADM
from exposan.bsm1 import bsm1 as bsm1_sys, A1, C1, RE as RWW, WAS
from exposan.interface import results_path

__all__ = ('create_system',)


def create_system():
    thermo_asm1 = qs.get_thermo() # ASM1 components loaded by the bsm1 module
    cmps_asm1 = thermo_asm1.chemicals
    
    # Subsequent units should be using ADM1 components
    cmps_adm1 = qs.processes.create_adm1_cmps()
    thermo_adm1 = qs.get_thermo()
    adm1 = qs.processes.ADM1()
    cmps_adm1.X_I.i_N = cmps_asm1.X_I.i_N
    
    
    J1 = ASMtoADM('J1', upstream=WAS, thermo=thermo_adm1, isdynamic=True, adm1_model=adm1) # WAS is C1.outs[2]
    temp = lambda t: 293.15
    AD1 = qs.sanunits.AnaerobicCSTR('AD1', ins=J1.outs[0], outs=('biogas', 'ad_eff'), isdynamic=True,
                                    model=adm1, exogenous_var=(temp,))
    J2 = ADMtoASM('J2', upstream=AD1-1, thermo=thermo_asm1, isdynamic=True, adm1_model=adm1)
    
    # Subsequent units should be using ASM1 components
    qs.set_thermo(thermo_asm1)
    RWW.disconnect_sink() # disconnect from A1 to avoid replacement warning
    M1 = qs.sanunits.Mixer('M1', ins=[RWW, J2.outs[0]], isdynamic=True)
    A1.ins[1] = M1.outs[0]
    
    sys = System(path=(*bsm1_sys.units, J1, AD1, J2, M1))
    sys.set_dynamic_tracker(A1, C1, J1, AD1, J2, M1)

    return sys


if __name__ == '__main__':
    t = 50
    t_step = 0.5
    method = 'BDF'
    sys = create_system()
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0, t),
        t_eval=np.arange(0, t+t_step, t_step),
        )
    # Just to test a random state
    states = ('S_su',)
    AD1 = sys.flowsheet.unit.AD1
    fig, ax = AD1.scope.plot_time_series(states)
    
    # # Output all states, #!!! seems to have problems
    # sys.scope.export(os.path.join(results_path, f'states_{t}_{method}.xlsx'))