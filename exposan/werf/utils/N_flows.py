# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''
import pandas as pd, numpy as np
from qsdsan.utils import ospath
from exposan.werf import results_path

__all__ = ('plantwide_N_mass_flows',
           'get_TN_flow',
           'get_orgN_flow',
           'get_NOx_flow',
           'get_NH4_flow',
           'get_SN2_flow')

get_TN_flow = lambda stream: stream.composite('N', flow=True, unit='kg/d') # already excludes dissolved N2
get_orgN_flow = lambda stream: stream.composite('N', flow=True, organic=True, unit='kg/d')

def get_NOx_flow(stream):
    if 'S_NO3' in stream.components: return stream.imass['S_NO3'] * 24
    else: return 0.

def get_NH4_flow(stream):
    if 'S_NH4' in stream.components: 
        return stream.composite('N', flow=True, subgroup=('S_NH4', 'X_struv'), unit='kg/d')
    else:
        return stream.composite('N', flow=True, subgroup=('S_IN', 'X_struv'), unit='kg/d')
        
def get_SN2_flow(stream):
    if 'S_N2' in stream.components: return stream.imass['S_N2'] * 24
    else: return 0.

def plantwide_N_mass_flows(system, save_as=None):
    '''
    Compiles the plant-wide nitrogen mass flow data, in kg-N/day.

    Parameters
    ----------
    system : :class:`qsdsan.System`
        Plant system object.
    save_as : str, optional
        File name to save the data frame as a csv. The default is None.

    Example
    -------
    >>> import os, numpy as np
    >>> from exposan.werf import create_system, data_path
    >>> from exposan.werf.utils import plantwide_N_mass_flows, load_state
    >>> sys = create_system('I3')
    >>> state_arr = np.load(os.path.join(data_path, f'ss/baseline_unopt/{sys.ID}.npy'))
    >>> load_state(sys, state_arr=state_arr)
    >>> sys.simulate(t_span=(0,300), method='BDF')
    >>> df = plantwide_N_mass_flows(sys)
    >>> df.head()
        stream source  sink       TN  NOx_N  NH4_N    org_N       N2      TKN
    0      RWW   None   ASR 1.51e+03      0    981      533      681 1.51e+03
    1      RAS     FC   ASR 1.15e+04   61.7   6.41 1.14e+04      399 1.14e+04
    2   reject     HD   ASR     58.6   2.52  0.261     55.8     16.3     56.1
    3  treated    ASR    FC 1.22e+04    156   16.2  1.2e+04 1.01e+03  1.2e+04
    4       SE     FC  None      167     92   9.57     65.1      595     74.7

    '''
    fs = system.flowsheet.stream
    data = []
    for s in fs:
        if not s.isempty():
            entry = []
            entry.append(s.ID)
            if not (s.source or s.sink): continue
            if s.source:
                if s.source not in system.path: continue
                entry.append(s.source.ID)
            else:
                entry.append(None)
            if s.sink: 
                if s.sink not in system.path: continue
                entry.append(s.sink.ID)
            else:
                entry.append(None)
            for getter in (get_TN_flow, get_NOx_flow, get_NH4_flow, get_orgN_flow, get_SN2_flow):
                entry.append(getter(s))
            data.append(entry)
    df = pd.DataFrame(data, 
                      columns=('stream', 'source', 'sink',
                               'TN', 'NOx_N', 'NH4_N', 'org_N', 'N2'))
    df['TKN'] = df.TN - df.NOx_N
    assert np.allclose(df.TN - df.NOx_N, df.NH4_N + df.org_N, atol=1e-4)
    df = df[df.TN > 0]
    add_rows = []
    for u in system.path:
        ins = sum(df.TN[df.sink == u.ID])
        outs = sum(df.TN[df.source == u.ID])
        degas = ins - outs
        if degas > 1e-6:
            add_rows.append(['vent', u.ID, None, ins-outs, 0,0,0, ins-outs, 0])
    add_rows = pd.DataFrame(add_rows, columns=df.columns)
    df = pd.concat([df, add_rows], ignore_index=True)
    if save_as:
        if not save_as.endswith('.csv'): save_as += '.csv'
        df.to_csv(ospath.join(results_path, save_as))
    else: 
        return df    