# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''
import pandas as pd
from qsdsan.utils import ospath
from exposan.werf import results_path

__all__ = ('plantwide_P_mass_flows',
           'get_TP_flow',
           'get_mineralP_flow',
           'get_SPO4_flow',
           'get_bioP_flow',
           'get_orgP_flow')

get_TP_flow = lambda stream: stream.composite('P', flow=True, unit='kg/d')
get_mineralP_flow = lambda stream: stream.composite(
    'P', flow=True, unit='kg/d',
    subgroup=('X_struv', 'X_newb', 'X_ACP', 'X_AlPO4', 'X_FePO4')
    )

def get_SPO4_flow(stream):
    if 'S_PO4' in stream.components: return stream.imass['S_PO4'] * 24
    else: return stream.imass['S_IP'] * 24

def get_bioP_flow(stream):
    if 'X_H' in stream.components: 
        return stream.composite('P', flow=True, unit='kg/d', 
                                subgroup=('X_H', 'X_PAO', 'X_PP', 'X_PHA', 'X_AUT'))
    else:
        return stream.composite('P', flow=True, unit='kg/d', 
                                subgroup=(
                                    'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 
                                    'X_ac', 'X_h2', 'X_PHA', 'X_PP', 'X_PAO'
                                    ))

def get_orgP_flow(stream): 
    if 'S_F' in stream.components:
        return stream.composite('P', flow=True, unit='kg/d',
                                subgroup=('S_F', 'S_I', 'X_S', 'X_I') )# excluding biomass
    else:
        return stream.composite('P', flow=True, unit='kg/d',
                                subgroup=('S_I', 'X_li', 'X_I') )# excluding biomass        


def plantwide_P_mass_flows(system, save_as=None):
    '''
    Compiles the plant-wide phosphorus mass flow data, in kg-P/day.

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
    >>> from exposan.werf.utils import plantwide_P_mass_flows, load_state
    >>> sys = create_system('H1')
    >>> state_arr = np.load(os.path.join(data_path, f'ss/baseline_unopt/{sys.ID}.npy'))
    >>> load_state(sys, state_arr=state_arr)
    >>> sys.simulate(t_span=(0,300), method='BDF')
    >>> df = plantwide_P_mass_flows(sys)
    >>> df.head()
        stream source sink   TP  PO4-P  mineral_P  bio_P  other_organic_P
     0     RWW   None   MD  266    189          0      0               77
     1  carbon   None  ASR    0      0          0      0                0
     2     ws1     MD   PC  266   97.6       91.7      0               77
     3  reject     HD   PC 62.7   31.1       23.9   1.96             5.72
     4      PE     PC  ASR  223    128       45.9  0.779             48.3
    
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
            for getter in (get_TP_flow, get_SPO4_flow, get_mineralP_flow, get_bioP_flow, get_orgP_flow):
                entry.append(getter(s))
            data.append(entry)
    df = pd.DataFrame(data, 
                      columns=('stream', 'source', 'sink',
                               'TP', 'PO4-P', 'mineral_P', 'bio_P', 'other_organic_P'))
    if save_as:
        if not save_as.endswith('.csv'): save_as += '.csv'
        df.to_csv(ospath.join(results_path, save_as))
    else: 
        return df    