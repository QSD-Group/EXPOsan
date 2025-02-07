# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from qsdsan.utils import AttrGetter
from exposan.werf.utils import plantwide_aeration_demand


__all__ = (
    'add_performance_metrics',
    )

def add_performance_metrics(model):
    
    metric = model.metric
    sys = model.system
    s = sys.flowsheet.stream
    cmps = s.RWW.components
    
    kwargs = dict(units='mg/L', element='Effluent')
    metric(getter=AttrGetter(s.SE, 'COD'), name='COD', **kwargs)
    metric(getter=AttrGetter(s.SE, 'BOD'), name='BOD', **kwargs)
    
    @metric(name='TSS', **kwargs)
    def get_TSS():
        return s.SE.get_TSS()
    
    metric(getter=AttrGetter(s.SE, 'TN'), name='TN', **kwargs)
   
    @metric(name='NH4_N', **kwargs)
    def get_NH4_N():
        return s.SE.iconc['S_NH4']
    
    metric(getter=AttrGetter(s.SE, 'TP'), name='TP', **kwargs)
    
    @metric(name='ortho_P', **kwargs)
    def get_orthoP():
        return s.SE.iconc['S_PO4']
    
    metric(getter=AttrGetter(s.SE, 'TOC'), name='TOC', **kwargs)
    
    @metric(name='CH4 production', units='kg/hr', element='Biogas')
    def get_CH4_production():
        if 'biogas' in s: return s.biogas.imass['S_ch4'] * s.biogas.components.S_ch4.i_mass
        else: return np.nan
    
    @metric(name='CH4 content', units='%', element='Biogas')
    def get_CH4_content():
        if 'biogas' in s: return s.biogas.imol['S_ch4']/sum(s.biogas.mol) * 100
        else: return np.nan
    
    @metric(name='sludge production', units='tonne/d', element='Sludge')
    def get_sludge_production():
        return sum(s.cake.mass * cmps.i_mass) * 24e-3
    
    _cached_aer = {}
    @metric(name='liquid aeration flowrate', units='m3/d', element='Aeration')
    def get_liquid_qair():
        _cached_aer.update(plantwide_aeration_demand(sys))
        qair = 0.
        for k,v in _cached_aer.items():
            if k == 'AED': continue
            else: qair += v
        return qair
    
    @metric(name='sludge aeration flowrate', units='m3/d', element='Aeration')
    def get_aed_qair():
        qair = _cached_aer.pop('AED', np.nan)
        _cached_aer.clear()
        return qair
    
    

    