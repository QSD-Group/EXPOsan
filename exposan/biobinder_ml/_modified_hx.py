# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 13:39:24 2025

@author: aliah
"""

from types import MethodType
import biosteam as bst

def configure_utilities():
    bst.HeatUtility.cooling_agents = [
        a for a in bst.HeatUtility.cooling_agents
        if a.ID in {'cooling_water','chilled_water','chilled_brine'}
    ]
    bst.HeatUtility.dT = max(bst.HeatUtility.dT, 10.0)
    for a in bst.HeatUtility.cooling_agents:
        if a.ID in {'chilled_water','chilled_brine'}:
            a.dT = max(a.dT, 3.0)

def limit_internal_hx(HTL, coldest_id='chilled_water', cap_margin_K=2.0, cooler_margin_K=3.0):
    agents = {a.ID: a for a in bst.HeatUtility.cooling_agents}
    a = agents.get(coldest_id, min(bst.HeatUtility.cooling_agents, key=lambda x: x.T + x.dT))
    cap = a.T + a.dT + cap_margin_K
    HTL.hx.T_lim1 = max(HTL.eff_T + cooler_margin_K, cap)

def set_eff_hx_temperature(HTL, clamp_to=None):
    HU = bst.HeatUtility
    eff = HTL.eff_hx
    original = eff.simulate_as_auxiliary_exchanger

    def wrapped(self, *args, **kwargs):
        kwargs.pop('duty', None)          # ignore duty
        self.T = HTL.eff_T                # temperature-target
        if clamp_to:
            order = [a.ID for a in HU.cooling_agents]
            cap_idx = order.index(clamp_to)
            old = HU.get_suitable_cooling_agent
            def picker(T_pinch):
                for a in HU.cooling_agents[:cap_idx+1]:
                    if T_pinch > a.T + a.dT: return a
                raise RuntimeError(f"Needs colder than {clamp_to}")
            try:
                HU.get_suitable_cooling_agent = staticmethod(picker)
                return original(*args, **kwargs)
            finally:
                HU.get_suitable_cooling_agent = old
        return original(*args, **kwargs)

    eff.simulate_as_auxiliary_exchanger = MethodType(wrapped, eff)
