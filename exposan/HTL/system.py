#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 15:31:22 2022

@author: yalinli_cabbi
"""

import qsdsan as qs
from qsdsan import sanunits as su
from exposan.HTL._components import create_components

cmps = create_components()

inf = qs.WasteStream('inf', Sludge=10, Water=10000)



with qs.System('sys') as sys:
    U1 = su.BeltThickener('U1', ins=inf)
    U2 = su.SludgeCentrifuge('U2', ins=U1-1)
    
sys.simulate()