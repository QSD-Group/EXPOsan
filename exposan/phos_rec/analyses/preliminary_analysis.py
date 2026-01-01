# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 15:18:10 2025

@author: 19481
"""

import matplotlib.pyplot as plt
from exposan.phos_rec import create_system

result = {}

for mass_flow in [100000, 500000, 1000000]:
    sys = create_system(water_mass_flow=mass_flow)
    result[mass_flow] = sys.LCA.get_total_impacts()['GlobalWarming']/(sys.flowsheet.product.F_mass*24*365*30)

plt.scatter(result.keys(), result.values())
plt.plot(result.keys(), result.values())