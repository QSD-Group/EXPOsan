#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 22:18:49 2024

@author: saumitrarai
"""

import qsdsan as qs
from qsdsan import processes as pc
import numpy as np

adm1_state_array = np.array([9.973e+00, 4.397e+00, 10.000e+00, 10.000e+00, 
                             10.000e+00, 10.000e+00, 4.723e+01, 10.000e+00, 
                             0.000e+00, 1.160e+03, 2.359e+01, 4.303e+02,
                             1.496e+01, 1.033e+04, 2.666e+04, 2.251e+04, 
                             10.000e+00, 10.000e+00, 10.000e+00, 10.000e+00, 
                             10.000e+00, 10.000e+00, 10.000e+00, 1.240e+04,
                             6.170e+01, 1.471e+01, 2.621e+01, 0.000e+00, 
                             0.000e+00, 10.000e+00, 10.000e+00, 0.000e+00, 
                             8.341e+01, 9.342e+05, 0.000e+00, 0.000e+00, 
                             4.254e-02, 1.397e+02, 3.081e+02])

adm1 = pc.create_adm1_p_extension_cmps()
qs.get_thermo()
ADM1 = pc.ADM1_p_extension()
rhos = pc.rhos_adm1_p_extension(adm1_state_array, ADM1.rate_function._params)