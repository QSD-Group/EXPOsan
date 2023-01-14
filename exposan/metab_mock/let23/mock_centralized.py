# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from qsdsan.utils import auom
import chaospy.distributions as shape, numpy as np
import chaospy

N_blower = 7
lhp = 3000
hhp = 5500


def E_per_rCOD(n_lhp, Q, BOD_in, rBOD_pri, BOD2COD_pre2nd, BOD_eff, BOD2COD_eff):
    avg_power = (n_lhp * lhp + (N_blower - n_lhp) * hhp) \
        * auom('hp').conversion_factor('kW')
    Q = Q * auom('mgd').conversion_factor('L/hr')
    pre2_BOD = BOD_in * (1-rBOD_pri)
    pre2_COD = pre2_BOD / BOD2COD_pre2nd
    eff_COD = BOD_eff / BOD2COD_eff
    rCOD = (pre2_COD - eff_COD) * auom('mg/L').conversion_factor('kg/L')
    return avg_power / (rCOD * Q)  # kWh/kg COD removal


def E_model():
    RVs = []
    RVs.append(shape.Binomial(N_blower, 0.5))  # distribution of the number of low-hp blowers
    RVs.append(shape.Triangle(170, 185, 250))  # distribution of influent Q [MGD]
    RVs.append(shape.Triangle(110, 190, 350))  # untreaded municipal WW BOD [mg/L]
    RVs.append(shape.Uniform(0.5, 0.6))        # primary treatment BOD removal rate
    RVs.append(shape.Uniform(0.4, 0.6))        # BOD/COD after primary treatment
    RVs.append(shape.Triangle(0, 5, 5))        # Effluent BOD 
    RVs.append(shape.Uniform(0.1, 0.3))        # BOD/COD after secondary treatment
    params = chaospy.J(*RVs)
    samples = params.sample(size=400, rule='latin_hypercube').T
    Es = np.array([E_per_rCOD(*row) for row in samples])
    return Es