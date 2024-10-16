# -*- coding: utf-8 -*-
"""
Created on Thu May 16 10:02:37 2024

@author: joy_c

Flores-Alsina, X., Solon, K., Kazadi Mbamba, C., Tait, S., Gernaey, K. v., 
Jeppsson, U., & Batstone, D. J. (2016). Appendix A in Modelling phosphorus (P), 
sulfur (S) and iron (Fe) interactions for dynamic simulations of anaerobic 
digestion processes. Water Research, 95, 370–382. 
https://doi.org/10.1016/J.WATRES.2016.03.012

"""
# Table 1.1 [mg/L]
inf_asm2d = dict(
    S_O2=0,
    S_F=26.44,
    S_A=17.66,
    S_I=27.23,
    S_NH4=18.58,
    S_N2=5.07,
    S_NO3=0.02,
    S_PO4=4.69,
    S_ALK=78.99,
    X_I=10964.41,
    X_S=19084.76,
    X_H=9479.39,
    X_PAO=3862.20,
    X_PP=450.87,
    X_PHA=24.64,
    X_AUT=333.79,
    S_K=19.79,
    S_Mg=189.87,
    S_Na=70,
    S_Cl=1035,
    S_Ca=300,
    S_SO4=0
    )

# Table 1.3 [kg/m3]
out_adm1p = dict(
    S_su=0.018,
    S_aa=0.008,
    S_ac=0.018,
    S_IC=0.021*12,
    S_IN=0.036*14,
    S_I=0.027,
    X_ch=8.020,
    X_pr=8.481,
    X_li=11.416,
    X_I=11.946,
    S_IP=0.006*31,
    X_PHA=0.025,
    # X_PP=0.015*31,
    X_PP=0.015,
    X_PAO=3.862,
    S_Na=0.003*23,
    S_K=0.001*39,
    S_Cl=0.029*35.5,
    S_Ca=0.007*40,
    S_Mg=0.008*24.3,
    S_N2=0.0004*14
    )

# Table 1.4 [kg/m3]
inf_adm1p = dict(
    S_su=0.013,
    S_aa=0.006,
    S_fa=0.116,
    S_va=0.012,
    S_bu=0.016,
    S_pro=0.019,
    S_ac=0.055,
    S_h2=2.65e-7,
    S_ch4=0.052,
    S_IC=0.059*12,
    S_IN=0.080*14,
    S_I=0.027,
    X_ch=1.441,
    X_pr=1.513,
    X_li=2.025,
    X_I=12.345,
    S_IP=0.007*31,
    X_PHA=0.252,
    X_PP=8.05e-6,
    X_biomass=3.600,
    S_Na=0.003*23,
    S_K=0.005*39,
    S_Cl=0.029*35.5,
    S_Ca=0.001*40,
    S_Mg=0.001*24.3,
    X_Ca3PO43=0.002, # kmol/m3
    X_MgNH4PO4=0.011, # kmol/m3
    S_N2=0.0004*14
    )

# Table 1.5 [mg/L]
out_asm2d = dict(
    S_F=134.43,
    S_A=353.82,
    S_I=27.23,
    S_NH4=1291.68,
    S_PO4=298.09,
    S_ALK=885.27, # S_IC
    X_I=12704.93,
    X_S=8218.94,
    S_K=208.84,
    S_Mg=28.29,
    S_Na=70,
    S_Cl=1035,
    S_Ca=20.45,
    X_Ca3PO43=722.17,
    X_MgNH4PO4=1578.52
    )