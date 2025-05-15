# -*- coding: utf-8 -*-

PF_CAS = dict(
    SRT = (3, 15),          # d
    F2M = (0.2, 0.4),       # kg BOD/kg MLVSS/d
    VolLoad = (0.3, 0.7),   # kg BOD/m3/d
    MLSS = (1000, 3000),    # mg/L
    tau = (4, 8)            # h
    )

mBardenpho = dict(
    SRT = (10, 20),
    MLSS = (3000, 4000), 
    tau = {
        'anae': (0.5, 1.5),
        'anox_1': (1, 3),
        'aero_1': (4, 12),
        'anox_2': (2, 4),
        'aero_2': (0.5, 1),
        },
    RAS = (50, 100),        # % of influent
    IR = (200, 400),        # internal recycle, % of influent
    )