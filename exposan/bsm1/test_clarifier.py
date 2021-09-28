# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 20:44:27 2021

@author: joy_c
"""
import pandas as pd
import numpy as np
from qsdsan.utils import load_data
tss = load_data("C:\\Users\\joy_c\\Dropbox\\PhD\\Research\\QSD\\codes_developing\\EXPOsan\\exposan\\bsm1\\TSS_gpsx_5d.csv")
cols = ['L%s'%i for i in range(1,11)]
dTSS = pd.DataFrame(columns=cols)

def _settling_flux(X, v_max, v_max_practical, X_min, rh, rp):
    X_star = max(X-X_min, 0)
    v = min(v_max_practical, v_max*(np.exp(-rh*X_star) - np.exp(-rp*X_star)))
    return X*max(v, 0)

X_t=3000
n=10
jf=3
hj=0.4
A=1500
vmax=474
vmaxp=250
rh=5.76e-4 
rp=2.86e-3
fns=2.28e-3

def dydt(Q_in, X_in, X, Q_s):
    Q_e = Q_in - Q_s
    X_min = X_in * fns
    Q_jout = np.array([Q_e if j < jf else Q_in if j == jf else Q_s for j in range(n)])
    flow_out = X*Q_jout
    flow_in = np.array([Q_e*X[j+1] if j < jf else Q_in*X_in if j == jf else Q_s*X[j-1] for j in range(n)])
    VX = [_settling_flux(xj, vmax, vmaxp, X_min, rh, rp) for xj in X]
    J = np.array([min(VX[j], VX[j+1]) for j in range(n-1)])
    # J = np.array([VX[j] if X[j+1] <= X_t and j < jf else min(VX[j], VX[j+1]) for j in range(n-1)])
    settle_out = np.append(J, 0)
    settle_in = np.insert(J, 0, 0)
    TSS_dot = ((flow_in - flow_out)/A + settle_in - settle_out)/hj        # (n,)
    return TSS_dot

X_in = 2350
Q_in = 36892
Q_s = 18831
X = np.array([10, 20, 40, 70, 200, 300, 350, 350, 2000, 4000])

for t, row in tss.iterrows():
    dtss = dydt(Q_in=36892, X_in=row[-1], X=row[:10], Q_s=18831)
    dTSS = dTSS.append(dict(zip(cols, dtss)), ignore_index=True)

dTSS.to_csv('dTSS_5d.csv')
