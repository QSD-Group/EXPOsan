# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 15:53:18 2023

@author: joy_c
"""
import numpy as np, pandas as pd
import qsdsan.processes as pc
from qsdsan.utils import auom, ospath, time_printer
from scipy.integrate import solve_ivp
from exposan.metab_mock import results_path, biomass_IDs, create_systems

Q = 5 # m3/d

cmps = pc.create_adm1_cmps()
bm_idx = cmps.indices(biomass_IDs)
S_idx = cmps.indices([i for i in cmps.IDs if i.startswith('S_')])
adm1 = pc.ADM1()
stoi_bk = adm1.stoichio_eval()
stoi_en = stoi_bk[:-3]  # no liquid-gas transfer
Rho_bk = adm1.rate_function

#%% functions & constants

from qsdsan.processes._adm1 import (
    mass2mol_conversion,
    T_correction_factor,
    acid_base_rxn,
    substr_inhibit,
    Hill_inhibit,
    non_compet_inhibit
    )
from scipy.optimize import brenth

rhos = np.zeros(19) # 19 kinetic processes, no gas-liquid transfer
Cs = np.empty(19)
def Rho_en(state_arr, params=Rho_bk._params, T_op=298.15):
    ks = params['rate_constants']
    Ks = params['half_sat_coeffs']
    cmps = params['components']
    pH_ULs = params['pH_ULs']
    pH_LLs = params['pH_LLs']
    KS_IN = params['KS_IN']
    KI_nh3 = params['KI_nh3']
    KIs_h2 = params['KIs_h2']
    Kab = params['Ka_base']
    Ka_dH = params['Ka_dH']
    T_base = params['T_base']

    Cs[:8] = state_arr[12:20]
    Cs[8:12] = state_arr[19:23]
    Cs[12:] = state_arr[16:23]

    substrates = state_arr[:8]

    S_va, S_bu, S_h2, S_IN = state_arr[[3,4,7,10]]
    unit_conversion = mass2mol_conversion(cmps)
    cmps_in_M = state_arr[:27] * unit_conversion
    weak_acids = cmps_in_M[[24, 25, 10, 9, 6, 5, 4, 3]]

    Kas = Kab * T_correction_factor(T_base, T_op, Ka_dH)

    rhos[:] = ks * Cs
    Monod = substr_inhibit(substrates, Ks)
    rhos[4:12] *= Monod
    if S_va > 0: rhos[7] *= 1/(1+S_bu/S_va)
    if S_bu > 0: rhos[8] *= 1/(1+S_va/S_bu)

    h = brenth(acid_base_rxn, 1e-14, 1.0,
            args=(weak_acids, Kas),
            xtol=1e-12, maxiter=100)
        
    nh3 = Kas[1] * weak_acids[2] / (Kas[1] + h)
    
    Iph = Hill_inhibit(h, pH_ULs, pH_LLs)
    Iin = substr_inhibit(S_IN, KS_IN)
    Ih2 = non_compet_inhibit(S_h2, KIs_h2)
    Inh3 = non_compet_inhibit(nh3, KI_nh3)
    rhos[4:12] *= Iph * Iin
    rhos[6:10] *= Ih2
    rhos[10] *= Inh3
    return rhos

n_cmps = len(cmps)
Diff = np.zeros(n_cmps)
Diff[:12] = [4.56e-6, 8.62e-6, 5.33e-6, 5.00e-6, 5.04e-6, 6.00e-6, 
             6.48e-6, 4.02e-4, 1.36e-4, 1.71e-4, 1.52e-4, 6.00e-6]

Diff[-3:-1] = 1.17e-4
# f_diff = 0.75       # bead-to-water diffusivity fraction
Cs_in = np.array([
    3.000e+03, 6.000e+02, 4.000e+02, 4.000e+02, 4.000e+02, 4.000e+02,
    4.000e+02, 5.000e-06, 5.000e-03, 4.804e+02, 1.401e+02, 2.000e+01,
    1.000e+02, 3.000e+02, 5.000e+02, 2.500e+02, 0.000e+00, 1.000e+00,
    1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 1.000e+00, 2.500e+01,
    4.000e+01, 2.000e+01, 9.900e+05
    ])*1e-3 # mg/L to kg/m3

_R = 8.3145e-2 # Universal gas constant, [bar/M/K]
_k_p = 5.0e4
_P_atm = 1.013
Pa_to_bar = auom('Pa').conversion_factor('bar')
p_vapor = lambda T: cmps.H2O.Psat(T)*Pa_to_bar
gas_mass2mol_conversion = (cmps.i_mass / cmps.chem_MW)[cmps.indices(('S_h2', 'S_ch4', 'S_IC'))]
# K_tss = 25/2        # TSS half saturation coefficient

def f_qgas(S_gas, T):
    p_gas = S_gas * _R * T
    P = sum(p_gas) + p_vapor(T) 
    q_gas = max(0, _k_p * (P - _P_atm))
    return q_gas
    
def compile_ode(r_beads=5e-3, l_bl=1e-5, f_void=0.1, 
                n_dz=10, f_diff=0.75, K_tss=11):
    dz = r_beads / n_dz
    zs = np.linspace(dz, r_beads, n_dz)
    dV = 4/3*np.pi*(zs)**3
    dV[1:] -= dV[:-1]
    V_bead = (4/3*np.pi*r_beads**3)
    D = Diff * f_diff   # Diffusivity in beads
    k = Diff/l_bl       # mass transfer coeffient through liquid boundary layer
    D_ov_dz2 = D[S_idx]/(dz**2)     # (n_soluble,)
    D_ov_dz = D[S_idx]/dz
    _1_ov_z = 1/zs                  # (n_dz,)
        
    def dydt(t, y, HRT, detach=False, retain=False):
        Q, T = y[-2:]
        V_liq = HRT * Q
        V_bed = V_liq / f_void
        V_beads = V_bed - V_liq
        V_gas = V_bed * 0.1
        A_beads = 3 * V_beads / r_beads # m2, total bead surface area
    
        S_gas = y[-5:-2]
        Cs_bk = y[-(n_cmps+5):-5]        # bulk liquid concentrations
        Cs_en = y[:n_dz*n_cmps].reshape((n_dz, n_cmps))  # each row is one control volume
        
        # Transformation
        rhos_en = np.apply_along_axis(Rho_en, 1, Cs_en, T_op=T)
        Rs_en = rhos_en @ stoi_en       # n_dz * n_cmps
        rhos_bk = Rho_bk(y[-(n_cmps+5):])
        Rs_bk = np.dot(stoi_bk.T, rhos_bk) # n_cmps+5
        q_gas = f_qgas(S_gas, T)
        gas_transfer = - q_gas*S_gas/V_gas + rhos_bk[-3:] * V_liq/V_gas * gas_mass2mol_conversion
        
        # Detachment -- particulates
        if detach:
            tss = np.sum(Cs_en * (cmps.x*cmps.i_mass), axis=1)
            x_net_growth = np.sum(Rs_en * cmps.x, axis=1)/np.sum(Cs_en * cmps.x, axis=1) # d^(-1), equivalent to k_de
            u_de = 1/(1+np.exp(K_tss-tss)) * np.maximum(x_net_growth, 0)
            de_en = np.diag(u_de) @ (Cs_en * cmps.x)
            tot_de = np.sum(np.diag(dV) @ de_en, axis=0) / V_bead  # detachment per unit volume of beads
        else:
            de_en = tot_de = 0
            
        # Mass transfer (centered differences) -- MATLAB/Simulink
        # S_en = Cs_en[:, S_idx]
        # dSdz2 = np.zeros_like(Cs_en)
        # dSdz2[0, S_idx] = (1+dz/(2*zs[0]))**2 * (S_en[1] - S_en[0])/dz
        # dSdz2[1:-1, S_idx] = (np.diag((1+dz/(2*zs[1:-1]))**2) @ (S_en[2:]-S_en[1:-1]) \
        #     - np.diag((1-dz/(2*zs[1:-1]))**2) @ (S_en[1:-1]-S_en[:-2]))/(dz**2)
        # C_lf = Cs_en[-1].copy()
        # S_en[-1] = C_lf[S_idx] = (D/dz*Cs_en[-2] + k*Cs_bk)[S_idx]/(D/dz+k)[S_idx]  # Only solubles based on outer b.c.
        # dSdz2[-1, S_idx] = ((1+dz/(2*zs[-1]))**2 * (S_en[-1]-S_en[-2]) \
        #     -(1-dz/(2*zs[-1]))**2 * (S_en[-2]-S_en[-3]))/(dz**2)
                
        # dCdz = Cs_en * 0.
        # d2Cdz2 = Cs_en * 0.
        # dCdz[1:-1] = (Cs_en[2:] - Cs_en[:-2])/(2*dz)
        # d2Cdz2[1:-1] = (Cs_en[2:] + Cs_en[:-2] - 2*Cs_en[1:-1])/(dz**2)
        # d2Cdz2[0,:] = (Cs_en[1] - Cs_en[0])*(1+dz/(2*zs[0]))**2/dz  # inner boundary condition, forward difference
        # C_lf = Cs_en[-1].copy()
        # C_lf[S_idx] = (D*Cs_en[-2] + k*Cs_bk*dz)[S_idx]/(D+k*dz)[S_idx]  # Only solubles based on outer b.c.
        # d2Cdz2[-1,:] = ((1+dz/(2*zs[-1]))**2*(C_lf-Cs_en[-2]) - (1-dz/(2*zs[-1]))**2*(Cs_en[-2]-Cs_en[-3]))/(dz**2)  # backward difference
        
        #!!! Mass transfer (centered differences) -- MOL solubles only
        C_lf = Cs_en[-1]
        J_lf = k*(Cs_bk - C_lf)
        S_en = Cs_en[:, S_idx]
        M_transfer = np.zeros_like(Cs_en)
        M_transfer[1:-1, S_idx] = D_ov_dz2 * (S_en[2:] - 2*S_en[1:-1] + S_en[:-2])\
            + D_ov_dz * (np.diag(_1_ov_z[1:-1]) @ (S_en[2:] - S_en[:-2]))
        M_transfer[0, S_idx] = 2 * D_ov_dz2 * (S_en[1] - S_en[0])
        M_transfer[-1, S_idx] = 2 * D_ov_dz2 * (S_en[-2] - S_en[-1])\
            + 2 * (1/dz + _1_ov_z[-1]) * J_lf[S_idx]
        
        # Mass balance
        # dCdt_en = D*dSdz2 + Rs_en - de_en
        # dCdt_en = D*(d2Cdz2 + np.diag(2/zs) @ dCdz) + Rs_en - de_en
        # dCdt_en = dJ_dz + Rs_en - de_en
        dCdt_en = M_transfer + Rs_en - de_en
        if retain:
            Cs_eff = Cs_bk.copy()
            Cs_eff[bm_idx] = 0.
            # dCdt_bk = (Cs_in-Cs_eff)/HRT - A_beads/V_liq*k*(Cs_bk-C_lf) + Rs_bk + V_beads*tot_de/V_liq
            dCdt_bk = (Cs_in-Cs_eff)/HRT - A_beads/V_liq*J_lf + Rs_bk + V_beads*tot_de/V_liq
        else:
            # dCdt_bk = (Cs_in-Cs_bk)/HRT - A_beads/V_liq*k*(Cs_bk-C_lf) + Rs_bk + V_beads*tot_de/V_liq
            dCdt_bk = (Cs_in-Cs_bk)/HRT - A_beads/V_liq*J_lf + Rs_bk + V_beads*tot_de/V_liq

        dy = y*0.
        dy[:n_dz*n_cmps] = dCdt_en.flatten()
        dy[-(n_cmps+5):-5] = dCdt_bk
        dy[-5:-2] = gas_transfer
        
        return dy
    
    return dydt

#%%
# HRT = 1
C0 = np.array([
    1.204e-02, 5.323e-03, 9.959e-02, 1.084e-02, 1.411e-02, 1.664e-02,
    4.592e-02, 2.409e-07, 7.665e-02, 5.693e-01, 1.830e-01, 2e-02,
    2.424e-01, 2.948e-02, 4.766e-02, 2.603e-02, 4.708e+00, 1.239e+00,
    4.838e-01, 1.423e+00, 8.978e-01, 2.959e+00, 1.467e+00, 4.924e-02,
    4.000e-02, 2.000e-02, 9.900e+02
    ])

y0_bulk = np.array([
    1.204e-02, 5.323e-03, 9.959e-02, 1.084e-02, 1.411e-02, 1.664e-02,
    4.592e-02, 2.409e-07, 7.665e-02, 5.693e-01, 1.830e-01, 3.212e-02,
    2.424e-01, 2.948e-02, 4.766e-02, 2.603e-02, 0, 0,
    0, 0, 0, 0, 0, 0,   # no biomass in bulk liquid initially
    4.000e-02, 2.000e-02, 9.900e+02, 3.922e-07, 2.228e-02, 1.732e-02,
    Q, 298.15 
    ])

y0_from_bulk = lambda n_dz: np.append(np.tile(C0, n_dz), y0_bulk)

#%% steady-state spatial profile
def y0_even(n_dz, TSS0):
    C0_even = C0.copy()
    # C0_even = Cs_in.copy()    
    C0_even[bm_idx] = TSS0/0.777/7
    return np.append(np.tile(C0_even, n_dz), y0_bulk)

def spatial_profiling(HRTs=[1, 0.5, 10/24, 8/24], TSS0=5, n_dz=10, 
                      detach=False, save_to='', **ode_kwargs):
    y0 = y0_even(n_dz, TSS0)
    dfs = {}
    dydt = compile_ode(n_dz=n_dz, **ode_kwargs)
    for tau in HRTs:
        print(f'HRT = {tau:.2f} d')
        sol = solve_ivp(dydt, t_span=(0, 400), y0=y0, method='BDF', args=(tau, detach))
        y_ss = sol.y.T[-1]
        C_ss = y_ss[:n_cmps*(n_dz+1)].reshape((n_dz+1, n_cmps))
        Xbio_ss = C_ss[:, bm_idx]
        df_c = pd.DataFrame(C_ss, columns=cmps.IDs, index=[*range(n_dz), 'bulk'])
        df_c['biomass_COD'] = np.sum(Xbio_ss, axis=1)
        dfs[f'{tau:.2f}'] = df_c
    
    file_name = save_to or f'TSS0_{TSS0}_{detach}.xlsx'
    path = ospath.join(results_path, file_name)
    with pd.ExcelWriter(path) as writer:
        for k, df in dfs.items():
            df.to_excel(writer, sheet_name=k)
    return dfs

#%% suspended vs. attached growth

@time_printer
def suspended_vs_attached(HRT, dydt, n_dz, r_beads, sys, f_void=0.1,
                          y0=None, detach=False, retain=False):
    print('\n')
    print('='*12)
    print(f'HRT: {HRT:.3f} d')
    u, = sys.units
    u.V_liq = HRT * 5
    u.V_gas = HRT * 0.5
    inf, = sys.feeds
    bg, eff = sys.products
    sys.simulate(state_reset_hook='reset_cache', t_span=(0, 400), method='BDF')
    sus_rcod = 1-eff.COD/inf.COD
    i = 0
    while sus_rcod < 0.7 and i <= 5:
        i += 1
        print('add VSS...')
        u._concs[bm_idx] *= 2
        sys.simulate(state_reset_hook='reset_cache', t_span=(0, 400), method='BDF')
        sus_rcod = 1-eff.COD/inf.COD
    sus_bm = sum(u._state[bm_idx] * cmps.i_mass[bm_idx])
    if y0 is None:
        y0_susp = np.append(u._state, 298.15)
        if not retain: y0_susp[bm_idx] = 0.
        y0 = np.append(np.tile(u._state[:n_cmps], n_dz), y0_susp)
    try: y = converge(dydt, y0, args=(HRT, detach, retain))
    except:
        try:
            y0 = y0_even(n_dz, TSS0=min(2.5/HRT, 22))
            y = converge(dydt, y0, args=(HRT, detach, retain))
        except:
            return sus_rcod, np.nan, sus_bm, np.nan, None        
    bk_Cs = y[-(n_cmps+5):-5].copy()
    if retain: bk_Cs[bm_idx] = 0.
    att_COD = sum(bk_Cs * cmps.i_COD * (1-cmps.g)) * 1e3
    att_rcod = 1-att_COD/inf.COD
    bk_bm, en_bm = biomass_CODs(y, n_dz, r_beads)
    att_bm = ((1-f_void)*en_bm + f_void*bk_bm) * 0.777
    print(bk_bm*0.777, en_bm*0.777)
    return [sus_rcod, att_rcod, sus_bm, att_bm, y]

def HRT_suspended_vs_encap(HRTs=[1, 0.5, 10/24, 8/24, 4/24, 2/24, 1/24], 
                           frac_retain=1.0, r_beads=1e-5, l_bl=1e-6, n_dz=8,
                           **ode_kwargs):
    record = []
    # rows = []
    dydt = compile_ode(r_beads=r_beads, l_bl=l_bl, n_dz=n_dz, **ode_kwargs)
    # r_beads = ode_kwargs.pop('r_beads', 5e-3)
    # n_dz = ode_kwargs.pop('n_dz', 10)
    
    sys, = create_systems(which='C')
    u, = sys.units
    u._f_retain = (u._f_retain > 0) * frac_retain
    u.encapsulate_concentration = 1e5
    y = y0_even(n_dz, TSS0=2.5)
    # y = None
    for tau in HRTs:
        out = suspended_vs_attached(tau, dydt, n_dz, r_beads, sys, y,
                                    detach=True, retain=True)
        y = out.pop(-1)
        record.append(out)
        # rows.append(tau)
    
    # df = pd.DataFrame(record, index=rows, 
    df = pd.DataFrame(record, index=HRTs, 
                      columns=pd.MultiIndex.from_product(
                          [['rCOD', 'biomass [g TSS/L]'], ['Suspended', 'Encapsulated']]
                          ))
    df.index.name = 'HRT [d]'
    df.to_excel(ospath.join(results_path, 'suspended_vs_encap_retain.xlsx'))

    return df

#%% encapsulated biomass, HRT
def biomass_CODs(y, n_dz, r_beads):
    en_bm = np.sum(y[:n_dz*n_cmps].reshape((n_dz, n_cmps))[:,bm_idx], axis=1)
    bk_bm = np.sum(y[n_dz*n_cmps: ((n_dz+1)*n_cmps)][bm_idx])
    dz = r_beads / n_dz
    zs = np.linspace(dz, r_beads, n_dz)
    dV = 4/3*np.pi*(zs)**3
    V_bead = dV[-1]
    dV[1:] -= dV[:-1]
    C_en_avg = np.dot(en_bm, dV)/V_bead
    return bk_bm, C_en_avg

def SRT(y, n_dz, HRT, r_beads, f_void):
    '''SRT given there's no retention of biomass in bulk liquid'''
    bk_bm, C_en_avg = biomass_CODs(y, n_dz, r_beads)
    return (1 + (1-f_void)/f_void*C_en_avg/bk_bm) * HRT

def converge(dydt, y, args, threshold=1e-4, once=False):
    dy_max = 1.
    while dy_max > threshold:
        sol = solve_ivp(dydt, t_span=(0, 400), y0=y, method='BDF', args=args)
        y = sol.y.T[-1]
        rcod = 1 - sum(y[-(n_cmps+5):-5] * cmps.i_COD)/6.801
        dy = dydt(0, y, *args)
        dy_max = np.max(np.abs(dy))
        print(f'rCOD {rcod:.3f}, dy_max {dy_max:.2e}')
        if once: break
    return y

@time_printer
def encap(TSS_init, HRT, dydt, n_dz, f_void=0.1, r_beads=5e-3, 
          detach=False, retain=False, y0=None, threshold=1e-4):
    msg = f'r_beads {r_beads*1e3:.2f}mm, HRT {HRT:.2f}d'
    print(f'\n{msg}\n{"="*len(msg)}')
    if y0 is None: y = y0_even(n_dz, TSS_init)
    else: y = y0
    try: 
        y = converge(dydt, y, args=(HRT, detach, retain), threshold=threshold)
    except: 
        try: 
            print('Try bulk y0...')
            # y = y0_even(n_dz, TSS_init*2)
            y = y0_from_bulk(n_dz)
            y = converge(dydt, y, args=(HRT, detach, retain), threshold=threshold)
        except: return [np.nan]*4
    bk_Cs = y[-(n_cmps+5):-5]
    en_Cs = y[-(2*n_cmps+5):-(n_cmps+5)]
    att_COD = sum(bk_Cs * cmps.i_COD * (1-cmps.g)) * 1e3
    att_rcod = 1-att_COD/6801
    att_bm = sum(en_Cs[bm_idx] * cmps.i_mass[bm_idx])
    srt = SRT(y, n_dz, HRT, r_beads, f_void)
    return att_rcod, att_bm, srt, y

def HRT_init_TSS(TSS=[1, 2, 5, 10, 30], HRTs=[1, 0.5, 10/24, 8/24, 4/24, 2/24, 1/24],
                 detach=False, retain=False, **ode_kwargs):
    rcod = []
    ss_bm = []
    dydt = compile_ode(detach=detach, **ode_kwargs)
    r_beads = ode_kwargs.pop('r_beads', 5e-3)
    f_void = ode_kwargs.pop('f_void', 0.1)
    n_dz = ode_kwargs.pop('n_dz', 10)
    for tss in TSS:
        l_rcod = []
        l_bm = []
        for tau in HRTs:
            out = encap(tss, tau, dydt, n_dz, f_void=f_void, r_beads=r_beads, 
                        detach=detach, retain=retain)
            l_rcod.append(out[0])
            l_bm.append(out[1])
        rcod.append(l_rcod)
        ss_bm.append(l_bm)
    
    df_rcod = pd.DataFrame(rcod).T
    df_bm = pd.DataFrame(ss_bm).T
    
    for i in (df_rcod, df_bm):
        i.columns = pd.MultiIndex.from_product([['initial TSS [g/L]'], TSS])
        i.index = HRTs
        i.index.name = 'HRT [d]'
    
    with pd.ExcelWriter(ospath.join(results_path, 'HRT_vs_encapTSS_detach.xlsx')) as writer:
        df_rcod.to_excel(writer, sheet_name='rCOD')
        df_bm.to_excel(writer, sheet_name='Biomass TSS')
        
#%% bead size

def bead_size_HRT(HRTs=[1, 0.5, 10/24, 8/24, 4/24, 2/24, 1/24], 
                  bead_size=[5e-3, 1e-3, 5e-4, 1e-4, 1e-5], 
                  voidage=[0.39, 0.39, 0.5, 0.6, 0.7],
                  detach=False, retain=False, **ode_kwargs):
    # ode_kwargs.pop('r_beads', None)
    rcod = []
    ss_bm = []
    srt = []
    # n_dz = ode_kwargs.pop('n_dz', 8)
    # f_void = ode_kwargs.pop('f_void', 0.1)
    for r_beads, f_void in zip(bead_size, voidage):
        l_bl = min(1e-5, r_beads/10)
        if r_beads > 1e-3: n_dz = 10
        elif r_beads < 1e-4: n_dz = 5
        else: n_dz = 8
        dydt = compile_ode(r_beads=r_beads, l_bl=l_bl, f_void=f_void, n_dz=n_dz, **ode_kwargs)
        l_rcod = []
        l_bm = []
        l_srt = []
        y0 = y0_from_bulk(n_dz)
        for tau in HRTs:
            tss = min(2.5/tau, 22)
            z1, z2, z3, y0 = encap(tss, tau, dydt, n_dz, f_void, 
                                   r_beads, detach, retain, y0, 5e-4)
            l_rcod.append(z1)
            l_bm.append(z2)
            l_srt.append(z3)
        rcod.append(l_rcod)
        ss_bm.append(l_bm)
        srt.append(l_srt)
    
    sys, = create_systems(which='C')
    u, = sys.units
    inf, = sys.feeds
    bg, eff = sys.products
    u.encapsulate_concentration = 1e5
    l_rcod = []
    for tau, tau_s in zip(HRTs, srt[-1]):
        u._f_retain[bm_idx] = 1-tau/tau_s
        u.V_liq = tau*5
        u.V_gas = tau*0.5
        sys.simulate(state_reset_hook='reset_cache',
                     t_span=(0, 400), method='BDF')
        l_rcod.append(1-eff.COD/inf.COD)
        
    df_rcod = pd.DataFrame(rcod).T
    df_bm = pd.DataFrame(ss_bm).T
    df_srt = pd.DataFrame(srt).T
    
    for i in (df_rcod, df_bm, df_srt):
        i.columns = pd.MultiIndex.from_product([['bead size [m]'], bead_size])
        i.index = HRTs
        i.index.name = 'HRT [d]'
    df_rcod[('Suspended','')] = l_rcod
    with pd.ExcelWriter(ospath.join(results_path, 
                                    'HRT_vs_rbeads_vs_suspended.xlsx')) as writer:
        df_rcod.to_excel(writer, sheet_name='rCOD')
        df_bm.to_excel(writer, sheet_name='Biomass TSS')
        df_srt.to_excel(writer, sheet_name='SRT')

    return df_rcod, df_bm, df_srt

#%%
def voidage_vs_hrt(voidage = np.linspace(0.39, 0.9, 10),
                   HRTs=np.array([1, 0.5, 10/24, 8/24, 4/24, 2/24, 1/24]),
                   r_beads=1e-5, l_bl=1e-6, n_dz=5, f_diff=0.75, K_tss=11):
    sys, = create_systems(which='C')
    u, = sys.units
    inf, = sys.feeds
    bg, eff = sys.products
    u.encapsulate_concentration = 1e5
    df_rcod = pd.DataFrame(index=np.int32(HRTs*24), 
                           columns=pd.MultiIndex.from_product(
                               [np.round(voidage,3), ['encap','suspend']]
                               ))
    df_rcod.index.name = 'HRT [h]'
    df_rcod.columns.names = ['voidage', 'biomass']
    df_bm = df_rcod.copy()
    df_srt = df_rcod.copy()
    # df_srt = pd.DataFrame(index=np.int32(HRTs*24), columns=np.round(voidage,3))
    # df_srt.index.name = df_rcod.index.name
    # df_srt.columns.name = 'voidage'
    for f_void in voidage:
        dydt = compile_ode(r_beads=r_beads, l_bl=l_bl, f_void=f_void, n_dz=n_dz,
                           f_diff=f_diff, K_tss=K_tss)
        y0 = y0_from_bulk(n_dz)
        for tau in HRTs:
            itau = int(tau*24)
            ivoid = round(f_void,3)
            msg = f'voidage {ivoid}, HRT {itau}h'
            print(f'\n{msg}\n{"="*len(msg)}')
            y0 = converge(dydt, y0, (tau, True, False), 
                          threshold=1e-4, once=False)           
            bk_Cs = y0[-(n_cmps+5):-5].copy()
            att_COD = sum(bk_Cs * cmps.i_COD * (1-cmps.g)) * 1e3
            df_rcod.loc[itau, (ivoid, 'encap')] = 1-att_COD/inf.COD
            bk_bm, en_bm = biomass_CODs(y0, n_dz, r_beads)
            df_bm.loc[itau, (ivoid, 'encap')] = ((1-f_void)*en_bm + f_void*bk_bm) * 0.777
            df_srt.loc[itau, (ivoid, 'encap')] = srt = SRT(y0, n_dz, tau, r_beads, f_void)
            
            u._f_retain[bm_idx] = 1-tau/srt
            u.V_liq = tau*5
            u.V_gas = tau*0.5
            sys.simulate(state_reset_hook='reset_cache',
                         t_span=(0, 400), method='BDF')
            df_rcod.loc[itau, (ivoid, 'suspend')] = 1-eff.COD/inf.COD
            df_bm.loc[itau, (ivoid, 'suspend')] = sp_bm = sum(u._state[bm_idx]) * 0.777
            df_srt.loc[itau, (ivoid, 'suspend')] = tau * sp_bm/(sum(eff.state[bm_idx])/1e3*0.777)
            
    with pd.ExcelWriter(ospath.join(results_path, 'void_vs_HRT.xlsx')) as writer:
        df_rcod.to_excel(writer, sheet_name='rCOD')
        df_bm.to_excel(writer, sheet_name='biomass TSS')
        df_srt.to_excel(writer, sheet_name='SRT')
    
    return df_rcod, df_bm, df_srt

#%%
@time_printer
def run(tau=0.5, TSS0=5, r_beads=5e-3, l_bl=1e-5, f_void=0.1, n_dz=10, 
        f_diff=0.75, K_tss=11, detach=False, retain=False, y0=None):
    if y0 is None: 
        # y0 = y0_even(n_dz, TSS0)
        y0 = y0_from_bulk(n_dz)
    dydt = compile_ode(r_beads, l_bl, f_void, n_dz, f_diff, K_tss)
    print(f'HRT = {tau:.2f} d')
    y_ss = converge(dydt, y0, args=(tau, detach, retain))
    C_ss = y_ss[:n_cmps*(n_dz+1)].reshape((n_dz+1, n_cmps))
    Xbio_ss = C_ss[:, bm_idx]
    df_c = pd.DataFrame(C_ss, columns=cmps.IDs, index=[*range(n_dz), 'bulk'])
    df_c['biomass_COD'] = np.sum(Xbio_ss, axis=1)
    bk = C_ss[-1].copy()
    if retain: bk[bm_idx] = 0.
    eff_cod = np.sum(bk * cmps.i_COD * (1-cmps.g))
    rcod = 1-eff_cod/6.801
    srt = SRT(y_ss, n_dz=n_dz, HRT=tau, r_beads=r_beads, f_void=f_void)
    return df_c, y_ss, rcod, srt
    

if __name__ == '__main__':
    # dfs = bead_size_HRT(detach=True)
    # dfs = bead_size_HRT(bead_size=[1e-4, 1e-5], detach=True, retain=False, 
    #                     n_dz=8, f_void=0.95)
    # de_rcod, de_bm = bead_size_HRT(detach=True, n_dz=20)

    # for detach in (True,):
    #     print(f'detach: {detach}')
    #     for r_beads in (5e-3, 1e-3, 5e-4, 1e-4, 1e-5):
    #         print(f'r_beads: {r_beads}\n{"="*15}')
    #         dfs = spatial_profiling(
    #             HRTs=[1, 0.5, 1/3], TSS0=6, 
    #             r_beads=r_beads,  detach=detach,
    #             save_to=f'spatial_r{r_beads}_{detach}_new.xlsx'
    #             )
    # df, y, rcod, srt = run(tau=1, TSS0=2.5, r_beads=1e-5, l_bl=1e-6, n_dz=5, 
    #                        f_void=0.5, detach=True, retain=False)
    # df = run(tau=1/6, TSS0=6, detach=True)
   
    # df = HRT_suspended_vs_encap()
    rcod, tss, srt = voidage_vs_hrt()
