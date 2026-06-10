# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

This module builds the Benchmark Simulation Model No. 2 (BSM2) and its
phosphorus-extended variant through a single ``create_system`` entry point:

    - ``kind='bsm2'``  : the classic BSM2 using ASM1 + ADM1 [1]_.
    - ``kind='bsm2p'`` : a phosphorus-removal variant using mASM2d + ADM1p.

``create_subsys`` builds ``bsm1p``, the activated-sludge sub-train of the
``bsm2p`` plant. It is a diagnostic system for isolating the bioreactor train,
not a published benchmark, and is intentionally not part of ``create_system``.

References
----------
.. [1] Alex, J.; Benedetti, L.; Copp, J. B.; Gernaey, K. V.; Jeppsson, U.;
    Nopens, I.; Pons, M. N.; Rosen, C.; Steyer, J. P.; Vanrolleghem, P. A.
    Benchmark Simulation Model No. 2 (BSM2).
    http://iwa-mia.org/wp-content/uploads/2022/09/TR3_BSM_TG_Tech_Report_no_3_BSM2_General_Description.pdf.

'''

import os, numpy as np, pandas as pd, qsdsan as qs
from qsdsan import (
    process_models as pc,
    unit_operations as su,
    WasteStream,
    )
from qsdsan.utils import time_printer, ospath, load_data, get_SRT
from exposan.bsm2 import data_path, figures_path, results_path
from exposan.bsm1 import data_path as bsm1_path

__all__ = ('create_system',)


# %%

# =============================================================================
# Initial-condition loaders (lazy, per biokinetic family)
# =============================================================================

def _load_asm1_init():
    '''Initial conditions for the ASM1/ADM1 "bsm2" plant (from CSVs).'''
    asinit = pd.read_csv(os.path.join(data_path, 'asm1init.csv'), index_col=0)
    asm1init = asinit.to_dict('index')
    settler1dinit = pd.read_csv(
        os.path.join(data_path, 'settler1dinit.csv'), index_col=0).to_dict('index')
    adm1init = pd.read_csv(
        os.path.join(data_path, 'adm1init.csv'), index_col=0).to_dict('index')
    return asinit, asm1init, settler1dinit, adm1init


def _load_masm2d_init():
    '''Initial conditions for the mASM2d/ADM1p "bsm2p"/"bsm1p" systems (from xlsx).'''
    dfs = load_data(ospath.join(data_path, 'bsm2p_init.xlsx'), sheet=None)
    inf_concs = dfs['asm'].iloc[0].to_dict()
    c1init = dfs['asm'].iloc[1].to_dict()
    asinit = dfs['asm'].iloc[2:]
    adinit = dfs['adm'].iloc[0].to_dict()
    c2init = dfs['settler'].to_dict('index')
    c2init['s'] = {k: v for k, v in c2init['s'].items() if v > 0}
    c2init['x'] = {k: v for k, v in c2init['x'].items() if v > 0}
    c2init['tss'] = [v for k, v in c2init['tss'].items() if v > 0]
    return inf_concs, c1init, asinit, adinit, c2init


# %%

# =============================================================================
# bsm2: classic BSM2 (ASM1 + ADM1)
# =============================================================================

def _create_bsm2(flowsheet, reactor_model, default_init_conds):
    Q = 20648.361 # influent flowrate [m3/d]
    Q_intr = 3 * Q # activated sludge process internal recycle [m3/d]
    Q_ras = Q # recycle sludge flowrate
    Q_was = 300 # sludge wastage flowrate
    Temp = 273.15+14.85808 # temperature [K]
    T_ad = 273.15+35
    V_an = 1500 # anoxic zone tank volume
    V_ae = 3000 # aerated zone tank volume
    SOSAT1 = 8 # O2 saturation concentration at 15 degC

    flowsheet = flowsheet or qs.Flowsheet('bsm2')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    asinit, asm1init, settler1dinit, adm1init = _load_asm1_init()
    default_inf_kwargs = {
        'concentrations': asm1init['inf'],
        'units': ('m3/d', 'mg/L'),
        }

    # ASM1 components and process model
    cmps_asm1 = pc.create_asm1_cmps()
    asm1 = pc.ASM1(components=cmps_asm1,
                   path=os.path.join(bsm1_path, '_asm1.tsv'))
    thermo_asm1 = qs.get_thermo()
    DO_ID = 'S_O'
    aer1 = aer2 = pc.DiffusedAeration('aer1', DO_ID, KLa=120, DOsat=SOSAT1, V=V_ae)
    aer3 = pc.DiffusedAeration('aer3', DO_ID, KLa=60, DOsat=SOSAT1, V=V_ae)

    # Influent
    inf = WasteStream('inf', T=Temp)
    inf.set_flow_by_concentration(Q, **default_inf_kwargs)
    carb = WasteStream('carbon', T=Temp)
    carb.set_flow_by_concentration(2, {'S_S':400}, units=('m3/d', 'kg/m3'))

    # Primary clarifier using the Otterpohl-Freund model
    C1 = su.PrimaryClarifierBSM2(
        'C1', ins=(inf, 'reject'),
        outs=('C1_eff', 'C1_underflow'),
        isdynamic=True,
        volume=900,
        f_corr=0.65,
        ratio_uf=0.007, # f_PS
        )

    # Activated-sludge train (BSM1 reactors)
    if reactor_model == 'CSTR':
        A1 = su.CSTR('A1', ins=[C1-0, 'RWW', 'RAS', carb], V_max=V_an,
                     aeration=None, suspended_growth_model=asm1)
        A2 = su.CSTR('A2', A1-0, V_max=V_an,
                     aeration=None, suspended_growth_model=asm1)
        O1 = su.CSTR('O1', A2-0, V_max=V_ae, aeration=aer1,
                     DO_ID=DO_ID, suspended_growth_model=asm1)
        O2 = su.CSTR('O2', O1-0, V_max=V_ae, aeration=aer2,
                     DO_ID=DO_ID, suspended_growth_model=asm1)
        O3 = su.CSTR('O3', O2-0, [1-A1, 'treated'],
                     split=[0.6, 0.4],
                     V_max=V_ae, aeration=aer3,
                     DO_ID=DO_ID, suspended_growth_model=asm1)
        # 10-layer one-dimensional settler model, Table 4
        C2 = su.FlatBottomCircularClarifier(
            'C2', O3-1, ['effluent', 2-A1, 'WAS'],
            underflow=Q_ras, wastage=Q_was,
            surface_area=1500, height=4, N_layer=10,
            feed_layer=5, # from top to bottom, 6th if from bottom to top
            X_threshold=3000, v_max=474, v_max_practical=250,
            rh=5.76e-4, rp=2.86e-3, fns=2.28e-3,
            )
        as_units = (A1, A2, O1, O2, O3)
        as_recycle = (O3-0,)
        as_track = (A1, O3)
    elif reactor_model == 'PFR':
        AS = su.PFR('AS', ins=[C1-0, 'RAS', carb], outs='treated',
                    V_tanks=[V_an]*2+[V_ae]*3,
                    influent_fractions=[[1,0,0,0,0]]*3,
                    internal_recycles=[(4,0,Q_intr)],
                    kLa=[0,0,120,120,60], DO_ID=DO_ID, DO_sat=SOSAT1,
                    suspended_growth_model=asm1)
        C2 = su.FlatBottomCircularClarifier(
            'C2', AS-0, ['effluent', 1-AS, 'WAS'],
            underflow=Q_ras, wastage=Q_was,
            surface_area=1500, height=4, N_layer=10,
            feed_layer=5,
            X_threshold=3000, v_max=474, v_max_practical=250,
            rh=5.76e-4, rp=2.86e-3, fns=2.28e-3,
            )
        as_units = (AS,)
        as_recycle = ()
        as_track = (AS,)
    else:
        raise ValueError('`reactor_model` can only be "CSTR" or "PFR", '
                         f'not {reactor_model}.')

    TC1 = su.Thickener('TC1', C2-2, outs=['thickened_sludge', ''],
                       thickening_perc=7, TSS_removal_perc=98)
    M1 = su.Mixer('M1', ins=(C1-1, TC1-0))

    # Switch to ADM1 components for the anaerobic digester
    cmps_adm1 = pc.create_adm1_cmps()
    thermo_adm1 = qs.get_thermo()
    adm1 = pc.ADM1()
    cmps_adm1.X_I.i_N = cmps_asm1.X_I.i_N # slight difference

    J1 = su.ASMtoADM('J1', upstream=M1-0, thermo=thermo_adm1, isdynamic=True,
                     adm1_model=adm1, T=T_ad, pH=7.2631) # WAS is C1.outs[2]
    AD1 = su.AnaerobicCSTR('AD1', ins=J1.outs[0], outs=('biogas', 'AD_eff'), isdynamic=True,
                           V_liq=3400, V_gas=300, T=T_ad, model=adm1,)
    AD1.algebraic_h2 = True
    # Switch back to ASM1 components
    J2 = su.ADMtoASM('J2', upstream=AD1-1, thermo=thermo_asm1, isdynamic=True, adm1_model=adm1)
    J2.bio_to_xs = 0.79
    qs.set_thermo(thermo_asm1)

    # Dewatering
    C3 = su.Centrifuge(ID='C3', ins=J2-0, outs=['digested_sludge', ''],
                       thickening_perc=28, TSS_removal_perc=96.29)

    M2 = su.Mixer('M2', ins=(TC1-1, C3-1), outs=1-C1)

    #!!! Should have a storage tank with HRT = 1,
    # where the outs should have a bypass stream and an out stream.
    # Now equivalent to 100% bypass.

    if default_init_conds:
        C1.set_init_conc(**asm1init['C1'])
        if reactor_model == 'CSTR':
            for i in ('A1', 'A2', 'O1', 'O2', 'O3'):
                getattr(flowsheet.unit, i).set_init_conc(**asm1init[i])
        else:
            AS.set_init_conc(concentrations=asinit.iloc[1:-1])
        C2.set_init_TSS(list(settler1dinit['C2'].values()))
        AD1.set_init_conc(**adm1init['AD1'])

    sys = qs.System('bsm2_sys',
                    path=(C1, *as_units, C2, TC1, M1, J1, AD1, J2, C3, M2),
                    recycle=(*as_recycle, C2-1, M2-0),
                    )
    sys.set_tolerance(mol=1e-5, rmol=1e-5)
    sys.maxiter = 5000
    sys.set_dynamic_tracker(C1, *as_track, C2, J1, AD1, J2, C3)

    return sys


# %%

# =============================================================================
# bsm2p: phosphorus-extended BSM2 (mASM2d + ADM1p)
# =============================================================================

# mmp kinetics ('Musvoto' or 'KM')
mmp = 'KM'

def _create_bsm2p(flowsheet, reactor_model, default_init_conds):
    Q = 20648.4 # influent flowrate [m3/d]
    Q_intr = 3 * Q # activated sludge process internal recycle [m3/d]
    Q_ras = 1 * Q # recycle sludge flowrate
    Q_was = 600 # sludge wastage flowrate
    Temp = 273.15+20 # temperature [K]
    T_ad = 273.15+35
    V_anae = 1000 # anaerobic zone tank volume
    V_anox = 1500 # anoxic zone tank volume
    V_ae = 3000 # aerated zone tank volume
    SOSAT1 = 8 # O2 saturation concentration at 15 degC

    flowsheet = flowsheet or qs.Flowsheet('bsm2p')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    inf_concs, c1init, asinit, adinit, c2init = _load_masm2d_init()
    default_inf_kwargs = {
        'flow_tot': Q,
        'concentrations': inf_concs,
        'units': ('m3/d', 'mg/L'),
        }

    cmps_asm = pc.create_masm2d_cmps()
    asm = pc.mASM2d(components=cmps_asm,
                    electron_acceptor_dependent_decay=True,)
    thermo_asm = qs.get_thermo()

    # Influent
    inf = WasteStream('inf', T=Temp)
    inf.set_flow_by_concentration(**default_inf_kwargs)
    carb = WasteStream('carbon', T=Temp)
    carb.set_flow_by_concentration(2, {'S_A':400}, units=('m3/d', 'kg/m3'))

    # Primary clarifier using the Otterpohl-Freund model
    C1 = su.PrimaryClarifierBSM2(
        'C1', ins=(inf, 'reject'),
        outs=('C1_eff', 'C1_underflow'),
        isdynamic=True,
        volume=900,
        f_corr=0.65,
        ratio_uf=0.007, # f_PS
        )

    DO_ID = 'S_O2'
    gstrip = True
    c2_kwargs = dict(
        underflow=Q_ras, wastage=Q_was,
        surface_area=1500, height=4, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=474, v_max_practical=250,
        rh=5.76e-4, rp=2.86e-3, fns=2.28e-3,
        )

    # Activated-sludge train (anaerobic + anoxic + aerated, for bio-P removal)
    if reactor_model == 'PFR':
        AS = su.PFR('AS', ins=[C1-0, 'RAS', carb], outs='treated',
                    N_tanks_in_series=7,
                    V_tanks=[V_anae]*2+[V_anox]*2+[V_ae]*3,
                    influent_fractions=[[1]+[0]*6]*2 + [[0,0,1,0,0,0,0]],
                    internal_recycles=[(6,2,Q_intr*1.1)],
                    kLa=[0]*4+[240,120,60],
                    DO_ID=DO_ID, DO_sat=SOSAT1,
                    suspended_growth_model=asm,
                    gas_stripping=gstrip)
        C2 = su.FlatBottomCircularClarifier(
            'C2', AS-0, ['effluent', 1-AS, 'WAS'],
            **c2_kwargs
            )
        as_units = (AS,)
        as_recycle = ()
    elif reactor_model == 'CSTR':
        aer1 = aer2 = pc.DiffusedAeration('aer1', DO_ID, KLa=120, DOsat=SOSAT1, V=V_ae)
        aer3 = pc.DiffusedAeration('aer3', DO_ID, KLa=60, DOsat=SOSAT1, V=V_ae)
        anae_kwargs = dict(V_max=V_anae, aeration=None, suspended_growth_model=asm, gas_stripping=gstrip)
        anox_kwargs = dict(V_max=V_anox, aeration=None, suspended_growth_model=asm, gas_stripping=gstrip)
        ae_kwargs = dict(V_max=V_ae, DO_ID=DO_ID, suspended_growth_model=asm, gas_stripping=gstrip)
        A1 = su.CSTR('A1', ins=[C1-0, 'RAS', carb], **anae_kwargs)
        A2 = su.CSTR('A2', A1-0, **anae_kwargs)
        A3 = su.CSTR('A3', [A2-0, 'RWW'], **anox_kwargs)
        A4 = su.CSTR('A4', A3-0, **anox_kwargs)
        O1 = su.CSTR('O1', A4-0, aeration=aer1, **ae_kwargs)
        O2 = su.CSTR('O2', O1-0, aeration=aer2, **ae_kwargs)
        O3 = su.CSTR('O3', O2-0, [1-A3, 'treated'], split=[Q_intr, Q+Q_ras],
                     aeration=aer3, **ae_kwargs)
        C2 = su.FlatBottomCircularClarifier(
            'C2', O3-1, ['effluent', 1-A1, 'WAS'],
            **c2_kwargs
            )
        as_units = (A1, A2, A3, A4, O1, O2, O3)
        as_recycle = (O3-0,)
    else:
        raise ValueError('`reactor_model` can only be "CSTR" or "PFR", '
                         f'not {reactor_model}.')

    TC1 = su.IdealClarifier('TC1', C2-2, outs=['', 'thickened_WAS'],
                            sludge_flow_rate=30.9,
                            sludge_MLSS=7.0e4,)
    M1 = su.Mixer('M1', ins=(C1-1, TC1-1))

    # Switch to ADM1p components for the anaerobic digester
    pc.create_adm1p_cmps()
    thermo_adm = qs.get_thermo()
    adm = pc.ADM1p(
        f_bu_su=0.1328, f_pro_su=0.2691, f_ac_su=0.4076,
        q_ch_hyd=0.3, q_pr_hyd=0.3, q_li_hyd=0.3,
        mmp_kinetics=mmp,
        )

    J1 = su.mASM2dtoADM1p('J1', upstream=M1-0, thermo=thermo_adm, isdynamic=True,
                          adm1_model=adm, asm2d_model=asm)
    AD = su.AnaerobicCSTR('AD', ins=J1.outs[0], outs=('biogas', 'AD_eff'), isdynamic=True,
                          V_liq=3400, V_gas=300, T=T_ad, model=adm,
                          pH_ctrl=7.0,)
    AD.algebraic_h2 = False
    J2 = su.ADM1ptomASM2d('J2', upstream=AD-1, thermo=thermo_asm, isdynamic=True,
                          adm1_model=adm, asm2d_model=asm)
    # Switch back to ASM components
    qs.set_thermo(thermo_asm)

    # Dewatering
    C3 = su.IdealClarifier('C3', J2-0, outs=['', 'digested_sludge'],
                           sludge_flow_rate=9.6,
                           sludge_MLSS=2.8e5,)
    M2 = su.Mixer('M2', ins=(TC1-0, C3-0))
    HD = su.HydraulicDelay('HD', ins=M2-0, outs=1-C1)

    if default_init_conds:
        C1.set_init_conc(**c1init)
        if reactor_model == 'PFR':
            AS.set_init_conc(concentrations=asinit)
        else:
            asdct = asinit.to_dict('index')
            for i in as_units:
                i.set_init_conc(**asdct[i.ID])
        C2.set_init_solubles(**c2init['s'])
        C2.set_init_sludge_solids(**c2init['x'])
        C2.set_init_TSS(c2init['tss'])
        AD.set_init_conc(**adinit)

    sys = qs.System('bsm2p',
                    path=(C1, *as_units, C2, TC1, M1, J1, AD, J2, C3, M2, HD),
                    recycle=(*as_recycle, C2-1, HD-0),
                    )
    sys.set_dynamic_tracker(AD, C2-0)

    return sys


# %%

# =============================================================================
# Public entry point
# =============================================================================

_default_reactor_model = {'bsm2': 'CSTR', 'bsm2p': 'PFR'}

# Default dynamic-simulation settings per kind, used by `load` (and the
# __main__ harness below). Keep these in sync with the values documented in
# the `if __name__ == '__main__'` block.
default_simulate_kwargs = {
    'bsm2':  dict(t_span=(0, 30),  method='RK23'),
    'bsm2p': dict(t_span=(0, 300), method='BDF'),
    }

def create_system(flowsheet=None, kind='bsm2', reactor_model=None,
                  default_init_conds=True):
    '''
    Create a Benchmark Simulation Model No. 2 system.

    Parameters
    ----------
    flowsheet : obj
        Flowsheet where this system will be created on.
    kind : str
        Which BSM2 configuration to build:

        - "bsm2" (default): the classic BSM2 using ASM1 + ADM1.
        - "bsm2p": a phosphorus-removal variant using mASM2d + ADM1p.
    reactor_model : str, optional
        "CSTR" to model each zone of the activated-sludge reactor as a CSTR
        with explicit internal recirculation, or "PFR" to model the entire
        reactor as a single unit with implicit internal recirculation.
        Defaults to "CSTR" for "bsm2" and "PFR" for "bsm2p".
    default_init_conds : bool
        Whether to use the default initial conditions shipped with the module.

    See Also
    --------
    `create_subsys` builds ``bsm1p``, the activated-sludge sub-train of the
    ``bsm2p`` plant (a diagnostic system, not a benchmark).
    '''
    kind = kind.lower()
    reactor_model = (reactor_model or _default_reactor_model.get(kind)) or ''
    reactor_model = reactor_model.upper()
    if kind == 'bsm2':
        return _create_bsm2(flowsheet, reactor_model, default_init_conds)
    elif kind == 'bsm2p':
        return _create_bsm2p(flowsheet, reactor_model, default_init_conds)
    raise ValueError(f'`kind` can only be "bsm2" or "bsm2p", not "{kind}".')


# %%

# =============================================================================
# bsm1p: activated-sludge sub-train of bsm2p (diagnostic, not a benchmark)
# =============================================================================

default_PE_concs = dict(
    S_N2=18, S_NH4=25, S_F=87.0, S_I=21.8, X_S=112.7+39.6, X_I=29.0, S_PO4=8.0,
    S_IC=75.6, S_Ca=140, S_Mg=50, S_K=28, S_Na=3.76*23, S_Cl=12*35.45,
    S_A=30
    )

def create_subsys():
    '''
    Build ``bsm1p``: the activated-sludge sub-train (mASM2d) of the ``bsm2p``
    plant, fed a primary-effluent-like influent. Used to study the bioreactor
    train in isolation; it is not a published benchmark and is not exposed
    through `create_system`.
    '''
    Q = 20648.4
    Q_intr = 3 * Q
    Q_ras = 1 * Q
    Q_was = 600
    Temp = 273.15+20
    V_anae = 1000
    V_anox = 1500
    V_ae = 3000
    SOSAT1 = 8

    flowsheet = qs.Flowsheet('bsm1p')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    _, _, asinit, _, c2init = _load_masm2d_init()

    cmps_asm = pc.create_masm2d_cmps()
    asm = pc.mASM2d(components=cmps_asm,
                    electron_acceptor_dependent_decay=True,
                    k_h=2.46, mu_H=4.23, q_fe=2.11, b_H=0.28, mu_PAO=0.82,
                    q_PP=1.23, q_PHA=2.46, b_PAO=0.14, b_PP=0.14, b_PHA=0.14,
                    mu_AUT=0.61, b_AUT=0.09)

    # Influent
    inf = WasteStream('inf', T=Temp)
    inf.set_flow_by_concentration(20446, default_PE_concs, units=('m3/d', 'mg/L'))

    DO_ID = 'S_O2'
    aer1 = aer2 = pc.DiffusedAeration('aer1', DO_ID, KLa=120, DOsat=SOSAT1, V=V_ae)
    aer3 = pc.DiffusedAeration('aer3', DO_ID, KLa=60, DOsat=SOSAT1, V=V_ae)
    gstrip = True
    anae_kwargs = dict(V_max=V_anae, aeration=None, suspended_growth_model=asm, gas_stripping=gstrip)
    anox_kwargs = dict(V_max=V_anox, aeration=None, suspended_growth_model=asm, gas_stripping=gstrip)
    ae_kwargs = dict(V_max=V_ae, DO_ID=DO_ID, suspended_growth_model=asm, gas_stripping=gstrip)
    c2_kwargs = dict(
        underflow=Q_ras, wastage=Q_was,
        surface_area=1500, height=4, N_layer=10, feed_layer=5,
        X_threshold=3000, v_max=474, v_max_practical=250,
        rh=5.76e-4, rp=2.86e-3, fns=2.28e-3
        )

    A1 = su.CSTR('A1', ins=[inf, 'RAS'], **anae_kwargs)
    A2 = su.CSTR('A2', A1-0, **anae_kwargs)
    A3 = su.CSTR('A3', [A2-0, 'RWW'], **anox_kwargs)
    A4 = su.CSTR('A4', A3-0, **anox_kwargs)
    O1 = su.CSTR('O1', A4-0, aeration=aer1, **ae_kwargs)
    O2 = su.CSTR('O2', O1-0, aeration=aer2, **ae_kwargs)
    O3 = su.CSTR('O3', O2-0, [1-A3, 'treated'], split=[Q_intr, Q+Q_ras],
                 aeration=aer3, **ae_kwargs)
    C2 = su.FlatBottomCircularClarifier(
        'C2', O3-1, ['effluent', 1-A1, 'WAS'],
        **c2_kwargs
        )

    asdct = asinit.to_dict('index')
    for i in (A1, A2, A3, A4, O1, O2, O3):
        i.set_init_conc(**asdct[i.ID])
    C2.set_init_solubles(**c2init['s'])
    C2.set_init_sludge_solids(**c2init['x'])
    C2.set_init_TSS(c2init['tss'])

    sub = qs.System('bsm1p',
                    path=(A1, A2, A3, A4, O1, O2, O3, C2),
                    recycle=(O3-0, C2-1,)
                    )
    sub.set_dynamic_tracker(A1, A3, O3, C2-0)

    return sub


# %%

@time_printer
def run(sys, t, t_step, method=None, **kwargs):
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}')
    print(f'Time span 0-{t}d \n')

    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0,t),
        t_eval=np.arange(0, t+t_step, t_step),
        method=method,
        print_t=True,
        # rtol=1e-2,
        # atol=1e-3,
        # export_state_to=f'results/sol_{t}d_{method}.xlsx',
        **kwargs)

#%%
if __name__ == '__main__':
    kind = 'bsm2'
    # kind = 'bsm2p'
    sys = create_system(kind=kind)
    dct = globals()
    dct.update(sys.flowsheet.to_dict())

    t = default_simulate_kwargs[kind]['t_span'][1]
    method = default_simulate_kwargs[kind]['method']
    # method options: 'RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA'
    t_step = 1 if kind == 'bsm2' else 0.1

    run(sys, t, t_step, method=method)
    sys.diagram()
    # sys.diagram(file=os.path.join(figures_path, 'bsm2_sys'), format='png')
