#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
This module is developed by:
    
    Aaron Marszewski <aaronpm3@illinois.edu>
    Rishabh Puri <rp34@illinois.edu>
    Siqi Tang <siqit@outlook.com>
    Yuyao Huang <yuyaoh2@illinois.edu>
    

This python file is used to perform uncertainty and sensitivity analysis for Enviroloo Clear Reinvented Toilet system.

TOGGLE: el.INCLUDE_STRUVITE (set in exposan/enviroloo/__init__.py)
    True  -> system includes urine diversion + MgCl2 dosing + Struvite
             Reactor (SR) + Struvite Redissolution (SRD) upstream of CT
    False -> baseline system, CT receives combined toilet waste directly
             (no urine diversion, no SR/SRD)

Urine diversion mass balance (when INCLUDE_STRUVITE=True), corrected
(Option B): mass is split by the literature excretion fraction
(f_urine_N, f_urine_P; Rose et al. 2015 via Lohman et al. 2020 SI Table S2),
then divided by each stream's OWN volume (Q_urine or Q_remaining) to
compute physically correct concentrations -- this produces a urine stream
that is correctly MORE concentrated than the combined stream, consistent
with Etter et al. 2011 (stored urine: 195-388 mg P/L) vs. this system's
fresh, diverted, non-stored urine (~65 mg P/L).

MgCl2 dosing target Mg:P = 1.1:1 (Etter et al. 2011).
SR eff_PO4_mgL = 6.55 mg/L expected (90% removal, Etter et al. 2011/
Maurer et al. 2006), Triangular(1.31, 6.55, 27.50) -- see _EL_SR.tsv.
SR precip_yield = 0.70 expected (EnviroLoo field reactor, 68.7% harvest
efficiency), Triangular(0.55, 0.70, 0.85) -- see _EL_SR.tsv.
'''

# %% 
import qsdsan as qs 
from qsdsan import (
      ImpactItem, 
     )
from qsdsan.utils import (
    ospath, 
    time_printer, 
    load_data, 
    get_SRT,
    )
from qsdsan import Components, processes as pc, sanunits as su

from exposan.enviroloo import _units as elu
from exposan.enviroloo import (
    data_path,
    _load_lca_data,
    ppl, baseline_ppl, scale_factor
    )
from exposan import enviroloo as el

folder = ospath.dirname(__file__)

__all__ = ('create_systemEL',)
#%%

'''
Name notes:

To make programming more convenient, we use the following names for the units in the system:

Control Room Housing: ELH 
Collection tank: CT
Primary clarifier: PC
Anoxic tank: AnoT
Aerobic tank: AerT
Membrane tank: MemT
Clear water tank: CWT
Pressure tank: PT
Primary clarifier return pump: P_PC_return
Glucose agitation pump: P_AnoT_agitation
Glucose dosing pump: P_AnoT_dosing
Anoxic mixing pump: P_AnoT_mixing.
PAC agitation pump: P_AerT_agitation
PAC dosing pump: P_AerT_dosing
Self-priming pump: P_MT_selfpriming
Lift pump: P_CT_lift
Aerobic blower: B_AerT
Membrane blower: B_MemT
Clear water pump: P_CWT
Ozone generator: O3_gen
Micro bubble pump (Ozone dosing pump): P_O3_dosing
Air dissolving pump: P_AirDissolved

Struvite Reactor: SR (only present if el.INCLUDE_STRUVITE = True)
Struvite Redissolution: SRD (only present if el.INCLUDE_STRUVITE = True)
MgCl2 Mixer: M_MgCl2 (only present if el.INCLUDE_STRUVITE = True)

'''

Temp = 273.15 + 20   # temperature [K]

# Hydraulic flows
Q_w   = 4.48 * scale_factor    # m3/day -- total wastewater flow
Q_ras = 2.24 * scale_factor    # m3/day -- nitrate return flow
Q_was = 0.05 * 24              # m3/day -- waste activated sludge flow

biomass_IDs = ('X_H', 'X_AUT', 'X_PAO')

# =============================================================================
# Field-measured/design-basis toilet waste concentrations
# Combined urine + feces + flush water stream (no diversion)
# =============================================================================
toilet_waste = {
    'S_NH4': 19,      # mg/L
    'S_NO3': 1.6,     # mg/L
    'S_PO4': 16.7,    # mg/L
    'S_O2':  3.0,     # mg/L
    'S_F':   65,      # mg/L
    'S_I':   25,      # mg/L
    'X_S':   40,      # mg/L
    'S_IC':  0.148,   # mg/L
    'S_K':   0.0694,  # mg/L
    'S_Mg':  0.00833, # mg/L
    'S_Ca':  0.0117,  # mg/L
    'S_Na':  0.775,   # mg/L
    'S_Cl':  0.687,   # mg/L
}

# =============================================================================
# Urine diversion parameters -- only used when el.INCLUDE_STRUVITE = True
# =============================================================================

# Per capita urine volume (Lohman et al. 2020, SI Table S16)
V_urine_Lpcpd = 1.4            # L/cap/day
Q_urine = ppl * V_urine_Lpcpd / 1000.0   # m3/day

# Urine/feces N and P split fractions (Lohman et al. 2020 SI Table S2,
# citing Rose et al. 2015)
f_urine_N = 0.88    # fraction of total N excreted in urine
f_urine_P = 0.61    # fraction of total P excreted in urine

Q_total     = Q_w * 0.9
Q_remaining = Q_total - Q_urine

# Components for which a literature excretion fraction is defined.
_excretion_fraction = {
    'S_NH4': f_urine_N,   # 88% of N excreted in urine (Rose et al. 2015)
    'S_PO4': f_urine_P,   # 61% of P excreted in urine (Rose et al. 2015)
}

def _split_concentration(component):
    """
    CORRECTED urine diversion mass balance (Option B).

    Split a component's concentration between the urine and remaining
    (feces + flush water) streams via mass balance: mass is split by the
    literature excretion fraction, then converted back to concentration
    using each stream's own volume. This fixes the prior error of applying
    the fraction directly to the combined-stream concentration, which
    incorrectly produced a urine stream LESS concentrated than the
    combined stream (physically backwards, since urine occupies a much
    smaller share of flow (~15.6%) than its share of excreted N/P mass).
    """
    C_total = toilet_waste.get(component, 0.0)
    M_total = C_total * Q_total

    f = _excretion_fraction.get(component)
    if f is not None:
        # Components with a defined excretion split (S_NH4, S_PO4)
        M_urine     = M_total * f
        M_remaining = M_total - M_urine
        C_urine     = M_urine / Q_urine
        C_remaining = M_remaining / Q_remaining
    else:
        # Components without a defined literature excretion split:
        #   S_O2, S_F, S_I, X_S: assumed entirely in feces/flush stream
        #     (zero in urine; primarily flush-water/feces-derived)
        #   S_IC, S_K, S_Mg, S_Ca, S_Na, S_Cl: retained at the combined
        #     concentration in BOTH streams (no excretion-fraction basis
        #     currently defined; conservative simplifying assumption,
        #     flagged for future refinement using e.g. Rose et al. 2015
        #     K/Mg/Ca excretion fractions if available)
        if component in ('S_O2', 'S_F', 'S_I', 'X_S'):
            C_urine     = 0.0
            C_remaining = C_total * Q_total / Q_remaining
        else:
            C_urine     = C_total
            C_remaining = C_total

    return max(C_urine, 0.0), max(C_remaining, 0.0)


def _build_urine_split():
    """Compute urine_waste and remaining_waste dicts. Only called when
    el.INCLUDE_STRUVITE = True."""
    urine_waste     = {}
    remaining_waste = {}
    for _component in toilet_waste:
        _c_urine, _c_remaining = _split_concentration(_component)
        urine_waste[_component]     = _c_urine
        remaining_waste[_component] = _c_remaining

    print(f"\n--- Corrected Urine Diversion Split (Option B) ---")
    print(f"Q_urine     : {Q_urine:.4f} m3/d  ({Q_urine/Q_total*100:.1f}% of flow)")
    print(f"Q_remaining : {Q_remaining:.4f} m3/d  ({Q_remaining/Q_total*100:.1f}% of flow)")
    print(f"S_PO4  -> urine: {urine_waste['S_PO4']:.3f} mg/L   remaining: {remaining_waste['S_PO4']:.3f} mg/L")
    print(f"S_NH4  -> urine: {urine_waste['S_NH4']:.3f} mg/L   remaining: {remaining_waste['S_NH4']:.3f} mg/L")
    print("---------------------------------------------------\n")

    return urine_waste, remaining_waste


def _compute_MgCl2_dose(urine_waste):
    """
    MgCl2 dosing parameters (Etter et al. 2011, Mg:P = 1.1:1 molar).
    Only used when el.INCLUDE_STRUVITE = True.
    """
    MW_P       = 30.97
    MW_Mg      = 24.31
    MW_MgCl2   = 95.21
    Mg_P_ratio = 1.1   # molar, Etter et al. 2011

    C_PO4_urine   = urine_waste['S_PO4']           # mg/L (corrected, ~65.5 mg/L)
    C_Mg_required = C_PO4_urine * (MW_Mg / MW_P) * Mg_P_ratio  # mg/L total Mg needed
    C_Mg_natural  = urine_waste['S_Mg']             # mg/L already in urine
    delta_Mg_mgL  = C_Mg_required - C_Mg_natural    # mg/L to add

    Mg_add_g_d   = delta_Mg_mgL * Q_urine                      # g/d
    Mg_add_kg_hr = Mg_add_g_d / 24.0 / 1000.0                  # kg/hr (for WasteStream)
    MgCl2_kg_d   = Mg_add_g_d * (MW_MgCl2 / MW_Mg) / 1000.0    # kg/d  (for cost)

    print(f"--- MgCl2 Dosing Parameters (corrected) ---")
    print(f"C_PO4 in urine stream  : {C_PO4_urine:.3f} mg/L")
    print(f"C_Mg required (1.1:1)  : {C_Mg_required:.3f} mg/L")
    print(f"C_Mg natural in urine  : {C_Mg_natural:.4f} mg/L")
    print(f"Delta Mg to add        : {delta_Mg_mgL:.3f} mg/L")
    print(f"Q_urine                : {Q_urine:.3f} m3/d")
    print(f"Mg mass flow added     : {Mg_add_g_d:.3f} g/d = {Mg_add_kg_hr:.6f} kg/hr")
    print(f"MgCl2 dose (cost acctg): {MgCl2_kg_d:.4f} kg/d")
    print(f"Mg:P molar ratio       : {Mg_P_ratio} : 1")
    print("-------------------------------\n")

    return Mg_add_kg_hr, MgCl2_kg_d, MW_P, MW_Mg


# %% Create Universal Units and Functions
def batch_init(sys, path, sheet):
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    for k in sys.units:
        if k.ID.startswith(('O', 'A', 'B')):
            k.set_init_conc(**dct[k.ID])


# %% Create EnviroLoo Clear system

def create_components(set_thermo=True):
    masm2d_cmps = pc.create_masm2d_cmps(set_thermo=False)
    cmps = Components([*masm2d_cmps])
    # cmps.compile(ignore_inaccurate_molar_weight=True)
    cmps.compile(skip_checks=True)
    
    if set_thermo:
        qs.set_thermo(cmps)
    return cmps


# %%

def create_systemEL(flowsheet=None, inf_kwargs={}, masm_kwargs={}, init_conds={},
                  aeration_processes=()):
    """
    EnviroLoo Clear system.

    When el.INCLUDE_STRUVITE = True, the system includes urine diversion,
    MgCl2 dosing, and a Struvite Reactor (SR) + Struvite Redissolution
    (SRD) upstream of CT:

        urine_ins --+
                     +-- M_MgCl2 (Mixer) --> SR --> SRD --+
        Mg_dose  ----+                                     |
                                                            v
        remaining_ins ------------------------------> CT (Mixer)
                                                            |
                                                           PC -> A1 -> O1
                                                            -> B1 -> ... -> ELH

    When el.INCLUDE_STRUVITE = False (baseline), CT receives the combined
    toilet waste stream directly, with no urine diversion or SR/SRD:

        toilet_ins ---> CT --> PC -> A1 -> O1 -> B1 -> ... -> ELH
    """

    _load_lca_data()
    cmps = create_components()

    masm2d = pc.mASM2d(**masm_kwargs)
    kwargs_O = dict(V_max=6.35, aeration=2,    DO_ID='S_O2',
                    suspended_growth_model=masm2d)
    kwargs_1 = dict(V_max=6.35, aeration=None, DO_ID=None,
                    suspended_growth_model=masm2d)
    kwargs_2 = dict(V_max=2.89, aeration=2,    DO_ID='S_O2',
                    suspended_growth_model=masm2d)

    if el.INCLUDE_STRUVITE:
        # =====================================================================
        # SR-INCLUSIVE PATH: urine diversion + MgCl2 dosing + SR + SRD
        # =====================================================================

        urine_waste, remaining_waste = _build_urine_split()
        Mg_add_kg_hr, MgCl2_kg_d, MW_P, MW_Mg = _compute_MgCl2_dose(urine_waste)

        # ---- 1. Urine stream ----
        urine_ins = qs.WasteStream('urine_ins', T=Temp)
        urine_ins.set_flow_by_concentration(
            Q_urine, concentrations=urine_waste, units=('m3/d', 'mg/L'))

        # ---- 2. MgCl2 dosing stream ----
        Mg_dose = qs.WasteStream('Mg_dose', T=Temp)
        Mg_dose.imass['S_Mg'] = Mg_add_kg_hr   # kg/hr
        item = ImpactItem.get_item('MgCl2_item').copy('Mg_dose_item', set_as_source=True)
        Mg_dose.stream_impact_item = item

        # ---- 3. Mixer: urine + MgCl2 dose, before SR ----
        M_MgCl2 = su.Mixer(
            'M_MgCl2',
            ins=(urine_ins, Mg_dose),
            outs=qs.WasteStream('urine_with_Mg', T=Temp),
        )

        # ---- 4. Remaining stream (feces + flush water) ----
        remaining_ins = qs.WasteStream('remaining_ins', T=Temp)
        remaining_ins.set_flow_by_concentration(
            Q_remaining, concentrations=remaining_waste, units=('m3/d', 'mg/L'))

        # ---- 6. Struvite Reactor (SR) ----
        # eff_PO4_mgL = 6.55 mg/L expected (90% removal, Etter et al. 2011/
        # Maurer et al. 2006), applied to corrected influent ~65.5 mg/L.
        # precip_yield = 0.70 expected (EnviroLoo field reactor, 68.7%
        # harvest efficiency). See _EL_SR.tsv for full uncertainty ranges.
        SR = elu.StruviteReactor(
            'SR',
            ins=M_MgCl2-0,
            outs=(
                qs.WasteStream('struvite_recovered', T=Temp),
                qs.WasteStream('loss',               T=Temp),
                qs.WasteStream('effluent_SR',         T=Temp),
            ),
            isdynamic=True,
            eff_PO4_mgL=6.55,
            precip_yield=0.70,
            dose_MgCl2_kg_d=MgCl2_kg_d,
        )

        # ---- 7. Struvite Redissolution (SRD) ----
        SRD = elu.StruviteRedissolution(
            'SRD',
            ins=SR-2,
            outs=qs.WasteStream('effluent_SRD', T=Temp),
            isdynamic=True,
            k_max=2.61,
            d_p=0.3,
            HRT_min=5.0,
        )

        # ---- 8. Collection Tank (CT): remaining_ins + SRD effluent ----
        effluent_CT = qs.WasteStream('effluent_CT', T=Temp)
        effluent_CT.set_flow_by_concentration(
            Q_w * 0.9, concentrations=toilet_waste, units=('m3/d', 'mg/L'))

        CT = elu.EL_CT(
            'CT',
            ins=(remaining_ins, SRD-0),
            outs=effluent_CT,
            isdynamic=True, V_max=10, aeration=None,
            suspended_growth_model=None,
            ppl=ppl, baseline_ppl=baseline_ppl,
        )

        upstream_units = (M_MgCl2, SR, SRD)

    else:
        # =====================================================================
        # BASELINE PATH: no urine diversion, no SR/SRD
        # =====================================================================

        toilet_ins = qs.WasteStream('toilet_waste', T=Temp)
        toilet_ins.set_flow_by_concentration(
            Q_w * 0.9, concentrations=toilet_waste, units=('m3/d', 'mg/L'))

        effluent_CT = qs.WasteStream('effluent_CT', T=Temp)
        effluent_CT.set_flow_by_concentration(
            Q_w * 0.9, concentrations=toilet_waste, units=('m3/d', 'mg/L'))

        CT = elu.EL_CT(
            'CT', ins=(toilet_ins), outs=(effluent_CT),
            isdynamic=True, V_max=10, aeration=None,
            suspended_growth_model=None,
            ppl=ppl, baseline_ppl=baseline_ppl,
        )

        upstream_units = ()

    # =========================================================================
    # Shared downstream units -- identical for both INCLUDE_STRUVITE branches
    # =========================================================================

    PC = elu.EL_PC(
        'PC', ins=(CT-0, 'RAS_PC'), outs=('effluent_PC_total', 'sludge_PC'),
        ppl=ppl, baseline_ppl=baseline_ppl,
        solids_removal_efficiency=0.85,
        isdynamic=True,
        sludge_flow_rate=Q_w * 0.01,
    )

    Glucose = qs.WasteStream('Glucose_Dose', T=Temp)

    A1 = elu.EL_Anoxic(
        'A1', ins=(PC-0, 'RAS_A1', Glucose), outs=('effluent_AnoxT',),
        isdynamic=True, ppl=ppl, baseline_ppl=baseline_ppl, **kwargs_1,
    )
    item = ImpactItem.get_item('Glucose_item').copy('A1_glucose_item', set_as_source=True)
    A1.ins[2].stream_impact_item = item

    PAC = qs.WasteStream('PAC_Dose', T=Temp)

    O1 = elu.EL_Aerobic(
        'O1', ins=(A1-0, PAC), outs=('effluent_AeroT',),
        isdynamic=True, ppl=ppl, baseline_ppl=baseline_ppl, **kwargs_O,
    )
    item = ImpactItem.get_item('PAC_item').copy('O1_PAC_item', set_as_source=True)
    O1.ins[1].stream_impact_item = item

    B1 = elu.EL_CMMBR(
        'B1', ins=O1-0, outs=('effluent_MembT', 'sludge_MembT'),
        isdynamic=True, ppl=ppl, baseline_ppl=baseline_ppl,
        pumped_flow=(Q_ras * 2 + Q_was),
        **kwargs_2,
    )

    S2 = su.Splitter('S2', ins=B1-1, outs=['RAS', 'WAS'],
                      split=Q_ras * 2 / (Q_ras * 2 + Q_was))
    S1 = su.Splitter('S1', ins=S2-0, outs=[1-PC, 1-A1], split=0.01)

    CWT = elu.EL_CWT(
        'CWT', ins=(B1-0), outs=('effluent_CWT'),
        isdynamic=True, V_max=12, aeration=None, suspended_growth_model=None,
        ppl=ppl, baseline_ppl=baseline_ppl,
    )

    PV = elu.EL_WindSolar('PV', ins=CWT-0, outs='effluent', isdynamic=True,
                          ppl=ppl, baseline_ppl=baseline_ppl)

    ELH = elu.EL_Housing('ELH', ins=PV-0, outs='effluent', isdynamic=True,
                         ppl=ppl, baseline_ppl=baseline_ppl)

    # =========================================================================
    # Build system -- path includes upstream_units (empty tuple if baseline)
    # =========================================================================

    sys = qs.System(
        'EL',
        path=(*upstream_units, CT, PC, A1, O1, B1, S2, S1, CWT, PV, ELH),
    )
    sys.set_dynamic_tracker(A1, O1, B1, B1-0, B1-1)

    batch_init(
        sys, ospath.join(data_path, "units_data/bsm2p_init.xlsx"), sheet='el')

    GWP        = qs.ImpactIndicator('GlobalWarming', alias='GWP', unit='kg CO2-eq')
    Ecosystems = qs.ImpactIndicator('H_Ecosystems',  alias='Ecosystems', unit='points')
    Health     = qs.ImpactIndicator('H_Health',       alias='Health', unit='points')
    Resources  = qs.ImpactIndicator('H_Resources',    alias='Resources', unit='points')

    tea = qs.TEA(system=sys, discount_rate=0.05, lifetime=10, simulate_system=False)
    lca = qs.LCA(system=sys, lifetime=10, lifetime_unit='yr',
                 indicators=(GWP, Ecosystems, Health, Resources), simulate_system=False)

    return sys


# %%

@time_printer
def run(t, method=None, **kwargs):
    sys = create_systemEL()
    batch_init(
        sys, ospath.join(data_path, "units_data/bsm2p_init.xlsx"), sheet='el')
    return sys


if __name__ == '__main__':
    t = 100
    method = 'RK23'
    msg = f'Method {method}'
    print(f'\n{msg}\n{"-"*len(msg)}')
    print(f'Time span 0-{t}d \n')
    print(f'INCLUDE_STRUVITE = {el.INCLUDE_STRUVITE}\n')

    sysEL = run(t, method=method)

    fs = sysEL.flowsheet.stream
    fu = sysEL.flowsheet.unit

    sysEL.simulate(
        t_span=(0, t),
        method=method,
    )
    sysEL.diagram()

    act_units = [
        u.ID for u in sysEL.units
        if isinstance(u, (elu.EL_Aerobic, elu.EL_Anoxic, elu.EL_CMMBR))
    ]
    print(f'act_units = {act_units}')

    srt = get_SRT(sysEL, biomass_IDs, wastage=[fs.WAS, fs.effluent_MembT],
                  active_unit_IDs=act_units)
    print(f'Estimated SRT assuming at steady state is {srt} days\n')

    fig, axis = fs.effluent_MembT.scope.plot_time_series(
        ('S_F', 'S_A', 'X_H', 'S_NH4', 'S_NO3', 'S_PO4', 'X_I', 'S_I', 'S_N2'))

    # ---- SR-specific diagnostics, only meaningful if INCLUDE_STRUVITE=True ----
    if el.INCLUDE_STRUVITE:
        sr_inf = fu.SR.ins[0]
        MW_P, MW_Mg = 30.97, 24.31
        C_Mg_sr  = float(sr_inf.iconc['S_Mg'])
        C_PO4_sr = float(sr_inf.iconc['S_PO4'])
        mol_Mg   = C_Mg_sr  / MW_Mg
        mol_PO4  = C_PO4_sr / MW_P
        print(f"\n--- SR Influent Verification ---")
        print(f"S_Mg  in SR influent : {C_Mg_sr:.3f} mg/L")
        print(f"S_PO4 in SR influent : {C_PO4_sr:.3f} mg/L")
        print(f"Mg:P molar ratio     : {mol_Mg/mol_PO4:.3f} : 1  (target 1.1:1)")
        print("--------------------------------\n")

        elu.StruviteReactor.report_recovery(fu.SR)

        print('--- SR Effluent (after struvite precipitation) ---')
        print(f"S_PO4 : {fs.effluent_SR.iconc['S_PO4']:.3f} mg/L")
        print(f"S_NH4 : {fs.effluent_SR.iconc['S_NH4']:.3f} mg/L")
        print(f"S_Mg  : {fs.effluent_SR.iconc['S_Mg']:.3f}  mg/L")

        print('\n--- SRD Effluent (after redissolution) ---')
        print(f"S_PO4 : {fs.effluent_SRD.iconc['S_PO4']:.3f} mg/L")
        print(f"S_NH4 : {fs.effluent_SRD.iconc['S_NH4']:.3f} mg/L")
        print(f"S_Mg  : {fs.effluent_SRD.iconc['S_Mg']:.3f}  mg/L")

    print('\n--- MBR Effluent ---')
    print(f"S_PO4 : {fs.effluent_MembT.iconc['S_PO4']:.3f} mg/L")
    print(f"S_NH4 : {fs.effluent_MembT.iconc['S_NH4']:.3f} mg/L")
    print(f"S_NO3 : {fs.effluent_MembT.iconc['S_NO3']:.3f} mg/L")

    qs.PowerUtility.price = 0

    tea1 = qs.TEA(system=sysEL, discount_rate=0.05, lifetime=10)
    el.get_TEA_metrics_breakdown(sysEL, include_breakdown=True)

    GWP        = qs.ImpactIndicator('GlobalWarming', alias='GWP', unit='kg CO2-eq')
    Ecosystems = qs.ImpactIndicator('H_Ecosystems',  alias='Ecosystems', unit='points/cap/yr')
    Health     = qs.ImpactIndicator('H_Health',       alias='Health', unit='points/cap/yr')
    Resources  = qs.ImpactIndicator('H_Resources',    alias='Resources', unit='points/cap/yr')

    lca1 = qs.LCA(system=sysEL, lifetime=10, lifetime_unit='yr',
                  indicators=(GWP, Ecosystems, Health, Resources))
    lca1.show()