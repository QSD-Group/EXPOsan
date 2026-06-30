#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
METAB TEA sweep, 1-stage UASB systems only (PASSIVE gas extraction)

Operating space:
    Q = 0.350 m3/d
    COD = 10 to 50 g/L
    HRT = 3 to 5 d

Outputs (base columns first, then normalization checks, then CAPEX line-items,
then annualized breakdown line-items, then totals):

Base:
    Q_m3_d
    COD_g_L
    HRT_d
    COD_removal_pct
    annualized_cost_USD_per_yr
    system_ID

Normalization checks:
    operating_hours_yr
    COD_inf_g_L
    COD_eff_g_L
    F_inf_m3_day
    F_eff_m3_day
    COD_removed_tonne_yr
    USD_per_tonne_COD_removed

CAPEX line-items (raw installed costs, USD):
    CAPEX__<UnitID>__<installed_cost_name>

Annualized breakdown line-items (USD/yr):
    ANNUAL__<category>
"""

import numpy as np
import pandas as pd
import qsdsan as qs

from exposan.metab import create_system
from exposan.metab.utils import categorize_cashflow


# ==========================
# USER INPUTS
# ==========================
OUTFILE = "METAB_TEA_UASB_1stage_passive_small_HRT.xlsx"

Q_M3_D = 0.350
COD_RANGE = np.linspace(10, 50, 9)   # g/L
HRT_RANGE = np.linspace(0.1, 3, 9)     # d

T_C = 22
DISCOUNT_RATE = 0.05
SYSTEM_LIFETIME = 30

# Simulation settings
T_SPAN = (0, 400)
METHOD = "BDF"

# Debug
DEBUG_FIRST_CASE = False


# ==========================
# Helpers
# ==========================
def crf(r, n):
    """Capital recovery factor."""
    r = float(r)
    n = float(n)
    if n <= 0:
        raise ValueError("Lifetime must be > 0.")
    if r == 0:
        return 1.0 / n
    return r * (1 + r) ** n / ((1 + r) ** n - 1)


def sanitize_for_column(name):
    """Make column names cleaner and more stable."""
    return (
        str(name)
        .replace(" ", "_")
        .replace(" - ", "_")
        .replace("/", "_")
        .replace(",", "")
        .replace("(", "")
        .replace(")", "")
    )


def get_scaled_inf_concs(target_cod_gL):
    """
    Scale default ADM1 influent to hit desired total COD.

    target_cod_gL is in g/L.
    s.inf.COD is in mg/L.

    Returns concentrations in mg/L suitable for:
    set_flow_by_concentration(..., units=("m3/d", "mg/L"))
    """
    tmp = create_system(
        n_stages=1,
        reactor_type="UASB",
        gas_extraction="P",
        Q=Q_M3_D,
        T=T_C,
        discount_rate=DISCOUNT_RATE,
        lifetime=SYSTEM_LIFETIME,
    )

    s = tmp.flowsheet.stream
    cmps = s.inf.components

    baseline_cod_mgL = float(s.inf.COD)
    if baseline_cod_mgL <= 0:
        qs.main_flowsheet.clear()
        raise ValueError(f"Baseline influent COD is non-positive: {baseline_cod_mgL} mg/L")

    target_cod_mgL = float(target_cod_gL) * 1000.0
    scale = target_cod_mgL / baseline_cod_mgL

    baseline = {cid: float(s.inf.iconc[cid]) for cid in cmps.IDs}

    # Scale only COD-contributing components
    scaled = baseline.copy()
    for cid, i_cod in zip(cmps.IDs, cmps.i_COD):
        if i_cod > 0:
            scaled[cid] = baseline[cid] * scale

    # Keep water unchanged
    if "H2O" in scaled:
        scaled["H2O"] = baseline["H2O"]

    # Scale nitrogen source with loading
    if "S_IN" in scaled:
        scaled["S_IN"] = baseline["S_IN"] * scale

    qs.main_flowsheet.clear()
    return scaled


def set_total_HRT_uasb(sys, tau_days):
    """
    Set UASB liquid volume from total HRT.

    V_liq = Q * HRT
    where:
        Q is from stream volumetric flow in m3/hr
        HRT is in d
        V_liq is in m3
    """
    u = sys.flowsheet.unit
    s = sys.flowsheet.stream

    V = float(s.inf.F_vol) * 24.0 * float(tau_days)
    u.R1.V_liq = V
    u.R1.V_gas = 0.1 * V

def extract_capex_line_items(sys):
    """
    Flat dict of installed cost line-items (raw USD).

    Keys: CAPEX__<UnitID>__<installed_cost_name>
    """
    out = {}
    for u in sys.cost_units:
        for name, cost in u.installed_costs.items():
            out[f"CAPEX__{u.ID}__{sanitize_for_column(name)}"] = float(cost)
    return out


# ==========================
# Single evaluation
# ==========================
def run_case(cod_gL, hrt_d, do_debug=False):

    sys = create_system(
        n_stages=1,
        reactor_type="UASB",
        gas_extraction="P",
        Q=Q_M3_D,
        T=T_C,
        discount_rate=DISCOUNT_RATE,
        lifetime=SYSTEM_LIFETIME,
    )

    # Set influent concentrations
    inf_concs = get_scaled_inf_concs(cod_gL)
    sys.flowsheet.stream.inf.set_flow_by_concentration(
        Q_M3_D,
        concentrations=inf_concs,
        units=("m3/d", "mg/L"),
    )

    # Set HRT
    set_total_HRT_uasb(sys, hrt_d)

    # Simulate
    sys.simulate(
        state_reset_hook="reset_cache",
        method=METHOD,
        t_span=T_SPAN,
    )

    if do_debug:
        print("\n--- Simulation complete. Entering debug mode ---")
        import pdb
        pdb.set_trace()

    s = sys.flowsheet.stream

    # Try common effluent names safely
    eff = None
    for eff_name in ("eff_dg", "eff", "effluent"):
        if hasattr(s, eff_name):
            eff = getattr(s, eff_name)
            break
    if eff is None:
        qs.main_flowsheet.clear()
        raise AttributeError("Could not find effluent stream. Checked: eff_dg, eff, effluent")

    COD_inf_mgL = float(s.inf.COD)
    COD_eff_mgL = float(eff.COD)

    F_inf_m3_hr = float(s.inf.F_vol)
    F_eff_m3_hr = float(eff.F_vol)

    # mg/L = g/m3, so kg/m3 = mg/L / 1000
    COD_inf_kg_m3 = COD_inf_mgL / 1000.0
    COD_eff_kg_m3 = COD_eff_mgL / 1000.0

    COD_in_kg_hr = COD_inf_kg_m3 * F_inf_m3_hr
    COD_out_kg_hr = COD_eff_kg_m3 * F_eff_m3_hr
    COD_removed_kg_hr = max(COD_in_kg_hr - COD_out_kg_hr, 0.0)

    operating_hours = float(sys.operating_hours)
    COD_removed_tonne_yr = (COD_removed_kg_hr * operating_hours) / 1000.0

    annualized_cost_USD_per_yr = -1.0 * float(sys.TEA.annualized_NPV)

    USD_per_tonne_COD_removed = (
        annualized_cost_USD_per_yr / COD_removed_tonne_yr
        if COD_removed_tonne_yr > 0 else np.nan
    )

    COD_removal_pct = (
        (1.0 - COD_eff_mgL / COD_inf_mgL) * 100.0
        if COD_inf_mgL > 0 else np.nan
    )

    results = {
        "system_ID": sys.ID,
        "Q_m3_d": float(Q_M3_D),
        "COD_g_L": float(cod_gL),
        "HRT_d": float(hrt_d),
        "COD_removal_pct": float(COD_removal_pct),
        "annualized_cost_USD_per_yr": annualized_cost_USD_per_yr,
        "operating_hours_yr": operating_hours,
        "COD_inf_g_L": COD_inf_mgL / 1000.0,
        "COD_eff_g_L": COD_eff_mgL / 1000.0,
        "F_inf_m3_day": F_inf_m3_hr * 24.0,
        "F_eff_m3_day": F_eff_m3_hr * 24.0,
        "COD_removed_tonne_yr": COD_removed_tonne_yr,
        "USD_per_tonne_COD_removed": USD_per_tonne_COD_removed,
    }

    # Annualized breakdown from PV categories
    cf_pv = categorize_cashflow(sys.TEA)
    r = float(getattr(sys.TEA, "discount_rate", DISCOUNT_RATE))
    n = float(getattr(sys.TEA, "lifetime", SYSTEM_LIFETIME))
    CRF = crf(r, n)

    for k_pv, pv_val in cf_pv.items():
        if k_pv == "cnpv":
            continue
        k_ann = f"ANNUAL__{sanitize_for_column(k_pv)}"
        pv_val = float(pv_val)
        results[k_ann] = pv_val * CRF if np.isfinite(pv_val) else np.nan

    # Raw CAPEX line items
    results.update(extract_capex_line_items(sys))

    qs.main_flowsheet.clear()
    return results


# ==========================
# Run full sweep
# ==========================
def main():
    rows = []

    total = len(COD_RANGE) * len(HRT_RANGE)
    counter = 0
    debug_fired = False

    for cod in COD_RANGE:
        for hrt in HRT_RANGE:
            counter += 1
            print(f"{counter}/{total} | UASB-1stage-PASSIVE | COD={cod:.1f} g/L | HRT={hrt:.2f} d")

            do_debug = DEBUG_FIRST_CASE and (not debug_fired)
            rows.append(run_case(cod, hrt, do_debug=do_debug))
            if do_debug:
                debug_fired = True

    df = pd.DataFrame(rows)

    base_cols = [
        "system_ID",
        "Q_m3_d",
        "COD_g_L",
        "HRT_d",
        "COD_removal_pct",
        "annualized_cost_USD_per_yr",
    ]

    norm_cols = [
        "operating_hours_yr",
        "COD_inf_g_L",
        "COD_eff_g_L",
        "F_inf_m3_day",
        "F_eff_m3_day",
        "COD_removed_tonne_yr",
        "USD_per_tonne_COD_removed",
    ]

    annual_breakdown_cols = sorted(
        [c for c in df.columns if c.startswith("ANNUAL__") and c != "ANNUAL__total_USD_per_yr"]
    )
    total_cols = [c for c in ["ANNUAL__total_USD_per_yr"] if c in df.columns]
    capex_cols = sorted([c for c in df.columns if c.startswith("CAPEX__")])

    ordered = [
        c for c in (base_cols + norm_cols + annual_breakdown_cols + total_cols + capex_cols)
        if c in df.columns
    ]
    remaining = [c for c in df.columns if c not in ordered]

    df = df.loc[:, ordered + remaining]
    df.to_excel(OUTFILE, index=False)

    print("\nSaved:", OUTFILE)


if __name__ == "__main__":
    main()