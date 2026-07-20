# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 14:56:05 2026

@author: aliah
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
import flexsolve as flx
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import qmc
import glob
import qsdsan as qs
from chaospy import distributions as shape
import concurrent.futures as cf
import multiprocessing as mp
from joblib import load
from shap import TreeExplainer
import time
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score

# ---- Your system imports (EXPOsan biobinder)
from exposan.biobinder_ml import (
    HTL_yields,
    feedstock_composition,
    tea_kwargs,
    central_dry_flowrate as default_central,
)

# ---- Your RF yield prediction utilities (EXPOsan biobinder_ml)
from exposan.biobinder_ml._feedstocks import get_feedstock_composition
from exposan.biobinder_ml._rf_htl_predictor import predict_htl_yields_from_rf
from exposan.biobinder_ml._ml_features import FeatureDefaults
from exposan.biobinder_ml.Dist_flex import create_system
# from exposan.biobinder_ml._irr_brent import solve_IRR_brent
import flexsolve as flx
from biosteam import TEA

RUN_MODE = "full" #full/shap_only
# ============================================================
# CONFIG
# ============================================================
OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M")

FEEDSTOCKS = ["fog", "sludge", "manure", "green", "food",] #

CONFIGS = [
    {"name": "CHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central,
                                            "decentralized_HTL": False, "decentralized_upgrading": False,
                                            "skip_EC": True, "generate_H2": False, "EC_config": None}},
    {"name": "DHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central,
                                            "decentralized_HTL": True, "decentralized_upgrading": False,
                                            "skip_EC": True, "generate_H2": False, "EC_config": None}},
]

N_SAMPLES = 500
SEED = 42

PREFER = "max_top_ratio"
RF_MODEL_PATH = os.path.join(OUTPUT_DIR, "rf_yield_model.joblib")

USE_FEEDSTOCK_COMPOSITION_UNCERTAINTY = True
USE_PROCESS_KNOB_UNCERTAINTY = True
USE_RF_YIELD_NOISE = False
USE_SYSTEM_PARAM_UNCERTAINTY = True

# ---- Targets in surrogate (IRR in % + GWP)
TARGETS = ["IRR_pct", "GWP"]
CAPS = {
        # "MSP": 50,
    "IRR_pct": 50.0,   # $/kg cap magnitude (±20) – you can also use [0, 20] if MSP is non-negative
    "GWP": 100.0,   # kgCO2e/kg cap magnitude (tune to your model scale)
    
}

# # ---- Cutoff ratio uncertainty
# CUTOFF_FRACS_BASE = [0.03, 0.61, 0.36]
# CUTOFF_LIGHT_FIXED = float(CUTOFF_FRACS_BASE[0])
# CUTOFF_RATIO_BASE = float(CUTOFF_FRACS_BASE[1] / CUTOFF_FRACS_BASE[2])  


# ============================================================
# Helpers
# ============================================================
def safe_float(x):
    try:
        return float(x)
    except Exception:
        return np.nan
    try:
        return float(x)
    except Exception:
        return np.nan

def clip_symmetric(x, cap):
    if cap is None:
        return x
    return np.clip(x, -cap, cap)

def normalize_dict(d, keys):
    s = sum(float(d[k]) for k in keys)
    if s <= 0:
        return d
    for k in keys:
        d[k] = float(d[k]) / s
    return d

def _lhs_samples(n, d, seed=42):
    sampler = qmc.LatinHypercube(d=d, seed=seed)
    return sampler.random(n=n)

def _ppf(dist, u):
    return float(dist.ppf(u))


def find_latest_meta(feedstock_id, config_name, prefer, output_dir="results"):
    pattern = os.path.join(
        output_dir,
        f"META_SHAP3_{feedstock_id}_{config_name}_prefer-{prefer}_*.xlsx"
    )
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No META files found for {feedstock_id} | {config_name}")
    return files[-1]

import numpy as np

def add_tea_capex_checks(out, sys, tea, top_n=5, min_cost=1.0):
    """
    Mirrors your ALL COST CHECKS sheet:
      - ISBL installed cost (tea.ISBL_installed_equipment_cost)
      - OSBL installed cost (tea.OSBL_installed_equipment_cost)
      - tea.TDC (your Total Direct Costs)
      - tea.FCI (your Fixed Capital Investment)
    Plus a compact ISBL top-unit list (one cell) to diagnose why ISBL jumps.
    """

    # --- Direct-cost pieces (what your sheet shows) ---
    out["CAPEX_ISBL_installed_$"] = float(getattr(tea, "ISBL_installed_equipment_cost", np.nan))
    out["CAPEX_OSBL_installed_$"] = float(getattr(tea, "OSBL_installed_equipment_cost", np.nan))

    # In BioSTEAM TEA, TDC is "Total depreciable capital",
    # BUT in your HTL TEA subclass you are using it as "Total Direct Costs (TDC)" in the sheet.
    out["CAPEX_TDC_$"] = float(getattr(tea, "TDC", np.nan))
    out["CAPEX_FCI_$"] = float(getattr(tea, "FCI", np.nan))
    out["CAPEX_TCI_$"] = float(getattr(tea, "TCI", np.nan))

    # --- Your extra "of ISBL" adders are % of ISBL installed ---
    # (these are fixed % in create_tea kwargs; record them and the implied $)
    isbl = out["CAPEX_ISBL_installed_$"]
    if np.isfinite(isbl):
        out["CAPEX_buildings_$"] = 0.04  * isbl
        out["CAPEX_site_dev_$"]  = 0.10  * isbl
        out["CAPEX_add_piping_$"]= 0.045 * isbl
    else:
        out["CAPEX_buildings_$"] = np.nan
        out["CAPEX_site_dev_$"]  = np.nan
        out["CAPEX_add_piping_$"]= np.nan

    # --- Consistency check: (ISBL + OSBL + buildings + site + piping) - tea.TDC ---
    # This should be ~0 if definitions align.
    if np.isfinite(out["CAPEX_TDC_$"]) and np.isfinite(out["CAPEX_OSBL_installed_$"]) and np.isfinite(isbl):
        implied_tdc = (
            isbl
            + out["CAPEX_OSBL_installed_$"]
            + out["CAPEX_buildings_$"]
            + out["CAPEX_site_dev_$"]
            + out["CAPEX_add_piping_$"]
        )
        out["CAPEX_TDC_check_$"] = float(implied_tdc - out["CAPEX_TDC_$"])
    else:
        out["CAPEX_TDC_check_$"] = np.nan

    # --- Compact diagnostic: which ISBL unit(s) drove ISBL installed cost ---
    OSBL = set(getattr(tea, "OSBL_units", []) or [])
    rows = []
    for u in sys.units:
        if u in OSBL:
            continue
        val = float(getattr(u, "installed_cost", 0.0) or 0.0)
        if val >= min_cost:
            rows.append((u.ID, val))
    rows.sort(key=lambda x: x[1], reverse=True)

    top = rows[:top_n]
    rest = rows[top_n:]
    out["CAPEX_ISBL_top_units"] = "; ".join([f"{uid}:{val:,.0f}" for uid, val in top]) if top else ""
    out["CAPEX_ISBL_rest_$"] = float(sum(v for _, v in rest)) if rest else 0.0

    # Useful: fraction of ISBL from the biggest unit (often CrudeHeavyDis)
    out["CAPEX_ISBL_top1_share"] = float(top[0][1] / isbl) if top and np.isfinite(isbl) and isbl > 0 else np.nan

def hx_duty_kJ_per_hr(hx):
    # HXutility: duty is energy change across the exchanger (outs[0] vs ins[0])
    try:
        return float(hx.outs[0].H - hx.ins[0].H)  # kJ/hr (BioSTEAM convention)
    except Exception:
        return np.nan

def safe_get_design(col, key):
    try:
        v = col.design_results.get(key, np.nan)
        # Some design entries can be strings like '100+'
        return float(v) if np.isfinite(v) else np.nan
    except Exception:
        try:
            return float(col.design_results[key])
        except Exception:
            return np.nan

def add_column_diagnostics(out, col, prefix="CHD_"):
    # --- design numbers ---
    out[prefix+"R"]     = safe_get_design(col, "Reflux")
    out[prefix+"Rmin"]  = safe_get_design(col, "Minimum reflux")

    # These can be '100+' or '?', store as string-safe too:
    out[prefix+"N_theoretical_raw"] = col.design_results.get("Theoretical stages", None)
    out[prefix+"N_feed_raw"]        = col.design_results.get("Theoretical feed stage", None)

    # If you want numeric versions when possible:
    try:
        out[prefix+"N_theoretical"] = float(col.design_results.get("Theoretical stages", np.nan))
    except Exception:
        out[prefix+"N_theoretical"] = np.nan

    # --- duties (kJ/hr) ---
    out[prefix+"Q_cond_kJhr"] = hx_duty_kJ_per_hr(col.condenser)
    out[prefix+"Q_reb_kJhr"]  = hx_duty_kJ_per_hr(col.reboiler)

    # --- utility costs (USD/hr) from each HXutility if available ---
    # These are more interpretable than SYS totals.
    try:
        out[prefix+"cond_utility_cost_USDhr"] = float(col.condenser.utility_cost)
    except Exception:
        out[prefix+"cond_utility_cost_USDhr"] = np.nan

    try:
        out[prefix+"reb_utility_cost_USDhr"]  = float(col.reboiler.utility_cost)
    except Exception:
        out[prefix+"reb_utility_cost_USDhr"]  = np.nan

    # --- mass flows (helps find cases with near-zero denom causing nonsense) ---
    try:
        out[prefix+"F_feed_kgph"] = float(col.ins[0].F_mass)
        out[prefix+"F_top_kgph"]  = float(col.outs[0].F_mass)
        out[prefix+"F_bot_kgph"]  = float(col.outs[1].F_mass)
    except Exception:
        out[prefix+"F_feed_kgph"] = out[prefix+"F_top_kgph"] = out[prefix+"F_bot_kgph"] = np.nan

def _NPV_at_r(r, CF, dur):
    r = float(r)
    if (not np.isfinite(r)) or (r <= -0.999999):
        return np.nan
    return (CF / (1.0 + r)**dur).sum()

def solve_IRR_aitken_limit(tea, maxiter=200, xtol=1e-6, ytol=10.0,  npv_tol_abs=1e3,):
    # Grab arrays
    CF = tea.cashflow_array
    dur = tea._get_duration_array()

    # 1) cashflows must be finite
    if (CF is None) or (dur is None) or (not np.isfinite(CF).all()) or (not np.isfinite(dur).all()):
        return np.nan, "CF_or_dur_nonfinite"

    args = (CF, dur)

    # 2) initial guess
    IRR_guess = float(getattr(tea, "IRR", 0.10) or 0.10)
    if (not np.isfinite(IRR_guess)) or (IRR_guess <= -0.999999):
        IRR_guess = 0.10

    x0 = IRR_guess
    x1 = 1.0001 * IRR_guess + 1e-3

    # 3) basic "can we evaluate f" sanity at start points
    y0 = _NPV_at_r(x0, *args)
    y1 = _NPV_at_r(x1, *args)
    if (not np.isfinite(y0)) or (not np.isfinite(y1)):
        return np.nan, "endpoint_nonfinite"

    # 4) run solver
    try:
        irr = flx.aitken_secant(
            _NPV_at_r,
            x0, x1,
            xtol=xtol, ytol=ytol, maxiter=maxiter,
            args=args,
            checkiter=False
        )
    except Exception as e:
        return np.nan, f"solver_exception:{type(e).__name__}"

    # 5) validate result
    if (not np.isfinite(irr)) or (irr <= -0.999999):
        return np.nan, "solver_returned_invalid"

    #require that NPV is ~0 at solution and reject fakes (since checkiter=False)
    y = _NPV_at_r(irr, *args)
    if not np.isfinite(y):
        return np.nan, "npv_at_solution_nonfinite"
    # ytol=10 USD
    if abs(float(y)) > float(npv_tol_abs):
        return np.nan, f"npv_not_close:{float(y):.3g}"

    return float(irr), "ok"


# def solve_IRR_bounded(tea, IRR_min=-0.9, IRR_max=5.0, maxiter=200):
#     CF = tea.cashflow_array
#     dur = tea._get_duration_array()

#     # Hard fail if cashflow is not finite
#     if not np.isfinite(CF).all():
#         return np.nan

#     def f(r):
#         # Disallow r <= -1 to avoid division by zero / negative base issues
#         if r <= -0.999999:
#             return np.nan
#         return (CF / (1.0 + r)**dur).sum()

#     # Sample endpoints
#     y_lo = f(IRR_min)
#     y_hi = f(IRR_max)
#     if (not np.isfinite(y_lo)) or (not np.isfinite(y_hi)):
#         return np.nan

#     # If no sign change, try to find a bracket by scanning
#     if np.sign(y_lo) == np.sign(y_hi):
#         grid = np.array([-0.5, 0.0, 0.1, 0.3, 0.6, 1.0, 2.0, 3.0, 5.0])
#         ys = np.array([f(r) for r in grid])
#         ok = np.isfinite(ys)
#         grid, ys = grid[ok], ys[ok]
#         # find first adjacent sign change
#         idx = np.where(np.sign(ys[:-1]) != np.sign(ys[1:]))[0]
#         if idx.size == 0:
#             return np.nan
#         a, b = grid[idx[0]], grid[idx[0] + 1]
#     else:
#         a, b = IRR_min, IRR_max

#     # Bracketing solver (robust)
#     try:
#         return flx.IQ_interpolation(f, a, b, xtol=1e-6, ytol=10.0, maxiter=maxiter, checkiter=False)
#     except Exception:
#         # last resort bisection
#         try:
#             return flx.bisect(f, a, b, xtol=1e-6, ytol=10.0, maxiter=maxiter, checkiter=False)
#         except Exception:
#             return np.nan

# def mark_runaway(out, prefix="CrudeHeavyDis_"):
#     R  = out.get(prefix+"R", np.nan)
#     Qc = out.get(prefix+"Q_cond_kJhr", np.nan)
#     Qr = out.get(prefix+"Q_reb_kJhr", np.nan)

#     # Thresholds: pick conservative; you can tighten after you look at distributions
#     runaway = False
#     if np.isfinite(R) and R > 50: runaway = True            # reflux ratio insane
#     if np.isfinite(Qc) and abs(Qc) > 1e10: runaway = True   # kJ/hr absurd
#     if np.isfinite(Qr) and abs(Qr) > 1e10: runaway = True

#     out["RUNAWAY_sep"] = bool(runaway)

# # after add_column_diagnostics(...)
# mark_runaway(out, prefix="CrudeHeavyDis_")

# def compute_cutoff_fracs_from_ratio(ratio_m_over_h, light=CUTOFF_LIGHT_FIXED):
#     """Reconstruct cutoff_fracs=[light,m,h] given ratio m/h, holding light fixed."""
#     r = float(ratio_m_over_h)
#     if not np.isfinite(r) or r <= 0:
#         raise ValueError(f"Invalid cutoff ratio (m/h): {ratio_m_over_h}")
#     rem = 1.0 - float(light)
#     h = rem / (1.0 + r)
#     m = rem - h
#     return [float(light), float(m), float(h)]


# def apply_cutoff_ratio(sys, ratio_m_over_h):
#     """Apply sampled cutoff ratio to BiocrudeSplitter._cutoff_fracs."""
#     splitter = sys.flowsheet.unit.BiocrudeSplitter
#     new_fracs = compute_cutoff_fracs_from_ratio(ratio_m_over_h, light=CUTOFF_LIGHT_FIXED)
#     splitter._cutoff_fracs = new_fracs

#     # best-effort refresh hooks
#     for maybe in ("_update", "_update_splits", "update", "_reset_cache"):
#         if hasattr(splitter, maybe) and callable(getattr(splitter, maybe)):
#             try:
#                 getattr(splitter, maybe)()
#                 break
#             except Exception:
#                 pass

#     return new_fracs


# def compute_irr_percent(tea, products):
#     """
#     IRR metric as %:
#       IRR_pct = tea.solve_IRR() * 100
#     """
#     if not hasattr(tea, "solve_IRR"):
#         raise AttributeError("TEA object does not have solve_IRR().")
#     irr = tea.solve_IRR(products)
#     return float(irr) * 100.0


def compute_gwp_per_kg_biobinder(sys):
    """
    GWP in kg CO2e/kg biobinder (operation_only, annual, allocated impacts).
    """
    lca = sys.LCA
    biobinder = sys.flowsheet.stream.biobinder
    impacts = lca.get_allocated_impacts(streams=(biobinder,), operation_only=True, annual=True)
    gwp = impacts["GWP"] / (biobinder.F_mass * lca.system.operating_hours)
    return float(gwp)


# ============================================================
# Define uncertain inputs (distributions)
# ============================================================
def build_uncertainty_spec(sys, config_name):
    spec = []

    # (A) Feedstock composition uncertainty
    if USE_FEEDSTOCK_COMPOSITION_UNCERTAINTY:
        spec += [
            {"name": "FS_factor_Water", "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Carb",  "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Prot",  "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Lipid", "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Ash",   "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
        ]

    # (B) HTL process knobs
    if USE_PROCESS_KNOB_UNCERTAINTY:
        spec += [
            {"name": "Temperature (C)",          "dist": shape.Uniform(260, 320), "group": "process"},
            {"name": "Residence Time",           "dist": shape.Uniform(10, 30),   "group": "process"},
            {"name": "Solid content (w/w) %",    "dist": shape.Uniform(15, 25),   "group": "process"},
            {"name": "Pre-processing",           "dist": shape.DiscreteUniform(0, 1), "group": "process_cat"},
            {"name": "Catalyst",                 "dist": shape.DiscreteUniform(1, 3), "group": "process_cat"},
            {"name": "Reactor Type",             "dist": shape.DiscreteUniform(0, 2), "group": "process_cat"},
            {"name": "Solvent",                  "dist": shape.DiscreteUniform(1, 3), "group": "process_cat"},
        ]

    # (B2) Cutoff ratio uncertainty (medium/heavy)
    # spec += [
    #     {"name": "Cutoff_ratio_M_over_H",
    #      "dist": shape.Uniform(0.5 * CUTOFF_RATIO_BASE, 1.5 * CUTOFF_RATIO_BASE),
    #      "group": "cutoff"}
    # ]

    # (C) System/TEA uncertainties (NO IRR sampling!)
    if USE_SYSTEM_PARAM_UNCERTAINTY:
        spec += [
            {"name": "Biofuel_price",          "dist": shape.Uniform(1.02, 1.25),            "group": "economics"},
            # {"name": "N_Fertilizer_Price",     "dist": shape.Triangle(0.85, 0.9, 1.19),      "group": "economics"},
            {"name": "Natural_Gas_Price",      "dist": shape.Normal(0.14, 0.0211),           "group": "economics"},
            {"name": "Electricity_Price",      "dist": shape.Triangle(0.03, 0.074, 0.1059),  "group": "economics"},
            {"name": "Uptime_Ratio",           "dist": shape.Uniform(0.8, 1.0),              "group": "operations"},
            {"name": "income_tax",             "dist": shape.Uniform(0.15, 0.35),            "group": "policy_finance"},
            # {"name": "IRR",                    "dist": shape.Uniform(0.05, 0.15),            "group": "policy_finance"}
        ]

        if "CHCU" in config_name:
            spec += [{"name": "Feedstock_Transportation_Cost", "dist": shape.Uniform(0.048, 0.059), "group": "logistics"}]
        if "DHCU" in config_name:
            spec += [{"name": "Biocrude_Transportation_Cost",  "dist": shape.Uniform(0.089, 0.109), "group": "logistics"}]

    # (D) RF yield noise (optional)
    # if USE_RF_YIELD_NOISE:
    #     spec += [
    #         {"name": "RF_noise_biocrude", "dist": shape.Normal(0.0, 0.02), "group": "yield_noise"},
    #         {"name": "RF_noise_aqueous",  "dist": shape.Normal(0.0, 0.02), "group": "yield_noise"},
    #         {"name": "RF_noise_gas",      "dist": shape.Normal(0.0, 0.02), "group": "yield_noise"},
    #         {"name": "RF_noise_char",     "dist": shape.Normal(0.0, 0.02), "group": "yield_noise"},
    #     ]

    return spec


def generate_lhs_samples(spec, n_samples, seed=42):
    d = len(spec)
    U = _lhs_samples(n_samples, d, seed=seed)
    samples = []
    for i in range(n_samples):
        row = {}
        for j, s in enumerate(spec):
            v = _ppf(s["dist"], U[i, j])
            if "DiscreteUniform" in str(type(s["dist"])) or s["name"] in ["Catalyst", "Reactor Type", "Solvent", "Pre-processing"]:
                v = int(np.round(v))
            row[s["name"]] = float(v)
        samples.append(row)
    return [s["name"] for s in spec], samples


# ============================================================
# Apply sampled scenario to system
# ============================================================
def apply_sample_to_system(sample, sys, tea, config_name, rf_model=None, base_wet=None):
    wet = dict(base_wet) if base_wet is not None else {}

    # (1) Feedstock composition
    if USE_FEEDSTOCK_COMPOSITION_UNCERTAINTY and base_wet is not None:
        comp_map = {
            "Water": "FS_factor_Water",
            "Carbohydrates": "FS_factor_Carb",
            "Proteins": "FS_factor_Prot",
            "Lipids": "FS_factor_Lipid",
            "Ash": "FS_factor_Ash",
        }
        new = {}
        for k, fac_key in comp_map.items():
            fac = float(sample.get(fac_key, 1.0))
            new[k] = max(0.0, wet[k] * fac)

        s = sum(new.values())
        if s > 0:
            for k in new:
                new[k] /= s
        wet = new

    for k in list(feedstock_composition.keys()):
        feedstock_composition[k] = wet.get(k, feedstock_composition[k])

    # ------------------------------------------------------------
    # (2) Process knobs (for RF) + ALSO push into flowsheet (HTL + Conditioning)
    # ------------------------------------------------------------
    process = {}
    T_C = None
    tau_min = None
    solid_pct = None

    if USE_PROCESS_KNOB_UNCERTAINTY:
        T_C = float(sample.get("Temperature (C)", 280.0))
        tau_min = float(sample.get("Residence Time", 15.0))          # minutes
        solid_pct = float(sample.get("Solid content (w/w) %", 20.0))  # %

        process = {
            "Temperature (C)": T_C,
            "Residence Time": tau_min,
            "Solid content (w/w) %": solid_pct,
            "Pre-processing": int(sample.get("Pre-processing", FeatureDefaults().pre_processing)),
            "Catalyst": int(sample.get("Catalyst", FeatureDefaults().catalyst)),
            "Reactor Type": int(sample.get("Reactor Type", FeatureDefaults().reactor_type)),
            "Solvent": int(sample.get("Solvent", FeatureDefaults().solvent)),
            "Reactor Volume (mL)": float(FeatureDefaults().reactor_volume_ml),
        }

        # ---- PUSH INTO FLOWSHEET ----
        # HTL operating conditions
        if hasattr(sys.flowsheet.unit, "HTL"):
            sys.flowsheet.unit.HTL.T = T_C + 273.15      # K
            sys.flowsheet.unit.HTL.tau = tau_min / 60.0  # hours (your HTL_kwargs uses tau=15/60)

        # Conditioning solid loading
        # (you showed FeedstockCond.target_HTL_solid_loading exists)
        if hasattr(sys.flowsheet.unit, "FeedstockCond"):
            sl_frac = solid_pct / 100.0                  # fraction 0.15–0.25
            sys.flowsheet.unit.FeedstockCond.target_HTL_solid_loading = sl_frac

        # Optional: lock values in case later specs reset them
        # (remove if you don’t need it)
        try:
            htl = sys.flowsheet.unit.HTL
            print("AFTER-SET HTL.T(C):", htl.T - 273.15)
            print("AFTER-SET HTL.tau(min):", htl.tau * 60)

            cond = sys.flowsheet.unit.FeedstockCond
            print("AFTER-SET solid loading (%):", cond.target_HTL_solid_loading * 100)
            
            def _lock_operating_conditions():
                htl.T = T_C + 273.15
                htl.tau = tau_min / 60.0
                cond.target_HTL_solid_loading = (solid_pct / 100.0)

            cond.add_specification(_lock_operating_conditions)
            cond.run_after_specifications = True
        except Exception:
            pass

    # (3) RF yields + noise
    yields_used = None
    if rf_model is not None and base_wet is not None:
        dry_frac = 1.0 - wet["Water"]
        if dry_frac <= 0:
            raise ValueError("Non-physical water fraction in sampled feedstock composition.")

        carb_wt    = wet["Carbohydrates"] / dry_frac * 100.0
        protein_wt = wet["Proteins"]      / dry_frac * 100.0
        lipids_wt  = wet["Lipids"]        / dry_frac * 100.0
        ash_wt     = wet["Ash"]           / dry_frac * 100.0

        if not process:
            fd = FeatureDefaults()
            process = {
                "Temperature (C)": 280.0,
                "Residence Time": 15.0,
                "Solid content (w/w) %": 20.0,
                "Pre-processing": fd.pre_processing,
                "Catalyst": fd.catalyst,
                "Reactor Type": fd.reactor_type,
                "Solvent": fd.solvent,
                "Reactor Volume (mL)": fd.reactor_volume_ml,
            }

        y = predict_htl_yields_from_rf(
            rf_model=rf_model,
            carb_wt=carb_wt,
            protein_wt=protein_wt,
            lipids_wt=lipids_wt,
            ash_wt=ash_wt,
            process=process,
        )

        if USE_RF_YIELD_NOISE:
            y_noisy = {
                "biocrude": max(0.0, y["biocrude"] + float(sample.get("RF_noise_biocrude", 0.0))),
                "aqueous":  max(0.0, y["aqueous"]  + float(sample.get("RF_noise_aqueous", 0.0))),
                "gas":      max(0.0, y["gas"]      + float(sample.get("RF_noise_gas", 0.0))),
                "char":     max(0.0, y["char"]     + float(sample.get("RF_noise_char", 0.0))),
            }
            y_noisy = normalize_dict(y_noisy, ["biocrude", "aqueous", "gas", "char"])
            yields_used = y_noisy
        else:
            yields_used = y
        
        sys.flowsheet.unit.HTL.dw_yields.clear()
        sys.flowsheet.unit.HTL.dw_yields.update(yields_used)
        print("RF process used:", process["Temperature (C)"], process["Residence Time"], process["Solid content (w/w) %"])
        print("FINAL HTL.dw_yields (used):", dict(sys.flowsheet.unit.HTL.dw_yields))

    # (4) system params
    if USE_SYSTEM_PARAM_UNCERTAINTY:
        if "Biofuel_price" in sample:
            sys.flowsheet.stream.biofuel.price = float(sample["Biofuel_price"])
        # if "N_Fertilizer_Price" in sample:
        #     sys.flowsheet.stream.recovered_N.price = float(sample["N_Fertilizer_Price"])
        if "Natural_Gas_Price" in sample:
            sys.flowsheet.stream.natural_gas.price = float(sample["Natural_Gas_Price"])
        if "Electricity_Price" in sample:
            qs.PowerUtility.price = float(sample["Electricity_Price"])
        if "Uptime_Ratio" in sample:
            sys.flowsheet.unit.HTL.Uptime_ratio = float(sample["Uptime_Ratio"])

        if "Feedstock_Transportation_Cost" in sample:
            sys.flowsheet.unit.FeedstockTrans.cost = float(sample["Feedstock_Transportation_Cost"])
        if "Biocrude_Transportation_Cost" in sample:
            sys.flowsheet.unit.BiocrudeTrans.cost = float(sample["Biocrude_Transportation_Cost"])

        tea.IRR = tea_kwargs["IRR"]

        if "income_tax" in sample:
            tea.income_tax = float(sample["income_tax"])
        else:
            tea.income_tax = tea_kwargs["income_tax"]

    # (5) cutoff ratio
    cutoff_applied = None
    # if "Cutoff_ratio_M_over_H" in sample:
    #     cutoff_applied = apply_cutoff_ratio(sys, sample["Cutoff_ratio_M_over_H"])

    return {
        "process": process,
        "yields_used": yields_used,
        "wet_comp_used": wet,
        "cutoff_fracs_applied": cutoff_applied,
    }


# ============================================================
# Run one full simulation with one sample
# ============================================================
def simulate_one(sample, config_name, config_kwargs, feedstock_id, rf_model):
    out = dict(sample)

    print(
        "SAMPLED:",
        sample.get("Temperature (C)"),
        sample.get("Residence Time"),
        sample.get("Solid content (w/w) %"),
    )

    sys = create_system(feedstock_id=feedstock_id, **config_kwargs)
    tea = sys.TEA

    base_wet = get_feedstock_composition(feedstock_id)

    details = apply_sample_to_system(
        sample=sample,
        sys=sys,
        tea=tea,
        config_name=config_name,
        rf_model=rf_model,
        base_wet=base_wet,
    )

    try:
        sys.simulate()
        hours = sys.operating_hours

        out["Feedstock_price_$/tonne"] = float(sys.flowsheet.stream.scaled_feedstock.price * 1000)

        col = sys.flowsheet.unit.CrudeHeavyDis
        feed = col.ins[0]
        out["COL_feed_F_mass"] = float(feed.F_mass)
        out["COL_feed_T"] = float(feed.T)
        out["COL_feed_H"] = float(feed.H)
        nz = [ID for ID, m in zip(feed.chemicals.IDs, feed.mol) if m > 0]
        out["COL_feed_nonzero_n"] = len(nz)

        flash = sys.flowsheet.unit.BiofuelFlash
        best = getattr(col, "_runtime_best", None)
        out["SEP_policy_top_ratio"] = float(best["top_ratio"]) if best else np.nan
        out["SEP_column_top_ratio"] = float(col.outs[0].F_mass / col.F_mass_out) if col.F_mass_out else np.nan
        out["SEP_flash_vapor_loss_frac_of_column_top"] = float(flash.outs[0].F_mass / col.outs[0].F_mass) if col.outs[0].F_mass else np.nan
        out["SEP_biofuel_liquid_mass"] = float(flash.outs[1].F_mass)

        out.update(
            {
                "SEP_LHK_light": col.LHK[0] if isinstance(col.LHK, tuple) else col.LHK,
                "SEP_LHK_heavy": col.LHK[1] if isinstance(col.LHK, tuple) else None,
                "SEP_Lr": float(col.Lr),
                "SEP_Hr": float(col.Hr),
            }
        )

        D = col.design_results
        out["CHD_N_theoretical"] = float(D.get("Theoretical stages", np.nan))
        out["CHD_R"] = float(D.get("Reflux", np.nan))
        out["CHD_Rmin"] = float(D.get("Minimum reflux", np.nan))
        out["CHD_Diameter_ft"] = float(D.get("Rectifier diameter", D.get("Diameter", np.nan)))
        out["CHD_Height_ft"] = float(D.get("Rectifier height", D.get("Height", np.nan)))
        out["CHD_twall_in"] = float(D.get("Rectifier wall thickness", D.get("Wall thickness", np.nan)))
        out["CHD_weight_lb"] = float(D.get("Rectifier weight", D.get("Weight", np.nan)))

        alpha_mean = col._estimate_mean_volatilities_relative_to_heavy_key()
        LK_index = col._LHK_vle_index[0]
        out["CHD_alpha_LK"] = float(alpha_mean[LK_index])

        # ---------------------------
        # Deep distillation diagnostics (dp/bp, K, alpha, Fenske/Underwood)
        # ---------------------------
        import numpy as _np

        col = sys.flowsheet.unit.CrudeHeavyDis
        feed = col.ins[0]
        Dstream, Bstream = col.outs

        out["COL_feed_P_Pa"] = float(getattr(feed, "P", _np.nan))
        try:
            out["COL_feed_F_mol"] = float(feed.F_mol)
        except Exception:
            out["COL_feed_F_mol"] = _np.nan

        best = getattr(col, "_runtime_best", None)
        out["SEP_policy_found"] = bool(best is not None)
        if best:
            out["SEP_policy_prefer"] = str(best.get("prefer", ""))
            out["SEP_policy_LHK"] = str(best.get("LHK", ""))
            out["SEP_policy_Lr"] = float(best.get("Lr", _np.nan))
            out["SEP_policy_Hr"] = float(best.get("Hr", _np.nan))
            out["SEP_policy_top_ratio"] = float(best.get("top_ratio", _np.nan))
            out["SEP_policy_err"] = float(best.get("err", _np.nan))

        try:
            out["SEP_LHK_light"] = col.LHK[0]
            out["SEP_LHK_heavy"] = col.LHK[1]
        except Exception:
            out["SEP_LHK_light"] = str(getattr(col, "LHK", ""))
            out["SEP_LHK_heavy"] = _np.nan

        out["SEP_Lr_final"] = float(getattr(col, "Lr", _np.nan))
        out["SEP_Hr_final"] = float(getattr(col, "Hr", _np.nan))

        IDs = getattr(col, "_IDs_vle", None)
        LHK_vle = getattr(col, "_LHK_vle_index", None)

        out["CHD_P_Pa"] = float(getattr(col, "P", _np.nan))
        out["CHD_partial_condenser"] = bool(getattr(col, "_partial_condenser", True))
        out["CHD_IDs_vle_n"] = len(IDs) if IDs else 0

        if IDs and LHK_vle is not None:
            LK_vle, HK_vle = int(LHK_vle[0]), int(LHK_vle[1])

            zD = Dstream.get_normalized_mol(IDs)
            zB = Bstream.get_normalized_mol(IDs)
            out["CHD_zD_LK"] = float(zD[LK_vle])
            out["CHD_zD_HK"] = float(zD[HK_vle])
            out["CHD_zB_LK"] = float(zB[LK_vle])
            out["CHD_zB_HK"] = float(zB[HK_vle])

            dew = getattr(col, "_dew_point", None)
            bub = getattr(col, "_bubble_point", None)

            try:
                if getattr(col, "_partial_condenser", True):
                    dp = dew(zD, P=col.P)
                else:
                    dp = bub(zD, P=col.P)
                bp = bub(zB, P=col.P)

                out["CHD_dp_T_K"] = float(dp.T)
                out["CHD_bp_T_K"] = float(bp.T)

                out["CHD_dp_min_x_raw"] = float(_np.min(dp.x))
                out["CHD_bp_min_x_raw"] = float(_np.min(bp.x))

                out["CHD_dp_x_LK"] = float(dp.x[LK_vle])
                out["CHD_dp_y_LK"] = float(dp.y[LK_vle])
                out["CHD_dp_x_HK"] = float(dp.x[HK_vle])
                out["CHD_dp_y_HK"] = float(dp.y[HK_vle])

                out["CHD_bp_x_LK"] = float(bp.x[LK_vle])
                out["CHD_bp_y_LK"] = float(bp.y[LK_vle])
                out["CHD_bp_x_HK"] = float(bp.x[HK_vle])
                out["CHD_bp_y_HK"] = float(bp.y[HK_vle])

                eps = 1e-16
                Kdist = dp.y / _np.maximum(dp.x, eps)
                Kbot = bp.y / _np.maximum(bp.x, eps)

                out["CHD_Kdist_LK"] = float(Kdist[LK_vle])
                out["CHD_Kdist_HK"] = float(Kdist[HK_vle])
                out["CHD_Kbot_LK"] = float(Kbot[LK_vle])
                out["CHD_Kbot_HK"] = float(Kbot[HK_vle])

                out["CHD_alpha_dist_LK_over_HK"] = float(Kdist[LK_vle] / max(Kdist[HK_vle], eps))
                out["CHD_alpha_bot_LK_over_HK"] = float(Kbot[LK_vle] / max(Kbot[HK_vle], eps))

                alpha_mean = col._estimate_mean_volatilities_relative_to_heavy_key()
                out["CHD_alpha_LK_used"] = float(alpha_mean[LK_vle])
                out["CHD_alpha_LK_le_1"] = bool(alpha_mean[LK_vle] <= 1.0)            
                # ratio = (LK_D/HK_D) * (HK_B/LK_B)
                try:
                    LK_i, HK_i = col._LHK_index
                    LK_D = float(Dstream.mol[LK_i])
                    HK_D = float(Dstream.mol[HK_i])
                    LK_B = float(Bstream.mol[LK_i])
                    HK_B = float(Bstream.mol[HK_i])
                    ratio = (LK_D / max(HK_D, 1e-30)) * (HK_B / max(LK_B, 1e-30))
                    out["CHD_Fenske_ratio"] = float(ratio)
                    out["CHD_log10_Fenske_ratio"] = float(_np.log10(max(ratio, 1e-30)))
                    out["CHD_log10_alpha_LK"] = float(_np.log10(max(alpha_mean[LK_vle], 1e-30)))
                except Exception:
                    out["CHD_Fenske_ratio"] = _np.nan
                    out["CHD_log10_Fenske_ratio"] = _np.nan
                    out["CHD_log10_alpha_LK"] = _np.nan
                #Underwood diagnostics (theta + Rm)
                try:
                    q = col.get_feed_quality()
                    z_f = col.ins[0].get_normalized_mol(IDs)
                    theta = col._solve_Underwood_constant(alpha_mean, alpha_mean[LK_vle])
                    z_d = Dstream.get_normalized_mol(IDs)
                    Rm = float((alpha_mean * z_d / (alpha_mean - theta)).sum() - 1.0)
                    out["CHD_q"] = float(q)
                    out["CHD_theta"] = float(theta)
                    out["CHD_Rm_calc"] = float(Rm)
                    out["CHD_Rm_bad"] = bool((not _np.isfinite(Rm)) or (Rm < 0))
                except Exception:
                    out["CHD_q"] = _np.nan
                    out["CHD_theta"] = _np.nan
                    out["CHD_Rm_calc"] = _np.nan
                    out["CHD_Rm_bad"] = _np.nan

            except Exception as _e:
                out["CHD_dp_bp_error"] = str(_e)
        else:
            out["CHD_dp_bp_error"] = "missing _IDs_vle or _LHK_vle_index"

        try:
            Po_psi = col.P * 0.000145078
            out["CHD_Po_psi"] = float(Po_psi)
            out["CHD_is_vacuum"] = bool(Po_psi < 14.68)
        except Exception:
            out["CHD_Po_psi"] = _np.nan
            out["CHD_is_vacuum"] = _np.nan

        # TEA results
        out["TEA_FCI"] = safe_float(tea.FCI)
        out["TEA_FOC"] = safe_float(tea.FOC)
        out["TEA_VOC"] = safe_float(tea.VOC)
        out["TEA_AOC"] = safe_float(tea.AOC)
        out["SYS_sales"] = safe_float(sys.sales)
        out["SYS_material_cost"] = safe_float(sys.material_cost)
        out["SYS_utility_cost"] = safe_float(sys.utility_cost)

        try:
            out["SYS_net_electricity_kWh"] = safe_float(
                sys.get_electricity_consumption() - sys.get_electricity_production()
            )
        except Exception:
            out["SYS_net_electricity_kWh"] = np.nan

        try:
            out["SYS_power_utility_cost"] = safe_float(sys.power_utility.cost * hours)
        except Exception:
            out["SYS_power_utility_cost"] = np.nan

        try:
            out["SYS_heating_duty_MJ_per_hr"] = safe_float(sys.get_heating_duty() / 1000.0)
        except Exception:
            out["SYS_heating_duty_MJ_per_hr"] = np.nan

        try:
            out["SYS_cooling_duty_MJ_per_hr"] = safe_float(sys.get_cooling_duty() / 1000.0)
        except Exception:
            out["SYS_cooling_duty_MJ_per_hr"] = np.nan

        try:
            out["CrudeHeavyDis_utility_cost"] = safe_float(col.utility_cost * hours)
        except Exception:
            out["CrudeHeavyDis_utility_cost"] = np.nan

        try:
            out["CHP_utility_cost"] = safe_float(sys.flowsheet.unit.CHP.utility_cost * hours)
        except Exception:
            out["U_CHP_utility_cost"] = np.nan

        try:
            out["TEA_NPV"] = safe_float(tea.NPV)
        except Exception:
            out["TEA_NPV"] = np.nan

        try:
            out["TEA_net_earnings"] = safe_float(tea.net_earnings)
        except Exception:
            out["TEA_net_earnings"] = np.nan

        biobinder = sys.flowsheet.stream.biobinder
        biofuel = sys.flowsheet.stream.biofuel
        out["Biobinder_MT"] = safe_float(biobinder.F_mass * hours / 1000)
        out["Biofuel_MT"] = safe_float(biofuel.F_mass * hours / 1000)

        sep_top = np.nan
        den = float(biofuel.F_mass + biobinder.F_mass)
        if den > 0:
            sep_top = biofuel.F_mass / den
        out["SEP_biofuel_mass_frac"] = float(sep_top) if np.isfinite(sep_top) else np.nan

        prod_ratio = np.nan
        if biofuel.F_mass and biofuel.F_mass > 0:
            prod_ratio = biobinder.F_mass / biofuel.F_mass

        biobinder.price = 0.1
        irr, irr_reason = solve_IRR_aitken_limit(tea)
        out["IRR_reason"] = irr_reason
        out["IRR_raw"] = safe_float(irr)
        IRR_pct = irr * 100.0 if np.isfinite(irr) else np.nan
        out["IRR_pct"] = safe_float(IRR_pct)
        try:
            irr_bs = tea.solve_IRR()
        except Exception as e:
            irr_bs = np.nan
            out["IRR_biosteam_error"] = str(e)
        out["IRR_biosteam_raw"] = safe_float(irr_bs)
        out["IRR_biosteam_pct"] = safe_float(irr_bs * 100.0) if np.isfinite(irr_bs) else np.nan

        # irr, irr_reason = solve_IRR_brent(tea)
        # IRR_pct=irr*100 if np.isfinite(irr) else np.nan
        # out["IRR_reason"] = irr_reason
        # IRR = flx.aitken_secant(
        #      lambda r, CF, dur: NPV_at_IRR(max(r, -0.999999), CF, dur),
        #      IRR, 1.0001*IRR + 1e-3,
        #      xtol=1e-6, ytol=10., maxiter=200,
        #      args=args, checkiter=False
        #  )
      
        # IRR_pct=tea.solve_IRR()*100
        add_tea_capex_checks(out, sys, tea, top_n=5)

        gwp = compute_gwp_per_kg_biobinder(sys)

        out.update(
            {
                "OK": True,
                "IRR_pct": safe_float(IRR_pct),
                "GWP": safe_float(gwp),
                "product_ratio_biobinder_over_biofuel": safe_float(prod_ratio),
            }
        )
      
    except Exception as e:
        out.update(
            {
                "OK": False,
                "IRR_pct": np.nan,
                "GWP": np.nan,
                "product_ratio_biobinder_over_biofuel": np.nan,
                "Error": str(e),
            }
        )

    if details["yields_used"] is not None:
        out.update(
            {
                "Y_biocrude": details["yields_used"]["biocrude"],
                "Y_aqueous": details["yields_used"]["aqueous"],
                "Y_gas": details["yields_used"]["gas"],
                "Y_char": details["yields_used"]["char"],
            }
        )

    try:
        splitter = sys.flowsheet.unit.BiocrudeSplitter
        cf = list(splitter.cutoff_fracs)
        out["Cutoff_frac_light"] = safe_float(cf[0])
        out["Cutoff_frac_medium"] = safe_float(cf[1])
        out["Cutoff_frac_heavy"] = safe_float(cf[2])
    except Exception:
        out["Cutoff_frac_light"] = np.nan
        out["Cutoff_frac_medium"] = np.nan
        out["Cutoff_frac_heavy"] = np.nan

    return out



# ============================================================
# Train surrogate + SHAP
# ============================================================
def train_surrogate_and_shap(df, feature_cols, target, cap=None, n_estimators=800):
    d = df.copy()
    d = d[d["OK"]].copy()
    d = d[np.isfinite(d[target])].copy()

    if len(d) < 30:
        raise RuntimeError(f"Too few OK rows for {target}: {len(d)}")
    
    y = d[target].astype(float)
    if cap is not None:
        keep = y.abs() <= float(cap)
        d = d.loc[keep].copy()
        y = d[target].astype(float)
        print(f"[cap filter] {target}: kept {len(d)} rows with |{target}| <= {cap}")
  
    X = d[feature_cols].copy()
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=SEED)

    rf = RandomForestRegressor(
        n_estimators=n_estimators,
        random_state=SEED,
        n_jobs=-1,
        min_samples_leaf=2,
    )
    rf.fit(Xtr, ytr)

    r2 = r2_score(yte, rf.predict(Xte))

    expl = TreeExplainer(rf)
    shap_vals = expl.shap_values(X)

    return rf, r2, X, d, np.array(shap_vals)


def save_shap_longform(X, d_ok, shap_vals, feature_cols, target, out_path, r2_test, extra_cols=None):
    shap_col = f"{target}_SHAP"
    val_col  = f"{target}_value"
    
    extra_cols = extra_cols or []
    d_ok["prefer"] = PREFER

    rows = []
    for i in range(X.shape[0]):
        sid = int(d_ok.iloc[i].get("Sample_ID", i))
        for j, feat in enumerate(feature_cols):
            rows.append({
                "Sample_ID": sid,
                "Target": target,
                "Target_R2_test": float(r2_test),
                "Feature": feat,
                "Feature value": float(X.iloc[i, j]),
                shap_col: float(shap_vals[i, j]),
                val_col: float(d_ok.iloc[i].get(f"{target}_capped", d_ok.iloc[i][target])),
                "Feedstock": d_ok.iloc[i].get("Feedstock", ""),
                "Scenario": d_ok.iloc[i].get("Scenario", ""),
                "product_ratio_biobinder_over_biofuel": d_ok.iloc[i].get("product_ratio_biobinder_over_biofuel", np.nan),
            })
    df_out = pd.DataFrame(rows)
    df_out.to_excel(out_path, index=False)
    print(f"📁 Saved SHAP longform: {out_path}")

def _run_one_task(task):
    """
    task is a plain dict with only picklable stuff.
    Returns a dict result (picklable) that you already produce.
    """
    # Import inside worker to avoid any weird fork/spawn import ordering issues
    from joblib import load

    sample = task["sample"]
    config_name = task["config_name"]
    config_kwargs = task["config_kwargs"]
    feedstock_id = task["feedstock_id"]
    rf_model_path = task["rf_model_path"]

    rf_model = load(rf_model_path) if (rf_model_path and os.path.exists(rf_model_path)) else None

    # IMPORTANT: make sure Sample_ID/Feedstock/Scenario are set here too
    # (in case you do it outside)
    return simulate_one(sample, config_name, config_kwargs, feedstock_id, rf_model)

# ============================================================
# MAIN
# ============================================================
def main():
    N_WORKERS = 10  # or mp.cpu_count() - 2 for safety
    CHUNK = 1       # keep 1 if runs vary a lot; can try 2-5 later
    rf_model = None
    if os.path.exists(RF_MODEL_PATH):
        rf_model = load(RF_MODEL_PATH)
        print(f"✅ RF yield model loaded: {RF_MODEL_PATH}")
    else:
        print(f"⚠️ RF model not found at {RF_MODEL_PATH}. Proceeding without RF yield propagation.")

    for feedstock_id in FEEDSTOCKS:
        for cfg in CONFIGS:
            config_name = cfg["name"]
            config_kwargs = cfg["config_kwargs"]

            print("\n" + "="*80)
            print(f"FEEDSTOCK={feedstock_id} | CONFIG={config_name}")
            print("="*80)

            # meta_path = find_latest_meta(feedstock_id, config_name, PREFER, OUTPUT_DIR)

            if RUN_MODE == "full":
                run_ts = datetime.now().strftime("%Y%m%d_%H%M")
                meta_path = os.path.join(
                                    OUTPUT_DIR,
                              f"META_SHAP3_{feedstock_id}_{config_name}_prefer-{PREFER}_{run_ts}.xlsx"
                 )
                # ---------- RUN SIMULATIONS ----------

                sys_tmp = create_system(feedstock_id=feedstock_id, **config_kwargs)
                spec = build_uncertainty_spec(sys_tmp, config_name)
                var_names = [s["name"] for s in spec]

                feature_names, samples = generate_lhs_samples(spec, N_SAMPLES, seed=SEED)

                tasks = []
                for i, sample in enumerate(samples, 1):
                    sample = dict(sample)  # ensure plain dict
                    sample["Sample_ID"] = i
                    sample["Feedstock"] = feedstock_id
                    sample["Scenario"] = config_name
                    tasks.append({
                        "sample": sample,
                        "config_name": config_name,
                        "config_kwargs": config_kwargs,   # dict -> picklable
                        "feedstock_id": feedstock_id,
                        "rf_model_path": RF_MODEL_PATH,
                    })
                ctx = mp.get_context("spawn")  # explicit for Windows stability
                results = []
                with cf.ProcessPoolExecutor(max_workers=N_WORKERS, mp_context=ctx) as ex:
                    futs = [ex.submit(_run_one_task, t) for t in tasks]
                    for k, fut in enumerate(cf.as_completed(futs), 1):
                        try:
                            res = fut.result()
                        except Exception as e:
                            res = {"OK": False, "Error": f"worker_exception:{type(e).__name__}:{e}"}
                        results.append(res)
                        if k % 10 == 0:
                           print(f"✅ completed {k}/{len(tasks)}")
                             

                for i, sample in enumerate(samples, 1):
                    sample["Sample_ID"] = i
                    sample["Feedstock"] = feedstock_id
                    sample["Scenario"] = config_name

                    print(f"▶️ {feedstock_id} | {config_name} | sample {i}/{N_SAMPLES}")
                    res = simulate_one(sample, config_name, config_kwargs, feedstock_id, rf_model)
                    results.append(res)
                    
                    
                results.sort(key=lambda d: d.get("Sample_ID", 10**9))
                df = pd.DataFrame(results)
                df.to_excel(meta_path, index=False)
                print(f"📁 Saved meta dataset: {meta_path}")
    #         else:
    #             # ---------- SHAP-ONLY OR FULL ----------
    #             meta_path = find_latest_meta(feedstock_id, config_name, PREFER, OUTPUT_DIR)
    #             #         print(f"📖 Loading meta dataset: {meta_path}")
    #             df = pd.read_excel(meta_path)
 

    #         feature_cols = [c for c in df.columns if c in [
    #             # uncertainty knobs
    #             "FS_factor_Water", "FS_factor_Carb", "FS_factor_Prot", "FS_factor_Lipid", "FS_factor_Ash",
    #             "Temperature (C)", "Residence Time", "Solid content (w/w) %",
    #             "Pre-processing", "Catalyst", "Reactor Type", "Solvent",
    #             "Biofuel_price", "Natural_Gas_Price", "Electricity_Price",
    #             "Uptime_Ratio", "income_tax",
    #             # intermediate + separation
    #             "product_ratio_biobinder_over_biofuel",
    #             "Y_biocrude",
    #         ]]

    #         for target in TARGETS:
    #             cap = CAPS.get(target, None)

    #             rf_surr, r2, X_ok, d_ok, shap_vals = train_surrogate_and_shap(
    #                 df=df,
    #                 feature_cols=feature_cols,
    #                 target=target,
    #                 cap=cap,
    #             )

    #             print(f"✅ {target} surrogate trained | R2(test) = {r2:.3f}")

    #             shap_path = os.path.join(
    #                 OUTPUT_DIR,
    #                 f"SHAP3_SURROGATE_{feedstock_id}_{config_name}_{target}_prefer-{PREFER}_{timestamp}.xlsx"
    #             )

    #             save_shap_longform(
    #                 X=X_ok,
    #                 d_ok=d_ok,
    #                 shap_vals=shap_vals,
    #                 feature_cols=feature_cols,
    #                 target=target,
    #                 out_path=shap_path,
    #                 r2_test=r2,
    #                 extra_cols=["Feedstock", "Scenario", "prefer"],
    #             )

    # print("\n✅ SHAP analysis complete.")


if __name__ == "__main__":
    main()
