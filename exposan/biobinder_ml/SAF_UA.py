# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

"""
UA (LHS) for SAF_FLEX-based system with rich CrudeHeavyDis diagnostics.

Key fixes vs old UA:
- Uses exposan.biobinder_ml.SAF_FLEX.create_system
- No global biobinder feedstock_composition dict
- Uses exposan.saf.tea_kwargs
- Correct SAF stream IDs for H2: H2_HC, H2_HT
- Correct transport cost attribute: transportation_unit_cost
- Removes biobinder/BiofuelFlash assumptions
- Keeps + extends CrudeHeavyDis diagnostics
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import time
import glob
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import qmc
from chaospy import distributions as shape
import concurrent.futures as cf
import multiprocessing as mp
from joblib import load

import qsdsan as qs
import flexsolve as flx

from shap import TreeExplainer
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score

# Feedstock composition (wet fractions)
from exposan.biobinder_ml._feedstocks import get_feedstock_composition, BIOBINDER_FEEDSTOCKS

# RF yield predictor (optional)
from exposan.biobinder_ml._rf_htl_predictor import predict_htl_yields_from_rf
from exposan.biobinder_ml._ml_features import FeatureDefaults

# SAF constants + default TEA kwargs (IMPORTANT: use SAF tea_kwargs)
from exposan.saf import _HHV_per_GGE, price_dct, tea_kwargs 
from exposan.biobinder_ml.SAF_FLEX import get_MFSP, get_GWP
# SAF_FLEX system factory (the one that runs)
from exposan.biobinder_ml.SAF_FLEX import create_system as create_saf_system

# ============================================================
# CONFIG
# ============================================================
RUN_MODE = "full"  # full / shap_only (shap_only block kept commented as in your script)
OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M")

FEEDSTOCKS = ["food", "sludge", "manure", "green"]


from exposan.biobinder_ml import central_dry_flowrate as default_central

CONFIGS = [
    {"name": "CHCU_No_EC",
     "config_kwargs": dict(
         flowsheet=None,
         central_dry_flowrate=default_central,
         decentralized_HTL=False,
         decentralized_upgrading=False,
         include_PSA=True,
         include_EC=False,
     )},
     {"name": "DHCU_No_EC",
     "config_kwargs": dict(
          flowsheet=None,
          central_dry_flowrate=default_central,
          decentralized_HTL=True,
          decentralized_upgrading=False,
          include_PSA=True,
          include_EC=False,
      )},
]

N_SAMPLES = 10000
SEED = 42
plant_scale_dtpd=110
PREFER = "max_top_ratio_saf"
RF_MODEL_PATH = os.path.join(OUTPUT_DIR, "rf_yield_model.joblib")

# economics mode:
# - "market": gasoline/jet/diesel prices drive sales; mixed_fuel is bookkeeping (price=0)
# - "mfsp": set gasoline/jet/diesel price=0; solve mixed_fuel price to hit breakeven
SAF_ECON_MODE = "market"  # "market" or "mfsp"

USE_FEEDSTOCK_COMPOSITION_UNCERTAINTY = True
USE_PROCESS_KNOB_UNCERTAINTY = True
USE_RF_YIELD_NOISE = False
USE_SYSTEM_PARAM_UNCERTAINTY = True

TARGETS = ["IRR_pct", "GWP"]
CAPS = {"IRR_pct": 100.0, "GWP": 100.0}


# ============================================================
# Tippping Fee Distribution
# ============================================================

import chaospy as cp
from scipy import stats
tea_indices = qs.utils.indices.tea_indices
_short_ton_to_metric_tonne = 0.907185
cost_year = 2020
PCE_indices = tea_indices['PCEPI']
FOG_year=2021
FOG_factor= PCE_indices[cost_year]/PCE_indices[FOG_year]
SnowdenSwan_year = 2016
SnowdenSwan_factor = PCE_indices[cost_year]/PCE_indices[SnowdenSwan_year]
Green_year= 2000
Green_factor= PCE_indices[cost_year]/PCE_indices[Green_year]
Badgett_year= 2019
Badgett_factor= PCE_indices[cost_year]/PCE_indices[Badgett_year]
# # -----------------------------
# # Custom fitted distributions
# # -----------------------------
# class SludgePriceDist(cp.UserDistribution):
#     def __init__(self):
#         # Fitted to sludge avoided-disposal / tipping-fee data
#         # Units: 2020$ / wet tonne
#         self.a, self.loc, self.scale = -42.3, 8.8, 77.1
#         super().__init__(lower=-160, upper=15)

#     def _cdf(self, x):
#         return stats.skewnorm.cdf(x, self.a, loc=self.loc, scale=self.scale)

#     def _ppf(self, q):
#         return stats.skewnorm.ppf(q, self.a, loc=self.loc, scale=self.scale)


# class ManurePriceDist(cp.UserDistribution):
#     def __init__(self):
#         # Fitted to manure avoided-disposal / tipping-fee data
#         # Units: 2020$ / wet tonne
#         self.a, self.loc, self.scale = -1000000.0, 7.5, 112.8
#         super().__init__(lower=-320, upper=10)

#     def _cdf(self, x):
#         return stats.skewnorm.cdf(x, self.a, loc=self.loc, scale=self.scale)

#     def _ppf(self, q):
#         return stats.skewnorm.ppf(q, self.a, loc=self.loc, scale=self.scale)


# sludge_dist = SludgePriceDist()
# manure_dist = ManurePriceDist()
import numpy as np

# -----------------------------
# MANURE (empirical CDF)
# -----------------------------
MANURE_CDF_P = np.array([
    0.00, 0.03, 0.04, 0.05, 0.05, 0.07, 0.10, 0.11, 0.13, 0.15,
    0.17, 0.20, 0.23, 0.27, 0.30, 0.33, 0.36, 0.39, 0.43, 0.46,
    0.50, 0.54, 0.58, 0.62, 0.66, 0.69, 0.73, 0.76, 0.80, 0.83,
    0.86, 0.89, 0.93, 0.96, 0.98, 1.00
], dtype=float)

MANURE_CDF_X = np.array([
    -310.7, -276.7, -256.6, -216.4, -173.6, -123.3, -119.5, -100.6,
     -90.6,  -88.1,  -73.0,  -73.0,  -67.9,  -65.4,  -64.2,  -57.9,
     -52.8,  -49.1,  -47.8,  -44.0,  -40.3,  -35.2,  -32.7,  -28.9,
     -28.9,  -27.7,  -22.6,  -18.9,  -18.9,  -15.1,  -15.1,  -13.8,
     -12.6,   -3.8,    6.3,    7.5
], dtype=float)

def manure_price_ppf(u):
    u = float(np.clip(u, 0.0, 1.0))
    return float(np.interp(u, MANURE_CDF_P, MANURE_CDF_X))


# -----------------------------
# SLUDGE (empirical CDF)
# -----------------------------
SLUDGE_CDF_P = np.array([
    0.000, 0.028, 0.088, 0.155, 0.218, 0.244, 0.281, 0.330, 0.382,
    0.428, 0.472, 0.535, 0.593, 0.667, 0.740, 0.813, 0.885, 0.950,
    1.000
], dtype=float)

SLUDGE_CDF_X = np.array([
    -152.9, -123.5, -120.6, -113.2, -102.9, -69.1, -58.8, -51.5, -48.5,
     -39.7,  -42.6,  -29.4,  -22.1,   -5.9,    0.0,    2.9,    2.9,    7.4,
       8.8
], dtype=float)

def sludge_price_ppf(u):
    u = float(np.clip(u, 0.0, 1.0))
    return float(np.interp(u, SLUDGE_CDF_P, SLUDGE_CDF_X))

# -----------------------------
# Central feedstock prices
# Units: 2020$ / wet tonne
# -----------------------------
# Replace these numbers if you want updated calibration later.
feedprice_tonne_dct = {
    "fog":    +390.0,    # yellow grease purchase cost, 2020$ / wet tonne
    "green":  -16/_short_ton_to_metric_tonne *Green_factor,   # Leaves & Yardwaste: $20/ton broome county NY; Atlantic County NJ 27$/ton
    "food":   -36/_short_ton_to_metric_tonne,  # 36 $/wet short ton -> 39.69 $/wet metric tonne, Li et al. 2026
    "sludge": -40.0*Badgett_factor,     # representative sludge tipping credit, Li et al 2026
    "manure": -53.0*Badgett_factor,     # representative manure tipping credit, Li et al. 2026
}

# -----------------------------
# Feedstock uncertainty distributions
# Units: 2020$ / wet tonne
# -----------------------------
feedprice_dist_dct = {
    "sludge": "empirical_sludge",
    "manure": "empirical_manure",

    "food": shape.Triangle(
        lower   = -55.0 / _short_ton_to_metric_tonne,
        midpoint= -36.0 / _short_ton_to_metric_tonne,
        upper   = -18.0 / _short_ton_to_metric_tonne,
    ),

    "green": shape.Triangle(
        lower   = -40.0 / _short_ton_to_metric_tonne,
        midpoint= -16.0 / _short_ton_to_metric_tonne * Green_factor,
        upper   = 0.0,
    ),

    "fog": shape.Triangle(
        lower   = 300.0,
        midpoint= 390.0,
        upper   = 600.0,
    ),
}

# ============================================================
# Helpers
# ============================================================
def _ppf(dist, u):
    if isinstance(dist, str):
        if dist == "empirical_manure":
            return manure_price_ppf(u)
        elif dist == "empirical_sludge":
            return sludge_price_ppf(u)
        else:
            raise ValueError(f"Unknown custom distribution tag: {dist}")
    return float(dist.ppf(u))
def safe_float(x):
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


def find_latest_meta(feedstock_id, config_name, prefer, output_dir="results"):
    pattern = os.path.join(
        output_dir,
        f"META_SHAP3_{feedstock_id}_{config_name}_prefer-{prefer}_*.xlsx"
    )
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No META files found for {feedstock_id} | {config_name}")
    return files[-1]

def hx_duty_kJ_per_hr(hx):
    try:
        return float(hx.outs[0].H - hx.ins[0].H)  # kJ/hr
    except Exception:
        return np.nan

def safe_get_design(unit, key):
    try:
        v = unit.design_results.get(key, np.nan)
        return float(v) if np.isfinite(v) else np.nan
    except Exception:
        try:
            return float(unit.design_results[key])
        except Exception:
            return np.nan

def _NPV_at_r(r, CF, dur):
    r = float(r)
    if (not np.isfinite(r)) or (r <= -0.999999):
        return np.nan
    return (CF / (1.0 + r)**dur).sum()

def solve_IRR_aitken_limit(tea, maxiter=200, xtol=1e-6, ytol=10.0, npv_tol_abs=1e3):
    CF = tea.cashflow_array
    dur = tea._get_duration_array()

    if (CF is None) or (dur is None) or (not np.isfinite(CF).all()) or (not np.isfinite(dur).all()):
        return np.nan, "CF_or_dur_nonfinite"

    args = (CF, dur)
    IRR_guess = float(getattr(tea, "IRR", 0.10) or 0.10)
    if (not np.isfinite(IRR_guess)) or (IRR_guess <= -0.999999):
        IRR_guess = 0.10

    x0 = IRR_guess
    x1 = 1.0001 * IRR_guess + 1e-3

    y0 = _NPV_at_r(x0, *args)
    y1 = _NPV_at_r(x1, *args)
    if (not np.isfinite(y0)) or (not np.isfinite(y1)):
        return np.nan, "endpoint_nonfinite"

    try:
        irr = flx.aitken_secant(
            _NPV_at_r,
            x0, x1,
            xtol=xtol, ytol=ytol, maxiter=maxiter,
            args=args, checkiter=False
        )
    except Exception as e:
        return np.nan, f"solver_exception:{type(e).__name__}"

    if (not np.isfinite(irr)) or (irr <= -0.999999):
        return np.nan, "solver_returned_invalid"

    y = _NPV_at_r(irr, *args)
    if not np.isfinite(y):
        return np.nan, "npv_at_solution_nonfinite"
    if abs(float(y)) > float(npv_tol_abs):
        return np.nan, f"npv_not_close:{float(y):.3g}"

    return float(irr), "ok"

# Gasoline gallon equivalent
get_GGE = lambda sys, fuel, annual=True: fuel.HHV/1e3/_HHV_per_GGE*max(1, bool(annual)*sys.operating_hours)

def get_GGE_hr(stream):
    try:
        return float(stream.HHV) / _HHV_per_GGE
    except Exception:
        return np.nan

def saf_blend_ratios(sys, basis="mass"):
    s = sys.flowsheet.stream
    g, j, d = s.gasoline, s.jet, s.diesel

    if basis == "mass":
        denom = float(g.F_mass + j.F_mass + d.F_mass)
        if denom <= 0:
            return (np.nan, np.nan, np.nan)
        return (float(g.F_mass/denom), float(j.F_mass/denom), float(d.F_mass/denom))

    if basis == "gge":
        gg, gj, gd = get_GGE_hr(g), get_GGE_hr(j), get_GGE_hr(d)
        denom = float(gg + gj + gd)
        if (not np.isfinite(denom)) or denom <= 0:
            return (np.nan, np.nan, np.nan)
        return (float(gg/denom), float(gj/denom), float(gd/denom))

    raise ValueError("basis must be 'mass' or 'gge'")



# def compute_irr_pct(sys):
#     """
#     IRR consistent with SAF_ECON_MODE:
#       - market: gasoline/jet/diesel priced, mixed_fuel price=0
#       - mfsp: gasoline/jet/diesel price=0, solve mixed_fuel price
#     """
#     tea = sys.TEA
#     s = sys.flowsheet.stream

#     if SAF_ECON_MODE == "mfsp":
#         for st in ("gasoline", "jet", "diesel"):
#             if hasattr(s, st):
#                 getattr(s, st).price = 0.0
#         if hasattr(s, "mixed_fuel"):
#             _ = get_MFSP(sys)
#     else:
#         if hasattr(s, "mixed_fuel"):
#             s.mixed_fuel.price = 4.87

#     irr, reason = solve_IRR_aitken_limit(tea)
#     if not np.isfinite(irr):
#         return np.nan, reason
#     return float(irr * 100.0), "ok"

def compute_gwp_per_gge_mixed_fuel(sys):
    # SAF_FLEX may not attach sys.LCA; keep safe.
    if (not hasattr(sys, "LCA")) or (sys.LCA is None):
        return np.nan
    lca = sys.LCA
    mf = sys.flowsheet.stream.mixed_fuel
    impacts = lca.get_allocated_impacts(streams=(mf,), operation_only=True, annual=True)
    gwp_annual = float(impacts["GWP"])
    gge_annual = float(get_GGE_hr(mf) * sys.operating_hours)
    if gge_annual <= 0:
        return np.nan
    return float(gwp_annual / gge_annual)

def wet_to_dry_ton_price(price_usd_per_kg_wet, moisture):
    """
    Convert $/kg wet feedstock to $/dry ton.
    """
    dry_frac = 1 - moisture
    if dry_frac <= 0:
        raise ValueError("Dry fraction must be > 0.")
    return price_usd_per_kg_wet * 1000 / dry_frac


def predict_mfsp_3param(scale_dtpd, y_biocrude, feed_cost_usd_per_dry_ton):
    """
    3-parameter MFSP correlation.
    
    scale_dtpd : dry ton/day
    y_biocrude : fraction (0-1)
    feed_cost_usd_per_dry_ton : USD/dry ton
    """
    if scale_dtpd <= 0:
        raise ValueError("scale_dtpd must be > 0")
    if y_biocrude <= 0:
        raise ValueError("y_biocrude must be > 0")

    return (
        6.607 * (scale_dtpd ** -0.6577) / (y_biocrude ** 1.195)
        + 0.00321 * feed_cost_usd_per_dry_ton / (y_biocrude ** 1.062)
        + 2.698
    )
def get_feed_cost_usd_per_dry_ton(feedstock, feedprice_dct, BIOBINDER_FEEDSTOCKS):
    moisture = BIOBINDER_FEEDSTOCKS[feedstock]['moisture']
    wet_price = feedprice_dct[feedstock]   # $/kg wet
    return wet_to_dry_ton_price(wet_price, moisture)

def add_tea_capex_checks(out, sys, tea, top_n=5, min_cost=1.0):
    out["CAPEX_ISBL_installed_$"] = float(getattr(tea, "ISBL_installed_equipment_cost", np.nan))
    out["CAPEX_OSBL_installed_$"] = float(getattr(tea, "OSBL_installed_equipment_cost", np.nan))
    out["CAPEX_TDC_$"] = float(getattr(tea, "TDC", np.nan))
    out["CAPEX_FCI_$"] = float(getattr(tea, "FCI", np.nan))
    out["CAPEX_TCI_$"] = float(getattr(tea, "TCI", np.nan))

    isbl = out["CAPEX_ISBL_installed_$"]
    if np.isfinite(isbl):
        out["CAPEX_buildings_$"]  = 0.04  * isbl
        out["CAPEX_site_dev_$"]   = 0.10  * isbl
        out["CAPEX_add_piping_$"] = 0.045 * isbl
    else:
        out["CAPEX_buildings_$"]  = np.nan
        out["CAPEX_site_dev_$"]   = np.nan
        out["CAPEX_add_piping_$"] = np.nan

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
    out["CAPEX_ISBL_top1_share"] = float(top[0][1] / isbl) if top and np.isfinite(isbl) and isbl > 0 else np.nan

# ============================================================
# Uncertainty spec
# ============================================================
def build_uncertainty_spec(sys, config_name, feedstock_id):
    spec = []

    # -------------------------------------------------
    # (A) Feedstock composition uncertainty
    # -------------------------------------------------
    if USE_FEEDSTOCK_COMPOSITION_UNCERTAINTY:
        spec += [
            {"name": "FS_factor_Water", "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Carb",  "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Prot",  "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Lipid", "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Ash",   "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
        ]

    # -------------------------------------------------
    # (B) HTL process knobs
    # -------------------------------------------------
    if USE_PROCESS_KNOB_UNCERTAINTY:
        spec += [
            {"name": "Temperature (C)",       "dist": shape.Uniform(260, 320), "group": "process"},
            {"name": "Residence Time",        "dist": shape.Uniform(10, 30),   "group": "process"},
            {"name": "Solid content (w/w) %", "dist": shape.Uniform(15, 25),   "group": "process"},
            {"name": "Pre-processing",        "dist": shape.DiscreteUniform(0, 1), "group": "process_cat"},
            {"name": "Catalyst",              "dist": shape.DiscreteUniform(1, 3), "group": "process_cat"},
            {"name": "Reactor Type",          "dist": shape.DiscreteUniform(0, 2), "group": "process_cat"},
            {"name": "Solvent",               "dist": shape.DiscreteUniform(1, 3), "group": "process_cat"},
        ]

    # -------------------------------------------------
    # (C) Downstream product-quality / allocation uncertainty
    #     Keep crude_fracs fixed; vary HT slate + HT oil yield
    # -------------------------------------------------
    spec += [
        {"name": "HT_factor_gasoline", "dist": shape.LogNormal(0.0, 0.10), "group": "product_quality"},
        {"name": "HT_factor_jet",      "dist": shape.LogNormal(0.0, 0.10), "group": "product_quality"},
        {"name": "HT_factor_diesel",   "dist": shape.LogNormal(0.0, 0.10), "group": "product_quality"},
        {"name": "HC_oil_yield_abs",         "dist": shape.Triangle(0.68, 0.7335, 0.78), "group": "upgrading"},
        {"name": "HC_hydrogen_rxn_factor",    "dist": shape.LogNormal(0.0, 0.12),  "group": "upgrading"},
        {"name": "HC_catalyst_life_factor",   "dist": shape.LogNormal(0.0, 0.15),  "group": "upgrading"},

        # base HT oil_yield in SAF_FLEX = 0.8637
        {"name": "HT_oil_yield_abs",   "dist": shape.Triangle(0.82, 0.8637, 0.90), "group": "upgrading"},
    ]

    # -------------------------------------------------
    # (D) System / TEA uncertainties
    #     Fuel prices fixed as requested
    # -------------------------------------------------
    if USE_SYSTEM_PARAM_UNCERTAINTY:
        spec += [
            # feedstock price from PROVIDED distributions, in $ / wet tonne
            {"name": "Feedstock_price_wet_tonne", "dist": feedprice_dist_dct[feedstock_id], "group": "economics"},

            {"name": "price_H2",         "dist": shape.Uniform(1.45, 1.77), "group": "economics"},
            {"name": "price_HCcatalyst", "dist": shape.Triangle(
                0.7 * price_dct["HCcatalyst"],
                price_dct["HCcatalyst"],
                1.3 * price_dct["HCcatalyst"]), "group": "economics"},
            {"name": "price_HTcatalyst", "dist": shape.Triangle(
                0.7 * price_dct["HTcatalyst"],
                price_dct["HTcatalyst"],
                1.3 * price_dct["HTcatalyst"]), "group": "economics"},

            {"name": "Natural_Gas_Price", "dist": shape.Normal(
                price_dct["natural_gas"],
                0.10 * price_dct["natural_gas"]), "group": "economics"},
            {"name": "Electricity_Price", "dist": shape.Triangle(0.03, 0.074, 0.1059), "group": "economics"},
            {"name": "Uptime_Ratio",      "dist": shape.Uniform(0.8, 1.0), "group": "operations"},
            {"name": "income_tax",        "dist": shape.Uniform(0.15, 0.35), "group": "policy_finance"},
        ]

        if "CHCU" in config_name:
            spec += [{"name": "Feedstock_Transportation_Cost", "dist": shape.Uniform(0.048, 0.059), "group": "logistics"}]
        if "DHCU" in config_name:
            spec += [{"name": "Biocrude_Transportation_Cost",  "dist": shape.Uniform(0.089, 0.109), "group": "logistics"}]

    return spec

def generate_lhs_samples(spec, n_samples, seed=42):
    d = len(spec)
    U = _lhs_samples(n_samples, d, seed=seed)
    samples = []
    for i in range(n_samples):
        row = {}
        for j, s in enumerate(spec):
            v = _ppf(s["dist"], U[i, j])
            if ("DiscreteUniform" in str(type(s["dist"]))) or (s["name"] in ["Catalyst", "Reactor Type", "Solvent", "Pre-processing"]):
                v = int(np.round(v))
            row[s["name"]] = float(v)
        samples.append(row)
    return [s["name"] for s in spec], samples

# ============================================================
# Apply sample to SAF_FLEX system
# ============================================================
def apply_sample_to_system(sample, sys, tea, config_name, rf_model=None, base_wet=None):
    wet = dict(base_wet) if base_wet is not None else {}

    # -------------------------------------------------
    # (1) Feedstock composition uncertainty
    # -------------------------------------------------
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

        ssum = sum(new.values())
        if ssum > 0:
            for k in new:
                new[k] /= ssum
        wet = new

        if hasattr(sys.flowsheet.unit, "FeedstockCond"):
            sys.flowsheet.unit.FeedstockCond.feedstock_composition = wet.copy()

    # -------------------------------------------------
    # (2) Process knobs
    # -------------------------------------------------
    process = {}
    if USE_PROCESS_KNOB_UNCERTAINTY:
        T_C = float(sample.get("Temperature (C)", 280.0))
        tau_min = float(sample.get("Residence Time", 15.0))
        solid_pct = float(sample.get("Solid content (w/w) %", 20.0))

        fd = FeatureDefaults()
        process = {
            "Temperature (C)": T_C,
            "Residence Time": tau_min,
            "Solid content (w/w) %": solid_pct,
            "Pre-processing": int(sample.get("Pre-processing", fd.pre_processing)),
            "Catalyst": int(sample.get("Catalyst", fd.catalyst)),
            "Reactor Type": int(sample.get("Reactor Type", fd.reactor_type)),
            "Solvent": int(sample.get("Solvent", fd.solvent)),
            "Reactor Volume (mL)": float(fd.reactor_volume_ml),
        }

        if hasattr(sys.flowsheet.unit, "HTL"):
            sys.flowsheet.unit.HTL.T = T_C + 273.15
            sys.flowsheet.unit.HTL.tau = tau_min / 60.0

        if hasattr(sys.flowsheet.unit, "FeedstockCond"):
            sys.flowsheet.unit.FeedstockCond.target_HTL_solid_loading = solid_pct / 100.0

    # -------------------------------------------------
    # (3) RF yields -> update HTL.dw_yields
    # -------------------------------------------------
    yields_used = None
    if (rf_model is not None) and (base_wet is not None):
        dry_frac = 1.0 - wet.get("Water", 0.0)
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
            yields_used = normalize_dict(y_noisy, ["biocrude", "aqueous", "gas", "char"])
        else:
            yields_used = y

        if hasattr(sys.flowsheet.unit, "HTL"):
            sys.flowsheet.unit.HTL.dw_yields.clear()
            sys.flowsheet.unit.HTL.dw_yields.update(yields_used)

    # -------------------------------------------------
    # (4) System / TEA uncertainties
    # -------------------------------------------------
    if USE_SYSTEM_PARAM_UNCERTAINTY:
        s = sys.flowsheet.stream
        u = sys.flowsheet.unit

        # feedstock price: sampled in $ / wet tonne -> convert to $ / kg wet
        if "Feedstock_price_wet_tonne" in sample and hasattr(s, "scaled_feedstock"):
            s.scaled_feedstock.price = float(sample["Feedstock_price_wet_tonne"]) / 1000.0

        # keep gasoline / jet / diesel prices FIXED
        # do not overwrite them here

        if "price_H2" in sample:
            if hasattr(s, "H2_HC"):
                s.H2_HC.price = float(sample["price_H2"])
            if hasattr(s, "H2_HT"):
                s.H2_HT.price = float(sample["price_H2"])

        if "price_HCcatalyst" in sample and hasattr(s, "HCcatalyst_in"):
            s.HCcatalyst_in.price = float(sample["price_HCcatalyst"])
        if "price_HTcatalyst" in sample and hasattr(s, "HTcatalyst_in"):
            s.HTcatalyst_in.price = float(sample["price_HTcatalyst"])

        if "Natural_Gas_Price" in sample and hasattr(s, "natural_gas"):
            s.natural_gas.price = float(sample["Natural_Gas_Price"])

        if "Electricity_Price" in sample:
            qs.PowerUtility.price = float(sample["Electricity_Price"])
        if "Uptime_Ratio" in sample:
            ur = float(sample["Uptime_Ratio"])
            sys.operating_hours = 365 * 24 * ur
            if hasattr(sys, "LCA") and sys.LCA is not None:
               sys.LCA.uptime_ratio = ur

        if "Feedstock_Transportation_Cost" in sample and hasattr(u, "FeedstockTrans"):
            u.FeedstockTrans.transportation_unit_cost = float(sample["Feedstock_Transportation_Cost"])
        if "Biocrude_Transportation_Cost" in sample and hasattr(u, "BiocrudeTrans"):
            u.BiocrudeTrans.transportation_unit_cost = float(sample["Biocrude_Transportation_Cost"])

        if "income_tax" in sample:
            tea.income_tax = float(sample["income_tax"])
        else:
            tea.income_tax = float(tea_kwargs.get("income_tax", tea.income_tax))
            
        if hasattr(u, "HC"):
            HC = u.HC
            if "HC_oil_yield_abs" in sample:
                y = float(sample["HC_oil_yield_abs"])
                y = np.clip(y, 0.65, 0.85)
                base_gas = 0.2665
                gas = min(base_gas, 0.999 - y)
                HC.oil_yield = y
                HC.gas_yield = gas
            # 2) Hydrogen consumption
            if "HC_hydrogen_rxn_factor" in sample:
                HC.hydrogen_rxned_to_inf_oil = 0.0111 * float(sample["HC_hydrogen_rxn_factor"])
            # 3) Catalyst lifetime
            if "HC_catalyst_life_factor" in sample:
                base = 5 * 365 * 24
                HC.catalyst_lifetime = base * float(sample["HC_catalyst_life_factor"])            
        # -------------------------------------------------
        # HT oil slate uncertainty
        # Base oil_fracs = [0.3455, 0.4479, 0.2066]
        # -------------------------------------------------
        if hasattr(u, "HT"):
            HT = u.HT

            base_oil_fracs = np.array([0.3455, 0.4479, 0.2066], dtype=float)
            facs = np.array([
                float(sample.get("HT_factor_gasoline", 1.0)),
                float(sample.get("HT_factor_jet",      1.0)),
                float(sample.get("HT_factor_diesel",   1.0)),
            ], dtype=float)

            oil_fracs = base_oil_fracs * facs
            oil_fracs = oil_fracs / oil_fracs.sum()

            HT.oil_composition = {
                "Gasoline": float(oil_fracs[0]),
                "Jet":      float(oil_fracs[1]),
                "Diesel":   float(oil_fracs[2]),
            }

            # -------------------------------------------------
            # HT oil yield uncertainty
            # Base HT oil_yield = 0.8637
            # Keep gas_yield <= 1 - oil_yield
            # -------------------------------------------------
            if "HT_oil_yield_abs" in sample:
                new_oil_yield = float(sample["HT_oil_yield_abs"])
                new_oil_yield = min(max(new_oil_yield, 0.75), 0.95)

                base_gas_yield = 0.2143
                max_allowed_gas = max(0.0, 0.999 - new_oil_yield)
                new_gas_yield = min(base_gas_yield, max_allowed_gas)

                HT.oil_yield = float(new_oil_yield)
                HT.gas_yield = float(new_gas_yield)

    return {
        "process": process,
        "yields_used": yields_used,
        "wet_comp_used": wet,
    }
# ============================================================
# CrudeHeavyDis deep diagnostics (max info)
# ============================================================
def add_crudeheavy_dis_diagnostics(out, col, prefix="CHD_"):
    import numpy as _np

    feed = col.ins[0]
    Dstream, Bstream = col.outs

    # Basic feed state
    out[prefix+"feed_F_mass_kgph"] = safe_float(feed.F_mass)
    out[prefix+"feed_F_mol_kmolph"] = safe_float(getattr(feed, "F_mol", np.nan))
    out[prefix+"feed_T_K"] = safe_float(feed.T)
    out[prefix+"feed_P_Pa"] = safe_float(getattr(feed, "P", np.nan))
    out[prefix+"feed_H_kJph"] = safe_float(getattr(feed, "H", np.nan))

    # Nonzero components in feed
    try:
        nz = [ID for ID, m in zip(feed.chemicals.IDs, feed.mol) if m > 0]
        out[prefix+"feed_nonzero_n"] = len(nz)
        out[prefix+"feed_nonzero_IDs_head"] = ",".join(nz[:30])  # keep cell manageable
    except Exception:
        out[prefix+"feed_nonzero_n"] = np.nan
        out[prefix+"feed_nonzero_IDs_head"] = ""

    # Column settings / results
    out[prefix+"P_Pa"] = safe_float(getattr(col, "P", np.nan))
    out[prefix+"P_psi"] = safe_float(getattr(col, "P", np.nan) * 0.000145078) if np.isfinite(getattr(col, "P", np.nan)) else np.nan
    out[prefix+"is_vacuum"] = bool(out[prefix+"P_psi"] < 14.68) if np.isfinite(out[prefix+"P_psi"]) else np.nan

    # policy / runtime best
    best = getattr(col, "_runtime_best", None)
    out[prefix+"policy_found"] = bool(best is not None)
    if best:
        out[prefix+"policy_prefer"] = str(best.get("prefer", ""))
        out[prefix+"policy_LHK"] = str(best.get("LHK", ""))
        out[prefix+"policy_Lr"] = safe_float(best.get("Lr", np.nan))
        out[prefix+"policy_Hr"] = safe_float(best.get("Hr", np.nan))
        out[prefix+"policy_top_ratio"] = safe_float(best.get("top_ratio", np.nan))
        out[prefix+"policy_err"] = safe_float(best.get("err", np.nan))

    # final LHK/Lr/Hr
    try:
        out[prefix+"LHK_light"] = col.LHK[0]
        out[prefix+"LHK_heavy"] = col.LHK[1]
    except Exception:
        out[prefix+"LHK_light"] = str(getattr(col, "LHK", ""))
        out[prefix+"LHK_heavy"] = np.nan
    out[prefix+"Lr_final"] = safe_float(getattr(col, "Lr", np.nan))
    out[prefix+"Hr_final"] = safe_float(getattr(col, "Hr", np.nan))

    # top ratio actually achieved
    try:
        denom = float(getattr(col, "F_mass_out", np.nan))
        out[prefix+"column_top_ratio"] = float(col.outs[0].F_mass / denom) if np.isfinite(denom) and denom > 0 else np.nan
    except Exception:
        out[prefix+"column_top_ratio"] = np.nan

    # Design results (grab a lot)
    D = getattr(col, "design_results", {}) or {}
    out[prefix+"N_theoretical_raw"] = D.get("Theoretical stages", None)
    out[prefix+"N_feed_raw"] = D.get("Theoretical feed stage", None)
    out[prefix+"Reflux"] = safe_get_design(col, "Reflux")
    out[prefix+"Rmin"] = safe_get_design(col, "Minimum reflux")

    out[prefix+"Diameter_ft"] = safe_float(D.get("Rectifier diameter", D.get("Diameter", np.nan)))
    out[prefix+"Height_ft"] = safe_float(D.get("Rectifier height", D.get("Height", np.nan)))
    out[prefix+"twall_in"] = safe_float(D.get("Rectifier wall thickness", D.get("Wall thickness", np.nan)))
    out[prefix+"weight_lb"] = safe_float(D.get("Rectifier weight", D.get("Weight", np.nan)))

    # condenser/reboiler duties + costs
    try:
        out[prefix+"Q_cond_kJhr"] = hx_duty_kJ_per_hr(col.condenser)
        out[prefix+"cond_utility_cost_USDhr"] = safe_float(col.condenser.utility_cost)
    except Exception:
        out[prefix+"Q_cond_kJhr"] = np.nan
        out[prefix+"cond_utility_cost_USDhr"] = np.nan

    try:
        out[prefix+"Q_reb_kJhr"] = hx_duty_kJ_per_hr(col.reboiler)
        out[prefix+"reb_utility_cost_USDhr"] = safe_float(col.reboiler.utility_cost)
    except Exception:
        out[prefix+"Q_reb_kJhr"] = np.nan
        out[prefix+"reb_utility_cost_USDhr"] = np.nan

    # VLE internals (IDs, LK/HK indices)
    IDs = getattr(col, "_IDs_vle", None)
    LHK_vle = getattr(col, "_LHK_vle_index", None)
    out[prefix+"partial_condenser"] = bool(getattr(col, "_partial_condenser", True))
    out[prefix+"IDs_vle_n"] = len(IDs) if IDs else 0

    if IDs and (LHK_vle is not None):
        LK_vle, HK_vle = int(LHK_vle[0]), int(LHK_vle[1])

        # distillate/bottom compositions for IDs
        try:
            zD = Dstream.get_normalized_mol(IDs)
            zB = Bstream.get_normalized_mol(IDs)
            out[prefix+"zD_LK"] = safe_float(zD[LK_vle])
            out[prefix+"zD_HK"] = safe_float(zD[HK_vle])
            out[prefix+"zB_LK"] = safe_float(zB[LK_vle])
            out[prefix+"zB_HK"] = safe_float(zB[HK_vle])
            out[prefix+"zD_min"] = safe_float(np.min(zD))
            out[prefix+"zB_min"] = safe_float(np.min(zB))
        except Exception as e:
            out[prefix+"z_comp_error"] = str(e)

        dew = getattr(col, "_dew_point", None)
        bub = getattr(col, "_bubble_point", None)

        try:
            if getattr(col, "_partial_condenser", True):
                dp = dew(zD, P=col.P)
            else:
                dp = bub(zD, P=col.P)
            bp = bub(zB, P=col.P)

            out[prefix+"dp_T_K"] = safe_float(dp.T)
            out[prefix+"bp_T_K"] = safe_float(bp.T)

            out[prefix+"dp_min_x_raw"] = safe_float(_np.min(dp.x))
            out[prefix+"bp_min_x_raw"] = safe_float(_np.min(bp.x))

            out[prefix+"dp_x_LK"] = safe_float(dp.x[LK_vle])
            out[prefix+"dp_y_LK"] = safe_float(dp.y[LK_vle])
            out[prefix+"dp_x_HK"] = safe_float(dp.x[HK_vle])
            out[prefix+"dp_y_HK"] = safe_float(dp.y[HK_vle])

            out[prefix+"bp_x_LK"] = safe_float(bp.x[LK_vle])
            out[prefix+"bp_y_LK"] = safe_float(bp.y[LK_vle])
            out[prefix+"bp_x_HK"] = safe_float(bp.x[HK_vle])
            out[prefix+"bp_y_HK"] = safe_float(bp.y[HK_vle])

            eps = 1e-16
            Kdist = dp.y / _np.maximum(dp.x, eps)
            Kbot  = bp.y / _np.maximum(bp.x, eps)

            out[prefix+"Kdist_LK"] = safe_float(Kdist[LK_vle])
            out[prefix+"Kdist_HK"] = safe_float(Kdist[HK_vle])
            out[prefix+"Kbot_LK"]  = safe_float(Kbot[LK_vle])
            out[prefix+"Kbot_HK"]  = safe_float(Kbot[HK_vle])

            out[prefix+"alpha_dist_LK_over_HK"] = safe_float(Kdist[LK_vle] / max(Kdist[HK_vle], eps))
            out[prefix+"alpha_bot_LK_over_HK"]  = safe_float(Kbot[LK_vle] / max(Kbot[HK_vle], eps))

            alpha_mean = col._estimate_mean_volatilities_relative_to_heavy_key()
            out[prefix+"alpha_LK_used"] = safe_float(alpha_mean[LK_vle])
            out[prefix+"alpha_LK_le_1"] = bool(alpha_mean[LK_vle] <= 1.0)

            # Fenske ratio using actual LK/HK mole flows (shortcut internal indices)
            try:
                LK_i, HK_i = col._LHK_index
                LK_D = float(Dstream.mol[LK_i])
                HK_D = float(Dstream.mol[HK_i])
                LK_B = float(Bstream.mol[LK_i])
                HK_B = float(Bstream.mol[HK_i])
                ratio = (LK_D / max(HK_D, 1e-30)) * (HK_B / max(LK_B, 1e-30))
                out[prefix+"Fenske_ratio"] = safe_float(ratio)
                out[prefix+"log10_Fenske_ratio"] = safe_float(_np.log10(max(ratio, 1e-30)))
                out[prefix+"log10_alpha_LK"] = safe_float(_np.log10(max(alpha_mean[LK_vle], 1e-30)))
            except Exception as e:
                out[prefix+"Fenske_error"] = str(e)

            # Underwood diagnostics
            try:
                q = col.get_feed_quality()
                theta = col._solve_Underwood_constant(alpha_mean, alpha_mean[LK_vle])
                z_d = Dstream.get_normalized_mol(IDs)
                Rm = float((alpha_mean * z_d / (alpha_mean - theta)).sum() - 1.0)
                out[prefix+"q"] = safe_float(q)
                out[prefix+"theta"] = safe_float(theta)
                out[prefix+"Rm_calc"] = safe_float(Rm)
                out[prefix+"Rm_bad"] = bool((not _np.isfinite(Rm)) or (Rm < 0))
            except Exception as e:
                out[prefix+"Underwood_error"] = str(e)

        except Exception as e:
            out[prefix+"dp_bp_error"] = str(e)
    else:
        out[prefix+"dp_bp_error"] = "missing _IDs_vle or _LHK_vle_index"

    # Mass flows of top/bottom
    try:
        out[prefix+"F_top_kgph"] = safe_float(Dstream.F_mass)
        out[prefix+"F_bot_kgph"] = safe_float(Bstream.F_mass)
    except Exception:
        out[prefix+"F_top_kgph"] = np.nan
        out[prefix+"F_bot_kgph"] = np.nan

# ============================================================
# Run one full simulation
# ============================================================
def simulate_one(sample, config_name, config_kwargs, feedstock_id, rf_model):
    t_run0 = time.perf_counter()
    out = dict(sample)

    sys = create_saf_system(feedstock_id=feedstock_id, **config_kwargs)
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
        hours = float(sys.operating_hours)

        s = sys.flowsheet.stream
        u = sys.flowsheet.unit

        # -------------------------------------------------
        # Realized feedstock price
        # scaled_feedstock.price is $/kg wet in the SAF system
        # -------------------------------------------------
        try:
            out["Feedstock_price_wet_$_per_kg"] = safe_float(s.scaled_feedstock.price)
            out["Feedstock_price_$/tonne"] = safe_float(s.scaled_feedstock.price * 1000.0)
        except Exception:
            out["Feedstock_price_wet_$_per_kg"] = np.nan
            out["Feedstock_price_$/tonne"] = np.nan

        # -------------------------------------------------
        # CrudeHeavyDis diagnostics
        # -------------------------------------------------
        col = u.CrudeHeavyDis
        add_crudeheavy_dis_diagnostics(out, col, prefix="CHD_")

        # -------------------------------------------------
        # Splitter cutoff fractions
        # -------------------------------------------------
        try:
            splitter = u.CrudeSplitter
            cf = list(splitter.cutoff_fracs)
            out["Cutoff_frac_light"] = safe_float(cf[0])
            out["Cutoff_frac_medium"] = safe_float(cf[1])
            out["Cutoff_frac_heavy"] = safe_float(cf[2])
        except Exception:
            out["Cutoff_frac_light"] = np.nan
            out["Cutoff_frac_medium"] = np.nan
            out["Cutoff_frac_heavy"] = np.nan

        # -------------------------------------------------
        # Realized HT oil yield + realized HT oil slate
        # -------------------------------------------------
        try:
            HT = u.HT
            out["HT_oil_yield_realized"] = safe_float(getattr(HT, "oil_yield", np.nan))
            out["HT_gas_yield_realized"] = safe_float(getattr(HT, "gas_yield", np.nan))
            out["HT_hydrogen_rxned_to_inf_oil_realized"] = safe_float(
                getattr(HT, "hydrogen_rxned_to_inf_oil", np.nan)
            )
            out["HT_catalyst_lifetime_hr_realized"] = safe_float(
                getattr(HT, "catalyst_lifetime", np.nan)
            )

            oil_comp = getattr(HT, "oil_composition", None)
            if oil_comp:
                g = safe_float(oil_comp.get("Gasoline", np.nan))
                j = safe_float(oil_comp.get("Jet", np.nan))
                d = safe_float(oil_comp.get("Diesel", np.nan))

                out["HT_oil_frac_gasoline"] = g
                out["HT_oil_frac_jet"] = j
                out["HT_oil_frac_diesel"] = d

                denom = g + j + d
                if np.isfinite(denom) and denom > 0:
                    out["HT_oil_frac_gasoline_norm"] = safe_float(g / denom)
                    out["HT_oil_frac_jet_norm"] = safe_float(j / denom)
                    out["HT_oil_frac_diesel_norm"] = safe_float(d / denom)
                else:
                    out["HT_oil_frac_gasoline_norm"] = np.nan
                    out["HT_oil_frac_jet_norm"] = np.nan
                    out["HT_oil_frac_diesel_norm"] = np.nan
            else:
                out["HT_oil_frac_gasoline"] = np.nan
                out["HT_oil_frac_jet"] = np.nan
                out["HT_oil_frac_diesel"] = np.nan
                out["HT_oil_frac_gasoline_norm"] = np.nan
                out["HT_oil_frac_jet_norm"] = np.nan
                out["HT_oil_frac_diesel_norm"] = np.nan

        except Exception:
            out["HT_oil_yield_realized"] = np.nan
            out["HT_gas_yield_realized"] = np.nan
            out["HT_hydrogen_rxned_to_inf_oil_realized"] = np.nan
            out["HT_catalyst_lifetime_hr_realized"] = np.nan
            out["HT_oil_frac_gasoline"] = np.nan
            out["HT_oil_frac_jet"] = np.nan
            out["HT_oil_frac_diesel"] = np.nan
            out["HT_oil_frac_gasoline_norm"] = np.nan
            out["HT_oil_frac_jet_norm"] = np.nan
            out["HT_oil_frac_diesel_norm"] = np.nan

        # -------------------------------------------------
        # TEA summary
        # -------------------------------------------------
        out["TEA_FCI"] = safe_float(tea.FCI)
        out["TEA_FOC"] = safe_float(getattr(tea, "FOC", np.nan))
        out["TEA_VOC"] = safe_float(getattr(tea, "VOC", np.nan))
        out["TEA_AOC"] = safe_float(getattr(tea, "AOC", np.nan))
        out["SYS_sales"] = safe_float(getattr(sys, "sales", np.nan))
        out["SYS_material_cost"] = safe_float(getattr(sys, "material_cost", np.nan))
        out["SYS_utility_cost"] = safe_float(getattr(sys, "utility_cost", np.nan))

        # -------------------------------------------------
        # Electricity / duties
        # -------------------------------------------------
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

        # -------------------------------------------------
        # NPV diagnostics
        # -------------------------------------------------
        try:
            out["TEA_NPV"] = safe_float(tea.NPV)
            CF = tea.cashflow_array
            dur = tea._get_duration_array()
            out["NPV_at_0pct"] = safe_float(_NPV_at_r(0.0, CF, dur))
            out["NPV_at_5pct"] = safe_float(_NPV_at_r(0.05, CF, dur))
            out["NPV_at_10pct"] = safe_float(_NPV_at_r(0.10, CF, dur))
            out["NPV_at_minus10pct"] = safe_float(_NPV_at_r(-0.10, CF, dur))
        except Exception as e:
            out["NPV_eval_error"] = str(e)

        try:
            out["TEA_net_earnings"] = safe_float(getattr(tea, "net_earnings", np.nan))
        except Exception:
            out["TEA_net_earnings"] = np.nan

        # -------------------------------------------------
        # MFSP / GWP reporting (optional, keep for diagnostics)
        # BUT NOT used as the IRR revenue basis
        # -------------------------------------------------
        try:
            out["MFSP_mixed_fuel_$/GGE"] = safe_float(get_MFSP(sys))
        except Exception:
            out["MFSP_mixed_fuel_$/GGE"] = np.nan

        try:
            out["GWP_mixed_fuel_kgCO2e_per_GGE"] = safe_float(get_GWP(sys))
        except Exception:
            out["GWP_mixed_fuel_kgCO2e_per_GGE"] = np.nan

        # -------------------------------------------------
        # SAF-cut-based IRR logic
        # Do NOT overwrite mixed_fuel.price
        # IRR comes from gasoline / jet / diesel sales already in the TEA
        # -------------------------------------------------
        irr_raw, irr_reason = solve_IRR_aitken_limit(tea)
        out["IRR_reason"] = irr_reason
        out["IRR_raw"] = safe_float(irr_raw)
        out["IRR_pct"] = safe_float(irr_raw * 100.0) if np.isfinite(irr_raw) else np.nan

        # -------------------------------------------------
        # Final product slate metrics
        # -------------------------------------------------
        gasoline = s.gasoline
        jet = s.jet
        diesel = s.diesel

        g_mass = safe_float(gasoline.F_mass)
        j_mass = safe_float(jet.F_mass)
        d_mass = safe_float(diesel.F_mass)

        out["gasoline_F_mass_kgph"] = g_mass
        out["jet_F_mass_kgph"] = j_mass
        out["diesel_F_mass_kgph"] = d_mass

        total_fuel_mass = g_mass + j_mass + d_mass if np.isfinite(g_mass + j_mass + d_mass) else np.nan
        if np.isfinite(total_fuel_mass) and total_fuel_mass > 0:
            out["blend_mass_gasoline"] = safe_float(g_mass / total_fuel_mass)
            out["blend_mass_jet"] = safe_float(j_mass / total_fuel_mass)
            out["blend_mass_diesel"] = safe_float(d_mass / total_fuel_mass)

            # SAF-side allocation / product-ratio metrics
            out["product_ratio_jet_over_gasoline_plus_diesel"] = safe_float(
                j_mass / (g_mass + d_mass)
            ) if (g_mass + d_mass) > 0 else np.nan

            out["jet_share_total_fuel_mass"] = safe_float(j_mass / total_fuel_mass)
            out["gasoline_share_total_fuel_mass"] = safe_float(g_mass / total_fuel_mass)
            out["diesel_share_total_fuel_mass"] = safe_float(d_mass / total_fuel_mass)
        else:
            out["blend_mass_gasoline"] = np.nan
            out["blend_mass_jet"] = np.nan
            out["blend_mass_diesel"] = np.nan
            out["product_ratio_jet_over_gasoline_plus_diesel"] = np.nan
            out["jet_share_total_fuel_mass"] = np.nan
            out["gasoline_share_total_fuel_mass"] = np.nan
            out["diesel_share_total_fuel_mass"] = np.nan

        # GGE-basis slate metrics
        g_g, j_g, d_g = saf_blend_ratios(sys, basis="gge")
        out["blend_gge_gasoline"] = safe_float(g_g)
        out["blend_gge_jet"] = safe_float(j_g)
        out["blend_gge_diesel"] = safe_float(d_g)
        out["jet_share_total_fuel_gge"] = safe_float(j_g)

        # -------------------------------------------------
        # 3-param MFSP comparison block
        # Use REALIZED feedstock price, not a separate dict
        # -------------------------------------------------
        try:
            y_biocrude = np.nan
            if details.get("yields_used", None) is not None:
                y_biocrude = safe_float(details["yields_used"].get("biocrude", np.nan))

            feed_cost_dry_ton = get_feed_cost_usd_per_dry_ton(
                feedstock=feedstock_id,
                feedprice_dct={feedstock_id: safe_float(s.scaled_feedstock.price)},  # $/kg wet realized
                BIOBINDER_FEEDSTOCKS=BIOBINDER_FEEDSTOCKS,
            )

            mfsp_pred_3param = predict_mfsp_3param(
                scale_dtpd=plant_scale_dtpd,
                y_biocrude=y_biocrude,
                feed_cost_usd_per_dry_ton=feed_cost_dry_ton,
            )

            out["Feedstock_cost_dry_$_per_ton"] = safe_float(feed_cost_dry_ton)
            out["MFSP_pred_3param_$/GGE"] = safe_float(mfsp_pred_3param)

            actual_mfsp = safe_float(out.get("MFSP_mixed_fuel_$/GGE", np.nan))
            out["MFSP_ratio_actual_to_3param"] = (
                safe_float(actual_mfsp / mfsp_pred_3param)
                if np.isfinite(actual_mfsp) and np.isfinite(mfsp_pred_3param) and mfsp_pred_3param != 0
                else np.nan
            )
            out["MFSP_diff_actual_minus_3param"] = (
                safe_float(actual_mfsp - mfsp_pred_3param)
                if np.isfinite(actual_mfsp) and np.isfinite(mfsp_pred_3param)
                else np.nan
            )
            out["MFSP_pct_error_vs_3param"] = (
                safe_float(100 * (actual_mfsp - mfsp_pred_3param) / mfsp_pred_3param)
                if np.isfinite(actual_mfsp) and np.isfinite(mfsp_pred_3param) and mfsp_pred_3param != 0
                else np.nan
            )
        except Exception as e:
            out["Feedstock_cost_dry_$_per_ton"] = np.nan
            out["MFSP_pred_3param_$/GGE"] = np.nan
            out["MFSP_ratio_actual_to_3param"] = np.nan
            out["MFSP_diff_actual_minus_3param"] = np.nan
            out["MFSP_pct_error_vs_3param"] = np.nan
            out["MFSP_3param_error"] = str(e)

        # -------------------------------------------------
        # CAPEX diagnostics
        # -------------------------------------------------
        add_tea_capex_checks(out, sys, tea, top_n=5)

        out.update({
            "OK": True,
            "GWP": safe_float(out["GWP_mixed_fuel_kgCO2e_per_GGE"]),
        })

    except Exception as e:
        out.update({
            "OK": False,
            "IRR_pct": np.nan,
            "GWP": np.nan,
            "Error": str(e),
        })

    # -------------------------------------------------
    # Record HTL yields used
    # -------------------------------------------------
    if details.get("yields_used", None) is not None:
        out.update({
            "Y_biocrude": safe_float(details["yields_used"].get("biocrude", np.nan)),
            "Y_aqueous":  safe_float(details["yields_used"].get("aqueous", np.nan)),
            "Y_gas":      safe_float(details["yields_used"].get("gas", np.nan)),
            "Y_char":     safe_float(details["yields_used"].get("char", np.nan)),
        })

    out["Run_time_s"] = float(time.perf_counter() - t_run0)
    return out

# ============================================================
# Surrogate + SHAP (kept for parity; optional)
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
    val_col = f"{target}_value"

    d_ok = d_ok.copy()
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
                val_col: float(d_ok.iloc[i].get(target, np.nan)),
                "Feedstock": d_ok.iloc[i].get("Feedstock", ""),
                "Scenario": d_ok.iloc[i].get("Scenario", ""),
            })
    pd.DataFrame(rows).to_excel(out_path, index=False)
    print(f"📁 Saved SHAP longform: {out_path}")

# ============================================================
# Parallel worker
# ============================================================
def _run_one_task(task):
    import os
    from joblib import load
    sample = task["sample"]
    config_name = task["config_name"]
    config_kwargs = task["config_kwargs"]
    feedstock_id = task["feedstock_id"]
    rf_model_path = task["rf_model_path"]
    rf_model = load(rf_model_path) if (rf_model_path and os.path.exists(rf_model_path)) else None
    return simulate_one(sample, config_name, config_kwargs, feedstock_id, rf_model)


# SAF_PRICE= 1.27
# ============================================================
# MAIN
# ============================================================
def main():
    N_WORKERS = 10

    if os.path.exists(RF_MODEL_PATH):
        print(f"✅ RF yield model found: {RF_MODEL_PATH}")
    else:
        print(f"⚠️ RF model not found at {RF_MODEL_PATH}. Proceeding without RF yield propagation.")

    for feedstock_id in FEEDSTOCKS:
        for cfg in CONFIGS:
            config_name = cfg["name"]
            config_kwargs = cfg["config_kwargs"]

            print("\n" + "=" * 80)
            print(f"FEEDSTOCK={feedstock_id} | CONFIG={config_name}")
            print("=" * 80)

            if RUN_MODE == "full":
                run_ts = datetime.now().strftime("%Y%m%d_%H%M")
                meta_path = os.path.join(
                    OUTPUT_DIR,
                    f"META_SHAP3_SAF_{feedstock_id}_{config_name}_prefer-{PREFER}_price-{SAF_ECON_MODE}_{run_ts}.xlsx",
                )

                # BUILD SAMPLES (IMPORTANT: create SAF system here)
                sys_tmp = create_saf_system(feedstock_id=feedstock_id, **config_kwargs)
                spec = build_uncertainty_spec(sys_tmp, config_name, feedstock_id)
                feature_names, samples = generate_lhs_samples(spec, N_SAMPLES, seed=SEED)

                # BUILD TASKS
                tasks = []
                for i, srow in enumerate(samples, 1):
                    srow = dict(srow)
                    srow["Sample_ID"] = i
                    srow["Feedstock"] = feedstock_id
                    srow["Scenario"] = config_name
                    tasks.append({
                        "sample": srow,
                        "config_name": config_name,
                        "config_kwargs": config_kwargs,
                        "feedstock_id": feedstock_id,
                        "rf_model_path": RF_MODEL_PATH,
                    })

                # RUN PARALLEL
                ctx = mp.get_context("spawn")
                results = []
                t0 = time.perf_counter()
                done = ok = fail = 0

                print(f"🚀 Launching {len(tasks)} simulations with {N_WORKERS} workers...")

                with cf.ProcessPoolExecutor(max_workers=N_WORKERS, mp_context=ctx) as ex:
                    future_to_sid = {}
                    for t in tasks:
                        sid = t["sample"]["Sample_ID"]
                        fut = ex.submit(_run_one_task, t)
                        future_to_sid[fut] = sid

                    for fut in cf.as_completed(future_to_sid):
                        sid = future_to_sid[fut]
                        done += 1
                        try:
                            res = fut.result()
                        except Exception as e:
                            res = {
                                "Sample_ID": sid,
                                "Feedstock": feedstock_id,
                                "Scenario": config_name,
                                "OK": False,
                                "Error": f"worker_exception:{type(e).__name__}:{e}",
                            }

                        results.append(res)
                        if res.get("OK", False):
                            ok += 1
                            status = "OK"
                        else:
                            fail += 1
                            status = "FAIL"

                        elapsed = time.perf_counter() - t0
                        rate = done / elapsed if elapsed > 0 else float("nan")
                        eta = (len(tasks) - done) / rate if rate and rate > 0 else float("nan")
                        rt = res.get("Run_time_s", None)
                        rt_txt = f", run={rt:.1f}s" if isinstance(rt, (int, float)) and np.isfinite(rt) else ""
                        print(f"✅ {done}/{len(tasks)} | Sample {sid} {status}{rt_txt} | OK={ok} FAIL={fail} | ETA~{eta/60:.1f} min")

                # SAVE
                results.sort(key=lambda d: d.get("Sample_ID", 10**9))
                df = pd.DataFrame(results)
                df.to_excel(meta_path, index=False)
                print(f"📁 Saved meta dataset: {meta_path}")

if __name__ == "__main__":
    main()
