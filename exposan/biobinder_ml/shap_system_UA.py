# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 14:56:05 2026

@author: aliah
"""

import os
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import qmc

import qsdsan as qs
from chaospy import distributions as shape

from joblib import load
from shap import TreeExplainer

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


# ============================================================
# CONFIG
# ============================================================
OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M")

FEEDSTOCKS = ["sludge", "manure", "green", "food", "fog"]

CONFIGS = [
    {"name": "CHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central,
                                            "decentralized_HTL": False, "decentralized_upgrading": False,
                                            "skip_EC": True, "generate_H2": False, "EC_config": None}},
    {"name": "DHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central,
                                            "decentralized_HTL": True, "decentralized_upgrading": False,
                                            "skip_EC": True, "generate_H2": False, "EC_config": None}},
]

N_SAMPLES = 100  
SEED = 42

PREFER = "max_top_ratio"
RF_MODEL_PATH = os.path.join(OUTPUT_DIR, "rf_yield_model.joblib")

USE_FEEDSTOCK_COMPOSITION_UNCERTAINTY = True
USE_PROCESS_KNOB_UNCERTAINTY = True
USE_RF_YIELD_NOISE = True
USE_SYSTEM_PARAM_UNCERTAINTY = True

# ---- Targets in surrogate (IRR in % + GWP)
TARGETS = ["IRR_pct", "GWP"]
CAPS = {
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


def compute_irr_percent(tea, products):
    """
    IRR metric as %:
      IRR_pct = tea.solve_IRR() * 100
    """
    if not hasattr(tea, "solve_IRR"):
        raise AttributeError("TEA object does not have solve_IRR().")
    irr = tea.solve_IRR(products)
    return float(irr) * 100.0


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
        col = sys.flowsheet.unit.CrudeHeavyDis
        
        out.update({
             "SEP_LHK_light": col.LHK[0] if isinstance(col.LHK, tuple) else col.LHK,
             "SEP_LHK_heavy": col.LHK[1] if isinstance(col.LHK, tuple) else None,
             "SEP_Lr": float(col.Lr),
             "SEP_Hr": float(col.Hr),
             })

        biobinder = sys.flowsheet.stream.biobinder
        biofuel   = sys.flowsheet.stream.biofuel

        # product ratio metric
        prod_ratio = np.nan
        if biofuel.F_mass and biofuel.F_mass > 0:
            prod_ratio = biobinder.F_mass / biofuel.F_mass

        # IRR output as %
        irr_pct = compute_irr_percent(tea, products=(biobinder,))

        # GWP output
        gwp = compute_gwp_per_kg_biobinder(sys)

        out.update({
            "OK": True,
            "IRR_pct": safe_float(irr_pct),
            "GWP": safe_float(gwp),
            "product_ratio_biobinder_over_biofuel": safe_float(prod_ratio),
        })

    except Exception as e:
        out.update({
            "OK": False,
            "IRR_pct": np.nan,
            "GWP": np.nan,
            "product_ratio_biobinder_over_biofuel": np.nan,
            "Error": str(e),
        })

    # yields used (optional)
    if details["yields_used"] is not None:
        out.update({
            "Y_biocrude": details["yields_used"]["biocrude"],
            "Y_aqueous":  details["yields_used"]["aqueous"],
            "Y_gas":      details["yields_used"]["gas"],
            "Y_char":     details["yields_used"]["char"],
        })

    # cutoff applied
    # if details["cutoff_fracs_applied"] is not None:
    #     cf = details["cutoff_fracs_applied"]
    #     out["Cutoff_frac_light"]  = safe_float(cf[0])
    #     out["Cutoff_frac_medium"] = safe_float(cf[1])
    #     out["Cutoff_frac_heavy"]  = safe_float(cf[2])
    #     out["Cutoff_ratio_M_over_H_used"] = safe_float(cf[1] / cf[2]) if cf[2] else np.nan
    try:
        splitter = sys.flowsheet.unit.BiocrudeSplitter
        cf = list(splitter.cutoff_fracs)  # public property
        out["Cutoff_frac_light"]  = safe_float(cf[0])
        out["Cutoff_frac_medium"] = safe_float(cf[1])
        out["Cutoff_frac_heavy"]  = safe_float(cf[2])
    except Exception as e:
        
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


# ============================================================
# MAIN
# ============================================================
def main():
    rf_model = None
    if os.path.exists(RF_MODEL_PATH):
        rf_model = load(RF_MODEL_PATH)
        print(f"✅ RF yield model loaded: {RF_MODEL_PATH}")
    else:
        print(f"⚠️ RF model not found at {RF_MODEL_PATH}. Proceeding without RF yield propagation.")
        print("   In that case, SHAP-3 will NOT attribute uncertainty to yield-related upstream knobs.")

    for feedstock_id in FEEDSTOCKS:
        for cfg in CONFIGS:
            config_name = cfg["name"]
            config_kwargs = cfg["config_kwargs"]

            print("\n" + "="*80)
            print(f"FEEDSTOCK={feedstock_id} | CONFIG={config_name} | N={N_SAMPLES}")
            print("="*80)

            sys_tmp = create_system(feedstock_id=feedstock_id, **config_kwargs)

            spec = build_uncertainty_spec(sys_tmp, config_name)
            var_names = [s["name"] for s in spec]

            feature_names, samples = generate_lhs_samples(spec, N_SAMPLES, seed=SEED)

            results = []
            for i, sample in enumerate(samples, 1):
                sample["Sample_ID"] = i
                sample["Feedstock"] = feedstock_id
                sample["Scenario"] = config_name

                print(f"▶️ {feedstock_id} | {config_name} | sample {i}/{N_SAMPLES}")
                res = simulate_one(sample, config_name, config_kwargs, feedstock_id, rf_model)
                results.append(res)

            df = pd.DataFrame(results)

            meta_path = os.path.join(OUTPUT_DIR, f"META_SHAP3_{feedstock_id}_{config_name}_prefer-{PREFER}_{timestamp}.xlsx")
            df.to_excel(meta_path, index=False)
            print(f"📁 Saved meta dataset: {meta_path}")
            print(f"df: {df.shape} | OK: {df['OK'].sum()} / {len(df)}")

            feature_cols = [c for c in var_names if c in df.columns]
            if "product_ratio_biobinder_over_biofuel" in df.columns:
                feature_cols.append("product_ratio_biobinder_over_biofuel")

            for target in TARGETS:
                cap = CAPS.get(target, None)
                rf_surr, r2, X_ok, d_ok, shap_vals = train_surrogate_and_shap(
                    df=df,
                    feature_cols=feature_cols,
                    target=target,
                    cap=cap,
                )
                cap_msg = f"(capped ±{cap})" if cap is not None else ""
                print(f"✅ {target} surrogate trained {cap_msg} | R2(test) = {r2:.3f}")

                shap_path = os.path.join(OUTPUT_DIR, f"SHAP3_SURROGATE_{feedstock_id}_{config_name}_{target}_prefer-{PREFER}_{timestamp}.xlsx")
                save_shap_longform(
                    X=X_ok,
                    d_ok=d_ok,
                    shap_vals=shap_vals,
                    feature_cols=feature_cols,
                    target=target,
                    out_path=shap_path,
                    r2_test=r2,
                    extra_cols=["Feedstock", "Scenario", "prefer"],
                )

    print("\n✅ SHAP-3 uncertainty propagation complete (IRR% + GWP).")


if __name__ == "__main__":
    main()
