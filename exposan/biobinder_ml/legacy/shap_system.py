# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 13:24:44 2026

@author: aliah
"""

# ============================================================
# SHAP-3: End-to-end uncertainty propagation (System sensitivity)
# ============================================================
# What this script does
# ---------------------
# 1) Samples uncertain "upstream" inputs (feedstock composition, process knobs,
#    optional RF-yield noise) + uncertain "system/TEA" parameters (prices, uptime,
#    IRR, taxes, EC params, transport cost, etc.)
# 2) Runs the full EXPOsan biobinder system to compute downstream outputs:
#       MSP ($/kg), IRR (fraction), GWP (kg CO2e/kg)
# 3) Trains a surrogate model for each output on successful runs
# 4) Uses SHAP on the surrogates to quantify which uncertainties dominate
#
# Notes
# -----
# • You do NOT need to change create_system() for this to work, provided that:
#   - your create_system already uses HTL_yields (global) and feedstock_composition (global),
#   - OR you can override HTL_yields/feedstock_composition BEFORE calling create_system.
# • If your create_system now computes RF yields internally from feedstock/process,
#   keep it ON, but then the sampled "RF noise" should be applied by perturbing
#   the yields AFTER prediction and BEFORE create_system, which this script supports.
#
# Requirements
# ------------
# pip/conda: shap, scikit-learn, openpyxl
#
# ============================================================

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
    # create_system,
    HTL_yields,                # global dict used by HTL unit in your flowsheet
    feedstock_composition,     # global dict in biobinder module (if used)
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

# Feedstocks to run SHAP-3 for (you can run one or many)
FEEDSTOCKS = ["sludge", "manure", "green", "food", "fog"]  

# Scenarios (same naming you used earlier)
CONFIGS = [
    {"name": "CHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central,
                                            "decentralized_HTL": False, "decentralized_upgrading": False,
                                            "skip_EC": True, "generate_H2": False, "EC_config": None}},
    # {"name": "CHCU_EC",    "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central,
    #                                         "decentralized_HTL": False, "decentralized_upgrading": False,
    #                                         "skip_EC": False, "generate_H2": False, "EC_config": None}},
    {"name": "DHCU_No_EC", "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central,
                                            "decentralized_HTL": True, "decentralized_upgrading": False,
                                            "skip_EC": True, "generate_H2": False, "EC_config": None}},
    # {"name": "DHCU_EC",    "config_kwargs": {"flowsheet": None, "central_dry_flowrate": default_central,
    #                                         "decentralized_HTL": True, "decentralized_upgrading": False,
    #                                         "skip_EC": False, "generate_H2": False, "EC_config": None}},
]

# Sampling
N_SAMPLES = 1000  
SEED = 42

PREFER = "max_top_ratio"

# RF model path (if you want yield propagation via RF → HTL yields)
RF_MODEL_PATH = os.path.join(OUTPUT_DIR, "rf_yield_model.joblib")

# Toggle these uncertainty blocks on/off
USE_FEEDSTOCK_COMPOSITION_UNCERTAINTY = True
USE_PROCESS_KNOB_UNCERTAINTY = True
USE_RF_YIELD_NOISE = True  # optional stochastic yield perturbation around RF prediction
USE_SYSTEM_PARAM_UNCERTAINTY = True

# Which downstream targets to build surrogates for
TARGETS = ["MSP", "GWP"]  # IRR will be your TEA.IRR if you sample it; otherwise baseline

# Caps to stabilize surrogate training (important!)
CAPS = {
    "MSP": 20.0,   # $/kg cap magnitude (±20) – you can also use [0, 20] if MSP is non-negative
    "GWP": 50.0,   # kgCO2e/kg cap magnitude (tune to your model scale)
    # IRR typically is within [0, 0.3] so no cap needed; you can add if desired
}

# ============================================================
# Helpers
# ============================================================
def clip_symmetric(x, cap):
    if cap is None:
        return x
    return np.clip(x, -cap, cap)

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return np.nan

def normalize_dict(d, keys):
    s = sum(float(d[k]) for k in keys)
    if s <= 0:
        return d
    for k in keys:
        d[k] = float(d[k]) / s
    return d

def _lhs_samples(n, d, seed=42):
    sampler = qmc.LatinHypercube(d=d, seed=seed)
    return sampler.random(n=n)  # U(0,1)

def _ppf(dist, u):
    # Chaospy distributions support .ppf()
    return float(dist.ppf(u))

def decode_feature_name(name):
    # Optional: keep as-is; you can expand for nicer labels
    return name

# ============================================================
# Define uncertain inputs (distributions)
# ============================================================
def build_uncertainty_spec(sys, config_name):
    """
    Return a list of uncertainty variables with:
      - name
      - apply(sample_value, sys, tea, lca)
      - dist (chaospy distribution)
      - group (for later grouping/plots)
    """
    spec = []

    # -----------------------------
    # (A) Feedstock composition uncertainty (wet fractions)
    # -----------------------------
    # We sample around the canonical feedstock composition from get_feedstock_composition().
    # Constraint: components must remain non-negative; then renormalize.
    # You can adjust CV/widths per component.
    if USE_FEEDSTOCK_COMPOSITION_UNCERTAINTY:
        # We will sample multiplicative factors for each component (lognormal-ish)
        # then renormalize to sum to 1.
        # Factors centered at 1.0 with ~10% spread (tune).
        spec += [
            {"name": "FS_factor_Water", "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Carb",  "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Prot",  "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Lipid", "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
            {"name": "FS_factor_Ash",   "dist": shape.LogNormal(0.0, 0.10), "group": "composition"},
        ]

    # -----------------------------
    # (B) HTL process knobs (these are your 7 “knobs”)
    # -----------------------------
    if USE_PROCESS_KNOB_UNCERTAINTY:
        spec += [
            {"name": "Temperature (C)",          "dist": shape.Uniform(260, 320), "group": "process"},
            {"name": "Residence Time",           "dist": shape.Uniform(10, 30),   "group": "process"},
            {"name": "Solid content (w/w) %",    "dist": shape.Uniform(15, 25),   "group": "process"},
            # Categorical knobs: sample as discrete integers
            {"name": "Pre-processing",           "dist": shape.DiscreteUniform(0, 1), "group": "process_cat"},
            {"name": "Catalyst",                 "dist": shape.DiscreteUniform(1, 3), "group": "process_cat"},
            {"name": "Reactor Type",             "dist": shape.DiscreteUniform(0, 2), "group": "process_cat"},
            {"name": "Solvent",                  "dist": shape.DiscreteUniform(1, 3), "group": "process_cat"},
        ]

    # -----------------------------
    # (C) System/TEA/LCA uncertainties (your earlier UA parameters)
    # -----------------------------
    if USE_SYSTEM_PARAM_UNCERTAINTY:
        # Stream prices, utilities, uptime, TEA policy
        spec += [
            {"name": "Biofuel_price",          "dist": shape.Uniform(1.02, 1.25),          "group": "economics"},
            {"name": "N_Fertilizer_Price",     "dist": shape.Triangle(0.85, 0.9, 1.19),    "group": "economics"},
            {"name": "Natural_Gas_Price",      "dist": shape.Normal(0.14, 0.0211),         "group": "economics"},
            {"name": "Electricity_Price",      "dist": shape.Triangle(0.03, 0.074, 0.1059),"group": "economics"},
            {"name": "Uptime_Ratio",           "dist": shape.Uniform(0.8, 1.0),            "group": "operations"},
            {"name": "IRR",                    "dist": shape.Uniform(0.05, 0.15),          "group": "policy_finance"},
            {"name": "income_tax",             "dist": shape.Uniform(0.15, 0.35),          "group": "policy_finance"},
        ]

        # Transport cost depends on CHCU vs DHCU
        if "CHCU" in config_name:
            spec += [{"name": "Feedstock_Transportation_Cost", "dist": shape.Uniform(0.048, 0.059), "group": "logistics"}]
        if "DHCU" in config_name:
            spec += [{"name": "Biocrude_Transportation_Cost",  "dist": shape.Uniform(0.089, 0.109), "group": "logistics"}]

        # EC params only if EC scenario
        if config_name in ["CHCU_EC", "DHCU_EC"]:
            for unit_name in ["HTL_EC", "Upgrading_EC"]:
                # Only include if unit exists in flowsheet
                try:
                    _ = getattr(sys.flowsheet.unit, unit_name)
                    spec += [
                        {"name": f"{unit_name}_EO_voltage",     "dist": shape.Uniform(2.5, 10),         "group": "EC"},
                        {"name": f"{unit_name}_ED_voltage",     "dist": shape.Uniform(2.5, 50),         "group": "EC"},
                        {"name": f"{unit_name}_electrode_cost","dist": shape.Triangle(225, 40000, 80000),"group": "EC"},
                        {"name": f"{unit_name}_AEM_cost",      "dist": shape.Triangle(0, 170, 250),     "group": "EC"},
                        {"name": f"{unit_name}_CEM_cost",      "dist": shape.Triangle(0, 190, 250),     "group": "EC"},
                        {"name": f"{unit_name}_COD_removal",   "dist": shape.Uniform(0.9, 1.0),         "group": "EC"},
                        {"name": f"{unit_name}_N_recovery",    "dist": shape.Uniform(0.72, 0.88),       "group": "EC"},
                    ]
                except AttributeError:
                    pass

    # -----------------------------
    # (D) RF-yield noise (optional)
    # -----------------------------
    # We add independent noise to each yield component in fraction space, then renormalize.
    # If you want correlated noise, switch to Dirichlet noise.
    if USE_RF_YIELD_NOISE:
        spec += [
            {"name": "RF_noise_biocrude", "dist": shape.Normal(0.0, 0.02), "group": "yield_noise"},
            {"name": "RF_noise_aqueous",  "dist": shape.Normal(0.0, 0.02), "group": "yield_noise"},
            {"name": "RF_noise_gas",      "dist": shape.Normal(0.0, 0.02), "group": "yield_noise"},
            {"name": "RF_noise_char",     "dist": shape.Normal(0.0, 0.02), "group": "yield_noise"},
        ]

    return spec


# ============================================================
# Sampling engine
# ============================================================
def generate_lhs_samples(spec, n_samples, seed=42):
    names = [s["name"] for s in spec]
    d = len(spec)
    U = _lhs_samples(n_samples, d, seed=seed)
    samples = []
    for i in range(n_samples):
        row = {}
        for j, s in enumerate(spec):
            v = _ppf(s["dist"], U[i, j])
            # For discrete uniforms coming back as float, round safely to int
            if "DiscreteUniform" in str(type(s["dist"])) or s["name"] in ["Catalyst", "Reactor Type", "Solvent", "Pre-processing"]:
                v = int(np.round(v))
            row[s["name"]] = float(v)
        samples.append(row)
    return names, samples


# ============================================================
# Apply a sampled scenario to system
# ============================================================
def apply_sample_to_system(sample, sys, tea, lca, config_name, rf_model=None, feedstock_id=None, base_wet=None):
    """
    Applies sampled uncertainties:
      - feedstock composition (global dict or module dict)
      - process knobs -> used to RF-predict yields
      - RF noise -> perturbs yields
      - TEA/LCA/system params -> sets attributes on flowsheet objects
    Returns a dict with:
      - feature snapshot (all sampled inputs)
      - predicted/used yields (biocrude/aqueous/gas/char)
    """
    # ---------- (1) Feedstock composition update (global) ----------
    wet = dict(base_wet) if base_wet is not None else {}

    if USE_FEEDSTOCK_COMPOSITION_UNCERTAINTY and base_wet is not None:
        # Multiply components by sampled factors; renormalize
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

        # Renormalize to sum 1
        s = sum(new.values())
        if s > 0:
            for k in new:
                new[k] /= s

        wet = new

    # Push to module-global feedstock_composition if your create_system reads it,
    # AND also pass feedstock_id into create_system via config_kwargs outside.
    # (You are already passing feedstock_id via create_system arguments in your newer scripts.)
    # Here we update the module-global dict just in case.
    for k in list(feedstock_composition.keys()):
        feedstock_composition[k] = wet.get(k, feedstock_composition[k])

    # ---------- (2) Build process conditions (knobs) ----------
    process = {}
    if USE_PROCESS_KNOB_UNCERTAINTY:
        process = {
            "Temperature (C)": float(sample.get("Temperature (C)", 280.0)),
            "Residence Time": float(sample.get("Residence Time", 15.0)),
            "Solid content (w/w) %": float(sample.get("Solid content (w/w) %", 20.0)),
            "Pre-processing": int(sample.get("Pre-processing", FeatureDefaults().pre_processing)),
            "Catalyst": int(sample.get("Catalyst", FeatureDefaults().catalyst)),
            "Reactor Type": int(sample.get("Reactor Type", FeatureDefaults().reactor_type)),
            "Solvent": int(sample.get("Solvent", FeatureDefaults().solvent)),
            "Reactor Volume (mL)": float(FeatureDefaults().reactor_volume_ml),
        }

    # ---------- (3) Predict yields via RF (based on sampled feedstock + sampled process) ----------
    yields_used = None
    if rf_model is not None and base_wet is not None:
        # Convert wet fractions -> dry wt%
        dry_frac = 1.0 - wet["Water"]
        if dry_frac <= 0:
            raise ValueError("Non-physical water fraction in sampled feedstock composition.")

        carb_wt    = wet["Carbohydrates"] / dry_frac * 100.0
        protein_wt = wet["Proteins"]      / dry_frac * 100.0
        lipids_wt  = wet["Lipids"]        / dry_frac * 100.0
        ash_wt     = wet["Ash"]           / dry_frac * 100.0

        # If process knobs are OFF, use defaults
        if not process:
            process = {
                "Temperature (C)": 280.0,
                "Residence Time": 15.0,
                "Solid content (w/w) %": 20.0,
                "Pre-processing": FeatureDefaults().pre_processing,
                "Catalyst": FeatureDefaults().catalyst,
                "Reactor Type": FeatureDefaults().reactor_type,
                "Solvent": FeatureDefaults().solvent,
                "Reactor Volume (mL)": FeatureDefaults().reactor_volume_ml,
            }

        y = predict_htl_yields_from_rf(
            rf_model=rf_model,
            carb_wt=carb_wt,
            protein_wt=protein_wt,
            lipids_wt=lipids_wt,
            ash_wt=ash_wt,
            process=process,
        )

        # Add optional RF yield noise (and renormalize)
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

        # Push yields to module-global HTL_yields so create_system uses them
        HTL_yields.update(yields_used)

    # ---------- (4) Apply system/TEA/LCA parameter samples ----------
    if USE_SYSTEM_PARAM_UNCERTAINTY:
        # Prices/utilities/uptime
        if "Biofuel_price" in sample:
            sys.flowsheet.stream.biofuel.price = float(sample["Biofuel_price"])
        if "N_Fertilizer_Price" in sample:
            sys.flowsheet.stream.recovered_N.price = float(sample["N_Fertilizer_Price"])
        if "Natural_Gas_Price" in sample:
            sys.flowsheet.stream.natural_gas.price = float(sample["Natural_Gas_Price"])
        if "Electricity_Price" in sample:
            qs.PowerUtility.price = float(sample["Electricity_Price"])
        if "Uptime_Ratio" in sample:
            sys.flowsheet.unit.HTL.Uptime_ratio = float(sample["Uptime_Ratio"])

        # Transport costs
        if "Feedstock_Transportation_Cost" in sample:
            sys.flowsheet.unit.FeedstockTrans.cost = float(sample["Feedstock_Transportation_Cost"])
        if "Biocrude_Transportation_Cost" in sample:
            sys.flowsheet.unit.BiocrudeTrans.cost = float(sample["Biocrude_Transportation_Cost"])

        # TEA policy
        if "IRR" in sample:
            tea.IRR = float(sample["IRR"])
        else:
            tea.IRR = tea_kwargs["IRR"]

        if "income_tax" in sample:
            tea.income_tax = float(sample["income_tax"])
        else:
            tea.income_tax = tea_kwargs["income_tax"]

        # EC params
        if config_name in ["CHCU_EC", "DHCU_EC"]:
            for unit_name in ["HTL_EC", "Upgrading_EC"]:
                unit = getattr(sys.flowsheet.unit, unit_name, None)
                if unit is None:
                    continue
                for suffix, attr in [
                    ("EO_voltage", "EO_voltage"),
                    ("ED_voltage", "ED_voltage"),
                    ("electrode_cost", "electrode_cost"),
                    ("AEM_cost", "Anion_exchange_membrane_cost"),
                    ("CEM_cost", "Cation_exchange_membrane_cost"),
                    ("COD_removal", "COD_removal"),
                    ("N_recovery", "N_recovery"),
                ]:
                    key = f"{unit_name}_{suffix}"
                    if key in sample:
                        setattr(unit, attr, float(sample[key]))

    # Return extra details for logging
    return {
        "process": process,
        "yields_used": yields_used,
        "wet_comp_used": wet,
    }


# ============================================================
# Run one full simulation with one sample
# ============================================================
def simulate_one(sample, config_name, config_kwargs, feedstock_id, rf_model):
    """
    Build system fresh for each sample (robust),
    apply sample, simulate, compute outputs.
    """
    out = dict(sample)  # record inputs

    # 1) Create system
    sys = create_system(feedstock_id=feedstock_id, **config_kwargs)
    tea = sys.TEA
    lca = sys.LCA

    # 2) Build base composition for this feedstock from your library
    base_wet = get_feedstock_composition(feedstock_id)

    # 3) Apply sampled uncertainties (updates globals + sets sys params)
    details = apply_sample_to_system(
        sample=sample,
        sys=sys,
        tea=tea,
        lca=lca,
        config_name=config_name,
        rf_model=rf_model,
        feedstock_id=feedstock_id,
        base_wet=base_wet,
    )

    # 4) Simulate + compute KPIs
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

        MSP = tea.solve_price(biobinder)
        # GWP in kg CO2e/kg product
        GWP = (
            lca.get_allocated_impacts(streams=(biobinder,), operation_only=True, annual=True)["GWP"]
            / (biobinder.F_mass * lca.system.operating_hours)
        )

        out.update({
            "OK": True,
            "MSP": safe_float(MSP),
            "GWP": safe_float(GWP),
            "IRR": safe_float(tea.IRR),
        })

    except Exception as e:
        out.update({
            "OK": False,
            "MSP": np.nan,
            "GWP": np.nan,
            "IRR": safe_float(tea.IRR),
            "Error": str(e),
        })

    # 5) log extra details (optional)
    if details["yields_used"] is not None:
        out.update({
            "Y_biocrude": details["yields_used"]["biocrude"],
            "Y_aqueous":  details["yields_used"]["aqueous"],
            "Y_gas":      details["yields_used"]["gas"],
            "Y_char":     details["yields_used"]["char"],
        })

    if details["process"]:
        for k, v in details["process"].items():
            out[f"KNOB::{k}"] = v

    if details["wet_comp_used"]:
        for k, v in details["wet_comp_used"].items():
            out[f"WET::{k}"] = v

    return out


# ============================================================
# Train surrogate + SHAP
# ============================================================
def train_surrogate_and_shap(df, feature_cols, target, cap=None, n_estimators=800):
    """
    Train RF surrogate for target, compute SHAP values on X_all.
    Returns rf_model, r2_test, shap_values (n_rows x n_features)
    """
    d = df.copy()
    d = d[d["OK"]].copy()
    d = d[np.isfinite(d[target])].copy()

    if len(d) < 30:
        raise RuntimeError(f"Too few OK rows for {target}: {len(d)}")

    y = d[target].astype(float).values
    if cap is not None:
        y = clip_symmetric(y, cap)

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
    shap_vals = expl.shap_values(X)  # (n_rows, n_features)

    return rf, r2, X, d, np.array(shap_vals)


def save_shap_longform(X, d_ok, shap_vals, feature_cols, target, r2_test, out_path, extra_cols=None):
    """
    Long-form SHAP table:
      one row per (sample, feature)
    """
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
                f"{target}_SHAP": float(shap_vals[i, j]),
                f"{target}_value": float(d_ok.iloc[i][target]),
                **{c: d_ok.iloc[i].get(c, np.nan) for c in extra_cols},
            })
    df_out = pd.DataFrame(rows)
    df_out.to_excel(out_path, index=False)
    print(f"📁 Saved SHAP longform: {out_path}")


# ============================================================
# MAIN
# ============================================================
def main():
    # Load RF yield model (optional but recommended for true end-to-end propagation)
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

            # Build a temporary system once to detect available units (EC params etc.)
            sys_tmp = create_system(feedstock_id=feedstock_id, **config_kwargs)

            # Build uncertainty spec for this scenario
            spec = build_uncertainty_spec(sys_tmp, config_name)
            var_names = [s["name"] for s in spec]

            # Generate LHS samples
            feature_names, samples = generate_lhs_samples(spec, N_SAMPLES, seed=SEED)

            # Run Monte Carlo (fresh system per sample)
            results = []
            for i, sample in enumerate(samples, 1):
                sample["Sample_ID"] = i
                sample["Feedstock"] = feedstock_id
                sample["Scenario"] = config_name

                print(f"▶️ {feedstock_id} | {config_name} | sample {i}/{N_SAMPLES}")
                res = simulate_one(sample, config_name, config_kwargs, feedstock_id, rf_model)
                results.append(res)

            df = pd.DataFrame(results)

            # Save raw meta dataset
            meta_path = os.path.join(OUTPUT_DIR, f"META_SHAP3_{feedstock_id}_{config_name}_{timestamp}.xlsx")
            df.to_excel(meta_path, index=False)
            print(f"📁 Saved meta dataset: {meta_path}")
            print(f"df: {df.shape} | OK: {df['OK'].sum()} / {len(df)}")

            # Build surrogate feature matrix columns:
            # Use only the sampled uncertainty variables (var_names) as features
            feature_cols = [c for c in var_names if c in df.columns]

            # Train + SHAP per target
            for target in TARGETS:
                cap = CAPS.get(target, None)

                # For IRR:
                # If you did NOT sample IRR, it will be constant and surrogate training will fail.
                # In that case, skip IRR.
                if target == "IRR" and ("IRR" not in df.columns or df["IRR"].nunique(dropna=True) <= 2):
                    print("⚠️ Skipping IRR surrogate (IRR not varying / too few unique values).")
                    continue

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
                    r2_test=r2,
                    out_path=shap_path,
                    extra_cols=["Feedstock", "Scenario", "prefer"],
                )

    print("\n✅ SHAP-3 end-to-end uncertainty propagation complete.")


if __name__ == "__main__":
    main()