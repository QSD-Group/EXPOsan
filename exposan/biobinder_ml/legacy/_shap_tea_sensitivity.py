# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 13:36:47 2026

@author: aliah
"""

import os
import numpy as np
import pandas as pd
from datetime import datetime

from joblib import load
from shap import TreeExplainer

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

from exposan.biobinder_ml._feedstocks import get_feedstock_composition

# Your feature adapter + decoding maps (same as your working yield SHAP script)
from exposan.biobinder_ml._ml_features import (
    FEATURE_COLS,
    make_feature_row,
    PREPROCESSING_MAP,
    CATALYST_MAP,
    SOLVENT_MAP,
    REACTOR_TYPE_MAP,
)

# Your system + globals (same structure you already use)
from exposan.biobinder_ml import (
    HTL_yields,          # global dict used by create_system
    tea_kwargs,
)
from exposan.biobinder_ml.Dist_flex_shap import create_system


# -----------------------------
# Decode helper (same pattern)
# -----------------------------
def decode_feature_value(feature, value):
    try:
        iv = int(value)
    except Exception:
        return value  # continuous feature

    if feature == "Pre-processing":
        return PREPROCESSING_MAP.get(iv, f"Unknown ({iv})")

    if feature == "Reactor Type":
        return REACTOR_TYPE_MAP.get(iv, f"Unknown ({iv})")

    if feature == "Catalyst":
        return CATALYST_MAP.get(iv, f"Unknown ({iv})")

    if feature == "Solvent":
        return SOLVENT_MAP.get(iv, f"Unknown ({iv})")

    return value


# -----------------------------
# Yield prediction helper
# -----------------------------
def predict_yields_wt_pct(rf_model, X_row: pd.DataFrame):
    """
    Returns (wt% array length 4) in the same output order your RF uses.
    Assumption: [Biocrude, Aqueous, Gas, Char] in wt%.
    """
    y_pred = rf_model.predict(X_row)[0]
    return np.array(y_pred, dtype=float)


def wt_pct_to_frac(y_wt_pct):
    y = np.array(y_wt_pct, dtype=float) / 100.0
    s = y.sum()
    if not np.isfinite(s) or s <= 0:
        raise ValueError(f"Non-physical RF yields: {y_wt_pct}")
    return y / s


# -----------------------------
# Main script
# -----------------------------
def main():
    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")

    # -------------------------------------------------
    # Load yield RF model (same as your working script)
    # -------------------------------------------------
    rf_model = load(os.path.join(output_dir, "rf_yield_model.joblib"))
    print("✅ RF yield model loaded")

    # -------------------------------------------------
    # Choose feedstock + composition -> dry wt%
    # -------------------------------------------------
    feedstock_id = "manure"  # change as needed
    wet_comp = get_feedstock_composition(feedstock_id)

    dry_frac = 1.0 - wet_comp["Water"]
    if dry_frac <= 0:
        raise ValueError(f"Bad water fraction for {feedstock_id}: {wet_comp['Water']}")

    carb_wt    = wet_comp["Carbohydrates"] / dry_frac * 100
    protein_wt = wet_comp["Proteins"]      / dry_frac * 100
    lipids_wt  = wet_comp["Lipids"]        / dry_frac * 100
    ash_wt     = wet_comp["Ash"]           / dry_frac * 100

    # -------------------------------------------------
    # Process conditions (keep this name)
    # You can vary these later; for now we’ll create samples
    # by changing categorical encodes + (optionally) temperature/time/solid%.
    # -------------------------------------------------
    base_process_conditions = {
        "Temperature (C)": 280,
        "Residence Time": 15,
        "Solid content (w/w) %": 20,
        # you can include these too if your FEATURE_COLS expects them:
        # "Pre-processing": 1,
        # "Catalyst": 1,
        # "Reactor Type": 0,
        # "Solvent": 1,
        # "Reactor Volume (mL)": ...
    }

    # -------------------------------------------------
    # Build a small DOE-style sample set (edit as you like)
    # Keep naming “process_conditions” per sample.
    # -------------------------------------------------
    # Categorical codes you gave:
    # Pre-processing: 0/1
    # Reactor Type: 0=batch, 1=PFR, 2=CSTR
    # Catalyst: many codes (use a subset for now)
    # Solvent: 1=no solvent, 2=50% ethanol, 3=seawater
    preprocess_codes = [0, 1]
    reactor_type_codes = [0, 1, 2]
    catalyst_codes = [1, 2, 3, 4]   # pick a few to start
    solvent_codes = [1, 2, 3]

    # Also vary continuous knobs lightly (optional)
    temperatures = [260, 280, 300]
    residence_times = [30, 45, 60]
    solid_loadings = [15, 20, 25]

    samples = []
    for pp in preprocess_codes:
        for rt in reactor_type_codes:
            for cat in catalyst_codes:
                for solv in solvent_codes:
                    for T in temperatures:
                        for tau in residence_times:
                            for sc in solid_loadings:
                                process_conditions = dict(base_process_conditions)
                                process_conditions.update({
                                    "Pre-processing": pp,
                                    "Reactor Type": rt,
                                    "Catalyst": cat,
                                    "Solvent": solv,
                                    "Temperature (C)": T,
                                    "Residence Time": tau,
                                    "Solid content (w/w) %": sc,
                                })
                                samples.append(process_conditions)

    print(f"🔧 Meta sampling plan: {len(samples)} samples")

    # -------------------------------------------------
    # Build meta-dataset: features -> predicted yields -> simulate -> MSP/GWP
    # -------------------------------------------------
    rows = []

    # Use ONE config (match your earlier scripts). Add more if you want.
    config_kwargs = {
        "flowsheet": None,
        "decentralized_HTL": False,
        "decentralized_upgrading": False,
        "skip_EC": True,
        "generate_H2": False,
        "EC_config": None,
        # include your flowrate args if your create_system needs them
        # "central_dry_flowrate": default_central,
    }

    for i, process_conditions in enumerate(samples, 1):
        # ---- feature row (same adapter as your yield SHAP script)
        X = make_feature_row(
            carb_wt=carb_wt,
            protein_wt=protein_wt,
            lipids_wt=lipids_wt,
            ash_wt=ash_wt,
            process=process_conditions,
        )
        print("DEBUG feature snapshot:",
                X[["Pre-processing","Reactor Type","Catalyst","Solvent",
               "Temperature (C)","Residence Time","Solid content (w/w) %"]].iloc[0].to_dict())

        # store raw features as dict
        feature_row = X.iloc[0].to_dict()

        try:
            # ---- predict yields (wt%) then to fractions
            y_wt = predict_yields_wt_pct(rf_model, X)
            y_frac = wt_pct_to_frac(y_wt)

            # ---- update global HTL_yields BEFORE create_system (your pattern)
            # Assumed order: [Biocrude, Aqueous, Gas, Char]
            HTL_yields.update({
                "biocrude": float(y_frac[0]),
                "aqueous":  float(y_frac[1]),
                "gas":      float(y_frac[2]),
                "char":     float(y_frac[3]),
            })
            if "feedstock_id" in config_kwargs:
               raise ValueError("feedstock_id should be passed via loop variable, not config_kwargs, to avoid confusion.")


            # ---- create + simulate
            print("Before create_system:", HTL_yields)
            sys = create_system(**config_kwargs,feedstock_id=feedstock_id, process_conditions=process_conditions)
            print("After create_system:", HTL_yields)
            tea = sys.TEA
            lca = sys.LCA

            tea.IRR = tea_kwargs["IRR"]
            tea.income_tax = tea_kwargs["income_tax"]

            sys.simulate()

            # ---- outputs (match your previous naming)
            biobinder = sys.flowsheet.stream.biobinder

            MSP = float(tea.solve_price(biobinder))

            # Keep your GWP definition consistent with your existing script
            impacts = lca.get_allocated_impacts(
                streams=(biobinder,),
                operation_only=True,
                annual=True,
            )
            GWP = float(impacts["GWP"]) / (biobinder.F_mass * lca.system.operating_hours)

            rows.append({
                "Feedstock": feedstock_id,
                "Timestamp": timestamp,
                "Sample_ID": i,
                **feature_row,

                "Pred_Biocrude_wt%": float(y_wt[0]),
                "Pred_Aqueous_wt%":  float(y_wt[1]),
                "Pred_Gas_wt%":      float(y_wt[2]),
                "Pred_Char_wt%":     float(y_wt[3]),

                "Pred_Biocrude_frac": float(y_frac[0]),
                "Pred_Aqueous_frac":  float(y_frac[1]),
                "Pred_Gas_frac":      float(y_frac[2]),
                "Pred_Char_frac":     float(y_frac[3]),

                "MSP": MSP,
                "GWP": GWP,
                "OK": True,
                "Error": "",
            })

        except Exception as e:
            rows.append({
                "Feedstock": feedstock_id,
                "Timestamp": timestamp,
                "Sample_ID": i,
                **feature_row,

                "MSP": np.nan,
                "GWP": np.nan,
                "OK": False,
                "Error": repr(e),
            })

        if i % 50 == 0:
            ok_ct = sum(r["OK"] for r in rows)
            print(f"  progress: {i}/{len(samples)} | OK={ok_ct}")

    df_meta = pd.DataFrame(rows)

    meta_path = os.path.join(output_dir, f"META_TEA_SHAP_{feedstock_id}_{timestamp}.xlsx")
    df_meta.to_excel(meta_path, index=False)
    print(f"📁 Saved meta dataset: {meta_path}")

    # ---- sanity
    ok = df_meta["OK"].astype(bool)
    print("df_meta:", df_meta.shape, "| OK:", ok.sum(), "/", len(df_meta))
    if ok.sum() == 0:
        print("Top errors:")
        print(df_meta["Error"].value_counts().head(15))
        raise RuntimeError("No successful simulations (OK=0). Fix errors above first.")

    # -------------------------------------------------
    # Train surrogates (MSP, GWP) on successful rows
    # -------------------------------------------------
    
    df_ok = df_meta[df_meta["OK"]].copy()
    
    # Only keep rows with finite targets
    df_ok = df_ok[np.isfinite(df_ok["MSP"]) & np.isfinite(df_ok["GWP"])].copy()
    MSP_CAP = 20.0
    df_ok["MSP_raw"] = df_ok["MSP"].astype(float)
    df_ok["MSP_cap"] = df_ok["MSP_raw"].clip(lower=-MSP_CAP, upper=MSP_CAP)
    
    if len(df_ok) < 10:
        raise RuntimeError(f"Too few successful rows to train surrogates: {len(df_ok)}")

    X_all = df_ok[FEATURE_COLS].copy()

    # ---- MSP surrogate
    y_msp = df_ok["MSP_cap"].values
    Xtr, Xte, ytr, yte = train_test_split(X_all, y_msp, test_size=0.2, random_state=42)
    rf_msp = RandomForestRegressor(n_estimators=800, random_state=42, n_jobs=-1)
    rf_msp.fit(Xtr, ytr)
    r2_msp = r2_score(yte, rf_msp.predict(Xte))
    print(f"✅ MSP surrogate trained (capped ±{MSP_CAP}) | R2(test) = {r2_msp:.3f}")

    # ---- GWP surrogate
    y_gwp = df_ok["GWP"].values
    Xtr, Xte, ytr, yte = train_test_split(X_all, y_gwp, test_size=0.2, random_state=42)
    rf_gwp = RandomForestRegressor(n_estimators=800, random_state=42, n_jobs=-1)
    rf_gwp.fit(Xtr, ytr)
    r2_gwp = r2_score(yte, rf_gwp.predict(Xte))
    print(f"✅ GWP surrogate trained | R2(test) = {r2_gwp:.3f}")

    # -------------------------------------------------
    # SHAP on surrogates
    # -------------------------------------------------
    # SHAP expects a matrix; we’ll compute per-row shap for df_ok
    expl_msp = TreeExplainer(rf_msp)
    shap_msp = expl_msp.shap_values(X_all)  # shape: (n_rows, n_features)

    expl_gwp = TreeExplainer(rf_gwp)
    shap_gwp = expl_gwp.shap_values(X_all)

    # -------------------------------------------------
    # Save SHAP tables (with feature values + decoded categorical)
    # -------------------------------------------------
    def save_surrogate_shap(shap_vals_2d, target_name, target_r2):
        # Build long-form table: one row per (sample, feature)
        out_rows = []
        for row_i in range(X_all.shape[0]):
            sample_id = int(df_ok.iloc[row_i]["Sample_ID"])
            feat_vals = X_all.iloc[row_i]
            for feat_j, feat in enumerate(FEATURE_COLS):
                raw_val = feat_vals[feat]
                if target_name == "MSP":
                    target_val = float(df_ok.iloc[row_i]["MSP_cap"])
                    target_raw = float(df_ok.iloc[row_i]["MSP_raw"])
                else:
                    target_val = float(df_ok.iloc[row_i][target_name])
                    target_raw = None

                row = {
                "Feedstock": feedstock_id,
                "Timestamp": timestamp,
                "Target": target_name,
                "Target_R2_test": float(target_r2),
                "Sample_ID": sample_id,
                "Feature": feat,
                "Feature value": raw_val,
                "Feature value (decoded)": decode_feature_value(feat, raw_val),
                f"{target_name}_SHAP": float(shap_vals_2d[row_i, feat_j]),
                f"{target_name}_value": target_val,
                  }
                if target_name == "MSP":
                    row["MSP_raw_value"] = target_raw
                    row["MSP_cap_value"] = target_val
                    row["MSP_cap_used"] = MSP_CAP
                out_rows.append(row)
                # out_rows.append({
                #     "Feedstock": feedstock_id,
                #     "Timestamp": timestamp,
                #     "Target": target_name,
                #     "Target_R2_test": float(target_r2),
                #     "Sample_ID": sample_id,
                #     "Feature": feat,
                #     "Feature value": raw_val,
                #     "Feature value (decoded)": decode_feature_value(feat, raw_val),
                #     f"{target_name}_SHAP": float(shap_vals_2d[row_i, feat_j]),
                #     # include actual target value for that sample
                #     f"{target_name}_value": float(df_ok.iloc[row_i][target_name]),
                # })

        df_out = pd.DataFrame(out_rows)

        path = os.path.join(output_dir, f"SURROGATE_SHAP_{feedstock_id}_{target_name}_{timestamp}.xlsx")
        df_out.to_excel(path, index=False)
        print(f"📁 Saved surrogate SHAP: {path}")

    save_surrogate_shap(np.array(shap_msp), "MSP", r2_msp)
    save_surrogate_shap(np.array(shap_gwp), "GWP", r2_gwp)

    print("✅ Done.")


if __name__ == "__main__":
    main()
