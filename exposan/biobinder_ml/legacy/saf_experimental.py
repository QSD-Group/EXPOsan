# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 13:41:21 2025

@author: aliah
"""

import os
import numpy as np
import pandas as pd
from datetime import datetime
from joblib import load
from qsdsan.utils import clear_lca_registries
from exposan.saf import (
    create_system,
    config_baseline,
    HTL_yields,
    tea_kwargs,
)
from exposan.biobinder_ml import (
    cutoff_fracs_from_feed,
    limit_internal_hx,
    configure_utilities,
    set_eff_hx_temperature,
    find_Lr_Hr
    )

# ============================================================
# 1️⃣ PATHS AND LOADING MODEL + DATA
# ============================================================
# Create a directory for results
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)
yield_table_path = os.path.join(output_dir, "HTL_all_yield_normalized.xlsx")
yt = pd.read_excel(yield_table_path, sheet_name="Sheet1")
yt.columns = yt.columns.str.strip()

required_cols = [ "Solids wt%", "Biocrude wt%", "Aqueous wt%", "Gas wt%"]
missing_req = [c for c in required_cols if c not in yt.columns]
if missing_req:
    raise ValueError(f"Missing required columns in HTL_yield_type.xlsx: {missing_req}")

# Drop rows fully empty in yield columns, convert wt% -> fractions
yt = yt.dropna(subset=["Solids wt%", "Biocrude wt%", "Aqueous wt%", "Gas wt%"], how="all").copy()
for col in ["Solids wt%", "Biocrude wt%", "Aqueous wt%", "Gas wt%"]:
    yt[col] = pd.to_numeric(yt[col], errors="coerce") / 100.0

# Sum check with ±5% tolerance (no normalization; skip invalid rows)
yt["Sum_fractions"] = yt[["Solids wt%", "Biocrude wt%", "Aqueous wt%", "Gas wt%"]].sum(axis=1)
tolerance = 0.05
valid_mask = (yt["Sum_fractions"] >= 1 - tolerance) & (yt["Sum_fractions"] <= 1 + tolerance)
valid_rows = yt.loc[valid_mask].copy()
skipped_rows = yt.loc[~valid_mask].copy()

# Prepare fraction columns for clarity
valid_rows["Solids_f"]   = valid_rows["Solids wt%"]
valid_rows["Biocrude_f"] = valid_rows["Biocrude wt%"]
valid_rows["Aqueous_f"]  = valid_rows["Aqueous wt%"]
valid_rows["Gas_f"]      = valid_rows["Gas wt%"]

# Assign unique index for traceability
valid_rows = valid_rows.reset_index(drop=True)
valid_rows["Sample_ID"] = valid_rows.index + 1  # start from 1 for readability

timestamp = datetime.now().strftime("%Y%m%d_%H%M")
# ============================================================
# 2️⃣ PREDICT YIELDS
# ============================================================
# pred_yields = rf_model.predict(X)
# pred_df = pd.DataFrame(pred_yields, columns=["Biocrude wt%", "Aqueous wt%", "Gas wt%", "Solids wt%"])
# pred_df["Sum wt%"] = pred_df.sum(axis=1)

# # Keep only valid (sum ≈ 100%)
# mask = (pred_df["Sum wt%"] >= 95) & (pred_df["Sum wt%"] <= 105)
# valid_df = pred_df[mask].copy()
# invalid_df = pred_df[~mask].copy()

# df_pred = pd.concat([df, pred_df], axis=1)
# print(f"✅ Predicted yields for {len(valid_df)} valid samples (±5%)")

# ============================================================
# 3️⃣ RUN SAF SYSTEM SIMULATIONS
# ============================================================
timestamp = datetime.now().strftime("%Y%m%d_%H%M")
all_results = []

print(f"\n⚙️ Running SAF system using predicted yields ({len(valid_rows)} samples)...")

for i, row in valid_rows.iterrows():
    try:
        # Explicit scalar conversion to avoid Series ambiguity
        Y_bio  = float(row["Biocrude wt%"]) / 100
        Y_aq   = float(row["Aqueous wt%"]) / 100
        Y_gas  = float(row["Gas wt%"]) / 100
        Y_char = float(row["Solids wt%"]) / 100

        # Update SAF HTL yields
        HTL_yields.update({
            "biocrude": Y_bio,
            "gas":      Y_gas,
            "char":     Y_char,
            "aqueous":  Y_aq,
        })

        # Clear registry & rebuild fresh system
        clear_lca_registries()
        sys = create_system(flowsheet=None, **config_baseline)
        tea = sys.TEA
        lca = sys.LCA
        tea.IRR = tea_kwargs["IRR"]
        tea.income_tax = tea_kwargs["income_tax"]

        sys.simulate()

        # Compute MFSP, IRR, GWP per GGE
        mixed_fuel = sys.flowsheet.stream.mixed_fuel
        MFSP = sys.TEA.solve_price(mixed_fuel)
        mixed_fuel.price= 4.87 #$/GGE
        IRR  = tea.solve_IRR() * 100
        all_impacts = lca.get_allocated_impacts(streams=(mixed_fuel,), operation_only=True, annual=True)
        total_GWP = all_impacts["GWP"]

        HHV_per_GGE = 120.21  # MJ per GGE
        GGE_annual = mixed_fuel.HHV / 1e3 / HHV_per_GGE * sys.operating_hours
        GWP_per_GGE = total_GWP / GGE_annual if GGE_annual else np.nan

        all_results.append({
            "Index": i,
            "Biocrude wt%": row["Biocrude wt%"],
            "Aqueous wt%": row["Aqueous wt%"],
            "Gas wt%": row["Gas wt%"],
            "Solids wt%": row["Solids wt%"],
            "MFSP ($/GGE)": MFSP,
            "IRR (%)": IRR,
            "GWP (kg CO2e/GGE)": GWP_per_GGE,
        })

        print(f"   ✅ Sample {i}: MFSP=${MFSP:.2f}/GGE | IRR={IRR:.1f}% | GWP={GWP_per_GGE:.3f}")

    except Exception as e:
        all_results.append({
            "Index": i,
            "Biocrude wt%": float(row["Biocrude wt%"]),
            "Aqueous wt%": float(row["Aqueous wt%"]),
            "Gas wt%": float(row["Gas wt%"]),
            "Solids wt%": float(row["Solids wt%"]),
            "MFSP ($/GGE)": np.nan,
            "IRR (%)": np.nan,
            "GWP (kg CO2e/GGE)": np.nan,
            "Error": str(e),
        })
        print(f"   ❌ Sample {i} failed: {e}")

# ============================================================
# 4️⃣ SAVE RESULTS
# ============================================================
out_path = os.path.join(output_dir, f"SAF_Yield_TEA_Results_{timestamp}.xlsx")
with pd.ExcelWriter(out_path) as w:
    pd.DataFrame(all_results).to_excel(w, sheet_name="results", index=False)
    skipped_rows.to_excel(w, sheet_name="skipped", index=False)

print(f"\n📁 All results saved to: {out_path}")
