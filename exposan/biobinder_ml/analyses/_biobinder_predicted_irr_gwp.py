# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import os
import numpy as np
import pandas as pd
from datetime import datetime
from joblib import load
import qsdsan as qs
from chaospy import distributions as shape
import numpy as np
import pandas as pd
from scipy.stats import qmc
import os
from datetime import datetime
from exposan.biobinder_ml import (
    # create_system,
    _HHV_per_GGE,
    _load_components,
    _load_process_settings,
    _units as u,
    BiobinderTEA,
    central_dry_flowrate as default_central,
    data_path,
    feedstock_composition,
    HTL_yields,
    pilot_dry_flowrate as default_pilot,
    price_dct,
    results_path,
    tea_kwargs,
    uptime_ratio,
    limit_internal_hx,
    configure_utilities,
    set_eff_hx_temperature,
    )
from exposan.biobinder_ml.Dist_flex import (
    create_system)

# Create a directory for results
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

# yield_table_path = os.path.join(output_dir, "HTL_all_yield_normalized.xlsx")  # adjust if needed
# yt = pd.read_excel(yield_table_path, sheet_name="Sheet1")
# yt.columns = yt.columns.str.strip()

model_path = os.path.join(output_dir, "rf_yield_model.joblib")
data_path  = os.path.join(output_dir, "HTL_all_yield_experimental.xlsx")

print(f"🔍 Loading RF model from: {model_path}")
rf_model = load(model_path)

# ============================================================
# 2️⃣ LOAD DATA
# ============================================================
df = pd.read_excel(data_path, sheet_name="Sheet1")
expected_features = rf_model.feature_names_in_

# Align columns and encode categorical variables
X = df.reindex(columns=expected_features, fill_value=0)
for col in X.columns:
    if X[col].dtype == "object":
        X[col] = pd.factorize(X[col])[0]
X = X.fillna(X.median())

# ============================================================
# 3️⃣ PREDICT YIELDS
# ============================================================
pred_yields = rf_model.predict(X)
pred_df = pd.DataFrame(pred_yields, columns=["Biocrude wt%", "Aqueous wt%", "Gas wt%", "Solids wt%"])
pred_df["Sum wt%"] = pred_df.sum(axis=1)

mask = (pred_df["Sum wt%"] >= 95) & (pred_df["Sum wt%"] <= 105)
valid_df = pred_df[mask].copy()
invalid_df = pred_df[~mask].copy()
df_pred = pd.concat([df, pred_df], axis=1)

print(f"✅ Predicted yields for {len(valid_df)} valid samples (±5%)")

# ============================================================
# 4️⃣ RUN SIMULATIONS FOR EACH CONFIGURATION
# ============================================================
configs = [
    {
        "name": "CHCU_No_EC",
        "config_kwargs": {
            "flowsheet": None,
            "central_dry_flowrate": default_central,
            "decentralized_HTL": False,
            "decentralized_upgrading": False,
            "skip_EC": True,
            "generate_H2": False,
            "EC_config": None,
        },
    },
]

timestamp = datetime.now().strftime("%Y%m%d_%H%M")
all_results = []

for config in configs:
    cfg_name = config["name"]
    print(f"\n⚙️ Running configuration: {cfg_name}")

    for i, row in valid_df.iterrows():
        try:
            # Safely extract scalar values from predicted yields
            Y_bio  = float(row["Biocrude wt%"]) / 100
            Y_aq   = float(row["Aqueous wt%"]) / 100
            Y_gas  = float(row["Gas wt%"]) / 100
            Y_char = float(row["Solids wt%"]) / 100

            # Update global yield dictionary
            HTL_yields.update({
                "biocrude": Y_bio,
                "gas":      Y_gas,
                "char":     Y_char,
                "aqueous":  Y_aq,
            })

            # Create and simulate the system
            sys = create_system(**config["config_kwargs"])
            tea = sys.TEA
            lca = sys.LCA
            tea.IRR = tea_kwargs["IRR"]
            tea.income_tax = tea_kwargs["income_tax"]

            sys.simulate()

            biobinder = sys.flowsheet.stream.biobinder
            biofuel   = sys.flowsheet.stream.biofuel
            MSP  = tea.solve_price(biobinder)
            IRR  = tea.solve_IRR() * 100
            GWP  = lca.get_allocated_impacts(
                streams=(biobinder,), operation_only=True, annual=True
            )["GWP"] / (biobinder.F_mass * lca.system.operating_hours)
            ratio = biobinder.F_mass / biofuel.F_mass if biofuel.F_mass else np.nan

            all_results.append({
                "Config": cfg_name,
                "Index": i,
                "Biocrude wt%": row["Biocrude wt%"],
                "Aqueous wt%": row["Aqueous wt%"],
                "Gas wt%": row["Gas wt%"],
                "Solids wt%": row["Solids wt%"],
                "MSP ($/kg)": MSP,
                "IRR (%)": IRR,
                "GWP (kg CO2e/kg)": GWP,
                "Product Ratio": ratio
            })

            print(f"   ✅ Sample {i}: MSP=${MSP:.2f}/kg | IRR={IRR:.1f}% | GWP={GWP:.3f}")

        except Exception as e:
            all_results.append({
                "Config": cfg_name,
                "Index": i,
                "Biocrude wt%": float(row['Biocrude wt%']),
                "Aqueous wt%": float(row['Aqueous wt%']),
                "Gas wt%": float(row['Gas wt%']),
                "Solids wt%": float(row['Solids wt%']),
                "MSP ($/kg)": np.nan,
                "IRR (%)": np.nan,
                "GWP (kg CO2e/kg)": np.nan,
                "Error": str(e)
            })
            print(f"   ❌ Sample {i} failed: {e}")

# ============================================================
# 5️⃣ SAVE RESULTS
# ============================================================
out_path = os.path.join(output_dir, f"RF_Yield_TEA_Results_{timestamp}.xlsx")
with pd.ExcelWriter(out_path) as w:
    pd.DataFrame(all_results).to_excel(w, sheet_name="results", index=False)
    invalid_df.to_excel(w, sheet_name="skipped", index=False)

print(f"\n📁 All results saved to: {out_path}")
