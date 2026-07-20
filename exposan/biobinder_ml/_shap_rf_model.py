# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import pandas as pd
import numpy as np
from shap import TreeExplainer
from joblib import load
from exposan.biobinder_ml._feedstocks import get_feedstock_composition
from exposan.biobinder_ml._ml_features import FEATURE_COLS, make_feature_row, PREPROCESSING_MAP
from exposan.biobinder_ml._ml_features import CATALYST_MAP, SOLVENT_MAP, REACTOR_TYPE_MAP
from datetime import datetime
import os
# -------------------------------------------------
# Load model
# -------------------------------------------------
rf_model = load("results/rf_yield_model.joblib")
explainer = TreeExplainer(rf_model)
print("RF model loaded")

# -------------------------------------------------
# Feedstock → dry-basis composition
# -------------------------------------------------
feedstock_id = "sludge"
wet_comp = get_feedstock_composition(feedstock_id)

dry_frac = 1.0 - wet_comp["Water"]

carb_wt    = wet_comp["Carbohydrates"] / dry_frac * 100
protein_wt = wet_comp["Proteins"]      / dry_frac * 100
lipids_wt  = wet_comp["Lipids"]        / dry_frac * 100
ash_wt     = wet_comp["Ash"]           / dry_frac * 100

# -------------------------------------------------
# Process conditions (must exist)
# -------------------------------------------------
process_conditions = {
    "Temperature (C)": 280,
    "Residence Time": 15,
    "Solid content (w/w) %": 20,
}

# -------------------------------------------------
# Feature row
# -------------------------------------------------
X = make_feature_row(
    carb_wt=carb_wt,
    protein_wt=protein_wt,
    lipids_wt=lipids_wt,
    ash_wt=ash_wt,
    process=process_conditions,
)


y = rf_model.predict(X)[0]   # wt%
y_pred = rf_model.predict(X)[0]

predicted_yields = {
    "Predicted biocrude wt%": float(y_pred[0]),
    "Predicted aqueous wt%":  float(y_pred[1]),
    "Predicted gas wt%":      float(y_pred[2]),
    "Predicted char wt%":     float(y_pred[3]),
}

# Optional: fractions
predicted_yields_frac = {
    "Predicted biocrude frac": predicted_yields["Predicted biocrude wt%"] / 100,
    "Predicted aqueous frac":  predicted_yields["Predicted aqueous wt%"] / 100,
    "Predicted gas frac":      predicted_yields["Predicted gas wt%"] / 100,
    "Predicted char frac":     predicted_yields["Predicted char wt%"] / 100,
}
# -------------------------------------------------
# SHAP
# -------------------------------------------------
shap_vals = explainer.shap_values(X)
feature_values = X.iloc[0]
shap_arr = shap_vals[0]  # shape: (n_features, 4)
timestamp = datetime.now().strftime("%Y%m%d_%H%M")
feedstock = feedstock_id  # already defined earlier

output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

output_names = [
    "Biocrude",
    "Aqueous",
    "Gas",
    "Char",
]

# -------------------------------------------------
# Build SHAP DataFrame
# -------------------------------------------------
shap_out = pd.DataFrame(
    shap_arr,
    index=FEATURE_COLS,
    columns=[f"{o}_SHAP" for o in output_names],
)

# -------------------------------------------------
# Insert feature name + raw feature value
# -------------------------------------------------
shap_out.insert(0, "Feature", shap_out.index)
shap_out.insert(1, "Feature value", feature_values.values)

# -------------------------------------------------
# Decode categorical feature values
# -------------------------------------------------
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

shap_out.insert(
    2,
    "Feature value (decoded)",
    [
        decode_feature_value(f, v)
        for f, v in zip(shap_out["Feature"], shap_out["Feature value"])
    ],
)

# -------------------------------------------------
# Optional rounding for readability
# -------------------------------------------------
shap_out = shap_out.round(4)

# -------------------------------------------------
# Add metadata
# -------------------------------------------------
shap_out.insert(0, "Feedstock", feedstock)
shap_out.insert(1, "Timestamp", timestamp)

for col, val in predicted_yields.items():
    shap_out[col] = val

for col, val in predicted_yields_frac.items():
    shap_out[col] = round(val, 4)

shap_out.reset_index(drop=True, inplace=True)

# -------------------------------------------------
# Save
# -------------------------------------------------
base = f"SHAP_RF_{feedstock}_{timestamp}"
# csv_path  = os.path.join(output_dir, f"{base}.csv")
xlsx_path = os.path.join(output_dir, f"{base}.xlsx")

# shap_out.to_csv(csv_path, index=False)
shap_out.to_excel(xlsx_path, index=False)

print("✅ SHAP + decoded feature values saved")
# print(f"   → {csv_path}")
print(f"   → {xlsx_path}")