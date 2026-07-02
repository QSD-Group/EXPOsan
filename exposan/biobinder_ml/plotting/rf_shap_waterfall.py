# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 13:25:13 2026

@author: aliah
"""


# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Ali Ahmad <aa3056@scarletmail.rutgers.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
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
feedstock_id = "food"
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

y = rf_model.predict(X)[0]
y_pred = rf_model.predict(X)[0]

predicted_yields = {
    "Predicted biocrude wt%": float(y_pred[0]),
    "Predicted aqueous wt%":  float(y_pred[1]),
    "Predicted gas wt%":      float(y_pred[2]),
    "Predicted char wt%":     float(y_pred[3]),
}

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
timestamp = datetime.now().strftime("%Y%m%d_%H%M")
feedstock = feedstock_id

output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

output_names = [
    "Biocrude",
    "Aqueous",
    "Gas",
    "Char",
]

# robust handling for SHAP output
if isinstance(shap_vals, list):
    # older shape: list[n_outputs] of arrays (n_samples, n_features)
    shap_arr = np.column_stack([np.array(v)[0, :] for v in shap_vals])
else:
    shap_vals_arr = np.array(shap_vals)
    if shap_vals_arr.ndim == 3:
        # modern shape: (n_samples, n_features, n_outputs)
        shap_arr = shap_vals_arr[0]
    else:
        raise ValueError(f"Unexpected SHAP shape: {shap_vals_arr.shape}")

# -------------------------------------------------
# Waterfall plot for one selected output
# -------------------------------------------------
output_idx = 0   # 0=Biocrude, 1=Aqueous, 2=Gas, 3=Char
output_name = output_names[output_idx]

local_shap = shap_arr[:, output_idx]
base_values = np.array(explainer.expected_value).reshape(-1)
base_value = float(base_values[output_idx])
pred_value = float(y_pred[output_idx])

pretty_feature_names = [
    "Carbohydrates wt%" if f == "Carbohydrates" else
    "Proteins wt%" if f == "Proteins" else
    "Lipids wt%" if f == "Lipids" else
    "Ash wt%" if f == "Ash" else
    "Temperature (°C)" if f == "Temperature (C)" else
    "Residence Time (min)" if f == "Residence Time" else
    "Solid content (w/w) %" if f == "Solid content (w/w) %" else
    "Pre-processing" if f == "Pre-processing" else
    "Catalyst" if f == "Catalyst" else
    "Solvent" if f == "Solvent" else
    "Reactor Type" if f == "Reactor Type" else
    "Reactor Volume (mL)" if f == "Reactor Volume (mL)" else
    "C%" if f == "C%" else
    "H%" if f == "H%" else
    "N%" if f == "N%" else
    "S%" if f == "S%" else
    "O%" if f == "O%" else
    f
    for f in FEATURE_COLS
]

def decode_feature_value(feature, value):
    try:
        iv = int(value)
        is_int_like = abs(float(value) - iv) < 1e-9
    except Exception:
        is_int_like = False
        iv = value

    if feature == "Pre-processing" and is_int_like:
        return {0: "No", 1: "Yes"}.get(iv, f"Unknown ({iv})")

    if feature == "Catalyst" and is_int_like:
        return CATALYST_MAP.get(iv, f"Unknown ({iv})")

    if feature == "Solvent" and is_int_like:
        return SOLVENT_MAP.get(iv, f"Unknown ({iv})")

    if feature == "Reactor Type" and is_int_like:
        return REACTOR_TYPE_MAP.get(iv, f"Unknown ({iv})")

    return value

def format_display_value(feature, value):
    decoded = decode_feature_value(feature, value)
    if isinstance(decoded, (int, float, np.integer, np.floating)):
        return f"{float(decoded):.2f}"
    return str(decoded)

decoded_feature_values = [
    format_display_value(f, v)
    for f, v in zip(FEATURE_COLS, feature_values.values)
]

exp = shap.Explanation(
    values=local_shap,
    base_values=base_value,
    data=feature_values.values,              # keep numeric for SHAP math
    display_data=np.array(decoded_feature_values, dtype=object),  # shown on plot
    feature_names=pretty_feature_names,      # name only, no value here
)

plt.figure(figsize=(11, 7))
shap.plots.waterfall(exp, max_display=10, show=False)
plt.title(f"{output_name} Waterfall ({feedstock} waste)", fontsize=16)
plt.tight_layout()
# plt.figure(figsize=(12.5, 7.5))
# shap.plots.waterfall(exp, max_display=10, show=False)

# ax = plt.gca()

# # -------------------------------------------------
# # Layout: more room on the left for long labels
# # -------------------------------------------------
# plt.subplots_adjust(left=0.42, right=0.98, top=0.90, bottom=0.12)

# # -------------------------------------------------
# # Bold title
# # -------------------------------------------------
# ax.set_title("Biocrude Waterfall", fontsize=20, fontweight="bold", pad=12)

# # -------------------------------------------------
# # Bold x-axis ticks
# # -------------------------------------------------
# ax.tick_params(axis='x', labelsize=16, width=1.5)
# for lab in ax.get_xticklabels():
#     lab.set_fontweight('bold')

# # -------------------------------------------------
# # Bold y-axis feature labels
# # -------------------------------------------------
# ax.tick_params(axis='y', labelsize=16)
# for lab in ax.get_yticklabels():
#     lab.set_fontweight('bold')

# # -------------------------------------------------
# # Move left-side feature text farther left and keep visible
# # This helps long values like Reactor Volume, O%, C%, etc.
# # -------------------------------------------------
# for txt in ax.texts:
#     s = txt.get_text()

#     # top-right prediction text
#     if "f(x)" in s:
#         txt.set_fontweight("bold")
#         txt.set_fontsize(18)

#     # bottom expected value text
#     elif "E[f(X)]" in s:
#         txt.set_fontweight("bold")
#         txt.set_fontsize(18)

#     # left-side feature/value labels
#     elif " = " in s:
#         txt.set_fontsize(18)
#         txt.set_fontweight("bold")
#         txt.set_clip_on(False)

#         # shift a bit further left so text sits outside bars
#         x, y = txt.get_position()
#         txt.set_x(x - 0.25)

# # -------------------------------------------------
# # Make contribution labels on bars bold too
# # -------------------------------------------------
# for txt in ax.texts:
#     s = txt.get_text().strip()
#     if s.startswith("+") or s.startswith("−") or s.startswith("-"):
#         txt.set_fontweight("bold")
#         txt.set_fontsize(16)

# plt.tight_layout()
waterfall_path = os.path.join(output_dir, f"WATERFALL_{feedstock}_{output_name}_{timestamp}.png")
plt.savefig(waterfall_path, dpi=300, bbox_inches="tight")
plt.close()

print("✅ Waterfall plot saved")
print(f"   → {waterfall_path}")
print(f"   Base value ({output_name}): {base_value:.4f}")
print(f"   Prediction ({output_name}): {pred_value:.4f}")
print(f"   Check base + shap sum: {(base_value + local_shap.sum()):.4f}")

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

# # -------------------------------------------------
# # Decode categorical feature values
# # -------------------------------------------------
# def decode_feature_value(feature, value):
#     # continuous features stay numeric
#     try:
#         iv = int(value)
#         is_int_like = abs(float(value) - iv) < 1e-9
#     except Exception:
#         is_int_like = False
#         iv = value

#     if feature == "Pre-processing" and is_int_like:
#         # adjust this mapping if your encoding is different
#         return {0: "No", 1: "Yes"}.get(iv, f"Unknown ({iv})")

#     if feature == "Reactor Type" and is_int_like:
#         return REACTOR_TYPE_MAP.get(iv, f"Unknown ({iv})")

#     if feature == "Catalyst" and is_int_like:
#         return CATALYST_MAP.get(iv, f"Unknown ({iv})")

#     if feature == "Solvent" and is_int_like:
#         return SOLVENT_MAP.get(iv, f"Unknown ({iv})")

#     return value

# shap_out.insert(
#     2,
#     "Feature value (decoded)",
#     [
#         decode_feature_value(f, v)
#         for f, v in zip(shap_out["Feature"], shap_out["Feature value"])
#     ],
# )

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
# base = f"SHAP_RF_{feedstock}_{timestamp}"
# xlsx_path = os.path.join(output_dir, f"{base}.xlsx")

# shap_out.to_excel(xlsx_path, index=False)

# print("✅ SHAP + decoded feature values saved")
# print(f"   → {xlsx_path}")