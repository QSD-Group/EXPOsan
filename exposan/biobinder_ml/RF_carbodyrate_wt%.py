# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 17:03:51 2025

@author: aliah
"""
import os as os
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import r2_score, mean_absolute_error

output_dir = "results"
os.makedirs(output_dir, exist_ok=True)
# Load your dataset
data = pd.read_excel(os.path.join(output_dir,"HTL_all_yield_with_estimated_protein.xlsx"))


# --- Derived ratios ---
data["H/C"] = data["H%"] / data["C%"]
data["O/C"] = data["O%"] / data["C%"]
data["N/C"] = data["N%"] / data["C%"]

# --- Features and target ---
features = ["C%", "H%", "O%", "N%", "S%", "H/C", "O/C", "N/C"]
target = "Carbohydrates wt%"

# Drop rows with missing target for training
train_df = data.dropna(subset=[target])[features + [target]].dropna()
X_train = train_df[features]
y_train = train_df[target]

# Train Random Forest (Xing-style)
rf = RandomForestRegressor(
    n_estimators=300, max_depth=12, random_state=42, n_jobs=-1
)
rf.fit(X_train, y_train)

# --- Predict for rows with missing Carbohydrates ---
mask_missing = data[target].isna() & data[features].notna().all(axis=1)
X_missing = data.loc[mask_missing, features]
predictions = rf.predict(X_missing)

# --- Fill the missing values ---
data.loc[mask_missing, target] = predictions

# --- Save updated dataset ---
data.to_excel("HTL_all_yield_with_RF_imputed_carbohydrates.xlsx", index=False)

print(f"✅ Imputed {len(predictions)} rows successfully.")
print("Saved as: HTL_all_yield_with_RF_imputed_carbohydrates.xlsx")
