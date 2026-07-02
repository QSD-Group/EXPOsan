# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 18:07:15 2025

@author: aliah
"""


import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)
import xgboost as xgb
# requirements: pip install xgboost shap scikit-learn pandas openpyxl numpy
import pandas as pd, numpy as np
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import ElasticNet
from sklearn.multioutput import MultiOutputRegressor
from xgboost import XGBRegressor
import shap
import os

# ---- load data ----


output_dir = "results"
os.makedirs(output_dir, exist_ok=True)
# Load your dataset
file_path = os.path.join(output_dir, "HTL_all_yield_normalized_split_stratified_10pct.xlsx")

# Load the split sheets
train_df = pd.read_excel(file_path, sheet_name="Train")
test_df  = pd.read_excel(file_path, sheet_name="Test")

# Define target columns
yield_cols = ["Biocrude wt%", "Aqueous wt%", "Gas wt%", "Solids wt%"]

# Separate features and targets
Xtr = train_df.drop(columns=yield_cols).select_dtypes(include=[np.number])
ytr = train_df[yield_cols]
Xte = test_df.drop(columns=yield_cols).select_dtypes(include=[np.number])
yte = test_df[yield_cols]

# Impute any residual missing values
imp = SimpleImputer(strategy="median")
Xtr = pd.DataFrame(imp.fit_transform(Xtr), columns=Xtr.columns)
Xte = pd.DataFrame(imp.transform(Xte), columns=Xte.columns)

# ===============================================================
# 2. Define helpers for metrics + coverage
# ===============================================================

def coverage(y_true, y_pred, thr):
    """Return per-target and overall % of samples within ±thr of true values."""
    abs_err = np.abs(y_true - y_pred)
    return (abs_err < thr).mean(axis=0) * 100.0, (abs_err < thr).mean() * 100.0


def metrics_table(name, y_true, y_pred, k_params):
    """Compute model metrics (R², RMSE, MAE, MAPE, AIC, Coverage)."""
    n = y_true.shape[0]
    rows = []
    
    # Compute coverage thresholds
    cov5_t, cov5 = coverage(y_true.values, y_pred, 5.0)
    cov10_t, cov10 = coverage(y_true.values, y_pred, 10.0)
    
    # Loop over each yield target
    for i, c in enumerate(y_true.columns):
        yt, yp = y_true.iloc[:, i].values, y_pred[:, i]
        mse = mean_squared_error(yt, yp)
        rows.append(dict(
            Model=name,
            Target=c,
            R2=r2_score(yt, yp),
            RMSE=np.sqrt(mse),
            MAE=mean_absolute_error(yt, yp),
            MAPE=np.nanmean(np.abs((yt - yp) / np.where(yt == 0, np.nan, yt))) * 100.0,
            AIC_approx=n * np.log(mse + 1e-12) + 2 * k_params,
            **{"Coverage_<5p": cov5_t[i], "Coverage_<10p": cov10_t[i]}
        ))

    dfm = pd.DataFrame(rows)
    
    # Macro-average row (overall)
    macro = dfm.mean(numeric_only=True).to_dict()
    macro.update(Model=name, Target="__macro__")
    macro["Coverage_<5p"] = cov5
    macro["Coverage_<10p"] = cov10

    return pd.concat([dfm, pd.DataFrame([macro])], ignore_index=True)


# ===============================================================
# 3. Random Forest
# ===============================================================
rf = RandomForestRegressor(
    n_estimators=600, max_depth=None,
    max_features="sqrt", random_state=42, n_jobs=-1
)
rf.fit(Xtr, ytr)
rf_pred = rf.predict(Xte)
rf_metrics = metrics_table("RandomForest", yte, rf_pred, k_params=Xtr.shape[1])

rf_importance = pd.Series(rf.feature_importances_, index=Xtr.columns)\
    .sort_values(ascending=False).reset_index()\
    .rename(columns={"index": "Feature", 0: "Importance"})

# ===============================================================
# 4. XGBoost
# ===============================================================
xgb = MultiOutputRegressor(XGBRegressor(
    n_estimators=700, learning_rate=0.05, max_depth=4,
    subsample=0.9, colsample_bytree=0.9, reg_lambda=1.0,
    random_state=42, n_jobs=-1, tree_method="hist", verbosity=0
))
xgb.fit(Xtr, ytr)
xgb_pred = xgb.predict(Xte)
xgb_metrics = metrics_table("XGBoost", yte, xgb_pred, k_params=Xtr.shape[1])

# Average feature importance across all 4 yield outputs
xgb_imp = np.array([est.feature_importances_ for est in xgb.estimators_]).mean(axis=0)
xgb_importance = pd.DataFrame({"Feature": Xtr.columns, "Importance": xgb_imp})\
    .sort_values("Importance", ascending=False)

# ===============================================================
# 5. Save results
# ===============================================================
with pd.ExcelWriter("Results_RF_XGB_TrainTest_2.xlsx", engine="openpyxl") as w:
    rf_metrics.to_excel(w, sheet_name="RF_Metrics", index=False)
    rf_importance.to_excel(w, sheet_name="RF_FeatureImportance", index=False)
    xgb_metrics.to_excel(w, sheet_name="XGB_Metrics", index=False)
    xgb_importance.to_excel(w, sheet_name="XGB_FeatureImportance", index=False)

print("✅ Training complete — results saved as Results_RF_XGB_TrainTest.xlsx")
