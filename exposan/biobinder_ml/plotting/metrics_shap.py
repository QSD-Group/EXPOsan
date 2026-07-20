# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 12:02:30 2025

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
df = pd.read_excel(os.path.join(output_dir,"HTL_all_yield_normalized.xlsx"))

yield_cols = ["Biocrude wt%","Aqueous wt%","Gas wt%","Solids wt%"]
X = df.drop(columns=yield_cols).select_dtypes(include=[np.number]).copy()
y = df[yield_cols].copy()

# Impute numeric
X = pd.DataFrame(SimpleImputer(strategy="median").fit_transform(X), columns=X.columns)

# Split
Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=42)

# ==============================================================
# Helper functions
# ==============================================================
def coverage(y_true, y_pred, thr):
    abs_err = np.abs(y_true - y_pred)
    return (abs_err < thr).mean(axis=0)*100.0, (abs_err < thr).mean()*100.0

def metrics_table(name, y_true, y_pred, k_params):
    n = y_true.shape[0]
    rows=[]
    for i,c in enumerate(y_true.columns):
        yt, yp = y_true.iloc[:,i].values, y_pred[:,i]
        mse = mean_squared_error(yt, yp)
        rows.append(dict(
            Model=name, Target=c,
            R2=r2_score(yt, yp),
            RMSE=np.sqrt(mse),
            MAE=mean_absolute_error(yt, yp),
            MAPE=np.nanmean(np.abs((yt-yp)/np.where(yt==0,np.nan,yt)))*100.0,
            AIC_approx=n*np.log(mse+1e-12)+2*k_params
        ))
    dfm = pd.DataFrame(rows)
    macro = dfm.mean(numeric_only=True).to_dict()
    macro.update(Model=name, Target="__macro__")
    return pd.concat([dfm, pd.DataFrame([macro])], ignore_index=True)

sheets = {}

# ==============================================================
# Random Forest
# ==============================================================
rf = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
rf.fit(Xtr, ytr)
rf_pred = rf.predict(Xte)

cov5_t, cov5 = coverage(yte.values, rf_pred, 5.0)
cov10_t, cov10 = coverage(yte.values, rf_pred, 10.0)

rf_metrics = metrics_table("RandomForest", yte, rf_pred, k_params=Xtr.shape[1])
rf_metrics.loc[rf_metrics["Target"]!="__macro__", "Coverage_<5%"] = cov5_t
rf_metrics.loc[rf_metrics["Target"]!="__macro__", "Coverage_<10%"] = cov10_t
rf_metrics.loc[rf_metrics["Target"]=="__macro__", ["Coverage_<5%","Coverage_<10%"]] = [cov5, cov10]

sheets["RF_Metrics"] = rf_metrics
sheets["RF_Feature_Importance"] = pd.Series(rf.feature_importances_, index=Xtr.columns)\
    .sort_values(ascending=False).reset_index().rename(columns={"index":"Feature",0:"Importance"})

# ==============================================================
# XGBoost
# ==============================================================
xgb = MultiOutputRegressor(XGBRegressor(
    n_estimators=600, learning_rate=0.05, max_depth=4,
    subsample=0.9, colsample_bytree=0.9, reg_lambda=1.0,
    random_state=42, n_jobs=-1, tree_method="hist", verbosity=0
))
xgb.fit(Xtr, ytr)
xgb_pred = xgb.predict(Xte)

cov5_t, cov5 = coverage(yte.values, xgb_pred, 5.0)
cov10_t, cov10 = coverage(yte.values, xgb_pred, 10.0)

xgb_metrics = metrics_table("XGBoost", yte, xgb_pred, k_params=Xtr.shape[1])
xgb_metrics.loc[xgb_metrics["Target"]!="__macro__", "Coverage_<5%"] = cov5_t
xgb_metrics.loc[xgb_metrics["Target"]!="__macro__", "Coverage_<10%"] = cov10_t
xgb_metrics.loc[xgb_metrics["Target"]=="__macro__", ["Coverage_<5%","Coverage_<10%"]] = [cov5, cov10]

sheets["XGB_Metrics"] = xgb_metrics

# Average feature importances
xgb_imp = np.array([est.feature_importances_ for est in xgb.estimators_]).mean(axis=0)
sheets["XGB_Feature_Importance"] = pd.DataFrame({"Feature":Xtr.columns,"Importance":xgb_imp})\
    .sort_values("Importance", ascending=False)

# ==============================================================
# Elastic Net
# ==============================================================
enet = MultiOutputRegressor(ElasticNet(alpha=1e-3, l1_ratio=0.5, random_state=42, max_iter=10000))
enet.fit(Xtr, ytr)
enet_pred = enet.predict(Xte)

cov5_t, cov5 = coverage(yte.values, enet_pred, 5.0)
cov10_t, cov10 = coverage(yte.values, enet_pred, 10.0)

enet_metrics = metrics_table("ElasticNet", yte, enet_pred, k_params=Xtr.shape[1])
enet_metrics.loc[enet_metrics["Target"]!="__macro__", "Coverage_<5%"] = cov5_t
enet_metrics.loc[enet_metrics["Target"]!="__macro__", "Coverage_<10%"] = cov10_t
enet_metrics.loc[enet_metrics["Target"]=="__macro__", ["Coverage_<5%","Coverage_<10%"]] = [cov5, cov10]

sheets["ENet_Metrics"] = enet_metrics

# ==============================================================
# SHAP Analysis
# ==============================================================
sub_n = min(len(Xte), 200)
Xte_sub = Xte.iloc[:sub_n, :]

# --- Random Forest SHAP (new API) ---
try:
    expl_rf = shap.Explainer(rf, Xtr)
    sv_rf = expl_rf(Xte_sub)
    shap_df_rf = pd.DataFrame({
        "Feature": Xte_sub.columns,
        "Mean|SHAP|": np.abs(sv_rf.values).mean(axis=0)
    }).sort_values("Mean|SHAP|", ascending=False).reset_index(drop=True)
    shap_df_rf.insert(0, "Target", "AllOutputs")
    sheets["RF_SHAP_Global"] = shap_df_rf
except Exception as e:
    print(f"[WARN] RF SHAP computation failed: {e}")
    sheets["RF_SHAP_Global"] = pd.DataFrame([{"Note": f"RF SHAP failed: {e}"}])


# --- XGB SHAP ---
xgb_shap_global = []
for i, est in enumerate(xgb.estimators_):
    try:
        expl_xgb = shap.Explainer(est, Xtr)
        sv_xgb = expl_xgb(Xte_sub)
        shap_df_xgb = pd.DataFrame({
            "Feature": Xte_sub.columns,
            "Mean|SHAP|": np.abs(sv_xgb.values).mean(axis=0)
        }).sort_values("Mean|SHAP|", ascending=False).reset_index(drop=True)
        shap_df_xgb.insert(0, "Target", y.columns[i])
        xgb_shap_global.append(shap_df_xgb)
    except Exception as e:
        print(f"[WARN] XGB SHAP for target {y.columns[i]} failed: {e}")
        xgb_shap_global.append(pd.DataFrame([{"Target": y.columns[i], "Note": str(e)}]))

sheets["XGB_SHAP_Global"] = pd.concat(xgb_shap_global, ignore_index=True)

# ==============================================================
# Summary Sheet
# ==============================================================
sheets["Summary"] = pd.concat([
    rf_metrics[rf_metrics["Target"]=="__macro__"],
    xgb_metrics[xgb_metrics["Target"]=="__macro__"],
    enet_metrics[enet_metrics["Target"]=="__macro__"]
], ignore_index=True)[["Model","R2","RMSE","MAE","MAPE","AIC_approx","Coverage_<5%","Coverage_<10%"]]

# ==============================================================
# Export Results
# ==============================================================
with pd.ExcelWriter("Results_ML_XGB_ElasticNet_SHAP.xlsx", engine="openpyxl") as w:
    for name, df_sheet in sheets.items():
        df_sheet.to_excel(w, sheet_name=name[:31], index=False)

print("✅ Saved Results_ML_XGB_ElasticNet_SHAP.xlsx successfully.")