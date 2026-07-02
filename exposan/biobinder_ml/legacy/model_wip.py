# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 12:41:01 2025

@author: aliah
"""

import pandas as pd, numpy as np
from sklearn.model_selection import KFold
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

# --- load normalized dataset ---
df = pd.read_excel("HTL_all_yield_consistent_normalized.xlsx")
yield_cols = ["Biocrude wt%","Aqueous wt%","Gas wt%","Solids wt%"]
X = df.drop(columns=yield_cols).select_dtypes(include=[np.number]).copy()
y = df[yield_cols].copy().values

# impute missing numeric features
X = SimpleImputer(strategy="median").fit_transform(X)

# helpers
def coverage(y_true, y_pred, thr):
    abs_err = np.abs(y_true - y_pred)
    per_target = (abs_err < thr).mean(axis=0)*100.0
    return per_target, per_target.mean()

def fold_metrics(y_true, y_pred, model_name, fold_id, targets):
    rows=[]
    for i,t in enumerate(targets):
        yt, yp = y_true[:,i], y_pred[:,i]
        mse = mean_squared_error(yt, yp)
        rows.append(dict(Model=model_name, Fold=fold_id, Target=t,
                         R2=r2_score(yt, yp),
                         RMSE=np.sqrt(mse),
                         MAE=mean_absolute_error(yt, yp),
                         MAPE=np.nanmean(np.abs((yt-yp)/np.where(yt==0,np.nan,yt)))*100.0))
    dfm = pd.DataFrame(rows)
    macro = dfm.mean(numeric_only=True).to_dict()
    macro.update(Model=model_name, Fold=fold_id, Target="__macro__")
    return pd.concat([dfm, pd.DataFrame([macro])], ignore_index=True)

# 5-fold CV
kf = KFold(n_splits=5, shuffle=True, random_state=42)
targets = yield_cols
all_dnn, all_rf, all_gbr, cov_rows = [], [], [], []
for f,(tr,va) in enumerate(kf.split(X), 1):
    Xtr, Xva = X[tr], X[va]
    ytr, yva = y[tr], y[va]

    # DNN (multi-output)
    dnn = Pipeline([
        ("scaler", StandardScaler()),
        ("mlp", MLPRegressor(hidden_layer_sizes=(64,128,128,64,64,64,64),
                             activation="relu", alpha=0.01,
                             learning_rate_init=5e-4, max_iter=600,
                             random_state=42, early_stopping=True,
                             n_iter_no_change=20, verbose=False))
    ])
    dnn.fit(Xtr, ytr); yhat_dnn = dnn.predict(Xva)
    all_dnn.append(fold_metrics(yva, yhat_dnn, "DNN(MLP)", f, targets))
    _, cov5 = coverage(yva, yhat_dnn, 5.0)
    _, cov10 = coverage(yva, yhat_dnn, 10.0)
    cov_rows.append({"Model":"DNN(MLP)","Fold":f,"Coverage_<5%":cov5,"Coverage_<10%":cov10})

    # Random Forest
    rf = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
    rf.fit(Xtr, ytr); yhat_rf = rf.predict(Xva)
    all_rf.append(fold_metrics(yva, yhat_rf, "RandomForest", f, targets))
    _, cov5 = coverage(yva, yhat_rf, 5.0)
    _, cov10 = coverage(yva, yhat_rf, 10.0)
    cov_rows.append({"Model":"RandomForest","Fold":f,"Coverage_<5%":cov5,"Coverage_<10%":cov10})

    # Gradient Boosting (per target)
    gbr = GradientBoostingRegressor(random_state=42, n_estimators=600, learning_rate=0.05, max_depth=3)
    yhat_gbr = np.column_stack([gbr.fit(Xtr, ytr[:,i]).predict(Xva) for i in range(y.shape[1])])
    all_gbr.append(fold_metrics(yva, yhat_gbr, "GradientBoosting", f, targets))
    _, cov5 = coverage(yva, yhat_gbr, 5.0)
    _, cov10 = coverage(yva, yhat_gbr, 10.0)
    cov_rows.append({"Model":"GradientBoosting","Fold":f,"Coverage_<5%":cov5,"Coverage_<10%":cov10})

# aggregate + save
dnn_df = pd.concat(all_dnn, ignore_index=True)
rf_df  = pd.concat(all_rf, ignore_index=True)
gbr_df = pd.concat(all_gbr, ignore_index=True)
cov_df = pd.DataFrame(cov_rows)

summary = pd.concat([
    dnn_df[dnn_df["Target"]=="__macro__"].groupby("Model").mean(numeric_only=True).reset_index(),
    rf_df[rf_df["Target"]=="__macro__"].groupby("Model").mean(numeric_only=True).reset_index(),
    gbr_df[gbr_df["Target"]=="__macro__"].groupby("Model").mean(numeric_only=True).reset_index(),
], ignore_index=True)

summary = summary.merge(cov_df.groupby("Model")[["Coverage_<5%","Coverage_<10%"]].mean().reset_index(),
                        on="Model", how="left")

with pd.ExcelWriter("Results_KFold_DNN_vs_Trees.xlsx") as w:
    dnn_df.to_excel(w, sheet_name="DNN_KFold_Metrics", index=False)
    rf_df.to_excel(w, sheet_name="RF_KFold_Metrics", index=False)
    gbr_df.to_excel(w, sheet_name="GBR_KFold_Metrics", index=False)
    cov_df.to_excel(w, sheet_name="Fold_Coverage", index=False)
    summary.to_excel(w, sheet_name="Summary", index=False)
print("Saved Results_KFold_DNN_vs_Trees.xlsx")
