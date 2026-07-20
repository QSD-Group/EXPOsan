# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 13:32:23 2025

@author: aliah
"""

import os, warnings
import numpy as np, pandas as pd
from sklearn.model_selection import KFold
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

file_path = os.path.join(output_dir, "HTL_all_yield_normalized_split_stratified_10pct.xlsx")

train_df = pd.read_excel(file_path, sheet_name="Train")
test_df  = pd.read_excel(file_path, sheet_name="Test")

yield_cols = ["Biocrude wt%", "Aqueous wt%", "Gas wt%", "Solids wt%"]
X_train_full = train_df.drop(columns=yield_cols).select_dtypes(include=[np.number])
y_train_full = train_df[yield_cols].values
X_test = test_df.drop(columns=yield_cols).select_dtypes(include=[np.number])
y_test = test_df[yield_cols].values

# ---------------------------------------------------------------
# 2. Helpers
# ---------------------------------------------------------------
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

# ---------------------------------------------------------------
# 3. 10-Fold CV on TRAIN
# ---------------------------------------------------------------
kf = KFold(n_splits=10, shuffle=True, random_state=42)
targets = yield_cols
dnn_all, rf_all, gbr_all, cov_rows = [], [], [], []

for f,(tr,va) in enumerate(kf.split(X_train_full), 1):
    print(f"Fold {f}...")

    Xtr, Xva = X_train_full.iloc[tr], X_train_full.iloc[va]
    ytr, yva = y_train_full[tr], y_train_full[va]

    imp = SimpleImputer(strategy="median")
    Xtr = pd.DataFrame(imp.fit_transform(Xtr), columns=X_train_full.columns)
    Xva = pd.DataFrame(imp.transform(Xva), columns=X_train_full.columns)

    # ---------------- DNN ----------------
    dnn = Pipeline([
        ("scaler", StandardScaler()),
        ("mlp", MLPRegressor(hidden_layer_sizes=(64,128,128,64,64,64,64),
                             activation="relu", alpha=0.01,
                             learning_rate_init=5e-4, max_iter=600,
                             random_state=42, early_stopping=True,
                             n_iter_no_change=20, verbose=False))
    ])
    dnn.fit(Xtr, ytr)
    yhat_dnn = dnn.predict(Xva)
    dnn_all.append(fold_metrics(yva, yhat_dnn, "DNN(MLP)", f, targets))
    _, c5 = coverage(yva, yhat_dnn, 5.0); _, c10 = coverage(yva, yhat_dnn, 10.0)
    cov_rows.append({"Model":"DNN(MLP)","Fold":f,"Coverage_<5%":c5,"Coverage_<10%":c10})

    # ---------------- RF ----------------
    rf = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
    rf.fit(Xtr, ytr)
    yhat_rf = rf.predict(Xva)
    rf_all.append(fold_metrics(yva, yhat_rf, "RandomForest", f, targets))
    _, c5 = coverage(yva, yhat_rf, 5.0); _, c10 = coverage(yva, yhat_rf, 10.0)
    cov_rows.append({"Model":"RandomForest","Fold":f,"Coverage_<5%":c5,"Coverage_<10%":c10})

    # ---------------- GBR ----------------
    yhat_gbr = np.zeros_like(yva)
    for i in range(y_train_full.shape[1]):
        gbr = GradientBoostingRegressor(random_state=42,
                                        n_estimators=600, learning_rate=0.05, max_depth=3)
        gbr.fit(Xtr, ytr[:,i])
        yhat_gbr[:,i] = gbr.predict(Xva)
    gbr_all.append(fold_metrics(yva, yhat_gbr, "GradientBoosting", f, targets))
    _, c5 = coverage(yva, yhat_gbr, 5.0); _, c10 = coverage(yva, yhat_gbr, 10.0)
    cov_rows.append({"Model":"GradientBoosting","Fold":f,"Coverage_<5%":c5,"Coverage_<10%":c10})

# ---------------------------------------------------------------
# 4. Aggregate CV results
# ---------------------------------------------------------------
dnn_df = pd.concat(dnn_all, ignore_index=True)
rf_df  = pd.concat(rf_all,  ignore_index=True)
gbr_df = pd.concat(gbr_all, ignore_index=True)
cov_df = pd.DataFrame(cov_rows)

summary = pd.concat([
    dnn_df[dnn_df["Target"]=="__macro__"],
    rf_df[rf_df["Target"]=="__macro__"],
    gbr_df[gbr_df["Target"]=="__macro__"]
], ignore_index=True).groupby("Model").mean(numeric_only=True).reset_index()

summary = summary.merge(
    cov_df.groupby("Model")[["Coverage_<5%","Coverage_<10%"]].mean().reset_index(),
    on="Model", how="left"
)

# ---------------------------------------------------------------
# 5. Refit on full TRAIN + evaluate on TEST
# ---------------------------------------------------------------
imp = SimpleImputer(strategy="median")
Xtr_all = pd.DataFrame(imp.fit_transform(X_train_full), columns=X_train_full.columns)
Xte_all = pd.DataFrame(imp.transform(X_test), columns=X_train_full.columns)

models_final = {}
metrics_final = []

# ---- DNN ----
dnn_final = Pipeline([
    ("scaler", StandardScaler()),
    ("mlp", MLPRegressor(hidden_layer_sizes=(64,128,128,64,64,64,64),
                         activation="relu", alpha=0.01,
                         learning_rate_init=5e-4, max_iter=600,
                         random_state=42, early_stopping=True,
                         n_iter_no_change=20, verbose=False))
])
dnn_final.fit(Xtr_all, y_train_full)
y_pred_dnn = dnn_final.predict(Xte_all)
metrics_final.append(fold_metrics(y_test, y_pred_dnn, "DNN(MLP)_Test", "HoldOut", targets))

# ---- Random Forest ----
rf_final = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
rf_final.fit(Xtr_all, y_train_full)
y_pred_rf = rf_final.predict(Xte_all)
metrics_final.append(fold_metrics(y_test, y_pred_rf, "RandomForest_Test", "HoldOut", targets))

# ---- Gradient Boosting ----
y_pred_gbr = np.zeros_like(y_test)
for i in range(y_train_full.shape[1]):
    gbr_final = GradientBoostingRegressor(random_state=42,
                                          n_estimators=600, learning_rate=0.05, max_depth=3)
    gbr_final.fit(Xtr_all, y_train_full[:,i])
    y_pred_gbr[:,i] = gbr_final.predict(Xte_all)
metrics_final.append(fold_metrics(y_test, y_pred_gbr, "GradientBoosting_Test", "HoldOut", targets))

test_metrics = pd.concat(metrics_final, ignore_index=True)

# ---------------------------------------------------------------
# 6. Save all
# ---------------------------------------------------------------
with pd.ExcelWriter(os.path.join(output_dir, "Results_10Fold_DNN_RF_GBR_TrainTest.xlsx")) as w:
    dnn_df.to_excel(w, sheet_name="DNN_10Fold", index=False)
    rf_df.to_excel(w,  sheet_name="RF_10Fold", index=False)
    gbr_df.to_excel(w, sheet_name="GBR_10Fold", index=False)
    cov_df.to_excel(w, sheet_name="Fold_Coverage", index=False)
    summary.to_excel(w, sheet_name="CV_Summary", index=False)
    test_metrics.to_excel(w, sheet_name="TestSet_Results", index=False)

print("✅ All done — results saved to 'Results_10Fold_DNN_RF_GBR_TrainTest.xlsx'")