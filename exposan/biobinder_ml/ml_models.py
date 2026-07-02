# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# ===============================================================
# DNN Upgrade: Physics-Aware Softmax Neural Network
# ===============================================================
import os, warnings
import numpy as np, pandas as pd
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from datetime import datetime
from sklearn.multioutput import MultiOutputRegressor
warnings.filterwarnings("ignore")
timestamp = datetime.now().strftime("%Y%m%d_%H%M")
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
# 3. Define improved DNN (Softmax + Huber)
# ---------------------------------------------------------------
def build_dnn(input_dim):
    inputs = keras.Input(shape=(input_dim,))
    x = layers.BatchNormalization()(inputs)
    x = layers.Dense(128, activation="relu", kernel_regularizer=keras.regularizers.l2(1e-2))(x)
    x = layers.Dense(128, activation="relu", kernel_regularizer=keras.regularizers.l2(1e-2))(x)
    x = layers.Dropout(0.2)(x)
    x = layers.Dense(64, activation="relu", kernel_regularizer=keras.regularizers.l2(1e-2))(x)
    logits = layers.Dense(4)(x)
    outputs = layers.Lambda(lambda t: 100.0 * tf.nn.softmax(t))(logits)
    model = keras.Model(inputs, outputs)
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=5e-4),
        loss=keras.losses.Huber(delta=2.0),
        metrics=[keras.metrics.MeanAbsoluteError(name="MAE")]
    )
    return model

# ---------------------------------------------------------------
# 4. 10-Fold CV on TRAIN
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

    scaler = StandardScaler()
    Xtr_s = scaler.fit_transform(Xtr)
    Xva_s = scaler.transform(Xva)

    # ---------------- DNN (Keras Softmax) ----------------
    dnn = build_dnn(Xtr_s.shape[1])
    cb = [
        keras.callbacks.EarlyStopping(patience=12, restore_best_weights=True),
        keras.callbacks.ReduceLROnPlateau(patience=6, factor=0.5, min_lr=1e-5)
    ]
    dnn.fit(Xtr_s, ytr, validation_split=0.15, epochs=400, batch_size=64, verbose=0, callbacks=cb)
    yhat_dnn = dnn.predict(Xva_s)
    dnn_all.append(fold_metrics(yva, yhat_dnn, "DNN(Softmax)", f, targets))
    _, c5 = coverage(yva, yhat_dnn, 5.0); _, c10 = coverage(yva, yhat_dnn, 10.0)
    cov_rows.append({"Model":"DNN(Softmax)","Fold":f,"Coverage_<5%":c5,"Coverage_<10%":c10})

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
# 5. Aggregate CV results
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
# 6. Refit on full TRAIN + evaluate on TEST
# ---------------------------------------------------------------
imp = SimpleImputer(strategy="median")
Xtr_all = pd.DataFrame(imp.fit_transform(X_train_full), columns=X_train_full.columns)
Xte_all = pd.DataFrame(imp.transform(X_test), columns=X_train_full.columns)

# scale for DNN
scaler = StandardScaler()
Xtr_s = scaler.fit_transform(Xtr_all)
Xte_s = scaler.transform(Xte_all)

models_final = {}
metrics_final = []

# ---- DNN (Softmax) ----
dnn_final = build_dnn(Xtr_s.shape[1])
cb = [
    keras.callbacks.EarlyStopping(patience=12, restore_best_weights=True),
    keras.callbacks.ReduceLROnPlateau(patience=6, factor=0.5, min_lr=1e-5)
]
dnn_final.fit(Xtr_s, y_train_full, validation_split=0.15, epochs=400, batch_size=64, verbose=0, callbacks=cb)
y_pred_dnn = dnn_final.predict(Xte_s)
metrics_final.append(fold_metrics(y_test, y_pred_dnn, "DNN(Softmax)_Test", "HoldOut", targets))

# ---- Random Forest ----
# rf_final = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
# rf_final.fit(Xtr_all, y_train_full)
# y_pred_rf = rf_final.predict(Xte_all)
# metrics_final.append(fold_metrics(y_test, y_pred_rf, "RandomForest_Test", "HoldOut", targets))
# ---- Random Forest + Grid Search ----
print("\n🔎 Tuning Random Forest with GridSearchCV...")

rf_base = RandomForestRegressor(
    random_state=42,
    n_jobs=-1
)

rf_model = MultiOutputRegressor(rf_base, n_jobs=-1)

rf_param_grid = {
    "estimator__n_estimators": [300, 500, 800],
    "estimator__max_depth": [None, 10, 20, 30],
    "estimator__min_samples_split": [2, 5, 10],
    "estimator__min_samples_leaf": [1, 2, 4],
    "estimator__max_features": ["sqrt", 0.5, 0.8]
}

rf_grid = GridSearchCV(
    estimator=rf_model,
    param_grid=rf_param_grid,
    scoring="neg_root_mean_squared_error",
    cv=5,
    n_jobs=-1,
    verbose=2,
    refit=True
)

rf_grid.fit(Xtr_all, y_train_full)

rf_final = rf_grid.best_estimator_
y_pred_rf = rf_final.predict(Xte_all)

print("✅ Best RF parameters found:")
print(rf_grid.best_params_)
print(f"✅ Best CV score (neg RMSE): {rf_grid.best_score_:.4f}")

metrics_final.append(
    fold_metrics(y_test, y_pred_rf, "RandomForest_Tuned_Test", "HoldOut", targets)
)
rf_best_params_df = pd.DataFrame([rf_grid.best_params_])
rf_best_score_df = pd.DataFrame([{
    "Best_CV_Score_neg_RMSE": rf_grid.best_score_
}])
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
# 7. Save all
# ---------------------------------------------------------------


with pd.ExcelWriter(os.path.join(output_dir, f"Results_10Fold_DNN_Softmax_RF_GBR_{datetime.now():%Y%m%d_%H%M}.xlsx")) as w:
    dnn_df.to_excel(w, sheet_name="DNN_10Fold", index=False)
    rf_df.to_excel(w,  sheet_name="RF_10Fold",  index=False)
    gbr_df.to_excel(w, sheet_name="GBR_10Fold", index=False)
    cov_df.to_excel(w, sheet_name="Fold_Coverage", index=False)
    summary.to_excel(w, sheet_name="CV_Summary", index=False)
    test_metrics.to_excel(w, sheet_name="TestSet_Results", index=False)
    rf_best_params_df.to_excel(w, sheet_name="RF_Best_Params", index=False)
    rf_best_score_df.to_excel(w, sheet_name="RF_Best_CV_Score", index=False)


# from joblib import dump, load
# dump(rf_final, "results/rf_yield_model.joblib")


# ---------------------------------------------------------------
# 7. Parity and Residual Plots for All Models × Yields
# ---------------------------------------------------------------
# import matplotlib.pyplot as plt
# import seaborn as sns

# plt.style.use("seaborn-v0_8-whitegrid")
# fig_dir = os.path.join(output_dir, f"plots_{timestamp}")
# os.makedirs(fig_dir, exist_ok=True)

# pred_dict = {
#     "DNN": y_pred_dnn,
#     "RF":  y_pred_rf,
#     "GBR": y_pred_gbr
# }

# # ---- LOOP THROUGH EACH YIELD ----
# for i, col in enumerate(yield_cols):
#     plt.figure(figsize=(14,6))

#     # === Parity plot ===
#     plt.subplot(1,2,1)
#     for model, preds in pred_dict.items():
#         sns.scatterplot(
#             x=y_test[:,i], y=preds[:,i],
#             s=60, alpha=0.7, edgecolor="none", label=model
#         )
#     lims = [min(y_test[:,i].min(), preds[:,i].min())-2,
#             max(y_test[:,i].max(), preds[:,i].max())+2]
#     plt.plot(lims, lims, "k--", lw=2)
#     plt.xlabel("Observed", fontsize=14, fontweight="bold")
#     plt.ylabel("Predicted", fontsize=14, fontweight="bold")
#     plt.title(f"Parity Plot – {col}", fontsize=16, fontweight="bold")
#     plt.legend(fontsize=12)

#     # === Residual plot ===
#     plt.subplot(1,2,2)
#     for model, preds in pred_dict.items():
#         residuals = preds[:,i] - y_test[:,i]
#         sns.scatterplot(
#             x=y_test[:,i], y=residuals,
#             s=60, alpha=0.7, edgecolor="none", label=model
#         )
#     plt.axhline(0, color="k", linestyle="--", lw=2)
#     plt.xlabel("Observed", fontsize=14, fontweight="bold")
#     plt.ylabel("Residual (Pred - Obs)", fontsize=14, fontweight="bold")
#     plt.title(f"Residual Plot – {col}", fontsize=16, fontweight="bold")
#     plt.legend(fontsize=12)

#     plt.tight_layout()
#     fname = f"ParityResidual_{col.replace('%','pct').replace(' ','_')}.png"
#     plt.savefig(os.path.join(fig_dir, fname), dpi=300, bbox_inches="tight")
#     plt.close()
#     print(f"✅ Saved {fname}")


# # ---------------------------------------------------------------
# # 6. SHAP analysis for model interpretability
# # ---------------------------------------------------------------
# # ---------------------------------------------------------------
# # 6. SHAP Analysis (DNN, RF, GBR)
# # ---------------------------------------------------------------
# import shap, matplotlib.pyplot as plt
# from sklearn.ensemble import GradientBoostingRegressor
# import numpy as np, os

# print("⚙️ Computing SHAP values (may take several minutes)...")

# # Background and test subsets for SHAP speed
# X_bg = Xtr_all.sample(n=min(200, len(Xtr_all)), random_state=42)
# X_sample = Xte_all.sample(n=min(300, len(Xte_all)), random_state=42)
# feature_names = list(Xtr_all.columns)

# # ===============================================================
# # Generic SHAP function for any model (multi/single-output)
# # ===============================================================
# def compute_shap_and_plot(model, X_bg, X_sample, feature_names, model_name):
#     print(f"\n🔍 SHAP analysis for {model_name}...")
#     explainer = shap.Explainer(model, X_bg)
#     shap_values = explainer(X_sample)

#     # --- Case A: list (1 per output)
#     if isinstance(shap_values, list):
#         for i, sv in enumerate(shap_values):
#             yield_name = yield_cols[i] if i < len(yield_cols) else f"Output_{i}"
#             print(f"  • Plotting {yield_name}...")
#             shap.summary_plot(sv, X_sample, feature_names=feature_names, show=False)
#             plt.title(f"{model_name} — {yield_name}")
#             plt.tight_layout()
#             fname = os.path.join(output_dir, f"{model_name}_SHAP_{yield_name.replace('%','pct').replace(' ','_')}.png")
#             plt.savefig(fname, dpi=300, bbox_inches="tight")
#             plt.close()
#             print(f"    ✅ Saved {fname}")

#     # --- Case B: 3D array (samples × features × targets)
#     elif hasattr(shap_values, "values") and shap_values.values.ndim == 3:
#         n_targets = shap_values.values.shape[2]
#         for i in range(n_targets):
#             yield_name = yield_cols[i] if i < len(yield_cols) else f"Output_{i}"
#             print(f"  • Plotting {yield_name}...")
#             vals_i = shap_values.values[:, :, i]
#             shap.summary_plot(vals_i, X_sample, feature_names=feature_names, show=False)
#             plt.title(f"{model_name} — {yield_name}")
#             plt.tight_layout()
#             fname = os.path.join(output_dir, f"{model_name}_SHAP_{yield_name.replace('%','pct').replace(' ','_')}.png")
#             plt.savefig(fname, dpi=300, bbox_inches="tight")
#             plt.close()
#             print(f"    ✅ Saved {fname}")

#     # --- Case C: single-output model
#     else:
#         shap.summary_plot(shap_values, X_sample, feature_names=feature_names, show=False)
#         plt.title(f"{model_name} — All Outputs")
#         plt.tight_layout()
#         fname = os.path.join(output_dir, f"{model_name}_SHAP_AllOutputs.png")
#         plt.savefig(fname, dpi=300, bbox_inches="tight")
#         plt.close()
#         print(f"    ✅ Saved {fname}")

#     return shap_values


# # ===============================================================
# # Run SHAP for each model
# # ===============================================================
# # DNN (Keras or sklearn-MLP — works both ways)
# shap_dnn = compute_shap_and_plot(dnn_final, X_bg, X_sample, feature_names, "DNN(MLP)")

# # Random Forest (multi-output)
# shap_rf = compute_shap_and_plot(rf_final, X_bg, X_sample, feature_names, "RandomForest")

# # Gradient Boosting — train & compute per yield
# print("\n🔹 GradientBoosting (per-yield SHAP)")
# shap_gbr_list = []
# for i, target in enumerate(yield_cols):
#     print(f"  • Training + SHAP for {target}...")
#     gbr = GradientBoostingRegressor(random_state=42,
#                                     n_estimators=600, learning_rate=0.05, max_depth=3)
#     gbr.fit(Xtr_all, y_train_full[:, i])
#     explainer = shap.Explainer(gbr, X_bg)
#     shap_values = explainer(X_sample)
#     shap_gbr_list.append(shap_values)

#     shap.summary_plot(shap_values, X_sample, feature_names=feature_names, show=False)
#     plt.title(f"GradientBoosting — {target}")
#     plt.tight_layout()
#     fname = os.path.join(output_dir, f"GradientBoosting_SHAP_{target.replace('%','pct').replace(' ','_')}.png")
#     plt.savefig(fname, dpi=300, bbox_inches="tight")
#     plt.close()
#     print(f"    ✅ Saved {fname}")

# shap_gbr = shap_gbr_list

# # ===============================================================
# # 7. Aggregate mean(|SHAP|) feature importances and export to Excel
# # ===============================================================
# def summarize_shap_to_table(shap_values, X_sample, model_name):
#     feature_names = X_sample.columns
#     records = []
#     if isinstance(shap_values, list):
#         for i, sv in enumerate(shap_values):
#             vals = np.abs(sv.values) if hasattr(sv, "values") else np.abs(sv)
#             mean_abs = np.mean(vals, axis=0)
#             for feat, imp in zip(feature_names, mean_abs):
#                 records.append({"Model": model_name, "Yield": yield_cols[i],
#                                 "Feature": feat, "Mean|SHAP|": imp})
#     elif hasattr(shap_values, "values") and shap_values.values.ndim == 3:
#         n_targets = shap_values.values.shape[2]
#         for i in range(n_targets):
#             mean_abs = np.mean(np.abs(shap_values.values[:, :, i]), axis=0)
#             for feat, imp in zip(feature_names, mean_abs):
#                 records.append({"Model": model_name, "Yield": yield_cols[i],
#                                 "Feature": feat, "Mean|SHAP|": imp})
#     else:
#         vals = shap_values.values if hasattr(shap_values, "values") else shap_values
#         mean_abs = np.mean(np.abs(vals), axis=0)
#         for feat, imp in zip(feature_names, mean_abs):
#             records.append({"Model": model_name, "Yield": "AllOutputs",
#                             "Feature": feat, "Mean|SHAP|": imp})
#     return pd.DataFrame(records)

# print("\n📊 Aggregating SHAP feature importances...")
# df_dnn = summarize_shap_to_table(shap_dnn, X_sample, "DNN(MLP)")
# df_rf  = summarize_shap_to_table(shap_rf,  X_sample, "RandomForest")
# df_gbr = summarize_shap_to_table(shap_gbr, X_sample, "GradientBoosting")

# df_all = pd.concat([df_dnn, df_rf, df_gbr], ignore_index=True)
# df_ranked = (df_all.groupby(["Model", "Yield", "Feature"], as_index=False)
#                     .mean(numeric_only=True)
#                     .sort_values(["Model", "Yield", "Mean|SHAP|"], ascending=[True, True, False]))

# shap_excel = os.path.join(output_dir, "Feature_Importance_SHAP.xlsx")
# with pd.ExcelWriter(shap_excel, engine="openpyxl") as w:
#     df_all.to_excel(w, sheet_name="Raw_SHAP", index=False)
#     df_ranked.to_excel(w, sheet_name="Ranked_SHAP", index=False)

# print(f"✅ SHAP feature importance tables saved → {shap_excel}")


# ---------------------------------------------------------------
# 6. SHAP Feature Importance & Summary Plots (per yield type)
# ---------------------------------------------------------------
# import shap, matplotlib.pyplot as plt

# print("⚙️ Computing SHAP values (may take several minutes)...")

# # Background (train subset) and test subset for speed
# X_bg = Xtr_all.sample(n=min(200, len(Xtr_all)), random_state=42)
# X_sample = Xte_all.sample(n=min(300, len(Xte_all)), random_state=42)
# feature_names = list(Xtr_all.columns)

# # --- Helper ---
# def compute_shap_and_plot(model, X_bg, X_sample, feature_names, model_name):
#     print(f"\n🔍 SHAP analysis for {model_name}...")
#     explainer = shap.Explainer(model, X_bg)
#     shap_values = explainer(X_sample)

#     # 1️⃣ Case A — SHAP returns list (1 per output)
#     if isinstance(shap_values, list):
#         for i, sv in enumerate(shap_values):
#             yield_name = yield_cols[i] if i < len(yield_cols) else f"Output_{i}"
#             print(f"  • Plotting {yield_name}...")
#             shap.summary_plot(sv, X_sample, feature_names=feature_names, show=False)
#             plt.title(f"{model_name} — {yield_name}")
#             plt.tight_layout()
#             fname = f"{model_name}_SHAP_{yield_name.replace('%','pct').replace(' ','_')}.png"
#             plt.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches="tight")
#             plt.close()
#             print(f"    ✅ Saved {fname}")

#     # 2️⃣ Case B — SHAP returns 3D array (samples × features × targets)
#     elif hasattr(shap_values, "values") and shap_values.values.ndim == 3:
#         n_targets = shap_values.values.shape[2]
#         for i in range(n_targets):
#             yield_name = yield_cols[i] if i < len(yield_cols) else f"Output_{i}"
#             print(f"  • Plotting {yield_name}...")
#             vals_i = shap_values.values[:, :, i]
#             shap.summary_plot(vals_i, X_sample, feature_names=feature_names, show=False)
#             plt.title(f"{model_name} — {yield_name}")
#             plt.tight_layout()
#             fname = f"{model_name}_SHAP_{yield_name.replace('%','pct').replace(' ','_')}.png"
#             plt.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches="tight")
#             plt.close()
#             print(f"    ✅ Saved {fname}")

#     # 3️⃣ Case C — Single-output model
#     else:
#         shap.summary_plot(shap_values, X_sample, feature_names=feature_names, show=False)
#         plt.title(f"{model_name} — All Outputs")
#         plt.tight_layout()
#         fname = f"{model_name}_SHAP_AllOutputs.png"
#         plt.savefig(os.path.join(output_dir, fname), dpi=300, bbox_inches="tight")
#         plt.close()
#         print(f"    ✅ Saved {fname}")

#     return shap_values

# # --- Run for each model ---
# shap_dnn = compute_shap_and_plot(dnn_final, X_bg, X_sample, feature_names, "DNN(MLP)")
# shap_rf  = compute_shap_and_plot(rf_final,  X_bg, X_sample, feature_names, "RandomForest")
# shap_gbr = compute_shap_and_plot(gbr_final, X_bg, X_sample, feature_names, "GradientBoosting")

# print("\n🎨 All SHAP summary plots (per yield) saved successfully!")

# ---------------------------------------------------------------
# 6. SHAP Feature Importance & Summary Plots
# ---------------------------------------------------------------
# import shap, matplotlib.pyplot as plt

# print("⚙️ Computing SHAP values (may take several minutes)...")

# # For speed & interpretability, use a random subset
# X_bg = Xtr_all.sample(n=min(200, len(Xtr_all)), random_state=42)
# X_sample = Xte_all.sample(n=min(300, len(Xte_all)), random_state=42)
# feature_names = list(Xtr_all.columns)

# # --- Helper to compute & save plots ---
# def compute_shap_and_plot(model, X_bg, X_sample, feature_names, model_name):
#     print(f"🔍 SHAP analysis for {model_name}...")
#     explainer = shap.Explainer(model, X_bg)
#     shap_values = explainer(X_sample)
#     shap.summary_plot(shap_values, X_sample, feature_names=feature_names, show=False)
#     plt.title(f"{model_name} — SHAP Summary")
#     plt.tight_layout()
#     plot_path = os.path.join(output_dir, f"{model_name}_SHAP_Summary.png")
#     plt.savefig(plot_path, dpi=300, bbox_inches="tight")
#     plt.close()
#     print(f"✅ Saved {plot_path}")
#     return shap_values

# # --- Run for each model ---
# shap_dnn = compute_shap_and_plot(dnn_final, X_bg, X_sample, feature_names, "DNN(MLP)")
# shap_rf  = compute_shap_and_plot(rf_final,  X_bg, X_sample, feature_names, "RandomForest")
# shap_gbr = compute_shap_and_plot(gbr_final, X_bg, X_sample, feature_names, "GradientBoosting")

# print("🎨 All SHAP summary plots saved.")

# # ===============================================================
# # 7. SHAP feature importance (robust version)
# # ===============================================================
# import shap
# from tqdm import tqdm
# # Ensure test matrix has correct column names
# Xte_all = pd.DataFrame(Xte_all, columns=X_train_full.columns)

# print("\n⚙️ Computing SHAP values (this may take a few minutes)...")

# # --- Helper function ---
# def summarize_shap(shap_values, features, model_name):
#     """
#     Handle both single-output and multi-output SHAP arrays.
#     shap_values: list or np.ndarray
#     features: feature names
#     """
#     # Case 1: list of arrays (multi-target)
#     if isinstance(shap_values, list):
#         # e.g. [array(n,features) for each target]
#         n_targets = len(shap_values)
#         abs_mean = np.stack([np.abs(sv).mean(axis=0) for sv in shap_values], axis=1)
#         cols = yield_cols[:n_targets]
#     else:
#         # Case 2: single array (n, features)
#         abs_mean = np.abs(shap_values).mean(axis=0).reshape(-1, 1)
#         cols = [yield_cols[0]]  # dummy for single target

#     df = pd.DataFrame(abs_mean, index=list(features)[:abs_mean.shape[0]], columns=cols)
#     df["MeanAbsSHAP"] = df.mean(axis=1)
#     df = df.sort_values("MeanAbsSHAP", ascending=False).reset_index()
#     df.rename(columns={"index":"Feature"}, inplace=True)
#     df.insert(0, "Model", model_name)
#     return df

# # --- Tree models ---
# explainer_rf  = shap.TreeExplainer(rf_final)
# shap_values_rf  = explainer_rf.shap_values(Xte_all)
# rf_shap_df  = summarize_shap(shap_values_rf,  X_train_full.columns, "RandomForest")

# # For GBR: compute per-target SHAPs (list)
# shap_values_gbr = []
# for i in range(y_train_full.shape[1]):
#     gbr_tmp = GradientBoostingRegressor(random_state=42,
#                                         n_estimators=600, learning_rate=0.05, max_depth=3)
#     gbr_tmp.fit(Xtr_all, y_train_full[:, i])
#     expl = shap.TreeExplainer(gbr_tmp)
#     shap_values_gbr.append(expl.shap_values(Xte_all))
# gbr_shap_df = summarize_shap(shap_values_gbr, X_train_full.columns, "GradientBoosting")

# # --- DNN (Softmax): use KernelExplainer ---
# print("→ DNN SHAP (approx, using sample subset)...")
# X_background = Xtr_s[np.random.choice(len(Xtr_s), min(100, len(Xtr_s)), replace=False)]
# explainer_dnn = shap.KernelExplainer(dnn_final.predict, X_background)
# X_sample = Xte_s[:200]
# shap_values_dnn = explainer_dnn.shap_values(X_sample, nsamples=200)
# dnn_shap_df = summarize_shap(np.array(shap_values_dnn), X_train_full.columns, "DNN(Softmax)")

# # --- Combine & Save ---
# shap_summary = pd.concat([rf_shap_df, gbr_shap_df, dnn_shap_df], ignore_index=True)
# shap_summary_path = os.path.join(output_dir, "SHAP_FeatureImportance.xlsx")

# with pd.ExcelWriter(shap_summary_path) as w:
#     rf_shap_df.to_excel(w, sheet_name="RF_SHAP", index=False)
#     gbr_shap_df.to_excel(w, sheet_name="GBR_SHAP", index=False)
#     dnn_shap_df.to_excel(w, sheet_name="DNN_SHAP", index=False)
#     shap_summary.to_excel(w, sheet_name="Summary", index=False)

# print(f"✅ SHAP analysis complete — saved to '{shap_summary_path}'")

# # --- Optional plot for GBR ---
# try:
#     shap.summary_plot(shap_values_gbr[0], Xte_all, feature_names=X_train_full.columns,
#                       show=False, max_display=15)
#     import matplotlib.pyplot as plt
#     plt.title("GBR SHAP Summary (Biocrude example)")
#     plt.tight_layout()
#     plt.savefig(os.path.join(output_dir, "GBR_SHAP_Summary.png"), dpi=300)
#     plt.close()
# except Exception as e:
#     print("⚠️ Plot skipped:", e)

