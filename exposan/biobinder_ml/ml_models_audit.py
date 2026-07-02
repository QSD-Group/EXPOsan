# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""
Extended HTL yield model benchmarking script

What this script does
---------------------
1. Loads predefined Train/Test sheets
2. Runs 10-fold CV on Train only
3. Trains/evaluates:
      - DNN (multi-output)
      - Random Forest (multi-output)
      - Gradient Boosting (4 separate single-output models)
4. Reports BOTH:
      A) Your style:
         - R2
         - RMSE
         - MAE
         - MAPE
      B) Other-user / HTL-style:
         - Median residual
         - Mean absolute residual
         - Median absolute residual
         - % < 5 wt%
         - % < 10 wt%
5. Refits on full Train and evaluates on hidden Test
6. Audits train/test overlap and potential leakage-like issues
7. Reports model complexity proxies
8. Saves everything to one Excel workbook

"""

import os
import warnings
import numpy as np
import pandas as pd

from sklearn.model_selection import KFold
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

warnings.filterwarnings("ignore")


# ===============================================================
# 1. Paths and settings
# ===============================================================
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

file_path = os.path.join(output_dir, "HTL_all_yield_normalized_split_stratified_10pct.xlsx")
out_xlsx = os.path.join(output_dir, "Results_10Fold_DNN_RF_GBR_Extended_tuned.xlsx")

yield_cols = ["Biocrude wt%", "Aqueous wt%", "Gas wt%", "Solids wt%"]

cv_random_state = 42
rf_random_state = 42
gbr_random_state = 42
dnn_random_state = 42

n_splits = 10


# ===============================================================
# 2. Load data
# ===============================================================
train_df = pd.read_excel(file_path, sheet_name="Train")
test_df = pd.read_excel(file_path, sheet_name="Test")

# Keep only numeric predictors, same as your original script
X_train_full = train_df.drop(columns=yield_cols).select_dtypes(include=[np.number]).copy()
y_train_full = train_df[yield_cols].to_numpy(dtype=float)

X_test = test_df.drop(columns=yield_cols).select_dtypes(include=[np.number]).copy()
y_test = test_df[yield_cols].to_numpy(dtype=float)

feature_cols = X_train_full.columns.tolist()
targets = yield_cols

print("=" * 70)
print("Data loaded")
print(f"Train shape: X={X_train_full.shape}, y={y_train_full.shape}")
print(f"Test  shape: X={X_test.shape}, y={y_test.shape}")
print("=" * 70)

# # ===============================================================
# # 2.5 Duplicate Feature Variability Analysis (Range + Histogram)
# # ===============================================================
# import numpy as np
# import matplotlib.pyplot as plt
# import os

# print("\n🔍 Running duplicate feature variability analysis...")

# # Separate features and targets
# X = train_df.drop(columns=yield_cols)
# y = train_df[yield_cols]

# # Group identical feature rows
# dup_groups = X.groupby(list(X.columns)).groups
# dup_indices = [idxs for idxs in dup_groups.values() if len(idxs) > 1]

# print(f"Total duplicate feature groups: {len(dup_indices)}")

# # ---------------------------------------------------------------
# # Compute variability stats
# # ---------------------------------------------------------------
# group_stats = []

# for i, idxs in enumerate(dup_indices):
#     vals = y.iloc[idxs]["Biocrude wt%"]
#     group_stats.append({
#         "Group": f"G{i}",
#         "std": np.std(vals),
#         "range": vals.max() - vals.min(),
#         "min": vals.min(),
#         "max": vals.max(),
#         "indices": idxs
#     })

# stats_df = pd.DataFrame(group_stats)

# # ---------------------------------------------------------------
# # Select top variable groups (cleaner visualization)
# # ---------------------------------------------------------------
# top_n = 12
# top_groups = stats_df.sort_values("range", ascending=False).head(top_n).reset_index(drop=True)

# # ---------------------------------------------------------------
# # 📊 RANGE PLOT (Best for slides)
# # ---------------------------------------------------------------
# plt.figure(figsize=(10, 5))

# for i, row in top_groups.iterrows():
#     plt.plot([i, i], [row["min"], row["max"]], linewidth=3)

# plt.ylabel("Biocrude yield (wt%)", fontsize=12)
# plt.title("Yield variability for identical inputs (Top groups)", fontsize=13)
# plt.xticks(range(len(top_groups)), top_groups["Group"], rotation=45)
# plt.grid(True, alpha=0.3)

# plt.tight_layout()

# range_path = os.path.join(output_dir, "duplicate_variability_range.png")
# plt.savefig(range_path, dpi=300)
# plt.show()

# print(f"✅ Range plot saved → {range_path}")

# # ---------------------------------------------------------------
# # 📊 HISTOGRAM (Best for paper)
# # ---------------------------------------------------------------
# ranges = stats_df["range"].values

# plt.figure(figsize=(6, 4))
# plt.hist(ranges, bins=12)
# plt.xlabel("Yield range (wt%)", fontsize=11)
# plt.ylabel("Frequency", fontsize=11)
# plt.title("Distribution of yield variability for identical inputs", fontsize=12)
# plt.grid(True, alpha=0.3)

# plt.tight_layout()

# hist_path = os.path.join(output_dir, "duplicate_variability_hist.png")
# plt.savefig(hist_path, dpi=300)
# plt.show()

# print(f"✅ Histogram saved → {hist_path}")
# ===============================================================
# 3. Helper functions
# ===============================================================
def safe_mape(y_true, y_pred):
    """
    MAPE in %, skipping zero denominators.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    denom = np.where(y_true == 0, np.nan, y_true)
    return np.nanmean(np.abs((y_true - y_pred) / denom)) * 100.0


def coverage(y_true, y_pred, thr):
    """
    Percentage of predictions with absolute error < threshold, per target and macro.
    """
    abs_err = np.abs(np.asarray(y_true) - np.asarray(y_pred))
    per_target = (abs_err < thr).mean(axis=0) * 100.0
    return per_target, float(np.mean(per_target))


def compute_target_metrics(yt, yp):
    """
    Metrics for one target.
    """
    yt = np.asarray(yt, dtype=float)
    yp = np.asarray(yp, dtype=float)

    resid = yt - yp
    abs_resid = np.abs(resid)

    mse = mean_squared_error(yt, yp)
    out = {
        "R2": r2_score(yt, yp),
        "RMSE": np.sqrt(mse),
        "MAE": mean_absolute_error(yt, yp),
        "MAPE": safe_mape(yt, yp),

        # HTL-style / other-user style
        "MedianResidual": np.nanmedian(resid),
        "MeanAbsResidual": np.nanmean(abs_resid),
        "MedianAbsResidual": np.nanmedian(abs_resid),
        "Pct_lt_5wt": (abs_resid < 5).mean() * 100.0,
        "Pct_lt_10wt": (abs_resid < 10).mean() * 100.0,
    }
    return out


def evaluate_multioutput(y_true, y_pred, model_name, split_name, fold_label, target_names):
    """
    Build a dataframe with per-target metrics + macro average row.
    """
    rows = []
    for i, t in enumerate(target_names):
        met = compute_target_metrics(y_true[:, i], y_pred[:, i])
        row = {
            "Model": model_name,
            "Split": split_name,
            "Fold": fold_label,
            "Target": t,
        }
        row.update(met)
        rows.append(row)

    df = pd.DataFrame(rows)

    macro = {
        "Model": model_name,
        "Split": split_name,
        "Fold": fold_label,
        "Target": "__macro__",
    }
    for col in ["R2", "RMSE", "MAE", "MAPE",
                "MedianResidual", "MeanAbsResidual", "MedianAbsResidual",
                "Pct_lt_5wt", "Pct_lt_10wt"]:
        macro[col] = df[col].mean()

    df = pd.concat([df, pd.DataFrame([macro])], ignore_index=True)
    return df


def summarize_macro(df):
    """
    Aggregate macro rows across folds or splits.
    """
    macro_df = df[df["Target"] == "__macro__"].copy()

    group_cols = ["Model", "Split"]
    metric_cols = ["R2", "RMSE", "MAE", "MAPE",
                   "MedianResidual", "MeanAbsResidual", "MedianAbsResidual",
                   "Pct_lt_5wt", "Pct_lt_10wt"]

    mean_df = macro_df.groupby(group_cols)[metric_cols].mean().reset_index()
    std_df = macro_df.groupby(group_cols)[metric_cols].std().reset_index()

    std_df = std_df.rename(columns={c: f"{c}_std" for c in metric_cols})
    out = mean_df.merge(std_df, on=group_cols, how="left")
    return out


def summarize_per_target(df):
    """
    Aggregate per-target rows across folds or splits.
    """
    df2 = df[df["Target"] != "__macro__"].copy()

    group_cols = ["Model", "Split", "Target"]
    metric_cols = ["R2", "RMSE", "MAE", "MAPE",
                   "MedianResidual", "MeanAbsResidual", "MedianAbsResidual",
                   "Pct_lt_5wt", "Pct_lt_10wt"]

    mean_df = df2.groupby(group_cols)[metric_cols].mean().reset_index()
    std_df = df2.groupby(group_cols)[metric_cols].std().reset_index()
    std_df = std_df.rename(columns={c: f"{c}_std" for c in metric_cols})

    out = mean_df.merge(std_df, on=group_cols, how="left")
    return out


def rf_leaf_complexity(rf_model):
    """
    Complexity proxy for RF: total terminal leaves across all trees.
    """
    return int(sum(tree.get_n_leaves() for tree in rf_model.estimators_))


def gbr_leaf_complexity(gbr_models):
    """
    Complexity proxy for GBR: total terminal leaves across all trees
    across all 4 separate GBR models.
    """
    total = 0
    for model in gbr_models:
        # GradientBoostingRegressor.estimators_ shape = (n_estimators, 1)
        for est in model.estimators_.ravel():
            total += est.get_n_leaves()
    return int(total)


def dnn_weight_complexity(dnn_pipeline):
    """
    Complexity proxy for DNN: total weights + biases.
    """
    mlp = dnn_pipeline.named_steps["mlp"]
    total = 0
    for w, b in zip(mlp.coefs_, mlp.intercepts_):
        total += w.size + b.size
    return int(total)


def leakage_audit(train_df, test_df, yield_cols):
    """
    Audit exact row overlaps and exact feature overlaps between Train and Test.
    """
    feature_cols_full = [c for c in train_df.columns if c not in yield_cols]

    # exact full-row overlap
    full_overlap = train_df.merge(test_df, how="inner", on=train_df.columns.tolist())

    # exact feature-only overlap
    train_feat = train_df[feature_cols_full].copy()
    test_feat = test_df[feature_cols_full].copy()
    feat_overlap = train_feat.merge(test_feat, how="inner", on=feature_cols_full)

    # duplicate counts within each set
    train_dup_full = train_df.duplicated().sum()
    test_dup_full = test_df.duplicated().sum()
    train_dup_feat = train_feat.duplicated().sum()
    test_dup_feat = test_feat.duplicated().sum()

    audit = pd.DataFrame([
        {"Check": "Train rows", "Value": len(train_df)},
        {"Check": "Test rows", "Value": len(test_df)},
        {"Check": "Train duplicate full rows", "Value": int(train_dup_full)},
        {"Check": "Test duplicate full rows", "Value": int(test_dup_full)},
        {"Check": "Train duplicate feature rows", "Value": int(train_dup_feat)},
        {"Check": "Test duplicate feature rows", "Value": int(test_dup_feat)},
        {"Check": "Exact duplicate full rows across Train/Test", "Value": int(len(full_overlap))},
        {"Check": "Exact duplicate feature rows across Train/Test", "Value": int(len(feat_overlap))},
    ])

    return audit, full_overlap, feat_overlap


def make_dnn():
    return Pipeline([
        ("scaler", StandardScaler()),
        ("mlp", MLPRegressor(
            hidden_layer_sizes=(64, 128, 128, 64, 64, 64, 64),
            activation="relu",
            alpha=0.01,
            learning_rate_init=5e-4,
            max_iter=600,
            random_state=dnn_random_state,
            early_stopping=True,
            n_iter_no_change=20,
            verbose=False
        ))
    ])


def make_rf():
    return RandomForestRegressor(
        n_estimators=500,
        random_state=rf_random_state,
        n_jobs=-1
    )


def make_gbr():
    return GradientBoostingRegressor(
        random_state=gbr_random_state,
        n_estimators=600,
        learning_rate=0.05,
        max_depth=3
    )
def make_rf_tuned():
    return RandomForestRegressor(
        n_estimators=500,
        max_depth=10,
        max_features=0.8,
        min_samples_split=2,
        min_samples_leaf=1,
        random_state=rf_random_state,
        n_jobs=-1
    )

# ===============================================================
# 4. 10-fold CV on TRAIN only
# ===============================================================
kf = KFold(n_splits=n_splits, shuffle=True, random_state=cv_random_state)

cv_all_rows = []
cv_complexity_rows = []

print("\nStarting 10-fold CV...\n")

for f, (tr_idx, va_idx) in enumerate(kf.split(X_train_full), start=1):
    print(f"Fold {f}/{n_splits}")

    Xtr = X_train_full.iloc[tr_idx].copy()
    Xva = X_train_full.iloc[va_idx].copy()
    ytr = y_train_full[tr_idx]
    yva = y_train_full[va_idx]

    # Leakage-safe imputation: fit only on fold-train
    imp = SimpleImputer(strategy="median")
    Xtr_imp = pd.DataFrame(imp.fit_transform(Xtr), columns=feature_cols)
    Xva_imp = pd.DataFrame(imp.transform(Xva), columns=feature_cols)

    # ---------------- DNN (multi-output) ----------------
    dnn = make_dnn()
    dnn.fit(Xtr_imp, ytr)
    yhat_dnn = dnn.predict(Xva_imp)

    dnn_eval = evaluate_multioutput(
        y_true=yva,
        y_pred=yhat_dnn,
        model_name="DNN(MLP)",
        split_name="CV",
        fold_label=f,
        target_names=targets
    )
    cv_all_rows.append(dnn_eval)

    cv_complexity_rows.append({
        "Model": "DNN(MLP)",
        "Split": "CV",
        "Fold": f,
        "Complexity_Definition": "weights_plus_biases",
        "Complexity_Value": dnn_weight_complexity(dnn)
    })

    # ---------------- RF (multi-output) ----------------
    rf = make_rf()
    rf.fit(Xtr_imp, ytr)
    yhat_rf = rf.predict(Xva_imp)

    rf_eval = evaluate_multioutput(
        y_true=yva,
        y_pred=yhat_rf,
        model_name="RandomForest",
        split_name="CV",
        fold_label=f,
        target_names=targets
    )
    cv_all_rows.append(rf_eval)

    cv_complexity_rows.append({
        "Model": "RandomForest",
        "Split": "CV",
        "Fold": f,
        "Complexity_Definition": "sum_terminal_leaves",
        "Complexity_Value": rf_leaf_complexity(rf)
    })

    # ---------------- GBR (4 separate single-output models) ----------------
    yhat_gbr = np.zeros_like(yva, dtype=float)
    gbr_models_fold = []

    for i in range(ytr.shape[1]):
        gbr = make_gbr()
        gbr.fit(Xtr_imp, ytr[:, i])
        yhat_gbr[:, i] = gbr.predict(Xva_imp)
        gbr_models_fold.append(gbr)

    gbr_eval = evaluate_multioutput(
        y_true=yva,
        y_pred=yhat_gbr,
        model_name="GradientBoosting",
        split_name="CV",
        fold_label=f,
        target_names=targets
    )
    cv_all_rows.append(gbr_eval)

    cv_complexity_rows.append({
        "Model": "GradientBoosting",
        "Split": "CV",
        "Fold": f,
        "Complexity_Definition": "sum_terminal_leaves_across_4_models",
        "Complexity_Value": gbr_leaf_complexity(gbr_models_fold)
    })

cv_results = pd.concat(cv_all_rows, ignore_index=True)
cv_complexity = pd.DataFrame(cv_complexity_rows)

cv_macro_summary = summarize_macro(cv_results)
cv_target_summary = summarize_per_target(cv_results)


# ===============================================================
# 5. Refit on full TRAIN and evaluate on hidden TEST
# ===============================================================
print("\nTraining final models on full Train and evaluating on Test...\n")

imp_final = SimpleImputer(strategy="median")
Xtr_all = pd.DataFrame(imp_final.fit_transform(X_train_full), columns=feature_cols)
Xte_all = pd.DataFrame(imp_final.transform(X_test), columns=feature_cols)

holdout_rows = []
holdout_complexity_rows = []

# ---------------- DNN final ----------------
dnn_final = make_dnn()
dnn_final.fit(Xtr_all, y_train_full)
y_pred_dnn = dnn_final.predict(Xte_all)

holdout_rows.append(
    evaluate_multioutput(
        y_true=y_test,
        y_pred=y_pred_dnn,
        model_name="DNN(MLP)",
        split_name="Holdout",
        fold_label="Holdout",
        target_names=targets
    )
)

holdout_complexity_rows.append({
    "Model": "DNN(MLP)",
    "Split": "Holdout",
    "Fold": "Holdout",
    "Complexity_Definition": "weights_plus_biases",
    "Complexity_Value": dnn_weight_complexity(dnn_final)
})

# ---------------- RF final (TUNED) ----------------
rf_final = make_rf_tuned()
rf_final.fit(Xtr_all, y_train_full)
y_pred_rf = rf_final.predict(Xte_all)

holdout_rows.append(
    evaluate_multioutput(
        y_true=y_test,
        y_pred=y_pred_rf,
        model_name="RandomForest",
        split_name="Holdout",
        fold_label="Holdout",
        target_names=targets
    )
)

holdout_complexity_rows.append({
    "Model": "RandomForest",
    "Split": "Holdout",
    "Fold": "Holdout",
    "Complexity_Definition": "sum_terminal_leaves",
    "Complexity_Value": rf_leaf_complexity(rf_final)
})

# ---------------- GBR final ----------------
y_pred_gbr = np.zeros_like(y_test, dtype=float)
gbr_models_final = []

for i in range(y_train_full.shape[1]):
    gbr_final = make_gbr()
    gbr_final.fit(Xtr_all, y_train_full[:, i])
    y_pred_gbr[:, i] = gbr_final.predict(Xte_all)
    gbr_models_final.append(gbr_final)

holdout_rows.append(
    evaluate_multioutput(
        y_true=y_test,
        y_pred=y_pred_gbr,
        model_name="GradientBoosting",
        split_name="Holdout",
        fold_label="Holdout",
        target_names=targets
    )
)

holdout_complexity_rows.append({
    "Model": "GradientBoosting",
    "Split": "Holdout",
    "Fold": "Holdout",
    "Complexity_Definition": "sum_terminal_leaves_across_4_models",
    "Complexity_Value": gbr_leaf_complexity(gbr_models_final)
})

holdout_results = pd.concat(holdout_rows, ignore_index=True)
holdout_complexity = pd.DataFrame(holdout_complexity_rows)

holdout_macro_summary = summarize_macro(holdout_results)
holdout_target_summary = summarize_per_target(holdout_results)


# ===============================================================
# 6. Combined summary tables
# ===============================================================
combined_results = pd.concat([cv_results, holdout_results], ignore_index=True)
combined_complexity = pd.concat([cv_complexity, holdout_complexity], ignore_index=True)

combined_macro_summary = pd.concat(
    [cv_macro_summary, holdout_macro_summary],
    ignore_index=True
)

combined_target_summary = pd.concat(
    [cv_target_summary, holdout_target_summary],
    ignore_index=True
)


# ===============================================================
# 7. Leakage audit
# ===============================================================
audit_df, full_overlap_df, feat_overlap_df = leakage_audit(train_df, test_df, yield_cols)


# ===============================================================
# 8. Optional: save final RF model for downstream use
# ===============================================================
# Uncomment if you want to overwrite or save a fresh RF joblib.
# from joblib import dump
# rf_joblib_path = os.path.join(output_dir, "rf_yield_model_retrained.joblib")
# dump(rf_final, rf_joblib_path)


# ===============================================================
# 9. Save all outputs
# ===============================================================
with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:

    # Raw CV / holdout results
    cv_results.to_excel(writer, sheet_name="CV_AllMetrics_Raw", index=False)
    holdout_results.to_excel(writer, sheet_name="Holdout_AllMetrics_Raw", index=False)
    combined_results.to_excel(writer, sheet_name="AllMetrics_Combined", index=False)

    # Summary tables
    cv_macro_summary.to_excel(writer, sheet_name="CV_Macro_Summary", index=False)
    cv_target_summary.to_excel(writer, sheet_name="CV_Target_Summary", index=False)

    holdout_macro_summary.to_excel(writer, sheet_name="Holdout_Macro_Summary", index=False)
    holdout_target_summary.to_excel(writer, sheet_name="Holdout_Target_Summary", index=False)

    combined_macro_summary.to_excel(writer, sheet_name="Combined_Macro_Summary", index=False)
    combined_target_summary.to_excel(writer, sheet_name="Combined_Target_Summary", index=False)

    # Complexity
    cv_complexity.to_excel(writer, sheet_name="CV_Model_Complexity", index=False)
    holdout_complexity.to_excel(writer, sheet_name="Holdout_Model_Complexity", index=False)
    combined_complexity.to_excel(writer, sheet_name="All_Model_Complexity", index=False)

    # Leakage / overlap audit
    audit_df.to_excel(writer, sheet_name="Leakage_Audit", index=False)
    full_overlap_df.to_excel(writer, sheet_name="Exact_Full_Overlap", index=False)
    feat_overlap_df.to_excel(writer, sheet_name="Exact_Feature_Overlap", index=False)

    # Metadata
    meta_df = pd.DataFrame([
        {"Key": "Input file", "Value": file_path},
        {"Key": "Output file", "Value": out_xlsx},
        {"Key": "Yield columns", "Value": ", ".join(yield_cols)},
        {"Key": "Numeric feature columns used", "Value": len(feature_cols)},
        {"Key": "Feature names", "Value": ", ".join(feature_cols)},
        {"Key": "CV splits", "Value": n_splits},
        {"Key": "CV random_state", "Value": cv_random_state},
        {"Key": "RF random_state", "Value": rf_random_state},
        {"Key": "GBR random_state", "Value": gbr_random_state},
        {"Key": "DNN random_state", "Value": dnn_random_state},
        {"Key": "RF output structure", "Value": "single multi-output model"},
        {"Key": "DNN output structure", "Value": "single multi-output model"},
        {"Key": "GBR output structure", "Value": "4 separate single-output models"},
    ])
    meta_df.to_excel(writer, sheet_name="Run_Metadata", index=False)

print("\n" + "=" * 70)
print(f"Done. Extended results saved to:\n{out_xlsx}")
print("=" * 70)


import pandas as pd
import matplotlib.pyplot as plt

# Load train data
train_df = pd.read_excel(file, sheet_name="Train")

# Separate features and targets
X = train_df.drop(columns=yield_cols)
y = train_df[yield_cols]

# Identify duplicate feature groups
dup_groups = X.groupby(list(X.columns)).groups

# Keep only groups with duplicates
dup_indices = [idxs for idxs in dup_groups.values() if len(idxs) > 1]

# Build plotting dataframe
rows = []
for i, idxs in enumerate(dup_indices):
    for idx in idxs:
        rows.append({
            "Group": f"G{i}",
            "Biocrude": y.iloc[idx]["Biocrude wt%"]
        })

plot_df = pd.DataFrame(rows)

# Plot
plt.figure()
plot_df.boxplot(column="Biocrude", by="Group")
plt.xticks(rotation=90)
plt.ylabel("Biocrude yield (wt%)")
plt.title("Variability in yield for identical input features")
plt.suptitle("")
plt.show()