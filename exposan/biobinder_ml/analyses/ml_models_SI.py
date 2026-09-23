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
      - DNN (Physics-Aware Softmax Keras Neural Network)
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
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

from sklearn.model_selection import KFold
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

warnings.filterwarnings("ignore")


# ===============================================================
# 1. Paths and settings
# ===============================================================
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

file_path = os.path.join(output_dir, "HTL_all_yield_normalized_split_stratified_10pct.xlsx")
out_xlsx = os.path.join(output_dir, "Results_10Fold_DNN_RF_GBR_Extended_tuned_SI_figure.xlsx")

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
        for est in model.estimators_.ravel():
            total += est.get_n_leaves()
    return int(total)


def dnn_weight_complexity(dnn_model):
    """
    Complexity proxy for Keras Softmax DNN: total trainable parameters.
    """
    return int(dnn_model.count_params())


def leakage_audit(train_df, test_df, yield_cols):
    """
    Audit exact row overlaps and exact feature overlaps between Train and Test.
    """
    feature_cols_full = [c for c in train_df.columns if c not in yield_cols]

    full_overlap = train_df.merge(test_df, how="inner", on=train_df.columns.tolist())

    train_feat = train_df[feature_cols_full].copy()
    test_feat = test_df[feature_cols_full].copy()
    feat_overlap = train_feat.merge(test_feat, how="inner", on=feature_cols_full)

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


def build_dnn(input_dim):
    """
    Physics-Aware Softmax + Huber Loss Keras Neural Network (from ml_models.py)
    """
    tf.random.set_seed(dnn_random_state)
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

    # ---------------- DNN (Physics-Aware Softmax) ----------------
    scaler_fold = StandardScaler()
    Xtr_s = scaler_fold.fit_transform(Xtr_imp)
    Xva_s = scaler_fold.transform(Xva_imp)

    dnn = build_dnn(Xtr_s.shape[1])
    cb = [
        keras.callbacks.EarlyStopping(patience=12, restore_best_weights=True),
        keras.callbacks.ReduceLROnPlateau(patience=6, factor=0.5, min_lr=1e-5)
    ]
    dnn.fit(Xtr_s, ytr, validation_split=0.15, epochs=400, batch_size=64, verbose=0, callbacks=cb)
    yhat_dnn = dnn.predict(Xva_s, verbose=0)

    dnn_eval = evaluate_multioutput(
        y_true=yva,
        y_pred=yhat_dnn,
        model_name="DNN(Softmax)",
        split_name="CV",
        fold_label=f,
        target_names=targets
    )
    cv_all_rows.append(dnn_eval)

    cv_complexity_rows.append({
        "Model": "DNN(Softmax)",
        "Split": "CV",
        "Fold": f,
        "Complexity_Definition": "total_trainable_parameters",
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

# ---------------- DNN final (Physics-Aware Softmax) ----------------
scaler_final = StandardScaler()
Xtr_all_s = scaler_final.fit_transform(Xtr_all)
Xte_all_s = scaler_final.transform(Xte_all)

dnn_final = build_dnn(Xtr_all_s.shape[1])
cb_final = [
    keras.callbacks.EarlyStopping(patience=12, restore_best_weights=True),
    keras.callbacks.ReduceLROnPlateau(patience=6, factor=0.5, min_lr=1e-5)
]
dnn_final.fit(Xtr_all_s, y_train_full, validation_split=0.15, epochs=400, batch_size=64, verbose=0, callbacks=cb_final)
y_pred_dnn = dnn_final.predict(Xte_all_s, verbose=0)

holdout_rows.append(
    evaluate_multioutput(
        y_true=y_test,
        y_pred=y_pred_dnn,
        model_name="DNN(Softmax)",
        split_name="Holdout",
        fold_label="Holdout",
        target_names=targets
    )
)

holdout_complexity_rows.append({
    "Model": "DNN(Softmax)",
    "Split": "Holdout",
    "Fold": "Holdout",
    "Complexity_Definition": "total_trainable_parameters",
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
# 5B. EXPORT HOLDOUT PREDICTIONS FOR FIGURE REPRODUCIBILITY
# ===============================================================

prediction_rows = []

for i, target in enumerate(yield_cols):
    for j in range(len(y_test)):
        prediction_rows.append({
            "Observation": j + 1,
            "Target": target,
            "Experimental": y_test[j, i],
            "RF_Predicted": y_pred_rf[j, i],
            "GBR_Predicted": y_pred_gbr[j, i],
            "DNN_Predicted": y_pred_dnn[j, i],
            "RF_Residual": y_pred_rf[j, i] - y_test[j, i],
            "GBR_Residual": y_pred_gbr[j, i] - y_test[j, i],
            "DNN_Residual": y_pred_dnn[j, i] - y_test[j, i],
        })

holdout_predictions = pd.DataFrame(prediction_rows)

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
# 8. Save all outputs
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
    holdout_predictions.to_excel(
        writer,
        sheet_name="Holdout_Predictions",
        index=False
    )

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
        {"Key": "DNN output structure", "Value": "single multi-output model (Keras Softmax)"},
        {"Key": "GBR output structure", "Value": "4 separate single-output models"},
    ])
    meta_df.to_excel(writer, sheet_name="Run_Metadata", index=False)

print("\n" + "=" * 70)
print(f"Done. Extended results saved to:\n{out_xlsx}")
print("=" * 70)

# ===============================================================
# 9. FIGURE S1 — PARITY PLOTS
#    Rows: RF, GBR, DNN
#    Columns: Biocrude, Aqueous, Gas, Solids
# ===============================================================

import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error

figure_dir = os.path.join(output_dir, "SI_Figures")
os.makedirs(figure_dir, exist_ok=True)

# Display names
product_names = ["Biocrude", "Aqueous", "Gas", "Solids"]

models_plot = {
    "RF": y_pred_rf,
    "GBR": y_pred_gbr,
    "DNN": y_pred_dnn,
}

# Publication-oriented defaults
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
})

fig, axes = plt.subplots(
    nrows=3,
    ncols=4,
    figsize=(10.5, 7.6),
    sharex=True,
    sharey=True
)

# All yields use the same physical scale
axis_min = 0
axis_max = 100

panel_letters = [
    "(a)", "(b)", "(c)", "(d)",
    "(e)", "(f)", "(g)", "(h)",
    "(i)", "(j)", "(k)", "(l)"
]

panel = 0

for row, (model_name, y_pred) in enumerate(models_plot.items()):

    for col, product_name in enumerate(product_names):

        ax = axes[row, col]

        y_obs = y_test[:, col]
        y_hat = y_pred[:, col]

        # Metrics
        r2 = r2_score(y_obs, y_hat)
        rmse = np.sqrt(mean_squared_error(y_obs, y_hat))

        # Experimental vs predicted observations
        ax.scatter(
            y_obs,
            y_hat,
            s=24,
            facecolors="none",
            edgecolors="black",
            linewidths=0.8,
            alpha=0.85
        )

        # Perfect-agreement line
        ax.plot(
            [axis_min, axis_max],
            [axis_min, axis_max],
            linestyle="--",
            linewidth=1.0,
            color="black"
        )

        ax.set_xlim(axis_min, axis_max)
        ax.set_ylim(axis_min, axis_max)

        ax.set_xticks(np.arange(0, 101, 20))
        ax.set_yticks(np.arange(0, 101, 20))

        # Keep panels square so slope = 1 visually
        ax.set_aspect("equal", adjustable="box")

        # Column headings
        if row == 0:
            ax.set_title(product_name, fontweight="bold")

        # Row labels
        if col == 0:
            ax.set_ylabel(
                f"{model_name}\nPredicted yield (wt%)",
                fontweight="bold"
            )

        # Bottom x-axis labels only
        if row == 2:
            ax.set_xlabel(
                "Experimental yield (wt%)",
                fontweight="bold"
            )

        # R2 and RMSE annotation
        ax.text(
            0.05,
            0.94,
            f"$R^2$ = {r2:.3f}\nRMSE = {rmse:.2f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8
        )

        # Panel letter
        ax.text(
            0.96,
            0.05,
            panel_letters[panel],
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=9,
            fontweight="bold"
        )

        # Clean publication appearance
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        panel += 1

fig.tight_layout(w_pad=1.0, h_pad=1.0)

# Vector format for manuscript/SI
fig.savefig(
    os.path.join(figure_dir, "Figure_S1_Parity_AllModels.pdf"),
    bbox_inches="tight"
)

# High-resolution raster copy
fig.savefig(
    os.path.join(figure_dir, "Figure_S1_Parity_AllModels.png"),
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)

print("Figure S1 parity plots saved.")

# ===============================================================
# 10. FIGURE S2 — RESIDUAL DIAGNOSTIC PLOTS
#     Residual = Predicted - Experimental
#     Rows: RF, GBR, DNN
#     Columns: Biocrude, Aqueous, Gas, Solids
# ===============================================================

fig, axes = plt.subplots(
    nrows=3,
    ncols=4,
    figsize=(10.5, 7.6),
    sharex=True
)

# Determine one symmetric residual range for all panels
all_residuals = []

for y_pred in models_plot.values():
    all_residuals.append((y_pred - y_test).ravel())

all_residuals = np.concatenate(all_residuals)

residual_limit = np.ceil(
    np.max(np.abs(all_residuals)) / 5
) * 5

panel = 0

for row, (model_name, y_pred) in enumerate(models_plot.items()):

    for col, product_name in enumerate(product_names):

        ax = axes[row, col]

        y_obs = y_test[:, col]
        residual = y_pred[:, col] - y_obs

        # Summary statistics
        mean_residual = np.mean(residual)
        mae = np.mean(np.abs(residual))

        ax.scatter(
            y_obs,
            residual,
            s=24,
            facecolors="none",
            edgecolors="black",
            linewidths=0.8,
            alpha=0.85
        )

        # Zero-error reference
        ax.axhline(
            0,
            linestyle="--",
            linewidth=1.0,
            color="black"
        )

        ax.set_xlim(0, 100)
        ax.set_ylim(-residual_limit, residual_limit)

        ax.set_xticks(np.arange(0, 101, 20))

        # Column headings
        if row == 0:
            ax.set_title(product_name, fontweight="bold")

        # Row-specific y labels
        if col == 0:
            ax.set_ylabel(
                f"{model_name}\nResidual (wt%)",
                fontweight="bold"
            )

        if row == 2:
            ax.set_xlabel(
                "Experimental yield (wt%)",
                fontweight="bold"
            )

        # Diagnostic annotation
        ax.text(
            0.05,
            0.94,
            f"Mean residual = {mean_residual:.2f}\nMAE = {mae:.2f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8
        )

        # Panel letter
        ax.text(
            0.96,
            0.05,
            panel_letters[panel],
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=9,
            fontweight="bold"
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        panel += 1

fig.tight_layout(w_pad=1.0, h_pad=1.0)

fig.savefig(
    os.path.join(figure_dir, "Figure_S2_Residuals_AllModels.pdf"),
    bbox_inches="tight"
)

fig.savefig(
    os.path.join(figure_dir, "Figure_S2_Residuals_AllModels.png"),
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)

print("Figure S2 residual plots saved.")