# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

Analyses:
    1. Pooled Spearman and Pearson correlation
    2. Pooled RF SHAP
    3. Per-case RF SHAP audit
    4. Comparison of per-case SHAP, pooled RF SHAP, and Spearman ranks

Final model specification:
    11 inputs
    800 RF trees
    min_samples_leaf = 2
    random_state = 42
    80/20 train/test split
    SHAP calculated on all retained observations

Output:
    RF_SHAP_avg_pooled_spearman_comparison_FINAL.xlsx
'''

import os
import re
import glob
import numpy as np
import pandas as pd

from scipy.stats import spearmanr
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

import shap


# ============================================================
# SETTINGS
# ============================================================

INPUT_DIR = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results"
OUTPUT_DIR = os.path.join(INPUT_DIR, "summary_outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTPUT_XLSX = os.path.join(
    OUTPUT_DIR,
    "RF_SHAP_avg_pooled_spearman_comparison_FINAL.xlsx"
)

META_PATTERN = os.path.join(
    INPUT_DIR,
    "META_SHAP3_10000_*_*HCU_No_EC*.xlsx"
)

TARGET = "IRR_pct"
OK_COL = "OK"

EXPECTED_FEEDSTOCKS = ["food", "sludge", "manure", "green"]
EXPECTED_CONFIGS = ["CHCU_No_EC", "DHCU_No_EC"]

FEATURES = [
    "Feedstock_price_$/tonne",
    "product_ratio_biobinder_over_biofuel",
    "Y_biocrude",
    "Solid content (w/w) %",
    "Temperature (C)",
    "income_tax",
    "Residence Time",
    "Natural_Gas_Price",
    "Electricity_Price",
    "Biofuel_price",
    "Uptime_Ratio",
]

DISPLAY_NAME = {
    "Feedstock_price_$/tonne": "Feedstock price",
    "product_ratio_biobinder_over_biofuel": "Biobinder to biofuel product ratio",
    "Y_biocrude": "Biocrude yield",
    "Solid content (w/w) %": "Solids content",
    "Temperature (C)": "Temperature",
    "income_tax": "Income tax",
    "Residence Time": "Residence time",
    "Natural_Gas_Price": "Natural gas price",
    "Electricity_Price": "Electricity price",
    "Biofuel_price": "Biofuel price",
    "Uptime_Ratio": "Uptime ratio",
}

RF_KWARGS = dict(
    n_estimators=800,
    random_state=42,
    min_samples_leaf=2,
    n_jobs=-1,
)

TEST_SIZE = 0.20
RANDOM_STATE = 42


# ============================================================
# HELPERS
# ============================================================

def parse_meta_filename(path):
    fname = os.path.basename(path)

    m = re.search(
        r"META_SHAP3_10000_(?P<feedstock>.+?)_(?P<config>CHCU_No_EC|DHCU_No_EC)",
        fname,
    )

    if m:
        return m.group("feedstock").lower(), m.group("config")

    return None, None


def clean_ok_true(df):
    if OK_COL not in df.columns:
        return df.copy()

    ok = df[OK_COL]

    if ok.dtype == bool:
        return df.loc[ok].copy()

    ok_str = ok.astype(str).str.strip().str.lower()

    return df.loc[
        ok_str.isin(["true", "1", "yes", "ok"])
    ].copy()


def find_latest_meta_files():
    files = glob.glob(META_PATTERN)

    if not files:
        raise FileNotFoundError(
            f"No META files found using pattern:\n{META_PATTERN}"
        )

    candidates = []

    for f in files:
        feedstock, config = parse_meta_filename(f)

        if feedstock not in EXPECTED_FEEDSTOCKS:
            continue

        if config not in EXPECTED_CONFIGS:
            continue

        candidates.append(
            {
                "path": f,
                "file": os.path.basename(f),
                "feedstock": feedstock,
                "config": config,
                "mtime": os.path.getmtime(f),
            }
        )

    if not candidates:
        raise ValueError(
            "META files were found, but none matched the expected "
            "feedstock/configuration cases."
        )

    candidate_df = pd.DataFrame(candidates)

    candidate_df = candidate_df.sort_values(
        ["feedstock", "config", "mtime"]
    )

    latest_df = (
        candidate_df
        .groupby(["feedstock", "config"], as_index=False)
        .tail(1)
        .copy()
    )

    latest_df = latest_df.sort_values(
        ["feedstock", "config"]
    ).reset_index(drop=True)

    expected_cases = {
        (feedstock, config)
        for feedstock in EXPECTED_FEEDSTOCKS
        for config in EXPECTED_CONFIGS
    }

    selected_cases = set(
        zip(latest_df["feedstock"], latest_df["config"])
    )

    missing_cases = sorted(expected_cases - selected_cases)

    if missing_cases:
        print("\nWARNING: Missing expected META cases:")
        for case in missing_cases:
            print("   ", case)

    if len(latest_df) != 8:
        print(
            f"\nWARNING: Expected 8 latest META files, "
            f"but found {len(latest_df)}."
        )

    print("\n" + "=" * 80)
    print("LATEST META FILES SELECTED")
    print("=" * 80)

    for _, row in latest_df.iterrows():
        timestamp = pd.to_datetime(
            row["mtime"],
            unit="s"
        ).strftime("%Y-%m-%d %H:%M:%S")

        print(
            f"{row['feedstock']:8s} | "
            f"{row['config']:11s} | "
            f"{timestamp} | "
            f"{row['file']}"
        )

    return latest_df, candidate_df


def load_latest_meta_files():
    latest_df, candidate_df = find_latest_meta_files()

    rows = []
    inventory = []

    for _, file_row in latest_df.iterrows():

        f = file_row["path"]
        feedstock = file_row["feedstock"]
        config = file_row["config"]

        print(f"\nLoading: {os.path.basename(f)}")

        df = pd.read_excel(f)
        n_raw = len(df)

        df = clean_ok_true(df)
        n_ok = len(df)

        required = FEATURES + [TARGET]

        missing = [
            c for c in required
            if c not in df.columns
        ]

        if missing:
            print("Missing columns:")
            for col in missing:
                print("   ", col)

            raise KeyError(
                f"{os.path.basename(f)} is missing required columns: "
                f"{missing}"
            )

        sub = df[required].copy()

        sub = sub.apply(
            pd.to_numeric,
            errors="coerce"
        )

        sub = sub.replace(
            [np.inf, -np.inf],
            np.nan
        )

        n_before_dropna = len(sub)

        sub = sub.dropna().copy()

        n_removed_nonfinite = n_before_dropna - len(sub)

        sub["feedstock"] = feedstock
        sub["config"] = config
        sub["case"] = f"{feedstock}_{config}"
        sub["source_file"] = os.path.basename(f)

        rows.append(sub)

        inventory.append(
            {
                "file": os.path.basename(f),
                "feedstock": feedstock,
                "config": config,
                "n_raw": n_raw,
                "n_ok": n_ok,
                "n_removed_nonfinite": n_removed_nonfinite,
                "n_used": len(sub),
                "modified_time": pd.to_datetime(
                    file_row["mtime"],
                    unit="s"
                ),
            }
        )

    if not rows:
        raise ValueError(
            "No usable META files remained after filtering."
        )

    pooled = pd.concat(
        rows,
        ignore_index=True
    )

    inventory_df = pd.DataFrame(inventory)

    print("\n" + "=" * 80)
    print("FINAL POOLED DATASET")
    print("=" * 80)

    print(f"Cases used: {pooled['case'].nunique()}")
    print(f"Rows used:  {len(pooled):,}")

    print("\nRows by case:")

    case_counts = (
        pooled
        .groupby(
            ["feedstock", "config"]
        )
        .size()
        .reset_index(name="n")
    )

    print(case_counts.to_string(index=False))

    return pooled, inventory_df, candidate_df


def corr_table(df, features, target):
    rows = []

    for feat in features:

        x = pd.to_numeric(
            df[feat],
            errors="coerce"
        )

        y = pd.to_numeric(
            df[target],
            errors="coerce"
        )

        valid = (
            x.notna()
            & y.notna()
            & np.isfinite(x)
            & np.isfinite(y)
        )

        x = x.loc[valid]
        y = y.loc[valid]

        pearson = np.nan
        spearman = np.nan
        p_spearman = np.nan

        if (
            len(x) > 2
            and x.nunique() > 1
            and y.nunique() > 1
        ):
            pearson = x.corr(
                y,
                method="pearson"
            )

            sp = spearmanr(x, y)

            spearman = sp.correlation
            p_spearman = sp.pvalue

        rows.append(
            {
                "Feature": feat,
                "Variable": DISPLAY_NAME.get(
                    feat,
                    feat
                ),
                "N": len(x),
                "Pearson_r": pearson,
                "abs_Pearson_r": (
                    abs(pearson)
                    if pd.notna(pearson)
                    else np.nan
                ),
                "Spearman_rho": spearman,
                "abs_Spearman_rho": (
                    abs(spearman)
                    if pd.notna(spearman)
                    else np.nan
                ),
                "Spearman_p": p_spearman,
            }
        )

    out = pd.DataFrame(rows)

    out["Pearson_rank"] = (
        out["abs_Pearson_r"]
        .rank(
            ascending=False,
            method="dense"
        )
        .astype("Int64")
    )

    out["Spearman_rank"] = (
        out["abs_Spearman_rho"]
        .rank(
            ascending=False,
            method="dense"
        )
        .astype("Int64")
    )

    return out


# def fit_rf_and_shap(
#     df,
#     features,
#     target,
#     label,
# ):
#     X = df[features].copy()
#     y = df[target].copy()

#     X_train, X_test, y_train, y_test = (
#         train_test_split(
#             X,
#             y,
#             test_size=TEST_SIZE,
#             random_state=RANDOM_STATE,
#         )
#     )

#     rf = RandomForestRegressor(
#         **RF_KWARGS
#     )

#     rf.fit(
#         X_train,
#         y_train
#     )

#     pred_train = rf.predict(X_train)
#     pred_test = rf.predict(X_test)

#     rmse_train = np.sqrt(
#         mean_squared_error(
#             y_train,
#             pred_train
#         )
#     )

#     rmse_test = np.sqrt(
#         mean_squared_error(
#             y_test,
#             pred_test
#         )
#     )

#     metrics = {
#         "model": label,
#         "n_total": len(df),
#         "n_train": len(X_train),
#         "n_test": len(X_test),
#         "R2_train": r2_score(
#             y_train,
#             pred_train
#         ),
#         "R2_test": r2_score(
#             y_test,
#             pred_test
#         ),
#         "RMSE_train": rmse_train,
#         "RMSE_test": rmse_test,
#         "MAE_train": mean_absolute_error(
#             y_train,
#             pred_train
#         ),
#         "MAE_test": mean_absolute_error(
#             y_test,
#             pred_test
#         ),
#     }

#     # --------------------------------------------------------
#     # SHAP on ALL retained observations
#     # --------------------------------------------------------

#     X_shap = X.copy()

#     explainer = shap.TreeExplainer(rf)

#     shap_values = explainer.shap_values(
#         X_shap
#     )

#     shap_values = np.asarray(
#         shap_values
#     )

#     shap_abs = np.abs(
#         shap_values
#     ).mean(axis=0)

#     shap_df = pd.DataFrame(
#         {
#             "Feature": features,
#             "Variable": [
#                 DISPLAY_NAME.get(f, f)
#                 for f in features
#             ],
#             f"{label}_mean_abs_SHAP": shap_abs,
#         }
#     )

#     shap_col = (
#         f"{label}_mean_abs_SHAP"
#     )

#     rank_col = (
#         f"{label}_SHAP_rank"
#     )

#     shap_df[rank_col] = (
#         shap_df[shap_col]
#         .rank(
#             ascending=False,
#             method="dense"
#         )
#         .astype(int)
#     )

#     shap_total = shap_df[
#         shap_col
#     ].sum()

#     shap_max = shap_df[
#         shap_col
#     ].max()

#     shap_df[
#         f"{label}_SHAP_share"
#     ] = (
#         shap_df[shap_col]
#         / shap_total
#     )

#     shap_df[
#         f"{label}_normalized_SHAP"
#     ] = (
#         shap_df[shap_col]
#         / shap_max
#     )

#     return (
#         rf,
#         shap_df,
#         pd.DataFrame([metrics]),
#     )

def fit_rf_and_shap(
    df,
    features,
    target,
    label,
):
    X = df[features].copy()
    y = df[target].copy()

    X_train, X_test, y_train, y_test = (
        train_test_split(
            X,
            y,
            test_size=TEST_SIZE,
            random_state=RANDOM_STATE,
        )
    )

    rf = RandomForestRegressor(
        **RF_KWARGS
    )

    rf.fit(
        X_train,
        y_train
    )

    pred_train = rf.predict(X_train)
    pred_test = rf.predict(X_test)

    rmse_train = np.sqrt(
        mean_squared_error(
            y_train,
            pred_train
        )
    )

    rmse_test = np.sqrt(
        mean_squared_error(
            y_test,
            pred_test
        )
    )

    metrics = {
        "model": label,
        "n_total": len(df),
        "n_train": len(X_train),
        "n_test": len(X_test),
        "R2_train": r2_score(
            y_train,
            pred_train
        ),
        "R2_test": r2_score(
            y_test,
            pred_test
        ),
        "RMSE_train": rmse_train,
        "RMSE_test": rmse_test,
        "MAE_train": mean_absolute_error(
            y_train,
            pred_train
        ),
        "MAE_test": mean_absolute_error(
            y_test,
            pred_test
        ),
    }

    # --------------------------------------------------------
    # SHAP on reproducible sample
    # --------------------------------------------------------

    SHAP_SAMPLE_N = 3000

    if len(X) > SHAP_SAMPLE_N:
        X_shap = X.sample(
            n=SHAP_SAMPLE_N,
            random_state=RANDOM_STATE
        )
    else:
        X_shap = X.copy()

    print(
        f"Calculating SHAP for {len(X_shap):,} "
        f"of {len(X):,} observations..."
    )

    explainer = shap.TreeExplainer(rf)

    shap_values = explainer.shap_values(
        X_shap,
        check_additivity=False
    )

    shap_values = np.asarray(
        shap_values
    )

    shap_abs = np.abs(
        shap_values
    ).mean(axis=0)

    shap_df = pd.DataFrame(
        {
            "Feature": features,
            "Variable": [
                DISPLAY_NAME.get(f, f)
                for f in features
            ],
            f"{label}_mean_abs_SHAP": shap_abs,
        }
    )

    shap_col = (
        f"{label}_mean_abs_SHAP"
    )

    rank_col = (
        f"{label}_SHAP_rank"
    )

    shap_df[rank_col] = (
        shap_df[shap_col]
        .rank(
            ascending=False,
            method="dense"
        )
        .astype(int)
    )

    shap_total = shap_df[
        shap_col
    ].sum()

    shap_max = shap_df[
        shap_col
    ].max()

    shap_df[
        f"{label}_SHAP_share"
    ] = (
        shap_df[shap_col]
        / shap_total
    )

    shap_df[
        f"{label}_normalized_SHAP"
    ] = (
        shap_df[shap_col]
        / shap_max
    )

    return (
        rf,
        shap_df,
        pd.DataFrame([metrics]),
    )
# ============================================================
# LOAD EXACTLY THE LATEST META FILE FOR EACH CASE
# ============================================================

pooled_df, inventory_df, candidate_df = (
    load_latest_meta_files()
)


# ============================================================
# 1. POOLED SPEARMAN / PEARSON
# ============================================================

print("\n" + "=" * 80)
print("POOLED CORRELATION ANALYSIS")
print("=" * 80)

corr_df = corr_table(
    pooled_df,
    FEATURES,
    TARGET
)

corr_df = corr_df.sort_values(
    "Spearman_rank"
).reset_index(drop=True)

print(
    corr_df[
        [
            "Variable",
            "N",
            "Spearman_rho",
            "Spearman_p",
            "Spearman_rank",
            "Pearson_r",
            "Pearson_rank",
        ]
    ].to_string(index=False)
)


# ============================================================
# 2. POOLED RF SHAP
# ============================================================

print("\n" + "=" * 80)
print("POOLED RF SHAP")
print("=" * 80)

(
    pooled_rf,
    pooled_shap_df,
    pooled_metrics_df,
) = fit_rf_and_shap(
    pooled_df,
    FEATURES,
    TARGET,
    label="pooled_RF",
)

pooled_shap_df = (
    pooled_shap_df
    .sort_values(
        "pooled_RF_SHAP_rank"
    )
    .reset_index(drop=True)
)

print(
    pooled_shap_df[
        [
            "Variable",
            "pooled_RF_mean_abs_SHAP",
            "pooled_RF_normalized_SHAP",
            "pooled_RF_SHAP_share",
            "pooled_RF_SHAP_rank",
        ]
    ].to_string(index=False)
)


# # ============================================================
# # 3. PER-CASE RF SHAP AUDIT
# # ============================================================

# print("\n" + "=" * 80)
# print("PER-CASE RF SHAP AUDIT")
# print("=" * 80)

# case_shap_rows = []
# case_metric_rows = []

# for case, g in pooled_df.groupby(
#     "case",
#     sort=True
# ):

#     if len(g) < 100:
#         print(
#             f"Skipping small case: "
#             f"{case}, n={len(g)}"
#         )
#         continue

#     print(
#         f"\nRunning {case}: "
#         f"n={len(g):,}"
#     )

#     (
#         rf,
#         shap_df,
#         metrics_df,
#     ) = fit_rf_and_shap(
#         g,
#         FEATURES,
#         TARGET,
#         label="case_RF",
#     )

#     feedstock = g[
#         "feedstock"
#     ].iloc[0]

#     config = g[
#         "config"
#     ].iloc[0]

#     shap_df["case"] = case
#     shap_df["feedstock"] = feedstock
#     shap_df["config"] = config

#     metrics_df["case"] = case
#     metrics_df["feedstock"] = feedstock
#     metrics_df["config"] = config

#     case_shap_rows.append(
#         shap_df
#     )

#     case_metric_rows.append(
#         metrics_df
#     )


# case_shap_df = pd.concat(
#     case_shap_rows,
#     ignore_index=True
# )

# case_metrics_df = pd.concat(
#     case_metric_rows,
#     ignore_index=True
# )


# # ============================================================
# # 4. AVERAGE PER-CASE SHAP
# # ============================================================

# avg_case_shap_df = (
#     case_shap_df
#     .groupby(
#         ["Feature", "Variable"],
#         as_index=False
#     )
#     .agg(
#         avg_case_RF_mean_abs_SHAP=(
#             "case_RF_mean_abs_SHAP",
#             "mean"
#         ),
#         avg_case_RF_mean_rank=(
#             "case_RF_SHAP_rank",
#             "mean"
#         ),
#         avg_case_RF_median_rank=(
#             "case_RF_SHAP_rank",
#             "median"
#         ),
#         n_cases=(
#             "case",
#             "nunique"
#         ),
#     )
# )

# avg_case_shap_df[
#     "avg_case_RF_SHAP_rank"
# ] = (
#     avg_case_shap_df[
#         "avg_case_RF_mean_abs_SHAP"
#     ]
#     .rank(
#         ascending=False,
#         method="dense"
#     )
#     .astype(int)
# )

# avg_case_max = (
#     avg_case_shap_df[
#         "avg_case_RF_mean_abs_SHAP"
#     ].max()
# )

# avg_case_total = (
#     avg_case_shap_df[
#         "avg_case_RF_mean_abs_SHAP"
#     ].sum()
# )

# avg_case_shap_df[
#     "avg_case_RF_normalized_SHAP"
# ] = (
#     avg_case_shap_df[
#         "avg_case_RF_mean_abs_SHAP"
#     ]
#     / avg_case_max
# )

# avg_case_shap_df[
#     "avg_case_RF_SHAP_share"
# ] = (
#     avg_case_shap_df[
#         "avg_case_RF_mean_abs_SHAP"
#     ]
#     / avg_case_total
# )

# avg_case_shap_df = (
#     avg_case_shap_df
#     .sort_values(
#         "avg_case_RF_SHAP_rank"
#     )
#     .reset_index(drop=True)
# )


# ============================================================
# 5. FINAL COMPARISON
# ============================================================

# final = (
#     avg_case_shap_df
#     .merge(
#         pooled_shap_df,
#         on=[
#             "Feature",
#             "Variable"
#         ],
#         how="left",
#     )
#     .merge(
#         corr_df[
#             [
#                 "Feature",
#                 "N",
#                 "Pearson_r",
#                 "abs_Pearson_r",
#                 "Pearson_rank",
#                 "Spearman_rho",
#                 "abs_Spearman_rho",
#                 "Spearman_p",
#                 "Spearman_rank",
#             ]
#         ],
#         on="Feature",
#         how="left",
#     )
# )

# final[
#     "rank_shift_avg_to_pooled_RF"
# ] = (
#     final[
#         "pooled_RF_SHAP_rank"
#     ]
#     - final[
#         "avg_case_RF_SHAP_rank"
#     ]
# )

# final[
#     "rank_shift_avg_to_spearman"
# ] = (
#     final[
#         "Spearman_rank"
#     ]
#     - final[
#         "avg_case_RF_SHAP_rank"
#     ]
# )

# final[
#     "rank_shift_pooled_RF_to_spearman"
# ] = (
#     final[
#         "Spearman_rank"
#     ]
#     - final[
#         "pooled_RF_SHAP_rank"
#     ]
# )

# final = (
#     final
#     .sort_values(
#         "avg_case_RF_SHAP_rank"
#     )
#     .reset_index(drop=True)
# )

# ============================================================
# FINAL COMPARISON: POOLED RF SHAP VS POOLED SPEARMAN
# ============================================================

final = (
    pooled_shap_df
    .merge(
        corr_df[
            [
                "Feature",
                "N",
                "Pearson_r",
                "abs_Pearson_r",
                "Pearson_rank",
                "Spearman_rho",
                "abs_Spearman_rho",
                "Spearman_p",
                "Spearman_rank",
            ]
        ],
        on="Feature",
        how="left",
    )
)

final[
    "rank_shift_pooled_RF_to_spearman"
] = (
    final["Spearman_rank"]
    - final["pooled_RF_SHAP_rank"]
)

final = (
    final
    .sort_values(
        "pooled_RF_SHAP_rank"
    )
    .reset_index(drop=True)
)
# ============================================================
# 6. SI TABLE S23
# ============================================================

table_S23 = (
    corr_df[
        [
            "Variable",
            "Spearman_rho",
            "abs_Spearman_rho",
            "Spearman_p",
            "Spearman_rank",
            "N",
        ]
    ]
    .copy()
)

table_S23.columns = [
    "Model input",
    "Spearman rho",
    "Absolute rho",
    "p value",
    "Rank",
    "N",
]

table_S23 = (
    table_S23
    .sort_values("Rank")
    .reset_index(drop=True)
)


# ============================================================
# 7. SI TABLE S24
# ============================================================

table_S24 = (
    final[
        [
            "Variable",
            "avg_case_RF_SHAP_rank",
            "pooled_RF_SHAP_rank",
            "Spearman_rank",
        ]
    ]
    .copy()
)

table_S24.columns = [
    "Model input",
    "Eight case SHAP rank",
    "Pooled RF SHAP rank",
    "Pooled Spearman rank",
]

table_S24 = (
    table_S24
    .sort_values(
        "Eight case SHAP rank"
    )
    .reset_index(drop=True)
)


# ============================================================
# 8. RANK AGREEMENT
# ============================================================

rank_cols = [
    "avg_case_RF_SHAP_rank",
    "pooled_RF_SHAP_rank",
    "Spearman_rank",
    "Pearson_rank",
]

agreement_rows = []

for i, a in enumerate(rank_cols):
    for b in rank_cols[i + 1:]:

        rho = spearmanr(
            final[a],
            final[b]
        ).correlation

        agreement_rows.append(
            {
                "rank_method_1": a,
                "rank_method_2": b,
                "Spearman_rank_agreement": rho,
            }
        )

agreement_df = pd.DataFrame(
    agreement_rows
)


# ============================================================
# 9. SAVE EXCEL
# ============================================================

with pd.ExcelWriter(
    OUTPUT_XLSX,
    engine="openpyxl"
) as writer:

    table_S23.to_excel(
        writer,
        sheet_name="Table_S23_Spearman",
        index=False
    )

    table_S24.to_excel(
        writer,
        sheet_name="Table_S24_Ranks",
        index=False
    )

    final.to_excel(
        writer,
        sheet_name="final_rank_comparison",
        index=False
    )

    corr_df.to_excel(
        writer,
        sheet_name="pooled_corr",
        index=False
    )

    pooled_shap_df.to_excel(
        writer,
        sheet_name="pooled_RF_SHAP",
        index=False
    )

    # avg_case_shap_df.to_excel(
    #     writer,
    #     sheet_name="avg_case_RF_SHAP",
    #     index=False
    # )

    # case_shap_df.to_excel(
    #     writer,
    #     sheet_name="per_case_RF_SHAP",
    #     index=False
    # )

    pooled_metrics_df.to_excel(
        writer,
        sheet_name="pooled_RF_metrics",
        index=False
    )

    # case_metrics_df.to_excel(
    #     writer,
    #     sheet_name="per_case_RF_metrics",
    #     index=False
    # )

    agreement_df.to_excel(
        writer,
        sheet_name="rank_agreement",
        index=False
    )

    inventory_df.to_excel(
        writer,
        sheet_name="selected_file_inventory",
        index=False
    )

    candidate_df.to_excel(
        writer,
        sheet_name="all_candidate_files",
        index=False
    )

    pooled_df[
        FEATURES
        + [
            TARGET,
            "feedstock",
            "config",
            "case",# SHAP on ALL retained observations
            "source_file",
        ]
    ].to_excel(
        writer,
        sheet_name="pooled_selected_data",
        index=False
    )


# ============================================================
# 10. PRINT FINAL RESULTS
# ============================================================

print("\n" + "=" * 80)
print("TABLE S23: POOLED SPEARMAN")
print("=" * 80)

print(
    table_S23.to_string(
        index=False
    )
)

print("\n" + "=" * 80)
print("TABLE S24: RANK COMPARISON")
print("=" * 80)

print(
    table_S24.to_string(
        index=False
    )
)

print("\n" + "=" * 80)
print("POOLED RF METRICS")
print("=" * 80)

print(
    pooled_metrics_df.to_string(
        index=False
    )
)

print("\n" + "=" * 80)
print("PER-CASE RF METRICS")
print("=" * 80)

# print(
#     case_metrics_df[
#         [
#             "case",
#             "n_total",
#             "R2_test",
#             "RMSE_test",
#             "MAE_test",
#         ]
#     ].to_string(
#         index=False
#     )
# )

print("\n" + "=" * 80)
print("CHECKS")
print("=" * 80)

print(
    f"Number of selected META files: "
    f"{len(inventory_df)}"
)

print(
    f"Number of unique cases: "
    f"{pooled_df['case'].nunique()}"
)

print(
    f"Number of pooled rows: "
    f"{len(pooled_df):,}"
)

print(
    f"Number of model inputs: "
    f"{len(FEATURES)}"
)

if len(FEATURES) != 11:
    raise RuntimeError(
        "Expected exactly 11 model inputs."
    )

if pooled_df["case"].nunique() != 8:
    raise RuntimeError(
        "Expected exactly 8 feedstock/configuration cases."
    )

print("\nAll final checks passed.")

print(
    f"\nSaved workbook:\n"
    f"{OUTPUT_XLSX}"
)