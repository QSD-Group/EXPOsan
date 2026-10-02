# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 17:56:55 2025

@author: aliah
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer
import os
# ===============================================================
# 1. Load dataset
# ===============================================================

output_dir = "results"
os.makedirs(output_dir, exist_ok=True)
# Load your dataset
df = pd.read_excel(os.path.join(output_dir,"HTL_all_yield_normalized.xlsx"))
# Define yield columns
yield_cols = ["Biocrude wt%", "Aqueous wt%", "Gas wt%", "Solids wt%"]
test_size=0.10

# ===============================================================
# 2. Median imputation for numeric columns
# ===============================================================
num_cols = df.select_dtypes(include=[np.number]).columns
imp = SimpleImputer(strategy="median")
df[num_cols] = imp.fit_transform(df[num_cols])

for col in yield_cols:
    n_unique = df[col].nunique()
    n_bins = min(4, max(3, n_unique // 10 if n_unique > 10 else n_unique))
    try:
        df[f"{col}_bin"] = pd.qcut(df[col], q=n_bins, labels=False, duplicates="drop")
    except ValueError:
        df[f"{col}_bin"] = pd.cut(df[col], bins=n_bins, labels=False, include_lowest=True)

# Composite stratification label
combo_col = "yield_bin_combo"
df[combo_col] = df[[f"{col}_bin" for col in yield_cols]].astype(str).agg("_".join, axis=1)

# ==============================
# Make stratification feasible
# (1) Merge classes with <2 rows into 'Other'
# (2) If classes > n_test_samples, iteratively merge smallest classes into 'Other'
# ==============================
def make_stratification_feasible(df, combo_col, test_size):
    N = len(df)
    n_test = int(np.floor(test_size * N))
    # Step 1: merge classes with <2 rows (need ≥1 for train and ≥1 for test)
    counts = df[combo_col].value_counts()
    rare = counts[counts < 2].index
    if len(rare) > 0:
        df.loc[df[combo_col].isin(rare), combo_col] = "Other"

    # Recompute counts
    counts = df[combo_col].value_counts()

    # Step 2: ensure number of classes ≤ n_test
    # If too many classes, merge smallest classes into 'Other' until feasible
    def ensure_class_count_ok(df, counts, combo_col, n_test):
        classes = counts.index.tolist()
        while len(classes) > n_test:
            # Find the smallest (non-'Other') class and merge it into 'Other'
            smallest = counts[counts.index != "Other"].sort_values().index[0]
            df.loc[df[combo_col] == smallest, combo_col] = "Other"
            counts = df[combo_col].value_counts()
            classes = counts.index.tolist()
        return df, counts

    df, counts = ensure_class_count_ok(df, counts, combo_col, n_test)

    # Also ensure after merging that 'Other' has at least 2 rows (it will)
    return df

df = make_stratification_feasible(df, combo_col, test_size)

# ==============================
# Split
# ==============================
train_df, test_df = train_test_split(
    df, test_size=test_size, random_state=42, stratify=df[combo_col]
)

# Cleanup helper columns
drop_cols = [f"{col}_bin" for col in yield_cols] + [combo_col]
train_df = train_df.drop(columns=drop_cols)
test_df  = test_df.drop(columns=drop_cols)

out_path = "HTL_all_yield_normalized_split_stratified_10pct.xlsx"
with pd.ExcelWriter(out_path, engine="openpyxl") as w:
    train_df.to_excel(w, "Train", index=False)
    test_df.to_excel(w,  "Test",  index=False)

print("✅ Saved:", out_path)
print(f"Train: {len(train_df)} rows, Test: {len(test_df)} rows")