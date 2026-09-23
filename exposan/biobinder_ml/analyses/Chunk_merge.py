# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 18:31:35 2026

@author: aliah
"""

import os
import glob
import re
import pandas as pd
from collections import defaultdict

OUTPUT_DIR = "results"
MASTER_N_SAMPLES = 10000

pattern = os.path.join(OUTPUT_DIR, f"META_SHAP3_{MASTER_N_SAMPLES}_*_chunk-*-of-*.xlsx")
files = sorted(glob.glob(pattern))

if not files:
    raise FileNotFoundError(f"No chunk files found for pattern:\n{pattern}")

groups = defaultdict(list)

rx = re.compile(
    rf"META_SHAP3_{MASTER_N_SAMPLES}_(?P<feedstock>.+?)_(?P<config>CHCU_No_EC|DHCU_No_EC)_prefer-(?P<prefer>.+?)_chunk-\d+-of-\d+_.*\.xlsx$"
)

for f in files:
    name = os.path.basename(f)
    m = rx.match(name)
    if m:
        key = (m.group("feedstock"), m.group("config"), m.group("prefer"))
        groups[key].append(f)
    else:
        print(f"⚠️ Skipping unmatched file: {name}")

for (feedstock, config, prefer), group_files in groups.items():
    print("\n" + "=" * 80)
    print(f"Merging: feedstock={feedstock} | config={config} | prefer={prefer}")
    print("=" * 80)

    dfs = []
    for f in sorted(group_files):
        try:
            df = pd.read_excel(f)
            df["source_file"] = os.path.basename(f)
            dfs.append(df)
        except Exception as e:
            print(f"⚠️ Failed to read {f}: {e}")

    if not dfs:
        print("⚠️ No readable files in this group, skipping.")
        continue

    df_all = pd.concat(dfs, ignore_index=True)

    if "Sample_ID" in df_all.columns:
        df_all.sort_values("Sample_ID", inplace=True)
        df_all.drop_duplicates(subset=["Sample_ID"], keep="last", inplace=True)

    merged_path = os.path.join(
        OUTPUT_DIR,
        f"META_SHAP3_{MASTER_N_SAMPLES}_{feedstock}_{config}_prefer-{prefer}_MERGED.xlsx"
    )

    df_all.to_excel(merged_path, index=False)

    print(f"✅ Saved: {merged_path}")
    print(f"Rows: {len(df_all)}")
    if "OK" in df_all.columns:
        ok = df_all["OK"].astype(bool).sum()
        print(f"OK={ok}, FAIL={len(df_all)-ok}")