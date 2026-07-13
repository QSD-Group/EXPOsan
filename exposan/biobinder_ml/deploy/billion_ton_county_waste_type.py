# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os
import pandas as pd
import numpy as np

# --- 1. SET CUSTOM INPUT/OUTPUT FILE PATHS ---
csv_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\Billion Ton\billionton_23_wastes_download20251002-021830.csv"
output_excel_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\Billion Ton\county_waste_scenarios_comparison_near_term.xlsx"

# Guarantee target directory layout structure exists
os.makedirs(os.path.dirname(output_excel_path), exist_ok=True)

# --- 2. LOAD DATASET ---
print("Loading Billion-Ton 2023 dataset...")
df = pd.read_csv(csv_path)

# Enforce clean string-based 5-digit formatting for FIPS codes
df['fips'] = df['fips'].astype(str).str.zfill(5)

# Filter for the standard optimization scenario baseline
df_filtered = df[df['scenario_name'] == 'near-term'].copy()

# --- 3. DEFINE COMPREHENSIVE CONFIGURATION DICTIONARY ---
disposition_rules = {
    "Food Wastes":                             {"moisture": 0.75,  "landfill_frac": 0.60},
    "Manure Wastes":                           {"moisture": 0.74,  "landfill_frac": 0.00},
    "Plastics":                                {"moisture": 0.00,  "landfill_frac": 0.75},
    "Paper, Cardboard & MSW Organics":         {"moisture": 0.055, "landfill_frac": 0.32},
    "MSW Organics other":                      {"moisture": 0.50,  "landfill_frac": 0.80},
    "Green, Yard, Wood & Agricultural Wastes": {"moisture": 0.55,  "landfill_frac": 0.28}
}

def assign_htl_group(resource_name):
    res = str(resource_name).strip()
    
    # Group 1: Food Wastes
    if res in ["Food waste, residential", "Food waste, nonresidential"]:
        return "Food Wastes"
    # Group 2: Manure Wastes
    elif res in ["Manure, beef", "Manure, dairy", "Manure, swine", "Manure, poultry"]:
        return "Manure Wastes"
    # Group 3: Plastics
    elif res == "Plastics":
        return "Plastics"
    # Group 4: Paper & Cardboard (Matches dataset string context)
    elif res in ["Paper", "Cardboard", "Paper and paperboard"]:
        return "Paper, Cardboard & MSW Organics"
    # Group 5: MSW Organics other (Captures text and leather elements)
    elif res in ["Rubber and leather", "Textiles"] or "MSW Organics" in res:
        return "MSW Organics other"
    # Group 6: Green, Yard, Wood & Agricultural Wastes
    elif any(kw in res.lower() for kw in ["yard waste", "wood", "residue", "stover", "straw"]):
        return "Green, Yard, Wood & Agricultural Wastes"
        
    return None

df_filtered['HTL_Group'] = df_filtered['resource'].apply(assign_htl_group)
df_valid = df_filtered[df_filtered['HTL_Group'].notna()].copy()

# --- 4. EXECUTE COHERENT COMPILATION LOOP ---
print("Aggregating annual inventory metrics...")
county_grouped = df_valid.groupby(['fips', 'county', 'state', 'HTL_Group'])['production_total'].sum().reset_index()

def compute_all_metrics(row):
    group = row['HTL_Group']
    dry_total = float(row['production_total'])
    
    rules = disposition_rules[group]
    m_frac = rules["moisture"]
    lf_frac = rules["landfill_frac"]
    
    # Calculate mass balances
    wet_total = dry_total / (1.0 - m_frac)
    dry_landfill = dry_total * lf_frac
    wet_landfill = wet_total * lf_frac
    
    return pd.Series([dry_total, wet_total, dry_landfill, wet_landfill])

metrics_matrix = county_grouped.apply(compute_all_metrics, axis=1)
metrics_matrix.columns = ['Dry_Total', 'Wet_Total', 'Dry_Landfill', 'Wet_Landfill']
county_grouped = pd.concat([county_grouped, metrics_matrix], axis=1)

# --- 5. GENERATE VARIATION 1: ALL WASTE AVAILABLE ---
print("Structuring Variation 1 sheet format...")
pivot_dry_all = county_grouped.pivot(index=['fips', 'county', 'state'], columns='HTL_Group', values='Dry_Total').fillna(0.0)
pivot_wet_all = county_grouped.pivot(index=['fips', 'county', 'state'], columns='HTL_Group', values='Wet_Total').fillna(0.0)

pivot_dry_all.columns = [f"{col}_Dry_Tons_Year" for col in pivot_dry_all.columns]
pivot_wet_all.columns = [f"{col}_Wet_Tons_Year" for col in pivot_wet_all.columns]

sheet1_df = pd.concat([pivot_dry_all, pivot_wet_all], axis=1).reset_index()
sheet1_df['Total_Dry_Tons_Year'] = sheet1_df[[c for c in sheet1_df.columns if "_Dry_Tons_Year" in c]].sum(axis=1)
sheet1_df['Total_Wet_Tons_Year'] = sheet1_df[[c for c in sheet1_df.columns if "_Wet_Tons_Year" in c]].sum(axis=1)

# --- 6. GENERATE VARIATION 2: ONLY LANDFILLED PORTION ---
print("Structuring Variation 2 sheet format...")
pivot_dry_lf = county_grouped.pivot(index=['fips', 'county', 'state'], columns='HTL_Group', values='Dry_Landfill').fillna(0.0)
pivot_wet_lf = county_grouped.pivot(index=['fips', 'county', 'state'], columns='HTL_Group', values='Wet_Landfill').fillna(0.0)

pivot_dry_lf.columns = [f"{col}_Landfilled_Dry_Tons_Year" for col in pivot_dry_lf.columns]
pivot_wet_lf.columns = [f"{col}_Landfilled_Wet_Tons_Year" for col in pivot_wet_lf.columns]

sheet2_df = pd.concat([pivot_dry_lf, pivot_wet_lf], axis=1).reset_index()
sheet2_df['Total_Landfilled_Dry_Tons_Year'] = sheet2_df[[c for c in sheet2_df.columns if "_Landfilled_Dry_Tons_Year" in c]].sum(axis=1)
sheet2_df['Total_Landfilled_Wet_Tons_Year'] = sheet2_df[[c for c in sheet2_df.columns if "_Landfilled_Wet_Tons_Year" in c]].sum(axis=1)

# --- 7. WRITE BOTH TO EXCEL WORKBOOK TABS ---
print("Compiling dataframes into final multi-tab Excel spreadsheet workbook...")
with pd.ExcelWriter(output_excel_path) as writer:
    sheet1_df.to_excel(writer, sheet_name="All Waste Available", index=False)
    sheet2_df.to_excel(writer, sheet_name="Only Landfilled Portion", index=False)

print(f"\nProcessing Success! Comprehensive file exported to:\n{output_excel_path}")