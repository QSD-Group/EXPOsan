# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import pandas as pd
import os

# --- 1. FILE PATHS AND CONFIGURATION ---
pop_csv_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\co-est2025-alldata.csv" #https://www2.census.gov/programs-surveys/popest/datasets/2020-2025/counties/totals/
# Pointing to the local copy you just saved:
gazetteer_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\2025_Gaz_counties_national.txt"   #https://www2.census.gov/geo/docs/maps-data/data/gazetteer/2025_Gazetteer/

output_csv_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\county_geo_demographics_summary.csv"

print("Starting demographic and geographic data integration...")

# --- 2. LOAD AND CLEAN LOCAL CENSUS POPULATION DATA ---
if not os.path.exists(pop_csv_path):
    raise FileNotFoundError(f"Could not locate the population file at: {pop_csv_path}")

df_pop = pd.read_csv(pop_csv_path, encoding='latin-1')

# Synthesize 5-digit county FIPS
df_pop['county_fips'] = df_pop['STATE'].astype(str).str.zfill(2) + df_pop['COUNTY'].astype(str).str.zfill(3)

# Filter for SUMLEV 050 (True Counties)
df_pop = df_pop[df_pop['SUMLEV'] == 50].copy()

# Determine active target population column
target_pop_col = 'POPESTIMATE2025' 
if target_pop_col not in df_pop.columns:
    available_est = [col for col in df_pop.columns if 'POPESTIMATE' in col]
    print(f"Warning: '{target_pop_col}' not found. Available estimates: {available_est}")
    target_pop_col = available_est[-1]
    print(f"Using fallback population column: {target_pop_col}")

df_pop_clean = df_pop[['county_fips', target_pop_col]].copy()
df_pop_clean = df_pop_clean.rename(columns={target_pop_col: 'Population'})


# --- 3. LOAD LOCAL CENSUS GAZETTEER DATA ---
if not os.path.exists(gazetteer_path):
    raise FileNotFoundError(f"Please download the Gazetteer file and place it at: {gazetteer_path}")

print(f"Loading local US Census Gazetteer data from: {gazetteer_path}")

# Added encoding='latin-1' to safely read text characters in county names
df_gaz = pd.read_csv(gazetteer_path, sep='|', encoding='latin-1')

# Strip white spaces from headers
df_gaz.columns = df_gaz.columns.str.strip()

# Enforce zero-padded 5-digit string alignment on GEOID
df_gaz['county_fips'] = df_gaz['GEOID'].astype(str).str.zfill(5)

# Extract spatial attributes
df_land = df_gaz[['county_fips', 'NAME', 'USPS', 'ALAND_SQMI']].copy()
df_land = df_land.rename(columns={'NAME': 'County_Name', 'USPS': 'State_Abbr', 'ALAND_SQMI': 'Land_Area_SqMi'})


# --- 4. COMPUTE GEOGRAPHIC JOIN AND POPULATION DENSITY ---
print("Merging datasets and calculating localized density profiles...")
df_geo_features = pd.merge(df_pop_clean, df_land, on='county_fips', how='inner')

# Guard against zero-division
df_geo_features = df_geo_features[df_geo_features['Land_Area_SqMi'] > 0].copy()

# Compute Density
df_geo_features['Population_Density'] = df_geo_features['Population'] / df_geo_features['Land_Area_SqMi']


# --- 5. DATA SANITY CHECK & EXPORT ---
print(f"\nProcessing Complete! Unified Data Matrix Summary:")
print(f"Total Unique Counties Monitored: {len(df_geo_features)}")
print(f"Unique states/territories mapped: {df_geo_features['State_Abbr'].nunique()}")

sc_check = df_geo_features[df_geo_features['State_Abbr'] == 'SC']
dc_check = df_geo_features[df_geo_features['State_Abbr'] == 'DC']
print(f" -> South Carolina Validation: Mapped {len(sc_check)} counties.")
print(f" -> District of Columbia Validation: Mapped {len(dc_check)} jurisdictions.")

df_geo_features.to_csv(output_csv_path, index=False)
print(f"\nMaster geo-demographic summary successfully exported to:\n{output_csv_path}")

print("\nSample Output View:")
print(df_geo_features.head())