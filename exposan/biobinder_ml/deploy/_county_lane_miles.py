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

# --- 1. CONFIGURATION AND FILE PATHS ---
data_dir = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data"

# Inputs
lane_miles_path = os.path.join(data_dir, "county_lane_miles_summary.csv")
demographics_path = os.path.join(data_dir, "county_geo_demographics_summary.csv")

# FAF5 Truck Data Folder
faf_dir = os.path.join(data_dir, "All_Experimental_Disaggregation_Factors") #https://www.bts.gov/faf, Experimental County-Level Estimates
truck_origin_path = os.path.join(faf_dir, "truck_origin_factors.csv")
truck_destination_path = os.path.join(faf_dir, "truck_destination_factors.csv")

# Output Master Matrix Path
master_output_path = os.path.join(data_dir, "master_infrastructure_model_matrix.csv")

print("Initializing Master Infrastructure Data Stream Merge...")

# --- 2. LOAD AND AGGREGATE FAF5 TRUCK FREIGHT FACTORS ---
print("Loading and consolidating FAF5 regional truck disaggregation factors...")
if not os.path.exists(truck_origin_path) or not os.path.exists(truck_destination_path):
    raise FileNotFoundError(f"Ensure truck factors exist inside: {faf_dir}")

df_orig = pd.read_csv(truck_origin_path)
df_dest = pd.read_csv(truck_destination_path)

# Ensure county FIPS codes are consistently zero-padded 5-digit strings
df_orig['county_fips'] = df_orig['dms_orig_cnty'].astype(str).str.zfill(5)
df_dest['county_fips'] = df_dest['dms_dest_cnty'].astype(str).str.zfill(5)

# Group by county_fips to find total regional freight footprint across all commodity groups
# This calculates the absolute county-level truck weight multiplier profile
df_orig_agg = df_orig.groupby('county_fips', as_index=False)['f_orig'].sum()
df_dest_agg = df_dest.groupby('county_fips', as_index=False)['f_dest'].sum()

# Merge Origin and Destination footprints into a single Structural Wear Index per county
df_faf = pd.merge(df_orig_agg, df_dest_agg, on='county_fips', how='outer').fillna(0)
df_faf['Structural_Freight_Intensity'] = df_faf['f_orig'] + df_faf['f_dest']
df_faf_clean = df_faf[['county_fips', 'Structural_Freight_Intensity']].copy()


# --- 3. LOAD CORE INFRASTRUCTURE AND DEMOGRAPHIC DATA ---
print("Loading lane miles summary and demographic data lookup...")
df_lanes = pd.read_csv(lane_miles_path)
df_demo = pd.read_csv(demographics_path)

# Clean FIPS columns for a bulletproof merge key connection
df_lanes['county_fips'] = df_lanes['county_fips'].astype(str).str.zfill(5)
df_demo['county_fips'] = df_demo['county_fips'].astype(str).str.zfill(5)


# --- 4. STEPWISE SEQUENTIAL MERGING ---
print("Merging structural spatial data streams...")
# First merge: Snap Demographics to your core physical Lane Mile footprint
df_master = pd.merge(df_lanes, df_demo, on='county_fips', how='inner')

# Second merge: Left-join the FAF5 Truck freight wear coefficients (fill unlisted with 0)
df_master = pd.merge(df_master, df_faf_clean, on='county_fips', how='left').fillna({'Structural_Freight_Intensity': 0})


# --- 5. APPLY FHWA SPRING 2025 LONG-TERM BASELINE VMT METRICS ---
print("Calculating localized growth rates via Spring 2025 FHWA baseline parameters...")

# Baseline annual compound growth vectors from 2023-2053 (30 Year Baseline Outlook, https://www.fhwa.dot.gov/policyinformation/tables/vmt/vmt_forecast_sum.cfm
LDV_growth_baseline = 0.005       # 0.5% annual growth for Light-Duty Vehicles
SingleUnit_truck_growth = 0.020   # 2.0% annual growth for Single-Unit Trucks
Combination_truck_growth = 0.009  # 0.9% annual growth for Long-haul Combination Trucks

# Blended Freight Compound Rate for Highways (approx. 1.45% average blend)
heavy_truck_baseline = (SingleUnit_truck_growth + Combination_truck_growth) / 2.0

def calculate_localized_growth(row):
    """
    Differentiates the system tiers based on functional classifications:
    Tiers 1, 2, 3: Heavy Freight Corridors (scaled by FAF5 Structural Weight Index)
    Tiers 4, 5, 6, 7: Local/Collector Roads (scaled by Census Population Density trends)
    """
    sys_tier = float(row['functional_system'])
    
    if sys_tier in [1.0, 2.0, 3.0]:
        # Interstates & Arterials: Baseline truck demand scaled up or down by local freight intensity
        # Adding 1.0 shifts the allocation fraction to a positive localized scalar multiplier
        return heavy_truck_baseline * (1.0 + row['Structural_Freight_Intensity'])
    else:
        # Collectors & Locals: Regulated by light commuter travel scaled by county density tier
        # Prevents negative growth if population density is low; caps variance safely
        density_mod = max(0.1, min(row['Population_Density'] / 100.0, 3.0))
        return LDV_growth_baseline * density_mod

# Map the localized baseline rate column onto the matrix
df_master['Assigned_Annual_Growth_Rate'] = df_master.apply(calculate_localized_growth, axis=1)


# --- 6. EXPORT MASTER MATRIX ---
df_master.to_csv(master_output_path, index=False)
print(f"\nIntegration Complete! Master model matrix written to:\n{master_output_path}")

print("\nIntegrated Matrix Sample Output (First 5 Rows):")
print(df_master[['state_abbr', 'county_fips', 'functional_system', 'lane_miles', 'Population_Density', 'Structural_Freight_Intensity', 'Assigned_Annual_Growth_Rate']].head())