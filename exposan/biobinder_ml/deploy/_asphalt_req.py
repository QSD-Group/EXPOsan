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
import numpy as np

# --- 1. CONFIGURATION AND FILE PATHS ---
data_dir = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data"
master_input_path = os.path.join(data_dir, "master_infrastructure_model_matrix.csv")
binder_output_path = os.path.join(data_dir, "county_asphalt_binder_demand_projections.csv")

print("Initializing Demand-Side Asphalt Binder Projection Model...")

if not os.path.exists(master_input_path):
    raise FileNotFoundError(f"Missing master infrastructure matrix file at: {master_input_path}. Run the merge script first.")

df = pd.read_csv(master_input_path)
df['county_fips'] = df['county_fips'].astype(str).str.zfill(5)


# --- 2. STRUCTURAL MODELING COEFFICIENTS & CLIMATE MODIFIERS ---
# [Ref A] Lane dimensions and geometry: AASHTO Policy on Geometric Design of Highways and Streets.
# [Ref B] Bulk compacted HMA density (~145 lbs/ft³): Asphalt Institute MS-2 Asphalt Mix Design Methods.
# [Ref C] Standard baseline binder content (5.0% total weight): NCHRP Report 815 / State DOT HMA Mix Specifications.
# [Ref D] Design thickness overlays & structural recurrence intervals (R): Federal Highway Administration (FHWA) 
#          Pavement Management Guide & Life-Cycle Cost Analysis (LCCA) standard pavement performance curves.
# [Ref E] Environmental/Climate Zone Pavement Life-Cycle Performance Modifiers:
#          Sourced from FHWA Long-Term Pavement Performance (LTPP) Climatic Zone Data Framework.
#          Official Documentation & Data Access Link: https://www.fhwa.dot.gov/pavement/pavementpolicy/linkages/design_pm/

LANE_WIDTH_FT = 12.0          # Standard US Traffic Lane Width (AASHTO)
LANE_MILE_LEN_FT = 5280.0     # Linear feet in one mile
LANE_MILE_AREA_SQFT = LANE_WIDTH_FT * LANE_MILE_LEN_FT  # 63,360 sq ft per lane-mile
HMA_DENSITY_LBS_CUFT = 145.0  # Compacted Bulk Density of Hot Mix Asphalt (Asphalt Institute)
BINDER_CONTENT_PCT = 0.05     # 5.0% AC Binder by weight of total mix (Standard Target Baseline)

# Mapping Functional Classifications (HPMS Tiers) to structural baseline pavement lifespan [Ref D]
# Replaced fixed "recurrence_R" with "base_R" to permit dynamic regional climate modifications.
structural_design_matrix = {
    1.0: {"name": "Interstates",              "thickness_in": 3.0, "base_R": 12.0}, 
    2.0: {"name": "Freeways & Expressways",   "thickness_in": 2.5, "base_R": 14.0}, 
    3.0: {"name": "Other Principal Arterials", "thickness_in": 2.5, "base_R": 14.0},
    4.0: {"name": "Minor Arterials",          "thickness_in": 2.0, "base_R": 16.5}, 
    5.0: {"name": "Major Collectors",         "thickness_in": 2.0, "base_R": 16.5},
    6.0: {"name": "Minor Collectors",         "thickness_in": 1.5, "base_R": 20.0}, 
    7.0: {"name": "Local Roads",              "thickness_in": 1.5, "base_R": 20.0}
}

# CLIMATE LIFE-CYCLE MODIFIERS [Ref E]/NCAT Report 18-02, 15% variance assumption
# Adjusts pavement lifespan based on weather distress (Ice/Potholes vs Dry/Heat oxidation)
# Multipliers < 1.0 accelerate failure intervals (boosting annual binder tonnage required)
CLIMATE_LIFESPAN_MODIFIER = {
    'Wet-Freeze': 0.83,      # Shortens lifespan by ~17% (Heavy frost heaving / spring thaws)
    'Dry-Freeze': 0.95,      # Shortens lifespan by ~5% (Thermal transverse cracking)
    'Dry-Non-Freeze': 1.05,  # Extends lifespan by ~5% (Optimal solid subgrade support, high oxidation)
    'Wet-Non-Freeze': 1.15   # Extends lifespan by ~15% (No freeze cycles, standard stripping vulnerability)
}


# --- 3. EXECUTE DEMAND-SIDE STRUCTURAL PAVEMENT MATHEMATICS ---
print("Executing layer volume-to-mass conversions and annualized climate-maintenance loops...")

def calculate_baseline_annual_binder_tons(row):
    f_sys = float(row['functional_system'])
    lane_miles = float(row['lane_miles'])
    
    if f_sys not in structural_design_matrix or lane_miles <= 0:
        return 0.0
    
    # Extract baseline structural geometry constraints [Ref D]
    design = structural_design_matrix[f_sys]
    thickness_ft = design["thickness_in"] / 12.0
    
    # DYNAMIC CLIMATE CALIBRATION LOOP [Ref E]
    # Checks for 'Climate_Zone' column in dataset, defaults safely to 'Wet-Freeze' if unlisted
    county_climate = row.get('Climate_Zone', 'Wet-Freeze')
    climate_factor = CLIMATE_LIFESPAN_MODIFIER.get(county_climate, 1.0)
    
    # Calculate the localized recurrence cycle adjusted by environmental strain
    R_years = design["base_R"] * climate_factor
    
    # Step A: Compute physical material volume per lane-mile (cubic feet)
    volume_per_lane_mile = LANE_MILE_AREA_SQFT * thickness_ft
    
    # Step B: Convert physical volume to mass tonnage (US short tons)
    total_hma_tons_per_lane_mile = (volume_per_lane_mile * HMA_DENSITY_LBS_CUFT) / 2000.0
    
    # Step C: Isolate liquid asphalt binder component weight
    binder_tons_per_lane_mile = total_hma_tons_per_lane_mile * BINDER_CONTENT_PCT
    
    # Step D: Apply annualized lifecycle maintenance recurrence interval (Annualized Baseline Demand)
    # The climate-adjusted 'R_years' now shifts material demand dynamically based on the region.
    annual_pavement_cycle_demand = (lane_miles * binder_tons_per_lane_mile) / R_years
    
    return annual_pavement_cycle_demand

# Compute baseline 2023 structural lifecycle target demand (Short Tons)
df['Binder_Demand_2023'] = df.apply(calculate_baseline_annual_binder_tons, axis=1)

# --- 4. RUN FORECAST HORIZONS (COMPONED 5-YEAR INTERVALS 2023-2053) ---
# References for Reviewers:
# [Ref E] Componed rate application: Spring 2025 FHWA Long-Term National VMT Forecasts (Baseline Outlook).
print("Projecting demand profiles across 5-year compounding forecast intervals (2023-2053)...")

forecast_years = [2028, 2033, 2038, 2043, 2048, 2053]

for year in forecast_years:
    time_delta = year - 2023
    # Compounding Equation: Demand_t = Demand_baseline * (1 + r)^t
    df[f'Binder_Demand_{year}'] = df['Binder_Demand_2023'] * (np.power(1.0 + df['Assigned_Annual_Growth_Rate'], time_delta))

# --- 5. AGGREGATE MATRICES UP TO UNIFIED COUNTY LEVEL ---
print("Aggregating functional system rows into unified county geographic profiles...")

# Group columns to drop functional layout detail and present final data purely at the county level
aggregation_dictionary = {
    'Population': 'first',
    'Land_Area_SqMi': 'first',
    'Population_Density': 'first',
    'lane_miles': 'sum',
    'Binder_Demand_2023': 'sum'
}
for year in forecast_years:
    aggregation_dictionary[f'Binder_Demand_{year}'] = 'sum'

df_county_demand = df.groupby(['state_abbr', 'county_fips', 'County_Name'], as_index=False).agg(aggregation_dictionary)

# Clean up column names for final deliverable reporting
df_county_demand = df_county_demand.rename(columns={'lane_miles': 'Total_County_Lane_Miles'})

# Round decimal strings to reasonable structural measurements (nearest tenth of a short ton)
tonnage_cols = ['Binder_Demand_2023'] + [f'Binder_Demand_{y}' for y in forecast_years]
df_county_demand[tonnage_cols] = df_county_demand[tonnage_cols].round(1)
df_county_demand['Total_County_Lane_Miles'] = df_county_demand['Total_County_Lane_Miles'].round(2)
df_county_demand['Population_Density'] = df_county_demand['Population_Density'].round(2)

# --- 6. FILE EXPORT AND INTEGRITY VERIFICATION ---
df_county_demand.to_csv(binder_output_path, index=False)
print(f"\nModel Execution Successful! Outputs compiled and exported to:\n{binder_output_path}")

print("\nReviewer Check - Sample Output Preview (First 5 Counties, Short Tons of Binder/Year):")
preview_cols = ['state_abbr', 'County_Name', 'Total_County_Lane_Miles', 'Binder_Demand_2023', 'Binder_Demand_2033', 'Binder_Demand_2043', 'Binder_Demand_2053']
print(df_county_demand[preview_cols].head())