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
import numpy as np
import pandas as pd


# ==================================================================================================
# 1. CONFIGURATION AND FILE PATHS
# ==================================================================================================

data_dir = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data"

# Inputs
lane_miles_path = os.path.join(data_dir, "county_lane_miles_structural_summary.csv")
demographics_path = os.path.join(data_dir, "county_geo_demographics_summary.csv")

# FAF5 Truck Data Folder
faf_dir = os.path.join(data_dir, "All_Experimental_Disaggregation_Factors")  # https://www.bts.gov/faf
truck_origin_path = os.path.join(faf_dir, "truck_origin_factors.csv")
truck_destination_path = os.path.join(faf_dir, "truck_destination_factors.csv")

# Output Master Matrix Path
master_output_path = os.path.join(data_dir, "master_infrastructure_model_matrix.csv")
missing_demographics_output_path = os.path.join(data_dir, "counties_missing_demographics.csv")

# Minimum county-functional-class HPMS surface coverage required before using the county fraction
MIN_COUNTY_SURFACE_COVERAGE = 0.50

print("Initializing Master Infrastructure Data Stream Merge...")


# ==================================================================================================
# 2. FILE EXISTENCE VERIFICATION
# ==================================================================================================

required_files = [lane_miles_path, demographics_path, truck_origin_path, truck_destination_path]
missing_files = [path for path in required_files if not os.path.exists(path)]

if missing_files:
    raise FileNotFoundError("The following required input files are missing:\n" + "\n".join(missing_files))


# ==================================================================================================
# 3. LOAD AND AGGREGATE FAF5 TRUCK FREIGHT FACTORS
# ==================================================================================================

print("Loading and consolidating FAF5 regional truck disaggregation factors...")

df_orig = pd.read_csv(truck_origin_path)
df_dest = pd.read_csv(truck_destination_path)

# Ensure county FIPS codes are consistently zero-padded five-digit strings
df_orig["county_fips"] = df_orig["dms_orig_cnty"].astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(5)
df_dest["county_fips"] = df_dest["dms_dest_cnty"].astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(5)

# Group by county FIPS to determine the total regional freight footprint across all commodity groups
df_orig_agg = df_orig.groupby("county_fips", as_index=False)["f_orig"].sum()
df_dest_agg = df_dest.groupby("county_fips", as_index=False)["f_dest"].sum()

# Merge origin and destination footprints into a single structural freight intensity index
df_faf = pd.merge(df_orig_agg, df_dest_agg, on="county_fips", how="outer").fillna(0)
df_faf["Structural_Freight_Intensity"] = df_faf["f_orig"] + df_faf["f_dest"]
df_faf_clean = df_faf[["county_fips", "Structural_Freight_Intensity"]].copy()


# ==================================================================================================
# 4. LOAD CORE INFRASTRUCTURE AND DEMOGRAPHIC DATA
# ==================================================================================================

print("Loading HPMS lane-mile structural summary and demographic data...")

df_lanes = pd.read_csv(lane_miles_path, dtype={"county_fips": str})
df_demo = pd.read_csv(demographics_path, dtype={"county_fips": str})

# Clean FIPS columns for a consistent merge key
df_lanes["county_fips"] = df_lanes["county_fips"].astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(5)
df_demo["county_fips"] = df_demo["county_fips"].astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(5)

# Check the fields required to calculate county, state, and national asphalt fractions
required_lane_columns = {
    "state_abbr",
    "county_fips",
    "functional_system",
    "lane_miles",
    "known_surface_lane_miles",
    "asphalt_surface_lane_miles",
    "asphalt_surface_fraction_reported",
    "surface_type_coverage",
}

missing_lane_columns = required_lane_columns.difference(df_lanes.columns)

if missing_lane_columns:
    raise KeyError(
        "The HPMS structural summary is missing the following required columns:\n"
        + "\n".join(sorted(missing_lane_columns))
    )

# Convert required structural fields to numeric
surface_numeric_columns = [
    "functional_system",
    "lane_miles",
    "known_surface_lane_miles",
    "asphalt_surface_lane_miles",
    "asphalt_surface_fraction_reported",
    "surface_type_coverage",
]

for column in surface_numeric_columns:
    df_lanes[column] = pd.to_numeric(df_lanes[column], errors="coerce")

# Remove unusable rows before calculating surface fractions
invalid_lane_rows = (
    df_lanes["county_fips"].isna()
    | df_lanes["state_abbr"].isna()
    | df_lanes["functional_system"].isna()
    | df_lanes["lane_miles"].isna()
    | (df_lanes["lane_miles"] <= 0)
)

if invalid_lane_rows.any():
    print(f"Removing {invalid_lane_rows.sum():,} HPMS rows with missing identifiers or nonpositive lane miles.")

df_lanes = df_lanes.loc[~invalid_lane_rows].copy()

# Identify valid HPMS counties that do not have a matching demographic record
lane_counties = df_lanes[["state_abbr", "county_fips"]].drop_duplicates()
demographic_counties = set(df_demo["county_fips"])

missing_demographic_counties = (
    lane_counties.loc[~lane_counties["county_fips"].isin(demographic_counties)]
    .sort_values(["state_abbr", "county_fips"])
)

print("\nCounties missing from demographic data:")

if missing_demographic_counties.empty:
    print("None")
else:
    print(missing_demographic_counties.groupby("state_abbr").size().sort_values(ascending=False))

missing_demographic_counties.to_csv(missing_demographics_output_path, index=False)
print(f"Missing-demographics audit written to:\n{missing_demographics_output_path}")


# ==================================================================================================
# 5. CALCULATE STATE-FUNCTIONAL-CLASS ASPHALT FRACTIONS
# ==================================================================================================

print("Calculating state-functional-class asphalt-surface fractions...")

state_surface = (
    df_lanes.groupby(["state_abbr", "functional_system"], as_index=False)
    .agg(
        State_Total_Lane_Miles=("lane_miles", "sum"),
        State_Known_Surface_Lane_Miles=("known_surface_lane_miles", "sum"),
        State_Asphalt_Surface_Lane_Miles=("asphalt_surface_lane_miles", "sum"),
    )
)

state_surface["State_Asphalt_Surface_Fraction"] = (
    state_surface["State_Asphalt_Surface_Lane_Miles"]
    / state_surface["State_Known_Surface_Lane_Miles"].replace(0, np.nan)
)

state_surface["State_Surface_Type_Coverage"] = (
    state_surface["State_Known_Surface_Lane_Miles"]
    / state_surface["State_Total_Lane_Miles"].replace(0, np.nan)
)

state_surface = state_surface[
    [
        "state_abbr",
        "functional_system",
        "State_Asphalt_Surface_Fraction",
        "State_Surface_Type_Coverage",
    ]
]


# ==================================================================================================
# 6. CALCULATE NATIONAL-FUNCTIONAL-CLASS ASPHALT FRACTIONS
# ==================================================================================================

print("Calculating national-functional-class asphalt-surface fractions...")

national_surface = (
    df_lanes.groupby("functional_system", as_index=False)
    .agg(
        National_Total_Lane_Miles=("lane_miles", "sum"),
        National_Known_Surface_Lane_Miles=("known_surface_lane_miles", "sum"),
        National_Asphalt_Surface_Lane_Miles=("asphalt_surface_lane_miles", "sum"),
    )
)

national_surface["National_Asphalt_Surface_Fraction"] = (
    national_surface["National_Asphalt_Surface_Lane_Miles"]
    / national_surface["National_Known_Surface_Lane_Miles"].replace(0, np.nan)
)

national_surface["National_Surface_Type_Coverage"] = (
    national_surface["National_Known_Surface_Lane_Miles"]
    / national_surface["National_Total_Lane_Miles"].replace(0, np.nan)
)

national_surface = national_surface[
    [
        "functional_system",
        "National_Asphalt_Surface_Fraction",
        "National_Surface_Type_Coverage",
    ]
]


# ==================================================================================================
# 7. MERGE STATE AND NATIONAL FALLBACK VALUES INTO THE HPMS COUNTY MATRIX
# ==================================================================================================

print("Assigning county, state, and national asphalt-surface fractions...")

df_lanes = pd.merge(
    df_lanes,
    state_surface,
    on=["state_abbr", "functional_system"],
    how="left",
    validate="many_to_one",
)

df_lanes = pd.merge(
    df_lanes,
    national_surface,
    on="functional_system",
    how="left",
    validate="many_to_one",
)

# County value is used only when county-functional-class surface coverage is at least 50%
county_fraction_available = (
    df_lanes["surface_type_coverage"].ge(MIN_COUNTY_SURFACE_COVERAGE)
    & df_lanes["asphalt_surface_fraction_reported"].notna()
)

# State value is used when the county value is not sufficiently supported
state_fraction_available = df_lanes["State_Asphalt_Surface_Fraction"].notna()

# National value is used only when the county and state values are unavailable
national_fraction_available = df_lanes["National_Asphalt_Surface_Fraction"].notna()

df_lanes["Asphalt_Surface_Fraction"] = np.nan
df_lanes["Asphalt_Surface_Basis"] = "Unresolved"
df_lanes["Asphalt_Surface_Coverage"] = np.nan

# Priority 1: County-functional-class HPMS value
df_lanes.loc[county_fraction_available, "Asphalt_Surface_Fraction"] = df_lanes.loc[
    county_fraction_available, "asphalt_surface_fraction_reported"
]

df_lanes.loc[county_fraction_available, "Asphalt_Surface_Basis"] = "County-functional class HPMS"

df_lanes.loc[county_fraction_available, "Asphalt_Surface_Coverage"] = df_lanes.loc[
    county_fraction_available, "surface_type_coverage"
]

# Priority 2: State-functional-class HPMS value
use_state_fallback = ~county_fraction_available & state_fraction_available

df_lanes.loc[use_state_fallback, "Asphalt_Surface_Fraction"] = df_lanes.loc[
    use_state_fallback, "State_Asphalt_Surface_Fraction"
]

df_lanes.loc[use_state_fallback, "Asphalt_Surface_Basis"] = "State-functional class HPMS fallback"

df_lanes.loc[use_state_fallback, "Asphalt_Surface_Coverage"] = df_lanes.loc[
    use_state_fallback, "State_Surface_Type_Coverage"
]

# Priority 3: National-functional-class HPMS value
use_national_fallback = ~county_fraction_available & ~state_fraction_available & national_fraction_available

df_lanes.loc[use_national_fallback, "Asphalt_Surface_Fraction"] = df_lanes.loc[
    use_national_fallback, "National_Asphalt_Surface_Fraction"
]

df_lanes.loc[use_national_fallback, "Asphalt_Surface_Basis"] = "National-functional class HPMS fallback"

df_lanes.loc[use_national_fallback, "Asphalt_Surface_Coverage"] = df_lanes.loc[
    use_national_fallback, "National_Surface_Type_Coverage"
]

# Ensure final columns are numeric
df_lanes["Asphalt_Surface_Fraction"] = pd.to_numeric(df_lanes["Asphalt_Surface_Fraction"], errors="coerce")
df_lanes["Asphalt_Surface_Coverage"] = pd.to_numeric(df_lanes["Asphalt_Surface_Coverage"], errors="coerce")

# Stop if any county-functional-class row remains unresolved
unresolved_surface = df_lanes["Asphalt_Surface_Fraction"].isna()

if unresolved_surface.any():
    unresolved_rows = df_lanes.loc[
        unresolved_surface,
        ["state_abbr", "county_fips", "functional_system"],
    ]

    raise ValueError(
        f"{len(unresolved_rows):,} county-functional-class rows have no asphalt-surface fraction:\n"
        + unresolved_rows.head(20).to_string(index=False)
    )

# Verify physical limits
invalid_surface_fraction = ~df_lanes["Asphalt_Surface_Fraction"].between(0.0, 1.0, inclusive="both")

if invalid_surface_fraction.any():
    invalid_rows = df_lanes.loc[
        invalid_surface_fraction,
        ["state_abbr", "county_fips", "functional_system", "Asphalt_Surface_Fraction"],
    ]

    raise ValueError(
        "Asphalt_Surface_Fraction contains values outside the physical range of 0–1:\n"
        + invalid_rows.head(20).to_string(index=False)
    )

# Report assignment results
print("\nAsphalt-surface assignment basis:")
print(df_lanes["Asphalt_Surface_Basis"].value_counts())

weighted_asphalt_fraction = (
    (df_lanes["lane_miles"] * df_lanes["Asphalt_Surface_Fraction"]).sum()
    / df_lanes["lane_miles"].sum()
)

print(f"\nFinal lane-mile-weighted asphalt-surface fraction: {weighted_asphalt_fraction:.1%}")

# Remove temporary state and national fallback calculation columns
df_lanes = df_lanes.drop(
    columns=[
        "State_Asphalt_Surface_Fraction",
        "State_Surface_Type_Coverage",
        "National_Asphalt_Surface_Fraction",
        "National_Surface_Type_Coverage",
    ]
)


# ==================================================================================================
# 8. STEPWISE SEQUENTIAL MERGING
# ==================================================================================================

print("\nMerging structural spatial data streams...")

hpms_counties_before_merge = df_lanes["county_fips"].nunique()
hpms_lane_miles_before_merge = df_lanes["lane_miles"].sum()

# First merge: attach county demographics to the county-functional-class HPMS matrix
df_master = pd.merge(
    df_lanes,
    df_demo,
    on="county_fips",
    how="left",
    validate="many_to_one",
    indicator="Demographic_Merge_Status",
)

df_master["Demographic_Data_Available"] = df_master["Demographic_Merge_Status"].eq("both")
df_master["Population_Density_Imputed"] = df_master["Population_Density"].isna()

# County_Name is used as a grouping key by the binder-demand script. Give unmatched
# counties a transparent FIPS-based label so pandas groupby does not remove them.
df_master["County_Name"] = df_master["County_Name"].fillna(
    df_master["state_abbr"] + " county FIPS " + df_master["county_fips"]
)

# Population density affects only the future growth rate for functional systems 4–7.
# The value is filled for file completeness, while the growth function below applies
# the neutral FHWA light-duty baseline whenever the demographic record is unavailable.
national_median_population_density = df_demo["Population_Density"].dropna().median()

if pd.isna(national_median_population_density):
    raise ValueError("No valid Population_Density values are available in the demographic file.")

df_master["Population_Density"] = df_master["Population_Density"].fillna(national_median_population_density)
df_master = df_master.drop(columns=["Demographic_Merge_Status"])

# Second merge: attach county-level FAF5 truck freight intensity
df_master = pd.merge(
    df_master,
    df_faf_clean,
    on="county_fips",
    how="left",
    validate="many_to_one",
)

df_master["Structural_Freight_Intensity"] = df_master["Structural_Freight_Intensity"].fillna(0.0)

print(f"Master matrix rows after demographic and FAF5 merging: {len(df_master):,}")
print(f"Unique master counties: {df_master['county_fips'].nunique():,}")
print(f"Counties using neutral demographic-growth fallback: {df_master.loc[df_master['Population_Density_Imputed'], 'county_fips'].nunique():,}")


# ==================================================================================================
# 9. APPLY FHWA SPRING 2025 LONG-TERM BASELINE VMT METRICS
# ==================================================================================================

print("Calculating localized growth rates using Spring 2025 FHWA baseline parameters...")

# Baseline annual compound growth vectors for the 2023–2053 outlook
# https://www.fhwa.dot.gov/policyinformation/tables/vmt/vmt_forecast_sum.cfm
LDV_growth_baseline = 0.005
SingleUnit_truck_growth = 0.020
Combination_truck_growth = 0.009

# Blended freight compound rate for major highways
heavy_truck_baseline = (SingleUnit_truck_growth + Combination_truck_growth) / 2.0


def calculate_localized_growth(row):
    """
    Differentiate growth rates by functional class.

    Functional systems 1–3:
        Major freight corridors scaled by FAF5 structural freight intensity.

    Functional systems 4–7:
        Collectors and local roads scaled by county population density.
    """

    sys_tier = float(row["functional_system"])

    if sys_tier in [1.0, 2.0, 3.0]:
        return heavy_truck_baseline * (1.0 + row["Structural_Freight_Intensity"])

    if row["Population_Density_Imputed"]:
        return LDV_growth_baseline

    density_mod = max(0.1, min(row["Population_Density"] / 100.0, 3.0))
    return LDV_growth_baseline * density_mod


df_master["Assigned_Annual_Growth_Rate"] = df_master.apply(calculate_localized_growth, axis=1)


# ==================================================================================================
# 10. FINAL INTEGRITY CHECKS
# ==================================================================================================

# Verify that each county-functional-class combination is unique
duplicate_keys = df_master.duplicated(
    subset=["state_abbr", "county_fips", "functional_system"],
    keep=False,
)

if duplicate_keys.any():
    duplicate_rows = df_master.loc[
        duplicate_keys,
        ["state_abbr", "county_fips", "functional_system"],
    ]

    raise ValueError(
        "Duplicate county-functional-class rows were created during merging:\n"
        + duplicate_rows.head(20).to_string(index=False)
    )

# Verify that the demographic and FAF5 merges did not remove any valid HPMS counties or lane miles
master_counties_after_merge = df_master["county_fips"].nunique()
master_lane_miles_after_merge = df_master["lane_miles"].sum()

if master_counties_after_merge != hpms_counties_before_merge:
    raise ValueError(
        f"County loss detected during merging: HPMS had {hpms_counties_before_merge:,} counties, "
        f"but the master matrix has {master_counties_after_merge:,}."
    )

if not np.isclose(master_lane_miles_after_merge, hpms_lane_miles_before_merge, rtol=0.0, atol=0.1):
    raise ValueError(
        f"Lane-mile loss detected during merging: HPMS had {hpms_lane_miles_before_merge:,.4f}, "
        f"but the master matrix has {master_lane_miles_after_merge:,.4f}."
    )

# Verify required model inputs
required_master_columns = [
    "lane_miles",
    "Asphalt_Surface_Fraction",
    "Population_Density",
    "Structural_Freight_Intensity",
    "Assigned_Annual_Growth_Rate",
]

missing_master_values = df_master[required_master_columns].isna().sum()
missing_master_values = missing_master_values[missing_master_values > 0]

if not missing_master_values.empty:
    raise ValueError(
        "Missing required model values remain in the master matrix:\n"
        + missing_master_values.to_string()
    )

# Round only the final presentation columns
df_master["lane_miles"] = df_master["lane_miles"].round(4)
df_master["Asphalt_Surface_Fraction"] = df_master["Asphalt_Surface_Fraction"].round(4)
df_master["Asphalt_Surface_Coverage"] = df_master["Asphalt_Surface_Coverage"].round(4)
df_master["Assigned_Annual_Growth_Rate"] = df_master["Assigned_Annual_Growth_Rate"].round(8)


# ==================================================================================================
# 11. EXPORT MASTER MATRIX
# ==================================================================================================

df_master.to_csv(master_output_path, index=False)

print(f"\nIntegration Complete! Master model matrix written to:\n{master_output_path}")

print("\nIntegrated Matrix Sample Output:")
preview_columns = [
    "state_abbr",
    "county_fips",
    "functional_system",
    "lane_miles",
    "Asphalt_Surface_Fraction",
    "Asphalt_Surface_Basis",
    "Asphalt_Surface_Coverage",
    "Population_Density",
    "Structural_Freight_Intensity",
    "Assigned_Annual_Growth_Rate",
]

print(df_master[preview_columns].head().to_string(index=False))

print("\nFinal asphalt-surface assignment counts:")
print(df_master["Asphalt_Surface_Basis"].value_counts())

print(
    "\nFinal master lane-mile-weighted asphalt fraction: "
    f"{(df_master['lane_miles'] * df_master['Asphalt_Surface_Fraction']).sum() / df_master['lane_miles'].sum():.1%}"
)

