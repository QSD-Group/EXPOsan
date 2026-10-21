# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to:
https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
"""

import os
import glob
from io import StringIO

import numpy as np
import pandas as pd
import pyogrio
import requests


# =============================================================================
# 1. FILE PATH CONFIGURATION
# =============================================================================

HPMS_DIR = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\HPMS2024"
OUTPUT_CSV = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\county_lane_miles_structural_summary.csv"
COUNTY_SUMMARY_CSV = os.path.join(HPMS_DIR, "hpms_county_summary_2024.csv")
CALIBRATION_AUDIT_CSV = os.path.join(HPMS_DIR, "hm60_state_functional_class_calibration_audit.csv")
EXTREME_CALIBRATION_CSV = os.path.join(HPMS_DIR, "hm60_extreme_calibration_factors.csv")
EXTRACTION_DIAGNOSTIC_CSV = os.path.join(HPMS_DIR, "hpms_state_class_extraction_diagnostics.csv")
HM60_URL = "https://www.fhwa.dot.gov/policyinformation/statistics/2024/hm60.cfm"
CALIBRATION_REVIEW_LOW = 0.50
CALIBRATION_REVIEW_HIGH = 2.00

# HPMS surface codes with an asphalt wearing surface.
ASPHALT_SURFACE_CODES = {2, 6, 7, 8}

# Concrete pavement without a current asphalt wearing surface.
CONCRETE_SURFACE_CODES = {3, 4, 5, 9, 10}

# All recognized HPMS surface-type codes.
KNOWN_SURFACE_CODES = set(range(1, 12))

STATE_ABBR_TO_FIPS = {
    "AL": "01", "AK": "02", "AZ": "04", "AR": "05", "CA": "06",
    "CO": "08", "CT": "09", "DE": "10", "DC": "11", "FL": "12",
    "GA": "13", "HI": "15", "ID": "16", "IL": "17", "IN": "18",
    "IA": "19", "KS": "20", "KY": "21", "LA": "22", "ME": "23",
    "MD": "24", "MA": "25", "MI": "26", "MN": "27", "MS": "28",
    "MO": "29", "MT": "30", "NE": "31", "NV": "32", "NH": "33",
    "NJ": "34", "NM": "35", "NY": "36", "NC": "37", "ND": "38",
    "OH": "39", "OK": "40", "OR": "41", "PA": "42", "RI": "44",
    "SC": "45", "SD": "46", "TN": "47", "TX": "48", "UT": "49",
    "VT": "50", "VA": "51", "WA": "53", "WV": "54", "WI": "55",
    "WY": "56", "PR": "72",
}

STATE_NAME_TO_ABBR = {
    "Alabama": "AL", "Alaska": "AK", "Arizona": "AZ", "Arkansas": "AR", "California": "CA",
    "Colorado": "CO", "Connecticut": "CT", "Delaware": "DE", "District of Columbia": "DC", "Florida": "FL",
    "Georgia": "GA", "Hawaii": "HI", "Idaho": "ID", "Illinois": "IL", "Indiana": "IN", "Iowa": "IA",
    "Kansas": "KS", "Kentucky": "KY", "Louisiana": "LA", "Maine": "ME", "Maryland": "MD",
    "Massachusetts": "MA", "Michigan": "MI", "Minnesota": "MN", "Mississippi": "MS", "Missouri": "MO",
    "Montana": "MT", "Nebraska": "NE", "Nevada": "NV", "New Hampshire": "NH", "New Jersey": "NJ",
    "New Mexico": "NM", "New York": "NY", "North Carolina": "NC", "North Dakota": "ND", "Ohio": "OH",
    "Oklahoma": "OK", "Oregon": "OR", "Pennsylvania": "PA", "Rhode Island": "RI", "South Carolina": "SC",
    "South Dakota": "SD", "Tennessee": "TN", "Texas": "TX", "Utah": "UT", "Vermont": "VT",
    "Virginia": "VA", "Washington": "WA", "West Virginia": "WV", "Wisconsin": "WI", "Wyoming": "WY",
    "Puerto Rico": "PR",
}

STATE_FIPS_TO_ABBR = {state_fips: state_abbr for state_abbr, state_fips in STATE_ABBR_TO_FIPS.items()}


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

def safe_divide(numerator, denominator):
    """Divide two Series and return NaN where the denominator is zero."""
    return numerator.div(denominator.replace(0, np.nan))


def get_available_columns(gdb_path, layer_name):
    """Return a lowercase-to-original field-name mapping for one layer."""
    layer_info = pyogrio.read_info(gdb_path, layer=layer_name)
    return {str(field).lower(): str(field) for field in layer_info["fields"]}


def build_county_fips(county_value, state_fips):
    """Build a five-digit county FIPS from either a county code or full FIPS."""
    if pd.isna(county_value):
        return None

    county_text = str(county_value).strip()
    if county_text.endswith(".0"):
        county_text = county_text[:-2]
    county_text = "".join(character for character in county_text if character.isdigit())

    if not county_text:
        return None
    if len(county_text) == 5:
        return county_text if county_text[:2] == state_fips and county_text[-3:] != "000" else None
    if len(county_text) == 4:
        padded_full_fips = county_text.zfill(5)
        if padded_full_fips[:2] == state_fips and padded_full_fips[-3:] != "000":
            return padded_full_fips
    if len(county_text) <= 3:
        county_code = county_text.zfill(3)
        return state_fips + county_code if county_code != "000" else None
    return None


def clean_hm60_number(value):
    """Convert an HM-60 HTML cell to a numeric lane-mile value."""
    text = str(value).replace(",", "").replace("\u00a0", " ").strip()
    text = text.replace("—", "").replace("–", "").replace("-", "")
    return pd.to_numeric(text, errors="coerce")


def load_hm60_state_functional_totals():
    """Download FHWA HM-60 and return rural-plus-urban lane miles by state and functional class."""
    print(f"Downloading official FHWA HM-60 control totals:\n{HM60_URL}")
    response = requests.get(HM60_URL, timeout=60)
    response.raise_for_status()
    tables = pd.read_html(StringIO(response.text))

    hm60_wide = None
    for table in tables:
        flattened = table.copy()
        flattened.columns = [" ".join(str(part) for part in column if str(part) != "nan").strip() if isinstance(column, tuple) else str(column).strip() for column in flattened.columns]
        first_column = flattened.iloc[:, 0].astype(str).str.strip()
        if first_column.isin(STATE_NAME_TO_ABBR).sum() >= 40 and flattened.shape[1] >= 17:
            hm60_wide = flattened
            break

    if hm60_wide is None:
        raise ValueError("Could not identify the FHWA HM-60 state lane-mile table. Check whether the webpage layout changed.")

    first_column = hm60_wide.columns[0]
    hm60_wide[first_column] = hm60_wide[first_column].astype(str).str.strip().str.replace(r"\s*\d+$", "", regex=True)
    hm60_wide = hm60_wide[hm60_wide[first_column].isin(STATE_NAME_TO_ABBR)].copy()
    hm60_wide["state_abbr"] = hm60_wide[first_column].map(STATE_NAME_TO_ABBR)

    # Expected HM-60 layout after the state column:
    # rural classes 1-7, rural total, urban classes 1-7, urban total, overall total.
    if hm60_wide.shape[1] < 18:
        raise ValueError(f"HM-60 table has only {hm60_wide.shape[1]} columns; at least 18 were expected.")

    records = []
    for _, row in hm60_wide.iterrows():
        for functional_system in range(1, 8):
            rural_value = clean_hm60_number(row.iloc[functional_system])
            urban_value = clean_hm60_number(row.iloc[functional_system + 8])
            official_lane_miles = np.nansum([rural_value, urban_value])
            records.append({"state_abbr": row["state_abbr"], "f_system": float(functional_system), "HM60_Official_Lane_Miles": official_lane_miles})

    hm60_long = pd.DataFrame(records)
    if hm60_long.duplicated(["state_abbr", "f_system"]).any():
        raise ValueError("Duplicate state-functional-class rows were produced while parsing HM-60.")
    return hm60_long


def load_hpms_county_summary():
    """Read official HPMS County Summary mileage for functional systems 6 and 7."""
    if not os.path.exists(COUNTY_SUMMARY_CSV):
        raise FileNotFoundError(
            "The official HPMS County Summary CSV is required for functional systems 6 and 7.\n"
            f"Expected file:\n{COUNTY_SUMMARY_CSV}\n\n"
            "Export/download the 2024 HPMS County Summary data with these fields:\n"
            "State_Code, County_Code, F_System, RMC_L_System_Length"
        )

    county_summary = pd.read_csv(COUNTY_SUMMARY_CSV)
    normalized_names = {str(column).strip().lower(): column for column in county_summary.columns}
    field_aliases = {
        "state_code": ["state_code", "stateid", "state_id"],
        "county_code": ["county_code", "countyid", "county_id"],
        "f_system": ["f_system", "fsystem", "functional_system"],
        "rmc_l_system_length": ["rmc_l_system_length", "systemlength", "system_length"],
    }

    rename_columns = {}
    missing = []
    for standard_name, aliases in field_aliases.items():
        matched_name = next((normalized_names[alias] for alias in aliases if alias in normalized_names), None)
        if matched_name is None:
            missing.append(f"{standard_name} (accepted aliases: {aliases})")
        else:
            rename_columns[matched_name] = standard_name

    if missing:
        raise KeyError(
            "HPMS County Summary is missing required fields:\n"
            + "\n".join(missing)
            + f"\n\nColumns found: {list(county_summary.columns)}"
        )

    county_summary = county_summary.rename(columns=rename_columns)
    required = list(field_aliases)
    for column in required:
        county_summary[column] = pd.to_numeric(
            county_summary[column].astype(str).str.replace(",", "", regex=False).str.strip(),
            errors="coerce",
        )

    county_summary = county_summary.loc[
        county_summary["state_code"].notna()
        & county_summary["county_code"].notna()
        & county_summary["f_system"].isin([6, 7])
        & county_summary["rmc_l_system_length"].notna()
        & (county_summary["rmc_l_system_length"] >= 0)
    ].copy()

    county_summary["state_fips"] = county_summary["state_code"].astype(int).astype(str).str.zfill(2)
    county_summary["state_abbr"] = county_summary["state_fips"].map(STATE_FIPS_TO_ABBR)
    county_summary["county_id"] = county_summary["state_fips"] + county_summary["county_code"].astype(int).astype(str).str.zfill(3)
    county_summary["f_system"] = county_summary["f_system"].astype(float)

    if county_summary["state_abbr"].isna().any():
        bad_codes = sorted(county_summary.loc[county_summary["state_abbr"].isna(), "state_fips"].unique())
        raise ValueError(f"Unrecognized state FIPS codes in HPMS County Summary: {bad_codes}")

    # FHWA calculates lane miles for these summary systems as road miles × 2.
    county_summary["lane_miles"] = county_summary["rmc_l_system_length"] * 2.0
    county_summary = county_summary.groupby(["state_abbr", "county_id", "f_system"], as_index=False).agg(
        lane_miles=("lane_miles", "sum"),
        total_section_miles=("rmc_l_system_length", "sum"),
    )
    return county_summary


# =============================================================================
# 3. MAIN HPMS EXTRACTION PIPELINE
# =============================================================================

def main():
    print("--- Initializing HPMS County Structural Pipeline ---")

    # =========================================================================
    # 1. LOCATE HPMS GEODATABASE AND STATE LAYERS
    # =========================================================================

    gdb_paths = glob.glob(os.path.join(HPMS_DIR, "*.gdb"))

    if not gdb_paths:
        raise FileNotFoundError(f"Could not locate a valid .gdb folder inside:\n{HPMS_DIR}")

    gdb_path = gdb_paths[0]
    all_layers = [layer[0] for layer in pyogrio.list_layers(gdb_path)]

    print(f"Connected to: {gdb_path}")
    print(f"Processing {len(all_layers)} HPMS state layers...")

    aggregated_results = []
    extraction_diagnostics = []

    required_columns = ["through_lanes", "sectionlength", "f_system", "county_id"]
    optional_columns = ["surface_type", "lane_width", "aadt", "aadt_combination", "last_overlay_thickness", "urban_code"]

    # =========================================================================
    # 2. PROCESS EACH STATE LAYER
    # =========================================================================

    for layer_name in all_layers:
        state_parts = layer_name.split("_")
        state_abbr = state_parts[2] if len(state_parts) >= 3 else "UNKNOWN"

        print(f"\nProcessing {state_abbr}: {layer_name}")

        try:
            # -----------------------------------------------------------------
            # Identify available fields
            # -----------------------------------------------------------------

            layer_info = pyogrio.read_info(gdb_path, layer=layer_name)
            available_columns = {str(field).lower(): str(field) for field in layer_info["fields"]}

            missing_required = [column for column in required_columns if column not in available_columns]

            if missing_required:
                print(f"  [Warning] Skipping layer because required fields are missing: {missing_required}")
                continue

            requested_columns = required_columns + [column for column in optional_columns if column in available_columns]
            original_column_names = [available_columns[column] for column in requested_columns]

            # -----------------------------------------------------------------
            # Read tabular fields without geometry
            # -----------------------------------------------------------------

            df = pyogrio.read_dataframe(
                gdb_path,
                layer=layer_name,
                columns=original_column_names,
                read_geometry=False,
            )

            # Normalize actual geodatabase field names to lowercase.
            df = df.rename(columns={available_columns[column]: column for column in requested_columns})

            # Create blank optional fields when they are absent from a layer.
            for column in optional_columns:
                if column not in df.columns:
                    df[column] = np.nan

            if df.empty:
                print("  [Warning] Empty layer.")
                continue

            # -----------------------------------------------------------------
            # Convert fields to numeric
            # -----------------------------------------------------------------

            numeric_columns = [
                "through_lanes",
                "sectionlength",
                "f_system",
                "county_id",
                "surface_type",
                "lane_width",
                "aadt",
                "aadt_combination",
                "last_overlay_thickness",
                "urban_code",
            ]

            for column in numeric_columns:
                df[column] = pd.to_numeric(df[column], errors="coerce")

            # -----------------------------------------------------------------
            # Standardize county FIPS before filtering and diagnose every
            # exclusion reason separately.
            # -----------------------------------------------------------------

            if state_abbr not in STATE_ABBR_TO_FIPS:
                raise KeyError(f"No state FIPS mapping found for {state_abbr}.")

            state_fips = STATE_ABBR_TO_FIPS[state_abbr]
            df["county_fips"] = df["county_id"].apply(lambda value: build_county_fips(value, state_fips))

            invalid_county = df["county_fips"].isna()
            invalid_functional_class = ~df["f_system"].isin(range(1, 8))
            invalid_section_length = df["sectionlength"].isna() | (df["sectionlength"] <= 0)
            invalid_through_lanes = df["through_lanes"].isna() | (df["through_lanes"] <= 0)

            print(
                f"  Raw segments: {len(df):,} | "
                f"Invalid county: {invalid_county.sum():,} | "
                f"Invalid functional class: {invalid_functional_class.sum():,} | "
                f"Invalid length: {invalid_section_length.sum():,} | "
                f"Missing/zero through lanes: {invalid_through_lanes.sum():,}"
            )

            # -----------------------------------------------------------------
            # Retain segments with missing through lanes using a two-lane
            # fallback. Track the affected lane mileage for quality control.
            # -----------------------------------------------------------------

            df["through_lanes_imputed"] = invalid_through_lanes
            df.loc[invalid_through_lanes, "through_lanes"] = 2.0

            # Exclude only records that cannot be assigned to a county/class or
            # have no valid section length.
            # State-by-class diagnostics are calculated before exclusions so an
            # extreme calibration factor can be traced to its source.
            diagnostic_source = df.copy()
            diagnostic_source["valid_county"] = ~invalid_county
            diagnostic_source["valid_length"] = ~invalid_section_length
            diagnostic_source["valid_lanes_original"] = ~invalid_through_lanes
            diagnostic_source["road_miles_before_filter"] = diagnostic_source["sectionlength"].where(~invalid_section_length, 0.0)
            diagnostic_source["lane_miles_before_filter"] = diagnostic_source["road_miles_before_filter"] * diagnostic_source["through_lanes"].where(~invalid_through_lanes, 2.0)

            class_diagnostic = diagnostic_source.loc[diagnostic_source["f_system"].isin(range(1, 8))].groupby("f_system", as_index=False).agg(
                Raw_Segment_Count=("sectionlength", "size"),
                Valid_County_Segment_Count=("valid_county", "sum"),
                Valid_Length_Segment_Count=("valid_length", "sum"),
                Original_Valid_Lane_Count_Segments=("valid_lanes_original", "sum"),
                HPMS_Road_Miles_Before_Filter=("road_miles_before_filter", "sum"),
                HPMS_Lane_Miles_Before_Filter=("lane_miles_before_filter", "sum"),
            )
            class_diagnostic["state_abbr"] = state_abbr
            extraction_diagnostics.append(class_diagnostic)

            invalid_segments = invalid_county | invalid_functional_class | invalid_section_length

            if invalid_segments.any():
                print(f"  Excluding {invalid_segments.sum():,} segments with missing county, class, or usable length.")

            df = df.loc[~invalid_segments].copy()

            if df.empty:
                print("  [Warning] No valid segments remained after cleaning.")
                continue

            # -----------------------------------------------------------------
            # Standardize county FIPS and functional system
            # -----------------------------------------------------------------

            df["county_id"] = df["county_fips"]
            df["f_system"] = df["f_system"].astype(float)

            # -----------------------------------------------------------------
            # Calculate segment lane miles
            # -----------------------------------------------------------------

            df["lane_miles"] = df["sectionlength"] * df["through_lanes"]
            df["section_observation_lane_miles"] = df["lane_miles"]
            df["section_observation_section_miles"] = df["sectionlength"]

            # In HPMS summary coding, 99999 represents rural. Class-6 section
            # mileage is retained as an official-mileage component only when
            # it is truly urban; all class-6/7 sections remain available as
            # observations for surface and other structural attributes.
            is_urban_minor_collector = (
                df["f_system"].eq(6.0)
                & df["urban_code"].notna()
                & df["urban_code"].ne(99999)
            )
            df["urban_class6_lane_miles"] = np.where(is_urban_minor_collector, df["lane_miles"], 0.0)
            df["urban_class6_section_miles"] = np.where(is_urban_minor_collector, df["sectionlength"], 0.0)

            df["imputed_lane_count_lane_miles"] = np.where(
                df["through_lanes_imputed"],
                df["lane_miles"],
                0.0,
            )
            df["urban_class6_imputed_lane_miles"] = np.where(
                is_urban_minor_collector & df["through_lanes_imputed"],
                df["lane_miles"],
                0.0,
            )

            # -----------------------------------------------------------------
            # Surface-type calculations
            #
            # HPMS codes:
            # 1     = Unpaved
            # 2     = Bituminous/asphalt
            # 3     = JPCP
            # 4     = JRCP
            # 5     = CRCP
            # 6     = AC overlay over AC
            # 7     = AC overlay over jointed concrete
            # 8     = AC overlay over CRCP
            # 9     = Unbonded PCC overlay
            # 10    = Bonded PCC overlay
            # 11    = Other
            #
            # Asphalt wearing surface: 2, 6, 7, 8
            # Concrete wearing surface: 3, 4, 5, 9, 10
            # -----------------------------------------------------------------

            df["surface_known"] = df["surface_type"].isin(KNOWN_SURFACE_CODES)
            df["surface_asphalt"] = df["surface_type"].isin(ASPHALT_SURFACE_CODES)
            df["surface_concrete"] = df["surface_type"].isin(CONCRETE_SURFACE_CODES)
            df["surface_unpaved"] = df["surface_type"].eq(1)
            df["surface_other"] = df["surface_type"].eq(11)

            df["known_surface_lane_miles"] = np.where(df["surface_known"], df["lane_miles"], 0.0)
            df["asphalt_surface_lane_miles"] = np.where(df["surface_asphalt"], df["lane_miles"], 0.0)
            df["concrete_surface_lane_miles"] = np.where(df["surface_concrete"], df["lane_miles"], 0.0)
            df["unpaved_surface_lane_miles"] = np.where(df["surface_unpaved"], df["lane_miles"], 0.0)
            df["other_surface_lane_miles"] = np.where(df["surface_other"], df["lane_miles"], 0.0)

            # -----------------------------------------------------------------
            # Lane-width calculations
            # -----------------------------------------------------------------

            valid_lane_width = df["lane_width"].between(6.0, 16.0, inclusive="both")

            df["valid_width_lane_miles"] = np.where(valid_lane_width, df["lane_miles"], 0.0)
            df["width_x_lane_miles"] = np.where(valid_lane_width, df["lane_width"] * df["lane_miles"], 0.0)

            # -----------------------------------------------------------------
            # AADT calculations
            #
            # AADT is weighted by section length rather than lane miles because
            # AADT represents roadway-section traffic, not traffic per lane.
            # -----------------------------------------------------------------

            valid_aadt = df["aadt"].notna() & (df["aadt"] >= 0)

            df["valid_aadt_section_miles"] = np.where(valid_aadt, df["sectionlength"], 0.0)
            df["aadt_x_section_miles"] = np.where(valid_aadt, df["aadt"] * df["sectionlength"], 0.0)

            # -----------------------------------------------------------------
            # Combination-truck AADT calculations
            # -----------------------------------------------------------------

            valid_truck_aadt = df["aadt_combination"].notna() & (df["aadt_combination"] >= 0)

            df["valid_truck_section_miles"] = np.where(valid_truck_aadt, df["sectionlength"], 0.0)
            df["truck_aadt_x_section_miles"] = np.where(
                valid_truck_aadt,
                df["aadt_combination"] * df["sectionlength"],
                0.0,
            )

            # -----------------------------------------------------------------
            # Last-overlay-thickness calculations
            #
            # A broad 0.1–12 inch validity range is retained for diagnostic
            # purposes. Do not insert this directly into the demand model until
            # units and coverage are verified.
            # -----------------------------------------------------------------

            valid_overlay = df["last_overlay_thickness"].between(0.1, 12.0, inclusive="both")

            df["valid_overlay_lane_miles"] = np.where(valid_overlay, df["lane_miles"], 0.0)
            df["overlay_x_lane_miles"] = np.where(
                valid_overlay,
                df["last_overlay_thickness"] * df["lane_miles"],
                0.0,
            )

            # -----------------------------------------------------------------
            # Aggregate within state to county × functional system
            # -----------------------------------------------------------------

            summary = df.groupby(["county_id", "f_system"], as_index=False).agg(
                lane_miles=("lane_miles", "sum"),
                total_section_miles=("sectionlength", "sum"),
                section_observation_lane_miles=("section_observation_lane_miles", "sum"),
                section_observation_section_miles=("section_observation_section_miles", "sum"),
                imputed_lane_count_lane_miles=("imputed_lane_count_lane_miles", "sum"),
                urban_class6_lane_miles=("urban_class6_lane_miles", "sum"),
                urban_class6_section_miles=("urban_class6_section_miles", "sum"),
                urban_class6_imputed_lane_miles=("urban_class6_imputed_lane_miles", "sum"),
                known_surface_lane_miles=("known_surface_lane_miles", "sum"),
                asphalt_surface_lane_miles=("asphalt_surface_lane_miles", "sum"),
                concrete_surface_lane_miles=("concrete_surface_lane_miles", "sum"),
                unpaved_surface_lane_miles=("unpaved_surface_lane_miles", "sum"),
                other_surface_lane_miles=("other_surface_lane_miles", "sum"),
                valid_width_lane_miles=("valid_width_lane_miles", "sum"),
                width_x_lane_miles=("width_x_lane_miles", "sum"),
                valid_aadt_section_miles=("valid_aadt_section_miles", "sum"),
                aadt_x_section_miles=("aadt_x_section_miles", "sum"),
                valid_truck_section_miles=("valid_truck_section_miles", "sum"),
                truck_aadt_x_section_miles=("truck_aadt_x_section_miles", "sum"),
                valid_overlay_lane_miles=("valid_overlay_lane_miles", "sum"),
                overlay_x_lane_miles=("overlay_x_lane_miles", "sum"),
            )

            summary["state_abbr"] = state_abbr
            aggregated_results.append(summary)

            state_lane_miles = summary["lane_miles"].sum()
            state_imputed_lane_miles = summary["imputed_lane_count_lane_miles"].sum()
            state_imputed_share = state_imputed_lane_miles / state_lane_miles if state_lane_miles > 0 else np.nan

            print(
                f"  Retained segments: {len(df):,} | "
                f"Lane miles: {state_lane_miles:,.1f} | "
                f"Lane miles using two-lane fallback: {state_imputed_share:.1%}"
            )

        except Exception as error:
            print(f"  [Warning] Skipped {layer_name}: {error}")
            continue

    # =========================================================================
    # 3. COMPILE NATIONAL RESULTS
    # =========================================================================

    if not aggregated_results:
        print("\nPipeline failed. No state-level results were extracted.")
        return

    print("\n--- Compiling National County-Level Structural Summary ---")

    final_df = pd.concat(aggregated_results, ignore_index=True)

    sum_columns = [
        "lane_miles",
        "total_section_miles",
        "section_observation_lane_miles",
        "section_observation_section_miles",
        "imputed_lane_count_lane_miles",
        "urban_class6_lane_miles",
        "urban_class6_section_miles",
        "urban_class6_imputed_lane_miles",
        "known_surface_lane_miles",
        "asphalt_surface_lane_miles",
        "concrete_surface_lane_miles",
        "unpaved_surface_lane_miles",
        "other_surface_lane_miles",
        "valid_width_lane_miles",
        "width_x_lane_miles",
        "valid_aadt_section_miles",
        "aadt_x_section_miles",
        "valid_truck_section_miles",
        "truck_aadt_x_section_miles",
        "valid_overlay_lane_miles",
        "overlay_x_lane_miles",
    ]

    final_df = final_df.groupby(
        ["state_abbr", "county_id", "f_system"],
        as_index=False,
    )[sum_columns].sum()

    # Preserve the full section-level aggregation as the independent source of
    # surface, width, traffic, truck-traffic, and overlay observations.
    section_observations = final_df.copy()

    # Construct official lower-system mileage independently. County Summary
    # supplies rural class 6 and all class 7 road miles using FHWA's two-lane
    # convention. True urban class-6 section mileage is then added.
    hpms_county_summary = load_hpms_county_summary()
    hpms_county_summary["imputed_lane_count_lane_miles"] = hpms_county_summary["lane_miles"]

    keys = ["state_abbr", "county_id", "f_system"]
    urban_class6_mileage = section_observations.loc[
        section_observations["f_system"].eq(6.0),
        keys + ["urban_class6_lane_miles", "urban_class6_section_miles", "urban_class6_imputed_lane_miles"],
    ].rename(
        columns={
            "urban_class6_lane_miles": "lane_miles",
            "urban_class6_section_miles": "total_section_miles",
            "urban_class6_imputed_lane_miles": "imputed_lane_count_lane_miles",
        }
    )

    lower_mileage = pd.concat(
        [
            hpms_county_summary[keys + ["lane_miles", "total_section_miles", "imputed_lane_count_lane_miles"]],
            urban_class6_mileage,
        ],
        ignore_index=True,
    ).groupby(keys, as_index=False)[["lane_miles", "total_section_miles", "imputed_lane_count_lane_miles"]].sum()

    internal_mileage_columns = {
        "lane_miles",
        "total_section_miles",
        "imputed_lane_count_lane_miles",
        "urban_class6_lane_miles",
        "urban_class6_section_miles",
        "urban_class6_imputed_lane_miles",
    }
    attribute_columns = [column for column in sum_columns if column not in internal_mileage_columns]
    lower_attributes = section_observations.loc[
        section_observations["f_system"].isin([6.0, 7.0]),
        keys + attribute_columns,
    ]
    lower_systems = lower_mileage.merge(lower_attributes, on=keys, how="left", validate="one_to_one")
    lower_systems[attribute_columns] = lower_systems[attribute_columns].fillna(0.0)

    upper_systems = section_observations.loc[section_observations["f_system"].isin([1.0, 2.0, 3.0, 4.0, 5.0])].copy()
    final_df = pd.concat([upper_systems, lower_systems], ignore_index=True, sort=False)

    print(
        "\nApplied official County Summary/HM-60 mileage to functional systems 6 and 7 while retaining "
        f"section-level structural observations: {len(hpms_county_summary):,} county-summary rows."
    )

    # =========================================================================
    # 4. CALIBRATE HPMS COUNTY SHARES TO OFFICIAL HM-60 STATE × CLASS TOTALS
    # =========================================================================

    hm60_lookup = load_hm60_state_functional_totals()

    state_hpms_totals = final_df.groupby(["state_abbr", "f_system"], as_index=False).agg(
        State_HPMS_Raw_Lane_Miles=("lane_miles", "sum")
    )

    calibration_lookup = state_hpms_totals.merge(
        hm60_lookup,
        on=["state_abbr", "f_system"],
        how="left",
        validate="one_to_one",
    )

    calibration_lookup["Lane_Mile_Calibration_Factor"] = safe_divide(
        calibration_lookup["HM60_Official_Lane_Miles"],
        calibration_lookup["State_HPMS_Raw_Lane_Miles"],
    )

    calibration_lookup["Lane_Mile_Calibration_Basis"] = np.where(
        calibration_lookup["Lane_Mile_Calibration_Factor"].notna(),
        "FHWA HM-60 state-functional class control total",
        "Uncalibrated: HM-60 total unavailable",
    )

    calibration_lookup.loc[
        calibration_lookup["f_system"].isin([6, 7])
        & calibration_lookup["Lane_Mile_Calibration_Factor"].notna(),
        "Lane_Mile_Calibration_Basis",
    ] = "HPMS County Summary × 2 lanes plus applicable urban sections; calibrated to FHWA HM-60"

    calibration_lookup["Needs_Calibration_Review"] = (
        (calibration_lookup["Lane_Mile_Calibration_Factor"] < CALIBRATION_REVIEW_LOW)
        | (calibration_lookup["Lane_Mile_Calibration_Factor"] > CALIBRATION_REVIEW_HIGH)
        | calibration_lookup["Lane_Mile_Calibration_Factor"].isna()
    )

    final_df = final_df.merge(
        calibration_lookup,
        on=["state_abbr", "f_system"],
        how="left",
        validate="many_to_one",
    )

    # Missing official values remain uncalibrated and are clearly labeled.
    final_df["Lane_Mile_Calibration_Factor"] = final_df["Lane_Mile_Calibration_Factor"].fillna(1.0)

    # HM-60 calibrates only the official mileage stream. Section-observation
    # quantities remain unchanged so pavement composition and weighted
    # structural attributes retain their original HPMS evidence.
    for column in ["lane_miles", "total_section_miles", "imputed_lane_count_lane_miles"]:
        final_df[column] = final_df[column] * final_df["Lane_Mile_Calibration_Factor"]

    # Export extraction diagnostics with HM-60 values for direct investigation.
    extraction_audit = pd.concat(extraction_diagnostics, ignore_index=True).rename(columns={"f_system": "functional_system"})
    extraction_audit = extraction_audit.merge(
        calibration_lookup.rename(columns={"f_system": "functional_system"}),
        on=["state_abbr", "functional_system"],
        how="left",
        validate="one_to_one",
    )
    extraction_audit["County_Filter_Retention_Fraction"] = safe_divide(
        extraction_audit["State_HPMS_Raw_Lane_Miles"],
        extraction_audit["HPMS_Lane_Miles_Before_Filter"],
    )
    extraction_audit.to_csv(EXTRACTION_DIAGNOSTIC_CSV, index=False)

    calibration_audit = calibration_lookup.rename(columns={"f_system": "functional_system"}).sort_values(
        "Lane_Mile_Calibration_Factor",
        ascending=False,
    )
    calibration_audit.to_csv(CALIBRATION_AUDIT_CSV, index=False)

    extreme_calibration = calibration_audit.loc[calibration_audit["Needs_Calibration_Review"]].copy()
    extreme_calibration.to_csv(EXTREME_CALIBRATION_CSV, index=False)

    print("\nLargest unique state × functional-class calibration factors:")
    print(calibration_audit.head(20).to_string(index=False))
    print("\nCalibration factors outside the 0.5–2.0 review range:")
    print(extreme_calibration.to_string(index=False) if not extreme_calibration.empty else "None")

    # Reconciliation: each calibrated state/class total should equal HM-60.
    reconciliation = final_df.groupby(["state_abbr", "f_system"], as_index=False).agg(
        Calibrated_Lane_Miles=("lane_miles", "sum"),
        HM60_Official_Lane_Miles=("HM60_Official_Lane_Miles", "first"),
    )
    reconciliation["Absolute_Difference"] = (
        reconciliation["Calibrated_Lane_Miles"] - reconciliation["HM60_Official_Lane_Miles"]
    ).abs()
    failed_reconciliation = reconciliation.loc[
        reconciliation["HM60_Official_Lane_Miles"].notna() & (reconciliation["Absolute_Difference"] > 0.1)
    ]
    if not failed_reconciliation.empty:
        raise ValueError("Some calibrated state-functional-class totals do not reconcile to HM-60 within 0.1 lane mile.")

    # =========================================================================
    # 5. DERIVED STRUCTURAL FACTORS AND COVERAGE
    # =========================================================================

    final_df["through_lanes_imputation_fraction"] = safe_divide(
        final_df["imputed_lane_count_lane_miles"],
        final_df["lane_miles"],
    ).fillna(0.0)

    final_df["asphalt_surface_fraction_reported"] = safe_divide(
        final_df["asphalt_surface_lane_miles"],
        final_df["known_surface_lane_miles"],
    )

    final_df["concrete_surface_fraction_reported"] = safe_divide(
        final_df["concrete_surface_lane_miles"],
        final_df["known_surface_lane_miles"],
    )

    final_df["unpaved_surface_fraction_reported"] = safe_divide(
        final_df["unpaved_surface_lane_miles"],
        final_df["known_surface_lane_miles"],
    )

    final_df["other_surface_fraction_reported"] = safe_divide(
        final_df["other_surface_lane_miles"],
        final_df["known_surface_lane_miles"],
    )

    final_df["surface_type_coverage"] = safe_divide(
        final_df["known_surface_lane_miles"],
        final_df["section_observation_lane_miles"],
    ).fillna(0.0).clip(lower=0.0, upper=1.0)

    final_df["average_lane_width_ft"] = safe_divide(
        final_df["width_x_lane_miles"],
        final_df["valid_width_lane_miles"],
    )

    final_df["lane_width_coverage"] = safe_divide(
        final_df["valid_width_lane_miles"],
        final_df["section_observation_lane_miles"],
    ).fillna(0.0).clip(lower=0.0, upper=1.0)

    final_df["average_aadt"] = safe_divide(
        final_df["aadt_x_section_miles"],
        final_df["valid_aadt_section_miles"],
    )

    final_df["aadt_coverage"] = safe_divide(
        final_df["valid_aadt_section_miles"],
        final_df["section_observation_section_miles"],
    ).fillna(0.0).clip(lower=0.0, upper=1.0)

    final_df["average_combination_truck_aadt"] = safe_divide(
        final_df["truck_aadt_x_section_miles"],
        final_df["valid_truck_section_miles"],
    )

    final_df["truck_aadt_coverage"] = safe_divide(
        final_df["valid_truck_section_miles"],
        final_df["section_observation_section_miles"],
    ).fillna(0.0).clip(lower=0.0, upper=1.0)

    final_df["average_last_overlay_thickness_in"] = safe_divide(
        final_df["overlay_x_lane_miles"],
        final_df["valid_overlay_lane_miles"],
    )

    final_df["overlay_thickness_coverage"] = safe_divide(
        final_df["valid_overlay_lane_miles"],
        final_df["section_observation_lane_miles"],
    ).fillna(0.0).clip(lower=0.0, upper=1.0)

    # =========================================================================
    # 6. FINAL COLUMN NAMES AND OUTPUT LAYOUT
    # =========================================================================

    final_df = final_df.rename(
        columns={
            "county_id": "county_fips",
            "f_system": "functional_system",
        }
    )

    output_columns = [
        "state_abbr",
        "county_fips",
        "functional_system",
        "lane_miles",
        "State_HPMS_Raw_Lane_Miles",
        "HM60_Official_Lane_Miles",
        "Lane_Mile_Calibration_Factor",
        "Lane_Mile_Calibration_Basis",
        "Needs_Calibration_Review",
        "through_lanes_imputation_fraction",
        "section_observation_lane_miles",
        "section_observation_section_miles",
        "known_surface_lane_miles",
        "asphalt_surface_lane_miles",
        "concrete_surface_lane_miles",
        "unpaved_surface_lane_miles",
        "other_surface_lane_miles",
        "asphalt_surface_fraction_reported",
        "concrete_surface_fraction_reported",
        "unpaved_surface_fraction_reported",
        "other_surface_fraction_reported",
        "surface_type_coverage",
        "average_lane_width_ft",
        "lane_width_coverage",
        "average_aadt",
        "aadt_coverage",
        "average_combination_truck_aadt",
        "truck_aadt_coverage",
        "average_last_overlay_thickness_in",
        "overlay_thickness_coverage",
    ]

    final_df = final_df[output_columns]
    final_df = final_df.sort_values(["state_abbr", "county_fips", "functional_system"])

    # =========================================================================
    # 7. ROUND VALUES
    # =========================================================================

    lane_mile_columns = [
        "lane_miles",
        "section_observation_lane_miles",
        "section_observation_section_miles",
        "known_surface_lane_miles",
        "asphalt_surface_lane_miles",
        "concrete_surface_lane_miles",
        "unpaved_surface_lane_miles",
        "other_surface_lane_miles",
    ]

    fraction_columns = [
        "through_lanes_imputation_fraction",
        "asphalt_surface_fraction_reported",
        "concrete_surface_fraction_reported",
        "unpaved_surface_fraction_reported",
        "other_surface_fraction_reported",
        "surface_type_coverage",
        "lane_width_coverage",
        "aadt_coverage",
        "truck_aadt_coverage",
        "overlay_thickness_coverage",
    ]

    final_df[lane_mile_columns] = final_df[lane_mile_columns].round(4)
    final_df[fraction_columns] = final_df[fraction_columns].round(4)
    final_df["average_lane_width_ft"] = final_df["average_lane_width_ft"].round(2)
    final_df["average_aadt"] = final_df["average_aadt"].round(0)
    final_df["average_combination_truck_aadt"] = final_df["average_combination_truck_aadt"].round(0)
    final_df["average_last_overlay_thickness_in"] = final_df["average_last_overlay_thickness_in"].round(2)

    # =========================================================================
    # 8. EXPORT AND NATIONAL INTEGRITY CHECKS
    # =========================================================================

    final_df.to_csv(OUTPUT_CSV, index=False)

    total_lane_miles = final_df["lane_miles"].sum()
    total_known_surface_lane_miles = final_df["known_surface_lane_miles"].sum()
    total_asphalt_lane_miles = final_df["asphalt_surface_lane_miles"].sum()
    total_section_observation_lane_miles = final_df["section_observation_lane_miles"].sum()
    total_imputed_lane_miles = (
        final_df["lane_miles"] * final_df["through_lanes_imputation_fraction"]
    ).sum()

    national_surface_coverage = (
        total_known_surface_lane_miles / total_section_observation_lane_miles
        if total_section_observation_lane_miles > 0
        else np.nan
    )

    national_asphalt_fraction_reported = (
        total_asphalt_lane_miles / total_known_surface_lane_miles
        if total_known_surface_lane_miles > 0
        else np.nan
    )

    national_lane_imputation_fraction = (
        total_imputed_lane_miles / total_lane_miles
        if total_lane_miles > 0
        else np.nan
    )
    unique_counties = final_df[
    [
        "state_abbr",
        "county_fips",
    ]
      ].drop_duplicates()

    print(f"\nSuccess! County structural matrix exported to:\n{OUTPUT_CSV}")
    print(f"Total output rows: {len(final_df):,}")
    print(f"Unique counties: "f"{len(unique_counties):,}")
    print(f"Total lane miles: {total_lane_miles:,.1f}")
    print(f"Lane miles using two-lane fallback: {national_lane_imputation_fraction:.1%}")
    print(f"Surface-type coverage within section observations: {national_surface_coverage:.1%}")
    print(f"Asphalt fraction among known surfaces: {national_asphalt_fraction_reported:.1%}")
    print(f"Rows with surface data: {(final_df['surface_type_coverage'] > 0).sum():,}")
    print(f"Rows with lane-width data: {(final_df['lane_width_coverage'] > 0).sum():,}")
    print(f"Rows with AADT data: {(final_df['aadt_coverage'] > 0).sum():,}")
    print(f"Rows with truck-AADT data: {(final_df['truck_aadt_coverage'] > 0).sum():,}")
    print(f"Rows with overlay-thickness data: {(final_df['overlay_thickness_coverage'] > 0).sum():,}")
    print(f"\nCalibration audit:\n{CALIBRATION_AUDIT_CSV}")
    print(f"Extreme calibration factors:\n{EXTREME_CALIBRATION_CSV}")
    print(f"Extraction diagnostics:\n{EXTRACTION_DIAGNOSTIC_CSV}")
    print("\nUnique state × functional-class calibration factor summary:")
    print(calibration_audit["Lane_Mile_Calibration_Factor"].describe().to_string())

    # Compare against the approximate FHWA 2024 national scale.
    if total_lane_miles < 8_000_000 or total_lane_miles > 10_000_000:
        print(
            "\n[Warning] National lane-mile total is outside the expected "
            "approximately 8–10 million range. Review exclusions and duplicate segments."
        )

    print("\nSample output:")
    print(final_df.head().to_string(index=False))


if __name__ == "__main__":
    main()
