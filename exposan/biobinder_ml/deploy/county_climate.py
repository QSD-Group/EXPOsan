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
Build county climate factors for the asphalt-binder demand model.

Priority:
1. Unanimous county climate classification from LTPP state reports.
2. Federal climate estimate using:
   - NASA POWER daily precipitation and temperature, 1991–2020.
   - U.S. Census TIGER/Line county boundaries.
3. State lane-mile-weighted fallback if a federal request fails.
4. Neutral factor 1.0 only if all other methods fail.

FHWA/LTPP classification:
    Wet:    annual precipitation > 508 mm
    Freeze: annual freezing index > 83.3 deg C-days

Freezing index:
    FI = sum(max(0, -daily_mean_temperature_C))

"""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import re
import time

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
from shapely.geometry import Point


# =============================================================================
# 1. CONFIGURATION
# =============================================================================

DATA_DIR = Path(
    r"C:\Work\Rutgers\QSDsan\EXPOsan"
    r"\exposan\biobinder_ml\deploy\Data"
)

LTPP_DIR = DATA_DIR / "State LTPP"

MASTER_INPUT = DATA_DIR / "master_infrastructure_model_matrix.csv"

MASTER_OUTPUT = (
    DATA_DIR
    / "master_infrastructure_model_matrix_with_climate.csv"
)

# Stores federal climate calculations so reruns do not repeat API requests.
FEDERAL_CACHE = DATA_DIR / "county_federal_climate_cache.csv"

# Optional diagnostic table with one row per county.
COUNTY_CLIMATE_AUDIT = DATA_DIR / "county_climate_assignment_audit.csv"

# Census TIGER/Line county boundary file.
CENSUS_COUNTY_URL = (
    "https://www2.census.gov/geo/tiger/TIGER2023/COUNTY/"
    "tl_2023_us_county.zip"
)

# 30-year climate-normal period.
CLIMATE_START_YEAR = 1991
CLIMATE_END_YEAR = 2020

# Representative points sampled inside each county.
#
# 1 = one internal county point; fastest.
# 3 = recommended for this national macro model.
# 5 = better spatial averaging but considerably slower.
COUNTY_SAMPLE_POINTS = 3

# NASA POWER requests should be kept modest.
MAX_WORKERS = 4
REQUEST_TIMEOUT_SECONDS = 120
REQUEST_RETRIES = 4

# FHWA/LTPP classification thresholds.
WET_THRESHOLD_MM = 508.0
FREEZE_THRESHOLD_C_DAYS = 83.3

# Lifespan modifiers used by the binder-demand model.
CLIMATE_LIFESPAN_MODIFIER = {
    "Wet-Freeze": 0.83,
    "Dry-Freeze": 0.95,
    "Dry-Non-Freeze": 1.05,
    "Wet-Non-Freeze": 1.15,
}


# =============================================================================
# 2. GENERAL UTILITIES
# =============================================================================

def normalize_fips(values):
    """Convert county FIPS values to five-character strings."""
    return (
        values.astype(str)
        .str.replace(r"\.0$", "", regex=True)
        .str.strip()
        .str.zfill(5)
    )


def normalize_county_name(value):
    """Normalize county-equivalent names for matching."""
    if pd.isna(value):
        return None

    value = str(value).upper().strip()

    endings = [
        " COUNTY",
        " PARISH",
        " BOROUGH",
        " CENSUS AREA",
        " MUNICIPALITY",
        " CITY AND BOROUGH",
        " MUNICIPIO",
    ]

    for ending in endings:
        if value.endswith(ending):
            value = value[:-len(ending)]

    value = re.sub(r"[^A-Z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def normalize_climate_zone(value):
    """Convert LTPP climate labels to one consistent format."""
    if pd.isna(value):
        return np.nan

    text = str(value).strip().lower()
    text = text.replace(",", "-")
    text = text.replace("_", "-")
    text = re.sub(r"\s+", "", text)

    mapping = {
        "wet-freeze": "Wet-Freeze",
        "dry-freeze": "Dry-Freeze",
        "wet-nonfreeze": "Wet-Non-Freeze",
        "wet-no-freeze": "Wet-Non-Freeze",
        "dry-nonfreeze": "Dry-Non-Freeze",
        "dry-no-freeze": "Dry-Non-Freeze",
    }

    return mapping.get(text, np.nan)


def classify_ltpp_zone(precipitation_mm, freezing_index_c_days):
    """Apply FHWA/LTPP wet/dry and freeze/non-freeze thresholds."""
    moisture = (
        "Wet"
        if precipitation_mm > WET_THRESHOLD_MM
        else "Dry"
    )

    temperature = (
        "Freeze"
        if freezing_index_c_days > FREEZE_THRESHOLD_C_DAYS
        else "Non-Freeze"
    )

    return f"{moisture}-{temperature}"


# =============================================================================
# 3. READ THE MASTER MATRIX
# =============================================================================

def read_master_file():
    if not MASTER_INPUT.exists():
        raise FileNotFoundError(
            f"Master matrix not found:\n{MASTER_INPUT}"
        )

    master = pd.read_csv(
        MASTER_INPUT,
        dtype={"county_fips": str},
    )

    required = {
        "state_abbr",
        "county_fips",
        "County_Name",
        "lane_miles",
    }

    missing = required.difference(master.columns)

    if missing:
        raise KeyError(
            f"Master matrix is missing required columns: {sorted(missing)}"
        )

    master["county_fips"] = normalize_fips(
        master["county_fips"]
    )

    return master


# =============================================================================
# 4. EXTRACT COUNTY CLIMATE LABELS FROM ALL LTPP REPORTS
# =============================================================================

def find_ltpp_workbooks():
    patterns = [
        "State_Summary_Report_*.xlsx",
        "**/State_Summary_Report_*.xlsx",
    ]

    paths = set()

    for pattern in patterns:
        paths.update(LTPP_DIR.glob(pattern))

    return sorted(paths)


def state_fips_from_filename(path):
    """
    Extract the state FIPS appearing after the state name.

    Example:
    State_Summary_Report_Alabama_1_Default_04_02_20.xlsx
                                      ^
                                      state FIPS = 1
    """
    match = re.search(
        r"_(\d{1,2})_Default_",
        path.name,
        flags=re.IGNORECASE,
    )

    if not match:
        return None

    return match.group(1).zfill(2)


def extract_ltpp_records(path):
    """Read county and Climatic Zone from General Information."""
    state_fips = state_fips_from_filename(path)

    if state_fips is None:
        print(f"Skipping; state FIPS not found in filename: {path.name}")
        return pd.DataFrame()

    try:
        table = pd.read_excel(
            path,
            sheet_name="General Information",
            header=1,
        )
    except Exception as error:
        print(f"Could not read {path.name}: {error}")
        return pd.DataFrame()

    required = {"County", "Climatic Zone"}

    if not required.issubset(table.columns):
        print(
            f"Skipping {path.name}; General Information does not "
            f"contain {sorted(required)}"
        )
        return pd.DataFrame()

    records = table[["County", "Climatic Zone"]].copy()

    records["state_fips"] = state_fips
    records["county_normalized"] = (
        records["County"].apply(normalize_county_name)
    )
    records["Climate_Zone"] = (
        records["Climatic Zone"].apply(normalize_climate_zone)
    )
    records["source_file"] = path.name

    records = records.dropna(
        subset=["county_normalized", "Climate_Zone"]
    )

    return records


def build_ltpp_county_lookup(master):
    """
    Match LTPP county names to master county FIPS.

    A county receives an LTPP assignment only if every LTPP section
    in that county has the same normalized climate zone.
    """
    workbooks = find_ltpp_workbooks()

    if not workbooks:
        raise FileNotFoundError(
            f"No State_Summary_Report_*.xlsx files found under:\n"
            f"{LTPP_DIR}"
        )

    extracted = [
        extract_ltpp_records(path)
        for path in workbooks
    ]

    extracted = [
        frame for frame in extracted
        if not frame.empty
    ]

    if not extracted:
        raise RuntimeError(
            "No usable General Information records were extracted."
        )

    ltpp = pd.concat(extracted, ignore_index=True)

    county_crosswalk = (
        master[
            ["county_fips", "County_Name"]
        ]
        .drop_duplicates()
        .copy()
    )

    county_crosswalk["state_fips"] = (
        county_crosswalk["county_fips"].str[:2]
    )

    county_crosswalk["county_normalized"] = (
        county_crosswalk["County_Name"]
        .apply(normalize_county_name)
    )

    matched = ltpp.merge(
        county_crosswalk,
        on=["state_fips", "county_normalized"],
        how="left",
        validate="many_to_one",
    )

    unmatched = matched["county_fips"].isna().sum()

    if unmatched:
        print(
            f"Warning: {unmatched:,} LTPP section records could not "
            f"be matched to a master county."
        )

    matched = matched.dropna(subset=["county_fips"])

    summary = (
        matched.groupby("county_fips", as_index=False)
        .agg(
            LTPP_Section_Count=("Climate_Zone", "size"),
            LTPP_Zone_Count=("Climate_Zone", "nunique"),
            LTPP_Zones=(
                "Climate_Zone",
                lambda x: " | ".join(sorted(set(x))),
            ),
        )
    )

    # Use LTPP only if all sections within the county agree.
    summary["LTPP_Climate_Zone"] = np.where(
        summary["LTPP_Zone_Count"].eq(1),
        summary["LTPP_Zones"],
        np.nan,
    )

    return summary


# =============================================================================
# 5. OBTAIN FEDERAL COUNTY BOUNDARIES
# =============================================================================

def load_census_counties(master_counties):
    print("Reading Census TIGER/Line county boundaries...")

    counties = gpd.read_file(CENSUS_COUNTY_URL)

    counties["county_fips"] = (
        counties["STATEFP"].astype(str).str.zfill(2)
        + counties["COUNTYFP"].astype(str).str.zfill(3)
    )

    counties = counties[
        counties["county_fips"].isin(
            master_counties["county_fips"]
        )
    ][["county_fips", "geometry"]].copy()

    if counties.empty:
        raise RuntimeError(
            "No Census county geometries matched the master FIPS codes."
        )

    return counties


# =============================================================================
# 6. SELECT REPRESENTATIVE POINTS WITHIN EACH COUNTY
# =============================================================================

def representative_points(geometry, number_of_points=3):
    """
    Generate internal points distributed across a county.

    The county is projected into a global equal-area CRS before the
    candidate grid is constructed. Points are returned as longitude,
    latitude pairs.

    This provides an approximate county-area average without querying
    every climate-grid cell.
    """
    one_county = gpd.GeoSeries(
        [geometry],
        crs="EPSG:4326",
    ).to_crs("EPSG:6933")

    projected = one_county.iloc[0]

    # Always include a guaranteed internal representative point.
    points = [projected.representative_point()]

    if number_of_points > 1:
        min_x, min_y, max_x, max_y = projected.bounds

        # Create more candidates than needed to improve coverage.
        grid_size = max(5, int(np.ceil(np.sqrt(number_of_points * 8))))

        x_values = np.linspace(min_x, max_x, grid_size + 2)[1:-1]
        y_values = np.linspace(min_y, max_y, grid_size + 2)[1:-1]

        candidates = [
            Point(x, y)
            for x in x_values
            for y in y_values
            if projected.contains(Point(x, y))
        ]

        if candidates:
            # Choose spatially distributed candidates.
            candidate_indices = np.linspace(
                0,
                len(candidates) - 1,
                min(number_of_points - 1, len(candidates)),
            ).round().astype(int)

            points.extend(
                candidates[index]
                for index in candidate_indices
            )

    point_series = gpd.GeoSeries(
        points[:number_of_points],
        crs="EPSG:6933",
    ).to_crs("EPSG:4326")

    return [
        (float(point.x), float(point.y))
        for point in point_series
    ]


# =============================================================================
# 7. NASA POWER DAILY CLIMATE REQUEST
# =============================================================================

def request_nasa_power_daily(longitude, latitude):
    """
    Retrieve daily NASA POWER data.

    T2M:
        Daily mean air temperature, degrees C.

    PRECTOTCORR:
        Bias-corrected daily precipitation, mm/day.
    """
    endpoint = (
        "https://power.larc.nasa.gov/api/temporal/daily/point"
    )

    parameters = {
        "parameters": "T2M,PRECTOTCORR",
        "community": "AG",
        "longitude": longitude,
        "latitude": latitude,
        "start": f"{CLIMATE_START_YEAR}0101",
        "end": f"{CLIMATE_END_YEAR}1231",
        "format": "JSON",
    }

    last_error = None

    for attempt in range(1, REQUEST_RETRIES + 1):
        try:
            response = requests.get(
                endpoint,
                params=parameters,
                timeout=REQUEST_TIMEOUT_SECONDS,
            )
            response.raise_for_status()
            payload = response.json()

            parameter_data = payload[
                "properties"
            ]["parameter"]

            temperature = pd.Series(
                parameter_data["T2M"],
                name="temperature_c",
                dtype=float,
            )

            precipitation = pd.Series(
                parameter_data["PRECTOTCORR"],
                name="precipitation_mm",
                dtype=float,
            )

            climate = pd.concat(
                [temperature, precipitation],
                axis=1,
            )

            climate.index = pd.to_datetime(
                climate.index,
                format="%Y%m%d",
            )

            # POWER missing-value flag.
            climate = climate.replace(-999.0, np.nan)
            climate["year"] = climate.index.year

            annual = climate.groupby("year").agg(
                annual_precipitation_mm=(
                    "precipitation_mm",
                    "sum",
                ),
                valid_temperature_days=(
                    "temperature_c",
                    "count",
                ),
                valid_precipitation_days=(
                    "precipitation_mm",
                    "count",
                ),
            )

            annual_freezing_index = (
                climate["temperature_c"]
                .clip(upper=0)
                .abs()
                .groupby(climate["year"])
                .sum(min_count=300)
            )

            annual["annual_freezing_index_c_days"] = (
                annual_freezing_index
            )

            # Retain reasonably complete years only.
            annual = annual[
                (annual["valid_temperature_days"] >= 350)
                & (annual["valid_precipitation_days"] >= 350)
            ]

            if len(annual) < 20:
                raise ValueError(
                    f"Only {len(annual)} complete climate years found."
                )

            return {
                "annual_precipitation_mm": float(
                    annual["annual_precipitation_mm"].mean()
                ),
                "freezing_index_c_days": float(
                    annual["annual_freezing_index_c_days"].mean()
                ),
                "valid_years": int(len(annual)),
            }

        except Exception as error:
            last_error = error

            if attempt < REQUEST_RETRIES:
                time.sleep(2 ** attempt)

    raise RuntimeError(
        f"NASA POWER request failed at "
        f"({latitude:.5f}, {longitude:.5f}): {last_error}"
    )


def calculate_county_federal_climate(record):
    """Average climate metrics across representative county points."""
    county_fips = record["county_fips"]
    geometry = record["geometry"]

    try:
        points = representative_points(
            geometry,
            COUNTY_SAMPLE_POINTS,
        )

        point_results = [
            request_nasa_power_daily(lon, lat)
            for lon, lat in points
        ]

        precipitation = np.mean([
            result["annual_precipitation_mm"]
            for result in point_results
        ])

        freezing_index = np.mean([
            result["freezing_index_c_days"]
            for result in point_results
        ])

        climate_zone = classify_ltpp_zone(
            precipitation,
            freezing_index,
        )

        return {
            "county_fips": county_fips,
            "Federal_Annual_Precipitation_mm": precipitation,
            "Federal_Freezing_Index_C_Days": freezing_index,
            "Federal_Climate_Zone": climate_zone,
            "Federal_Sample_Points": len(points),
            "Federal_Valid_Years": min(
                result["valid_years"]
                for result in point_results
            ),
            "Federal_Status": "Success",
        }

    except Exception as error:
        return {
            "county_fips": county_fips,
            "Federal_Annual_Precipitation_mm": np.nan,
            "Federal_Freezing_Index_C_Days": np.nan,
            "Federal_Climate_Zone": np.nan,
            "Federal_Sample_Points": 0,
            "Federal_Valid_Years": 0,
            "Federal_Status": f"Failed: {error}",
        }


def build_federal_climate_lookup(counties_needed, census_counties):
    """
    Retrieve federal climate values only for counties that do not have
    a unanimous LTPP classification.
    """
    if FEDERAL_CACHE.exists():
        cache = pd.read_csv(
            FEDERAL_CACHE,
            dtype={"county_fips": str},
        )
        cache["county_fips"] = normalize_fips(
            cache["county_fips"]
        )
    else:
        cache = pd.DataFrame(
            columns=["county_fips"]
        )

    successful_cached = set(
        cache.loc[
            cache.get(
                "Federal_Status",
                pd.Series(index=cache.index, dtype=str),
            ).eq("Success"),
            "county_fips",
        ]
    )

    needed_fips = set(counties_needed["county_fips"])
    missing_fips = needed_fips.difference(successful_cached)

    to_process = census_counties[
        census_counties["county_fips"].isin(missing_fips)
    ]

    print(
        f"Federal climate cache contains "
        f"{len(successful_cached):,} successful counties."
    )
    print(
        f"Federal climate requests required for "
        f"{len(to_process):,} counties."
    )

    new_results = []

    if not to_process.empty:
        records = to_process.to_dict("records")

        with ThreadPoolExecutor(
            max_workers=MAX_WORKERS
        ) as executor:

            future_lookup = {
                executor.submit(
                    calculate_county_federal_climate,
                    record,
                ): record["county_fips"]
                for record in records
            }

            completed = 0

            for future in as_completed(future_lookup):
                result = future.result()
                new_results.append(result)
                completed += 1

                if completed % 25 == 0:
                    print(
                        f"Processed {completed:,} of "
                        f"{len(records):,} federal counties..."
                    )

    if new_results:
        new_frame = pd.DataFrame(new_results)

        cache = pd.concat(
            [cache, new_frame],
            ignore_index=True,
        )

        # Retain the most recent attempt for each county.
        cache = cache.drop_duplicates(
            subset="county_fips",
            keep="last",
        )

        cache.to_csv(FEDERAL_CACHE, index=False)

    return cache[
        cache["county_fips"].isin(needed_fips)
    ].copy()


# =============================================================================
# 8. ASSIGN THE FINAL CLIMATE ZONE AND FACTOR
# =============================================================================

def assign_final_climate(master, ltpp_lookup, federal_lookup):
    """
    Assign a final climate zone and lifespan factor to every county.

    Assignment priority:
        1. Unanimous LTPP county classification.
        2. Federal NASA POWER county classification.
        3. State lane-mile-weighted climate-factor fallback.
        4. Neutral national factor of 1.0.

    This version is compatible with older pandas versions because it
    does not use the include_groups argument in groupby.apply().
    """

    # -------------------------------------------------------------------------
    # 1. Aggregate the master matrix to one row per county.
    # -------------------------------------------------------------------------
    county_master = (
        master[
            [
                "state_abbr",
                "county_fips",
                "County_Name",
                "lane_miles",
            ]
        ]
        .groupby(
            [
                "state_abbr",
                "county_fips",
                "County_Name",
            ],
            as_index=False,
        )
        .agg(
            Total_County_Lane_Miles=(
                "lane_miles",
                "sum",
            )
        )
    )

    # -------------------------------------------------------------------------
    # 2. Add the LTPP county classifications.
    # -------------------------------------------------------------------------
    county = county_master.merge(
        ltpp_lookup,
        on="county_fips",
        how="left",
        validate="one_to_one",
    )

    # -------------------------------------------------------------------------
    # 3. Add the federal county climate estimates.
    # -------------------------------------------------------------------------
    county = county.merge(
        federal_lookup,
        on="county_fips",
        how="left",
        validate="one_to_one",
    )

    # -------------------------------------------------------------------------
    # 4. Select the final climate zone.
    #
    # Use LTPP when all sampled sections in the county agree.
    # Otherwise, use the federal county estimate.
    # -------------------------------------------------------------------------
    county["Final_Climate_Zone"] = (
        county["LTPP_Climate_Zone"]
        .fillna(county["Federal_Climate_Zone"])
    )

    county["Climate_Zone_Factor"] = (
        county["Final_Climate_Zone"]
        .map(CLIMATE_LIFESPAN_MODIFIER)
    )

    # -------------------------------------------------------------------------
    # 5. Document the source of each county assignment.
    # -------------------------------------------------------------------------
    ltpp_basis = (
        "LTPP General Information | "
        + county["Final_Climate_Zone"].fillna("")
    )

    federal_basis = (
        "Federal NASA POWER 1991-2020; "
        + county["Federal_Sample_Points"]
        .fillna(0)
        .astype(int)
        .astype(str)
        + " county points; precip="
        + county["Federal_Annual_Precipitation_mm"]
        .round(1)
        .astype(str)
        + " mm/yr; FI="
        + county["Federal_Freezing_Index_C_Days"]
        .round(1)
        .astype(str)
        + " C-days | "
        + county["Final_Climate_Zone"].fillna("")
    )

    county["Climate_Basis"] = np.select(
        [
            county["LTPP_Climate_Zone"].notna(),
            county["Federal_Climate_Zone"].notna(),
        ],
        [
            ltpp_basis,
            federal_basis,
        ],
        default="Unresolved",
    )

    # -------------------------------------------------------------------------
    # 6. Calculate state lane-mile-weighted fallback factors.
    #
    # This is used only when both the LTPP and federal county
    # classifications are unavailable.
    # -------------------------------------------------------------------------
    valid = county["Climate_Zone_Factor"].notna()

    state_fallback_data = county.loc[
        valid,
        [
            "state_abbr",
            "Climate_Zone_Factor",
            "Total_County_Lane_Miles",
        ],
    ].copy()

    # Make sure lane-mile weights are numeric.
    state_fallback_data["Total_County_Lane_Miles"] = (
        pd.to_numeric(
            state_fallback_data[
                "Total_County_Lane_Miles"
            ],
            errors="coerce",
        )
        .fillna(0.0)
    )

    # Avoid a zero-weight division.
    state_fallback_data["Climate_Weight"] = (
        state_fallback_data[
            "Total_County_Lane_Miles"
        ]
        .clip(lower=0.000001)
    )

    state_fallback_data["Weighted_Climate_Factor"] = (
        state_fallback_data[
            "Climate_Zone_Factor"
        ]
        * state_fallback_data[
            "Climate_Weight"
        ]
    )

    state_fallback = (
        state_fallback_data
        .groupby(
            "state_abbr",
            as_index=False,
        )
        .agg(
            Weighted_Factor_Sum=(
                "Weighted_Climate_Factor",
                "sum",
            ),
            Climate_Weight_Sum=(
                "Climate_Weight",
                "sum",
            ),
        )
    )

    state_fallback["State_Climate_Fallback"] = (
        state_fallback["Weighted_Factor_Sum"]
        / state_fallback["Climate_Weight_Sum"]
    )

    state_fallback = state_fallback[
        [
            "state_abbr",
            "State_Climate_Fallback",
        ]
    ]

    # -------------------------------------------------------------------------
    # 7. Merge state fallback factors into the county table.
    # -------------------------------------------------------------------------
    county = county.merge(
        state_fallback,
        on="state_abbr",
        how="left",
        validate="many_to_one",
    )

    use_state_fallback = (
        county["Climate_Zone_Factor"].isna()
        & county["State_Climate_Fallback"].notna()
    )

    county.loc[
        use_state_fallback,
        "Climate_Zone_Factor",
    ] = county.loc[
        use_state_fallback,
        "State_Climate_Fallback",
    ]

    county.loc[
        use_state_fallback,
        "Climate_Basis",
    ] = (
        "State lane-mile-weighted fallback; "
        "LTPP and federal county climate unavailable"
    )

    # -------------------------------------------------------------------------
    # 8. Apply Connecticut legacy-county climate fallback.
    #
    # Current Census TIGER data use Connecticut planning regions, whereas
    # HPMS retains the eight legacy county FIPS codes. Apply Wet-Freeze only
    # when LTPP, federal county climate, and state fallback are unavailable.
    # -------------------------------------------------------------------------
    use_connecticut_fallback = (
        county["state_abbr"].eq("CT")
        & county["Climate_Zone_Factor"].isna()
    )

    county.loc[
        use_connecticut_fallback,
        "Final_Climate_Zone",
    ] = "Wet-Freeze"

    county.loc[
        use_connecticut_fallback,
        "Climate_Zone_Factor",
    ] = CLIMATE_LIFESPAN_MODIFIER["Wet-Freeze"]

    county.loc[
        use_connecticut_fallback,
        "Climate_Basis",
    ] = (
        "Connecticut Wet-Freeze state fallback; legacy county FIPS "
        "unavailable in current Census TIGER boundaries"
    )

    # -------------------------------------------------------------------------
    # 9. Apply a neutral factor only when no county or state estimate exists.
    # -------------------------------------------------------------------------
    unresolved = county[
        "Climate_Zone_Factor"
    ].isna()

    county.loc[
        unresolved,
        "Climate_Zone_Factor",
    ] = 1.0

    county.loc[
        unresolved,
        "Climate_Basis",
    ] = (
        "National neutral fallback; "
        "LTPP, federal county climate, and state fallback unavailable"
    )

    # -------------------------------------------------------------------------
    # 10. Final integrity checks.
    # -------------------------------------------------------------------------
    if county["county_fips"].duplicated().any():
        duplicate_fips = county.loc[
            county["county_fips"].duplicated(
                keep=False
            ),
            "county_fips",
        ].unique()

        raise RuntimeError(
            "Duplicate county FIPS found after climate assignment: "
            + ", ".join(duplicate_fips[:20])
        )

    if county["Climate_Zone_Factor"].isna().any():
        raise RuntimeError(
            "Climate assignment finished with missing factors."
        )

    return county

# =============================================================================
# 9. MERGE THE TWO REQUESTED COLUMNS INTO THE MASTER MATRIX
# =============================================================================

def main():
    print("Reading master infrastructure matrix...")
    master = read_master_file()

    print("Extracting LTPP county climate classifications...")
    ltpp_lookup = build_ltpp_county_lookup(master)

    unique_counties = (
        master[
            ["state_abbr", "county_fips", "County_Name"]
        ]
        .drop_duplicates()
    )

    counties_without_unanimous_ltpp = (
        unique_counties.merge(
            ltpp_lookup[
                ["county_fips", "LTPP_Climate_Zone"]
            ],
            on="county_fips",
            how="left",
        )
    )

    counties_without_unanimous_ltpp = (
        counties_without_unanimous_ltpp[
            counties_without_unanimous_ltpp[
                "LTPP_Climate_Zone"
            ].isna()
        ]
    )

    print(
        f"Counties requiring federal climate assignment: "
        f"{len(counties_without_unanimous_ltpp):,}"
    )

    census_counties = load_census_counties(
        unique_counties
    )

    federal_lookup = build_federal_climate_lookup(
        counties_without_unanimous_ltpp,
        census_counties,
    )

    county_assignment = assign_final_climate(
        master,
        ltpp_lookup,
        federal_lookup,
    )

    # Write the detailed one-row-per-county audit file.
    county_assignment.to_csv(
        COUNTY_CLIMATE_AUDIT,
        index=False,
    )

    # Add exactly the two requested columns to the master matrix.
    merge_columns = county_assignment[
        [
            "county_fips",
            "Climate_Zone_Factor",
            "Climate_Basis",
        ]
    ]

    output = master.drop(
        columns=[
            "Climate_Zone_Factor",
            "Climate_Basis",
        ],
        errors="ignore",
    )

    output = output.merge(
        merge_columns,
        on="county_fips",
        how="left",
        validate="many_to_one",
    )

    output.to_csv(
        MASTER_OUTPUT,
        index=False,
    )

    # -------------------------------------------------------------------------
    # Integrity checks
    # -------------------------------------------------------------------------
    if len(output) != len(master):
        raise RuntimeError(
            "Row count changed during the county climate merge."
        )

    if output["Climate_Zone_Factor"].isna().any():
        raise RuntimeError(
            "Some master rows still have missing climate factors."
        )

    print("\nClimate assignment completed.")
    print(f"Master output:\n{MASTER_OUTPUT}")
    print(f"County audit:\n{COUNTY_CLIMATE_AUDIT}")
    print(f"Federal cache:\n{FEDERAL_CACHE}")

    print("\nAssignment basis:")
    print(
        county_assignment["Climate_Basis"]
        .str.split(" | ", regex=False)
        .str[0]
        .value_counts(dropna=False)
    )

    print("\nFinal climate-zone counts:")
    print(
        county_assignment[
            "Final_Climate_Zone"
        ].value_counts(dropna=False)
    )

    print("\nLTPP counties with conflicting section labels:")
    conflicts = county_assignment[
        county_assignment["LTPP_Zone_Count"].fillna(0) > 1
    ][
        [
            "state_abbr",
            "County_Name",
            "LTPP_Zones",
            "Federal_Climate_Zone",
        ]
    ]

    print(
        conflicts.to_string(index=False)
        if not conflicts.empty
        else "None"
    )


if __name__ == "__main__":
    main()
