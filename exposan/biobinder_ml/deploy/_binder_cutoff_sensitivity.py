# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# !!! Temporarily ignoring warnings
import warnings
warnings.filterwarnings('ignore')
from openpyxl.cell.cell import (
    ILLEGAL_CHARACTERS_RE,
)
import os
import sys
import importlib.util
from datetime import datetime

import numpy as np
import pandas as pd


# ============================================================
# 1. USER SETTINGS
# ============================================================

DISTILLATION_POLICY = "max_Lr_then_max_Hr"
# Valid choices:
# "max_top_ratio"
# "min_err"
# "max_Hr_then_max_Lr"
# "max_Lr_then_max_Hr"

LIGHT_FRACTION = 0.01

# Medium fraction values:
# 0.30, 0.40, 0.50, 0.60, 0.70
MEDIUM_FRACTIONS = np.arange(0.30, 0.71, 0.10)

feedstocks = ["food", "green", "manure"]

config_kwargs = {
    "decentralized_HTL": True,
    "decentralized_upgrading": False,
    "skip_EC": True,
    "generate_H2": False,
    "EC_config": None,
}


# ============================================================
# 2. PATHS
# ============================================================

dist_flex_path = (
    r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy"
    r"\Data\Billion Ton\Dist_flex.py"
)

input_excel = (
    r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy"
    r"\Data\Billion Ton"
    r"\county_waste_scenarios_comparison_near_term.xlsx"
)

output_directory = (
    r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy"
    r"\Data\Billion Ton"
)

timestamp_file = datetime.now().strftime("%Y%m%d_%H%M")

output_excel = os.path.join(
    output_directory,
    (
        "county_biobinder_cutoff_sensitivity_"
        f"{DISTILLATION_POLICY}_{timestamp_file}.xlsx"
    ),
)

os.makedirs(output_directory, exist_ok=True)


# ============================================================
# 3. IMPORT DIST_FLEX
# ============================================================

dist_dir = os.path.dirname(dist_flex_path)

if dist_dir not in sys.path:
    sys.path.insert(0, dist_dir)

spec = importlib.util.spec_from_file_location(
    "Dist_flex",
    dist_flex_path,
)

if spec is None or spec.loader is None:
    raise ImportError(
        f"Could not create import specification for:\n{dist_flex_path}"
    )

Dist_flex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(Dist_flex)

create_system = Dist_flex.create_system


# ============================================================
# 4. CREATE CUTOFF SCENARIOS
# ============================================================

cutoff_scenarios = []

for medium_fraction in MEDIUM_FRACTIONS:

    medium_fraction = round(float(medium_fraction), 2)

    heavy_fraction = round(
        1.0 - LIGHT_FRACTION - medium_fraction,
        2,
    )

    cutoff_fracs = [
        LIGHT_FRACTION,
        medium_fraction,
        heavy_fraction,
    ]

    if not np.isclose(sum(cutoff_fracs), 1.0, atol=1e-8):
        raise ValueError(
            f"Invalid cutoff fractions: {cutoff_fracs}; "
            f"sum={sum(cutoff_fracs)}"
        )

    cutoff_scenarios.append(
        {
            "Cutoff_scenario": (
                f"L{int(round(LIGHT_FRACTION * 100)):02d}_"
                f"M{int(round(medium_fraction * 100)):02d}_"
                f"H{int(round(heavy_fraction * 100)):02d}"
            ),
            "Cutoff_frac_light": LIGHT_FRACTION,
            "Cutoff_frac_medium": medium_fraction,
            "Cutoff_frac_heavy": heavy_fraction,
            "Cutoff_fracs": cutoff_fracs,
        }
    )


# ============================================================
# 5. GENERAL HELPERS
# ============================================================

def safe_float(value):
    """Convert value to finite float or return NaN."""
    try:
        value = float(value)
        return value if np.isfinite(value) else np.nan
    except Exception:
        return np.nan


def safe_bool(value):
    """Convert value to bool where possible."""
    try:
        return bool(value)
    except Exception:
        return False


def get_stream(flowsheet, names):
    """Return the first matching stream from a list of IDs."""
    for name in names:
        if hasattr(flowsheet.stream, name):
            return getattr(flowsheet.stream, name)
    return None


def safe_water_mass(stream):
    """Return water mass flow in kg/hr."""
    try:
        return float(stream.imass["Water"])
    except Exception:
        return 0.0


def annual_kg(stream, system):
    """Convert stream flow from kg/hr to kg/yr."""
    if stream is None:
        return np.nan

    return safe_float(
        stream.F_mass * system.operating_hours
    )


def get_yield_value(yields_dict, *keys):
    """Return the first available yield value."""
    for key in keys:
        if key in yields_dict:
            return safe_float(yields_dict[key])
    return np.nan

def clean_excel_value(value):
    """
    Remove control characters that cannot be written
    to an XLSX cell.
    """
    if isinstance(value, str):
        return ILLEGAL_CHARACTERS_RE.sub(
            "",
            value,
        )

    return value


def clean_dataframe_for_excel(df):
    """
    Clean all string cells before writing with openpyxl.
    """
    if df is None or df.empty:
        return df

    cleaned = df.copy()

    object_columns = cleaned.select_dtypes(
        include=["object"]
    ).columns

    for column in object_columns:
        cleaned[column] = cleaned[column].map(
            clean_excel_value
        )

    return cleaned

# ============================================================
# 6. RICH SIMULATION AND DISTILLATION DIAGNOSTICS
# ============================================================

def collect_diagnostics(
    flowsheet,
    system,
    feedstock,
    cutoff_info,
    distillation_policy,
):
    """
    Collect feed, yield, product, separation-policy, VLE,
    design, and conversion-factor diagnostics.
    """

    hours = safe_float(
        system.operating_hours
    )

    # ========================================================
    # STREAMS
    # ========================================================

    feed = get_stream(
        flowsheet,
        [
            "scaled_feedstock",
            "feedstock",
            "Feedstock",
            "HTL_feed",
        ],
    )

    biocrude = get_stream(
        flowsheet,
        [
            "scaled_biocrude",
            "splitted_crude",
            "crude_to_dist",
            "transported_biocrude",
            "biocrude",
        ],
    )

    biobinder = get_stream(
        flowsheet,
        [
            "biobinder",
            "cooled_biobinder",
            "hot_biobinder",
        ],
    )

    biofuel = get_stream(
        flowsheet,
        [
            "biofuel",
            "cooled_biofuel",
            "hot_biofuel",
        ],
    )

    if feed is None:
        raise RuntimeError(
            f"Could not find feed stream for {feedstock}."
        )

    if biocrude is None:
        raise RuntimeError(
            f"Could not find scaled biocrude stream for {feedstock}."
        )

    if biobinder is None:
        raise RuntimeError(
            f"Could not find biobinder stream for {feedstock}."
        )

    if biofuel is None:
        raise RuntimeError(
            f"Could not find biofuel stream for {feedstock}."
        )

    # ========================================================
    # FEED AND PRODUCT FLOWS
    # ========================================================

    water_kg_hr = safe_water_mass(
        feed
    )

    dry_feed_kg_hr = safe_float(
        feed.F_mass - water_kg_hr
    )

    dry_feed_kg_yr = (
        dry_feed_kg_hr * hours
    )

    biocrude_kg_yr = annual_kg(
        biocrude,
        system,
    )

    biobinder_kg_yr = annual_kg(
        biobinder,
        system,
    )

    biofuel_kg_yr = annual_kg(
        biofuel,
        system,
    )

    product_sum_kg_yr = (
        biobinder_kg_yr
        + biofuel_kg_yr
    )

    # ========================================================
    # HTL YIELDS
    # ========================================================

    htl_yields_used = dict(
        Dist_flex.HTL_yields
    )

    Y_biocrude = get_yield_value(
        htl_yields_used,
        "biocrude",
    )

    Y_aqueous = get_yield_value(
        htl_yields_used,
        "aqueous",
    )

    Y_gas = get_yield_value(
        htl_yields_used,
        "gas",
    )

    Y_char = get_yield_value(
        htl_yields_used,
        "char",
        "solid",
        "solids",
    )

    # ========================================================
    # DISTILLATION UNITS AND POLICY
    # ========================================================

    col = flowsheet.unit.CrudeHeavyDis
    splitter = flowsheet.unit.BiocrudeSplitter
    flash = flowsheet.unit.BiofuelFlash

    runtime_best = (
        getattr(
            col,
            "_runtime_best",
            None,
        )
        or {}
    )

    actual_cutoffs = list(
        splitter.cutoff_fracs
    )

    column_top_ratio = (
        safe_float(
            col.outs[0].F_mass
            / col.F_mass_out
        )
        if col.F_mass_out
        else np.nan
    )

    column_bottom_ratio = (
        1.0 - column_top_ratio
        if np.isfinite(
            column_top_ratio
        )
        else np.nan
    )

    flash_vapor_loss_fraction = (
        safe_float(
            flash.outs[0].F_mass
            / col.outs[0].F_mass
        )
        if col.outs[0].F_mass
        else np.nan
    )

    policy_target_ratio = safe_float(
        runtime_best.get(
            "target_ratio"
        )
    )

    policy_search_ratio = safe_float(
        runtime_best.get(
            "search_top_ratio"
        )
    )

    policy_search_error = safe_float(
        runtime_best.get(
            "search_error"
        )
    )

    policy_final_ratio = safe_float(
        runtime_best.get(
            "top_ratio"
        )
    )

    policy_final_error = safe_float(
        runtime_best.get(
            "err"
        )
    )

    policy_tolerance = safe_float(
        runtime_best.get(
            "tolerance"
        )
    )

    policy_stream_ratio_difference = (
        abs(
            policy_final_ratio
            - column_top_ratio
        )
        if (
            np.isfinite(
                policy_final_ratio
            )
            and np.isfinite(
                column_top_ratio
            )
        )
        else np.nan
    )

    policy_final_ratio_matches_stream = (
        bool(
            policy_stream_ratio_difference
            <= 1e-6
        )
        if np.isfinite(
            policy_stream_ratio_difference
        )
        else False
    )

    policy_target_satisfied = safe_bool(
        runtime_best.get(
            "target_satisfied",
            (
                policy_final_error
                <= policy_tolerance
                if (
                    np.isfinite(
                        policy_final_error
                    )
                    and np.isfinite(
                        policy_tolerance
                    )
                )
                else False
            ),
        )
    )

    # ========================================================
    # MAIN OUTPUT ROW
    # ========================================================

    row = {
        # --------------------------------------------
        # Run identifiers
        # --------------------------------------------

        "Cutoff_scenario":
            cutoff_info[
                "Cutoff_scenario"
            ],

        "Feedstock":
            feedstock,

        "Simulation_OK":
            True,

        "Simulation_error":
            "",

        "Distillation_policy_requested":
            distillation_policy,

        "Distillation_policy_selected":
            runtime_best.get(
                "prefer",
                "",
            ),

        # --------------------------------------------
        # Requested cutoff assumptions
        # --------------------------------------------

        "Requested_cutoff_frac_light":
            cutoff_info[
                "Cutoff_frac_light"
            ],

        "Requested_cutoff_frac_medium":
            cutoff_info[
                "Cutoff_frac_medium"
            ],

        "Requested_cutoff_frac_heavy":
            cutoff_info[
                "Cutoff_frac_heavy"
            ],

        "Actual_cutoff_frac_light":
            safe_float(
                actual_cutoffs[0]
            ),

        "Actual_cutoff_frac_medium":
            safe_float(
                actual_cutoffs[1]
            ),

        "Actual_cutoff_frac_heavy":
            safe_float(
                actual_cutoffs[2]
            ),

        "Cutoff_temperature_light_C":
            150.0,

        "Cutoff_temperature_heavy_C":
            343.0,

        # --------------------------------------------
        # Feed and operating scale
        # --------------------------------------------

        "Operating_hours_per_year":
            hours,

        "Total_feed_kg_per_hr":
            safe_float(
                feed.F_mass
            ),

        "Feed_water_kg_per_hr":
            water_kg_hr,

        "Dry_feed_kg_per_hr":
            dry_feed_kg_hr,

        "Dry_feed_kg_per_year":
            dry_feed_kg_yr,

        "Dry_feed_metric_ton_per_year":
            dry_feed_kg_yr / 1000,

        "Dry_feed_metric_ton_per_day":
            (
                dry_feed_kg_yr
                / 1000
                / 365
            ),

        # --------------------------------------------
        # RF/HTL yields
        # --------------------------------------------

        "Y_biocrude":
            Y_biocrude,

        "Y_aqueous":
            Y_aqueous,

        "Y_gas":
            Y_gas,

        "Y_char":
            Y_char,

        "Yield_sum":
            np.nansum(
                [
                    Y_biocrude,
                    Y_aqueous,
                    Y_gas,
                    Y_char,
                ]
            ),

        # --------------------------------------------
        # Product flows
        # --------------------------------------------

        "Scaled_biocrude_kg_per_hr":
            safe_float(
                biocrude.F_mass
            ),

        "Scaled_biocrude_kg_per_year":
            biocrude_kg_yr,

        "Biobinder_kg_per_hr":
            safe_float(
                biobinder.F_mass
            ),

        "Biobinder_kg_per_year":
            biobinder_kg_yr,

        "Biofuel_kg_per_hr":
            safe_float(
                biofuel.F_mass
            ),

        "Biofuel_kg_per_year":
            biofuel_kg_yr,

        # --------------------------------------------
        # Conversion factors
        # --------------------------------------------

        "Stream_biocrude_kg_per_kg_dry_feed":
            (
                biocrude_kg_yr
                / dry_feed_kg_yr
                if dry_feed_kg_yr > 0
                else np.nan
            ),

        "Biobinder_kg_per_kg_dry_feed":
            (
                biobinder_kg_yr
                / dry_feed_kg_yr
                if dry_feed_kg_yr > 0
                else np.nan
            ),

        "Biofuel_kg_per_kg_dry_feed":
            (
                biofuel_kg_yr
                / dry_feed_kg_yr
                if dry_feed_kg_yr > 0
                else np.nan
            ),

        "Biobinder_fraction_of_scaled_biocrude":
            (
                biobinder_kg_yr
                / biocrude_kg_yr
                if biocrude_kg_yr > 0
                else np.nan
            ),

        "Biofuel_fraction_of_scaled_biocrude":
            (
                biofuel_kg_yr
                / biocrude_kg_yr
                if biocrude_kg_yr > 0
                else np.nan
            ),

        "Biobinder_fraction_of_final_products":
            (
                biobinder_kg_yr
                / product_sum_kg_yr
                if product_sum_kg_yr > 0
                else np.nan
            ),

        "Biofuel_fraction_of_final_products":
            (
                biofuel_kg_yr
                / product_sum_kg_yr
                if product_sum_kg_yr > 0
                else np.nan
            ),

        "Product_mass_recovery_from_biocrude":
            (
                product_sum_kg_yr
                / biocrude_kg_yr
                if biocrude_kg_yr > 0
                else np.nan
            ),

        # --------------------------------------------
        # Runtime separation policy
        # --------------------------------------------

        "SEP_policy_found":
            safe_bool(
                runtime_best
            ),

        "SEP_policy_LHK":
            str(
                runtime_best.get(
                    "LHK",
                    "",
                )
            ),

        "SEP_policy_Lr":
            safe_float(
                runtime_best.get(
                    "Lr"
                )
            ),

        "SEP_policy_Hr":
            safe_float(
                runtime_best.get(
                    "Hr"
                )
            ),

        "SEP_policy_target_ratio":
            policy_target_ratio,

        "SEP_policy_tolerance":
            policy_tolerance,

        "SEP_policy_search_top_ratio":
            policy_search_ratio,

        "SEP_policy_search_error":
            policy_search_error,

        "SEP_policy_final_top_ratio":
            policy_final_ratio,

        "SEP_policy_final_error":
            policy_final_error,

        "SEP_policy_target_satisfied":
            policy_target_satisfied,

        "SEP_policy_candidate_count":
            safe_float(
                runtime_best.get(
                    "candidate_count"
                )
            ),

        "SEP_policy_LHK_candidate_count":
            safe_float(
                runtime_best.get(
                    "LHK_candidate_count"
                )
            ),

        "SEP_policy_stream_ratio_difference":
            policy_stream_ratio_difference,

        "SEP_policy_final_ratio_matches_stream":
            policy_final_ratio_matches_stream,

        # Backward-compatible fields
        "SEP_policy_top_ratio":
            policy_final_ratio,

        "SEP_policy_error":
            policy_final_error,

        # --------------------------------------------
        # Final applied separation settings
        # --------------------------------------------

        "SEP_LHK_light":
            (
                col.LHK[0]
                if isinstance(
                    col.LHK,
                    tuple,
                )
                else str(
                    col.LHK
                )
            ),

        "SEP_LHK_heavy":
            (
                col.LHK[1]
                if isinstance(
                    col.LHK,
                    tuple,
                )
                else ""
            ),

        "SEP_Lr_final":
            safe_float(
                col.Lr
            ),

        "SEP_Hr_final":
            safe_float(
                col.Hr
            ),

        "SEP_column_top_ratio":
            column_top_ratio,

        "SEP_column_bottom_ratio":
            column_bottom_ratio,

        "SEP_flash_vapor_loss_fraction":
            flash_vapor_loss_fraction,

        "Column_feed_kg_per_hr":
            safe_float(
                col.ins[0].F_mass
            ),

        "Column_distillate_kg_per_hr":
            safe_float(
                col.outs[0].F_mass
            ),

        "Column_bottoms_kg_per_hr":
            safe_float(
                col.outs[1].F_mass
            ),

        "Column_mass_balance_error_kg_per_hr":
            safe_float(
                col.ins[0].F_mass
                - col.outs[0].F_mass
                - col.outs[1].F_mass
            ),
    }

    # ========================================================
    # COLUMN DESIGN RESULTS
    # ========================================================

    design = (
        col.design_results
        or {}
    )

    row.update(
        {
            "CHD_theoretical_stages":
                safe_float(
                    design.get(
                        "Theoretical stages",
                        np.nan,
                    )
                ),

            "CHD_reflux":
                safe_float(
                    design.get(
                        "Reflux",
                        np.nan,
                    )
                ),

            "CHD_minimum_reflux":
                safe_float(
                    design.get(
                        "Minimum reflux",
                        np.nan,
                    )
                ),

            "CHD_rectifier_diameter_ft":
                safe_float(
                    design.get(
                        "Rectifier diameter",
                        design.get(
                            "Diameter",
                            np.nan,
                        ),
                    )
                ),

            "CHD_rectifier_height_ft":
                safe_float(
                    design.get(
                        "Rectifier height",
                        design.get(
                            "Height",
                            np.nan,
                        ),
                    )
                ),

            "CHD_wall_thickness_in":
                safe_float(
                    design.get(
                        "Rectifier wall thickness",
                        design.get(
                            "Wall thickness",
                            np.nan,
                        ),
                    )
                ),

            "CHD_weight_lb":
                safe_float(
                    design.get(
                        "Rectifier weight",
                        design.get(
                            "Weight",
                            np.nan,
                        ),
                    )
                ),

            "CHD_purchase_cost":
                safe_float(
                    np.nansum(
                        list(
                            (
                                col.baseline_purchase_costs
                                or {}
                            ).values()
                        )
                    )
                ),

            "CHD_utility_cost_per_hr":
                safe_float(
                    getattr(
                        col,
                        "utility_cost",
                        np.nan,
                    )
                ),
        }
    )

    # ========================================================
    # VLE / RELATIVE-VOLATILITY DIAGNOSTICS
    # ========================================================

    try:
        IDs = getattr(
            col,
            "_IDs_vle",
            None,
        )

        LHK_vle = getattr(
            col,
            "_LHK_vle_index",
            None,
        )

        if not IDs or LHK_vle is None:
            raise RuntimeError(
                "Missing _IDs_vle or _LHK_vle_index."
            )

        LK_vle = int(
            LHK_vle[0]
        )

        HK_vle = int(
            LHK_vle[1]
        )

        Dstream, Bstream = col.outs

        zD = Dstream.get_normalized_mol(
            IDs
        )

        zB = Bstream.get_normalized_mol(
            IDs
        )

        dew = getattr(
            col,
            "_dew_point",
            None,
        )

        bubble = getattr(
            col,
            "_bubble_point",
            None,
        )

        if dew is None or bubble is None:
            raise RuntimeError(
                "Missing dew- or bubble-point solver."
            )

        if getattr(
            col,
            "_partial_condenser",
            True,
        ):
            dp = dew(
                zD,
                P=col.P,
            )
        else:
            dp = bubble(
                zD,
                P=col.P,
            )

        bp = bubble(
            zB,
            P=col.P,
        )

        eps = 1e-16

        Kdist = (
            dp.y
            / np.maximum(
                dp.x,
                eps,
            )
        )

        Kbot = (
            bp.y
            / np.maximum(
                bp.x,
                eps,
            )
        )

        alpha_mean = (
            col
            ._estimate_mean_volatilities_relative_to_heavy_key()
        )

        row.update(
            {
                "CHD_pressure_Pa":
                    safe_float(
                        col.P
                    ),

                "CHD_partial_condenser":
                    safe_bool(
                        getattr(
                            col,
                            "_partial_condenser",
                            True,
                        )
                    ),

                "CHD_VLE_component_count":
                    len(IDs),

                "CHD_distillate_dewpoint_K":
                    safe_float(
                        dp.T
                    ),

                "CHD_bottoms_bubblepoint_K":
                    safe_float(
                        bp.T
                    ),

                "CHD_distillate_min_x":
                    safe_float(
                        np.min(
                            dp.x
                        )
                    ),

                "CHD_bottoms_min_x":
                    safe_float(
                        np.min(
                            bp.x
                        )
                    ),

                "CHD_zD_LK":
                    safe_float(
                        zD[LK_vle]
                    ),

                "CHD_zD_HK":
                    safe_float(
                        zD[HK_vle]
                    ),

                "CHD_zB_LK":
                    safe_float(
                        zB[LK_vle]
                    ),

                "CHD_zB_HK":
                    safe_float(
                        zB[HK_vle]
                    ),

                "CHD_Kdist_LK":
                    safe_float(
                        Kdist[LK_vle]
                    ),

                "CHD_Kdist_HK":
                    safe_float(
                        Kdist[HK_vle]
                    ),

                "CHD_Kbottom_LK":
                    safe_float(
                        Kbot[LK_vle]
                    ),

                "CHD_Kbottom_HK":
                    safe_float(
                        Kbot[HK_vle]
                    ),

                "CHD_alpha_dist_LK_over_HK":
                    safe_float(
                        Kdist[LK_vle]
                        / max(
                            Kdist[HK_vle],
                            eps,
                        )
                    ),

                "CHD_alpha_bottom_LK_over_HK":
                    safe_float(
                        Kbot[LK_vle]
                        / max(
                            Kbot[HK_vle],
                            eps,
                        )
                    ),

                "CHD_alpha_LK_used":
                    safe_float(
                        alpha_mean[LK_vle]
                    ),

                "CHD_alpha_LK_le_1":
                    safe_bool(
                        alpha_mean[LK_vle]
                        <= 1.0
                    ),

                "CHD_max_K_observed":
                    safe_float(
                        max(
                            np.nanmax(
                                Kdist
                            ),
                            np.nanmax(
                                Kbot
                            ),
                        )
                    ),

                "CHD_VLE_error":
                    "",
            }
        )

        # ====================================================
        # FENSKE DIAGNOSTICS
        # ====================================================

        try:
            LK_index, HK_index = (
                col._LHK_index
            )

            LK_D = safe_float(
                Dstream.mol[
                    LK_index
                ]
            )

            HK_D = safe_float(
                Dstream.mol[
                    HK_index
                ]
            )

            LK_B = safe_float(
                Bstream.mol[
                    LK_index
                ]
            )

            HK_B = safe_float(
                Bstream.mol[
                    HK_index
                ]
            )

            fenske_ratio = (
                LK_D
                / max(
                    HK_D,
                    1e-30,
                )
                * HK_B
                / max(
                    LK_B,
                    1e-30,
                )
            )

            row.update(
                {
                    "CHD_Fenske_ratio":
                        safe_float(
                            fenske_ratio
                        ),

                    "CHD_log10_Fenske_ratio":
                        safe_float(
                            np.log10(
                                max(
                                    fenske_ratio,
                                    1e-30,
                                )
                            )
                        ),

                    "CHD_log10_alpha_LK":
                        safe_float(
                            np.log10(
                                max(
                                    alpha_mean[
                                        LK_vle
                                    ],
                                    1e-30,
                                )
                            )
                        ),
                }
            )

        except Exception as error:
            row[
                "CHD_Fenske_error"
            ] = str(error)

        # ====================================================
        # UNDERWOOD DIAGNOSTICS
        # ====================================================

        try:
            feed_quality = (
                col.get_feed_quality()
            )

            theta = (
                col
                ._solve_Underwood_constant(
                    alpha_mean,
                    alpha_mean[
                        LK_vle
                    ],
                )
            )

            z_distillate = (
                Dstream
                .get_normalized_mol(
                    IDs
                )
            )

            minimum_reflux = float(
                (
                    alpha_mean
                    * z_distillate
                    / (
                        alpha_mean
                        - theta
                    )
                ).sum()
                - 1.0
            )

            minimum_reflux_bad = (
                not np.isfinite(
                    minimum_reflux
                )
                or minimum_reflux < 0
            )

            row.update(
                {
                    "CHD_feed_quality_q":
                        safe_float(
                            feed_quality
                        ),

                    "CHD_Underwood_theta":
                        safe_float(
                            theta
                        ),

                    "CHD_minimum_reflux_calculated":
                        safe_float(
                            minimum_reflux
                        ),

                    "CHD_minimum_reflux_bad":
                        safe_bool(
                            minimum_reflux_bad
                        ),
                }
            )

        except Exception as error:
            row[
                "CHD_Underwood_error"
            ] = str(error)

            row[
                "CHD_minimum_reflux_bad"
            ] = True

    except Exception as error:
        row.update(
            {
                "CHD_VLE_error":
                    str(error),

                "CHD_alpha_LK_used":
                    np.nan,

                "CHD_alpha_LK_le_1":
                    np.nan,

                "CHD_minimum_reflux_bad":
                    True,
            }
        )

    # ========================================================
    # OVERALL DISTILLATION-DIAGNOSTIC FLAG
    # ========================================================

    row[
        "Distillation_diagnostics_OK"
    ] = bool(
        row[
            "SEP_policy_found"
        ]
        and row[
            "SEP_policy_final_ratio_matches_stream"
        ]
        and not safe_bool(
            row.get(
                "CHD_alpha_LK_le_1",
                True,
            )
        )
        and not safe_bool(
            row.get(
                "CHD_minimum_reflux_bad",
                True,
            )
        )
        and (
            row.get(
                "CHD_VLE_error",
                "",
            )
            == ""
        )
    )

    return row

# ============================================================
# 7. RUN ALL CUTOFF × FEEDSTOCK SIMULATIONS
# ============================================================

diagnostic_rows = []
conversion_rows = []
failed_rows = []

print("\n" + "=" * 80)
print("RUNNING CUTOFF-FRACTION SENSITIVITY")
print("=" * 80)

for cutoff_info in cutoff_scenarios:

    cutoff_fracs = cutoff_info[
        "Cutoff_fracs"
    ]

    scenario = cutoff_info[
        "Cutoff_scenario"
    ]

    print("\n" + "#" * 80)
    print(
        f"SCENARIO: {scenario}\n"
        f"Light={cutoff_fracs[0]:.2f}, "
        f"Medium={cutoff_fracs[1]:.2f}, "
        f"Heavy={cutoff_fracs[2]:.2f}"
    )
    print("#" * 80)

    for feedstock in feedstocks:

        print(
            f"\nRunning {feedstock} | "
            f"{scenario} | "
            f"policy={DISTILLATION_POLICY}"
        )

        try:
            system = create_system(
                feedstock_id=feedstock,
                cutoff_fracs=cutoff_fracs,
                distillation_policy=(
                    DISTILLATION_POLICY
                ),
                **config_kwargs,
            )

            system.simulate()

            flowsheet = system.flowsheet

            row = collect_diagnostics(
                flowsheet=flowsheet,
                system=system,
                feedstock=feedstock,
                cutoff_info=cutoff_info,
                distillation_policy=(
                    DISTILLATION_POLICY
                ),
            )

            diagnostic_rows.append(row)

            conversion_rows.append(
                {
                    "Cutoff_scenario":
                        scenario,

                    "Cutoff_frac_light":
                        cutoff_fracs[0],

                    "Cutoff_frac_medium":
                        cutoff_fracs[1],

                    "Cutoff_frac_heavy":
                        cutoff_fracs[2],

                    "Feedstock":
                        feedstock,

                    "Distillation_policy":
                        DISTILLATION_POLICY,

                    "Y_biocrude":
                        row["Y_biocrude"],

                    "Stream_biocrude_factor":
                        row[
                            "Stream_biocrude_kg_per_kg_dry_feed"
                        ],

                    "Biobinder_conversion_factor":
                        row[
                            "Biobinder_kg_per_kg_dry_feed"
                        ],

                    "Biofuel_conversion_factor":
                        row[
                            "Biofuel_kg_per_kg_dry_feed"
                        ],

                    "Biobinder_fraction_of_biocrude":
                        row[
                            "Biobinder_fraction_of_scaled_biocrude"
                        ],

                    "Biofuel_fraction_of_biocrude":
                        row[
                            "Biofuel_fraction_of_scaled_biocrude"
                        ],

                    "Selected_LHK_light":
                        row["SEP_LHK_light"],

                    "Selected_LHK_heavy":
                        row["SEP_LHK_heavy"],

                    "Selected_Lr":
                        row["SEP_Lr_final"],

                    "Selected_Hr":
                        row["SEP_Hr_final"],

                    "Selected_column_top_ratio":
                        row[
                            "SEP_column_top_ratio"
                        ],

                    "CHD_alpha_LK":
                        row.get(
                            "CHD_alpha_LK_used",
                            np.nan,
                        ),

                    "CHD_theoretical_stages":
                        row.get(
                            "CHD_theoretical_stages",
                            np.nan,
                        ),

                    "CHD_reflux":
                        row.get(
                            "CHD_reflux",
                            np.nan,
                        ),

                    "Simulation_OK":
                        True,

                    "Simulation_error":
                        "",
                }
            )

            print(
                f"  Binder factor: "
                f"{row['Biobinder_kg_per_kg_dry_feed']:.6f}"
            )

            print(
                f"  Binder/biocrude: "
                f"{row['Biobinder_fraction_of_scaled_biocrude']:.6f}"
            )

        except Exception as error:

            error_text = str(error)

            print(
                f"  FAILED: {error_text}"
            )

            failed_row = {
                "Cutoff_scenario":
                    scenario,

                "Cutoff_frac_light":
                    cutoff_fracs[0],

                "Cutoff_frac_medium":
                    cutoff_fracs[1],

                "Cutoff_frac_heavy":
                    cutoff_fracs[2],

                "Feedstock":
                    feedstock,

                "Distillation_policy":
                    DISTILLATION_POLICY,

                "Simulation_OK":
                    False,

                "Simulation_error":
                    error_text,
            }

            failed_rows.append(failed_row)

            conversion_rows.append(
                {
                    **failed_row,

                    "Y_biocrude":
                        np.nan,

                    "Stream_biocrude_factor":
                        np.nan,

                    "Biobinder_conversion_factor":
                        np.nan,

                    "Biofuel_conversion_factor":
                        np.nan,

                    "Biobinder_fraction_of_biocrude":
                        np.nan,

                    "Biofuel_fraction_of_biocrude":
                        np.nan,

                    "Selected_LHK_light":
                        "",

                    "Selected_LHK_heavy":
                        "",

                    "Selected_Lr":
                        np.nan,

                    "Selected_Hr":
                        np.nan,

                    "Selected_column_top_ratio":
                        np.nan,

                    "CHD_alpha_LK":
                        np.nan,

                    "CHD_theoretical_stages":
                        np.nan,

                    "CHD_reflux":
                        np.nan,
                }
            )


conversion_df = pd.DataFrame(
    conversion_rows
)

diagnostic_df = pd.DataFrame(
    diagnostic_rows
)

failed_df = pd.DataFrame(
    failed_rows
)


# ============================================================
# 8. COUNTY INPUT-PROCESSING FUNCTION
# ============================================================

def process_county_sheet(
    sheet_name,
    suffix,
    conversion_factors,
):
    """
    Apply cutoff- and feedstock-specific binder factors
    to the county dry-waste inventory.
    """

    county_df = pd.read_excel(
        input_excel,
        sheet_name=sheet_name,
    )

    base_columns = [
        "fips",
        "county",
        "state",
    ]

    food_column = (
        f"Food Wastes{suffix}"
        "_Dry_Tons_Year"
    )

    green_column = (
        "Green, Yard, Wood & Agricultural Wastes"
        f"{suffix}_Dry_Tons_Year"
    )

    manure_column = (
        f"Manure Wastes{suffix}"
        "_Dry_Tons_Year"
    )

    required_columns = (
        base_columns
        + [
            food_column,
            green_column,
            manure_column,
        ]
    )

    missing_columns = [
        column
        for column in required_columns
        if column not in county_df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing columns in sheet "
            f"'{sheet_name}':\n"
            f"{missing_columns}\n\n"
            f"Available columns:\n"
            f"{list(county_df.columns)}"
        )

    output = county_df[
        required_columns
    ].copy()

    output = output.rename(
        columns={
            food_column:
                "Food_Dry_Tons_Year",

            green_column:
                "Green_Dry_Tons_Year",

            manure_column:
                "Manure_Dry_Tons_Year",
        }
    )

    output[
        "Food_Biobinder_Dry_Tons_Year"
    ] = (
        output["Food_Dry_Tons_Year"]
        * conversion_factors["food"]
    )

    output[
        "Green_Biobinder_Dry_Tons_Year"
    ] = (
        output["Green_Dry_Tons_Year"]
        * conversion_factors["green"]
    )

    output[
        "Manure_Biobinder_Dry_Tons_Year"
    ] = (
        output["Manure_Dry_Tons_Year"]
        * conversion_factors["manure"]
    )

    output[
        "Total_Biobinder_Dry_Tons_Year"
    ] = (
        output[
            "Food_Biobinder_Dry_Tons_Year"
        ]
        + output[
            "Green_Biobinder_Dry_Tons_Year"
        ]
        + output[
            "Manure_Biobinder_Dry_Tons_Year"
        ]
    )

    return output


# ============================================================
# 9. APPLY CUTOFF-SPECIFIC FACTORS TO COUNTY DATA
# ============================================================

county_results = []
national_summary_rows = []

availability_cases = [
    {
        "Sheet_name":
            "All Waste Available",

        "Suffix":
            "",

        "Availability_case":
            "All Waste Available",
    },
    {
        "Sheet_name":
            "Only Landfilled Portion",

        "Suffix":
            "_Landfilled",

        "Availability_case":
            "Only Landfilled Portion",
    },
]

for cutoff_info in cutoff_scenarios:

    scenario = cutoff_info[
        "Cutoff_scenario"
    ]

    scenario_conversion = conversion_df[
        (
            conversion_df[
                "Cutoff_scenario"
            ]
            == scenario
        )
        & (
            conversion_df[
                "Simulation_OK"
            ]
            == True
        )
    ].copy()

    factors = (
        scenario_conversion
        .set_index("Feedstock")[
            "Biobinder_conversion_factor"
        ]
        .to_dict()
    )

    required_feedstocks = set(
        feedstocks
    )

    available_feedstocks = set(
        factors.keys()
    )

    missing_feedstocks = (
        required_feedstocks
        - available_feedstocks
    )

    if missing_feedstocks:
        print(
            f"\nSkipping county propagation for "
            f"{scenario}; missing valid factors for: "
            f"{sorted(missing_feedstocks)}"
        )
        continue

    for case in availability_cases:

        county_output = process_county_sheet(
            sheet_name=case[
                "Sheet_name"
            ],
            suffix=case["Suffix"],
            conversion_factors=factors,
        )

        county_output.insert(
            0,
            "Distillation_policy",
            DISTILLATION_POLICY,
        )

        county_output.insert(
            0,
            "Cutoff_frac_heavy",
            cutoff_info[
                "Cutoff_frac_heavy"
            ],
        )

        county_output.insert(
            0,
            "Cutoff_frac_medium",
            cutoff_info[
                "Cutoff_frac_medium"
            ],
        )

        county_output.insert(
            0,
            "Cutoff_frac_light",
            cutoff_info[
                "Cutoff_frac_light"
            ],
        )

        county_output.insert(
            0,
            "Availability_case",
            case["Availability_case"],
        )

        county_output.insert(
            0,
            "Cutoff_scenario",
            scenario,
        )

        county_results.append(
            county_output
        )

        national_summary_rows.append(
            {
                "Availability_case":
                    case["Availability_case"],

                "Cutoff_scenario":
                    scenario,

                "Cutoff_frac_light":
                    cutoff_info[
                        "Cutoff_frac_light"
                    ],

                "Cutoff_frac_medium":
                    cutoff_info[
                        "Cutoff_frac_medium"
                    ],

                "Cutoff_frac_heavy":
                    cutoff_info[
                        "Cutoff_frac_heavy"
                    ],

                "Distillation_policy":
                    DISTILLATION_POLICY,

                "Food_conversion_factor":
                    factors["food"],

                "Green_conversion_factor":
                    factors["green"],

                "Manure_conversion_factor":
                    factors["manure"],

                "Food_waste_dry_ton_per_year":
                    county_output[
                        "Food_Dry_Tons_Year"
                    ].sum(),

                "Green_waste_dry_ton_per_year":
                    county_output[
                        "Green_Dry_Tons_Year"
                    ].sum(),

                "Manure_waste_dry_ton_per_year":
                    county_output[
                        "Manure_Dry_Tons_Year"
                    ].sum(),

                "Food_binder_dry_ton_per_year":
                    county_output[
                        "Food_Biobinder_Dry_Tons_Year"
                    ].sum(),

                "Green_binder_dry_ton_per_year":
                    county_output[
                        "Green_Biobinder_Dry_Tons_Year"
                    ].sum(),

                "Manure_binder_dry_ton_per_year":
                    county_output[
                        "Manure_Biobinder_Dry_Tons_Year"
                    ].sum(),

                "Total_binder_dry_ton_per_year":
                    county_output[
                        "Total_Biobinder_Dry_Tons_Year"
                    ].sum(),
            }
        )


if county_results:
    county_cut_df = pd.concat(
        county_results,
        ignore_index=True,
    )
else:
    county_cut_df = pd.DataFrame()


national_summary_df = pd.DataFrame(
    national_summary_rows
)


# ============================================================
# 10. WIDE CONVERSION-FACTOR TABLE
# ============================================================

successful_conversion_df = conversion_df[
    conversion_df["Simulation_OK"] == True
].copy()

if not successful_conversion_df.empty:

    conversion_factor_wide = (
        successful_conversion_df
        .pivot_table(
            index=[
                "Cutoff_scenario",
                "Cutoff_frac_light",
                "Cutoff_frac_medium",
                "Cutoff_frac_heavy",
                "Distillation_policy",
            ],
            columns="Feedstock",
            values=(
                "Biobinder_conversion_factor"
            ),
            aggfunc="first",
        )
        .reset_index()
    )

    conversion_factor_wide.columns.name = None

    conversion_factor_wide = (
        conversion_factor_wide.rename(
            columns={
                "food":
                    "Food_Binder_Factor",

                "green":
                    "Green_Binder_Factor",

                "manure":
                    "Manure_Binder_Factor",
            }
        )
    )

else:
    conversion_factor_wide = pd.DataFrame()


# ============================================================
# 11. EXPORT EXCEL WORKBOOK
# ============================================================

conversion_df = clean_dataframe_for_excel(
    conversion_df
)

conversion_factor_wide = (
    clean_dataframe_for_excel(
        conversion_factor_wide
    )
)

national_summary_df = (
    clean_dataframe_for_excel(
        national_summary_df
    )
)

diagnostic_df = clean_dataframe_for_excel(
    diagnostic_df
)

failed_df = clean_dataframe_for_excel(
    failed_df
)

county_cut_df = clean_dataframe_for_excel(
    county_cut_df
)


# ============================================================
# RUN METADATA
# ============================================================

run_timestamp = datetime.now().strftime(
    "%Y-%m-%d %H:%M:%S"
)

metadata_df = pd.DataFrame(
    {
        "Parameter": [
            "Run timestamp",
            "Distillation policy",
            "Fixed light fraction",
            "Medium fraction minimum",
            "Medium fraction maximum",
            "Medium fraction increment",
            "Heavy fraction calculation",
            "Heavy-column target calculation",
            "Target-ratio tolerance",
            "Tolerance interpretation",
            "Light cutoff temperature",
            "Heavy cutoff temperature",
            "Feedstocks",
            "HTL configuration",
            "HTL yield method",
            "Input county workbook",
            "Dist_flex source",
        ],
        "Value": [
            run_timestamp,
            DISTILLATION_POLICY,
            LIGHT_FRACTION,
            float(
                np.min(
                    MEDIUM_FRACTIONS
                )
            ),
            float(
                np.max(
                    MEDIUM_FRACTIONS
                )
            ),
            0.10,
            (
                "1 - light fraction "
                "- medium fraction"
            ),
            (
                "medium fraction / "
                "(medium fraction + "
                "heavy fraction)"
            ),
            0.10,
            (
                "Diagnostic only; the closest "
                "physically feasible candidate "
                "is retained even when the "
                "target error exceeds tolerance"
            ),
            "150 °C",
            "343 °C",
            ", ".join(feedstocks),
            "DHCU | No EC",
            "RF-predicted HTL yields",
            input_excel,
            dist_flex_path,
        ],
    }
)

metadata_df = clean_dataframe_for_excel(
    metadata_df
)


# ============================================================
# WRITE WORKBOOK
# ============================================================

with pd.ExcelWriter(
    output_excel,
    engine="openpyxl",
) as writer:

    conversion_df.to_excel(
        writer,
        sheet_name="Cut Conversion Factors",
        index=False,
    )

    conversion_factor_wide.to_excel(
        writer,
        sheet_name="Conversion Factors Wide",
        index=False,
    )

    national_summary_df.to_excel(
        writer,
        sheet_name="National Cut Summary",
        index=False,
    )

    diagnostic_df.to_excel(
        writer,
        sheet_name="Distillation Diagnostics",
        index=False,
    )

    failed_df.to_excel(
        writer,
        sheet_name="Failed Simulations",
        index=False,
    )

    county_cut_df.to_excel(
        writer,
        sheet_name="County Cut Potential",
        index=False,
    )

    metadata_df.to_excel(
        writer,
        sheet_name="Run Metadata",
        index=False,
    )


# ============================================================
# FINAL CONSOLE SUMMARY
# ============================================================

print("\n" + "=" * 80)
print("EXECUTION COMPLETE")
print("=" * 80)

print(
    f"Successful simulations: "
    f"{len(diagnostic_df)}"
)

print(
    f"Failed simulations: "
    f"{len(failed_df)}"
)

print(
    f"County scenario rows: "
    f"{len(county_cut_df):,}"
)

if (
    not diagnostic_df.empty
    and "SEP_policy_target_satisfied"
    in diagnostic_df.columns
):
    target_met_count = int(
        diagnostic_df[
            "SEP_policy_target_satisfied"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    target_not_met_count = (
        len(diagnostic_df)
        - target_met_count
    )

    print(
        f"Target ratio satisfied: "
        f"{target_met_count}"
    )

    print(
        f"Target ratio not satisfied: "
        f"{target_not_met_count}"
    )

if (
    not diagnostic_df.empty
    and "Distillation_diagnostics_OK"
    in diagnostic_df.columns
):
    diagnostics_ok_count = int(
        diagnostic_df[
            "Distillation_diagnostics_OK"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    diagnostics_bad_count = (
        len(diagnostic_df)
        - diagnostics_ok_count
    )

    print(
        f"Distillation diagnostics OK: "
        f"{diagnostics_ok_count}"
    )

    print(
        f"Distillation diagnostics flagged: "
        f"{diagnostics_bad_count}"
    )

print(
    f"\nWorkbook saved to:\n"
    f"{output_excel}"
)