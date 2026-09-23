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
RF-predicted HTL yields → full biobinder TEA/LCA (EcoInvent 3.11)

Outputs two prediction treatments:
1. Original: retain unmodified predictions with total yield of 95–105 wt%.
2. Normalized: retain nonnegative predictions with total yield of 90–110 wt%,
   then normalize the four yields to exactly 100 wt%.

Experimental and predicted yields are run through the same full-LCA system.
"""
import warnings
warnings.filterwarnings('ignore')
import os
import flexsolve as flx
import numpy as np
import pandas as pd
from joblib import load
import qsdsan as qs
from chaospy import distributions as shape
from scipy.stats import qmc
from qsdsan.utils import clear_lca_registries
from datetime import datetime
from exposan.biobinder_ml import (
    # create_system,
    _HHV_per_GGE,
    _load_components,
    _load_process_settings,
    _units as u,
    BiobinderTEA,
    central_dry_flowrate as default_central,
    data_path,
    feedstock_composition,
    HTL_yields,
    pilot_dry_flowrate as default_pilot,
    price_dct,
    results_path,
    tea_kwargs,
    uptime_ratio,
    limit_internal_hx,
    configure_utilities,
    set_eff_hx_temperature,
    TRACI_CATEGORIES,
    )

from exposan.biobinder_ml.lca_biobinder import create_system

# ============================================================
# SETTINGS
# ============================================================

OUTPUT_DIR = "results"
MODEL_FILE = os.path.join(OUTPUT_DIR, "rf_yield_model.joblib")
DATA_FILE = os.path.join(OUTPUT_DIR, "HTL_all_yield_experimental.xlsx")
OUTPUT_FILE = os.path.join(
    OUTPUT_DIR,
    f"RF_Yield_Full_TRACI_No_Credit_{datetime.now():%Y%m%d_%H%M}.xlsx",
)

YIELD_COLS = [
    "Biocrude wt%",
    "Aqueous wt%",
    "Gas wt%",
    "Solids wt%",
]

# Change only this mapping if the experimental-yield columns in
# HTL_all_yield_experimental.xlsx use different names.
EXPERIMENTAL_YIELD_COLS = {
    "Biocrude wt%": "Biocrude wt%",
    "Aqueous wt%": "Aqueous wt%",
    "Gas wt%": "Gas wt%",
    "Solids wt%": "Solids wt%",
}

CONFIGS = {
    "CHCU_No_EC": {
        "flowsheet": None,
        "central_dry_flowrate": default_central,
        "decentralized_HTL": False,
        "decentralized_upgrading": False,
        "skip_EC": True,
        "generate_H2": False,
        "EC_config": None,
    },
}

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def prepare_rf_inputs(data, model):
    """Align RF features and retain the original preprocessing method."""
    X = data.reindex(columns=model.feature_names_in_, fill_value=0).copy()

    for col in X.select_dtypes(include="object"):
        X[col] = pd.factorize(X[col])[0]

    X = X.apply(pd.to_numeric, errors="coerce")
    return X.fillna(X.median(numeric_only=True)).fillna(0)


def make_experimental_yields(data):
    """Extract experimental yields using the mapping above."""
    missing = [
        source for source in EXPERIMENTAL_YIELD_COLS.values()
        if source not in data.columns
    ]

    if missing:
        raise KeyError(
            "Experimental yield columns were not found:\n"
            f"{missing}\n\n"
            "Update EXPERIMENTAL_YIELD_COLS at the top of the script."
        )

    out = pd.DataFrame(index=data.index)

    for standard_name, source_name in EXPERIMENTAL_YIELD_COLS.items():
        out[standard_name] = pd.to_numeric(data[source_name], errors="coerce")

    out["Raw Sum wt%"] = out[YIELD_COLS].sum(axis=1)
    return out
def _NPV_at_r(r, CF, dur):
    r = float(r)
    if (not np.isfinite(r)) or (r <= -0.999999):
        return np.nan
    return (CF / (1.0 + r)**dur).sum()

def solve_IRR_aitken_limit(tea, maxiter=200, xtol=1e-6, ytol=10.0, npv_tol_abs=1e3):
    CF = tea.cashflow_array
    dur = tea._get_duration_array()

    if (CF is None) or (dur is None) or (not np.isfinite(CF).all()) or (not np.isfinite(dur).all()):
        return np.nan, "CF_or_dur_nonfinite"

    args = (CF, dur)
    IRR_guess = float(getattr(tea, "IRR", 0.10) or 0.10)
    if (not np.isfinite(IRR_guess)) or (IRR_guess <= -0.999999):
        IRR_guess = 0.10

    x0 = IRR_guess
    x1 = 1.0001 * IRR_guess + 1e-3

    y0 = _NPV_at_r(x0, *args)
    y1 = _NPV_at_r(x1, *args)
    if (not np.isfinite(y0)) or (not np.isfinite(y1)):
        return np.nan, "endpoint_nonfinite"

    try:
        irr = flx.aitken_secant(
            _NPV_at_r,
            x0, x1,
            xtol=xtol, ytol=ytol, maxiter=maxiter,
            args=args, checkiter=False
        )
    except Exception as e:
        return np.nan, f"solver_exception:{type(e).__name__}"

    if (not np.isfinite(irr)) or (irr <= -0.999999):
        return np.nan, "solver_returned_invalid"

    y = _NPV_at_r(irr, *args)
    if not np.isfinite(y):
        return np.nan, "npv_at_solution_nonfinite"
    if abs(float(y)) > float(npv_tol_abs):
        return np.nan, f"npv_not_close:{float(y):.3g}"

    return float(irr), "ok"

def make_prediction_cases(predicted):
    """Create original and normalized prediction datasets."""
    pred = predicted.copy()
    pred["Raw Sum wt%"] = pred[YIELD_COLS].sum(axis=1)
    pred["Mass Balance Error pp"] = abs(pred["Raw Sum wt%"] - 100)

    # Original rule: unchanged predictions totaling 95–105 wt%.
    original_mask = (
        pred["Raw Sum wt%"].between(95, 105)
        & pred[YIELD_COLS].notna().all(axis=1)
    )
    original_valid = pred.loc[original_mask].copy()
    original_skipped = pred.loc[~original_mask].copy()

    # Normalized rule:
    # nonnegative predictions, positive total, and raw total of 90–110 wt%.
    normalized_mask = (
        pred["Raw Sum wt%"].between(90, 110)
        & pred[YIELD_COLS].notna().all(axis=1)
        & pred[YIELD_COLS].ge(0).all(axis=1)
        & pred["Raw Sum wt%"].gt(0)
    )

    normalized_valid = pred.loc[normalized_mask].copy()
    normalized_skipped = pred.loc[~normalized_mask].copy()

    normalized_valid[YIELD_COLS] = (
        normalized_valid[YIELD_COLS]
        .div(normalized_valid["Raw Sum wt%"], axis=0)
        .mul(100)
    )
    normalized_valid["Normalized Sum wt%"] = (
        normalized_valid[YIELD_COLS].sum(axis=1)
    )

    return {
        "original": (original_valid, original_skipped),
        "normalized": (normalized_valid, normalized_skipped),
    }


def impact_column(category):
    """Return a readable per-kg impact column name."""
    unit = TRACI_CATEGORIES[category].get("unit", "")
    return f"{category} ({unit}/kg biobinder)"


def run_full_lca(yield_data, yield_source):
    """Run TEA and all TRACI categories for every yield row."""
    records = []

    for config_name, config_kwargs in CONFIGS.items():
        print(f"\n{'=' * 68}\n{yield_source}: {config_name}\n{'=' * 68}")

        for index, row in yield_data.iterrows():

            # ====================================================
            # YIELD INPUTS
            # ====================================================

            biocrude_wt = float(row["Biocrude wt%"])
            aqueous_wt = float(row["Aqueous wt%"])
            gas_wt = float(row["Gas wt%"])
            solids_wt = float(row["Solids wt%"])

            non_solids_wt = (
                biocrude_wt
                + aqueous_wt
                + gas_wt
            )

            if (
                np.isfinite(non_solids_wt)
                and non_solids_wt > 0
            ):
                biocrude_sf = (
                    biocrude_wt
                    / non_solids_wt
                    * 100.0
                )
                aqueous_sf = (
                    aqueous_wt
                    / non_solids_wt
                    * 100.0
                )
                gas_sf = (
                    gas_wt
                    / non_solids_wt
                    * 100.0
                )

            else:
                biocrude_sf = np.nan
                aqueous_sf = np.nan
                gas_sf = np.nan

            # ====================================================
            # INITIAL RESULT ROW
            # ====================================================

            result = {
                "Config": config_name,
                "Index": index,
                "Yield Source": yield_source,

                # Original yields
                "Biocrude wt%": biocrude_wt,
                "Aqueous wt%": aqueous_wt,
                "Gas wt%": gas_wt,
                "Solids wt%": solids_wt,

                # Solids-free yields
                "Biocrude solids-free wt%": biocrude_sf,
                "Aqueous solids-free wt%": aqueous_sf,
                "Gas solids-free wt%": gas_sf,
                "Solids-free Sum wt%": (
                    biocrude_sf
                    + aqueous_sf
                    + gas_sf
                    if np.isfinite(biocrude_sf)
                    else np.nan
                ),

                # Yield-quality diagnostics
                "Raw Sum wt%": row.get(
                    "Raw Sum wt%",
                    np.nan,
                ),
                "Mass Balance Error pp": row.get(
                    "Mass Balance Error pp",
                    np.nan,
                ),
                "Normalized Sum wt%": row.get(
                    "Normalized Sum wt%",
                    np.nan,
                ),

                # Distillation diagnostics initialized for
                # every row, including failed simulations.
                "Selected LHK": "",
                "LHK Role": "",
                "Fallback Used": np.nan,
                "LHK Candidate Count": np.nan,
                "LHK Attempt Count": np.nan,
                "Verification Failure Count": np.nan,

                "Selected Lr": np.nan,
                "Selected Hr": np.nan,

                "Target Top Ratio": np.nan,
                "Target Tolerance": np.nan,
                "Search Top Ratio": np.nan,
                "Search Error": np.nan,
                "Selected Top Ratio": np.nan,
                "Distillation Error": np.nan,
                "Target Satisfied": np.nan,

                # Process and TEA/LCA outputs initialized
                "Biobinder kg/hr": np.nan,
                "Biofuel kg/hr": np.nan,
                "Annual Biobinder kg/yr": np.nan,
                "MSP ($/kg)": np.nan,
                "IRR raw": np.nan,
                "IRR (%)": np.nan,
                "IRR Solver Status": "",
                "Product Ratio": np.nan,

                **{
                    impact_column(cat): np.nan
                    for cat in TRACI_CATEGORIES
                },

                "Status": "",
                "Error": "",
            }

            try:
                # =================================================
                # INSERT THIS SAMPLE'S YIELDS
                # =================================================

                HTL_yields.update({
                    "biocrude": biocrude_wt / 100.0,
                    "aqueous": aqueous_wt / 100.0,
                    "gas": gas_wt / 100.0,
                    "char": solids_wt / 100.0,
                })

                # Required because lca_biobinder registers the same
                # indicator and impact-item IDs for each new system.
                clear_lca_registries()

                # =================================================
                # CREATE AND SIMULATE SYSTEM
                # =================================================

                system = create_system(
                    **config_kwargs
                )

                system.simulate()

                tea = system.TEA
                lca = system.LCA

                tea.IRR = tea_kwargs["IRR"]
                tea.income_tax = (
                    tea_kwargs["income_tax"]
                )

                biobinder = (
                    system.flowsheet.stream.biobinder
                )
                biofuel = (
                    system.flowsheet.stream.biofuel
                )

                # =================================================
                # DISTILLATION LOGGING
                # =================================================

                dist = (
                    system.flowsheet.unit.CrudeHeavyDis
                )

                dist_info = getattr(
                    dist,
                    "_runtime_best",
                    {},
                ) or {}

                selected_lhk = dist_info.get(
                    "LHK",
                    (),
                )

                result.update({
                    "Selected LHK": (
                        str(tuple(selected_lhk))
                        if selected_lhk
                        else ""
                    ),

                    "LHK Role": dist_info.get(
                        "LHK_role",
                        "",
                    ),

                    "Fallback Used": dist_info.get(
                        "fallback_used",
                        np.nan,
                    ),

                    "LHK Candidate Count": (
                        dist_info.get(
                            "LHK_candidate_count",
                            np.nan,
                        )
                    ),

                    "LHK Attempt Count": (
                        dist_info.get(
                            "LHK_attempt_count",
                            np.nan,
                        )
                    ),

                    "Verification Failure Count": (
                        dist_info.get(
                            "verification_failure_count",
                            np.nan,
                        )
                    ),

                    "Selected Lr": dist_info.get(
                        "Lr",
                        np.nan,
                    ),

                    "Selected Hr": dist_info.get(
                        "Hr",
                        np.nan,
                    ),

                    "Target Top Ratio": (
                        dist_info.get(
                            "target_ratio",
                            np.nan,
                        )
                    ),

                    "Target Tolerance": (
                        dist_info.get(
                            "tolerance",
                            np.nan,
                        )
                    ),

                    "Search Top Ratio": (
                        dist_info.get(
                            "search_top_ratio",
                            np.nan,
                        )
                    ),

                    "Search Error": (
                        dist_info.get(
                            "search_error",
                            np.nan,
                        )
                    ),

                    "Selected Top Ratio": (
                        dist_info.get(
                            "top_ratio",
                            np.nan,
                        )
                    ),

                    "Distillation Error": (
                        dist_info.get(
                            "err",
                            np.nan,
                        )
                    ),

                    "Target Satisfied": (
                        dist_info.get(
                            "target_satisfied",
                            np.nan,
                        )
                    ),
                })

                # A successful simulation using the runtime
                # specification should have _runtime_best.
                if not dist_info:
                    raise RuntimeError(
                        "Distillation simulation completed "
                        "without _runtime_best diagnostics."
                    )

                # =================================================
                # PRODUCT FLOW
                # =================================================

                annual_mass = (
                    biobinder.F_mass
                    * system.operating_hours
                )

                if (
                    not np.isfinite(annual_mass)
                    or annual_mass <= 0
                ):
                    raise ValueError(
                        "Annual biobinder production is "
                        "nonfinite, zero, or negative."
                    )

                # Optional safeguard against numerically successful
                # runs with effectively zero biobinder output.
                #
                # Choose a threshold appropriate for your plant.
                #
                # if annual_mass < 1_000:
                #     raise ValueError(
                #         "Annual biobinder production is too "
                #         f"small ({annual_mass:.6g} kg/yr)."
                #     )

                product_ratio = (
                    biobinder.F_mass
                    / biofuel.F_mass
                    if (
                        np.isfinite(biofuel.F_mass)
                        and biofuel.F_mass > 0
                    )
                    else np.nan
                )

                # =================================================
                # LCA
                # =================================================

                impacts = lca.get_allocated_impacts(
                    streams=(biobinder,),
                    operation_only=True,
                    annual=True,
                )

                impact_results = {
                    impact_column(cat):
                    impacts.get(cat, np.nan)
                    / annual_mass
                    for cat in TRACI_CATEGORIES
                }

                # =================================================
                # TEA
                # =================================================

                msp = tea.solve_price(
                    biobinder
                )

                irr_raw, irr_status = (
                    solve_IRR_aitken_limit(
                        tea,
                        maxiter=200,
                        xtol=1e-6,
                        ytol=10.0,
                        npv_tol_abs=1e3,
                    )
                )

                irr_percent = (
                    irr_raw * 100.0
                    if np.isfinite(irr_raw)
                    else np.nan
                )

                # =================================================
                # SAVE SUCCESSFUL RESULT
                # =================================================

                result.update({
                    "Biobinder kg/hr":
                        biobinder.F_mass,

                    "Biofuel kg/hr":
                        biofuel.F_mass,

                    "Annual Biobinder kg/yr":
                        annual_mass,

                    "MSP ($/kg)":
                        msp,

                    "IRR raw":
                        irr_raw,

                    "IRR (%)":
                        irr_percent,

                    "IRR Solver Status":
                        irr_status,

                    "Product Ratio":
                        product_ratio,

                    **impact_results,

                    "Status":
                        "Successful",

                    "Error":
                        "",
                })

                # =================================================
                # CONSOLE OUTPUT
                # =================================================

                irr_text = (
                    f"{irr_percent:.2f}%"
                    if np.isfinite(irr_percent)
                    else f"NaN [{irr_status}]"
                )

                lhk_text = (
                    str(tuple(selected_lhk))
                    if selected_lhk
                    else "not logged"
                )

                print(
                    f"Sample {index}: "
                    f"SF Bio={biocrude_sf:.2f}% | "
                    f"SF Aq={aqueous_sf:.2f}% | "
                    f"SF Gas={gas_sf:.2f}% | "
                    f"LHK={lhk_text} | "
                    f"Lr={result['Selected Lr']:.3f} | "
                    f"Hr={result['Selected Hr']:.3f} | "
                    f"IRR={irr_text} | "
                    f"GlobalWarming="
                    f"{result[impact_column('GlobalWarming')]:.4g}"
                )

            except Exception as error:
                # Distillation diagnostics may have been partly
                # populated before a later TEA or LCA failure.
                result.update({
                    "Biobinder kg/hr": np.nan,
                    "Biofuel kg/hr": np.nan,
                    "Annual Biobinder kg/yr": np.nan,
                    "MSP ($/kg)": np.nan,
                    "IRR raw": np.nan,
                    "IRR (%)": np.nan,
                    "IRR Solver Status":
                        "simulation_failed",
                    "Product Ratio": np.nan,

                    **{
                        impact_column(cat): np.nan
                        for cat in TRACI_CATEGORIES
                    },

                    "Status": "Failed",
                    "Error": str(error),
                })

                print(
                    f"Sample {index} failed: {error}"
                )

            records.append(result)

    return pd.DataFrame(records)


def make_comparison(experimental_results, predicted_results, policy):
    """Create a wide experimental-versus-predicted comparison table."""
    metrics = [
        "MSP ($/kg)",
        "IRR (%)",
        "Product Ratio",
        *[impact_column(cat) for cat in TRACI_CATEGORIES],
    ]

    exp = experimental_results[
        ["Config", "Index", *metrics]
    ].rename(columns={
        col: f"Experimental {col}" for col in metrics
    })

    pred = predicted_results[
        ["Config", "Index", *metrics]
    ].rename(columns={
        col: f"Predicted {col}" for col in metrics
    })

    comparison = exp.merge(
        pred,
        on=["Config", "Index"],
        how="inner",
        validate="one_to_one",
    )

    comparison.insert(2, "Prediction Rule", policy)

    for category in TRACI_CATEGORIES:
        unit_col = impact_column(category)
        exp_col = f"Experimental {unit_col}"
        pred_col = f"Predicted {unit_col}"

        comparison[f"Residual {category}"] = (
            comparison[pred_col] - comparison[exp_col]
        )

    return comparison


# ============================================================
# LOAD DATA AND PREDICT YIELDS
# ============================================================

print(f"Loading RF model:\n{MODEL_FILE}")
rf_model = load(MODEL_FILE)

print(f"\nLoading experimental data:\n{DATA_FILE}")
df = pd.read_excel(DATA_FILE, sheet_name="Sheet1")
df.columns = df.columns.str.strip()

X = prepare_rf_inputs(df, rf_model)

predicted_yields = pd.DataFrame(
    rf_model.predict(X),
    columns=YIELD_COLS,
    index=df.index,
)

experimental_yields = make_experimental_yields(df)
prediction_cases = make_prediction_cases(predicted_yields)

original_valid, original_skipped = prediction_cases["original"]
normalized_valid, normalized_skipped = prediction_cases["normalized"]

print(f"\nTotal samples: {len(df)}")
print(
    f"Original 95–105 rule: "
    f"{len(original_valid)} retained, "
    f"{len(original_skipped)} skipped"
)
print(
    f"Normalized 90–110 rule: "
    f"{len(normalized_valid)} retained, "
    f"{len(normalized_skipped)} skipped"
)


# ============================================================
# RUN EXPERIMENTAL AND PREDICTED FULL-LCA CASES
# ============================================================

experimental_results = run_full_lca(
    experimental_yields,
    "Experimental yields",
)

original_results = run_full_lca(
    original_valid,
    "RF predicted: original 95–105 rule",
)

# normalized_results = run_full_lca(
#     normalized_valid,
#     "RF predicted: normalized 90–110 rule",
# )

original_comparison = make_comparison(
    experimental_results,
    original_results,
    "Original 95–105 wt%",
)

# normalized_comparison = make_comparison(
#     experimental_results,
#     normalized_results,
#     "Normalized 90–110 wt%",
# )


# ============================================================
# METADATA AND SAVE
# ============================================================

metadata = pd.DataFrame([
    {
        "Key": cat,
        "Alias": info.get("alias", cat),
        "Method": info.get("method", ""),
        "Category": info.get("category", ""),
        "Indicator Unit": info.get("unit", ""),
        "Reported Unit": (
            f"{info.get('unit', '')}/kg biobinder"
        ),
    }
    for cat, info in TRACI_CATEGORIES.items()
])

all_predictions = pd.concat(
    [
        predicted_yields,
        predicted_yields[YIELD_COLS]
        .sum(axis=1)
        .rename("Raw Sum wt%"),
    ],
    axis=1,
)
import re

_ILLEGAL_EXCEL_CHARS = re.compile(r"[\x00-\x08\x0B\x0C\x0E-\x1F]")

def clean_for_excel(df):
    df = df.copy()

    for col in df.select_dtypes(include=["object"]).columns:
        df[col] = df[col].map(
            lambda x: _ILLEGAL_EXCEL_CHARS.sub("", x)
            if isinstance(x, str) else x
        )

    return df

experimental_results = clean_for_excel(experimental_results)
original_results = clean_for_excel(original_results)
# normalized_results = clean_for_excel(normalized_results)

original_comparison = clean_for_excel(original_comparison)
# normalized_comparison = clean_for_excel(normalized_comparison)

original_skipped = clean_for_excel(original_skipped)
# normalized_skipped = clean_for_excel(normalized_skipped)

all_predictions = clean_for_excel(all_predictions)
metadata = clean_for_excel(metadata)

with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
    experimental_results.to_excel(
        writer, sheet_name="experimental_LCA", index=False
    )

    original_results.to_excel(
        writer, sheet_name="original_predicted_LCA", index=False
    )
    original_comparison.to_excel(
        writer, sheet_name="original_comparison", index=False
    )
    original_skipped.to_excel(
        writer, sheet_name="original_skipped", index=True
    )

    # normalized_results.to_excel(
    #     writer, sheet_name="normalized_pred_LCA", index=False
    # )
    # normalized_comparison.to_excel(
    #     writer, sheet_name="normalized_comparison", index=False
    # )
    # normalized_skipped.to_excel(
    #     writer, sheet_name="normalized_skipped", index=True
    # )

    all_predictions.to_excel(
        writer, sheet_name="all_RF_predictions", index=True
    )
    metadata.to_excel(
        writer, sheet_name="TRACI_metadata", index=False
    )

print(f"\nFull RF–LCA analysis saved to:\n{OUTPUT_FILE}")