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
import sys
import importlib.util
import pandas as pd
import qsdsan as qs
from datetime import datetime
# ============================================================
# PATHS
# ============================================================
DISTILLATION_POLICY = "min_err"  # max_top_ratio, min_err, max_Hr_then_max_Lr, max_Lr_then_max_Hr
timestamp = datetime.now().strftime("%Y%m%d_%H%M")
dist_flex_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\Billion Ton\Dist_flex.py"

input_excel = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\Billion Ton\county_waste_scenarios_comparison_near_term.xlsx"

output_excel = rf"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\Billion Ton\county_biobinder_potential_food_green_manure_htl_{DISTILLATION_POLICY}_{timestamp}.xlsx"

os.makedirs(os.path.dirname(output_excel), exist_ok=True)

# ============================================================
# IMPORT DIST_FLEX
# ============================================================

dist_dir = os.path.dirname(dist_flex_path)
if dist_dir not in sys.path:
    sys.path.insert(0, dist_dir)

spec = importlib.util.spec_from_file_location("Dist_flex", dist_flex_path)
Dist_flex = importlib.util.module_from_spec(spec)
Dist_flex.DISTILLATION_POLICY = DISTILLATION_POLICY
spec.loader.exec_module(Dist_flex)

create_system = Dist_flex.create_system

# ============================================================
# SETTINGS
# ============================================================

feedstocks = ["food", "green", "manure"]

config_kwargs = dict(
    decentralized_HTL=True,
    decentralized_upgrading=False,
    skip_EC=True,
    generate_H2=False,
    EC_config=None,
)

# ============================================================
# HELPERS
# ============================================================

def get_stream(fs, names):
    for name in names:
        if hasattr(fs.stream, name):
            return getattr(fs.stream, name)
    return None


def safe_water_mass(stream):
    try:
        return stream.imass["Water"]
    except Exception:
        return 0.0


def annual_kg(stream, sys_obj):
    return stream.F_mass * sys_obj.operating_hours


def print_stream_diagnostic(fs, sys_obj, feedstock):
    print("\n" + "=" * 80)
    print(f"DIAGNOSTIC FOR FEEDSTOCK: {feedstock}")
    print("=" * 80)

    print(f"Operating hours: {sys_obj.operating_hours:,.0f} hr/yr")

    feed = get_stream(fs, ["scaled_feedstock", "feedstock", "Feedstock", "HTL_feed"])
    if feed is None:
        raise RuntimeError("Could not find feed stream.")

    water_kg_hr = safe_water_mass(feed)
    dry_feed_kg_hr = feed.F_mass - water_kg_hr
    dry_feed_kg_yr = dry_feed_kg_hr * sys_obj.operating_hours

    print("\n--- FEED ---")
    print(f"Feed stream ID: {feed.ID}")
    print(f"Total feed: {feed.F_mass:,.3f} kg/hr")
    print(f"Water:      {water_kg_hr:,.3f} kg/hr")
    print(f"Dry feed:   {dry_feed_kg_hr:,.3f} kg/hr")
    print(f"Dry feed:   {dry_feed_kg_yr:,.3f} kg/yr")
    print(f"Dry feed:   {dry_feed_kg_yr / 1000:,.3f} metric ton/yr")
    print(f"Dry feed:   {dry_feed_kg_yr / 1000 / 365:,.3f} metric ton/day")

    stream_groups = {
        "biocrude": [
            "scaled_biocrude", "biocrude_trans_surrogate",
            "splitted_crude", "crude_to_dist",
            "transported_biocrude", "HTLbiocrude", "biocrude"
        ],
        "biobinder": [
            "biobinder", "hot_biobinder", "cooled_biobinder",
            "Biobinder", "bio_binder"
        ],
        "biofuel": [
            "biofuel", "hot_biofuel", "cooled_biofuel",
            "biofuel_additive", "fuel", "mixed_fuel"
        ],
    }

    found = {}

    print("\n--- MAIN PRODUCT STREAMS ---")
    for label, names in stream_groups.items():
        s = get_stream(fs, names)
        found[label] = s

        if s is None:
            print(f"{label}: NOT FOUND")
            continue

        s_kg_yr = annual_kg(s, sys_obj)

        print(f"{label}: {s.ID}")
        print(f"  {s.F_mass:,.3f} kg/hr")
        print(f"  {s_kg_yr:,.3f} kg/yr")
        print(f"  {s_kg_yr / 1000:,.3f} metric ton/yr")

        if dry_feed_kg_yr > 0:
            print(f"  yield vs dry feed = {s_kg_yr / dry_feed_kg_yr:.6f} kg/kg dry feed")

    # ---- additional diagnostics ----
    print("\n--- HTL YIELD + DISTILLATION DIAGNOSTICS ---")

    htl_yields_used = dict(Dist_flex.HTL_yields)
    print("HTL_yields used:")
    for k, v in htl_yields_used.items():
        print(f"  {k}: {v}")

    biocrude = found["biocrude"]
    biobinder = found["biobinder"]
    biofuel = found["biofuel"]

    biocrude_kg_yr = annual_kg(biocrude, sys_obj) if biocrude else None
    biobinder_kg_yr = annual_kg(biobinder, sys_obj) if biobinder else None
    biofuel_kg_yr = annual_kg(biofuel, sys_obj) if biofuel else None

    if biocrude_kg_yr and dry_feed_kg_yr > 0:
        print(f"Stream-based biocrude yield = {biocrude_kg_yr / dry_feed_kg_yr:.6f} kg/kg dry feed")

    if biobinder_kg_yr and biocrude_kg_yr:
        print(f"Biobinder / scaled biocrude = {biobinder_kg_yr / biocrude_kg_yr:.6f}")

    if biofuel_kg_yr and biocrude_kg_yr:
        print(f"Biofuel / scaled biocrude = {biofuel_kg_yr / biocrude_kg_yr:.6f}")

    if biobinder_kg_yr and biofuel_kg_yr:
        product_sum = biobinder_kg_yr + biofuel_kg_yr
        print(f"Biobinder fraction of biofuel+biobinder = {biobinder_kg_yr / product_sum:.6f}")
        print(f"Biofuel fraction of biofuel+biobinder = {biofuel_kg_yr / product_sum:.6f}")

    if hasattr(fs.unit, "CrudeHeavyDis"):
        col = fs.unit.CrudeHeavyDis
        print(f"CrudeHeavyDis Lr = {col.Lr}")
        print(f"CrudeHeavyDis Hr = {col.Hr}")
        print(f"CrudeHeavyDis LHK = {col.LHK}")

    print("\n--- STREAM IDS CONTAINING crude / binder / fuel / bio ---")
    for s in fs.stream:
        sid = s.ID.lower()
        if any(k in sid for k in ["crude", "binder", "fuel", "bio"]):
            print(
                f"{s.ID:40s}"
                f"{s.F_mass:15,.3f} kg/hr"
                f"{s.F_mass * sys_obj.operating_hours:18,.3f} kg/yr"
            )

    print("=" * 80 + "\n")

    return {
        "feed": feed,
        "dry_feed_kg_hr": dry_feed_kg_hr,
        "dry_feed_kg_yr": dry_feed_kg_yr,
        "HTL_yields_used": htl_yields_used,
        "Y_biocrude_input": htl_yields_used.get("biocrude"),
        "biocrude": found["biocrude"],
        "biobinder": found["biobinder"],
        "biofuel": found["biofuel"],
    }


# ============================================================
# GET SIMULATED CONVERSION FACTORS
# ============================================================

conversion_factors = {}
diagnostic_rows = []

print("\nSimulating systems to extract biobinder conversion factors...\n")

for feedstock in feedstocks:
    print(f"Running {feedstock}...")

    try:
        # Uncomment only if needed:
        # qs.clear_lca_registries()

        sys_obj = create_system(feedstock_id=feedstock, **config_kwargs)
        sys_obj.simulate()

        fs = sys_obj.flowsheet

        diag = print_stream_diagnostic(fs, sys_obj, feedstock)

        binder = diag["biobinder"]
        if binder is None:
            raise RuntimeError("Could not find biobinder stream.")

        dry_feed_kg_yr = diag["dry_feed_kg_yr"]
        binder_kg_yr = binder.F_mass * sys_obj.operating_hours

        if dry_feed_kg_yr <= 0:
            raise ValueError("Annual dry feed flow is zero or negative.")

        factor = binder_kg_yr / dry_feed_kg_yr
        conversion_factors[feedstock] = factor

        biocrude_kg_yr = (
            diag["biocrude"].F_mass * sys_obj.operating_hours
            if diag["biocrude"] is not None else None
        )

        biofuel_kg_yr = (
            diag["biofuel"].F_mass * sys_obj.operating_hours
            if diag["biofuel"] is not None else None
        )

        diagnostic_rows.append({
            "Feedstock": feedstock,
            "Operating_hours_per_year": sys_obj.operating_hours,
            "Dry_feed_kg_per_year": dry_feed_kg_yr,
            "Dry_feed_metric_ton_per_year": dry_feed_kg_yr / 1000,
            "Dry_feed_metric_ton_per_day": dry_feed_kg_yr / 1000 / 365,
            "Biocrude_kg_per_year": biocrude_kg_yr,
            "Biobinder_kg_per_year": binder_kg_yr,
            "Biofuel_kg_per_year": biofuel_kg_yr,
            "Biobinder_kg_per_kg_dry_feed": factor,
            "Distillation_policy": DISTILLATION_POLICY,
        })

        print(f"FINAL FACTOR for {feedstock}: {factor:.6f} kg binder/kg dry feed\n")

    except Exception as e:
        print(f"ERROR for {feedstock}: {e}")
        raise

print("Final simulated conversion factors:")
print(conversion_factors)

# ============================================================
# APPLY TO COUNTY WASTE SCENARIOS
# ============================================================

def process_sheet(sheet_name, suffix):
    df = pd.read_excel(input_excel, sheet_name=sheet_name)

    base_cols = ["fips", "county", "state"]

    food_col = f"Food Wastes{suffix}_Dry_Tons_Year"
    green_col = f"Green, Yard, Wood & Agricultural Wastes{suffix}_Dry_Tons_Year"
    manure_col = f"Manure Wastes{suffix}_Dry_Tons_Year"

    required = base_cols + [food_col, green_col, manure_col]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(
            f"Missing columns in sheet '{sheet_name}':\n{missing}\n\n"
            f"Available columns:\n{list(df.columns)}"
        )

    out = df[required].copy()

    out = out.rename(columns={
        food_col: "Food_Dry_Tons_Year",
        green_col: "Green_Dry_Tons_Year",
        manure_col: "Manure_Dry_Tons_Year",
    })

    out["Food_Biobinder_Dry_Tons_Year"] = (
        out["Food_Dry_Tons_Year"] * conversion_factors["food"]
    )

    out["Green_Biobinder_Dry_Tons_Year"] = (
        out["Green_Dry_Tons_Year"] * conversion_factors["green"]
    )

    out["Manure_Biobinder_Dry_Tons_Year"] = (
        out["Manure_Dry_Tons_Year"] * conversion_factors["manure"]
    )

    out["Total_Biobinder_Dry_Tons_Year"] = (
        out["Food_Biobinder_Dry_Tons_Year"]
        + out["Green_Biobinder_Dry_Tons_Year"]
        + out["Manure_Biobinder_Dry_Tons_Year"]
    )

    return out


print("\nApplying factors to county dry-ton waste data...\n")

all_waste = process_sheet(
    sheet_name="All Waste Available",
    suffix=""
)

landfilled = process_sheet(
    sheet_name="Only Landfilled Portion",
    suffix="_Landfilled"
)

# ============================================================
# EXPORT
# ============================================================

run_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:

    # --------------------------------------------------------
    # Conversion factors
    # --------------------------------------------------------
    pd.DataFrame({
        "Feedstock": list(conversion_factors.keys()),
        "Biobinder_DryTon_per_DryTonFeed": list(conversion_factors.values()),
    }).to_excel(
        writer,
        sheet_name="Conversion Factors",
        index=False,
    )

    # --------------------------------------------------------
    # Simulation diagnostics
    # --------------------------------------------------------
    pd.DataFrame(diagnostic_rows).to_excel(
        writer,
        sheet_name="Simulation Diagnostics",
        index=False,
    )

    # --------------------------------------------------------
    # Run metadata
    # --------------------------------------------------------
    pd.DataFrame({
        "Parameter": [
            "Run Timestamp",
            "Distillation Policy",
            "Feedstocks",
            "HTL Configuration",
        ],
        "Value": [
            run_timestamp,
            DISTILLATION_POLICY,
            ", ".join(feedstocks),
            "DHCU | No EC | RF HTL yields",
        ],
    }).to_excel(
        writer,
        sheet_name="Run Metadata",
        index=False,
    )

    # --------------------------------------------------------
    # County outputs
    # --------------------------------------------------------
    all_waste.to_excel(
        writer,
        sheet_name="All Waste Available",
        index=False,
    )

    landfilled.to_excel(
        writer,
        sheet_name="Only Landfilled Portion",
        index=False,
    )

print("\nDone.")
print(f"Saved to:\n{output_excel}")