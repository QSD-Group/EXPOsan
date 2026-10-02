# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 13:52:33 2025

@author: aliah
"""

# ============================================================
# 0. Imports
# ============================================================

from joblib import load
import qsdsan as qs

from exposan.biobinder_ml.feedstocks import get_feedstock_composition
from exposan.biobinder_ml.rf_yields import predict_htl_yields_from_rf
from exposan.biobinder_ml.features import assert_model_feature_compatibility

from exposan.biobinder_ml import (
    _load_components,
    _load_process_settings,
    _units as u,
    BiobinderTEA,
    central_dry_flowrate as default_central,
    pilot_dry_flowrate as default_pilot,
    price_dct,
)

from exposan.biobinder_ml.Dist_flex import create_system


# ============================================================
# 1. Load RF model ONCE
# ============================================================

rf_model = load("results/rf_yield_model.joblib")
assert_model_feature_compatibility(rf_model)
print("✅ RF model loaded")


# ============================================================
# 2. Choose feedstock
# ============================================================

feedstock_id = "fog"   # 'sludge', 'food', 'fog', 'green', 'manure'

wet_comp = get_feedstock_composition(feedstock_id)

print(f"Feedstock used is {feedstock_id}")
print("Corresponding composition (wet basis):")
for k, v in wet_comp.items():
    print(f"  {k}: {v:.2f}")


# ============================================================
# 3. Convert wet → dry basis (FOR RF ONLY)
# ============================================================

dry_frac = 1.0 - wet_comp["Water"]

carb_wt    = wet_comp["Carbohydrates"] / dry_frac * 100.0
protein_wt = wet_comp["Proteins"]      / dry_frac * 100.0
lipids_wt  = wet_comp["Lipids"]        / dry_frac * 100.0
ash_wt     = wet_comp["Ash"]           / dry_frac * 100.0


# ============================================================
# 4. Predict HTL yields using RF (dry-basis)
# ============================================================

HTL_yields = predict_htl_yields_from_rf(
    rf_model=rf_model,
    carb_wt=carb_wt,
    protein_wt=protein_wt,
    lipids_wt=lipids_wt,
    ash_wt=ash_wt,
    process={
        "Temperature (C)": 280,
        "Residence Time": 15,          # minutes
        "Solid content (w/w) %": 20,   # matches Conditioning
    },
)

print("RF-predicted HTL yields:")
for k, v in HTL_yields.items():
    print(f"  {k}: {v:.3f}")


# ============================================================
# 5. Create system (UNCHANGED structure)
# ============================================================

sys = create_system(
    flowsheet=None,
    central_dry_flowrate=default_central,
    decentralized_HTL=False,
    decentralized_upgrading=False,
    skip_EC=True,
    generate_H2=False,
    EC_config=None,
    # HTL_yields picked up inside create_system
)


# ============================================================
# 6. Run simulation
# ============================================================

sys.simulate()

tea = sys.TEA
lca = sys.LCA

biobinder = sys.flowsheet.stream.biobinder
biofuel   = sys.flowsheet.stream.biofuel


# ============================================================
# 7. Results
# ============================================================

MSP = tea.solve_price(biobinder)
IRR = tea.solve_IRR() * 100

GWP = lca.get_allocated_impacts(
    streams=(biobinder,),
    operation_only=True,
    annual=True,
)["GWP"] / (biobinder.F_mass * lca.system.operating_hours)

ratio = biobinder.F_mass / biofuel.F_mass if biofuel.F_mass else float("nan")

print(f"\nInternal rate of return of the system is {IRR:.2f}%")
print(f"OPEX intensity: ${tea.annual_operating_cost / sys.feedstock.F_mass:.2f} per ton-feed")
print(f"CAPEX intensity: ${tea.installed_equipment_cost / sys.feedstock.F_mass:.2f} per ton-feed")
print(f"Global warming potential of the biobinder is {GWP:.4f} kg CO2e/kg.")
print(f"Product ratio (biobinder/biofuel): {ratio:.3f}")
