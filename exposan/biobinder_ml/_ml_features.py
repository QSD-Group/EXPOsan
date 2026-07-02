
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Ali Ahmad <aa3056@scarletmail.rutgers.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Any, List

import numpy as np
import pandas as pd


# Feng et al.
LIPID_2_C   = 0.750
PROTEIN_2_C = 0.545
CARBO_2_C   = 0.400

LIPID_2_H   = 0.125
PROTEIN_2_H = 0.0682  
CARBO_2_H   = 0.0667

PROTEIN_2_N = 0.159



FEATURE_COLS: List[str] = [
    "Pre-processing",
    "Carbohydrates wt%",
    "Protein wt%",
    "Lipids wt%",
    "Ash wt%",
    "C%",
    "H%",
    "O%",
    "N%",
    "HHV Biomass",
    "Catalyst",
    "Reactor Type",
    "Reactor Volume (mL)",
    "Solid content (w/w) %",
    "Residence Time",
    "Temperature (C)",
    "Solvent",
]


PREPROCESSING_MAP = {
    0: "No",
    1: "Yes",
}

REACTOR_TYPE_MAP = {
    0: "Batch / Autoclave",
    1: "Plug Flow Reactor (PFR)",
    2: "CSTR",
}

CATALYST_MAP = {
    0: "Unknown / Not specified",
    1: "No catalyst",
    2: "K2CO3 (5 wt%)",
    3: "Na2CO3",
    4: "KOH (0.06 g)",
    5: "HZSM-5",
    6: "ZSM-5",
    7: "Ca(OH)2 (1 wt%)",
    8: "Fe, Ni",
    9: "KOH",
    10: "K2CO3 (1 wt%)",
    11: "KOH (1 wt%)",
    12: "Fe",
    13: "Na2CO3 (1 wt%)",
    14: "Ni",
    15: "Fe,Ni/Al2O3",
    16: "Ni/Al2O3",
    17: "Bentonite clay",
    18: "Alumina",
    19: "Kaolinite",
    20: "Talc",
    21: "Montmorillonite",
    22: "Attapulgite",
    23: "Phlogopite",
    24: "Vermiculite",
    25: "NiMo-20",
    26: "Diatomaceous earth",
    27: "NiMo-15",
    28: "NiMo-10",
    29: "NiMo-0",
    30: "H2SO4 + O2",
    31: "H2SO4 + H2",
    32: "H2SO4",
    33: "KOH + O2",
    34: "KOH + H2",
    35: "Ni/SiO2-Al2O3",
    36: "Ce/HZSM-5",
    37: "Pd/C",
    38: "ZnO",
    39: "ZMS-5",
    40: "Meixnerite",
}

SOLVENT_MAP = {
    0: "Unknown / Not specified",
    1: "No solvent",
    2: "50% ethanol",
    3: "Seawater",
}

@dataclass(frozen=True)
class FeatureDefaults:
    """defaults from biobinder module"""
    pre_processing: int = 1
    catalyst: int = 1
    reactor_type: int = 1
    solvent: int = 1
    reactor_volume_ml: float = 10459.7*1000 # pilot scale 28.88 L x (4157.93/11.46)
    solid_content_w_w_pct: float = 20
    residence_time: float = 30.0
    temperature_c: float = 280 # °C
    #if given
    hhv_biomass: Optional[float] = None


def compute_CHON_from_biochem_dry_wt_pct(
    carb_wt: float, protein_wt: float, lipids_wt: float, ash_wt: float
) -> Dict[str, float]:
    """
    Input: dry-basis wt% including ash (carb+protein+lipids+ash = 100)
    Output: C/H/O/N as wt% (dry basis), with O by difference incl. ash.
    """
    C = CARBO_2_C * carb_wt + PROTEIN_2_C * protein_wt + LIPID_2_C * lipids_wt
    H = CARBO_2_H * carb_wt + PROTEIN_2_H * protein_wt + LIPID_2_H * lipids_wt
    N = PROTEIN_2_N * protein_wt
    O = 100.0 - C - H - N - ash_wt
    return {"C%": C, "H%": H, "N%": N, "O%": O}


def compute_HHV_from_CHON_wt_pct(C, H, O, N):
    """
    Higher Heating Value (HHV) in MJ/kg
    C, H, O, N in wt% (dry basis)

    HHV = 0.3516*C + 1.16225*H - 0.1109*O + 0.0628*N
    """
    return (
        0.3516 * C
        + 1.16225 * H
        - 0.1109 * O
        + 0.0628 * N
    )

def make_feature_row(
    *,
    carb_wt: float,
    protein_wt: float,
    lipids_wt: float,
    ash_wt: float,
    process: Optional[Dict[str, Any]] = None,
    defaults: FeatureDefaults = FeatureDefaults(),
) -> pd.DataFrame:
    """
    Build a single-row DataFrame with FEATURE_COLS.

    `process` can include any of:
      - "Pre-processing", "Catalyst", "Reactor Type", "Solvent"  (int codes)
      - "Reactor Volume (mL)" (float)
      - "Solid content (w/w) %" (float)
      - "Residence Time" (float)
      - "Temperature (C)" (float)
      - "HHV Biomass" (float)  [optional override]
    """
    process = process or {}

    #check
    s = carb_wt + protein_wt + lipids_wt + ash_wt
    if not (95.0 <= s <= 105.0):
        raise ValueError(
            f"Dry composition should sum ~100 wt%. Got {s:.3f} "
            f"(carb={carb_wt}, protein={protein_wt}, lipids={lipids_wt}, ash={ash_wt})."
        )

    chon = compute_CHON_from_biochem_dry_wt_pct(carb_wt, protein_wt, lipids_wt, ash_wt)
    hhv = process.get("HHV Biomass", None)
    if hhv is None:
        if defaults.hhv_biomass is not None:
            hhv = float(defaults.hhv_biomass)
        else:
            hhv = compute_HHV_from_CHON_wt_pct(
            chon["C%"],
            chon["H%"],
            chon["O%"],
            chon["N%"],
           )

    row = {
        "Pre-processing": int(process.get("Pre-processing", defaults.pre_processing)),
        "Carbohydrates wt%": float(carb_wt),
        "Protein wt%": float(protein_wt),
        "Lipids wt%": float(lipids_wt),
        "Ash wt%": float(ash_wt),
        "C%": float(chon["C%"]),
        "H%": float(chon["H%"]),
        "O%": float(chon["O%"]),
        "N%": float(chon["N%"]),
        "HHV Biomass": float(hhv),
        "Catalyst": int(process.get("Catalyst", defaults.catalyst)),
        "Reactor Type": int(process.get("Reactor Type", defaults.reactor_type)),
        "Reactor Volume (mL)": float(process.get("Reactor Volume (mL)", defaults.reactor_volume_ml)),
        "Solid content (w/w) %": float(process.get("Solid content (w/w) %", defaults.solid_content_w_w_pct)),
        "Residence Time": float(process.get("Residence Time", defaults.residence_time)),
        "Temperature (C)": float(process.get("Temperature (C)", defaults.temperature_c)),
        "Solvent": int(process.get("Solvent", defaults.solvent)),
    }

    X = pd.DataFrame([row], columns=FEATURE_COLS)
    return X


def assert_model_feature_compatibility(rf_model) -> None:
    """
    Ensure the trained model expects the same feature columns/order.
    """
    if hasattr(rf_model, "feature_names_in_"):
        expected = list(rf_model.feature_names_in_)
        if expected != FEATURE_COLS:
            raise ValueError(
                "Feature mismatch vs rf_model.feature_names_in_.\n"
                f"Model expects: {expected}\n"
                f"Adapter provides: {FEATURE_COLS}\n"
            )
