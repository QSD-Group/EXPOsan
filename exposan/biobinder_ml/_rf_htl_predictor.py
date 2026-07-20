# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from exposan.biobinder_ml._ml_features import (
    make_feature_row,
    assert_model_feature_compatibility,
)
from typing import Dict
import numpy as np


def predict_htl_yields_from_rf(
    *,
    rf_model,
    carb_wt: float,
    protein_wt: float,
    lipids_wt: float,
    ash_wt: float,
    process: dict,
) -> Dict[str, float]:
    """
    Predict HTL yields (dry-weight fractions) using RF model.

    Returns
    -------
    dict with keys:
      - 'biocrude'
      - 'aqueous'
      - 'gas'
      - 'char'
    """

    # Build feature row
    X = make_feature_row(
        carb_wt=carb_wt,
        protein_wt=protein_wt,
        lipids_wt=lipids_wt,
        ash_wt=ash_wt,
        process=process,
    )

    # Optional guard
    assert_model_feature_compatibility(rf_model)

    # Predict (wt%)
    y = rf_model.predict(X)[0]

    # Map to HTL yields (convert to fractions)
    yields = {
        "biocrude": float(y[0]) / 100.0,
        "aqueous":  float(y[1]) / 100.0,
        "gas":      float(y[2]) / 100.0,
        "char":     float(y[3]) / 100.0,
    }

    # Defensive normalization
    s = sum(yields.values())
    if s <= 0:
        raise ValueError("RF predicted non-physical HTL yields.")

    for k in yields:
        yields[k] /= s

    return yields