
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Ali Ahmad <aa3056@scarletmail.rutgers.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
__all__ = (
    'BIOBINDER_FEEDSTOCKS',
    'SUPPORTED_FEEDSTOCKS',
    'get_feedstock_composition',
)

# ---------------------------------------------------------------------
# Canonical feedstock definitions (baseline values from exposan.htl)
# ---------------------------------------------------------------------


BIOBINDER_FEEDSTOCKS = {

    'sludge': dict(
        moisture=0.70,
        ash_dw=0.257,
        lipid_afdw=0.204,
        protein_afdw=0.463,
    ),

    'food': dict(
        moisture=0.74,
        ash_dw=0.0679,
        lipid_afdw=0.22,
        protein_afdw=0.20,
    ),

    'fog': dict(
        moisture=0.35,
        ash_dw=0.01865,
        lipid_afdw=0.987,
        protein_afdw=0.002,
    ),

    'green': dict(
        moisture=0.342,
        ash_dw=0.134,
        lipid_afdw=0.018,
        protein_afdw=0.049,
    ),

    'manure': dict(
        moisture=0.6634,
        ash_dw=0.3056,
        lipid_afdw=0.092325,
        protein_afdw=0.216375,
    ),
}

SUPPORTED_FEEDSTOCKS = tuple(BIOBINDER_FEEDSTOCKS.keys())




def get_feedstock_composition(feedstock: str):
    """
    Return feedstock composition formatted for the biobinder
    Conditioning unit.

    Parameters
    ----------
    feedstock : str
        One of: 'sludge', 'food', 'fog', 'green', 'manure'

    Returns
    -------
    As-received feedstock composition (dry solids + inherent moisture), before conditioning
    """
    if feedstock not in BIOBINDER_FEEDSTOCKS:
        raise ValueError(
            f"Unknown feedstock '{feedstock}'. "
            f"Choose from {SUPPORTED_FEEDSTOCKS}."
        )

    fs = BIOBINDER_FEEDSTOCKS[feedstock]

    moisture = fs['moisture']
    ash_dw = fs['ash_dw']
    lipid_afdw = fs['lipid_afdw']
    protein_afdw = fs['protein_afdw']

    # Sanity checks (fail fast)
    if not (0 < moisture < 1):
        raise ValueError(f"Invalid moisture for '{feedstock}': {moisture}")

    if lipid_afdw + protein_afdw > 1:
        raise ValueError(
            f"Lipid + protein exceed AFDW for '{feedstock}'."
        )

    # Dry mass fractions
    dry_frac = 1 - moisture
    ash = dry_frac * ash_dw

    lipid = dry_frac * lipid_afdw
    protein = dry_frac * protein_afdw
    carbohydrate = dry_frac * (1 - lipid_afdw - protein_afdw)

    # Final composition (wet basis)
    composition = {
        'Water': moisture,
        'Lipids': lipid,
        'Proteins': protein,
        'Carbohydrates': carbohydrate,
        'Ash': ash,
    }

    # Optional: enforce closure (floating point safety)
    total = sum(composition.values())
    if abs(total - 1) > 1e-6:
        composition['Water'] += (1 - total)

    return composition



