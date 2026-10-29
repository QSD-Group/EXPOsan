
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
    'get_recovery_composition',
    'AD_BIODEGRADABILITY_DEFAULTS',
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

# https://doi.org/10.1016/j.biortech.2010.10.035 Table 4, g methane-COD/g feed-COD.
# Named-feed defaults are literature proxies, not measurements of these feeds.
AD_BIODEGRADABILITY_DEFAULTS = {
    'food': dict(value=.60, substrate='Potatoes; proxy for mixed food waste',
                 source='https://doi.org/10.1016/j.biortech.2010.10.035',
                 basis='Table 4, potatoes, methane-COD/total COD'),
    'fog': dict(value=(.44+.46)/2, substrate='Suspended and settled FOG midpoint',
                source='https://doi.org/10.1016/j.biortech.2010.10.035',
                basis='Table 4 midpoint of 0.44 and 0.46, methane-COD/total COD'),
    'green': dict(value=.24, substrate='Switchgrass; proxy for green waste',
                  source='https://doi.org/10.1016/j.biortech.2010.10.035',
                  basis='Table 4, switchgrass, methane-COD/total COD'),
    'manure': dict(value=.55, substrate='Raw dairy manure',
                   source='https://doi.org/10.1016/j.biortech.2010.10.035',
                   basis='Table 4, dairy manure, methane-COD/total COD'),
    # https://doi.org/10.1016/j.watres.2013.07.019 Mixed municipal sludge:
    # 188-214 mL CH4/g CODfed; range midpoint / 350 mL CH4/g methane-COD.
    # This is methane-COD recovery, not the reported 58-65% COD removal.
    'sludge': dict(value=((188+214)/2)/350,
                   substrate='Seven mixed municipal sewage sludges; range-midpoint proxy',
                   source='https://doi.org/10.1016/j.watres.2013.07.019',
                   basis='(188+214)/2/350, methane-COD/total COD'),
}
for _name, _data in AD_BIODEGRADABILITY_DEFAULTS.items():
    BIOBINDER_FEEDSTOCKS[_name].update(
        ad_biodegradability=_data['value'],
        ad_biodegradability_source=_data['source'],
        ad_biodegradability_substrate=_data['substrate'],
        ad_biodegradability_basis=_data['basis'],
        ad_biodegradability_status='literature proxy',
    )
del _name, _data

# Additional measurements can be stored in each BIOBINDER_FEEDSTOCKS record:
# cellulose_afdw, hemicellulose_afdw, lignin_afdw (kg/kg ash-free dry matter),
# ad_biodegradability (degradable COD/total COD), ad_reactor_conversion
# (achieved methane/corrected BMP at specified operation), characterization_source.
# AD literature proxies above remain overridable by feed-specific measurements.


def get_recovery_composition(feedstock, *, characterization=None, resolve_fibers=False):
    """Return the HTL wet-basis composition, fractionated for recovery models.

    Preserve the existing HTL mapping: lipid/protein fractions multiply the
    nominal dry fraction; ash displaces water when closing the wet balance.
    This historical convention differs from a strict analytical AFDW mapping.
    All technologies must use the same effective composition for comparison.
    Carbohydrates include cellulose, hemicellulose and lignin. Optional fiber
    details are resolved only when resolve_fibers=True; they partition the
    umbrella total and are never added on top of it.
    ``characterization`` overrides a copy, never the canonical feedstock data.
    """
    from math import isfinite
    if isinstance(feedstock, str):
        if feedstock not in BIOBINDER_FEEDSTOCKS:
            raise ValueError(f'Unknown feedstock {feedstock!r}; choose {SUPPORTED_FEEDSTOCKS}')
        fs = dict(BIOBINDER_FEEDSTOCKS[feedstock])
    else:
        fs = dict(feedstock)
    fs.update(characterization or {})
    for key in ('moisture', 'ash_dw', 'lipid_afdw', 'protein_afdw'):
        value = fs[key]
        if not isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f'{key} must be a finite fraction')
    organic = 1-fs['moisture']
    remainder = 1-fs['lipid_afdw']-fs['protein_afdw']
    if remainder < -1e-12:
        raise ValueError('Protein and lipid fractions exceed ash-free dry matter')
    fiber_keys = ('cellulose_afdw', 'hemicellulose_afdw', 'lignin_afdw')
    fibers = [fs.get(k) for k in fiber_keys]
    if any(v is not None for v in fibers) and not all(v is not None for v in fibers):
        raise ValueError('Supply cellulose_afdw, hemicellulose_afdw and lignin_afdw together')
    if all(v is not None for v in fibers):
        if any(not isfinite(v) or v < 0 for v in fibers) or sum(fibers) > remainder+1e-12:
            raise ValueError('Fiber fractions must be nonnegative and fit the AFDW remainder')
    ash = (1-fs['moisture'])*fs['ash_dw']
    result = dict(Water=fs['moisture']-ash, Ash=ash,
                  Lipids=organic*fs['lipid_afdw'], Proteins=organic*fs['protein_afdw'],
                  Carbohydrates=organic*max(remainder, 0.))
    if result['Water'] < 0:
        raise ValueError('The legacy HTL composition mapping produces negative water for this characterization')
    if resolve_fibers and all(v is not None for v in fibers):
        result.update(Cellulose=organic*fibers[0], Hemicellulose=organic*fibers[1],
                      Lignin=organic*fibers[2])
        result['Carbohydrates'] = max(result['Carbohydrates']-organic*sum(fibers),0.)
    return result




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



