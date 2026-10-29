# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Buswell methane potential with the Labatut et al. biodegradability correction."""
from math import isfinite

SOURCE = 'https://doi.org/10.1016/j.biortech.2010.10.035'
# Labatut et al. Table of molecular formulae; kmol atoms per kmol compound.
# Fiber repeat-unit approximations extend the lumped Cornell characterization.
# Lignin formula: https://github.com/nomenazo/PaCTEA/blob/main/rcomposting.m
# Its theoretical COD potential is included BEFORE the whole-feed measured
# COD correction; excluding it would count recalcitrance twice. This is an
# aggregate gas balance, not a claim of biological lignin conversion.
FORMULAS = dict(Carbohydrates=(6, 12, 6, 0), Proteins=(16, 24, 5, 4),
                Lipids=(8, 16, 1, 0), Cellulose=(6, 10, 5, 0),
                Hemicellulose=(5, 8, 4, 0), Lignin=(20,22,7,0))
ATOMIC = (12.011, 1.00794, 15.9994, 14.0067)


def predict(composition, biodegradability, reactor_conversion):
    """Return yields per kg wet feed, gas at 0 C and 1 atm.

    biodegradability = degradable COD/total COD, determined independently of
    this prediction (Labatut et al. equation 2). reactor_conversion is the
    fraction of corrected BMP achieved by the specified operating reactor.
    This function accepts explicit corrections. The system factory resolves
    named-feed literature proxies and reports their provenance separately.
    """
    for name, value in (('ad_biodegradability', biodegradability),
                        ('reactor_conversion', reactor_conversion)):
        if value is None or not isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f'Supply {name} in [0, 1] with an independent measurement or cited scenario')
    ch4 = co2 = nh3 = water = converted = carbon = 0.
    correction = biodegradability*reactor_conversion
    for key, mass in composition.items():
        if key in ('Water', 'Ash'):
            continue
        formula = FORMULAS.get(key)
        if formula is None:
            raise ValueError(f'Unsupported AD biochemical fraction: {key}')
        a,b,c,d = formula
        mw = sum(x*y for x,y in zip(formula, ATOMIC))
        carbon += mass*a*ATOMIC[0]/mw
        mol = mass/mw
        ch4 += mol*(a/2+b/8-c/4-3*d/8)
        co2 += mol*(a/2-b/8+c/4+3*d/8)
        nh3 += mol*d
        water += mol*(a-b/4-c/2+3*d/4)*18.01528
        converted += mass
    return dict(theoretical_methane_nm3_kg_wet=ch4*22.414,
                methane_nm3_kg_wet=ch4*22.414*correction,
                co2_nm3_kg_wet=co2*22.414*correction,
                ammonia_kg_kg_wet=nh3*17.03052*correction,
                reaction_water_kg_kg_wet=water*correction,
                converted_organics_kg_kg_wet=converted*correction,
                carbon_kg_kg_wet=carbon, biodegradability=biodegradability,
                reactor_conversion=reactor_conversion, source=SOURCE)
