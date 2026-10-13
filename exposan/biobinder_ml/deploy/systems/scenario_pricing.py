# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import math
PNNL_USD_GGE = 2.23
# Historical USGC FOB assessment 2022-05-16, not a current quotation.
VLSFO_USD_TONNE = 856.75
VLSFO_LHV_MJ_KG = 40.5  # editable engineering reference assumption
MARINE_PRICE_YEAR_FACTOR = 1.0 # explicitly historical nominal 2022;
MARINE_QUALITY_DISCOUNT = 0.0 # scenario
ESTERIFICATION_USD_KG = 0.08 # engineering screening netback;

def marine_price(stream):
    if stream.F_mass <= 0: return 0.0
    lhv = float(stream.LHV)/1000/float(stream.F_mass)
    if not math.isfinite(lhv) or lhv <= 0:
        raise ValueError('Marine stream requires finite positive LHV')
    return VLSFO_USD_TONNE/1000*MARINE_PRICE_YEAR_FACTOR*lhv/VLSFO_LHV_MJ_KG*(1-MARINE_QUALITY_DISCOUNT)

def refinery_price(stream, hhv_per_gge):
    if stream.F_mass <= 0: return 0.0
    # Same HHV convention as the inherited Binder price calculation.
    return PNNL_USD_GGE*float(stream.HHV)/1000/hhv_per_gge/stream.F_mass

def apply_price(system, case, hhv_per_gge):
    fs=system.flowsheet.stream
    details={}
    if case.startswith('Refinery_'):
        s=fs.biofuel
        s.price=refinery_price(s,hhv_per_gge)
        details.update(Price_basis='PNNL-37046 refinery BEV transferred to lighter cut; HHV GGE convention',PNNL_USD_GGE=PNNL_USD_GGE,HHV_MJ_GGE=hhv_per_gge)
    elif case=='Binder_esterification':
        s=fs.biofuel
        gross=float(s.price)
        if gross < ESTERIFICATION_USD_KG: raise ValueError('Netback exceeds gross selling price')
        s.price=gross-ESTERIFICATION_USD_KG
        details.update(Gross_biofuel_USD_kg=gross,Esterification_USD_kg=ESTERIFICATION_USD_KG,External_processing_deduction_USD_year=s.F_mass*system.operating_hours*ESTERIFICATION_USD_KG)
    elif case=='Joint_marine':
        s=fs.intermediate;s.price=marine_price(s)
        details.update(VLSFO_USD_tonne=VLSFO_USD_TONNE,VLSFO_reference_date='2022-05-16',VLSFO_LHV_MJ_kg=VLSFO_LHV_MJ_KG,Marine_price_year_factor=MARINE_PRICE_YEAR_FACTOR,Marine_quality_discount=MARINE_QUALITY_DISCOUNT)
    else: return details
    details.update(Priced_stream=s.ID,Sale_price_USD_kg=s.price,Sale_mass_kg_year=s.F_mass*system.operating_hours)
    return details
