#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Sat Dec 30 08:01:06 2023

@author: jiananfeng

Note the word 'sludge' in this file refers to either sludge or biosolids.
'''

import qsdsan as qs
from chaospy import distributions as shape
from qsdsan.utils import auom, DictAttrSetter
from exposan.htl import create_geospatial_system

# for parameters/numbers not explained in the model, see notes in the system for references

__all__ = ('create_geospatial_model',)

# kg/m3
biocrude_density = 980
# kg/m3, this is for sludge with a moisture content higher than 80%,
# google 'Design of wastewater treatment sludge thickeners Iowa State University'
sludge_density = 1000
# kg/L
water_density = 1
# 42-47 MJ/kg, https://world-nuclear.org/information-library/facts-and-figures/heat-values-of-various-fuels
crude_oil_HHV = 44.5

_m3perh_to_MGD = auom('m3/h').conversion_factor('MGD')
_MMgal_to_L = auom('gal').conversion_factor('L')*1000000
_mile_to_km = auom('mile').conversion_factor('km')
_lb_to_kg = auom('lb').conversion_factor('kg')
_m3_to_ft3 = auom('m3').conversion_factor('ft3')
_oil_barrel_to_m3 = auom('oil_barrel').conversion_factor('m3')
_oil_barrel_to_L = auom('oil_barrel').conversion_factor('L')

# GDPCTPI (Gross Domestic Product: Chain-type Price Index), [2]
GDPCTPI = {2007: 86.352,
           2008: 87.977,
           2009: 88.557,
           2010: 89.619,
           2011: 91.466,
           2012: 93.176,
           2013: 94.786,
           2014: 96.436,
           2015: 97.277,
           2016: 98.208,
           2017: 100.000,
           2018: 102.290,
           2019: 104.008,
           2020: 105.407,
           2021: 110.220,
           2022: 117.995,
           2023: 122.284}

# TODO: why cannot run more than one time of model = create_geospatial_model(system=sys, test_run=True)

# TODO: add uncertainty for new units (parameters, cost, CI, etc.)
def create_geospatial_model(system=None,
                            test_run=False,
                            sludge_ash=[],
                            sludge_lipid=[],
                            sludge_protein=[],
                            include_check=True):
    '''
    Create a model based on the given system
    (or create the system based on the given configuration).
    
    Parameters
    ----------
    system : obj or str
        Can be either a :class:`System` objective for which the model will be created,
        or one of the allowed configurations ("baseline", "no_P", "PSA").
    '''
    
    # TODO: remove [0]
    sys = create_geospatial_system() if not system else system
    flowsheet = sys.flowsheet
    unit = flowsheet.unit
    stream = flowsheet.stream
    model = qs.Model(sys)
    param = model.parameter
    metric = model.metric
    
    WWTP = unit.WWTP
    H1 = unit.H1
    HTL = unit.HTL
    CHG = unit.CHG
    CHP = unit.CHP
    
    raw_wastewater = stream.raw_wastewater
    H2SO4 = stream.H2SO4
    # TODO: update this one in other HTL models, which might have CHG_catalyst_in = CHG.ins[1]
    virgin_CHG_catalyst = stream.virgin_CHG_catalyst
    NaOH = stream.NaOH
    biocrude = stream.biocrude
    natural_gas = stream.natural_gas
    solid_ash = stream.solid_ash
    cooling_tower_chemicals = stream.cooling_tower_chemicals
    
    tea = sys.TEA
    lca = sys.LCA
    
    # =========================================================================
    # WWTP
    # =========================================================================
    if test_run:        
        dist = shape.Triangle(0.174,0.257,0.414)
        @param(name='sludge_dw_ash',
               element=WWTP,
               kind='coupled',
               units='-',
               baseline=0.257,
               distribution=dist)
        def set_sludge_dw_ash(i):
            WWTP.sludge_dw_ash=i
        
        dist = shape.Triangle(0.08,0.204,0.308)
        @param(name='sludge_afdw_lipid',
               element=WWTP,
               kind='coupled',
               units='-',
               baseline=0.204,
               distribution=dist)
        def set_sludge_afdw_lipid(i):
            WWTP.sludge_afdw_lipid=i
        
        dist = shape.Triangle(0.38,0.463,0.51)
        @param(name='sludge_afdw_protein',
               element=WWTP,
               kind='coupled',
               units='-',
               baseline=0.463,
               distribution=dist)
        def set_sludge_afdw_protein(i):
            WWTP.sludge_afdw_protein=i
    else:
        # triangle for no_digestion (since we have enough datapoints to generate triangular distributions)
        if sludge_ash[-1] == 'no_digestion':
            dist = shape.Triangle(sludge_ash[0],sludge_ash[1],sludge_ash[2])
            @param(name='sludge_dw_ash',
                   element=WWTP,
                   kind='coupled',
                   units='-',
                   baseline=sludge_ash[1],
                   distribution=dist)
            def set_sludge_dw_ash(i):
                WWTP.sludge_dw_ash=i
            
            dist = shape.Triangle(sludge_lipid[0],sludge_lipid[1],sludge_lipid[2])
            @param(name='sludge_afdw_lipid',
                   element=WWTP,
                   kind='coupled',
                   units='-',
                   baseline=sludge_lipid[1],
                   distribution=dist)
            def set_sludge_afdw_lipid(i):
                WWTP.sludge_afdw_lipid=i
            
            dist = shape.Triangle(sludge_protein[0],sludge_protein[1],sludge_protein[2])
            @param(name='sludge_afdw_protein',
                   element=WWTP,
                   kind='coupled',
                   units='-',
                   baseline=sludge_protein[1],
                   distribution=dist)
            def set_sludge_afdw_protein(i):
                WWTP.sludge_afdw_protein=i
        # uniform for digestion (since we just have one datapoint for the cases of digestion)
        else:
            dist = shape.Uniform(sludge_ash[0],sludge_ash[2])
            @param(name='sludge_dw_ash',
                   element=WWTP,
                   kind='coupled',
                   units='-',
                   baseline=sludge_ash[1],
                   distribution=dist)
            def set_sludge_dw_ash(i):
                WWTP.sludge_dw_ash=i
            
            dist = shape.Uniform(sludge_lipid[0],sludge_lipid[2])
            @param(name='sludge_afdw_lipid',
                   element=WWTP,
                   kind='coupled',
                   units='-',
                   baseline=sludge_lipid[1],
                   distribution=dist)
            def set_sludge_afdw_lipid(i):
                WWTP.sludge_afdw_lipid=i
            
            dist = shape.Uniform(sludge_protein[0],sludge_protein[2])
            @param(name='sludge_afdw_protein',
                   element=WWTP,
                   kind='coupled',
                   units='-',
                   baseline=sludge_protein[1],
                   distribution=dist)
            def set_sludge_afdw_protein(i):
                WWTP.sludge_afdw_protein=i
    
    # the uncertainty range here covers both sludge and biosolids
    dist = shape.Triangle(0.1944,0.3927,0.5556)
    @param(name='N_2_P',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.3927,
           distribution=dist)
    def set_N_2_P(i):
        WWTP.N_2_P=i

    dist = shape.Uniform(0.675,0.825)
    @param(name='lipid_2_C',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.75,
           distribution=dist)
    def set_lipid_2_C(i):
        WWTP.lipid_2_C=i
    
    dist = shape.Uniform(0.4905,0.5995)
    @param(name='protein_2_C',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.545,
           distribution=dist)
    def set_protein_2_C(i):
        WWTP.protein_2_C=i
    
    dist = shape.Uniform(0.36,0.44)
    @param(name='carbo_2_C',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.4,
           distribution=dist)
    def set_carbo_2_C(i):
        WWTP.carbo_2_C=i
    
    dist = shape.Uniform(0.1125,0.1375)
    @param(name='lipid_2_H',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.125,
           distribution=dist)
    def set_lipid_2_H(i):
        WWTP.lipid_2_H=i
    
    dist = shape.Uniform(0.0614,0.075)
    @param(name='protein_2_H',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.0682,
           distribution=dist)
    def set_protein_2_H(i):
        WWTP.protein_2_H=i
    
    dist = shape.Uniform(0.06,0.0733)
    @param(name='carbo_2_H',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.0667,
           distribution=dist)
    def set_carbo_2_H(i):
        WWTP.carbo_2_H=i

    dist = shape.Uniform(0.1432,0.1750)
    @param(name='protein_2_N',
           element=WWTP,
           kind='coupled',
           units='-',
           baseline=0.159,
           distribution=dist)
    def set_protein_2_N(i):
        WWTP.protein_2_N=i
    
    # =========================================================================
    # HTL
    # =========================================================================
    dist = shape.Uniform(0.0170348,0.0227131)
    @param(name='enforced heating transfer coefficient',
           element=H1,
           kind='coupled',
           units='kW/m2/K',
           baseline=0.0198739,
           distribution=dist)
    def set_U(i):
        H1.U=i
    
    dist = shape.Triangle(0.692,0.846,1)
    @param(name='lipid_2_biocrude',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.846,
           distribution=dist)
    def set_lipid_2_biocrude(i):
        HTL.lipid_2_biocrude=i
    
    dist = shape.Normal(0.445,0.030)
    @param(name='protein_2_biocrude',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.445,
           distribution=dist)
    def set_protein_2_biocrude(i):
        HTL.protein_2_biocrude=i
    
    dist = shape.Normal(0.205,0.050)
    @param(name='carbo_2_biocrude',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.205,
           distribution=dist)
    def set_carbo_2_biocrude(i):
        HTL.carbo_2_biocrude=i
    
    dist = shape.Normal(0.074,0.020)
    @param(name='protein_2_gas',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.074,
           distribution=dist)
    def set_protein_2_gas(i):
        HTL.protein_2_gas=i
    
    dist = shape.Normal(0.418,0.030)
    @param(name='carbo_2_gas',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.418,
           distribution=dist)
    def set_carbo_2_gas(i):
        HTL.carbo_2_gas=i
        
    dist = shape.Normal(-8.370,0.939)
    @param(name='biocrude_C_slope',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=-8.370,
           distribution=dist)
    def set_biocrude_C_slope(i):
        HTL.biocrude_C_slope=i
        
    dist = shape.Normal(68.55,0.367)
    @param(name='biocrude_C_intercept',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=68.55,
           distribution=dist)
    def set_biocrude_C_intercept(i):
        HTL.biocrude_C_intercept=i
        
    dist = shape.Normal(0.133,0.005)
    @param(name='biocrude_N_slope',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.133,
           distribution=dist)
    def set_biocrude_N_slope(i):
        HTL.biocrude_N_slope=i
        
    dist = shape.Normal(-2.610,0.352)
    @param(name='biocrude_H_slope',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=-2.610,
           distribution=dist)
    def set_biocrude_H_slope(i):
        HTL.biocrude_H_slope=i
    
    dist = shape.Normal(8.200,0.138)
    @param(name='biocrude_H_intercept',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=8.200,
           distribution=dist)
    def set_biocrude_H_intercept(i):
        HTL.biocrude_H_intercept=i
        
    dist = shape.Normal(478,18.878)
    @param(name='HTLaqueous_C_slope',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=478,
           distribution=dist)
    def set_HTLaqueous_C_slope(i):
        HTL.HTLaqueous_C_slope=i
        
    dist = shape.Triangle(0.715,0.764,0.813)
    @param(name='TOC_TC',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.764,
           distribution=dist)
    def set_TOC_TC(i):
        HTL.TOC_TC=i
    
    dist = shape.Normal(1.750,0.122)
    @param(name='hydrochar_C_slope',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=1.750,
           distribution=dist)
    def set_hydrochar_C_slope(i):
        HTL.hydrochar_C_slope=i
    
    dist = shape.Triangle(0.035,0.063,0.102)
    @param(name='biocrude_moisture_content',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.063,
           distribution=dist)
    def set_biocrude_moisture_content(i):
        HTL.biocrude_moisture_content=i
    
    dist = shape.Uniform(0.84,0.88)
    @param(name='hydrochar_P_recovery_ratio',
           element=HTL,
           kind='coupled',
           units='-',
           baseline=0.86,
           distribution=dist)
    def set_hydrochar_P_recovery_ratio(i):
        HTL.hydrochar_P_recovery_ratio=i
    
    # =========================================================================
    # CHG
    # =========================================================================
    dist = shape.Triangle(2.86,3.562,3.99)
    @param(name='WHSV',
           element=CHG,
           kind='coupled',
           units='kg/hr/kg',
           baseline=3.562,
           distribution=dist)
    def set_CHG_WHSV(i):
        CHG.WHSV=i
    
    dist = shape.Triangle(3960,7920,15840)
    @param(name='catalyst_lifetime',
           element=CHG,
           kind='coupled',
           units='hr',
           baseline=7920,
           distribution=dist)
    def set_CHG_catalyst_lifetime(i):
        CHG.catalyst_lifetime=i
    
    dist = shape.Triangle(0.1893,0.5981,0.7798)
    @param(name='gas_C_2_total_C',
           element=CHG,
           kind='coupled',
           units='-',
           baseline=0.5981,
           distribution=dist)
    def set_gas_C_2_total_C(i):
        CHG.gas_C_2_total_C=i
    
    # =========================================================================
    # TEA
    # =========================================================================
    dist = shape.Triangle(0.6,1,1.4)
    @param(name='HTL_TIC_factor',
           element='TEA',
           kind='isolated',
           units='-',
           baseline=1,
           distribution=dist)
    def set_HTL_TIC_factor(i):
        HTL.TIC_factor=i
        
    dist = shape.Triangle(0.6,1,1.4)
    @param(name='CHG_TIC_factor',
           element='TEA',
           kind='isolated',
           units='-',
           baseline=1,
           distribution=dist)
    def set_CHG_TIC_factor(i):
        CHG.TIC_factor=i
    
    dist = shape.Uniform(980,1470)
    @param(name='CHP_unit_TIC',
           element='TEA',
           kind='isolated',
           units='-',
           baseline=1225,
           distribution=dist)
    def set_CHP_unit_TIC(i):
        CHP.unit_TIC=i
    
    dist = shape.Triangle(0,0.03,0.05)
    @param(name='IRR',
           element='TEA',
           kind='isolated',
           units='-',
           baseline=0.03,
           distribution=dist)
    def set_IRR(i):
        tea.IRR=i
    
    # the sources of chemical prices can be found in geospatial_HTL_systems.py
    H2SO4_price = (0.043*1+0.0002*(93/5-1))/(93/5)/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    dist = shape.Uniform(H2SO4_price*0.9,H2SO4_price*1.1)
    @param(name='5% H2SO4 price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=H2SO4_price,
           distribution=dist)
    def set_H2SO4_price(i):
        H2SO4.price=i
    
    NaOH_price = 0.2384/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    dist = shape.Uniform(NaOH_price*0.9,NaOH_price*1.1)
    @param(name='NaOH price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=NaOH_price,
           distribution=dist)
    def set_NaOH_price(i):
        NaOH.price=i
    
    CHG_catalyst_price = 60/_lb_to_kg/GDPCTPI[2011]*GDPCTPI[2022]
    dist = shape.Triangle(CHG_catalyst_price*0.5,CHG_catalyst_price,CHG_catalyst_price*2)
    @param(name='CHG catalyst price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=CHG_catalyst_price,
           distribution=dist)
    def set_catalyst_price(i):
        virgin_CHG_catalyst.price=i
    
    biocrude_price_min = 71.59/_oil_barrel_to_m3/biocrude_density/crude_oil_HHV*HTL.biocrude_HHV
    biocrude_price_max = 123.70/_oil_barrel_to_m3/biocrude_density/crude_oil_HHV*HTL.biocrude_HHV
    biocrude_price_ave = 94.53/_oil_barrel_to_m3/biocrude_density/crude_oil_HHV*HTL.biocrude_HHV
    dist = shape.Triangle(biocrude_price_min,biocrude_price_ave,biocrude_price_max)
    @param(name='biocrude price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=biocrude_price_ave,
           distribution=dist)
    def set_biocrude_price(i):
        biocrude.price=i
    
    # TODO: update in writing
    # from _heat_utility.py (biosteam): 3.49672 $/kmol
    natural_gas_price_min = 0.218545*0.9
    natural_gas_price_max = 0.218545*1.1
    natural_gas_price_ave = 0.218545
    dist = shape.Triangle(natural_gas_price_min,natural_gas_price_ave,natural_gas_price_max)
    @param(name='natural gas price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=natural_gas_price_ave,
           distribution=dist)
    def set_CH4_price(i):
        natural_gas.price=i
    
    ash_disposal_price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2022]
    dist = shape.Uniform(ash_disposal_price*1.1,ash_disposal_price*0.9)
    @param(name='ash disposal price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=ash_disposal_price,
           distribution=dist)
    def set_ash_disposal_price(i):
        solid_ash.price=i
    
    cooling_chemicals_price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    dist = shape.Uniform(cooling_chemicals_price*0.9,cooling_chemicals_price*1.1)
    @param(name='cooling tower chemicals price',
           element='TEA',
           kind='isolated',
           units='$/kg',
           baseline=cooling_chemicals_price,
           distribution=dist)
    def set_cooling_tower_chemicals_price(i):
        cooling_tower_chemicals.price=i
    
    sludge_transportation_price = WWTP.ww_2_dry_sludge*\
                                    (4.56/sludge_density*1000/0.2+0.072/_mile_to_km/sludge_density*1000/0.2*WWTP.sludge_distance)/\
                                        GDPCTPI[2015]*GDPCTPI[2022]/3.79/(10**6)/WWTP.sludge_distance
    dist = shape.Uniform(sludge_transportation_price*0.9,sludge_transportation_price*1.1)
    @param(name='sludge transportation price',
           element='TEA',
           kind='isolated',
           units='$/kg/km',
           baseline=sludge_transportation_price,
           distribution=dist)
    def set_sludge_transportation_price(i):
        qs.ImpactItem.get_all_items()['Sludge_trucking'].price=i
    
    
    biocrude_transportation_price = (5.67/biocrude_density+0.07/biocrude_density*WWTP.biocrude_distance)/GDPCTPI[2008]*GDPCTPI[2022]/WWTP.biocrude_distance
    dist = shape.Uniform(biocrude_transportation_price*0.9,biocrude_transportation_price*1.1)
    @param(name='biocrude transportation price',
           element='TEA',
           kind='isolated',
           units='$/kg/km',
           baseline=biocrude_transportation_price,
           distribution=dist)
    def set_biocrude_transportation_price(i):
        qs.ImpactItem.get_all_items()['Biocrude_trucking'].price=i
    
    # =========================================================================
    # anhydrous NH3
    # =========================================================================
    # TODO: add parameters
        
    # =========================================================================
    # urea
    # =========================================================================
    # TODO: add parameters  
    try:
        makeup_MEA = stream.makeup_MEA
        MEA_price_min = 1.93
        MEA_price_max = 2.31
        MEA_price_ave = 2.13
        dist = shape.Triangle(MEA_price_min,MEA_price_ave,MEA_price_max)
        @param(name='MEA price',
               element='TEA',
               kind='isolated',
               units='$/kg',
               baseline=MEA_price_ave,
               distribution=dist)
        def set_MEA_price(i):
            makeup_MEA.price=i
    except AttributeError:
        pass
        
    # =========================================================================
    # UAN
    # =========================================================================
    # TODO: add parameters
    try:
        HNO3 = stream.HNO3
        HNO3_price_min = 0.43*0.7 + 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]*0.3
        HNO3_price_max = 0.53*0.7 + 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]*0.3
        HNO3_price_ave = 0.497*0.7 + 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]*0.3
        dist = shape.Triangle(HNO3_price_min,HNO3_price_ave,HNO3_price_max)
        @param(name='HNO3 price',
               element='TEA',
               kind='isolated',
               units='$/kg',
               baseline=HNO3_price_ave,
               distribution=dist)
        def set_HNO3_price(i):
            HNO3.price=i
    except AttributeError:
        pass
    
    # =========================================================================
    # LCA (unifrom ± 10%)
    # =========================================================================
    CF = 'GlobalWarming'
    # do not get joint distribution for multiple times, since the baselines for LCA will change
    for item in qs.ImpactItem.get_all_items().keys():
        if qs.ImpactItem.get_item(item).CFs and qs.ImpactItem.get_item(item).CFs[CF] != 0:
            abs_small = 0.9*qs.ImpactItem.get_item(item).CFs[CF]
            abs_large = 1.1*qs.ImpactItem.get_item(item).CFs[CF]
            dist = shape.Uniform(min(abs_small,abs_large),max(abs_small,abs_large))
            @param(name=f'{item}_{CF}',
                   setter=DictAttrSetter(qs.ImpactItem.get_item(item),'CFs',CF),
                   element='LCA',
                   kind='isolated',
                   units=qs.ImpactIndicator.get_indicator(CF).unit,
                   baseline=qs.ImpactItem.get_item(item).CFs[CF],
                   distribution=dist)
            def set_LCA(i):
                qs.ImpactItem.get_item(item).CFs[CF]=i
    
    # =========================================================================
    # metrics
    # =========================================================================
    @metric(name='sludge_management_cost', units='$/tonne dry sludge', element='geospatial')
    def get_sludge_management_cost():
        return -tea.solve_price(raw_wastewater)*water_density*_MMgal_to_L/WWTP.ww_2_dry_sludge
    
    @metric(name='sludge_CI', units='kg CO2/tonne dry sludge', element='geospatial')
    def get_sludge_CI():
        return lca.get_total_impacts(operation_only=True,
                                     exclude=(raw_wastewater,),
                                     annual=True)['GlobalWarming']/raw_wastewater.F_vol/_m3perh_to_MGD/WWTP.ww_2_dry_sludge/(sys.operating_hours/24)
    
    @metric(name='biocrude_production', units='BPD', element='geospatial')
    def get_biocrude_production():
        return biocrude.F_mass/biocrude_density*1000/_oil_barrel_to_L*24
    
    # TODO: add N production/offset and P production/offset
    
    if include_check:
        @metric(name='sludge_afdw_carbohydrate', units='-', element='test')
        def get_sludge_afdw_carbohydrate():
            return WWTP.sludge_afdw_carbo
    
    return model