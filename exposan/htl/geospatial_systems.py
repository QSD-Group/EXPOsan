#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Mon Jun 5 08:46:28 2023

@author: jiananfeng

Note the word 'sludge' in this file refers to either sludge or biosolids.

References:
[1] Snowden-Swan, L. J.; Li, S.; Thorson, M. R.; Schmidt, A. J.; Cronin, D. J.;
    Zhu, Y.; Hart, T. R.; Santosa, D. M.; Fox, S. P.; Lemmon, T. L.; Swita, M. S.
    Wet Waste Hydrothermal Liquefaction and Biocrude Upgrading to Hydrocarbon Fuels:
    2022 State of Technology; PNNL-33622; Pacific Northwest National Lab. (PNNL),
    Richland, WA (United States), 2022. https://doi.org/10.2172/1897670.
[2] https://hpoilgas.in/Page/Detail/naturalgasoverview (accessed 2025-02-03).
[3] https://world-nuclear.org/information-library/facts-and-figures/heat-values-of-various-fuels
    (accessed 2025-02-03).
[4] https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20).
[5] https://data.bls.gov/cgi-bin/srgate (accessed 2024-08-06).
[6] https://mymarketnews.ams.usda.gov/report/list?market_types%5B0%5D=214
    (accessed 2025-02-01).
[7] https://www.sludgeprocessing.com/sludge-dewatering/sludge-drying-beds-lagoons/
    (accessed 2024-08-03).
[8] Metcalf & Eddy (2003) Wastewater Engineering: Treatment and Reuse. 4th Edition,
    McGraw-Hill, New York.
[9] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.;
    Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.;
    Snowden-Swan, L. J.; Davis, R.; Kinchin, C. Process Design and Economics for
    the Conversion of Algal Biomass to Hydrocarbons: Whole Algae Hydrothermal
    Liquefaction and Upgrading; PNNL--23227, 1126336; 2014; p PNNL--23227, 1126336.
    https://doi.org/10.2172/1126336.
[10] https://www.macrotrends.net/1369/crude-oil-price-history-chart
    (accessed 2024-08-03).
[11] Davis, R. E.; Grundl, N. J.; Tao, L.; Biddy, M. J.; Tan, E. C.; Beckham, G. T.;
    Humbird, D.; Thompson, D. N.; Roni, M. S. Process Design and Economics for the
    Conversion of Lignocellulosic Biomass to Hydrocarbon Fuels and Coproducts:
    2018 Biochemical Design Case Update; Biochemical Deconstruction and Conversion
    of Biomass to Fuels and Products via Integrated Biorefinery Pathways;
    NREL/TP--5100-71949, 1483234; 2018; p NREL/TP--5100-71949, 1483234.
    https://doi.org/10.2172/1483234.
[12] https://businessanalytiq.com/procurementanalytics/index/monoethanolamine-price-index/
    (accessed 2025-02-03)
[13] https://businessanalytiq.com/procurementanalytics/index/nitric-acid-price-index/
    (accessed 2025-02-03).
[14] Marufuzzaman, M.; Ekşioğlu, S. D.; Hernandez, R. Truck versus Pipeline
     Transportation Cost Analysis of Wastewater Sludge. Transportation Research Part A:
     Policy and Practice 2015, 74, 14–30. https://doi.org/10.1016/j.tra.2015.02.001.
[15] Pootakham, T.; Kumar, A. Bio-Oil Transport by Pipeline: A Techno-Economic
     Assessment. Bioresource Technology 2010, 101 (18), 7137–7143.
     https://doi.org/10.1016/j.biortech.2010.03.136.
[16] Snowden-Swan, L. J.; Zhu, Y.; Bearden, M. D.; Seiple, T. E.; Jones, S. B.;
     Schmidt, A. J.; Billing, J. M.; Hallen, R. T.; Hart, T. R.; Liu, J.;
     Albrecht, K. O.; Fox, S. P.; Maupin, G. D.; Elliott, D. C.
     Conceptual Biorefinery Design and Research Targeted for 2022:
     Hydrothermal Liquefacation Processing of Wet Waste to Fuels; PNNL-27186;
     Pacific Northwest National Lab. (PNNL), Richland, WA (United States), 2017.
     https://doi.org/10.2172/1415710.
[17] Stewart, D. W.; Cortés-Peña, Y. R.; Li, Y.; Stillwell, A. S.; Khanna, M.;
     Guest, J. S. Implications of Biorefinery Policy Incentives and
     Location-SpecificEconomic Parameters for the Financial Viability of Biofuels.
     Environ. Sci. Technol. 2023. https://doi.org/10.1021/acs.est.2c07936.
'''

# TODO: changes to implement:
# TODO 1: add new units and contextualized fertilizers prices (https://mymarketnews.ams.usda.gov/report/list?market_types%5B0%5D=214) in writing

import os, qsdsan as qs, biosteam as bst, pandas as pd
from qsdsan import sanunits as qsu
from qsdsan.utils import auom, clear_lca_registries
from exposan.htl import _load_components, create_tea, state_income_tax_rate_2022, _sanunits as su
from biosteam.units import IsenthalpicValve, Stripper, MolecularSieve, IsothermalCompressor, Splitter
from biosteam import settings

__all__ = ('create_geospatial_system','biocrude_density')

# kg/m3, [1]
biocrude_density = 980
# kg/m3, this is for sludge with a moisture content higher than 80%,
# google 'Design of wastewater treatment sludge thickeners Iowa State University'
sludge_density = 1000
# kg/m3, [2]
natural_gas_density = 0.65
# 42-47 MJ/kg, [3]
crude_oil_HHV = 44.5

_mile_to_km = auom('mile').conversion_factor('km')
_lb_to_kg = auom('lb').conversion_factor('kg')
_m3_to_ft3 = auom('m3').conversion_factor('ft3')
_oil_barrel_to_m3 = auom('oil_barrel').conversion_factor('m3')

# GDPCTPI (Gross Domestic Product: Chain-type Price Index), [4]
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

# the labor index can be found in [5] with the series id CEU3232500008,
# remember to select 'include annual average'
labor_index = {2014: 21.49,
               2015: 21.76,
               2016: 22.71,
               2017: 24.29,
               2018: 25.46,
               2019: 25.46,
               2020: 26.03,
               2021: 26.69,
               2022: 27.36,
               2023: 29.77}

# for parameters, unless otherwise stated, refer to the original HTL system model
def create_geospatial_system(# MGD
                             size=10,
                             # 0: no; 1: yes
                             sludge_transportation=0,
                             # in km, this is the slduge transportation total
                             # distance (normalized to total sludge amount)
                             sludge_distance=100,
                             # km
                             biocrude_distance=100,
                             # average values below are for sludge aggregation analyses
                             average_sludge_dw_ash=None,
                             average_sludge_afdw_lipid=None,
                             average_sludge_afdw_protein=None,
                             # 0: no; 1: yes
                             anaerobic_digestion=0,
                             # 0: no; 1: yes
                             aerobic_digestion=0,
                             # dry tonne sludge/day/MGD raw wastewater
                             ww_2_dry_sludge_ratio=1,
                             state='IL',
                             nitrogen_fertilizer='NH3',
                             # TODO: update the balancing-area-level grid CI in geospatial_analysis.py
                             # use balancing-area-level in the analysis
                             # kg CO2 eq/kWh
                             elec_GHG=0.40,
                             # use county-level in the analysis
                             wage_adjustment=1
                             ):
    
    if nitrogen_fertilizer not in ('NH3','urea','UAN'):
        raise ValueError("nitrogen_fertilizer can only be 'NH3', 'urea', or 'UAN'.")
    
    flowsheet_ID = 'htl_geospatial'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    folder = os.path.dirname(__file__)
    
    # electricity price in 2022$/kWh
    elec_price = pd.read_excel(folder + '/data/state_elec_price_2022.xlsx', 'elec_price_2022')
    bst.PowerUtility.price = elec_price[elec_price['state']==state]['price'].iloc[0]/100
    
    # fertilizer price in $/US-ton, [6]
    # UAN30 price includes UAN28/30/30-32
    if state == 'AL':
        DAP_price = 972.5
        NH3_price = 1407.5
        urea_price = 850
        UAN30_price = 705
    
    elif state == 'IA':
        DAP_price = 1000
        NH3_price = 1386.5
        urea_price = 867
        UAN30_price = 658.5
        
    elif state == 'IL':
        DAP_price = 962.5
        NH3_price = 1417.5
        urea_price = 875
        UAN30_price = 702.5
        
    elif state == 'NC':
        DAP_price = 995
        NH3_price = 1407.5
        urea_price = 855
        UAN30_price = 605
        
    elif state == 'OK':
        DAP_price = 1022
        NH3_price = 1312.5
        urea_price = 882.5
        UAN30_price = 607.5
        
    elif state == 'SC':
        DAP_price = 1002.2
        NH3_price = 1407.5
        urea_price = 895
        UAN30_price = 587.5
        
    else:
        DAP_price = 984.5
        NH3_price = 1407.5
        urea_price = 855
        UAN30_price = 605
    
    # convert price to $/kg
    DAP_price = DAP_price/0.907185/1000
    NH3_price = NH3_price/0.907185/1000
    urea_price = urea_price/0.907185/1000
    UAN30_price = UAN30_price/0.907185/1000
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # =========================================================================
    # pretreatment
    # =========================================================================
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)
    
    # assume the moisture content of sludge is 80% in all cases
    # for lagoon, the sludge will dry at the base of the lagoon (to an assumed 80% moisture content, see [7])
    # from [1]:
    # ash (dw%) of undigested sludge: 0.266, 0.192, 0.237, 0.174, 0.206, 0.308 (average: 0.231)
    # protein (afdw%) of undigested sludge: 0.464, 0.454, 0.38, 0.467, 0.485, 0.484 (average: 0.456)
    # lipid (afdw%) of undigested sludge: 0.308, 0.08, 0.27, 0.176, 0.178, 0.225 (average: 0.206)
    # ash (dw%), protein (afdw%), and lipid (afdw%) for anaerobically digested sludge are 0.414, 0.510, and 0.193, respectively
    # for aerobically digested sludge, assume the same biological compositions as anaerobically digested sludge,
    # but with a higher ash content
    # mass reduction assumption (from IEDO work, originally from [8]):
    # 42.5% VSS reduction for anaerobic digestion, 47.5% VSS reduction for aerobic digestion
    # assume X in sludge-to-be-digested (not necessarily having the same biochemical compositions as sludge)
    # is ash (cannot be digested), Y is VSS (can be digested)
    # for anaerobic digestion: X/(X+Y*(1-0.425)) = 0.414
    # we can calculate that for aerobic digestion: X/(X+Y*(1-0.475)) = 0.436
    if sludge_transportation == 0:
        if aerobic_digestion == 1:
            WWTP = su.WWTP(ID='WWTP',
                           ins=raw_wastewater,
                           outs=('sludge','treated_water'),
                           ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                           sludge_moisture=0.8,
                           sludge_dw_ash=0.436,
                           sludge_afdw_lipid=0.193,
                           sludge_afdw_protein=0.510,
                           operation_hours=8760,
                           sludge_distance=sludge_distance,
                           biocrude_distance=biocrude_distance)
        elif anaerobic_digestion == 1:
            WWTP = su.WWTP(ID='WWTP',
                           ins=raw_wastewater,
                           outs=('sludge','treated_water'),
                           ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                           sludge_moisture=0.8,
                           sludge_dw_ash=0.414,
                           sludge_afdw_lipid=0.193,
                           sludge_afdw_protein=0.510,
                           operation_hours=8760,
                           sludge_distance=sludge_distance,
                           biocrude_distance=biocrude_distance)
        else:
            WWTP = su.WWTP(ID='WWTP',
                           ins=raw_wastewater, outs=('sludge','treated_water'),
                           ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                           sludge_moisture=0.8,
                           sludge_dw_ash=0.231,
                           sludge_afdw_lipid=0.206,
                           sludge_afdw_protein=0.456,
                           operation_hours=8760,
                           sludge_distance=sludge_distance,
                           biocrude_distance=biocrude_distance)
    else:
        assert average_sludge_dw_ash != None, 'set average_sludge_dw_ash manually'
        assert average_sludge_afdw_lipid != None, 'set average_sludge_afdw_lipid manually'
        assert average_sludge_afdw_protein != None, 'set average_sludge_afdw_protein manually'
        
        WWTP = su.WWTP(ID='WWTP',
                       ins=raw_wastewater,
                       outs=('sludge','treated_water'),
                       ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                       sludge_moisture=0.8,
                       sludge_dw_ash=average_sludge_dw_ash,
                       sludge_afdw_lipid=average_sludge_afdw_lipid,
                       sludge_afdw_protein=average_sludge_afdw_protein,
                       operation_hours=8760,
                       sludge_distance=sludge_distance,
                       biocrude_distance=biocrude_distance)
    
    P1 = qsu.SludgePump(ID='P1',
                        ins=WWTP-0,
                        outs='pressed_sludge',
                        P=3049.7*6894.76,
                        init_with='Stream')
    
    # =========================================================================
    # HTL
    # =========================================================================
    H1 = qsu.HXutility(ID='H1',
                       ins=P1-0,
                       outs='heated_sludge',
                       T=351+273.15,
                       U=0.0198739,
                       init_with='Stream',
                       rigorous=True)
    
    # assume no value/cost and no environmental benefit/impact associated with hydrochar
    HTL = qsu.HydrothermalLiquefaction(ID='HTL',
                                       ins=H1-0,
                                       outs=('hydrochar','HTL_aqueous_undefined',
                                             'biocrude_to_be_stored','offgas_HTL'),
                                       mositure_adjustment_exist_in_the_system=False)
    
    # =========================================================================
    # CHG
    # =========================================================================
    HTLaqueous = su.HTLaqueous(ID='HTLaqueous',
                               ins=HTL-1,
                               outs='HTL_aqueous_defined')
    
    # TODO: add the following note to other HTL systems
    # 7.8%_Ru/C and 7.8%_Ru/C_out are not valid names; they are not in sys.flowsheet.stream but the price and LCA are included
    # although this have no impact on the results, still change the names to 'virgin_CHG_catalyst and 'used_CHG_catalyst'
    CHG = qsu.CatalyticHydrothermalGasification(ID='CHG',
                                                ins=(HTLaqueous-0, 'virgin_CHG_catalyst'),
                                                outs=('CHG_out','used_CHG_catalyst'))
    # CHG price: [9]
    CHG.ins[1].price = 60/_lb_to_kg/GDPCTPI[2011]*GDPCTPI[2022]
    
    S2WS1 = su.StreamTypeConverter(ID='S2WS1',
                                   ins=CHG-1,
                                   outs='CHG_catalyst_out',
                                   init_with='WasteStream')
    
    V1 = IsenthalpicValve(ID='V1',
                          ins=CHG-0,
                          outs='depressed_cooled_CHG',
                          P=50*6894.76,
                          vle=True)
    
    F1 = qsu.Flash(ID='F1',
                   ins=V1-0,
                   outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15,
                   P=50*6894.76,
                   thermo=settings.thermo.ideal())
    
    # =========================================================================
    # Storage and disposal
    # =========================================================================
    # store for 3 days based on [9]
    BiocrudeTank = qsu.StorageTank(ID='BiocrudeTank',
                                   ins=HTL-2,
                                   outs=('biocrude'),
                                   tau=3*24,
                                   init_with='WasteStream',
                                   vessel_material='Carbon steel')
    # TODO: update biocrude cost and CI calculation in writing
    # assume biocrude replace crude oil of the same amount of energy
    # 2022 average closing price for crude oil: 94.53 $/oil barrel, [10]
    # assume biocrude has an HHV of 35 MJ/kg
    # in the model, biocrude HHV will be calculated as HTL.biocrude_HHV
    BiocrudeTank.outs[0].price = 94.53/_oil_barrel_to_m3/biocrude_density/crude_oil_HHV*35
    
    GasMixer = qsu.Mixer(ID='GasMixer',
                         ins=(HTL-3, F1-0),
                         outs=('fuel_gas'),
                         init_with='Stream')
    
    # =========================================================================
    # nutrient recovery - part 1
    # =========================================================================
    
    H2SO4_Tank = qsu.StorageTank(ID='H2SO4_Tank',
                                 ins='H2SO4',
                                 outs='H2SO4_out',
                                 init_with='WasteStream',
                                 tau=24,
                                 vessel_material='Stainless steel')
    # 0.5 M H2SO4: ~5%
    # based on 93% H2SO4 and fresh water (dilute onsite to 5%) prices in [11]
    H2SO4_Tank.ins[0].price = (0.043*1+0.0002*(93/5-1))/(93/5)/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    
    AcidEx = su.AcidExtraction(ID='AcidEx',
                               ins=(HTL-0, H2SO4_Tank-0),
                               outs=('residual','extracted'))
    
    PreStripper = su.PreStripper(ID='PreStripper',
                                 ins=(F1-1,'NaOH'),
                                 outs='NH3_solution')
    # 0.2384 2016$/lb, [11]
    PreStripper.ins[1].price = 0.2384/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    
    # use steam could reduce the required air:liquid volume ratio from 3000:1 to ~300:1
    # set the flow here as size*2 results in a volume ratio of ~350:1
    water_steam = qs.Stream(ID='water_steam', H2O=size*2, phase='l', T=25+273.15)
    water_steam.price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    
    boiler = qsu.HXutility(ID='boiler',
                           ins=water_steam,
                           outs='steam',
                           T=390,
                           init_with='Stream',
                           rigorous=True)
    
    NH3Stripper = Stripper(ID='NH3Stripper',
                           N_stages=2,
                           ins=(PreStripper-0, boiler-0),
                           outs=('vapor','liquid'),
                           solute='NH3')
    
    # use 3 molecular sieves to ensure the removal of water
    # assume 2% (1%-5%) of NH3 loss for each pass
    NH3MS1 = MolecularSieve(ID='NH3MS1',
                            ins=NH3Stripper-0,
                            outs=('NH3_rich_1','water_rich_1'),
                            split=dict(Water=0.16,
                                       NH3=0.98))
    
    NH3MS2 = MolecularSieve(ID='NH3MS2',
                            ins=NH3MS1-0,
                            outs=('NH3_rich_2','water_rich_2'),
                            split=dict(Water=0.08,
                                       NH3=0.98))
    
    NH3MS3 = MolecularSieve(ID='NH3MS3',
                            ins=NH3MS2-0,
                            outs=('NH3_rich_3','water_rich_3'),
                            split=dict(Water=0,
                                       NH3=0.98))
    
    DAPSyn = su.DAPSynthesis(ID='DAPSyn',
                             ins=(AcidEx-1, NH3MS3-0),
                             outs=('DAP','excess_NH3','DAPSyn_effluent'))
    DAPSyn.outs[0].price = DAP_price
    
    if nitrogen_fertilizer == 'NH3':
        NH3Compressor = IsothermalCompressor(ID='NH3Compressor',
                                             ins=DAPSyn-1,
                                             outs='anhydrous_ammonia_gas',
                                             P=2e6,
                                             eta=1,
                                             vle=True)
        
        NH3Cooler = qsu.HXutility(ID='NH3Cooler',
                                  ins=NH3Compressor-0,
                                  outs='anhydrous_ammonia',
                                  T=25+273.15,
                                  init_with='Stream',
                                  rigorous=True)
        
        S2WS2 = su.StreamTypeConverter(ID='S2WS2',
                                       ins=NH3Cooler-0,
                                       outs='anhydrous_ammonia_LCA',
                                       init_with='WasteStream')
        S2WS2.outs[0].price = NH3_price
    
    # =========================================================================
    # facilities
    # =========================================================================
    # TODO: update in writing
    # previously used 86 C with self-defined heat utility
    # now use natural gas for heating, and use default T_min_app
    qsu.HeatExchangerNetwork(ID='HXN',
                             force_ideal_thermo=True)
    
    # assume no value/cost and no environmental benefit/impact associated with emission, since they are all captured, utilitzed, or included in the natural_gas item
    # the CHP here can ususally meet the heat requirement but not the electricity
    # buying additional natural gas to produce electricity does not provide benefit; therefore, set supplement_power_utility=False
    CHP = qsu.CombinedHeatPower(ID='CHP',
                                ins=(GasMixer-0, 'natural_gas', 'air'),
                                outs=('emission','solid_ash'),
                                init_with='WasteStream',
                                supplement_power_utility=False)
    # from _heat_utility.py (biosteam): 3.49672 $/kmol
    CHP.ins[1].price = 0.218545
    # 1.41 MM 2016$/year for 4270/4279 kg/hr ash, 7880 annual operating hours, from [11]
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2022]
    
    # TODO: consider adding CT and its TEA (price for cooling_tower_chemicals) and LCA items (CT_chemicals in the 'Other' category) for other systems (HTL, HTL-PFAS)
    # construction cost for CT is based on the flow rate of cooling_tower_chemicals in the current version of BioSTEAM
    CT = bst.facilities.CoolingTower(ID='CT')
    # cooling_tower_chemicals: 1.7842 2016$/lb, [11]
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    
    # CWP uses electricity to generate chilled water
    # the water cost can be ignored since the water can be recirculated
    CWP = bst.ChilledWaterPackage(ID='CWP')
    
    # =========================================================================
    # nutrient recovery - part 2
    # =========================================================================
    
    if nitrogen_fertilizer != 'NH3':
        # TODO: add cost and CI for ins
        # to be conservative, capture all CO2 from CHP emission (more than enough for urea synthesis and UAN synthesis)
        # if no natural gas directly burned in the CHP (rather than as heat utilities for other units):
        # then carbon in urea and carbon in UAN are biogenic, can reduce their CI
        # if there is natural gas directly burned in the CHP:
        # include the emission from natural gas combustion in LCA and reduce the CI of urea synthesis and UAN synthesis
        # AmineAbsorption includes a stripper based on the description
        CC = su.AmineAbsorption(ID='CC',
                                ins=(CHP-0, 'makeup_MEA', 'makeup_water'),
                                outs=('vent','CO2'),
                                CO2_recovery=0.9)
        
        # TODO: use these values in geospatial_models.py
        # min: 1.93, max: 2.31, average: 2.13, [12]
        CC.ins[1].price = 2.13
        CC.ins[2].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
        
        if nitrogen_fertilizer == 'urea':
            UreaSyn = su.UreaSynthesis(ID='UreaSyn',
                                       ins=(DAPSyn-1, CC-1, 'additional_carbon_dioxide'),
                                       outs=('urea','urea_vapor','urea_waste','excess_carbon_dioxide'))
            UreaSyn.outs[0].price = urea_price
        else:
            UANSyn = su.UANSynthesis(ID='UANSyn',
                                     ins=(DAPSyn-1, CC-1, 'additional_carbon_dioxide', 'HNO3', 'UAN_water'),
                                     outs=('UAN30','UAN_vapor','UAN_waste','UAN_excess_carbon_dioxide'),
                                     UAN_concentration=30)
            # TODO: use these values in geospatial_models.py
            # min: 0.43, max: 0.53, average: 0.497, [13]
            # calculate the price of 70 wt/wt% HNO3 solution by adding water
            UANSyn.ins[3].price = 0.497*0.7 + 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]*0.3
            UANSyn.ins[4].price = 0.0002/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
            UANSyn.outs[0].price = UAN30_price
    
    sys = qs.System.from_units(ID='sys_geospatial',
                               units=list(flowsheet.unit),
                               operating_hours=WWTP.operation_hours)
    # simulate twice to simulate CC after CHP
    sys.simulate()
    sys.simulate()
    
    # add cost and CI for carbon dioxide if there is additional CO2 usage in urea and/or UAN synthesis
    try:
        assert sys.flowsheet.additional_carbon_dioxide.F_mass == 0
    except AttributeError:
        pass
    
    # =========================================================================
    # LCA
    # =========================================================================
    # LCA data from ecoinvent 3.8 (cutoff model and TRACI method) unless otherwise stated
    
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.48748859)
    
    # deionized water
    CT_chemicals = qs.ImpactItem('CT_chemicals', functional_unit='kg')
    CT_chemicals.add_indicator(GlobalWarming, 0.00042012744)
    
    Sludge_trucking = qs.ImpactItem('Sludge_trucking', functional_unit='kg*km')
    # assume the transported sludge has 80% moisture content
    # 1 gal water = 3.79 kg water
    # 'market for transport, freight, lorry, unspecified' (0.13004958 kg CO2 eq/metric ton/km, including empty return trips)
    # assume the sludge/biosolids decomposition is minimal during transportation since our transportation distance is not long
    Sludge_trucking.add_indicator(GlobalWarming, WWTP.ww_2_dry_sludge*0.13004958/0.2/3.79/(10**6))
    # 4.56 $/m3, 0.072 $/m3/mile ([14], likely 2015$)
    Sludge_trucking.price = WWTP.ww_2_dry_sludge*\
        (4.56/sludge_density*1000/0.2+0.072/_mile_to_km/sludge_density*1000/0.2*WWTP.sludge_distance)/\
            GDPCTPI[2015]*GDPCTPI[2022]/3.79/(10**6)/WWTP.sludge_distance
    
    Sludge_transportation = qs.Transportation('Sludge_trucking',
                                              linked_unit=WWTP,
                                              item=Sludge_trucking,
                                              load_type='mass',
                                              load=stream.raw_wastewater.F_mass,
                                              load_unit='kg',
                                              distance=sludge_transportation*WWTP.sludge_distance,
                                              distance_unit='km',
                                              # set to 1 h since load = kg/h
                                              interval='1',
                                              interval_unit='h')
    WWTP.transportation = Sludge_transportation
    
    Biocrude_trucking = qs.ImpactItem('Biocrude_trucking', functional_unit='kg*km')
    # 0.13004958 kg CO2 eq/metric ton/km ('market for transport, freight, lorry, unspecified'))
    Biocrude_trucking.add_indicator(GlobalWarming, 0.13004958/1000)
    # transportation cost: 5.67 2008$/m3 (fixed cost) and 0.07 2008$/m3/km (variable cost), [15]
    Biocrude_trucking.price = (5.67/biocrude_density+0.07/biocrude_density*WWTP.biocrude_distance)/GDPCTPI[2008]*GDPCTPI[2022]/WWTP.biocrude_distance
        
    Biocrude_transportation = qs.Transportation('Biocrude_trucking',
                                                linked_unit=BiocrudeTank,
                                                item=Biocrude_trucking,
                                                load_type='mass',
                                                load=stream.biocrude.F_mass,
                                                load_unit='kg',
                                                distance=WWTP.biocrude_distance,
                                                distance_unit='km',
                                                # set to 1 h since load = kg/h
                                                interval='1',
                                                interval_unit='h')
    BiocrudeTank.transportation = Biocrude_transportation
    
    impact_items = {'CHG_catalyst': [stream.CHG_catalyst_out, 471.098936962268],
                    'H2SO4':        [stream.H2SO4, 0.005529872568],
                    'NaOH':         [stream.NaOH, 1.2497984],
                    'DAP':          [stream.DAP, -1.456692],
                    # TODO: do not use market or market group for products for other systems (i.e., CO2 sorbent, HTL-PFAS)
                    # include emission
                    'natural_gas':  [stream.natural_gas, 0.47016123/natural_gas_density+44/16],
                    # use market or market group for biocrude since we want to offset transportation and then add our own transportation part
                    # 0.22290007 kg CO2 eq/kg petroleum ('market for petroleum')
                    'biocrude':     [stream.biocrude, -0.22290007/crude_oil_HHV*HTL.biocrude_HHV],
                    'ash_disposal': [stream.solid_ash, 0.0082744841]}
    
    if nitrogen_fertilizer == 'NH3':
        impact_items['anhydrous_ammonia_LCA'] = [stream.anhydrous_ammonia_LCA, -2.4833472]
    elif nitrogen_fertilizer == 'urea':
        impact_items['makeup_MEA'] = [stream.makeup_MEA, 3.0923397]
        impact_items['makeup_water'] = [stream.makeup_water, 0.00042012744]
        # for every kg of urea produced, 44.009/60.06 kg of biogenic CO2 is used
        impact_items['urea'] = [stream.urea, -1.2510711 - 44.009/60.06]
    else:
        impact_items['makeup_MEA'] = [stream.makeup_MEA, 3.0923397]
        impact_items['makeup_water'] = [stream.makeup_water, 0.00042012744]
        impact_items['HNO3'] = [stream.HNO3, 1.624862648]
        impact_items['UAN_water'] = [stream.UAN_water, 0.00042012744]
        # the formula of UAN is CH6N4O4, for every kg of UAN30 produced, 1*0.3/14.0067/4*44.009 of biogenic CO2 is used
        impact_items['UAN30'] = [stream.UAN30, -1.6799471 - 1*0.3/14.0067/4*44.009]
    
    for item in impact_items.items():
        qs.StreamImpactItem(ID=item[0], linked_stream=item[1][0], GlobalWarming=item[1][1])
    
    qs.LCA(system=sys, lifetime=30, lifetime_unit='yr',
           # 0.48748859 is the GHG level with the Electricity item from ecoinvent,
           # we cannot list electricity GHG one state by one state,
           # but we can adjust the electricity amount to reflect different GHG of electricity at different states
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30/0.48748859*elec_GHG,
           # note LCA for CT_chemicals was included in the 'Other' category while it should be in the 'Stream' category
           # the effect is minimal since (i) this part of LCA is negligible and (ii) we do not use LCA breakdown results in the HTL geospatial analysis
           CT_chemicals=lambda:CT.ins[2].F_mass*sys.flowsheet.WWTP.operation_hours*30)
    
    # =========================================================================
    # TEA
    # =========================================================================
    # TODO: update in other HTL systems as well (the original system, HTL-PFAS)
    # based on the labor cost for the HTL plant from [16], 2014 level:
    # 1 plant manager (0.15 MM$/year)
    # 1 plant engineer (0.07 MM$/year)
    # 1 maintenance supervisor (0.06 MM$year)
    # 1 lab manager (0.06 MM$year)
    # variable cost (proportional to the sludge amount, the following is for a
    # plant of 110 dry ton [100 dry metric tonne] sludge per day):
    # 3 shift supervisors (0.14 MM$/year)
    # 1 lab technican (0.04 MM$/year)
    # 1 maintenance technician (0.04 MM$/year)
    # 4 shift operators (0.19 MM$/year)
    # 1 yard employee (0.03 MM$/year)
    # 1 clerk & secretary (0.04 MM$/year)
    # annual wage [$/year]
    wage = (0.34/labor_index[2014]*labor_index[2022]+\
            0.48/labor_index[2014]*labor_index[2022]*size*ww_2_dry_sludge_ratio/100)*10**6
    
    wage *= wage_adjustment
    
    # set income_tax_value = 0 to calculate net income
    create_tea(system=sys,
               IRR_value=0.03,
               income_tax_value=0,
               finance_interest_value=0.03,
               labor_cost_value=wage)
    
    # do not include tax credit as tax credit like 45Z is usually for final products (e.g., diesel)
    # and the allocation of the credit needs to be negociate with oil refineries
    federal_income_tax_rate_value = 0.21
    
    if state == 'US':
        # use the mode state income tax from [17]
        # from this citation: state income tax: [min, mode, max]: [0%, 6.5%, 12%]
        state_income_tax_rate_value = 0.065
    else:
        state_income_tax_rate_value = state_income_tax_rate_2022(state=state,
                                                                 sales=sys.sales,
                                                                 net_income=sys.TEA.net_earnings)
    
    income_tax_rate = federal_income_tax_rate_value + state_income_tax_rate_value
    
    create_tea(sys, IRR_value=0.03,
               income_tax_value=income_tax_rate,
               finance_interest_value=0.03,
               labor_cost_value=wage)
    
    return sys