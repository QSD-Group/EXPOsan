#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Mon Jun 5 08:46:28 2023

@author: jiananfeng

[1] Snowden-Swan, L. J.; Li, S.; Thorson, M. R.; Schmidt, A. J.; Cronin, D. J.;
    Zhu, Y.; Hart, T. R.; Santosa, D. M.; Fox, S. P.; Lemmon, T. L.; Swita, M. S.
    Wet Waste Hydrothermal Liquefaction and Biocrude Upgrading to Hydrocarbon Fuels:
    2022 State of Technology; PNNL-33622; Pacific Northwest National Lab. (PNNL),
    Richland, WA (United States), 2022. https://doi.org/10.2172/1897670.
[2] https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20)
[3] https://data.bls.gov/cgi-bin/srgate (accessed 2024-08-06)
[4] Marufuzzaman, M.; Ekşioğlu, S. D.; Hernandez, R. Truck versus Pipeline
    Transportation Cost Analysis of Wastewater Sludge. Transportation Research Part A:
    Policy and Practice 2015, 74, 14–30. https://doi.org/10.1016/j.tra.2015.02.001.
[5] https://www.sludgeprocessing.com/sludge-dewatering/sludge-drying-beds-lagoons/
    (accessed 2024-08-03)
[6] Davis, R. E.; Grundl, N. J.; Tao, L.; Biddy, M. J.; Tan, E. C.; Beckham, G. T.;
    Humbird, D.; Thompson, D. N.; Roni, M. S. Process Design and Economics for the
    Conversion of Lignocellulosic Biomass to Hydrocarbon Fuels and Coproducts:
    2018 Biochemical Design Case Update; Biochemical Deconstruction and Conversion
    of Biomass to Fuels and Products via Integrated Biorefinery Pathways;
    NREL/TP--5100-71949, 1483234; 2018; p NREL/TP--5100-71949, 1483234.
    https://doi.org/10.2172/1483234.
[7] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.;
    Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.;
    Snowden-Swan, L. J.; Davis, R.; Kinchin, C. Process Design and Economics for
    the Conversion of Algal Biomass to Hydrocarbons: Whole Algae Hydrothermal
    Liquefaction and Upgrading; PNNL--23227, 1126336; 2014; p PNNL--23227, 1126336.
    https://doi.org/10.2172/1126336.
[8] Al-Obaidani, S.; Curcio, E.; Macedonio, F.; Di Profio, G.; Al-Hinai, H.;
    Drioli, E. Potential of Membrane Distillation in Seawater Desalination:
    Thermal Efficiency, Sensitivity Study and Cost Estimation. Journal of Membrane
    Science 2008, 323 (1), 85–98. https://doi.org/10.1016/j.memsci.2008.06.006.
[9] https://businessanalytiq.com/procurementanalytics/index/ammonium-sulfate-index/
    (accessed 2024-08-03)
[10] https://www.macrotrends.net/1369/crude-oil-price-history-chart
    (accessed 2024-08-03)
[11] Pootakham, T.; Kumar, A. Bio-Oil Transport by Pipeline: A Techno-Economic
    Assessment. Bioresource Technology 2010, 101 (18), 7137–7143.
    https://doi.org/10.1016/j.biortech.2010.03.136.
[12] https://en.wikipedia.org/wiki/Methane (accessed 2024-08-03)
[13] https://www.plinovodi.si/en/transmission-system/environment-and-safety/about-natural-gas/
    (accessed 2024-08-08)
[14] https://www.eia.gov/dnav/ng/hist/n3035us3A.htm accessed (2024-08-03)
[15] Pootakham, T.; Kumar, A. A Comparison of Pipeline versus Truck Transport of
     Bio-Oil. Bioresource Technology 2010, 101 (1), 414–421.
     https://doi.org/10.1016/j.biortech.2009.07.077.
[16] Stewart, D. W.; Cortés-Peña, Y. R.; Li, Y.; Stillwell, A. S.; Khanna, M.;
     Guest, J. S. Implications of Biorefinery Policy Incentives and
     Location-SpecificEconomic Parameters for the Financial Viability of Biofuels.
     Environ. Sci. Technol. 2023. https://doi.org/10.1021/acs.est.2c07936.
[17] Snowden-Swan, L. J.; Zhu, Y.; Bearden, M. D.; Seiple, T. E.; Jones, S. B.;
     Schmidt, A. J.; Billing, J. M.; Hallen, R. T.; Hart, T. R.; Liu, J.;
     Albrecht, K. O.; Fox, S. P.; Maupin, G. D.; Elliott, D. C.
     Conceptual Biorefinery Design and Research Targeted for 2022:
     Hydrothermal Liquefacation Processing of Wet Waste to Fuels; PNNL-27186;
     Pacific Northwest National Lab. (PNNL), Richland, WA (United States), 2017.
     https://doi.org/10.2172/1415710.
'''

import qsdsan as qs, biosteam as bst, pandas as pd
from qsdsan import sanunits as qsu
from qsdsan.utils import auom, clear_lca_registries
from exposan.htl import (_load_components, create_tea, _oil_barrel_to_L,
                         state_income_tax_rate_2022, _sanunits as su)
from biosteam.units import IsenthalpicValve
from biosteam import settings

# TODO: for LCA, use ecoinvent 3.8, cutoff, TRACI, update this in the manuscript as well
# TODO: refactor the code wherever necessary
# TODO: do we want to include TEA and LCA for CT?

__all__ = ('create_geospatial_system','biocrude_density')

# kg/m3, [1]
biocrude_density = 980
# kg/m3, this is for sludge with a moisture content higher than 80%,
# google 'Design of wastewater treatment sludge thickeners Iowa State University'
sludge_density = 1000

_mile_to_km = auom('mile').conversion_factor('km')
_lb_to_kg = auom('lb').conversion_factor('kg')
_m3_to_ft3 = auom('m3').conversion_factor('ft3')
_oil_barrel_to_m3 = auom('oil_barrel').conversion_factor('m3')

# GDPCTPI (Gross Domestic Product: Chain-type Price Index), [2]
GDPCTPI = {2008: 87.977,
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

# the labor index can be found in [3] with the series id CEU3232500008,
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

def _load_process_settings(location='IL'):    
    DPO_chem = qs.Chemical(ID='DPO_chem', search_ID='101-84-8')
    DPO = qs.Component.from_chemical(ID='DPO', chemical=DPO_chem,
                                     particle_size='Soluble',
                                     degradability='Slowly',
                                     organic=True)
    
    BIP_chem = qs.Chemical(ID='BIP_chem', search_ID='92-52-4')
    BIP = qs.Component.from_chemical(ID='BIP', chemical=BIP_chem,
                                     particle_size='Soluble',
                                     degradability='Slowly',
                                     organic=True)
    
    HTF_thermo = bst.Thermo((DPO, BIP,))
    
    HTF = bst.UtilityAgent(ID='HTF', DPO=0.735, BIP=0.265, T=673.15, P=951477,
                           phase='g', thermo=HTF_thermo, regeneration_price=1)
    
    bst.HeatUtility.heating_agents.append(HTF)
    
    bst.CE = qs.CEPCI_by_year[2022]
    
    # TODO: for electricity, we aim to use balancing area for GHG; need to decide the price - does price also follows balancing area?
    folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/'
    elec = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')
    
    # TODO: note 'average' is used for the scenario 3 of sludge aggregation. Change it to the balancing area avarage value
    if location == 'average':
        bst.PowerUtility.price = elec['price (10-year median)'].mean()/100
    else:
        bst.PowerUtility.price = elec[elec['state']==location]['price (10-year median)'].iloc[0]/100

# for parameters, unless otherwise stated, refer to the original HTL system model
def create_geospatial_system(size=6.3, # MGD
                             distance=100, # km
                             sludge_transportation=0, # 0: no; 1: yes
                             # TODO: remove sludge_distance, only add it in the analysis when necessay
                              sludge_distance=100, # in km, this is the slduge transportation total distance (normalized to total sludge amount)
                             # average values below are for sludge aggregation analyses
                             average_sludge_dw_ash=None,
                             average_sludge_afdw_lipid=None,
                             average_sludge_afdw_protein=None,
                             anaerobic_digestion=0, # 0: no; 1: yes
                             aerobic_digestion=0, # 0: no; 1: yes
                             ww_2_dry_sludge_ratio=1, # dry tonne sludge/day/MGD raw wastewater
                             # TODO: replace 'state' with balancing are (balnc_area in short)
                             state='IL',
                             # TODO: use balancing area based electricity CI
                             elec_GHG=0.37 # kg CO2 eq/kWh, use state-avarage values
                             ):
    
    flowsheet_ID = 'htl_geospatial'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    _load_process_settings(location=state)
    
    # =========================================================================
    # pretreatment (Area 000)
    # =========================================================================
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)
    
    # assume the moisture content of sludge is 80% in all cases
    # for lagoon, the sludge will dry at the base of the lagoon (to an assumed 80% moisture content, see [5])
    # TODO: check this note, update if necessary
    # see 12/25/2023 notes for the compositions of different types of sludge
    if sludge_transportation == 0:
        if anaerobic_digestion == 1:
            WWTP = su.WWTP(ID='S000',
                           ins=raw_wastewater,
                           outs=('sludge','treated_water'),
                           ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                           sludge_moisture=0.8,
                           sludge_dw_ash=0.414,
                           sludge_afdw_lipid=0.193,
                           sludge_afdw_protein=0.510,
                           operation_hours=8760)
        elif aerobic_digestion == 1:
            WWTP = su.WWTP(ID='S000',
                           ins=raw_wastewater,
                           outs=('sludge','treated_water'),
                           ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                           sludge_moisture=0.8,
                           sludge_dw_ash=0.468,
                           sludge_afdw_lipid=0.193,
                           sludge_afdw_protein=0.510,
                           operation_hours=8760)
        else:
            WWTP = su.WWTP(ID='S000',
                           ins=raw_wastewater, outs=('sludge','treated_water'),
                           ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                           sludge_moisture=0.8,
                           sludge_dw_ash=0.231,
                           sludge_afdw_lipid=0.206,
                           sludge_afdw_protein=0.456,
                           operation_hours=8760)
    else:
        assert average_sludge_dw_ash != None, 'set average_sludge_dw_ash manually'
        assert average_sludge_afdw_lipid != None, 'set average_sludge_afdw_lipid manually'
        assert average_sludge_afdw_protein != None, 'set average_sludge_afdw_protein manually'
        
        WWTP = su.WWTP(ID='S000',
                       ins=raw_wastewater,
                       outs=('sludge','treated_water'),
                       ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                       sludge_moisture=0.8,
                       sludge_dw_ash=average_sludge_dw_ash,
                       sludge_afdw_lipid=average_sludge_afdw_lipid,
                       sludge_afdw_protein=average_sludge_afdw_protein,
                       operation_hours=8760)
    WWTP.register_alias('WWTP')
    
    P1 = qsu.SludgePump(ID='A100',
                        ins=WWTP-0,
                        outs='pressed_sludge',
                        P=3049.7*6894.76,
                        init_with='Stream')
    P1.include_construction = False
    P1.register_alias('P1')
    
    # TODO: only when sludge_transportation != 0
    # TODO: add sludge transportation cost in the analysis
    # assume the transported sludge has 80% moisture content
    # 4.56 $/m3, 0.072 $/m3/mile ([4], likely 2015$)
    # 1 gal water = 3.79 kg water
    # raw_wastewater.price = WWTP.ww_2_dry_sludge*\
    #     sludge_transportation*\
    #      (4.56/GDPCTPI[2015]*GDPCTPI[2022]/sludge_density*1000/0.2+\
    #       0.072/GDPCTPI[2015]*GDPCTPI[2022]/_mile_to_km/sludge_density*1000/0.2*sludge_distance)/\
    #         3.79/(10**6)
    
    # =========================================================================
    # HTL (Area 100)
    # =========================================================================
    H1 = qsu.HXutility(ID='A110',
                       include_construction=False,
                       ins=P1-0,
                       outs='heated_sludge',
                       T=351+273.15,
                       U=0.0198739,
                       init_with='Stream',
                       rigorous=True)
    H1.register_alias('H1')
    
    # assume no value/cost and no environmental benefit/impact associated with hydrochar
    HTL = qsu.HydrothermalLiquefaction(ID='A120',
                                       ins=H1-0,
                                       outs=('hydrochar','HTL_aqueous',
                                             'biocrude_to_be_stored','offgas_HTL'),
                                       mositure_adjustment_exist_in_the_system=False)
    HTL.register_alias('HTL')
    
    # =========================================================================
    # CHG (Area 200)
    # =========================================================================
    H2SO4_Tank = qsu.StorageTank(ID='T200',
                                 ins='H2SO4',
                                 outs=('H2SO4_out'),
                                 init_with='WasteStream',
                                 tau=24,
                                 vessel_material='Stainless steel')
    # 0.5 M H2SO4: ~5%
    # based on 93% H2SO4 and fresh water (dilute to 5%) prices in [6]
    H2SO4_Tank.ins[0].price = (0.043*1+0.0002*(93/5-1))/(93/5)/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    H2SO4_Tank.register_alias('H2SO4_Tank')
    
    M1 = su.HTLmixer(ID='A210',
                     ins=(HTL-1, ''),
                     outs=('mixture',))
    M1.register_alias('M1')
    
    # TODO: add the following note to other HTL systems
    # 7.8%_Ru/C and 7.8%_Ru/C_out are not valid names; they are not in sys.flowsheet.stream but the price and LCA are included
    # although this have no impact on the results, still change the names to 'virgin_CHG_catalyst and 'used_CHG_catalyst'
    CHG = qsu.CatalyticHydrothermalGasification(ID='A230',
                                                ins=(M1-0, 'virgin_CHG_catalyst'),
                                                outs=('CHG_out', 'used_CHG_catalyst'))
    # CHG price: [7]
    CHG.ins[1].price = 60/_lb_to_kg/GDPCTPI[2011]*GDPCTPI[2022]
    CHG.register_alias('CHG')
    
    V1 = IsenthalpicValve(ID='A240',
                          ins=CHG-0,
                          outs='depressed_cooled_CHG',
                          P=50*6894.76,
                          vle=True)
    V1.register_alias('V1')
    
    F1 = qsu.Flash(ID='A250',
                   ins=V1-0,
                   outs=('CHG_fuel_gas','N_riched_aqueous'),
                   T=60+273.15,
                   P=50*6894.76,
                   thermo=settings.thermo.ideal())
    F1.register_alias('F1')
    F1.include_construction = False
    
    # assume no value/cost and no environmental benefit/impact associated with MemDis_ww (the system is in a WWTP) and solution
    MemDis = qsu.MembraneDistillation(ID='A260',
                                      ins=(F1-1, H2SO4_Tank-0, 'NaOH', 'Membrane_in'),
                                      outs=('ammonium_sulfate','MemDis_ww','Membrane_out','solution'),
                                      init_with='WasteStream')
    MemDis.ins[2].price = 0.2384/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    # RO membrane price: [8] (likely 2008$)
    MemDis.ins[3].price = 90/GDPCTPI[2008]*GDPCTPI[2022]
    # ammonium sulfate (2022 average): [9]
    MemDis.outs[0].price = 0.41
    MemDis.include_construction = False
    MemDis.register_alias('MemDis')
    
    # =========================================================================
    # Storage, and disposal (Area 300)
    # =========================================================================
    BiocrudeTank = qsu.StorageTank(ID='T300',
                                   ins=HTL-2,
                                   outs=('biocrude'),
                                   tau=3*24,
                                   init_with='WasteStream',
                                   vessel_material='Carbon steel')
    # TODO: in the main manuscript or the SI, specifically say we displace crude oil by biocrude
    # TODO: (compromise a little bit, like saying this blending may need further research, with citations, but also saying it is not impossible, also with citations)
    # assume the biocrude has the same price as crude oil (for LCA: we assume the crude oil can be displaced by the biocrude)
    # 2022 average closing price for crude oil: 94.53 $/oil barrel, [10]
    # transportation cost: 5.67 2008$/m3 (fixed cost) and 0.07 2008$/m3/km (variable cost), [11]
    BiocrudeTank.outs[0].price = 94.53/_oil_barrel_to_m3/biocrude_density-\
        (5.67/biocrude_density+0.07/biocrude_density*distance)/GDPCTPI[2008]*GDPCTPI[2022]
    # store for 3 days based on [7]
    BiocrudeTank.register_alias('BiocrudeTank')
    
    PC1 = qsu.PhaseChanger(ID='S300',
                           ins=CHG-1,
                           outs='CHG_catalyst_out',
                           phase='s')
    PC1.register_alias('PC1')
    
    GasMixer = qsu.Mixer(ID='S310',
                         ins=(HTL-3, F1-0),
                         outs=('fuel_gas'),
                         init_with='Stream')
    GasMixer.register_alias('GasMixer')
    
    # =========================================================================
    # facilities
    # =========================================================================
    qsu.HeatExchangerNetwork(ID='HXN',
                             T_min_app=86,
                             force_ideal_thermo=True)
    
    # assume no value/cost and no environmental benefit/impact associated with emission (they are not treated and are biogenic; only CO2 from natural gas is non-biogenic but we have included the environmental impact for it)
    # use CHP to produce electricity does not provide benefit; therefore, set supplement_power_utility=False
    CHP = qsu.CombinedHeatPower(ID='CHP',
                                include_construction=False,
                                ins=(GasMixer-0, 'natural_gas', 'air'),
                                outs=('emission','solid_ash'),
                                init_with='WasteStream',
                                supplement_power_utility=False)
    # assume natural gas is CH4 (density is 0.657 kg/m3, [12])
    # the density of natural gas (0.68 kg/m3, [13]) is larger than that of CH4, therefore, the calculation here is conservative
    # 2022 US NG industrial price: 7.66 $/1000 ft3, [14]
    CHP.ins[1].price = 7.66/1000*_m3_to_ft3/0.657
    # 1.41 MM 2016$/year for 4270/4279 kg/hr ash, 7880 annual operating hours, from [6]
    CHP.outs[1].price = -1.41*10**6/7880/4270/GDPCTPI[2016]*GDPCTPI[2022]
    CHP.register_alias('CHP')
    
    # not sure if the construction cost for CT is right or not (based on return cooling water flow rate or makeup water flow rate or others)
    # however, the effect of this is small (2%) and can be ignored
    # the cost and CI of CT streams are also small and can be ignored
    # note there is actually a default 3 $/kg for cooling_tower_chemicals
    # we change the price to 1.7842 $/lb from [6]
    # but the flow rate is so small, so there is actually no significant difference
    # TODO: consider adding CT and its TEA and LCA items (if any) for other systems (HTL, HTL-PFAS)
    CT = bst.facilities.CoolingTower(ID='CT')
    CT.ins[2].price = 1.7842/_lb_to_kg/GDPCTPI[2016]*GDPCTPI[2022]
    CT.register_alias('CT')
    
    sys = qs.System.from_units(ID='sys_geospatial',
                               units=list(flowsheet.unit),
                               operating_hours=WWTP.operation_hours)
    sys.register_alias('sys')

    # =========================================================================
    # add stream impact items
    # =========================================================================
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='TRACI',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential.')
    
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.48748859)
    
    # TODO: update parameters for Sludge_trucking
    Sludge_trucking = qs.ImpactItem('Sludge_trucking', functional_unit='kg*km')
    # TODO: update GWP and cost
    Sludge_trucking.add_indicator(GlobalWarming, WWTP.ww_2_dry_sludge*\
                      1*\
                      0.13004958/0.2/3.79/(10**6))
    # TODO: transportation was not included in TEA
    Sludge_trucking.price = 1000000000
    
    impact_items = {# TODO: only when sludge_transportation != 0
                    # TODO: add sludge transportation GHG in the analysis
                    # assume the transported sludge has 80% moisture content
                    # 1 gal water = 3.79 kg water
                    # 'market for transport, freight, lorry, unspecified' (0.13004958 kg CO2 eq/metric ton/km)
                    # assume the sludge/biosolids decomposition is minimal during transportation since our transportation distance is not long
                    # 'sludge':       [stream.raw_wastewater,
                    #                  WWTP.ww_2_dry_sludge*\
                    #                  sludge_transportation*\
                    #                  0.13004958/0.2*sludge_distance/3.79/(10**6)],
                    'CHG_catalyst': [stream.CHG_catalyst_out, 471.098936962268],
                    'H2SO4':        [stream.H2SO4, 0.005529872568],
                    'NaOH':         [stream.NaOH, 1.2497984],
                    'RO_membrane':  [stream.Membrane_in, 2.2709135],
                    # TODO: address this problem (i.e., do not use market or market group for products) for other systems (i.e., CO2 sorbent, HTL-PFAS)
                    # use 'ammonium sulfate production' for 'NH42SO4' instead of market or market group since we don't offset things like transportation
                    'NH42SO4':      [stream.ammonium_sulfate, -1.1139067],
                    'natural_gas':  [stream.natural_gas, 2.5199928],
                    # use market or market group for biocrude since we want to offset transportation and then add our own transportation part
                    # 0.22290007 kg CO2 eq/kg petroleum ('market for petroleum')
                    # TODO: we don't use the biocrude transportation CI from the literature; instead, we use 'market for transport, freight, lorry, unspecified' (0.13004958 kg CO2 eq/metric ton/km) from ecoinvent instead; update this in the manuscript 
                    # 89 g CO2/m3/km: carbon intensity of truck transportation, [15]
                    'biocrude':     [stream.biocrude, 0.13004958/1000*distance-0.22290007],
                    'ash_disposal': [stream.solid_ash, 0.0082744841]}
    
    for item in impact_items.items():
        qs.StreamImpactItem(ID=item[0], linked_stream=item[1][0], GlobalWarming=item[1][1])
    
    # TODO: do we also want to separate biocrude transportation LCA and TEA?
    # TODO: update parameters for Sludge_transportation
    Sludge_transportation = qs.Transportation('Sludge_trucking', linked_unit=WWTP, item=Sludge_trucking,
                                              load_type='mass', load=stream.raw_wastewater.F_mass, load_unit='kg',
                                              distance=sludge_distance, distance_unit='km',
                                              interval='1', interval_unit='h')
    WWTP.transportation = Sludge_transportation
    
    # simulate first to enable the calculation of the income tax rate
    sys.simulate()
    
    # TODO: do we want to include tax credit?
    federal_income_tax_rate_value = 0.21
    
    # TODO: ignore annualized depreciation, add this in the manuscript or in the SI
    
    # breakpoint()
    
    # TODO: check here
    # consistent with the TEA parameters
    wage = (0.34/labor_index[2014]*labor_index[2022]+\
                      0.48/labor_index[2014]*labor_index[2022]*size*ww_2_dry_sludge_ratio/100)*\
                     10**6
    
    TIC = sys.installed_cost
    
    TDC = TIC*(1 + 0.04+0.09+0.045)
    
    FCI = TDC * (1+0.1+0.2+0.4+0.1+0.1)
    
    FOC = wage + 0.9*wage + 0.03*TIC + 0.007*FCI
    
    VOC = sys.material_cost + sys.utility_cost
    
    annual_net_income = sys.sales - FOC - VOC
    
    # TODO: make sure BioSTEAM has already consider the situration when annual_net_income is less than 0
    assert annual_net_income
    
    # TODO: mention this in the main manuscript or the SI and add citation
    if state == 'average':
        # use the mode state income tax from [16]
        # from this citation: state income tax: [min, mode, max]: [0%, 6.5%, 12%]
        state_income_tax_rate_value = 0.065
    else:
        state_income_tax_rate_value = state_income_tax_rate_2022(state=state,
                                                                 sales=sys.sales,
                                                                 net_income=annual_net_income)
    
    income_tax_rate = federal_income_tax_rate_value + state_income_tax_rate_value
    
    # TODO: update in other HTL systems as well (the original system, HTL-PFAS)
    # based on the labor cost for the HTL plant from [17], 2014 level:
    # 1 plant manager (0.15 MM$/year)
    # 1 plant engineer (0.07 MM$/year)
    # 1 maintenance supervisor (0.06 MM$year)
    # 1 lab manager (0.06 MM$year)
    # variable cost (proportional to the sludge amount, the following is for a plant of 110 dry ton [100 dry metric tonne] sludge per day)
    # 3 shift supervisors (0.14 MM$/year)
    # 1 lab technican (0.04 MM$/year)
    # 1 maintenance technician (0.04 MM$/year)
    # 4 shift operators (0.19 MM$/year)
    # 1 yard employee (0.03 MM$/year)
    # 1 clerk & secretary (0.04 MM$/year)
    create_tea(sys, IRR_value=0.03,
               income_tax_value=income_tax_rate,
               finance_interest_value=0.03,
               labor_cost_value=wage)
    
    assert sys.get_heating_duty() == 0
    assert sys.get_cooling_duty() == 0
    
    # 0.48748859 is the GHG level with the Electricity item from ecoinvent,
    # we cannot list electricity GHG one state by one state,
    # but we can adjust the electricity amount to reflect different GHG of electricity at different states
    qs.LCA(system=sys, lifetime=30, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30/0.48748859*elec_GHG)
    
    # biocrude production in BPD (barrel per day)
    biocrude_barrel = BiocrudeTank.outs[0].F_mass/biocrude_density*1000/_oil_barrel_to_L*24
    
    return sys, biocrude_barrel