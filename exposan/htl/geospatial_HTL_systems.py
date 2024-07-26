#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Mon Jun 5 08:46:28 2023

@author: jiananfeng
'''

import os, qsdsan as qs, biosteam as bst, pandas as pd
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
# TODO: we are using 2022$, should we use 2022 tax rate for state_income_tax_rate?
from exposan.htl import _load_components, create_tea, _oil_barrel_to_L, state_income_tax_rate, _sanunits as su
from biosteam.units import IsenthalpicValve
from biosteam import settings

# TODO: change all TEA and LCA data to the U.S.-based values
# if impossible, especially for LCA: use US, then RER, then RoW, then GLO
# also, remember to use LCA data collected for the 'cutoff' model and 'IPCC' method
# TODO: confirm IPCC 2013 (no LT) is ok

__all__ = ('create_geospatial_system','biocrude_density')

biocrude_density = 980 # kg/m3, Snowden-Swan et al. 2022 SOT, PNNL
# TODO: sludge_density is not used, do we still want to keep it here?
sludge_density = 1000 # kg/m3, this is for sludge with a moisture content higher than 80%, google 'Design of wastewater treatment sludge thickeners Iowa State University'

# GDPCTPI (Gross Domestic Product: Chain-type Price Index)
# https://fred.stlouisfed.org/series/GDPCTPI (accessed 2024-05-20)
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
    
    # TODO: update all costs to 2022$
    bst.CE = qs.CEPCI_by_year[2022]
    
    # TODO: for electricity, we aim to use balancing area for GHG. Need to decide the price - does price also follows balancing area?
    folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/'
    elec = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')
    
    # TODO: note 'average' is used for the scenario 3 of sludge aggregation. Change it to the balancing area avarage value
    if location == 'average':
        bst.PowerUtility.price = elec['price (10-year median)'].mean()/100
    else:
        bst.PowerUtility.price = elec[elec['state']==location]['price (10-year median)'].iloc[0]/100

# TODO: check the units of waste_cost and waste_GHG
# for unit parameters, unless otherwise states, refer to the original HTL system model
def create_geospatial_system(waste_cost=450, # based on the share of sludge management methods of each facilities
                             waste_GHG=200, # based on the share of sludge management methods of each facilities
                             size=8, # MGD
                             distance=100, # km
                             sludge_transportation=0, # 0: no; 1: yes
                             # TODO: pay attention to 'normalized to total sludge amount'
                             sludge_distance=0, # in km, this is the slduge transportation total distance (normalized to total sludge amount)
                             # TODO: average values below are for sludge aggregation analyses
                             average_sludge_dw_ash=None,
                             average_sludge_afdw_lipid=None,
                             average_sludge_afdw_protein=None,
                             anaerobic_digestion=0, # 0: no; 1: yes
                             aerobic_digestion=0, # 0: no; 1: yes
                             ww_2_dry_sludge_ratio=1, # tonne sludge/MGD raw wastewater
                             # TODO: replace 'state' with balancing are (balnc_area in short)
                             state='IL',
                             # TODO: use balancing area based electricity CI
                             elec_GHG=0.37 # use state-avarage values
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
    
    folder = os.path.dirname(__file__)
    qs.ImpactIndicator.load_from_file(os.path.join(folder, 'data/geo_impact_indicators.csv'))
    # TODO: update LCA items to be US-based, see comments on the top
    qs.ImpactItem.load_from_file(os.path.join(folder, 'data/geo_impact_items.xlsx'))
    
    # =============================================================================
    # pretreatment (Area 000)
    # =============================================================================
    # raw wastewater into a WRRF, in MGD
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)
    
    # TODO: check this note
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
    # TODO: include construction for cost only, but excluding it from LCA (although it is not a big part)
    P1.include_construction = True
    P1.register_alias('P1')
    
    # TODO: check the units and calculations, also update it to 2022$ if necessary
    # TODO: adjust the format
    # 1 gal water = 3.79 kg water
    # 5.06 $/m3, 0.05 $/m3/km (Marufuzzaman et al. Transportation Research Part A. 2015, converted to 2020 price)
    raw_wastewater.price = -WWTP.ww_2_dry_sludge*(waste_cost-sludge_transportation*(5.06*5+0.05*5*sludge_distance))/3.79/(10**6)
    
    # =============================================================================
    # HTL (Area 100)
    # =============================================================================
    H1 = qsu.HXutility(ID='A110',
                       include_construction=True,
                       ins=P1-0,
                       outs='heated_sludge',
                       T=351+273.15,
                       U=0.0198739,
                       init_with='Stream',
                       rigorous=True)
    H1.register_alias('H1')
    
    # TODO: what is the implication of this?
    # !!! we assume we don't need to adjust moisture content (all = 80%)
    # for lagoon, the sludge will dry at the base of the lagoon (see https://www.sludgeprocessing.com/sludge-dewatering/sludge-drying-beds-lagoons/)
    
    HTL = qsu.HydrothermalLiquefaction(ID='A120',
                                       ins=H1-0,
                                       outs=('hydrochar','HTL_aqueous','biocrude_to_be_stored','offgas_HTL'),
                                       mositure_adjustment_exist_in_the_system=False)
    HTL.register_alias('HTL')
    
    # =============================================================================
    # CHG (Area 200)
    # =============================================================================
    H2SO4_Tank = qsu.StorageTank(ID='T200',
                                 ins='H2SO4',
                                 outs=('H2SO4_out'),
                                 init_with='WasteStream',
                                 tau=24,
                                 vessel_material='Stainless steel')
    # TODO: update the cost to 2022$/kg; find a better price if possible
    # based on 93% H2SO4 and fresh water (dilute to 5%) price found in Davis 2020$/kg
    H2SO4_Tank.ins[0].price = 0.00658
    H2SO4_Tank.register_alias('H2SO4_Tank')
    
    SP1 = qsu.ReversedSplitter(ID='S200',
                               ins=H2SO4_Tank-0,
                               outs=('H2SO4_P','H2SO4_N'),
                               init_with='Stream')
    SP1.register_alias('SP1')

    # TODO: what is this? Can we just put '' in M1
    M1_outs1 = ''
    
    M1 = su.HTLmixer(ID='A210',
                     ins=(HTL-1, M1_outs1),
                     outs=('mixture',))
    M1.register_alias('M1')
    
    CHG = qsu.CatalyticHydrothermalGasification(ID='A230',
                                                ins=(M1-0, '7.8%_Ru/C'),
                                                outs=('CHG_out', '7.8%_Ru/C_out'))
    # TODO: update the cost to 2022$/kg
    CHG.ins[1].price = 134.53
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
    F1.include_construction = True
    
    MemDis = qsu.MembraneDistillation(ID='A260',
                                      ins=(F1-1, SP1-1, 'NaOH', 'Membrane_in'),
                                      outs=('ammonium_sulfate','MemDis_ww', 'Membrane_out','solution'),
                                      init_with='WasteStream')
    # TODO: update the cost to 2022$/kg
    # TODO: is membrane considered as part of construction? (tend to consider it as material instead)
    MemDis.ins[2].price = 0.5256
    MemDis.outs[0].price = 0.3236
    MemDis.register_alias('MemDis')
    
    # =============================================================================
    # Storage, and disposal (Area 300)
    # =============================================================================
    BiocrudeTank = qsu.StorageTank(ID='T300',
                                   ins=HTL-2,
                                   outs=('biocrude'),
                                   tau=3*24,
                                   init_with='WasteStream',
                                   vessel_material='Carbon steel')
    # store for 3 days based on Jones 2014
    BiocrudeTank.register_alias('CrudeOilTank')
    
    # TODO: update the price to 2022$ using GDPCTPI
    # check the units and calculations here
    # 5.67: fixed cost, 0.07: variable cost ($/m3, Pootakham et al. Bio-oil transport by pipeline: A techno-economic assessment. Bioresource Technology, 2010)
    # note cost in this paper was in 2008$, we need to convert it to 2020$ based on Gross Domestic Product chain-type price index
    # 1 2008$ = 1.20 $2020 based on https://fred.stlouisfed.org/series/GDPCTPI
    BiocrudeTank.outs[0].price = -5.67*1.20/biocrude_density-0.07*1.20/biocrude_density*distance+0.3847
    
    PC1 = qsu.PhaseChanger(ID='S300',
                           ins=CHG-1,
                           outs='CHG_catalyst_out',
                           phase='s')
    PC1.register_alias('PC1')
    
    GasMixer = qsu.Mixer(ID='S310',
                         ins=(HTL-3, F1-0,),
                         outs=('fuel_gas'),
                         init_with='Stream')
    GasMixer.register_alias('GasMixer')
    
    # =============================================================================
    # facilities
    # =============================================================================
    qsu.HeatExchangerNetwork(ID='HXN',
                             T_min_app=86,
                             force_ideal_thermo=True)
    
    # TODO: compare with the CHP in the original HTL system, why supplement_power_utility=False?
    CHP = qsu.CombinedHeatPower(ID='CHP',
                                include_construction=True,
                                ins=(GasMixer-0, 'natural_gas', 'air'),
                                outs=('emission','solid_ash'),
                                init_with='WasteStream',
                                supplement_power_utility=False)
    # TODO: make sure the signs are correct for these two items
    CHP.ins[1].price = bst.stream_prices['Natural gas']
    CHP.outs[1].price = bst.stream_prices['Ash disposal']
    
    # TODO: add LCA for the CT
    # TODO: consider adding CT for the system in other places (original model, HTL-PFAS)
    bst.facilities.CoolingTower(ID='CT')
    
    sys = qs.System.from_units(ID='sys_geospatial',
                               units=list(flowsheet.unit),
                               operating_hours=WWTP.operation_hours)
    sys.register_alias('sys')

    # =============================================================================
    # add stream impact items
    # =============================================================================
    # TODO: check the notes in impact_items
    impact_items = {# 1 gal water = 3.79 kg water
                    # 0.719 kg CO2 eq/dry tonne/km (Zhao et al. Journal of Environmental Sciences. 2023)
                    'sludge':       [stream.sludge_assumed_in_wastewater,
                                     -WWTP.ww_2_dry_sludge*\
                                     (waste_GHG-sludge_transportation*\
                                     (0.719*sludge_distance))/3.79/(10**6)],
                    'CHG_catalyst': [stream.CHG_catalyst_out, 484.7862509],
                    'H2SO4':        [stream.H2SO4, 0.008205666],
                    'NaOH':         [stream.NaOH, 1.2514],
                    'RO':           [stream.Membrane_in, 2.2663],
                    'natural_gas':  [stream.natural_gas, 1.584234288],
                    'NH42SO4':      [stream.ammonium_sulfate, -1.2499],
                    # crude oil/petroleum (transportation is included in crude oil item, we offset that first, but we need to add our own transportation)
                    # 89 g CO2/m3/km: carbon intensity of truck transportation (Pootakham et al. A comparison of pipeline versus truck transport of bio-oil. Bioresource Technology, 2010)
                    # others are scaled based on transportation data from EcoInvent and GWP from the paper above
                    'biocrude':     [stream.biocrude, 89/1000/biocrude_density*distance-0.22304]}
    
    for item in impact_items.items():
        qs.StreamImpactItem(ID=item[0], linked_stream=item[1][0], GlobalWarming=item[1][1])
    
    # TODO: make sure doing the following does not affect the results
    # simulate first to enable the calculation of the income tax rate
    sys.simulate()
    
    # TODO: do we want to include tax credit?
    federal_income_tax_rate_value = 0.21
    
    # TODO: check the calculation here
    annual_sales = BiocrudeTank.outs[0].F_mass*365*24*0.3847 + MemDis.outs[0].F_mass*365*24*0.3236
    annual_material_cost = sum([s.cost for s in sys.feeds[1:] if ((s.price > 0) & (s.ID != 'sludge_assumed_in_wastewater'))]) * sys.operating_hours
    annual_utility_cost = sum([u.utility_cost for u in sys.cost_units]) * sys.operating_hours
    annual_net_income = annual_sales - annual_material_cost - annual_utility_cost
    
    if state == 'average':
        # TODO: check this citation
        # use the mode state income tax from Steward et al. ES&T 2023
        state_income_tax_rate_value = 0.065
    else:
        state_income_tax_rate_value = state_income_tax_rate(state=state, sales=annual_sales, net_income=annual_net_income)
    
    income_tax_rate = federal_income_tax_rate_value + state_income_tax_rate_value
    
    # TODO: check the parameter values for TEA, especially labor costs
    # TODO: update in other HTL systems as well (the original system, HTL-PFAS)
    # TODO: do we want to use 2022 salary level?
    # for labor cost (2020 salary level)

    # 1 plant manager
    # 1 plant engineer
    # 1 maintenance supervisor
    # 1 lab manager
    # variable cost (proportional to the sludge amount, the following is for a plant of 100 dry metric tonne per day)
    # 3 shift supervisors
    # 1 lab technican
    # 1 maintenance technician
    # 4 shift operators
    # 1 yard employee
    # 1 clerk & secretary
    
    # the labor index can be found in https://data.bls.gov/cgi-bin/srgate with the series id CEU3232500008, remember to select 'include annual average'
    # the labor cost would be considered as the same for both the systems in the HTL model paper (including hydroprocessing and struvite recovery) and in the HTL geospatial paper
    # this is reasonable based on the labor cost data from
    # Snowden-Swan et al. 2017. Conceptual Biorefinery Design and Research Targeted for 2022: Hydrothermal Liquefaction Processing of Wet Waste to Fuels
    # (having separate tables for the HTL plant and the hydroprocessing plant serving 10 HTL plants) and
    # Jones et al. 2014. Process Design and Economics for the Conversion of Algal Biomass to Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading
    # (including hydroprocessing and CHG)
    create_tea(sys, IRR_value=0.03,
               income_tax_value=income_tax_rate,
               finance_interest_value=0.03,
               labor_cost_value=(0.41+0.59*size*ww_2_dry_sludge_ratio/100)*10**6)
    
    # TODO: since we have CT, do we still need Cooling here?
    # TODO: check the following lines of notes
    # 0.67848 is the GHG level with the Electricity item from ecoinvent,
    # we cannot list electricity GHG one state by one state,
    # but we can adjust the electricity amount to reflect different GHG of electricity at different states
    qs.LCA(system=sys,
           lifetime=30,
           lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30/0.67848*elec_GHG,
           Cooling=lambda:sys.get_cooling_duty()/1000*30)
    
    # biocrude production in BPD (barrel per day)
    biocrude_barrel = BiocrudeTank.outs[0].F_mass/biocrude_density*1000/_oil_barrel_to_L*24
    
    return sys, biocrude_barrel