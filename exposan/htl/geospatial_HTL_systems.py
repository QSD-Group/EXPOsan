#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Mon Jun 5 08:46:28 2023

@author: jiananfeng

References:
    
(1) Lang, C.; Lee, B. Heat Transfer Fluid Life Time Analysis of Diphenyl
    Oxide/Biphenyl Grades for Concentrated Solar Power Plants. Energy Procedia
    2015, 69, 672–680. https://doi.org/10.1016/j.egypro.2015.03.077.
'''

import os, qsdsan as qs, biosteam as bst, pandas as pd
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import _load_components, create_tea, _oil_barrel_to_L, state_income_tax_rate, _sanunits as su
from biosteam.units import IsenthalpicValve
from biosteam import settings

# TODO: change all TEA and LCA data to the U.S.-based values
# if impossible, especially for LCA: use US, then RER, then RoW, then GLO
# also, remember to use LCA data collected for the 'cutoff' model and 'IPCC' method
# TODO: confirm IPCC 2013 (no LT) is ok

__all__ = ('create_geospatial_system','biocrude_density')

biocrude_density = 980 # kg/m3, Snowden-Swan et al. 2022 SOT, PNNL
sludge_density = 1000 # kg/m3, this is for sludge with a moisture content higher than 80%, google 'Design of wastewater treatment sludge thickeners Iowa State University'

def _load_process_settings(location='IL'):

    # see _process_settings.py in htl module for more information
    DPO_chem = qs.Chemical('DPO_chem', search_ID='101-84-8')
    BIP_chem = qs.Chemical('BIP_chem', search_ID='92-52-4')
    
    DPO = qs.Component.from_chemical('DPO', chemical=DPO_chem, particle_size='Soluble',
                                     degradability='Slowly', organic=True)
    
    BIP = qs.Component.from_chemical('BIP', chemical=BIP_chem, particle_size='Soluble',
                                     degradability='Slowly', organic=True)
    
    HTF_thermo = bst.Thermo((DPO, BIP,))
    
    HTF = bst.UtilityAgent('HTF', DPO=0.735, BIP=0.265, T=673.15, P=951477,
                           phase='g', thermo=HTF_thermo, regeneration_price=1)
    
    bst.HeatUtility.heating_agents.append(HTF)
    
    # use 2020$ to match up with latest PNNL report
    bst.CE = qs.CEPCI_by_year[2020]
    
    folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/'
    
    elec = pd.read_excel(folder + 'state_elec_price_GHG.xlsx', 'summary')
    
    if location == 'average':
        bst.PowerUtility.price = elec['price (10-year median)'].mean()/100
    else:
        bst.PowerUtility.price = elec[elec['state']==location]['price (10-year median)'].iloc[0]/100

def create_geospatial_system(waste_cost=450, # based on the share of sludge management methods of each facilities
                             waste_GHG=200, # based on the share of sludge management methods of each facilities
                             size=8, # in MGD
                             distance=100, # in km, using Google Maps API
                             sludge_transportation=0, # if 0, there is no sludge transportation; if 1, there is sludge transportation
                             sludge_distance=0, # in km, this is the slduge transportation total distance (normalized to total sludge amount)
                             average_sludge_dw_ash=None,
                             average_sludge_afdw_lipid=None,
                             average_sludge_afdw_protein=None,
                             anaerobic_digestion=0,
                             aerobic_digestion=0,
                             ww_2_dry_sludge_ratio=1,
                             state='IL',
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
    
    folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/'
    qs.ImpactIndicator.load_from_file(os.path.join(folder, 'data/impact_indicators.csv'))
    qs.ImpactItem.load_from_file(os.path.join(folder, 'data/impact_items.xlsx'))
    
    raw_wastewater = qs.WasteStream('sludge_assumed_in_wastewater', H2O=size, units='MGD', T=25+273.15)
    # set H2O equal to the total raw wastewater into the WWTP
    
    # =============================================================================
    # pretreatment (Area 000)
    # =============================================================================
    # see 12/25/2023 notes for the compositions of different types of sludge
    if sludge_transportation == 0:
        if anaerobic_digestion == 1:
            WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                           ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                           # how much metric tonne/day sludge can be produced by 1 MGD of ww
                           sludge_moisture=0.8,
                           sludge_dw_ash=0.414,
                           sludge_afdw_lipid=0.193,
                           sludge_afdw_protein=0.510,
                           operation_hours=8760)
            
        elif aerobic_digestion == 1:
            WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                           ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                           # how much metric tonne/day sludge can be produced by 1 MGD of ww
                           sludge_moisture=0.8,
                           sludge_dw_ash=0.468,
                           sludge_afdw_lipid=0.193,
                           sludge_afdw_protein=0.510,
                           operation_hours=8760)
            
        else:
            WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                           ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                           # how much metric tonne/day sludge can be produced by 1 MGD of ww
                           sludge_moisture=0.8,
                           sludge_dw_ash=0.231,
                           sludge_afdw_lipid=0.206,
                           sludge_afdw_protein=0.456,
                           operation_hours=8760)
            
    else:
        WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                       ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                       # how much metric tonne/day sludge can be produced by 1 MGD of ww
                       sludge_moisture=0.8,
                       sludge_dw_ash=average_sludge_dw_ash,
                       sludge_afdw_lipid=average_sludge_afdw_lipid,
                       sludge_afdw_protein=average_sludge_afdw_protein,
                       operation_hours=8760)
        
    P1 = qsu.SludgePump('A100', ins=WWTP-0, outs='pressed_sludge', P=3049.7*6894.76,
              init_with='Stream')
    P1.include_construction = True
        
    WWTP.register_alias('WWTP')
    
    P1.register_alias('P1')
    # Jones 2014: 3049.7 psia
    
    raw_wastewater.price = -WWTP.ww_2_dry_sludge*(waste_cost - sludge_transportation*(5.06*5 + 0.05*5*sludge_distance))/3.79/(10**6)
    # 1 gal water = 3.79 kg water
    # 5.06 $/m3, 0.05 $/m3/km (Marufuzzaman et al. Transportation Research Part A. 2015, converted to 2020 price)

    # =============================================================================
    # HTL (Area 100)
    # =============================================================================
    
    H1 = qsu.HXutility('A110', include_construction=True,
                       ins=P1-0, outs='heated_sludge', T=351+273.15,
                       U=0.0198739, init_with='Stream', rigorous=True)
    # feed T is low, thus high viscosity and low U (case B in Knorr 2013)
    # U: 3, 3.5, 4 BTU/hr/ft2/F as minimum, baseline, and maximum
    # U: 0.0170348, 0.0198739, 0.0227131 kW/m2/K
    # but not in other heat exchangers (low viscosity, don't need U to enforce total heat transfer efficiency)
    # unit conversion: https://www.unitsconverters.com/en/Btu(It)/Hmft2mdegf-To-W/M2mk/Utu-4404-4398
    
    H1.register_alias('H1')
    
    # !!! we assume we don't need to adjust moisture content (all = 80%)
    # for lagoon, the sludge will dry at the base of the lagoon (see https://www.sludgeprocessing.com/sludge-dewatering/sludge-drying-beds-lagoons/)
    
    HTL = qsu.HydrothermalLiquefaction('A120', ins=H1-0,
                                       outs=('hydrochar','HTL_aqueous','biocrude_to_be_stored','offgas_HTL'),
                                       mositure_adjustment_exist_in_the_system=False)
    HTL.register_alias('HTL')
    
    # =============================================================================
    # CHG (Area 200)
    # =============================================================================
    
    H2SO4_Tank = qsu.StorageTank('T200', ins='H2SO4', outs=('H2SO4_out'),
                             init_with='WasteStream', tau=24, vessel_material='Stainless steel')
    H2SO4_Tank.ins[0].price = 0.00658 # based on 93% H2SO4 and fresh water (dilute to 5%) price found in Davis 2020$/kg
    H2SO4_Tank.register_alias('H2SO4_Tank')
    
    SP1 = qsu.ReversedSplitter('S200',ins=H2SO4_Tank-0, outs=('H2SO4_P','H2SO4_N'),
                               init_with='Stream')
    SP1.register_alias('SP1')

    M1_outs1 = ''
    
    M1 = su.HTLmixer('A210', ins=(HTL-1, M1_outs1), outs=('mixture',))
    M1.register_alias('M1')
    
    CHG = qsu.CatalyticHydrothermalGasification('A230',
                                                ins=(M1-0, '7.8%_Ru/C'),
                                                outs=('CHG_out', '7.8%_Ru/C_out'))
    CHG.ins[1].price = 134.53
    CHG.register_alias('CHG')
    
    V1 = IsenthalpicValve('A240', ins=CHG-0, outs='depressed_cooled_CHG', P=50*6894.76, vle=True)
    V1.register_alias('V1')
    
    F1 = qsu.Flash('A250', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                      T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    F1.register_alias('F1')
    F1.include_construction = True
    
    MemDis = qsu.MembraneDistillation('A260', ins=(F1-1, SP1-1, 'NaOH', 'Membrane_in'),
                                  outs=('ammonium_sulfate','MemDis_ww', 'Membrane_out','solution'),
                                  init_with='WasteStream')
    MemDis.ins[2].price = 0.5256
    MemDis.outs[0].price = 0.3236
    MemDis.register_alias('MemDis')
    
    # =============================================================================
    # Storage, and disposal (Area 300)
    # =============================================================================
    
    BiocrudeTank = qsu.StorageTank('T300', ins=HTL-2, outs=('biocrude'),
                                    tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')
    # store for 3 days based on Jones 2014
    BiocrudeTank.register_alias('CrudeOilTank')
    
    BiocrudeTank.outs[0].price = -5.67*1.20/biocrude_density - 0.07*1.20/biocrude_density*distance + 0.3847
    # 5.67: fixed cost, 0.07: variable cost ($/m3, Pootakham et al. Bio-oil transport by pipeline: A techno-economic assessment. Bioresource Technology, 2010)
    # note cost in this paper was in 2008$, we need to convert it to 2020$ based on Gross Domestic Product chain-type price index
    # 1 2008$ = 1.20 $2020 based on https://fred.stlouisfed.org/series/GDPCTPI
    
    PC1 = qsu.PhaseChanger('S300', ins=CHG-1, outs='CHG_catalyst_out', phase='s')
    PC1.register_alias('PC1')
    
    GasMixer = qsu.Mixer('S310', ins=(HTL-3, F1-0,),
                          outs=('fuel_gas'), init_with='Stream')
    GasMixer.register_alias('GasMixer')
    
    # =============================================================================
    # facilities
    # =============================================================================
    
    qsu.HeatExchangerNetwork('HXN', T_min_app=86, force_ideal_thermo=True)
    # 86 K: Jones et al. PNNL, 2014
    
    CHP = qsu.CombinedHeatPower('CHP', include_construction=True,
                                ins=(GasMixer-0, 'natural_gas', 'air'),
                                outs=('emission','solid_ash'), init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.ins[1].price = 0.1685
    
    bst.facilities.CoolingTower('CT')
    
    sys = qs.System.from_units('sys_geospatial',
                               units=list(flowsheet.unit),
                               operating_hours=WWTP.operation_hours)
    sys.register_alias('sys')

# =============================================================================
# add stream impact items
# =============================================================================
   
    # only GlobalWarming can be used, since the values for other CFs are 0 for the sludge_item

    # add impact for waste sludge
    qs.StreamImpactItem(ID='sludge_in_wastewater_item',
                        linked_stream=stream.sludge_assumed_in_wastewater,
                        Acidification=0,
                        Ecotoxicity=0,
                        Eutrophication=0,
                        GlobalWarming=-WWTP.ww_2_dry_sludge*(waste_GHG - sludge_transportation*(0.719*sludge_distance))/3.79/(10**6),
                        # 1 gal water = 3.79 kg water
                        # 0.719 kg CO2 eq/dry tonne/km (Zhao et al. Journal of Environmental Sciences. 2023)
                        OzoneDepletion=0,
                        PhotochemicalOxidation=0,
                        Carcinogenics=0,
                        NonCarcinogenics=0,
                        RespiratoryEffects=0)
    
    # CHG catalyst
    qs.StreamImpactItem(ID='CHG_catalyst_item',
                        linked_stream=stream.CHG_catalyst_out,
                        Acidification=991.6544196,
                        Ecotoxicity=15371.08292,
                        Eutrophication=0.45019348,
                        GlobalWarming=484.7862509,
                        OzoneDepletion=2.23437E-05,
                        PhotochemicalOxidation=6.735405072,
                        Carcinogenics=1.616793132,
                        NonCarcinogenics=27306.37232,
                        RespiratoryEffects=3.517184526)
    
    # membrane distillation and acid extraction
    qs.StreamImpactItem(ID='H2SO4_item',
                        linked_stream=stream.H2SO4,
                        Acidification=0.019678922,
                        Ecotoxicity=0.069909345,
                        Eutrophication=4.05E-06,
                        GlobalWarming=0.008205666,
                        OzoneDepletion=8.94E-10,
                        PhotochemicalOxidation=5.04E-05,
                        Carcinogenics=1.74E-03,
                        NonCarcinogenics=1.68237815,
                        RespiratoryEffects=9.41E-05)
    
    # membrane distillation
    qs.StreamImpactItem(ID='NaOH_item',
                        linked_stream=stream.NaOH,
                        Acidification=0.33656,
                        Ecotoxicity=0.77272,
                        Eutrophication=0.00032908,
                        GlobalWarming=1.2514,
                        OzoneDepletion=7.89E-07,
                        PhotochemicalOxidation=0.0033971,
                        Carcinogenics=0.0070044,
                        NonCarcinogenics=13.228,
                        RespiratoryEffects=0.0024543)
    
    qs.StreamImpactItem(ID='RO_item',
                        linked_stream=stream.Membrane_in,
                        Acidification=0.53533,
                        Ecotoxicity=0.90848,
                        Eutrophication=0.0028322,
                        GlobalWarming=2.2663,
                        OzoneDepletion=0.00000025541,
                        PhotochemicalOxidation=0.0089068,
                        Carcinogenics=0.034791,
                        NonCarcinogenics=31.8,
                        RespiratoryEffects=0.0028778)
    
    # natural gas
    qs.StreamImpactItem(ID='natural_gas_item',
                        linked_stream=stream.natural_gas,
                        Acidification=0.086654753,
                        Ecotoxicity=0.104609119,
                        Eutrophication=7.60E-05,
                        GlobalWarming=1.584234288,
                        OzoneDepletion=1.36924E-07,
                        PhotochemicalOxidation=0.001115155,
                        Carcinogenics=0.00144524,
                        NonCarcinogenics=3.506685032,
                        RespiratoryEffects=0.000326529)
    
    # ammonium sulfate
    qs.StreamImpactItem(ID='NH42SO4_item',
                        linked_stream=stream.ammonium_sulfate,
                        Acidification=-0.72917,
                        Ecotoxicity=-3.4746,
                        Eutrophication=-0.0024633,
                        GlobalWarming=-1.2499,
                        OzoneDepletion=-6.12E-08,
                        PhotochemicalOxidation=-0.0044519,
                        Carcinogenics=-0.036742,
                        NonCarcinogenics=-62.932,
                        RespiratoryEffects=-0.0031315)
    
    # crude oil/petroleum (transportation is included in crude oil item, we offset that first, but we need to add our own transportation)
    qs.StreamImpactItem(ID='biocrude_item',
                        linked_stream=stream.biocrude,
                        Acidification=89/1000/biocrude_density/(0.12917/1000)*0.12698/1000*distance-0.1617,
                        Ecotoxicity=89/1000/biocrude_density/(0.12917/1000)*0.25445/1000*distance-0.10666,
                        Eutrophication=89/1000/biocrude_density/(0.12917/1000)*0.00024901/1000*distance-0.00096886,
                        GlobalWarming=89/1000/biocrude_density*distance-0.22304,
                        OzoneDepletion=89/1000/biocrude_density/(0.12917/1000)*0.000000016986/1000*distance-0.00000060605,
                        PhotochemicalOxidation=89/1000/biocrude_density/(0.12917/1000)*0.001655/1000*distance-0.0013914,
                        Carcinogenics=89/1000/biocrude_density/(0.12917/1000)*0.00046431/1000*distance-0.00030447,
                        NonCarcinogenics=89/1000/biocrude_density/(0.12917/1000)*1.9859/1000*distance-1.0441,
                        RespiratoryEffects=89/1000/biocrude_density/(0.12917/1000)*0.00022076/1000*distance-0.00068606)
    # 89 g CO2/m3/km: carbon intensity of truck transportation (Pootakham et al. A comparison of pipeline versus truck transport of bio-oil. Bioresource Technology, 2010)
    # others are scaled based on transportation data from EcoInvent and GWP from the paper above
    
    sys.simulate() # simulate first to enable the calculation of the income tax rate
    
    federal_income_tax_rate_value = 0.21
    
    annual_sales = BiocrudeTank.outs[0].F_mass*365*24*0.3847 + MemDis.outs[0].F_mass*365*24*0.3236
    annual_material_cost = sum([s.cost for s in sys.feeds[1:] if ((s.price > 0) & (s.ID != 'sludge_assumed_in_wastewater'))]) * sys.operating_hours
    annual_utility_cost = sum([u.utility_cost for u in sys.cost_units]) * sys.operating_hours
    annual_net_income = annual_sales - annual_material_cost - annual_utility_cost
    
    if state == 'average':
        state_income_tax_rate_value = 0.065 # use the mode state income tax from Steward et al. ES&T 2023
    else:
        state_income_tax_rate_value = state_income_tax_rate(state=state, sales=annual_sales, net_income=annual_net_income)
    
    income_tax_rate = federal_income_tax_rate_value + state_income_tax_rate_value
    
    create_tea(sys, IRR_value=0.03, income_tax_value=income_tax_rate, finance_interest_value=0.03, labor_cost_value=(0.41+0.59*size*ww_2_dry_sludge_ratio/100)*10**6)
    
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
    
    qs.LCA(system=sys, lifetime=30, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30/0.67848*elec_GHG,
           # 0.67848 is the GHG level with the Electricity item from ecoinvent,
           # we cannot list electricity GHG one state by one state,
           # but we can adjust the electricity amount to reflect different GHG of electricity at different states
           Cooling=lambda:sys.get_cooling_duty()/1000*30)
    
    biocrude_barrel = BiocrudeTank.outs[0].F_mass/biocrude_density*1000/_oil_barrel_to_L*24 # in BPD (barrel per day)
    
    return sys, biocrude_barrel