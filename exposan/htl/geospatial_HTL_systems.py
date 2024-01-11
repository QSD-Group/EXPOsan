#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 5 08:46:28 2023

@author: jiananfeng
"""

import os
import qsdsan as qs, biosteam as bst
import pandas as pd
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import (_load_components, create_tea,)
from exposan.htl import _sanunits as su
from biosteam.units import IsenthalpicValve
from biosteam import settings

__all__ = ('create_spatial_system',)

density_biocrude = 980 # kg/m3 Snowden-Swan et al. 2022 SOT, PNNL

def _load_process_settings(location='IL'):
# =============================================================================
# add a heating agent
# =============================================================================
    # use DOWTHERM(TM) A Heat Transfer Fluid (HTF) as the heating agent
    # DOWTHERM(TM) A HTF = 73.5% diphenyl oxide (DPO) + 26.5% Biphenyl (BIP)
    # critical temperature for HTF: 497 C
    # critical pressure for HTF: 313.4 kPa
    # https://www.dow.com/en-us/pdp.dowtherm-a-heat-transfer-fluid.238000z.\
    # html#tech-content (accessed 11-16-2022)
    
    DPO_chem = qs.Chemical('DPO_chem', search_ID='101-84-8')
    BIP_chem = qs.Chemical('BIP_chem', search_ID='92-52-4')
    
    DPO = qs.Component.from_chemical('DPO', chemical=DPO_chem, particle_size='Soluble',
                                     degradability='Slowly', organic=True)
    
    BIP = qs.Component.from_chemical('BIP', chemical=BIP_chem, particle_size='Soluble',
                                     degradability='Slowly', organic=True)
    
    HTF_thermo = bst.Thermo((DPO, BIP,))
    
    HTF = bst.UtilityAgent('HTF', DPO=0.735, BIP=0.265, T=673.15, P=951477, phase='g',
                           # 400 C (673.15 K) and 138 psig (951477 pa) are max temp and pressure for HTF
                           thermo=HTF_thermo,
                           # T_limit = 495 F (530.372 K) is the highest temp that vapor can exist
                           regeneration_price=1) # Lang
                           # use default heat transfer efficiency (1)
    # Temperature and pressure: https://www.dow.com/content/dam/dcc/documents/\
    # en-us/app-tech-guide/176/176-01334-01-dowtherm-heat-transfer-fluids-\
    # engineering-manual.pdf?iframe=true (accessed on 11-16-2022)
    bst.HeatUtility.heating_agents.append(HTF)

    bst.CE = qs.CEPCI_by_year[2020] # use 2020$ to match up with latest PNNL report
    
# =============================================================================
# set utility prices
# =============================================================================

    folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/'
    
    elec_input = pd.read_excel(folder + 'HTL_geospatial_model_input.xlsx')

    bst.PowerUtility.price = elec_input[elec_input['state']==location]['price (10-year median)'].iloc[0]/100

# =============================================================================
# create the HTL spatial system
# =============================================================================

def create_spatial_system(waste_cost=400, # assumed to be 400 for all WRRFs
                          size=100, # in MGD
                          distance=30, # in km, using Google Maps API
                          solid_fate=1, # from Seiple et al. 2020
                          ww_2_dry_sludge_ratio=1, # use real solid amount/WRRFs influent
                          state='IL',
                          elec_GHG=0.365593393303875):

    flowsheet_ID = 'htl_geospatial'
    
    # Clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    _load_process_settings(location=state)
    
    # folder = os.path.dirname(__file__)
    folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/'
    qs.ImpactIndicator.load_from_file(os.path.join(folder, 'data/impact_indicators.csv'))
    qs.ImpactItem.load_from_file(os.path.join(folder, 'data/impact_items.xlsx'))
    
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)
    # set H2O equal to the total raw wastewater into the WWTP
    
    # =============================================================================
    # pretreatment (Area 000)
    # =============================================================================
    
    if solid_fate == 8: # lagoon and constructed wetland
        WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                       ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                       # how much metric ton/day sludge can be produced by 1 MGD of ww
                       sludge_moisture=0.99, # dewatering needed
                       sludge_dw_ash=0.257, # w/o digestion
                       sludge_afdw_lipid=0.204, # w/o digestion
                       sludge_afdw_protein=0.463, # w/o digestion
                       operation_hours=8760)
        
        SluC = qsu.SludgeCentrifuge('A000', ins=WWTP-0,
                                outs=('supernatant','compressed_sludge'),
                                init_with='Stream',
                                solids=('Sludge_lipid','Sludge_protein',
                                        'Sludge_carbo','Sludge_ash'),
                                sludge_moisture=0.8)
        SluC.register_alias('SluC')
        
        P1 = qsu.SludgePump('A100', ins=SluC-1, outs='pressed_sludge', P=3049.7*6894.76,
                  init_with='Stream')
    
    elif solid_fate in (1, 2, 4): # anaerobic/aerobic digestion
        WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                       ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                       # how much metric ton/day sludge can be produced by 1 MGD of ww
                       sludge_moisture=0.8, # dewatering not needed
                       sludge_dw_ash=0.454166666666667, # w/ digestion
                       sludge_afdw_lipid=0.1693, # w/ digestion
                       sludge_afdw_protein=0.5185, # w/ digestion
                       operation_hours=8760)
        
        P1 = qsu.SludgePump('A100', ins=WWTP-0, outs='pressed_sludge', P=3049.7*6894.76,
                  init_with='Stream')
        
    else:
        WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                       ww_2_dry_sludge=ww_2_dry_sludge_ratio,
                       # how much metric ton/day sludge can be produced by 1 MGD of ww
                       sludge_moisture=0.8, # dewatering not needed
                       sludge_dw_ash=0.257, # w/o digestion
                       sludge_afdw_lipid=0.204, # w/o digestion
                       sludge_afdw_protein=0.463, # w/o digestion
                       operation_hours=8760)
        
        P1 = qsu.SludgePump('A100', ins=WWTP-0, outs='pressed_sludge', P=3049.7*6894.76,
                  init_with='Stream')
        
    WWTP.register_alias('WWTP')
    
    P1.register_alias('P1')
    # Jones 2014: 3049.7 psia
    
    raw_wastewater.price = -WWTP.ww_2_dry_sludge*waste_cost/3.79/(10**6) # 1 gal water = 3.79 kg water

    # =============================================================================
    # HTL (Area 100)
    # =============================================================================
    
    H1 = qsu.HXutility('A110', include_construction=True,
                       ins=P1-0, outs='heated_sludge', T=351+273.15,
                       U=0.0795, init_with='Stream', rigorous=True)
    # feed T is low, thus high viscosity and low U (case B in Knorr 2013)
    # U: 3, 14, 15 BTU/hr/ft2/F as minimum, baseline, and maximum
    # U: 0.0170348, 0.0794957, 0.085174 kW/m2/K
    # H1: SS PNNL 2020: 50 (17-76) Btu/hr/ft2/F ~ U = 0.284 (0.096-0.4313) kW/m2/K
    # but not in other pumps (low viscosity, don't need U to enforce total heat transfer efficiency)
    # unit conversion: https://www.unitsconverters.com/en/Btu(It)/Hmft2mdegf-To-W/M2mk/Utu-4404-4398
    H1.register_alias('H1')
    
    if solid_fate == 8:
        HTL = qsu.HydrothermalLiquefaction('A120', ins=H1-0,
                                           outs=('hydrochar','HTL_aqueous','biocrude','offgas_HTL'),
                                           mositure_adjustment_exist_in_the_system=True)
    else:
        HTL = qsu.HydrothermalLiquefaction('A120', ins=H1-0,
                                           outs=('hydrochar','HTL_aqueous','biocrude','offgas_HTL'),
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

    AcidEx = su.AcidExtraction('A200', ins=(HTL-0, SP1-0),
                               outs=('residual','extracted'))
    AcidEx.register_alias('AcidEx')
    # AcidEx.outs[0].price = -0.055 # SS 2021 SOT PNNL report page 24 Table 9
    # not include residual for TEA and LCA for now
    
    M1_outs1 = AcidEx.outs[1]
    M1 = su.HTLmixer('A210', ins=(HTL-1, M1_outs1), outs=('mixture',))
    M1.register_alias('M1')
    
    StruPre = su.StruvitePrecipitation('A220', ins=(M1-0,'MgCl2','NH4Cl','MgO'),
                                       outs=('struvite','CHG_feed'))
    StruPre.ins[1].price = 0.5452
    StruPre.ins[2].price = 0.13
    StruPre.ins[3].price = 0.2
    StruPre.outs[0].price = 0.661
    StruPre.register_alias('StruPre')
    
    CHG = qsu.CatalyticHydrothermalGasification('A230',
                                                ins=(StruPre-1, '7.8%_Ru/C'),
                                                outs=('CHG_out', '7.8%_Ru/C_out'))
    CHG.ins[1].price = 134.53
    CHG.register_alias('CHG')
    
    V1 = IsenthalpicValve('A240', ins=CHG-0, outs='depressed_cooled_CHG', P=50*6894.76, vle=True)
    V1.register_alias('V1')
    
    F1 = qsu.Flash('A250', ins=V1-0, outs=('CHG_fuel_gas','N_riched_aqueous'),
                      T=60+273.15, P=50*6894.76, thermo=settings.thermo.ideal())
    F1.register_alias('F1')
    
    MemDis = qsu.MembraneDistillation('A260', ins=(F1-1, SP1-1, 'NaOH', 'Membrane_in'),
                                  outs=('ammonium_sulfate','MemDis_ww', 'Membrane_out','solution'),
                                  init_with='WasteStream')
    MemDis.ins[2].price = 0.5256
    MemDis.outs[0].price = 0.3236
    MemDis.register_alias('MemDis')
    
    # =============================================================================
    # Storage, and disposal (Area 300)
    # =============================================================================
    
    CrudeOilTank = qsu.StorageTank('T300', ins=HTL-2, outs=('crude_oil'),
                                    tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')
    # store for 3 days based on Jones 2014
    CrudeOilTank.register_alias('CrudeOilTank')
    
    CrudeOilTank.outs[0].price = -5.67/density_biocrude - 0.07/density_biocrude*distance + 0.3847
    # 5.67: fixed cost, 0.07: variable cost (Pootakham et al. Bio-oil transport by pipeline: A techno-economic assessment. Bioresource Technology, 2010)
    
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
    
    sys = qs.System.from_units('sys_geospatial',
                               units=list(flowsheet.unit),
                               operating_hours=WWTP.operation_hours)
    sys.register_alias('sys')

# =============================================================================
# add stream impact items
# =============================================================================
   
    # transportation - crude oil   
    qs.StreamImpactItem(ID='transportation_item',
                        linked_stream=stream.crude_oil,
                        Acidification=89/1000/density_biocrude/(0.12917/1000)*0.12698/1000*distance-0.1617,
                        Ecotoxicity=89/1000/density_biocrude/(0.12917/1000)*0.25445/1000*distance-0.10666,
                        Eutrophication=89/1000/density_biocrude/(0.12917/1000)*0.00024901/1000*distance-0.00096886,
                        GlobalWarming=89/1000/density_biocrude*distance-0.22304,
                        OzoneDepletion=89/1000/density_biocrude/(0.12917/1000)*0.000000016986/1000*distance-0.00000060605,
                        PhotochemicalOxidation=89/1000/density_biocrude/(0.12917/1000)*0.001655/1000*distance-0.0013914,
                        Carcinogenics=89/1000/density_biocrude/(0.12917/1000)*0.00046431/1000*distance-0.00030447,
                        NonCarcinogenics=89/1000/density_biocrude/(0.12917/1000)*1.9859/1000*distance-1.0441,
                        RespiratoryEffects=89/1000/density_biocrude/(0.12917/1000)*0.00022076/1000*distance-0.00068606)
    # 89 g CO2/m3/km: carbon intensity of truck transportation (Pootakham et al. A comparison of pipeline versus truck transport of bio-oil. Bioresource Technology, 2010)
    # others are scaled based on transportation data from EcoInvent and GWP from the paper above
    
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
    
    # Membrane distillation and acid extraction
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
    
    # Membrane distillation
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
    
    # Struvite precipitation
    qs.StreamImpactItem(ID='MgCl2_item',
                        linked_stream=stream.MgCl2,
                        Acidification=0.77016,
                        Ecotoxicity=0.97878,
                        Eutrophication=0.00039767,
                        GlobalWarming=2.8779,
                        OzoneDepletion=4.94E-08,
                        PhotochemicalOxidation=0.0072306,
                        Carcinogenics=0.0050938,
                        NonCarcinogenics=8.6916,
                        RespiratoryEffects=0.004385)

    qs.StreamImpactItem(ID='NH4Cl_item',
                        linked_stream=stream.NH4Cl,
                        Acidification=0.34682,
                        Ecotoxicity=0.90305, 
                        Eutrophication=0.0047381,
                        GlobalWarming=1.525,
                        OzoneDepletion=9.22E-08,
                        PhotochemicalOxidation=0.0030017,
                        Carcinogenics=0.010029,
                        NonCarcinogenics=14.85,
                        RespiratoryEffects=0.0018387)
    
    qs.StreamImpactItem(ID='MgO_item',
                        linked_stream=stream.MgO,
                        Acidification=0.12584,
                        Ecotoxicity=2.7949,
                        Eutrophication=0.00063607,
                        GlobalWarming=1.1606,
                        OzoneDepletion=1.54E-08,
                        PhotochemicalOxidation=0.0017137,
                        Carcinogenics=0.018607,
                        NonCarcinogenics=461.54,
                        RespiratoryEffects=0.0008755)
    
    # Heating and power utilities
    qs.StreamImpactItem(ID='natural_gas_item',
                        linked_stream=stream.natural_gas,
                        Acidification=0.083822558,
                        Ecotoxicity=0.063446198,
                        Eutrophication=7.25E-05,
                        GlobalWarming=1.584234288,
                        OzoneDepletion=1.23383E-07,
                        PhotochemicalOxidation=0.000973731,
                        Carcinogenics=0.000666424,
                        NonCarcinogenics=3.63204,
                        RespiratoryEffects=0.000350917)
    
    # Struvite
    qs.StreamImpactItem(ID='struvite_item',
                        linked_stream=stream.struvite,
                        Acidification=-0.122829597,
                        Ecotoxicity=-0.269606396,
                        Eutrophication=-0.000174952,
                        GlobalWarming=-0.420850152,
                        OzoneDepletion=-2.29549E-08,
                        PhotochemicalOxidation=-0.001044087,
                        Carcinogenics=-0.002983018,
                        NonCarcinogenics=-4.496533528,
                        RespiratoryEffects=-0.00061764)
    
    # Ammonium sulfate
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
    
    create_tea(sys, IRR_value=0.03, finance_interest_value=0.03)
    
    qs.LCA(system=sys, lifetime=30, lifetime_unit='yr',
           Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30/0.67848*elec_GHG,
           # 0.67848 is the GHG level with the Electricity item,
           # we can adjust the electricity amount to reflect different GHG of electricity at different locations
           Cooling=lambda:sys.get_cooling_duty()/1000*30)
    
    sys.simulate()
    
    biocrude_barrel = CrudeOilTank.outs[0].F_mass/0.98/3.78541/42*24 # in BPD (barrel per day)
    
    return sys, biocrude_barrel