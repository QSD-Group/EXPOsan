#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 08:46:28 2023

@author: jiananfeng
"""

import os, qsdsan as qs, biosteam as bst, pandas as pd
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import (
    _load_components,
    create_tea,
    )
from exposan.htl import _sanunits as su

__all__ = ('create_spatial_system',)

def _load_process_settings(location='IL'):
# =============================================================================
#     add a heating agent
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
#     set utility prices
# =============================================================================

    folder = '/Users/jiananfeng/Desktop/PhD_CEE/NSF_PFAS/HTL_geospatial/'
    
    elec_input = pd.read_excel(folder + 'HTL_spatial_model_input.xlsx')

    bst.PowerUtility.price = elec_input[elec_input['state']==location]['price (10-year median)'].iloc[0]/100
    
    # # These utilities are provided by CHP thus cost already considered
    # # setting the regeneration price to 0 or not will not affect the final results
    # # as the utility cost will be positive for the unit that consumes it
    # # but negative for HXN/CHP as they produce it
    # for adj in ('low', 'medium', 'high'):
    #     steam = bst.HeatUtility.get_agent(f'{adj}_pressure_steam')
    #     steam.heat_transfer_price = steam.regeneration_price = 0.

def create_spatial_system(waste_price=400, waste_GHG=800, size=50, distance=30, solid_fate=1, ww_2_dry_sludge_ratio=1, solid_reduction=0, state='IL', elec_GHG=0.365593393303875):
# 400: Seiple et al. 2020, use 400 for now,
# later, we may want to give different prices according to the original disposal methods.
# size in MGD
# distance in km

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
    
    # Construction here, StreamImpactItem after TEA
    folder = os.path.dirname(__file__)
    qs.ImpactIndicator.load_from_file(os.path.join(folder, 'data/impact_indicators.csv'))
    qs.ImpactItem.load_from_file(os.path.join(folder, 'data/impact_items.xlsx'))
    
    raw_wastewater = qs.WasteStream('raw_wastewater', H2O=size, units='MGD', T=25+273.15)
    # Jones baseline: 1276.6 MGD, 1.066e-4 $/kg ww
    # set H2O equal to the total raw wastewater into the WWTP
    
    # =============================================================================
    # pretreatment (Area 000)
    # =============================================================================
    
    if solid_fate == 8:
        WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                       ww_2_dry_sludge=ww_2_dry_sludge_ratio*(1-solid_reduction), # use real sludge data
                       # how much metric ton/day sludge can be produced by 1 MGD of ww
                       sludge_moisture=0.99, sludge_dw_ash=0.257, 
                       sludge_afdw_lipid=0.204, sludge_afdw_protein=0.463, operation_hours=7920)
        
        SluC = qsu.SludgeCentrifuge('A000', ins=WWTP-0,
                                outs=('supernatant','compressed_sludge'),
                                init_with='Stream',
                                solids=('Sludge_lipid','Sludge_protein',
                                        'Sludge_carbo','Sludge_ash'),
                                sludge_moisture=0.8)
        SluC.register_alias('SluC')
        
        P1 = qsu.SludgePump('A100', ins=SluC-1, outs='pressed_sludge', P=3049.7*6894.76,
                  init_with='Stream')
    
    elif solid_fate in (1,2,4):
        WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                       ww_2_dry_sludge=ww_2_dry_sludge_ratio*(1-solid_reduction), # use real sludge data
                       # how much metric ton/day sludge can be produced by 1 MGD of ww
                       sludge_moisture=0.8, sludge_dw_ash=0.454166666666667, 
                       sludge_afdw_lipid=0.1693, sludge_afdw_protein=0.5185, operation_hours=7920)
        
        P1 = qsu.SludgePump('A100', ins=WWTP-0, outs='pressed_sludge', P=3049.7*6894.76,
                  init_with='Stream')
        
    else:
        WWTP = su.WWTP('S000', ins=raw_wastewater, outs=('sludge','treated_water'),
                       ww_2_dry_sludge=ww_2_dry_sludge_ratio*(1-solid_reduction), # use real sludge data
                       # how much metric ton/day sludge can be produced by 1 MGD of ww
                       sludge_moisture=0.8, sludge_dw_ash=0.257, 
                       sludge_afdw_lipid=0.204, sludge_afdw_protein=0.463, operation_hours=7920)
        
        P1 = qsu.SludgePump('A100', ins=WWTP-0, outs='pressed_sludge', P=3049.7*6894.76,
                  init_with='Stream')
        
    WWTP.register_alias('WWTP')
    
    raw_wastewater.price = -WWTP.ww_2_dry_sludge*waste_price/3.79/(10**6)

    # =============================================================================
    # HTL (Area 100)
    # =============================================================================
    
    P1.register_alias('P1')
    # Jones 2014: 3049.7 psia
    
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
        HTL = qsu.HydrothermalLiquefaction('A120', ins=H1-0, outs=('biochar','HTL_aqueous','biocrude','offgas_HTL'), dewatered_unit_exist_in_the_system=True)
    else:
        HTL = qsu.HydrothermalLiquefaction('A120', ins=H1-0, outs=('biochar','HTL_aqueous','biocrude','offgas_HTL'), dewatered_unit_exist_in_the_system=False)
    HTL.register_alias('HTL')
    
    # =============================================================================
    # Storage, and disposal (Area 500)
    # =============================================================================
    
    CrudeOilTank = qsu.StorageTank('T500', ins=HTL-2, outs=('crude_oil'),
                                    tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')
    # store for 3 days based on Jones 2014
    CrudeOilTank.register_alias('CrudeOilTank')
    
    CrudeOilTank.outs[0].price = -0.000075*distance + 0.3847

    # =============================================================================
    # facilities
    # =============================================================================
    
    qsu.HeatExchangerNetwork('HXN', force_ideal_thermo=True)
    
    sys = qs.System.from_units(
        'sys_geospatial',
        units=list(flowsheet.unit), 
        operating_hours=WWTP.operation_hours, # 7920 hr Jones
        )
    sys.register_alias('sys')

    ##### Add stream impact items #####

    # add impact for waste sludge
    qs.StreamImpactItem(ID='transportation_item',
                        linked_stream=stream.crude_oil,
                        Acidification=0.12698/1000*distance-0.1617,
                        Ecotoxicity=0.25445/1000*distance-0.10666,
                        Eutrophication=0.00024901/1000*distance-0.00096886,
                        GlobalWarming=0.12917/1000*distance-0.22304,
                        OzoneDepletion=0.000000016986/1000*distance-0.00000060605,
                        PhotochemicalOxidation=0.001655/1000*distance-0.0013914,
                        Carcinogenics=0.00046431/1000*distance-0.00030447,
                        NonCarcinogenics=1.9859/1000*distance-1.0441,
                        RespiratoryEffects=0.00022076/1000*distance-0.00068606)
    
    create_tea(sys)
    qs.LCA(
        system=sys, lifetime=30, lifetime_unit='yr',
        Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())/0.67848*elec_GHG*30, # 0.67848 is the GHG level with the Electricity item, we can adjust the electricity amount to reflect different GHG of electricity at different locations
        Cooling=lambda:sys.get_cooling_duty()/1000*30,
        )
    sys.simulate()
    
    biocrude_barrel = CrudeOilTank.outs[0].F_mass/0.98/3.78541/42*7920/365 # in BPD (barrel per day)
    
    return sys, biocrude_barrel