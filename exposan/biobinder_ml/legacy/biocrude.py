# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 11:53:17 2025

@author: aliah
"""

# !!! Temporarily ignoring warnings
# import warnings
# warnings.filterwarnings('ignore')

import os, biosteam as bst, qsdsan as qs
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import create_tea
from exposan.biobinder import (
    _HHV_per_GGE,
    _load_components,
    _load_process_settings,
    _units as u,
    BiobinderTEA,
    central_dry_flowrate as default_central,
    data_path,
    feedstock_composition,
    HTL_yields,
    pilot_dry_flowrate as default_pilot,
    price_dct,
    results_path,
    tea_kwargs,
    uptime_ratio,
    )

_psi_to_Pa = 6894.76


# %%

__all__ = ('create_system',)

def create_system(
        flowsheet=None,
        central_dry_flowrate=None,
        pilot_dry_flowrate=None,
        decentralized_HTL=False,
        decentralized_upgrading=False,
        skip_EC=False,
        generate_H2=False,
        EC_config=None,
        crude_fracs=None,
        oil_fracs=None,
        ):
    qs.main_flowsheet.clear()
    print(f"Received flowsheet: {flowsheet}")
    
    if central_dry_flowrate is None:
        central_dry_flowrate = default_central
    if pilot_dry_flowrate is None:
        pilot_dry_flowrate = default_pilot
    
    if decentralized_HTL is False:
        if decentralized_upgrading is False:
            flowsheet_ID = 'bb_CHCU'
            N_HTL = N_upgrading = 1
        else:
            raise ValueError('Centralized HTL, decentralized upgrading is not a valid configuration.')
    else:
        if decentralized_upgrading is False:
            flowsheet_ID = 'bb_DHCU'
            N_HTL = round(central_dry_flowrate/pilot_dry_flowrate)
            N_upgrading = 1
            pilot_dry_flowrate = central_dry_flowrate/N_HTL
        else:
            flowsheet_ID = 'bb_DHDU'
            N_HTL = N_upgrading = 1
            central_dry_flowrate = pilot_dry_flowrate
        
    if skip_EC is True:
        if generate_H2 is True:
            raise ValueError('Cannot generate H2 without EC.')
        if EC_config is not None:
            raise ValueError('Cannot set EC configurations while `skil_EC` is True.')
    
    if flowsheet is None:
        print(f"Creating new flowsheet with ID: {flowsheet_ID}")
        flowsheet = qs.Flowsheet(flowsheet_ID)
        qs.main_flowsheet.set_flowsheet(flowsheet)
    else:
        print(f"Using provided flowsheet with ID: {flowsheet.ID}")
    print(f"Active flowsheet set to: {qs.main_flowsheet.ID}")
    
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()  
    
    _load_process_settings()
    _load_components()
    
    scaled_feedstock = qs.WasteStream('scaled_feedstock', price=price_dct['tipping'])
    FeedstockScaler = u.Scaler(
        'FeedstockScaler', ins=scaled_feedstock, outs='feedstock',
        scaling_factor=N_HTL, reverse=True,
        )
       
    FeedstockTrans = u.Transportation(
        'FeedstockTrans',
        ins=(FeedstockScaler-0, 'feedstock_trans_surrogate'),
        outs=('transported_feedstock',),
        N_unit=N_HTL,
        copy_ins_from_outs=True,
        transportation_unit_cost=0, # will be adjusted later
        transportation_distance=1, # 25 km ref [1]
        )

    # Price accounted for in PWC
    scaled_process_water = qs.WasteStream('scaled_process_water')
    ProcessWaterScaler = u.Scaler(
        'ProcessWaterScaler', ins=scaled_process_water, outs='htl_process_water',
        scaling_factor=N_HTL, reverse=True,
        )

    FeedstockCond = u.Conditioning(
        'FeedstockCond', ins=(FeedstockTrans-0, ProcessWaterScaler.outs[0]),
        outs='conditioned_feedstock',
        feedstock_composition=feedstock_composition,
        feedstock_dry_flowrate=central_dry_flowrate if decentralized_HTL is False else pilot_dry_flowrate,
        target_HTL_solid_loading=0.2,
        )
    FeedstockCond.N_unit = N_HTL # init doesn't take this property

    aqueous_composition = {
        'N': 0.48/100,
        }
    aqueous_composition['HTLaqueous'] = 1 - sum(aqueous_composition.values())
    HTL_kwargs = dict(
        ID='HTL',
        ins=FeedstockCond.outs[0],
        outs=('gas','HTL_aqueous','biocrude','hydrochar'),
        T=280+273.15,
        P=12.4e6, # may lead to HXN error when HXN is included
        # P=101325, # setting P to ambient pressure not practical, but it has minimum effects on the results (several cents)
        tau=15/60,
        dw_yields=HTL_yields,
        gas_composition={'CO2': 1,},
        aqueous_composition=aqueous_composition,
        biocrude_composition={'HTLbiocrude': 1,},
        char_composition={'HTLchar': 1},
        internal_heat_exchanging=True,
        eff_T=60+273.15, # 140.7°F
        eff_P=30*_psi_to_Pa,
        )
    if decentralized_HTL is False:
        HTL_unit = u.CentralizedHTL
    else:
        HTL_unit = u.PilotHTL
        HTL_kwargs['N_unit'] = N_HTL
    HTL = HTL_unit(**HTL_kwargs)

    HTLgasScaler = u.Scaler(
        'HTLgasScaler', ins=HTL-0, outs='scaled_HTLgas',
        scaling_factor=N_HTL, reverse=False,
        )
    
    HTLcharScaler = u.Scaler(
        'HTLcharScaler', ins=HTL.outs[-1], outs='scaled_HTLchar',
        scaling_factor=N_HTL, reverse=False,
        )
    default_EO_voltage = 5
    default_ED_voltage = 30
    default_electrode_cost = 40000
    default_anion_exchange_membrane_cost = 170
    default_cation_exchange_membrane_cost = 190
    
   
    # Only need separate ECs for decentralized HTL, centralized upgrading
    if N_HTL != N_upgrading:
        HTL_EC = u.Electrochemical(
            'HTL_EC',
            ins=(HTL-1, 'HTL_EC_replacement_surrogate'),
            outs=('HTL_EC_gas', 'HTL_EC_H2', 'HTL_EC_N', 'HTL_EC_P', 'HTL_EC_K', 'HTLww_to_disposal'),
            include_PSA=generate_H2,
            EO_voltage=EC_config.get('EO_voltage', default_EO_voltage) if EC_config else default_EO_voltage,
            ED_voltage=EC_config.get('ED_voltage', default_ED_voltage) if EC_config else default_ED_voltage,
            electrode_cost=EC_config.get('electrode_cost', default_electrode_cost) if EC_config else default_electrode_cost,
            anion_exchange_membrane_cost=EC_config.get('anion_exchange_membrane_cost', default_anion_exchange_membrane_cost) if EC_config else default_anion_exchange_membrane_cost,
            cation_exchange_membrane_cost=EC_config.get('cation_exchange_membrane_cost', default_cation_exchange_membrane_cost) if EC_config else default_cation_exchange_membrane_cost,
        )
        HTL_EC.N_unit = N_HTL
        HTL_EC.skip = skip_EC
        
        HTL_ECgasScaler = u.Scaler(
            'HTL_ECgasScaler', ins=HTL_EC.outs[0], outs='scaled_HTLgas',
            scaling_factor=N_HTL, reverse=False,
            )
        
        HTL_ECH2Scaler = u.Scaler(
            'HTL_ECH2Scaler', ins=HTL_EC.outs[1], outs='scaled_HTLH2',
            scaling_factor=N_HTL, reverse=False,
            )
        
        HTL_ECNScaler = u.Scaler(
            'HTL_ECNScaler', ins=HTL_EC.outs[2], outs='scaled_HTLN',
            scaling_factor=N_HTL, reverse=False,
            )
        
        HTL_ECPScaler = u.Scaler(
            'HTL_ECPScaler', ins=HTL_EC.outs[3], outs='scaled_HTLP',
            scaling_factor=N_HTL, reverse=False,
            )
        
        HTL_ECKScaler = u.Scaler(
            'HTL_ECKScaler', ins=HTL_EC.outs[4], outs='scaled_HTLK',
            scaling_factor=N_HTL, reverse=False,
            )

        HTL_ECwwScaler = u.Scaler(
            'HTL_ECscaler', ins=HTL_EC.outs[5], outs='scaled_HTLww_to_disposal',
            scaling_factor=N_HTL, reverse=False,
            )
        
        streams_to_upgrading_EC_lst = []
        streams_to_CHP_lst = []
        H2_streams = [HTL_ECH2Scaler-0]
        N_streams = [HTL_ECNScaler-0]
        P_streams = [HTL_ECPScaler-0]
        K_streams = [HTL_ECKScaler-0]
        liquids_to_disposal_lst = [HTL_ECwwScaler-0]
        solids_to_disposal_lst = [HTLcharScaler-0]
    else:
        streams_to_upgrading_EC_lst = [HTL-1]
        streams_to_CHP_lst = [HTLgasScaler-0, HTLcharScaler-0]
        H2_streams = []
        N_streams = []
        P_streams = []
        K_streams = []
        liquids_to_disposal_lst = []
        solids_to_disposal_lst = []

    BiocrudeTrans = u.Transportation(
        'BiocrudeTrans',
        ins=(HTL-2, 'biocrude_trans_surrogate'),
        outs=('transported_biocrude',),
        N_unit=N_HTL,
        transportation_unit_cost=0, # will be adjusted later
        transportation_distance=1,
        )
    # BiocrudeScaler = u.Scaler(
    #     'BiocrudeScaler', ins=BiocrudeTrans-0, outs='scaled_biocrude',
    #     scaling_factor=N_HTL, reverse=False,
    #     )
    # Biocrude Price $/kg 
    price_dct['biocrude'] = 1
    biocrude = qs.WasteStream('biocrude', price=price_dct['biocrude'])
    
    BiocrudeScaler = u.Scaler(
        'BiocrudeScaler', ins=BiocrudeTrans-0, outs='scaled_biocrude',
        scaling_factor=N_HTL, reverse=False,
        )
    
    BiocrudeStorage = qsu.StorageTank(
        'BiocrudeStorage', BiocrudeScaler-0, outs=biocrude,
        tau=24*3, vessel_material='Stainless steel',
        include_construction=False,
        )
        
    def adjust_prices():
        FeedstockTrans._run()        
        # Centralized HTL and upgrading, transport feedstock
        if decentralized_HTL is False:
            dw_price = price_dct['trans_feedstock'] # $/dry mass
            factor = 1 - FeedstockTrans.ins[0].imass['Water']/FeedstockTrans.ins[0].F_mass
            FeedstockTrans.transportation_unit_cost = dw_price * factor
            BiocrudeTrans.transportation_unit_cost = 0
        # Decentralized HTL, centralized upgrading, transport biocrude
        elif decentralized_upgrading is False:
            FeedstockTrans.transportation_unit_cost = 0
            GGE_price = price_dct['trans_biocrude'] # $/GGE
            # 1e3 to convert from kJ/hr to MJ/hr, 264.172 is m3/hr to gal/hr
            factor = BiocrudeTrans.ins[0].HHV/1e3/(BiocrudeTrans.ins[0].F_vol*264.172)/_HHV_per_GGE
            BiocrudeTrans.transportation_unit_cost = GGE_price * factor
        # Decentralized HTL and upgrading, no transportation needed
        else:
            FeedstockTrans.transportation_unit_cost = BiocrudeTrans.transportation_unit_cost = 0
        
            
    BiocrudeStorage.add_specification(adjust_prices)
    #Other products
    # recovered_H2 = qs.WasteStream('recovered_H2', price=price_dct['H2'])
    # H2mixer = qsu.Mixer('H2mixer', ins=H2_streams, outs=recovered_H2)

    # recovered_N = qs.WasteStream('recovered_N', price=price_dct['N'])
    # Nmixer = qsu.Mixer('Nmixer', ins=N_streams, outs=recovered_N)

    # recovered_P = qs.WasteStream('recovered_P', price=price_dct['P'])
    # Pmixer = qsu.Mixer('Pmixer', ins=P_streams, outs=recovered_P)
    
    # recovered_K = qs.WasteStream('recovered_K', price=price_dct['K'])
    # Kmixer = qsu.Mixer('Kmixer', ins=K_streams, outs=recovered_K)

    natural_gas = qs.WasteStream('natural_gas', CH4=1, price=price_dct['natural_gas'])
    
    CHPmixer = qsu.Mixer('CHPmixer', ins=streams_to_CHP_lst,)
    CHP = qsu.CombinedHeatPower('CHP', 
                                ins=(CHPmixer-0, natural_gas, 'air'),
                                outs=('gas_emissions', 'CHP_solids_to_disposal'),
                                init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.outs[1].price = 0
    solids_to_disposal_lst.append(CHP.outs[-1])
    
    solids_to_disposal = qs.WasteStream('solids_to_disposal', price=price_dct['solids'])
    SolidsDisposalMixer = qsu.Mixer('SolidsDisposalMixer',
                                    ins=solids_to_disposal_lst,
                                    outs=solids_to_disposal)
            
    PWC = u.ProcessWaterCenter(
            'ProcessWaterCenter',
            process_water_streams=scaled_process_water,
            process_water_price=price_dct['process_water']
            )
    PWC.register_alias('PWC')
    @PWC.add_specification
    def run_scalers():
        FeedstockScaler._run()
        ProcessWaterScaler._run()
        PWC._run()
        
    
    ww_to_disposal = qs.WasteStream('ww_to_disposal')
    WWStorage = qsu.StorageTank(
        'WWStorage',liquids_to_disposal_lst, outs=ww_to_disposal, #use liquids_to_disposal_lst as ins for DHCU
        tau=24*1, vessel_material='Stainless steel',
        include_construction=False,
        )
    # @WWStorage.add_specification
    def ww_price():
        WWStorage._run()
        COD_mass_content = HTL.outs[1].COD*HTL.outs[1].F_vol/1e3 # mg/L*m3/hr to kg/hr
        factor = COD_mass_content/HTL.outs[1].F_mass
        ww_to_disposal.price = price_dct['COD']*factor
    
    WWStorage.add_specification(ww_price)
    # WWStorage.run_after_specifications = True
    
    sys = qs.System.from_units(
            'sys',
        units=list(flowsheet.unit),
        operating_hours=365 * 24 * uptime_ratio,
          )

    for unit in sys.units: unit.include_construction = False
        
    gwp_dict = {
          'feedstock': 0,
          'landfill': 400/1e3, # nearly 400 kg CO2e/tonne, Nordahl et al., 2020
          'composting': -41/1e3, # -41 kg CO2e/tonne, Nordahl et al., 2020
          'anaerobic_digestion': (-36-2)/2, # -36 to -2 kg CO2e/tonne, Nordahl et al., 2020
          'trans_feedstock': 0.011856, # 78 km, https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/9393/impact_assessment, Snowden-Swan PNNL 32731
          'trans_biocrude': 0.024472, # 100 miles,https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/9393/impact_assessment, Snowden-Swan PNNL 32731
          'H2': -9.3368, # Compressed gaseous H2 produced from NA NG for SAF Prodcution, GREET 2024
          'natural_gas': 0.3674+(44/16), # NA NG from Shale and Conventional Recovery, include combustion, GREET 2024
          'process_water': 0,
          'electricity': 0.3855 , # kg CO2e/kWh, Non-distributed U.S. Mix, GREET 2024
          'steam': 85.4330/1e3, # g CO2e/MJ, Mix: Natural Gas and Still Gas, GREET 2024
          'cooling_water': 0.068359242, # https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/14408/impact_assessment
          'chilled_water': 0.068359242, # https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/14408/impact_assessment
          'cooling': 0.066033, # kg CO2e/MJ, Feng et al., 2024
          'diesel': -0.6456 , # kg CO2e/kg, Conventional Diesel from Crude Oil for US Refineries, GREET 2024
          'N': -3.5022, # kg CO2e/kg N, Mix: Nitrogen Average, GREET 2024
          'P': -1.1506*(98/31), # H3PO4, Production of Phosphoric Acid, GREET 2024
          'K': -1.6145*(56/39), # KOH, Production, GREET 2024
          'COD': 1.7, # Li et al., 2023
          'wastewater': 0.2629/1e3, # kg CO2e/m3, Industrial Wastewater Treatment, GREET 2024
          'ethylene': 1.0159 #US Average Production, GREET 2024
          }
    
    price_dct['electricity'] = 0.074 #2020$, https://www.energy.gov/sites/default/files/2024-12/hydrogen-shot-water-electrolysis-technology-assessment.pdf
    price_dct['natural_gas'] = 0.14  # $/kg, 2020$, https://www.eia.gov/dnav/ng/hist/n3035us3m.htm 
    qs.PowerUtility.price = EC_config.get('electricity_price', price_dct['electricity']) if EC_config else price_dct['electricity']
    electricity_GHG = EC_config.get('electricity_GHG', gwp_dict['electricity']) if EC_config else gwp_dict['electricity']

    tea = create_tea(sys, cls=BiobinderTEA, **tea_kwargs)
    land_factor = 115000/24 #115000 gal/day, PNNL 32731
    tea.land = lambda: 90000+(tea.system.flowsheet.unit.BiocrudeTrans.ins[0].F_vol*264.172)/land_factor*4.692*1e6 # 6% of 78.2M$, PNNL 32731

# Load impact indicators and items
    # clear_lca_registries()
    # qs.ImpactIndicator.load_from_file(os.path.join(data_path, 'impact_indicators.csv'))
    # qs.ImpactItem.load_from_file(os.path.join(data_path, 'impact_items.xlsx'))
    # print("Loaded Impact Items:")
   
    GWP = qs.ImpactIndicator('GWP',
                              alias='GlobalWarmingPotential',
                              method='GREET',
                              category='environmental impact',
                              unit='kg CO2-eq',)
           
    feedstock_item = qs.StreamImpactItem(
        ID='feedstock_item',
        linked_stream=scaled_feedstock,
        GWP=-gwp_dict['landfill'],
        )
    trans_feedstock_item = qs.StreamImpactItem(
        ID='feedstock_trans_surrogate_item',
        linked_stream=FeedstockTrans.ins[1],
        GWP=gwp_dict['trans_feedstock'],
        )
    process_water_item = qs.StreamImpactItem(
        ID='scaled_process_water_item',
        linked_stream=ProcessWaterScaler.ins[0],
        GWP=gwp_dict['process_water'],
        )
    trans_biocrude_item=qs.StreamImpactItem(
        ID='biocrude_trans_surrogate_item',
        linked_stream=BiocrudeTrans.ins[1],
        GWP=gwp_dict['trans_biocrude'],
        )
    natural_gas_item = qs.StreamImpactItem(
        ID='natural_gas_item',
        linked_stream=natural_gas,
        GWP=gwp_dict['natural_gas'],
        )
    ww_to_disposal_item = qs.StreamImpactItem(
        ID='ww_to_disposal_item',
        linked_stream=ww_to_disposal,
        GWP=gwp_dict['wastewater'], 
        )
    solids_to_disposal_item = qs.StreamImpactItem(
        ID='solids_to_disposal_item',
        linked_stream=solids_to_disposal,
        GWP=gwp_dict['trans_feedstock'],
        )
    e_item = qs.ImpactItem(
        ID='e_item',
        GWP=electricity_GHG,
        )
    steam_item = qs.ImpactItem(
        ID='steam_item',
        GWP=gwp_dict['steam'],
        )
    cooling_item = qs.ImpactItem(
        ID='cooling_item',
        GWP=gwp_dict['cooling'],
        )
    # LCA adjustment based on system configuration
    fake_stream= qs.SanStream('Nothing', price=0)
    nothing_item = qs.StreamImpactItem(
        ID='nothing_item',
        linked_stream=fake_stream,
        GWP=0,
        )
    if decentralized_HTL is False:
    # Centralized HTL, Centralized upgrading
        trans_feedstock_item.linked_stream = FeedstockTrans.ins[1]  # feedstock transportation stream
        trans_biocrude_item.linked_stream = fake_stream             # No biocrude transportation
    elif decentralized_upgrading is False:
    # Decentralized HTL, centralized upgrading
        trans_feedstock_item.linked_stream = fake_stream            # No feedstock transportation
        trans_biocrude_item.linked_stream = BiocrudeTrans.ins[1]    #biocrude tranportation stream
    else:
    # Fully decentralized (no transportation needed)
        trans_feedstock_item.linked_stream = fake_stream            # No feedstock transportation
        trans_biocrude_item.linked_stream = fake_stream             # 
    
    tea = create_tea(sys, cls=BiobinderTEA, **tea_kwargs)
    lifetime = tea_kwargs['duration'][1] - tea_kwargs['duration'][0]
    lca = qs.LCA(
            system=sys,
            lifetime=lifetime,
            simulate_system=False,
            uptime_ratio=sys.operating_hours / (365 * 24),
            e_item=lambda: (sys.get_electricity_consumption() - sys.get_electricity_production()) * lifetime,
            steam_item=lambda: sys.get_heating_duty() / 1000 * lifetime,
            cooling_item=lambda: sys.get_cooling_duty() / 1000 * lifetime, # this function will run during LCA
            )
    return sys

def simulate_and_print(sys, save_report=False):
    sys.simulate()
    tea = sys.TEA
    lca = sys.LCA
    biocrude = sys.flowsheet.stream.biocrude

    # Solve for the minimum selling price (MSP) of biocrude
    biocrude.price = MSP = tea.solve_price(biocrude)
    print(f'Minimum selling price of the biocrude is ${MSP:.2f}/kg.')
    all_impacts = lca.get_allocated_impacts(streams=(biocrude,), operation_only=True, annual=True)
    GWP = all_impacts['GWP']/(biocrude.F_mass*lca.system.operating_hours)
    print(f'Global warming potential of the biocrude is {GWP:.4f} kg CO2e/kg.')
    biocrude.price = 0.50   #$/kg
    IRR= tea.solve_IRR()
    print(f'Internal rate of return of the process is {IRR * 100:.2f}%')
      
if __name__ == '__main__':    
  
    config_kwargs = dict(
      flowsheet=None, 
      central_dry_flowrate=None,
      pilot_dry_flowrate=None,
      EC_config=None
      )
    # What to do with HTL-AP
    config_kwargs.update(dict(skip_EC=True, generate_H2=False, EC_config=None)) # no EC
    # config_kwargs.update(dict(skip_EC=False, generate_H2=False, EC_config=None)) # EC, recover nutrients only
    # config_kwargs.update(dict(skip_EC=False, generate_H2=True, EC_config=None)) # EC, recover nutrients and generate H2
    # config_kwargs.update(dict(skip_EC=False, generate_H2=True, EC_config=EC_future_config)) # EC, recovery nutrients, generate H2, optimistic assumptions
    
    # Decentralized vs. centralized configuration
    # config_kwargs.update(dict(decentralized_HTL=False, decentralized_upgrading=False)) # CHCU
    config_kwargs.update(dict(decentralized_HTL=True, decentralized_upgrading=False)) # DHCU
    
    # Distillation column cost calculation doesn't scale down well, so the cost is very high now.
    # But maybe don't need to do the DHDU scenario, if DHCU isn't too different from CHCU
    # However, maybe the elimination of transportation completely will make a difference
    # config_kwargs.update(dict(decentralized_HTL=True, decentralized_upgrading=True)) # DHDU
    
    sys = create_system(**config_kwargs)
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    # tea = sys.TEA
    # lca = sys.LCA
    
    simulate_and_print(sys)