# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References
[1] Snowden-Swan et al., Wet Waste Hydrothermal Liquefaction and Biocrude Upgrading to Hydrocarbon Fuels:
    2021 State of Technology; PNNL-32731; Pacific Northwest National Lab. (PNNL), Richland, WA (United States), 2022.
    https://doi.org/10.2172/1863608.
'''

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
        
    if skip_EC is True and generate_H2 is True:
        raise ValueError('Cannot generate H2 without EC.')
    
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

    # Only need separate ECs for decentralized HTL, centralized upgrading
    if N_HTL != N_upgrading:
        HTL_EC = u.Electrochemical(
            'HTL_EC',
            ins=(HTL-1, 'HTL_EC_replacement_surrogate'),
            outs=('HTL_EC_gas', 'HTL_EC_H2', 'HTL_EC_N', 'HTL_EC_P', 'HTL_EC_K', 'HTLww_to_disposal'),
            include_PSA=generate_H2,
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
    
    BiocrudeScaler = u.Scaler(
        'BiocrudeScaler', ins=BiocrudeTrans-0, outs='scaled_biocrude',
        scaling_factor=N_HTL, reverse=False,
        )
    
    crude_fracs = [0.0339, 0.8104+0.1557]
    oil_fracs = [0.5316, 0.4684]
    BiocrudeSplitter = u.BiocrudeSplitter(
        'BiocrudeSplitter', ins=BiocrudeScaler-0, outs='splitted_crude',
        biocrude_IDs=('HTLbiocrude'),
        cutoff_fracs=[crude_fracs[0], crude_fracs[1]*oil_fracs[0], crude_fracs[1]*oil_fracs[1]], # light (water): medium/heavy (biocrude/char)
        cutoff_Tbs=(150+273.15, 343+273.15,),
        )
    
    CrudePump = qsu.Pump('CrudePump', init_with='Stream', 
                         ins=BiocrudeSplitter-0, outs='crude_to_dist',)
    
    # Separate water from organics (bp<150°C)
    CrudeLightDis = qsu.ShortcutColumn(
        'CrudeLightDis', ins=CrudePump-0,
        outs=('crude_light','crude_medium_heavy'),
        LHK=BiocrudeSplitter.keys[0],
        P=50*_psi_to_Pa,
        Lr=0.87,
        Hr=0.98,
        k=2, is_divided=True)
    
    CrudeLightFlash = qsu.Flash('CrudeLightFlash', ins=CrudeLightDis-0, T=298.15, P=101325,)
    streams_to_CHP_lst.append(CrudeLightFlash.outs[0])
    streams_to_upgrading_EC_lst.append(CrudeLightFlash.outs[1])

    # Separate fuel from biobinder
    CrudeHeavyDis = qsu.ShortcutColumn(
        'CrudeHeavyDis', ins=CrudeLightDis-1,
        outs=('hot_biofuel','hot_biobinder'),
        LHK=BiocrudeSplitter.keys[1],
        P=50*_psi_to_Pa,
        Lr=0.75,
        Hr=0.89, # 0.85 and 0.89 are more better
        k=2, is_divided=True)
    
    # import numpy as np
    # from exposan.saf.utils import find_Lr_Hr
    # oil_fracs = [0.5316, 0.4684]
    # Lr_range = np.arange(0.5, 1, 0.05)
    # Hr_range = np.arange(0.75, 1, 0.05)
    # results = find_Lr_Hr(CrudeHeavyDis, Lr_trial_range=Lr_range, Hr_trial_range=Hr_range)
    # # results = find_Lr_Hr(CrudeHeavyDis, target_light_frac=oil_fracs[0], Lr_trial_range=Lr_range, Hr_trial_range=Hr_range)
    # results_df, Lr, Hr = results
    
    CrudeHeavyDis_run = CrudeHeavyDis._run
    CrudeHeavyDis_design = CrudeHeavyDis._design
    CrudeHeavyDis_cost = CrudeHeavyDis._cost
    def run_design_cost():
        CrudeHeavyDis_run()
        try:
            CrudeHeavyDis_design()
            CrudeHeavyDis_cost()
            if all([v>0 for v in CrudeHeavyDis.baseline_purchase_costs.values()]):
                # Save for later debugging
                # print('design')
                # print(CrudeHeavyDis.design_results)
                # print('cost')
                # print(CrudeHeavyDis.baseline_purchase_costs)
                # print(CrudeHeavyDis.installed_costs) # this will be empty
                return
        except: pass
        raise RuntimeError('`CrudeHeavyDis` simulation failed.')

    # Simulation may converge at multiple points, filter out unsuitable ones
    def screen_results():
        ratio0 = oil_fracs[0]
        lb, ub = round(ratio0,2)-0.05, round(ratio0,2)+0.05
        try: 
            run_design_cost()
            status = True
        except: 
            status = False
        def get_ratio():
            if CrudeHeavyDis.F_mass_out > 0:
                return CrudeHeavyDis.outs[0].F_mass/CrudeHeavyDis.F_mass_out
            return 0
        n = 0
        ratio = get_ratio()
        while (status is False) or (ratio<lb) or (ratio>ub):
            try: 
                run_design_cost()
                status = True
            except: 
                status = False
            ratio = get_ratio()
            n += 1
            if n >= 100:
                status = False
                raise RuntimeError(f'No suitable solution for `CrudeHeavyDis` within {n} simulation.')
    CrudeHeavyDis._run = screen_results

    def do_nothing(): pass
    CrudeHeavyDis._design = CrudeHeavyDis._cost = do_nothing

    BiofuelFlash = qsu.Flash('BiofuelFlash', ins=CrudeHeavyDis-0, outs=('', 'cooled_biofuel',),
                              T=298.15, P=101325)
    streams_to_upgrading_EC_lst.append(BiofuelFlash.outs[0])
    
    BiobinderHX = qsu.HXutility('BiobinderHX', ins=CrudeHeavyDis-1, outs=('cooled_biobinder',),
                                T=298.15)
    
    UpgradingECmixer = qsu.Mixer('UpgradingECmixer', ins=streams_to_upgrading_EC_lst, outs='ww_to_upgrading_EC',)
    Upgrading_EC = u.Electrochemical(
        'Upgrading_EC',
        ins=(UpgradingECmixer-0, 'Upgrading_EC_replacement_surrogate'),
        outs=('Upgrading_EC_gas', 'Upgrading_EC_H2', 'Upgrading_EC_N', 'Upgrading_EC_P', 'Upgrading_EC_K', 'ECww_to_disposal'),
        include_PSA=generate_H2,
        )
    Upgrading_EC.N_unit = N_upgrading
    Upgrading_EC.skip = skip_EC
    H2_streams.append(Upgrading_EC-1)
    N_streams.append(Upgrading_EC-2)
    P_streams.append(Upgrading_EC-3)
    K_streams.append(Upgrading_EC-4)
    
    liquids_to_disposal_lst.append(Upgrading_EC.outs[-1])
    
    ww_to_disposal = qs.WasteStream('ww_to_disposal')
    WWdisposalMixer = qsu.Mixer('WWdisposalMixer', ins=liquids_to_disposal_lst, outs=ww_to_disposal)
    @WWdisposalMixer.add_specification        
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
        
        # Wastewater
        WWdisposalMixer._run()
        COD_mass_content = ww_to_disposal.COD*ww_to_disposal.F_vol/1e3 # mg/L*m3/hr to kg/hr
        factor = COD_mass_content/ww_to_disposal.F_mass
        ww_to_disposal.price = price_dct['COD']*factor

    # 3-day storage time as in the SAF module
    biofuel = qs.WasteStream('biofuel', price=price_dct['diesel'])
    BiofuelStorage = qsu.StorageTank(
        'BiofuelStorage',
        BiofuelFlash-1, outs=biofuel,
        tau=24*3, vessel_material='Stainless steel',
        include_construction=False,
        )
    def adjust_biofuel_price():
        BiofuelStorage._run()
        GGE = biofuel.HHV/1e3/_HHV_per_GGE # MJ/gal
        GGE_per_kg = GGE/biofuel.F_mass
        biofuel.price = price_dct['diesel'] * GGE_per_kg
    BiofuelStorage.add_specification(adjust_biofuel_price)
    
    biobinder = qs.WasteStream('biobinder', price=price_dct['biobinder'])
    BiobinderStorage = qsu.StorageTank(
        'HeavyFracStorage', BiobinderHX-0, outs=biobinder,
        tau=24*3, vessel_material='Stainless steel',
        include_construction=False,
        )
    
    # Other co-products
    recovered_H2 = qs.WasteStream('recovered_H2', price=price_dct['H2'])
    H2mixer = qsu.Mixer('H2mixer', ins=H2_streams, outs=recovered_H2)

    recovered_N = qs.WasteStream('recovered_N', price=price_dct['N'])
    Nmixer = qsu.Mixer('Nmixer', ins=N_streams, outs=recovered_N)

    recovered_P = qs.WasteStream('recovered_P', price=price_dct['P'])
    Pmixer = qsu.Mixer('Pmixer', ins=P_streams, outs=recovered_P)
    
    recovered_K = qs.WasteStream('recovered_K', price=price_dct['K'])
    Kmixer = qsu.Mixer('Kmixer', ins=K_streams, outs=recovered_K)

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
    
    # Potentially recycle the water from aqueous filtration (will be ins[2])
    # Can ignore this if the feedstock moisture if close to desired range
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

    sys = qs.System.from_units(
        'sys',
    units=list(flowsheet.unit),
    operating_hours=365 * 24 * uptime_ratio,
     )

    for unit in sys.units:
     unit.include_construction = False


    streams_with_impacts = [
    i for i in sys.feeds + sys.products
    if i.isempty() is False and i.imass['Water'] != i.F_mass and 'surrogate' not in i.ID
    ]

    tea = create_tea(sys, **tea_kwargs)


    clear_lca_registries()
    # qs.ImpactIndicator.load_from_file(os.path.join(data_path, 'impact_indicators.csv'))
    # qs.ImpactItem.load_from_file(os.path.join(data_path, 'impact_items.xlsx'))
    # print("Loaded Impact Items:")
  
    gwp_dict = {
        'feedstock': 0,
        'landfill': 400/1e3, # nearly 400 kg CO2e/tonne, Nordahl et al., 2020
        'composting': -41/1e3, # -41 kg CO2e/tonne, Nordahl et al., 2020
        'anaerobic_digestion': (-36-2)/2, # -36 to -2 kg CO2e/tonne, Nordahl et al., 2020
        'trans_feedstock': 0.011856, # 78 km, https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/9393/impact_assessment, Snowden-Swan PNNL 32731
        'trans_biocrude': 0.024472, # 100 miles,https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/9393/impact_assessment, Snowden-Swan PNNL 32731
        'H2': -10.71017675, # https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/24913/impact_assessment
        'natural_gas': 0.780926344, # https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/4866/impact_assessment
        'process_water': 0,
        'electricity': 0.465474829, # https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/13670/impact_assessment
        'steam': 0.126312684, # kg CO2e/MJ, https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/7479/impact_assessment
        'cooling': 0.068359242, # https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/14408/impact_assessment
        'diesel': -0.801163967, # kg CO2e/kg, https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/13381/impact_assessment
        'N': -0.441913058, #liquid, https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/11489/impact_assessment
        'P': -1.344*(98/31), # H3PO4, https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/8421/impact_assessment
        'K': -4.669210326*(56/39), # KOH, https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/5111/impact_assessment
        'COD': 1.7, # Li et al., 2023
        'wastewater': 0.477724554/1e3, # kg CO2e/m3, https://ecoquery.ecoinvent.org/3.10/cutoff/dataset/26546/impact_assessment
        }

    GWP = qs.ImpactIndicator('GWP',
                              alias='GlobalWarmingPotential',
                              method='Ecoinvent',
                              category='environmental impact',
                              unit='kg CO2-eq',)
           
    feedstock_item = qs.StreamImpactItem(
        ID='feedstock_item',
        linked_stream=scaled_feedstock,
        GWP=gwp_dict['feedstock'],
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
    recovered_N_item = qs.StreamImpactItem(
        ID='recovered_N_item',
        linked_stream=recovered_N,
        GWP=gwp_dict['N'],
        )
    recovered_P_item = qs.StreamImpactItem(
        ID='recovered_P_item',
        linked_stream=recovered_P,
        GWP=gwp_dict['P'],
        )
    recovered_K_item = qs.StreamImpactItem(
        ID='recovered_K_item',
        linked_stream=recovered_K,
        GWP=gwp_dict['K'],
        )
    recovered_H2_item = qs.StreamImpactItem(
        ID='recovered_H2_item',
        linked_stream=recovered_H2,
        GWP=gwp_dict['H2'],
        )
    biofuel_item = qs.StreamImpactItem(
        ID='biofuel_item',
        linked_stream=biofuel,
        GWP=gwp_dict['diesel'],
        )
    e_item = qs.ImpactItem(
        ID='e_item',
        GWP=gwp_dict['electricity'],
        )
    steam_item = qs.ImpactItem(
        ID='steam_item',
        GWP=gwp_dict['steam'],
        )
    cooling_item = qs.ImpactItem(
        ID='cooling_item',
        GWP=gwp_dict['cooling'],
        )
    # for item in qs.ImpactItem.registry:
    #   print(f"- ID: {item.ID}, Functional Unit: {item.functional_unit}")
    fake_stream= qs.SanStream('Nothing', price=0)
    nothing_item = qs.StreamImpactItem(
        ID='nothing_item',
        linked_stream=fake_stream,
        GWP=0,
        )
    def adjust_lca_items():
   
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
           trans_biocrude_item.linked_stream = fake_stream             # No biocrude transportation

    adjust_lca_items()
    lifetime = tea_kwargs['duration'][1] - tea_kwargs['duration'][0]
            
    lca = qs.LCA(
    system=sys,
    lifetime=lifetime,
    simulate_system=False,
    uptime_ratio=sys.operating_hours / (365 * 24),
    e_item=lambda: (sys.get_electricity_consumption() - sys.get_electricity_production()) * lifetime,
    steam_item=lambda: sys.get_heating_duty() / 1000 * lifetime,
    cooling_item=lambda: sys.get_cooling_duty() / 1000 * lifetime,
        )
  
    return sys


def simulate_and_print(sys, save_report=False):
    sys.simulate()
    tea = sys.TEA
    # lca = sys.LCA
    biobinder = sys.flowsheet.stream.biobinder

    biobinder.price = MSP = tea.solve_price(biobinder)
    print(f'Minimum selling price of the biobinder is ${MSP:.2f}/kg.')
        
    all_impacts = lca.get_allocated_impacts(streams=(biobinder,), operation_only=True, annual=True)
    GWP = all_impacts['GWP']/(biobinder.F_mass*lca.system.operating_hours)

    print(f'Global warming potential of the biobinder is {GWP:.4f} kg CO2e/kg.')
    if save_report:
        # Use `results_path` and the `join` func can make sure the path works for all users
        sys.save_report(file=os.path.join(results_path, f'{sys.ID}.xlsx'))

if __name__ == '__main__':
    config_kwargs = dict(
        flowsheet=None, 
        central_dry_flowrate=None,
        pilot_dry_flowrate=None,
        )

    # What to do with HTL-AP
    config_kwargs.update(dict(skip_EC=False, generate_H2=False,))
    # config_kwargs.update(dict(skip_EC=False, generate_H2=True,))
    # config_kwargs.update(dict(skip_EC=True, generate_H2=False,))
    
    # Decentralized vs. centralized configuration
    # config_kwargs.update(dict(decentralized_HTL=False, decentralized_upgrading=False))
    config_kwargs.update(dict(decentralized_HTL=True, decentralized_upgrading=False))

    #config_kwargs.update(dict(decentralized_HTL=True, decentralized_upgrading=True))

    
    # Distillation column cost calculation doesn't scale down well, so the cost is very high now.
    # But maybe don't need to do the DHDU scenario, if DHCU isn't too different from CHCU
    # However, maybe the elimination of transportation completely will make a difference
    # config_kwargs.update(dict(decentralized_HTL=True, decentralized_upgrading=True))
    
    sys = create_system(**config_kwargs)
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    tea = sys.TEA
    lca = sys.LCA    
    
    simulate_and_print(sys)
