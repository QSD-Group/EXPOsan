# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:46:11 2024

@author: Yalin Li
"""

def create_CHCU_system(
        flowsheet=None,
        total_dry_flowrate=central_dry_flowrate, # 110 tpd
        ):
    '''
    Centralized HTL and upgrading.
    '''
    N_HTL = N_upgrading = 1
    _load_process_settings()
    flowsheet = qs.Flowsheet('bb_CHCU')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    _load_components()
    
    feedstock = qs.WasteStream('feedstock', price=price_dct['tipping'])

    FeedstockTrans = safu.Transportation(
        'FeedstockTrans',
        ins=(feedstock, 'feedstock_trans_surrogate'),
        outs=('transported_feedstock',),
        N_unit=N_HTL,
        copy_ins_from_outs=True,
        transportation_unit_cost=0, # will be adjusted later
        transportation_distance=1, # 25 km ref [1]; #!!! changed, distance already included in the unit cost
        )

    FeedstockCond = safu.Conditioning(
        'FeedstockCond', ins=(FeedstockTrans-0, 'htl_process_water'),
        outs='conditioned_feedstock',
        feedstock_composition=feedstock_composition,
        feedstock_dry_flowrate=total_dry_flowrate,
        )
    FeedstockCond.N_unit = N_HTL # init doesn't take this property

    HTL_kwargs = dict(
        ID='HTL',
        ins=FeedstockCond.outs[0],
        outs=('gas','HTL_aqueous','biocrude','hydrochar'),
        dw_yields=HTL_yields,
        gas_composition={
            'CH4':0.050,
            'C2H6':0.032,
            'CO2':0.918
            },
        aqueous_composition={'HTLaqueous': 1},
        biocrude_composition={
            'Water': 0.063,
            'HTLbiocrude': 1-0.063,
            },
        char_composition={'HTLchar': 1},
        internal_heat_exchanging=True,
        )
    
    HTL_unit = u.CentralizedHTL
    HTL = HTL_unit(**HTL_kwargs)

    BiocrudeTrans = safu.Transportation(
        'BiocrudeTrans',
        ins=(HTL-2, 'biocrude_trans_surrogate'),
        outs=('transported_biocrude',),
        N_unit=N_HTL,
        transportation_unit_cost=0, # no need to transport biocrude
        transportation_distance=1, # cost considered average transportation distance
        )
    
    BiocrudeScaler = u.Scaler(
        'BiocrudeScaler', ins=BiocrudeTrans-0, outs='scaled_biocrude',
        scaling_factor=N_HTL, reverse=False,
        )
    
    BiocrudeSplitter = safu.BiocrudeSplitter(
        'BiocrudeSplitter', ins=BiocrudeScaler-0, outs='splitted_biocrude',
        cutoff_Tb=343+273.15, light_frac=0.5316)
    
    CrudePump = qsu.Pump('CrudePump', init_with='Stream', 
                         ins=BiocrudeSplitter-0, outs='crude_to_dist',
                         P=1530.0*_psi_to_Pa,)
    # Jones 2014: 1530.0 psia
    
    CrudeSplitter = safu.BiocrudeSplitter(
        'CrudeSplitter', ins=CrudePump-0, outs='splitted_crude',
        biocrude_IDs=('HTLbiocrude'),
        cutoff_fracs=[0.0339, 0.8104+0.1557], # light (water): medium/heavy (biocrude/char)
        cutoff_Tbs=(150+273.15, ),
        )
    
    # Separate water from organics (bp<150°C)
    CrudeLightDis = qsu.ShortcutColumn(
        'CrudeLightDis', ins=CrudeSplitter-0,
        outs=('crude_light','crude_medium_heavy'),
        LHK=CrudeSplitter.keys[0],
        P=50*_psi_to_Pa,
        Lr=0.87,
        Hr=0.98,
        k=2, is_divided=True)
    
    CrudeLightFlash = qsu.Flash('CrudeLightFlash', ins=CrudeLightDis-0, T=298.15, P=101325,)
    
    # Separate biocrude from biobinder
    CrudeHeavyDis = qsu.ShortcutColumn(
        'CrudeHeavyDis', ins=CrudeLightDis-1,
        outs=('hot_biofuel','hot_biobinder'),
        LHK=('Biofuel', 'Biobinder'),
        P=50*_psi_to_Pa,
        Lr=0.89,
        Hr=0.85,
        k=2, is_divided=True)

    BiofuelFlash = qsu.Flash('BiofuelFlash', ins=CrudeHeavyDis-0, outs=('', 'cooled_biofuel',),
                              T=298.15, P=101325)
    
    BiobinderHX = qsu.HXutility('BiobinderHX', ins=CrudeHeavyDis-1, outs=('cooled_biobinder',),
                                T=298.15)
    
    ww_to_disposal = qs.WasteStream('ww_to_disposal')
    WWmixer = qsu.Mixer('WWmixer',
                        ins=(HTL.outs[1], CrudeLightFlash.outs[1], BiofuelFlash.outs[0]),
                        outs='HTL_aq',)
    @WWmixer.add_specification
    def adjust_prices():
        # Centralized HTL and upgrading, transport feedstock
        dw_price = price_dct['trans_feedstock']
        factor = 1 - FeedstockTrans.ins[0].imass['Water']/FeedstockTrans.ins[0].F_mass
        FeedstockTrans.transportation_unit_cost = dw_price * factor
        BiocrudeTrans.transportation_unit_cost = 0
        
        # Wastewater
        COD_mass_content = sum(ww_to_disposal.imass[i.ID]*i.i_COD for i in ww_to_disposal.components)
        factor = COD_mass_content/ww_to_disposal.F_mass
        ww_to_disposal.price = min(price_dct['wastewater'], price_dct['COD']*factor)
    WWmixer.add_specification(adjust_prices)

    HTL_EC = safu.Electrochemical(
        'HTL_EC',
        ins=(WWmixer-0, 'HTL_EC_replacement_surrogate'),
        outs=('HTL_EC_gas', 'HTL_EC_H2', 'HTL_EC_N', 'HTL_EC_P', 'HTL_EC_K', ww_to_disposal),
        )
    HTL_EC.N_unit = N_upgrading
    
    # 3-day storage time as in the SAF module
    biofuel = qs.WasteStream('biofuel', price=price_dct['diesel'])
    BiofuelStorage = qsu.StorageTank(
        'BiofuelStorage',
        BiofuelFlash-1, outs=biofuel,
        tau=24*3, vessel_material='Stainless steel',
        include_construction=False,
        )
    
    biobinder = qs.WasteStream('biobinder', price=price_dct['biobinder'])
    BiobinderStorage = qsu.StorageTank(
        'HeavyFracStorage', BiobinderHX-0, outs=biobinder,
        tau=24*3, vessel_material='Stainless steel',
        include_construction=False,
        )
    
    natural_gas = qs.WasteStream('natural_gas', CH4=1, price=price_dct['natural_gas'])
    solids_to_disposal = qs.WasteStream('solids_to_disposal', price=price_dct['solids'])
    CHPMixer = qsu.Mixer('CHPMixer',
                         ins=(HTL.outs[0], CrudeLightFlash.outs[0], HTL.outs[-1]),
                         )
    CHP = qsu.CombinedHeatPower('CHP', 
                                ins=CHPMixer.outs[0],
                                outs=('gas_emissions', solids_to_disposal),
                                init_with='WasteStream',
                                supplement_power_utility=False)
    
    # Potentially recycle the water from aqueous filtration (will be ins[2])
    # Can ignore this if the feedstock moisture if close to desired range
    PWC = safu.ProcessWaterCenter(
        'ProcessWaterCenter',
        process_water_streams=FeedstockCond.ins[0],
        process_water_price=price_dct['process_water']
        )
    PWC.register_alias('PWC')

    sys = qs.System.from_units(
        'sys',
        units=list(flowsheet.unit),
        operating_hours=365*24*uptime_ratio,
        )
    for unit in sys.units: unit.include_construction = False
    
    tea = create_tea(sys, **tea_kwargs)
    
    return sys
        
def create_DHCU_system(
        flowsheet=None,
        total_dry_flowrate=central_dry_flowrate, # 110 tpd
        ):
    '''
    Decentralized HTL, centralized upgrading.
    '''
    N_HTL = round(total_dry_flowrate/pilot_dry_flowrate)
    N_upgrading = 1
    _load_process_settings()
    flowsheet = qs.Flowsheet('bb_DHCU')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    _load_components()
    
    feedstock = qs.WasteStream('feedstock', price=price_dct['tipping'])

    FeedstockTrans = safu.Transportation(
        'FeedstockTrans',
        ins=(feedstock, 'feedstock_trans_surrogate'),
        outs=('transported_feedstock',),
        N_unit=N_HTL,
        copy_ins_from_outs=True,
        transportation_unit_cost=0, # will be adjusted later
        transportation_distance=1, # 25 km ref [1]; #!!! changed, distance already included in the unit cost
        )

    FeedstockCond = safu.Conditioning(
        'FeedstockCond', ins=(FeedstockTrans-0, 'htl_process_water'),
        outs='conditioned_feedstock',
        feedstock_composition=feedstock_composition,
        feedstock_dry_flowrate=total_dry_flowrate,
        )
    FeedstockCond.N_unit = N_HTL # init doesn't take this property

    HTL_kwargs = dict(
        ID='HTL',
        ins=FeedstockCond.outs[0],
        outs=('gas','HTL_aqueous','biocrude','hydrochar'),
        dw_yields=HTL_yields,
        gas_composition={
            'CH4':0.050,
            'C2H6':0.032,
            'CO2':0.918
            },
        aqueous_composition={'HTLaqueous': 1},
        biocrude_composition={
            'Water': 0.063,
            'HTLbiocrude': 1-0.063,
            },
        char_composition={'HTLchar': 1},
        internal_heat_exchanging=True,
        )
    
    HTL_unit = u.PilotHTL
    HTL = HTL_unit(**HTL_kwargs)

    HTL_EC = safu.Electrochemical(
        'HTL_EC',
        ins=(HTL-1, 'HTL_EC_replacement_surrogate'),
        outs=('HTL_EC_gas', 'HTL_EC_H2', 'HTL_EC_N', 'HTL_EC_P', 'HTL_EC_K', 'HTL_EC_WW'),
        )
    HTL_EC.N_unit = N_HTL

    BiocrudeTrans = safu.Transportation(
        'BiocrudeTrans',
        ins=(HTL-2, 'biocrude_trans_surrogate'),
        outs=('transported_biocrude',),
        N_unit=N_HTL,
        transportation_unit_cost=0, # no need to transport biocrude
        transportation_distance=1, # cost considered average transportation distance
        )
    
    BiocrudeScaler = u.Scaler(
        'BiocrudeScaler', ins=BiocrudeTrans-0, outs='scaled_biocrude',
        scaling_factor=N_HTL, reverse=False,
        )
    
    BiocrudeSplitter = safu.BiocrudeSplitter(
        'BiocrudeSplitter', ins=BiocrudeScaler-0, outs='splitted_biocrude',
        cutoff_Tb=343+273.15, light_frac=0.5316)
    
    CrudePump = qsu.Pump('CrudePump', init_with='Stream', 
                         ins=BiocrudeSplitter-0, outs='crude_to_dist',
                         P=1530.0*_psi_to_Pa,)
    # Jones 2014: 1530.0 psia
    
    CrudeSplitter = safu.BiocrudeSplitter(
        'CrudeSplitter', ins=CrudePump-0, outs='splitted_crude',
        biocrude_IDs=('HTLbiocrude'),
        cutoff_fracs=[0.0339, 0.8104+0.1557], # light (water): medium/heavy (biocrude/char)
        cutoff_Tbs=(150+273.15, ),
        )
    
    # Separate water from organics (bp<150°C)
    CrudeLightDis = qsu.ShortcutColumn(
        'CrudeLightDis', ins=CrudeSplitter-0,
        outs=('crude_light','crude_medium_heavy'),
        LHK=CrudeSplitter.keys[0],
        P=50*_psi_to_Pa,
        Lr=0.87,
        Hr=0.98,
        k=2, is_divided=True)
    
    CrudeLightFlash = qsu.Flash('CrudeLightFlash', ins=CrudeLightDis-0, T=298.15, P=101325,)
    
    # Separate biocrude from biobinder
    CrudeHeavyDis = qsu.ShortcutColumn(
        'CrudeHeavyDis', ins=CrudeLightDis-1,
        outs=('hot_biofuel','hot_biobinder'),
        LHK=('Biofuel', 'Biobinder'),
        P=50*_psi_to_Pa,
        Lr=0.89,
        Hr=0.85,
        k=2, is_divided=True)

    BiofuelFlash = qsu.Flash('BiofuelFlash', ins=CrudeHeavyDis-0, outs=('', 'cooled_biofuel',),
                              T=298.15, P=101325)
    
    BiobinderHX = qsu.HXutility('BiobinderHX', ins=CrudeHeavyDis-1, outs=('cooled_biobinder',),
                                T=298.15)
    
    WWmixer = qsu.Mixer('WWmixer',
                        ins=(CrudeLightFlash-1, BiofuelFlash-0,),
                        outs='HTL_aq',
                        )
    
    Upgrading_EC = safu.Electrochemical(
        'Upgrading_EC',
        ins=(WWmixer-0, 'HTL_EC_replacement_surrogate'),
        outs=('Upgrading_EC_gas', 'Upgrading_EC_H2', 'Upgrading_EC_N', 'Upgrading_EC_P', 'Upgrading_EC_K', 'Upgrading_EC_WW'),
        )
    HTL_EC.N_unit = N_upgrading
    
    # 3-day storage time as in the SAF module
    biofuel = qs.WasteStream('biofuel', price=price_dct['diesel'])
    BiofuelStorage = qsu.StorageTank(
        'BiofuelStorage',
        BiofuelFlash-1, outs=biofuel,
        tau=24*3, vessel_material='Stainless steel',
        include_construction=False,
        )
    
    biobinder = qs.WasteStream('biobinder', price=price_dct['biobinder'])
    BiobinderStorage = qsu.StorageTank(
        'HeavyFracStorage', BiobinderHX-0, outs=biobinder,
        tau=24*3, vessel_material='Stainless steel',
        include_construction=False,
        )
    
    natural_gas = qs.WasteStream('natural_gas', CH4=1, price=price_dct['natural_gas'])
    CHPMixer = qsu.Mixer('CHPMixer', ins=CrudeLightFlash.outs[0])
    CHP = qsu.CombinedHeatPower('CHP', 
                                ins=CHPMixer.outs[0],
                                outs=('gas_emissions', 'CHP_solids_to_disposal'),
                                init_with='WasteStream',
                                supplement_power_utility=False)
    CHP.outs[1].price = 0 # ash disposal will be accounted for in SolidsDisposal
    
    # Potentially recycle the water from aqueous filtration (will be ins[2])
    # Can ignore this if the feedstock moisture if close to desired range
    PWC = safu.ProcessWaterCenter(
        'ProcessWaterCenter',
        process_water_streams=FeedstockCond.ins[0],
        process_water_price=price_dct['process_water']
        )
    PWC.register_alias('PWC')
    
    solids_to_disposal = qs.WasteStream('solids_to_disposal', price=price_dct['solids'])
    SolidsDisposal = qsu.Mixer('SolidsDisposal',
                               ins=(HTL.outs[-1], CHP.outs[1]),
                               outs=solids_to_disposal,)
    
    ww_to_disposal = qs.WasteStream('ww_to_disposal')
    WWdisposalMixer = qsu.Mixer('WWdisposalMixer',
                                ins=(HTL_EC.outs[-1], Upgrading_EC.outs[-1]),
                                outs=ww_to_disposal)
    @WWdisposalMixer.add_specification
    def adjust_prices():
        # Decentralized HTL, centralized upgrading, transport biocrude
        FeedstockTrans.transportation_unit_cost = 0
        GGE_price = price_dct['trans_biocrude'] # $/GGE
        factor = BiocrudeTrans.ins[0].HHV/_HHV_per_GGE/BiocrudeTrans.ins[0].F_mass
        BiocrudeTrans.transportation_unit_cost = GGE_price * factor #!!! need to check the calculation
        
        # Wastewater
        WWdisposalMixer._run()
        COD_mass_content = sum(ww_to_disposal.imass[i.ID]*i.i_COD for i in ww_to_disposal.components)
        factor = COD_mass_content/ww_to_disposal.F_mass
        ww_to_disposal.price = min(price_dct['wastewater'], price_dct['COD']*factor)
    WWdisposalMixer.add_specification(adjust_prices)
    
    sys = qs.System.from_units(
        'sys',
        units=list(flowsheet.unit),
        operating_hours=365*24*uptime_ratio,
        )
    for unit in sys.units: unit.include_construction = False
    
    tea = create_tea(sys, **tea_kwargs)
    
    return sys

def create_DHDU_system(
        flowsheet=None,
        total_dry_flowrate=pilot_dry_flowrate, # 11.46 dry kg/hr
        ):
    sys = create_CHCU_system(
        flowsheet=flowsheet,
        total_dry_flowrate=pilot_dry_flowrate,
        )
    
    return sys

def create_system(
        flowsheet=None,
        total_dry_flowrate=central_dry_flowrate, # 110 tpd
        decentralized_HTL=False,
        decentralized_upgrading=False,        
        ):
    kwargs = dict(
        flowsheet=None,
        total_dry_flowrate=central_dry_flowrate,
        )
    if decentralized_HTL is False:
        if decentralized_upgrading is False:
            f = create_CHCU_system
        else:
            raise ValueError('Centralized HTL, decentralized upgrading is not a valid configuration.')
    else:
        if decentralized_upgrading is False:
            f = create_DHCU_system
        else:
            f = create_DHDU_system
    sys = f(**kwargs)
    return sys