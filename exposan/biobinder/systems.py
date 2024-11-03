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
import warnings
warnings.filterwarnings('ignore')

import os, biosteam as bst, qsdsan as qs
# from biosteam.units import IsenthalpicValve
# from biosteam import settings
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import create_tea
from exposan.saf import _units as safu
from exposan.biobinder import (
    _HHV_per_GGE,
    _load_components,
    _load_process_settings,
    _units as u,
    central_dry_flowrate,
    data_path,
    feedstock_composition,
    HTL_yields,
    pilot_dry_flowrate,
    price_dct,
    results_path,
    tea_kwargs,
    uptime_ratio,
    )

_psi_to_Pa = 6894.76


# %%

__all__ = ('create_system',)

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
        ins=FeedstockCond.ins[0],
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
    
    gas_streams = [HTL.outs[0]]
    ww_streams = [HTL.outs[1]]
    solids_streams = [HTL.outs[-1]]

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
    ww_streams.append(CrudeLightFlash.outs[0])
    
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
    
    HTLaqMixer = qsu.Mixer('HTLaqMixer', ins=ww_streams, outs='HTL_aq')
    
    HTL_EC = safu.Electrochemical(
        'HTL_EC',
        ins=(HTLaqMixer-0, 'HTL_EC_replacement_surrogate'),
        outs=('HTL_EC_gas', 'HTL_EC_H2', 'HTL_EC_N', 'HTL_EC_P', 'HTL_EC_K', 'HTL_EC_WW'),
        )
    HTL_EC.N_unit = N_upgrading
    
    # 3-day storage time as in the SAF module
    biofuel = qs.WasteStream('biofuel', price=price_dct['diesel'])
    BiofuelStorage = qsu.StorageTank(
        'BiofuelStorage',
        BiofuelFlash-1, outs=biofuel,
        tau=24*3, vessel_material='Stainless steel')
    
    biobinder = qs.WasteStream('biobinder', price=price_dct['biobinder'])
    BiobinderStorage = qsu.StorageTank(
        'HeavyFracStorage', BiobinderHX-0, outs=biobinder,
        tau=24*3, vessel_material='Stainless steel')
    
    natural_gas = qs.WasteStream('natural_gas', CH4=1, price=price_dct['natural_gas'])
    solids_to_disposal = qs.WasteStream('solids_to_disposal', price=price_dct['solids'])
    CHPMixer = qsu.Mixer('CHPMixer', ins=gas_streams+solids_streams)
    CHP = qsu.CombinedHeatPower('CHP', 
                                ins=gas_streams+solids_streams,
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
    
    solids_to_disposal = qs.WasteStream('solids_to_disposal')
    SolidsDisposal = qsu.Mixer('SolidsDisposal', ins=solids_streams, outs=(solids_to_disposal, ),)
    
    ww_to_disposal = qs.WasteStream('ww_to_disposal')
    WWmixer = qsu.Mixer('WWmixer', ins=HTL_EC.outs[-1], outs=ww_to_disposal)
    @WWmixer.add_specification
    def adjust_prices():
        # Centralized HTL and upgrading, transport feedstock
        dw_price = price_dct['trans_feedstock']
        factor = 1 - FeedstockTrans.ins[0].imass['Water']/FeedstockTrans.ins[0].F_mass
        FeedstockTrans.transportation_unit_cost = dw_price * factor
        BiocrudeTrans.transportation_unit_cost = 0
        
        # Wastewater
        WWmixer._run()
        COD_mass_content = sum(ww_to_disposal.imass[i.ID]*i.i_COD for i in ww_to_disposal.components)
        factor = COD_mass_content/ww_to_disposal.F_mass
        ww_to_disposal.price = min(price_dct['wastewater'], price_dct['COD']*factor)
    ww_to_disposal.source.add_specification(adjust_prices)
    
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
        ins=FeedstockCond.ins[0],
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
    
    gas_streams = [HTL.outs[0]]
    solids_streams = [HTL.outs[-1]]

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
    
    HTLaqMixer = qsu.Mixer('HTLaqMixer', ins=CrudeLightFlash.outs[0], outs='HTL_aq')
    
    Upgrading_EC = safu.Electrochemical(
        'HTL_EC',
        ins=(HTLaqMixer-0, 'HTL_EC_replacement_surrogate'),
        outs=('Upgrading_EC_gas', 'Upgrading_EC_H2', 'Upgrading_EC_N', 'Upgrading_EC_P', 'Upgrading_EC_K', 'Upgrading_EC_WW'),
        )
    HTL_EC.N_unit = N_upgrading
    
    # 3-day storage time as in the SAF module
    biofuel = qs.WasteStream('biofuel', price=price_dct['diesel'])
    BiofuelStorage = qsu.StorageTank(
        'BiofuelStorage',
        BiofuelFlash-1, outs=biofuel,
        tau=24*3, vessel_material='Stainless steel')
    
    biobinder = qs.WasteStream('biobinder', price=price_dct['biobinder'])
    BiobinderStorage = qsu.StorageTank(
        'HeavyFracStorage', BiobinderHX-0, outs=biobinder,
        tau=24*3, vessel_material='Stainless steel')
    
    natural_gas = qs.WasteStream('natural_gas', CH4=1, price=price_dct['natural_gas'])
    solids_to_disposal = qs.WasteStream('solids_to_disposal', price=price_dct['solids'])
    CHPMixer = qsu.Mixer('CHPMixer', ins=gas_streams+solids_streams)
    CHP = qsu.CombinedHeatPower('CHP', 
                                ins=gas_streams+solids_streams,
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
    
    solids_to_disposal = qs.WasteStream('solids_to_disposal')
    SolidsDisposal = qsu.Mixer('SolidsDisposal', ins=solids_streams, outs=(solids_to_disposal, ),)
    
    ww_to_disposal = qs.WasteStream('ww_to_disposal')
    WWmixer = qsu.Mixer('WWmixer', ins=HTL_EC.outs[-1], outs=ww_to_disposal)
    @WWmixer.add_specification
    def adjust_prices():
        # Decentralized HTL, centralized upgrading, transport biocrude
        FeedstockTrans.transportation_unit_cost = 0
        GGE_price = price_dct['trans_biocrude'] # $/GGE
        factor = BiocrudeTrans.ins[0].HHV/_HHV_per_GGE/BiocrudeTrans.ins[0].F_mass
        BiocrudeTrans.transportation_unit_cost = GGE_price * factor #!!! need to check the calculation
        
        # Wastewater
        WWmixer._run()
        COD_mass_content = sum(ww_to_disposal.imass[i.ID]*i.i_COD for i in ww_to_disposal.components)
        factor = COD_mass_content/ww_to_disposal.F_mass
        ww_to_disposal.price = min(price_dct['wastewater'], price_dct['COD']*factor)
    ww_to_disposal.source.add_specification(adjust_prices)
    
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


def simulate_and_print(sys, save_report=False):
    sys.simulate()
    tea = sys.TEA
    lca = sys.LCA
    biobinder = sys.flowsheet.stream.biobinder

    biobinder.price = MSP = tea.solve_price(biobinder)
    print(f'Minimum selling price of the biobinder is ${MSP:.2f}/kg.')
        
    all_impacts = lca.get_allocated_impacts(streams=(biobinder,), operation_only=True, annual=True)
    GWP = all_impacts['GlobalWarming']/(biobinder.F_mass*lca.system.operating_hours)
    print(f'Global warming potential of the biobinder is {GWP:.4f} kg CO2e/kg.')
    if save_report:
        # Use `results_path` and the `join` func can make sure the path works for all users
        sys.save_report(file=os.path.join(results_path, f'{sys.ID}.xlsx'))

if __name__ == '__main__':
    # sys = create_system(decentralized_HTL=False, decentralized_upgrading=False)
    sys = create_system(decentralized_HTL=True, decentralized_upgrading=False)
    # sys = create_system(decentralized_HTL=True, decentralized_upgrading=True)
    # simulate_and_print(sys)

    # # Separate biocrude from biobinder
    # CrudeHeavyDis = qsu.ShortcutColumn(
    #     'CrudeHeavyDis', ins=CrudeLightDis-1,
    #     outs=('crude_medium','char'),
    #     LHK=CrudeSplitter.keys[1],
    #     P=50*_psi_to_Pa,
    #     Lr=0.89,
    #     Hr=0.85,
    #     k=2, is_divided=True)

    # CrudeHeavyDis_run = CrudeHeavyDis._run
    # CrudeHeavyDis_design = CrudeHeavyDis._design
    # CrudeHeavyDis_cost = CrudeHeavyDis._cost
    # def run_design_cost():
    #     CrudeHeavyDis_run()
    #     try:
    #         CrudeHeavyDis_design()
    #         CrudeHeavyDis_cost()
    #         if all([v>0 for v in CrudeHeavyDis.baseline_purchase_costs.values()]):
    #             # Save for later debugging
    #             # print('design')
    #             # print(CrudeHeavyDis.design_results)
    #             # print('cost')
    #             # print(CrudeHeavyDis.baseline_purchase_costs)
    #             # print(CrudeHeavyDis.installed_costs) # this will be empty
    #             return
    #     except: pass
    #     raise RuntimeError('`CrudeHeavyDis` simulation failed.')

    
    # # Simulation may converge at multiple points, filter out unsuitable ones
    # def screen_results():
    #     ratio0 = CrudeSplitter.cutoff_fracs[1]/sum(CrudeSplitter.cutoff_fracs[1:])
    #     lb, ub = round(ratio0,2)-0.02, round(ratio0,2)+0.02
    #     try: 
    #         run_design_cost()
    #         status = True
    #     except: 
    #         status = False
    #     def get_ratio():
    #         if CrudeHeavyDis.F_mass_out > 0: 
    #             return CrudeHeavyDis.outs[0].F_mass/CrudeHeavyDis.F_mass_out
    #         return 0
    #     n = 0
    #     ratio = get_ratio()
    #     while (status is False) or (ratio<lb) or (ratio>ub):
    #         try: 
    #             run_design_cost()
    #             status = True
    #         except: 
    #             status = False
    #         ratio = get_ratio()
    #         n += 1
    #         if n >= 10:
    #             status = False
    #             raise RuntimeError(f'No suitable solution for `CrudeHeavyDis` within {n} simulation.')
    # CrudeHeavyDis._run = screen_results

    # def do_nothing(): pass
    # CrudeHeavyDis._design = CrudeHeavyDis._cost = do_nothing
    
    # # Lr_range = Hr_range = np.arange(0.05, 1, 0.05)
    # # results = find_Lr_Hr(CrudeHeavyDis, target_light_frac=crude_char_fracs[0], Lr_trial_range=Lr_range, Hr_trial_range=Hr_range)
    # # results_df, Lr, Hr = results

    # # Shortcut column uses the Fenske-Underwood-Gilliland method,
    # # better for hydrocarbons according to the tutorial
    # # https://biosteam.readthedocs.io/en/latest/API/units/distillation.html
    # FracDist = qsu.ShortcutColumn(
    #     'FracDist', ins=BiocrudeSplitter-0,
    #     outs=('biocrude_light','biocrude_heavy'),
    #     LHK=('Biofuel', 'Biobinder'), # will be updated later
    #     P=50*6894.76, # outflow P, 50 psig
    #     # Lr=0.1, Hr=0.5,
    #     y_top=188/253, x_bot=53/162,
    #     k=2, is_divided=True)
    # @FracDist.add_specification
    # def adjust_LHK():
    #     FracDist.LHK = BiocrudeSplitter.keys[0]
    #     FracDist._run()
        
    # Lr_range = Hr_range = np.linspace(0.05, 0.95, 19)
    # Lr_range = Hr_range = np.linspace(0.01, 0.2, 20)
    # results = find_Lr_Hr(FracDist, Lr_trial_range=Lr_range, Hr_trial_range=Hr_range)
    # results_df, Lr, Hr = results
    
    # def adjust_prices():
    #     # Centralized HTL and upgrading, transport feedstock
    #     if decentralized_HTL is False:
    #         dw_price = price_dct['trans_feedstock']
    #         factor = 1 - FeedstockTrans.ins[0].imass['Water']/FeedstockTrans.ins[0].F_mass
    #         FeedstockTrans.transportation_unit_cost = dw_price * factor
    #         BiocrudeTrans.transportation_unit_cost = 0
    #     # Decentralized HTL, centralized upgrading, transport biocrude
    #     elif decentralized_upgrading is False:
    #         FeedstockTrans.transportation_unit_cost = 0
    #         GGE_price = price_dct['trans_biocrude'] # $/GGE
    #         factor = BiocrudeTrans.ins[0].HHV/_HHV_per_GGE/BiocrudeTrans.ins[0].F_mass
    #         BiocrudeTrans.transportation_unit_cost = GGE_price * factor #!!! need to check the calculation
    #     # Decentralized HTL and upgrading, no transportation needed
    #     else:
    #         FeedstockTrans.transportation_unit_cost = BiocrudeTrans.transportation_unit_cost = 0