

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
[2] Feng et al, Characterizing the Opportunity Space for Sustainable
    Hydrothermal Valorization of Wet Organic Wastes.
    Environ. Sci. Technol. 2024, 58 (5), 2528–2541.
    https://doi.org/10.1021/acs.est.3c07394.
[3] Nordahl et al., Life-Cycle Greenhouse Gas Emissions and Human Health Trade-Offs
    of Organic Waste Management Strategies.
    Environ. Sci. Technol. 2020, 54 (15), 9200–9209.
    https://doi.org/10.1021/acs.est.0c00364.
[4] Lopez-Ruiz et al., Electrocatalytic Valorization into H2 and Hydrocarbons
    of an Aqueous Stream Derived from Hydrothermal Liquefaction.
    J Appl Electrochem 2021, 51 (1), 107–118.
    https://doi.org/10.1007/s10800-020-01452-x.
[5] Lopez-Ruiz et al., Low-Temperature Electrochemical Wastewater Oxidation;
    PNNL-35535; 2023. https://doi.org/10.2172/2332860.


'''

# !!! Temporarily ignoring warnings
import warnings
warnings.filterwarnings('ignore')

import os, numpy as np, biosteam as bst, qsdsan as qs
from biosteam import IsenthalpicValve
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import create_tea
from exposan.saf import (
    _HHV_per_GGE,
    _load_components,
    _load_process_settings,
    _units as u,
    data_path,
    dry_flowrate,
    # feedstock_composition,
    find_Lr_Hr,
    get_mass_energy_balance,
    gwp_dct,
    HTL_yields,
    price_dct,
    results_path,
    tea_kwargs,
    uptime_ratio,
    )
from exposan.biobinder_ml._feedstocks import get_feedstock_composition
from exposan.biobinder_ml import (
    cutoff_fracs_from_feed,
    limit_internal_hx,
    configure_utilities,
    set_eff_hx_temperature,
    )
from exposan.biobinder_ml._distill_utils import add_runtime_LHK_LrHr_spec
_psi_to_Pa = 6894.76
_m3_to_gal = 264.172


# %%

__all__ = (
    'config_baseline',
    'config_EC',
    'config_EC_future',
    'create_system',
    'get_GWP',
    'get_MFSP',
    'simulate_and_print',
    )

def create_system(
        flowsheet=None,
        include_PSA=True,
        include_EC=False,
        feedstock_id='fog',
        dry_flowrate=dry_flowrate,
        # feedstock_composition=feedstock_composition,
        electricitry_price=None,
        electricitry_GHG=None,
        ):
    _load_process_settings()
    _load_components()

    if not flowsheet:
        if not include_PSA:
            flowsheet_ID = 'no_PSA'
        elif include_EC is False: flowsheet_ID = 'baseline'
        elif include_EC is True: flowsheet_ID = 'EC'
        elif type(include_EC) is dict: flowsheet_ID = 'EC_future'
        else: raise ValueError('Invalid system configuration.')
        flowsheet = qs.Flowsheet(flowsheet_ID)
        qs.main_flowsheet.set_flowsheet(flowsheet)
    else:
        qs.main_flowsheet.set_flowsheet(flowsheet)
        
    feedstock_composition = get_feedstock_composition(feedstock_id)
    print(f'Feedstock used is {feedstock_id}')
    print("Corresponding composition:")
    for k, v in feedstock_composition.items():
         print(f"  {k}: {v:.2f}")
    
    feedstock = qs.WasteStream('feedstock', price=price_dct['tipping'])
    feedstock.imass[list(feedstock_composition.keys())] = list(feedstock_composition.values())
    feedstock.F_mass = dry_flowrate / (1-feedstock_composition['Water'])
    
    feedstock_water = qs.Stream('feedstock_water', Water=1)
    
    FeedstockTrans = u.Transportation(
        'FeedstockTrans',
        ins=(feedstock, 'transportation_surrogate'),
        outs=('transported_feedstock',),
        copy_ins_from_outs=False,
        transportation_unit_cost=1, # already considered distance, will be adjusted later
        transportation_distance=1,
        N_unit=1,
        )
    
    FeedstockWaterPump = qsu.Pump('FeedstockWaterPump', ins=feedstock_water)

    FeedstockCond = u.Conditioning(
        'FeedstockCond', ins=(FeedstockTrans-0, FeedstockWaterPump-0),
        outs='conditioned_feedstock',
        feedstock_composition=None,
        feedstock_dry_flowrate=dry_flowrate,
        target_HTL_solid_loading=0.2,
        )
    @FeedstockCond.add_specification
    def adjust_feedstock_composition():
        FeedstockCond._run()
        FeedstockWaterPump._run()
    
    MixedFeedstockPump = qsu.Pump('MixedFeedstockPump', ins=FeedstockCond-0)
    
    # =========================================================================
    # Hydrothermal Liquefaction (HTL)
    # =========================================================================
    aqueous_composition = {
        'N': 0.48/100,
        'P': 0.60/100,
        'K': 0.72/100,
        }
    aqueous_composition['HTLaqueous'] = 1 - sum(aqueous_composition.values())
    HTL = u.HydrothermalLiquefaction(
        'HTL', ins=MixedFeedstockPump-0,
        outs=('', '', 'HTL_crude', 'ash'),
        T=280+273.15,
        P=12.4e6, # may lead to HXN error when HXN is included
        # P=101325, # setting P to ambient pressure not practical, but it has minimum effects on the results (several cents)
        tau=15/60,
        dw_yields=HTL_yields,
        gas_composition={'CO2': 1}, 
        aqueous_composition=aqueous_composition,
        biocrude_composition={'Biocrude': 1},
        char_composition={'HTLchar': 1},
        internal_heat_exchanging=True,
        eff_T=60+273.15, # 140.7°F
        eff_P=30*_psi_to_Pa,
        use_decorated_cost=True,
        )
    HTL.register_alias('HydrothermalLiquefaction')
    limit_internal_hx(HTL, coldest_id='chilled_water', cap_margin_K=2.0, cooler_margin_K=3.0)
    set_eff_hx_temperature(HTL, clamp_to='cooling_water')  # or 'chilled_water'
    
    CrudePump = qsu.Pump('CrudePump', ins=HTL-2, outs='crude_to_dist',
              init_with='Stream')
    
    # Light (water): medium (biocrude): heavy (char)
    original_crude_fracs = [0.01, 0.69, 0.30] # to account for the non-volatile crude fracs (inorganics)
    ratio = (HTL_yields['biocrude']+HTL_yields['char'])/HTL_yields['biocrude']
    crude_fracs = [i*ratio for i in original_crude_fracs[:2]]
    crude_fracs.append(1-sum(crude_fracs))
    # cutoff_fracs = cutoff_fracs_from_feed(feedstock_composition, HTL_yields, dewatered=True)
    
    CrudeSplitter = u.BiocrudeSplitter(
        'CrudeSplitter', ins=CrudePump-0, outs='splitted_crude',
        biocrude_IDs=('HTLbiocrude'),
        cutoff_fracs=crude_fracs,
        cutoff_Tbs=(150+273.15, 300+273.15,),
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
    
    CrudeLightFlash = qsu.Flash('CrudeLightFlash', ins=CrudeLightDis-0,
                                T=298.15, P=101325,)
                                # thermo=settings.thermo.ideal())
    HTLaqMixer = qsu.Mixer('HTLaqMixer', ins=(HTL-1, CrudeLightFlash-1), outs='HTL_aq')
    
    # # Separate fuel from char
    # CrudeHeavyDis = qsu.ShortcutColumn(
    #     'CrudeHeavyDis', ins=CrudeLightDis-1,
    #     outs=('crude_medium','char'),
    #     LHK=('4M-PHYNO', 'INDOLE'),
    #     P=50*_psi_to_Pa,
    #     Lr=0.98, # 0.57, 0.65, 0.97
    #     Hr=0.89,
    #     # y_top=0.95,
    #     # x_bot=0.05,
    #     # product_specification_format='Composition',
    #     k=2,
    #     is_divided=True)
    
    # target_light_frac = cutoff_fracs[1]   # 

    # # Reasonable trial ranges
    # Lr_range = np.linspace(0.6, 1.0, 9)
    # Hr_range = np.linspace(0.8, 1.0, 9)

    # def optimize_CrudeHeavyDis():
    #     """Find and apply best Lr/Hr for CrudeHeavyDis to match predicted cutoff fraction."""
    #     try:
    #       print(f"\n[Optimization] Target biofuel fraction: {target_light_frac:.3f}")
    #       results_df, best_Lr, best_Hr, max_yield = find_Lr_Hr(
    #            CrudeHeavyDis,
    #            Lr_trial_range=Lr_range,
    #            Hr_trial_range=Hr_range,
    #           )
    #       if best_Lr and best_Hr:
    #             CrudeHeavyDis.Lr = best_Lr
    #             CrudeHeavyDis._Lr = best_Lr
    #             CrudeHeavyDis.Hr = best_Hr
    #             CrudeHeavyDis._Hr = best_Hr
    #             print(f"[Optimization] Applied optimal Lr={best_Lr:.3f}, Hr={best_Hr:.3f} "
    #                     f"→ Biofuel yield={max_yield:.3f}")
    #       else:
    #             print("[Optimization] No convergence found; keeping default Lr/Hr.")
    #     except Exception as e:
    #      print(f"[Optimization] Distillation optimization failed: {e}")
    
    # CrudeHeavyDis.check_LHK = False
    # CrudeHeavyDis.add_specification(optimize_CrudeHeavyDis)
    # CrudeHeavyDis.run_after_specifications = True
    
    #Separate biocrude from char
    CrudeHeavyDis = qsu.ShortcutColumn(
        "CrudeHeavyDis",
        ins=CrudeLightDis - 1,
        outs=("crude_medium", "char"),
        LHK=("4M-PHYNO", "INDOLE"),
        P=50 * _psi_to_Pa,
        Lr=0.8,
        Hr=0.85,
        k=2,
        is_divided=True,
    )

    ratio0 = CrudeSplitter.cutoff_fracs[1]/sum(CrudeSplitter.cutoff_fracs[1:])
    add_runtime_LHK_LrHr_spec(
        splitter=CrudeSplitter,
        col=CrudeHeavyDis,
        idx=1,
        fallback_LHK=("4M-PHYNO", "INDOLE"),
        target_ratio=ratio0,
        tol=0.05,              # matches SAF lb/ub = ±0.05
        prefer="min_err",      # closest to original SAF screen-to-window
        require_design=True,
        require_cost=False,
        vle_screen=True,
        design_sanity=True,
    )

    # def get_ratio():
    #     if CrudeHeavyDis.F_mass_out > 0:
    #         return CrudeHeavyDis.outs[0].F_mass/CrudeHeavyDis.F_mass_out
    #     return 0

    # # Simulation may converge at multiple points, filter out unsuitable ones
    # def screen_results():
    #     n = 0
    #     status = False
    #     while (status is False) and (n<20):
    #         try:
    #             CrudeHeavyDis._run()
    #             ratio = get_ratio()
    #             assert(lb<=ratio<=ub)
    #             CrudeHeavyDis._design()
    #             CrudeHeavyDis._cost()
    #             assert(all([v>0 for v in CrudeHeavyDis.baseline_purchase_costs.values()]))
    #             status = True
    #         except:
    #             n += 1
    #             status = False
    #     if n >= 20:
    #         raise RuntimeError(f'No suitable solution for `CrudeHeavyDis` within {n} simulation.')
    # CrudeHeavyDis.add_specification(screen_results)
    # CrudeHeavyDis.run_after_specifications = True
    
    # import numpy as np
    # Lr_range = Hr_range = np.arange(0.85, 0.96, 0.01)
    # ratio0 = CrudeSplitter.cutoff_fracs[1]/sum(CrudeSplitter.cutoff_fracs[1:])
    # results = find_Lr_Hr(CrudeHeavyDis, target_light_frac=ratio0, Lr_trial_range=Lr_range, Hr_trial_range=Hr_range)
    # results_df, Lr, Hr = results
    
    # =========================================================================
    # Hydrocracking
    # =========================================================================
    
    # 10 wt% Fe-ZSM
    HCcatalyst_in = qs.WasteStream('HCcatalyst_in', HCcatalyst=1, price=price_dct['HCcatalyst'])
    
    HC = u.Hydroprocessing(
        'HC',
        ins=(CrudeHeavyDis-0, 'H2_HC', HCcatalyst_in),
        outs=('HC_out','HCcatalyst_out'),
        T=400+273.15,
        P=1500*_psi_to_Pa,
        WHSV=0.625,
        catalyst_ID='HCcatalyst',
        catalyst_lifetime=5*uptime_ratio*365*24, # 5 years [1]
        hydrogen_rxned_to_inf_oil=0.0111, # data from expt
        hydrogen_ratio=5.556,
        include_PSA=include_PSA,
        gas_yield=0.2665,
        oil_yield=0.7335,
        gas_composition={ # [1] after the first hydroprocessing
            'CO2': 1-0.08809, # 0.08809 is the sum of all other gases
            'CH4':0.02280, 'C2H6':0.02923,
            'C3H8':0.01650, 'C4H10':0.00870,
            'TWOMBUTAN':0.00408, 'NPENTAN':0.00678,
            },
        oil_composition={
           'TWOMPENTA':0.00408, 'HEXANE':0.00408,
           'TWOMHEXAN':0.00408, 'HEPTANE':0.00408,
           'CC6METH':0.01020, 'PIPERDIN':0.00408,
           'TOLUENE':0.01020, 'THREEMHEPTA':0.01020,
           'OCTANE':0.01020, 'ETHCYC6':0.00408,
           'ETHYLBEN':0.02040, 'OXYLENE':0.01020,
           'C9H20':0.00408, 'PROCYC6':0.00408,
           'C3BENZ':0.01020, 'FOURMONAN':0,
           'C10H22':0.00203, 'C4BENZ':0.01223,
           'C11H24':0.02040, 'C10H12':0.02040,
           'C12H26':0.02040, 'OTTFNA':0.01020,
           'C6BENZ':0.02040, 'OTTFSN':0.02040,
           'C7BENZ':0.02040, 'C8BENZ':0.02040,
           'C10H16O4':0.01837, 'C15H32':0.06120,
           'C16H34':0.18360, 'C17H36':0.08160, 
           'C18H38':0.04080, 'C19H40':0.04080,
           'C20H42':0.10200, 'C21H44':0.04080,
           'TRICOSANE':0.04080, 'C24H38O4':0.00817,
           'C26H42O4':0.01020, 'C30H62':0.00203,
           },
        aqueous_composition={'Water':1},
        internal_heat_exchanging=False,
        use_decorated_cost='Hydrocracker',
        tau=15/60, # set to the same as HTL
        V_wf=0.4, # Towler
        length_to_diameter=2, diameter=None,
        N=None, V=None, auxiliary=False,
        mixing_intensity=None, kW_per_m3=0,
        wall_thickness_factor=1.5,
        vessel_material='Stainless steel 316',
        vessel_type='Vertical',
        )
    HC.register_alias('Hydrocracking')
    # In [1], HC is costed for a multi-stage HC, but commented that the cost could be
    # $10-70 MM (originally $25 MM for a 6500 bpd system),
    # HC.cost_items['Hydrocracker'].cost = 10e6
    
    HC_HX = qsu.HXutility(
        'HC_HX', ins=HC-0, outs='cooled_HC_eff', T=60+273.15,
        init_with='Stream', rigorous=True)
    
    # To depressurize products
    HC_IV = IsenthalpicValve('HC_IV', ins=HC_HX-0, outs='cooled_depressed_HC_eff', P=30*6894.76, vle=True)
    
    # To separate products, can be adjusted to minimize fuel chemicals in the gas phase
    HCflash = qsu.Flash('HCflash', ins=HC_IV-0, outs=('HC_fuel_gas','HC_liquid'),
                        T=60.2+273.15, P=30*_psi_to_Pa,)
    
    HCpump = qsu.Pump('HCpump', ins=HCflash-1, init_with='Stream')
    
    # Separate water from oil
    HCliquidSplitter = qsu.Splitter('HCliquidSplitter', ins=HCpump-0,
                                    outs=('HC_ww','HC_oil'),
                                    split={'H2O':1}, init_with='Stream')
    class DummyHX:
          ins = [None]
          outs = [None]

    HC.eff_hx = DummyHX()
    

    HC.eff_hx = None
    # =========================================================================
    # Hydrotreating
    # =========================================================================
    
    # Pd/Al2O3
    HTcatalyst_in = qs.WasteStream('HTcatalyst_in', HTcatalyst=1, price=price_dct['HTcatalyst'])
    
    # Light (gasoline, <C8): medium (jet, C8-C14): heavy (diesel, >C14)
    oil_fracs = [0.3455, 0.4479, 0.2066]
    HT = u.Hydroprocessing(
        'HT',
        ins=(HCliquidSplitter-1, 'H2_HT', HTcatalyst_in),
        outs=('HTout','HTcatalyst_out'),
        WHSV=0.625,
        catalyst_lifetime=2*uptime_ratio*365*24, # 2 years [1]
        catalyst_ID='HTcatalyst',
        T=300+273.15,
        P=1500*_psi_to_Pa,
        hydrogen_rxned_to_inf_oil=0.0207, # data from expt
        hydrogen_ratio=3,
        include_PSA=include_PSA,
        gas_yield=0.2143,
        oil_yield=0.8637,
        gas_composition={'CO2':0.03880, 'CH4':0.00630,}, # [1] after the second hydroprocessing
        oil_composition={
            'Gasoline': oil_fracs[0],
            'Jet': oil_fracs[1],
            'Diesel': oil_fracs[2],
            },
        aqueous_composition={'Water':1},
        internal_heat_exchanging=False,
        use_decorated_cost='Hydrotreater',
        tau=0.5, V_wf=0.4, # Towler
        length_to_diameter=2, diameter=None,
        N=None, V=None, auxiliary=False,
        mixing_intensity=None, kW_per_m3=0,
        wall_thickness_factor=1,
        vessel_material='Stainless steel 316',
        vessel_type='Vertical',
        )
    HT.register_alias('Hydrotreating')
    
    HT_HX = qsu.HXutility('HT_HX',ins=HT-0, outs='cooled_HT_eff', T=60+273.15,
                          init_with='Stream', rigorous=True)
    
    HT_IV = IsenthalpicValve('HT_IV', ins=HT_HX-0, outs='cooled_depressed_HT_eff',
                             P=717.4*_psi_to_Pa, vle=True)
    
    # To separate products, can be adjusted to minimize fuel chemicals in the gas phase
    HTflash = qsu.Flash('HTflash', ins=HT_IV-0, outs=('HT_fuel_gas','HT_liquid'),
                        T=43+273.15, P=55*_psi_to_Pa)
    
    HTpump = qsu.Pump('HTpump', ins=HTflash-1, init_with='Stream')
    
    # Separate water from oil, if any
    HTliquidSplitter = qsu.Splitter('HTliquidSplitter', ins=HTpump-0,
                                    outs=('HT_ww','HT_oil'),
                                    split={'H2O':1}, init_with='Stream')
    
    # Separate gasoline from jet and diesel
    GasolineDis = qsu.ShortcutColumn(
        'OilLightDis', ins=HTliquidSplitter-1,
        outs=('hot_gasoline','jet_diesel'),
        LHK=('Gasoline', 'Jet'),
        Lr=0.99,
        Hr=0.99,
        k=2, is_divided=True)
    # Lr_range = Hr_range = np.linspace(0.05, 0.95, 19)
    # Lr_range = Hr_range = np.linspace(0.01, 0.2, 20)
    # results = find_Lr_Hr(GasolineDis, Lr_trial_range=Lr_range, Hr_trial_range=Hr_range, target_light_frac=oil_fracs[0])
    # results_df, Lr, Hr = results
    
    GasolineFlash = qsu.Flash('GasolineFlash', ins=GasolineDis-0, outs=('', 'cooled_gasoline',),
                              T=298.15, P=101325)
    
    # Separate jet from diesel
    JetDis = qsu.ShortcutColumn(
        'JetDis', ins=GasolineDis-1,
        outs=('hot_jet','hot_diesel'),
        LHK=('Jet', 'Diesel'),
        Lr=0.99,
        Hr=0.99,
        k=2, is_divided=True)
    # Lr_range = Hr_range = np.linspace(0.05, 0.95, 19)
    # Lr_range = Hr_range = np.linspace(0.01, 0.2, 20)
    # results = find_Lr_Hr(JetDis, Lr_trial_range=Lr_range, Hr_trial_range=Hr_range, target_light_frac=oil_fracs[1]/(1-oil_fracs[0]))
    # results_df, Lr, Hr = results
    
    JetFlash = qsu.Flash('JetFlash', ins=JetDis-0, outs=('', 'cooled_jet',), T=298.15, P=101325)
    
    DieselHX = qsu.HXutility('DieselHX',ins=JetDis-1, outs='cooled_diesel', T=298.15,
                             init_with='Stream', rigorous=True)
    
    HT.eff_hx = DummyHX()
    
    # =========================================================================
    # Products and Wastes
    # =========================================================================
    
    GasolinePC = qsu.PhaseChanger('GasolinePC', ins=GasolineFlash-1)
    gasoline = qs.WasteStream('gasoline', Gasoline=1)
    # gasoline.price = price_dct['gasoline']/(gasoline.rho/_m3_to_gal)
    # Storage time assumed to be 3 days per [1]
    GasolineTank = qsu.StorageTank('GasolineTank', ins=GasolinePC-0, outs=(gasoline),
                                    tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')
    
    JetPC = qsu.PhaseChanger('JetPC', ins=JetFlash-1)
    jet = qs.WasteStream('jet', Jet=1)
    # jet.price = price_dct['jet']/(jet.rho/_m3_to_gal)
    JetTank = qsu.StorageTank('JetTank', ins=JetPC-0, outs=(jet,),
                              tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')
    
    DieselPC = qsu.PhaseChanger('DieselPC', ins=DieselHX-0)
    diesel = qs.WasteStream('diesel', Jet=1)
    # diesel.price = price_dct['diesel']/(diesel.rho/_m3_to_gal)
    DieselTank = qsu.StorageTank('DieselTank', ins=DieselPC-0, outs=(diesel,),
                                 tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')
    
    # Combine all fuel to get a one fuel selling price
    mixed_fuel = qs.WasteStream('mixed_fuel')
    FuelMixer = qsu.Mixer('FuelMixer', ins=(GasolineTank-0, JetTank-0, DieselTank-0), outs=mixed_fuel)
    
    # =========================================================================
    # Electrochemical Unit
    # =========================================================================
    # All wastewater streams
    ww_streams = [HTLaqMixer-0, HCliquidSplitter-0, HTliquidSplitter-0]
    # Wastewater sent to municipal wastewater treatment plant
    ww_to_disposal = qs.WasteStream('ww_to_disposal')

    WWmixer = qsu.Mixer('WWmixer', ins=ww_streams)
    
    fuel_gases = [
        HTL-0, CrudeLightFlash-0, # HTL gases
        HCflash-0, HTflash-0, # post-hydroprocessing gases
        GasolineFlash-0, JetFlash-0, # final distillation fuel gases
        ]
       
    recovered_N = qs.WasteStream('recovered_N', price=price_dct['N'])
    recovered_P = qs.WasteStream('recovered_P', price=price_dct['P'])
    recovered_K = qs.WasteStream('recovered_K', price=price_dct['K'])

    EC = u.SAFElectrochemical(
        'EC',
        ins=(WWmixer-0, 'EC_replacement_surrogate'),
        outs=('EC_gas', 'EC_H2', recovered_N, recovered_P, recovered_K, ww_to_disposal),
        COD_removal=0.95, # assumed
        N_recovery=0.8,
        P_recovery=0.99,
        K_recovery=0.8,
        include_PSA=include_PSA,
        PSA_efficiency=0.95,
        )
    EC.register_alias('Electrochemical')
    EC.include_PSA_cost = False # HC/HT has PSA
    fuel_gases.append(EC-0)
    if include_EC is False:
        EC.skip = True
    else:
        EC.skip = False
        if type(include_EC) is dict:
            for attr, val in include_EC.items():
                setattr(EC, attr, val)

    def adjust_prices():
        # Transportation
        dry_price = price_dct['trans_feedstock']
        factor = 1 - FeedstockTrans.ins[0].imass['Water']/FeedstockTrans.ins[0].F_mass
        FeedstockTrans.transportation_unit_cost = dry_price * factor
        # Wastewater
        ww_to_disposal.source._run()
        COD_mass_content = ww_to_disposal.COD*ww_to_disposal.F_vol/1e3 # mg/L*m3/hr to kg/hr
        factor = COD_mass_content/ww_to_disposal.F_mass
        ww_to_disposal.price = price_dct['COD']*factor
        ww_to_disposal_item = qs.ImpactItem.get_item('ww_to_disposal_item')
        try: ww_to_disposal_item.CFs['GWP'] = gwp_dct['COD']*factor
        except: pass
    ww_to_disposal.source.add_specification(adjust_prices)

    GasMixer = qsu.Mixer('GasMixer', ins=fuel_gases, outs=('waste_gases'))

    # =========================================================================
    # Facilities
    # =========================================================================
    
    # Adding HXN only saves cents/GGE with HTL internal HX, eliminate for simpler system
    # HXN = qsu.HeatExchangerNetwork('HXN', T_min_app=86, force_ideal_thermo=True)
    # 86 K: Jones et al. PNNL, 2014
    
    natural_gas = qs.WasteStream('natural_gas', CH4=1, price=price_dct['natural_gas'])
    solids_to_disposal = qs.WasteStream('solids_to_disposal', price=price_dct['solids'])
    CHPmixer = qsu.Mixer('CHPmixer', ins=(GasMixer-0, CrudeHeavyDis-1, HTL-3))
    CHP = qsu.CombinedHeatPower('CHP', 
                                ins=(CHPmixer-0, natural_gas, 'air'),
                                outs=('gas_emissions', solids_to_disposal),
                                init_with='WasteStream',
                                supplement_power_utility=False)
    
    H2C = u.HydrogenCenter(
        'H2C',
        process_H2_streams=(HC.ins[1], HT.ins[1]),
        recycled_H2_streams=EC-1,
        )
    H2C.register_alias('HydrogenCenter')
    H2C.makeup_H2_price = H2C.excess_H2_price = price_dct['H2'] # expected H2 price
    # H2C.makeup_H2_price = H2C.excess_H2_price = 33.4 # current H2 price
    
    PWC = u.ProcessWaterCenter('PWC', process_water_streams=[feedstock_water],)
    PWC.register_alias('ProcessWaterCenter')
    PWC.process_water_price = price_dct['process_water']
    
    # =========================================================================
    # System, TEA, LCA
    # =========================================================================   
    sys = qs.System.from_units(
        'sys',
        units=list(flowsheet.unit),
        operating_hours=365*24*uptime_ratio,
        )
    for unit in sys.units: unit.include_construction = False
    
    tea = create_tea(sys, **tea_kwargs)  
    

    # Add characterization factors for each impact item
    clear_lca_registries()
    GWP = qs.ImpactIndicator('GWP',
                             alias='GlobalWarmingPotential',
                             method='GREET',
                             category='environmental impact',
                             unit='kg CO2-eq',)
    feedstock_item = qs.StreamImpactItem(
        ID='feedstock_item',
        linked_stream=feedstock,
        # feedstock, landfill, composting, anaerobic_digestion
        # may or may not be good to assume landfill offsetting
        GWP=-gwp_dct['landfill'],
        )
    trans_feedstock_item = qs.StreamImpactItem(
        ID='trans_feedstock_item',
        linked_stream=FeedstockTrans.ins[1],
        GWP=gwp_dct['trans_feedstock'],
        )
    makeup_H2_item = qs.StreamImpactItem(
        ID='makeup_H2_item',
        linked_stream=H2C.ins[0],
        GWP=gwp_dct['H2'],
        )
    excess_H2_item = qs.StreamImpactItem(
        ID='excess_H2_item',
        linked_stream=H2C.outs[1],
        GWP=-gwp_dct['H2'],
        )
    HCcatalyst_item = qs.StreamImpactItem(
        ID='HCcatalyst_item',
        linked_stream=HC.ins[-1],
        GWP=gwp_dct['HCcatalyst'],
        )
    HTcatalyst_item = qs.StreamImpactItem(
        ID='HTcatalyst_item',
        linked_stream=HT.ins[-1],
        GWP=gwp_dct['HTcatalyst'],
        )
    natural_gas_item = qs.StreamImpactItem(
        ID='natural_gas_item',
        linked_stream=natural_gas,
        GWP=gwp_dct['natural_gas'],
        )
    # Assume no impacts from process water
    # process_water_item = qs.StreamImpactItem(
    #     ID='process_water_item',
    #     linked_stream=PWC.ins[-1],
    #     GWP=gwp_dct['process_water'],
    #     )
    ww_to_disposal_item = qs.StreamImpactItem(
        ID='ww_to_disposal_item',
        linked_stream=ww_to_disposal,
        GWP=gwp_dct['COD'], # will be updated based on COD content
        )
    solids_to_disposal_item = qs.StreamImpactItem(
        ID='solids_to_disposal_item',
        linked_stream=CHP.outs[1],
        GWP=gwp_dct['solids'],
        )
    qs.PowerUtility.price = electricitry_price if electricitry_price is not None else qs.PowerUtility.price
    e_item = qs.ImpactItem(
        ID='e_item',
        GWP=electricitry_GHG if electricitry_GHG is not None else gwp_dct['electricity'],
        )
    steam_item = qs.ImpactItem(
        ID='steam_item',
        GWP=gwp_dct['steam'],
        )
    cooling_item = qs.ImpactItem(
        ID='cooling_item',
        GWP=gwp_dct['cooling'],
        )
    recovered_N_item = qs.StreamImpactItem(
        ID='recovered_N_item',
        linked_stream=recovered_N,
        GWP=gwp_dct['N'],
        )
    recovered_P_item = qs.StreamImpactItem(
        ID='recovered_P_item',
        linked_stream=recovered_P,
        GWP=gwp_dct['P'],
        )
    recovered_K_item = qs.StreamImpactItem(
        ID='recovered_K_item',
        linked_stream=recovered_K,
        GWP=gwp_dct['K'],
        )

    lifetime = tea.duration[1]-tea.duration[0]
    lca = qs.LCA(
        system=sys,
        lifetime=lifetime,
        uptime_ratio=uptime_ratio,
        simulate_system=False,
        e_item=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
        steam_item=lambda:sys.get_heating_duty()/1000*lifetime, # kJ/yr to MJ/yr, include natural gas, but all offset in CHP
        cooling_item=lambda:sys.get_cooling_duty()/1000*lifetime, # kJ/yr to MJ/yr
        )
    
    return sys

# %%

# =========================================================================
# Result outputting
# =========================================================================
def get_GGE(sys, fuel, annual=True):
    """Gasoline gallon equivalent for a fuel stream (GGE/hr if annual=False)."""
    # HHV in kJ/hr -> MJ/hr, then / (MJ/GGE)
    gge_per_hr = fuel.HHV / 1e3 / _HHV_per_GGE
    return gge_per_hr * (sys.operating_hours if annual else 1)


def fuel_blend_ratios(sys, basis="mass"):
    """
    Returns (gasoline_ratio, jet_ratio, diesel_ratio) into mixed_fuel.
    basis: "mass" or "gge"
    """
    s = sys.flowsheet.stream
    gasoline, jet, diesel = s.gasoline, s.jet, s.diesel

    if basis.lower() == "mass":
        g = float(gasoline.F_mass)
        j = float(jet.F_mass)
        d = float(diesel.F_mass)
        tot = g + j + d
        if tot <= 0:
            return (np.nan, np.nan, np.nan)
        return (g / tot, j / tot, d / tot)

    if basis.lower() == "gge":
        gG = float(gasoline.HHV / 1e3 / _HHV_per_GGE)
        jG = float(jet.HHV / 1e3 / _HHV_per_GGE)
        dG = float(diesel.HHV / 1e3 / _HHV_per_GGE)
        tot = gG + jG + dG
        if tot <= 0:
            return (np.nan, np.nan, np.nan)
        return (gG / tot, jG / tot, dG / tot)

    raise ValueError("basis must be 'mass' or 'gge'")
def get_MFSP(sys, print_msg=False):
    mixed_fuel = sys.flowsheet.stream.mixed_fuel
    mixed_fuel.price = sys.TEA.solve_price(mixed_fuel)
    MFSP = mixed_fuel.cost / get_GGE(sys, mixed_fuel, annual=False)
    if print_msg:
        print(f"Minimum selling price of all fuel is ${MFSP:.2f}/GGE.")
    return MFSP


def simulate_and_print(sys, save_report=False):
    sys.simulate()

    stream = sys.flowsheet.stream
    tea = sys.TEA

    fuels = (gasoline, jet, diesel) = (stream.gasoline, stream.jet, stream.diesel)

    # ADDED: blend ratios into mixed_fuel
    rm = fuel_blend_ratios(sys, basis="mass")
    re = fuel_blend_ratios(sys, basis="gge")
    print("\nMixed fuel blend ratios into 'mixed_fuel'")
    print("----------------------------------------")
    print(f"Mass basis: gasoline={rm[0]:.3f}, jet={rm[1]:.3f}, diesel={rm[2]:.3f}")
    print(f"GGE basis : gasoline={re[0]:.3f}, jet={re[1]:.3f}, diesel={re[2]:.3f}")

    MFSP = get_MFSP(sys, print_msg=True)

    table = tea.get_cashflow_table()
    c = qs.currency
    for attr in ("NPV", "AOC", "sales", "net_earnings"):
        uom = c if attr in ("NPV", "CAPEX") else (c + "/yr")
        print(f"{attr} is {getattr(tea, attr):,.0f} {uom}")

    if save_report:
        sys.save_report(file=os.path.join(os.getcwd(), f"sys_{sys.flowsheet.ID}.xlsx"))

    return MFSP


# =========================
# Run
# =========================

if __name__ == "__main__":
    config = dict(include_PSA=True, include_EC=False)

    sys = create_system(flowsheet=None, **config)
    dct = globals()
    dct.update(sys.flowsheet.to_dict())

    simulate_and_print(sys)