
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>
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

# -*- coding: utf-8 -*-
# !!! Temporarily ignoring warnings
import warnings
warnings.filterwarnings('ignore')

import os
import numpy as np
import biosteam as bst
import qsdsan as qs
from biosteam import IsenthalpicValve
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import create_tea

# ============================
# UPSTREAM: use biobinder_ml
# ============================
from exposan.saf import (
    _HHV_per_GGE,
    _load_components,
    _load_process_settings,
    _units as u_saf,
    data_path,
    dry_flowrate,
    find_Lr_Hr,
    get_mass_energy_balance,
    gwp_dct,
    # HTL_yields,
    price_dct,
    results_path,
    tea_kwargs,
    uptime_ratio,
)

from exposan.biobinder_ml import (
    central_dry_flowrate as default_central,
    pilot_dry_flowrate as default_pilot,
    feedprice_dct as feedprice_dct,
    HTL_yields,
    limit_internal_hx,
    configure_utilities,
    set_eff_hx_temperature,
    BiobinderTEA,
)
from exposan.biobinder_ml._feedstocks import get_feedstock_composition
from exposan.biobinder_ml import _units as u_bb
from exposan.biobinder_ml.Dist_flex import gwp_dict
# ============================
# DOWNSTREAM: SAF helpers
# ============================
from exposan.saf import price_dct, gwp_dct  # keep SAF CFs if you want SAF GWP dictionary

from exposan.biobinder_ml._distill_utils import add_runtime_LHK_LrHr_spec

_psi_to_Pa = 6894.76
_m3_to_gal = 264.172
SAF_ECON_MODE = "market"  # "market" or "mfsp"
configure_utilities()

__all__ = ("create_system", "get_MFSP", "simulate_and_print", "fuel_blend_ratios")

class _DummyEffHX:
    """
    Minimal object that satisfies exposan.saf._units.Hydroprocessing._design()
    as 'eff_hx'. It prevents thermo/chemical mismatch during mix_from(self.outs).
    """
    def __init__(self, ID, thermo):
        self.ID = ID
        # IMPORTANT: both ins and outs must be real streams (not None)
        self.ins = [bst.Stream(f"{ID}_in", thermo=thermo)]
        self.outs = [bst.Stream(f"{ID}_out", thermo=thermo)]
        # Hydroprocessing._design uses eff_hx.Hnet
        self.Hnet = 0.0

    def simulate_as_auxiliary_exchanger(self, ins, outs, duty=None):
        self._last_duty = duty
        return

def create_system(
    flowsheet=None,
    feedstock_id="fog",

    # biobinder-style knobs
    central_dry_flowrate=None,
    pilot_dry_flowrate=None,
    decentralized_HTL=False,
    decentralized_upgrading=False,

    # SAF toggles
    include_PSA=True,
    include_EC=False,

    # electricity overrides
    electricity_price=None,
    electricity_GHG=None,

    # optional override of HTL yields (e.g., RF predicted)
    HTL_yields_override=None,

    # runtime distillation tuning
    prefer="max_top_ratio",
    tol=0.05,

    # catch-all so UA doesn’t crash
    **kwargs,
):
    # ---------------------------------------------------------------------
    # 0) Clear registries like biobinder does (prevents stale items)
    # ---------------------------------------------------------------------
    qs.main_flowsheet.clear()
    clear_lca_registries()

    if central_dry_flowrate is None:
        central_dry_flowrate = default_central
    if pilot_dry_flowrate is None:
        pilot_dry_flowrate = default_pilot

    # ---------------------------------------------------------------------
    # 1) Match biobinder CHCU/DHCU/DHDU logic exactly
    # ---------------------------------------------------------------------
    if decentralized_HTL is False:
        if decentralized_upgrading is False:
            flowsheet_ID = "saf_CHCU"
            N_HTL = N_upgrading = 1
        else:
            raise ValueError("Centralized HTL + decentralized upgrading is not a valid configuration.")
    else:
        if decentralized_upgrading is False:
            flowsheet_ID = "saf_DHCU"
            N_HTL = round(central_dry_flowrate / pilot_dry_flowrate)
            N_upgrading = 1
            pilot_dry_flowrate = central_dry_flowrate / N_HTL
        else:
            flowsheet_ID = "saf_DHDU"
            N_HTL = N_upgrading = 1
            central_dry_flowrate = pilot_dry_flowrate

    # ---------------------------------------------------------------------
    # 2) Flowsheet init
    # ---------------------------------------------------------------------
    if flowsheet is None:
        flowsheet = qs.Flowsheet(flowsheet_ID)
        qs.main_flowsheet.set_flowsheet(flowsheet)
    else:
        qs.main_flowsheet.set_flowsheet(flowsheet)

    _load_process_settings()
    _load_components()

    # -----------------------------------------
    # 3) Feedstock scaling + pricing
    # -----------------------------------------
    scaled_feedstock = qs.WasteStream("scaled_feedstock", price=0.0)
    try:
        scaled_feedstock.price = float(feedprice_dct[feedstock_id])
    except KeyError:
        raise KeyError(f"No feedstock price defined for '{feedstock_id}' in feedprice_dct")

    FeedstockScaler = u_bb.Scaler(
        "FeedstockScaler", ins=scaled_feedstock, outs="feedstock",
        scaling_factor=N_HTL, reverse=True,
    )

    FeedstockTrans = u_bb.Transportation(
        "FeedstockTrans",
        ins=(FeedstockScaler - 0, "feedstock_trans_surrogate"),
        outs=("transported_feedstock",),
        N_unit=N_HTL,
        copy_ins_from_outs=True,
        transportation_unit_cost=0,   # set later if you want
        transportation_distance=1,
    )

    # Process water (biobinder style)
    scaled_process_water = qs.WasteStream("scaled_process_water")
    ProcessWaterScaler = u_bb.Scaler(
        "ProcessWaterScaler", ins=scaled_process_water, outs="htl_process_water",
        scaling_factor=N_HTL, reverse=True,
    )

    feedstock_composition = get_feedstock_composition(feedstock_id)
    print(f"Feedstock used is {feedstock_id}")
    print("Corresponding composition:")
    for k, v in feedstock_composition.items():
        print(f"  {k}: {v:.4f}")

    FeedstockCond = u_bb.Conditioning(
        "FeedstockCond",
        ins=(FeedstockTrans - 0, ProcessWaterScaler.outs[0]),
        outs="conditioned_feedstock",
        feedstock_composition=feedstock_composition,
        feedstock_dry_flowrate=central_dry_flowrate if decentralized_HTL is False else pilot_dry_flowrate,
        target_HTL_solid_loading=0.2,
    )
    FeedstockCond.N_unit = N_HTL

    # ---------------------------------------------------------------------
    # 4) HTL (biobinder CentralizedHTL / PilotHTL)
    # ---------------------------------------------------------------------
    aqueous_composition = {"N": 0.48 / 100}
    aqueous_composition["HTLaqueous"] = 1 - sum(aqueous_composition.values())

    dw_yields = HTL_yields_override if HTL_yields_override is not None else HTL_yields

    HTL_kwargs = dict(
        ID="HTL",
        ins=FeedstockCond.outs[0],
        outs=("gas", "HTL_aqueous", "biocrude", "hydrochar"),
        T=280 + 273.15,
        P=12.4e6,
        tau=15 / 60,
        dw_yields=dw_yields,
        gas_composition={"CO2": 1},
        aqueous_composition=aqueous_composition,
        biocrude_composition={"HTLbiocrude": 1},
        char_composition={"HTLchar": 1},
        internal_heat_exchanging=True,
        eff_T=60 + 273.15,
        eff_P=30 * _psi_to_Pa,
    )

    if decentralized_HTL is False:
        HTL_unit = u_bb.CentralizedHTL
    else:
        HTL_unit = u_bb.PilotHTL
        HTL_kwargs["N_unit"] = N_HTL

    HTL = HTL_unit(**HTL_kwargs)

    # scale gas/char to central totals (biobinder pattern)
    HTLgasScaler = u_bb.Scaler("HTLgasScaler", ins=HTL - 0, outs="scaled_HTLgas", scaling_factor=N_HTL, reverse=False)
    HTLcharScaler = u_bb.Scaler("HTLcharScaler", ins=HTL - 1 + 2, outs="scaled_HTLchar", scaling_factor=N_HTL, reverse=False) \
        if False else u_bb.Scaler("HTLcharScaler", ins=HTL.outs[-1], outs="scaled_HTLchar", scaling_factor=N_HTL, reverse=False)

    # internal HX guardrails (same helpers you use in biobinder)
    limit_internal_hx(HTL, coldest_id="chilled_water", cap_margin_K=2.0, cooler_margin_K=3.0)
    set_eff_hx_temperature(HTL, clamp_to="cooling_water")

    # ---------------------------------------------------------------------
    # 5) Biocrude transport + scaling to central (biobinder pattern)
    #     - This ensures downstream upgrading sees ONE stream representing the total system
    # ---------------------------------------------------------------------
    BiocrudeTrans = u_bb.Transportation(
        "BiocrudeTrans",
        ins=(HTL - 2, "biocrude_trans_surrogate"),
        outs=("transported_biocrude",),
        N_unit=N_HTL,
        transportation_unit_cost=0,
        transportation_distance=1,
    )

    BiocrudeScaler = u_bb.Scaler(
        "BiocrudeScaler",
        ins=BiocrudeTrans - 0,
        outs="scaled_biocrude",
        scaling_factor=N_HTL,
        reverse=False,
    )

    # =========================================================================
    # ================== DOWNSTREAM SAF SECTION STARTS HERE ====================
    # =========================================================================

    # --- crude split ---
    cutoff_fracs = [0.01, 0.69, 0.30]  # keep your SAF baseline split for now

    CrudeSplitter = u_bb.BiocrudeSplitter(
        "CrudeSplitter",
        ins=BiocrudeScaler - 0,
        outs="splitted_crude",
        biocrude_IDs=("HTLbiocrude",),
        cutoff_fracs=cutoff_fracs,
        cutoff_Tbs=(150 + 273.15, 300 + 273.15),
    )

    CrudePump = qsu.Pump("CrudePump", ins=CrudeSplitter - 0, outs="crude_to_dist", init_with="Stream")

    CrudeLightDis = qsu.ShortcutColumn(
        "CrudeLightDis",
        ins=CrudePump - 0,
        outs=("crude_light", "crude_medium_heavy"),
        LHK=CrudeSplitter.keys[0],
        P=50 * _psi_to_Pa,
        Lr=0.87,
        Hr=0.98,
        k=2,
        is_divided=True,
    )

    CrudeLightFlash = qsu.Flash(
        "CrudeLightFlash",
        ins=CrudeLightDis - 0,
        outs=("crude_light_gas", "crude_light_liq"),
        T=298.15,
        P=101325,
    )
    
    HTLaqMixer = qsu.Mixer('HTLaqMixer', ins=(HTL-1, CrudeLightFlash-1), outs='HTL_aq')

    # --- separate medium from char ---
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

    ratio0 = CrudeSplitter.cutoff_fracs[1] / sum(CrudeSplitter.cutoff_fracs[1:])
    add_runtime_LHK_LrHr_spec(
        splitter=CrudeSplitter,
        col=CrudeHeavyDis,
        idx=1,
        fallback_LHK=("4M-PHYNO", "INDOLE"),
        target_ratio=ratio0,
        tol=float(tol),
        prefer=str(prefer),
        require_design=True,
        require_cost=False,
        vle_screen=True,
        design_sanity=True,
    )

    # =========================================================================
    # Hydrocracking (keep SAF downstream logic)
    # =========================================================================
    HCcatalyst_in = qs.WasteStream("HCcatalyst_in", HCcatalyst=1, price=price_dct["HCcatalyst"])

    HC = u_saf.Hydroprocessing(
        "HC",
        ins=(CrudeHeavyDis - 0, "H2_HC", HCcatalyst_in),
        outs=("HC_out", "HCcatalyst_out"),
        T=400 + 273.15,
        P=1500 * _psi_to_Pa,
        WHSV=0.625,
        catalyst_ID="HCcatalyst",
        catalyst_lifetime=5 * uptime_ratio * 365 * 24,
        hydrogen_rxned_to_inf_oil=0.0111,
        hydrogen_ratio=5.556,
        include_PSA=include_PSA,
        gas_yield=0.2665,
        oil_yield=0.7335,
        gas_composition={
            "CO2": 1 - 0.08809,
            "CH4": 0.02280, "C2H6": 0.02923,
            "C3H8": 0.01650, "C4H10": 0.00870,
            "TWOMBUTAN": 0.00408, "NPENTAN": 0.00678,
        },
        oil_composition={
            'TWOMPENTA': 0.00408, 'HEXANE': 0.00408,
            'TWOMHEXAN': 0.00408, 'HEPTANE': 0.00408,
            'CC6METH': 0.01020, 'PIPERDIN': 0.00408,
            'TOLUENE': 0.01020, 'THREEMHEPTA': 0.01020,
            'OCTANE': 0.01020, 'ETHCYC6': 0.00408,
            'ETHYLBEN': 0.02040, 'OXYLENE': 0.01020,
            'C9H20': 0.00408, 'PROCYC6': 0.00408,
            'C3BENZ': 0.01020, 'FOURMONAN': 0,
            'C10H22': 0.00203, 'C4BENZ': 0.01223,
            'C11H24': 0.02040, 'C10H12': 0.02040,
            'C12H26': 0.02040, 'OTTFNA': 0.01020,
            'C6BENZ': 0.02040, 'OTTFSN': 0.02040,
            'C7BENZ': 0.02040, 'C8BENZ': 0.02040,
            'C10H16O4': 0.01837, 'C15H32': 0.06120,
            'C16H34': 0.18360, 'C17H36': 0.08160,
            'C18H38': 0.04080, 'C19H40': 0.04080,
            'C20H42': 0.10200, 'C21H44': 0.04080,
            'TRICOSANE': 0.04080, 'C24H38O4': 0.00817,
            'C26H42O4': 0.01020, 'C30H62': 0.00203,
        },
        aqueous_composition={"Water": 1},
        internal_heat_exchanging=False,
        use_decorated_cost="Hydrocracker",
        tau=15/60,
        V_wf=0.4,
        length_to_diameter=2,
        vessel_material="Stainless steel 316",
        vessel_type="Vertical",
    )

    HC_HX = qsu.HXutility("HC_HX", ins=HC - 0, outs="cooled_HC_eff", T=60 + 273.15, init_with="Stream", rigorous=True)
    HC.eff_hx = _DummyEffHX("HC_eff_hx", bst.settings.thermo)
    HC_IV = IsenthalpicValve("HC_IV", ins=HC_HX - 0, outs="cooled_depressed_HC_eff", P=30 * _psi_to_Pa, vle=True)
    HCflash = qsu.Flash("HCflash", ins=HC_IV - 0, outs=("HC_fuel_gas", "HC_liquid"), T=60.2 + 273.15, P=30 * _psi_to_Pa)
    HCpump = qsu.Pump("HCpump", ins=HCflash - 1, init_with="Stream")

    HCliquidSplitter = qsu.Splitter(
        "HCliquidSplitter",
        ins=HCpump - 0,
        outs=("HC_ww", "HC_oil"),
        split={"H2O": 1},
        init_with="Stream",
    )

    # =========================================================================
    # Hydrotreating
    # =========================================================================
    HTcatalyst_in = qs.WasteStream("HTcatalyst_in", HTcatalyst=1, price=price_dct["HTcatalyst"])
    oil_fracs = [0.3455, 0.4479, 0.2066]

    HT = u_saf.Hydroprocessing(
        "HT",
        ins=(HCliquidSplitter - 1, "H2_HT", HTcatalyst_in),
        outs=("HTout", "HTcatalyst_out"),
        WHSV=0.625,
        catalyst_lifetime=2 * uptime_ratio * 365 * 24,
        catalyst_ID="HTcatalyst",
        T=300 + 273.15,
        P=1500 * _psi_to_Pa,
        hydrogen_rxned_to_inf_oil=0.0207,
        hydrogen_ratio=3,
        include_PSA=include_PSA,
        gas_yield=0.2143,
        oil_yield=0.8637,
        gas_composition={"CO2": 0.03880, "CH4": 0.00630},
        oil_composition={"Gasoline": oil_fracs[0], "Jet": oil_fracs[1], "Diesel": oil_fracs[2]},
        aqueous_composition={"Water": 1},
        internal_heat_exchanging=False,
        use_decorated_cost="Hydrotreater",
        tau=0.5,
        V_wf=0.4,
        length_to_diameter=2,
        vessel_material="Stainless steel 316",
        vessel_type="Vertical",
    )

    HT_HX = qsu.HXutility("HT_HX", ins=HT - 0, outs="cooled_HT_eff", T=60 + 273.15, init_with="Stream", rigorous=True)
    HT.eff_hx = _DummyEffHX("HT_eff_hx", bst.settings.thermo)
    HT_IV = IsenthalpicValve("HT_IV", ins=HT_HX - 0, outs="cooled_depressed_HT_eff", P=717.4 * _psi_to_Pa, vle=True)
    HTflash = qsu.Flash("HTflash", ins=HT_IV - 0, outs=("HT_fuel_gas", "HT_liquid"), T=43 + 273.15, P=55 * _psi_to_Pa)
    HTpump = qsu.Pump("HTpump", ins=HTflash - 1, init_with="Stream")

    HTliquidSplitter = qsu.Splitter(
        "HTliquidSplitter",
        ins=HTpump - 0,
        outs=("HT_ww", "HT_oil"),
        split={"H2O": 1},
        init_with="Stream",
    )

    GasolineDis = qsu.ShortcutColumn(
        "OilLightDis",
        ins=HTliquidSplitter - 1,
        outs=("hot_gasoline", "jet_diesel"),
        LHK=("Gasoline", "Jet"),
        Lr=0.99,
        Hr=0.99,
        k=2,
        is_divided=True,
    )
    GasolineFlash = qsu.Flash("GasolineFlash", ins=GasolineDis - 0, outs=("", "cooled_gasoline"), T=298.15, P=101325)

    JetDis = qsu.ShortcutColumn(
        "JetDis",
        ins=GasolineDis - 1,
        outs=("hot_jet", "hot_diesel"),
        LHK=("Jet", "Diesel"),
        Lr=0.99,
        Hr=0.99,
        k=2,
        is_divided=True,
    )
    JetFlash = qsu.Flash("JetFlash", ins=JetDis - 0, outs=("", "cooled_jet"), T=298.15, P=101325)
    DieselHX = qsu.HXutility("DieselHX", ins=JetDis - 1, outs="cooled_diesel", T=298.15, init_with="Stream", rigorous=True)

    # =========================================================================
    # Products
    # =========================================================================
    def _set_price_from_gal(stream, price_per_gal):
        """Convert $/gal to $/kg using the stream density."""
        rho = getattr(stream, "rho", None)
        if rho is None or (not np.isfinite(rho)) or rho <= 0:
            fallback_rho = {
            "gasoline": 740.0,   # kg/m3
            "jet": 800.0,
            "diesel": 832.0,
             }
            rho = fallback_rho.get(stream.ID, 800.0)
        rho_kg_per_gal = rho / 264.172
        stream.price = float(price_per_gal) / float(rho_kg_per_gal)

        if rho_kg_per_gal <= 0:
           raise ValueError(f"Non-positive density for stream {stream.ID}")
        stream.price = float(price_per_gal) / float(rho_kg_per_gal)
    GasolinePC = qsu.PhaseChanger("GasolinePC", ins=GasolineFlash - 1)
    gasoline = qs.WasteStream("gasoline", Gasoline=1)
    GasolineTank = qsu.StorageTank("GasolineTank", ins=GasolinePC - 0, outs=(gasoline,),
                                   tau=3*24, init_with="WasteStream", vessel_material="Carbon steel")

    JetPC = qsu.PhaseChanger("JetPC", ins=JetFlash - 1)
    jet = qs.WasteStream("jet", Jet=1)
    JetTank = qsu.StorageTank("JetTank", ins=JetPC - 0, outs=(jet,),
                              tau=3*24, init_with="WasteStream", vessel_material="Carbon steel")

    DieselPC = qsu.PhaseChanger("DieselPC", ins=DieselHX - 0)
    diesel = qs.WasteStream("diesel", Diesel=1)
    DieselTank = qsu.StorageTank("DieselTank", ins=DieselPC - 0, outs=(diesel,),
                                 tau=3*24, init_with="WasteStream", vessel_material="Carbon steel")

    mixed_fuel = qs.WasteStream("mixed_fuel", price=0.0)
    GasolineTank._run()
    JetTank._run()
    DieselTank._run()
    def _update_product_prices_and_mixed_fuel():        
      
        # MARKET / CUTS MODE: gasoline, jet, diesel are sold products
        if SAF_ECON_MODE == "market":
           _set_price_from_gal(gasoline, price_dct["gasoline"])
           _set_price_from_gal(jet,      price_dct["jet"])
           _set_price_from_gal(diesel,   price_dct["diesel"])
           mixed_fuel.price = 0.0

        # MFSP MODE: cuts are bookkeeping only, mixed fuel is the economic product
        elif SAF_ECON_MODE == "mfsp":
           gasoline.price = 0.0
           jet.price = 0.0
           diesel.price = 0.0
           mixed_fuel.price = 0.0

        else:
           raise ValueError(f"Unknown SAF_ECON_MODE: {SAF_ECON_MODE}")    

        # Always keep a bookkeeping mixed-fuel stream for reporting
        mixed_fuel.mix_from((gasoline, jet, diesel))

    # attach to last product unit so it updates each simulation
    DieselTank.add_specification(_update_product_prices_and_mixed_fuel)
    DieselTank.run_after_specifications = True
     

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

    EC = u_saf.SAFElectrochemical(
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
    
    H2C = u_saf.HydrogenCenter(
        'H2C',
        process_H2_streams=(HC.ins[1], HT.ins[1]),
        recycled_H2_streams=EC-1,
        )
    H2C.register_alias('HydrogenCenter')
    H2C.makeup_H2_price = H2C.excess_H2_price = price_dct['H2'] # expected H2 price
    # H2C.makeup_H2_price = H2C.excess_H2_price = 33.4 # current H2 price
    PWC = u_bb.ProcessWaterCenter(
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
    
    
    # =========================================================================
    # System, TEA, LCA
    # =========================================================================   
    sys = qs.System.from_units(
        'sys',
        units=list(flowsheet.unit),
        operating_hours=365*24*uptime_ratio,
        )
    for unit in sys.units: unit.include_construction = False
    tea = create_tea(sys, cls=BiobinderTEA, **tea_kwargs)
    land_factor = 115000/24 #115000 gal/day, PNNL 32731
    tea.land = lambda: 90000+(tea.system.flowsheet.unit.BiocrudeTrans.ins[0].F_vol*264.172)/land_factor*4.692*1e6 # 6% of 78.2M$, PNNL 32731
    

    # Add characterization factors for each impact item
    clear_lca_registries()
    GWP = qs.ImpactIndicator('GWP',
                             alias='GlobalWarmingPotential',
                             method='GREET',
                             category='environmental impact',
                             unit='kg CO2-eq',)

    feedstock_item = qs.StreamImpactItem(
        ID='feedstock_item',
        linked_stream=scaled_feedstock,
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
    trans_biocrude_item=qs.StreamImpactItem(
        ID='biocrude_trans_surrogate_item',
        linked_stream=BiocrudeTrans.ins[1],
        GWP=gwp_dict['trans_biocrude'],
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
    qs.PowerUtility.price = electricity_price if electricity_price is not None else qs.PowerUtility.price
    e_item = qs.ImpactItem(
        ID='e_item',
        GWP=electricity_GHG if electricity_GHG is not None else gwp_dct['electricity'],
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
        trans_biocrude_item.linked_stream = fake_stream             # No biocrude transportation

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


# ============================
# Result outputting
# ============================

# Gasoline gallon equivalent
get_GGE = lambda sys, fuel, annual=True: fuel.HHV/1e3/_HHV_per_GGE*max(1, bool(annual)*sys.operating_hours)


def fuel_blend_ratios(sys, basis="mass"):
    s = sys.flowsheet.stream
    gasoline, jet, diesel = s.gasoline, s.jet, s.diesel

    if basis.lower() == "mass":
        g, j, d = float(gasoline.F_mass), float(jet.F_mass), float(diesel.F_mass)
        tot = g + j + d
        return (np.nan, np.nan, np.nan) if tot <= 0 else (g / tot, j / tot, d / tot)

    if basis.lower() == "gge":
        gG = float(gasoline.HHV / 1e3 / _HHV_per_GGE)
        jG = float(jet.HHV / 1e3 / _HHV_per_GGE)
        dG = float(diesel.HHV / 1e3 / _HHV_per_GGE)
        tot = gG + jG + dG
        return (np.nan, np.nan, np.nan) if tot <= 0 else (gG / tot, jG / tot, dG / tot)

    raise ValueError("basis must be 'mass' or 'gge'")


def get_MFSP(sys, print_msg=False, as_reporter=True):
    """
    Return mixed-fuel MFSP in $/GGE.

    Dual-mode behavior
    ------------------
    market mode:
        gasoline / jet / diesel remain the economic products for TEA/IRR.
        If as_reporter=True, this function computes a reporting-only MFSP
        for the bookkeeping mixed_fuel stream, then restores original prices.

    mfsp mode:
        gasoline / jet / diesel should already be zero-priced, and mixed_fuel
        is treated as the economic product. This function solves mixed_fuel.price
        directly and leaves it set on the stream.

    Parameters
    ----------
    sys : qs.System
    print_msg : bool
        Whether to print the MFSP.
    as_reporter : bool
        In market mode, temporarily zero cut prices and solve mixed_fuel.price
        only for reporting. Original prices are restored afterward.

    Returns
    -------
    float
        MFSP of mixed_fuel in $/GGE.
    """
    s = sys.flowsheet.stream
    mixed_fuel = s.mixed_fuel

    if mixed_fuel.F_mass <= 0 or mixed_fuel.HHV <= 0:
        return np.nan

    # -----------------------------
    # MFSP economic mode
    # -----------------------------
    if SAF_ECON_MODE == "mfsp":
        mixed_fuel.price = sys.TEA.solve_price(mixed_fuel)
        mfsp = mixed_fuel.cost / get_GGE(sys, mixed_fuel, annual=False)
        if print_msg:
            print(f"Minimum selling price of all fuel is ${mfsp:.2f}/GGE.")
        return float(mfsp)

    # -----------------------------
    # Market/cuts mode
    # -----------------------------
    if SAF_ECON_MODE == "market":
        if not as_reporter:
            # not meaningful to solve MFSP while cuts remain the sold products
            return np.nan

        # Save original prices
        orig_prices = {
            "gasoline": getattr(s.gasoline, "price", 0.0),
            "jet": getattr(s.jet, "price", 0.0),
            "diesel": getattr(s.diesel, "price", 0.0),
            "mixed_fuel": getattr(mixed_fuel, "price", 0.0),
        }

        try:
            # Temporarily convert to MFSP-style pricing basis
            s.gasoline.price = 0.0
            s.jet.price = 0.0
            s.diesel.price = 0.0
            mixed_fuel.price = sys.TEA.solve_price(mixed_fuel)

            mfsp = mixed_fuel.cost / get_GGE(sys, mixed_fuel, annual=False)

            if print_msg:
                print(f"Reporting-only mixed-fuel MFSP is ${mfsp:.2f}/GGE.")

            return float(mfsp)

        finally:
            # Restore market-mode prices so IRR logic remains cut-based
            s.gasoline.price = orig_prices["gasoline"]
            s.jet.price = orig_prices["jet"]
            s.diesel.price = orig_prices["diesel"]
            mixed_fuel.price = orig_prices["mixed_fuel"]

    raise ValueError(f"Unknown SAF_ECON_MODE: {SAF_ECON_MODE}")

# In kg CO2e/GGE
def get_GWP(sys, print_msg=False):
    """
    Return mixed-fuel GWP in kg CO2e / GGE.

    Dual-mode behavior
    ------------------
    market mode:
        gasoline / jet / diesel are the economic products, but mixed_fuel
        is still maintained as a bookkeeping blend stream for reporting.
        This function reports aggregate fuel-pool GWP on a GGE basis.

    mfsp mode:
        mixed_fuel is also the economic aggregate product, so the same
        calculation is still valid.

    Parameters
    ----------
    sys : qs.System
    print_msg : bool
        Whether to print the GWP.

    Returns
    -------
    float
        GWP of mixed_fuel in kg CO2e / GGE.
    """
    if (not hasattr(sys, "LCA")) or (sys.LCA is None):
        return np.nan

    s = sys.flowsheet.stream
    if not hasattr(s, "mixed_fuel"):
        return np.nan

    mixed_fuel = s.mixed_fuel

    if mixed_fuel.F_mass <= 0 or mixed_fuel.HHV <= 0:
        return np.nan

    try:
        all_impacts = sys.LCA.get_allocated_impacts(
            streams=(mixed_fuel,),
            operation_only=True,
            annual=True,
        )
        gge_annual = get_GGE(sys, mixed_fuel, annual=True)
        if not np.isfinite(gge_annual) or gge_annual <= 0:
            return np.nan

        gwp = float(all_impacts["GWP"]) / float(gge_annual)

        if print_msg:
            print(f"Global warming potential of all fuel is {gwp:.2f} kg CO2e/GGE.")

        return float(gwp)

    except Exception:
        return np.nan
def get_fuel_properties(sys, fuel):
    HHV = fuel.HHV/fuel.F_mass/1e3 # MJ/kg
    rho = fuel.rho/_m3_to_gal # kg/gal
    return HHV, rho, get_GGE(sys, fuel, annual=False)

def simulate_and_print(sys, save_report=False):
    sys.simulate()
    rm = fuel_blend_ratios(sys, basis="mass")
    re = fuel_blend_ratios(sys, basis="gge")
    print("\nMixed fuel blend ratios into 'mixed_fuel'")
    print("----------------------------------------")
    print(f"Mass basis: gasoline={rm[0]:.3f}, jet={rm[1]:.3f}, diesel={rm[2]:.3f}")
    print(f"GGE basis : gasoline={re[0]:.3f}, jet={re[1]:.3f}, diesel={re[2]:.3f}")
  

    MFSP = get_MFSP(sys, print_msg=True)
    # mixed_fuel.price = 1.27 #$/kg, 1.27$/kg at fossil jet price with 45Z (1$/gallon or 0.31/kg transfereable)
     # WA $2.09 case ($0.75/L or $0.96/kg fossil jet A Yao et al. 2025, $0.31/kg 45Z, $0.62 state SAF credit, upto $0.20/kg LCFS)  
     # (IL $1.73/kg,  $0.46/kg state SAF credit), Green Premiums $0.31-0.62/kg excluded
    mf = sys.flowsheet.stream.mixed_fuel
    print("F_mass kg/hr:", mf.F_mass)
    print("HHV kJ/hr:", mf.HHV)
    print("price $/kg:", mf.price)
    print("cost $/hr:", mf.cost)
    print("GGE/hr:", get_GGE(sys, mf, annual=False))
    print("MFSP $/GGE:", mf.cost / get_GGE(sys, mf, annual=False))
    print("Implied GGE/kg:", get_GGE(sys, mf, annual=False) / mf.F_mass if mf.F_mass > 0 else np.nan)

    tea = sys.TEA
    c = qs.currency
    for attr in ("NPV", "AOC", "sales", "net_earnings"):
        uom = c if attr in ("NPV", "CAPEX") else (c + "/yr")
        print(f"{attr} is {getattr(tea, attr):,.0f} {uom}")
    
    GWP = get_GWP(sys, print_msg=True)

    if save_report:
        sys.save_report(file=os.path.join(os.getcwd(), f"sys_{sys.flowsheet.ID}.xlsx"))

    return MFSP


if __name__ == "__main__":
    sys = create_system(
        feedstock_id="food",
        decentralized_HTL=True,        
        decentralized_upgrading=False,
        include_PSA=True,
        include_EC=False,
    )
    dct = globals()
    dct.update(sys.flowsheet.to_dict())
    simulate_and_print(sys)