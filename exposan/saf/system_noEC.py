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

'''

# !!! Temporarily ignoring warnings
# import warnings
# warnings.filterwarnings('ignore')
import os, numpy as np, biosteam as bst, qsdsan as qs
from biosteam.units import IsenthalpicValve
from biosteam import settings
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import (
    create_tea,    
    )
from exposan.biobinder import _units as bbu, find_Lr_Hr
from exposan.saf import (
    create_components,
    # data_path,
    # results_path,
    # _load_components,
    # _load_process_settings,
    # create_tea,
    _units as u,
    )

_psi_to_Pa = 6894.76
_m3_to_gal = 264.172

# All in 2020 $/kg unless otherwise noted
price_dct = {
    'H2': 1.61, # Feng et al.
    'HCcatalyst': 3.52, # Fe-ZSM5, CatCost modified from ZSM5
    'HTcatalyst': 75.18, # Pd/Al2O3, CatCost modified from 2% Pt/TiO2
    'nature_gas': 0.1685,
    'process_water': 0,
    'gasoline': 2.5, # 2024$/gal
    'jet': 3.53, # 2024$/gal
    'diesel': 3.45, # 2024$/gal
    'solids': -0,
    'wastewater': -0, #!!! need to update
    }

flowsheet_ID = 'saf_noEC'
flowsheet = qs.Flowsheet(flowsheet_ID)
qs.main_flowsheet.set_flowsheet(flowsheet)
saf_cmps = create_components(set_thermo=True)

feedstock = qs.Stream('feedstock')
feedstock_water = qs.Stream('feedstock_water', Water=1)

FeedstockTrans = bbu.Transportation(
    'FeedstockTrans',
    ins=(feedstock, 'transportation_surrogate'),
    outs=('transported_feedstock',),
    N_unit=1,
    copy_ins_from_outs=True,
    transportation_distance=78, # km ref [1]
    )

FeedstockWaterPump = qsu.Pump('FeedstockWaterPump', ins=feedstock_water)

#!!! Need to update the composition (moisture/ash)
moisture = 0.7566
feedstock_composition = {
    'Water': moisture,
    'Lipids': (1-moisture)*0.5315,
    'Proteins': (1-moisture)*0.0255,
    'Carbohydrates': (1-moisture)*0.3816,
    'Ash': (1-moisture)*0.0614,
    }
FeedstockCond = bbu.Conditioning(
    'FeedstockCond', ins=(FeedstockTrans-0, FeedstockWaterPump-0),
    outs='conditioned_feedstock',
    feedstock_composition=feedstock_composition,
    feedstock_dry_flowrate=110*907.185/(24*0.9), # 110 dry sludge tpd ref [1]; 90% upfactor
    N_unit=1,
    )
@FeedstockCond.add_specification
def adjust_feedstock_composition():
    FeedstockCond._run()
    FeedstockTrans._run()
    FeedstockWaterPump._run()

MixedFeedstockPump = qsu.Pump('MixedFeedstockPump', ins=FeedstockCond-0)

# =============================================================================
# Hydrothermal Liquefaction (HTL)
# =============================================================================
HTLpreheater = qsu.HXutility(
    'HTLpreheater',
    ins=MixedFeedstockPump-0, outs='heated_feedstock', T=280+273.15,
    U=0.0198739, init_with='Stream', rigorous=True
    )

HTL = u.HydrothermalLiquefaction(
    'HTL', ins=HTLpreheater-0,
    outs=('','','HTL_crude','HTL_char'),
    dw_yields={
        'gas': 0.006,
        'aqueous': 0.192,
        'biocrude': 0.802,
        'char': 0,
        },
    gas_composition={'CO2': 1},
    aqueous_composition={'HTLaqueous': 1},
    biocrude_composition={'Biocrude': 1},
    char_composition={'HTLchar': 1},
    )
HTL.register_alias('HydrothermalLiquefaction')

CrudePump = qsu.Pump('CrudePump', ins=HTL-2, outs='crude_to_dist', P=1530.0*_psi_to_Pa,
          init_with='Stream')
# Jones 2014: 1530.0 psia

# Light (water): medium (biocrude): heavy (char)
# Split off the light compounds (bp<150°C)
crude_fracs = (0.0339, 0.8104, 0.1557)
CrudeSplitter = bbu.BiocrudeSplitter(
    'CrudeSplitter', ins=CrudePump-0, outs='splitted_crude',
    biocrude_IDs=('HTLbiocrude'),
    cutoff_Tbs=(150+273.15, 300+273.15,),
    cutoff_fracs=crude_fracs,
    )

# Separate water from organics
CrudeLightDis = qsu.ShortcutColumn(
    'CrudeDis', ins=CrudeSplitter-0,
    outs=('crude_light','crude_medium_heavy'),
    LHK=CrudeSplitter.keys[0],
    P=50*_psi_to_Pa,
    Lr=0.87,
    Hr=0.98,
    k=2, is_divided=True)
# results_df, Lr, Hr = find_Lr_Hr(CrudeLightDis, target_light_frac=crude_fracs[0])

CrudeLightFlash = qsu.Flash('CrudeLightFlash', ins=CrudeLightDis-0,
                            T=298.15, P=101325,)
                            # thermo=settings.thermo.ideal())
HTLgasMixer = qsu.Mixer('HTLgasMixer', ins=(HTL-0, CrudeLightFlash-0), outs='HTL_gas')
HTLaqMixer = qsu.Mixer('HTLaqMixer', ins=(HTL-1, CrudeLightFlash-1), outs='HTL_aq')

# Separate biocrude from char
CrudeHeavyDis = qsu.ShortcutColumn(
    'CrudeHeavyDis', ins=CrudeLightDis-1,
    outs=('crude_medium','char'),
    LHK=CrudeSplitter.keys[1],
    P=50*_psi_to_Pa,
    Lr=0.89,
    Hr=0.85,
    k=2, is_divided=True)
# results_df, Lr, Hr = find_Lr_Hr(CrudeHeavyDis, target_light_frac=crude_fracs[1]/(1-crude_fracs[0]))
    

# =============================================================================
# Hydrocracking
# =============================================================================

# include_PSA = False # want to compare with vs. w/o PSA

# External H2, will be updated after HT and HC
H2 = qs.WasteStream('H2', H2=1, price=price_dct['H2'])
H2splitter= qsu.ReversedSplitter('H2splitter', ins='H2', outs=('HC_H2', 'HT_H2'),
                            init_with='WasteStream')

# 10 wt% Fe-ZSM
HCcatalyst_in = qs.WasteStream('HCcatalyst_in', HCcatalyst=1, price=price_dct['HCcatalyst'])

# Original as in Feng et al., adjust to 73.35
HC = u.Hydrocracking(
    'HC',
    ins=(CrudeHeavyDis-0, H2splitter-0, HCcatalyst_in),
    outs=('HC_out','HCcatalyst_out'),
    catalyst_ID='HCcatalyst',
    HCin_T=None,
    HCrxn_T=400+273.15,
    gas_yield=0.2665,
    oil_yield=0.7335,
    hydrogen_P=1500*_psi_to_Pa,
    hydrogen_rxned_to_inf_oil=0.0111,
    )
HC.register_alias('Hydrocracking')

HC_HX = qsu.HXutility(
    'HC_HX', ins=HC-0, outs='cooled_HC_eff', T=60+273.15,
    init_with='Stream', rigorous=True)

# To depressurize products
HC_IV = IsenthalpicValve('HC_IV', ins=HC_HX-0, outs='cooled_depressed_HC_eff', P=30*6894.76, vle=True)

# To separate products
HCflash = qsu.Flash('HC_Flash', ins=HC_IV-0, outs=('HC_fuel_gas','HC_liquid'),
                    T=60.2+273.15, P=30*_psi_to_Pa,)

HCpump = qsu.Pump('HCpump', ins=HCflash-1, init_with='Stream')

# Separate water from oil
HCliquidSplitter = qsu.Splitter('HCliquidSplitter', ins=HCpump-0,
                                outs=('HC_ww','HC_oil'),
                                split={'H2O':1}, init_with='Stream')


# =============================================================================
# Hydrotreating
# =============================================================================

# Pd/Al2O3
HTcatalyst_in = qs.WasteStream('HCcatalyst_in', HTcatalyst=1, price=price_dct['HTcatalyst'])

HT = u.Hydrotreating(
    'HT',
    ins=(HCliquidSplitter-1, H2splitter-1, HTcatalyst_in),
    outs=('HTout','HTcatalyst_out'),
    catalyst_ID='HTcatalyst',
    HTin_T=None,
    HTrxn_T=300+273.15,
    gas_yield=0.2143,
    oil_yield=0.8637,
    hydrogen_P=1500*_psi_to_Pa,
    hydrogen_rxned_to_inf_oil=0.0207,
    )
HT.register_alias('Hydrotreating')


HT_HX = qsu.HXutility('HT_HX',ins=HT-0, outs='cooled_HT_eff', T=60+273.15,
                      init_with='Stream', rigorous=True)

HT_IV = IsenthalpicValve('HT_IV', ins=HT_HX-0, outs='cooled_depressed_HT_eff',
                         P=717.4*_psi_to_Pa, vle=True)

HTflash = qsu.Flash('HTflash', ins=HT_IV-0, outs=('HT_fuel_gas','HT_oil'),
                    T=43+273.15, P=55*_psi_to_Pa)

HTpump = qsu.Pump('HTpump', ins=HTflash-1, init_with='Stream')

# Separate water from oil
HTliquidSplitter = qsu.Splitter('HTliquidSplitter', ins=HTpump-0,
                                outs=('HT_ww','HT_oil'),
                                split={'H2O':1}, init_with='Stream')

# Light (gasoline, <C8): medium (jet, C8-C14): heavy (diesel, >C14)
oil_fracs = (0.2143, 0.5638, 0.2066)

# Separate gasoline from jet and diesel
GasolineDis = qsu.ShortcutColumn(
    'OilLightDis', ins=HTliquidSplitter-1,
    outs=('hot_gasoline','jet_diesel'),
    LHK=('OXYLENE', 'C9H20'),
    Lr=0.99,
    Hr=0.99,
    k=2, is_divided=True)
# Lr_range = Hr_range = np.linspace(0.05, 0.95, 19)
# Lr_range = Hr_range = np.linspace(0.01, 0.2, 20)
# results_df, Lr, Hr = find_Lr_Hr(GasolineDis, Lr_trial_range=Lr_range, Hr_trial_range=Hr_range, target_light_frac=oil_fracs[0])

GasolineFlash = qsu.Flash('GasolineFlash', ins=GasolineDis-0, outs=('', 'cooled_gasoline',),
                          T=298.15, P=101325)

# Separate jet from diesel
JetDis = qsu.ShortcutColumn(
    'JetDis', ins=GasolineDis-1,
    outs=('hot_jet','hot_diesel'),
    LHK=('C14H22', 'C15H32'),
    Lr=0.75,
    Hr=0.2,
    k=2, is_divided=True)
# Lr_range = Hr_range = np.linspace(0.05, 0.95, 19)
# Lr_range = Hr_range = np.linspace(0.01, 0.2, 20)
# results_df, Lr, Hr = find_Lr_Hr(JetDis, Lr_trial_range=Lr_range, Hr_trial_range=Hr_range, target_light_frac=oil_fracs[1]/(1-oil_fracs[0]))

JetFlash = qsu.Flash('JetFlash', ins=JetDis-0, outs=('', 'cooled_jet',), T=298.15, P=101325)

DieselHX = qsu.HXutility('DieselHX',ins=JetDis-1, outs='cooled_diesel', T=298.15,
                         init_with='Stream', rigorous=True)


# =============================================================================
# Electrochemical Units
# =============================================================================


# =============================================================================
# Products and Wastes
# =============================================================================

GasolinePC = qsu.PhaseChanger('GasolinePC', ins=GasolineFlash-1)
gasoline = qs.WasteStream('gasoline', Gasoline=1)
gasoline.price = price_dct['gasoline']/(gasoline.rho/_m3_to_gal)
# Storage time assumed to be 3 days per [1]
GasolineTank = qsu.StorageTank('GasolineTank', ins=GasolinePC-0, outs=(gasoline),
                                tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')

JetPC = qsu.PhaseChanger('JetPC', ins=JetFlash-1)
jet = qs.WasteStream('jet', Jet=1)
jet.price = price_dct['jet']/(jet.rho/_m3_to_gal)
JetTank = qsu.StorageTank('JetTank', ins=JetPC-0, outs=(jet,),
                          tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')

DieselPC = qsu.PhaseChanger('DieselPC', ins=DieselHX-0)
diesel = qs.WasteStream('diesel', Jet=1)
diesel.price = price_dct['diesel']/(diesel.rho/_m3_to_gal)
DieselTank = qsu.StorageTank('DieselTank', ins=DieselPC-0, outs=(diesel,),
                             tau=3*24, init_with='WasteStream', vessel_material='Carbon steel')

# Gas emissions
WasteGasMixer = qsu.Mixer('WasteGasMixer', ins=(HTLgasMixer-0,),
                      outs=('gas_emissions'), init_with='Stream')

# All fuel gases sent to CHP for heat generation
FuelGasMixer = qsu.Mixer('FuelGasMixer',
                         ins=(HCflash-0, HTflash-0, GasolineFlash-0, JetFlash-0,),
                         outs=('fuel_gas'), init_with='Stream')
# Run this toward the end to make sure H2 flowrate can be updated
@FuelGasMixer.add_specification
def update_H2_flow():
    H2splitter._run()
    FuelGasMixer._run()

# All wastewater, assumed to be sent to municipal wastewater treatment plant
wastewater = qs.WasteStream('wastewater', price=price_dct['wastewater'])
WWmixer = qsu.Mixer('WWmixer',
                    ins=(HTLaqMixer-0, HCliquidSplitter-0, HTliquidSplitter-0),
                    outs=wastewater, init_with='Stream')

# All solids, assumed to be disposed to landfill
disposed_solids = qs.WasteStream('solids', price=price_dct['solids'])
SolidsMixer = qsu.Mixer('SolidsMixer', ins=CrudeHeavyDis-1,
                    outs=disposed_solids, init_with='Stream')

# =============================================================================
# Facilities
# =============================================================================

qsu.HeatExchangerNetwork('HXN', T_min_app=86, force_ideal_thermo=True)
# 86 K: Jones et al. PNNL, 2014

nature_gas = qs.WasteStream('nature_gas', CH4=1, price=price_dct['nature_gas'])
CHP = qsu.CombinedHeatPower('CHP', 
                            ins=(FuelGasMixer-0, 'natural_gas', 'air'),
                            outs=('emission','solid_ash'), init_with='WasteStream',
                            supplement_power_utility=False)


# %%

sys = qs.System.from_units(
    'sys_noEC',
    units=list(flowsheet.unit),
    operating_hours=7920, # 90% uptime
    )
for unit in sys.units: unit.include_construction = False

tea = create_tea(sys, IRR_value=0.1, income_tax_value=0.21, finance_interest_value=0.08,
                  labor_cost_value=1.81*10**6)

# lca = qs.LCA(
#     system=sys,
#     lifetime=lifetime,
#     uptime_ratio=sys.operating_hours/(365*24),
#     Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
#     # Heating=lambda:sys.get_heating_duty()/1000*lifetime,
#     Cooling=lambda:sys.get_cooling_duty()/1000*lifetime,
#     )

def simulate_and_print():
    sys.simulate()
    jet.price = tea.solve_price(jet)
    jet_price = jet.price*(jet.rho/_m3_to_gal)
    print(f'Minimum selling price of the jet is ${jet_price:.2f}/kg.')
    
    c = qs.currency
    for attr in ('NPV','AOC', 'sales', 'net_earnings'):
        uom = c if attr in ('NPV', 'CAPEX') else (c+('/yr'))
        print(f'{attr} is {getattr(tea, attr):,.0f} {uom}')
        
    # all_impacts = lca.get_allocated_impacts(streams=(biobinder,), operation_only=True, annual=True)
    # GWP = all_impacts['GlobalWarming']/(biobinder.F_mass*lca.system.operating_hours)
    # print(f'Global warming potential of the biobinder is {GWP:.4f} kg CO2e/kg.')

if __name__ == '__main__':
    simulate_and_print()