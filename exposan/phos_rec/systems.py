#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Xuan Wang <easonbiubiu99@gmail.com>
    
    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import biosteam as bst, qsdsan as qs
from qsdsan.utils import clear_lca_registries, auom
from qsdsan import sanunits as qsu
from exposan.phos_rec import _load_components, create_tea
from exposan.phos_rec import _sanunits as su

# unit & basis conversion
_C_to_K = 273.15
_MMBTU_to_MJ = auom('BTU').conversion_factor('MJ') * 10**6
_ton_to_kg = auom('ton').conversion_factor('kg')

# year
lifetime = 30

# hr/year
operation_hours = 365*24

__all__ = (
    'create_system',
)

def create_system(water_mass_flow=1000000):
    
    flowsheet_ID = 'phos_rec'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]

    flowsheet = qs.Flowsheet(flowsheet_ID)
    unit = flowsheet.unit
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    # TODO: instead, assign concentrations to components except water
    # TODO: consider adding a cost for fe_sludge as credits (e.g., avoid landfilling, etc.)
    # Inert represents other inorganics; now VSS/TSS = 0.74, which is reasonable
    fe_sludge = qs.WasteStream(ID='fe_sludge', Fe3=180, Org=5000, PO4=300, Water=water_mass_flow,
                               Ca2=150, Mg2=100, Inert=1000, units='kg/d')
    
    P1 = qsu.SludgePump(ID='P1', ins=fe_sludge, outs='sludge_to_ST')
    
    # the same as AF: 4 days of reaction+ 3 hr for cleaning and unloading
    ST = qsu.StorageTank(ID='ST', ins=P1-0, outs='sludge_to_P2', tau=4*24+3)
    
    P2 = qsu.SludgePump(ID='P2', ins=ST-0, outs='sludge_to_AF')
    
    AF = su.AcidogenicFermenter(ID='AF', ins=(P2-0, 'food_waste'), outs=('fermentate', 'fermentation_gas'),
                                food_sludge_ratio=1)
    
    #FC = qsu.SludgeCentrifuge(ID='FC', ins=AF-0, outs=('fermentation_supernatant', 'residue'),
                              #sludge_moisture=0.85, solids=('Inert','Residue'))
    FC = su.SludgeCentrifugeWithElementFlow(
    ID='FC',
    ins=AF-0,
    outs=('fermentation_supernatant', 'residue'),
    sludge_moisture=0.85,
    solids=('Inert', 'Residue')
)
    
    
    
    SP = su.SelectivePrecipitation(ID='SP', ins=(FC-0, 'acid', 'oxidant'), outs='slurry')
    
    # TODO: Org seems low, but conservative: less organics for heat generation, and more go back to WRRF headworks
    #PC = qsu.SludgeCentrifuge(ID='PC', ins=SP-0, outs=('precipitation_supernatant', 'precipitate'), 
                              #sludge_moisture=0.92, solids=('FePO4_2H2O',))
    PC = su.SludgeCentrifugeWithElementFlow(
    ID='PC',
    ins=SP-0,
    outs=('precipitation_supernatant', 'precipitate'),
    sludge_moisture=0.92,
    solids=('FePO4_2H2O',)
)                         
    
    HD = su.HeatDrying(ID='HD', ins=(PC-1, 'heat_drying_natural_gas'), outs=('dried_precipitate', 'heat_drying_vapor'),
                       T=105+_C_to_K)
    
    # feed and product are solids, therefore, no internal heat exchangers for SI
    # vapor is 700 °C, but its heat is not recovered to supply AF or SP, which
    # operate at much lower temperatures (37 °C and 40 °C, respectively)
    # when ΔT is large, there can be problems such as scale mismatch and operational complexity
    SI = su.Sintering(ID='SI', ins=(HD-0, 'sintering_natural_gas', 'air'), outs=('product', 'sintering_vapor'))
    # ---- for recovery calculation in SI @property
    SI.feedstock = fe_sludge
    SI.food_waste = AF.ins[1]   # 可选：如果你想把 food_waste 计入 C 的分母
    
    sys = qs.System.from_units(ID='phos_rec', units=list(flowsheet.unit), operating_hours=operation_hours)
    
    sys.simulate()
    
    # TODO: may need update
    # fe_sludge.price = - 0.5 * (fe_sludge.F_mass - fe_sludge.imass['Water']) / fe_sludge.F_mass
    sludge_credit_used_per_ton_dry = 62.28 #$/ton dry solids in https://www.wasteoptima.com/blog/eref-2024-landfill-tipping-fees?utm_source
    dry_kg_d = fe_sludge.F_mass - fe_sludge.imass['Water']
    fe_sludge.price = - sludge_credit_used_per_ton_dry * (dry_kg_d/1000) / fe_sludge.F_mass
    
    # price of H2SO4 is 0.08$/kg, from Everbatt2023
    SP.ins[1].price = 0.08 * stream.acid.imass['H2SO4'] / stream.acid.F_mass
    
    # price of H2O2 is 1.46$/kg, from Everbatt2023
    SP.ins[2].price = 1.46 * stream.oxidant.imass['H2O2'] / stream.oxidant.F_mass
    
    # TODO: biosteam has a default price
    # price of natural gas is 7.7$/MMBTU, from Everbatt2023, LHV_ng = 50.1 MJ/kg
    # natural gas consumption was calculated using the lower heating value (LHV = 50.1 MJ/kg), assuming that the latent heat of water vapor in the flue gas was not recovered
    HD.ins[1].price = 7.7 / _MMBTU_to_MJ * 50.1
    
    SI.ins[1].price = 7.7 / _MMBTU_to_MJ * 50.1
    
    # the landfill cost (tipping fee) can be found in https://www.wasteoptima.com/blog/eref-2024-landfill-tipping-fees?utm_source with the $62.28 per ton
    FC.outs[1].price = - 62.28 / _ton_to_kg
    # TODO: labor cost needs to be considered
    # TODO: not consider avoided costs for sludge landfilling
    # assume food waste is free, labor cost is 0 $/kg
    # electricity and steam use default prices
 
    GlobalWarming = qs.ImpactIndicator(ID='GlobalWarming',
                                       method='IPCC',
                                       category='environmental impact',
                                       unit='kg CO2-eq',
                                       description='Global Warming Potential')
    
    # market group for electricity, medium voltage, United States of America (US)
    Electricity = qs.ImpactItem('Electricity', functional_unit='kWh')
    Electricity.add_indicator(GlobalWarming, 0.484477212926751)
    
    # market for heat, from steam, in chemical industry, Europe (RER)
    Steam = qs.ImpactItem('Steam', functional_unit='MJ')
    Steam.add_indicator(GlobalWarming, 0.110608480188443)
    
    # market for sulfuric acid, Europe (RER):0.122270486489691
    # water is igonred
    qs.StreamImpactItem(
        ID='H2SO4_item',
        linked_stream=stream.acid,
        GlobalWarming=0.122270486489691*stream.acid.imass['H2SO4']/stream.acid.F_mass
    )
    
    # market for hydrogen peroxide, without water, in 50% solution state, Europe (RER): 1.7460300053961
    # water is igonred
    qs.StreamImpactItem(
        ID='H2O2_item',
        linked_stream=stream.oxidant,
        GlobalWarming=1.7460300053961*stream.oxidant.imass['H2O2']/stream.oxidant.F_mass
    )
    
    # TODO: may need update: 30 is a rough number to consider sludge landfilling at 80% moisture content
    # market for process-specific burdens, sanitary landfill, Rest-of-World (RoW): 0.00582207178953914
    qs.StreamImpactItem(
        ID='fe_sludge',
        linked_stream=stream.fe_sludge,
        GlobalWarming=-0.00582207178953914/30
    )
    
    # market for process-specific burdens, sanitary landfill, Rest-of-World (RoW): 0.00582207178953914
    qs.StreamImpactItem(
        ID='residue',
        linked_stream=stream.residue,
        GlobalWarming=0.00582207178953914
    )
    
    # assume biogenic CO2
    qs.StreamImpactItem(
        ID='fermentation_gas',
        linked_stream=stream.fermentation_gas,
        GlobalWarming=0
    )
    
    # market for natural gas, low pressure, United States of America (US)
    qs.StreamImpactItem(
        ID='heat_drying_natural_gas',
        linked_stream=stream.heat_drying_natural_gas,
        GlobalWarming=0.573511217751451 * stream.heat_drying_natural_gas.F_vol / stream.heat_drying_natural_gas.F_mass
    )
    
    # market for natural gas, low pressure, United States of America (US)
    qs.StreamImpactItem(
        ID='sintering_natural_gas',
        linked_stream=stream.sintering_natural_gas,
        GlobalWarming=0.573511217751451 * stream.sintering_natural_gas.F_vol / stream.sintering_natural_gas.F_mass
    )

    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=0) 
    
    qs.LCA(
        system=sys, lifetime=lifetime, lifetime_unit='yr',
        Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30,
        Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
    )
    
    return sys