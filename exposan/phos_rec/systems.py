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
from qsdsan import unit_operations as qsu
from exposan.phos_rec import _load_components, create_tea
from exposan.phos_rec import _sanunits as su

# unit & basis conversion
_C_to_K = 273.15
_MMBTU_to_MJ = auom('BTU').conversion_factor('MJ')*10**6
_ton_to_kg = auom('ton').conversion_factor('kg')

# MW
MW_CH4 = 16

# year
lifetime = 30

# hr/year
operation_hours = 365*24

__all__ = (
    'get_precipitation_supernatant_price',
    'create_system'
)

def get_precipitation_supernatant_price(food_sludge_ratio, fermentation_time):
    VFA_conc_dict = {
        (0, 0): 229.836,
        (0, 12): 229.836,
        (0, 24): 229.836,
        (0, 36): 229.836,
        (0, 48): 229.836,
        (0, 60): 229.836,
        (0, 72): 229.836,
        (0, 84): 315.014,
        (0, 96): 324.748,
        (0, 108): 474.744,
        (0, 120): 503.468,
        (0, 132): 574.002,
        (1/3, 0): 229.836,
        (1/3, 12): 318.997,
        (1/3, 24): 891.924,
        (1/3, 36): 930.369,
        (1/3, 48): 942.929,
        (1/3, 60): 1019.904,
        (1/3, 72): 1030.104,
        (1/3, 84): 1058.436,
        (1/3, 96): 1065.347,
        (1/3, 108): 1235.25,
        (1/3, 120): 1295.057,
        (1/3, 132): 1486.757,
        (2/3, 0): 229.836,
        (2/3, 12): 619.055,
        (2/3, 24): 1342.854,
        (2/3, 36): 1975.393,
        (2/3, 48): 2195.785,
        (2/3, 60): 2540.621,
        (2/3, 72): 2677.336,
        (2/3, 84): 2740.366,
        (2/3, 96): 2754.408,
        (2/3, 108): 2835.845,
        (2/3, 120): 2882.61,
        (2/3, 132): 3123.319,
        (3/3, 0): 229.836,
        (3/3, 12): 755.134,
        (3/3, 24): 1534.219,
        (3/3, 36): 2006.594,
        (3/3, 48): 2624.859,
        (3/3, 60): 3049.32,
        (3/3, 72): 3413.96,
        (3/3, 84): 3716.996,
        (3/3, 96): 3763.005,
        (3/3, 108): 3953.696,
        (3/3, 120): 4114.03,
        (3/3, 132): 4500,
        (4/3, 132): 4650
      }
        
    vfa_conc = VFA_conc_dict[(food_sludge_ratio), fermentation_time]
    price = 0.002*(vfa_conc/4500)
    
    return price

# TODO: add in the manuscript: except when X-axis is size, the system input is assumed to be 100 dry tonne solids·day-1
# sludge_cost_credit: $/tonne dry solids (https://pubs.acs.org/doi/10.1021/es505329q)
def create_system(dry_solids_tonne_per_day=100, food_sludge_ratio=1, fermentation_time=132, sludge_cost_credit=300, sludge_CI_credit=1300, perspective='FePO4'):
    
    if perspective not in ['FePO4','sludge']:
        raise KeyError('perspective can only be FePO4 or sludge')
    
    flowsheet_ID = 'phos_rec'
    
    # clear flowsheet and registry for reloading
    if hasattr(qs.main_flowsheet.flowsheet, flowsheet_ID):
        getattr(qs.main_flowsheet.flowsheet, flowsheet_ID).clear()
        clear_lca_registries()
    
    bst.CE = qs.CEPCI_by_year[2023]

    flowsheet = qs.Flowsheet(flowsheet_ID)
    stream = flowsheet.stream
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    _load_components()
    
    conversion_factor = dry_solids_tonne_per_day/6730
    
    # Inert represents other inorganics; now VSS/TSS = 0.74, which is reasonable
    fe_sludge = qs.WasteStream(ID='fe_sludge', Fe3=180*conversion_factor, Org=5000*conversion_factor,
                               PO4=300*conversion_factor, Water=1000000*conversion_factor, Ca2=150*conversion_factor,
                               Mg2=100*conversion_factor, Inert=1000*conversion_factor, units='tonne/d')
    
    P1 = qsu.SludgePump(ID='P1', ins=fe_sludge, outs='sludge_to_ST')
    
    # the same as AF: 5 days of reaction + 3 hr for cleaning and unloading
    ST = qsu.StorageTank(ID='ST', ins=P1-0, outs='sludge_to_P2', tau=5*24+3)
    
    P2 = qsu.SludgePump(ID='P2', ins=ST-0, outs='sludge_to_AF')
    
    # TODO: include this assumption in the manuscript or the SI
    # assume food waste is neutral, e.g., it has no cost, CI or related credits
    AF = su.AcidogenicFermenter(ID='AF', ins=(P2-0, 'food_waste'), outs=('fermentate', 'fermentation_gas'),
                                food_sludge_ratio=food_sludge_ratio, fermentation_time=fermentation_time)
    
    FC = su.SludgeCentrifugeWithElementFlow(ID='FC', ins=AF-0, outs=('fermentation_supernatant', 'residue'),
                                            sludge_moisture=0.85, solids=('Inert','Residue'))
    
    SP = su.SelectivePrecipitation(ID='SP', ins=(FC-0, 'acid', 'oxidant'), outs='slurry')
    
    PC = su.SludgeCentrifugeWithElementFlow(ID='PC', ins=SP-0, outs=('precipitation_supernatant', 'precipitate'),
                                            sludge_moisture=0.935, solids=('FePO4_2H2O',))
    
    HD = su.HeatDrying(ID='HD', ins=(PC-1, 'heat_drying_natural_gas'), outs=('dried_precipitate', 'heat_drying_vapor'),
                       T=105+_C_to_K)
    
    # feed and product are solids, therefore, no internal heat exchangers for SI
    # vapor is 700 °C, but its heat is not recovered to supply AF or SP, which
    # operate at much lower temperatures (37 °C and 40 °C, respectively)
    # when ΔT is large, there can be problems such as scale mismatch and operational complexity
    SI = su.Sintering(ID='SI', ins=(HD-0, 'sintering_natural_gas', 'air'), outs=('product', 'sintering_vapor'))
    # for recovery calculation in SI @property
    SI.feedstock = fe_sludge
    
    sys = qs.System.from_units(ID='phos_rec', units=list(flowsheet.unit), operating_hours=operation_hours)
    
    sys.simulate()
    
    sludge_dry_kg_hr = fe_sludge.F_mass - fe_sludge.imass['Water']
    fe_sludge.price = -sludge_cost_credit*(sludge_dry_kg_hr/1000)/fe_sludge.F_mass
    
    # 0.08 $/kg, from Everbatt2023
    SP.ins[1].price = 0.08 * stream.acid.imass['H2SO4']/stream.acid.F_mass
    
    # 1.46 $/kg, from Everbatt2023
    SP.ins[2].price = 1.46 * stream.oxidant.imass['H2O2']/stream.oxidant.F_mass
    
    # 3.49672 $/kmol, from BioSTEAM
    HD.ins[1].price = 3.49672/MW_CH4
    
    # TODO: any reference?
    # 15000 rmb/ton -> 2.14 $/kg
    SI.outs[0].price = 2.14
    
    SI.ins[1].price = 3.49672/MW_CH4
    
    # 62.28 $/ton, https://www.wasteoptima.com/blog/eref-2024-landfill-tipping-fees?utm_source
    FC.outs[1].price = - 62.28/_ton_to_kg
    
    PC.outs[0].price = get_precipitation_supernatant_price(AF.food_sludge_ratio, AF.fermentation_time)
    
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
    
    qs.StreamImpactItem(
        ID='fe_sludge',
        linked_stream=stream.fe_sludge,
        GlobalWarming=-sludge_CI_credit*(sludge_dry_kg_hr/1000)/fe_sludge.F_mass
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
        GlobalWarming=0.573511217751451*stream.heat_drying_natural_gas.F_vol/stream.heat_drying_natural_gas.F_mass
    )
    
    # market for natural gas, low pressure, United States of America (US)
    qs.StreamImpactItem(
        ID='sintering_natural_gas',
        linked_stream=stream.sintering_natural_gas,
        GlobalWarming=0.573511217751451*stream.sintering_natural_gas.F_vol/stream.sintering_natural_gas.F_mass
    )
    
    # the CI of FePO4: https://www.sciencedirect.com/science/article/pii/S1383586624004258#s0085  2-6 kg CO2-eq/kg FePO4
    qs.StreamImpactItem(
    ID='FePO4_credit',
    linked_stream=stream.product,
    GlobalWarming=-4
    )
    
    # TODO: any reference?
    # TODO: confirm -0.002 is for for AF.food_sludge_ratio = 1 not 4/3
    # -0.002 for AF.food_sludge_ratio = 1 and AF.fermentation_time = 132
    qs.StreamImpactItem(
    ID='VFA_credit',
    linked_stream=stream.precipitation_supernatant,
    GlobalWarming=-0.002/4500*get_precipitation_supernatant_price(AF.food_sludge_ratio, AF.fermentation_time)
    )
    
    # TODO: include this assumption in the manuscript or the SI
    # assume no additional labor cost
    if perspective == 'FePO4':
        create_tea(sys,
                   duration=(2023, 2023+lifetime),
                   IRR_value=0.1,
                   income_tax_value=0.3,
                   finance_interest_value=0.08,
                   labor_cost_value=0)
    else:
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