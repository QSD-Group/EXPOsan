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

# MW
MW_CH4 = 16

# year
lifetime = 30

# hr/year
operation_hours = 365*24

__all__ = (
    'create_system',
)

# TODO: use temp_ratio to scale fe_sludge temporarily
def create_system(temp_ratio=1, food_sludge_ratio=1, HRT=132, sludge_credit=300, sludge_CI_credit=300):
    
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
    
    # TODO: instead, assign concentrations to components except water
    # Inert represents other inorganics; now VSS/TSS = 0.74, which is reasonable
    fe_sludge = qs.WasteStream(ID='fe_sludge', Fe3=180*temp_ratio, Org=5000*temp_ratio,
                               PO4=300*temp_ratio, Water=1000000*temp_ratio, Ca2=150*temp_ratio,
                               Mg2=100*temp_ratio, Inert=1000*temp_ratio, units='kg/d')
    
    P1 = qsu.SludgePump(ID='P1', ins=fe_sludge, outs='sludge_to_ST')
    
    # the same as AF: 5 days of reaction + 3 hr for cleaning and unloading
    ST = qsu.StorageTank(ID='ST', ins=P1-0, outs='sludge_to_P2', tau=5*24+3)
    
    P2 = qsu.SludgePump(ID='P2', ins=ST-0, outs='sludge_to_AF')
    
    AF = su.AcidogenicFermenter(ID='AF', ins=(P2-0, 'food_waste'), outs=('fermentate', 'fermentation_gas'),
                                food_sludge_ratio=food_sludge_ratio, HRT=HRT)
    
    FC = su.SludgeCentrifugeWithElementFlow(ID='FC', ins=AF-0, outs=('fermentation_supernatant', 'residue'),
                                            sludge_moisture=0.85, solids=('Inert','Residue'))
    
    SP = su.SelectivePrecipitation(ID='SP', ins=(FC-0, 'acid', 'oxidant'), outs='slurry')
    
    # TODO: may need to remove the scenario where food_sludge_ratio = 0 due to its very high costs and GHG emissions, also due to its different moisture content set here (only 0.99 works for the scenario where food_sludge_ratio = 0)
    # if food_sludge_ratio != 0:
    #     PC = su.SludgeCentrifugeWithElementFlow(ID='PC', ins=SP-0, outs=('precipitation_supernatant', 'precipitate'),
    #                                             sludge_moisture=0.92, solids=('FePO4_2H2O',))
    # else:
    #     PC = su.SludgeCentrifugeWithElementFlow(ID='PC', ins=SP-0, outs=('precipitation_supernatant', 'precipitate'),
    #                                             sludge_moisture=0.99, solids=('FePO4_2H2O',))
        
    # precip_mass = (SP-0).imass['FePO4_2H2O']
    
    # if precip_mass < 1e-6:
    #     pc_moisture = 0.99
    # else:
    #     pc_moisture = 0.80
        
    if food_sludge_ratio == 0:
        PC = su.SludgeCentrifugeWithElementFlow(ID='PC', ins=SP-0, outs=('precipitation_supernatant', 'precipitate'),
                                            sludge_moisture=0.999, solids=('FePO4_2H2O',))
    else:
        PC = su.SludgeCentrifugeWithElementFlow(ID='PC', ins=SP-0, outs=('precipitation_supernatant', 'precipitate'),
                                            sludge_moisture=0.92, solids=('FePO4_2H2O',))
        
        
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
    
    # the value：https://pubs.acs.org/doi/10.1021/es505329q
    sludge_credit_used_per_ton_dry = sludge_credit # $/ton dry solids
    temp_ratio = temp_ratio # size of WRRF
    sludge_dry_kg_d = fe_sludge.F_mass - fe_sludge.imass['Water']
    fe_sludge.price = -sludge_credit_used_per_ton_dry*(sludge_dry_kg_d/1000)/fe_sludge.F_mass
    
    # the value: https://www.sciencedirect.com/science/article/pii/S0959652621022757
    food_waste_credit_used_per_ton_dry = 75.8 # $/ton dry solids
    food_waste_dry_kg_d = AF.ins[1].F_mass - AF.ins[1].imass['Water']
    AF.ins[1].price = -food_waste_credit_used_per_ton_dry*(food_waste_dry_kg_d/1000)/fe_sludge.F_mass
    
    # price of H2SO4 is 0.08$/kg, from Everbatt2023
    SP.ins[1].price = 0.08 * stream.acid.imass['H2SO4'] / stream.acid.F_mass
    
    # price of H2O2 is 1.46$/kg, from Everbatt2023
    SP.ins[2].price = 1.46 * stream.oxidant.imass['H2O2'] / stream.oxidant.F_mass
    
    HD.ins[1].price = bst.HeatUtility.heating_agents[-1].regeneration_price/MW_CH4
    
    # TODO: this only applies to the sludge management perspective; the value needs to be updated: the price of FePO4 is 15000 rmb/ton, means 2.14 $/kg
    SI.outs[0].price = 2.14
    
    SI.ins[1].price = bst.HeatUtility.heating_agents[-1].regeneration_price/MW_CH4
    
    # the landfill cost (tipping fee) can be found in https://www.wasteoptima.com/blog/eref-2024-landfill-tipping-fees?utm_source with the $62.28 per ton
    FC.outs[1].price = - 62.28 / _ton_to_kg
    
    # TODO: food waste now is emission-free, but it could provide credits as well
    # TODO: now, electricity and steam use default prices
    
    # TODO: food_sludge_ratio and HRT also affect the concentration of VFAs, and therefore their price
    # TODO: this is a guess from Xuan: the outs is rich in VFAs at 4000-5000 mg-VFAs/L (COD = 1 mg VFA)= 4.45-5.5 kg-COD/m3 =5 kg/m3 value=5*0.4=2$/m3 =0.002 $/kg
    
    def get_precipitation_supernatant_price(food_sludge_ratio, HRT):
        VFA_conc_dict = {
            (0, 0):229.836,    #mg/L
            (0, 12):229.836,
            (0, 24):229.836,
            (0, 36):229.836,
            (0, 48):229.836,
            (0, 60):229.836,
            (0, 72):229.836,
            (0, 84):315.014,
            (0, 96):324.748,
            (0, 108):474.744,
            (0, 120):503.468,
            (0, 132):574.002,
            (1/3, 0):229.836,    #mg/L
            (1/3, 12):318.997,
            (1/3, 24):891.924,
            (1/3, 36):930.369,
            (1/3, 48):942.929,
            (1/3, 60):1019.904,
            (1/3, 72):1030.104,
            (1/3, 84):1058.436,
            (1/3, 96):1065.347,
            (1/3, 108):1235.25,
            (1/3, 120):1295.057,
            (1/3, 132):1486.757,
            (2/3, 0):229.836,
            (2/3, 12):619.055,
            (2/3, 24):1342.854,
            (2/3, 36):1975.393,
            (2/3, 48):2195.785,
            (2/3, 60):2540.621,
            (2/3, 72):2677.336,
            (2/3, 84):2740.366,
            (2/3, 96):2754.408,
            (2/3, 108):2835.845,
            (2/3, 120):2882.61,
            (2/3, 132):3123.319,
            (3/3, 0):229.836,    #mg/L
            (3/3, 12):755.134,
            (3/3, 24):1534.219,
            (3/3, 36):2006.594,
            (3/3, 48):2624.859,
            (3/3, 60):3049.32,
            (3/3, 72):3413.96,
            (3/3, 84):3716.996,
            (3/3, 96):3763.005,
            (3/3, 108):3953.696,
            (3/3, 120):4114.03,
            (3/3, 132):4500,
            (4/3, 132):4650,
          }
            
        vfa_conc = VFA_conc_dict[(food_sludge_ratio),HRT]  # mg/L
        price = 0.002 * (vfa_conc / 4500)
        
        return price
    
    ratio = AF.food_sludge_ratio
    HRT = AF.HRT
    PC.outs[0].price = get_precipitation_supernatant_price(ratio, HRT)
 
    # PC.outs[0].price = 0.002
 
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
        GlobalWarming=-0.00582207178953914/30*sludge_CI_credit
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
    
    # the CI of FePO4: https://www.sciencedirect.com/science/article/pii/S1383586624004258#s0085  2-3.5-6 kg CO2-eq/kg FePO4
    qs.StreamImpactItem(
    ID='FePO4_credit',
    linked_stream=stream.product,
    GlobalWarming = -4
    )
    # kg CO2/kg FePO4
    
    qs.StreamImpactItem(
    ID='VFA_credit',
    linked_stream=stream.precipitation_supernatant,
    GlobalWarming = -0.002
    )
    
    # TODO: add labor costs
    # labor cost=1 FTE*8000h/year*25$-wage/hour=200,000 $/year
    # wage:U.S. Bureau of Labor Statistics. Producer Price Index–Metals and Metal Products: Mid–Atlantic Information Office: U.S. Bureau of Labor Statistics. Available online: https://www.bls.gov/regions/mid-atlantic/data/producerpriceindexmetals_us_table.htm (accessed on 28 June 2022).
    # FTE amount: https://www.mheducation.com/highered/product/wastewater-engineering-treatment-resource-recovery-metcalf-eddy/M9780073401188.html (highly automated treatment processes require minimal operator attention)
    # other refer-1: the labor cost is €462.25 per month:https://www.sciencedirect.com/science/article/pii/S0959652617305504
    # other refer-2: operation labor in need (number of operator·hours per year):https://pubs.acs.org/doi/10.1021/acssuschemeng.0c05189
    create_tea(sys,
               duration=(2023, 2023+lifetime),
               IRR_value=0.03,
               income_tax_value=0.3,
               finance_interest_value=0.03,
               labor_cost_value=200000) 
    
    qs.LCA(
        system=sys, lifetime=lifetime, lifetime_unit='yr',
        Electricity=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*30,
        Steam=lambda:abs(sum(i.duty/1000*operation_hours*lifetime for i in sys.heat_utilities if 'steam' in i.ID)),
    )
    
    return sys