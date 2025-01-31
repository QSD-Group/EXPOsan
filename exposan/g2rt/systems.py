#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from qsdsan import (
    Flowsheet, 
    main_flowsheet,
    WasteStream,
    sanunits as qsu,
    ImpactItem,
    System, TEA, LCA,
    )
from exposan.g2rt import _sanunits as su
from exposan.surt import _sanunits as ssu
from exposan.g2rt._components import create_components
from exposan.g2rt import (
    create_components,
    _load_lca_data,
    _load_components,
    update_resource_recovery_settings,
    )
import biosteam as bst
from qsdsan.utils import clear_lca_registries
from exposan.utils import (
    add_fugitive_items,
    get_decay_k,
    )

__all__ = ('create_system',)

#%%

# from exposan.g2rt import (
#     _load_components,
#     _load_lca_data,
#     discount_rate,
#     get_toilet_user,
#     max_CH4_emission,
#     operator_daily_wage,
#     ppl,
#     price_dct,
#     )

# __all__= ('create_system',)

#%%
# =============================================================================
# Universal units and functions
# =============================================================================

def batch_create_streams(prefix, phases=('liq', 'sol')):
    item = ImpactItem.get_item('CH4_item').copy(f'{prefix}_CH4_item', set_as_source=True)
    WasteStream('CH4', phase='g', stream_impact_item=item)

    item = ImpactItem.get_item('N2O_item').copy(f'{prefix}_N2O_item', set_as_source=True)
    WasteStream('N2O', phase='g', stream_impact_item=item)
    
    item = ImpactItem.get_item('NO_item').copy(f'{prefix}_NO_item', set_as_source=True)
    WasteStream('NO', phase='g', stream_impact_item=item)
    
    item = ImpactItem.get_item('SO2_item').copy(f'{prefix}_SO2_item', set_as_source=True)
    WasteStream('SO2', phase='g', stream_impact_item=item)
    
    item = ImpactItem.get_item('NH3_item').copy(f'{prefix}_NH3_item', set_as_source=True)
    WasteStream('NH3', phase='g', stream_impact_item=item)
    
    price_dct = update_resource_recovery_settings()[0]
    item = ImpactItem.get_item('H2O_item').copy(f'{prefix}_H2O_item', set_as_source=True)
    WasteStream('H2O', phase='l', price=price_dct.get('H2O'), stream_impact_item=item)

    for nutrient in ('N', 'P', 'K'):
        for phase in phases:
            original = ImpactItem.get_item(f'{nutrient}_item')
            new = original.copy(f'{phase}_{nutrient}_item', set_as_source=True)
            WasteStream(f'{phase}_{nutrient}', phase='l',
                        price=price_dct[nutrient], stream_impact_item=new)

    def create_stream_with_impact_item(stream_ID='', item_ID='', dct_key=''):
        item_ID = item_ID or stream_ID+'_item'
        dct_key = dct_key or stream_ID
        item = ImpactItem.get_item(item_ID).copy(f'{prefix}_{item_ID}', set_as_source=True)
        WasteStream(f'{stream_ID}', phase='l',
                    price=price_dct.get(dct_key) or 0., stream_impact_item=item)


    # create_stream_with_impact_item(stream_ID='reverse_osmosis_treated',
    #                                dct_key = 'H2O')
    # create_stream_with_impact_item('H2O')
    WasteStream('H2O_vapor', phase='g')
    WasteStream('H2O_vapor1', phase='g')
    WasteStream('air', phase = 'g')
    # create_stream_with_impact_item(stream_ID='NH3_gas', dct_key='NH3',item_ID = 'NH3_item')
    
    WasteStream('NH3_gas', phase='g')
    WasteStream('NH3_gas1', phase='g')

# %%

# =========================================================================
# Unit parameters
# =========================================================================
default_ppl = 6
default_lifetime = 10 #10 years
discount_rate = 0.05

# Time take for full degradation, [yr]
tau_deg = 2
# Log reduction at full degradation
log_deg = 3

max_CH4_emission = 0.25

emptying_fee = 0.15

# Nutrient loss during application
app_loss = dict.fromkeys(('NH3', 'NonNH3', 'P', 'K', 'Mg', 'Ca'), 0.02)
app_loss['NH3'] = 0.05


# %%
# =============================================================================
# System A
# Volume reduction toilet based on
# https://patentimages.storage.googleapis.com/73/c1/74/aa1f6a78b89957/WO2023288326A1.pdf

# =============================================================================

def create_systemA(flowsheet=None, ppl=default_ppl, lifetime=default_lifetime):
    # TODO: Set flowsheet to avoid stream replacement warnings
    flowsheet = flowsheet or main_flowsheet
    batch_create_streams('A')
    streamA = flowsheet.stream

    #### Human Inputs ####
    A1 = su.Excretion('A1', outs=('urine','feces'),)
    
    # recycle_fw = WasteStream('reverse_osmosis_treated')
    # A2 = su.SURT('A2',
    #              ins= (A1-0, A1-1, 'toilet_paper', 'flushing_water', recycle_fw),
    #              outs = ('mixed_waste','A2_CH4','A2_N2O'),
    #              N_user = 6, N_tot_user=default_ppl, lifetime = 10, 
    #              if_include_front_end=True, if_toilet_paper=True,
    #              if_flushing=True, if_cleansing=False,
    #              if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
    #              CAPEX=500*max(1, default_ppl/100), OPEX_over_CAPEX=0.06,
    #              decay_k_COD=get_decay_k(),
    #              decay_k_N=get_decay_k(),
    #              max_CH4_emission=max_CH4_emission)
    mixer = su.FWMixer('Mixer',
                      ins = ('tap_water', streamA['H2O']),
                      outs = ('flushing_water'),
                      N_tot_user = ppl
                      )
    A2 = su.SURT('A2',
                  ins= (A1-0, A1-1, 'toilet_paper'),
                  outs = ('mixed_waste','A2_CH4','A2_N2O'),
                  N_user = 6, N_tot_user=ppl, lifetime = lifetime, 
                  if_include_front_end=True, if_toilet_paper=True,
                  if_flushing=True, if_cleansing=False,
                  if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                  CAPEX=500*max(1, default_ppl/100), OPEX_over_CAPEX=0.06,
                  decay_k_COD=get_decay_k(),
                  decay_k_N=get_decay_k(),
                  max_CH4_emission=max_CH4_emission)
    
    A3 = su.G2RTSolidsSeparation('A3',
                                ins = (A2-0,mixer-0),
                                outs = ('A3_liquid','A3_solid'),
                                )
    
    # make a recycle loop
    # recycle_uf = WasteStream('ultrafiltration_reject')
    # A4 = su.G2RTBeltSeparation('A4', 
    #                       ins = (A3-0,A3-1,recycle_uf),
    #                       outs = ('A4_liquid','A4_solid'),
    #                       )
    UFmixer = su.UFMixer('UFMixer',
                         ins = (A3-0,'ultrafiltration_reject'),
                         outs = 'UFmixed_liquid'
                         )
    A4 = su.G2RTBeltSeparation('A4', 
                          ins = (UFmixer-0,A3-1),
                          outs = ('A4_liquid','A4_solid'),
                          )
    A12 = su.G2RTSolidsTank('A12',
                             ins = A4-1,
                             outs = 'A12_solids'
                             )
    
    A13 = su.G2RTLiquidsTank('A13',
                             ins = (A4-0,'filter_press_liquid','concentrator_liquid'),
                             outs = 'A12_liquids'
                             )
    # A5 = su.G2RTUltrafiltration('A5',
    #                                   ins = A13-0,
    #                                   outs = ('A5_treated',recycle_uf)
    #                                   )
    A5 = su.G2RTUltrafiltration('A5',
                                      ins = A13-0,
                                      outs = ('A5_treated',1-UFmixer)
                                      )
    
    # A6 = su.G2RTReverseOsmosis('A6', 
    #                         ins = A5-0,
    #                         outs = (recycle_fw,'A6_brine')
    #                         )
    A6 = su.G2RTReverseOsmosis('A6', 
                            ins = A5-0,
                            outs = (1-mixer,'A6_brine')
                            )
    
    # item = ImpactItem.get_item('H2O_item').copy('A_H2O_item', set_as_source=True)
    # A6.outs[0].stream_impact_item = item

    A7 = su.VRConcentrator('A7',
                           ins = A6-1,
                           outs = ('A7_consensed_waste',2-A13,'A7_N2O','A7_CH4',
                                   'A7_NH3_gas',streamA['H2O_vapor1'])
                           )

    A8 = su.G2RThomogenizer('A8', 
                            ins = A12-0, 
                            outs = 'A8_grinded')
    
    A9 = su.VRpasteurization('A9', 
                                 ins = A8-0,
                                 outs = 'A9_pasteurized'
                                 )
    A10 = su.VolumeReductionFilterPress('A10',
                                        ins = A9-0,
                                        outs = (1-A13,'A10_pressed_solid_cake')
                                        )
    A11 = su.VRdryingtunnel('A11',
                           ins = (A7-0,A10-1),
                           outs = ('A11_solid_cakes','A11_N2O', 'A11_CH4',
                                   'A11_NH3_gas',streamA['H2O_vapor'])
                           )
    A14 = su.VolumeReductionCombustor('A14',
                                      ins = (A11-0,'Wood_pellets','Air'),
                                      outs = ("Ash","hot_gas","A14_CH4", "A14_N2O", 
                                              streamA['NO'],streamA['SO2'],'A14_NH3_gas'),
                                      if_sludge_service = True, lifetime = lifetime
                                      )
    
    # CH4 emissions from MURT, concentrator, drying tunnel, and combustor
    A15 = qsu.Mixer('A15', ins=(A2-1, A7-3, A11-2, A14-2), outs=streamA['CH4'])
    A15.add_specification(lambda: add_fugitive_items(A15, 'CH4_item'))
    A15.line = 'fugitive CH4 mixer'
    # N2O emissions from MURT, concentrator, drying tunnel, and combustor
    A16 = qsu.Mixer('A16', ins=(A2-2, A7-2, A11-1, A14-3), outs=streamA['N2O'])
    A16.add_specification(lambda: add_fugitive_items(A16, 'N2O_item'))
    A16.line = 'fugitive N2O mixer'
    # NH3 emissions from concentrator, drying tunnel, and combustor
    A17 = qsu.Mixer('A17', ins=(A7-4, A11-3, A14-6), outs=streamA['NH3'])
    A17.add_specification(lambda: add_fugitive_items(A17, 'NH3_item'))
    A17.line = 'fugitive NH3 mixer'
    
    # N, P, K split for recovery
    A18 = qsu.ComponentSplitter('A18', ins=A14-0,
                               outs=(streamA['sol_N'], streamA['sol_P'], streamA['sol_K'],
                                     'A_sol_non_fertilizers'),
                               split_keys=('N', 'P', 'K'))
    A19 = su.G2RTControls('A19', ins = A11-4,
                          outs= 'A19_out')
    A20 = su.G2RTHousing('A20', ins = A19-0, outs = "hot_gas")
    qsu.HeatExchangerNetwork('HXN')
    
    sysA_1 = System('sysA_1',
                    path = (A1,mixer,A2,A3,UFmixer,A4,A12,A13,A5),
                    recycle = A5-1
                    )
    # sysA = System('sysA',
    #               path = (sysA_1,A6,A7,A8,A9,A10,A11,A14,A15,A16),
    #               recycle = A10-0
    #               )
    sysA_2 = System('sysA_2',
                  path = (sysA_1,A6),
                  recycle = A6-0
                  )
    sysA_3 = System('sysA_3',
                    path = (sysA_2,A7,A8,A9),
                    recycle = A7-1
                    )
    sysA = System('sysA',
                  path=(sysA_3,A10,A11,A14,A15,A16,A17,A18,A19,A20),
                  recycle = A10-0
                  )
    
    teaA = TEA(system=sysA, discount_rate=discount_rate,
               start_year=2020, lifetime=lifetime, uptime_ratio=1,
               lang_factor=None, annual_maintenance=0,
               annual_labor=0)
    
    # batch_create_streams('A')
    get_powerA = sum([u.power_utility.rate for u in sysA.units]) * (24 * 365 * teaA.lifetime)
    
    LCA(system=sysA, lifetime=lifetime, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerA)
    # sysA.simulate()
    return sysA
    

# %%
# =============================================================================
# System B
# micro Supercritical Water Oxidation toilet based on
# https://patentimages.storage.googleapis.com/57/6a/81/72a168a92be44c/WO2023288331A1.pdf
# =============================================================================
def create_systemB(flowsheet=None, ppl=default_ppl, lifetime= default_lifetime):
    # TODO: Set flowsheet to avoid stream replacement warnings
    flowsheet = flowsheet or main_flowsheet
    batch_create_streams('B')
    streamB = flowsheet.stream

    #### Human Inputs ####
    B1 = su.Excretion('B1', outs=('urine','feces'),)

    mixer = su.FWMixer('Mixer',
                      ins = ('tap_water', streamB['H2O']),
                      outs = ('flushing_water'),
                      N_tot_user = ppl
                      )
    B2 = su.SURT('B2',
                  ins= (B1-0, B1-1, 'toilet_paper'),
                  outs = ('mixed_waste','B2_CH4','B2_N2O'),
                  N_user = 6, N_tot_user=ppl, lifetime = lifetime, 
                  if_include_front_end=True, if_toilet_paper=True,
                  if_flushing=True, if_cleansing=False,
                  if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                  CAPEX=500*max(1, ppl/100), OPEX_over_CAPEX=0.06,
                  decay_k_COD=get_decay_k(),
                  decay_k_N=get_decay_k(),
                  max_CH4_emission=max_CH4_emission)
    
    B3 = su.G2RTSolidsSeparation('B3',
                                ins = (B2-0,mixer-0),
                                outs = ('B3_liquid','B3_solid'),
                                )
    
    # make a recycle loop
    # recycle_uf = WasteStream('ultrafiltration_reject')
    # A4 = su.G2RTBeltSeparation('A4', 
    #                       ins = (A3-0,A3-1,recycle_uf),
    #                       outs = ('A4_liquid','A4_solid'),
    #                       )
    UFmixer = su.UFMixer('UFMixer',
                         ins = (B3-0,'ultrafiltration_reject'),
                         outs = 'UFmixed_liquid'
                         )
    B4 = su.G2RTBeltSeparation('B4', 
                          ins = (UFmixer-0,B3-1),
                          outs = ('B4_liquid','B4_solid'),
                          )
    B5 = su.G2RTSolidsTank('B5',
                             ins = B4-1,
                             outs = 'B5_solids'
                             )
    
    B6 = su.G2RThomogenizer('B6', 
                            ins = B5-0, 
                            outs = 'B6_grinded')
    
    B7 = su.G2RTLiquidsTank('B7',
                             ins = (B4-0,'concentrator_liquid'),
                             outs = 'B7_liquids'
                             )

    B8 = su.G2RTUltrafiltration('B8', 
                                ins = B7-0, 
                                outs = ('B8_treated',1-UFmixer))
    
    # A6 = su.G2RTReverseOsmosis('A6', 
    #                         ins = A5-0,
    #                         outs = (recycle_fw,'A6_brine')
    #                         )
    B9 = su.G2RTReverseOsmosis('B9', 
                            ins = B8-0,
                            outs = (1-mixer,'B9_brine') #TODO: sys.flowsheet.unit.A6.outs[0].stream_impact_item
                            )
    
    # item = ImpactItem.get_item('H2O_item').copy('A_H2O_item', set_as_source=True)
    # A6.outs[0].stream_impact_item = item
    
    B10 = su.mSCWOGasModule('B10',
                           ins = streamB['air'],
                           outs = 'B10_compressed_air',
                           )
    
    B11 = su.mSCWOReactorModule('B11',
                                ins = (B10-0,B6-0),
                                outs = ('B11_gas','B11_N2O','B11_liquid_effluent','B11_ash_waste'),
                                )

    B12 = su.mSCWOConcentratorModule('B12',
                           ins = (B9-1,B11-2,B11-3),
                           outs = ('B12_consensed_waste',1-B7,'B12_N2O','B12_CH4',
                                   'B12_NH3_gas',streamB['H2O_vapor1'])
                           )

    # B9 = su.VRpasteurization('B9', 
    #                              ins = B8-0,
    #                              outs = 'B9_pasteurized'
    #                              )
    # B10 = su.VolumeReductionFilterPress('B10',
    #                                     ins = B9-0,
    #                                     outs = (1-B13,'B10_pressed_solid_cake')
    #                                     )
    
    B13 = su.VRdryingtunnel('B13',
                           ins = B12-0,
                           outs = ('B13_solid_cakes','B13_N2O', 'B13_CH4',
                                   'B13_NH3_gas',streamB['H2O_vapor'])
                           )
    
    # B14 = su.VolumeReductionCombustor('B14',
    #                                   ins = (B11-0,'Wood_pellets'),
    #                                   outs = ("Ash","hot_gas","A14_CH4", "A14_N2O", 
    #                                           streamA['NO'],streamA['SO2'],'A14_NH3_gas'),
    #                                   if_sludge_service = True
    #                                   )
    
    # CH4 emissions from MURT, concentrator, drying tunnel
    B14 = qsu.Mixer('B14', ins=(B2-1, B12-3, B13-2), outs=streamB['CH4'])
    B14.add_specification(lambda: add_fugitive_items(B14, 'CH4_item'))
    B14.line = 'fugitive CH4 mixer'
    # N2O emissions from MURT, mSCWO reactor, concentrator, drying tunnel
    B15 = qsu.Mixer('B15', ins=(B2-2, B11-1, B12-2, B13-1), outs=streamB['N2O'])
    B15.add_specification(lambda: add_fugitive_items(B15, 'N2O_item'))
    B15.line = 'fugitive N2O mixer'
    
    # NH3 emissions from concentrator, drying tunnel
    B16 = qsu.Mixer('B16', ins=(B12-4, B13-3), outs=streamB['NH3'])
    B16.add_specification(lambda: add_fugitive_items(B16, 'NH3_item'))
    B16.line = 'fugitive NH3 mixer'
    
    # N, P, K split for recovery
    B17 = qsu.ComponentSplitter('B17', ins=B13-0,
                               outs=(streamB['sol_N'], streamB['sol_P'], streamB['sol_K'],
                                     'B_sol_non_fertilizers'),
                               split_keys=('N', 'P', 'K'))
    B18 = su.G2RTControls('B18', ins = B13-4,
                          outs= 'B18_out')
    B19 = su.G2RTHousing('B19', ins = B18-0, outs = "hot_gas")
    
    sysB_1 = System('sysB_1',
                    path = (B1,mixer,B2,B3,UFmixer,B4,B5,B6,B7,B8),
                    recycle = B8-1
                    )
    sysB_2 = System('sysB_2',
                  path = (sysB_1,B9,B10,B11),
                  recycle = B9-0
                  )

    sysB = System('sysB',
                  path=(sysB_2,B12,B13,B14,B15,B16,B17,B18,B19),
                  recycle = B12-1
                  )
    
    teaB = TEA(system=sysB, discount_rate=discount_rate,
               start_year=2020, lifetime=lifetime, uptime_ratio=1,
               lang_factor=None, annual_maintenance=0,
               annual_labor=0)
    
    # batch_create_streams('A')
    get_powerB = sum([u.power_utility.rate for u in sysB.units]) * (24 * 365 * teaB.lifetime)
    
    LCA(system=sysB, lifetime=lifetime, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerB)
    # sysA.simulate()
    return sysB
#%%
# =============================================================================
# System C
# Volume reduction toilet with combustion energy recovery, based on
# https://patentimages.storage.googleapis.com/73/c1/74/aa1f6a78b89957/WO2023288326A1.pdf
# =============================================================================

def create_systemC(flowsheet=None, ppl=default_ppl, lifetime= default_lifetime):
    # TODO: Set flowsheet to avoid stream replacement warnings
    flowsheet = flowsheet or main_flowsheet
    batch_create_streams('C')
    streamA = flowsheet.stream

    #### Human Inputs ####
    A1 = su.Excretion('A1', outs=('urine','feces'),)
    
    # recycle_fw = WasteStream('reverse_osmosis_treated')
    # A2 = su.SURT('A2',
    #              ins= (A1-0, A1-1, 'toilet_paper', 'flushing_water', recycle_fw),
    #              outs = ('mixed_waste','A2_CH4','A2_N2O'),
    #              N_user = 6, N_tot_user=default_ppl, lifetime = 10, 
    #              if_include_front_end=True, if_toilet_paper=True,
    #              if_flushing=True, if_cleansing=False,
    #              if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
    #              CAPEX=500*max(1, default_ppl/100), OPEX_over_CAPEX=0.06,
    #              decay_k_COD=get_decay_k(),
    #              decay_k_N=get_decay_k(),
    #              max_CH4_emission=max_CH4_emission)
    mixer = su.FWMixer('Mixer',
                      ins = ('tap_water', streamA['H2O']),
                      outs = ('flushing_water'),
                      N_tot_user = 6
                      )
    A2 = su.SURT('A2',
                  ins= (A1-0, A1-1, 'toilet_paper'),
                  outs = ('mixed_waste','A2_CH4','A2_N2O'),
                  N_user = 6, N_tot_user=ppl, lifetime = lifetime, 
                  if_include_front_end=True, if_toilet_paper=True,
                  if_flushing=True, if_cleansing=False,
                  if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                  CAPEX=500*max(1, ppl/100), OPEX_over_CAPEX=0.06,
                  decay_k_COD=get_decay_k(),
                  decay_k_N=get_decay_k(),
                  max_CH4_emission=max_CH4_emission)
    
    A3 = su.G2RTSolidsSeparation('A3',
                                ins = (A2-0,mixer-0),
                                outs = ('A3_liquid','A3_solid'),
                                )
    
    # make a recycle loop
    # recycle_uf = WasteStream('ultrafiltration_reject')
    # A4 = su.G2RTBeltSeparation('A4', 
    #                       ins = (A3-0,A3-1,recycle_uf),
    #                       outs = ('A4_liquid','A4_solid'),
    #                       )
    UFmixer = su.UFMixer('UFMixer',
                         ins = (A3-0,'ultrafiltration_reject'),
                         outs = 'UFmixed_liquid'
                         )
    A4 = su.G2RTBeltSeparation('A4', 
                          ins = (UFmixer-0,A3-1),
                          outs = ('A4_liquid','A4_solid'),
                          )
    A12 = su.G2RTSolidsTank('A12',
                             ins = A4-1,
                             outs = 'A12_solids'
                             )
    
    A13 = su.G2RTLiquidsTank('A13',
                             ins = (A4-0,'filter_press_liquid','concentrator_liquid'),
                             outs = 'A12_liquids'
                             )
    # A5 = su.G2RTUltrafiltration('A5',
    #                                   ins = A13-0,
    #                                   outs = ('A5_treated',recycle_uf)
    #                                   )
    A5 = su.G2RTUltrafiltration('A5',
                                      ins = A13-0,
                                      outs = ('A5_treated',1-UFmixer)
                                      )
    
    # A6 = su.G2RTReverseOsmosis('A6', 
    #                         ins = A5-0,
    #                         outs = (recycle_fw,'A6_brine')
    #                         )
    A6 = su.G2RTReverseOsmosis('A6', 
                            ins = A5-0,
                            outs = (1-mixer,'A6_brine')
                            )
    
    # item = ImpactItem.get_item('H2O_item').copy('A_H2O_item', set_as_source=True)
    # A6.outs[0].stream_impact_item = item

    A7 = su.VRConcentrator('A7',
                           ins = A6-1,
                           outs = ('A7_consensed_waste',2-A13,'A7_N2O','A7_CH4',
                                   'A7_NH3_gas',streamA['H2O_vapor1'])
                           )

    A8 = su.G2RThomogenizer('A8', 
                            ins = A12-0, 
                            outs = 'A8_grinded')
    
    A9 = su.VRpasteurization('A9', 
                                 ins = A8-0,
                                 outs = 'A9_pasteurized'
                                 )
    A10 = su.VolumeReductionFilterPress('A10',
                                        ins = A9-0,
                                        outs = (1-A13,'A10_pressed_solid_cake')
                                        )
    A11 = su.VRdryingtunnel('A11',
                           ins = (A7-0,A10-1),
                           outs = ('A11_solid_cakes','A11_N2O', 'A11_CH4',
                                   'A11_NH3_gas',streamA['H2O_vapor'])
                           )
    A14 = ssu.BiomassCombustion('A14',
                                      ins = (A11-0,'Wood_pellets','Air'),
                                      outs = ("Ash","hot_gas","A14_CH4", "A14_N2O", 
                                              streamA['NO'],streamA['SO2'],'A14_NH3_gas'),
                                      if_sludge_service = False, lifetime = lifetime
                                      )
    
    # CH4 emissions from MURT, concentrator, drying tunnel, and combustor
    A15 = qsu.Mixer('A15', ins=(A2-1, A7-3, A11-2, A14-2), outs=streamA['CH4'])
    A15.add_specification(lambda: add_fugitive_items(A15, 'CH4_item'))
    A15.line = 'fugitive CH4 mixer'
    # N2O emissions from MURT, concentrator, drying tunnel, and combustor
    A16 = qsu.Mixer('A16', ins=(A2-2, A7-2, A11-1, A14-3), outs=streamA['N2O'])
    A16.add_specification(lambda: add_fugitive_items(A16, 'N2O_item'))
    A16.line = 'fugitive N2O mixer'
    # NH3 emissions from concentrator, drying tunnel, and combustor
    A17 = qsu.Mixer('A17', ins=(A7-4, A11-3, A14-6), outs=streamA['NH3'])
    A17.add_specification(lambda: add_fugitive_items(A17, 'NH3_item'))
    A17.line = 'fugitive NH3 mixer'
    
    # N, P, K split for recovery
    A18 = qsu.ComponentSplitter('A18', ins=A14-0,
                               outs=(streamA['sol_N'], streamA['sol_P'], streamA['sol_K'],
                                     'A_sol_non_fertilizers'),
                               split_keys=('N', 'P', 'K'))
    A19 = su.G2RTControls('A19', ins = A11-4,
                          outs= 'A19_out')
    A20 = su.G2RTHousing('A20', ins = A19-0, outs = "hot_gas")
    qsu.HeatExchangerNetwork('HXN')
    
    sysA_1 = System('sysA_1',
                    path = (A1,mixer,A2,A3,UFmixer,A4,A12,A13,A5),
                    recycle = A5-1
                    )
    # sysA = System('sysA',
    #               path = (sysA_1,A6,A7,A8,A9,A10,A11,A14,A15,A16),
    #               recycle = A10-0
    #               )
    sysA_2 = System('sysA_2',
                  path = (sysA_1,A6),
                  recycle = A6-0
                  )
    sysA_3 = System('sysA_3',
                    path = (sysA_2,A7,A8,A9),
                    recycle = A7-1
                    )
    sysA = System('sysA',
                  path=(sysA_3,A10,A11,A14,A15,A16,A17,A18,A19,A20),
                  recycle = A10-0
                  )
    
    
    teaA = TEA(system=sysA, discount_rate=discount_rate,
               start_year=2020, lifetime=lifetime, uptime_ratio=1,
               lang_factor=None, annual_maintenance=0,
               annual_labor=0)
    
    get_powerA = sum([u.power_utility.rate for u in sysA.units]) * (24 * 365 * teaA.lifetime)
    
    LCA(system=sysA, lifetime=lifetime, lifetime_unit='yr', uptime_ratio=1, e_item=get_powerA)
    # sysA.simulate()
    return sysA
    



#%%

# =============================================================================
# create system
# =============================================================================

def create_system(system_ID='A', flowsheet=None, ppl=default_ppl, lifetime = default_lifetime):
    ID = system_ID.lower().lstrip('sys').upper()  # so that it'll work for "sysA"/"A"
    reload_lca = False

    # Set flowsheet to avoid stream replacement warnings
    if flowsheet is None:
        flowsheet_ID = f'g2rt{ID}' #TODO: what is this ID based on?
        if hasattr(main_flowsheet.flowsheet, flowsheet_ID): # clear flowsheet
            getattr(main_flowsheet.flowsheet, flowsheet_ID).clear()
            clear_lca_registries()
            reload_lca = True
        flowsheet = Flowsheet(flowsheet_ID)
        main_flowsheet.set_flowsheet(flowsheet)

    _load_components()
    _load_lca_data(reload_lca)

    if system_ID == 'A': f = create_systemA
    elif system_ID == 'B': f = create_systemB
    elif system_ID == 'C': f = create_systemC
    else: raise ValueError(f'`system_ID` can only be "A" or "B" or "C", not "{ID}".')

    try: system = f(flowsheet, ppl=ppl, lifetime = lifetime)
    except:
        _load_components(reload=True)
        system = f(flowsheet, ppl=ppl, lifetime = lifetime)

    return system




    
    
    

