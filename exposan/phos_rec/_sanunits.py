# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Xuan Wang <easonbiubiu99@gmail.com>
    
    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

import qsdsan as qs
from qsdsan import SanUnit
from qsdsan.sanunits import HXutility
from qsdsan.equipments import VerticalMixer
from biosteam import Stream
from biosteam.units.decorators import cost
from thermosteam.reaction import ParallelReaction

_C_to_K = 273.15

__all__ = (
    'AcidogenicFermenter',
    'SelectivePrecipitation',
    'HeatDrying',
    'Sintering'
)

# =============================================================================
# AcidogenicFermenter
# =============================================================================

# TODO: may need to update F_BM for verticle_mixer
class AcidogenicFermenter(SanUnit):
    '''
    Fermentation of sludge and food waste.
    
    Parameters
    ----------
    ins : iterable
        fe_sludge, food_waste.
    outs : iterable
        gas, fermentate.
    food_sludge_ratio : float
        Mass ratio of organics bewteen food waste and sludge, [-].
    food_waste_moisture : float
        Mositure content of food waste, [-].
    org_to_gas : float
        Mass ratio of total organics converted to gas, [-].
    org_to_vfa : float
        Mass ratio of total organics converted to volatile fatty acids (VFAs), [-].
    org_to_ethanol : float
        Mass ratio of total organics converted to ethanol, [-].
    org_to_residue : float
        Mass ratio of total organics in the residue after acidogenic fermentation, [-].
    VFA_ratio : float
        Mass fractions of individual VFAs, [-].
    fe_reduction : float
        Fe3+ reduction ratio during acidogenic fermentation, [-].
    # TODO: confirm the unit of HRT
    HRT : float
        Hydraulic retention time, [d].
    # TODO: what is the temperature for slurry; is cooling necessary?
    T : float
        Required temperature for fermentation, [K].
    '''
    _N_ins = 2
    _N_outs = 2
    
    # TODO: may need to update F_BM for heat_exchanger
    auxiliary_unit_names=('heat_exchanger',)
    
    _units = {
        'Flow rate': 'm3/d',
        'HRT': 'd',
        'Reactor volume': 'm3',
        'Required mixing power': 'kW'
    }
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 food_sludge_ratio=1,
                 food_waste_moisture=0.2,
                 org_to_gas=0.05,
                 org_to_vfa=0.65, 
                 org_to_ethanol=0.02,
                 org_to_residue=0.25,
                 VFA_ratio={'Ac': 0.5, 'Pr': 0.24, 'Bu': 0.23, 'Va': 0.02, 'Lac': 0.01},
                 fe_reduction=0.98,
                 HRT=3,
                 T=37+_C_to_K): 
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.food_sludge_ratio = food_sludge_ratio
        self.food_waste_moisture = food_waste_moisture
        self.org_to_gas = org_to_gas
        self.org_to_vfa = org_to_vfa
        self.org_to_ethanol = org_to_ethanol
        self.org_to_residue = org_to_residue
        self.VFA_ratio = VFA_ratio
        self.fe_reduction = fe_reduction
        self.HRT = HRT
        self.T = T
        hx_in = Stream(ID=f'{ID}_hx_in')
        hx_out = Stream(ID=f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out)
        self.verticle_mixer = VerticalMixer(ID='verticle_mixer', linked_unit=self)
        self.equipments = (self.verticle_mixer,)
        
        
    def _run(self):
        fe_sludge, food_waste = self.ins
        gas, fermentate = self.outs
        
        fe_sludge.phase = 'l'
        food_waste.phase = 'l'
        gas.phase = 'g'
        fermentate.phase = 'l'
        
        if self.food_sludge_ratio not in [0, 1/3, 2/3, 1, 4/3]:
            raise RuntimeError('food_sludge_ratio must be one of the follow: 0, 1/3, 2/3, 1, 4/3.')
        
        food_waste.imass['Org'] = fe_sludge.imass['Org']*self.food_sludge_ratio
        food_waste.imass['H2O'] = food_waste.imass['Org']/(1-self.food_waste_moisture) * self.food_waste_moisture
        
        fermentate.mix_from(self.ins)
        
        org_total = fermentate.imass['Org']
        
        # gas production
        gas.imass['CO2'] = org_total*self.org_to_gas
        
        # VFAs and inorganics in fermentate
        vfa_mass = org_total * self.org_to_vfa
        
        fermentate.imass['Acetic_acid'] = vfa_mass * self.VFA_ratio['Ac']
        fermentate.imass['Propionic_acid'] = vfa_mass * self.VFA_ratio['Pr']
        fermentate.imass['Butyric_acid'] = vfa_mass * self.VFA_ratio['Bu']
        fermentate.imass['Valeric_acid'] = vfa_mass * self.VFA_ratio['Va']
        fermentate.imass['Lactic_acid'] = vfa_mass * self.VFA_ratio['Lac']
        
        fermentate.imass['Ethanol'] = org_total*self.org_to_ethanol
        
        fermentate.imass['Residue'] = org_total*self.org_to_residue
        
        if self.org_to_gas + self.org_to_vfa + self.org_to_ethanol + self.org_to_residue > 1:
            raise RuntimeError('org cannot be balanced')
        
        fermentate.imass['Org'] *= (1 - self.org_to_gas - self.org_to_vfa - self.org_to_ethanol - self.org_to_residue)

        fermentate.imass['Fe2'] = fermentate.imass['Fe3'] * self.fe_reduction
        fermentate.imass['Fe3'] -= fermentate.imass['Fe2']
        
        metal_P_release = {
            0: {'metal': 2.124263, 'P': 11.0693},
            1/3: {'metal': 21.13664, 'P': 44.31886},
            2/3: {'metal': 60.67351, 'P': 71.22673},
            1: {'metal': 83.07559, 'P': 82.30645},
            4/3: {'metal': 85, 'P': 83}
        }
        
        metal_P_to_residue = 0
        
        for metal in ['Fe2','Fe3','Ca2','Mg2']:
            metal_P_to_residue += fermentate.imass[metal]*(1 - metal_P_release[self.food_sludge_ratio]['metal']/100)
            fermentate.imass[metal] *= metal_P_release[self.food_sludge_ratio]['metal']/100
        
        metal_P_to_residue += fermentate.imass['PO4']*(1 - metal_P_release[self.food_sludge_ratio]['P']/100)
        fermentate.imass['PO4'] *= metal_P_release[self.food_sludge_ratio]['P']/100
        
        fermentate.imass['Residue'] += metal_P_to_residue
        
    def _design(self):
        fe_sludge, food_waste = self.ins
        
        Q = fe_sludge.F_vol + food_waste.F_vol
        
        D = self.design_results
        
        D['Flow rate'] = Q
        D['HRT'] = self.HRT
        # TODO: any headspace? (e.g., multiply it by 1.2)
        D['Reactor volume'] = Q * self.HRT
        
        self.add_equipment_design()
        
        D.update(self.verticle_mixer.design_results)
        
        # kW
        self.add_power_utility(D['Required mixing power'])
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from(self.ins)
        hx_outs0.copy_flow(hx_ins0)
        hx_ins0.T = self.ins[0].T
        hx_outs0.T = self.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)

# TODO: check    
#     def _cost(self):
#         V = self.design_result['Reactor volume']
#         self.baseline_purchase_costs['Acidogenic fermenter'] = 3000 * V ** 0.6

# =============================================================================
# SelectivePrecipitation
# =============================================================================

class SelectivePrecipitation(SanUnit):
    '''
    Selective precipitation of FePO4 by adding acid and oxidant.
    
    Parameters
    ----------
    ins : iterable
        supernatant, acid, oxidant.
    outs : iterable
        slurry.
    # TODO: this is not used
    target_pH : float
        Target pH for Fe3+ and PO43- precipitation, [-].
    acid_dose : float
        Mass ratio of acid (measured as pure H2SO4) to the supernatant after acidogenic fermentation, [-].
    P_recovery : float
        Phosphorus recovery ratio via of Fe3+ and PO43- precipitation, [-].
    # TODO: what is the temperature for slurry; is cooling necessary?
    T : float
        Required temperature for FePO4 precipitation, [K].
    '''
    _N_ins = 3
    _N_outs = 1
    
    # TODO: may need to update F_BM for heat_exchanger
    auxiliary_unit_names=('heat_exchanger',)
    
    _units = {
        'Required mixing power': 'kW'
    }
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                target_pH=2.0,
                acid_dose=0.003,
                P_recovery=0.80,
                T=40+_C_to_K):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target_pH = target_pH
        self.acid_dose = acid_dose
        self.P_recovery = P_recovery
        self.T = T
        hx_in = Stream(ID=f'{ID}_hx_in')
        hx_out = Stream(ID=f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out)
        self.verticle_mixer = VerticalMixer(ID='verticle_mixer', linked_unit=self)
        self.equipments = (self.verticle_mixer,)
        
    def _run(self):
        supernatant, acid, oxidant = self.ins
        slurry = self.outs[0]
        
        supernatant.phase = 'l'
        acid.phase = 'l'
        oxidant.phase = 'l'
        slurry.phase = 'l'
        
        # acid dosing while adjusting the pH
        acid.imass['H2SO4'] = supernatant.F_mass * self.acid_dose
        # assume the volumetic ratio between H2SO4 and H2O is 1:1
        acid.ivol['H2O'] = acid.ivol['H2SO4']
        
        # H2O2 comsumption (H2O2 + 2Fe2+ -> 2Fe3+ + 2OH-)
        # double H2O2 to ensure complete oxidation
        oxidant.imass['H2O2'] = supernatant.imass['Fe2'] / 56 / 2 * 34 * 2
        # assume the volumetic ratio between H2O2 and H2O is 3:7
        oxidant.ivol['H2O'] = oxidant.ivol['H2O2'] / 3 * 7
        
        slurry.mix_from(self.ins)
        
        # Fe2+ oxidation to Fe3+
        slurry.imass['Fe3'] += slurry.imass['Fe2']
        slurry.imass['Fe2'] = 0
        
        # FePO4 precipitation
        P_mol = slurry.imass['PO4'] / 95
        P_precip_mol = P_mol * self.P_recovery
        
        Fe_mol = slurry.imass['Fe3'] / 56
        Fe_reacted_mol = min(Fe_mol, P_precip_mol)
        
        slurry.imass['PO4'] -= Fe_reacted_mol * 95
        slurry.imass['Fe3'] -= Fe_reacted_mol * 56
        slurry.imass['FePO4_2H2O'] = Fe_reacted_mol * 187
        
        slurry.imass['H2O2'] = 0
        slurry.imass['H2O'] = 0
        slurry.imass['H2O'] = supernatant.F_mass + acid.F_mass + oxidant.F_mass - slurry.F_mass
        
        slurry.T = self.T

    def _design(self):
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from(self.ins)
        hx_outs0.copy_flow(hx_ins0)
        hx_ins0.T = self.ins[0].T
        hx_outs0.T = self.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        
        # D = self.design_results
        
        # TODO: add reactor volume
        # D['Reactor volume'] = 
        
        # self.add_equipment_design()
        
        # D.update(self.verticle_mixer.design_results)
        
        # kW
        # self.add_power_utility(D['Required mixing power'])
    
    def _cost(self):
        pass

# =============================================================================
# HeatDrying
# =============================================================================

# TODO: confirm: unit_heat is higher than needed heat of water evaporation, so this is kind of like an efficiency
# TODO: in Sintering, also use the unit_heat value here and also assume a natural gas combustion efficiency of 80%?
# n: 0.7, assumed
# BM: 3.17, from biosteam/units/heat_exchange.py
@cost(ID='Dryer 1', basis='Half dry solids flow', units='tonne/day',
      cost=220000*80/2/3.17, S=80, CE=qs.CEPCI_by_year[2018], n=0.7, BM=3.17)
@cost(ID='Dryer 2', basis='Half dry solids flow', units='tonne/day',
      cost=220000*80/2/3.17, S=80, CE=qs.CEPCI_by_year[2018], n=0.7, BM=3.17)
class HeatDrying(SanUnit):
    '''
    Heat drying of sludge or biosolids.
    
    from [1]: 220000 2018$ CAPEX/dry tonne solids/day for a 80 dry tonne/day.
    # assume CAPEX / installed cost = 2 and installed cost / purchase cost (BM) = 3.17
    
    Parameters
    ----------
    ins : iterable
        input_sludge.
    outs : iterable
        dried_solids, vapor.
    target_moisture : float
        Targeted moisture content, [-].
    T : float
        Heat drying temperature, [K].
    unit_heat : float
        Energy for removing unit water from solids, [GJ/tonne water].
    unit_electricity : float
        Electricity for heat drying, [kWh/dry tonne solids].
    natural_gas_HHV : float
        Higher heating value of natural gas, [MJ/m3].
    
    See Also
    --------
    :class:`exposan.htl.landscape_sanunits.HeatDrying` under `pfas` branch
    
    References
    ----------
    [1] Hao, X.; Chen, Q.; van Loosdrecht, M. C. M.; Li, J.; Jiang, H.
        Sustainable Disposal of Excess Sludge: Incineration without
        Anaerobic Digestion. Water Research 2020, 170, 115298.
        https://doi.org/10.1016/j.watres.2019.115298.
    '''
    _N_ins = 2
    _N_outs = 2
    
    _units = {'Half dry solids flow': 'tonne/day'}
    
    # TODO: confirm parameter values
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', lifetime=25, target_moisture=0,
                 T=90 + _C_to_K, unit_heat=4.5, natural_gas_HHV=39,
                 unit_electricity=214):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, lifetime=lifetime)
        self.target_moisture = target_moisture
        self.T = T
        self.unit_heat = unit_heat
        self.unit_electricity = unit_electricity
        self.natural_gas_HHV = natural_gas_HHV
    
    def _run(self):
        input_sludge, natural_gas = self.ins
        dried_solids, vapor = self.outs
        
        input_sludge.phase = 'l'
        natural_gas.phase = 'g'
        dried_solids.phase = 's'
        vapor.phase = 'g'
        
        # dry tonne/h
        self.dry_solids = (input_sludge.F_mass - input_sludge.imass['H2O'])/1000
        
        # tonne/h
        dried_solids_mass_flow = self.dry_solids/(1 - self.target_moisture)
        
        dried_solids.copy_like(input_sludge)
        
        dried_solids.imass['H2O'] = dried_solids_mass_flow*1000*self.target_moisture
        
        # TODO: heat drying may remove VFAs and ethanol
        vapor.imass['H2O'] = input_sludge.F_mass - dried_solids.F_mass
        
        # use natural gas for heat drying base on the BEAM*2024 model
        natural_gas.ivol['CH4'] = vapor.imass['H2O']/1000*self.unit_heat*1000/self.natural_gas_HHV
        
        dried_solids.T = vapor.T = self.T
    
    def _design(self):
        self.design_results['Half dry solids flow'] = self.dry_solids*24/2
        
        # kW
        self.add_power_utility(self.dry_solids*self.unit_electricity)

# =============================================================================
# Sintering
# =============================================================================

class Sintering(SanUnit):
    '''
    Fermentation of sludge and food waste.
    
    Parameters
    ----------
    ins : iterable
        feed, natural_gas, air.
    outs : iterable
        gas, product.
    T : float
        Sintering precipitation, [K].
    combustion_eff : float
        Combustion efficiency to convert the energy embedded in feed to heat.
    natural_gas_HHV : float
        Higher heating value of natural gas, [MJ/m3].
    '''
    _N_ins = 2
    _N_outs = 2
    
    # TODO: add later
    _units = {}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 T=700 + _C_to_K, combustion_eff=0.8, natural_gas_HHV=39): 
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.T = T
        self.combustion_eff = combustion_eff
        self.natural_gas_HHV = natural_gas_HHV
        
    def _run(self):
        # TODO: air was ignored for now
        feed, natural_gas = self.ins
        gas, product = self.outs
        
        feed.phase = 's'
        natural_gas.phase = 'g'
        product.phase = 's'
        gas.phase = 'g'
        
        product.T = gas.T = self.T
        
        cmps = self.components
        rxns = []
        for cmp in [cmp for cmp in cmps if feed.imass[cmp.ID] > 0]:
            if cmp.locked_state in ('l', 's') and (not cmp.organic or cmp.degradability=='Undegradable'):
                continue
            
            rxn = cmp.get_combustion_reaction()
            
            if rxn:
                rxns.append(rxn)
        
        combustion_rxns = self.combustion_reactions = ParallelReaction(rxns)
        
        gas.copy_flow(feed)
        combustion_rxns.force_reaction(gas.mol)
        
        product.copy_flow(gas, IDs=('Fe3', 'PO4', 'Ca2', 'Mg2', 'FePO4_2H2O'), remove=True)
        
        # TODO: 4.5 GJ/tonne is from HeatDrying
        # kJ/hr
        NG_energy_flow = ((gas.H + gas.HHV + product.H + product.HHV + 4.5*1000*product.imass['FePO4_2H2O']/187*36) - (feed.H + feed.HHV))/self.combustion_eff
        
        natural_gas.ivol['CH4'] = NG_energy_flow/1000/self.natural_gas_HHV
        
        product.imass['FePO4'] = product.imass['FePO4_2H2O']/187*151
        gas.imass['H2O'] += product.imass['FePO4_2H2O']/187*36
        
        product.imass['FePO4_2H2O'] = 0
        gas.imass['O2'] = 0
    
    def _design(self):
        pass
    
    def _cost(self):
        pass