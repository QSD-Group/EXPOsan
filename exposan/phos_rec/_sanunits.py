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
from qsdsan import SanUnit
import qsdsan as qs

class AcidogenicFermenter(SanUnit):
    _N_ins = 2
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 sludge_food_ratio=1,
                 food_waste_moisture=0.2,
                 org_to_gas=0.1,
                 org_to_vfa=0.6,               
                 gas_split = {'CO2': 0.85,
                              'H2': 0.1,
                              'CH4': 0.05},
                 HRT=3,
                 fe_reduction = 0.98): 
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.sludge_food_ratio = sludge_food_ratio
        self.food_waste_moisture = food_waste_moisture
        self.org_to_gas = org_to_gas
        self.org_to_vfa = org_to_vfa
        self.gas_split = gas_split
        self.HRT = HRT
        self.fe_reduction = fe_reduction
        
        
        
    def _run(self):
        
        fe_sludge, food_waste = self.ins
        gas, fermentate = self.outs
        cmps = qs.get_components()
        
        # TODO: sludge:food waste ratio, plus cmps = qs.get_componments()?
        
        food_waste.imass['Org']=fe_sludge.imass['Org']/self.sludge_food_ratio
        food_waste.imass['H2O']=food_waste.imass['Org']/(1-self.food_waste_moisture) * self.food_waste_moisture
        
        org_total = fe_sludge.imass['Org']+food_waste.imass['Org']
        
# =============================================================================
#         Gas production
# =============================================================================

        C_mol = org_total * cmps['Org'].i_C / 12
        gas_total = C_mol * self.org_to_gas
        
        gas.imass['CO2']=gas_total*self.gas_split['CO2'] * 44 #the same to the upon code
        gas.imass['H2']=gas_total*self.gas_split['H2']*2
        gas.imass['CH4']=gas_total*self.gas_split['CH4']*16
        
# =============================================================================
#         VFAs and inorganics in fermentate
# =============================================================================
        # TODO: use split to define the ratio of VFAS

        vfa_mass = org_total * self.org_to_vfa
        fermentate.imass['Acetic_acid'] = vfa_mass * 0.6
        fermentate.imass['Propionic_acid'] = vfa_mass * 0.15
        fermentate.imass['Butyric_acid'] = vfa_mass * 0.15
        fermentate.imass['Valeric_acid'] = vfa_mass * 0.05
        fermentate.imass['Lactic_acid'] = vfa_mass * 0.1
        fermentate.imass['Ethanol'] = vfa_mass * 0.1
        fermentate.imass['Org'] = org_total * (1 - self.org_to_vfa - self.org_to_gas)*0.65             # *0.65
        
        for cmp in ('H2O', 'PO4', 'Fe3', 'Ca2', 'Mg2'):
            fermentate.imass[cmp] = (fe_sludge.imass[cmp] + food_waste.imass[cmp])

        fermentate.imass['Fe2'] = fermentate.imass['Fe3'] * self.fe_reduction
        fermentate.imass['Fe3'] -= fermentate.imass['Fe2']
        
        
        
        fermentate.imass['PO4'] *= 0.85
        fermentate.imass['Fe3'] *= 0.8
        fermentate.imass['Fe2'] *= 0.95
        fermentate.imass['Ca2'] *= 0.75
        fermentate.imass['Mg2'] *= 0.75

                
        fermentate.imass['Residue'] = fe_sludge.F_mass + food_waste.F_mass - gas.F_mass - fermentate.F_mass
        
# class SelectivePrecipitation(SanUnit):
#     _N_ins = 3
#     _N_outs = 2
    
    
#     def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
#                 target_pH=2.0,
#                 acid_dose=0.02,
#                 fe2_oxidation=1,
#                 P_recovery=0.85,
#                 T=40):
#         SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
#         self.target_pH = target_pH
#         self.acid_dose = acid_dose
#         self.fe2_oxidation = fe2_oxidation
#         self.P_recovery = P_recovery
#         self.T = T
        
    
    
    
    # def _design(self):
    #     fe_sludge, food_waste = self.ins
    #     Q = fe_sludge.F_vol + food_waste.F_vol
        
    #     self.design_resluts['Flow rate (m3/d)'] = Q
    #     self.design_results['HRT(d)']=self.HRT
    #     self.design_resluts['Reactor volume(m3)'] = Q * self.HRT
        
    #     pass
    
    # def _cost(self):
    #     V = self.design_result['Reactor volume (m3)']
    #     self.baseline_purchase_costs['Acidogenic fermenter'] = 3000 * V **0.6
        
    #     pass
    
    
        
    
    

