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

class AcidogenicFermenter(SanUnit):
    _N_ins = 2
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 food_sludge_ratio=1,
                 food_waste_moisture=0.2,
                 org_to_gas=0.05,
                 org_to_vfa=0.65, 
                 org_to_ethanol=0.02,
                 org_to_residue=0.25,
                 VFA_ratio={'Ac': 0.5, 'Pr': 0.24, 'Bu': 0.23, 'Va': 0.02, 'Lac': 0.01},
                 metal_release=0.8,
                 P_release=0.82,
                 HRT=3,
                 fe_reduction = 0.98): 
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.food_sludge_ratio = food_sludge_ratio
        self.food_waste_moisture = food_waste_moisture
        self.org_to_gas = org_to_gas
        self.org_to_vfa = org_to_vfa
        self.org_to_ethanol = org_to_ethanol
        self.org_to_residue = org_to_residue
        self.VFA_ratio = VFA_ratio
        self.metal_release = metal_release
        self.P_release = P_release
        self.HRT = HRT
        self.fe_reduction = fe_reduction
        
    def _run(self):
        
        fe_sludge, food_waste = self.ins
        gas, fermentate = self.outs
        
        if self.food_sludge_ratio not in [0, 1/3, 2/3, 1, 4/3]:
            raise ValueError('food_sludge_ratio must be one of the follow: 0, 1/3, 2/3, 1, 4/3.')
        
        food_waste.imass['Org'] = fe_sludge.imass['Org']*self.food_sludge_ratio
        food_waste.imass['H2O'] = food_waste.imass['Org']/(1-self.food_waste_moisture) * self.food_waste_moisture
        
        fermentate.mix_from((fe_sludge, food_waste))
        
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
            raise Warning('org cannot be balanced')
        
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
        
        metal_P_to_residual = 0
        
        for metal in ['Fe2','Fe3','Ca2','Mg2']:
            metal_P_to_residual += fermentate.imass[metal]*(1 - metal_P_release[self.food_sludge_ratio]['metal']/100)
            fermentate.imass[metal] *= metal_P_release[self.food_sludge_ratio]['metal']/100
        
        metal_P_to_residual += fermentate.imass['PO4']*(1 - metal_P_release[self.food_sludge_ratio]['P']/100)
        fermentate.imass['PO4'] *= metal_P_release[self.food_sludge_ratio]['P']/100
        
        fermentate.imass['Residue'] += metal_P_to_residual
        
#     def _design(self):
#         fe_sludge, food_waste = self.ins
#         Q = fe_sludge.F_vol + food_waste.F_vol
        
#         self.design_resluts['Flow rate (m3/d)'] = Q
#         self.design_results['HRT(d)']=self.HRT
#         self.design_resluts['Reactor volume(m3)'] = Q * self.HRT
        
#         pass
    
#     def _cost(self):
#         V = self.design_result['Reactor volume (m3)']
#         self.baseline_purchase_costs['Acidogenic fermenter'] = 3000 * V **0.6
        
#         pass



        
# class SelectivePrecipitation(SanUnit):
#     _N_ins = 3
#     _N_outs = 1
    
    
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
        
        
#     def _run(self):
#         supernatant, H2SO4, H2O2 = self.ins
#         Slurry = self.outs
        
        
    
    
    
   
    
    
        
    
    

