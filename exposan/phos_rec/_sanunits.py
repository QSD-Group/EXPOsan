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
    
    '''
    
    food_sludge_ratio : the mass ratio of organics bewteen food wastes and sludge, [-]
    food_waste_moisture ： the mositure content of the input food wastes, [-]
    org_to_gas ： the mass ratio of total organics converted to gas, [-]
    org_to_vfa ：the mass ratio of total organics converted to VFAs, [-]
    org_to_ethanol ：the mass ratio of total organics converted to ethanol, [-]
    org_to_residue ：the mass ratio of total organics in the residue after acidogenic fermentation, [-]
    VFA_ratio ：the mass fractions of individual VFAs (Acetic_acid,Propionic_acid, Butyric_acid, Valeric_acid, and Lactic_acid), [-]
    metal_release ：the release ratio pf metal ions during acidogenic fermentation (mainly the Fe2, Fe3, Ca2, Mg2), [-]
    P_release ：the phosphorus release ratio during acidogenic fermentation, [-]
    fe_reduction ：the Fe reduction and release ratio during acidogenic fermentation, [-]
    fermentate ：the slurry inclding both supernatant and solid residue after acidogenic fermentation, [-]
    metal_P_release ：the release ratio of metals and phosphorus after acidogenic fermentation, [-]
    metal_P_to_residue ：the ratio of metal elements retained in the sludge (residue), [-]
    
    '''
   
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



        
class SelectivePrecipitation(SanUnit):
    
    
    
    '''
    
    target_pH : the target pH (pH=2) for Fe3+ and PO43- precipitation, [-]
    acid_dose : the mass ratio of acid (H2SO4) to the supernatant after acidogenic fermentation, [-]
    P_recovery : the phosphorus recovery ratio via of Fe3+ and PO43- precipitation, [-]
    T : the required temperature for the FePO4 precipitation process, [oC]
    supernatant : the liquid phase obtained after acidogenic fermentation and centrifugation, [-]
    acid : H2SO4 mixed with water at a 1:1 volume ratio, [-]
    oxidant : H2O2 mixed with water at a 3:7 volume ratio, [-]
    slurry : the mixture containing Fe-P precipitate, supernatant, acid, and oxidant after selective precipitation, [-]    
    P_mol : the molar amount of PO43- in the supernatant after acidogenic fermentaion and oxidation by H2O2, [-]
    P_precip_mol : the molar amount of PO43- precipitated during selective precipitation, [-]
    Fe_reacted_mol : the molar amount of Fe3+ consumed during selective precipitation, [-]
    FePO4_2H2O : the precipitated product formed from Fe3+ and PO43- at pH 2 and temperature of 40 degree Celsius, [-]
        
    '''  
    
    _N_ins = 3
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                target_pH=2.0,
                acid_dose=0.02, # H2SO4-supernatant mass ratio
                P_recovery=0.80,
                T=40):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.target_pH = target_pH
        self.acid_dose = acid_dose
        self.P_recovery = P_recovery
        self.T = T
        
    def _run(self):
        supernatant, acid, oxidant = self.ins
        slurry = self.outs[0]
        
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
        
        slurry.imass['H2O'] = 0
        slurry.imass['H2O'] = supernatant.F_mass + acid.F_mass + oxidant.F_mass - slurry.F_mass
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
   
    
    
        
    
    

