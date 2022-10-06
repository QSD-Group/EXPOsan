#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Jianan Feng <jiananf2@illinois.edu>
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References:
    
(1) Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J. 
    Quantitative Evaluation of an Integrated System for Valorization of Wastewater 
    Algae as Bio-Oil, Fuel Gas, and Fertilizer Products. 
    Environ. Sci. Technol. 2018, 52 (21), 12717–12727. 
    https://doi.org/10.1021/acs.est.8b04035.
    
(2) Leow, S.; Witter, J. R.; Vardon, D. R.; Sharma, B. K.;
    Guest, J. S.; Strathmann, T. J. Prediction of Microalgae Hydrothermal
    Liquefaction Products from Feedstock Biochemical Composition.
    Green Chem. 2015, 17 (6), 3584–3599. https://doi.org/10.1039/C5GC00574D.

    
(3) Snowden-Swan, L. J.; Zhu, Y.; Jones, S. B.; Elliott, D. C.; Schmidt, A. J.; 
    Hallen, R. T.; Billing, J. M.; Hart, T. R.; Fox, S. P.; Maupin, G. D. 
    Hydrothermal Liquefaction and Upgrading of Municipal Wastewater Treatment 
    Plant Sludge: A Preliminary Techno-Economic Analysis; 
    PNNL--25464, 1258731; 2016; https://doi.org/10.2172/1258731.


(4) Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
    Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
    Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
    Process Design and Economics for the Conversion of Algal Biomass to Hydrocarbons: 
    Whole Algae Hydrothermal Liquefaction and Upgrading; PNNL--23227, 1126336; 2014; 
    https://doi.org/10.2172/1126336.
'''

from qsdsan import SanUnit

__all__ = (
    'HTL',
    'HT',
    'AcidExtraction',
    'StruvitePrecipitation',
    'CHG',
    'MembraneDistillation')

class HTL(SanUnit):
    
    '''
    HTL converts dewatered sludge to biocrude, aqueous, off-gas, and biochar under
    elevated temperature (350°C) and pressure. The products percentage (wt%) can be evaluatad
    using revised MCA model (Li et al., 2017, Leow et al., 2018) with known sludge
    composition (protein%, lipid%, and carbohydrate%)
    
    Parameters
    ----------
    ins: Iterable (stream)
        Dewatered biosolid or sludge
    outs: Iterable (stream)
        biochar, others
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',
                 biochar_C_N_ratio=15.5,biochar_C_P_ratio=2.163,COD_ratio=0.602,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.biochar_C_N_ratio=biochar_C_N_ratio
        self.biochar_C_P_ratio=biochar_C_P_ratio
        self.COD_ratio=COD_ratio

    _N_ins = 1
    _N_outs = 2
        
    def _run(self):
        
        cmps=self.components
        
        dewatered_sludge=self.ins[0]
        biochar,others=self.outs
        
        dewatered_sludge_afdw = dewatered_sludge.imass['Sludge_lipid']\
            + dewatered_sludge.imass['Sludge_protein']\
            + dewatered_sludge.imass['Sludge_carbo']
                               
        lipid_ratio = dewatered_sludge.imass['Sludge_lipid']/dewatered_sludge_afdw
        protein_ratio = dewatered_sludge.imass['Sludge_protein']/dewatered_sludge_afdw
        carbo_ratio = dewatered_sludge.imass['Sludge_carbo']/dewatered_sludge_afdw
        
        sludge_C_mass = dewatered_sludge_afdw*cmps.Sludge_carbo.i_C
        sludge_N_mass = dewatered_sludge_afdw*cmps.Sludge_carbo.i_N
        sludge_P_mass = dewatered_sludge_afdw*cmps.Sludge_carbo.i_P
        
        # if abs(lipid_ratio+protein_ratio+carbo_ratio-1)>0.05:
        #     warn('The sum of the sludge composition is not 1')
        
        others.imass['H2O'] = dewatered_sludge.imass['H2O']
        
        Biocrude_mass = (0.846*lipid_ratio + 0.445*protein_ratio\
            + 0.205*carbo_ratio) * dewatered_sludge_afdw
        others.imass['C_c']=Biocrude_mass*0.721
        others.imass['N_c']=Biocrude_mass*protein_ratio*0.133
        others.imass['Others_c']=Biocrude_mass-others.imass['C_c']-others.imass['N_c']
        
        biochar_mass = 0.377*carbo_ratio * dewatered_sludge_afdw
        biochar.imass['C_s']=min(1.75*carbo_ratio, 0.65)*biochar_mass
        biochar.imass['N_s']=biochar.imass['C_s']/self.biochar_C_N_ratio
        biochar.imass['P_s']=min(sludge_P_mass,biochar.imass['C_s']/self.biochar_C_P_ratio)
        biochar.imass['Others_s']=biochar_mass-biochar.imass['C_s']-biochar.imass['N_s']-\
            biochar.imass['P_s']
            
        others.imass['CO2'] = (0.074*protein_ratio + 0.418*carbo_ratio)\
            * dewatered_sludge_afdw
            
        HTLaqueous_mass = (0.154*lipid_ratio + 0.481*protein_ratio)\
            * dewatered_sludge_afdw
        cmps = self.components
        cmps.C_l.i_COD=self.COD_ratio
        others.imass['C_l']=sludge_C_mass-biochar.imass['C_s']-others.imass['C_c']-others.imass['CO2']*12/44
        others.imass['N_l']=sludge_N_mass-biochar.imass['N_s']-others.imass['N_c']
        others.imass['P_l']=sludge_P_mass-biochar.imass['P_s']
        others.imass['Others_l']=HTLaqueous_mass-others.imass['C_l']-others.imass['N_l']-others.imass['P_l']

    def _design(self):
        pass
    
    def _cost(self):
        pass

class HT(SanUnit):
    
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature (405°C) and pressure
    to produce upgraded biooil. Co-products include fuel gas and char. The amount of
    biooil and fuel gas can be estimated using values from Li et al., 2018.
    The amount of char can be calculated based on mass closure.
    
    Parameters
    ----------
    ins: Iterable (stream)
        biocurde
    outs: Iterable (stream)
        char, others
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',\
                 biooil_ratio=0.581,gas_ratio=0.047,c5h12_ratio=0.09,co_ratio=0.128,\
                     co2_ratio=0.007,c2h6_ratio=0.188,c3h8_ratio=0.107,biooil_C_ratio=0.855,
                     biooil_N_ratio=0.01,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.biooil_ratio=biooil_ratio
        self.gas_ratio=gas_ratio
        self.c5h12_ratio=c5h12_ratio
        self.co_ratio=co_ratio
        self.co2_ratio=co2_ratio
        self.c2h6_ratio=c2h6_ratio
        self.c3h8_ratio=c3h8_ratio
        self.biooil_C_ratio=biooil_C_ratio
        self.biooil_N_ratio=biooil_N_ratio

    _N_ins = 1
    _N_outs = 2
        
    def _run(self):
        
        biocrude = self.ins[0]
        char,others = self.outs
        
        massin = biocrude.imass['C_c']+biocrude.imass['N_c']+biocrude.imass['Others_c']
        biooil_mass = massin * self.biooil_ratio
        gas_mass = massin * self.gas_ratio
        char_mass = massin * (1-self.biooil_ratio-self.gas_ratio)
        
        others.imass['C_o']=biooil_mass*self.biooil_C_ratio
        others.imass['N_o']=biooil_mass*self.biooil_N_ratio
        others.imass['Others_o']=biooil_mass*(1-self.biooil_C_ratio-self.biooil_N_ratio)
        
        for i in (self.c5h12_ratio,self.co_ratio,self.co2_ratio,self.c2h6_ratio,self.c3h8_ratio,
                    1-self.c5h12_ratio-self.co_ratio-self.co2_ratio-self.c2h6_ratio-\
                    self.c3h8_ratio):
            if i<0: i=0
    
        fuelgas_HT_composition = {
            'C5H12':self.c5h12_ratio,
            'CO':self.co_ratio,
            'CO2':self.co2_ratio,
            'C2H6':self.c2h6_ratio,
            'C3H8':self.c3h8_ratio,
            'CH4':1-self.c5h12_ratio-self.co_ratio-self.co2_ratio-self.c2h6_ratio-\
                self.c3h8_ratio
            }
        for name,ratio in fuelgas_HT_composition.items():
            others.imass[name] = gas_mass * ratio
        
        char.imass['C_s'] = biocrude.imass['C_c']-others.imass['C_o']-others.imass['CH4']*12/16-\
            others.imass['CO']*12/28-others.imass['CO2']*12/44-others.imass['C2H6']*24/30-\
            others.imass['C3H8']*36/44-others.imass['C5H12']*60/72
        char.imass['N_s'] = biocrude.imass['N_c']-others.imass['N_o']
        char.imass['Others_s'] = char_mass-char.imass['C_s']-char.imass['N_s']
    
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
   
class AcidExtraction(SanUnit):
    
    '''
    H2SO4 is added to biochar from HTL to extract P. 
    
    Parameters
    ----------
    ins: Iterable (stream)
        biochar, acid_P
    outs: Iterable (stream)
        residual, extracted_P
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',P_recovery_rate=0.95,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.P_recovery_rate=P_recovery_rate

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        biochar, acid = self.ins
        residual, extracted = self.outs
        
        biochar_mass = biochar.imass['C_s']+biochar.imass['N_s']+biochar.imass['P_s']+biochar.imass['Others_s']
        
        acid.imass['H2SO4']=biochar_mass*10*0.5*98/1000 #0.5 M H2SO4 10 mL/1 g Biochar
        acid.imass['H2O']=biochar_mass*10*1.05-acid.imass['H2SO4'] #0.5 M H2SO4 density: 1.05 kg/L 
        #https://www.fishersci.com/shop/products/sulfuric-acid-1n-0-5m-standard-solution-thermo-
        #scientific/AC124240010 (accessed 10-6-2022)
        
        residual.copy_like(biochar)
        residual.imass['P_s']=biochar.imass['P_s']*(1-self.P_recovery_rate)
        
        extracted.imass['H2O']=acid.imass['H2O']
        extracted.imass['Others_l']=acid.imass['H2SO4']
        extracted.imass['P_l']=biochar.imass['P_s']*self.P_recovery_rate
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
    
class HTLmixer(SanUnit):
    '''
    
    Parameters
    ----------
    ins: Iterable (stream)
        aqueous, extracted_P
    outs: Iterable (stream)
        mixture
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',
                  **kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)

    _N_ins = 2
    _N_outs = 1
        
    def _run(self):
        
        aqueous, extracted = self.ins
        mixture = self.outs[0]
        
        mixture.copy_like(extracted)
        mixture.imass['C_l']+=aqueous.imass['C_l']
        mixture.imass['N_l']+=aqueous.imass['N_l']
        mixture.imass['P_l']+=aqueous.imass['P_l']
        mixture.imass['Others_l']+=aqueous.imass['Others_l']
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
    
class StruvitePrecipitation(SanUnit):
    '''
    extracted_P and HTL aqueous are mixed together (Mixer) before adding MgCl2 and struvite precipitation.
    
    Parameters
    ----------
    ins: Iterable (stream)
        mixture, supply_MgCl2
    outs: Iterable (stream)
        struvite, effluent
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',Mg_P_ratio=1,
                 P_recovery_rate=0.95,P_in_struvite=0.127,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.Mg_P_ratio=Mg_P_ratio
        self.P_recovery_rate=P_recovery_rate
        self.P_in_struvite=P_in_struvite

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        mixture, supply_MgCl2 = self.ins
        struvite, effluent = self.outs
        
        supply_MgCl2.imass['MgCl2']=mixture.imass['P_l']/31*95.211*self.Mg_P_ratio
        struvite.imass['Struvite']=mixture.imass['P_l']*self.P_recovery_rate/self.P_in_struvite
        
        effluent.copy_like(mixture)
        effluent.imass['P_l'] -= struvite.imass['Struvite']*self.P_in_struvite
        effluent.imass['N_l'] -= struvite.imass['Struvite']*self.P_in_struvite/31*14
        effluent.imass['Others_l'] += (supply_MgCl2.imass['MgCl2']-\
                                       struvite.imass['Struvite']*self.P_in_struvite*(1+14/31))
          
    def _design(self):
        pass
    
    def _cost(self):
        pass
    
class CHG(SanUnit):
   
    '''
    CHG serves to reduce the COD content in the aqueous phase and produce fuel gas 
    under elevated temperature (350°C) and pressure.
    
    Parameters
    ----------
    ins: Iterable (stream)
        feed
    outs: Iterable (stream)
        fuelgas, effluent
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',ch4_ratio=0.244,
                 co_ratio=0.029,co2_ratio=0.150,c2h6_ratio=0.043,toc_tc_ratio=0.764,
                 toc_to_gas_c_ratio=0.262,COD_removal_rate=0.975,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.ch4_ratio=ch4_ratio
        self.co_ratio=co_ratio
        self.co2_ratio=co2_ratio
        self.c2h6_ratio=c2h6_ratio
        self.toc_tc_ratio=toc_tc_ratio
        self.toc_to_gas_c_ratio=toc_to_gas_c_ratio
        self.COD_removal_rate=COD_removal_rate
        
    _N_ins = 1
    _N_outs = 2
        
    def _run(self):
        
        feed = self.ins[0]
        fuelgas, effluent = self.outs
        
        for i in (self.ch4_ratio,self.co_ratio,self.co2_ratio,self.c2h6_ratio,
                  1-self.ch4_ratio-self.co_ratio-self.co2_ratio-self.c2h6_ratio):
            if i < 0: i = 0
        
        fuel_gas_mass = feed.imass['C_l']*self.toc_tc_ratio*self.toc_to_gas_c_ratio/(self.\
            ch4_ratio*12/16+self.co_ratio*12/28+self.co2_ratio*12/44+self.c2h6_ratio*24/30)

        fuelgas_CHG_composition = {
            'CH4':self.ch4_ratio,
            'CO':self.co_ratio,
            'CO2':self.co2_ratio,
            'C2H6':self.c2h6_ratio,
            'H2':1-self.ch4_ratio-self.co_ratio-self.co2_ratio-self.c2h6_ratio
            }
        
        for name,ratio in fuelgas_CHG_composition.items():
            fuelgas.imass[name] = fuel_gas_mass * ratio
        
        cmps=self.components
        cmps.C_l_after_COD_removal.i_COD=cmps.C_l.i_COD
        cmps.C_l_after_COD_removal.i_COD*=(1-self.COD_removal_rate)
        
        effluent.copy_like(feed)
        effluent.imass['C_l']=0
        effluent.imass['C_l_after_COD_removal']=feed.imass['C_l']*(1-self.toc_tc_ratio*self.toc_to_gas_c_ratio)
        effluent
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
    
class MembraneDistillation(SanUnit):
    
    '''
    Membrane distillation will be modeled using vapor-liquid-liquid equilibrium later.
    Here is simplified.
    
    Parameters
    ----------
    ins: Iterable (stream)
        influent, acid
    outs: Iterable (stream)
        ammoniasulfate, ww
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',N_S_ratio=2,
                 N_recovery_rate=0.9,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.N_S_ratio=N_S_ratio
        self.N_recovery_rate=N_recovery_rate

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        influent, acid = self.ins
        ammoniasulfate, ww = self.outs
        
        acid.imass['H2SO4'] = influent.imass['N_l']/14/self.N_S_ratio*98
        acid.imass['H2O'] = acid.imass['H2SO4']*1000/98/0.5*1.05-acid.imass['H2SO4']
        ammoniasulfate.imass['NH42SO4']=influent.imass['N_l']*self.N_recovery_rate/14*132
        ammoniasulfate.imass['H2O']=acid.imass['H2O']
        ww.copy_like(influent)
        ww.imass['Others_l']+=ww.imass['N_l']*self.N_recovery_rate+acid.imass['H2SO4']-\
            ammoniasulfate.imass['NH42SO4']
        ww.imass['N_l']*=(1-self.N_recovery_rate)
        
    def _design(self):
        pass
    
    def _cost(self):
        pass