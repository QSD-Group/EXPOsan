#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Jianan Feng <jiananf2@illinois.edu>
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

import biosteam as bst
from qsdsan import SanUnit
from qsdsan.sanunits import HXutility
from biosteam.units import Flash as FL
from biosteam import PowerUtility
from thermosteam import MultiStream, separations
from math import pi
import numpy as np
from biosteam.units import design_tools as design

__all__ = (
    'SludgeLab',
    'HTL',
    'HT',
    'HC',
    'AcidExtraction',
    'HTLmixer',
    'StruvitePrecipitation',
    'CHG',
    'MembraneDistillation'
    )

# =============================================================================
# Sludge Lab
# =============================================================================

class SludgeLab(SanUnit):
    '''
    SludgeLab is a fake unit that can set up sludge biochemical compositions and 
    calculate sludge elemental compositions.
    
    Model method: just _run, no _design or _cost.
    
    Parameters
    ----------
    ins: Iterable (stream)
        fake_sludge
    outs: Iterable (stream)
        real_sludge
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream', 
                 sludge_moisture=0.99, sludge_dw_protein=0.341,
                 sludge_dw_lipid=0.226, sludge_dw_carbo=0.167, sludge_P_ratio = 0.019,
                 #data are from SS PNNL 2021
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.sludge_moisture = sludge_moisture
        self.sludge_dw_protein = sludge_dw_protein
        self.sludge_dw_carbo = sludge_dw_carbo
        self.sludge_dw_lipid = sludge_dw_lipid
        self.sludge_dw_ash = 1 - sludge_dw_protein - sludge_dw_carbo - sludge_dw_lipid
        self.sludge_P_ratio = sludge_P_ratio
        #set P as an independent variable, assume S is 0
        
    _N_ins = 1
    _N_outs = 1
    
    def _run(self):
        
        fake_sludge = self.ins[0]
        real_sludge = self.outs[0]
        
        real_sludge.imass['H2O'] = fake_sludge.F_mass * self.sludge_moisture
        sludge_dw = fake_sludge.F_mass*(1 - self.sludge_moisture)
        real_sludge.imass['Sludge_protein'] = sludge_dw*self.sludge_dw_protein
        real_sludge.imass['Sludge_carbo'] = sludge_dw*self.sludge_dw_carbo
        real_sludge.imass['Sludge_lipid'] = sludge_dw*self.sludge_dw_lipid
        real_sludge.imass['Sludge_ash'] = sludge_dw*self.sludge_dw_ash
          
        
    #all sludge elemental analysis are based on empirical equation
    @property
    def sludge_C_ratio(self):
       return self.sludge_dw_carbo*0.44 + self.sludge_dw_lipid*0.75 + self.sludge_dw_protein*0.53
    #https://pubmed.ncbi.nlm.nih.gov/2061559/ (accessed 2022-10-27)
    #https://encyclopedi/a2.thefreedictionary.com/Proteins (accessed 2022-10-27)
    
    @property
    def sludge_H_ratio(self):
       return self.sludge_C_ratio/7
    #based on SS PNNL 2021 data, H ~ C/7
   
    @property
    def sludge_N_ratio(self):
       return self.sludge_dw_protein*0.16
    #https://www.fao.org/3/y5022e/y5022e03.htm#:~:text=On%20the%20basis%20of%20
    #early,is%20confounded%20by%20two%20considerations (accessed 2022-10-27)
   
    @property
    def sludge_P_ratio(self):
       return self._sludge_P_ratio
    #set P as an indepedent variable since hard to find any association with sludge biochemical compositions
    
    @sludge_P_ratio.setter
    def sludge_P_ratio(self, i):
        if not 0 <= i <= 1:
            raise AttributeError('`sludge_P` must be within [0, 1], '
                                f'the provided value {i} is outside this range.')
        self._sludge_P_ratio = i
    
    @property
    def sludge_O_ratio(self):
       return 1 - self.sludge_C_ratio - self.sludge_H_ratio - self.sludge_N_ratio -\
           self.sludge_P_ratio - self.sludge_dw_ash*0.75
    #sludge_O is calculated based on mass balance closure
    #asd * 0.75 since double count some elements. 0.75 is based on SS PNNL 2021.
    
    @property
    def AOSc(self):
       return (3*self.sludge_N_ratio/14 + 2*self.sludge_O_ratio/16 -\
               self.sludge_H_ratio/1)/(self.sludge_C_ratio/12)
   
    def _design(self):
        pass
    
    def _cost(self):
        pass

# =============================================================================
# HTL (ignore three phase separator for now, ask Yalin)
# =============================================================================

class HTL(SanUnit):
    
    '''
    HTL converts dewatered sludge to biocrude, aqueous, off-gas, and biochar under
    elevated temperature (350°C) and pressure. The products percentage (wt%) can be evaluated
    using revised MCA model (Li et al., 2017, Leow et al., 2018) with known sludge
    composition (protein%, lipid%, and carbohydrate%, all afdw%).
    
    Model method: empirical model (based on MCA model and experimental data).
    
    Parameters
    ----------
    ins: Iterable (stream)
        dewatered_sludge
    outs: Iterable (stream)
        biochar, HTLaqueous, biocrude, offgas
    '''
    auxiliary_unit_names=('heat_exchanger',)

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream',
                 biocrude_moisture_content=0.056, #Jones PNNL 2014
                 biochar_C_N_ratio=15.5, biochar_C_P_ratio=2.163, #based on Yalin's data
                 ch4_ratio=0.05, c2h6_ratio = 0.03,
                 biochar_pre=3029.7*6894.76, #Jones 2014: 3029.7 psia
                 HTLaqueous_pre=30*6894.76, #Jones 2014: 30 psia
                 biocrude_pre=30*6894.76,
                 offgas_pre=30*6894.76,
                 eff_T=60+273.15, #Jones 2014
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.biocrude_moisture_content = biocrude_moisture_content

        self.biochar_C_N_ratio = biochar_C_N_ratio
        self.biochar_C_P_ratio = biochar_C_P_ratio
        
        self.ch4_ratio = ch4_ratio
        self.c2h6_ratio = c2h6_ratio
        self.co2_ratio = 1 - self.ch4_ratio - self.c2h6_ratio
        
        self.biochar_pre = biochar_pre
        self.HTLaqueous_pre = HTLaqueous_pre
        self.biocrude_pre = biocrude_pre
        self.offgas_pre = offgas_pre
        self.eff_T = eff_T
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'{ID}_hx', ins=hx_in, outs=hx_out)

    _N_ins = 1
    _N_outs = 4
    
    def _run(self):
        
        dewatered_sludge = self.ins[0]
        biochar, HTLaqueous, biocrude, offgas = self.outs
        
        dewatered_sludge_afdw = dewatered_sludge.imass['Sludge_lipid'] +\
                                dewatered_sludge.imass['Sludge_protein'] +\
                                dewatered_sludge.imass['Sludge_carbo']
                               
        lipid_ratio = dewatered_sludge.imass['Sludge_lipid']/dewatered_sludge_afdw
        protein_ratio = dewatered_sludge.imass['Sludge_protein']/dewatered_sludge_afdw
        carbo_ratio = dewatered_sludge.imass['Sludge_carbo']/dewatered_sludge_afdw

        #the following calculations are based on revised MCA model
        biochar.imass['Biochar'] = 0.377*carbo_ratio*dewatered_sludge_afdw     
        HTLaqueous.imass['HTLaqueous'] = (0.154*lipid_ratio + 0.481*protein_ratio)\
            *dewatered_sludge_afdw #HTLaqueous is TDS in aqueous phase
            
        for ratio in (self.ch4_ratio, self.c2h6_ratio, self.co2_ratio):
            if ratio < 0: ratio = 0
         
        offgas_imass = (0.074*protein_ratio + 0.418*carbo_ratio)\
            *dewatered_sludge_afdw
            
        offgas.imass['CH4'] = offgas_imass*self.ch4_ratio
        offgas.imass['C2H6'] = offgas_imass*self.c2h6_ratio
        offgas.imass['CO2'] = offgas_imass*self.co2_ratio
            
        biocrude.imass['Biocrude'] = (0.846*lipid_ratio + 0.445*protein_ratio\
            + 0.205*carbo_ratio)*dewatered_sludge_afdw
        biocrude.imass['H2O'] = biocrude.imass['Biocrude']/(1 - self.biocrude_moisture_content) -\
            biocrude.imass['Biocrude']
        HTLaqueous.imass['H2O'] = dewatered_sludge.imass['H2O'] - biocrude.imass['H2O'] +\
            dewatered_sludge.imass['Sludge_ash'] #assume ash goes to water
        
        biochar.phase = 's'
        offgas.phase = 'g'
        
        biochar.P = self.biochar_pre
        HTLaqueous.P = self.HTLaqueous_pre
        biocrude.P = self.biocrude_pre
        offgas.P = self.offgas_pre
        
        for stream in self.outs: stream.T = self.eff_T
        
        self.sludgelab = self.ins[0]._source.ins[0]._source.ins[0]._source.ins[0]._source.ins[0]._source

    @property
    def biochar_C_ratio(self):
        return min(1.75*self.sludgelab.sludge_dw_carbo, 0.65) #revised MCA model
    
    @property
    def biochar_N_ratio(self):
        return self.biochar_C_ratio/self.biochar_C_N_ratio 
    
    @property
    def biochar_P_ratio(self):
        return self.biochar_C_ratio/self.biochar_C_P_ratio

    @property
    def biochar_C(self):
        return self.outs[0].F_mass*self.biochar_C_ratio

    @property
    def biochar_N(self):
        return self.outs[0].F_mass*self.biochar_N_ratio

    @property
    def biochar_P(self):
        return min((self.ins[0].F_mass - self.ins[0].imass['H2O'])*self.sludgelab.sludge_P_ratio,
                   self.outs[0].F_mass*self.biochar_P_ratio) #make sure biochar P smaller than total P

    @property
    def biocrude_C_ratio(self):
        return (self.sludgelab.AOSc*(-8.37) + 68.55)/100 #revised MCA model

    @property
    def biocrude_N_ratio(self):
        return 0.133*self.sludgelab.sludge_dw_protein #revised MCA model

    @property
    def biocrude_C(self):
        return self.outs[2].F_mass*self.biocrude_C_ratio

    @property
    def biocrude_N(self):
        return self.outs[2].F_mass*self.biocrude_N_ratio

    @property
    def offgas_C(self):
        return self.outs[3].imass['CO2']*12/44 + self.outs[3].imass['CH4']*12/16 +\
            self.outs[3].imass['C2H6']*24/30

    #C, N, and P in aqueous phase are calculated base on mass balance closure
    @property
    def HTLaqueous_C(self):
        return (self.ins[0].F_mass - self.ins[0].imass['H2O'])*self.sludgelab.sludge_C_ratio -\
            self.biochar_C - self.biocrude_C - self.offgas_C

    @property
    def HTLaqueous_N(self):
        return (self.ins[0].F_mass - self.ins[0].imass['H2O'])*self.sludgelab.sludge_N_ratio -\
            self.biochar_N - self.biocrude_N

    @property
    def HTLaqueous_P(self):
        return (self.ins[0].F_mass - self.ins[0].imass['H2O'])*self.sludgelab.sludge_P_ratio -\
            self.biochar_P

    def _design(self):
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from(self.outs)
        hx_outs0.mix_from(self.outs)
        hx_ins0.T = self.ins[0].T
        hx.T = hx_outs0.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
    
    def _cost(self):
        pass

# =============================================================================
# Acid Extraction
# =============================================================================

class AcidExtraction(SanUnit):
    
    '''
    H2SO4 is added to biochar from HTL to extract P. 
    
    Model method: assume P recovery ratio, add filters in _design and _cost.
    
    Parameters
    ----------
    ins: Iterable (stream)
        biochar, acid
    outs: Iterable (stream)
        residual, extracted
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream', acid_vol=10,
                 P_acid_recovery_ratio=0.95,
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.acid_vol = acid_vol
        self.P_acid_recovery_ratio = P_acid_recovery_ratio

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        biochar, acid = self.ins
        residual, extracted = self.outs
        
        acid.imass['H2SO4'] = biochar.F_mass*self.acid_vol*0.5*98/1000 #0.5 M H2SO4 acid_vol (10 mL/1 g) Biochar
        acid.imass['H2O'] = biochar.F_mass*self.acid_vol*1.05 - acid.imass['H2SO4'] #0.5 M H2SO4 density: 1.05 kg/L 
        #https://www.fishersci.com/shop/products/sulfuric-acid-1n-0-5m-standard-solution-thermo-
        #scientific/AC124240010 (accessed 10-6-2022)
        
        residual.imass['Residual'] = biochar.F_mass*\
            (1 - self.ins[0]._source.biochar_P_ratio*self.P_acid_recovery_ratio)
        
        extracted.copy_like(acid)
        extracted.imass['P'] = biochar.F_mass - residual.F_mass #assume just P can be extracted
        
        residual.phase = 's'
        
        residual.T = extracted.T = biochar.T
        #H2SO4 reacts with biochar to release heat and temperature will be increased mixture's temperature
        
    @property
    def residual_C(self):
        return self.ins[0]._source.biochar_C
    
    @property
    def residual_N(self):
        return self.ins[0]._source.biochar_N

    @property
    def residual_P(self):
        return self.ins[0]._source.biochar_P - self.outs[1].imass['P']
        
    def _design(self):
        pass
    
    def _cost(self):
        pass
    
# =============================================================================
# HTL mixer
# =============================================================================

class HTLmixer(SanUnit):
    '''
    
    Model method: elements separation, need _design and _cost.
    
    Parameters
    ----------
    ins: Iterable (stream)
        HTLaqueous, extracted
    outs: Iterable (stream)
        mixture
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream',
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

    _N_ins = 2
    _N_outs = 1
        
    def _run(self):
        
        HTLaqueous, extracted = self.ins
        mixture = self.outs[0]
        
        mixture.imass['C'] = self.ins[0]._source.HTLaqueous_C
        mixture.imass['N'] = self.ins[0]._source.HTLaqueous_N
        mixture.imass['P'] = self.ins[0]._source.HTLaqueous_P + extracted.imass['P']
        mixture.imass['H2O'] = HTLaqueous.F_mass + extracted.F_mass -\
            mixture.imass['C'] - mixture.imass['N'] - mixture.imass['P'] #Represented by H2O except C, N, P
        
        mixture.T = extracted.T
        mixture.P = HTLaqueous.P
        
    def _design(self):
        pass
    
    def _cost(self):
        pass

# =============================================================================
# Struvite Precipitation
# =============================================================================

class StruvitePrecipitation(SanUnit):
    '''
    extracted_P and HTL aqueous are mixed together (Mixer) before adding MgCl2 and struvite precipitation.
    
    Model method: P recovery rate with uncertainty from literature data. If mol(N)<mol(P), add NH4Cl to mol(N):mol(P)=1:1
    
    Parameters
    ----------
    ins: Iterable (stream)
        mixture, supply_MgCl2, supply_NH4Cl
    outs: Iterable (stream)
        struvite, effluent
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream', Mg_P_ratio=1,
                 P_pre_recovery_ratio=0.95, P_in_struvite=0.127,
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.Mg_P_ratio = Mg_P_ratio
        self.P_pre_recovery_ratio = P_pre_recovery_ratio
        self.P_in_struvite = P_in_struvite

    _N_ins = 3
    _N_outs = 2
        
    def _run(self):
        
        mixture, supply_MgCl2, supply_NH4Cl = self.ins
        struvite, effluent = self.outs
        
        if mixture.imass['P']/31 > mixture.imass['N']/14:
            supply_NH4Cl.imass['NH4Cl'] = (mixture.imass['P']/31 - mixture.imass['N']/14)*53.5 #make sure N:P >= 1:1
        
        supply_MgCl2.imass['MgCl2'] = mixture.imass['P']/31*95.211*self.Mg_P_ratio #Mg:P = 1:1
        struvite.imass['Struvite'] = mixture.imass['P']*self.P_pre_recovery_ratio/self.P_in_struvite
        
        supply_MgCl2.phase = 's'
        
        effluent.copy_like(mixture)
        effluent.imass['P'] -= struvite.imass['Struvite']*self.P_in_struvite
        effluent.imass['N'] += supply_NH4Cl.imass['NH4Cl']*14/53.5 - struvite.imass['Struvite']*self.P_in_struvite/31*14
        effluent.imass['H2O'] += (supply_MgCl2.imass['MgCl2'] + supply_NH4Cl.imass['NH4Cl'] -\
                                  struvite.imass['Struvite']*(1 - self.P_in_struvite*(1+14/31)))
        struvite.phase = 's'    
            
        struvite.T = mixture.T
        effluent.T = mixture.T
        
    @property
    def struvite_P(self):
        return self.outs[0].imass['Struvite']*self.P_in_struvite

    @property
    def struvite_N(self):
        return self.struvite_P*14/31

    def _design(self):
        pass
    
    def _cost(self):
        pass

# =============================================================================
# CHG
# =============================================================================

class CHG(SanUnit):
   
    '''
    CHG serves to reduce the COD content in the aqueous phase and produce fuel gas 
    under elevated temperature (350°C) and pressure. The outlet will be cooled down
    and separated by a flash unit.
    
    Model method: use experimental data, assume no NH3 loss for now.
    
    Parameters
    ----------
    ins: Iterable (stream)
        CHGfeed
    outs: Iterable (stream)
        CHGfuelgas, effluent
    '''
    
    auxiliary_unit_names=('heat_exchanger',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream',
                 ch4_ratio=0.244, co_ratio=0.029, co2_ratio=0.150,
                 c2h6_ratio=0.043, #fuel gas ratio are from Li 2018
                 toc_tc_ratio=0.764,
                 toc_to_gas_c_ratio=0.262,
                 CHGout_pre = 3065.7*6894.76, #Jones 2014: pressure before flash
                 eff_T=60+273.15, #Jones 2014: temperature after cooler
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.ch4_ratio = ch4_ratio
        self.co_ratio = co_ratio
        self.co2_ratio = co2_ratio
        self.c2h6_ratio = c2h6_ratio
        self.toc_tc_ratio = toc_tc_ratio
        self.toc_to_gas_c_ratio = toc_to_gas_c_ratio
        self.CHGout_pre = CHGout_pre
        self.eff_T = eff_T
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'{ID}_hx', ins=hx_in, outs=hx_out)
        
    _N_ins = 1
    _N_outs = 1
        
    def _run(self):
        
        CHGfeed = self.ins[0]
        CHGout = self.outs[0]
        
        for i in (self.ch4_ratio, self.co_ratio, self.co2_ratio, self.c2h6_ratio,
                  1 - self.ch4_ratio - self.co_ratio - self.co2_ratio - self.c2h6_ratio):
            if i < 0: i = 0
        
        CHGfuel_gas_mass = CHGfeed.imass['C']*self.toc_tc_ratio*self.toc_to_gas_c_ratio/(self.\
            ch4_ratio*12/16 + self.co_ratio*12/28 + self.co2_ratio*12/44 + self.c2h6_ratio*24/30)

        CHGfuelgas_composition = {
            'CH4':self.ch4_ratio,
            'CO':self.co_ratio,
            'CO2':self.co2_ratio,
            'C2H6':self.c2h6_ratio,
            'H2':1-self.ch4_ratio-self.co_ratio-self.co2_ratio-self.c2h6_ratio
            }
        
        for name,ratio in CHGfuelgas_composition.items():
            CHGout.imass[name] = CHGfuel_gas_mass*ratio
        
        
        # CHGout.imass['C'] *= (1 - self.toc_tc_ratio*self.toc_to_gas_c_ratio)
        # CHGout.imass['H2O'] = CHGfeed.F_mass - CHGfuel_gas_mass - CHGout.imass['C'] -\
        #     CHGout.imass['N'] - CHGout.imass['P']
            
        CHGout.imass['H2O'] = CHGfeed.F_mass - CHGfuel_gas_mass
        
        # CHGout.phase = 'l'
            
        CHGout.P = self.CHGout_pre
        
        for stream in self.outs: stream.T = self.eff_T
        

        
    def _design(self):
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from(self.outs)
        hx_outs0.mix_from(self.outs)
        hx_ins0.T = self.ins[0].T
        hx.T = hx_outs0.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
    
    def _cost(self):
        pass
    
# =============================================================================
# Membrane Distillation
# =============================================================================

class MembraneDistillation(SanUnit):
    
    '''
    Membrane distillation will be modeled using vapor-liquid-liquid equilibrium later.
    Here is simplified.
    
    Model method: 
        1. Feed pH = 10, permeate pH =1.5 (0.5 M H2SO4)
        2. All N in the feed are NH4+/NH3 (Jones PNNL 2014)
        3. 95% NH3 in feed can be transfered to permeate (assume 95% for now)
        4. All NH3 in permeate can form (NH4)2SO4 (which makes sense since just water evaporates)
        5. _design and _cost refer to A.A. et al., Membrane distillation: A comprehensive review
    
    
# =============================================================================
#     ############# add NaOH somewhere there or in other units
# =============================================================================
    
    
    
    Parameters
    ----------
    ins: Iterable (stream)
        influent, acid
    outs: Iterable (stream)
        ammoniasulfate, ww
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream',
                 N_S_ratio=2, pH=10, ammonia_transfer_ratio=0.95,
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.N_S_ratio = N_S_ratio
        self.pH = pH
        self.ammonia_transfer_ratio = ammonia_transfer_ratio

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        influent, acid = self.ins
        ammoniasulfate, ww = self.outs
        
        acid.imass['H2SO4'] = influent.imass['N']/14/self.N_S_ratio*98
        acid.imass['H2O'] = acid.imass['H2SO4']*1000/98/0.5*1.05 - acid.imass['H2SO4']
        
        
        pKa = 9.26 #ammonia pKa
        ammonia_to_ammonium = 10**(-pKa)/10**(-self.pH)
        ammonia_in_feed = influent.imass['N']/14*ammonia_to_ammonium/(1 + ammonia_to_ammonium)*17

        ammoniasulfate.imass['NH42SO4'] = ammonia_in_feed*self.ammonia_transfer_ratio/34*132
        ammoniasulfate.imass['H2O'] = acid.imass['H2O']
        ammoniasulfate.imass['H2SO4'] = acid.F_mass + ammoniasulfate.imass['NH42SO4']/132*26 -\
            ammoniasulfate.imass['NH42SO4'] - ammoniasulfate.imass['H2O']
        ww.copy_like(influent)
        ww.imass['N'] -= influent.imass['N']*self.ammonia_transfer_ratio*ammonia_to_ammonium/(1 + ammonia_to_ammonium)
        ammoniasulfate.T = acid.T
        ww.T = acid.T
        ww.P = acid.P
        
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
    
    Model method: use experimental data.
    
    Parameters
    ----------
    ins: Iterable (stream)
    biocrude,hydrogen_gas
    outs: Iterable (stream)
    HTaqueous, fuel_gas, gasoline, diesel, heavy_oil
    '''
    
    auxiliary_unit_names=('heat_exchanger',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream',
                 biooil_ratio=0.77, gas_ratio=0.07,#Jones et al., 2014
                 biocrude_h2_ratio = 18.14,
                 gasoline_ratio = 0.103, heavy_oil_ratio = 0.113, #see spreadsheet for calculation
                 co_ratio=0.128, co2_ratio=0.007, c2h6_ratio=0.188,
                 c3h8_ratio=0.107, c4h10_ratio=0.09,
                 biooil_C_ratio=0.855, biooil_N_ratio=0.01,
                 HTaqueous_pre = 55*6894.76, #Jones 2014: 55 psia
                 fuel_gas_pre = 20*6894.76,
                 gasoline_pre = 25*6894.76,
                 diesel_pre = 18.7*6894.76,
                 heavy_oil_pre = 18.7*6894.76,
                 HTaqueous_T=47+273.15,
                 fuel_gas_T=44+273.15, #Jones 2014
                 gasoline_T=109+273.15,
                 diesel_T=265+273.15,
                 heavy_oil_T=381+273.15,
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.biooil_ratio = biooil_ratio
        self.gas_ratio = gas_ratio
        self.biocrude_h2_ratio = biocrude_h2_ratio
        self.gasoline_ratio = gasoline_ratio
        self.heavy_oil_ratio = heavy_oil_ratio  
        self.c4h10_ratio = c4h10_ratio
        self.co_ratio = co_ratio
        self.co2_ratio = co2_ratio
        self.c2h6_ratio = c2h6_ratio
        self.c3h8_ratio = c3h8_ratio
        self.biooil_C_ratio = biooil_C_ratio
        self.biooil_N_ratio = biooil_N_ratio
        self.HTaqueous_pre = HTaqueous_pre
        self.fuel_gas_pre = fuel_gas_pre
        self.gasoline_pre = gasoline_pre
        self.diesel_pre = diesel_pre
        self.heavy_oil_pre = heavy_oil_pre
        self.HTaqueous_T = HTaqueous_T
        self.fuel_gas_T = fuel_gas_T
        self.gasoline_T = gasoline_T
        self.diesel_T = diesel_T
        self.heavy_oil_T = heavy_oil_T
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'{ID}_hx', ins=hx_in, outs=hx_out)

    _N_ins = 2
    _N_outs = 5
        
    def _run(self):
        
        biocrude,hydrogen_gas = self.ins
        HTaqueous, fuel_gas, gasoline, diesel, heavy_oil = self.outs
        
        # hydrogen_gas.imass['H2'] = biocrude.F_mass/self.biocrude_h2_ratio
        hydrogen_gas.phase = 'g'
        
        biooil_total_mass = biocrude.F_mass*self.biooil_ratio
        
        for i in (self.gasoline_ratio, self.heavy_oil_ratio):
            if i < 0: i = 0
            
        gasoline.imass['Gasoline'] = biooil_total_mass*self.gasoline_ratio
        diesel.imass['Diesel'] = biooil_total_mass*(1 - self.gasoline_ratio - self.heavy_oil_ratio)
        heavy_oil.imass['Heavy_oil'] = biooil_total_mass*self.heavy_oil_ratio
        
        gas_mass = biocrude.F_mass*self.gas_ratio
        
        for i in (self.c4h10_ratio, self.co_ratio, self.co2_ratio, self.c2h6_ratio, self.c3h8_ratio,
                    1 - self.c4h10_ratio - self.co_ratio - self.co2_ratio - self.c2h6_ratio -\
                    self.c3h8_ratio):
            if i < 0: i = 0
    
        fuelgas_HT_composition = {
            'CO':self.co_ratio,
            'CO2':self.co2_ratio,
            'C2H6':self.c2h6_ratio,
            'C3H8':self.c3h8_ratio,
            'C4H10':self.c4h10_ratio,
            'CH4':1 - self.c4h10_ratio - self.co_ratio - self.co2_ratio - self.c2h6_ratio -\
                self.c3h8_ratio
            }
            
        for name,ratio in fuelgas_HT_composition.items():
            fuel_gas.imass[name] = gas_mass*ratio
        
        HTaqueous.imass['HTaqueous'] = biocrude.F_mass + hydrogen_gas.F_mass - biooil_total_mass - gas_mass
        #HTaqueous is liquid waste from HT
        
        fuel_gas.phase = 'g'
        
        HTaqueous.P = self.HTaqueous_pre
        fuel_gas.P = self.fuel_gas_pre
        gasoline.P = self.gasoline_pre
        diesel.P = self.diesel_pre
        heavy_oil.P = self.heavy_oil_pre
        
        HTaqueous.T = self.HTaqueous_T
        fuel_gas.T = self.fuel_gas_T
        gasoline.T = self.gasoline_T
        diesel.T = self.diesel_T
        heavy_oil.T = self.heavy_oil_T
        
    @property
    def biooil_C(self):
        return min(self.outs[2].F_mass*self.biooil_C_ratio,self.ins[0]._source.biocrude_C)
        
    @property
    def HTfuel_gas_C(self):
        HTfuel_gas_C = 0
        fuelgas_carbo_ratio = {
            'C4H10':60/72,
            'CO':12/28,
            'CO2':12/44,
            'C2H6':24/30,
            'C3H8':36/44,
            'CH4':12/16
              }
        for name, ratio in fuelgas_carbo_ratio.items():
            HTfuel_gas_C+=self.outs[1].imass[name]*ratio
        return min(HTfuel_gas_C, self.ins[0]._source.biocrude_C-self.biooil_C)
    
    @property
    def biooil_N(self):
        return self.outs[2].F_mass*self.biooil_N_ratio

    biocrude_C_ratio = HTL.biocrude_C_ratio

    @property
    def HTaqueous_C(self):
        return max(0,self.ins[0]._source.biocrude_C - self.HTfuel_gas_C - self.biooil_C)

    @property
    def HTaqueous_N(self):
        return self.ins[0]._source.biocrude_N - self.biooil_N

    def _design(self):
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from(self.outs)
        hx_outs0.mix_from(self.outs)
        hx_ins0.T = self.ins[0].T
        hx.T = hx_outs0.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
    
    def _cost(self):
        pass
   
   
class HC(SanUnit):
    
    '''
    Hydrocracking further cracks down heavy part in HT biooil to diesel and gasoline.
    
    Model method: use experimental data.
    
    Parameters
    ----------
    ins: Iterable (stream)
    heavy_oil, hydrogen_gas
    outs: Iterable (stream)
    gasoline, diesel, off_gas
    '''
    
    auxiliary_unit_names=('heat_exchanger',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream',
                 heavy_oil_h2_ratio=45.17,
                 hydrogen_gas_T=128+273.15, #Jones 2014 #H2 heating and pressurization?
                 hydrogen_gas_P=1039.7*6894.76,
                 gasoline_hc_ratio=0.303,off_gas_ratio=0.0394,
                 gasoline_pre = 20*6894.76,
                 diesel_pre = 24*6894.76,
                 off_gas_pre = 20*6894.76, #Jones 2014
                 gasoline_T=60+273.15,
                 diesel_T=220+273.15,
                 off_gas_T=45+273.15, #Jones 2014
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo,init_with)
        self.heavy_oil_h2_ratio = heavy_oil_h2_ratio
        self.hydrogen_gas_T = hydrogen_gas_T
        self.hydrogen_gas_P = hydrogen_gas_P
        self.gasoline_hc_ratio = gasoline_hc_ratio #see spreadsheet for calculation
        self.off_gas_ratio = off_gas_ratio
        self.gasoline_pre = gasoline_pre
        self.diesel_pre = diesel_pre
        self.off_gas_pre = off_gas_pre  
        self.gasoline_T = gasoline_T
        self.diesel_T = diesel_T
        self.off_gas_T = off_gas_T
        
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'{ID}_hx', ins=hx_in, outs=hx_out)
        
    _N_ins = 2
    _N_outs = 3
        
    def _run(self):
        
        heavy_oil, hydrogen_gas = self.ins
        gasoline, diesel, off_gas = self.outs
        
        hydrogen_gas.imass['H2'] = heavy_oil.F_mass/self.heavy_oil_h2_ratio
        hydrogen_gas.phase = 'g'
        hydrogen_gas.T = self.hydrogen_gas_T
        hydrogen_gas.P = self.hydrogen_gas_P
        
        for i in [self.gasoline_hc_ratio, self.off_gas_ratio]:
            if i < 0: i = 0
        
        gasoline.imass['Gasoline'] = (heavy_oil.F_mass + hydrogen_gas.F_mass)*self.gasoline_hc_ratio

        diesel.imass['Diesel'] = (heavy_oil.F_mass + hydrogen_gas.F_mass)*\
            (1 - self.gasoline_hc_ratio - self.off_gas_ratio)
        
        off_gas.imass['CO2'] = (heavy_oil.F_mass + hydrogen_gas.F_mass)*self.off_gas_ratio
        
        off_gas.phase = 'g'

        gasoline.P = self.gasoline_pre
        diesel.P = self.diesel_pre
        off_gas.P = self.off_gas_pre
        
        gasoline.T = self.gasoline_T
        diesel.T = self.diesel_T
        off_gas.T = self.off_gas_T
        
    def _design(self):
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from(self.outs)
        hx_outs0.mix_from(self.outs)
        hx_ins0.T = self.ins[0].T
        hx.T = hx_outs0.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
    
    def _cost(self):
        pass
    
    

class Acidsplitter(SanUnit):
    
    '''
    Hydrocracking further cracks down heavy part in HT biooil to diesel and gasoline.
    
    Model method: use experimental data.
    
    Parameters
    ----------
    ins: Iterable (stream)
    acid_in
    outs: Iterable (stream)
    acid_for_N, acid_for_P
    '''
    
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='Stream',
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo,init_with)

    _N_ins = 1
    _N_outs = 2
        
    def _run(self):
        
        acid_in = self.ins[0]
        acid_for_N, acid_for_P = self.outs
        
        acid_in.mix_from((acid_for_N, acid_for_P))
    
    def _design(self):
        pass
    
    def _cost(self):
        pass
        







class Flash(SanUnit, FL):

    auxiliary_unit_names = FL.auxiliary_unit_names
    _units = FL._units
    _max_agile_design = FL._max_agile_design
    _F_BM_default = FL._F_BM_default
    _graphics = FL._graphics
    _N_outs = FL._N_outs

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', *,
                 V=None, T=None, Q=None, P=None, y=None, x=None,
                 vessel_material='Carbon steel',
                 vacuum_system_preference='Liquid-ring pump',
                 has_glycol_groups=False,
                 has_amine_groups=False,
                 vessel_type=None,
                 holdup_time=15,
                 surge_time=7.5,
                 has_mist_eliminator=False):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with=init_with)
        self._load_components()
        
        #: Enforced molar vapor fraction
        self.V = V
        
        #: Enforced operating temperature (K)
        self.T = T
        
        #: [array_like] Molar composition of vapor (for binary mixture)
        self.y = y
        
        #: [array_like] Molar composition of liquid (for binary mixture)
        self.x = x
        
        #: Enforced duty (kJ/hr)
        self.Q = Q
        
        #: Operating pressure (Pa)
        self.P = P
        
        #: [str] Vessel construction material
        self.vessel_material = vessel_material

        #: [str] If a vacuum system is needed, it will choose one according to this preference.
        self.vacuum_system_preference = vacuum_system_preference
        
        #: [bool] True if glycol groups are present in the mixture
        self.has_glycol_groups = has_glycol_groups
        
        #: [bool] True if amine groups are present in the mixture
        self.has_amine_groups = has_amine_groups
        
        #: [str] 'Horizontal', 'Vertical', or 'Default'
        self.vessel_type = vessel_type
        
        #: [float] Time it takes to raise liquid to half full (min)
        self.holdup_time = holdup_time
        
        #: [float] Time it takes to reach from normal to maximum liquied level (min)
        self.surge_time = surge_time
        
        #: [bool] True if using a mist eliminator pad
        self.has_mist_eliminator = has_mist_eliminator
        
    def _load_components(self):
        self._multi_stream = ms = MultiStream(None, thermo=self.thermo)
        self.heat_exchanger = HXutility(None, (None,), ms, thermo=self.thermo) 
        
    def reset_cache(self, isdynamic=None):
        self._multi_stream.reset_cache()
        self.heat_exchanger.reset_cache()
        
    @property
    def P(self):
        """Operating pressure (Pa)."""
        return self._P
    @P.setter
    def P(self, P):
        if P and P < 101325 and not self.power_utility:
            self.power_utility = PowerUtility()
        self._P = P

    @property
    def vapor(self):
        """Outlet vapor stream (equivalent to outs[0])."""
        return self._outs[0]
    @vapor.setter
    def vapor(self, vapor):
        self._outs[0] = vapor
    
    @property
    def liquid(self):
        """Outlet liquid stream (equivalent to outs[1])."""
        return self._outs[1]
    @liquid.setter
    def liquid(self, liquid):
        self._outs[1] = liquid

    def _default_vessel_type(self):
        vap, liq = self.outs
        F_mass_vap = vap.F_mass
        F_mass_liq = liq.F_mass 
        return 'Vertical' if F_mass_vap / F_mass_liq > 0.1 else 'Horizontal'

    def _run(self):
        separations.vle(self.ins[0], *self.outs, self.T, self.P, self.V, 
                        self.Q, self.x, self.y, self._multi_stream)
            
    def _design(self):
        vap, liq = self.outs
        self.no_vessel_needed = vap.isempty() or liq.isempty()
        if self.no_vessel_needed:
            self.design_results.clear()
        else:
            vessel_type = self.vessel_type
            if vessel_type == 'Vertical': 
                args = self._vertical_vessel_pressure_diameter_and_length()
            elif vessel_type == 'Horizontal': 
                args = self._horizontal_vessel_pressure_diameter_and_length()
            else: raise RuntimeError('unknown vessel type') # pragma: no cover
            self.design_results.update(
                self._vessel_design(*args)
            )
        if self.Q == 0.:
            self.heat_exchanger._setup() # Removes results
        else:
            self.heat_exchanger.simulate_as_auxiliary_exchanger(self.ins, [self._multi_stream])

    def _cost(self):
        D = self.design_results
        if not self.no_vessel_needed:
            self.baseline_purchase_costs.update(
                self._vessel_purchase_cost(D['Weight'], D['Diameter'], D['Length'])
            )
            self._cost_vacuum()

    def _cost_vacuum(self):
        P = self.P
        if not P or P > 101320: 
            self.vacuum_system = None
        else:
            Design = self.design_results
            R = Design['Diameter'] * 0.5
            volume = 0.02832 * np.pi * Design['Length'] * R * R # Volume ft3 to m3
            self.vacuum_system = bst.VacuumSystem(
                self, self.vacuum_system_preference, vessel_volume=volume,
            )

    def _design_parameters(self):
        # Retrieve run_args and properties
        vap, liq = self._outs
        rhov = vap.get_property('rho', 'lb/ft3')
        rhol = liq.get_property('rho', 'lb/ft3')
        P = liq.get_property('P', 'psi')  # Pressure (psi)

        vessel_type = self.vessel_type
        Th = self.holdup_time
        Ts = self.surge_time
        has_mist_eliminator = self.has_mist_eliminator

        # Calculate the volumetric flowrate
        Qv = vap.get_total_flow('ft^3 / s')
        Qll = liq.get_total_flow('ft^3 / min')

        # Calculate Ut and set Uv
        K = design.compute_Stokes_law_York_Demister_K_value(P)

        # Adjust K value
        if not has_mist_eliminator and vessel_type == 'Vertical': K /= 2

        # Adjust for amine or glycol groups:
        if self.has_glycol_groups: K *= 0.6
        elif self.has_amine_groups: K *= 0.8

        Ut = K*((rhol - rhov) / rhov)**0.5
        Uv = 0.75*Ut

        # Calculate Holdup and Surge volume
        Vh = Th*Qll
        Vs = Ts*Qll
        return rhov, rhol, P, Th, Ts, has_mist_eliminator, Qv, Qll, Ut, Uv, Vh, Vs

    def _vertical_vessel_pressure_diameter_and_length(self):
        rhov, rhol, P, Th, Ts, has_mist_eliminator, Qv, Qll, Ut, Uv, Vh, Vs = self._design_parameters()

        # Calculate internal diameter, Dvd
        Dvd = (4.0*Qv/(pi*Uv))**0.5
        if has_mist_eliminator:
            D = design.ceil_half_step(Dvd + 0.4)
        else:
            D = design.ceil_half_step(Dvd)

        # Obtaining low liquid level height, Hlll
        Hlll = 0.5
        if P < 300:
            Hlll = 1.25

        # Calculate the height from Hlll to Normal liq level, Hnll
        Hh = Vh/(pi/4.0*Dvd**2)
        if Hh < 1.0: Hh = 1.0

        # Calculate the height from Hnll to  High liq level, Hhll
        Hs = Vs/(pi/4.0*Dvd**2)
        if Hs < 0.5: Hs = 0.5

        # Calculate dN
        Qm = Qll + Qv
        lamda = Qll/Qm
        rhoM = rhol*lamda + rhov*(1-lamda)
        dN = (4*Qm/(pi*60.0/(rhoM**0.5)))**0.5
        dN = design.ceil_half_step(dN)

        # Calculate Hlin, assume with inlet diverter
        Hlin = 1.0 + dN

        # Calculate the vapor disengagement height
        Hv = 0.5*Dvd
        Hv2 = (2.0 if has_mist_eliminator else 3.0) + dN/2.0 # pragma: no cover
        if Hv2 < Hv: Hv = Hv2
        Hv = Hv

        # Calculate total height, Ht
        Hme = 1.5 if has_mist_eliminator else 0.0
        Ht = Hlll + Hh + Hs + Hlin + Hv + Hme
        Ht = design.ceil_half_step(Ht)

        # Find maximum and normal liquid level
        # Hhll = Hs + Hh + Hlll
        # Hnll = Hh + Hlll

        return P, D, Ht
        
    def _horizontal_vessel_pressure_diameter_and_length(self):
        rhov, rhol, P, Th, Ts, has_mist_eliminator, Qv, Qll, Ut, Uv, Vh, Vs = self._design_parameters()

        # Initialize LD
        if P > 0 and P <= 264.7:
            LD = 1.5/250.0*(P-14.7)+1.5
        elif P > 264.7 and P <= 514.7: # pragma: no cover
            LD = 1.0/250.0*(P-14.7)+2.0
        elif P > 514.7: # pragma: no cover
            LD = 5.0

        D = (4.0*(Vh+Vs)/(0.6*pi*LD))**(1.0/3.0)
        if D <= 4.0:
            D = 4.0
        else:
            D = design.ceil_half_step(D)

        for outerIter in range(50):
            At = pi*(D**2)/4.0 # Total area

            # Calculate Lower Liquid Area
            Hlll = round(0.5*D + 7.0)  
            Hlll = Hlll/12.0 # D is in ft but Hlll is in inches
            X = Hlll/D
            Y = design.HNATable(1, X)
            Alll = Y*At

            # Calculate the Vapor disengagement area, Av
            Hv = 0.2*D
            if has_mist_eliminator and Hv <= 2.0: Hv = 2.0
            elif Hv <= 1.0: Hv = 1.0
            else: Hv = design.ceil_half_step(Hv)
            Av = design.HNATable(1, Hv/D)*At
            
            # Calculate minimum length for surge and holdup
            L = (Vh + Vs)/(At - Av - Alll)
            # Calculate liquid dropout
            Phi = Hv/Uv
            # Calculate actual vapor velocity
            Uva = Qv/Av
            # Calculate minimum length for vapor disengagement
            Lmin = Uva*Phi
            Li = L
            
            for innerIter in range(50):
                if L < 0.8*Lmin: Hv += 0.5
                elif L > 1.2*Lmin:
                    if has_mist_eliminator and Hv <= 2.0: Hv = 2.0
                    elif not has_mist_eliminator and Hv <= 1.0: Hv = 1.0
                    else: Hv -= 0.5
                else: break
                Av = design.HNATable(1, Hv/D)*At
                Alll = design.HNATable(1, Hlll/D)*At
                Li = (Vh + Vs)/(At - Av - Alll)
                Phi = Hv/Uv
                Uva = Qv/Av
                Lmin = Uva*Phi
            
            L = Li
            LD = L/D
            # Check LD
            if LD < 1.2:
                if D <= 4.0: break
                else: D -= 0.5

            if LD > 7.2: # pragma: no cover
                D += 0.5
            else: break

        # Recalculate LD so it lies between 1.5 - 6.0
        while True:
            LD = L / D
            if (LD < 1.5) and D <= 4.0: L += 0.5
            elif LD < 1.5: D -= 0.5
            elif (LD > 6.0): D += 0.5
            else: break

        # # To check minimum Hv value
        # if int(has_mist_eliminator) == 1 and Hv <= 2.0:
        #     Hv = 2.0
        # if int(has_mist_eliminator) == 0 and Hv <= 1.0:
        #     Hv = 1.0

        # Calculate normal liquid level and High liquid level
        # Hhll = D - Hv
        # if (Hhll < 0.0):
        #     Hhll = 0.0
        # Anll = Alll + Vh/L
        # X = Anll/At
        # Y = HNATable(2, X)
        # Hnll = Y*D
        
        return P, D, L













