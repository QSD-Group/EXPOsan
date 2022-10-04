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

from warnings import warn
from qsdsan import SanUnit
import exposan

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
    elevated temperature and pressure. The products percentage (wt%) can be evaluatad
    using revised MCA model (Li et al., 2017, Leow et al., 2018) with known sludge
    composition (protein%, lipid%, and carbohydrate%)
    
    Parameters
    ----------
    ins: Iterable (stream)
        Dewatered biosolid or sludge
    outs: Iterable (stream)
        biocrude, aqueous, offgas, biochar
    T: Float
        Temperature in HTL (k)
    P: Float
        Pressure at room temperature (N2, psi)
        #pressure at high temperature can be found in Jones et al., 2014 and SS et al., 2016
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',
                 T=300+273.15,P=14.7,**kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.T = T
        self.P = P

    _N_ins = 1
    _N_outs = 4
        
    def _run(self):
        
        dewatered_sludge=self.ins[0]
        biocrude,aqueous,offgas,biochar=self.outs
        
        dewatered_sludge_afdw = dewatered_sludge.imass['Sludge_lipid']\
            + dewatered_sludge.imass['Sludge_protein']\
            + dewatered_sludge.imass['Sludge_carbo']
                               
        lipid_ratio = dewatered_sludge.imass['Sludge_lipid']/dewatered_sludge_afdw
        protein_ratio = dewatered_sludge.imass['Sludge_protein']/dewatered_sludge_afdw
        carbo_ratio = dewatered_sludge.imass['Sludge_carbo']/dewatered_sludge_afdw
        
        if lipid_ratio+protein_ratio+carbo_ratio != 1:
            warn('The sum of the sludge composition is not 1')
        
        #Revised MCA model
        biocrude.imass['Biocrude'] = (0.846*lipid_ratio + 0.445*protein_ratio\
            + 0.205*carbo_ratio) * dewatered_sludge_afdw
        aqueous.imass['HTLaqueous'] = (0.154*lipid_ratio + 0.481*protein_ratio)\
            * dewatered_sludge_afdw
        offgas.imass['CO2'] = (0.074*protein_ratio + 0.418*carbo_ratio)\
            * dewatered_sludge_afdw
        offgas.phase='g'
        biochar.imass['Biochar'] = 0.377*carbo_ratio * dewatered_sludge_afdw

    def _design(self):
        pass
    
    def _cost(self):
        pass

class HT(SanUnit):
    
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature and pressure
    to produce upgraded biooil. Co-products include fuel gas and char. The amount of
    biooil and fuel gas can be estimated using values from Li et al., 2018.
    The amount of char can be calculated based on mass closure.
    
    Parameters
    ----------
    ins: Iterable (stream)
        biocurde
    outs: Iterable (stream)
        biooil, fuelgas_HT, char
    T: Float
        Temperature in HT (k)
    P: Float
        Pressure at room temperature (H2, psi)
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',
                 T = 405+273.15, P=960, **kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.T = T
        self.P = P

    _N_ins = 1
    _N_outs = 3
        
    def _run(self):
        
        biocrude = self.ins[0]
        biooil, fuelgas_HT, char = self.outs
        
        massin = biocrude.imass['Biocrude']
        biooil.imass['Biooil'] = massin * 0.581
        fuelgas_HT_composition = {
            'CH4':0.48,
            'CO':0.128,
            'CO2':0.007,
            'C2H6':0.188,
            'C3H8':0.107,
            'C5H12':0.09
            }
        for name,ratio in fuelgas_HT_composition.items():
            fuelgas_HT.imass[name] = massin * 0.047 * ratio
        char.imass['Char'] = massin - biooil.imass['Biooil'] - massin * 0.047
        fuelgas_HT.phase='g'
        
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
    T: Float
        Temperature in AcidExtraction (k)
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',
                  T = 80+273.15, **kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.T = T

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        biochar, acid_P = self.ins
        residual, extracted_P = self.outs
        
        cmps = exposan.htl._components.create_components()
        residual.imass['Residual'] = biochar.imass['Biochar']-biochar.imass['Biochar']*\
            cmps.Biochar.i_P*0.95
        extracted_P.imass['H3PO4'] = biochar.imass['Biochar']*cmps.Biochar.i_P*0.95/31*98
        
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
        
        aqueous, extracted_P = self.ins
        mixture = self.outs[0]
        
        mixture.imass['Mixture'] = aqueous.imass['HTLaqueous']+extracted_P.imass['H3PO4']
          
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
        struvite, chgfeed
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',
                  **kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        mixture, supply_MgCl2 = self.ins
        struvite, chgfeed = self.outs
        
        cmps = exposan.htl._components.create_components()
        struvite.imass['Struvite'] = mixture.imass['Mixture']*cmps.Mixture.i_P*0.95/0.127
        struvite.phase='s'
        chgfeed.imass['CHGfeed'] = mixture.imass['Mixture']+supply_MgCl2.imass['MgCl2']\
            - struvite.imass['Struvite']
          
    def _design(self):
        pass
    
    def _cost(self):
        pass
    
class CHG(SanUnit):
   
    '''
    CHG serves to reduce the COD content in the aqueous phase and produce fuel gas 
    under elevated temperature and pressure.
    
    Parameters
    ----------
    ins: Iterable (stream)
        chgfeed
    outs: Iterable (stream)
        fuelgas_CHG, chgeffluent
    T: Float
        Temperature in CHG (k)  
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',
                  T = 350+273.15, **kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)
        self.T = T

    _N_ins = 1
    _N_outs = 2
        
    def _run(self):
        
        chgfeed = self.ins[0]
        fuelgas_CHG, chgeffluent = self.outs
        
        cmps = exposan.htl._components.create_components()
        fuelgas_CHG_composition = {
            'CH4':0.244,
            'CO':0.029,
            'CO2':0.150,
            'C2H6':0.043,
            'H2':0.534
            }
        for name,ratio in fuelgas_CHG_composition.items():
            fuelgas_CHG.imass[name] = chgfeed.imass['CHGfeed']*cmps.CHGfeed.i_C*0.764*0.262/0.3 * ratio
        
        chgeffluent.imass['CHGeffluent'] = chgfeed.imass['CHGfeed'] -\
            chgfeed.imass['CHGfeed']*cmps.CHGfeed.i_C*0.764*0.262/0.3
       
        
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
        chgeffluent, acid_N
    outs: Iterable (stream)
        ammoniasulfate, ww
    '''
    
    def __init__(self,ID='',ins=None,outs=(),thermo=None,init_with='WasteStream',
                  **kwargs):
        SanUnit.__init__(self,ID,ins,outs,thermo,init_with)

    _N_ins = 2
    _N_outs = 2
        
    def _run(self):
        
        chgeffluent, acid_N = self.ins
        ammoniasulfate, ww = self.outs
        
        cmps = exposan.htl._components.create_components()
        ammoniasulfate.imass['NH42SO4'] = chgeffluent.imass['CHGeffluent']*\
            cmps.CHGeffluent.i_N*0.9/14*132
        ww.imass['WW'] = chgeffluent.imass['CHGeffluent']+acid_N.imass['H2SO4']-\
            ammoniasulfate.imass['NH42SO4']
        
    def _design(self):
        pass
    
    def _cost(self):
        pass