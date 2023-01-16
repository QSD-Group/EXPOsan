# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 14:09:53 2022

@author: joy_c
"""
from qsdsan.utils import auom
from chemicals.elements import molecular_weight as get_mw
from math import pi

__all__ = ('encap_items',
           'encap_material_cost',
           'E_heat_temed',
           'encap_material_input')

#%%
encap_items = {
    'EG': ['ethylene glycol production', 'kg'],
    'MAA': ['methacrylic acid production', 'kg'],
    'PAM': ['polyacrylamide production', 'kg'],
    'FMD': ['oxidation of methanol', 'kg'],
    'CuSO4': ['copper sulfate production', 'kg'],
    'H2SO4': ['sulfuric acid production', 'kg'],
    'DMA': ['dimethylamine production', 'kg'],
    'EDCl': ['ethylene dichloride production', 'kg'],
    'KOH': ['potassium hydroxide production', 'kg'],
    'NaPS': ['sodium persulfate production', 'kg'],
    'GAC': ['activated carbon production, granular from hard coal', 'kg']
    }

# d_bead = 0.4*auom('inch').conversion_factor('m')  # bead diameter
# v_bead = 4/3 * pi * (d_bead/2)**3 /2   # half sphere, in m3

# # per batch beads
# ingredients = {
#     'PEGDMA_1000': 4.5/1000, # g to kg
#     'BIS': 225/1e6,          # mg to kg
#     'TEMED': 60/1e6,         # uL to L
#     'APS': 30/1e6,           # mg to kg
#     'PAC': 0.3/1000          # g to kg
#     }

unit_price = {
    'PEGDMA_1000': 1017.00,            # USD/kg; https://www.polysciences.com/default/polyethylene-glycol-dimethacrylate-pegdma-1000
    'BIS': 137/0.5,                    # USD/kg; https://www.sigmaaldrich.com/US/en/product/sial/146072
    'TEMED': 169.00,                   # USD/L;  https://www.sigmaaldrich.com/US/en/product/mm/808742
    'APS': 237/2.5,                    # USD/kg; https://www.sigmaaldrich.com/US/en/product/sigald/215589    
    'PAC': 393/5                       # USD/kg; https://www.sigmaaldrich.com/US/en/product/sigald/161551
    }

# USD/m3 beads
def encap_material_cost(V_beads, n_bead=45, d_bead=0.4, PEGDMA_1000=4.5e-3, 
                        BIS=2.25e-4, TEMED=6e-5, APS=3e-5, PAC=3e-4, 
                        unit_prices=[]):
    '''
    Calculate material cost of the specified amount of encapsulation beads, 
    based on `recipe <https://docs.google.com/document/d/1LfMMFil2SWCBZlhBxlHdP86P05nl_xenQlGQY-Ld1IE/edit>`.

    Parameters
    ----------
    V_beads : float
        Volume of encapsulation beads to produce, in m3.
    n_bead : int, optional
        Number of beads per batch according to the recipe. The default is 45.
    d_bead : float, optional
        Bead diameter [inch], half spheres according to the recipe. The default is 0.4.
    PEGDMA_1000 : float, optional
        Amount of PEGDMA_1000 input per batch, in kg. The default is 4.5e-3.
    BIS : float, optional
        Amount of BIS input per batch, in kg. The default is 2.25e-4.
    TEMED : float, optional
        Amount of TEMED input per batch, in kg. The default is 6e-5.
    APS : float, optional
        Amount of APS input per batch, in kg. The default is 3e-5.
    PAC : float, optional
        Amount of PAC input per batch, in kg. The default is 3e-4.
    unit_price : array-like, optional
        Unit price [USD] of the materials. Must in same order as previous inputs.
    '''
    d_bead *= auom('inch').conversion_factor('m')
    v_bead = 4/3 * pi * (d_bead/2)**3 /2
    n_batch = V_beads / (n_bead*v_bead)
    price = unit_prices or unit_price.values()
    return sum(x * p for x, p in zip([PEGDMA_1000, BIS, TEMED, APS, PAC], price)) * n_batch

#%%
mw_pegdma = get_mw(dict(C=46, H=86, O=22))
mw_eg = get_mw(dict(C=2, H=6, O=2))
mw_maa = get_mw(dict(C=4, H=6, O=2))

def make_pegdma(m):
    '''
    Estimate the mass [kg] of ethylene glycol and mathecrylic acid needed to produce 
    specified mass [kg] of PEGDMA_1000. Assume overall 90% yield. Stoichiometry:  
    [n] C2H6O2 + [2] C4H6O2 -> (C2H4O)n-C8H10O3 + [n+1] H2O.
    
    PEGDMA_1000 has molecular weight of roughly 1000. So n = 19, formular = C46H86O22
    '''
    n_pegdma = m/mw_pegdma
    m_eg = n_pegdma * 19 * mw_eg
    m_maa = n_pegdma * 2 * mw_maa
    return m_eg/0.9, m_maa/0.9

mw_bis = get_mw(dict(C=7, H=10, N=2, O=2))
mw_am = get_mw(dict(C=3, H=5, N=1, O=1))
mw_fmd = get_mw(dict(C=1, H=2, O=1))
mw_CuCl2 = get_mw(dict(Cu=1, Cl=2))
mw_CuSO4 = get_mw(dict(Cu=1, S=1, O=4))
mw_H2SO4 = get_mw(dict(H=2, S=1, O=4))

def make_bis(m):
    '''
    Estimate the mass of materials needed to produce specified mass [kg] of BIS. 
    Stoichiometry and yield based on example 5 in patent US2475846A.
    '''
    n_batch = m/mw_bis/0.15      # Per batch yields 0.25 kmol * 60% = 0.15 kmol of BIS
    # input per batch in kg
    m_pam = 0.5*mw_am
    m_fmd = 0.25*mw_fmd
    #!!! assume CuSO4 and H2SO4 reused for 500 batches on average
    n_reuse = 500
    m_CuSO4 = 0.4/mw_CuCl2*mw_CuSO4/n_reuse
    m_H2SO4 = 18.9/1.18*3*mw_H2SO4/n_reuse  # 6N sulfuric acid (i.e., 3M H2SO4) has a density of 1.18 kg/L
    return m_pam*n_batch, m_fmd*n_batch, m_CuSO4*n_batch, m_H2SO4*n_batch

mw_edcl = get_mw(dict(C=2, H=4, Cl=2))
mw_dma = get_mw(dict(C=2, H=7, N=1))
mw_temed = get_mw(dict(C=6, H=16, N=2))
mw_koh = get_mw(dict(K=1, O=1, H=1))

# Per hour temed production (https://patentimages.storage.googleapis.com/0f/02/b8/c757cf84169e05/US4053516.pdf) 
Q_water = 700 * 4182 * (165 - 20)  # in J, assume ambient temperature = 20 C; assume water and heat is recycled, so 100% heat efficiency
Q_edcl = 250 * 1297.5 * (165 - 20) # in J
m_temed = 250/mw_edcl*mw_temed

def make_temed(m):
    '''
    Estimate the mass of materials and heat energy needed to produce specified mass [kg]
    of TEMED based on example 1 in patent US4053516A.
    
    Assume (with recycling) 20% loss of dimethylamine, 100% conversion of ethylene dichloride,
    complete neutralization of HCl with KOH
    '''
    m_edcl = 250 * m/m_temed
    m_dma = m/mw_temed * 2 / (1-0.2) * mw_dma
    m_koh = m/mw_temed * 2 * mw_koh
    return m_dma, m_edcl, m_koh

E_heat_temed = lambda m: m/m_temed * (Q_water+Q_edcl) * 1e-6

mw_APS = get_mw(dict(N=2, H=8, S=2, O=8))
mw_NaPS = get_mw(dict(Na=2, S=2, O=8))

def sub_APS(m):
    '''1-to-1 (by mol) substitution of ammonium persulfate (APS) with sodium persulfate.'''
    return m/mw_APS*mw_NaPS

def encap_material_input(V_beads, n_bead=45, d_bead=0.4, PEGDMA_1000=4.5e-3, 
                         BIS=2.25e-4, TEMED=6e-5, APS=3e-5, PAC=3e-4):
    '''
    Calculate amount of raw materials needed to produce specified amount of 
    encapsulation beads, based on `recipe <https://docs.google.com/document/d/1LfMMFil2SWCBZlhBxlHdP86P05nl_xenQlGQY-Ld1IE/edit>`.

    Parameters
    ----------
    V_beads : float
        Volume of encapsulation beads to produce, in m3.
    n_bead : int, optional
        Number of beads per batch according to the recipe. The default is 45.
    d_bead : float, optional
        Bead diameter [inch], half spheres according to the recipe. The default is 0.4.
    PEGDMA_1000 : float, optional
        Amount of PEGDMA_1000 input per batch, in kg. The default is 4.5e-3.
    BIS : float, optional
        Amount of BIS input per batch, in kg. The default is 2.25e-4.
    TEMED : float, optional
        Amount of TEMED input per batch, in kg. The default is 6e-5.
    APS : float, optional
        Amount of APS input per batch, in kg. The default is 3e-5.
    PAC : float, optional
        Amount of PAC input per batch, in kg. The default is 3e-4.
    '''
    d_bead *= auom('inch').conversion_factor('m')
    v_bead = 4/3 * pi * (d_bead/2)**3 /2
    n_batch = V_beads / (n_bead*v_bead)
    EG, MAA = make_pegdma(n_batch*PEGDMA_1000)
    PAM, FMD, CuSO4, H2SO4 = make_bis(n_batch*BIS)
    DMA, EDCl, KOH = make_temed(n_batch*TEMED)
    NaPS = sub_APS(n_batch*APS)
    GAC = n_batch*PAC
    return dict(zip(encap_items.keys(), [EG, MAA, PAM, FMD, CuSO4, H2SO4, DMA, EDCl, KOH, NaPS, GAC]))



