# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import os

folder = os.path.dirname(__file__)
data_path = os.path.join(folder, 'data')
results_path = os.path.join(folder, 'results')
figures_path = os.path.join(folder, 'figures')
# To save simulation results and generated figures
if not os.path.isdir(results_path): os.mkdir(results_path)
if not os.path.isdir(figures_path): os.mkdir(figures_path)

from . import blower
from .blower import *

from . import units
from .units import *

__all__ = (
    'folder',
    'data_path',
    'results_path',
    'figures_path',
    'create_hap_cmps',
    *blower.__all__,
    *units.__all__,
	)

#%%
import qsdsan as qs
from qsdsan import Component, Components
from biorefineries.cane import create_sugarcane_chemicals

def create_hap_cmps(set_thermo=True, industrial_yeast_production=True):
    cmps_df = Components.load_default()
    chems = create_sugarcane_chemicals()
    Yeast = Component.from_chemical('Yeast', chemical=chems.Yeast,
                                    particle_size='Particulate', organic=True,
                                    degradability='Slowly', f_Vmass_Totmass=0.872)
    NH3 = Component('NH3', particle_size='Dissolved_gas', 
                    degradability='Undegradable', organic=False)
    CO2 = Component('CO2', particle_size='Dissolved_gas', 
                    degradability='Undegradable', organic=False)
    HAP = Component('HAP', search_ID='hydroxyapatite', 
                    particle_size='particulate', 
                    degradability='Undegradable', organic=False)
    ig_kwargs = dict(particle_size='Soluble', degradability='Undegradable', organic=False)
    org_kwargs = dict(particle_size='Soluble', degradability='Readily', organic=True)
    CaCl2 = Component('CaCl2', **ig_kwargs)
    
    # common urine consitituents
    urea = Component('Urea', **ig_kwargs)
    creatinine = Component('Creatinine', **org_kwargs)
    Hhip = Component('Hippuric_acid', search_ID='Hippuric acid', **org_kwargs)
    Hcit = Component('Citric_acid', search_ID='Citric acid', **org_kwargs)
    Hglu = Component.from_chemical('Glucuronic_acid', chemical='C6H10O7', **org_kwargs)
    Huric = Component('Uric_acid', search_ID='Uric acid', **org_kwargs)
    other_COD = cmps_df.S_F.copy('Other_COD')
    other_COD.i_N = other_COD.i_P = 0
    
    chloride = Component('Cl', search_ID='Cl-', **ig_kwargs)
    sodium = Component('Na', search_ID='Na+', **ig_kwargs)
    potassium = Component('K', search_ID='K+', **ig_kwargs)
    IS = Component('IS', search_ID='SO4-2', measured_as='S', **ig_kwargs)
    IP = Component('IP', search_ID='PO4-3', measured_as='P', **ig_kwargs)
    
    ash = Component.from_chemical('Ash', chemical=chems.Ash, **ig_kwargs)
    other_SS = ash.copy('other_SS')
    
    N2 = cmps_df.S_N2.copy('N2')
    O2 = cmps_df.S_O2.copy('O2')
    
    if industrial_yeast_production:
        cmps_boulardii = create_boulardii_cmps(False)
    else: cmps_boulardii = ()
    cmps = Components([cmps_df.H2O, N2, O2, Yeast, NH3, CO2, HAP, CaCl2, 
                       urea, creatinine, Hhip, Hcit, Hglu, Huric, other_COD,
                       chloride, sodium, potassium, IS, IP, ash, other_SS, 
                       *cmps_boulardii])
    cmps.default_compile()
    if set_thermo: qs.set_thermo(cmps)
    return cmps

def create_boulardii_cmps(default_compile=False):
    org_kwargs = dict(particle_size='Soluble', degradability='Readily', organic=True)
    ig_kwargs = dict(particle_size='Soluble', degradability='Undegradable', organic=False)
    glucose = Component('Glucose', **org_kwargs)
    molasses = Component('Molasses', search_ID='glucose', **org_kwargs)
    ammoium_sulfate = Component('Ammonium_sulfate', search_ID='(NH4)2SO4', **ig_kwargs)
    H3PO4 = Component('H3PO4', **ig_kwargs)
    MgSO4 = Component('MgSO4', **ig_kwargs)
    NaCl = Component('NaCl', **ig_kwargs)
    vb1 = Component('VB1', search_ID='Thiamine', **org_kwargs)
    vb2 = Component('VB2', search_ID='Riboflavin', **org_kwargs)
    vb5 = Component('VB5', search_ID='C18H32CaN2O10', **org_kwargs)
    vb6 = Component('VB6', search_ID='Pyridoxine', **org_kwargs)
    vb7 = Component('VB7', search_ID='Biotin', **org_kwargs)
    ethanol = Component('Ethanol', **org_kwargs)
    cmps = Components([glucose, molasses, ammoium_sulfate, H3PO4, MgSO4, NaCl,
                       vb1, vb2, vb5, vb6, vb7, ethanol])
    if default_compile:
        cmps_df = Components.load_default()
        chems = create_sugarcane_chemicals()
        Yeast = Component.from_chemical('Yeast', chemical=chems.Yeast,
                                        particle_size='Particulate', organic=True,
                                        degradability='Slowly')
        cmps = Components([*cmps, cmps_df.H2O, Yeast])
        cmps.default_compile()
    return cmps
    