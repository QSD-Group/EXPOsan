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

from . import units
from .units import *

__all__ = (
    'folder',
    'data_path',
    'results_path',
    'figures_path',
    'create_hap_cmps',
    *units.__all__,
	)

#%%
import qsdsan as qs
from qsdsan import Component, Components
from biorefineries.cane import create_sugarcane_chemicals

def create_hap_cmps(set_thermo=True):
    cmps_df = Components.load_default()
    chems = create_sugarcane_chemicals()
    Yeast = Component.from_chemical('Yeast', chemical=chems.Yeast,
                                    particle_size='Particulate', organic=True,
                                    degradability='Slowly')
    NH3 = Component('NH3', particle_size='Dissolved_gas', 
                    degradability='Undegradable', organic=False)
    CO2 = Component('CO2', particle_size='Dissolved_gas', 
                    degradability='Undegradable', organic=False)
    HAP = Component('HAP', search_ID='hydroxyapatite', 
                    particle_size='particulate', 
                    degradability='Undegradable', organic=False)
    CaCl2 = Component('CaCl2', particle_size='Soluble', 
                    degradability='Undegradable', organic=False)
    
    # common urine consitituents
    org_kwargs = dict(particle_size='Soluble', degradability='Readily', organic=True)
    urea = Component('Urea', **org_kwargs)
    creatinine = Component('Creatinine', **org_kwargs)
    Hhip = Component('Hippuric_acid', search_ID='Hippuric acid', **org_kwargs)
    Hcit = Component('Citric_acid', search_ID='Citric acid', **org_kwargs)
    Hglu = Component.from_chemical('Glucuronic_acid', chemical='C6H10O7', **org_kwargs)
    Huric = Component('Uric_acid', search_ID='Uric acid', **org_kwargs)
    other_COD = cmps_df.S_F.copy('Other_COD')
    
    ig_kwargs = dict(particle_size='Soluble', degradability='Undegradable', organic=False)
    chloride = Component('Cl', search_ID='Cl-', **ig_kwargs)
    sodium = Component('Na', search_ID='Na+', **ig_kwargs)
    potassium = Component('K', search_ID='K+', **ig_kwargs)
    IS = Component('IS', search_ID='SO4-2', measured_as='S', **ig_kwargs)
    IP = Component('IP', search_ID='PO4-3', measured_as='P', **ig_kwargs)
    
    ash = Component.from_chemical('Ash', chemical=chems.Ash, **ig_kwargs)
    other_SS = ash.copy('other_SS')
    
    cmps = Components([cmps_df.H2O, Yeast, NH3, CO2, HAP, CaCl2, 
                       urea, creatinine, Hhip, Hcit, Hglu, Huric, other_COD,
                       chloride, sodium, potassium, IS, IP, ash, other_SS])
    cmps.default_compile()
    if set_thermo: qs.set_thermo(cmps)
    return cmps
