#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

# Run uncertainty analysis and Spearman without country-specific settings
import os, numpy as np
from exposan.g2rt.models import (create_model,
                                      run_uncertainty,
                                      )
from exposan import g2rt as g2rt
# from g2rt import (default_ppl,
#                   default_lifetime)

__all__ = ('run_ppl_lifetime',
           # 'run_ppl',
           'run')

def run(model_IDs, seed=None, N=10000, city_specific=False, note=None, **model_kwargs):
    # Make it possible to run with one or more models
    if isinstance(model_IDs, str): model_IDs = (model_IDs, )
    for ID in model_IDs:
        model = create_model(ID, city_specific=city_specific, **model_kwargs)
        g2rt.INCLUDE_RESOURCE_RECOVERY = False
        run_uncertainty(model, note=note, seed=seed, N=N)

if __name__ == '__main__': #checks whether the script is being run directly (rather than imported as a module in another script). 
#This is often used to separate code that should only execute when the file is run as a standalone program.
    run(('A'), N=50) # running systems A and B for contextual analysis with DMsan

#%%
def run_N (model_IDs, Ns=10000, seed=None, city_specific=False, note=None, **model_kwargs):
    for ID in model_IDs:
        for N in Ns:
            model = create_model(ID, city_specific=city_specific, **model_kwargs)
            g2rt.INCLUDE_RESOURCE_RECOVERY = False
            run_uncertainty(model, note=f"{N}", seed=seed, N=N)
#%%
def run_resource_recovery(model_IDs, seed=None, N=10000, city_specific=False, **model_kwargs):
    '''
    Run uncertainty analysis for resource recovery True and False.

    Parameters:
    ----------
    model_IDs : str or list of str
        One or more model identifiers.
    seed : int, optional
        Random seed for reproducibility.
    N : int, default=10000
        Number of uncertainty samples.
    city_specific : bool, default=False
        Whether to use city-specific data.
    model_kwargs : dict
        Additional arguments to pass to `create_model`.
    
    Returns:
    --------
    None

    '''
    if isinstance(model_IDs, (str, np.str_)):  # Handle both Python and NumPy strings
        model_IDs = [model_IDs]

    g2rt.INCLUDE_RESOURCE_RECOVERY = False
    for ID in model_IDs:
        model = create_model(ID, city_specific=city_specific,**model_kwargs)
        if g2rt.INCLUDE_RESOURCE_RECOVERY == False:
            run_uncertainty(model, note=f"no_RR_default_{N}" ,seed=seed, N=N)
        else:run_uncertainty(model, note=f"RR_default_{N}" ,seed=seed, N=N)
    g2rt.INCLUDE_RESOURCE_RECOVERY = True
    for ID in model_IDs:
        model = create_model(ID, city_specific=city_specific,**model_kwargs)
        if g2rt.INCLUDE_RESOURCE_RECOVERY == False:
            run_uncertainty(model, note=f"no_RR_default_{N}" ,seed=seed, N=N)
        else:run_uncertainty(model, note=f"RR_default_{N}" ,seed=seed, N=N)
    
#%%
def run_ppl_lifetime(model_IDs, lifetimes=10, ppls=6, seed=None, N=10000, city_specific=False, **model_kwargs):
    '''
    Run uncertainty analysis for models with different lifetimes.

    Parameters:
    ----------
    model_IDs : str or list of str
        One or more model identifiers.
    lifetimes : list of int or float
        Lifetimes to test.
    ppls : list of int or float
        user numbers to test.
    seed : int, optional
        Random seed for reproducibility.
    N : int, default=10000
        Number of uncertainty samples.
    city_specific : bool, default=False
        Whether to use city-specific data.
    model_kwargs : dict
        Additional arguments to pass to `create_model`.
    
    Returns:
    --------
    None

    '''
    # Validate inputs
    if isinstance(model_IDs, (str, np.str_)):  # Handle both Python and NumPy strings
        model_IDs = [model_IDs]
    if not model_IDs:
        raise ValueError("`model_IDs` must be non-empty.")
        
    g2rt.INCLUDE_RESOURCE_RECOVERY = False
    for ID in model_IDs:
        for lifetime in lifetimes:
            for ppl in ppls:
                model = create_model(ID, city_specific=city_specific,lifetime= lifetime, ppl=ppl, **model_kwargs)
                run_uncertainty(model, note=f"ppl_{ppl}_lifetime_{lifetime}_N_{N}" ,seed=seed, N=N)
#%% 
def run_learning_curve (model_IDs, seed=None, N=10000, city_specific=False, **model_kwargs):
    '''
    Run uncertainty analysis for costs with pessimistic and optimistic learning curve assumptions.
    # Learning curve assumptions
    percent_CAPEX_to_scale = 0.65
    number_of_units = 100000

    if OPTIMISTIC_LEARNING_CURVE:
        percent_limit = 0.01  # Optimistic learning curve
        learning_curve_percent = 0.9  # Optimistic learning curve
    else:
        percent_limit = 0.03  # Pessimistic learning curve
        learning_curve_percent = 0.95  # Pessimistic learning curve
    Parameters:
    ----------
    model_IDs : str or list of str
        One or more model identifiers.
    seed : int, optional
        Random seed for reproducibility.
    N : int, default=10000
        Number of uncertainty samples.
    city_specific : bool, default=False
        Whether to use city-specific data.
    model_kwargs : dict
        Additional arguments to pass to `create_model`.
    
    Returns:
    --------
    None

    '''
    if isinstance(model_IDs, (str, np.str_)):  # Handle both Python and NumPy strings
        model_IDs = [model_IDs]

    g2rt.OPTIMISTIC_LEARNING_CURVE = False
    for ID in model_IDs:
        model = create_model(ID, city_specific=city_specific,**model_kwargs)
        if g2rt.OPTIMISTIC_LEARNING_CURVE == False:
            run_uncertainty(model, note=f"pessimistic_default_{N}" ,seed=seed, N=N)
        else:run_uncertainty(model, note=f"optimistic_default_{N}" ,seed=seed, N=N)
    g2rt.OPTIMISTIC_LEARNING_CURVE = True
    for ID in model_IDs:
        model = create_model(ID, city_specific=city_specific,**model_kwargs)
        if g2rt.OPTIMISTIC_LEARNING_CURVE == False:
            run_uncertainty(model, note=f"pessimistic_default_{N}" ,seed=seed, N=N)
        else:run_uncertainty(model, note=f"optimistic_default_{N}" ,seed=seed, N=N)

#%%
def run_mSCWO_replacement_cost(costs, model_ID='B', seed=None, N=10000, country_specific=False, **model_kwargs):
    '''
    Run uncertainty analysis for mSCWO with different material replacement costs.

    Parameters:
    ----------
    model_IDs : str or list of str
        One or more model identifiers.
    costs : list of float
        fraction of CAPEX for annual material replacement cost, 0-1
    seed : int, optional
        Random seed for reproducibility.
    N : int, default=10000
        Number of uncertainty samples.
    country_specific : bool, default=False
        Whether to use country-specific data.
    model_kwargs : dict
        Additional arguments to pass to `create_model`.
    
    Returns:
    --------
    None

    '''
    # Validate inputs
    if not costs:
        raise ValueError("`costs` must be non-empty.")
    g2rt.INCLUDE_RESOURCE_RECOVERY = False
    for cost in costs:
        model = create_model(model_ID, country_specific=country_specific, 
                             mscwo_replacement_cost=cost, **model_kwargs)
        run_uncertainty(model, note=f"mscwo_replace_cost{cost}_{N}" ,seed=seed, N=N)

#%%
def run_mSCWO_replacement_capital_cost(replacement_costs,system_costs, model_ID='B', seed=None, N=10000, city_specific=False, **model_kwargs):
    '''
    Run uncertainty analysis for mSCWO with different material replacement costs and system costs.

    Parameters:
    ----------
    model_IDs : str or list of str
        One or more model identifiers.
    costs : list of float
        fraction of CAPEX for annual material replacement cost, 0-1
    seed : int, optional
        Random seed for reproducibility.
    N : int, default=10000
        Number of uncertainty samples.
    country_specific : bool, default=False
        Whether to use country-specific data.
    model_kwargs : dict
        Additional arguments to pass to `create_model`.
    
    Returns:
    --------
    None

    '''
    # Validate inputs
    if not replacement_costs or not system_costs:
        raise ValueError("`replacement_costs` and `system_costs` must be non-empty.")
    g2rt.INCLUDE_RESOURCE_RECOVERY = False
    for replacement_cost in replacement_costs:
        for system_cost in system_costs:
            model = create_model(model_ID, city_specific=city_specific, 
                                 mscwo_replacement_cost=replacement_cost, 
                                 mscwo_equipment_cost=system_cost,
                                 **model_kwargs)
            run_uncertainty(model, note=f"mscwo_replace_{replacement_cost}_system_{system_cost}_{N}" ,seed=seed, N=N)
#%%
def run_e_CF(model_IDs, e_CFs, seed=None, N=10000, country_specific=False, **model_kwargs):
    '''
    Run uncertainty analysis for models with different user numbers.

    Parameters:
    ----------
    model_IDs : str or list of str
        One or more model identifiers.
    e_CFs : list of float
        The CO2 characterization factors of electricity in unit of kg CO2eq·kWh-1
        range: 0-2
    seed : int, optional
        Random seed for reproducibility.
    N : int, default=10000
        Number of uncertainty samples.
    country_specific : bool, default=False
        Whether to use country-specific data.
    model_kwargs : dict
        Additional arguments to pass to `create_model`.
    
    Returns:
    --------
    None

    '''
    # Validate inputs
    if isinstance(model_IDs, (str, np.str_)):  # Handle both Python and NumPy strings
        model_IDs = [model_IDs]
    if not model_IDs or not e_CFs:
        raise ValueError("Both `model_IDs` and `e_CFs` must be non-empty.")
        
    g2rt.INCLUDE_RESOURCE_RECOVERY = False
    for ID in model_IDs:
        for e_CF in e_CFs:
            model = create_model(ID, country_specific=country_specific,e_CF=e_CF, **model_kwargs)
            run_uncertainty(model, note=f"e_CF_{e_CF}_{N}" ,seed=seed, N=N)
#%%
def run_flushwater(model_IDs, flush_waters, seed=None, N=10000, city_specific=False, **model_kwargs):
    '''
    Run uncertainty analysis for models with different user numbers.

    Parameters:
    ----------
    model_IDs : str or list of str
        One or more model identifiers.
    flush_waters : list of float
        The flushing water in unit of kg/cap/hour
        range: 0-4.2
    seed : int, optional
        Random seed for reproducibility.
    N : int, default=10000
        Number of uncertainty samples.
    country_specific : bool, default=False
        Whether to use country-specific data.
    model_kwargs : dict
        Additional arguments to pass to `create_model`.
    
    Returns:
    --------
    None

    '''
    # Validate inputs
    if isinstance(model_IDs, (str, np.str_)):  # Handle both Python and NumPy strings
        model_IDs = [model_IDs]
    if not model_IDs or not flush_waters:
        raise ValueError("Both `model_IDs` and `flush_waters` must be non-empty.")
        
    g2rt.INCLUDE_RESOURCE_RECOVERY = False
    for ID in model_IDs:
        for flush_water in flush_waters:
            model = create_model(ID, city_specific=city_specific,flush_water=flush_water, **model_kwargs)
            run_uncertainty(model, note=f"flush_water_{flush_water}_{N}" ,seed=seed, N=N)
          
#%%
def run_vr_combustion_methane_EF(EFs, model_ID='A', seed=None, N=10000, country_specific=False, **model_kwargs):
    '''
    Run uncertainty analysis for mSCWO with different material replacement costs.

    Parameters:
    ----------
    model_ID : str or list of str
        A for volume reduction unit.
    EFs : list of float
        combustion CH4 emission factor, g-CH4/g-feces solids
    seed : int, optional
        Random seed for reproducibility.
    N : int, default=10000
        Number of uncertainty samples.
    country_specific : bool, default=False
        Whether to use country-specific data.
    model_kwargs : dict
        Additional arguments to pass to `create_model`.
    
    Returns:
    --------
    None

    '''
    # Validate inputs
    if not EFs:
        raise ValueError("`EFs` must be non-empty.")
    g2rt.INCLUDE_RESOURCE_RECOVERY = False
    for EF in EFs:
        model = create_model(model_ID, country_specific=country_specific, 
                             combustion_CH4_EF=EF, **model_kwargs)
        run_uncertainty(model, note=f"vr_combustion_CH4_EF{EF}_{N}" ,seed=seed, N=N)