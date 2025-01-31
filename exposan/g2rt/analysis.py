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

__all__ = ('run_lifetime',
           'run_ppl',)

def run(model_IDs, seed=None, N=1000, country_specific=False, **model_kwargs):
    # Make it possible to run with one or more models
    if isinstance(model_IDs, str): model_IDs = (model_IDs, )
    for ID in model_IDs:
        model = create_model(ID, country_specific=country_specific, **model_kwargs)
        g2rt.INCLUDE_RESOURCE_RECOVERY = True
        run_uncertainty(model, seed=seed, N=N)

if __name__ == '__main__': #checks whether the script is being run directly (rather than imported as a module in another script). 
#This is often used to separate code that should only execute when the file is run as a standalone program.
    run(('A'), N=50) # running systems A and B for contextual analysis with DMsan

#%%
def run_lifetime(model_IDs, lifetimes, seed=None, N=10000, country_specific=False, **model_kwargs):
    '''
    Run uncertainty analysis for models with different lifetimes.

    Parameters:
    ----------
    model_IDs : str or list of str
        One or more model identifiers.
    lifetimes : list of int or float
        Lifetimes to test.
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
    if not model_IDs or not lifetimes:
        raise ValueError("Both `model_IDs` and `lifetimes` must be non-empty.")
        
    g2rt.INCLUDE_RESOURCE_RECOVERY = False
    for ID in model_IDs:
        for lifetime in lifetimes:
            model = create_model(ID, country_specific=country_specific,lifetime= lifetime, **model_kwargs)
            run_uncertainty(model, note=f"lifetime_{lifetime}_N_{N}" ,seed=seed, N=N)

#%%
def run_ppl(model_IDs, ppls, seed=None, N=10000, country_specific=False, **model_kwargs):
    '''
    Run uncertainty analysis for models with different user numbers.

    Parameters:
    ----------
    model_IDs : str or list of str
        One or more model identifiers.
    ppls : list of int or float
        user numbers to test.
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
    if not model_IDs or not ppls:
        raise ValueError("Both `model_IDs` and `lifetimes` must be non-empty.")
        
    g2rt.INCLUDE_RESOURCE_RECOVERY = False
    for ID in model_IDs:
        for ppl in ppls:
            model = create_model(ID, country_specific=country_specific,ppl=ppl, **model_kwargs)
            run_uncertainty(model, note=f"ppl_{ppl}_N_{N}" ,seed=seed, N=N)