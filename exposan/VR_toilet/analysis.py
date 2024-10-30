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
from exposan.VR_toilet.models import (INCLUDE_RESOURCE_RECOVERY,
                                      create_model,
                                      run_uncertainty,
                                      )

        
def run(model_IDs, seed=None, N=1000, country_specific=False, **model_kwargs):
    # Make it possible to run with one or more models
    if isinstance(model_IDs, str): model_IDs = (model_IDs, )
    for ID in model_IDs:
        model = create_model(ID, country_specific=country_specific, **model_kwargs)
        run_uncertainty(model, seed=seed, N=N)

if __name__ == '__main__': #checks whether the script is being run directly (rather than imported as a module in another script). 
#This is often used to separate code that should only execute when the file is run as a standalone program.
    INCLUDE_RESOURCE_RECOVERY = True
    run(('A'), seed=5, N=500) # running systems A and B for contextual analysis with DMsan
        
