#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

import pandas as pd
from scipy import stats

data = pd.read_csv('XXX.csv')

for i in ['heat','CHGCAT','CHGYIELD','HTCAT','HTYIELD','HCCAT','HCYIELD','HTLCAP','CHGCAP', 
          'HTCAP','OIL','AQUEOUS','BIOCHAR','GAS']:
    if stats.kstest(data[i].dropna(), data[f'{i}_H'].dropna()).pvalue < 0.05:
        print(i)
        print(stats.kstest(data[i].dropna(), data[f'{i}_H'].dropna()))