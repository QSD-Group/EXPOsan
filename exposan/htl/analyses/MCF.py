#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 14:27:26 2023

@author: jiananfeng
"""

import pandas as pd
from scipy import stats

data = pd.read_csv('XXX.csv')

for i in ['heat','CHGCAT','CHGYIELD','HTCAT','HTYIELD','HCCAT','HCYIELD','HTLCAP','CHGCAP', 
          'HTCAP','OIL','AQUEOUS','BIOCHAR','GAS']:
    if stats.kstest(data[i].dropna(), data[f'{i}_H'].dropna()).pvalue < 0.05:
        print(i)
        print(stats.kstest(data[i].dropna(), data[f'{i}_H'].dropna()))