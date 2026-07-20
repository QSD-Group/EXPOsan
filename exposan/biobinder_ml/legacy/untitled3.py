# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 13:20:06 2026

@author: aliah
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

meta_path= r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\META_SHAP3_manure_CHCU_No_EC_prefer-max_top_ratio_20260128_1607.xlsx"
df = pd.read_excel(meta_path)
d = df[df["OK"]].copy()

d["IRR_pct"].describe()
d.nlargest(10, "IRR_pct")[[
    "IRR_pct",
    "product_ratio_biobinder_over_biofuel",
    "Biofuel_price",
    "Y_biocrude",
    "Y_aqueous",
    "Y_gas",
    "Y_char",
]]
