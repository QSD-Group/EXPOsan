#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 08:16:41 2023

@author: saumitrarai
"""
# libraries
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline  
# set width of bars
barWidth = 0.25

# set heights of bars
bars1 = [93.73390558, 86.66666667, 97.5, 84.5]
bars2 = [85,	68.79271071, 54.05405405,	83.81877023]
 
# Set position of bar on X axis
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]

 
# Make the plot 
plt.bar(r1, bars1, color='#ED586F', width=barWidth, edgecolor='white', label='Columbia model')
plt.bar(r2, bars2, color='#79bf82', width=barWidth, edgecolor='white', label='QSDsan model')
 
# Add xticks on the middle of the group bars
plt.xlabel('Water quality parameter', fontweight='bold')
plt.ylabel('Percentage removal (%)', fontweight='bold')
plt.xticks([r + barWidth for r in range(len(bars1))], ['COD','TN','TP','TSS'])
 
# Create legend & Show graphic
# plt.legend()
plt.show()

