#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 21:19:51 2023

@author: lane
"""
    
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

# Calculate dry basis biochar yield (% mass)
def get_yield(temp, AC): 
    f_AC_dec = AC/100 #converts % ash content of feedstock to decimal
    
    # predictive equation for % biochar dry basis (db) yield - derived via Excel solver
    db_yield = 100 * (1.18 * f_AC_dec ** 0.843 + (1 - f_AC_dec) * 2.106 * np.exp(-0.0066 * temp))
    return db_yield
   
# Calculate carbon sequestration potential of biochar produced (% mass)
def get_CS(temp, AC): 
    db_yield = get_yield(temp, AC)  
    
    # predictive equation for % biochar dry ash-free (daf) yield - Neves et al. 2011
    daf_yield = 100 * (0.106 + 2.43 * np.exp(-0.0066 * temp))
    
    # predictive equation for % biochar fixed carbon - derived via Excel solver
    FC_biochar = 87.786 * daf_yield ** -0.483

    # calculate biochar volatile matter and ash content using calculated values from above eqns
    AC_biochar = (db_yield - daf_yield) * 100 / db_yield
    VM_biochar = 100 - AC_biochar - FC_biochar
    
    # calculations for carbon sequestration (CS) potential [% mass C/mass biochar] 
    C_biochar = (0.474 * VM_biochar + 0.963 * FC_biochar + 0.067 * AC_biochar) / 100 # Klasson 2017
    Cafb = (0.474 * VM_biochar + 0.963 * FC_biochar + 0.067 * AC_biochar) / (100 - AC_biochar) # Klasson 2017
    C_feedstock = -0.50 * AC + 54.51 # Krueger et al. 2021
    R50 = 0.17 * Cafb + 0.00479 # Klasson 2017
    CS = db_yield * (C_biochar*100) * R50 / C_feedstock # Zhao et al. 2013
    return CS  

def get_mass_C(temp, AC):
    CS = get_CS(temp, AC)
    C_feedstock = -0.50 * AC + 54.51 # Krueger et al. 2021
    daily_in = 18*24 # Biogenic Refinery loading rate of 18 kg feedstock/h (dry basis) * 24 h/d
    mass_C = daily_in * (C_feedstock/100) * (CS/100)
    return mass_C
    
def heatmap(X, Y, Z, x_label="", y_label="", cbarlabel="", filename="", no_labels=False, **kwargs):
    """
    Create a heatmap from a numpy meshgrid, array-like, and two lists of labels.

    Parameters
    ----------
    X, Y
        A numpy meshgrid of x and y variables.
    Z
        An array or meshgrid of z-values to be plotted.
    x_label
        The label for the x-axis. Optional.
    y_label
        The label for the y-axis. Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    no_labels
        Removes labeling and returns only the colormap and colorbar to allow for formatting of the figure outside Python.
    **kwags
        Defults to pcolormesh. 
    """

    # Plot the heatmap
    fig , ax = plt.subplots(figsize=[6,5])
    im = ax.pcolormesh(X, Y, Z, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, location='right')
    
    if no_labels:
        # removes ticks, borders, labels
        ax.set_xticks([]) 
        ax.set_yticks([]) 
        ax.spines[:].set_visible(False)

    else:    
        # set colorbar label properties
        # cbar.ax.set_ylabel("\n"+cbarlabel, va="top", 
        #                    #fontsize=16
        #                    )
        cbar.minorticks_on()
        # cbar.fontsize = 16

        # set plot title
        ax.set_title(cbarlabel)

        # show ticks at intervals and label them
        ax.xaxis.set_major_locator(MultipleLocator(100))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        
        ax.tick_params(which='both', direction='inout')
        
        # label x and y axes
        ax.xaxis.set_label_text(x_label, 
                                #fontsize=16
                                )
        ax.yaxis.set_label_text(y_label, 
                                #fontsize=16
                                )
    
        ax.spines[:].set_lw(0.75)
        ax.spines.bottom.set_bounds(X.min(), X.max())
        ax.spines.top.set_bounds(X.min(), X.max())
        ax.spines.left.set_bounds(Y.min(), Y.max())
        ax.spines.right.set_bounds(Y.min(), Y.max())

    # fig.set_size_inches((15,13))
    plt.savefig('/Users/lanet/Desktop/GitHub/EXPOsan/exposan/biogenic_refinery/results/'+filename+'.pdf', format='pdf')
    return im, cbar

# create  x and y variables
x = np.arange(300, 901)
y = np.linspace(15.0, 75.0, num=600)
X, Y = np.meshgrid(x, y)

# generate data for yield plot
Z = get_yield(X, Y)

# plot yields heatmap
yieldplot = heatmap(X, Y, Z, 
                    x_label="Pyrolysis temperature [\N{DEGREE CELSIUS}]", 
                    y_label="Feedstock ash content [% db]", 
                    cmap='YlGnBu', 
                    cbarlabel="Biochar yield [% db]", 
                    filename="YieldPlot")

# generate data for carbon sequestration plot
Z = get_CS(X, Y)

# plot CS heatmap
CSplot = heatmap(X, Y, Z, 
                    x_label="Pyrolysis temperature [\N{DEGREE CELSIUS}]", 
                    y_label="Feedstock ash content [% db]", 
                    cmap='YlGnBu', 
                    cbarlabel="Carbon sequestration potential [%]",
                    filename="CSplot")

# # generate data for carbon mass plot
# Z = get_mass_C(X, Y)

# # plot mass C heatmap
# massCplot = heatmap(X, Y, Z, 
#                     x_label="Pyrolysis temperature [\N{DEGREE CELSIUS}]", 
#                     y_label="Feedstock ash content [% db]", 
#                     cmap='YlGnBu', 
#                     cbarlabel="Carbon sequestered [kg\N{DOT OPERATOR}d\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}]",
#                     filename="massCplot")
