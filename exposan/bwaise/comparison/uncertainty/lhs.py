#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 14:17:57 2017

@author: John Trimmer
"""

import numpy as np
from pyDOE import *
from scipy import stats

# LHS function to generate uniform distribution
def lhs_uniform(xmin, xmax, nsample):
    lhd = lhs(1, samples = nsample) #Generates n samples for each variable
    start = xmin  # defines the starting location of the distribution, at xmin
    width = xmax - xmin  # defines the total width of the distribution, from xmin to xmax
    lhd = stats.uniform.ppf(lhd, loc = start, scale = width)  # scales lhd to correspond to starting location and width
    return lhd

# LHS function to generate triangular distribution
def lhs_triangle(xmin, probable, xmax, nsample):
    lhd = lhs(1, samples = nsample) #Generates n samples for each variable
    start = xmin  # defines the starting location of the distribution, at xmin
    width = xmax - xmin  # defines the total width of the distribution, from xmin to xmax
    mode_loc = (probable - xmin) / (xmax - xmin)  # defines the location of the most probable value
    lhd = stats.triang.ppf(lhd, mode_loc, loc = start, scale = width)  # scales lhd to generate triangular distribution
    return lhd

# LHS function to generate normal distribution
#   optionally, xmin and xmax can be specified to truncate the distribution;
#   if they are not specified, the distribution will not be truncated
def lhs_normal(mean, stdev, nsample, xmin = -np.inf, xmax = np.inf):
    lhd = lhs(1, samples = nsample) #Generates n samples for each variable
    if np.isnan(xmin):
        xmin = -np.inf  # if xmin is not a number, set equal to default value (-inf)
    if np.isnan(xmax):
        xmax = np.inf  # if xmax is not a number, set equal to default value (inf)
    a = (xmin - mean) / stdev  # defines the location of the minimum value relative to standard normal distribution
    b = (xmax - mean) / stdev  # defines the location of the maximum value relative to standard normal distribution
    lhd = stats.truncnorm.ppf(lhd, a, b, loc = mean, scale = stdev)
    return lhd

# for an input series that contains the desired distribution, identify and perform
# the correct LHS function
def lhs_distribution(input_series, correlation_distributions, correlation_parameters, nsample):
    distribution = input_series.distribution
    if distribution == 'constant':
        lhd = np.full((nsample, 1), input_series.expected)
    elif distribution == 'uniform':
        lhd = lhs_uniform(input_series.low, input_series.high, nsample)
    elif distribution == 'triangular':
        lhd = lhs_triangle(input_series.low, input_series.expected, input_series.high, nsample)
    elif distribution == 'normal':
        lhd = lhs_normal(input_series.expected, input_series.standard_deviation, nsample,
                         xmin = input_series.low, xmax = input_series.high)

    # if correlations will be assessed for the parameter, add it to the matrix
    correlation = input_series.correlation
    if correlation == 'yes':
        # identify the first unused column
        i = 0
        while not(np.isnan(correlation_distributions[0,i])):
            i = i + 1
        # add distribution to that column
        correlation_distributions[:,i] = lhd[:,0]
        correlation_parameters[i] = input_series.name

    return lhd, correlation_distributions, correlation_parameters