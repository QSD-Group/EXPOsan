#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 09:42:58 2024

@author: jiananfeng
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.linalg import cholesky

# generate samples
np.random.seed(3221)
a_samples_known = np.random.triangular(0, 2, 10, 1000)
b_samples_known = np.random.exponential(2, 1000)

# step 1: apply Copula (introduce quantile dependency)
normal_samples = np.random.randn(len(a_samples_known), 2)

# define the correlation matrix
rho = 0
cov_matrix = np.array([[1, rho], [rho, 1]])

# apply Cholesky decomposition to introduce dependency
L = cholesky(cov_matrix)
dependent_samples = normal_samples @ L.T

# map the dependent standard normal samples to uniform [0, 1]
dependent_uniform = norm.cdf(dependent_samples)
a_dependent_uniform = dependent_uniform[:, 0]
b_dependent_uniform = dependent_uniform[:, 1]

# step 2: map the dependent uniform samples back to the original distributions
a_new_samples = np.quantile(a_samples_known, a_dependent_uniform)
b_new_samples = np.quantile(b_samples_known, b_dependent_uniform)

sum_results = a_new_samples + b_new_samples

# a
plt.figure(figsize=(8, 6))
plt.hist(a_samples_known, bins=30, color='skyblue', edgecolor='black', alpha=0.7, density=True)
plt.title('Distribution of a_samples_known (Normal Distribution)', fontsize=14)
plt.xlabel('Value', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.grid(alpha=0.3)
plt.show()

# b
plt.figure(figsize=(8, 6))
plt.hist(b_samples_known, bins=30, color='skyblue', edgecolor='black', alpha=0.7, density=True)
plt.title('Distribution of a_samples_known (Normal Distribution)', fontsize=14)
plt.xlabel('Value', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.grid(alpha=0.3)
plt.show()

# a & b correlation
plt.figure(figsize=(10, 6))
plt.scatter(a_new_samples, b_new_samples, alpha=0.6, color='blue', label='Dependent Samples')
plt.title('Scatter plot of dependent a and b (using known samples)')
plt.xlabel('a_new_samples')
plt.ylabel('b_new_samples')
plt.grid(True)
plt.legend()
plt.show()

# a+b
plt.hist(sum_results, bins=30, color='skyblue', edgecolor='black', alpha=0.7, density=True)
plt.title('Distribution of a_samples_known (Normal Distribution)', fontsize=14)
plt.xlabel('Value', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.grid(alpha=0.3)
plt.show()